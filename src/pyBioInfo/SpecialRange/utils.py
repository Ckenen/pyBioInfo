import re
import logging
import pysam
from pyBioInfo.IO.File import BedFile
from pyBioInfo.Utils import BlockTools
from .events import SoftClipEvent, InsertionEvent, DeletionEvent, MismatchEvent
from .sites import MismatchSite

class MismatchEventFactory(object):
    @classmethod
    def integrate_decision_maker_single(cls, event, mate):
        keep = True
        if mate.start <= event.start < mate.end:
            try:
                qualities = mate.get_aligned_qualities()
                score = qualities[mate.index(event.start, strandness=False)]
                if score >= event.score:
                    keep = False
            except ValueError:
                pass
        if keep:
            obj = event.copy()
            return obj
        return None

    @classmethod
    def integrate_decision_maker_double(cls, event1, event2):
        if event1.rbase == event2.rbase:
            if event1.tbase == event2.tbase:
                event = event1
                obj = event.copy()
                obj.score = max(event1.score, event2.score)
                obj.distance = max(event1.distance, event2.distance)
                return obj
            else:
                event = None
                if event1.score > event2.score:
                    event = event1
                elif event1.score < event2.score:
                    event = event2
                else:
                    if event1.distance > event2.distance:
                        event = event1
                    elif event1.distance < event2.distance:
                        event = event2
                    else:
                        return None
                        # if event1.score >= 30:
                        #     print(event1.chrom, event1.start, event1.name, event1.score, event1.rbase, event1.tbase)
                        #     print(event2.chrom, event2.start, event2.name, event2.score, event2.rbase, event2.tbase)
                        #     raise RuntimeError("Unsure.")
                        # else:
                        #     return None
                obj = event.copy()
                return obj
        else:
            raise RuntimeError("")
        return None

    @classmethod
    def integrate_fragment_mismatch_event(cls, fragment, array1, array2):
        i1 = 0
        i2 = 0
        results = []
        while True:
            if i1 >= len(array1):
                if i2 >= len(array2):
                    break
                else:
                    results.append(cls.integrate_decision_maker_single(array2[i2], fragment.mate1))
                    i2 += 1
            else:
                if i2 >= len(array2):
                    results.append(cls.integrate_decision_maker_single(array1[i1], fragment.mate2))
                    i1 += 1
                else:
                    event1, event2 = array1[i1], array2[i2]
                    if event1.start < event2.start:
                        results.append(cls.integrate_decision_maker_single(event1, fragment.mate2))
                        i1 += 1
                    elif event1.start > event2.start:
                        results.append(cls.integrate_decision_maker_single(event2, fragment.mate1))
                        i2 += 1
                    else:
                        results.append(cls.integrate_decision_maker_double(event1, event2))
                        i1 += 1
                        i2 += 1
        return list(filter(lambda item: item is not None, results))

    @classmethod
    def from_fragment(cls, fragment, fasta=None):
        array1 = cls.from_alignment(fragment.mate1, fasta)
        array2 = cls.from_alignment(fragment.mate2, fasta)
        results = cls.integrate_fragment_mismatch_event(fragment, array1, array2)
        return results

    @classmethod
    def mark_snp(cls, array, snps):
        i = 0
        for event in array:
            j = 0
            while True:
                k= i + j
                if k >= len(snps):
                    break
                else:
                    snp = snps[k]
                    if snp.start < event.start:
                        i += 1
                    elif snp.start == event.start:
                        if snp.rbase == event.rbase:
                            if snp.tbase == event.tbase:
                                event.is_snp = True
                                break
                            else:
                                j += 1
                        elif snp.rbase == MismatchEvent.BASE_MAPPER[event.rbase]:
                            if snp.tbase == MismatchEvent.BASE_MAPPER[event.tbase]:
                                event.is_snp = True
                                break
                            else:
                                j += 1
                        else:
                            raise RuntimeError("The reference base of SNP does not match the reference base of EVENT.")
                    else:
                        break
        return array

    @classmethod
    def from_alignment(cls, alignment, fasta=None):
        if fasta is None:
            return cls.from_alignment_by_tag(alignment)
        else:
            return cls.from_alignment_by_sequence_comparison(alignment, fasta)
    
    @classmethod
    def from_alignment_by_sequence_comparison(cls, alignment, fasta):
        raise NotImplementedError()
       

    @classmethod
    def from_alignment_by_tag(cls, alignment):
        mismatch_events = []

        chrom = alignment.chrom
        name = alignment.name
        strand = alignment.strand
        segment = alignment.segment

        # Extract query sequences and qualities
        query_sequence = segment.query_sequence
        query_qualities = segment.query_qualities
        qseq = []
        qual = []
        offset = 0
        for flag, count in segment.cigartuples:
            if flag == pysam.CHARD_CLIP:
                pass
            elif flag == pysam.CSOFT_CLIP: # Soft clip event
                offset2 = offset + count
                # soft_clip_events.append([offset, query_sequence[offset: offset2], query_qualities[offset: offset2]])
                offset = offset2
            elif flag == pysam.CREF_SKIP:
                pass
            elif flag == pysam.CDEL:
                qseq.append("-" * count)
                qual.extend([-1] * count) 
                pass
            elif flag == pysam.CINS:    # Insertion event
                offset2 = offset + count
                # insertion_events.append([offset, query_sequence[offset: offset2], query_qualities[offset: offset2]])
                offset = offset2
            elif flag == pysam.CMATCH:
                qseq.append(query_sequence[offset:offset + count])
                qual.extend(query_qualities[offset:offset + count])
                offset += count
            else:
                raise Exception()
        qseq = "".join(qseq)

        mdstr = segment.get_tag("MD")
        i = 0
        while len(mdstr) > 0:
            res = re.match("^[0-9]+", mdstr) # Mapped
            if res is not None:
                x, y = res.span()
                assert x == 0
                num = int(mdstr[0:y])
                mdstr = mdstr[y:]
                i += num
                continue
            res = re.match("^\^[A-Z]+", mdstr) # Deletion event
            if res is not None:
                x, y = res.span()
                assert x == 0
                i2 = i + y - 1
                # deletion_events.append([i, i2, mdstr[1:y]])
                mdstr = mdstr[y:]
                i = i2
                continue
            res = re.match("^[A-Z]", mdstr) # Mismatch event
            if res is not None:
                x, y = res.span()
                assert x == 0 and y == 1
                ref = mdstr[0]
                mismatch_events.append([i, ref, qseq[i], qual[i]])
                mdstr = mdstr[1:]
                i += 1
                continue
        
        blocks = alignment.blocks
        temp = []
        for idx, rbase, tbase, score in mismatch_events:
            position = BlockTools.get_position(blocks, idx, check=False, length=len(alignment))
            mate = 0
            if segment.is_read1:
                mate = 1
            elif segment.is_read2:
                mate = 2
            obj = MismatchEvent(chrom=chrom, 
                                name=name, 
                                start=position, 
                                strand=strand, 
                                mate=mate,
                                rbase=rbase, 
                                tbase=tbase, 
                                score=score, 
                                offset=idx, 
                                length=len(alignment))
            temp.append(obj)
        mismatch_events = temp

        return mismatch_events

class MismatchSiteFactory(object):
    @classmethod
    def from_bed_file(cls, path):
        with BedFile(path) as f:
            for record in f:
                rbase, tbase = record.name.split("-")
                obj = MismatchSite(chrom=record.chrom, start=record.start, strand=record.strand, name=record.name)
                obj.rbase = rbase
                obj.tbase = tbase
                yield obj


class EventFactory(object):
    @classmethod
    def from_alignment_by_tag(cls, alignment):
        soft_clip_events = []
        insertion_events = []
        deletion_events = []
        mismatch_events = []

        chrom = alignment.chrom
        name = alignment.name
        strand = alignment.strand
        segment = alignment.segment

        # Extract query sequences and qualities
        query_sequence = segment.query_sequence
        query_qualities = segment.query_qualities
        qseq = []
        qual = []
        offset = 0
        for flag, count in segment.cigartuples:
            if flag == pysam.CHARD_CLIP:
                pass
            elif flag == pysam.CSOFT_CLIP: # Soft clip event
                offset2 = offset + count
                soft_clip_events.append([offset, query_sequence[offset: offset2], query_qualities[offset: offset2]])
                offset = offset2
            elif flag == pysam.CREF_SKIP:
                pass
            elif flag == pysam.CDEL:
                qseq.append("-" * count)
                qual.extend([-1] * count) 
                pass
            elif flag == pysam.CINS:    # Insertion event
                offset2 = offset + count
                insertion_events.append([offset, query_sequence[offset: offset2], query_qualities[offset: offset2]])
                offset = offset2
            elif flag == pysam.CMATCH:
                qseq.append(query_sequence[offset:offset + count])
                qual.extend(query_qualities[offset:offset + count])
                offset += count
            else:
                raise Exception()
        qseq = "".join(qseq)

        mdstr = segment.get_tag("MD")
        i = 0
        while len(mdstr) > 0:
            res = re.match("^[0-9]+", mdstr) # Mapped
            if res is not None:
                x, y = res.span()
                assert x == 0
                num = int(mdstr[0:y])
                mdstr = mdstr[y:]
                i += num
                continue
            res = re.match("^\^[A-Z]+", mdstr) # Deletion event
            if res is not None:
                x, y = res.span()
                assert x == 0
                i2 = i + y - 1
                deletion_events.append([i, i2, mdstr[1:y]])
                mdstr = mdstr[y:]
                i = i2
                continue
            res = re.match("^[A-Z]", mdstr) # Mismatch event
            if res is not None:
                x, y = res.span()
                assert x == 0 and y == 1
                ref = mdstr[0]
                mismatch_events.append([i, ref, qseq[i], qual[i]])
                mdstr = mdstr[1:]
                i += 1
                continue
        
        blocks = alignment.blocks
        # soft_clip_events = []
        # insertion_events = []
        # deletion_events = []
        temp = []
        for idx, rbase, tbase, score in mismatch_events:
            position = BlockTools.get_position(blocks, idx, check=False, length=len(alignment))
            mate = 0
            if segment.is_read1:
                mate = 1
            elif segment.is_read2:
                mate = 2
            obj = MismatchEvent(chrom=chrom, 
                                name=name, 
                                start=position, 
                                strand=strand, 
                                mate=mate,
                                rbase=rbase, 
                                tbase=tbase, 
                                score=score, 
                                offset=idx, 
                                length=len(alignment))
            temp.append(obj)
        mismatch_events = temp

        return soft_clip_events, insertion_events, deletion_events, mismatch_events

    @classmethod
    def from_alignment_by_sequence_comparison(cls, alignment, fasta):
        raise NotImplementedError()

    @classmethod
    def from_alignment(cls, alignment, fasta=None):
        if fasta is None:
            return cls.from_alignment_by_tag(alignment)
        else:
            return cls.from_alignment_by_sequence_comparison(alignment, fasta)

    @classmethod
    def integrate_mismatch_event(cls, array):
        events = []
        for event in sorted(array):
            if len(events) > 0:
                last = events[-1]
                assert event.name == last.name
                if event.start == last.start:
                    assert event.rbase == last.rbase
                    if event.tbase == last.tbase:
                        if event.score > last.score:
                            last.score = event.score
                        elif event.distance > last.distance:
                            last.distance = event.distance
                    else:
                        if event.score > last.score:
                            last.tbase = event.tbase
                            last.score = event.tbase
                        elif event.distance > last.distance:
                            last.tbase = event.tbase
                            last.distance = event.distance
                        else:
                            obj = MismatchEvent(chrom=event.chrom, 
                                        start=event.start, 
                                        name=event.name, 
                                        rbase=event.rbase, 
                                        tbase=event.tbase, 
                                        strand=event.strand, 
                                        score=event.score, 
                                        distance=event.distance)
                            events.append(obj)
                elif event.start > last.start:
                    obj = MismatchEvent(chrom=event.chrom, 
                                        start=event.start, 
                                        name=event.name, 
                                        rbase=event.rbase, 
                                        tbase=event.tbase, 
                                        strand=event.strand, 
                                        score=event.score, 
                                        distance=event.distance)
                    events.append(obj)
                else:
                    raise RuntimeError("Unsorted array!")
            else:
                obj = MismatchEvent(chrom=event.chrom, 
                                    start=event.start, 
                                    name=event.name, 
                                    rbase=event.rbase, 
                                    tbase=event.tbase, 
                                    strand=event.strand, 
                                    score=event.score, 
                                    distance=event.distance)
                events.append(obj)
        return events

    @classmethod
    def from_fragment(cls, fragment, fasta=None):
        values1 = cls.from_alignment(fragment.mate1, fasta)
        values2 = cls.from_alignment(fragment.mate2, fasta)
        soft_clip_events = values1[0] + values2[0]
        insertion_events = values1[1] + values2[1]
        deletion_events = values1[2] + values2[2]
        mismatch_events = cls.integrate_mismatch_event(list(sorted(values1[3] + values2[3])))
        return soft_clip_events, insertion_events, deletion_events, mismatch_events