import re
import pysam
from pyBioInfo.Range import GRange
from pyBioInfo.IO.File import BedFile


# from pyBioInfo.Utils import ShiftLoader


class MismatchSite(GRange):
    MAPPER = {"A": "T", "C": "G", "T": "A", "G": "C", "N": 'N'}
    
    def __init__(self, chrom, start, end=None, name=None, strand=None, rbase=None, tbase=None, count=None, background=None, ratio=None):
        super(MismatchSite, self).__init__(chrom=chrom, start=start, end=start + 1, name=name, strand=strand)
        self.rbase = rbase
        self.tbase = tbase
        self.count = count
        self.background = background
        self.ratio = ratio

    def format(self, fmt=None, **kwargs):
        if fmt is None:
            items = [self.chrom, self.start, self.end, self.strand, self.name, self.ref, self.alt]
            line = "\t".join(map(str, items))
            return line
        return super(MismatchSite, self).format(fmt, **kwargs)
    
    def reverse_complement(self):
        self.rbase = self.MAPPER[self.rbase]
        self.tbase = self.MAPPER[self.tbase]
        

class MismatchSiteFactory(object):
    MAPPER = {"A": "T", "C": "G", "T": "A", "G": "C", "N": 'N'}

    # REGULAR_PATTERN1 = "[0-9]+[A-Z]"
    # REGULAR_PATTERN2 = "\^[A-Z]+"
    # REGULAR_PATTERN3 = "[0-9]+\^"
    # REGULAR_PATTERN4 = "[0-9]+$"

    @classmethod
    def _from_alignment_by_tag(cls, alignment):
        chrom = alignment.chrom
        name = alignment.name
        strand = alignment.strand
        segment = alignment.segment
        mate = None
        if segment.is_read1:
            mate = 1
        elif segment.is_read2:
            mate = 2

        query_sequence = segment.query_sequence
        query_qualities = segment.query_qualities
        query_alignment_sequence = []
        query_alignment_qualities = []
        offset = 0
        for flag, count in segment.cigartuples:
            if flag == pysam.CHARD_CLIP:
                pass
            elif flag == pysam.CSOFT_CLIP:
                offset += count
            elif flag == pysam.CREF_SKIP:
                pass
            elif flag == pysam.CDEL:
                # print("D")
                # offset += count
                pass
            elif flag == pysam.CINS:
                offset += count
            elif flag == pysam.CMATCH:
                query_alignment_sequence.append(query_sequence[offset:offset + count])
                query_alignment_qualities.extend(query_qualities[offset:offset + count])
                offset += count
            else:
                raise Exception()

        qseq = "".join(query_alignment_sequence)
        qual = query_alignment_qualities
        value = segment.get_tag("MD")

        i = 0
        j = 0
        sites = []
        # print(value)
        while len(value) > 0:
            # print(value)
            res = re.match("^[0-9]+", value) # Mapped
            if res is not None:
                x, y = res.span()
                num = int(value[x:y])
                value = value[y:]
                i += num
                j += num
                continue
            res = re.match("^\^[A-Z]+", value) # Deletion
            if res is not None:
                x, y = res.span()
                value = value[y:]
                i += y - x - 1
                continue
            res = re.match("^[A-Z]", value) # Mismatch
            if res is not None:
                x, y = res.span()
                ref = value[x:y]
                # print(ref)
                sites.append([i, ref, qseq[j], qual[j]])
                value = value[y:]
                i += y
                j += y
                continue

        results = []
        with alignment.ignore_strand():
            for index, ref, alt, score in sites:
                start = alignment.position(index)
                obj = MismatchSite(chrom=chrom, start=start, end=start + 1, strand=strand, name=name)
                obj["rbase"] = ref
                obj["tbase"] = alt
                obj["index"] = index
                obj["width"] = len(qseq)
                obj["score"] = score
                obj["mate"] = mate
                obj["alignment"] = alignment
                results.append(obj)
        # print("A", len(sites))
        return results

    @classmethod
    def _from_alignment_by_sequence_comparison(cls, alignment, fasta):
        # Extract mismatch site by base comparison mode.
        assert False
        results = []
        segment = alignment.segment
        chrom = segment.reference_name
        start = segment.reference_start
        name = segment.query_name
        strand = alignment.strand
        mate = None
        if segment.is_read1:
            mate = 1
        elif segment.is_read2:
            mate = 2
        query_sequence = segment.query_sequence
        query_qualities = segment.query_qualities
        width = len(query_sequence)
        offset = 0
        start0 = start
        for cigar, count in segment.cigartuples:
            if cigar == pysam.CMATCH:
                rseq = fasta.fetch(chrom=chrom, start=start0, end=start0 + count).upper()
                for k in range(count):
                    index = offset + k
                    qbase = query_sequence[index]
                    rbase = rseq[k]
                    if qbase != rbase:
                        x = start0 + k
                        y = x + 1
                        quality = None
                        if query_qualities is not None:
                            quality = query_qualities[index]
                        obj = MismatchSite(chrom=chrom, start=x, end=y, strand=strand, name=name, )
                        obj["source"] = "Alignment"
                        obj["score"] = quality
                        obj["rbase"] = rbase
                        obj["tbase"] = rbase
                        obj["index"] = index
                        obj["width"] = width
                        obj["mate"] = mate
                        obj["alignment"] = alignment
                        results.append(obj)
                offset += count
                start0 += count
            elif cigar == pysam.CINS:
                offset += count
            elif cigar == pysam.CDEL:
                start0 += count
            elif cigar == pysam.CREF_SKIP:
                start0 += count
            elif cigar == pysam.CSOFT_CLIP:
                offset += count
            elif cigar == pysam.CHARD_CLIP:
                pass
            else:
                raise RuntimeError()
        return results

    @classmethod
    def _remove_blacklist_site(cls, sites, blacklist):
        results1 = []
        results2 = []
        i = 0
        for site in sites:
            remove = False
            j = 0
            while True:
                k = i + j
                if k >= len(blacklist):
                    break
                item = blacklist[k]
                if item.start >= site.end:
                    break
                if item.end <= site.start:
                    i += 1
                    continue
                assert item.start == site.start
                if site.rbase == item.rbase:
                    if site.tbase == item.tbase:
                        remove = True
                        break
                    else:
                        j += 1
                        continue
                elif site.rbase == cls.MAPPER[item.rbase]:
                    if site.tbase == cls.MAPPER[item.tbase]:
                        remove = True
                        break
                    else:
                        j += 1
                        continue
                else:
                    # print(site.chrom, site.start, site.name, site.rbase, site.tbase)
                    # print(item.chrom, item.start, item.name, item.rbase, item.tbase)
                    raise RuntimeError()
            if remove:
                results2.append(site)
            else:
                results1.append(site)
        return results1, results2

    @classmethod
    def from_alignment(cls, alignment, fasta=None, relative=None, blacklist=None):
        results1 = None  # Purge
        results2 = None  # Hit blacklist

        if fasta is None:
            results = cls._from_alignment_by_tag(alignment)
        else:
            results = cls._from_alignment_by_sequence_comparison(alignment, fasta)

        need_reverse = False
        if relative is not None:
            if relative:
                if alignment.reverse:
                    need_reverse = True
            else:
                if alignment.forward:
                    need_reverse = True

        if need_reverse:
            for site in results:
                rbase, tbase = site.rbase, site.tbase
                rbase, tbase = cls.MAPPER[rbase], cls.MAPPER[tbase]
                site["rbase"] = rbase
                site["tbase"] = tbase

        if blacklist is None:
            results1 = results
            results2 = []
        else:
            results1, results2 = cls._remove_blacklist_site(results, blacklist)

        return results1, results2

    @classmethod
    def _merge_sites(cls, sites):
        results = []
        last = None
        for site in sites:
            if last is None:
                last = site
            else:
                if site.start == last.start and site.tbase == last.tbase:
                    assert site.rbase == last.rbase
                    if site["score"] > last["score"]:
                        last = site
                else:
                    results.append(last)
                    last = site
        if last is not None:
            results.append(last)
        return results

    @classmethod
    def from_fragment(cls, fragment, fasta=None, relative=True, blacklist=None):
        mate1, mate2 = fragment.sub_genomic_ranges
        sites1, sites1_blacklist = cls.from_alignment(alignment=mate1, fasta=fasta, relative=None, blacklist=blacklist)
        sites2, sites2_blacklist = cls.from_alignment(alignment=mate2, fasta=fasta, relative=None, blacklist=blacklist)
        sites = list(sorted(sites1 + sites2))
        sites_blacklist = list(sorted(sites1_blacklist + sites2_blacklist))

        sites = cls._merge_sites(sites)
        sites_blacklist = cls._merge_sites(sites_blacklist)

        need_reverse = False
        if relative is not None:
            if relative:
                if fragment.reverse:
                    need_reverse = True
            else:
                if fragment.forward:
                    need_reverse = True

        if need_reverse:
            for site in sites:
                rbase, tbase = site.rbase, site.tbase
                rbase, tbase = cls.MAPPER[rbase], cls.MAPPER[tbase]
                site["rbase"] = rbase
                site["tbase"] = tbase
            for site in sites_blacklist:
                rbase, tbase = site.rbase, site.tbase
                rbase, tbase = cls.MAPPER[rbase], cls.MAPPER[tbase]
                site["rbase"] = rbase
                site["tbase"] = tbase

        for site in sites:
            site.strand = fragment.strand
        for site in sites_blacklist:
            site.strand = fragment.strand

        return sites, sites_blacklist

    @classmethod
    def from_bed_file(cls, path):
        with BedFile(path) as f:
            for record in f:
                rbase, tbase = record.name.split("-")
                obj = MismatchSite(chrom=record.chrom, start=record.start, strand=record.strand, name=record.name)
                obj.rbase = rbase
                obj.tbase = tbase
                yield obj
