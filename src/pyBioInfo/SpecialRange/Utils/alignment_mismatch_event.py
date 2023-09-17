import re
import pysam
from pyBioInfo.Utils import BlockTools
from pyBioInfo.SpecialRange import MismatchSite


class AlignmentMismatchEvent(MismatchSite):
    def __init__(self, chrom, start, end, strand, name, rbase, tbase, score, distance, alignment):
        super(AlignmentMismatchEvent, self).__init__(chrom=chrom, start=start,
                                                     end=end, strand=strand, name=name,
                                                     rbase=rbase, tbase=tbase)
        self.score = score
        self.distance = distance
        self.alignment = alignment

    def __str__(self):
        return "\t".join(map(str, [self.chrom, self.start, self.end, self.name,
                                   self.score, self.strand, self.distance, self.rbase,
                                   self.tbase, self.is_snp]))


class AlignmentMismatchEventFactory(object):

    @classmethod
    def _from_alignment_by_sequence_comparison(cls, alignment, fasta):
        raise NotImplementedError()

    @classmethod
    def _from_alignment_by_tag(cls, alignment):
        array = []  # mismatch event
        chrom = alignment.chrom
        name = alignment.name
        strand = alignment.strand
        segment = alignment.segment
        aligned_query_sequence, aligned_query_qualities = alignment.get_aligned_sequences()
        mdstr = segment.get_tag("MD")
        i = 0
        while len(mdstr) > 0:
            res = re.match("^[0-9]+", mdstr)  # Mapped
            if res is not None:
                x, y = res.span()
                assert x == 0
                num = int(mdstr[0:y])
                mdstr = mdstr[y:]
                i += num
                continue
            res = re.match("^\^[A-Z]+", mdstr)  # Deletion event
            if res is not None:
                x, y = res.span()
                assert x == 0
                i2 = i + y - 1
                mdstr = mdstr[y:]
                i = i2
                continue
            res = re.match("^[A-Z]", mdstr)  # Mismatch event
            if res is not None:
                x, y = res.span()
                assert x == 0 and y == 1
                ref = mdstr[0]
                array.append([i, ref, aligned_query_sequence[i],
                              aligned_query_qualities[i]])
                mdstr = mdstr[1:]
                i += 1
                continue
        blocks = alignment.blocks
        events = []
        length = len(alignment)
        for idx, rbase, tbase, score in array:
            # 用BlockTools算position是因为不需要考虑strand信息
            position = alignment.position(idx, strandness=False)
            distance = min(idx, length - idx - 1)  # 距离两端最近的距离
            obj = AlignmentMismatchEvent(chrom=chrom, start=position, end=position + 1, name=name,
                                         strand=strand, rbase=rbase, tbase=tbase, score=score,
                                         distance=distance, alignment=alignment)
            events.append(obj)
        return events

    @classmethod
    def from_alignment(cls, alignment, fasta=None):
        if fasta is None:
            return cls._from_alignment_by_tag(alignment)
        else:
            return cls._from_alignment_by_sequence_comparison(alignment, fasta)