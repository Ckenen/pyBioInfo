# Class: MismatchSite
from .special_range import SpecialRange


class MismatchSite(SpecialRange):
    BASE_MAPPER = {"A": "T",
                   "C": "G",
                   "T": "A",
                   "G": "C",
                   "N": 'N'}

    def __init__(self, chrom, start, end, name, strand, rbase, tbase, is_snp=False):
        if end:
            assert end - start == 1
        super(MismatchSite, self).__init__(chrom=chrom, start=start, end=start + 1,
                                           name=name, strand=strand, rbase=rbase,
                                           tbase=tbase, is_snp=is_snp)

    def reverse_complement(self):
        self.rbase = self.BASE_MAPPER[self.rbase]
        self.tbase = self.BASE_MAPPER[self.tbase]

    def format(self, fmt=None, **kwargs):
        if fmt is None:
            items = [self.chrom, self.start, self.end,
                     self.strand, self.name, self.rbase, self.tbase]
            line = "\t".join(map(str, items))
            return line
        return super(MismatchSite, self).format(fmt, **kwargs)
