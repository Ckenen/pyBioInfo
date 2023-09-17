from pyBioInfo.Range import CRange

class SpecialRange(CRange):
    def __init__(self, chrom, start, end, strand, name, rbase, tbase, is_snp=False):
        super(SpecialRange, self).__init__(chrom=chrom, start=start, end=end, strand=strand, name=name)
        self.rbase = rbase
        self.tbase = tbase
        self.is_snp = is_snp
    
    def reverse_complement(self):
        raise NotImplementedError()
    
    def __str__(self):
        return "%s: %d-%d (%s) [%s -> %s] %s" % (self.chrom, self.start, self.end, self.strand, self.rbase, self.tbase, self.name)
    