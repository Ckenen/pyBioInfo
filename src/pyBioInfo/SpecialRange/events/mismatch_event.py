from pyBioInfo.Range import GRange


class MismatchEvent(GRange):
    BASE_MAPPER = {"A": "T", 
                   "C": "G", 
                   "T": "A", 
                   "G": "C", 
                   "N": 'N'}

    def __init__(self, chrom, start, strand=None, name=None, mate=None, rbase=None, tbase=None, score=None, offset=None, length=None, distance=None, is_snp=False):
        super(MismatchEvent, self).__init__(chrom=chrom, start=start, end=start + 1, name=name, strand=strand)
        self.mate = mate
        self.rbase = rbase
        self.tbase = tbase
        self.score = score
        self.offset = offset
        self.length = length
        self.distance = distance
        self.is_snp = is_snp
        if distance is None and self.offset is not None and self.length is not None:
            self.distance = min(self.offset, self.length - 1 - self.offset)

    def reverse_complement(self):
        self.rbase = self.BASE_MAPPER[self.rbase]
        self.tbase = self.BASE_MAPPER[self.tbase]

    def copy(self):
        obj = MismatchEvent(chrom=self.chrom, start=self.start, strand=self.strand, name=self.name,
                            mate=self.mate, rbase=self.rbase, tbase=self.tbase, score=self.score, 
                            offset=self.offset, length=self.length, distance=self.distance, is_snp=self.is_snp)
        return obj


    



    
        
    