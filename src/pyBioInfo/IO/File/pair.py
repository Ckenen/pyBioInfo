import os, logging
from pyBioInfo.Range import GRange
from .base import BaseFile
import gzip


class PairRange(GRange):
    def __init__(self, name, chrom1, start1, end1, strand1, chrom2, start2, end2, strand2, score1=None, score2=None):
        super(PairRange, self).__init__(chrom=chrom1, start=start1, end=end1, strand=strand1, name=name)
        self.chrom1 = chrom1
        self.start1 = start1
        self.end1 = end1
        self.strand1 = strand1
        self.chrom2 = chrom2
        self.start2 = start2
        self.end2 = end2
        self.strand2 = strand2
        self.score1 = score1
        self.score2 = score2
        self.additions = None
        
    def __str__(self):
        return "%s:%d-%d(%s)|%s:%d-%d(%s)" % (
            self.chrom1, self.start1, self.end1, self.strand1, 
            self.chrom2, self.start2, self.end2, self.strand2)
        
    def format(self, fmt="pair"):
        if fmt == "pair":
            items = [
                self.name, 
                self.chrom1, 
                self.start1 + 1, 
                self.chrom2, 
                self.start2 + 1,
                self.strand1, 
                self.strand2, 
                "." if self.score1 is None else self.score1, 
                "." if self.score2 is None else self.score2,
                self.end1 - self.start1,
                self.end2 - self.start2]
            if self.additions is not None:
                items = items + self.additions
            return "\t".join(map(str, items))
        raise ValueError()
        
    
class PairFile(BaseFile):
    def __init__(self, path, mode="rb", random=None):
        if mode[0] == "r" and random is None:
            if path.endswith(".pairs.gz"):
                if os.path.exists(path + ".px2"):
                    logging.warning("Index file exists, automatic set random=True")
                    random = True
                else:
                    logging.warning("Index file does not exists, automatic set random=False")
                    random = False
            else:
                random = False
        
        if random:
            assert path.endswith(".pairs.gz")
            assert os.path.exists(path)
            assert os.path.exists(path + ".px2")
        assert mode == "rb" or mode == "wb"
        self._references = None
        self._random = random
        super(PairFile, self).__init__(path, mode)
        self.open()
        
    def open(self):
        if self._handle is None:
            if self.mode == "rb":
                self._handle = gzip.open(self.path, "rt")
            elif self.mode == "wb":
                assert False
            else:
                raise RuntimeError()
    
    def close(self):
        if self._handle:
            self._handle.close()
            self._handle = None
            
    def fetch(self, chrom1=None, start1=None, end1=None, chrom2=None, start2=None, end2=None):
        if self._random:
            query = None
            if chrom1 is None:
                if chrom2 is None:
                    for line in self.handle:
                        yield self.parse_line(line)
                else:
                    if start2 is None:
                        start2 = 0
                    if end2 is None:
                        end2 = 999999999
                    query = "|%s:%d-%d" % (chrom2, start2, end2)
            else:
                if chrom2 is None:
                    if start1 is None:
                        start1 = 0
                    if end1 is None:
                        end1 = 999999
                    query = "%s:%d-%d|%" % (chrom1, start1, end1)
                else:
                    if start1 is None:
                        start1 = 0
                    if end1 is None:
                        end1 = 999999
                    if start2 is None:
                        start2 = 0
                    if end2 is None:
                        end2 = 999999999
                    query = "%s:%d-%d|%s:%d-%d" % (chrom1, start1, end1, chrom2, start2, end2)
            if query is not None:
                raise NotImplementedError()
        else:
            if chrom1 is None and start1 is None and end1 is None \
                and chrom2 is None and start2 is None and end2 is None:
                for line in self.handle:
                    yield self.parse_line(line)
            else:
                raise RuntimeError("Random access is not supported! Please set random=True.")
    
    @classmethod
    def parse_line(cls, x):
        if isinstance(x, str):
            row = x.strip().split("\t")
        elif isinstance(x, list):
            row = x
        else:
            raise RuntimeError()
        name = row[0]
        chrom1, start1 = row[1], int(row[2]) - 1
        chrom2, start2 = row[3], int(row[4]) - 1
        strand1, strand2 = row[5], row[6]
        score1 = None if row[7] == "." else int(row[7])
        score2 = None if row[8] == "." else int(row[8])
        length1, length2 = int(row[9]), int(row[10])
        end1, end2 = start1 + length1, start2 + length2
        pr = PairRange(name=name, 
                       chrom1=chrom1, start1=start1, end1=end1, strand1=strand1, score1=score1, 
                       chrom2=chrom2, start2=start2, end2=end2, strand2=strand2, score2=score2)
        if len(row) > 11:
            pr.additions = row[11:]
        return pr
        