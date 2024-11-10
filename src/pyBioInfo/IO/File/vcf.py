import os
import gzip
import pysam
from pyBioInfo.Range import GRange
from .base import BaseFile


REVERSE_COMPLEMENT_MAPPER = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "N": "N",
}


class VcfRecord(GRange):
    def __init__(self, chrom, pos, name, ref, alts, qual, filt, info):
        start = pos - 1
        end = pos
        super(VcfRecord, self).__init__(chrom=chrom,
                                        start=start, 
                                        end=end, 
                                        name=name, 
                                        strand="+")
        assert len(ref) == 1
        self.ref = ref
        if isinstance(alts, str):
            alts = [alts]
        self.alts = alts
        self.qual = qual
        self.filt = filt
        self.info = info

    def induce(self, ref, alt):
        if ref == self.ref:
            return alt in self.alts
        elif REVERSE_COMPLEMENT_MAPPER[ref] == self.ref:
            return REVERSE_COMPLEMENT_MAPPER[alt] in self.alts
        else:
            print(ref, alt, self.ref, self.alts)
            return False


class VcfFile(BaseFile):
    def __init__(self, path, mode="r", template=None, header=None):
        self._template = template
        self._header = header
        self._comments = []
        super(VcfFile, self).__init__(path, mode)
        self.open()
        
    def open(self):
        if self._handle is None:
            if self._mode == "r":
                if self._path.endswith(".gz"):
                    self._handle = gzip.open(self._path, "rt")
                else:
                    self._handle = open(self._path)
            elif self._mode == "w":
                raise NotImplementedError()
                
    def close(self):
        if self._handle is not None:
            self._handle.close()
            self._handle = None

    def fetch(self, chrom=None, start=None, end=None):
        if chrom is None:
            for line in self._handle:
                line = line.strip("\n")
                if line == "":
                    continue
                if line.startswith("#"):
                    self._comments.append(line)
                else:
                    yield self.parse_vcf_string(line)
        else:
            raise RuntimeError() 

    def __iter__(self):
        for x in self.fetch():
            yield x

    def write(self, record):
        assert "w" in self.mode
        if self.handle is None:
            if self.path.endswith(".gz"):
                self.mode = "wb"
                self.handle = gzip.open(self.path, self.mode)
            else:
                self.mode = "w+"
                self.handle = open(self.path, self.mode)
            for comment in self.comments:
                self.handle.write(comment + "\n")
        columns = [record.chrom, record.start + 1, record.name,
                   record.ref, ",".join(record.alts), record.qual, record.filt]
        info = {"SB": record.strand}
        for key, value in record.info.items():
            info[key] = value
        array = []
        for key, value in info.items():
            array.append("%s=%s" % (key, value))
        columns.append(";".join(array))
        self.handle.write("\t".join(map(str, columns)) + "\n")

    @classmethod
    def parse_vcf_string(cls, line):
        columns = line.split("\t")
        chrom, pos, name, ref, alts, qual, filt, info = columns
        pos = int(pos)
        alts = alts.split(",")
        qual = int(qual)
        temp = dict()
        if info != ".":
            for item in info.split(";"):            
                if item == "":
                    continue
                key, value = item.split("=")
                temp[key] = value
        info = temp
        record = VcfRecord(chrom, pos=pos, name=name, ref=ref,
            alts=alts, qual=qual, filt=filt, info=info)
        return record 
    

def VcfFileRandom(VcfFile):
    pass