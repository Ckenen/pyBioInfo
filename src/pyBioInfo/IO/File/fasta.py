import os
import gzip
from .base import BaseFile


REVERSE_COMLEMENT_MAPPER = {
    "A": "T", "C": "G", "G": "C", "T": "A", "N": "N",
    "a": "t", "c": "g", "g": "c", "t": "a", "n": "n"}


class FastaRecord(object):
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
    
    @property
    def seq(self):
        return self.sequence

    def __str__(self):
        return ">%s\n%s\n" % (self.name, self.sequence)

    def reverse_comlement(self):
        s = "".join([REVERSE_COMLEMENT_MAPPER[b] for b in self.sequence[::-1]])
        return FastaRecord(self.name, s)
    

class FastaFile(BaseFile):
    def __init__(self, path, mode="r", random=None):
        if mode == "r" and random is None:
            random = os.path.exists(path + ".fai")
        
        if random:
            assert mode == "r"
            assert os.path.exists(path)
            assert os.path.exists(path + ".fai")
        super(FastaFile, self).__init__(path, mode)
        self._random = random
        self._index = None
        self.open()
        
    def _load_index(self):
        self._index = dict()
        with open(self._path + ".fai") as f:
            for line in f:
                row = line.strip("\n").split("\t")
                seqname = row[0]
                length = int(row[1])
                offset = int(row[2])
                width1 = int(row[3])
                width2 = int(row[4])
                self._index[seqname] = [length, offset, width1, width2]
                
    def open(self):
        if self._handle is None:
            if self._random:
                assert self._mode == "r"
                self._handle = open(self._path)
                self._load_index()
            else:
                if self._mode == "r":
                    if self._path.endswith(".gz"):
                        self._handle = gzip.open(self._path, "rt")
                    else:
                        self._handle = open(self._path)
                elif self._mode == "w":
                    if self._path.endswith(".gz"):
                        self._handle = gzip.open(self._path, "wt")
                    else:
                        self._handle = open(self._path, "w+")
                else:
                    raise ValueError()

    def close(self):
        if self._handle:
            self._handle.close()
            self._handle = None
            
    def get_chrom_length(self, chrom):
        return self._index[chrom][0]
            
    def _get_sequence_by_coordinates(self, chrom, start=None, end=None):
        assert self._random
        h = self._handle
        length, offset, width1, width2 = self._index[chrom]
        if start is None:
            start = 0
        if end is None:
            end = length
        assert start >= 0 and end <= length
        n1 = int(start / width1) * width2 + start % width1
        n2 = int((end - 1) / width1) * width2 + (end - 1) % width1 + 1
        h.seek(offset + n1, 0)
        return h.read(n2 - n1).replace("\n", "")
    
    def _get_sequence_by_blocks(self, chrom, blocks):
        assert self._random
        ss = []
        for start, end in blocks:
            s = self._get_sequence_by_coordinates(chrom, start, end)
            ss.append(s)
        return "".join(ss)

    def fetch(self, chrom=None, start=None, end=None, strand="+", obj=None, strandness=True):
        if self._random:
            if obj:
                strand = obj.strand
                s = self._get_sequence_by_blocks(obj.chrom, obj.blocks)
                if strand == "-" and strandness:
                    s = self.reverse_complement(s)
                yield FastaRecord("unamed", s)
            else:
                if chrom is None:
                    for chrom2 in list(sorted(self._index)):
                        s = self._get_sequence_by_coordinates(chrom2)
                        if strand == "-" and strandness:
                            s = self.reverse_complement(s)
                        yield FastaRecord(chrom2, s)
                else:
                    if start is None:
                        start = 0
                    if end is None:
                        end = self.get_chrom_length(chrom)
                    s = self._get_sequence_by_coordinates(chrom, start, end)
                    if strand == "-" and strandness:
                        s = self.reverse_complement(s)
                    yield FastaRecord("%s:%d-%s" % (chrom, start, end), s)
        else:
            if chrom is None:
                name = None
                rows = None
                for line in self.handle:
                    line = line.strip("\n").strip()
                    if line.startswith(">"):
                        if name is not None:
                            seq = "".join(rows)
                            yield FastaRecord(name, seq)
                        name = line[1:]
                        rows = []
                    else:
                        if line != "":
                            rows.append(line)
                if name is not None:
                    seq = "".join(rows)
                    yield FastaRecord(name, seq)
            else:
                raise RuntimeError()
            
    def fast_fetch(self, chrom=None, start=None, end=None, strand="+", obj=None, strandness=True):
        if self._random:
            if obj:
                strand = obj.strand
                s = self._get_sequence_by_blocks(obj.chrom, obj.blocks)
                if strand == "-" and strandness:
                    s = self.reverse_complement(s)
                return s
            else:
                if chrom is None:
                    raise RuntimeError("Must provided chrom.")
                else:
                    if start is None:
                        start = 0
                    if end is None:
                        end = self.get_chrom_length(chrom)
                    s = self._get_sequence_by_coordinates(chrom, start, end)
                    if strand == "-" and strandness:
                        s = self.reverse_complement(s)
                    return s
        else:
            raise RuntimeError("Must random access.")

    def write(self, record=None, name=None, sequence=None, width=50):
        if record:
            name = record.name
            sequence = record.sequence
        self._handle.write(">%s\n" % name)
        if width is None:
            self._handle.write(">%s\n" % sequence)
        else:
            for i1 in range(0, len(sequence), width):
                i2 = min(i1 + width, len(sequence))
                seq2 = sequence[i1:i2]
                self.handle.write(seq2 + "\n")

    @classmethod
    def reverse_complement(cls, seq):
        return "".join([REVERSE_COMLEMENT_MAPPER[b] for b in seq[::-1]])
