from doctest import set_unittest_reportflags
import os
import gzip
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class FastaFile(object):
    def __init__(self, path, mode="r"):
        self.path = path
        self.mode = mode
        self.handle = None

        if self.mode == "r":
            if os.path.exists(self.path + ".fai"):
                self.handle = pysam.FastaFile(self.path)
            else:
                if self.path.endswith(".gz"):
                    self.handle = gzip.open(self.path, "rt")
                else:
                    self.handle = open(self.path)
        elif self.mode == "w":
            if self.path.endswith(".gz"):
                self.handle = gzip.open(self.path, "wt")
            else:
                self.handle = open(self.path, "w+")
        else:
            raise RuntimeError()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        if self.handle:
            self.handle.close()
            self.handle = None

    def __iter__(self):
        name = None
        rows = None
        for line in self.handle:
            line = line.strip("\n").strip()
            if line.startswith(">"):
                if name is not None:
                    sequence = "".join(rows)
                    obj = self._construct_seq_record(name, sequence)
                    yield obj
                name = line[1:]
                rows = []
            else:
                if line != "":
                    rows.append(line)
        if name is not None:
            sequence = "".join(rows)
            obj = self._construct_seq_record(name, sequence)
            yield obj

    def write(self, record):
        self.handle.write(">%s\n" % record.name)
        sequence = record.seq
        for i1 in range(0, len(sequence), 50):
            i2 = min(i1 + 50, len(sequence))
            part = sequence[i1:i2]
            self.handle.write(part + "\n")

    @classmethod
    def _construct_seq_record(cls, name, sequence):
        return SeqRecord(id=name, name=name, description=name, seq=sequence)

    def _fetch_by_chrom_start_end(self, chrom, start, end):
        sequence = self.handle.fetch(chrom, start, end)
        return sequence

    def _fetch_by_genomic_range(self, obj):
        sequences = []
        for block_start, block_end in obj.blocks:
            sequence = self._fetch_by_chrom_start_end(
                obj.chrom, block_start, block_end)
            sequences.append(sequence)
        sequence = "".join(sequences)
        return sequence

    def fetch(self, chrom=None, start=None, end=None, obj=None, strand="+", strandness=True, relative=True):
        if obj:
            sequence = self._fetch_by_genomic_range(obj)
            strand = obj.strand
        else:
            sequence = self._fetch_by_chrom_start_end(chrom, start, end)

        reverse = False
        if strandness:
            if strand == "+":
                if not relative:
                    reverse = True
            else:
                if relative:
                    reverse = True

        if reverse:
            sequence = str(Seq(sequence).reverse_complement())

        return sequence

#


class FastaRecord(object):
    MAPPER = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "N": "N",
        "a": "t",
        "c": "g",
        "g": "c",
        "t": "a",
        "n": "n"
    }

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def __str__(self):
        return ">%s\n%s\n" % (self.name, self.sequence)

    def reverse_comlement(self):
        s = "".join([self.MAPPER[b] for b in self.sequence[::-1]])
        return FastaRecord(self.name, s)


class FastaReader(object):
    def __init__(self, path):
        self.path = path
        self.handle = open(self.path)

    def __iter__(self):
        name = None
        tmp = None
        for line in self.handle:
            line = line.strip("\n")
            if line.startswith(">"):
                if tmp:
                    yield FastaRecord(name, "".join(tmp))
                name = line[1:].split()[0]
                tmp = []
            else:
                tmp.append(line)
        if tmp:
            yield FastaRecord(name, "".join(tmp))

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.handle:
            self.handle.close()
            self.handle = None


class FastaWriter(object):
    def __init__(self, path):
        self.path = path
        self.handle = open(self.path, "w+")

    def write(self, record=None, name=None, sequence=None, length=None):
        if record:
            name = record.name
            sequence = record.sequence
        assert name
        assert sequence
        self.handle.write(">")
        self.handle.write(name)
        self.handle.write("\n")
        if length:
            for m in range(0, len(sequence), length):
                n = min(m + length, len(sequence))
                self.handle.write(sequence[m:n])
                self.handle.write("\n")
        else:
            self.handle.write(sequence)
            self.handle.write("\n")

    def close(self):
        if self.handle:
            self.handle.close()
            self.handle = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class FastaRandomReader(object):
    MAPPER = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "N": "N",
        "a": "t",
        "c": "g",
        "g": "c",
        "t": "a",
        "n": "n"
    }

    def __init__(self, path):
        self.path = path
        self.handle = open(self.path)

        self.index = dict()
        with open(self.path + ".fai") as f:
            for line in f:
                row = line.strip("\n").split("\t")
                length = int(row[1])
                offset = int(row[2])
                width1 = int(row[3])
                width2 = int(row[4])
                self.index[row[0]] = [length, offset, width1, width2]

    def fetch(self, obj=None, chrom=None, start=None, end=None, strand="+"):
        if obj:
            strand = obj.strand
            ss = []
            for bs, be in obj.blocks():
                ss.append(self.fetch(chrom=obj.chrom, start=bs, end=be))
            s = "".join(ss)
        else:
            length, offset, width1, width2 = self.index[chrom]
            if start is None:
                start = 0
            if end is None:
                end = length
            assert start >= 0 and end <= length
            n1 = int(start / width1) * width2 + start % width1
            n2 = int((end - 1) / width1) * width2 + (end - 1) % width1 + 1
            self.handle.seek(offset + n1, 0)
            s = self.read(n2 - n1).replace("\n", "")
        if strand == "-":
            s = "".join([self.MAPPER[b] for b in s[::-1]])
        return s

    def close(self):
        if self.handle:
            self.handle.close()
            self.handle = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
