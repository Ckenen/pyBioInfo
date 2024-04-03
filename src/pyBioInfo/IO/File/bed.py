import os
import pysam
from .file import BaseFile
from pygz import PigzFile
from pyBioInfo.Range import TRange


class BedRecord(TRange):
    pass


class BedFile(BaseFile):
    def __init__(self, path, mode="r", ncol=12):
        super(BedFile, self).__init__(path, mode)
        self._ncol = ncol
        self.open()

    @property
    def ncol(self):
        return self._ncol

    def open(self):
        if self._mode == "r":
            if self._path.endswith(".gz"):
                self._handle = PigzFile(self.path, "rt")
            else:
                self._handle = open(self.path)
        elif self.mode == "w":
            if self.path.endswith(".gz"):
                self.handle = PigzFile(self.path, "wt")
            else:
                self.handle = open(self.path, "w+")
        else:
            raise ValueError()

    def close(self):
        if self._handle:
            self._handle.close()
            self._handle = None

    def fetch(self, chrom=None, start=None, end=None):
        if chrom is None:
            for line in self._handle:
                line = line.strip("\n").strip()
                if line.startswith("#"):
                    continue
                elif line == "":
                    continue
                else:
                    yield self.parse_bed_string(line, self._ncol)
        else:
            raise RuntimeError("Random access is not supported! Please try BedFileRandom.")

    def __iter__(self):
        for x in self.fetch():
            yield x
    
    def write(self, record):
        line = self.to_bed_string(record, self._ncol)
        self.handle.write(line + "\n")

    @classmethod
    def parse_bed_string(cls, line, column=12):
        row = line.split("\t")
        column = min(len(row), column)
        assert column >= 3

        chrom = row[0]
        start = int(row[1])
        end = int(row[2])

        name = None
        if column >= 4:
            name = row[3]

        strand = None
        if column >= 6:
            strand = row[5]

        thick = None
        if column >= 8:
            thick_start = int(row[6])
            thick_end = int(row[7])
            if thick_start < thick_end:
                thick = (thick_start, thick_end)

        blocks = None
        if column >= 12:
            # print(row)
            block_sizes = map(int, filter(lambda p: p != "", row[10].split(",")))
            block_offsets = map(int, filter(lambda p: p != "", row[11].split(",")))
            # print(block_sizes, block_offsets)
            temps = []
            for offset, size in zip(block_offsets, block_sizes):
                x = start + offset
                y = x + size
                temps.append((x, y))
            blocks = temps

        obj = TRange(chrom=chrom, start=start, end=end, name=name, strand=strand, thick=thick, blocks=blocks)

        if column >= 5:
            score = row[4]
            if score != ".":
                score = float(score)
                obj.score = score

        if column >= 9:
            color = row[8]
            obj.color = color

        return obj

    @classmethod
    def to_bed_string(cls, record, column=12):
        assert column >= 3

        columns = []
        if column >= 1:
            columns.append(record.chrom)
        if column >= 2:
            columns.append(record.start)
        if column >= 3:
            columns.append(record.end)
        if column >= 4:
            name = record.name
            if name is None:
                name = "None"
            columns.append(name)
        if column >= 5:
            try:
                score = record.score
            except AttributeError:
                score = "."
            if score is None:
                score = "."
            columns.append(score)
        if column >= 6:
            columns.append(record.strand)
        if column >= 7:
            thick_start = record.thick[0] if record.thick is not None else record.start
            columns.append(thick_start)
        if column >= 8:
            thick_end = record.thick[1] if record.thick is not None else record.start
            columns.append(thick_end)
        if column >= 9:
            try:
                color = record.color
            except AttributeError:
                color = "0,0,255" if record.is_forward else "255,0,0"
            columns.append(color)
        if column >= 10:
            columns.append(len(record.blocks))
        if column >= 11:
            block_sizes = [block_end - block_start for block_start, block_end in record.blocks]
            block_sizes_string = ",".join(map(str, block_sizes)) + ","
            columns.append(block_sizes_string)
        if column >= 12:
            block_offsets = [block_start - record.start for block_start, block_end in record.blocks]
            block_offsets_string = ",".join(map(str, block_offsets)) + ","
            columns.append(block_offsets_string)
        line = "\t".join(map(str, columns))
        return line
 

class BedFileRandom(BedFile):
    def __init__(self, path, ncol=12):
        assert path.endswith(".bed.gz")
        assert os.path.exists(path)
        assert os.path.exists(path + ".tbi")
        super(BedFileRandom, self).__init__(path, "r", ncol)

    def fetch(self, chrom=None, start=None, end=None):
        for line in self._handle.fetch(chrom, start, end):
            yield BedFile.parse_bed_string(line, self._ncol)

    def open(self):
        if self._handle is None:
            self._handle = pysam.TabixFile(self._path)

    def close(self):
        if self._handle is not None:
            self._handle.close()
            self._handle = None

