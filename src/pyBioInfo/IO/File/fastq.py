import gzip
from .base import BaseFile

class FastqRecord(object):
    def __init__(self, name, sequence, quality):
        self.name = name
        self.sequence = sequence
        self.quality = quality
        assert len(self.sequence) == len(self.quality)

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        s = "@%s\n%s\n%s\n%s\n" % (self.name, self.sequence, "+", self.quality)
        return s


class FastqFile(BaseFile):
    def __init__(self, path, mode="r"):
        assert mode == "r" or mode == "w"
        super(FastqFile, self).__init__(path, mode)
        self.open()
        
    def open(self):
        if self._handle is None:
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
        if self._handle is not None:
            self._handle.close()
            self._handle = None
    
    def fetch(self):
        name, sequence = None, None
        for i, line in enumerate(self._handle):
            j = i % 4
            if j == 0:
                assert line[0] == "@"
                name = line[1:-1]
            elif j == 1:
                sequence = line[:-1]
            elif j == 2:
                continue
            else:
                quality = line[:-1]
                yield FastqRecord(name, sequence, quality)

    def write(self, record):
        if isinstance(record, FastqRecord):
            self._handle.write("@%s\n" % record.name)
            self._handle.write("%s\n" % record.sequence)
            self._handle.write("+\n")
            self._handle.write("%s\n" % record.quality)
        else:
            raise TypeError()