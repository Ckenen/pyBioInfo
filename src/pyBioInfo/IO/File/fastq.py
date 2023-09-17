# import gzip
from pygz import PigzFile
# from .file_stream import TextFile


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


class FastqFile(object):
    def __init__(self, path, mode="r"):
        self.path = path
        self.mode = mode
        self._rtype = None
        self._wtype = None
        self._handle = None
        self._iterator = None

        if self.mode == "r":
            if self.path.endswith(".gz"):
                self._handle = PigzFile(self.path, "rt")
            else:
                self._handle = open(self.path)
        elif self.mode == "w":
            if self.path.endswith(".gz"):
                self._handle = PigzFile(self.path, "wt")
            else:
                self._handle = open(self.path, "w+")
        else:
            raise RuntimeError()

    def __iter__(self):
        name = None
        sequence = None
        for i, line in enumerate(self._handle):
            line = line.strip("\n")
            n = i % 4
            if n == 0:
                assert line[0] == "@"
                name = line[1:]
            elif n == 1:
                sequence = line
            elif n == 2:
                continue
            else:
                quality = line
                yield FastqRecord(name, sequence, quality)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def close(self):
        if self._handle:
            self._handle.close()
            self._handle = None

    def write(self, record):
        if isinstance(record, FastqRecord):
            self._handle.write("@%s\n" % record.name)
            self._handle.write("%s\n" % record.sequence)
            self._handle.write("+\n")
            self._handle.write("%s\n" % record.quality)
        else:
            raise TypeError()

    def __next__(self):
        if self._iterator is None:
            self._iterator = self.__iter__()
        return next(self._iterator)

