import abc
from abc import ABC
import gzip



class BioInfoFile(object):

    # __metaclass__ == abc.ABCMeta

    def __init__(self, path=None, mode=None):
        self._path = path
        self._mode = path
        self._status = None
        self._handle = None
        self._iterator = None
        self._start_iterate = False
        self._comment_char = "#"
        self._comment_lines = list()
        self._is_gz = False
        self._is_index = False

    @property
    def path(self):
        return self._path

    @property
    def mode(self):
        return self._mode

    @property
    def status(self):
        return self._status

    @property
    def handle(self):
        return self._handle

    @abc.abstractmethod
    def open(self):
        pass

    def close(self):
        if self._handle:
            self._handle.close()
        self._handle = None

    @abc.abstractmethod
    def read(self):
        pass

    # def reads(self, count):
    #     pass

    @abc.abstractmethod
    def fetch(self, *args, **kwargs):
        pass

    @abc.abstractmethod
    def write(self, obj):
        pass

    def writes(self, objs):
        for obj in objs:
            self.write(obj)

    def __iter__(self):
        assert self._start_iterate is False
        self._start_iterate = True
        for obj in self.read():
            yield obj

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class FileStream(object):
    __slots__ = ["_path", "_mode", "_status",
                 "_handle", "_iterator", "_has_start_iterate"]

    STATUS_INITED = 0
    STATUS_OPENED = 1
    STATUS_CLOSED = 2

    def __init__(self, source, mode="r"):
        self._path = source
        self._mode = mode
        self._status = self.STATUS_INITED
        self._handle = None
        self._iterator = None
        self._has_start_iterate = False
        self.open()

    @property
    def path(self):
        return self._path

    @property
    def mode(self):
        return self._mode

    @property
    def status(self):
        return self._status

    @property
    def handle(self):
        return self._handle

    def open(self):
        if self._status != self.STATUS_INITED:
            raise RuntimeError("This instance has been opened previously.")
        self._init_handle()
        self._status = self.STATUS_OPENED

    def _init_handle(self):
        raise RuntimeError("This method must be overridden in the subclass.")

    def close(self):
        if self._status != self.STATUS_OPENED:
            raise RuntimeError("This instance has not yet been opened.")
        if self._handle:
            self._handle.close()
        self._status = self.STATUS_CLOSED

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        if self._status == self.STATUS_OPENED:
            self.close()

    def _start_iterate(self):
        if "r" not in self.mode:
            raise RuntimeError("File can not be read!")
        if self._has_start_iterate:
            raise RuntimeError("Can only iterate once.")
        self._has_start_iterate = True

    def __iter__(self):
        self._start_iterate()
        for record in self._read_record():
            yield record

    def _read_record(self):
        raise RuntimeError("This method must be overridden in the subclass.")

    def __next__(self):
        if self._iterator is None:
            self._iterator = self.__iter__()
        return next(self._iterator)

    def write(self, record):
        if "w" not in self.mode:
            raise RuntimeError("File can not be wroted!")
        self._write_record(record)

    def _write_record(self, record):
        raise RuntimeError("This method must be overridden in the subclass.")


class TextFile(FileStream, ABC):

    def __init__(self, path, mode):
        super(TextFile, self).__init__(path, mode)
        self.comments = []
        self.compressed = False

    def _init_handle(self):
        self.compressed = self.path.endswith(".gz")
        if "r" in self.mode:
            if self.compressed:
                self._handle = gzip.open(self.path, self.mode)
            else:
                self._handle = open(self.path)
        elif "w" in self.mode:
            if self.compressed:
                self._handle = gzip.open(self.path, "wt")
            else:
                self._handle = open(self.path, "w+")
        else:
            raise ValueError("Unknown mode[%s]." % self.mode)
