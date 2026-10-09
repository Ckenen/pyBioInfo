# Class: FastqFile (reader and writer of FASTQ format files)

import gzip
import sys
from collections.abc import Iterator
from typing import TextIO, cast

if sys.version_info >= (3, 12):
    from typing import override
else:
    from typing import Callable  # pyright: ignore[reportUnreachable]

    def override(func: Callable) -> Callable:
        """No-op fallback for Python < 3.12 compatibility."""
        return func

from .base import BaseFile


class FastqRecord(object):
    """
    FastqRecord is a class for a single record of FASTQ files.

    It holds the read ``name`` (without the "@" prefix), the
    ``sequence`` and the Phred ``quality`` string.
    """

    def __init__(self, name: str, sequence: str, quality: str) -> None:
        """
        Initialize a FASTQ record.

        :param name: Name of the read.
        :param sequence: Read sequence string.
        :param quality: Quality string of the read.
        :raises AssertionError: If sequence and quality lengths differ.
        """
        self.name: str = name
        self.sequence: str = sequence
        self.quality: str = quality
        assert len(self.sequence) == len(self.quality)

    def __len__(self) -> int:
        """Length of the read sequence."""
        return len(self.sequence)

    @override
    def __str__(self) -> str:
        """Format the record as a four-line FASTQ text."""
        return "@%s\n%s\n%s\n%s\n" % (self.name, self.sequence, "+", self.quality)


class FastqFile(BaseFile):
    """
    FastqFile is a class for reading and writing FASTQ files.

    Plain text FASTQ and gzip-compressed FASTQ (``.fastq.gz``) are
    supported. Records are read sequentially (iteration) in groups of
    four lines and parsed into :py:class:`FastqRecord` objects; records
    are written by :py:meth:`write`.
    """

    # Attributes assigned by the untyped BaseFile base class; declared
    # here to keep their types known.
    _path: str
    _mode: str
    _handle: "TextIO | None"

    def __init__(self, path: str, mode: str = "r") -> None:
        """
        Initialize and open a FASTQ file object.

        :param path: Path to the FASTQ (or compressed .fastq.gz) file.
        :param mode: File mode; "r" for reading and "w" for writing.
        :raises AssertionError: If the mode is unsupported.
        """
        assert mode == "r" or mode == "w", "Mode %s is not supported!" % mode
        super(FastqFile, self).__init__(path, mode)
        self.open()

    @override
    def open(self) -> None:
        """
        Open the underlying file handle.

        Calling open on an already-opened file is a no-op.

        :raises ValueError: If the mode is unsupported.
        """
        if self._handle is None:
            if self._mode == "r":
                if self._path.endswith(".gz"):
                    self._handle = cast("TextIO", gzip.open(self._path, "rt"))
                else:
                    self._handle = open(self._path)
            elif self._mode == "w":
                if self._path.endswith(".gz"):
                    self._handle = cast("TextIO", gzip.open(self._path, "wt"))
                else:
                    self._handle = open(self._path, "w+")
            else:
                raise ValueError("Mode %s is not supported!" % self._mode)

    @override
    def close(self) -> None:
        """
        Close the underlying file handle.

        Calling close on an already-closed file is a no-op, so that
        explicit close followed by ``__del__`` does not raise.
        """
        if self._handle is not None:
            self._handle.close()
            # Note: the narrowed handle type is intentionally incompatible
            # with the heterogeneous "IO[Any] | None" of the base class.
            self._handle = None  # pyright: ignore[reportIncompatibleVariableOverride]

    @override
    def fetch(self) -> "Iterator[FastqRecord]":
        """
        Fetch the records of the file sequentially.

        Records are read in groups of four lines (name, sequence,
        "+" line and quality).

        :raises RuntimeError: If the file is not opened.
        :return: Iterator of FastqRecord objects.
        """
        handle = self._handle
        assert handle is not None, "File is not opened!"
        name, sequence = None, None
        for i, line in enumerate(handle):
            j = i % 4
            if j == 0:
                assert line[0] == "@", "Invalid FASTQ record name line!"
                name = line[1:-1]
            elif j == 1:
                sequence = line[:-1]
            elif j == 2:
                continue
            else:
                quality = line[:-1]
                assert name is not None and sequence is not None
                yield FastqRecord(name, sequence, quality)

    def write(self, record: FastqRecord) -> None:
        """
        Write a FASTQ record to the file.

        :param record: The FastqRecord to write.
        :raises RuntimeError: If the file is not opened.
        """
        handle = self._handle
        assert handle is not None, "File is not opened!"
        _ = handle.write("@%s\n" % record.name)
        _ = handle.write("%s\n" % record.sequence)
        _ = handle.write("+\n")
        _ = handle.write("%s\n" % record.quality)
