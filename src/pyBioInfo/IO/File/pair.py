# Class: PairFile (reader of PAIR format files)

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
from ...Range import CRange


class PairRange(CRange):
    """
    PairRange is a class for a single record of PAIR files.

    A PAIR record describes the two ends of a contact (e.g. a Hi-C
    read pair) with their coordinates, strands, mapqs and counts.
    """

    def __init__(self, chrom1: str, start1: int, end1: int, strand1: str,
                 chrom2: str, start2: int, end2: int, strand2: str,
                 mapq1: int, mapq2: int, count: int = 1) -> None:
        """
        Initialize a pair record.

        :param chrom1: Chromosome name of the first end.
        :param start1: Start of the first end (0-based, inclusive).
        :param end1: End of the first end (0-based, exclusive).
        :param strand1: Strand of the first end ("+" or "-").
        :param chrom2: Chromosome name of the second end.
        :param start2: Start of the second end (0-based, inclusive).
        :param end2: End of the second end (0-based, exclusive).
        :param strand2: Strand of the second end ("+" or "-").
        :param mapq1: Mapping quality of the first end.
        :param mapq2: Mapping quality of the second end.
        :param count: Number of contacts represented by the record.
        """
        super(PairRange, self).__init__(chrom=chrom1, start=start1, end=end1,
                                        strand=strand1)
        self.chrom2: str = chrom2
        self.start2: int = start2
        self.end2: int = end2
        self.strand2: str = strand2
        self.mapq1: int = mapq1
        self.mapq2: int = mapq2
        self.count: int = count

    @override
    def format(self, fmt: str = "pair") -> str:  # pyright: ignore[reportIncompatibleMethodOverride]
        """
        Format the record as a line of the given text format.

        :param fmt: Format name, case-insensitive; only "pair" is
                    supported.
        :return: Tab-separated PAIR line.
        :raises NotImplementedError: If the format is unsupported.
        """
        fmt = fmt.upper()
        if fmt == "PAIR":
            return "%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d\n" % (
                self.chrom, self.start, self.end, self.strand,
                self.chrom2, self.start2, self.end2, self.strand2,
                self.mapq1, self.mapq2, self.count)
        else:
            raise NotImplementedError("Format %s is not supported!" % fmt)

    @override
    def __str__(self) -> str:
        """Format the record as a PAIR line."""
        return self.format(fmt="PAIR")


class PairFile(BaseFile):
    """
    PairFile is a class for reading PAIR files.

    Plain text PAIR and gzip-compressed PAIR (``.pair.gz``) are
    supported. Records are read sequentially (iteration), optionally
    filtered by a mapping-quality threshold.
    """

    # Attributes assigned by the untyped BaseFile base class; declared
    # here to keep their types known.
    _path: str
    _mode: str
    _handle: "TextIO | None"

    def __init__(self, path: str, mode: str = "rb",
                 mapq: int = 255) -> None:
        """
        Initialize and open a PAIR file object.

        :param path: Path to the PAIR (or compressed .pair.gz) file.
        :param mode: File mode; only "rb" (reading) is supported.
        :param mapq: Minimum mapping quality of both ends for a record
                     to be yielded by :py:meth:`fetch`.
        :raises NotImplementedError: If the mode is unsupported.
        :raises AssertionError: If the mode is invalid.
        """
        assert mode == "rb" or mode == "wb", "Mode %s is not supported!" % mode
        if mode == "wb":
            raise NotImplementedError("Writing is not supported!")
        super(PairFile, self).__init__(path, mode)
        self._mapq: int = mapq
        self.open()

    @property
    def mapq(self) -> int:
        """Minimum mapping-quality threshold."""
        return self._mapq

    @override
    def open(self) -> None:
        """
        Open the underlying file handle.

        Calling open on an already-opened file is a no-op.

        :raises ValueError: If the mode is unsupported.
        """
        if self._handle is None:
            if self._mode == "rb":
                if self._path.endswith(".gz"):
                    self._handle = cast("TextIO", gzip.open(self._path, "rt"))
                else:
                    self._handle = open(self._path)
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
    def fetch(self) -> "Iterator[PairRange]":
        """
        Fetch the records of the file sequentially.

        Records whose both ends have a mapping quality below the
        :py:attr:`mapq` threshold are skipped.

        :raises RuntimeError: If the file is not opened.
        :return: Iterator of PairRange objects.
        """
        handle = self._handle
        assert handle is not None, "File is not opened!"
        for line in handle:
            row = line.strip("\n").split("\t")
            assert len(row) == 11
            if int(row[8]) < self._mapq and int(row[9]) < self._mapq:
                continue
            chrom1 = row[0]
            start1 = int(row[1])
            end1 = int(row[2])
            strand1 = row[3]
            chrom2 = row[4]
            start2 = int(row[5])
            end2 = int(row[6])
            strand2 = row[7]
            mapq1 = int(row[8])
            mapq2 = int(row[9])
            count = int(row[10])
            record = PairRange(chrom1, start1, end1, strand1,
                               chrom2, start2, end2, strand2,
                               mapq1, mapq2, count)
            yield record

    def write(self, record: PairRange) -> None:  # pyright: ignore[reportUnusedParameter]
        """
        Write a pair record to the file.

        :param record: The PairRange to write.
        :raises NotImplementedError: Always; writing is not supported.
        """
        raise NotImplementedError("Writing is not supported!")
