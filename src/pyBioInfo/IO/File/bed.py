# Class: BedFile (reader and writer of BED format files)

import gzip
import os
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

import pysam

from .base import BaseFile
from ...Range import CRange, TRange

"""
BED format (only the first three columns are mandatory):
1   chrom
2   start (0-base, inclusive)
3   end (0-base, exclusive)
4   name
5   score
6   strand
7   thick_start (0-base)
8   thick_end (0-base, not included)
9   item_rgb
10  block_count
11  block_sizes
12  block_starts
"""


class BedRecord(TRange):
    """
    BedRecord is a class for a single record (line) of BED files.

    It is a :py:class:`TRange` carrying the BED columns: the optional
    name, score and strand columns and, for BED12 records, the thick
    interval (``thick``) and the block structure (``blocks``). The
    itemRgb column is kept in the ``color`` attribute.
    """

    pass


class BedFile(BaseFile):
    """
    BedFile is a class for reading and writing BED files.

    Plain text BED, gzip-compressed BED (``.bed.gz``) and block-
    compressed BED with a tabix index (``.bed.gz`` + ``.tbi``) are
    supported. For indexed files random access via :py:meth:`fetch`
    is available; otherwise records can only be read sequentially
    (iteration). Each line is parsed into a :py:class:`BedRecord` by
    :py:meth:`parse_bed_string`; records are written by :py:meth:`write`.
    """

    # Attributes assigned by the untyped BaseFile base class; declared
    # here to keep their types known.
    _path: str
    _mode: str
    _handle: "TextIO | pysam.TabixFile | None"

    _ncol: int
    _random: bool

    def __init__(self, path: str, mode: str = "r", ncol: int = 12,
                 random: "bool | None" = None) -> None:
        """
        Initialize and open a BED file object.

        :param path: Path to the BED (or compressed .bed.gz) file.
        :param mode: File mode; "r" for reading and "w" for writing.
        :param ncol: Number of BED columns to read/write, 3 to 12.
        :param random: Enable random access. If None (default), it is
                       enabled automatically for indexed .bed.gz files.
        :raises ValueError: If the mode is unsupported.
        :raises AssertionError: If the random-access preconditions
                                are not met.
        """
        # Enable random access automatically for indexed files.
        if random is None:
            random = mode == "r" and path.endswith(".bed.gz") \
                and os.path.exists(path + ".tbi")
        if random:
            assert mode == "r", "Random access is only supported in read mode!"
            assert path.endswith(".bed.gz"), \
                "Random access requires a .bed.gz file!"
            assert os.path.exists(path), "File %s does not exist!" % path
            assert os.path.exists(path + ".tbi"), \
                "The tabix index (.tbi) file is missing!"
        super(BedFile, self).__init__(path, mode)
        self._ncol = ncol
        self._random = bool(random)
        self.open()

    @property
    def ncol(self) -> int:
        """Number of BED columns to read/write."""
        return self._ncol

    @property
    def random(self) -> bool:
        """Return True if random access is enabled for this file."""
        return self._random

    @override
    def open(self) -> None:
        """
        Open the underlying file handle.

        Random-access files are opened with pysam (tabix); plain and
        gzip-compressed files are opened as text streams. Calling open
        on an already-opened file is a no-op.

        :raises ValueError: If the mode is unsupported.
        """
        if self._handle is None:
            if self.random:
                self._handle = pysam.TabixFile(self.path)
            elif self.mode == "r":
                if self.path.endswith(".gz"):
                    self._handle = cast("TextIO", gzip.open(self.path, "rt"))
                else:
                    self._handle = open(self.path)
            elif self.mode == "w":
                if self.path.endswith(".gz"):
                    self._handle = cast("TextIO", gzip.open(self.path, "wt"))
                else:
                    self._handle = open(self.path, "w+")
            else:
                raise ValueError("Mode %s is not supported!" % self.mode)

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
    def fetch(self, chrom: "str | None" = None, start: "int | None" = None,
              end: "int | None" = None) -> "Iterator[BedRecord]":
        """
        Fetch records of the given chromosome (region).

        For random-access files records of the given region are
        fetched via tabix, and chrom=None iterates over all contigs;
        for sequential files only chrom=None (the whole file) is
        allowed. Comment lines and empty lines are skipped.

        :param chrom: Chromosome name; None iterates over all contigs.
        :param start: Region start (inclusive); None for the whole chromosome.
        :param end: Region end (exclusive); None for the whole chromosome.
        :raises RuntimeError: If random access is not enabled or the
                              file is not opened.
        :return: Iterator of BedRecord objects.
        """
        if self._random:
            handle = self.handle
            assert handle is not None
            assert isinstance(handle, pysam.TabixFile)
            if chrom is None:
                for contig in sorted(handle.contigs):
                    for line in cast("Iterator[str]", handle.fetch(contig)):
                        yield self.parse_bed_string(line, self._ncol)
            else:
                for line in cast("Iterator[str]", handle.fetch(chrom, start, end)):
                    yield self.parse_bed_string(line, self._ncol)
        else:
            if chrom is not None:
                raise RuntimeError(
                    "Random access is not supported! Please set random=True.")
            handle = self.handle
            assert handle is not None
            for line in cast("TextIO", handle):
                line = line.strip()
                if line.startswith("#") or line == "":
                    continue
                yield self.parse_bed_string(line, self.ncol)

    def write(self, record: CRange) -> None:
        """
        Write a record to the file as a BED line.

        The record is formatted with the first ``ncol`` BED columns.

        :param record: The range record to write.
        :raises RuntimeError: If the file is not opened.
        """
        handle = self.handle
        assert handle is not None
        assert not isinstance(handle, pysam.TabixFile)
        _ = handle.write(record.to_bed_string(ncol=self.ncol) + "\n")

    @classmethod
    def parse_bed_string(cls, line: str, column: int = 12) -> BedRecord:
        """
        Parse a single BED line into a BedRecord.

        :param line: One line of a BED file.
        :param column: Maximum number of BED columns to consider,
                       3 to 12.
        :return: The parsed BedRecord instance.
        :raises AssertionError: If the line has fewer than 3 columns.
        """
        row = line.split("\t")
        column = min(len(row), column)
        assert column >= 3

        chrom = row[0]
        start = int(row[1])
        end = int(row[2])

        name = row[3] if column >= 4 else None
        strand = row[5] if column >= 6 else None

        thick: "tuple[int, int] | None" = None
        if column >= 8:
            thick_start = int(row[6])
            thick_end = int(row[7])
            if thick_start < thick_end:
                thick = (thick_start, thick_end)

        blocks: "list[tuple[int, int]] | None" = None
        if column >= 12:
            block_sizes = [int(s) for s in row[10].split(",") if s != ""]
            block_offsets = [int(s) for s in row[11].split(",") if s != ""]
            blocks = [(start + offset, start + offset + size)
                      for offset, size in zip(block_offsets, block_sizes)]

        record = BedRecord(chrom=chrom, start=start, end=end, name=name,
                           strand=strand, thick=thick, blocks=blocks)

        if column >= 5 and row[4] != ".":
            record.score = float(row[4])
        if column >= 9:
            record.color = row[8]

        return record

    @classmethod
    def to_bed_string(cls, record: CRange, column: int = 12) -> str:
        """
        Format a record as a BED line with the first ``column`` columns.

        A thin wrapper around the record's own ``to_bed_string``
        method, kept for backward compatibility.

        :param record: The range record to format.
        :param column: Number of BED columns, 3 to 12.
        :return: Tab-separated BED text line.
        """
        return record.to_bed_string(ncol=column)
