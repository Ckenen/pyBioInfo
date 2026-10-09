# Class: VcfFile (reader and writer of VCF format files)

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
from ...Range import GRange


class VcfRecord(GRange):
    """
    VcfRecord is a class for a single record (line) of VCF files.

    It is a :py:class:`GRange` whose single block spans the reference
    allele, carrying the site ID, quality, filter status, alternate
    alleles and the INFO column string.
    """

    def __init__(self, chrom: str, start: int, end: int, strand: str,
                 id: "str | None", qual: "float | None",
                 filter: str, ref: str, alts: "list[str]",
                 info: "str | None") -> None:
        """
        Initialize a VCF record.

        :param chrom: Chromosome name.
        :param start: Start coordinate (0-based, inclusive).
        :param end: End coordinate (0-based, exclusive).
        :param strand: Strand: "+" (VCF is strandless, "+" by convention).
        :param id: Site identifier, or None for ".".
        :param qual: Quality score, or None.
        :param filter: Filter status ("PASS", "." or a list of filters).
        :param ref: Reference allele string.
        :param alts: List of alternate allele strings.
        :param info: INFO column string, or None.
        """
        super(VcfRecord, self).__init__(chrom=chrom,
                                        start=start,
                                        end=end,
                                        strand=strand,
                                        blocks=[(start, end)])
        self.id: "str | None" = id
        self.qual: "float | None" = qual
        self.filter: str = filter
        self.ref: str = ref
        self.alts: "list[str]" = alts
        self.info: "str | None" = info

    def induce(self, ref: str, alt: str) -> None:
        """
        Set the reference and alternate alleles of the record.

        :param ref: Reference allele string.
        :param alt: Alternate allele string.
        :raises AssertionError: If the alleles do not match the record.
        """
        assert ref == self.ref
        assert alt in self.alts
        self.ref = ref
        self.alts = [alt]


class VcfFile(BaseFile):
    """
    VcfFile is a class for reading and writing VCF files.

    Plain text VCF and gzip-compressed VCF (``.vcf.gz``) are supported.
    The header comment lines ("#...") are collected into
    :py:attr:`comments` while reading; each data line is parsed into a
    :py:class:`VcfRecord`.
    """

    # Attributes assigned by the untyped BaseFile base class; declared
    # here to keep their types known.
    _path: str
    _mode: str
    _handle: "TextIO | None"

    def __init__(self, path: str, mode: str = "r",
                 template: "object | None" = None,
                 header: "object | None" = None) -> None:
        """
        Initialize a VCF file object.

        The file is opened lazily: the underlying handle is created by
        :py:meth:`open`, either explicitly or on first access via
        :py:attr:`handle`.

        :param path: Path to the VCF (or compressed .vcf.gz) file.
        :param mode: File mode; "r" for reading and "w" for writing.
        :param template: Not used; kept for signature compatibility.
        :param header: Not used; kept for signature compatibility.
        :raises AssertionError: If the mode is invalid.
        """
        assert mode == "r" or mode == "w", "Mode %s is not supported!" % mode
        super(VcfFile, self).__init__(path, mode)
        self._template: "object | None" = template
        self._header: "object | None" = header
        # Comment lines ("#...") collected while reading.
        self._comments: "list[str]" = []

    @property
    def comments(self) -> "list[str]":
        """Comment lines ("#...") collected while reading."""
        return self._comments

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
    def fetch(self) -> "Iterator[VcfRecord]":
        """
        Fetch the records of the file sequentially.

        Comment lines are collected into :py:attr:`comments` and empty
        lines are skipped.

        :raises RuntimeError: If the file is not opened.
        :return: Iterator of VcfRecord objects.
        """
        handle = self._handle
        assert handle is not None, "File is not opened!"
        for line in handle:
            if line.startswith("#"):
                self._comments.append(line)
                continue
            line = line.strip("\n")
            if line == "":
                continue
            yield self.parse_vcf_string(line)

    @classmethod
    def parse_vcf_string(cls, line: str) -> VcfRecord:
        """
        Parse a single VCF line into a VcfRecord.

        :param line: One line of a VCF file.
        :return: The parsed VcfRecord instance.
        :raises AssertionError: If the line is malformed.
        """
        row = line.split("\t")
        assert len(row) == 10, "Invalid VCF record: %s" % line
        chrom = row[0]
        start = int(row[1]) - 1
        id = row[2] if row[2] != "." else None
        ref = row[3]
        end = start + len(ref)
        alts = row[4].split(",")
        qual = float(row[5]) if row[5] != "." else None
        filter = row[6]
        info = row[7] if row[7] != "." else None
        record = VcfRecord(chrom=chrom,
                           start=start,
                           end=end,
                           strand="+",
                           id=id,
                           qual=qual,
                           filter=filter,
                           ref=ref,
                           alts=alts,
                           info=info)
        return record

    def write(self, record: VcfRecord) -> None:
        """
        Write a VCF record to the file.

        The header comments collected while reading (or the ones added
        manually to :py:attr:`comments`) are written before the first
        record.

        :param record: The VcfRecord to write.
        :raises RuntimeError: If the file is not opened.
        """
        if self._mode != "w":
            raise RuntimeError("File is not opened in write mode!")
        self.open()
        handle = self._handle
        assert handle is not None, "File is not opened!"
        if len(self._comments) > 0:
            comments, self._comments = self._comments, []
            for comment in comments:
                _ = handle.write(comment)
        columns = [record.chrom,
                   str(record.start + 1),
                   record.id if record.id is not None else ".",
                   record.ref,
                   ",".join(record.alts),
                   str(record.qual) if record.qual is not None else ".",
                   record.filter,
                   record.info if record.info is not None else ".",
                   "GT",
                   "0|0"]
        _ = handle.write("\t".join(columns) + "\n")
