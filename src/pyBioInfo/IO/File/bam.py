# Class: BamFile (reader and writer of BAM format files)

from os import stat
import logging
import os
import sys
from collections.abc import Iterator
from typing import Any, cast

if sys.version_info >= (3, 12):
    from typing import override
else:
    from typing import Callable  # pyright: ignore[reportUnreachable]

    def override(func: Callable) -> Callable:
        """No-op fallback for Python < 3.12 compatibility."""
        return func

import pysam

from .base import BaseFile
from ...Range import GRange
from ...Utils import SegmentTools

# Note: the parsed-CIGAR type ("list[tuple[str, tuple[int, int],
# tuple[int, int], tuple[int, int]]]") and the parsed-MD-tag type
# ("list[tuple[str, tuple[int, int], str | None]]") are spelled out
# inline in the annotations below.


def _header_to_dict(header: "dict[str, Any] | pysam.AlignmentHeader") -> "dict[str, Any]":  # pyright: ignore[reportExplicitAny]
    """
    Copy a BAM header (dict or pysam AlignmentHeader) into a plain dict.

    :param header: A header dict or a pysam AlignmentHeader object.
    :return: A copied header dict.
    """
    if isinstance(header, dict):
        return header.copy()
    # AlignmentHeader.as_dict exists at runtime but is missing from the
    # bundled pysam stubs, hence the cast.
    return cast("dict[str, Any]", header.as_dict())  # pyright: ignore[reportAttributeAccessIssue, reportUnknownMemberType, reportExplicitAny]


class Alignment(GRange):
    """
    Alignment is a class for a single alignment record of BAM files.

    It is a :py:class:`GRange` (blocked genomic range) wrapping a
    ``pysam.AlignedSegment``: the alignment blocks are copied into the
    range while the segment itself is kept in the ``segment`` property.
    Alignment-derived data (aligned sequence, parsed CIGAR and MD tag)
    are computed lazily on first access and cached.
    """

    def __init__(self, segment: pysam.AlignedSegment) -> None:
        """
        Build an Alignment from a pysam segment.

        :param segment: An aligned segment of a BAM file.
        """
        chrom = segment.reference_name
        if chrom is None:
            raise RuntimeError("Segment is not mapped to a reference!")
        super(Alignment, self).__init__(
            chrom=chrom,
            name=segment.query_name,
            blocks=SegmentTools.get_blocks(segment),
            strand="-" if segment.is_reverse else "+")
        self._segment: pysam.AlignedSegment = segment
        self._aligned_sequences: str | None = None
        self._aligned_qualities: list[int] | None = None
        self._parsed_cigar: list[tuple[str, tuple[int, int], tuple[int, int], tuple[int, int]]] | None = None
        self._parsed_md_tag: list[tuple[str, tuple[int, int], str | None]] | None = None

    @property
    def segment(self) -> pysam.AlignedSegment:
        """The wrapped pysam segment."""
        return self._segment

    def get_aligned_sequences(self) -> "tuple[str, list[int]]":
        """
        Return the read sequence and qualities aligned to the reference.

        Deletion positions are filled with "-" and quality -1;
        insertions and clips are removed. The result is computed once
        and cached.

        :return: (aligned_sequence, aligned_qualities) tuple.
        """
        if self._aligned_sequences is None or self._aligned_qualities is None:
            self._aligned_sequences, self._aligned_qualities = \
                SegmentTools.get_aligned_sequence(self.segment)
        return self._aligned_sequences, self._aligned_qualities

    def get_parsed_cigar(self) -> "list[tuple[str, tuple[int, int], tuple[int, int], tuple[int, int]]]":
        """
        Return the parsed CIGAR string of the segment.

        The result is computed once and cached.

        :return: List of (operation, mapped_indexes, read_indexes,
                 reference_indexes) tuples.
        """
        if self._parsed_cigar is None:
            self._parsed_cigar = SegmentTools.parse_cigar(self.segment)
        return self._parsed_cigar

    def get_parsed_md_tag(self) -> "list[tuple[str, tuple[int, int], str | None]]":
        """
        Return the parsed MD tag of the segment.

        The result is computed once and cached.

        :return: List of (kind, mapped_indexes, bases) tuples.
        """
        if self._parsed_md_tag is None:
            self._parsed_md_tag = SegmentTools.parse_md_tag(self.segment)
        return self._parsed_md_tag

    def get_pared_md_tag(self) -> "list[tuple[str, tuple[int, int], str | None]]":
        """Deprecated alias of :py:meth:`get_parsed_md_tag`."""
        return self.get_parsed_md_tag()

    def get_query_base(self, position: int) -> "str | None":
        """
        Return the read base aligned at a reference position.

        :param position: Reference position (0-based).
        :return: The aligned base, or None if the position is not
                 aligned to a read base (deletion or outside blocks).
        """
        return SegmentTools.get_query_base(self.segment, position,
                                           self.get_parsed_cigar())


class BamFile(BaseFile):
    """
    BamFile is a class for reading and writing BAM files.

    Reading (mode "rb") iterates over the mapped records, either
    sequentially or, for indexed files, via random access with
    :py:meth:`fetch`. Writing (mode "wb") requires a header, given
    directly or via a ``template`` (a BamFile or pysam.AlignmentFile).
    Each record is exposed as an :py:class:`Alignment`.
    """

    # Attributes assigned by the untyped BaseFile base class; declared
    # here to keep their types known.
    _path: str
    _mode: str
    _handle: "pysam.AlignmentFile | None"
    _references: dict[str, int] | None
    _header: dict[str, Any] | None # pyright: ignore[reportExplicitAny]
    _random: bool

    def __init__(self, path: str, mode: str = "rb",
                 template: "BamFile | pysam.AlignmentFile | None" = None,
                 header: "dict[str, Any] | pysam.AlignmentHeader | None" = None,  # pyright: ignore[reportExplicitAny]
                 random: "bool | None" = None) -> None:
        """
        Initialize and open a BAM file object.

        :param path: Path to the BAM file.
        :param mode: File mode; "rb" for reading and "wb" for writing.
        :param template: Template providing the header for writing.
        :param header: Header dict (or pysam AlignmentHeader) for writing.
        :param random: Enable random access. If None (default), it is
                       enabled automatically for indexed .bam files.
        :raises ValueError: If the mode is unsupported.
        :raises AssertionError: If the random-access preconditions
                                are not met.
        """
        # Enable random access automatically for indexed files.
        if mode[0] == "r" and random is None:
            if path.endswith(".bam") and os.path.exists(path + ".bai"):
                logging.warning("Index file exists, automatic set random=True")
                random = True
            else:
                logging.warning(
                    "Index file does not exists, automatic set random=False")
                random = False

        if random:
            if not path.endswith(".bam"):
                raise RuntimeError("Random access requires a .bam file!")
            if not os.path.exists(path):
                raise RuntimeError("File %s does not exist!" % path)
            if not os.path.exists(path + ".bai"):
                raise RuntimeError("The index (.bai) file is missing!")
        
        if mode not in ("rb", "wb"):
            raise ValueError("Mode %s is not supported!" % mode)

        self._template: "BamFile | pysam.AlignmentFile | None" = template
        header_dict = _header_to_dict(header) if header is not None else None
        self._header = header_dict  
        self._references = None
        self._random = bool(random)
        super(BamFile, self).__init__(path, mode)
        self.open()

    @property
    def header(self) -> "dict[str, Any] | None":  # pyright: ignore[reportExplicitAny]
        """Header of the BAM file as a dict, or None if unavailable."""
        return self._header

    @property
    def references(self) -> "dict[str, int] | None":
        """Mapping of reference names to lengths, or None if unopened."""
        return self._references

    @property
    def random(self) -> bool:
        """Return True if random access is enabled for this file."""
        return self._random

    @override
    def open(self) -> None:
        """
        Open the underlying file handle.

        For reading, the header is loaded and exposed by the ``header``
        and ``references`` properties. For writing, the header is taken
        from the ``header`` argument or from the ``template``. Calling
        open on an already-opened file is a no-op.

        :raises RuntimeError: If neither header nor template is given
                              for writing.
        """
        if self._handle is None:
            if self._mode == "rb":
                self._handle = pysam.AlignmentFile(self._path)
                self._header = _header_to_dict(self._handle.header)
            elif self._mode == "wb":
                header = self._header
                if header is None:
                    template = self._template
                    if isinstance(template, BamFile):
                        handle = template._handle
                        assert handle is not None, "Template file is not opened!"
                        header = _header_to_dict(handle.header)
                    elif isinstance(template, pysam.AlignmentFile):
                        header = _header_to_dict(template.header)
                    else:
                        raise RuntimeError(
                            "Writing a BAM file requires a header or a template!")
                self._handle = pysam.AlignmentFile(self._path, "wb", header=header)
                self._header = header
            else:
                raise RuntimeError("Mode %s is not supported!" % self._mode)
            handle = self._handle
            header = self._header
            assert handle is not None and header is not None
            references: "dict[str, int]" = dict()
            for item in header["SQ"]:  # pyright: ignore[reportAny]
                references[item["SN"]] = item["LN"]
            self._references = references

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
              end: "int | None" = None) -> "Iterator[Alignment]":
        """
        Fetch the alignment records of the given chromosome (region).

        For random-access files records of the given region are fetched
        via the index, and chrom=None iterates over all references; for
        sequential files only chrom=None (the whole file) is allowed.
        Unmapped records are skipped.

        :param chrom: Chromosome name; None iterates over all references.
        :param start: Region start (inclusive); None for the whole chromosome.
        :param end: Region end (exclusive); None for the whole chromosome.
        :raises RuntimeError: If random access is not enabled or the
                              file is not opened.
        :return: Iterator of Alignment objects.
        """
        handle = self._handle
        assert handle is not None, "File is not opened!"
        if self._random:
            references = self._references
            assert references is not None, "File is not opened!"
            if chrom is None:
                for contig in sorted(references):
                    for segment in handle.fetch(contig=contig):
                        if segment.is_unmapped:
                            continue
                        yield Alignment(segment)
            else:
                for segment in handle.fetch(contig=chrom, start=start, stop=end):
                    if segment.is_unmapped:
                        continue
                    yield Alignment(segment)
        else:
            if chrom is not None:
                raise RuntimeError(
                    "Random access is not supported! Please set random=True.")
            for segment in handle.fetch(until_eof=True):
                if segment.is_unmapped:
                    continue
                yield Alignment(segment)

    def write(self, obj: object) -> None:
        """
        Write an alignment record to the file.

        :param obj: An Alignment (its wrapped segment is written) or a
                    raw pysam AlignedSegment.
        :raises TypeError: If the record type is unsupported.
        """
        handle = self._handle
        assert handle is not None, "File is not opened!"
        if isinstance(obj, Alignment):
            _ = handle.write(obj.segment)
        elif isinstance(obj, pysam.AlignedSegment):
            _ = handle.write(obj)
        else:
            raise TypeError("Unsupported record type: %s!" % type(obj))

    @staticmethod
    def detect_layout(bamfile: str) -> str:
        """
        Detect the layout of a BAM file.
        """
        with pysam.AlignmentFile(bamfile) as f:
            for segment in f.fetch(until_eof=True):
                if segment.is_paired:
                    return "PE"
                else:
                    return "SE"
        return "SE"