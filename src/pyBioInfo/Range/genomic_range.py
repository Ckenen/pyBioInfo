# Class: GRange (genomic range)

import sys
from collections.abc import Iterable

if sys.version_info >= (3, 12):
    from typing import override
else:
    from typing import Callable  # pyright: ignore[reportUnreachable]

    def override(func: Callable) -> Callable:
        """No-op fallback for Python < 3.12 compatibility."""
        return func

# Re-exported by the package (see pyBioInfo.Utils.__all__).
from pyBioInfo.Utils import BlockTools

from .chromosome_range import CRange


class GRange(CRange):
    """
    GRange is a class for genomic range with multiple blocks.

    A genomic range is a half-open interval [start, end) located on a
    chromosome, carrying an optional feature name and strand, whose
    interval may be split into several disjoint "blocks" (e.g. exons
    of a transcript). The blocks are stored as a sorted tuple of
    half-open intervals [start, end); the overall range spans from the
    first block start to the last block end, while len() returns the
    total length of all blocks (gaps between blocks excluded).

    The read-only ``start``/``end`` span coordinates, the chromosome
    name ``chrom``, the validated ``strand`` (with ``is_forward`` /
    ``is_reverse``), the 2-element (start, end) sequence protocol
    (iteration and [0]/[1] indexing), the total ordering and the
    distance relationships (overlap / adjacency / contact, span-based
    and chromosome-aware) are inherited from the base classes (IRange
    and CRange). On top of that, GRange provides the block-level
    containment tests (``is_contain`` / ``is_coincide``), while
    ``__str__`` and ``__len__`` are overridden to show the block count
    and to return the total block length (gaps excluded) instead of
    the span length.

    Strand-aware index/position conversions follow the block layout:
    indexes run over the concatenated blocks, skipping the gaps. On
    the reverse strand the internal indexing runs from the end
    towards the start, like a sequence read on the minus strand.

    Ordering rule: ranges are grouped by chromosome first (compared
    alphabetically), then ordered by start and end coordinates within
    the same chromosome; the ordering implementation is inherited from
    CRange (chromosome-aware) and IRange (coordinates).

    Examples:
        >>> obj = GRange("chr1", blocks=[(10, 15), (20, 25)])
        >>> obj.start, obj.end       # overall span
        (10, 25)
        >>> len(obj)                 # blocks only, gap excluded
        10
        >>> start, end = obj         # iterable as (start, end)
        (10, 25)
    """

    # Optional annotation columns used when formatting; subclasses may
    # override them. None means "not available" and is rendered as ".".

    # Common options
    score: "int | float | str | None" = None

    # GFF and GTF specific
    source: "str | None" = None
    feature: "str | None" = None
    frame: "int | str | None" = None
    
    # BED specific
    color: "str | None" = None


    def __init__(self, chrom: str, start: "int | None" = None, end: "int | None" = None,
                 name: "str | None" = None, strand: "str | None" = None,
                 blocks: "Iterable[Iterable[int]] | None" = None) -> None:
        """
        Initialize a genomic range.

        Either pass ``blocks`` (start/end are then optional and must be
        consistent with them), or pass ``start``/``end`` to create a
        single-block range.

        :param chrom: Chromosome name, e.g. "chr1".
        :param start: Overall start (inclusive); ignored if blocks given.
        :param end: Overall end (exclusive); defaults to start + 1 when
                    omitted (with start given); ignored if blocks given.
        :param name: Optional feature name.
        :param strand: Optional strand: None, ".", "+" or "-".
        :param blocks: Optional iterable of (start, end) block intervals.
        :raises ValueError: If coordinates are inconsistent or invalid.
        """

        blocks_tuple: "tuple[tuple[int, int], ...]"
        if blocks is None:
            assert start is not None
            if end is None:
                end = start + 1
            blocks_tuple = ((start, end),)
        else:
            blocks_tuple = tuple((block_start, block_end) for block_start, block_end in blocks)
            if start is not None and start != blocks_tuple[0][0]:
                raise ValueError("Inconsistent value of start and blocks.")
            if end is not None and end != blocks_tuple[-1][1]:
                raise ValueError("Inconsistent value of end and blocks.")

        # Validate the blocks (start >= 0, start < end, sorted and
        # non-overlapping) in both construction paths, then set the
        # chromosome, the overall span coordinates and the strand via
        # the CRange constructor.
        BlockTools.check_blocks(blocks_tuple)
        super().__init__(chrom, blocks_tuple[0][0], blocks_tuple[-1][1], strand=strand)

        self.name: "str | None" = name
        self._blocks: "tuple[tuple[int, int], ...]" = blocks_tuple

        # Cache the block total length so that len() is O(1).
        self._cache_length: int = sum(y - x for x, y in blocks_tuple)


    # Access to protected attributes.

    @property
    def blocks(self) -> "tuple[tuple[int, int], ...]":
        """Blocks of the range as a sorted tuple of (start, end) intervals."""
        return self._blocks

    # string format and length.

    @override
    def __str__(self) -> str:
        """Return a string with span, strand and block count."""
        return "%s: %d-%d [strand: %s, blocks: %d]" % (
            self.__class__.__name__,
            self._start, self._end,
            self.strand, len(self._blocks))

    @override
    def __len__(self) -> int:
        """Return the total length of all blocks (gaps excluded)."""
        return self._cache_length

    # index-position relevant operations.

    @override
    def get_position(self, index: int, strandness: bool = True) -> int:
        """
        Convert a block-aware 0-based index to an absolute coordinate.

        Indexes run over the concatenated blocks, skipping the gaps. On
        the reverse strand (with strandness enabled) index 0 maps to the
        last position of the last block.

        :param index: 0-based index, negative values count from the end.
        :param strandness: Apply strand-aware mapping if on reverse strand.
        :return: Absolute coordinate corresponding to the index.
        :raises ValueError: If the index lies outside [-length, length).
        """
        length = len(self)
        if -length <= index < length:
            if strandness and self.is_reverse:
                if index >= 0:
                    index = length - 1 - index
                else:
                    index = -length - 1 - index
            return BlockTools.get_position(self.blocks, index, length=length)
        else:
            raise ValueError("Index (%d) is out of range [%d-%d]." % (index, -length, length))

    @override
    def get_index(self, position: int, strandness: bool = True) -> int:
        """
        Convert an absolute coordinate to a block-aware 0-based index.

        Indexes run over the concatenated blocks, skipping the gaps. On
        the reverse strand (with strandness enabled) the index counts
        from the end towards the start.

        :param position: Absolute coordinate on the chromosome.
        :param strandness: Apply strand-aware indexing if on reverse strand.
        :return: 0-based index of the position inside the blocks.
        :raises ValueError: If the position lies outside the blocks.
        """
        return BlockTools.get_index(self.blocks, position, reverse=self.is_reverse and strandness)

    # Distance relationships are inherited from the base classes:
    # IRange provides the coordinate-only distance (is_overlap /
    # is_adjacent / is_away / is_contact on the outer span) and CRange
    # makes it chromosome-aware.

    def is_contain(self, other: "GRange") -> bool:
        """
        Return True if this range's blocks fully cover `other`'s blocks.

        :param other: Another GRange instance.
        """
        return BlockTools.is_contain(self.blocks, other.blocks)

    def is_coincide(self, other: "GRange") -> bool:
        """
        Return True if `other`'s blocks coincide with this range's blocks.

        Block boundaries must align, allowing `other` to be trimmed at
        the first block start and/or the last block end.

        :param other: Another GRange instance.
        """
        return BlockTools.is_coincide(self.blocks, other.blocks)

    # def clip(self, start, end):
    #     blocks = BlockTools.clip(blocks=self.blocks, start=start, end=end, extend=False)
    #     obj = None
    #     if len(blocks) > 0:
    #         obj = GRange(chrom=self.chrom, blocks=blocks, name=self.name, strand=self.strand)
    #     return obj
