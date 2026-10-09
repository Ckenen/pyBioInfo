# Class: TRange (thick range)

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

from .genomic_range import GRange


class TRange(GRange):
    """
    TRange is a class for genomic range with a thick interval.

    It extends GRange with an optional **thick** interval marking the
    coding part of the range (like the thick columns of a BED12
    record, e.g. the CDS interval of a transcript). The thick interval
    must lie inside the blocks: its start and its last position
    (end - 1) must both be covered by some block.

    Examples:
        >>> obj = TRange("chr1", blocks=[(10, 20), (40, 50)], thick=(40, 50))
        >>> obj.thick
        (40, 50)
    """

    _thick: tuple[int, int] | None
    _cache_forward_indexes: tuple[int, 'int | None', 'int | None', int] | None
    _cache_reverse_indexes: tuple[int, 'int | None', 'int | None', int] | None

    def __init__(self, chrom: str, start: "int | None" = None, end: "int | None" = None,
                 name: "str | None" = None, strand: "str | None" = None,
                 blocks: "Iterable[Iterable[int]] | None" = None,
                 thick: "tuple[int, int] | None" = None) -> None:
        """
        Initialize a thick range.

        :param chrom: Chromosome name, e.g. "chr1".
        :param start: Overall start (inclusive); ignored if blocks given.
        :param end: Overall end (exclusive); defaults to start + 1 when
                    omitted (with start given); ignored if blocks given.
        :param name: Optional feature name.
        :param strand: Optional strand: None, ".", "+" or "-".
        :param blocks: Optional iterable of (start, end) block intervals.
        :param thick: Optional thick interval (start, end), must be
                      included in the blocks.
        :raises ValueError: If coordinates are inconsistent or invalid,
                            or if the thick interval is out of the blocks.
        """
        super().__init__(chrom=chrom, start=start, end=end, name=name, strand=strand, blocks=blocks)

        self._thick = None

        if thick is not None:
            thick_start, thick_end = thick
            if thick_start < self.start:
                raise ValueError("The value of thick start[%d] can not smaller than start[%d]." % (thick_start, self.start))
            if thick_end > self.end:
                raise ValueError("The value of thick end[%d] can not larger than end[%d]." % (thick_end, self.end))
            if not BlockTools.is_include(self.blocks, thick_start):
                raise ValueError("The value of thick start[%d] is not included in blocks." % thick_start)
            if not BlockTools.is_include(self.blocks, thick_end - 1):
                raise ValueError("The value of thick end[%d] is not included in blocks." % thick_end)
            self._thick = (thick_start, thick_end)

        # Cached index ranges of the thick interval (see self.indexes()).
        self._cache_forward_indexes = None
        self._cache_reverse_indexes = None

    # Access to protected attributes.

    @property
    def thick(self) -> "tuple[int, int] | None":
        """Thick interval (start, end) of the range, or None."""
        return self._thick

    # string format and thick indexes.

    @override
    def __str__(self) -> str:
        """Return a string with span, strand, block count and thick interval."""
        text = "%s: %d-%d [strand: %s, blocks: %d" % (
            self.__class__.__name__,
            self._start, self._end,
            self.strand, len(self._blocks))
        if self.thick is not None:
            text += ", thick: %d-%d" % (self.thick[0], self.thick[1])
        return text + "]"

    def indexes(self, strandness: bool = True) -> "tuple[int, 'int | None', 'int | None', int]":
        """
        Return the (start, thick_start, thick_end, end) index tuple.

        The indexes are relative to the concatenated blocks. On the
        reverse strand (with strandness enabled) they run from the end
        towards the start. The thick indexes are None if there is no
        thick interval. Results are cached per strand.

        :param strandness: Use reverse-strand indexes if on reverse strand.
        :return: Tuple of (0, thick_start_index, thick_end_index, length).
        """

        i1 = 0
        i2: "int | None" = None
        i3: "int | None" = None
        i4 = len(self)

        thick = self.thick

        if thick is None:
            return (i1, i2, i3, i4)

        if self.is_reverse and strandness:
            if self._cache_reverse_indexes is None:
                thick_start, thick_end = thick
                i2 = self.get_index(thick_start)
                i3 = self.get_index(thick_end - 1)
                i2, i3 = min(i2, i3), max(i2, i3) + 1
                self._cache_reverse_indexes = (i1, i2, i3, i4)
            indexes = self._cache_reverse_indexes
        else:
            if self._cache_forward_indexes is None:
                thick_start, thick_end = thick
                i2 = self.get_index(thick_start, strandness=False)
                i3 = self.get_index(thick_end - 1, strandness=False)
                i2, i3 = min(i2, i3), max(i2, i3) + 1
                self._cache_forward_indexes = (i1, i2, i3, i4)
            indexes = self._cache_forward_indexes
        return indexes
