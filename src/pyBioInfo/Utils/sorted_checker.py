# SortedChecker: an iterator wrapper asserting the input is sorted.

from collections.abc import Generator, Iterable
from typing import Protocol

class _SortedItem(Protocol):
    """Structural type of the items checked by SortedChecker."""

    chrom: str
    start: int
    end: int


class SortedChecker(object):
    """
    Iterate over a container while asserting genomic sort order.

    Items must provide ``chrom``, ``start`` and ``end`` attributes and
    be sorted by chrom (lexicographic), then start, then end. As each
    item is yielded, it is compared against its predecessor and a
    ``RuntimeError`` is raised on the first ordering violation.
    """

    # Iterable of items to check.
    _container: "Iterable[_SortedItem]"
    # The previously yielded item (None before the first item).
    _last: "_SortedItem | None"
    # Iterator created lazily by __next__.
    _iterator: "Generator[_SortedItem, None, None] | None"
    # Whether chrom ordering is checked.
    _check_chrom: bool
    # Whether start ordering is checked.
    _check_start: bool
    # Whether end ordering is checked.
    _check_end: bool

    def __init__(self, container: "Iterable[_SortedItem]", check_chrom: bool = True,
                 check_start: bool = True, check_end: bool = True) -> None:
        """
        Initialize a SortedChecker.

        :param container: Iterable of items sorted by chrom, start, end.
        :param check_chrom: Whether to check chrom ordering.
        :param check_start: Whether to check start ordering.
        :param check_end: Whether to check end ordering.
        :raises AssertionError: If chrom or start checking is disabled.
        """
        self._container = container
        self._last = None
        self._iterator = None
        self._check_chrom = check_chrom
        self._check_start = check_start
        self._check_end = check_end
        assert self._check_chrom
        assert self._check_start

    def __iter__(self) -> "Generator[_SortedItem, None, None]":
        """
        Yield items one by one, asserting sort order against the
        previous item.

        :raises RuntimeError: If the chrom/start/end ordering is
                              violated between two consecutive items.
        """
        for current in self._container:
            if self._last is not None and self._check_chrom:
                if self._last.chrom > current.chrom:
                    raise RuntimeError(
                        "Error in sorted chrom [last: %s, current: %s]" % (self._last.chrom, current.chrom))
                elif self._last.chrom == current.chrom and self._check_start:
                    if self._last.start > current.start:
                        raise RuntimeError(
                            "Error in sorted start [chrom: %s, last: %d, current: %d]" % (
                                current.chrom, self._last.start, current.start))
                    elif self._last.start == current.start and self._check_end:
                        if self._last.end > current.end:
                            raise RuntimeError(
                                "Error in sorted end [chrom: %s, start: %d, last: %d, current: %d]" % (
                                    current.chrom, current.start, self._last.end, current.end))

                # if current < self._last:
                #     raise RuntimeError("The items is not sorted! [last: %s, current: %s]" % (self._last, current))
            yield current
            self._last = current

    def __next__(self) -> "_SortedItem":
        """Advance the lazily created iterator and return the next item."""
        if self._iterator is None:
            self._iterator = self.__iter__()
        return next(self._iterator)
