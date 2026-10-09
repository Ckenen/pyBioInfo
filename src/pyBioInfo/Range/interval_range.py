# Class: IRange (interval range)

import sys
from collections.abc import Iterator

if sys.version_info >= (3, 12):
    from typing import override
else:
    from typing import Callable  # pyright: ignore[reportUnreachable]

    def override(func: Callable) -> Callable:
        """No-op fallback for Python < 3.12 compatibility."""
        return func


class IRange(object):
    """
    IRange is a minimal class for a one-dimensional interval range.

    A half-open interval [start, end) on the integer axis: the start
    coordinate is inclusive while the end coordinate is exclusive, so
    the interval always has a positive length of (end - start). The
    class exposes the two read-only attributes ``start`` and ``end``
    plus the index/position conversions between absolute coordinates
    and 0-based offsets within the interval; it is the plain-coordinate
    base of the Range hierarchy.

    Instances also implement a fixed 2-element sequence protocol:
    ``len()`` returns the interval length, iteration and ``[0]``/``[1]``
    indexing yield the start and end coordinates (item assignment is
    forbidden), and ``str()`` renders as ``"<ClassName>: <start>-<end>"``.
    Total ordering is provided by ``_compare`` (by ``start`` then
    ``end``), which subclasses extend with higher-level keys, and the
    distance relationships (``is_overlap`` / ``is_adjacent`` /
    ``is_away`` / ``is_contact``) are decided by ``get_distance`` on the
    coordinates alone.

    Examples:
        >>> obj = IRange(5, 15)
        >>> obj.start, obj.end
        (5, 15)
        >>> len(obj)                 # interval length
        10
        >>> start, end = obj         # iterable as (start, end)
        (5, 15)
        >>> str(obj)
        'IRange: 5-15'
        >>> obj.get_index(8)         # coordinate -> 0-based index
        3
        >>> obj.get_position(-1)     # index -> coordinate, negative ok
        14
    """

    # Start and end coordinates.
    _start: int
    _end: int

    def __init__(self, start: int, end: int) -> None:
        """
        Initialize an interval range.

        :param start: Start coordinate (inclusive), must be >= 0.
        :param end: End coordinate (exclusive), must be > start.
        :raises ValueError: If start is negative or start >= end.
        """
        if start < 0:
            raise ValueError("The value of start (%d) must be larger or equal to 0." % start)
        if start >= end:
            raise ValueError("The value of start (%d) must be larger than that of end (%d)." % (start, end))
        self._start = start
        self._end = end

    @property
    def start(self) -> int:
        """Start coordinate (inclusive) of the range."""
        return self._start

    @property
    def end(self) -> int:
        """End coordinate (exclusive) of the range."""
        return self._end

    # Length, text form and sequence protocol.

    @override
    def __str__(self) -> str:
        """Return a string in the form "<ClassName>: <start>-<end>"."""
        return "%s: %d-%d" % (self.__class__.__name__, self.start, self.end)

    def __len__(self) -> int:
        """Return the length of the range, i.e. (end - start)."""
        return self.end - self.start

    def __iter__(self) -> "Iterator[int]":
        """Iterate over the range as a 2-element (start, end) sequence."""
        yield self.start
        yield self.end

    def __getitem__(self, item: int) -> int:
        """
        Return the start (index 0) or the end (index 1) coordinate.

        :param item: Index, 0 for start and 1 for end.
        :return: The requested coordinate.
        :raises RuntimeError: If the index is not 0 or 1.
        """
        if item == 0:
            return self.start
        elif item == 1:
            return self.end
        else:
            raise RuntimeError("Index (%d) is out of range." % item)

    def __setitem__(self, item: int, value: int) -> None:
        """
        Item assignment is not allowed; the coordinates are read-only.

        :param item: Index, 0 for start and 1 for end (unused).
        :param value: The new coordinate (unused).
        :raises RuntimeError: Always, since the range is immutable.
        """
        # Disable the default behaviors.
        raise RuntimeError("Not settable value.")

    # Comparison and ordering (by start and end coordinates only).

    def _compare(self, other: "IRange") -> int:
        """
        Compare this range with another one by coordinates only.

        Ordering rule: by ``start`` first, then by ``end``. Subclasses
        may override this to add higher-level keys (e.g. the chromosome
        name in CRange) before the coordinates.

        :param other: Another IRange instance.
        :return: Negative if self < other, 0 if equal, positive if self > other.
        """
        if self.start < other.start:
            return -1
        elif self.start == other.start:
            if self.end < other.end:
                return -1
            elif self.end == other.end:
                return 0
        return 1

    def __lt__(self, other: "IRange") -> bool:
        """Return True if this range sorts before `other`."""
        return self._compare(other) < 0

    def __le__(self, other: "IRange") -> bool:
        """Return True if this range sorts before or equals `other`."""
        return self._compare(other) <= 0

    def __gt__(self, other: "IRange") -> bool:
        """Return True if this range sorts after `other`."""
        return self._compare(other) > 0

    def __ge__(self, other: "IRange") -> bool:
        """Return True if this range sorts after or equals `other`."""
        return self._compare(other) >= 0

    @override
    def __eq__(self, other: object) -> bool:
        """Return True if `other` is an IRange with the same sort key."""
        if not isinstance(other, IRange):
            return NotImplemented
        return self._compare(other) == 0

    # Distance relationships (by start and end coordinates only).

    def get_distance(self, other: "IRange") -> int:
        """
        Return the signed distance between this range and `other`.

        The distance is measured between the two closest ends:
        negative (minus the overlap length) when the ranges overlap,
        zero when they are directly adjacent, and positive (the gap
        length) when there is a gap between them. Subclasses may
        override this to add higher-level separation (e.g. ranges on
        different chromosomes in CRange are always apart).

        :param other: Another IRange instance.
        :return: The signed distance between the ranges.
        """
        return max(other.start, self.start) - min(other.end, self.end)

    def is_overlap(self, other: "IRange") -> bool:
        """Return True if this range overlaps with `other`."""
        return self.get_distance(other) < 0

    def is_adjacent(self, other: "IRange") -> bool:
        """Return True if `other` starts exactly where this range ends."""
        return self.get_distance(other) == 0

    def is_away(self, other: "IRange") -> bool:
        """Return True if there is a positive gap between this range and `other`."""
        return self.get_distance(other) > 0

    def is_contact(self, other: "IRange") -> bool:
        """Return True if the ranges overlap or are adjacent (no gap)."""
        return self.get_distance(other) <= 0

    # Index-position relevant operations.

    def get_index(self, position: int) -> int:
        """
        Convert an absolute coordinate to a 0-based index within the range.

        :param position: Absolute coordinate, must lie within
                         [start, end).
        :return: The 0-based index of the position.
        :raises ValueError: If the position is out of the range.
        """
        if self.start <= position < self.end:
            return position - self.start
        raise ValueError("Position (%d) is out of range [%d, %d)." % (position, self.start, self.end))

    def get_position(self, index: int) -> int:
        """
        Convert a 0-based index within the range to an absolute coordinate.

        Negative indexes count from the end of the range like Python
        sequence indexes.

        :param index: 0-based index within the range; negative values
                      count backwards from the end.
        :return: The absolute coordinate of the index.
        :raises ValueError: If the index is out of the range.
        """
        length = self.end - self.start
        if -length <= index < 0:
            return self.end + index
        elif 0 <= index < length:
            return self.start + index
        raise ValueError("Index (%d) is out of range (%d, %d)." % (index, -length, length))
