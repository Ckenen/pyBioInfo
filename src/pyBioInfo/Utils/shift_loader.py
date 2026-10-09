# ShiftLoader: a buffered stream of sorted items supporting windowed fetches.
from collections.abc import Generator, Iterable
from typing import Generic, Protocol, TypeVar


class _ShiftItem(Protocol):
    """Structural type of the items handled by ShiftLoader."""

    @property
    def chrom(self) -> str: ...

    @property
    def start(self) -> int: ...

    @property
    def end(self) -> int: ...


_T = TypeVar("_T", bound="_ShiftItem")


class ShiftLoader(Generic[_T]):
    """
    A buffered loader over a stream sorted by chrom, start and end.

    Items are consumed lazily from the input iterable into an internal
    buffer. ``fetch`` yields the items overlapping a query window
    (or the window of a query object) and drops from the buffer the
    items that lie entirely before it. Because yielded items stay in
    the buffer, consecutive overlapping queries do not re-read the
    input.
    """

    # Iterable of sorted items to load from.
    _items: "Iterable[_T]"
    # Iterator pulling from _items (created once in __init__).
    _iterator: "Generator[_T, None, None]"
    # Buffered items not yet passed by the current query window.
    _buffer: "list[_T]"
    # Whether the input iterator is exhausted.
    _reach_tail: bool
    # Whether every loaded item is recorded in histories.
    keep_history: bool
    # All loaded items (only populated when history is True).
    histories: "list[_T]"
    # Number of items loaded so far.
    size: int
    # Chrom of the previous query window.
    _last_chrom: "str | None"
    # Start of the previous query window.
    _last_start: "int | None"
    # End of the previous query window.
    _last_end: "int | None"

    def __init__(self, items: "Iterable[_T]", keep_history: bool = False) -> None:
        """
        Initialize a ShiftLoader.

        :param items: Iterable of items sorted by chrom, start, end.
        :param history: If True, every loaded item is appended to
                        ``self.histories``.
        """
        self._items = items
        self._iterator = self.__iter__()
        self._buffer = []
        self._reach_tail = False
        self.keep_history = keep_history
        self.histories = []
        self.size = 0
        self._last_chrom = None
        self._last_start = None
        self._last_end = None

    def __iter__(self) -> "Generator[_T, None, None]":
        """Yield the items of the input iterable, one by one."""
        for obj in self._items:
            yield obj

    def _load_next_item(self) -> "_T | None":
        """
        Pull the next item from the input iterator into the buffer.

        :return: The next item, or None if the input is exhausted.
        """
        obj: "_T | None" = None
        if not self._reach_tail:
            try:
                obj = next(self._iterator)
                self.size += 1
                if self.keep_history:
                    self.histories.append(obj)
            except StopIteration:
                self._reach_tail = True
        return obj

    def _check_chrom_start_end(self, chrom: str, start: int, end: int) -> None:
        """
        Assert that query windows are requested in sorted order.

        :param chrom: Chrom of the query window.
        :param start: Start of the query window.
        :param end: End of the query window.
        :raises RuntimeError: If the window moves backwards.
        """
        if self._last_chrom is not None:
            if chrom < self._last_chrom:
                raise RuntimeError("Error in sorted chrom [last: %s, current: %s]" % (self._last_chrom, chrom))
            if chrom == self._last_chrom and self._last_start is not None and start < self._last_start:
                raise RuntimeError(
                    "Error in sorted start [chrom: %s, last: %d-%d, current: %d-%d]" % (
                        chrom, self._last_start, self._last_end or -1, start, end))
        self._last_chrom = chrom
        self._last_start = start
        self._last_end = end

    def flush(self) -> None:
        """Consume the remaining input (recorded in histories only when enabled)."""
        while True:
            obj = self._load_next_item()
            if obj is None:
                break

    def fetch(self, chrom: "str | None" = None, start: "int | None" = None,
              end: "int | None" = None, obj: "_ShiftItem | None" = None) -> "Generator[_T, None, None]":
        """
        Yield the buffered items overlapping a query window.

        Either provide the window explicitly (``chrom``, ``start``,
        ``end``) or via a query ``obj`` carrying the same attributes.
        Items ending before the window start are dropped from the
        buffer; the iteration stops at the first item starting at or
        after the window end.

        :param chrom: Chrom of the query window.
        :param start: Start of the query window (0-based, inclusive).
        :param end: End of the query window (exclusive).
        :param obj: A query object providing chrom/start/end.
        :raises RuntimeError: If the query window moves backwards.
        """

        if obj is not None:
            chrom = obj.chrom
            start = obj.start
            end = obj.end

        if chrom is None or start is None or end is None:
            raise ValueError("Either provide the window explicitly (chrom, start, end) or via a query obj")

        self._check_chrom_start_end(chrom, start, end)
        
        i = 0
        while True:
            if i >= len(self._buffer):
                if self._reach_tail:
                    break
                else:
                    new_obj = self._load_next_item()
                    if new_obj is None:
                        break
                    else:
                        self._buffer.append(new_obj)
            else:
                item = self._buffer[i]
                if item.chrom < chrom:
                    _ = self._buffer.pop(i)
                elif item.chrom == chrom:
                    if item.end <= start:
                        _ = self._buffer.pop(i)
                    elif item.start >= end:
                        break
                    else:
                        yield item
                        i += 1
                else:
                    break
