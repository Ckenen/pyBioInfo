# Class: CRange (chromosome range)

import sys
from collections.abc import Iterable
from typing import ClassVar

if sys.version_info >= (3, 12):
    from typing import override
else:
    from typing import Callable  # pyright: ignore[reportUnreachable]

    def override(func: Callable) -> Callable:
        """No-op fallback for Python < 3.12 compatibility."""
        return func

from .interval_range import IRange


class CRange(IRange):
    """
    CRange is a minimal class for a range located on a chromosome.

    Attributes carried in addition to IRange: the chromosome name
    ``chrom`` (a plain read/write attribute), the optional strand
    ``strand`` (None, ".", "+" or "-", validated on assignment, with
    the ``is_forward`` / ``is_reverse`` conveniences), strand-aware
    index/position conversions (on the reverse strand the internal
    indexing runs from the end towards the start), a
    chromosome-aware ordering: instances are grouped by ``chrom``
    in lexicographic order first, then ordered by the inherited
    ``start``/``end`` coordinates, and the text output methods
    (``format`` / ``to_bed_string`` / ``to_gff_string`` /
    ``to_tbl_string``) whose optional columns are resolved via
    attribute lookup.

    Examples:
        >>> obj = CRange("chr1", 5, 15)
        >>> obj.chrom, obj.start, obj.end
        ('chr1', 5, 15)
        >>> obj = CRange("chr1", 5, 15, strand="+")
        >>> obj.strand, obj.is_forward
        ('+', True)
        >>> CRange("chr1", 10, 20) < CRange("chr2", 5, 15)
        True
    """

    # Chromosome name.
    chrom: str

    # Allowed strand values: no strand, unknown, forward and reverse.
    VALID_STRAND_VALUE_LIST: ClassVar["tuple[str | None, str, str, str]"] = (None, ".", "+", "-")

    def __init__(self, chrom: str, start: int, end: int, strand: "str | None" = None) -> None:
        """
        Initialize a chromosome range.

        :param chrom: Chromosome name, e.g. "chr1".
        :param start: Start coordinate (inclusive), must be >= 0.
        :param end: End coordinate (exclusive), must be > start.
        :param strand: Optional strand: None, ".", "+" or "-".
        :raises ValueError: If the coordinates are invalid (see IRange)
                            or if the strand is not a valid value.
        """
        super().__init__(start, end)
        self.chrom = chrom
        self._strand: "str | None" = None
        self.strand = strand  # validated by the setter

    # Access to protected attributes.

    @property
    def strand(self) -> "str | None":
        """Strand of the range: None, ".", "+" (forward) or "-" (reverse)."""
        return self._strand

    @strand.setter
    def strand(self, value: "str | None") -> None:
        """
        Set the strand after validating its value.

        :param value: One of None, ".", "+", "-".
        :raises ValueError: If the value is not a valid strand.
        """
        if value not in self.VALID_STRAND_VALUE_LIST:  # Check strand
            raise ValueError("Invalid strand \"%s\". Valid value [None, \"+\", \"-\", \".\"]" % value)
        self._strand = value

    @property
    def is_forward(self) -> bool:
        """Return True if the range is on the forward ("+") strand."""
        return self.strand == "+"

    @property
    def is_reverse(self) -> bool:
        """Return True if the range is on the reverse ("-") strand."""
        return self.strand == "-"

    # Index-position relevant operations (strand-aware).

    @override
    def get_index(self, position: int, strandness: bool = True) -> int:
        """
        Convert an absolute coordinate to a 0-based index within the range.

        On the reverse strand (with strandness enabled) the index runs
        from the end towards the start, like a sequence read on the
        minus strand.

        :param position: Absolute coordinate, must lie within
                         [start, end).
        :param strandness: Use reverse-strand indexing if on reverse strand.
        :return: The 0-based index of the position.
        :raises ValueError: If the position is out of the range.
        """
        index = super().get_index(position)
        if self.is_reverse and strandness:
            index = len(self) - 1 - index
        return index

    @override
    def get_position(self, index: int, strandness: bool = True) -> int:
        """
        Convert a 0-based index within the range to an absolute coordinate.

        On the reverse strand (with strandness enabled) the index runs
        from the end towards the start, like a sequence read on the
        minus strand. Negative indexes count from the strand-aware end
        of the range like Python sequence indexes.

        :param index: 0-based index within the range; negative values
                      count backwards from the strand-aware end.
        :param strandness: Use reverse-strand indexing if on reverse strand.
        :return: The absolute coordinate of the index.
        :raises ValueError: If the index is out of the range.
        """
        if self.is_reverse and strandness:
            length = len(self)
            if -length <= index < 0:
                index = length + index
            elif not 0 <= index < length:
                raise ValueError("Index (%d) is out of range (%d, %d)." % (index, -length, length))
            index = length - 1 - index
        return super().get_position(index)

    # Formatted output.

    def to_bed_string(self, ncol: int = 6) -> str:
        """
        Format the range as a BED line with the first ``ncol`` columns.

        Optional columns (name, score, thick, color, blocks) are
        resolved via attribute lookup: attributes not provided by the
        class are treated as None, and a range without a block
        structure is rendered as a single block covering the whole
        span.

        :param ncol: Number of columns to emit, 3 to 12.
        :return: Tab-separated BED text line.
        """

        assert ncol >= 3
        bed_cols: "list[object]" = [self.chrom, self.start, self.end]

        blocks = getattr(self, "blocks", None)
        if blocks is None:
            # No block structure: treat the whole span as one block.
            blocks = ((self.start, self.end),)

        if ncol >= 4:
            name = getattr(self, "name", None)
            if name is None:
                name = "None"
            bed_cols.append(name)
        if ncol >= 5:
            score = getattr(self, "score", None)
            if score is None:
                score = "."
            bed_cols.append(score)
        if ncol >= 6:
            bed_cols.append(self.strand)
        if ncol >= 7:
            thick = getattr(self, "thick", None)
            bed_cols.append(thick[0] if thick is not None else self.start)
        if ncol >= 8:
            thick = getattr(self, "thick", None)
            bed_cols.append(thick[1] if thick is not None else self.start)
        if ncol >= 9:
            color = getattr(self, "color", None)
            if color is None:
                color = "0,0,255" if self.is_forward else "255,0,0"
            bed_cols.append(color)
        if ncol >= 10:
            bed_cols.append(len(blocks))
        if ncol >= 11:
            block_sizes = [v2 - v1 for v1, v2 in blocks]
            bed_cols.append(",".join(map(str, block_sizes)) + ",")
        if ncol >= 12:
            block_offsets = [v1 - self.start for v1, _ in blocks]
            bed_cols.append(",".join(map(str, block_offsets)) + ",")
        return "\t".join(map(str, bed_cols))

    def to_gff_string(self, attris: "list[str] | None" = None) -> str:
        """
        Format the range as a GFF line.

        Optional columns (source, feature, score, frame, name) are
        resolved via attribute lookup: attributes not provided by the
        class are rendered as "." (or None for the name attribute).

        :param attris: Optional list of extra attribute keys appended
                       after the mandatory "name" attribute.
        :return: Tab-separated GFF text line.
        """
        gff_cols: "list[object]" = [self.chrom]
        source = getattr(self, "source", None)
        if source is None:
            source = "."
        gff_cols.append(source)
        feature = getattr(self, "feature", None)
        if feature is None:
            feature = "."
        gff_cols.append(feature)

        gff_cols.append(self.start + 1)
        gff_cols.append(self.end)
        score = getattr(self, "score", None)
        if score is None:
            score = "."
        gff_cols.append(score)

        gff_cols.append(self.strand)
        frame = getattr(self, "frame", None)
        if frame is None:
            frame = "."
        gff_cols.append(frame)
        items: "list[str]"
        if attris is None:
            items = ["name"]
        else:
            items = [str(item) for item in attris]
        if "name" not in items:
            items.insert(0, "name")
        array: "list[str]" = []
        for key in items:
            if key == "name":
                value: object = getattr(self, "name", None)
            else:
                # Note: the index protocol only supports 0/1; str keys
                # are kept for compatibility with the original design.
                value = self[key]  # pyright: ignore[reportArgumentType]
            array.append("%s=%s" % (key, str(value)))
        gff_cols.append("; ".join(array))
        return "\t".join(map(str, gff_cols))

    def to_tbl_string(self, header: "Iterable[str] | None") -> str:
        """
        Format the range as a TBL line with the columns named by ``header``.

        Supported special keys: ``"blocks"`` (rendered as a ";"
        joined list of "start-end" block intervals, falling back to
        the whole span as a single block when the class provides no
        block structure) and ``"thick"`` (rendered as "start-end",
        or an empty cell when absent). All other keys are resolved
        via attribute lookup, falling back to the index protocol.

        :param header: Iterable of column keys.
        :return: Tab-separated TBL text line.
        :raises RuntimeError: If the header is None.
        """
        if header is None:
            raise RuntimeError()
        tbl_cols: "list[object]" = []
        for key in [str(item) for item in header]:
            if key == "blocks":
                blocks = getattr(self, "blocks", None)
                if blocks is None:
                    # No block structure: treat the whole span as one block.
                    blocks = ((self.start, self.end),)
                temps: "list[str]" = []
                for x, y in blocks:
                    temps.append("%d-%d" % (x, y))
                tbl_cols.append(";".join(temps))
            elif key == "thick":
                thick = getattr(self, "thick", None)
                if thick is None:
                    tbl_cols.append("")
                else:
                    thick_start, thick_end = thick
                    tbl_cols.append("%d-%d" % (thick_start, thick_end))
            else:
                try:
                    value = self.__getattribute__(key)
                except AttributeError:
                    # Fall back to the index protocol; it only supports
                    # 0/1, kept for compatibility with the original design.
                    value = self[key]  # pyright: ignore[reportArgumentType]
                tbl_cols.append(value)
        return "\t".join(map(str, tbl_cols))

    def format(self, fmt: str = "bed", ncol: int = 12, attris: "list[str] | None" = None, header: "Iterable[str] | None" = None) -> str:
        """
        Format the range as a line of the given text format.

        :param fmt: Format name, case-insensitive: "bed", "gff" or "tbl".
        :param ncol: Number of columns for the BED format (default 12).
        :param attris: Extra attribute keys for the GFF format.
        :param header: Column keys for the TBL format.
        :return: Tab-separated text line.
        :raises RuntimeError: If the format is unsupported or the header
                              is missing.
        """
        fmt = fmt.upper()
        if fmt == "BED":
            return self.to_bed_string(ncol=ncol)
        elif fmt == "GFF":
            return self.to_gff_string(attris=attris)
        elif fmt == "TBL":
            return self.to_tbl_string(header=header)
        raise RuntimeError("Unsupported format %s" % fmt)

    # Comparison and ordering (chromosome-aware).

    @override
    def _compare(self, other: "IRange") -> int:
        """
        Compare this range with another one, chromosome first.

        Ordering rule: by chromosome name first (compared
        alphabetically), then by the inherited ``start``/``end``
        coordinates within the same chromosome. Plain IRange operands
        (without a chromosome) fall back to the coordinate-only
        comparison.

        :param other: Another IRange or CRange instance.
        :return: Negative if self < other, 0 if equal, positive if self > other.
        """
        if isinstance(other, CRange):
            if self.chrom < other.chrom:
                return -1
            elif self.chrom == other.chrom:
                return super()._compare(other)
            return 1
        return super()._compare(other)

    # Distance relationships (chromosome-aware).

    @override
    def get_distance(self, other: "IRange") -> int:
        """
        Return the signed distance between this range and `other`.

        Ranges located on different chromosomes can never overlap or
        touch, so they are reported as maximally apart (a positive
        sentinel distance): they are always "away". Ranges on the same
        chromosome (and plain IRange operands without a chromosome)
        fall back to the coordinate-only distance of the base class.

        :param other: Another IRange or CRange instance.
        :return: The signed distance between the ranges.
        """
        if isinstance(other, CRange) and self.chrom != other.chrom:
            return sys.maxsize
        return super().get_distance(other)
