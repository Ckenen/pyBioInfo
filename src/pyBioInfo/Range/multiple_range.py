# Class: MRange (multiple genomic range())
from typing import Iterable, Sequence

from .genomic_range import GRange
from pyBioInfo.Utils import BlockTools


class MRange(GRange):
    """
    MRange is a genomic range composed of multiple groups of blocks.

    Each group of blocks (e.g. the blocks of each mate of a paired-end
    fragment) is stored in ``blocks_array``; the overall ``blocks`` of
    the range are the sorted and sutured union of all groups.

    Examples:
        >>> obj = MRange("chr1", blocks_array=[[(10, 20), (50, 60)],
        ...                                    [(30, 40)]])
        >>> obj.blocks                    # sutured union, gap <= 0
        ((10, 20), (30, 40), (50, 60))
    """

    # Immutable array of block groups; each element is a group of
    # (start, end) intervals (e.g. the blocks of one alignment mate).
    _blocks_array: "tuple[tuple[tuple[int, int], ...], ...]"

    def __init__(self, chrom: str, start: "int | None" = None, end: "int | None" = None,
                 name: "str | None" = None, strand: "str | None" = None,
                 blocks: "Sequence[Sequence[int]] | None" = None,
                 blocks_array: "Iterable[Sequence[Sequence[int]]] | None" = None) -> None:
        """
        Initialize a multiple genomic range.

        Either pass ``blocks_array`` (multiple groups of blocks; the
        overall blocks are computed as the sorted and sutured union of
        all groups), or pass ``start``/``end`` (a single trivial group
        is then created). ``blocks`` is accepted for signature
        compatibility with GRange but currently ignored.

        :param chrom: Chromosome name, e.g. "chr1".
        :param start: Overall start (inclusive); only used when
                      ``blocks_array`` is not given.
        :param end: Overall end (exclusive); only used when
                    ``blocks_array`` is not given.
        :param name: Optional feature name.
        :param strand: Optional strand: None, ".", "+" or "-".
        :param blocks: Ignored; kept for base-class compatibility.
        :param blocks_array: Iterable of block groups, each being an
                             iterable of (start, end) intervals.
        :raises ValueError: If neither ``blocks_array`` nor a valid pair
                            of ``start``/``end`` is given.
        :raises ValueError: If any group of blocks is invalid.
        """
        blocks_tuple: "tuple[tuple[int, int], ...]"
        if blocks_array is None:
            if start is None:
                raise ValueError()
            if end is None:
                raise ValueError()
            blocks_tuple = ((start, end),)
            blocks_array = ((start, end),)
        else:
            # Normalize each group into a tuple of int tuples and collect
            # all intervals for the overall (union) blocks.
            groups: "list[tuple[int, int]]" = []
            normalized_groups: "list[tuple[tuple[int, int], ...]]" = []
            for group in blocks_array:
                BlockTools.check_blocks(group)
                group_tuple = tuple((block_start, block_end) for block_start, block_end in group)
                normalized_groups.append(group_tuple)
                groups.extend(group_tuple)
            # Sort and suture the union of all groups; gaps of size <= 0
            # (i.e. touching or overlapping intervals) are merged.
            blocks_tuple = tuple((block_start, block_end) for block_start, block_end
                                 in BlockTools.suture(BlockTools.sorted(groups)))
            blocks_array = tuple(normalized_groups)
        assert blocks_tuple
        super(MRange, self).__init__(chrom=chrom, start=start, end=end,
                                     name=name, strand=strand, blocks=blocks_tuple)
        self._blocks_array = blocks_array

    @property
    def blocks_array(self) -> "tuple[tuple[tuple[int, int], ...], ...]":
        """Groups of blocks of the range, each as a tuple of intervals."""
        return self._blocks_array

    def coincide_from(self, other: "MRange") -> bool:
        """
        Test whether every group of ``self`` coincides with ``other``.

        A group of ``self`` coincides with ``other`` when all of its
        block boundaries are aligned to blocks of ``other`` (see
        ``BlockTools.is_coincide`` for the detailed directional rules).

        :param other: Another multiple genomic range.
        :return: True if all groups of self coincide with other.
        """
        for group in self._blocks_array:
            if not BlockTools.is_coincide(other.blocks, group):
                return False
        return True
