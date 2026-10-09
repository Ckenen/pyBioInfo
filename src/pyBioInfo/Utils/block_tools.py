import functools
from collections.abc import Iterable, Sequence


class BlockTools(object):
    """
    BlockTools is a collection of static utilities for block lists.

    A "block" is a half-open interval [start, end) represented as a
    two-element sequence; a block list is a sequence of such blocks,
    usually sorted by coordinates and non-overlapping (see
    :py:meth:`check_blocks`).
    """

    @classmethod
    def check_blocks(cls, blocks: "Sequence[Sequence[int]]") -> None:
        """Ensure the value of the coordinate of blocks are valid.
        
        Rules:
        1. start >= 0
        2. start >= previous block end
        3. start < end

        Args:
            blocks (list): list of blocks.

        Raises:
            RuntimeError: Empty block list.
            ValueError: The value of start or end not in compliance with the rules.
        """
        if len(blocks) == 0:
            raise RuntimeError("Empty blocks.")
        y0 = 0
        for x, y in blocks:
            if x < 0: # x >= 0.
                raise ValueError("Not in compliance with the rules: start (%d) >= 0." % x)
            if x < y0: # Allow zero-length gap, otherwise x <= y0.
                raise ValueError("Not in compliance with the rules: start (%d) >= previous block end (%d)." % (x, y0))
            if x >= y: # Not allow zero-length block, otherwise y < x.
                raise ValueError("Not in compliance with the rules: start (%d) < end (%d)." % (x, y))
            y0 = y
            
    @classmethod
    def is_valid_blocks(cls, blocks: "Sequence[Sequence[int]]") -> bool:
        """Return whether the blocks are valid.

        Args:
            blocks (list): list of blocks.

        Returns:
            bool: `True` if the blocks are valid.
        """
        if len(blocks) == 0:
            return False
        y0 = 0
        for x, y in blocks:
            if x < 0 or x < y0 or x >= y:
                return False
            y0 = y
        return True
    
    @classmethod
    def check_sorted(cls, blocks: "Sequence[Sequence[int]]") -> None:
        """Check the blocks are sorted by coordinates.

        Args:
            blocks (list): list of blocks.

        Raises:
            ValueError: The blocks are unsorted.
        """
        x0: "int | None" = None
        y0: "int | None" = None
        for x, y in blocks:
            if x0 is not None and y0 is not None:
                if (x < x0) or (x == x0 and y < y0):
                    raise ValueError("The blocks are unsorted!")
            x0, y0 = x, y
            
    @classmethod
    def is_sorted(cls, blocks: "Sequence[Sequence[int]]") -> bool:
        """Return whether the blocks are sorted by coordinates.

        Args:
            blocks (list): list of blocks.

        Returns:
            bool: `True` if the blocks are sorted.
        """
        x0: "int | None" = None
        y0: "int | None" = None
        for x, y in blocks:
            if x0 is not None and y0 is not None:
                if (x < x0) or (x == x0 and y < y0):
                    return False
            x0, y0 = x, y
        return True
    
    @classmethod
    def _full_block_comparison(cls, block1: "Sequence[int]", block2: "Sequence[int]") -> int:
        """Comparison between blocks by both start and end coordinates.

        Args:
            block1 (list): list of start and end.
            block2 (list): list of start and end

        Returns:
            int: return -1 (<), 0 (=) or 1 (>).
        """
        x1, y1 = block1
        x2, y2 = block2
        if x1 < x2:
            return -1
        elif x1 == x2:
            if y1 < y2:
                return -1
            elif y1 == y2:
                return 0
        return 1

    @classmethod
    def sorted(cls, blocks: "Iterable[Sequence[int]]") -> "list[Sequence[int]]":
        """Sorting blocks.

        Args:
            blocks (list): list of blocks.

        Returns:
            list: new list of sorted blocks.
        """
        return sorted(blocks, key=functools.cmp_to_key(cls._full_block_comparison))

    @classmethod
    def get_length(cls, blocks: "Iterable[Sequence[int]]") -> int:
        """Get the summation length of all blocks.

        Args:
            blocks (list): list of blocks.

        Returns:
            int: Summation length of all blocks.
        """
        return sum([y - x for x, y in blocks])

    @classmethod
    def get_index(cls, blocks: "Iterable[Sequence[int]]", position: int,
                  reverse: bool = False, length: "int | None" = None) -> int:
        """Get the index of certain position in blocks.

        Indexes run over the concatenated blocks, skipping the gaps
        between them.

        Args:
            blocks (list): list of blocks.
            position (int): value of position.
            reverse (bool, optional): return the reverse index. Defaults to False.
            length (int, optional): the summation length of all blocks. Defaults to None.

        Raises:
            ValueError: The position is not covered by any block.

        Returns:
            int: Index of certain position.
        """
        index: "int | None" = None
        offset = 0
        for x, y in blocks:
            if position >= y:
                offset += y - x
            elif x <= position < y:
                index = position - x + offset
                break
            elif position < x:
                break
        if index is None:
            raise ValueError("Position (%d) not in blocks." % position)
        if reverse:
            if length is None:
                length = cls.get_length(blocks)
            index = length - 1 - index
        return index
    
    @classmethod
    def get_position(cls, blocks: "Iterable[Sequence[int]]", index: int,
                     reverse: bool = False, length: "int | None" = None) -> int:
        """Get the position of certain index in blocks.

        Indexes run over the concatenated blocks, skipping the gaps
        between them; negative indexes count from the end.

        Args:
            blocks (list): list of blocks.
            index (int): 0-based index within the concatenated blocks.
            reverse (bool, optional): treat the index as on the reverse
                strand. Defaults to False.
            length (int, optional): the summation length of all blocks.
                Defaults to None.

        Raises:
            ValueError: The index is out of the blocks range.
            RuntimeError: The index is not covered by any block (invalid blocks).

        Returns:
            int: Position of certain index.
        """
        position: "int | None" = None
        if length is None:
            length = cls.get_length(blocks)
        if index >= length or index < -length:
            raise ValueError("Index (%d) out of blocks range (%d <= index < %d)." % (index, -length, length))
        if index < 0:
            index += length
        if reverse:
            index = length - 1 - index
        i1, i2 = 0, 0
        for x, y in blocks:
            i1, i2 = i2, i2 + y - x
            if index >= i2:
                continue
            elif i1 <= index < i2:
                position = index - i1 + x
                break
            else:
                raise RuntimeError("Unknown error.")
        if position is None:
            # Only reachable with invalid blocks (total length mismatch).
            raise RuntimeError("Index (%d) not found in blocks." % index)
        return position
    
    @classmethod
    def is_include(cls, blocks: "Iterable[Sequence[int]]", position: int) -> bool:
        """Determine if `position` is included in blocks.

        Args:
            blocks (list): List of blocks.
            position (int): Value of position.

        Returns:
            bool: `True` if `position` is included in blocks, else `False`.
        """
        for x, y in blocks:
            if position >= y:
                continue
            elif x <= position < y:
                return True
            elif position < x:
                return False
        return False
    
    @classmethod
    def move(cls, blocks: "Sequence[Sequence[int]]", step: int = 0) -> "list[tuple[int, int]]":
        """Move blocks.

        Args:
            blocks (list): list of blocks (start, end).
            step (int, optional): steps to move. Defaults to 0.

        Raises:
            RuntimeError: Invalid blocks after moving.

        Returns:
            list: new list of blocks.
        """
        blocks2: "list[tuple[int, int]]" = []
        if blocks[0][0] + step < 0:
            raise RuntimeError("Invalid blocks after moving.")
        for x, y in blocks:
            blocks2.append((x + step, y + step))
        return blocks2

    @classmethod
    def gaps(cls, blocks: "Iterable[Sequence[int]]") -> "list[tuple[int, int]]":
        """Get the gaps between adjacent blocks.

        Args:
            blocks (list): list of blocks.

        Returns:
            list: list of gaps.
        """
        blocks2: "list[tuple[int, int]]" = []
        previous_end: "int | None" = None
        for x, y in blocks:
            if previous_end is not None:
                blocks2.append((previous_end, x))
            previous_end = y
        return blocks2

    @classmethod
    def suture(cls, blocks: "Iterable[Sequence[int]]", gap: int = 0) -> "list[tuple[int, int]]":
        """Filling the gaps between blocks.

        Adjacent blocks whose gap size is not larger than `gap` are
        merged into a single block.

        Args:
            blocks (list): list of blocks.
            gap (int, optional): maximum gap size. Defaults to 0.

        Returns:
            list: list of blocks after suturing.
        """
        blocks2: "list[tuple[int, int]]" = []
        iterator = iter(blocks)
        first = next(iterator, None)
        if first is None:
            return blocks2
        x0, y0 = first
        for x, y in iterator:
            if x - y0 <= gap:
                x0, y0 = min(x0, x), max(y0, y)
            else:
                blocks2.append((x0, y0))
                x0, y0 = x, y
        blocks2.append((x0, y0))
        return blocks2

    @classmethod
    def fusion(cls, blocks1: "Iterable[Sequence[int]]", blocks2: "Iterable[Sequence[int]]") -> "list[tuple[int, int]]":
        """Merge two lists of blocks into one.

        Args:
            blocks1 (list): list of blocks.
            blocks2 (list): another list of blocks.

        Returns:
            list: list of fusion blocks.
        """
        blocks: "list[Sequence[int]]" = list(blocks1) + list(blocks2)
        blocks = sorted(blocks, key=lambda item: item[0])
        return cls.suture(blocks, 0)
        
    @classmethod
    def is_contain(cls, blocks1: "Sequence[Sequence[int]]", blocks2: "Sequence[Sequence[int]]") -> bool:
        """Return whether `blocks1` fully covers `blocks2`.

        Every block (or fragment) of `blocks2` must be covered by the
        blocks of `blocks1`.

        Args:
            blocks1 (list): list of covering blocks.
            blocks2 (list): list of covered blocks.

        Returns:
            bool: `True` if `blocks1` contains `blocks2`.
        """
        is_contain = True
        length1 = len(blocks1)
        length2 = len(blocks2)
        index1 = 0
        index2 = 0
        while index2 < length2:
            block2 = blocks2[index2]
            start2, end2 = block2
            while True:
                if index1 >= length1:
                    is_contain = False
                    break
                block1 = blocks1[index1]
                start1, end1 = block1
                if start2 < start1:
                    is_contain = False
                    break
                else:
                    if end1 >= end2:
                        break
                    else:
                        if end1 <= start2:
                            index1 += 1
                            continue
                        else:
                            start2 = end1
                            index1 += 1
                            continue
            if not is_contain:
                break
            index2 += 1
        return is_contain

    @classmethod
    def is_coincide(cls, 
        blocks_ref: "Sequence[Sequence[int]]", 
        blocks_que: "Sequence[Sequence[int]]", 
        tolerance: int = 0) -> bool:
        """Return whether `blocks_que` coincides with `blocks_ref`.

        `blocks_que` must have no more blocks than `blocks_ref`. The
        blocks of `blocks_que` must match a consecutive run of the
        blocks of `blocks_ref` (no reference block may be skipped). The
        first block of `blocks_que` may start later than the reference
        block it aligns to, and its last block may end earlier;
        intermediate block boundaries may differ by up to `tolerance`.
        The check is directional.

        Args:
            blocks_ref (list): list of reference blocks.
            blocks_que (list): list of query blocks.
            tolerance (int, optional): allowed boundary deviation for
                intermediate blocks. Defaults to 0.

        Returns:
            bool: `True` if `blocks_que` coincides with `blocks_ref`.
        """

        count_ref = len(blocks_ref)
        count_que = len(blocks_que)
        if count_ref < count_que:
            return False

        for offset in range(count_ref - count_que + 1):
            matched = True
            for i_que, (start_que, end_que) in enumerate(blocks_que):
                i_ref = offset + i_que
                start_ref, end_ref = blocks_ref[i_ref]

                if i_que == 0:
                    if i_ref == 0:
                        matched = start_que < end_ref
                    else:
                        matched = start_que < end_ref \
                            and start_que >= start_ref - tolerance \
                            and start_que > blocks_ref[i_ref - 1][1]
                else:
                    if i_ref == 0:
                        raise RuntimeError("Invalid blocks after moving.")
                    else:
                        matched = start_que < end_ref \
                            and abs(start_que - start_ref) <= tolerance \
                            and start_que > blocks_ref[i_ref - 1][1]
                if not matched:
                    break

                if i_que == count_que - 1:
                    if i_ref == count_ref - 1:
                        matched = end_que > start_ref
                    else:
                        matched = end_que > start_ref \
                            and end_que <= end_ref + tolerance \
                            and end_que < blocks_ref[i_ref + 1][0]
                else:
                    if i_ref == count_ref - 1:
                        raise RuntimeError("Invalid blocks after moving.")
                    else:
                        matched = end_que > start_ref \
                            and abs(end_que - end_ref) <= tolerance \
                            and end_que < blocks_ref[i_ref + 1][0]
                if not matched:
                    break
                
            if matched:
                return True

        return False

    @classmethod
    def get_boundary_bias(cls, blocks_ref: "Sequence[Sequence[int]]", blocks_que: "Sequence[Sequence[int]]") -> "tuple[int, int]":
        """Get the boundary bias of `blocks_que` relative to `blocks_ref`.

        The boundary of a block list is defined by the smallest block
        start and the largest block end (both lists are assumed sorted).
        The difference is measured per boundary:

        - ``start``: how far ``blocks_que``'s start lies from
          ``blocks_ref``'s start. Negative when ``blocks_que`` starts
          earlier (outside on the left), positive when it starts later
          (inside).
        - ``end``: how far ``blocks_que``'s end lies from
          ``blocks_ref``'s end. Positive when ``blocks_que`` ends later
          (outside on the right), negative when it ends earlier
          (inside).

        When a boundary lies outside ``blocks_ref`` the plain coordinate
        difference is used; when it lies inside, the block-aware index
        (gaps skipped) is used instead.

        Args:
            blocks_ref (list): list of reference blocks.
            blocks_que (list): list of query blocks.

        Raises:
            ValueError: If either block list is invalid, or if a query
                boundary falls inside a gap of ``blocks_ref``.

        Returns:
            tuple: (start_difference, end_difference).
        """
        cls.check_blocks(blocks_ref)
        cls.check_blocks(blocks_que)

        start1=  blocks_ref[0][0]
        end1 = blocks_ref[-1][1]
        start2 = blocks_que[0][0]
        end2 = blocks_que[-1][1]

        if start2 <= start1:
            difference_start = start2 - start1
        else:
            try:
                difference_start = cls.get_index(blocks_ref, start2)
            except ValueError:
                raise ValueError(f"Query start {start2} falls inside a gap of reference blocks.")
                
        if end2 >= end1:
            difference_end = end2 - end1
        else:
            try:
                difference_end = cls.get_index(blocks_ref, end2 - 1) - cls.get_length(blocks_ref)
            except ValueError:
                raise ValueError(f"Query end {end2} falls inside a gap of reference blocks.")

        return difference_start, difference_end

    @classmethod
    def is_matched_gap(cls, blocks_ref: "Sequence[Sequence[int]]", blocks_que: "Sequence[Sequence[int]]",
                       tolerance: int = 0) -> bool:
        """Return whether the gaps of `blocks_que` match a subset of the gaps of `blocks_ref`.

        `blocks_que` must have no more blocks than `blocks_ref`. The
        first gap of `blocks_que` may match any gap of `blocks_ref`;
        once matched, the remaining gaps of `blocks_que` must match the
        following gaps of `blocks_ref` consecutively. A gap matches when
        both its start and its end differ by no more than `tolerance`.

        Args:
            blocks_ref (list): list of reference blocks.
            blocks_que (list): list of query blocks.
            tolerance (int, optional): allowed coordinate deviation for
                each gap boundary. Defaults to 0.

        Raises:
            ValueError: If either block list is invalid.

        Returns:
            bool: `True` if the gaps of `blocks_que` match a consecutive
                run of the gaps of `blocks_ref`.
        """
        cls.check_blocks(blocks_ref)
        cls.check_blocks(blocks_que)
        if len(blocks_que) > len(blocks_ref):
            return False
        gaps_ref = cls.gaps(blocks_ref)
        gaps_que = cls.gaps(blocks_que)
        if not gaps_que:
            return True
        # The first query gap may match any reference gap; once matched,
        # the remaining query gaps must match consecutively.
        for start in range(len(gaps_ref) - len(gaps_que) + 1):
            matched = True
            for i, (qs, qe) in enumerate(gaps_que):
                rs, re = gaps_ref[start + i]
                if abs(rs - qs) > tolerance or abs(re - qe) > tolerance:
                    matched = False
                    break
            if matched:
                return True
        return False

    @classmethod
    def clip(cls, blocks: "Sequence[Sequence[int]]", start: int, end: int,
             template: "Sequence[Sequence[int]] | None" = None,
             extend: bool = False, extend_left: bool = False, extend_right: bool = False,
             supply: bool = False, supply_left: bool = False, supply_right: bool = False,
             limited: "tuple[int, int] | None" = None, suture: bool = False) -> "list[tuple[int, int]]":
        """Clip the blocks to the window [start, end).

        Only the blocks (or block fragments) inside the window are kept.
        Optionally, the clipped result may be extended towards the
        template blocks and supplied with plain filler blocks.

        Args:
            blocks (list): list of blocks to clip.
            start (int): window start (inclusive), must be >= 0.
            end (int): window end (exclusive), must be > start.
            template (list, optional): template blocks limiting the
                extend regions. Defaults to None.
            extend (bool, optional): enable both-side extension. Defaults to False.
            extend_left (bool, optional): enable left-side extension. Defaults to False.
            extend_right (bool, optional): enable right-side extension. Defaults to False.
            supply (bool, optional): enable both-side supply. Defaults to False.
            supply_left (bool, optional): enable left-side supply. Defaults to False.
            supply_right (bool, optional): enable right-side supply. Defaults to False.
            limited (tuple, optional): valid range (start, end) for the
                supply regions. Defaults to None.
            suture (bool, optional): suture the result blocks. Defaults to False.

        Raises:
            ValueError: Invalid window or template coordinates.

        Returns:
            list: list of clipped blocks.
        """

        # check the parameter.
        if start < 0:
            raise ValueError()

        if start >= end:
            raise ValueError()

        target_start = start
        target_end = end

        origin_start = blocks[0][0]
        origin_end = blocks[-1][1]

        # The limited range for extend.
        template_start = 0
        template_end: "int | None" = None
        if template is not None:
            template_start = template[0][0]
            template_end = template[-1][1]
            if template_start > origin_start:
                raise ValueError()
            if template_end < origin_end:
                raise ValueError()

        # The limited range for supply.
        limited_start = 0
        limited_end = 0
        if limited is not None:
            limited_start, limited_end = limited
        if limited_start:
            assert limited_start <= origin_start
            assert template_start <= limited_end
            assert limited_start <= target_start
        if limited_end:
            assert limited_end >= origin_end
            if template_end is not None:
                assert limited_end >= template_end
            assert limited_end >= target_end

        array_core: "list[tuple[int, int]]" = []
        array_extend_left: "list[tuple[int, int]]" = []
        array_extend_right: "list[tuple[int, int]]" = []
        array_supply_left: "list[tuple[int, int]]" = []
        array_supple_right: "list[tuple[int, int]]" = []

        for block_start, block_end in blocks:
            if block_end <= target_start:
                continue
            elif block_start >= target_end:
                break
            block_start1 = max(target_start, block_start)
            block_end1 = min(target_end, block_end)
            if block_start1 < block_end1:
                array_core.append((block_start1, block_end1))

        if extend:
            extend_left = True
            extend_right = True

        if supply:
            supply_left = True
            supply_right = True

        if target_start < origin_start and extend_left:
            if template is None:
                array_extend_left = [(target_start, origin_start)]
            else:
                array_extend_left = cls.clip(blocks=template,
                                             start=max(target_start, template_start),
                                             end=origin_start)
            if (target_start < template_start) and supply_left:
                array_supply_left.append((target_start, template_start))

        if target_end > origin_end and extend_right:
            if template is None:
                array_extend_right = [(origin_end, target_end)]
            else:
                # Bound to a local so the assert below does not narrow
                # the Optional type of template_end itself.
                template_end1 = template_end
                assert template_end1 is not None
                array_extend_right = cls.clip(blocks=template,
                                              start=origin_end,
                                              end=min(template_end1, target_end))
            if (template_end is not None) and (target_end > template_end) and supply_right:
                array_supple_right.append((template_end, target_end))

        array = array_supply_left + array_extend_left + array_core + array_extend_right + array_supple_right
        if suture:
            array = cls.suture(array)
        return array

    @classmethod
    def extend(cls, blocks: "Sequence[Sequence[int]]", left: int = 0, right: int = 0,  # pyright: ignore[reportUnusedParameter]
               template: "Sequence[Sequence[int]] | None" = None, extend: bool = False, supply: bool = False,  # pyright: ignore[reportUnusedParameter]
               trim: bool = False, valid: "tuple[int, int] | None" = None, suture: bool = False) -> "list[tuple[int, int]]":  # pyright: ignore[reportUnusedParameter]
        """Extend the blocks towards left and/or right. (Not implemented yet.)

        Args:
            blocks (list): list of blocks.
            left (int, optional): steps to extend towards left. Defaults to 0.
            right (int, optional): steps to extend towards right. Defaults to 0.
            template (list, optional): template blocks for extension. Defaults to None.
            extend (bool, optional): extend with template blocks if possible. Defaults to False.
            supply (bool, optional): supply with plain blocks if template is exhausted. Defaults to False.
            trim (bool, optional): trim to the valid range instead of raising. Defaults to False.
            valid (tuple, optional): valid range (start, end) for the extension. Defaults to None.
            suture (bool, optional): suture the result blocks. Defaults to False.

        Raises:
            NotImplementedError: Always.

        Returns:
            list: new list of extended blocks.
        """
        raise NotImplementedError()

        # assert isinstance(blocks, list)
        # length = cls.length(blocks)
        # assert length + left + right > 0
        # from pyBioInfo.Range import IRange
        #
        # blocks = [block.copy() for block in blocks]
        # if template is not None:
        #     assert isinstance(template, list)
        # valid_start = 0
        # valid_end = None
        # if valid is not None:
        #     # assert isinstance(valid, IntervalRange)
        #     valid_start = valid.start
        #     valid_end = valid.end
        # if valid_start is not None:
        #     assert valid_start <= blocks[0].start
        # if valid_end is not None:
        #     assert valid_end >= blocks[-1].end
        #
        # i = -1
        # if template is not None:
        #     i = len(template) - 1
        # while left != 0 and len(blocks) > 0:
        #     first = blocks[0]
        #     if left < 0:
        #         length = len(first)
        #         if left + length <= 0:
        #             blocks.pop(0)
        #             left += length
        #         else:
        #             first.start = first.start - left
        #             # left = 0
        #             break
        #     else:
        #         if extend:
        #             if i >= 0:
        #                 temp = template[i]
        #                 start1 = temp.start
        #                 end1 = min(first.start, temp.end)
        #                 length = end1 - start1
        #                 if length <= 0:
        #                     i -= 1
        #                     continue
        #                 start1 = max(start1, end1 - left)
        #                 flag = False
        #                 if start1 < valid_start:
        #                     flag = True
        #                     if trim:
        #                         start1 = valid_start
        #                     else:
        #                         raise ValueError()
        #                 blocks.insert(0, IRange(start1, end1))
        #                 if flag:
        #                     break
        #             elif supply:
        #                 start1 = first.start - left
        #                 end1 = first.end
        #                 flag = False
        #                 if start1 < valid_start:
        #                     flag = True
        #                     if trim:
        #                         start1 = valid_start
        #                     else:
        #                         raise ValueError()
        #                 blocks.insert(0, IRange(start1, end1))
        #                 left = 0
        #                 if flag:
        #                     break
        #
        # i = 0
        # if template is not None:
        #     i = len(template)
        # while right != 0 and len(blocks) > 0:
        #     last = blocks[-1]
        #     if right < 0:
        #         length = len(last)
        #         if length + right <= 0:
        #             blocks.pop(-1)
        #             right += length
        #         else:
        #             last.end = last.end + right
        #             # right = 0
        #             break
        #     else:
        #         if extend:
        #             if i > 0:
        #                 temp = template[-i]
        #                 start1 = max(temp.start, last.end)
        #                 end1 = temp.end
        #                 length = end1 - start1
        #                 if length <= 0:
        #                     i -= 1
        #                     continue
        #                 end1 = min(end1, start1 + right)
        #                 flag = False
        #                 if valid_end is not None and end1 > valid_end:
        #                     flag = True
        #                     if trim:
        #                         end1 = valid_end
        #                     else:
        #                         raise ValueError()
        #                 blocks.append(IRange(start1, end1))
        #                 if flag:
        #                     break
        #             elif supply:
        #                 start1 = last.end
        #                 end1 = start1 + right
        #                 flag = False
        #                 if valid_end is not None and end1 > valid_end:
        #                     flag = True
        #                     if trim:
        #                         end1 = valid_end
        #                     else:
        #                         raise ValueError()
        #                 blocks.append(IRange(start1, end1))
        #                 right = 0
        #                 if flag:
        #                     break
        # if suture:
        #     blocks = cls.suture(blocks, gap=0)
        # return blocks
