import functools

class BlockTools(object):
    
    @classmethod
    def check_blocks(cls, blocks):
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
    def is_valid_blocks(cls, blocks):
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
    def check_sorted(cls, blocks):
        """Check the blocks are sorted by coordinates.

        Args:
            blocks (list): list of blocks.

        Raises:
            ValueError: The blocks are unsorted.
        """
        x0, y0 = None, None
        for x, y in blocks:
            if x0 is not None:
                if (x < x0) or (x == x0 and y < y0):
                    raise ValueError("The blocks are unsorted!")
            x0, y0 = x, y
            
    @classmethod
    def is_sorted(cls, blocks):
        """Return whether the blocks are sorted by coordinates.

        Args:
            blocks (list): list of blocks.

        Returns:
            bool: `True` if the blocks are sorted.
        """
        x0, y0 = None, None
        for x, y in blocks:
            if x0 is not None:
                if (x < x0) or (x == x0 and y < y0):
                    return False
            x0, y0 = x, y
        return True
    
    @classmethod
    def _full_block_comparison(cls, block1, block2):
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
    def sorted(cls, blocks):
        """Sorting blocks.

        Args:
            blocks (list): list of blocks.

        Returns:
            list: new list of sorted blocks.
        """
        return list(sorted(blocks, key=functools.cmp_to_key(cls._full_block_comparison)))

    @classmethod
    def get_length(cls, blocks):
        """Get the summation length of all blocks.

        Args:
            blocks (list): list of blocks.

        Returns:
            int: Summation length of all blocks.
        """
        return sum([y - x for x, y in blocks])

    @classmethod
    def get_index(cls, blocks, position, reverse=False, length=None):
        """Get the index of certain position in blocks.

        Args:
            blocks (list): list of blocks.
            position (int): value of position.
            reverse (bool, optional): return the reverse index. Defaults to False.
            length (int, optional): the summation length of all blocks. Defaults to None.

        Raises:
            ValueError: _description_

        Returns:
            int: Index of certain position.
        """
        index = None
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
    def get_position(cls, blocks, index, reverse=False, length=None):
        """Get the position of certain index in blocks.

        Args:
            blocks (list): _description_
            index (int): _description_
            reverse (bool, optional): _description_. Defaults to False.
            length (int, optional): _description_. Defaults to None.

        Raises:
            ValueError: _description_
            RuntimeError: _description_

        Returns:
            int: Position of certain index.
        """
        position = None
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
        return position
    
    @classmethod
    def is_include(cls, blocks, position):
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
    def move(cls, blocks, step=0):
        """Move blocks.

        Args:
            blocks (list): list of blocks (start, end).
            step (int, optional): steps to move. Defaults to 0.

        Raises:
            RuntimeError: Invalid blocks after moving.

        Returns:
            list: new list of blocks.
        """
        blocks2 = []
        if blocks[0][0] + step < 0:
            raise RuntimeError("Invalid blocks after moving.")
        for x, y in blocks:
            blocks2.append((x + step, y + step))
        return blocks2

    @classmethod
    def gaps(cls, blocks):
        """_summary_

        Args:
            blocks (list): list of blocks.

        Returns:
            list: list of gaps.
        """
        blocks2 = []
        x0, y0 = None, None
        for x, y in blocks:
            if x0 is not None:
                blocks2.append((y0, x))
            x0 = x
            y0 = y
        return blocks2

    @classmethod
    def suture(cls, blocks, gap=0):
        """Filling the gaps between blocks.

        Args:
            blocks (list): list of blocks.
            gap (int, optional): maximum gap size. Defaults to 0.

        Returns:
            list: list of blocks after suturing.
        """
        blocks2 = []
        x0, y0 = None, None
        for x, y in blocks:
            if x0 is None:
                x0, y0 = x, y
            else:
                if x - y0 <= gap:
                    x0, y0 = min(x0, x), max(y0, y)
                else:
                    blocks2.append((x0, y0))
                    x0, y0 = x, y
        if x0 is not None:
            blocks2.append((x0, y0))
        return blocks2

    @classmethod
    def fusion(cls, blocks1, blocks2):
        """Merge two lists of blocks into one.

        Args:
            blocks1 (list): list of blocks.
            blocks2 (list): another list of blocks.

        Returns:
            list: list of fusion blocks.
        """
        blocks = []
        for block in blocks1:
            blocks.append(block)
        for block in blocks2:
            blocks.append(block)
        blocks = list(sorted(blocks, key=lambda item: item[0]))
        return cls.suture(blocks, 0)
        
    @classmethod
    def is_contain(cls, blocks1, blocks2):
        """_summary_

        Args:
            blocks1 (list): _description_
            blocks2 (list): _description_

        Returns:
            bool: _description_
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
    def is_coincide(cls, blocks1, blocks2):
        """_summary_

        Args:
            blocks1 (list): _description_
            blocks2 (list): _description_

        Returns:
            bool: _description_
        """
        is_coincide = True
        length1 = len(blocks1)
        length2 = len(blocks2)
        index1 = 0
        index2 = 0
        while index2 < length2:
            start2, end2 = blocks2[index2]
            while True:
                if index1 >= length1:
                    is_coincide = False
                    break
                start1, end1 = blocks1[index1]
                if start1 < start2:
                    if start2 >= end1:
                        index1 += 1
                        continue
                    else:
                        if end1 < end2:
                            is_coincide = False
                            break
                        elif end1 == end2:
                            if index2 == 0:
                                index1 += 1
                                break
                            else:
                                is_coincide = False
                                break
                        else:
                            if index2 == 0 and length2 == 1:
                                break
                            else:
                                is_coincide = False
                                break

                elif start1 == start2:
                    if end1 < end2:
                        is_coincide = False
                        break
                    elif end1 == end2:
                        index1 += 1
                        break
                    else:
                        if index2 == length2 - 1:
                            break
                        else:
                            is_coincide = False
                            break
                else:
                    is_coincide = False
                    break
            if not is_coincide:
                break
            index2 += 1

        return is_coincide

    @classmethod
    def clip(cls, blocks, start, end, template=None,
             extend=False, extend_left=False, extend_right=False,
             supply=False, supply_left=False, supply_right=False,
             limited=None, suture=False):
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
        template_end = None
        if template is not None:
            template_start = template[0][0]
            template_end = template[-1][1]
            if template_start > origin_start:
                raise ValueError()
            if template_end < origin_end:
                raise ValueError()

        # The limited range for supply.
        limited_start = 0
        limited_end = None
        if limited:
            limited_start, limited_end = limited
        if limited_start:
            assert limited_start <= origin_start
            assert limited_end <= template_start
            assert limited_start <= target_start
        if limited_end:
            assert limited_end >= origin_end
            if template_end:
                assert limited_end >= template_end
            assert limited_end >= target_end

        array_core = []
        array_extend_left = []
        array_extend_right = []
        array_supply_left = []
        array_supple_right = []

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
                array_extend_right = cls.clip(blocks=template,
                                              start=origin_end,
                                              end=min(template_end, target_end))
            if (template_end is not None) and (target_end > template_end) and supply_right:
                array_supple_right.append((template_end, target_end))

        array = array_supply_left + array_extend_left + array_core + array_extend_right + array_supple_right
        if suture:
            array = cls.suture(array)
        return array

    @classmethod
    def extend(cls, blocks, left=0, right=0, template=None, extend=False, supply=False, trim=False,
               valid=None, suture=False):
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
