import functools


class BlockTools(object):
    @classmethod
    def length(cls, blocks, check=True):
        if check:
            cls.check_valid_blocks(blocks)
        return sum([end - start for start, end in blocks])

    @classmethod
    def move(cls, blocks, step=0, check=True):
        if check:
            cls.check_valid_blocks(blocks)
        array = []
        for start, end in blocks:
            array.append((start + step, end + step))
        return array

    @classmethod
    def index(cls, blocks, position, check=True):
        if check:
            cls.check_valid_blocks(blocks)
        offset = 0
        for start, end in blocks:
            if end <= position:
                offset += end - start
            elif start <= position < end:
                return position - start + offset
            elif position < start:
                break
        raise ValueError("Position (%d) not in blocks." % position)

    @classmethod
    def is_include(cls, blocks, position, check=True):
        if check:
            cls.check_valid_blocks(blocks)
        for start, end in blocks:
            if position >= end:
                continue
            elif start <= position < end:
                return True
            elif position < start:
                break
        return False

    @classmethod
    def position(cls, blocks, index, check=True, length=None):
        if check:
            cls.check_valid_blocks(blocks)
        if length is None:
            length = cls.length(blocks, check=False)
        if index >= length or index < -length:
            raise ValueError("Index (%d) out of blocks range (%d <= index < %d)!" % (index, -length, length))
        if index < 0:
            index += length
        idx2 = 0
        for start, end in blocks:
            idx1 = idx2
            idx2 = idx1 + end - start
            if idx1 <= index < idx2:
                return index - idx1 + start
        raise RuntimeError("Unknown error.")

    @classmethod
    def gaps(cls, blocks, check=True):
        if check:
            cls.check_valid_blocks(blocks)
        array = []
        last_start = None
        last_end = None
        for start, end in blocks:
            if last_start is not None:
                array.append((last_end, start))
            last_start = start
            last_end = end
        return array

    @classmethod
    def suture(cls, blocks, gap=0, check=True):
        if check:
            cls.check_sorted(blocks)
        array = []
        last_start = None
        last_end = None
        for start, end in blocks:
            if last_start is None:
                last_start = start
                last_end = end
            else:
                if start - last_end <= gap:
                    last_start = min(last_start, start)
                    last_end = max(last_end, end)
                else:
                    array.append((last_start, last_end))
                    last_start = start
                    last_end = end
        if last_start is not None:
            array.append((last_start, last_end))
        return array

    @classmethod
    def fusion(cls, blocks1, blocks2, check=True):
        blocks = []
        for block in blocks1:
            blocks.append(block)
        for block in blocks2:
            blocks.append(block)
        blocks = list(sorted(blocks, key=lambda item: item[0]))
        return cls.suture(blocks, 0, check)
        

    @classmethod
    def contain(cls, blocks1, blocks2, check=True):
        if check:
            cls.check_valid_blocks(blocks1)
            cls.check_valid_blocks(blocks2)
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
    def coincide(cls, blocks1, blocks2, check=True):
        if check:
            cls.check_valid_blocks(blocks1)
            cls.check_valid_blocks(blocks2)
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

        # flag = True
        #
        # i = 0
        # j = 0
        # while i < len(blocks2):
        #     x2, y2 = blocks2[i]
        #
        #     while True:
        #         if j >= len(blocks1):
        #             flag = False
        #             break
        #         x1, y1 = blocks1[j]
        #
        #         if x1 < x2:
        #             if x2 >= y1:
        #                 j += 1
        #                 continue
        #             else:
        #                 if y1 < y2:
        #                     flag = False
        #                     break
        #                 elif y1 == y2:
        #                     if i == 0:
        #                         j += 1
        #                         break
        #                     else:
        #                         flag = False
        #                         break
        #                 else:
        #                     if i == 0 and len(blocks2) == 1:
        #                         break
        #                     else:
        #                         flag = False
        #                         break
        #         elif x1 == x2:
        #             if y1 < y2:
        #                 flag = False
        #                 break
        #             elif y1 == y2:
        #                 j += 1
        #                 break
        #             else:
        #                 if i == len(blocks2) - 1:
        #                     break
        #                 else:
        #                     flag = False
        #                     break
        #         else:
        #             flag = False
        #             break
        #
        #     if not flag:
        #         break
        #
        #     i += 1
        #
        # return flag

    @classmethod
    def clip(cls, blocks, start, end, template=None, check=True,
             extend=False, extend_left=False, extend_right=False,
             supply=False, supply_left=False, supply_right=False,
             limited=None, suture=False):
        # check the parameter.
        if start < 0:
            raise ValueError()

        if start >= end:
            raise ValueError()

        if check:
            cls.check_valid_blocks(blocks)
            if template is not None:
                cls.check_valid_blocks(blocks)

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
    def _sorted_compare(cls, block1, block2):
        block_start1, block_end1 = block1
        block_start2, block_end2 = block2
        if block_start1 < block_start2:
            return -1
        elif block_start1 == block_start2:
            if block_end1 < block_end2:
                return -1
            elif block_end1 == block_end2:
                return 0
        return 1

    @classmethod
    def sorted(cls, blocks):
        return list(sorted(blocks, key=functools.cmp_to_key(cls._sorted_compare)))

    @classmethod
    def is_sorted(cls, blocks):
        last_start = None
        last_end = None
        for start, end in blocks:
            if last_start is not None:
                if (start < last_start) or (start == last_start and end < last_end):
                    return False
            last_start = start
            last_end = end
        return True

    @classmethod
    def check_sorted(cls, blocks):
        last_start = None
        last_end = None
        for start, end in blocks:
            if last_start is not None:
                if (start < last_start) or (start == last_start and end < last_end):
                    raise ValueError("The blocks is unsorted!")
            last_start = start
            last_end = end

    @classmethod
    def check_valid_blocks(cls, blocks):
        last_end = 0
        for start, end in blocks:
            if start < 0:
                raise ValueError("The value of start (%d) must be larger or equal to 0." % (start))
            if start >= end:
                raise ValueError("The value of start (%d) must be lower than that of the end (%d)." % (start, end))
            if start < last_end:
                raise ValueError("The value of start (%d) must be larger or equal to the end (%d) of the previous block." % (start, last_end))
            last_end = end

    @classmethod
    def is_valid_blocks(cls, blocks):
        try:
            cls.check_valid_blocks(blocks)
            return True
        except ValueError:
            return False

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
