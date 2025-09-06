# import numpy as np
from collections import defaultdict
from functools import cmp_to_key
import pysam
# from .bundle_builder import Bundle, BundleBuilder


class SegmentTools(object):
    @classmethod
    def cmp_for_segments(cls, s1, s2):
        chrom1, chrom2 = s1.reference_name, s2.reference_name
        if chrom1 < chrom2:
            return -1
        elif chrom1 == chrom2:
            start1, start2 = s1.reference_start, s2.reference_start
            if start1 < start2:
                return -1
            elif start1 == start2:
                end1, end2 = s1.reference_end, s2.reference_end
                if end1 < end2:
                    return -1
                elif end1 == end2:
                    name1, name2 = s1.query_name, s2.query_name
                    if name1 < name2:
                        return -1
                    elif name1 == name2:
                        return 0
        return 1
    
    @classmethod
    def sort_segments(cls, segments):
        return sorted(segments, key=cmp_to_key(cls.cmp_for_segments))
    
    @classmethod
    def get_blocks(cls, segment, fill_deletion=True):
        return cls.get_block_from_segment(segment, fill_deletion)    
        # blocks = []
        # x = None  # block start
        # y = None  # block end
        # y0 = segment.reference_start  # last end
        # cigartuples = segment.cigartuples
        # for i in range(len(cigartuples)):
        #     flag, count = cigartuples[i]
        #     if flag == pysam.CMATCH:  # M
        #         x, y = y0, y0 + count
        #         y0 = y
        #         if len(blocks) > 0 and x == blocks[-1][1]:
        #             blocks[-1][1] = y
        #         else:
        #             blocks.append([x, y])
        #     elif flag == pysam.CDEL:  # D
        #         x, y = y0, y0 + count
        #         y0 = y
        #         if fill_deletion:
        #             # if len(blocks) == 0:
        #             #     blocks.append([x, y])
        #             # elif x == blocks[-1][1]:
        #             #     blocks[-1][1] = y
        #             if len(blocks) > 0 and x == blocks[-1][1]:
        #                 blocks[-1][1] = y
        #             else:
        #                 blocks.append([x, y])
        #     elif flag == pysam.CREF_SKIP:  # N
        #         y0 = y0 + count
        #     elif flag == pysam.CINS:  # I
        #         continue
        #     elif flag == pysam.CSOFT_CLIP:  # S
        #         continue
        #     elif flag == pysam.CHARD_CLIP:  # H
        #         continue
        #     else:
        #         raise RuntimeError("Unknown cigar flag: %s" % flag)
        # return blocks

    @classmethod
    def get_mapped_length(cls, segment, fill_deletion=True):
        length = 0
        for flag, count in segment.cigartuples:
            if flag == pysam.CMATCH:  # M
                length += count
            elif flag == pysam.CDEL:  # D
                if fill_deletion:
                    length += count
        return length

    @classmethod
    def get_events(cls, segment, parsed_cigar=None, parsed_md_tag=None, mapped_length=None, report_del=True, report_ins=True):
        """_summary_

        Args:
            segment (_type_): _description_
            parsed_cigar (_type_, optional): _description_. Defaults to None.
            parsed_md_tag (_type_, optional): _description_. Defaults to None.
            mapped_length (_type_, optional): _description_. Defaults to None.

        Returns:
            list: [pos, ref, alt, qua, dis]
        """
        if parsed_cigar is None:
            parsed_cigar = cls.parse_cigar(segment)
        if parsed_md_tag is None:
            parsed_md_tag = cls.parse_md_tag(segment)
        if mapped_length is None:
            mapped_length = cls.get_mapped_length(segment)
        
        sequence = segment.query_sequence
        qualities = segment.query_qualities
        events = []
        # x1, y1 = None, None # mapped index
        # x2, y2 = None, None # read index
        # x3, y3 = None, None # reference index
        # pos = None # position
        # ref = None # reference bases
        # alt = None # alternative bases
        # qua = None # event quality
        for tag in parsed_md_tag:
            if tag[0] == "X": # mismatch
                x1, y1 = tag[1]
                assert y1 - x1 == 1
                dis = min(x1, mapped_length - x1 - 1) # the closest distance to the boundary
                ref = tag[2]
                alt = None
                qua = None
                pos = None
                found = False
                for cigar in parsed_cigar:
                    if cigar[0] == "M" and cigar[1][0] <= x1 and y1 <= cigar[1][1]:
                        offset = x1 - cigar[1][0]
                        x2, y2 = cigar[2]
                        x3, y3 = cigar[3]
                        alt = sequence[x2 + offset]
                        qua = qualities[x2 + offset]
                        pos = x3 + offset
                        events.append(
                            [pos, ref, alt, qua, dis])
                        found = True
                        break
                assert found
            elif tag[0] == "D" and report_del: # deletion
                x1, y1 = tag[1]
                ref = tag[2]
                found = False
                for cigar in parsed_cigar:
                    if cigar[0] == "D" and cigar[1][0] <= x1 and y1 <= cigar[1][1]:
                        assert cigar[1][0] == x1 and cigar[1][1] == y1
                        dis = min(x1, mapped_length - y1)
                        pos = cigar[3][0]
                        events.append([pos, ref, "-", 0, dis])
                        found = True
                        break
                assert found
        if report_ins:
            for cigar in parsed_cigar:
                if cigar[0] == "I": # insertion
                    x1, y1 = cigar[1]
                    x2, y2 = cigar[2]
                    x3, y3 = cigar[3]
                    assert x3 == y3
                    pos = x3
                    ref = "-"
                    dis = min(x1, mapped_length - x1)
                    alt = sequence[x2:y2]
                    if len(alt) == 1:
                        qua = qualities[x2]
                    else:
                        qua = list(qualities[x2:y2])
                    events.append([pos, ref, alt, qua, dis])
        events.sort()
        return events

    @classmethod
    def get_reference_sequence(cls, segment, parsed_cigar=None, parsed_md_tag=None):
        seq = segment.get_reference_sequence
        return seq
        
    @classmethod
    def get_reference_sequence_base(cls, segment):
        counter = defaultdict(int)
        for base in segment.get_reference_sequence():
            counter[base.upper()] += 1
        return counter
    
    @classmethod
    def get_block_from_segment(cls, segment, fill_gap=True):
        blocks = []
        x = None  # block start
        y = None  # block end
        y0 = segment.reference_start  # last end
        for m, n in segment.cigartuples:
            if m == pysam.CMATCH:  # M
                x, y = y0, y0 + n
                y0 = y
                if len(blocks) > 0 and x == blocks[-1][1]:  # concat
                    blocks[-1][1] = y
                else:
                    blocks.append([x, y])
            elif m == pysam.CDEL:  # D
                x, y = y0, y0 + n
                y0 = y
                if fill_gap:
                    if len(blocks) > 0 and x == blocks[-1][1]:
                        blocks[-1][1] = y
                    else:
                        blocks.append([x, y])
            elif m == pysam.CREF_SKIP:  # N
                y0 = y0 + n
            elif m == pysam.CINS:  # I
                continue
            elif m == pysam.CSOFT_CLIP:  # S
                continue
            elif m == pysam.CHARD_CLIP:  # H
                continue
            else:
                raise RuntimeError("Unknown cigar type: %s" % m)
        return blocks

    @classmethod
    def get_aligned_sequence(cls, segment):
        # 把deletion的地方用横杠和-1的质量补上，insertion的地方去除掉
        query_sequence = segment.query_sequence
        query_qualities = segment.query_qualities
        aligned_query_sequence = []
        aligned_query_qualities = []

        offset = 0
        for flag, count in segment.cigartuples:
            if flag == pysam.CHARD_CLIP:
                continue
            elif flag == pysam.CSOFT_CLIP:
                offset = offset + count
            elif flag == pysam.CREF_SKIP:
                continue
            elif flag == pysam.CDEL:
                aligned_query_sequence.append("-" * count)
                aligned_query_qualities.extend([-1] * count)
            elif flag == pysam.CINS:
                offset += count
            elif flag == pysam.CMATCH:
                tmp = offset + count
                aligned_query_sequence.append(query_sequence[offset:tmp])
                aligned_query_qualities.extend(query_qualities[offset:tmp])
                offset = tmp
            else:
                raise RuntimeError(
                    "Unknown cigar type: %s [%d]" % (flag, count))
        aligned_query_sequence = "".join(aligned_query_sequence)
        return aligned_query_sequence, aligned_query_qualities

    @classmethod
    def merge_deletion_gap(cls, items):
        items = items.copy()
        while True:
            found = False
            i = 0
            for item in items:
                if item[0] == pysam.CDEL:
                    found = True
                    break
                i += 1
            if found:
                if i == 0:
                    if i == len(items) - 1:
                        raise ValueError("Only one DEL")
                    else:
                        item = items[i]
                        items.pop(i)
                        if items[i][0] == pysam.CMATCH and items[i][1] == item[2]:
                            items[i][1] = item[1]
                else:
                    if i == len(items) - 1:
                        item = items[i]
                        items.pop(i)
                        if items[-1][0] == pysam.CMATCH and items[-1][2] == item[1]:
                            items[-1][2] = item[2]
                    else:
                        item = items[i]
                        items.pop(i)
                        if items[i - 1][0] == pysam.CMATCH:
                            if items[i][0] == pysam.CMATCH:
                                if items[i - 1][2] == item[1]:
                                    if items[i][1] == item[2]:
                                        item1 = items[i]
                                        items.pop(i)
                                        items[i - 1][2] = item1[2]
                                    else:
                                        items[i - 1][2] = item[2]
                                else:
                                    if items[i][1] == item[2]:
                                        items[i][1] = item[1]
                                    else:
                                        continue
                            else:
                                if items[i - 1][2] == item[1]:
                                    items[i - 1][2] = item[2]
                        else:
                            if items[i][0] == pysam.CMATCH:
                                if items[i][1] == item[2]:
                                    items[i][2] = item[1]
                            else:
                                continue
            else:
                break
        return items

    @classmethod
    def merge_insertion_gap(cls, items):
        items = items.copy()
        return items

    @classmethod
    def parse_cigar(cls, segment):
        """
        return: [
            (
                operation, 
                (mapped_index_start, mapped_index_end), 
                (read_index_start, read_index_end), 
                (ref_position_start, ref_position_end)
            ),
            ...
        ]
        """
        results = []
        rps = segment.reference_start  # ref_position_start
        rpe = rps  # ref_position_end
        mis = 0  # mapped_index_start
        mie = 0  # mapped_index_end
        ris = 0  # read_index_start
        rie = 0  # read_index_end
        for operation, length in segment.cigartuples:
            if operation == pysam.CMATCH:
                rpe = rps + length
                mie = mis + length
                rie = ris + length
                results.append(("M", (mis, mie), (ris, rie), (rps, rpe)))
                rps = rpe
                mis = mie
                ris = rie
            elif operation == pysam.CINS:
                rie = ris + length
                results.append(("I", (mis, mis), (ris, rie), (rps, rpe)))
                ris = rie
            elif operation == pysam.CDEL:
                rpe = rps + length
                mie = mis + length
                results.append(("D", (mis, mie), (ris, ris), (rps, rpe)))
                rps = rpe
                mis = mie
            elif operation == pysam.CREF_SKIP:
                rpe = rps + length
                results.append(("N", (mis, mis), (ris, ris), (rps, rpe)))
                rps = rpe
            elif operation == pysam.CSOFT_CLIP:
                rie = ris + length
                results.append(("S", (mis, mis), (ris, rie), (rps, rpe)))
                ris = rie
            elif operation == pysam.CHARD_CLIP:
                rie = ris + length
                results.append(("H", (mis, mis), (ris, rie), (rps, rpe)))
                ris = rie
                # raise RuntimeError()
            elif operation == pysam.CPAD:
                raise RuntimeError()
            elif operation == pysam.CEQUAL:
                raise RuntimeError()
            elif operation == pysam.CDIFF:
                raise RuntimeError()
            elif operation == pysam.CBACK:
                raise RuntimeError()
            else:
                raise RuntimeError()
        return results

    @classmethod
    def parse_md_tag(cls, segment):
        """_summary_

        Args:
            segment (_type_): _description_

        Raises:
            RuntimeError: _description_

        Returns:
            _type_: _description_
        """
        results = []
        mis = 0  # mapped_index_start
        mie = 0  # mapped_index_end
        md = segment.get_tag("MD")
        i = 0
        while i < len(md):
            if "0" <= md[i] <= "9":
                j = i + 1
                while j < len(md) and "0" <= md[j] <= "9":
                    j += 1
                mie = mis + int(md[i:j])
                results.append(("=", (mis, mie), None))
                mis = mie
                i = j
            elif "A" <= md[i] <= "Z":
                j = i + 1
                while j < len(md) and "A" <= md[j] <= "Z":
                    j += 1
                mie = mis + j - i
                results.append(("X", (mis, mie), md[i:j]))
                mis = mie
                i = j
            elif md[i] == "^":
                j = i + 1
                while j < len(md) and "A" <= md[j] <= "Z":
                    j += 1
                mie = mis + j - i - 1
                results.append(("D", (mis, mie), md[i + 1:j]))
                mis = mie
                i = j
            else:
                raise RuntimeError()
        return results

    @classmethod
    def fetch_snvs(cls, segment=None, parsed_cigar=None, parsed_md_tag=None):
        """[summary]

        Args:
            segment ([type], optional): [description]. Defaults to None.
            parsed_cigar ([type], optional): [description]. Defaults to None.
            parsed_md_tag ([type], optional): [description]. Defaults to None.

        Raises:
            RuntimeError: [description]
            RuntimeError: [description]
            RuntimeError: [description]
            RuntimeError: [description]
            RuntimeError: [description]
        """
        if parsed_cigar is None:
            parsed_cigar = cls.parse_cigar(segment)
        if parsed_md_tag is None:
            parsed_md_tag = cls.parse_md_tag(segment)
        array1 = parsed_cigar
        array2 = parsed_md_tag
        sequence = segment.query_sequence
        qualities = segment.query_qualities
        results = []
        for item1 in array1:
            if item1[0] == "I":
                alt_seq = sequence[item1[2][0]:item1[2][1]]
                scores = qualities[item1[2][0]:item1[2][1]]
                results.append(("I",
                               item1[1], item1[2], item1[3],
                               None, alt_seq, scores))
        sequence = segment.query_sequence
        qualities = segment.query_qualities
        i1 = 0
        i2 = 0
        while i2 < len(array2):
            item2 = array2[i2]
            if item2[0] == "=":
                i2 += 1
                continue
            elif item2[0] == "X":
                assert item2[1][1] - item2[1][0] == 1
                found = False
                while i1 < len(array1):
                    item1 = array1[i1]
                    if item1[1][1] <= item2[1][0]:
                        i1 += 1
                        continue
                    elif item1[1][0] >= item2[1][1]:
                        raise RuntimeError()
                    else:
                        if item1[0] == "M":
                            offset = item2[1][0] - item1[1][0]
                            length = item2[1][1] - item2[1][0]
                            ris = item1[2][0] + offset
                            rie = ris + length
                            rps = item1[3][0] + offset
                            rpe = rps + length
                            alt_seq = sequence[ris:rie]
                            scores = qualities[ris:rie]
                            results.append(("X",
                                           item2[1], (ris, rie), (rps, rpe),
                                           item2[2], alt_seq, scores))
                            found = True
                            break
                        elif item1[0] == "I":
                            i1 += 1
                            continue
                        else:
                            raise RuntimeError()
                assert found
                i2 += 1
            elif item2[0] == "D":
                found = False
                while i1 < len(array1):
                    item1 = array1[i1]
                    if item1[1][1] <= item2[1][0]:
                        i1 += 1
                        continue
                    elif item1[1][0] >= item2[1][1]:
                        raise RuntimeError()
                    else:
                        if item1[0] == "D":
                            assert item1[1][0] == item2[1][0]
                            assert item1[1][1] == item2[1][1]
                            results.append(("D",
                                           item2[1], item1[2], item1[3],
                                           item2[2], None, None))
                            i1 += 1
                            found = True
                            break
                        else:
                            raise RuntimeError()
                assert found
                i2 += 1
            else:
                raise RuntimeError()
        results = list(sorted(results, key=lambda item: item[1]))
        return results

    @classmethod
    def get_mapped_length(cls, segment):
        length = 0
        for cigar in segment.cigartuples:
            if cigar[0] == pysam.CMATCH or cigar[0] == pysam.CDEL:
                length += cigar[1]
        return length
    
    @classmethod
    def get_clipped(cls, segment):
        clip1, clip2 = 0, 0
        cigars = segment.cigartuples
        if cigars[0][0] == pysam.CHARD_CLIP or cigars[0][0] == pysam.CSOFT_CLIP:
            clip1 = cigars[0][1]
        if cigars[-1][0] == pysam.CHARD_CLIP or cigars[-1][0] == pysam.CSOFT_CLIP:
            clip2 = cigars[1][1]
        return clip1, clip2
    
    @classmethod
    def get_query_base(cls, segment, position, parsed_cigar=None):
        if parsed_cigar is None:
            parsed_cigar = cls.parse_cigar(segment)
        base = None
        for k, mapped_idxs, read_idxs, contig_idxs in parsed_cigar:
            if contig_idxs[1] <= position:
                continue
            elif contig_idxs[0] > position:
                break
            else:
                if k == "M":
                    read_idx = position - contig_idxs[0] + read_idxs[0]
                    if parsed_cigar[0][0] == "H":
                        i = read_idx - parsed_cigar[0][2][1]
                    else:
                        i = read_idx
                    base = segment.query_sequence[i]
                elif k == "D":
                    base = "-"
                elif k == "N":
                    base = None
                else:
                    assert False
                break

        # for block in parsed_cigar:
        #     if block[3][0] <= position < block[3][1]:
        #         if block[0] == "M":
        #             offset = position - block[3][0]
        #             read_idx = block[2][0] + offset
        #             base = segment.query_sequence[read_idx]
        #         elif block[0] == "D":
        #             base = "-"
        #         else:
        #             assert False
        #         break

        return base
    
    
    
    
class SegmentPair(object):
    """
    An object that contains segment1 and segment2.
    """

    def __init__(self, mate1, mate2):
        # assert mate1.reference_name == mate2.reference_name
        # assert mate1.is_read1
        # assert mate2.is_read2
        self.chrom = mate1.reference_name
        self.start = min(mate1.reference_start, mate2.reference_start)
        self.end = max(mate1.reference_end, mate2.reference_end)
        self.mate1 = mate1
        self.mate2 = mate2

    def __lt__(self, other):
        if self.chrom < other.chrom:
            return True
        elif self.chrom == other.chrom:
            if self.start < other.start:
                return True
            elif self.start == other.start:
                if self.end < other.end:
                    return True
        return False


class SegmentPairBuilder(object):
    """
    Create SegmentPaired objects.
    """

    def __init__(self, segments):
        self.segments = segments

    @classmethod
    def _build_pair(cls, segments):
        # [s1, s2, s3, s4, ...]
        array1 = []
        # [[s1, s2, s3, s4, ...], [s5, s6, s7, s8, ...]]
        array2 = []
        # [pair1, pair2, ...]
        array3 = []

        # step 1: grouped by query name
        for segment in sorted(segments, key=lambda obj: obj.query_name):
            if len(array1) == 0:
                array1.append(segment)
            else:
                if segment.query_name == array1[0].query_name:
                    array1.append(segment)
                else:
                    array2.append(array1)
                    array1 = [segment]
        if len(array1) > 0:
            array2.append(array1)

        # step 2: build pair
        for array1 in array2:
            array1 = list(sorted(array1, key=lambda obj: obj.reference_start))
            while len(array1) >= 2:
                s1 = array1.pop(0)  # segment 1
                s2 = None  # segment 2
                idx = 0
                for i, tmp in enumerate(array1):
                    if tmp.reference_start < s1.next_reference_start:
                        continue
                    elif tmp.reference_start == s1.next_reference_start:
                        if tmp.next_reference_start == s1.reference_start:
                            idx = i
                            s2 = tmp
                            break
                    else:
                        break
                if s2:
                    array1.pop(idx)
                    if s1.is_read2:
                        s1, s2 = s2, s1
                    pair = SegmentPair(s1, s2)
                    array3.append(pair)

        # step 3: sorting pairs
        array3 = list(sorted(array3))
        for item in array3:
            yield item

    def __iter__(self):
        chrom, require_start = None, None
        count = None
        array = None
        last = None
        for segment in self.segments:
            assert segment.is_proper_pair
            chrom1 = segment.reference_name
            start1 = segment.reference_start
            start2 = segment.next_reference_start
            if chrom is None:
                chrom = chrom1
                require_start = max(start1, start2)
                array = [segment]
                count = 1
            else:
                if chrom1 == chrom:
                    assert start1 >= last.reference_start
                    if count >= 32 and start1 > require_start:
                        for pair in self._build_pair(array):
                            yield pair
                        require_start = max(start1, start2)
                        array = [segment]
                        count = 1
                    else:
                        require_start = max(require_start, start1, start2)
                        array.append(segment)
                        count += 1
                elif chrom1 > chrom:
                    for pair in self._build_pair(array):
                        yield pair
                    chrom = chrom1
                    require_start = max(start1, start2)
                    array = [segment]
                    count = 1
                else:
                    raise RuntimeError()
            last = segment
        if chrom:
            for pair in self._build_pair(array):
                yield pair
