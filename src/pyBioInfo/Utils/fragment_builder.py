import functools
from pyBioInfo.Utils import SortedChecker


class FragmentBuilder(object):
    def __init__(self, alignments):
        self.count = 100000  # Read count per round.
        self.alignments = alignments
        self.reads1 = []
        self.reads2 = []
        self.results = []

    @classmethod
    def sort_read1(cls, alignment1, alignment2):
        ret = -1
        chrom1 = alignment1.chrom
        chrom2 = alignment2.chrom
        if chrom1 > chrom2:
            ret = 1
        elif chrom1 == chrom2:
            start1 = alignment1.start
            start2 = alignment2.start
            if start1 > start2:
                ret = 1
            elif start1 == start2:
                name1 = alignment1.name
                name2 = alignment2.name
                if name1 > name2:
                    ret = 1
                elif name1 == name2:
                    ret = 0
        return ret

    @classmethod
    def sort_read2(cls, alignment1, alignment2):
        ret = -1
        segment1 = alignment1.segment
        segment2 = alignment2.segment
        chrom1 = segment1.next_reference_name
        chrom2 = segment2.next_reference_name
        if chrom1 > chrom2:
            ret = 1
        elif chrom1 == chrom2:
            start1 = segment1.next_reference_start
            start2 = segment2.next_reference_start
            if start1 > start2:
                ret = 1
            elif start1 == start2:
                name1 = alignment1.name
                name2 = alignment2.name
                if name1 > name2:
                    ret = 1
                elif name1 == name2:
                    ret = 0
        return ret

    @classmethod
    def _compare_for_mate2(cls, alignment1, alignment2):
        chrom1 = alignment1.chrom
        chrom2 = alignment2.segment.next_reference_name
        if chrom1 < chrom2:
            return -2
        elif chrom1 > chrom2:
            return 2
        else:
            start1 = alignment1.start
            start2 = alignment2.segment.next_reference_start
            if start1 < start2:
                return -2
            elif start1 > start2:
                return 2
            else:
                name1 = alignment1.name
                name2 = alignment2.name
                if name1 < name2:
                    return -2
                elif name1 > name2:
                    return 2
                else:
                    if alignment1.segment.next_reference_name == alignment2.chrom \
                            and alignment1.segment.next_reference_start == alignment2.start:
                        return 0
                    else:
                        return 1

    def export_results(self):
        pass

    def execute_round(self):
        self.reads1 = list(sorted(self.reads1, key=functools.cmp_to_key(self.sort_read1)))
        self.reads2 = list(sorted(self.reads2, key=functools.cmp_to_key(self.sort_read2)))

    def fragments(self):
        count = 0
        for alignment in SortedChecker(self.alignments, check_end=False):
            if alignment.segment.is_proper_pair:
                if alignment.segment.is_read1:
                    self.reads1.append(alignment)
                else:
                    self.reads2.append(alignment)
                count += 1
                if count >= self.count:
                    self.execute_round()
                    count = 0
                    for fragment in self.export_results():
                        yield fragment
        if count > 0:
            self.execute_round()
            count = 0
            for fragment in self.export_results():
                yield fragment



    # def fragments1(self):
    #     buffer1 = None
    #     buffer2 = None
    #     index1 = None
    #     index2 = None
    #     chrom = None
    #     for alignment in SortedChecker(self.alignments, check_end=False):
    #         segment = alignment.segment
    #         if not segment.is_proper_pair:
    #             continue
    #         if chrom is None:
    #             chrom = alignment.chrom
    #             buffer1 = []
    #             buffer2 = []
    #             index1 = 0
    #             index2 = 0
    #             if segment.is_read1:
    #                 buffer1.append(alignment)
    #             else:
    #                 buffer2.append(alignment)
    #         elif alignment.chrom == chrom:
    #             if segment.reference_start > segment.next_reference_start:
    #                 found = False
    #                 if segment.is_read1:
    #                     if len(buffer2) > 0 and buffer2[0].start <= segment.next_reference_start <= buffer2[-1].start:
    #                         for i, obj in enumerate(buffer2):
    #                             segment1 = obj.segment
    #                             if obj.start == segment.next_reference_start:
    #                                 if segment1.next_reference_start == segment.reference_start and obj.name == alignment.name and segment1.is_read2:
    #                                     found = True
    #                                     buffer2.pop(i)
    #                                     pair = [alignment, obj]
    #                                     yield pair
    #                                     break
    #                             elif obj.start > segment.next_reference_start:
    #                                 break
    #                 else:
    #                     if len(buffer1) > 0 and buffer1[0].start <= segment.next_reference_start <= buffer1[-1].start:
    #                         for i, obj in enumerate(buffer1):
    #                             segment1 = obj.segment
    #                             if obj.start == segment.next_reference_start:
    #                                 if segment1.next_reference_start == segment.reference_start and obj.name == alignment.name and segment1.is_read1:
    #                                     found = True
    #                                     buffer1.pop(i)
    #                                     pair = [obj, alignment]
    #                                     yield pair
    #                                     break
    #                             elif obj.start > segment.next_reference_start:
    #                                 break
    #             elif segment.reference_start == segment.next_reference_start:
    #                 found = False
    #                 if segment.is_read1:
    #                     if len(buffer2) > 0 and buffer2[0].start <= segment.next_reference_start <= buffer2[-1].start:
    #                         for i, obj in enumerate(buffer2):
    #                             segment1 = obj.segment
    #                             if obj.start == segment.next_reference_start:
    #                                 if segment1.next_reference_start == segment.reference_start and obj.name == alignment.name and segment1.is_read2:
    #                                     found = True
    #                                     buffer2.pop(i)
    #                                     pair = [alignment, obj]
    #                                     yield pair
    #                                     break
    #                             elif obj.start > segment.next_reference_start:
    #                                 break
    #                 else:
    #                     if len(buffer1) > 0 and buffer1[0].start <= segment.next_reference_start <= buffer1[-1].start:
    #                         for i, obj in enumerate(buffer1):
    #                             segment1 = obj.segment
    #                             if obj.start == segment.next_reference_start:
    #                                 if segment1.next_reference_start == segment.reference_start and obj.name == alignment.name and segment1.is_read1:
    #                                     found = True
    #                                     buffer1.pop(i)
    #                                     pair = [obj, alignment]
    #                                     yield pair
    #                                     break
    #                             elif obj.start > segment.next_reference_start:
    #                                 break
    #                 if not found:
    #                     if segment.is_read1:
    #                         buffer1.append(alignment)
    #                     else:
    #                         buffer2.append(alignment)
    #             else:
    #                 if segment.is_read1:
    #                     buffer1.append(alignment)
    #                 else:
    #                     buffer2.append(alignment)
    #         else:
    #             chrom = alignment.chrom
    #             buffer1 = []
    #             buffer2 = []
    #             index1 = 0
    #             index2 = 0
    #             if segment.is_read1:
    #                 buffer1.append(alignment)
    #             else:
    #                 buffer2.append(alignment)

    def __iter__(self):
        return self.fragments()


if __name__ == '__main__':
    from pyBioInfo.IO.File import BamFile

    path = "/home/chenzonggui/developing/pyBioInfo/911_0H_Input.bam"

    # buffer = []
    # with BamFile(path) as bam:
    #     for i, alignment in enumerate(bam):
    #         if i < 200000:
    #             buffer.append(alignment)
    #         elif i == 200000:
    #             print("Finish!")

    with BamFile(path) as bam:
        i = 0
        for fragment in FragmentBuilder(bam):
            read1, read2 = fragment
            i += 1
            print(i)
