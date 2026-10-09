# Class: SegmentPair, SegmentPairBuilder

from collections import defaultdict
from collections.abc import Iterable, Generator
import pysam
from .segment_tools import SegmentTools

class SegmentPair(object):
    """
    An object that contains segment1 and segment2.
    """

    # Chromosome name of both mates.
    chrom: str
    # Start of the earlier mate (0-based, inclusive).
    start: int
    # End of the later mate (exclusive).
    end: int
    # The read1 segment.
    mate1: pysam.AlignedSegment
    # The read2 segment.
    mate2: pysam.AlignedSegment

    def __init__(self, mate1: pysam.AlignedSegment, mate2: pysam.AlignedSegment) -> None:
        """
        Initialize a SegmentPair.

        :param mate1: The read1 segment.
        :param mate2: The read2 segment.
        """
        if mate1.reference_name is None:
            raise RuntimeError(f"Segment {mate1.query_name} has no chromosome")
        if mate2.reference_name is None:
            raise RuntimeError(f"Segment {mate2.query_name} has no chromosome")
        if mate1.reference_name != mate2.reference_name:
            raise ValueError(f"Segments {mate1.query_name} and {mate2.query_name} are not on the same chromosome")
        if mate1.query_name != mate2.query_name:
            raise ValueError(f"Segments {mate1.query_name} and {mate2.query_name} are not from the same template")
        if not mate1.is_read1 or not mate2.is_read2:
            raise ValueError(f"Segments {mate1.query_name} and {mate2.query_name} are not read1 and read2")
        if mate1.reference_end is None:
            raise RuntimeError(f"Segment {mate1.query_name} has no end")
        if mate2.reference_end is None:
            raise RuntimeError(f"Segment {mate2.query_name} has no end")
        self.chrom = mate1.reference_name
        self.start = min(mate1.reference_start, mate2.reference_start)
        self.end = max(mate1.reference_end, mate2.reference_end)
        self.mate1 = mate1
        self.mate2 = mate2

    def __lt__(self, other: "SegmentPair") -> bool:
        """Order pairs by chrom, start and end."""
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

    # Whether the input stream is sorted by genomic coordinates.
    is_sorted_by_coords: bool
    # Whether the input stream is sorted by query name.
    is_sorted_by_name: bool
    # Iterable of input segments.
    segments: "Iterable[pysam.AlignedSegment]"
    # Number of unmapped segments seen.
    number_of_unmapped: int
    # Total number of segments seen.
    number_of_reads: int
    # Number of read1 segments seen.
    number_of_read1: int
    # Number of read2 segments seen.
    number_of_read2: int
    # Number of pairs produced.
    number_of_pairs: int

    def __init__(self, segments: "Iterable[pysam.AlignedSegment]",
                 is_sorted_by_coords: bool = False, is_sorted_by_name: bool = False) -> None:
        """
        Initialize a SegmentPairBuilder.

        :param segments: Iterable of input segments.
        :param is_sorted_by_coords: The input is sorted by genomic
                                    coordinates.
        :param is_sorted_by_name: The input is sorted by query name.
        :raises RuntimeError: If both sorting flags are set.
        """
        if is_sorted_by_coords and is_sorted_by_name:
            raise RuntimeError("Input cannot be sorted by both coordinates and name")
        self.is_sorted_by_coords = is_sorted_by_coords
        self.is_sorted_by_name = is_sorted_by_name
        self.segments = segments

        self.number_of_unmapped = 0
        self.number_of_reads = 0
        self.number_of_read1 = 0
        self.number_of_read2 = 0
        self.number_of_pairs = 0

    def _build_pairs_from_coords_sorted_segments(self) -> "Generator[SegmentPair, None, None]":
        """
        Build pairs from segments sorted by genomic coordinates.

        Segments are buffered until the stream moves past the current
        window, then paired within the buffer.
        """
        assert self.is_sorted_by_coords
        assert not self.is_sorted_by_name

        chrom_last: str | None = None
        chrom_current: str | None = None
        start_current: int = -1
        start_max: int = -1
        start_require: int = -1
        segments_tmp: list[pysam.AlignedSegment] = []
        segment_last: pysam.AlignedSegment | None = None

        for segment in self.segments:
            if not segment.is_proper_pair:
                raise ValueError(f"Segment {segment.query_name} is not a proper pair")
            chrom_current = segment.reference_name
            if chrom_current is None:
                raise RuntimeError(f"Segment {segment.query_name} has no chromosome")
            start_current = segment.reference_start
            start_max = max(start_current, segment.next_reference_start)

            if chrom_last is not None:
                assert chrom_current >= chrom_last
                if chrom_current == chrom_last:
                    assert segment_last is not None
                    assert start_current >= segment_last.reference_start

                if chrom_current != chrom_last or start_current > start_require:
                    for pair in sorted(self._build_pairs_from_unsorted_segments(segments_tmp)):
                        yield pair
                    segments_tmp = []
                    start_require = -1
                    segment_last = None

            chrom_last = chrom_current
            start_require = max(start_require, start_max)
            segments_tmp.append(segment)
            segment_last = segment

        if chrom_last is not None:
            for pair in sorted(self._build_pairs_from_unsorted_segments(segments_tmp)):
                yield pair

    def _build_pairs_from_name_sorted_segments(self) -> "Generator[SegmentPair, None, None]":
        """Build pairs from segments sorted by query name."""
        assert self.is_sorted_by_name
        assert not self.is_sorted_by_coords

        for segments2 in SegmentTools.load_same_name_segments(self.segments):
            for pair in self._build_pairs_from_same_name_segments(segments2):
                yield pair

    def _build_pairs_from_unsorted_segments(self, segments: "Iterable[pysam.AlignedSegment] | None" = None) -> "Generator[SegmentPair, None, None]":
        """
        Build pairs from a stream without any sorting guarantee.

        Segments are grouped by query name; each group is paired
        individually and the groups are processed in query-name order.
        """
        
        if segments is None:
            assert not self.is_sorted_by_name
            assert not self.is_sorted_by_coords
            segments = self.segments

        count = 0
        data: "defaultdict[str, list[pysam.AlignedSegment]]" = defaultdict(list)
        for segment in segments:
            name = segment.query_name
            if name is None:
                raise RuntimeError(f"Segment {segment.query_name} has no query name")
            count += 1
            data[name].append(segment)
        for _, items in sorted(data.items()):
            for pair in self._build_pairs_from_same_name_segments(items):
                yield pair
        # print(count)

    def _build_pairs_from_same_name_segments(self, segments: "Iterable[pysam.AlignedSegment]") -> "Generator[SegmentPair, None, None]":
        """
        Build pairs from a list of segments sharing one query name.

        read1 segments are matched to read2 segments via ``is_mate``.
        """
        pairs: "list[SegmentPair]" = []
        reads1: "list[pysam.AlignedSegment]" = []
        reads2: "list[pysam.AlignedSegment]" = []
        query_name: "str | None"= None

        for s in segments:
            if query_name is None:
                query_name = s.query_name
            else:
                assert s.query_name == query_name

            self.number_of_reads += 1
            if s.is_unmapped:
                self.number_of_unmapped += 1
                continue
            if s.is_read1:
                self.number_of_read1 += 1
                reads1.append(s)
            elif s.is_read2:
                self.number_of_read2 += 1
                reads2.append(s)
            else:
                raise RuntimeError()

        for s1 in reads1:
            s2: "pysam.AlignedSegment | None" = None
            for i2 in range(len(reads2)):
                if SegmentTools.is_mate(s1, reads2[i2]):
                    s2 = reads2.pop(i2)
                    break

            if s2 is None:
                raise RuntimeError("Can't find mate 2 for segment: %s" % str(s1))

            pairs.append(SegmentPair(s1, s2))
            self.number_of_pairs += 1

        if len(reads2) > 0:
            raise RuntimeError("Can't find mate 1 for segment: %s" % str(reads2[0]))

        for pair in sorted(pairs):
            yield pair
            

    def __iter__(self) -> "Generator[SegmentPair, None, None]":
        """Iterate over the pairs produced by the configured mode."""
        iterator: "Generator[SegmentPair, None, None]"

        if self.is_sorted_by_coords:
            if self.is_sorted_by_name:
                raise RuntimeError()
            else:
                iterator = self._build_pairs_from_coords_sorted_segments()
        else:
            if self.is_sorted_by_name:
                iterator = self._build_pairs_from_name_sorted_segments()
            else:
                iterator = self._build_pairs_from_unsorted_segments()

        for pair in iterator:
            yield pair
