# Class: CRange (chromosome range)

from .interval_range import IRange

class CRange(IRange):
    VALID_STRAND_VALUE_LIST = [None, ".", "+", "-"]
    __slots__ = ["chrom", "name", "_strand"]    # Set slot to reduce memory usage.

    def __init__(self, chrom, start, end, name=None, strand=None):
        super(CRange, self).__init__(start, end)
        self.chrom = chrom
        self.name = name
        self.strand = strand
    
    @property
    def strand(self):
        return self._strand
    
    @strand.setter
    def strand(self, value):
        if value not in self.VALID_STRAND_VALUE_LIST:  # Check strand
            raise ValueError("Invalid strand \"%s\". Valid value [None, \"+\", \"-\", \".\"]" % value)
        self._strand = value

    @property
    def is_forward(self):
        return self.strand == "+"

    @property
    def is_reverse(self):
        return self.strand == "-"

    def get_index(self, position, strandness=True):
        if self.start <= position < self.end:
            if strandness and self.is_reverse:
                return self.end - position - 1
            else:
                return position - self.start
        raise ValueError("Position %d is out of range %d-%d." % (position, self.start, self.end))

    def get_position(self, index, strandness=True):
        length = len(self)
        if -length <= index < 0:
            if strandness and self.is_reverse:
                return self.start - index - 1
            else:
                return self.end + index
        elif 0 <= index < length:
            if strandness and self.is_reverse:
                return self.end - index - 1
            else:
                return self.start + index
        raise ValueError("Index %d is out of range %d-%d" % (index, -length, length))

    # overwrite comparison.

    def _compare(self, other):
        if self.chrom < other.chrom:
            return -1
        elif self.chrom == other.chrom:
            return super(CRange, self)._compare(other)
        return 1
    
    def format(self, fmt="BED"):
        fmt = fmt.upper()
        if fmt == "BED":
            return "\t".join(map(str, [self.chrom, self.start, self.end, self.name, ".", self.strand]))
        raise ValueError("Not support format %s." % fmt)
