# Class: TRange (thick range)

from .genomic_range import GRange
from pyBioInfo.Utils import BlockTools

class TRange(GRange):
    def __init__(self, thick=None, *args, **kwargs):
        super(TRange, self).__init__(*args, **kwargs)
        if thick is not None:
            thick_start, thick_end = thick
            if thick_start < self.start:
                raise ValueError("The value of thick start[%d] can not smaller than start[%d]." % (thick_start, self.start))
            if thick_end > self.end:
                raise ValueError("The value of thick end[%d] can not larger than end[%d]." % (thick_end, self.end))
            if not BlockTools.is_include(self.blocks, thick_start):
                raise ValueError("The value of thick start[%d] is not included in blocks." % thick_start)
            if not BlockTools.is_include(self.blocks, thick_end - 1):
                raise ValueError("The value of thick end[%d] is not included in blocks." % thick_end)
            thick = (thick_start, thick_end)
        self._thick = thick
        self._cache_forward_indexes = None
        self._cache_reverse_indexes = None

    @property
    def thick(self):
        return self._thick

    def __str__(self):
        if self._thick:
            tstr = "%d-%d" % (self._thick[0], self._thick[1])
        else:
            tstr = "None"
        return "%s: %d-%d [strand: %s, block count: %d, thick: %s]" % (self.__class__.__name__, self._start, self._end, self.strand, len(self._blocks), tstr)

    def indexes(self, strandness=True):
        # indexes = None
        if self.is_reverse and strandness:
            if self._cache_reverse_indexes is None:
                i1 = 0
                i2 = None
                i3 = None
                i4 = len(self)
                thick = self.thick
                if thick is not None:
                    thick_start, thick_end = thick
                    i2 = self.get_index(thick_start)
                    i3 = self.get_index(thick_end - 1)
                    i2, i3 = min(i2, i3), max(i2, i3) + 1
                self._cache_reverse_indexes = (i1, i2, i3, i4)
            indexes = self._cache_reverse_indexes
        else:
            if self._cache_forward_indexes is None:
                i1 = 0
                i2 = None
                i3 = None
                i4 = len(self)
                thick = self.thick
                if thick is not None:
                    thick_start, thick_end = thick
                    i2 = self.get_index(thick_start, strandness=False)
                    i3 = self.get_index(thick_end - 1, strandness=False)
                    i2, i3 = min(i2, i3), max(i2, i3) + 1
                self._cache_forward_indexes = (i1, i2, i3, i4)
            indexes = self._cache_forward_indexes
        return indexes

    def format(self, fmt="bed", **kwargs):
        fmt = fmt.lower()
        line = None
        if fmt == "bed":
            ncol = kwargs.get("ncol", 12)
            cols = []
            if ncol >= 1:
                cols.append(self.chrom)
            if ncol >= 2:
                cols.append(self.start)
            if ncol >= 3:
                cols.append(self.end)
            if ncol >= 4:
                name = self.name
                if name is None:
                    name = "None"
                cols.append(name)
            if ncol >= 5:
                try:
                    score = self.score
                except AttributeError:
                    score = "."
                if score is None:
                    score = "."
                cols.append(score)
            if ncol >= 6:
                cols.append(self.strand)
            if ncol >= 7:
                thick_start = self.thick[0] if self.thick is not None else self.start
                cols.append(thick_start)
            if ncol >= 8:
                thick_end = self.thick[1] if self.thick is not None else self.start
                cols.append(thick_end)
            if ncol >= 9:
                try:
                    color = self.color
                except AttributeError:
                    color = "0,0,255" if self.is_forward else "255,0,0"
                cols.append(color)
            if ncol >= 10:
                cols.append(len(self.blocks))
            if ncol >= 11:
                block_sizes = [block_end - block_start for block_start, block_end in self.blocks]
                block_sizes_string = ",".join(map(str, block_sizes)) + ","
                cols.append(block_sizes_string)
            if ncol >= 12:
                block_offsets = [block_start - self.start for block_start, block_end in self.blocks]
                block_offsets_string = ",".join(map(str, block_offsets)) + ","
                cols.append(block_offsets_string)
            line = "\t".join(map(str, cols))
        elif fmt == "gff":
            cols = [self.chrom]
            try:
                source = self.source
            except AttributeError:
                source = "."
            cols.append(source)
            try:
                feature = self.feature
            except AttributeError:
                feature = "."
            cols.append(feature)

            cols.append(self.start + 1)
            cols.append(self.end)
            try:
                score = self.score
            except AttributeError:
                score = "."
            cols.append(score)

            cols.append(self.strand)
            try:
                frame = self.frame
            except AttributeError:
                frame = "."
            cols.append(frame)
            items = kwargs.get("attris")
            if items is None:
                items = ["name"]
            else:
                items = [item for item in items]
            if "name" not in items:
                items.insert(0, "name")
            array = []
            for key in items:
                if key == "name":
                    value = self.name
                else:
                    value = self[key]
                array.append("%s=%s" % (key, str(value)))
            cols.append("; ".join(array))
            line = "\t".join(map(str, cols))
        elif fmt == "tbl":
            header = kwargs.get("header")
            if header is None:
                raise RuntimeError()
            cols = []
            for key in header:
                if key == "thick":
                    thick = self.thick
                    if thick is None:
                        cols.append("")
                    else:
                        thick_start, thick_end = thick
                        cols.append("%d-%d" % (thick_start, thick_end))
                elif key == "blocks":
                    temps = []
                    for x, y in self.blocks:
                        temps.append("%d-%d" % (x, y))
                    cols.append(";".join(temps))
                else:
                    try:
                        value = self.__getattribute__(key)
                    except AttributeError:
                        value = self[key]
                    cols.append(value)
            line = "\t".join(map(str, cols))
        else:
            raise RuntimeError("Unsupported format %s" % fmt)
        return line
