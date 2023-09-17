# Class: IRange (interval range)

class IRange(object):
    # Set slot to reduce memory usage.
    __slots__ = ["_start", "_end"] 

    def __init__(self, start, end):
        if start < 0:
            raise ValueError("The value of start (%d) must be larger or equal to 0." % start)
        if start >= end:
            raise ValueError("The value of start (%d) must be larger than that of end (%d)." % (start, end))
        self._start = start
        self._end = end

    # Read-only access to start and end.

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    def __iter__(self):
        """
        IRange is an iterator.
        The first element is start and the second is end. 
        """
        yield self.start
        yield self.end

    def __getitem__(self, item):
        """
        IRange is an list.
        start = obj[0]
        end = obj[1]
        """
        if item == 0:
            return self.start
        elif item == 1:
            return self.end
        raise RuntimeError("Index 0 for start, 1 for end.")

    def __setitem__(self, key, value):
        raise RuntimeError("Unsupported operations!")

    # string format and length.

    def __str__(self):
        return "%s: %d-%d" % (self.__class__.__name__, self.start, self.end)

    def __len__(self):
        return self.end - self.start
    
    # Position-index system.

    def index(self, position):
        if self.start <= position < self.end:
            return position - self.start
        raise ValueError("Position %d is out of range %d-%d." % (position, self.start, self.end))

    def position(self, index):
        length = len(self)
        if -length <= index < 0:
            return self.end + index
        elif 0 <= index < length:
            return self.start + index
        raise ValueError("Index %d is out of range %d-%d" % (index, -length, length))

    # Comparison

    def _compare(self, other):
        if self.start < other.start:
            return -1
        elif self.start == other.start:
            if self.end < other.end:
                return -1
            elif self.end == other.end:
                return 0
        return 1

    def __lt__(self, other):
        return self._compare(other) < 0

    def __le__(self, other):
        return self._compare(other) <= 0

    def __eq__(self, other):
        return self._compare(other) == 0

    def __gt__(self, other):
        return self._compare(other) > 0

    def __ge__(self, other):
        return self._compare(other) >= 0
    
    # Distance relationship

    def _distance(self, other):
        start = max(self.start, other.start)
        end = min(self.end, other.end)
        return start - end

    def overlap(self, other):
        return self._distance(other) < 0

    def adjacent(self, other):
        return self._distance(other) == 0

    def away(self, other):
        return self._distance(other) > 0

    def contact(self, other):
        return self._distance(other) <= 0

    def contain(self, other):
        return self.start <= other.start and self.end >= other.end

    def coincide(self, other):
        return self.contain(other)

    # Others

    def __and__(self, other):
        raise NotImplementedError()

    def __or__(self, other):
        raise NotImplementedError()

    def __xor__(self, other):
        raise NotImplementedError()

    def __add__(self, other):
        raise NotImplementedError()

    def __sub__(self, other):
        raise NotImplementedError()

    def __left__(self, step):
        raise NotImplementedError()

    def _move(self, step):
        start = self.start + step
        end = self.end + step
        return IRange(start, end)

    def __lshift__(self, step):
        return self._move(-step)
    
    def __rshift__(self, step):
        return self._move(step)
        
    
