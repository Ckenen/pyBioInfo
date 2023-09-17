class ShiftLoader(object):
    def __init__(self, items, history=False):
        self._items = items
        self._iterator = self.__iter__()
        self._buffer = []
        self._reach_tail = False
        self.history = history
        self.histories = []
        self.size = 0
        self._last_chrom = None
        self._last_start = None
        self._last_end = None

    def __iter__(self):
        for obj in self._items:
            yield obj

    def _load_next_item(self):
        obj = None
        if not self._reach_tail:
            try:
                obj = next(self._iterator)
                self.size += 1
                if self.history:
                    self.histories.append(obj)
            except StopIteration:
                self._reach_tail = True
        return obj

    def _check_chrom_start_end(self, chrom, start, end):
        if self._last_chrom is not None:
            if chrom < self._last_chrom:
                raise RuntimeError("%s < %s" % (chrom, self._last_chrom))
            elif chrom == self._last_chrom:
                if start < self._last_start:
                    print(self._last_chrom, self._last_start, self._last_end)
                    print(chrom, start, end)
                    raise RuntimeError()
        self._last_chrom = chrom
        self._last_start = start
        self._last_end = end

    def flush(self):
        while True:
            obj = self._load_next_item()
            if obj is None:
                break


    def fetch(self, chrom=None, start=None, end=None, obj=None):
        if obj is not None:
            chrom = obj.chrom
            start = obj.start
            end = obj.end
        self._check_chrom_start_end(chrom, start, end)
        i = 0
        while True:
            if i >= len(self._buffer):
                if self._reach_tail:
                    break
                else:
                    obj = self._load_next_item()
                    if obj is None:
                        break
                    else:
                        self._buffer.append(obj)
            else:
                obj = self._buffer[i]
                if obj.chrom < chrom:
                    self._buffer.pop(i)
                elif obj.chrom == chrom:
                    if obj.end <= start:
                        self._buffer.pop(i)
                    elif obj.start >= end:
                        break
                    else:
                        yield obj
                        i += 1
                else:
                    break
