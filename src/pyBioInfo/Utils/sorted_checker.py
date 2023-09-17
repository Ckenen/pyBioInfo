class SortedChecker(object):
    def __init__(self, container, check_chrom=True, check_start=True, check_end=True):
        self._container = container
        self._last = None
        self._iterator = None
        self._check_chrom = check_chrom
        self._check_start = check_start
        self._check_end = check_end
        assert self._check_chrom
        assert self._check_start

    def __iter__(self):
        for current in self._container:
            if self._last is not None and self._check_chrom:
                if self._last.chrom > current.chrom:
                    raise RuntimeError(
                        "Error in sorted chrom [last: %s, current: %d]" % (self._last.chrom, current.chrom))
                elif self._last.chrom == current.chrom and self._check_start:
                    if self._last.start > current.start:
                        raise RuntimeError(
                            "Error in sorted start [last: %s, current: %d]" % (self._last.chrom, current.chrom))
                    elif self._last.start == current.start and self._check_end:
                        if self._last.end > current.end:
                            raise RuntimeError(
                                "Error in sorted end [last: %s, current: %d]" % (self._last.chrom, current.chrom))

                # if current < self._last:
                #     raise RuntimeError("The items is not sorted! [last: %s, current: %s]" % (self._last, current))
            yield current
            self._last = current

    def __next__(self):
        if self._iterator is None:
            self._iterator = self.__iter__()
        return next(self._iterator)
