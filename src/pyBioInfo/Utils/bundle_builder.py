# Bundle: a group of neighboring genomic regions.
from typing import TYPE_CHECKING, ClassVar

import sys
from collections.abc import Generator, Iterable

if TYPE_CHECKING:
    # Imported lazily to break the import cycle
    # pyBioInfo.Range.genomic_range -> pyBioInfo.Utils -> this module.
    from pyBioInfo.Range import GRange

if sys.version_info >= (3, 12):
    from typing import override
else:
    def override(func):
        return func


class Bundle(object):
    """An aggregated group of neighboring genomic regions (produced by BundleBuilder)."""

    # Chromosome name.
    chrom: str
    # Smallest start within the group (inclusive).
    start_min: int
    # Largest start within the group (inclusive).
    start_max: int
    # Smallest end within the group (exclusive).
    end_min: int
    # Largest end within the group (exclusive).
    end_max: int
    # Mode that produced this Bundle (BundleBuilder.MODE_*).
    mode: int
    # Number of regions within the group.
    count: int
    # Original regions of the group; kept only when keep=True,
    # otherwise None.
    data: "list[GRange] | None"

    def __init__(self, chrom: str, start_min: int, start_max: int,
                 end_min: int, end_max: int, mode: int, count: int,
                 data: "list[GRange] | None") -> None:
        """
        Initialize a Bundle.

        :param chrom: Chromosome name.
        :param start_min: Smallest start within the group.
        :param start_max: Largest start within the group.
        :param end_min: Smallest end within the group.
        :param end_max: Largest end within the group.
        :param mode: Aggregation mode (BundleBuilder.MODE_START_ONLY
                     or MODE_START_AND_END).
        :param count: Number of regions within the group.
        :param data: Original regions of the group (None when keep=False).
        """
        self.chrom = chrom
        self.start_min = start_min
        self.start_max = start_max
        self.end_min = end_min
        self.end_max = end_max
        self.mode = mode
        self.count = count
        self.data = data

    def as_dict(self) -> "dict[str, str | int]":
        """Convert the Bundle to a dictionary (excluding data)."""
        d: "dict[str, str | int]" = dict()
        d["chrom"] = self.chrom
        d["start_min"] = self.start_min
        d["start_max"] = self.start_max
        d["end_min"] = self.end_min
        d["end_max"] = self.end_max
        d["mode"] = self.mode
        d["count"] = self.count
        return d

    @override
    def __str__(self) -> str:
        """Return a human-readable description string."""
        return "chrom: %s, start_min: %d, start_max: %d, end_min: %d, end_max: %d, count: %d" % (
            self.chrom,
            self.start_min,
            self.start_max,
            self.end_min,
            self.end_max,
            self.count
        )


class BundleBuilder(object):
    """Aggregate a stream of regions (sorted by chrom and start) into neighboring Bundles."""

    # Aggregate by start only: a new bundle is opened when the start of
    # an incoming region exceeds the largest start within the group by
    # more than min_spacing.
    MODE_START_ONLY: ClassVar[int] = 0
    # Aggregate by start and end: a new bundle is opened when the start
    # of an incoming region is at least end_max + min_spacing (regions
    # within a bundle are thus guaranteed not to overlap).
    MODE_START_AND_END: ClassVar[int] = 1

    # Iterable of regions to aggregate.
    _objs: "Iterable[GRange]"
    # Minimum number of regions to form a bundle.
    _min_capacity: int
    # Spacing threshold between neighboring regions.
    _min_spacing: int
    # Aggregation mode.
    _mode: int
    # Whether to keep the original region objects.
    _keep: bool

    def __init__(self, objs: "Iterable[GRange]", min_capacity: int = 1,
                 min_spacing: int = 0, mode: int = MODE_START_AND_END,
                 keep: bool = False) -> None:
        """
        Initialize a BundleBuilder.

        :param objs: Iterable of regions sorted by chrom and start.
        :param min_capacity: Minimum number of regions to form a bundle;
                             smaller remainders are merged into the previous bundle.
        :param min_spacing: Spacing threshold between neighboring regions;
                            a larger spacing opens a new bundle.
        :param mode: Aggregation mode (MODE_START_ONLY or MODE_START_AND_END).
        :param keep: Whether to keep the original region objects in Bundle.data.
        :raises ValueError: If mode is not a known mode.
        """
        self._objs = objs
        self._min_capacity = min_capacity
        self._min_spacing = min_spacing
        self._mode = mode
        self._keep = keep  # whether to keep the original objects

    def _build_by_start_only_mode(self) -> "Generator[Bundle, None, None]":
        """Aggregate in MODE_START_ONLY mode: new bundles are opened relative to the largest start of the group."""
        chrom: "str | None" = None
        start_min: int = 0
        start_max: int = 0
        end_min: int = 0
        end_max: int = 0
        count: int = 0
        data: "list[GRange] | None" = None

        last: "GRange | None" = None

        for obj in self._objs:

            if chrom is None or last is None:
                # The first region
                chrom = obj.chrom
                start_min = obj.start
                start_max = obj.start
                end_min = obj.end
                end_max = obj.end
                count = 1

                if self._keep:
                    data = [obj]

            elif obj.chrom == chrom:
                # Equal chromosome

                if obj.start < last.start:
                    raise RuntimeError("The input regions must be sorted by start (%s: %d-%d < %s: %d)" % (
                        obj.chrom, obj.start, obj.end, last.chrom, last.start))

                if count >= self._min_capacity and obj.start > start_max + self._min_spacing:
                    bundle = Bundle(
                        chrom=chrom,
                        start_min=start_min,
                        start_max=start_max,
                        end_min=end_min,
                        end_max=end_max,
                        mode=self._mode,
                        count=count,
                        data=data,
                    )
                    yield bundle

                    chrom = obj.chrom
                    start_min = obj.start
                    start_max = obj.start
                    end_min = obj.end
                    end_max = obj.end
                    count = 1

                    if self._keep:
                        data = [obj]

                else:
                    start_max = obj.start
                    end_min = min(end_min, obj.end)
                    end_max = max(end_max, obj.end)
                    count += 1

                    if self._keep and data is not None:
                        data.append(obj)

            elif obj.chrom > chrom:
                bundle = Bundle(
                    chrom=chrom,
                    start_min=start_min,
                    start_max=start_max,
                    end_min=end_min,
                    end_max=end_max,
                    mode=self._mode,
                    count=count,
                    data=data,
                )
                yield bundle

                chrom = obj.chrom
                start_min = obj.start
                start_max = obj.start
                end_min = obj.end
                end_max = obj.end
                count = 1

                if self._keep:
                    data = [obj]

            else:
                raise RuntimeError("The input regions must be sorted by chrom and start (%s: %d-%d < %s)" % (
                    obj.chrom, obj.start, obj.end, chrom))

            last = obj

        if chrom is not None:
            bundle = Bundle(
                chrom=chrom,
                start_min=start_min,
                start_max=start_max,
                end_min=end_min,
                end_max=end_max,
                mode=self._mode,
                count=count,
                data=data,
            )
            yield bundle

    def _build_by_start_and_end_mode(self) -> "Generator[Bundle, None, None]":
        """Aggregate in MODE_START_AND_END mode: new bundles are opened relative to the largest end of the group."""
        chrom: "str | None" = None
        start_min: int = 0
        start_max: int = 0
        end_min: int = 0
        end_max: int = 0
        count: int = 0
        data: "list[GRange] | None" = None

        last: "GRange | None" = None

        for obj in self._objs:
            if chrom is None or last is None:
                # The first region

                chrom = obj.chrom
                start_min = obj.start
                start_max = obj.start
                end_min = obj.end
                end_max = obj.end
                count = 1
                if self._keep:
                    data = [obj]

            elif obj.chrom == chrom:
                # Equal chromosome

                if obj.start < last.start:
                    raise RuntimeError("The input regions must be sorted by start (%s: %d-%d < %s: %d)" % (
                        obj.chrom, obj.start, obj.end, last.chrom, last.start))

                if count >= self._min_capacity and obj.start >= end_max + self._min_spacing:
                    bundle = Bundle(
                        chrom=chrom,
                        start_min=start_min,
                        start_max=start_max,
                        end_min=end_min,
                        end_max=end_max,
                        mode=self._mode,
                        count=count,
                        data=data,
                    )
                    yield bundle

                    chrom = obj.chrom
                    start_min = obj.start
                    start_max = obj.start
                    end_min = obj.end
                    end_max = obj.end
                    count = 1

                    if self._keep:
                        data = [obj]

                else:
                    start_max = obj.start
                    end_min = min(end_min, obj.end)
                    end_max = max(end_max, obj.end)
                    count += 1

                    if self._keep and data is not None:
                        data.append(obj)

            elif obj.chrom > chrom:
                bundle = Bundle(
                    chrom=chrom,
                    start_min=start_min,
                    start_max=start_max,
                    end_min=end_min,
                    end_max=end_max,
                    mode=self._mode,
                    count=count,
                    data=data,
                )
                yield bundle

                chrom = obj.chrom
                start_min = obj.start
                start_max = obj.start
                end_min = obj.end
                end_max = obj.end
                count = 1

                if self._keep:
                    data = [obj]

            else:
                raise RuntimeError("The input regions must be sorted by chrom and start (%s: %d-%d < %s)" % (
                    obj.chrom, obj.start, obj.end, chrom))

            last = obj

        if chrom is not None:
            bundle = Bundle(
                chrom=chrom,
                start_min=start_min,
                start_max=start_max,
                end_min=end_min,
                end_max=end_max,
                mode=self._mode,
                count=count,
                data=data,
            )
            yield bundle

    def __iter__(self) -> "Generator[Bundle, None, None]":
        """Iterate over Bundles produced by the configured mode."""
        if self._mode == self.MODE_START_ONLY:
            builder = self._build_by_start_only_mode()
        elif self._mode == self.MODE_START_AND_END:
            builder = self._build_by_start_and_end_mode()
        else:
            raise RuntimeError("Unknown mode (%s)" % self._mode)
        for bundle in builder:
            yield bundle
