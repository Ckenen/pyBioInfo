
class Bundle(object):
    def __init__(self, chrom, start_min, start_max, end_min, end_max, mode, count, data):
        self.chrom = chrom
        self.start_min = start_min
        self.start_max = start_max
        self.end_min = end_min
        self.end_max = end_max
        self.mode = mode
        self.count = count
        self.data = data

    def as_dict(self):
        d = dict()
        d["chrom"] = self.chrom
        d["start_min"] = self.start_min
        d["start_max"] = self.start_max
        d["end_min"] = self.end_min
        d["end_max"] = self.end_max
        d["mode"] = self.mode
        d["count"] = self.count
        return d

    def __str__(self):
        # return "\t".join(map(str, [self.chrom,
        #                            self.start_min,
        #                            self.start_max,
        #                            self.end_min,
        #                            self.end_max,
        #                            self.count]))
        return "chrom: %s, start_min: %d, start_max: %d, end_min: %d, end_max: %d, count: %d" % (
            self.chrom,
            self.start_min,
            self.start_max,
            self.end_min,
            self.end_max,
            self.count
        )


class BundleBuilder(object):
    MODE_START_ONLY = 0
    MODE_START_AND_END = 1

    def __init__(self, objs, min_capacity=1, min_spacing=0, mode=MODE_START_AND_END, keep=False):
        self._objs = objs
        self._min_capacity = min_capacity
        self._min_spacing = min_spacing
        self._mode = mode
        self._keep = keep  #

    def _build_by_start_only_mode(self):
        chrom = None
        start_min = None
        start_max = None
        end_min = None
        end_max = None
        count = None
        data = None
        
        last = None
        
        for obj in self._objs:
            
            if chrom is None:
                chrom = obj.chrom 
                start_min = obj.start
                start_max = obj.start
                end_min = obj.end
                end_max = obj.end
                count = 1
                
                if self._keep:
                    data = [obj]
                    
            elif obj.chrom == chrom:
                
                if obj.start < last.start:
                    raise RuntimeError("The input regions must be sorted by start (%s: %d-%d < %s: %d)" % (
                        obj.start, obj.start, obj.end, start_min))
                    
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
                    
                    if self._keep:
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
            
        if chrom:
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

    def _build_by_start_and_end_mode(self):
        chrom = None
        start_min = None
        start_max = None
        end_min = None
        end_max = None
        count = None
        data = None
        
        last = None
        
        for obj in self._objs:
            if chrom is None:
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
                        obj.start, obj.start, obj.end, start_min))
                    
                if count >= self._min_capacity and obj.start >= end_max + self._min_spacing:
                    bundle = Bundle(
                        chrom=chrom, 
                        start_min=start_min, 
                        start_max=start_max, 
                        end_min=end_min, 
                        end_max=end_max, 
                        mode=self._mode, 
                        count=count, 
                        data =data,
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
                    
                    if self._keep:
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
                    data =data,
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
                
        if chrom:
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

    def __iter__(self):
        if self._mode == self.MODE_START_ONLY:
            iter = self._build_by_start_only_mode()
        elif self._mode == self.MODE_START_AND_END:
            iter = self._build_by_start_and_end_mode()
        else:
            raise RuntimeError("Unknown mode (%s)" % self._mode)
        for bundle in iter:
            yield bundle
