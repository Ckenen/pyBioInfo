from pyBioInfo.SpecialRange import MismatchSite
from .alignment_mismatch_event import AlignmentMismatchEventFactory


class FragmentMismatchEvent(MismatchSite):
    def __init__(self, chrom, start, end, name, strand, score, distance, rbase, tbase, fragment):
        super(FragmentMismatchEvent, self).__init__(chrom=chrom, 
                                                    start=start, 
                                                    end=end,
                                                    name=name, 
                                                    strand=strand,
                                                    rbase=rbase, 
                                                    tbase=tbase)
        self.score = score
        self.distance = distance
        self.fragment = fragment


class FragmentMismatchEventFactory(object):
    @classmethod
    def _integrate_decision_maker_single(cls, fragment, event, mate):
        keep = True
        if mate.start <= event.start < mate.end:
            try:
                qualities = mate.get_aligned_sequences()[1]
                score = qualities[mate.index(event.start, strandness=False)]
                keep = event.score > score
            except ValueError:
                # 突变的位点在另外一条read的start和end之内，但不在block里面
                keep = True
                # print(fragment.format("BED"))
                # raise RuntimeError("Unknown error!")
        if keep:
            obj = FragmentMismatchEvent(
                chrom=event.chrom, start=event.start, end=event.end, name=fragment.name,
                strand=fragment.strand, score=event.score, distance=event.distance,
                rbase=event.rbase, tbase=event.tbase, fragment=fragment)
            return obj
        return None

    @classmethod
    def _integrate_decision_maker_double(cls, fragment, event1, event2):
        if event1.rbase == event2.rbase:
            if event1.tbase == event2.tbase:
                obj = FragmentMismatchEvent(chrom=event1.chrom,
                                            start=event1.start,
                                            end=event1.end,
                                            name=fragment.name,
                                            strand=fragment.strand, 
                                            score=max(event1.score, event2.score),
                                            distance=max(event1.distance, event2.distance),
                                            rbase=event1.rbase, 
                                            tbase=event1.tbase,
                                            fragment=fragment)
                return obj
            else:
                event = None
                if event1.score > event2.score:
                    event = event1
                elif event1.score < event2.score:
                    event = event2
                else:
                    if event1.distance > event2.distance:
                        event = event1
                    elif event1.distance < event2.distance:
                        event = event2
                    else:
                        return None
                obj = FragmentMismatchEvent(chrom=event.chrom,
                                            start=event.start,
                                            end=event.end,
                                            name=fragment.name,
                                            strand=fragment.strand, 
                                            score=event.score,
                                            distance=event.distance,
                                            rbase=event.rbase, 
                                            tbase=event.tbase,
                                            fragment=fragment)
                return obj
        else:
            raise RuntimeError("Unknown error!")
        # return None

    @classmethod
    def _integrate_fragment_mismatch_event(cls, fragment, array1, array2):
        i1 = 0
        i2 = 0
        results = []
        while True:
            if i1 >= len(array1):
                if i2 >= len(array2):
                    break
                else:
                    results.append(cls._integrate_decision_maker_single(
                        fragment, array2[i2], fragment.mate1))
                    i2 += 1
            else:
                if i2 >= len(array2):
                    results.append(cls._integrate_decision_maker_single(
                        fragment, array1[i1], fragment.mate2))
                    i1 += 1
                else:
                    event1, event2 = array1[i1], array2[i2]
                    if event1.start < event2.start:
                        results.append(cls._integrate_decision_maker_single(
                            fragment, event1, fragment.mate2))
                        i1 += 1
                    elif event1.start > event2.start:
                        results.append(cls._integrate_decision_maker_single(
                            fragment, event2, fragment.mate1))
                        i2 += 1
                    else:
                        results.append(
                            cls._integrate_decision_maker_double(fragment, event1, event2))
                        i1 += 1
                        i2 += 1
        return list(filter(lambda item: item is not None, results))

    @classmethod
    def from_fragment(cls, fragment, fasta=None):
        array1 = AlignmentMismatchEventFactory.from_alignment(
            fragment.mate1, fasta)
        array2 = AlignmentMismatchEventFactory.from_alignment(
            fragment.mate2, fasta)
        results = cls._integrate_fragment_mismatch_event(
            fragment, array1, array2)
        return results