from . import BlockTools
## categories, subcategories
# FSM, 
#   reference match, 
#   alternative 3'end, 
#   alternative 5'end, 
#   alternative 3'end and 5'end.
# ISM, 5'fragment, 3'fragment, internal fragment, intron retention, mono-exon
# NNC, 
# NIC, 
# Genic intron, 
# Genic genomic

SQANTI3_NOT_OVERLAP = 0
SQANTI3_OVERLAP = 1
SQANTI3_FSM = 2
SQANTI3_ISM = 3
SQANTI3_NIC = 4
SQANTI3_NNC = 5

SQANTI3_FSM_REFERENCE_MATCH = 0
SQANTI3_FSM_ALTERNATIVE_3_END = 1
SQANTI3_FSM_ALTERNATIVE_5_END = 2
SQANTI3_FSM_ALTERNATIVE_3_AND_5_END = 3

def is_fsm():
    pass

def compare_sqanti3(blocks1, blocks2):
    """An implementation of SQANTI3 algorithm.

    Args:
        blocks1 (list): reference list of blocks.
        blocks2 (list): query list of blocks.

    Returns:
        tuple: (categories, subcategories)
    """
    assert len(blocks1) >= 1
    assert len(blocks2) >= 1
    splice_bias = 0
    edge_bias = 50
    category, subcategory = None, None
    if len(blocks1) == 1:
        if len(blocks2) == 1:
            x1, y1 = blocks1[0][0], blocks1[0][1]
            x2, y2 = blocks2[0][0], blocks2[0][1]
            x3, y3 = max(x1, x2), min(y1, y2)
            if x3 < y3:
                category = SQANTI3_OVERLAP
            else:
                category = SQANTI3_NOT_OVERLAP
        else:
            x1, y1 = blocks1[0][0], blocks1[0][1]
            x2, y2 = blocks2[0][0], blocks2[-1][1]
            x3, y3 = max(x1, x2), min(y1, y2)
            if x3 < y3:
                category = SQANTI3_OVERLAP
            else:
                category = SQANTI3_NOT_OVERLAP
    else:
        if len(blocks2) == 1:
            pass
        else:
            # FSM
            gaps1 = BlockTools.gaps(blocks1)
            gaps2 = BlockTools.gaps(blocks2)
            if len(gaps1) == len(gaps2):
                is_fsm = True
                for gap1, gap2 in zip(gaps1, gaps2):
                    x1, y1 = gap1
                    x2, y2 = gap2
                    if abs(x1 - x2) > splice_bias or abs(y1 - y2) > splice_bias:
                        is_fsm = False
                        break
                if is_fsm:
                    category = SQANTI3_FSM
                    x1, y1 = blocks1[0][0], blocks1[-1][1]
                    x2, y2 = blocks2[0][0], blocks2[-1][1]
                    diff5 = abs(x1 - x2)
                    diff3 = abs(y1 - y2)
                    if diff5 > edge_bias:
                        if diff3 > edge_bias:
                            subcategory = SQANTI3_FSM_ALTERNATIVE_3_AND_5_END
                        else:
                            subcategory = SQANTI3_FSM_ALTERNATIVE_5_END
                    else:
                        if diff3 > edge_bias:
                            subcategory = SQANTI3_FSM_ALTERNATIVE_3_END
                        else:
                            subcategory = SQANTI3_FSM_REFERENCE_MATCH
                else:
                    pass
            elif len(gaps1) > len(gaps2):
                pass
            else:
                pass
    
    
    return category, subcategory
    

def compare_gffcompare(ref, query):
    pass

def compare(ref, query, engine="sqanti3"):
    if engine == "sqanti3":
        return compare_sqanti3(ref, query)
    elif engine == "gffcompare":
        return compare_gffcompare(ref, query)
    else:
        raise RuntimeError("Unsupported engine %s." % engine)