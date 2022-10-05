import bisect
import math
import sys
import util


def max_loc_length(chr_tranges):
    '''
    takes chr_tranges which is one element of tranges (tranges[i])
    returns maximum segment/locus length in chr_tranges
    '''

    max_len = 0

    for trange in chr_tranges:
        len_i = trange[1] - trange[0]
        if len_i > max_len:
            max_len = len_i

    return max_len


def derive_starts(chr_tranges):
    
    starts = []     ## start of trange
    stops = []      ## stop of trange
    indices = []    ## transcript index

    for trange in chr_tranges:
        starts.append(trange[0])
        stops.append(trange[1])
        indices.append(trange[3])

    return starts, stops, indices


def exons_overlap(exon1, exon2, p_exon_overlap):

    if exon1[0] != exon2[0]:
        return False        ## chromosome mismatch

    if exon1[1] != exon2[1]:
        return False        ## strand mismatch

    ## ensure start1 <= start2, just so following checks work:
    if exon1[2] > exon2[2]:
        exon1, exon2 = exon2, exon1

    if exon1[3] < exon2[2]:
        return False        ## exon1 ends before exon2 begins

    if exon1[2] > exon2[3]:
        return False        ## exon1 begins after exon2 ends

    len1 = 1 + exon1[3] - exon1[2]
    len2 = 1 + exon2[3] - exon2[2]

    if exon1[3] < exon2[3]: ## partial overlap
        olap = exon1[3] - exon2[2]
    else:
        olap = len2         ## exon1 contains exon2

    if (olap / min(len1, len2)) >= p_exon_overlap:
        return True
    else:
        return False


def same_gene(exons1, exons2, p_exon_overlap, p_exons_overlap):

    n_exons_min = min(len(exons1), len(exons2))
    n_exons_overlap_min = math.ceil(p_exons_overlap * n_exons_min)
    n_exons_overlap = 0

    for exon1 in exons1:
        for exon2 in exons2:
            if exons_overlap(exon1, exon2, p_exon_overlap):
                n_exons_overlap += 1
                if n_exons_overlap >= n_exons_overlap_min:
                    return True

    return False


def find_overlaps_nonfusion(dat, params):
    '''
    Populates dat['olaps'] w/ overlapping segments found in dat['tranges']
    '''

    olaps = dat['olaps']
    is_fusion = dat['is_fusion']
    transcripts = dat['transcripts']
    old_ids = dat['old_ids']
    p_exon_overlap = params.p_exon_overlap
    p_exons_overlap = params.p_exons_overlap

    for idx_chr, chr_tranges in enumerate(dat['tranges']):

        print(f"{util.elapsed(params)}: processing {dat['chrs'][idx_chr]}")

        ## indices into dat['transcripts']:
        starts, stops, indices = derive_starts(chr_tranges)
        max_len = max_loc_length(chr_tranges)

        for idx_trange, trange in enumerate(chr_tranges):

            i1 = trange[3]                 ## index into dat['transcripts']
            if is_fusion[i1]:
                continue

            if i1 not in olaps:
                olaps[i1] = {i1}

            start1 = trange[0]
            stop1 = trange[1]

            idx = bisect.bisect_left(starts, start1)
            idx_end = bisect.bisect_right(starts, stop1 + max_len)

            if idx_end >= len(starts):
                idx_end = len(starts) - 1

            while idx <= idx_end:

                i2 = indices[idx]   ## index into dat['transcripts']

                if i1 == i2:
                    pass            ## self
                elif stop1 < starts[idx]:
                    break           ## range1 ends before range2 begins
                elif start1 > stops[idx]:
                    pass            ## range1 begins after range2 ends
                elif idx_trange > idx and not is_fusion[i2]:
                    pass            ## only link later tranges to earlier, unless fusion
                elif not same_gene(
                    transcripts[i1],
                    transcripts[i2],
                    p_exon_overlap,
                    p_exons_overlap
                ):
                    pass
                else:
                    if i2 not in olaps:
                        olaps[i2] = set()
                    olaps[i2] |= olaps[i1]

                idx += 1


def find_overlaps_fusion(dat, params):
    '''
    Populates dat['olaps'] w/ overlapping segments found in dat['tranges']
    '''

    olaps = dat['olaps']
    is_fusion = dat['is_fusion']
    transcripts = dat['transcripts']
    p_exon_overlap = params.p_exon_overlap
    p_exons_overlap = params.p_exons_overlap

    for idx_chr, chr_tranges in enumerate(dat['tranges']):

        print(f"{util.elapsed(params)}: processing {dat['chrs'][idx_chr]}")

        ## indices into dat['transcripts']:
        starts, stops, indices = derive_starts(chr_tranges)
        max_len = max_loc_length(chr_tranges)

        for idx_trange, trange in enumerate(chr_tranges):

            i1 = trange[3]                 ## index into dat['transcripts']
            if not is_fusion[i1]:
                continue

            start1 = trange[0]
            stop1 = trange[1]
            left_end = max(start1 - max_len, 0)
            idx = bisect.bisect_left(starts, left_end)
            idx_end = bisect.bisect_right(starts, stop1)
            if idx_end >= len(starts):
                idx_end = len(starts) - 1

            while idx <= idx_end:

                i2 = indices[idx]     ## index into dat['transcripts']
                if is_fusion[i2]:     ## covers if i1 == i2
                    pass
                elif stop1 < starts[idx]:
                    break             ## range1 ends before range2 begins
                elif start1 > stops[idx]:
                    pass              ## range1 begins after range2 ends
                elif not same_gene(
                    transcripts[i1],
                    transcripts[i2],
                    p_exon_overlap,
                    p_exons_overlap
                ):
                    pass
                else:
                    if i1 not in olaps:
                        olaps[i1] = set()
                    olaps[i1] |= olaps[i2]

                idx += 1

            if i1 not in olaps:
                olaps[i1] = {i1}


def find_overlaps(dat, params):
    '''
    entry point
    Populates dat['olaps'] w/ overlapping segments found in dat['tranges']
    '''

    print(f"{util.elapsed(params)}: finding overlaps between non-fusions")
    find_overlaps_nonfusion(dat, params)

    print(f"{util.elapsed(params)}: finding overlaps between fusions")
    find_overlaps_fusion(dat, params)


