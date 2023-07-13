## system:
import bisect

## local:
import util


def derive_starts(tranges, transcripts):
    '''
    starts and indices hierarchies in register w/ tranges
      starts are starts of first exon of each transcript
      indices are corresponding transcript_idxs
    '''

    starts = []
    indices = []

    for chr_trange in tranges:
        chr_starts = []
        chr_indices = []
        for trange in chr_trange:
            ## start of first exon for transcript_idx trange[3]:
            chr_starts.append(transcripts[trange[3]][0][2])
            chr_indices.append(trange[3])
        starts.append(chr_starts)
        indices.append(chr_indices)

    return starts, indices


def match_starts(transcript, starts, indices, params):
    '''
    returns list of indices of transcripts with starts
      matching start of transcript
    '''

    idx_chr, strand, start, end = transcript[0]    ## first exon

    ## '-': tts; '+': tss:
    tol = (params.tol_tts, params.tol_tss)[strand]

    ## where in tranges starts might line up well enough:
    starts_chr = starts[idx_chr]
    idx = bisect.bisect_left(starts_chr, start - tol)
    idx_end = bisect.bisect_right(starts_chr, start + tol)
    if idx_end >= len(starts_chr):
        idx_end = len(starts_chr) - 1

    ## transcript_maybe = dat.transcripts[maybe[i]]:
    maybe = []
    while idx <= idx_end:
        maybe.append(indices[idx_chr][idx])
        idx += 1

    return maybe


## xref_transcripts() innermost loop:

def transcripts_match(transcript1, transcript2, tol_tss, tol_sj, tol_tts):
    '''
    returns True/False depending on whether transcript1 and transcript2
      match to w/i stringency set by tol_*.
    '''

    if len(transcript1) != len(transcript2):
        return False             ## differing number of exons

    idx_last = len(transcript1) - 1

    for idx, (exon1, exon2) in enumerate(zip(transcript1, transcript2)):
    
        if exon1[0] != exon2[0]:
            return False         ## chromosome mismatch
        
        if exon1[1] != exon2[1]:
            return False         ## strand mismatch

        if idx_last == 0:        ## single exon
            if exon1[1]:         ## '+' strand
                tol5 = tol_tss
                tol3 = tol_tts
            else:                ## '-' strand
                tol5 = tol_tts
                tol3 = tol_tss
        elif idx == 0:           ## first of multiple exons
            tol3 = tol_sj
            if exon1[1]:
                tol5 = tol_tss
            else:
                tol5 = tol_tts
        elif idx == idx_last:    ## last of multiple exons
            tol5 = tol_sj
            if exon1[1]:
                tol3 = tol_tts
            else:
                tol3 = tol_tss
        else:
            tol5 = tol3 = tol_sj ## internal exon

        if abs(exon1[2] - exon2[2]) > tol5:
            return False

        if abs(exon1[3] - exon2[3]) > tol3:
            return False

    return True


def xref_transcripts(dat, params):
    '''
    populates dat['xrefs'] 
      based on dat['tranges'] and dat['transcripts'];
      params.tol_* determines stringency
    '''

    ## starts are starts of first exon of each transcript;
    ## indices are corresponding transcript_idxs;
    ## everything in tranges hierarchy and order:
    print(f"  {util.elapsed(params)}: deriving starts")
    starts, indices = derive_starts(dat['tranges'], dat['transcripts'])

    print(f"  {util.elapsed(params)}: matching transcripts")
    for idx1, transcript1 in enumerate(dat['transcripts']):

        if idx1 in dat['xrefs']:
            continue                     ## already merged

        maybe_list = match_starts(transcript1, starts, indices, params)

        for idx2 in maybe_list:

            if idx1 == idx2:             ## self
                continue

            if idx2 in dat['xrefs']:     ## already merged
                continue

            if transcripts_match(
              transcript1, 
              dat['transcripts'][idx2], 
              params.tol_tss,
              params.tol_sj,
              params.tol_tts
            ):
                ## extend (retained) transcript1 to reach at least end of transcript2,
                ##   transcripts_match() already ensures exon numbers match;
                ##   sort ensures transcript1 start <= transcript2 start;
                ##   this ensures transcript1 end >= transcript2 end as well:
                idx_last = len(dat['transcripts'][idx1]) - 1
                end1 = dat['transcripts'][idx1][idx_last][3]   ## 3'end of transcript1 
                end2 = dat['transcripts'][idx2][idx_last][3]   ## 3'end of transcript2
                if end1 < end2:
                    dat['transcripts'][idx1][idx_last][3] = end2

                ## cross-referencing xrefs[merged_idx] -> kept_idx:
                dat['xrefs'][idx2] = idx1
                dat['transcripts'][idx2] = None
                idx_last = len(dat['transcripts'][idx1]) - 1


def resolve_xrefs(dat):
    '''
    ensures each val in dat['xrefs'] corresponds to a retained transcript
    '''

    for id2, id1 in dat['xrefs'].items():
        while id1 in dat['xrefs']:
            id1 = dat['xrefs'][id1]
        dat['xrefs'][id2] = id1



