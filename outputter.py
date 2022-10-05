def write_xref_records(dat, fh):

    old_ids = dat['old_ids']
    old_genes = dat['old_genes']
    new_ids = dat['new_ids']
    new_genes = dat['new_genes']
    xrefs = dat['xrefs']
    is_fusion = dat['is_fusion']

    line = '\t'.join([
        'old_transcript', 
        'old_gene', 
        'new_transcript', 
        'new_gene', 
        'rearranged', 
        'exemplar'
    ])
    fh.write(f"{line}\n")

    for idx, vals in enumerate(zip(old_ids, old_genes, new_ids, new_genes, is_fusion)):

        old_id, old_gene, new_id, new_gene, is_fusion = vals

        if new_id is None:
            kept = 'False'
            idx_kept = xrefs[idx]
            new_id = new_ids[idx_kept]
            new_gene = new_genes[idx_kept]
        else:
            kept = 'True'

        line = '\t'.join([old_id, old_gene, new_id, new_gene, str(is_fusion), kept])
        fh.write(f"{line}\n")


def write_xref_file(dat, params):

    try:
        file_out = f"{params.output_prefix}.xrefs.tsv"
        with open(file_out, 'w') as fh:
            write_xref_records(dat, fh)
    except Exception as e:
        raise Exception(f"write_xref_file: for file_out {file_out}: {e}")


def write_exon_recs(transcript_idx, rec, dat, fh):

    rec[2] = 'exon'
    chrs = dat['chrs']
    transcript = dat['transcripts'][transcript_idx]

    for exon in transcript:

        rec[0] = chrs[exon[0]]        ## chromosome
        rec[3] = str(exon[2])         ## start
        rec[4] = str(exon[3])         ## end
        rec[6] = ['-', '+'][exon[1]]  ## strand

        line = '\t'.join(rec)
        fh.write(f"{line}\n")

    ## number of exons written:
    return len(transcript)


def make_transcript_rec(transcript_idx, dat):

    transcript = dat['transcripts'][transcript_idx]
    new_gene_id = dat['new_genes'][transcript_idx]
    new_transcript_id = dat['new_ids'][transcript_idx]

    attributes = (
        f'gene_id "{new_gene_id}"; '
        f'transcript_id "{new_transcript_id}";'
    )

    exon1 = transcript[0]       ## 1st exon
    idx_chr1 = exon1[0]
    end = exon1[3]

    for exon in transcript:
        if exon[0] != idx_chr1:
            continue
        if exon[3] > end:
            end = exon[3]

    rec = [
        dat['chrs'][idx_chr1],  ## chromosome for transcript from 1st exon
        'gtfmerge',             ## source (constant)
        'transcript',           ## feature
        str(exon1[2]),          ## start
        str(end),               ## end
        '.',                    ## score (constant)
        ['-', '+'][exon1[1]],   ## strand
        '.',                    ## frame (constant)
        attributes              ## attributes (constant)
    ]

    return rec


def write_transcript_records(transcript_idx, dat, fh):

    rec = make_transcript_rec(transcript_idx, dat)
    line = '\t'.join(rec)
    fh.write(f"{line}\n")

    ## return number of exon records:
    return write_exon_recs(transcript_idx, rec, dat, fh)


def write_gtf_records(dat, fh):

    tranges_chr = dat['tranges']
    transcripts = dat['transcripts']
    done = set()

    if tranges_chr is None:
        return 0, 0

    n_exons = 0

    for trange_chr in tranges_chr:     ## in original chromosome order
        for trange in trange_chr:      ## in (start, end) order
            transcript_idx = trange[3]
            if transcript_idx in done:
                continue
            n_exons += write_transcript_records(transcript_idx, dat, fh) 
            done.add(transcript_idx)

    ## n_transcripts = len(done):
    return len(done), n_exons


def write_gtf_file(dat, params):

    try:
        file_out = f"{params.output_prefix}.gtf"
        with open(file_out, 'w') as fh:
            n_transcripts, n_exons = write_gtf_records(dat, fh)
    except Exception as e:
        raise Exception(f"write_gtf_file: for file_out {file_out}: {e}")

    return n_transcripts, n_exons


def report_inputs(dat):

    ids = dat['old_ids']
    genes = dat['old_genes']
    is_fusion = dat['is_fusion']
    
    n_transcripts = 0
    n_fusions = 0
    done_genes = set()
    done_ids = set()
 
    for idx, id1 in enumerate(ids):
        n_transcripts += 1
        done_ids.add(ids[idx])
        done_genes.add(genes[idx])
        if is_fusion[idx]:
            n_fusions += 1

    print(
        f"number of input transcripts: {n_transcripts}\n"
        f"number of input fusion transcripts: {n_fusions}\n"
        f"number of unique input transcript ids: {len(done_ids)}\n"
        f"number of unique input gene ids: {len(done_genes)}"
    )


def report_outputs(dat):

    ids = dat['new_ids']
    genes = dat['new_genes']
    is_fusion = dat['is_fusion']

    n_transcripts = 0
    n_fusions = 0
    done_genes = set()
    done_ids = set()

    for idx, id1 in enumerate(ids):
        if id1 is None:
            continue
        n_transcripts += 1
        done_ids.add(ids[idx])
        done_genes.add(genes[idx])
        if is_fusion[idx]:
            n_fusions += 1

    print(
        f"number of output transcripts: {n_transcripts}\n"
        f"number of output fusion transcripts: {n_fusions}\n"
        f"number of unique output transcript ids: {len(done_ids)}\n"
        f"number of unique output gene ids: {len(done_genes)}"
    )


def report(dat):
    print(f"number of chromosomes: {len(dat['chrs'])}")
    report_inputs(dat)
    report_outputs(dat)

