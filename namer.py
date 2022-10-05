import sys

## local:

import util

def name_transcripts(dat):
    '''
    populates dat['new_ids'] using dat['olaps']
    '''

    olaps = dat['olaps']
    tranges_chr = dat['tranges']
    new_genes = dat['new_genes']
    dat['new_ids'] = [None] * len(dat['transcripts'])
    new_ids = dat['new_ids']
    done = set()
    counts = {}

    for trange_chr in tranges_chr:     ## in original chromosome order
        for trange in trange_chr:      ## in (start, end) order
            idx = trange[3]            ## transcript = dat['transcripts'][idx]
            if idx in done:
                continue
            gene = new_genes[idx]
            if gene not in counts:
                counts[gene] = 0
            counts[gene] += 1
            new_ids[idx] = f"{gene}.{counts[gene]}"
            done.add(idx)


def name_genes_i(dat, params, name_fusions=False, n_genes=0):
    '''
    populates dat['new_genes'] using dat['olaps']
    '''

    olaps = dat['olaps']
    is_fusion = dat['is_fusion']
    tranges_chr = dat['tranges']
    new_genes = dat['new_genes']
    prefix = params.gene_prefix

    for trange_chr in tranges_chr:     ## in original chromosome order
        for trange in trange_chr:      ## in (start, end) order
            idx = trange[3]            ## transcript = dat['transcripts'][idx]
            if new_genes[idx] is not None:
                continue
            if is_fusion[idx] != name_fusions:
                continue
            links = olaps.get(idx)
            if idx in links:
                n_genes += 1
                new_genes[idx] = f"{prefix}{n_genes}"
            elif len(links) == 1:
                link_idx = list(links)[0]
                link_gene = new_genes[link_idx]
                if link_gene is not None:
                    new_genes[idx] = link_gene
            elif is_fusion[idx]:
                nom = 'fusion'
                for i in links:
                    nom += f':{new_genes[i]}'
                new_genes[idx] = nom
            else:
                ## this is a strange situation:
                nom = None
                for i in list(olaps[idx]):
                    if new_genes[i] is not None:
                        nom = new_genes[i]
                        break
                if nom is not None:
                    new_genes[idx] = nom
                else:
                    n_genes += 1
                    new_genes[idx] = f"{prefix}{n_genes}"

    return n_genes


def none_gene_found(dat, check_fusions=False):

    new_genes = dat['new_genes']
    is_fusion = dat['is_fusion']

    for tranges_chr in dat['tranges']:
        for trange in tranges_chr:
            transcript_idx = trange[3]
            if is_fusion[transcript_idx] != check_fusions:
                continue
            if new_genes[transcript_idx] is None:
                return True

    return False


def name_genes(dat, params):
    '''
    populates dat['new_genes'] using dat['olaps']
    '''

    dat['new_genes'] = [None] * len(dat['transcripts'])

    n_genes = 0

    while none_gene_found(dat, False):
        print(f"{util.elapsed(params)}: naming non-fusions;")
        n_genes = name_genes_i(dat, params, False, n_genes)

    while none_gene_found(dat, True):
        print(f"{util.elapsed(params)}: naming fusions;")
        n_genes = name_genes_i(dat, params, True, n_genes)


def xref_table(dat, params):
    '''
    repopulates dat['xref'] 
    '''

    xref_table = []
    old_ids = dat['old_ids']
    old_genes = dat['old_genes']
    new_ids = dat['new_ids']
    new_genes = dat['new_genes']
    xrefs = dat['xrefs']

    for idx, vals in enumerate(zip(old_ids, old_genes, new_ids, new_genes)):
        old_id, old_gene, new_id, new_gene = vals
        if new_id is None:
            kept = False
            idx_kept = xrefs[idx]
            new_id = new_ids[idx_kept]
            new_gene = new_genes[idx_kept]
        else:
            kept = True

        xref_table.append([old_id, old_gene, new_id, new_gene, kept])

    dat['xrefs'] = xref_table


