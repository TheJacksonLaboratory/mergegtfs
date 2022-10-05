#!/usr/bin/env python3

## system:
import sys
import time

## local:
import initializer
import inputter
import namer
import outputter
import overlapper
import resolver
import util

'''

dat = {
  ## label is from gtf_list_file; gtf_chr_id is from col 1 of gtf file;
  ## chr_id: f'{label}:{gtf_chr_id}'
  chrs: [chr_id0, chr_id1, chr_id2, ... ]
  chr2idx: {chr_id0: 0, chr_id1: 1, chr_id2: 2, ...}

  ## strand == '-': 0; strand == '+': 1;
  ## exon = [seq, strand, start, end]
  ## transcript = [exon1, exon2, ...]
  transcripts: [transcript0, transcript1, transcript2, ...]
  ## in register w/ transcripts; True only if exons on different 
  ##   chromosomes or if: 
  ##     abs(exon[i+1].start - exon[i].end) > params.max_intron_size
  is_fusion: [True, False, True, ...]
  new_ids: [new_transcript_id0, new_transcript_id1, ...]
  new_genes: [new_gene0, new_gene1, new_gene2, ...]

  ## transcript_id: f'{label}:{gtf_transcript_id}';
  ## ids in register w/ transcripts:
  old_ids: [transcript_id0, transcript_id1, transcript_id2, ... ]

  ## gene_ids: f'{label}:{gtf_gene_id}'
  old_genes: [gene_id0, gene_id1, gene_id2, ... ]

  ## tranges[i] corresponds to chrs[i]; transcript is transcripts[idx]:
  tranges: [[[start, end, strand, idx], ...], ...] 

  ## initially, *_transcript = transcripts[*_idx]:
  xrefs = { collapsed_idx: kept_idx, ... }  
  ## later, kept is True or False; others are strings:
  xrefs = [[old_id1, old_gene1, new_id1, new_gene1, kept1], ...]

  ## transcript1 = transcripts[idx1]:
  olaps = { idx1: {idx2, idx3, ...}, idx2: {idx1, idx3, ...}, ... }
}
'''

###############################################################################
## main:

if __name__ == '__main__':

    params = initializer.initialize()

    print(f"{util.elapsed(params)}: parsing gtf_list_file:")
    try:
        id2gtf = inputter.ingest_gtf_list_file(params.gtf_list_file)
    except Exception as e:
        sys.stderr.write(f"ERROR:31: {e}\n")
        sys.exit(31)

    dat = {
        ## seq_name: f'{label}:{seqid}'
        'chrs': [],
        'chr2idx': {},
        ## exon = [seq, strand, start, end]
        ## transcript = [exon1, exon2, ...]
        'transcripts': [],
        'is_fusion': [],
        'old_ids': [],
        'old_genes': [],
        'new_ids': [],
        'new_genes': [],
        'tranges': [],
        'xrefs': {},
        'olaps': {}
    }

    try:
        for label, gtf in id2gtf.items():
            print(f"{util.elapsed(params)}: ingesting {label}: {gtf}")
            inputter.ingest_gtf_file(label, gtf, dat, params)
    except Exception as e:
        sys.stderr.write(f"ERROR:33: {e}\n")
        sys.exit(33)

    if params.sort_exons:
        print(f"{util.elapsed(params)}: sorting exons")
        inputter.sort_exons(dat)
    elif params.rev_neg_exons:
        print(f"{util.elapsed(params)}: reversing negative strand exon order")
        inputter.rev_neg_exons(dat)

    print(f"{util.elapsed(params)}: identifying fusions")
    inputter.identify_fusions(dat, params)

    print(f"{util.elapsed(params)}: populating tranges")
    inputter.populate_tranges(dat, params)

    print(f"{util.elapsed(params)}: cross-referencing transcripts")
    resolver.xref_transcripts(dat, params)

    print(f"{util.elapsed(params)}: resolving xrefs")
    resolver.resolve_xrefs(dat)

    print(f"{util.elapsed(params)}: repopulating tranges")
    dat['tranges'] = []
    inputter.populate_tranges(dat, params)

    print(f"{util.elapsed(params)}: finding transcript overlaps")
    overlapper.find_overlaps(dat, params,)

    print(f"{util.elapsed(params)}: naming genes")
    namer.name_genes(dat, params)

    print(f"{util.elapsed(params)}: naming transcripts")
    namer.name_transcripts(dat)

    print(f"{util.elapsed(params)}: writing xref file")
    try:
        outputter.write_xref_file(dat, params)
    except Exception as e:
        sys.stderr.write(f"ERROR:41: {e}\n")
        sys.exit(41)

    print(f"{util.elapsed(params)}: writing gtf file")
    try:
        n_transcripts, n_exons = outputter.write_gtf_file(dat, params)
    except Exception as e:
        sys.stderr.write(f"ERROR:53: {e}\n")
        sys.exit(53)

    outputter.report(dat)

    print(f"{util.elapsed(params)}: completed")
    print(f"finished: {util.time_stamp()}")
    sys.exit(0)

