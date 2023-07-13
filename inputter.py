#!/usr/bin/env python

"""
parses gtf file; assumes input gtf has exon features that have 
  attributes transcript_id, and gene_id; 
"""

## system:
import re
import sys

## local:
import util


###############################################################################
## ingest gtf file stuff:

def parse_attributes(tok9, prefix='None'):
    '''
    Args:
      gene_id "PB.38"; transcript_id "PB.38.1";
    Returns: 
      tuple (gene_id, transcript_id)
    '''

    if prefix is None:
        prefix = ''

    toks = re.split('; |;| ', tok9)
    toks = [tok.strip('"') for tok in toks]

    transcript_id = None
    gene_id = None
    idx_last = len(toks) - 1
    idx = 0

    while(idx < idx_last):
        if toks[idx] == 'transcript_id':
            if transcript_id:
                raise Exception(f"parse_attributes: multiple transcript_ids in {tok9}")
            else:
                transcript_id = f"{prefix}:{toks[idx + 1]}"
                idx += 2
        elif toks[idx] == 'gene_id':
            if gene_id:
                raise Exception(f"parse_attributes: multiple gene_ids in {tok9}")
            else:
                gene_id = f"{prefix}:{toks[idx + 1]}"
                idx += 2
        else:
            idx += 1

    if transcript_id is None:
        raise Exception(f"no transcript_id found in attributes {tok9}")

    if gene_id is None:
        raise Exception(f"no gene_id found in attributes {tok9}")

    return transcript_id, gene_id


def parse_toks(toks, label):

    try:
        transcript_id, gene_id = parse_attributes(toks[8], label)
        start = int(toks[3])
        end = int(toks[4])
    except Exception as e:
        raise Exception(f"parse_toks:1: {e}")

    if toks[6] == '+':
        strand = 1
    elif toks[6] == '-':
        strand = 0
    else:
        raise Exception(
            f"parse_toks:2: unexpected strand '{toks[6]}'"
        )

    return start, end, strand, transcript_id, gene_id


def ingest_exon(label, toks, id2idx, dat):
    '''
    Args:
      label: label for gtf file to be prepended to sequence ids
      toks: tokens from exon line in gtf file
      id2idx:
      dat: object to be updated
    Returns: None
    Side-effect: updates dat

    KI270721.1    .    exon    50359   50976  40  -   .   gene_id "PB.38"; transcript_id "PB.38.1";
    '''

    try: 
        start, end, strand, transcript_id, gene_id = parse_toks(toks, label)
    except Exception as e:
        raise Exception(f"ingest_exon: for label {label}; toks {toks}: {e}")

    ## register chromosome:
    chr_idx = dat['chr2idx'].get(toks[0])
    if chr_idx is None:
        dat['chrs'].append(toks[0])
        chr_idx = len(dat['chrs']) - 1
        dat['chr2idx'][toks[0]] = chr_idx

    ## register transcript w/ first exon:
    rna_idx = id2idx.get(transcript_id)

    if rna_idx is None:
        dat['transcripts'].append([])
        dat['old_ids'].append(transcript_id)
        dat['old_genes'].append(gene_id)
        rna_idx = len(dat['transcripts']) - 1
        id2idx[transcript_id] = rna_idx

    dat['transcripts'][rna_idx].append([chr_idx, strand, start, end])


def sort_exons(dat):
    for exons in dat['transcripts']:
        exons.sort(key=lambda exon: (exon[2], exon[3]))


def rev_neg_exons(dat):
    for idx, exons in enumerate(dat['transcripts']):
        if exons[0][1] == 0:   ## first exon on negative strand
            exons.reverse()


def ingest_gtf_file(label, gtf, dat, params):
    '''
    entry point
    Args:
      label: label for gtf file to be prepended to sequence ids
      gtf: path to gtf file
      dat: data structure to be updated
      params: run-time configuration parameters
    Returns: None
    Side-effect: updates dat
    '''

    try:
        with open(gtf, 'r') as fh:
            id2idx = {}
            for line in fh:
                toks = line.strip().split('\t')
                if not toks:
                    continue
                if toks[0].startswith('#'):
                    continue
                if len(toks) != 9:
                    raise Exception(f"ncolumns != 9 on line: {line}")
                if toks[2] != 'exon':
                    continue
                ingest_exon(label, toks, id2idx, dat)
    except Exception as e:
        raise Exception(f"ingest_gtf_file: for {label}; {gtf}: {e}")


###############################################################################
## misc entrypoints:

def ingest_gtf_list_file(gtf_list_file):
    '''
    entry point
    input: gtf_list_file: a tab-delimited text file w/ one path to gtf file per line
    output: dict of labels and gtf file paths {label0: path0, label1: path1, ...}
    '''

    id2gtf = {}
    try:
        with open(gtf_list_file, 'r') as fh:
            for line in fh:
                line = line.strip()
                if not len(line):
                    continue
                if line.startswith('#'):
                    continue
                toks = line.split('\t')
                if len(toks) != 2:
                    raise Exception(f"len(toks) != 2 ({len(toks)}) on line: {line}")
                if toks[0] in id2gtf:
                    raise Exception(f"identifier '{toks[0]}' used more than once.")
                id2gtf[toks[0]] = toks[1]
    except Exception as e:
        raise Exception(f"ingest_gtf_list_file: for gtf_list_file {gtf_list_file}: {e}")

    return id2gtf


def identify_fusions(dat, params):
    '''
    entrypoint
    updates dat['is_fusion'] based on dat['transcripts']
    '''

    transcripts = dat['transcripts']
    is_fusion = dat['is_fusion']

    for transcript in transcripts:

        fusion = False
        chrom = transcript[0][0]
        strand = transcript[0][1]
        last_end = transcript[0][2]   ## first exon start

        for exon in transcript:

            if chrom != exon[0]:
                fusion = True         ## exons on different chromosomes
            elif strand != exon[1]:
                fusion = True         ## exons on different strands
            elif exon[2] < last_end: 
                fusion = True         ## exons in wrong order
            elif (exon[2] - last_end) > params.max_intron_length:
                fusion = True         ## intron between successive exons too big
            else:
                last_end = exon[3]    ## update end of last exon

            if fusion:
                break

        is_fusion.append(fusion)


def populate_tranges(dat, params):
    '''
    populate dat['tranges'] using info from dat['transcripts']
    '''

    tranges = dat['tranges']

    for chr_i in dat['chrs']:
        tranges.append([])

    for i_transcript, transcript in enumerate(dat['transcripts']):

        if not transcript:
            continue

        ## each entry corresponds to a genomic segment of transcript:
        chrstrands = {}         ## (chr, strand): [start, end]
        for exon in transcript:
            k = (exon[0], exon[1])
            if k not in chrstrands:
                chrstrands[k] = [sys.maxsize, -1]
            if chrstrands[k][0] > exon[2]:
                chrstrands[k][0] = exon[2]
            if chrstrands[k][1] < exon[3]:
                chrstrands[k][1] = exon[3]

        for k, v in chrstrands.items():
            i_chr, i_strand = k
            start, end = v
            tranges[i_chr].append([start, end, i_strand, i_transcript])

    print(f"  {util.elapsed(params)}: sorting tranges")
    for chr_i in tranges:
        chr_i.sort(key=lambda item: tuple(item[:]))


