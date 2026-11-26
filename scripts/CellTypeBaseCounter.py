#!/usr/bin/env python
"""
CellTypeBaseCounter.py - single-pass, cell-type-aware replacement for

    split_by_celltype-parallel.py  ->  samtools index  ->  BaseCellCounter.py (x N cell types)
    ->  MergeBaseCellCounts.py

It reads the ORIGINAL (unsplit) BAM once, assigns every read to a cell type
through its CB tag, and writes the merged "AllCellTypes" table that
SComatic's BaseCellCalling.step1.py reads, unchanged.

Per cell type, the counting rules copy BaseCellCounter.run_interval exactly
(same pysam pileup filters, same EasyReadPileup allele coding, same DP/NC
thresholds, same INFO string), so at each site the per-cell-type columns are
what the split-BAM route would have produced. Two things differ by design:

  1. --min_total_alt N (default 1): a site is only processed and written if
     at least N reads (summed over ALL reads in the column) carry a non-reference
     allele. With N=1 only sites with zero alt reads are dropped. Such sites can
     never appear in BaseCellCalling.step1 output (step1 only writes sites with
     >=1 significant alt allele), so step1 + step2 calls are unchanged. They
     DO change step1's *.coverage_cell_count.report.tsv (callable-site
     denominators). Use --min_total_alt 0 to keep every site.
     N=3 (= step1 --min_ac_reads) is faster still, but can change the step2
     "Clustered" flag, because step1 can report sites with 1-2 alt reads that
     step2 counts as neighbours.

  2. --max_dp is applied to the whole BAM, not per cell-type BAM. The default
     here is 0 (no cap). With a cap, the reads kept at sites deeper than the cap
     differ from the split route.

Barcodes are matched to the metadata 'Index' column with the same batch
prefix/suffix logic as split_by_celltype-parallel.py. If the exact key is not
found, the key is retried with the CB's "-N" suffix removed, so metadata with
or without "-1" both work.
"""

import argparse
import multiprocessing as mp
import os
import sys
import time

import pandas as pd
import pysam

ALLELES_OUT = ['A', 'C', 'T', 'G', 'I', 'D']          # order of the count fields
COUNT_KEYS = ('A', 'C', 'T', 'G', 'D', 'I', 'N', 'O')  # what BaseCellCounter counts
BASES = frozenset('ACGTN')
INFO_FIELD = '|'.join(['DP', 'NC', 'CC', 'BC', 'BQ', 'BCf', 'BCr'])

HEADER_CONCEPTS = (
    '##INFO=DP,Description="Depth of coverage">\n'
    '##INFO=NC,Description="Number of different cells">\n'
    '##INFO=CC,Description="Cell counts [A:C:T:G:I:D:N:O], where D means deletion, I insertion and O other type of character">\n'
    '##INFO=BC,Description="Base counts [A:C:T:G:I:D:N:O], where D means deletion, I insertion and O other type of character">\n'
    '##INFO=BQ,Description="Base quality sums [A:C:T:G:I:D:N:O], where D means deletion, I insertion and O other type of character">\n'
    '##INFO=BCf,Description="Base counts in forward reads [A:C:T:G:I:D:N:O], where D means deletion, I insertion and O other type of character">\n'
    '##INFO=BCr,Description="Base counts in reverse reads [A:C:T:G:I:D:N:O], where D means deletion, I insertion and O other type of character">'
)


# ----------------------------------------------------------------------------
# Barcode -> cell type
# ----------------------------------------------------------------------------
def resolve_batch_id(bam, sra_meta, batch_column, batch_id):
    """Same rules as split_by_celltype-parallel.py."""
    bam_name = os.path.basename(bam).split('.')[0]
    if batch_id is not None:
        return bam_name, batch_id
    if sra_meta is None:
        return bam_name, bam_name.split('_')[0]
    df = pd.read_csv(sra_meta)
    col = batch_column if batch_column else df.columns[-1]
    hit = df[df[df.columns[0]] == bam_name]
    if hit.empty:
        sys.exit(f'ERROR: {bam_name} not found in first column of {sra_meta}')
    return bam_name, str(hit[col].values[0])


def load_barcode_map(meta, celltypes, batch_id, batch_pos):
    md = pd.read_csv(meta, sep='\t', dtype=str)
    ctcol = 'Cell_type' if 'Cell_type' in md.columns else md.columns[1]
    if celltypes:
        md = md[md[ctcol].isin(celltypes)]
        ct_list = list(celltypes)
    else:
        if batch_pos == 'prefix':
            md_b = md[md['Index'].str.startswith(batch_id + '_')]
        elif batch_pos == 'suffix':
            md_b = md[md['Index'].str.endswith('_' + batch_id)]
        else:
            md_b = md
        ct_list = sorted(md_b[ctcol].unique())
        md = md[md[ctcol].isin(ct_list)]
    ct_index = {ct: i for i, ct in enumerate(ct_list)}
    cb2ct = {k: ct_index[v] for k, v in zip(md['Index'], md[ctcol])}
    return cb2ct, ct_list


class BarcodeResolver:
    """CB tag -> cell-type index (or -1). Results are cached per raw CB."""

    def __init__(self, cb2ct, batch_id, batch_pos):
        self.cb2ct = cb2ct
        self.cache = {}
        if batch_pos == 'prefix':
            self.fmt = lambda cb: f'{batch_id}_{cb}'
        elif batch_pos == 'suffix':
            self.fmt = lambda cb: f'{cb}_{batch_id}'
        else:
            self.fmt = lambda cb: cb

    def __call__(self, cb):
        ct = self.cache.get(cb)
        if ct is None:
            ct = self.cb2ct.get(self.fmt(cb))
            if ct is None and '-' in cb:
                ct = self.cb2ct.get(self.fmt(cb.rsplit('-', 1)[0]))
            ct = -1 if ct is None else ct
            self.cache[cb] = ct
        return ct


# ----------------------------------------------------------------------------
# Pileup helpers (copied logic from SComatic BaseCellCounter.EasyReadPileup)
# ----------------------------------------------------------------------------
def easy_read_pileup(lst, ref):
    new = []
    for x in lst:
        up = x.upper()
        if up in BASES:
            new.append(up)
        elif len(x) > 1 and x[1] == '-':
            new.append('D')
        elif len(x) > 1 and x[1] == '+':
            new.append('I')
        elif x == '*':
            new.append('O')
        else:
            new.append('NA')
    return new


def total_alt_count(seqs, ref):
    """Alt count over all reads, as EasyReadPileup's AC (ACGTN != ref, I, D)."""
    joined = ''.join(seqs)
    if len(joined) == len(seqs):          # no indel strings: vectorised count
        u = joined.upper()
        return (u.count('A') + u.count('C') + u.count('G') + u.count('T')
                + u.count('N') - u.count(ref))
    ac = 0
    for x in seqs:
        if len(x) > 1 and (x[1] == '-' or x[1] == '+'):
            ac += 1
        else:
            up = x.upper()
            if up in BASES and up != ref:
                ac += 1
    return ac


# ----------------------------------------------------------------------------
# Worker
# ----------------------------------------------------------------------------
G = {}


def init_worker(bam, fasta, cb2ct, batch_id, batch_pos, n_ct, opts):
    G['bam'] = pysam.AlignmentFile(bam, 'rb')
    G['fa'] = pysam.FastaFile(fasta)
    G['resolve'] = BarcodeResolver(cb2ct, batch_id, batch_pos)
    G['n_ct'] = n_ct
    G['o'] = opts


def run_window(win):
    chrom, start, end = win
    bam, fa, resolve, n_ct, o = G['bam'], G['fa'], G['resolve'], G['n_ct'], G['o']
    min_cov, min_cc, min_ac_ct, min_tot = o['min_dp'], o['min_cc'], o['min_ac'], o['min_total_alt']

    ref_seq = fa.fetch(chrom, start, end).upper()
    out = []
    n_cols = n_cand = 0

    for p in bam.pileup(chrom, start, end, min_base_quality=o['min_bq'],
                        min_mapping_quality=o['min_mq'], max_depth=o['max_dp'],
                        ignore_overlaps=False):
        pos = p.reference_pos
        if pos < start or pos >= end:
            continue
        ref = ref_seq[pos - start]
        if ref == 'N':
            continue
        # every cell type needs DP >= min_dp on its own, so the column must too
        if p.get_num_aligned() < min_cov:
            continue
        n_cols += 1

        seqs = p.get_query_sequences(mark_matches=True, add_indels=True)
        if min_tot > 0 and total_alt_count(seqs, ref) < min_tot:
            continue
        n_cand += 1

        new = easy_read_pileup(seqs, ref)
        quals = p.get_query_qualities()
        reads = p.pileups

        # ---- assign each read in the column to a cell type ----
        ct_idx = [-1] * len(reads)
        aln_cache = []
        cb_cache = [None] * len(reads)
        dp_ct = [0] * n_ct
        callable_ct = [0] * n_ct
        ac_ct = [0] * n_ct
        for i, pr in enumerate(reads):
            a = pr.alignment
            aln_cache.append(a)
            try:
                cb = a.get_tag('CB')
            except KeyError:
                continue
            c = resolve(cb)
            if c < 0:
                continue
            ct_idx[i] = c
            dp_ct[c] += 1
            cb_cache[i] = cb
            b = new[i]
            if b != 'NA':
                callable_ct[c] += 1
                if b != ref and b != 'O':
                    # EasyReadPileup AC: ACGTN != ref, plus D and I (never 'O'/'NA')
                    ac_ct[c] += 1

        # ---- per cell type, apply BaseCellCounter's thresholds ----
        active = [c for c in range(n_ct)
                  if dp_ct[c] >= min_cov and callable_ct[c] >= min_cov and ac_ct[c] >= min_ac_ct]
        if not active:
            continue

        bc = {c: dict.fromkeys(COUNT_KEYS, 0) for c in active}
        bq = {c: dict.fromkeys(COUNT_KEYS, 0) for c in active}
        bcf = {c: dict.fromkeys(COUNT_KEYS, 0) for c in active}
        bcr = {c: dict.fromkeys(COUNT_KEYS, 0) for c in active}
        cells = {c: {k: set() for k in COUNT_KEYS} for c in active}
        allcells = {c: set() for c in active}
        count = {c: 0 for c in active}
        active_set = set(active)

        for i, c in enumerate(ct_idx):
            if c not in active_set:
                continue
            a = aln_cache[i]
            if a.is_secondary or a.is_duplicate or a.is_supplementary:
                continue
            b = new[i]
            if b == 'NA':
                continue
            cb = cb_cache[i].split('-')[0]
            count[c] += 1
            bc[c][b] += 1
            bq[c][b] += quals[i]
            if a.is_reverse:
                bcr[c][b] += 1
            else:
                bcf[c][b] += 1
            cells[c][b].add(cb)
            allcells[c].add(cb)

        cols = ['NA'] * n_ct
        any_pass = False
        for c in active:
            if count[c] < min_cov:
                continue
            nc = len(allcells[c])
            if nc < min_cc:
                continue
            any_pass = True
            cols[c] = '|'.join([
                str(count[c]), str(nc),
                ':'.join(str(len(cells[c][x])) for x in ALLELES_OUT),
                ':'.join(str(bc[c][x]) for x in ALLELES_OUT),
                ':'.join(str(bq[c][x]) for x in ALLELES_OUT),
                ':'.join(str(bcf[c][x]) for x in ALLELES_OUT),
                ':'.join(str(bcr[c][x]) for x in ALLELES_OUT),
            ])
        if any_pass:
            p1 = str(pos + 1)
            out.append('\t'.join([chrom, p1, p1, ref, INFO_FIELD] + cols))

    return win, out, n_cols, n_cand


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def make_windows(bam, chrom, bed, bin_size):
    if bed:
        regions = []
        with open(bed) as f:
            for line in f:
                if line.startswith(('#', 'track', 'browser')) or not line.strip():
                    continue
                c, s, e = line.split('\t')[:3]
                if chrom in (None, 'all') or c == chrom:
                    regions.append((c, int(s), int(e)))
    else:
        with pysam.AlignmentFile(bam) as b:
            chroms = b.references if chrom in (None, 'all') else [chrom]
            regions = [(c, 0, b.get_reference_length(c)) for c in chroms]
    wins = []
    for c, s, e in regions:
        for ws in range(s, e, bin_size):
            wins.append((c, ws, min(ws + bin_size, e)))
    return wins


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bam', required=True, help='Original (unsplit), coordinate-sorted, indexed BAM')
    ap.add_argument('--ref', required=True, help='Reference FASTA (indexed)')
    ap.add_argument('--chrom', default='all', help='Chromosome to analyse, or "all"')
    ap.add_argument('--bed', default=None, help='Optional BED of regions to analyse')
    ap.add_argument('--meta', required=True, help='TSV with Index (barcode) and Cell_type columns')
    ap.add_argument('--celltypes', default=None, help='Comma-separated cell types (default: all for this batch)')
    ap.add_argument('--batch_id', default=None, help='Batch ID (overrides --sra_meta)')
    ap.add_argument('--sra_meta', default=None, help='CSV mapping BAM name (col 1) to batch ID')
    ap.add_argument('--batch_column', default=None, help='Batch column in --sra_meta (default: last)')
    ap.add_argument('--batch_pos', choices=['prefix', 'suffix', 'None'], default='suffix')
    ap.add_argument('--outfile', required=True, help='Merged AllCellTypes TSV (step1 input)')
    ap.add_argument('--colname', default='{sample}_{celltype}_{chrom}',
                    help='Column name template for each cell type. Default matches the '
                         'SomaticMutationPipeline split-BAM file names.')
    ap.add_argument('--nprocs', type=int, default=1)
    ap.add_argument('--bin', type=int, default=200000, help='Window size per task [200000]')
    ap.add_argument('--min_dp', type=int, default=5)
    ap.add_argument('--min_cc', type=int, default=5)
    ap.add_argument('--min_ac', type=int, default=0, help='Per-cell-type min alt reads (BaseCellCounter --min_ac)')
    ap.add_argument('--min_bq', type=int, default=20)
    ap.add_argument('--min_mq', type=int, default=255)
    ap.add_argument('--max_dp', type=int, default=0, help='Pileup depth cap over the whole BAM; 0 = none')
    ap.add_argument('--min_total_alt', type=int, default=1,
                    help='Skip sites with fewer alt reads than this over all reads; 0 = keep all [1]')
    ap.add_argument('--keep_empty_celltypes', action='store_true',
                    help='Write all-NA columns for cell types with no rows (default: omit, as the split route does)')
    args = ap.parse_args()

    t0 = time.time()
    sample, batch_id = resolve_batch_id(args.bam, args.sra_meta, args.batch_column, args.batch_id)
    celltypes = args.celltypes.split(',') if args.celltypes else None
    cb2ct, ct_list = load_barcode_map(args.meta, celltypes, batch_id, args.batch_pos)
    print(f'Sample {sample}, batch {batch_id}, {len(cb2ct)} barcodes, '
          f'{len(ct_list)} cell types: {",".join(ct_list)}', file=sys.stderr)

    wins = make_windows(args.bam, args.chrom, args.bed, args.bin)
    # NB: pysam treats max_depth=0 as its default cap (~8000), not "unlimited".
    max_dp = args.max_dp if args.max_dp > 0 else 2_000_000_000
    opts = dict(min_dp=args.min_dp, min_cc=args.min_cc, min_ac=args.min_ac, min_bq=args.min_bq,
                min_mq=args.min_mq, max_dp=max_dp, min_total_alt=args.min_total_alt)
    init_args = (args.bam, args.ref, cb2ct, batch_id, args.batch_pos, len(ct_list), opts)

    chrom_label = args.chrom
    colnames = [args.colname.format(sample=sample, celltype=ct, chrom=chrom_label) for ct in ct_list]

    os.makedirs(os.path.dirname(os.path.abspath(args.outfile)), exist_ok=True)
    body = args.outfile + '.body.partial'
    n_rows = n_cols = n_cand = 0
    seen = [False] * len(ct_list)   # cell types with >= 1 non-NA row

    def emit(out, lines):
        for ln in lines:
            v = ln.split('\t', 5 + len(ct_list))[5:]
            for i, x in enumerate(v):
                if not seen[i] and x != 'NA':
                    seen[i] = True
        out.write('\n'.join(lines) + '\n')

    with open(body, 'w') as out:
        if args.nprocs > 1:
            with mp.Pool(args.nprocs, initializer=init_worker, initargs=init_args) as pool:
                results = pool.imap(run_window, wins, chunksize=1)  # ordered
                for _, lines, nc, nd in results:
                    if lines:
                        emit(out, lines)
                    n_rows += len(lines); n_cols += nc; n_cand += nd
        else:
            init_worker(*init_args)
            for w in wins:
                _, lines, nc, nd = run_window(w)
                if lines:
                    emit(out, lines)
                n_rows += len(lines); n_cols += nc; n_cand += nd

    # Like BaseCellCounter + MergeBaseCellCounts, leave out cell types with no rows
    # (BaseCellCounter writes no TSV for them, so the merge never sees them).
    keep = [i for i, s in enumerate(seen) if s] if not args.keep_empty_celltypes else list(range(len(ct_list)))
    dropped = [ct_list[i] for i in range(len(ct_list)) if i not in keep]
    if dropped:
        print(f'Cell types with no rows (column omitted): {",".join(dropped)}', file=sys.stderr)
    tmp_out = args.outfile + '.partial'
    with open(tmp_out, 'w') as out, open(body) as inp:
        out.write('##fileDate=%s\n' % time.strftime('%d/%m/%Y'))
        out.write(HEADER_CONCEPTS + '\n')
        out.write('\t'.join(['#CHROM', 'Start', 'End', 'REF', 'INFO'] + [colnames[i] for i in keep]) + '\n')
        if len(keep) == len(ct_list):
            for ln in inp:
                out.write(ln)
        else:
            idx = list(range(5)) + [5 + i for i in keep]
            for ln in inp:
                el = ln.rstrip('\n').split('\t')
                out.write('\t'.join(el[j] for j in idx) + '\n')
    os.remove(body)
    os.replace(tmp_out, args.outfile)

    print(f'{len(wins)} windows; {n_cols:,} columns with depth >= {args.min_dp}; '
          f'{n_cand:,} passed --min_total_alt {args.min_total_alt}; {n_rows:,} rows written; '
          f'{time.time() - t0:.1f} s', file=sys.stderr)


if __name__ == '__main__':
    main()
