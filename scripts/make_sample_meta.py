#!/usr/bin/env python
"""Write a per-sample barcode -> cell type table whose barcodes match the CB tags
in the ORIGINAL (unsplit) BAM, for SComatic's SingleCellGenotype.py.

SingleCellGenotype.py strips everything from the first '-' of both the metadata
Index and the BAM CB tag before matching. With batch-tagged metadata
(e.g. "B4_c2_AAACCTGAGATCGGGT-1") that only matched the split BAMs, whose CB tags
had been rewritten by split_by_celltype. This script keeps the rows for one batch
and removes the batch prefix/suffix, so the original BAM can be genotyped directly.

usage: make_sample_meta.py --meta ALL.tsv --bam SAMPLE.bam --sra_meta sra.csv \
           --batch_column batchID --batch_pos prefix [--celltypes A,B] --out SAMPLE.meta.tsv
"""
import argparse
import os

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--meta', required=True)
    ap.add_argument('--bam', required=True)
    ap.add_argument('--sra_meta', default=None)
    ap.add_argument('--batch_column', default=None)
    ap.add_argument('--batch_id', default=None)
    ap.add_argument('--batch_pos', choices=['prefix', 'suffix', 'None'], default='suffix')
    ap.add_argument('--celltypes', default=None)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    bam_name = os.path.basename(a.bam).split('.')[0]
    batch = a.batch_id
    if batch is None and a.sra_meta:
        df = pd.read_csv(a.sra_meta)
        col = a.batch_column or df.columns[-1]
        batch = str(df.loc[df[df.columns[0]] == bam_name, col].values[0])
    elif batch is None:
        batch = bam_name.split('_')[0]

    md = pd.read_csv(a.meta, sep='\t', dtype=str)
    ctcol = 'Cell_type' if 'Cell_type' in md.columns else md.columns[1]
    if a.batch_pos == 'prefix':
        md = md[md['Index'].str.startswith(batch + '_')].copy()
        md['Index'] = md['Index'].str[len(batch) + 1:]
    elif a.batch_pos == 'suffix':
        md = md[md['Index'].str.endswith('_' + batch)].copy()
        md['Index'] = md['Index'].str[:-(len(batch) + 1)]
    if a.celltypes:
        md = md[md[ctcol].isin(a.celltypes.split(','))]
    md[['Index', ctcol]].rename(columns={ctcol: 'Cell_type'}).to_csv(a.out, sep='\t', index=False)
    print(f'{bam_name}: batch {batch}, {len(md)} barcodes written to {a.out}')


if __name__ == '__main__':
    main()
