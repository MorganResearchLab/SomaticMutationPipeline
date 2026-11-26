#!/usr/bin/env python
"""Compare two SComatic merged AllCellTypes tables (or step1/step2 outputs).

usage:
  compare_merged.py merged CLASSIC.tsv NEW.tsv   # merged count tables
  compare_merged.py calls  CLASSIC.tsv NEW.tsv   # step1 / step2 call tables
"""
import sys

ALLELES = ['A', 'C', 'T', 'G', 'I', 'D']


def load(path, first_ct_col):
    rows, cols = {}, None
    with open(path) as f:
        for line in f:
            if line.startswith('##'):
                continue
            el = line.rstrip('\n').split('\t')
            if line.startswith('#'):
                cols = el
                continue
            rows[(el[0], int(el[1]))] = el
    return cols, rows


def as_dict(cols, el, first):
    """Cell-type columns keyed by name (column order may differ between files)."""
    d = {cols[i]: el[i] for i in range(first, len(cols))}
    return el[:first], d


def total_alt(cols, el, first):
    ref = el[3]
    tot = 0
    for v in el[first:]:
        if v.startswith('NA'):
            continue
        bc = v.split('|')[3].split(':')
        tot += sum(int(n) for a, n in zip(ALLELES, bc) if a != ref)
    return tot


def main():
    mode, a, b = sys.argv[1:4]
    first = 5
    ca, ra = load(a, first)
    cb, rb = load(b, first)
    if mode == 'calls':
        # step1/step2: col 4 is the calling INFO block, cell-type columns start at 6
        first = 6
    print(f'classic rows: {len(ra):,}   new rows: {len(rb):,}')
    if set(ca[first:]) != set(cb[first:]):
        print('WARNING: cell-type columns differ:', set(ca[first:]) ^ set(cb[first:]))

    only_new = [k for k in rb if k not in ra]
    only_classic = [k for k in ra if k not in rb]
    both = [k for k in ra if k in rb]
    diff = []
    for k in both:
        fa, da = as_dict(ca, ra[k], first)
        fb, db = as_dict(cb, rb[k], first)
        if fa != fb or da != db:
            diff.append(k)
    print(f'shared rows: {len(both):,}   identical: {len(both) - len(diff):,}   differing: {len(diff):,}')
    print(f'rows only in new: {len(only_new):,}')
    print(f'rows only in classic: {len(only_classic):,}')
    if mode == 'merged' and only_classic:
        unexplained = [k for k in only_classic if total_alt(ca, ra[k], first) > 0]
        print(f'  of which have alt reads in some cell type (NOT explained by --min_total_alt): {len(unexplained):,}')
        for k in unexplained[:5]:
            print('   ', '\t'.join(ra[k])[:300])
    for k in (diff[:5]):
        print('DIFF', k)
        fa, da = as_dict(ca, ra[k], first)
        fb, db = as_dict(cb, rb[k], first)
        if fa != fb:
            print('   fixed cols  classic:', fa, '\n               new:    ', fb)
        for c in da:
            if da[c] != db.get(c):
                print(f'   {c}: classic {da[c]}  new {db.get(c)}')
    for k in only_new[:5]:
        print('ONLY_NEW', '\t'.join(rb[k])[:300])
    if mode == 'calls':
        for label, rows in (('classic', ra), ('new', rb)):
            p = sum(1 for el in rows.values() if el[5] == 'PASS')
            print(f'{label} PASS calls: {p}')


if __name__ == '__main__':
    main()
