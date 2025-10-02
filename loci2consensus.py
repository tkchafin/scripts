#!/usr/bin/env python3
import sys
import argparse

IUPAC_TO_SET = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'R': {'A','G'}, 'Y': {'C','T'}, 'S': {'G','C'}, 'W': {'A','T'},
    'K': {'G','T'}, 'M': {'A','C'},
    'B': {'C','G','T'}, 'D': {'A','G','T'}, 'H': {'A','C','T'}, 'V': {'A','C','G'},
    'N': set(), '-': set(), '.': set(), '?': set()
}
SET_TO_IUPAC = {frozenset(v): k for k, v in IUPAC_TO_SET.items() if v}
SET_TO_IUPAC[frozenset({'A','C','G','T'})] = 'N'

def parse_blocks(path):
    """Yield lists of sequences per locus from ipyrad-style .loci."""
    seqs = []
    fin = sys.stdin if path == '-' else open(path, 'r')
    try:
        for raw in fin:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('//'):
                if seqs:
                    yield seqs
                    seqs = []
                continue
            if line.startswith('>'):
                continue
            parts = line.split(None, 1)
            if len(parts) < 2:
                continue
            seq = parts[1].strip().replace(' ', '').rstrip('>')
            if seq:
                seqs.append(seq.upper())
        if seqs:
            yield seqs
    finally:
        if fin is not sys.stdin:
            fin.close()

def consensus(seqs):
    """Consensus ignoring N/gaps; if all missing, N."""
    if not seqs:
        return ''
    L = max(len(s) for s in seqs)
    out = []
    for i in range(L):
        bases = set()
        any_nonmissing = False
        for s in seqs:
            ch = s[i] if i < len(s) else 'N'
            aset = IUPAC_TO_SET.get(ch, set())
            if aset:
                any_nonmissing = True
                bases |= aset
        if not any_nonmissing:
            out.append('N')
        else:
            out.append(SET_TO_IUPAC.get(frozenset(bases), 'N'))
    return ''.join(out)

def main():
    ap = argparse.ArgumentParser(description="ipyrad .loci → FASTA with per-locus consensus (IUPAC).")
    ap.add_argument('input', help="Input .loci file (or '-' for stdin)")
    ap.add_argument('output', help="Output FASTA path (or '-' for stdout)")
    args = ap.parse_args()

    out = sys.stdout if args.output == '-' else open(args.output, 'w')
    try:
        for idx, seqs in enumerate(parse_blocks(args.input), start=1):
            cons = consensus(seqs)
            if cons:
                out.write(f">LOCUS_{idx}\n{cons}\n")
    finally:
        if out is not sys.stdout:
            out.close()

if __name__ == "__main__":
    main()

