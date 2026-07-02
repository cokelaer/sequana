#!/usr/bin/env python3
"""Detect G-quadruplex forming motifs.

A G-quadruplex (G4) forms where four or more runs of guanines (*G-islands*,
each at least ``min_rep`` bp) sit close together (separated by loops of at most
``max_spacer`` bp). This reproduces the *G_Quadruplex_Motif* output of the
non-B gfa / nBMST tool (``gfa_GQ.gff``), ported faithfully from its
``findGQ.c`` / ``getGislands``.

Both strands are scanned: runs of ``G`` give plus-strand motifs, runs of ``C``
give minus-strand motifs (the reported sequence is then the reverse complement,
i.e. G-rich). For each motif gfa reports:

- ``islands`` : number of merged G-islands (``conIls``),
- ``runs``    : how many minimal (``min_rep``) runs fit across the islands
  (``npos``, must be >= 4),
- ``max``     : the largest run size for which at least 4 runs still fit.

gfa defaults (from ``gfa.c``): min_rep = 3, max_spacer = 7.

Validated against gfa on a 30+ Mb *Leishmania* assembly: identical coordinates,
strand, islands/runs/max attributes. As for Z-DNA, the only field that differs
is ``composition`` because gfa's printer counts bases over ``[start, end-1)``
(an off-by-one dropping the last base); this module counts the full motif, so
its A/C/G/T counts are correct and match the reported ``sequence``.
"""
import numpy as np
import pandas as pd
from tqdm import tqdm

from sequana import FastA

_G, _C = 71, 67  # uppercase byte codes
_COMP = bytes.maketrans(b"ACGTN", b"TGCAN")


def _islands_py(seq, target, min_rep):
    """Return (start_1based, length) for every run of ``target`` >= min_rep."""
    strt, lng = [], []
    n = len(seq)
    run = 0
    for i in range(n + 1):
        c = seq[i] if i < n else None
        if c == target:
            run += 1
        else:
            if run >= min_rep:
                strt.append(i - run + 1)  # 1-based, as in gfa
                lng.append(run)
            run = 0
    return strt, lng


def _findgq_py(strt, lng, min_rep, max_spacer):
    """Faithful port of findGQ for one strand; islands are 1-based starts."""
    nIls = len(strt)
    starts, ends, nums, subs, lens = [], [], [], [], []
    i = 0
    while i < nIls:
        conIls = 1
        npos = (lng[i] + 1) // (min_rep + 1)
        i2 = i + 1
        while i2 < nIls and (strt[i2] - (strt[i2 - 1] + lng[i2 - 1])) <= max_spacer:
            conIls += 1
            npos += (lng[i2] + 1) // (min_rep + 1)
            i2 += 1
        if npos >= 4:
            maxGQ = min_rep
            for j in range(i, i2):
                k = lng[j]
                while k > maxGQ:
                    nposMax = (lng[j] + 1) // (k + 1)
                    m = j + 1
                    while m < i2:
                        nposMax += (lng[m] + 1) // (k + 1)
                        if nposMax >= 4:
                            maxGQ = k
                            break
                        if (lng[m] + 1) // (k + 1) == 0:
                            sm1 = strt[m + 1] if m + 1 < nIls else 0
                            if sm1 > (strt[m - 1] + lng[m - 1] + max_spacer):
                                break
                        m += 1
                    k -= 1
            starts.append(strt[i])
            ends.append(strt[i2 - 1] + lng[i2 - 1] - 1)
            nums.append(npos)
            subs.append(conIls)
            lens.append(maxGQ)
        i = i + conIls  # skip the islands merged into this motif
    return starts, ends, nums, subs, lens


try:
    from numba import njit

    @njit(cache=True)
    def _islands_numba(arr, target, min_rep, strt, lng, fill):
        n = arr.shape[0]
        ndx = 0
        run = 0
        for i in range(n + 1):
            c = arr[i] if i < n else 0
            if c == target:
                run += 1
            else:
                if run >= min_rep:
                    if fill:
                        strt[ndx] = i - run + 1
                        lng[ndx] = run
                    ndx += 1
                run = 0
        return ndx

    @njit(cache=True)
    def _findgq_numba(strt, lng, nIls, min_rep, max_spacer, ostart, oend, onum, osub, olen, fill):
        ndx = 0
        i = 0
        while i < nIls:
            conIls = 1
            npos = (lng[i] + 1) // (min_rep + 1)
            i2 = i + 1
            while i2 < nIls and (strt[i2] - (strt[i2 - 1] + lng[i2 - 1])) <= max_spacer:
                conIls += 1
                npos += (lng[i2] + 1) // (min_rep + 1)
                i2 += 1
            if npos >= 4:
                maxGQ = min_rep
                for j in range(i, i2):
                    k = lng[j]
                    while k > maxGQ:
                        nposMax = (lng[j] + 1) // (k + 1)
                        m = j + 1
                        while m < i2:
                            nposMax += (lng[m] + 1) // (k + 1)
                            if nposMax >= 4:
                                maxGQ = k
                                break
                            if (lng[m] + 1) // (k + 1) == 0:
                                sm1 = strt[m + 1] if m + 1 < nIls else 0
                                if sm1 > (strt[m - 1] + lng[m - 1] + max_spacer):
                                    break
                            m += 1
                        k -= 1
                if fill:
                    ostart[ndx] = strt[i]
                    oend[ndx] = strt[i2 - 1] + lng[i2 - 1] - 1
                    onum[ndx] = npos
                    osub[ndx] = conIls
                    olen[ndx] = maxGQ
                ndx += 1
            i = i + conIls
        return ndx

    _HAS_NUMBA = True
except ImportError:  # pragma: no cover
    _islands_numba = None
    _findgq_numba = None
    _HAS_NUMBA = False


def _scan_strand_numba(arr, target, min_rep, max_spacer):
    n = arr.shape[0]
    cap = n // (min_rep + 1) + 1
    strt = np.empty(cap, np.int64)
    lng = np.empty(cap, np.int64)
    nIls = _islands_numba(arr, target, min_rep, strt, lng, 1)
    strt = strt[:nIls]
    lng = lng[:nIls]
    z = np.empty(0, np.int64)
    count = _findgq_numba(strt, lng, nIls, min_rep, max_spacer, z, z, z, z, z, 0)
    ostart = np.empty(count, np.int64)
    oend = np.empty(count, np.int64)
    onum = np.empty(count, np.int64)
    osub = np.empty(count, np.int64)
    olen = np.empty(count, np.int64)
    _findgq_numba(strt, lng, nIls, min_rep, max_spacer, ostart, oend, onum, osub, olen, 1)
    return ostart, oend, onum, osub, olen


def _scan_strand_python(seq, target, min_rep, max_spacer):
    strt, lng = _islands_py(seq, target, min_rep)
    s, e, num, sub, ln = _findgq_py(strt, lng, min_rep, max_spacer)
    dt = np.int64
    return tuple(np.array(x, dt) for x in (s, e, num, sub, ln))


class GQuadruplex:
    """Detect G-quadruplex forming motifs like gfa ``gfa_GQ.gff``.

    :param fasta_file: input FASTA.
    :param min_rep: minimum guanine run length / G-island size (gfa default 3).
    :param max_spacer: maximum loop length between islands (gfa default 7).
    """

    def __init__(self, fasta_file, min_rep=3, max_spacer=7):
        self.fasta_file = fasta_file
        self.min_rep = min_rep
        self.max_spacer = max_spacer
        self.df = pd.DataFrame(columns=["seqid", "start", "end", "strand", "islands", "runs", "max", "sequence"])

    def run(self, progress=True):
        fa = FastA(self.fasta_file)
        frames = []
        for seqid in tqdm(fa.names, disable=not progress, desc="G-quadruplex", unit="seq"):
            seq_b = fa.sequences[fa.names.index(seqid)].upper().encode()
            # plus strand (G runs) then minus strand (C runs), matching gfa order
            for target, strand in ((_G, "+"), (_C, "-")):
                if _HAS_NUMBA:
                    arr = np.frombuffer(seq_b, dtype=np.uint8)
                    s, e, num, sub, ln = _scan_strand_numba(arr, target, self.min_rep, self.max_spacer)
                else:
                    s, e, num, sub, ln = _scan_strand_python(seq_b, target, self.min_rep, self.max_spacer)
                if len(s) == 0:
                    continue
                seqs = []
                for gs, ge in zip(s.tolist(), e.tolist()):
                    fwd = seq_b[gs - 1 : ge]  # gs is 1-based start, ge 1-based inclusive
                    if strand == "-":
                        fwd = fwd.translate(_COMP)[::-1]
                    seqs.append(fwd.decode())
                frames.append(
                    pd.DataFrame(
                        {
                            "seqid": seqid,
                            "start": s - 1,  # 0-based
                            "end": e,  # 0-based exclusive (== 1-based inclusive)
                            "strand": strand,
                            "islands": sub,
                            "runs": num,
                            "max": ln,
                            "sequence": seqs,
                        }
                    )
                )
        cols = ["seqid", "start", "end", "strand", "islands", "runs", "max", "sequence"]
        self.df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=cols)

    def to_bed(self, output_file, mode="w"):
        if self.df.empty:
            raise ValueError("No results. Run `.run()` first.")
        bed = self.df[["seqid", "start", "end"]].copy()
        bed["name"] = self.df["sequence"]
        bed["score"] = 0
        bed["strand"] = self.df["strand"]
        bed.to_csv(output_file, sep="\t", header=False, index=False, mode=mode)

    def to_gff(self, output_file, source="sequana"):
        """Write a gff3 comparable to gfa ``gfa_GQ.gff`` (1-based, inclusive)."""
        if self.df.empty:
            raise ValueError("No results. Run `.run()` first.")
        with open(output_file, "w") as fout:
            fout.write("##gff-version 3\n")
            for row in self.df.itertuples(index=False):
                start = row.start + 1  # gff is 1-based
                s = row.sequence.lower()
                comp = f"{s.count('a')}A/{s.count('c')}C/{s.count('g')}G/{s.count('t')}T"
                attrs = (
                    f"ID={row.seqid}_{start}_{row.end}_GQ;islands={row.islands};"
                    f"runs={row.runs};max={row.max};composition={comp};sequence={s}"
                )
                fout.write(
                    f"{row.seqid}\t{source}\tG_Quadruplex_Motif\t{start}\t{row.end}\t.\t{row.strand}\t.\t{attrs}\n"
                )
