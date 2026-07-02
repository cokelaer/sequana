#!/usr/bin/env python3
"""Detect Z-DNA forming motifs (alternating purine/pyrimidine runs).

Z-DNA is favoured by runs of alternating purine/pyrimidine dinucleotides,
specifically ``CG/GC`` (strong) and ``CA/TG/AC/GT`` (weak); ``TA/AT`` are
excluded. This reproduces the *Z_DNA_Motif* output of the non-B gfa / nBMST
tool (``gfa_Z.gff``), ported faithfully from its ``findZDNA.c``.

Each valid dinucleotide contributes to a Kadane/Vasquez (KV) style score
(``CG/GC`` = 25, the weak ones = 3). A run of at least ``min_z`` consecutive
alternating bases is reported with ``score = sum_of_dinucleotide_scores // 2``.
A motif is flagged ``subset=1`` when its score reaches ``min_kvscore`` (33).

gfa defaults (from ``gfa.c``): min_z = 10, min_kvscore = 33 (the comment in
gfa's ``is_subset.c`` saying 35 is stale; ``gfa.c`` uses 33).

Validated against gfa on a 30+ Mb *Leishmania* assembly: identical coordinates,
length, score and subset for all 60076 motifs. The only field that differs is
``composition`` because gfa counts bases over ``[start, end-1)`` (an off-by-one
that drops the last base); this module counts the full motif, so its A/C/G/T
counts are the correct ones and match the reported ``sequence``.
"""
import numpy as np
import pandas as pd
from tqdm import tqdm

from sequana import FastA


def _pupy_py(a, b):
    if a == "A":
        return 3 if b == "C" else 0
    if a == "T":
        return 3 if b == "G" else 0
    if a == "C":
        if b == "G":
            return 25
        return 3 if b == "A" else 0
    if a == "G":
        if b == "C":
            return 25
        return 3 if b == "T" else 0
    return 0


def _scan_python(arr, min_z):
    n = arr.shape[0]
    seq = arr.tobytes().decode()
    starts, ends, lens, scores = [], [], [], []
    npy = 1
    kvsum = 0
    i = 0
    while i < n - min_z:
        tmp = _pupy_py(seq[i], seq[i + 1])
        if tmp > 0:
            npy += 1
            kvsum += tmp
        else:
            if npy >= min_z:
                starts.append(i - npy + 1)
                ends.append(i + 1)
                lens.append(npy)
                scores.append(kvsum // 2)
            npy = 1
            kvsum = 0
        i += 1
    dt = np.int64
    return tuple(np.array(x, dt) for x in (starts, ends, lens, scores))


try:
    from numba import njit

    @njit(cache=True)
    def _pupy(a, b):
        if a == 65:  # A
            return 3 if b == 67 else 0  # AC
        if a == 84:  # T
            return 3 if b == 71 else 0  # TG
        if a == 67:  # C
            if b == 71:  # CG
                return 25
            return 3 if b == 65 else 0  # CA
        if a == 71:  # G
            if b == 67:  # GC
                return 25
            return 3 if b == 84 else 0  # GT
        return 0

    @njit(cache=True)
    def _scan_numba(arr, min_z):
        n = arr.shape[0]
        starts = np.empty(0, np.int64)
        ends = np.empty(0, np.int64)
        lens = np.empty(0, np.int64)
        scores = np.empty(0, np.int64)
        for fill in range(2):
            ndx = 0
            npy = 1
            kvsum = 0
            i = 0
            while i < n - min_z:
                tmp = _pupy(arr[i], arr[i + 1])
                if tmp > 0:
                    npy += 1
                    kvsum += tmp
                else:
                    if npy >= min_z:
                        if fill == 1:
                            starts[ndx] = i - npy + 1
                            ends[ndx] = i + 1
                            lens[ndx] = npy
                            scores[ndx] = kvsum // 2
                        ndx += 1
                    npy = 1
                    kvsum = 0
                i += 1
            if fill == 0:
                starts = np.empty(ndx, np.int64)
                ends = np.empty(ndx, np.int64)
                lens = np.empty(ndx, np.int64)
                scores = np.empty(ndx, np.int64)
        return starts, ends, lens, scores

    _HAS_NUMBA = True
except ImportError:  # pragma: no cover
    _scan_numba = None
    _HAS_NUMBA = False


class ZDNA:
    """Detect Z-DNA forming motifs like gfa ``gfa_Z.gff``.

    :param fasta_file: input FASTA.
    :param min_z: minimum run length in bp (gfa default 10).
    :param min_kvscore: KV score threshold for the ``subset`` flag (gfa
        default 33); a motif with ``score >= min_kvscore`` gets ``subset=1``.
    """

    def __init__(self, fasta_file, min_z=10, min_kvscore=33):
        self.fasta_file = fasta_file
        self.min_z = min_z
        self.min_kvscore = min_kvscore
        self.df = pd.DataFrame(columns=["seqid", "start", "end", "length", "score", "subset", "sequence"])

    def run(self, progress=True):
        fa = FastA(self.fasta_file)
        scan = _scan_numba if _HAS_NUMBA else _scan_python
        frames = []
        for seqid in tqdm(fa.names, disable=not progress, desc="Z-DNA", unit="seq"):
            seq_b = fa.sequences[fa.names.index(seqid)].upper().encode()
            arr = np.frombuffer(seq_b, dtype=np.uint8)
            starts, ends, lens, scores = scan(arr, self.min_z)
            if len(starts) == 0:
                continue
            subset = (scores >= self.min_kvscore).astype(np.int64)
            frames.append(
                pd.DataFrame(
                    {
                        "seqid": seqid,
                        "start": starts,
                        "end": ends,
                        "length": lens,
                        "score": scores,
                        "subset": subset,
                        "sequence": [seq_b[s:e].decode() for s, e in zip(starts.tolist(), ends.tolist())],
                    }
                )
            )
        cols = ["seqid", "start", "end", "length", "score", "subset", "sequence"]
        self.df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=cols)

    def to_bed(self, output_file, mode="w"):
        if self.df.empty:
            raise ValueError("No results. Run `.run()` first.")
        bed = self.df[["seqid", "start", "end"]].copy()
        bed["name"] = self.df["sequence"]
        bed["score"] = self.df["score"]
        bed["strand"] = "+"
        bed.to_csv(output_file, sep="\t", header=False, index=False, mode=mode)

    def to_gff(self, output_file, source="sequana"):
        """Write a gff3 comparable to gfa ``gfa_Z.gff`` (1-based, inclusive)."""
        if self.df.empty:
            raise ValueError("No results. Run `.run()` first.")
        with open(output_file, "w") as fout:
            fout.write("##gff-version 3\n")
            for row in self.df.itertuples(index=False):
                start = row.start + 1  # gff is 1-based
                s = row.sequence.lower()
                comp = f"{s.count('a')}A/{s.count('c')}C/{s.count('g')}G/{s.count('t')}T"
                attrs = (
                    f"ID={row.seqid}_{start}_{row.end}_ZDNA;length={row.length};"
                    f"score={row.score};composition={comp};sequence={s};subset={row.subset}"
                )
                fout.write(f"{row.seqid}\t{source}\tZ_DNA_Motif\t{start}\t{row.end}\t.\t+\t.\t{attrs}\n")
