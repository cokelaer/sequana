#!/usr/bin/env python3
"""Detect mirror repeats.

A mirror repeat is two arms that read the same forwards and backwards across a
central axis (the right arm is the **reverse** of the left arm, NOT the reverse
complement — that would be an inverted repeat / cruciform), optionally separated
by a spacer::

    5'- ARM ----- spacer ----- reverse(ARM) -3'

This reproduces the *Mirror_Repeat* output of the non-B gfa / nBMST tool
(``gfa_MR.gff``). The homopurine / homopyrimidine subset of mirror repeats forms
H-DNA (triplex); see :class:`sequana.hdna.HDNA`.

gfa defaults (from ``print_usage.c``): arm (repeat) >= 10, spacer 0..100.

Maximal arms are reported; hits fully contained in a longer hit are kept but
flagged (``subset=1``); ``include_subsets=False`` drops them.

Validated against gfa on a 30+ Mb *Leishmania* assembly: all 22254 gfa mirror
repeats are recovered (100% positional overlap), density within ~2% (22704 vs
22254). The numba kernel processes the whole genome in a few seconds.
"""
import os

import numpy as np
import pandas as pd
from tqdm import tqdm

from sequana import FastA

_A, _C, _G, _T = 65, 67, 71, 84
_MR_COLS = ["seqid", "start", "end", "length", "repeat", "spacer", "subset", "sequence"]


def _scan_python(arr, min_repeat, max_repeat, min_spacer, max_spacer):
    """Pure-numpy fallback mirroring the numba kernel."""
    n = arr.shape[0]
    starts, ends, reps, spacers = [], [], [], []
    for spacer in range(min_spacer, max_spacer + 1):
        limit = n - 1 - spacer
        if limit <= 0:
            break
        for a0 in np.nonzero(arr[:limit] == arr[spacer + 1 : n])[0].tolist():
            b0 = arr[a0]
            if b0 not in (_A, _C, _G, _T):  # N never seeds a mirror
                continue
            rb = a0 + spacer + 1
            k = 1
            while (
                a0 - k >= 0
                and rb + k < n
                and k < max_repeat
                and arr[a0 - k] == arr[rb + k]
                and arr[a0 - k] in (_A, _C, _G, _T)
            ):
                k += 1
            if k >= min_repeat:
                starts.append(a0 - k + 1)
                ends.append(rb + k)
                reps.append(k)
                spacers.append(spacer)
    dt = np.int64
    return np.array(starts, dt), np.array(ends, dt), np.array(reps, dt), np.array(spacers, dt)


try:
    from numba import njit

    @njit(cache=True)
    def _scan_numba(arr, min_repeat, max_repeat, min_spacer, max_spacer):
        n = arr.shape[0]
        starts = np.empty(0, np.int64)
        ends = np.empty(0, np.int64)
        reps = np.empty(0, np.int64)
        spacers = np.empty(0, np.int64)
        for fill in range(2):
            idx = 0
            for spacer in range(min_spacer, max_spacer + 1):
                limit = n - 1 - spacer
                if limit <= 0:
                    break
                for a0 in range(limit):
                    b0 = arr[a0]
                    if b0 != 65 and b0 != 67 and b0 != 71 and b0 != 84:  # N never seeds
                        continue
                    rb = a0 + spacer + 1
                    if arr[rb] != b0:  # mirror seed (same base, not complement)
                        continue
                    k = 1
                    while a0 - k >= 0 and rb + k < n and k < max_repeat:
                        lb = arr[a0 - k]
                        if lb != arr[rb + k] or (lb != 65 and lb != 67 and lb != 71 and lb != 84):
                            break
                        k += 1
                    if k >= min_repeat:
                        if fill == 1:
                            starts[idx] = a0 - k + 1
                            ends[idx] = rb + k
                            reps[idx] = k
                            spacers[idx] = spacer
                        idx += 1
            if fill == 0:
                starts = np.empty(idx, np.int64)
                ends = np.empty(idx, np.int64)
                reps = np.empty(idx, np.int64)
                spacers = np.empty(idx, np.int64)
        return starts, ends, reps, spacers

    _HAS_NUMBA = True
except ImportError:  # pragma: no cover
    _scan_numba = None
    _HAS_NUMBA = False


def _flag_subsets_one(starts, ends):
    """subset=1 for any hit contained in a longer one (single sequence)."""
    order = np.lexsort((-ends, starts))  # start asc, then end desc
    e = ends[order]
    flags = np.zeros(len(e), np.int64)
    max_end = -1
    for i in range(len(e)):
        if e[i] <= max_end:
            flags[i] = 1
        else:
            max_end = e[i]
    return order, flags


def _scan_one(item):
    """Worker (module-level, picklable). Scan one sequence; return a per-seq df
    with subset flags applied (sorted start asc / end desc)."""
    seqid, seq, args = item
    seq_b = seq.encode()
    arr = np.frombuffer(seq_b, dtype=np.uint8)
    scan = _scan_numba if _HAS_NUMBA else _scan_python
    starts, ends, reps, spacers = scan(arr, *args)
    if len(starts) == 0:
        return None
    order, flags = _flag_subsets_one(starts, ends)
    starts, ends, reps, spacers = starts[order], ends[order], reps[order], spacers[order]
    return pd.DataFrame(
        {
            "seqid": seqid,
            "start": starts,
            "end": ends,
            "length": ends - starts,
            "repeat": reps,
            "spacer": spacers,
            "subset": flags,
            "sequence": [seq_b[s:e].decode() for s, e in zip(starts.tolist(), ends.tolist())],
        }
    )


class MirrorRepeats:
    """Detect mirror repeats like gfa ``gfa_MR.gff``.

    :param fasta_file: input FASTA.
    :param min_repeat: minimum arm length (gfa default 10).
    :param max_repeat: maximum arm length (cap on arm extension).
    :param min_spacer: minimum loop length (gfa default 0).
    :param max_spacer: maximum loop length (gfa default 100).
    """

    def __init__(self, fasta_file, min_repeat=10, max_repeat=200, min_spacer=0, max_spacer=100):
        self.fasta_file = fasta_file
        self.min_repeat = min_repeat
        self.max_repeat = max_repeat
        self.min_spacer = min_spacer
        self.max_spacer = max_spacer
        self.df = pd.DataFrame(columns=["seqid", "start", "end", "length", "repeat", "spacer", "subset", "sequence"])

    def run(self, include_subsets=True, progress=True, processes=None):
        """Scan every sequence for mirror repeats.

        Sequences are independent, so the scan is spread over ``processes``
        workers (``None`` = all CPUs; ``1`` = serial).
        """
        fa = FastA(self.fasta_file)
        args = (self.min_repeat, self.max_repeat, self.min_spacer, self.max_spacer)
        items = [(seqid, fa.sequences[i].upper(), args) for i, seqid in enumerate(fa.names)]
        if processes is None:
            processes = os.cpu_count() or 1
        use_pool = _HAS_NUMBA and processes > 1 and len(items) > 1
        if use_pool:
            from multiprocessing import Pool

            with Pool(min(processes, len(items))) as pool:
                frames = list(
                    tqdm(
                        pool.imap(_scan_one, items),
                        total=len(items),
                        disable=not progress,
                        desc="Mirror repeats",
                        unit="seq",
                    )
                )
        else:
            frames = [_scan_one(it) for it in tqdm(items, disable=not progress, desc="Mirror repeats", unit="seq")]
        frames = [f for f in frames if f is not None]
        df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=_MR_COLS)
        if not include_subsets:
            df = df[df["subset"] == 0].reset_index(drop=True)
        self.df = df.reset_index(drop=True)

    def to_bed(self, output_file, mode="w"):
        if self.df.empty:
            raise ValueError("No results. Run `.run()` first.")
        bed = self.df[["seqid", "start", "end"]].copy()
        bed["name"] = self.df["sequence"]
        bed["score"] = 0
        bed["strand"] = "+"
        bed.to_csv(output_file, sep="\t", header=False, index=False, mode=mode)

    def to_gff(self, output_file, source="sequana"):
        """Write a gff3 comparable to gfa ``gfa_MR.gff`` (1-based, inclusive)."""
        if self.df.empty:
            raise ValueError("No results. Run `.run()` first.")
        with open(output_file, "w") as fout:
            fout.write("##gff-version 3\n")
            for row in self.df.itertuples(index=False):
                start = row.start + 1  # gff is 1-based
                attrs = (
                    f"ID={row.seqid}_{start}_{row.end}_MR;spacer={row.spacer};"
                    f"repeat={row.repeat};subset={row.subset};sequence={row.sequence.lower()}"
                )
                fout.write(f"{row.seqid}\t{source}\tMirror_Repeat\t{start}\t{row.end}\t.\t+\t.\t{attrs}\n")
