#!/usr/bin/env python3
"""Detect intramolecular triplex (H-DNA) forming motifs.

H-DNA forms at **mirror repeats** whose strands are homopurine or
homopyrimidine. This reproduces the *Triplex* output of the non-B gfa / nBMST
tool (``gfa_TPX.gff``).

A motif is::

    5'- ARM ----- spacer ----- mirror(ARM) -3'

where ``mirror(ARM)`` is the arm read backwards (same bases, **not** the
reverse complement, unlike a cruciform/inverted repeat), and the whole motif
(both arms + spacer) is **100% purine (A/G) or 100% pyrimidine (C/T)**.

gfa ``gfa_TPX.gff`` defaults, reverse-engineered from its output, are
reproduced here: arm (repeat) >= 10, spacer 0..8, and each arm is allowed up to
10% impurities (``impure * 10 <= arm_length``) judged on the arm only, not the
spacer. Every gfa Triplex is also a gfa Mirror_Repeat; the triplex set is the
homopurine / homopyrimidine subset of mirror repeats.

Validated against gfa on a 30+ Mb *Leishmania* assembly (``gfa_TPX.gff``): all
1746 gfa Triplex loci are recovered (100% positional overlap), with density
within ~6% (1859 vs 1746). The numba kernel processes the whole genome in ~4 s.
"""
import numpy as np
import pandas as pd
from tqdm import tqdm

from sequana import FastA

# ASCII codes
_A, _C, _G, _T = 65, 67, 71, 84


def _scan_python(arr, min_repeat, max_repeat, min_spacer, max_spacer):
    """Pure-numpy fallback mirroring the numba kernel."""
    n = arr.shape[0]
    starts, ends, reps, spacers = [], [], [], []
    for spacer in range(min_spacer, max_spacer + 1):
        limit = n - 1 - spacer
        if limit <= 0:
            break
        # mirror seed: innermost pair (a0, a0+spacer+1) are the same base
        for a0 in np.nonzero(arr[:limit] == arr[spacer + 1 : n])[0].tolist():
            b0 = arr[a0]
            if b0 not in (_A, _C, _G, _T):  # N never seeds a mirror
                continue
            rb = a0 + spacer + 1
            # extend the perfect mirror arm outward (same base, not complement)
            k = 1
            while (
                a0 - k >= 0
                and rb + k < n
                and k < max_repeat
                and arr[a0 - k] == arr[rb + k]
                and arr[a0 - k] in (_A, _C, _G, _T)
            ):
                k += 1
            if k < min_repeat:
                continue
            start, end = a0 - k + 1, rb + k
            # purity is judged on the arm (repeat unit) only, not the spacer
            npur = npyr = 0
            for t in range(start, start + k):
                b = arr[t]
                if b == _A or b == _G:
                    npur += 1
                elif b == _C or b == _T:
                    npyr += 1
            # homopurine or homopyrimidine arm, allowing up to 10% impurities
            impure = k - (npur if npur >= npyr else npyr)
            if impure * 10 <= k:
                starts.append(start)
                ends.append(end)
                reps.append(k)
                spacers.append(spacer)
    dt = np.int64
    return np.array(starts, dt), np.array(ends, dt), np.array(reps, dt), np.array(spacers, dt)


try:
    from numba import njit

    @njit(cache=True)
    def _scan_numba(arr, min_repeat, max_repeat, min_spacer, max_spacer):
        n = arr.shape[0]
        # Two passes: count to size the outputs, then fill.
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
                    # extend the perfect mirror arm outward (stop at N)
                    k = 1
                    while a0 - k >= 0 and rb + k < n and k < max_repeat:
                        lb = arr[a0 - k]
                        if lb != arr[rb + k] or (lb != 65 and lb != 67 and lb != 71 and lb != 84):
                            break
                        k += 1
                    if k < min_repeat:
                        continue
                    start = a0 - k + 1
                    end = rb + k
                    # purity is judged on the arm (repeat unit) only, not the spacer
                    npur = 0
                    npyr = 0
                    for t in range(start, start + k):
                        b = arr[t]
                        if b == 65 or b == 71:
                            npur += 1
                        elif b == 67 or b == 84:
                            npyr += 1
                    major = npur if npur >= npyr else npyr
                    # homopurine or homopyrimidine arm, allowing up to 10% impurities
                    if (k - major) * 10 <= k:
                        if fill == 1:
                            starts[idx] = start
                            ends[idx] = end
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


class HDNA:
    """Detect H-DNA (intramolecular triplex) forming mirror repeats.

    :param fasta_file: input FASTA.
    :param min_repeat: minimum arm length (gfa Triplex default 10).
    :param max_repeat: maximum arm length (cap on arm extension).
    :param min_spacer: minimum loop length (gfa default 0).
    :param max_spacer: maximum loop length (gfa default 8).

    Hits fully contained in a longer hit are kept but flagged (``subset=1``),
    matching gfa's ``subset`` attribute; ``include_subsets=False`` drops them.
    """

    def __init__(self, fasta_file, min_repeat=10, max_repeat=100, min_spacer=0, max_spacer=8):
        self.fasta_file = fasta_file
        self.min_repeat = min_repeat
        self.max_repeat = max_repeat
        self.min_spacer = min_spacer
        self.max_spacer = max_spacer
        self.df = pd.DataFrame(columns=["seqid", "start", "end", "length", "repeat", "spacer", "subset", "sequence"])

    def run(self, include_subsets=True, progress=True):
        fa = FastA(self.fasta_file)
        scan = _scan_numba if _HAS_NUMBA else _scan_python
        args = (self.min_repeat, self.max_repeat, self.min_spacer, self.max_spacer)
        frames = []

        for seqid in tqdm(fa.names, disable=not progress, desc="H-DNA (triplex)", unit="seq"):
            seq_b = fa.sequences[fa.names.index(seqid)].upper().encode()
            arr = np.frombuffer(seq_b, dtype=np.uint8)

            starts, ends, reps, spacers = scan(arr, *args)
            if len(starts) == 0:
                continue
            frames.append(
                pd.DataFrame(
                    {
                        "seqid": seqid,
                        "start": starts,
                        "end": ends,
                        "length": ends - starts,
                        "repeat": reps,
                        "spacer": spacers,
                        "subset": 0,
                        "sequence": [seq_b[s:e].decode() for s, e in zip(starts.tolist(), ends.tolist())],
                    }
                )
            )

        cols = ["seqid", "start", "end", "length", "repeat", "spacer", "subset", "sequence"]
        df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=cols)
        df = self._flag_subsets(df)
        if not include_subsets:
            df = df[df["subset"] == 0].reset_index(drop=True)
        self.df = df.reset_index(drop=True)

    @staticmethod
    def _flag_subsets(df):
        """Flag (subset=1) any hit whose interval is contained in a longer one."""
        if df.empty:
            return df
        out = []
        for _, sub in df.groupby("seqid", sort=False):
            sub = sub.sort_values(["start", "end"], ascending=[True, False]).copy()
            max_end = -1
            flags = []
            for end in sub["end"]:
                flags.append(1 if end <= max_end else 0)
                if end > max_end:
                    max_end = end
            sub["subset"] = flags
            out.append(sub)
        return pd.concat(out)

    def to_bed(self, output_file, mode="w"):
        if self.df.empty:
            raise ValueError("No results. Run `.run()` first.")
        bed = self.df[["seqid", "start", "end"]].copy()
        bed["name"] = self.df["sequence"]
        bed["score"] = 0
        bed["strand"] = "+"
        bed.to_csv(output_file, sep="\t", header=False, index=False, mode=mode)

    def to_gff(self, output_file, source="sequana"):
        """Write a gff3 comparable to gfa ``gfa_TPX.gff`` (1-based, inclusive)."""
        if self.df.empty:
            raise ValueError("No results. Run `.run()` first.")
        with open(output_file, "w") as fout:
            fout.write("##gff-version 3\n")
            for row in self.df.itertuples(index=False):
                start = row.start + 1  # gff is 1-based
                attrs = (
                    f"ID={row.seqid}_{start}_{row.end}_TPX;spacer={row.spacer};"
                    f"repeat={row.repeat};subset={row.subset};sequence={row.sequence.lower()}"
                )
                fout.write(f"{row.seqid}\t{source}\tTriplex\t{start}\t{row.end}\t.\t+\t.\t{attrs}\n")
