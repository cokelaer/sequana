#!/usr/bin/env python3
"""Detect short tandem repeats (STR / microsatellites).

A short tandem repeat is a short unit (period 1..9 bp) repeated head-to-tail at
least ``min_reps`` times, spanning at least ``min_span`` bp (including a partial
trailing copy). This reproduces the *Short_Tandem_Repeat* output of the non-B
gfa / nBMST tool (``gfa_STR.gff``), ported faithfully from its ``findSTR.c``.

gfa defaults (from ``gfa.c``): period 1..9, min_reps 3, min_span 10. For each
start position the smallest period giving a qualifying repeat is taken, then the
scan jumps past it; a new STR is kept only if it ends beyond the previous one.

The ``type`` attribute is gfa's ``nonBstr`` structural code, a 4-bit value
``8*isComplementary + 4*isSymmetric + 2*isAlternatingRY + 1*isEven`` describing
which non-B structure the unit can adopt.

Validated against gfa on a 30+ Mb *Leishmania* assembly: identical coordinates
and identical length/num/remainder/type attributes (exact match).
"""
import numpy as np
import pandas as pd
from tqdm import tqdm

from sequana import FastA


def _nonbstr_py(unit):
    L = len(unit)
    is_even = L % 2 == 0
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    if L >= 2:
        is_sym = True
        is_comp = is_even
        is_pupy = True
        for ii in range(L // 2):
            if unit[ii] != unit[L - 1 - ii]:
                is_sym = False
            if is_even and comp.get(unit[ii]) != unit[L - 1 - ii]:
                is_comp = False
        for k in range(1, L):
            a, b = unit[k], unit[k - 1]
            if a in "AG" and b in "AG":
                is_pupy = False
            if a in "TC" and b in "TC":
                is_pupy = False
    else:
        is_sym, is_comp, is_pupy = True, False, False
    return (1 if is_even else 0) + (2 if is_pupy else 0) + (4 if is_sym else 0) + (8 if is_comp else 0)


def _scan_python(arr, min_period, max_period, min_span, min_reps):
    n = arr.shape[0]
    seq = arr.tobytes().decode()
    starts, ends, periods, nums, subs, types = [], [], [], [], [], []
    last_end = -1
    i = 0
    while i < n - min_span:
        if seq[i] == "N":
            i += 1
            continue
        for rpsz in range(min_period, max_period + 1):
            reps = 1
            j = i + rpsz
            while j + rpsz <= n and seq[i : i + rpsz] == seq[j : j + rpsz]:
                reps += 1
                j += rpsz
                if j + rpsz >= n:
                    break
            if reps >= min_reps:
                remainder = 0
                rs, re = i, j
                while re < n and seq[rs] == seq[re]:
                    remainder += 1
                    rs += 1
                    re += 1
                if reps * rpsz + remainder >= min_span and (not starts or last_end < re):
                    starts.append(i)
                    ends.append(re)
                    periods.append(rpsz)
                    nums.append(reps)
                    subs.append(remainder)
                    types.append(_nonbstr_py(seq[i : i + rpsz]))
                    last_end = re
                    i = re - min_span + 1
                    break
        i += 1
    dt = np.int64
    return tuple(np.array(x, dt) for x in (starts, ends, periods, nums, subs, types))


try:
    from numba import njit

    @njit(cache=True)
    def _nonbstr(arr, start, length):
        is_even = length % 2 == 0
        if length >= 2:
            is_sym = True
            is_comp = is_even
            is_pupy = True
            for ii in range(length // 2):
                a = arr[start + ii]
                b = arr[start + length - 1 - ii]
                if a != b:
                    is_sym = False
                if is_even:
                    c = 0
                    if a == 65:
                        c = 84
                    elif a == 84:
                        c = 65
                    elif a == 67:
                        c = 71
                    elif a == 71:
                        c = 67
                    if c != b:
                        is_comp = False
            for k in range(1, length):
                a = arr[start + k]
                b = arr[start + k - 1]
                if (a == 65 or a == 71) and (b == 65 or b == 71):
                    is_pupy = False
                if (a == 84 or a == 67) and (b == 84 or b == 67):
                    is_pupy = False
        else:
            is_sym, is_comp, is_pupy = True, False, False
        code = 0
        if is_even:
            code += 1
        if is_pupy:
            code += 2
        if is_sym:
            code += 4
        if is_comp:
            code += 8
        return code

    @njit(cache=True)
    def _scan_numba(arr, min_period, max_period, min_span, min_reps):
        n = arr.shape[0]
        N = 78
        starts = np.empty(0, np.int64)
        ends = np.empty(0, np.int64)
        periods = np.empty(0, np.int64)
        nums = np.empty(0, np.int64)
        subs = np.empty(0, np.int64)
        types = np.empty(0, np.int64)
        for fill in range(2):
            ndx = 0
            last_end = -1
            i = 0
            while i < n - min_span:
                if arr[i] == N:
                    i += 1
                    continue
                rpsz = min_period
                while rpsz <= max_period:
                    reps = 1
                    j = i + rpsz
                    while j + rpsz <= n:
                        eq = True
                        for t in range(rpsz):
                            if arr[i + t] != arr[j + t]:
                                eq = False
                                break
                        if not eq:
                            break
                        reps += 1
                        j += rpsz
                        if j + rpsz >= n:
                            break
                    if reps >= min_reps:
                        remainder = 0
                        rs = i
                        re = j
                        while re < n and arr[rs] == arr[re]:
                            remainder += 1
                            rs += 1
                            re += 1
                        if reps * rpsz + remainder >= min_span and (ndx == 0 or last_end < re):
                            if fill == 1:
                                starts[ndx] = i
                                ends[ndx] = re
                                periods[ndx] = rpsz
                                nums[ndx] = reps
                                subs[ndx] = remainder
                                types[ndx] = _nonbstr(arr, i, rpsz)
                            ndx += 1
                            last_end = re
                            i = re - min_span + 1
                            break
                    rpsz += 1
                i += 1
            if fill == 0:
                starts = np.empty(ndx, np.int64)
                ends = np.empty(ndx, np.int64)
                periods = np.empty(ndx, np.int64)
                nums = np.empty(ndx, np.int64)
                subs = np.empty(ndx, np.int64)
                types = np.empty(ndx, np.int64)
        return starts, ends, periods, nums, subs, types

    _HAS_NUMBA = True
except ImportError:  # pragma: no cover
    _scan_numba = None
    _HAS_NUMBA = False


class ShortTandemRepeats:
    """Detect short tandem repeats (microsatellites) like gfa ``gfa_STR.gff``.

    :param fasta_file: input FASTA.
    :param min_period: minimum repeat-unit length (gfa default 1).
    :param max_period: maximum repeat-unit length (gfa default 9).
    :param min_span: minimum total span in bp (gfa default 10).
    :param min_reps: minimum number of full copies (gfa default 3).
    """

    def __init__(self, fasta_file, min_period=1, max_period=9, min_span=10, min_reps=3):
        self.fasta_file = fasta_file
        self.min_period = min_period
        self.max_period = max_period
        self.min_span = min_span
        self.min_reps = min_reps
        self.df = pd.DataFrame(columns=["seqid", "start", "end", "period", "num", "remainder", "type", "sequence"])

    def run(self, progress=True):
        fa = FastA(self.fasta_file)
        scan = _scan_numba if _HAS_NUMBA else _scan_python
        args = (self.min_period, self.max_period, self.min_span, self.min_reps)
        frames = []
        for seqid in tqdm(fa.names, disable=not progress, desc="Short tandem repeats", unit="seq"):
            seq_b = fa.sequences[fa.names.index(seqid)].upper().encode()
            arr = np.frombuffer(seq_b, dtype=np.uint8)
            starts, ends, periods, nums, subs, types = scan(arr, *args)
            if len(starts) == 0:
                continue
            frames.append(
                pd.DataFrame(
                    {
                        "seqid": seqid,
                        "start": starts,
                        "end": ends,
                        "period": periods,
                        "num": nums,
                        "remainder": subs,
                        "type": types,
                        "sequence": [seq_b[s:e].decode() for s, e in zip(starts.tolist(), ends.tolist())],
                    }
                )
            )
        cols = ["seqid", "start", "end", "period", "num", "remainder", "type", "sequence"]
        self.df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=cols)

    def to_bed(self, output_file, mode="w"):
        if self.df.empty:
            raise ValueError("No results. Run `.run()` first.")
        bed = self.df[["seqid", "start", "end"]].copy()
        bed["name"] = self.df["sequence"]
        bed["score"] = 0
        bed["strand"] = "+"
        bed.to_csv(output_file, sep="\t", header=False, index=False, mode=mode)

    def to_gff(self, output_file, source="sequana"):
        """Write a gff3 comparable to gfa ``gfa_STR.gff`` (1-based, inclusive)."""
        if self.df.empty:
            raise ValueError("No results. Run `.run()` first.")
        with open(output_file, "w") as fout:
            fout.write("##gff-version 3\n")
            for row in self.df.itertuples(index=False):
                start = row.start + 1  # gff is 1-based
                attrs = (
                    f"ID={row.seqid}_{start}_{row.end}_STR;length={row.period};"
                    f"x{row.num}+{row.remainder};type={row.type};sequence={row.sequence.lower()}"
                )
                fout.write(f"{row.seqid}\t{source}\tShort_Tandem_Repeat\t{start}\t{row.end}\t.\t+\t.\t{attrs}\n")
