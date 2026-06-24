#!/usr/bin/env python3
"""Detect direct (and slipped / tandem) repeats.

A direct repeat is two identical arms in the same orientation, optionally
separated by a spacer::

    5'- ARM ----- spacer ----- ARM -3'

When the spacer is zero the arm is extended into a big tandem repeat (BTR).
This reproduces the *Direct_Repeat* output of the non-B gfa / nBMST tool
(``gfa_DR.gff``), ported faithfully from its ``findDR.c``.

gfa defaults (from ``gfa.c``): min arm 10, max arm 300, spacer 0..10 (the usage
text saying spacer 100 is stale; ``gfa.c`` uses 10). The search is greedy: for
each start it tries the largest arm first and the smallest qualifying spacer,
records the first full match, and skips ahead so successive repeats must extend
beyond the previous one (the ``end`` escape variable in gfa).

A motif is flagged ``subset=1`` (gfa's "slipped" set) when its spacer is 0.

Validated against the gfa binary compiled from source on a 30+ Mb *Leishmania*
assembly: identical coordinates and attributes.
"""
import os

import numpy as np
import pandas as pd
from tqdm import tqdm

from sequana import FastA

_N = 78  # uppercase 'N'


def _scan_python(arr, mindir, maxdir, dspacer):
    n = arr.shape[0]
    starts, ends, sizes, spacers, nums, rems = [], [], [], [], [], []
    lasti = n - mindir * 2
    end = 0
    strti = 0
    while strti <= lasti:
        while strti < n and arr[strti] == _N:
            strti += 1
        if strti >= lasti:
            break
        size = maxdir
        while size >= mindir:
            if (size * 2 + dspacer) <= (end - strti):
                break  # even the largest config cannot escape previous DR
            sp_min = max(0, (end - strti) - (size * 2) + 2)
            sp_max = min(dspacer, lasti - strti)
            found = False
            sp = sp_min
            while sp <= sp_max:
                j = strti + size + sp
                i = strti
                k = 0
                while j < n and arr[i] == arr[j] and k < size and arr[i] != _N:
                    k += 1
                    j += 1
                    i += 1
                if k == size:
                    totlen = k
                    if sp == 0:
                        while j < n and arr[i] == arr[j]:
                            totlen += 1
                            j += 1
                            i += 1
                    starts.append(strti)
                    ends.append(j)
                    sizes.append(size)
                    spacers.append(sp)
                    nums.append(totlen // size)
                    rems.append(totlen % size)
                    end = j - 1
                    found = True
                    break
                sp += 1
            if found:
                break
            size -= 1
        strti += 1
    dt = np.int64
    return tuple(np.array(x, dt) for x in (starts, ends, sizes, spacers, nums, rems))


try:
    from numba import njit

    @njit(cache=True)
    def _scan_numba(arr, mindir, maxdir, dspacer):
        n = arr.shape[0]
        lasti = n - mindir * 2
        starts = np.empty(0, np.int64)
        ends = np.empty(0, np.int64)
        sizes = np.empty(0, np.int64)
        spacers = np.empty(0, np.int64)
        nums = np.empty(0, np.int64)
        rems = np.empty(0, np.int64)
        for fill in range(2):
            ndx = 0
            end = 0
            strti = 0
            while strti <= lasti:
                while strti < n and arr[strti] == _N:
                    strti += 1
                if strti >= lasti:
                    break
                # common path: no active escape (end behind), so sp_min is always 0
                # and the size-2*size break never fires -> hoist the arithmetic out
                escape = end > strti
                sp_max0 = dspacer
                if lasti - strti < sp_max0:
                    sp_max0 = lasti - strti
                size = maxdir
                while size >= mindir:
                    if escape:
                        if (size * 2 + dspacer) <= (end - strti):
                            break
                        sp_min = (end - strti) - (size * 2) + 2
                        if sp_min < 0:
                            sp_min = 0
                    else:
                        sp_min = 0
                    found = False
                    sp = sp_min
                    while sp <= sp_max0:
                        j = strti + size + sp
                        i = strti
                        k = 0
                        while j < n and arr[i] == arr[j] and k < size and arr[i] != _N:
                            k += 1
                            j += 1
                            i += 1
                        if k == size:
                            totlen = k
                            if sp == 0:
                                while j < n and arr[i] == arr[j]:
                                    totlen += 1
                                    j += 1
                                    i += 1
                            if fill == 1:
                                starts[ndx] = strti
                                ends[ndx] = j
                                sizes[ndx] = size
                                spacers[ndx] = sp
                                nums[ndx] = totlen // size
                                rems[ndx] = totlen % size
                            ndx += 1
                            end = j - 1
                            found = True
                            break
                        sp += 1
                    if found:
                        break
                    size -= 1
                strti += 1
            if fill == 0:
                starts = np.empty(ndx, np.int64)
                ends = np.empty(ndx, np.int64)
                sizes = np.empty(ndx, np.int64)
                spacers = np.empty(ndx, np.int64)
                nums = np.empty(ndx, np.int64)
                rems = np.empty(ndx, np.int64)
        return starts, ends, sizes, spacers, nums, rems

    _HAS_NUMBA = True
except ImportError:  # pragma: no cover
    _scan_numba = None
    _HAS_NUMBA = False


def _scan_one(item):
    """Worker (module-level so it can be pickled for multiprocessing).

    ``item`` = (seqid, upper-cased sequence string, (min_repeat, max_repeat,
    max_spacer)). Returns a per-sequence dict ready to become a DataFrame.
    """
    seqid, seq, args = item
    seq_b = seq.encode()
    arr = np.frombuffer(seq_b, dtype=np.uint8)
    scan = _scan_numba if _HAS_NUMBA else _scan_python
    starts, ends, sizes, spacers, nums, rems = scan(arr, *args)
    seqs = [seq_b[s:e].decode() for s, e in zip(starts.tolist(), ends.tolist())]
    return seqid, starts, ends, sizes, spacers, nums, rems, seqs


class DirectRepeats:
    """Detect direct/slipped repeats like gfa ``gfa_DR.gff``.

    :param fasta_file: input FASTA.
    :param min_repeat: minimum arm length (gfa default 10).
    :param max_repeat: maximum arm length (gfa default 300).
    :param max_spacer: maximum spacer between arms (gfa default 10).
    """

    def __init__(self, fasta_file, min_repeat=10, max_repeat=300, max_spacer=10):
        self.fasta_file = fasta_file
        self.min_repeat = min_repeat
        self.max_repeat = max_repeat
        self.max_spacer = max_spacer
        self.df = pd.DataFrame(
            columns=["seqid", "start", "end", "repeat", "spacer", "num", "remainder", "subset", "sequence"]
        )

    def run(self, progress=True, processes=None):
        """Scan every sequence for direct repeats.

        :param processes: number of worker processes. ``None`` (default) uses
            all CPUs; sequences are independent so this is an exact speed-up.
            Pass ``1`` to force a single in-process (serial) scan.
        """
        fa = FastA(self.fasta_file)
        args = (self.min_repeat, self.max_repeat, self.max_spacer)
        items = [(seqid, fa.sequences[i].upper(), args) for i, seqid in enumerate(fa.names)]
        if processes is None:
            processes = os.cpu_count() or 1
        # multiprocessing only pays off with numba and several sequences; a single
        # in-process scan keeps the numba/fallback monkeypatch honest for tests
        use_pool = _HAS_NUMBA and processes > 1 and len(items) > 1
        if use_pool:
            from multiprocessing import Pool

            with Pool(min(processes, len(items))) as pool:
                results = list(
                    tqdm(
                        pool.imap(_scan_one, items),
                        total=len(items),
                        disable=not progress,
                        desc="Direct repeats",
                        unit="seq",
                    )
                )
        else:
            results = [_scan_one(it) for it in tqdm(items, disable=not progress, desc="Direct repeats", unit="seq")]
        frames = []
        for seqid, starts, ends, sizes, spacers, nums, rems, seqs in results:
            if len(starts) == 0:
                continue
            frames.append(
                pd.DataFrame(
                    {
                        "seqid": seqid,
                        "start": starts,
                        "end": ends,
                        "repeat": sizes,
                        "spacer": spacers,
                        "num": nums,
                        "remainder": rems,
                        "subset": (spacers == 0).astype(np.int64),
                        "sequence": seqs,
                    }
                )
            )
        cols = ["seqid", "start", "end", "repeat", "spacer", "num", "remainder", "subset", "sequence"]
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
        """Write a gff3 comparable to gfa ``gfa_DR.gff`` (1-based, inclusive).

        Like gfa, ``composition`` counts only the first arm (``repeat`` bases)
        while ``sequence`` is the full motif span.
        """
        if self.df.empty:
            raise ValueError("No results. Run `.run()` first.")
        with open(output_file, "w") as fout:
            fout.write("##gff-version 3\n")
            for row in self.df.itertuples(index=False):
                start = row.start + 1  # gff is 1-based
                s = row.sequence.lower()
                arm = s[: row.repeat]
                comp = f"{arm.count('a')}A/{arm.count('c')}C/{arm.count('g')}G/{arm.count('t')}T"
                attrs = (
                    f"ID={row.seqid}_{start}_{row.end}_DR;spacer={row.spacer};"
                    f"repeat={row.repeat};x{row.num}+{row.remainder};"
                    f"composition={comp};sequence={s};subset={row.subset}"
                )
                fout.write(f"{row.seqid}\t{source}\tDirect_Repeat\t{start}\t{row.end}\t.\t+\t.\t{attrs}\n")
