import os

import numpy as np
import pandas as pd
from tqdm import tqdm

from sequana import FastA
from sequana.tools import reverse_complement

_IR_COLS = ["seqid", "start", "end", "length", "stem", "spacer", "subset", "sequence"]

# Single-base complement. Bases not in this table (e.g. N) never pair, so a
# window containing them is never reported as an inverted repeat.
_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

# Byte translation table: complement ACGT (upper), everything else -> 0x00.
# A complemented base 0x00 never equals a real base, so N and other ambiguous
# bases can never pair.
_COMP_TABLE = bytearray(256)
for _a, _b in zip(b"ACGT", b"TGCA"):
    _COMP_TABLE[_a] = _b
_COMP_TABLE = bytes(_COMP_TABLE)


def _scan_python(arr, carr, min_spacer, max_spacer, short_spacer_max, min_stem, min_stem_long, max_stem):
    """Pure-numpy fallback. Vectorise the seed test per spacer, extend seeds."""
    n = arr.shape[0]
    starts, ends, stems, spacers = [], [], [], []
    for spacer in range(min_spacer, max_spacer + 1):
        limit = n - 1 - spacer
        if limit <= 0:
            break
        ms = min_stem if spacer <= short_spacer_max else min_stem_long
        for a0 in np.nonzero(carr[:limit] == arr[spacer + 1 : n])[0].tolist():
            k = 1
            rb = a0 + spacer + 1
            while a0 - k >= 0 and rb + k < n and k < max_stem and carr[a0 - k] == arr[rb + k]:
                k += 1
            if k >= ms:
                starts.append(a0 - k + 1)
                ends.append(rb + k)
                stems.append(k)
                spacers.append(spacer)
    dt = np.int64
    return np.array(starts, dt), np.array(ends, dt), np.array(stems, dt), np.array(spacers, dt)


try:
    from numba import njit

    @njit(cache=True)
    def _scan_numba(arr, carr, min_spacer, max_spacer, short_spacer_max, min_stem, min_stem_long, max_stem):
        n = arr.shape[0]
        # Two passes: count hits to size the output, then fill. Avoids growing
        # arrays inside nopython code.
        starts = np.empty(0, np.int64)
        ends = np.empty(0, np.int64)
        stems = np.empty(0, np.int64)
        spacers = np.empty(0, np.int64)
        for fill in range(2):
            idx = 0
            for spacer in range(min_spacer, max_spacer + 1):
                limit = n - 1 - spacer
                if limit <= 0:
                    break
                ms = min_stem if spacer <= short_spacer_max else min_stem_long
                for a0 in range(limit):
                    if carr[a0] != arr[a0 + spacer + 1]:
                        continue
                    k = 1
                    rb = a0 + spacer + 1
                    while a0 - k >= 0 and rb + k < n and k < max_stem and carr[a0 - k] == arr[rb + k]:
                        k += 1
                    if k >= ms:
                        if fill == 1:
                            starts[idx] = a0 - k + 1
                            ends[idx] = rb + k
                            stems[idx] = k
                            spacers[idx] = spacer
                        idx += 1
            if fill == 0:
                starts = np.empty(idx, np.int64)
                ends = np.empty(idx, np.int64)
                stems = np.empty(idx, np.int64)
                spacers = np.empty(idx, np.int64)
        return starts, ends, stems, spacers

    _HAS_NUMBA = True
except ImportError:  # pragma: no cover
    _scan_numba = None
    _HAS_NUMBA = False


def _flag_subsets_one(starts, ends):
    """subset=1 for any hit contained in a longer one (single sequence).

    Sort by start ascending, end descending; sweep the running max end.
    Returns the sort order (indices) and the subset flags in that order.
    """
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
    """Worker (module-level, picklable). Scan one sequence, return a per-seq df
    with subset flags already applied (sorted start asc / end desc)."""
    seqid, seq, args = item
    seq_b = seq.encode()
    arr = np.frombuffer(seq_b, dtype=np.uint8)
    carr = np.frombuffer(seq_b.translate(_COMP_TABLE), dtype=np.uint8)  # non-ACGT -> 0x00
    scan = _scan_numba if _HAS_NUMBA else _scan_python
    starts, ends, stems, spacers = scan(arr, carr, *args)
    if len(starts) == 0:
        return None
    order, flags = _flag_subsets_one(starts, ends)
    starts, ends, stems, spacers = starts[order], ends[order], stems[order], spacers[order]
    return pd.DataFrame(
        {
            "seqid": seqid,
            "start": starts,
            "end": ends,
            "length": ends - starts,
            "stem": stems,
            "spacer": spacers,
            "subset": flags,
            "sequence": [seq_b[s:e].decode() for s, e in zip(starts.tolist(), ends.tolist())],
        }
    )


class Cruciforms:
    """Detect inverted repeats (cruciform-forming motifs).

    An inverted repeat is two arms that are reverse-complements of each other,
    optionally separated by a spacer (loop)::

        5'- ARM ----- spacer ----- revcomp(ARM) -3'

    This mirrors the *Inverted_Repeat* output of the non-B gfa / nBMST tool
    (``gfa_IR.gff``). Differences from a naive adjacent-arm scan:

    - a **spacer** (loop) between the two arms is allowed (``min_spacer`` ..
      ``max_spacer``); gfa default is 0..100,
    - for every (gap, spacer) center the **maximal** stem is reported instead
      of one hit per stem length, avoiding nested-register duplicates,
    - hits fully contained in a longer hit are kept but flagged (``subset=1``),
      matching gfa's ``subset`` attribute. Set ``include_subsets=False`` to drop
      them.

    gfa couples the loop and arm lengths: a short loop tolerates a short arm,
    a long loop requires a long arm. The two-tier default is reproduced here:
    spacer <= ``short_spacer_max`` (4) needs stem >= ``min_stem_len`` (6),
    otherwise stem >= ``min_stem_len_long`` (10).

    Validated against gfa on a 30+ Mb *Leishmania* assembly: 57,702 non-subset
    inverted repeats vs gfa's 57,740 (0.07% difference), with 100% reciprocal
    positional overlap on a test slice. With the numba kernel the whole genome
    is processed in ~25 s.

    :param fasta_file: input FASTA.
    :param min_stem_len: minimum arm length for short loops (gfa default 6).
    :param max_stem_len: maximum arm length (cap on stem extension).
    :param min_spacer: minimum loop length (gfa default 0).
    :param max_spacer: maximum loop length (gfa default 100).
    :param short_spacer_max: largest loop still allowed with ``min_stem_len``.
    :param min_stem_len_long: minimum arm length once spacer > short_spacer_max.

    .. note:: cost scales as ``len(genome) * (max_spacer + 1)``. For a whole
        genome lower ``max_spacer`` or run per-chromosome if too slow.
    """

    def __init__(
        self,
        fasta_file,
        min_stem_len=6,
        max_stem_len=115,
        min_spacer=0,
        max_spacer=100,
        short_spacer_max=4,
        min_stem_len_long=10,
    ):
        self.fasta_file = fasta_file
        self.min_stem_len = min_stem_len
        self.max_stem_len = max_stem_len
        self.min_spacer = min_spacer
        self.max_spacer = max_spacer
        self.short_spacer_max = short_spacer_max
        self.min_stem_len_long = min_stem_len_long
        self.df = pd.DataFrame(columns=["seqid", "start", "end", "length", "stem", "spacer", "subset", "sequence"])

    def is_cruciform(self, left, right):
        return left == reverse_complement(right)

    def run(self, include_subsets=True, progress=True, processes=None):
        """Scan every sequence for inverted repeats.

        Uses a numba-compiled kernel when numba is installed, otherwise a
        numpy-vectorised fallback. Sequences are independent, so the scan is
        spread over ``processes`` workers (``None`` = all CPUs; ``1`` = serial).
        """
        fa = FastA(self.fasta_file)
        args = (
            self.min_spacer,
            self.max_spacer,
            self.short_spacer_max,
            self.min_stem_len,
            self.min_stem_len_long,
            self.max_stem_len,
        )
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
                        desc="Inverted repeats",
                        unit="seq",
                    )
                )
        else:
            frames = [_scan_one(it) for it in tqdm(items, disable=not progress, desc="Inverted repeats", unit="seq")]
        frames = [f for f in frames if f is not None]
        df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=_IR_COLS)
        if not include_subsets:
            df = df[df["subset"] == 0].reset_index(drop=True)
        self.df = df.reset_index(drop=True)

    def to_bed(self, output_file):
        if self.df.empty:
            raise ValueError("No results. Run `.run()` first.")
        bed = self.df[["seqid", "start", "end"]].copy()
        bed["name"] = self.df["sequence"]
        bed["score"] = 0
        bed["strand"] = "+"
        bed.to_csv(output_file, sep="\t", header=False, index=False)

    def to_gff(self, output_file, source="sequana"):
        """Write a gff3 comparable to gfa ``gfa_IR.gff`` (1-based, inclusive)."""
        if self.df.empty:
            raise ValueError("No results. Run `.run()` first.")
        with open(output_file, "w") as fout:
            fout.write("##gff-version 3\n")
            for row in self.df.itertuples(index=False):
                start = row.start + 1  # gff is 1-based
                attrs = (
                    f"ID={row.seqid}_{start}_{row.end}_IR;spacer={row.spacer};"
                    f"repeat={row.stem};subset={row.subset};sequence={row.sequence.lower()}"
                )
                fout.write(f"{row.seqid}\t{source}\tInverted_Repeat\t{start}\t{row.end}\t.\t+\t.\t{attrs}\n")
