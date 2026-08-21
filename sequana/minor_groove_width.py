#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2026 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Prediction of the DNA minor groove width (MGW) from sequence.

The minor groove width is a DNA shape feature that strongly influences
protein-DNA recognition (electrostatic potential in the minor groove,
recognition of AT tracts, etc.). It cannot be read directly from the sequence
but can be predicted with the pentamer *sliding-window* model introduced by the
Rohs laboratory:

- Zhou et al. (2013) *DNAshape: a method for the high-throughput prediction of
  DNA structural features on a genomic scale*, Nucleic Acids Research 41,
  W56-W62.
- Chiu et al. (2016) *DNAshapeR: an R/Bioconductor package for DNA shape
  prediction and feature encoding*, Bioinformatics 32, 1211-1213.

Each of the 512 unique pentamers (1024 pentamers folded by reverse-complement
symmetry) is associated with an average MGW value in Angstroms, derived from
all-atom Monte Carlo simulations. To predict the MGW along a sequence, a
5-nucleotide window is slid base by base and the tabulated value is assigned to
the **central** nucleotide. The two nucleotides at each end therefore have no
predicted value (returned as ``NaN``).
"""
from collections import defaultdict

import colorlog

from sequana.datatools import sequana_data
from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana.tools import reverse_complement

logger = colorlog.getLogger(__name__)


__all__ = ["MinorGrooveWidth", "minor_groove_width"]


class MinorGrooveWidth:
    """Predict the minor groove width (MGW) along a DNA sequence.

    The prediction relies on the DNAshape pentamer query table (see module
    documentation). Values are expressed in Angstroms.

    ::

        from sequana.minor_groove_width import MinorGrooveWidth
        mgw = MinorGrooveWidth()
        values = mgw.predict("GATTACAGATTACA")
        mgw.plot("GATTACAGATTACA")

    The first two and last two positions are ``NaN`` since a full pentamer
    centred on those positions is not available. Windows containing a base
    outside ``A/C/G/T`` (e.g. ``N``) also yield ``NaN``.
    """

    #: length of the sliding window (pentamer)
    WINDOW = 5

    def __init__(self, table=None):
        """
        :param table: optional path to a custom pentamer query table (CSV with
            two columns ``pentamer,mgw``). If not provided, the bundled
            DNAshape table is used.
        """
        self._table = self._load_table(table)

    def _load_table(self, table=None):
        if table is None:
            table = sequana_data("dnashape_mgw_pentamers.csv")
        df = pd.read_csv(table)
        df["pentamer"] = df["pentamer"].str.upper()
        lookup = dict(zip(df["pentamer"], df["mgw"].astype(float)))

        # Fold by reverse-complement symmetry so that any of the 1024 pentamers
        # can be looked up even when only one strand is tabulated.
        full = defaultdict(lambda: np.nan)
        for pentamer, value in lookup.items():
            full[pentamer] = value
            rc = reverse_complement(pentamer)
            if rc not in lookup:
                full[rc] = value
        logger.debug(f"Loaded {len(lookup)} pentamers ({len(full)} with reverse complements)")
        return full

    def pentamer_mgw(self, pentamer):
        """Return the tabulated MGW (Angstroms) of a single pentamer.

        Reverse-complement symmetry is used, so a pentamer and its reverse
        complement share the same value. Returns ``NaN`` when the pentamer is
        not exactly 5 ``A/C/G/T`` bases.

        ::

            >>> from sequana.minor_groove_width import MinorGrooveWidth
            >>> MinorGrooveWidth().pentamer_mgw("AAAAA")
            3.38
        """
        pentamer = pentamer.upper()
        if len(pentamer) != self.WINDOW or set(pentamer) - set("ACGT"):
            return np.nan
        return self._table[pentamer]

    def predict(self, sequence):
        """Predict the MGW profile along a sequence.

        :param str sequence: a DNA sequence.
        :return: a numpy array of the same length as ``sequence``. The two
            flanking positions on each side are ``NaN`` and so is any position
            whose centred pentamer contains a non-``ACGT`` base.
        :rtype: numpy.ndarray
        """
        sequence = sequence.upper()
        N = len(sequence)
        values = np.full(N, np.nan)
        half = self.WINDOW // 2
        for i in range(N - self.WINDOW + 1):
            values[i + half] = self.pentamer_mgw(sequence[i : i + self.WINDOW])
        return values

    def predict_fasta(self, filename):
        """Predict the MGW profile for every sequence of a FastA file.

        :param str filename: path to a FastA file.
        :return: dictionary mapping each sequence name to its MGW numpy array.
        :rtype: dict
        """
        from sequana.fasta import FastA

        fasta = FastA(filename)
        return {name: self.predict(seq) for name, seq in zip(fasta.names, fasta.sequences)}

    def plot(self, sequence, ax=None, marker="o", color="C0", **kwargs):
        """Plot the MGW profile of a sequence.

        :param str sequence: a DNA sequence.
        :param ax: matplotlib axes. If None, the current axes is used.
        :param str marker: matplotlib marker.
        :param color: line/marker color.
        :return: the numpy array of MGW values.
        """
        values = self.predict(sequence)
        if ax is None:
            ax = pylab.gca()
        positions = np.arange(len(sequence))
        ax.plot(positions, values, marker=marker, color=color, **kwargs)
        ax.set_xlabel("Position (bp)")
        ax.set_ylabel("Minor groove width (Å)")
        ax.grid(True)
        return values


def minor_groove_width(sequence, table=None):
    """Convenience wrapper to predict the MGW profile of a sequence.

    :param str sequence: a DNA sequence.
    :param table: optional path to a custom pentamer query table.
    :return: numpy array of MGW values (Angstroms), with ``NaN`` at the two
        flanking positions on each side.

    ::

        >>> from sequana.minor_groove_width import minor_groove_width
        >>> values = minor_groove_width("AAAAA")
        >>> round(float(values[2]), 2)
        3.38
    """
    return MinorGrooveWidth(table=table).predict(sequence)
