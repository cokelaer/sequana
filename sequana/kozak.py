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
import functools
import os
import random
import re
from collections import Counter, defaultdict
from contextlib import contextmanager

import colorlog
from tqdm import tqdm

from sequana.fasta import FastA
from sequana.gff3 import GFF3
from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana.tools import reverse_complement

logger = colorlog.getLogger(__name__)


class Motif:
    def __init__(self, motif):
        if isinstance(motif, str):
            self.df = pd.read_csv(motif, index_col=0)
        else:
            self.df = motif

    def _get_entropy(self, base=2):
        """Compute Shannon entropy for each position in the motif.

        Shannon Entropy H(p) = -sum(p * log_base(p)) for DNA (4-letter alphabet).

        :param base: the base of the logarithm. Common values: 2 (bits, default),
            e (nats), 10 (bans).
        :return: numpy array of entropy values per position.
        """
        from scipy.stats import entropy

        freq = self.df[["A", "C", "G", "T"]].to_numpy()
        freq = freq / freq.sum(axis=1, keepdims=True)
        H = np.array([entropy(row, base=base) for row in freq])
        return H

    entropy = property(_get_entropy)

    def _get_information_content(self, base=2):
        """Compute information content (bits) for each position in the motif.

        Information Content (IC) = max_entropy - Shannon_entropy = 2 - H(p)
        for DNA (4-letter alphabet).

        :param base: the base of the logarithm. Common values: 2 (bits, default),
            e (nats), 10 (bans).
        :return: numpy array of information content values per position.
        """
        return 2 - self.entropy

    information_content = property(_get_information_content)

    def plot_entropy(self):
        """Plot Shannon entropy per position."""
        from pylab import axvline, plot, ylabel, ylim

        axvline(0.5, c="k")
        plot(self.df.index, self.entropy)
        ylabel("Bits")
        ylim([0, 2.1])

    def plot_information_content(self):
        """Plot information content (IC) per position."""
        from pylab import axvline, plot, ylabel, ylim

        axvline(0.5, c="k")
        plot(self.df.index, self.information_content)
        ylabel("Bits")
        ylim([0, 2.1])


class Kozak:
    """

        k = KOSAC("fasta", "gff", feature="gene")

    Filter only to keep ATG since others seems to ncRNA

    - raw Kozak sequence names and counts
    - a Kozak is e.g GGCRGG  . first position is the less important

    - for the enumeration of kmers, get of the rid of the Ns

    - odds ratio have 4 cases depending on the on enumeration:
        use entire genome
        use chromosome by chromosome
        use of gene on genome
        use gene on chromosomes


    Table of counts of Kozak sequences without dna ambiguities.
    - across the entire genome
    - by chromosomes

        counts = k.get_all_kmer_counts()
        counts_chroms = k.get_all_kmer_counts_by_chromosome()
        counts_genes = k.get_all_kmer_counts_genes_only()

        # proportions of kmer in genes:
        sum(list(counts_genes.values())) / length_genome

        # counts in chroms should equal counts in genomes:
        Sgenes = sum([sum(list(counts_chroms[x].values())) for x in counts_chroms.keys()])
        Sgenome = sum(list(counts.values()))

    ::

        k = Kozak("ecoli_MG1655.fa", "ecoli_MG1655.gff", "gene", "ID")
        df = k.get_data()
        k.plot_logo(df.query("start_codon=='ATG'"))


    """

    def __init__(self, fasta, gff, genetic_type="gene", attribute="ID", light=True):
        self.fasta = FastA(fasta)
        if isinstance(gff, GFF3):
            self.gff = gff
        else:
            self.gff = GFF3(gff, light=light)
        logger.info("scanning GFF")
        self._valid_genetic_types = self.gff.features  # some overhead but required
        self.genetic_type = genetic_type
        self.attribute = attribute
        self.set_context()

        # internal boundaries to precompute data
        self.RIGHT = 100
        self.LEFT = 100

        # place holder for various metrics
        self.metrics = {}
        self._background_method = "context"  # default

    @property
    def genetic_type(self):
        return getattr(self, "_genetic_type", None)

    @genetic_type.setter
    def genetic_type(self, value):
        if value not in self._valid_genetic_types:
            raise ValueError(f"Invalid genetic_type '{value}'. " f"Must be one of {sorted(self._valid_genetic_types)}")
        self._genetic_type = value

    def set_context(
        self,
        left_kozak=6,
        right_kozak=6,
        keep_ATG_only=True,
        include_start_codon=False,
        background_method="context",
        collapse_first_cds=True,
    ):
        """Configure context windows and feature-row collapsing.

        :param int left_kozak: number of nucleotides to keep upstream of the
            start codon.
        :param int right_kozak: number of nucleotides to keep downstream of
            the start codon.
        :param bool keep_ATG_only: if True, restrict downstream analyses to
            rows whose start codon is ``ATG``.
        :param bool include_start_codon: include the start codon itself in
            the Kozak window when True.
        :param str background_method: one of ``"context"``, ``"shuffled"``,
            or ``"uniform"``.
        :param bool collapse_first_cds: when True (default), collapse
            multi-exon CDS rows to one row per transcript (the 5'-most CDS,
            which is the only CDS row corresponding to a real start codon).
            See :meth:`_collapse_to_first_cds` for rationale. Set to False
            to recover the legacy behaviour where every CDS row is treated
            as a separate start (useful for benchmarking the bug fix).
        """
        assert left_kozak > 0
        assert right_kozak > 0
        self._left_kozak = left_kozak
        self._right_kozak = right_kozak
        self._keep_ATG_only = keep_ATG_only
        self._include_start_codon = include_start_codon
        self._background_method = background_method
        self._collapse_first_cds = collapse_first_cds

        _valid_methods = ["context", "shuffled", "uniform"]
        if background_method not in _valid_methods:
            raise ValueError(f"background_method must be one of {_valid_methods}")
        self._get_full_background.cache_clear()
        self._get_annotated_starts_by_chrom.cache_clear()
        try:
            del self._cached_df
        except AttributeError:
            pass

    def _get_left_kozak(self):
        return self._left_kozak

    left_kozak = property(_get_left_kozak)

    def _get_right_kozak(self):
        return self._right_kozak

    right_kozak = property(_get_right_kozak)

    def _get_keep_ATG_only(self):
        return self._keep_ATG_only

    def _set_keep_ATG_only(self, value):
        assert value in [True, False]
        self._keep_ATG_only = value

    keep_ATG_only = property(_get_keep_ATG_only, _set_keep_ATG_only)

    def _get_include_start_codon(self):
        return self._include_start_codon

    def _set_include_start_codon(self, value):
        assert value in [True, False]
        self._include_start_codon = value

    include_start_codon = property(_get_include_start_codon, _set_include_start_codon)

    def _get_collapse_first_cds(self):
        return self._collapse_first_cds

    def _set_collapse_first_cds(self, value):
        assert value in [True, False]
        self._collapse_first_cds = value
        try:
            del self._cached_df
        except AttributeError:
            pass

    collapse_first_cds = property(_get_collapse_first_cds, _set_collapse_first_cds)

    def _collapse_to_first_cds(self, df):
        """Collapse multi-exon CDS rows to one row per transcript.

        For multi-exon eukaryotic genes, GFF3 files contain one CDS row per
        coding exon. Only the 5'-most CDS row corresponds to the actual
        start codon; subsequent rows correspond to internal exon boundaries
        (splice junctions), which are NOT biological start codons.
        Including them pollutes the Kozak signal with what amounts to
        random in-frame codons and dilutes the ATG fraction toward
        ``1 / mean_exons_per_gene``.

        This helper groups rows by their parent transcript (taken from the
        ``Parent`` GFF attribute) and keeps the 5'-most CDS per group:

        * for ``+`` strand features, the row with the smallest ``start``;
        * for ``-`` strand features, the row with the largest ``stop``.

        If the ``Parent`` attribute is missing for every row (typical when
        ``genetic_type='gene'`` since gene features have no parent, or for
        flat prokaryotic GFFs without parent/child relationships), the
        dataframe is returned unchanged. Rows whose ``Parent`` is missing
        while others have one are also kept unchanged.

        :param pandas.DataFrame df: GFF3 dataframe pre-filtered to the
            selected ``genetic_type`` rows.
        :return: dataframe with at most one CDS row per transcript.
        :rtype: pandas.DataFrame
        """
        if len(df) == 0 or "attributes" not in df.columns:
            return df

        def _get_parent(attrs):
            if isinstance(attrs, dict):
                p = attrs.get("Parent")
                if isinstance(p, list):
                    return p[0] if p else None
                return p
            return None

        parents = df["attributes"].map(_get_parent)
        if parents.isna().all():
            logger.info("No 'Parent' attribute found on selected features; " "skipping multi-exon CDS collapse.")
            return df

        df = df.assign(_parent=parents)
        with_parent = df[df["_parent"].notna()]
        without_parent = df[df["_parent"].isna()]

        plus = (
            with_parent[with_parent["strand"] == "+"]
            .sort_values("start", ascending=True)
            .drop_duplicates("_parent", keep="first")
        )
        minus = (
            with_parent[with_parent["strand"] == "-"]
            .sort_values("stop", ascending=False)
            .drop_duplicates("_parent", keep="first")
        )
        out = pd.concat([plus, minus, without_parent]).drop(columns="_parent")
        logger.info(
            f"Collapsed multi-exon {self.genetic_type} rows: " f"{len(df)} -> {len(out)} (one row per transcript)."
        )
        return out

    def _compute(self):
        # Store all kozak sequences with generous left and right values
        if hasattr(self, "_cached_df"):
            return self._cached_df

        genetic_type = self.genetic_type
        LEFT = self.LEFT
        RIGHT = self.RIGHT

        logger.info("1. Checking consistency between FastA and GFF files")
        if len(set(self.gff.contig_names).intersection(self.fasta.names)) != len(self.fasta.names):
            logger.warning("GFF and FASTA sequence identifiers have different length")

        # trick to not load all data in memory but only selected genetic type
        self.gff.skip_types = [x for x in self._valid_genetic_types if x != self.genetic_type]
        self.gff._df = None

        # Reading GFF
        gff = self.gff.df

        # For multi-exon transcripts, GFF lists one CDS row per exon but only
        # the 5'-most CDS corresponds to a real start codon. Collapse here
        # before iterating, otherwise internal exon boundaries are mistaken
        # for start codons and ATG_ratio drops to ~1/mean_exons_per_gene.
        if self._collapse_first_cds:
            gff = self._collapse_to_first_cds(gff)

        # we split by chrom to get the sequence one by one.
        data = []
        warning_message = None
        for chrom, subdf in gff.groupby("seqid"):
            sequence = self.fasta[chrom].upper()
            for row in subdf.itertuples(index=False):
                try:
                    ID = getattr(row, self.attribute)
                except AttributeError:
                    ID = ""
                    if warning_message is None:  # show message only once
                        warning_message = (
                            f"Attribute {self.attribute} not found in GFF file. Setting ID to empty string."
                        )
                        logger.warning(warning_message)

                strand = row.strand
                start = row.start
                end = row.stop

                if strand == "+":
                    codon = sequence[start - 1 : start - 1 + 3]
                    kozak_left = sequence[start - 1 - LEFT : start - 1]
                    kozak_right = sequence[start + 3 : start + 3 + RIGHT]
                    data.append([chrom, ID, strand, start, end, codon, kozak_left, kozak_right])
                elif strand == "-":
                    codon = sequence[end - 3 : end]
                    codon = reverse_complement(codon)
                    koz = sequence[end : end + LEFT]
                    kozak_left = reverse_complement(koz)
                    koz = sequence[end - RIGHT - 3 : end - 3]
                    kozak_right = reverse_complement(koz)
                    data.append([chrom, ID, strand, start, end, codon, kozak_left, kozak_right])

        _cols = ["chrom", "ID", "strand", "start", "end", "start_codon", "kozak_left", "kozak_right"]
        df = pd.DataFrame(data, columns=_cols) if data else pd.DataFrame(columns=_cols)

        # Cleanup region where selected left and right sub sequence do not
        # have the correct length. This could be a gene right at the border of contig.
        N0 = len(df)
        if N0 == 0:
            logger.warning("No data found after parsing the GFF/FASTA files.")
            self.metrics["feature"] = self.genetic_type
            self.metrics["ATG_ratio"] = 0.0
            self._cached_df = df
            return df

        df = df[[len(x) == LEFT for x in df["kozak_left"].values]]
        df = df[[len(x) == RIGHT for x in df["kozak_right"].values]]
        ratio = len(df) / N0
        self.metrics["feature"] = self.genetic_type
        logger.info(
            f"Filtered to keep only rows with correct left and right Kozak lengths: {len(df)} / {N0} ({ratio:.1%})"
        )
        self.metrics["ATG_ratio"] = sum(df["start_codon"] == "ATG") / len(df) if len(df) > 0 else 0.0
        self._cached_df = df
        return df

    def get_data(self):
        # do not touch the original
        df = self._compute().copy()

        N = len(df)
        if N == 0:
            self.atg_contribution = 0.0
            return df

        if self.keep_ATG_only:
            df = df[df["start_codon"] == "ATG"]
            logger.info(f"Filtered to keep only ATG start codons: {len(df)} / {N} ({len(df)/N:.1%})")
            self.atg_contribution = 100.0 * len(df) / N
        else:
            n = len(df[df["start_codon"] == "ATG"])
            self.atg_contribution = 100 * float(n) / N
        df["kozak_left"] = [x[-self._left_kozak :] for x in df["kozak_left"]]
        df["kozak_right"] = [x[: self._right_kozak] for x in df["kozak_right"]]
        df["sequence"] = [x + y + z for x, y, z in zip(df["kozak_left"], df["start_codon"], df["kozak_right"])]
        freqs = Counter(df["sequence"])
        df["frequency"] = [freqs[x] / len(df) for x in df["sequence"].values]

        return df

    def filter_dataframe(self, df, strand=None, query=None, genes_set=None, attribute=None):
        if strand == "strand+":
            df = df.query("strand == '+'")
        elif strand == "strand-":
            df = df.query("strand == '-'")

        if query:
            df = df.query(query)

        if genes_set is not None:
            df = df[df[attribute].isin(genes_set)]

        return df

    @functools.lru_cache(maxsize=1)  # None means no limit on cache size
    def get_gc_per_chromosome(self, quiet=True):
        GCs = []
        chrom_names = []
        # we use a list rather than a dictionary to keep same order as in the
        # fasta
        from sequana.tools import fast_gc_content

        for chrom in tqdm(self.fasta.names, disable=quiet):
            sequence = self.fasta.sequences[self.fasta.names.index(chrom)]
            GCs.append(fast_gc_content(sequence))
            chrom_names.append(chrom)
        del sequence
        return chrom_names, GCs

    def plot_GC_per_chromosome(self, ylim=[0, 100]):
        chrom_names, GCs = self.get_gc_per_chromosome()
        pylab.plot(range(0, len(GCs)), [100 * x for x in GCs], "o-")
        if len(GCs) < 100:
            pylab.xticks(range(0, len(GCs)), chrom_names, rotation=90)
        pylab.xlabel("Chromosome")
        pylab.grid()
        pylab.ylim(ylim)
        pylab.ylabel("GC (%)")

    def _get_logo_data(
        self,
        df=None,
    ):

        if df is None:
            df = self.get_data()

        try:
            Nl = len(df["kozak_left"].iloc[0])
            Nr = len(df["kozak_right"].iloc[0])
        except KeyError:
            N = self.left_kozak
            return {"status": "Warning", "msg": "No data found"}

        if self.include_start_codon:
            pos = [defaultdict(int) for x in range(Nl + 3 + Nr)]
            for seql, start, seqr in zip(df.kozak_left.values, df.start_codon.values, df.kozak_right.values):
                # left sequence
                for i in range(Nl):
                    pos[i][seql[i]] += 1
                # start codon
                for i in range(3):
                    pos[Nl + i][start[i]] += 1

                # right sequence
                for i in range(Nr):
                    pos[Nl + 3 + i][seqr[i]] += 1
        else:
            pos = [defaultdict(int) for x in range(Nl + Nr)]
            for seql, seqr in zip(df.kozak_left.values, df.kozak_right.values):
                for i in range(Nl):
                    pos[i][seql[i]] += 1

                for i in range(Nr):
                    pos[Nl + i][seqr[i]] += 1

        logo_data = pd.DataFrame(pos).fillna(0)
        if "N" in logo_data.columns:
            del logo_data["N"]

        # ignore all non ACGT bases
        logo_data = logo_data[["A", "C", "G", "T"]]

        logo_data = logo_data.divide(logo_data.sum(axis=1), axis=0)

        L, R = self.left_kozak, self.right_kozak
        if self.include_start_codon:
            logo_data.index = list(range(-L, 0)) + list(range(1, R + 4))
        else:
            logo_data.index = list(range(-L, 0)) + list(range(4, R + 4))

        return logo_data

    _DNA_COLOR_SCHEMES = {
        "colorblind": {"A": "#0072B2", "C": "#E69F00", "G": "#009E73", "T": "#D55E00"},
        "classic": {"A": "green", "C": "blue", "T": "red", "G": "orange"},
    }

    # Maximum information content per position for a 4-letter (DNA) alphabet (log2(4) = 2 bits)
    _MAX_BITS = 2

    def _get_color_scheme(self, color_scheme):
        assert color_scheme in self._DNA_COLOR_SCHEMES, "color_scheme must be colorblind or classic"
        return self._DNA_COLOR_SCHEMES[color_scheme]

    def plot_logo(self, df=None, ax=None, color_scheme="colorblind"):
        logo_data = self._get_logo_data(df)
        self._plot_logo(logo_data, ax=ax, color_scheme=self._get_color_scheme(color_scheme))

        return logo_data

    def plot_logo_bits(self, df=None, ax=None, color_scheme="colorblind"):
        """Plot sequence logo with letter heights scaled by information content (bits).

        Unlike :meth:`plot_logo` which shows relative nucleotide frequencies with
        uniform column heights, this method scales each column by its information
        content (IC = 2 - Shannon entropy), so highly conserved positions appear
        tall (up to 2 bits) and variable positions appear short.

        :param df: output of :meth:`get_data`. If None, :meth:`get_data` is called.
        :param ax: matplotlib axes object. If None, the current axes is used.
        :param color_scheme: ``"colorblind"`` (default) or ``"classic"``.
        :return: the probability logo_data DataFrame (same as :meth:`plot_logo`).
        """
        import logomaker

        logo_data = self._get_logo_data(df)

        data = logo_data.copy()
        indices = data.index
        data.reset_index(inplace=True, drop=True)

        # Convert from probability matrix to information-content (bits) matrix
        info_data = logomaker.transform_matrix(data, from_type="probability", to_type="information")

        logo = logomaker.Logo(info_data, ax=ax, color_scheme=self._get_color_scheme(color_scheme), stack_order="fixed")
        pylab.axvline(self._left_kozak - 0.5, color="k", lw=2)
        if self.include_start_codon:
            pylab.axvline(self._left_kozak + 2.5, color="k", lw=2)
        pylab.xticks(range(len(data)), indices)
        pylab.ylabel("Bits")
        pylab.ylim([0, self._MAX_BITS])

        return logo_data

    def _add_purine_pyrimidine(self, df):
        df["R"] = df["A"] + df["G"]
        df["Y"] = df["C"] + df["T"]
        return df

    def plot_logo_purine_pyrimidine(self, df=None, ax=None):
        """
        df is the output of :meth:`get_data`
        """
        logo_data = self._get_logo_data(df)
        logo_data = self._add_purine_pyrimidine(logo_data)
        self._plot_logo(logo_data[["R", "Y"]], ax=ax, color_scheme={"Y": "purple", "R": "#78bc00"})
        return logo_data

    def get_information_content(self, motif):
        return Motif(motif).information_content

    @contextmanager
    def temporary_lr(self, left=None, right=None):
        old_left = self.left_kozak
        old_right = self.right_kozak
        try:
            if left is not None:
                self._left_kozak = left
            if right is not None:
                self._right_kozak = right
            yield
        finally:
            self._left_kozak = old_left
            self._right_kozak = old_right

    def _get_KL_data(self, df=None, n_boot=500, ci=95, left=None, right=None):

        if left or right:
            self._get_full_background.cache_clear()
            with self.temporary_lr(left=left, right=right):
                logo_data = self._get_logo_data(df=df)
                background = self.get_random_contexts()
                mu, low, high = self.bootstrap(
                    logo_data,
                    background,
                    n_boot=n_boot,
                    ci=ci,
                )
        else:
            logo_data = self._get_logo_data(df=df)
            background = self.get_random_contexts()
            mu, low, high = self.bootstrap(
                logo_data,
                background,
                n_boot=n_boot,
                ci=ci,
            )

        df = pd.DataFrame({"position": logo_data.index, "KL_divergence": mu, "ci_low": low, "ci_high": high})
        return df

    def plot_KL_divergence(self, df=None, ax=None, n_boot=500, ci=95):
        KL = self._get_KL_data(df, n_boot, ci)

        if ax is None:
            fig, ax = pylab.subplots(figsize=(12, 4))
        ax.plot(KL["position"], KL["KL_divergence"], marker="o")
        ax.fill_between(KL["position"], KL["ci_low"], KL["ci_high"], color="gray", alpha=0.3)
        ax.set_xlabel("Position relative to start codon")
        ax.set_ylabel("KL divergence (bits)")
        ax.set_title("Positional KL divergence vs background distribution")
        return KL

    def _plot_logo(self, data, color_scheme=None, ax=None):
        import logomaker

        data = data.copy()
        indices = data.index
        data.reset_index(inplace=True, drop=True)

        logo = logomaker.Logo(data, ax=ax, color_scheme=color_scheme, stack_order="fixed")
        for x in [0.25, 0.5, 0.75]:
            pylab.axhline(x, color="grey", zorder=-1, ls="--")
        pylab.axvline(self._left_kozak - 0.5, color="k", lw=2)

        if self.include_start_codon:
            pylab.axvline(self._left_kozak + 2.5, color="k", lw=2)
        pylab.xticks(range(len(data)), indices)
        _ = pylab.yticks([])

    def plot_kozak_chi2(self, motif=None, GC_mode=None, fontsize=10, noplot=False):
        """
        for GC computation, 2 options:
            1. genomic ATG background excluding annotated starts. Good for GC bias, codon bias, local sequence structure.
            2. genome wide base composition. simple and fast but inflates signal in GC biases genomes. does not control for ATG specific context. This is a simple a genome-wide background is estimated assuming strand symmetry:
                G = C = GC / 2
                A = T = AT / 2
        """
        if motif is None:
            motif = self._get_logo_data()
        df = motif  # just an alias

        if GC_mode is None:
            if self._background_method == "context":
                GC_mode = "context"
            else:
                # shuffled and uniform
                GC_mode = "genome"

        if GC_mode == "context":
            atg = self.get_random_contexts()
            seq = atg["kozak_left"].str.cat(atg["kozak_right"])
            GC = seq.str.count("[GC]").sum() / seq.str.len().sum()
        elif GC_mode in ["shuffled", "genome", "uniform"]:
            if "GC" in self.metrics:
                GC = self.metrics["GC"]
            else:
                GC = self.metrics.get("GC", self.fasta.GC_content()) / 100
                self.metrics["GC"] = GC
        else:
            raise ValueError("GC_mode must be context or genome")
        assert GC <= 1 and GC >= 0
        # computes chi2
        from scipy.stats import chisquare

        expected = {"A": 100 * (1 - GC) / 2.0, "C": 100 * GC / 2.0, "G": 100 * GC / 2.0, "T": 100 * (1 - GC) / 2}
        observed = 100 * df[["A", "C", "G", "T"]].copy()

        chi2 = []
        for _, row in observed.iterrows():
            chi2_stat, p_value = chisquare(
                [row["A"], row["C"], row["G"], row["T"]], f_exp=[expected[n] for n in ["A", "C", "G", "T"]]
            )
            chi2.append(p_value)
        observed["chi2"] = chi2
        if noplot:
            return observed

        df = observed
        fig, ax = pylab.subplots(figsize=(12, 10))
        # Plotting values
        ax.matshow(df[["A", "C", "G", "T"]].T, cmap="viridis", vmin=0, vmax=50)

        # Adding text annotations
        df = df[["A", "C", "G", "T"]].T
        for i in range(df.shape[0]):
            for j in range(df.shape[1]):
                ax.text(j, i, f"{df.values[i, j]:.1f}", ha="center", va="center", color="red", fontsize=fontsize)
        pylab.yticks([0, 1, 2, 3], ["A", "C", "G", "T"])
        stars = []
        for p in observed["chi2"]:
            if p < 0.01:
                stars.append("**")
            elif p < 0.05:
                stars.append("*")
            else:
                stars.append(" ")

        pylab.xticks(range(len(observed)), stars, fontsize=fontsize)
        pylab.axvline(self.left_kozak - 0.5, lw=2, color="k")
        if self.include_start_codon:
            pylab.axvline(self._left_kozak + 2.5, color="k", lw=2)

        return observed

    @functools.lru_cache(maxsize=1)
    def _get_annotated_starts_by_chrom(self):
        """
        Returns:
            dict[str, set[int]]: chromosome -> set of start positions (1-based)
        """

        self.gff.skip_types = [x for x in self.gff.features if x != self.genetic_type]
        df = self.gff.df

        starts = defaultdict(set)
        for chrom, start in zip(df.seqid.values, df.start.values):
            starts[chrom].add(int(start))

        return starts

    def _is_annotated_start(self, chromosome, position):
        return position in self._get_annotated_starts_by_chrom().get(chromosome, set())

    def get_random_contexts(self, Nmax=None, quiet=True):
        """
        Return a background distribution of Kozak contexts.
        """
        # We need Nmax at the end, but the core generation can be cached for a large number
        # unless Nmax is very small.
        # However, it's safer to just cache the full background and then sample if needed.
        bg = self._get_full_background(quiet=quiet)

        df = self.get_data()
        if Nmax is None:
            Nmax = len(df)

        if len(bg) == Nmax:
            return bg.copy()
        elif len(bg) > Nmax:
            return bg.iloc[:Nmax].copy()
        else:
            return bg.sample(n=Nmax, replace=True).copy()

    @functools.lru_cache(maxsize=None)
    def _get_full_background(self, quiet=True):
        df = self.get_data()
        Nmax = len(df)

        if self._background_method == "shuffled":
            return self._get_shuffled_background(df, Nmax)
        elif self._background_method == "uniform":
            return self._get_uniform_background(df, Nmax)
        elif self._background_method == "context":
            # context uses genomic sampling for background generation
            return self._get_genomic_background(df, Nmax, quiet)
        else:
            raise ValueError(f"Unsupported background method: {self._background_method}")

    def _get_uniform_background(self, df, Nmax):
        """
        Generate background by independently sampling nucleotides at each position
        based on overall genomic GC content.
        """
        if "GC" in self.metrics:
            GC = self.metrics["GC"]
        else:
            GC = self.fasta.GC_content() / 100
            self.metrics["GC"] = GC

        # Probabilities: [A, C, G, T]
        p = [(1 - GC) / 2, GC / 2, GC / 2, (1 - GC) / 2]
        bases = ["A", "C", "G", "T"]

        L = self.left_kozak
        R = self.right_kozak

        # Sample proportions for start codons from observed data
        codon_counts = df["start_codon"].value_counts(normalize=True)
        target_codons = codon_counts.index.tolist()
        proportions = codon_counts.values

        lefts = []
        rights = []
        starts = []

        for _ in range(Nmax):
            left = "".join(np.random.choice(bases, size=L, p=p))
            right = "".join(np.random.choice(bases, size=R, p=p))
            start = np.random.choice(target_codons, p=proportions)
            lefts.append(left)
            rights.append(right)
            starts.append(start)

        return pd.DataFrame({"kozak_left": lefts, "start_codon": starts, "kozak_right": rights})

    def _get_shuffled_background(self, df, Nmax):
        """
        Generate background by shuffling observed Kozak sequences.
        """
        # We sample with replacement if Nmax > len(df)
        if Nmax > len(df):
            df_sampled = df.sample(n=Nmax, replace=True)
        else:
            df_sampled = df.sample(n=Nmax)

        lefts = df_sampled["kozak_left"].tolist()
        rights = df_sampled["kozak_right"].tolist()
        starts = df_sampled["start_codon"].tolist()

        # Shuffle each sequence independently
        new_lefts = []
        new_rights = []
        for l, r in zip(lefts, rights):
            l_list = list(l)
            r_list = list(r)
            random.shuffle(l_list)
            random.shuffle(r_list)
            new_lefts.append("".join(l_list))
            new_rights.append("".join(r_list))
        lefts, rights = new_lefts, new_rights

        return pd.DataFrame({"kozak_left": lefts, "start_codon": starts, "kozak_right": rights})

    def _get_genomic_background(self, df, Nmax, quiet=True):
        """
        Generate background by sampling random genomic locations with matching start codons.
        """
        # Calculate start codon proportions from observed data
        codon_counts = df["start_codon"].value_counts(normalize=True)
        target_codons = codon_counts.index.tolist()
        proportions = codon_counts.values

        chrom_lengths = [l for l in self.fasta.lengths]
        cum_lengths = np.cumsum(chrom_lengths)
        total_length = cum_lengths[-1]

        lefts = []
        rights = []
        starts = []

        annotated = self._get_annotated_starts_by_chrom()

        while len(lefts) < Nmax:
            # Random global position
            gpos = random.randint(0, total_length - 4)

            # Map global position to chromosome
            chrom_idx = np.searchsorted(cum_lengths, gpos, side="right")
            local_pos = gpos if chrom_idx == 0 else gpos - cum_lengths[chrom_idx - 1]

            try:
                seq = self.fasta.sequences[chrom_idx].upper()
            except IndexError:
                continue

            chrom = self.fasta.names[chrom_idx]

            # Sample a target start codon based on observed proportions
            target = np.random.choice(target_codons, p=proportions)

            # Find next occurrence of this target codon
            local_pos = seq.find(target, local_pos)
            if local_pos == -1:
                continue

            # Context check (not an annotated start)
            if local_pos + 1 in annotated.get(chrom, set()):
                continue

            left = seq[local_pos - self.left_kozak : local_pos]
            right = seq[local_pos + 3 : local_pos + 3 + self.right_kozak]

            if len(left) != self.left_kozak or len(right) != self.right_kozak:
                continue

            lefts.append(left)
            rights.append(right)
            starts.append(target)

        return pd.DataFrame({"kozak_left": lefts, "start_codon": starts, "kozak_right": rights})

    def export_meme(self, filename, name="Kozak"):
        """PWM compatible with standard motif scanners"""
        df = self.get_data()
        pwm = self._get_logo_data(df)
        with open(filename, "w") as f:
            f.write("MEME version 4\n\n")
            f.write("ALPHABET= ACGT\n\n")
            f.write("strands: + -\n\n")
            f.write(f"MOTIF {name}\n")
            f.write(f"letter-probability matrix: alength= 4 w= {len(pwm)}\n")
            for _, row in pwm.iterrows():
                f.write(" ".join(f"{row[x]:.6f}" for x in "ACGT") + "\n")

    def kl_vs_random_atg(self, motif_df, random_df):
        """
        Compute position-wise divergence using Kullback–Leibler (KL) divergence
        between Kozak contexts and random (non-annotated) ATG contexts.

        This quantifies how specific the Kozak signal is compared
        to generic ATG neighborhoods.

        Compute positional KL divergence between the observed Kozak motif
        and a background nucleotide distribution.

        This method computes, for each position i:

            D_KL(P_i || Q) = sum_i (P_i x log2(P_i / Q_i))

        where P_i is the observed nucleotide frequency distribution
        at position i, and Q is a fixed background distribution.

        Notes
        -----
        - This quantity is mathematically related to Shannon entropy.
        - When Q is uniform, this is equivalent to classical
          sequence logo information content.
        - The result is deterministic for a given motif_df and  random_df input pair

        Parameters
        ----------
        motif_df : pandas.DataFrame
            DataFrame with columns ['A', 'C', 'G', 'T'] containing
            nucleotide frequencies per position.

        random_df : pandas.DataFrame
            Background distribution.

        Returns
        -------
        numpy.ndarray
            KL divergence (bits) for each motif position.
        """
        p = motif_df[["A", "C", "G", "T"]].values.astype(float) + 1e-12
        q = random_df[["A", "C", "G", "T"]].values.astype(float) + 1e-12

        # KL divergence: sum(p * log2(p/q))
        return np.sum(p * np.log2(p / q), axis=1)

    def bootstrap(self, df, contexts, n_boot=500, ci=95, sample_size=200):
        """
        Perform bootstrapping to compute confidence intervals for KL divergence.
        """
        n = len(contexts)
        if self.include_start_codon:
            seqs = contexts["kozak_left"].str.cat(contexts["start_codon"]).str.cat(contexts["kozak_right"]).values
            L = self.left_kozak + 3 + self.right_kozak
        else:
            seqs = contexts["kozak_left"].str.cat(contexts["kozak_right"]).values
            L = self.left_kozak + self.right_kozak

        # Fast one-hot encoding of all sequences
        # Mapping: A=0, C=1, G=2, T=3
        mapping = {"A": 0, "C": 1, "G": 2, "T": 3}

        # Convert all sequences to a numeric (N, L) array
        # This is much faster than strings in a loop
        numeric_seqs = np.zeros((n, L), dtype=int)
        for i, s in enumerate(seqs):
            # Skip non-ACGT characters if any (default to 0/A or we could handle Ns)
            numeric_seqs[i] = [mapping.get(c, 0) for c in s]

        # Pre-allocate boot array
        boot = np.zeros((n_boot, L))

        # Pre-compute P values for KL divergence calculation
        # p is (L, 4)
        p = df[["A", "C", "G", "T"]].values.astype(float) + 1e-12

        sample_size = min(n, sample_size)

        for b in tqdm(range(n_boot)):
            # Randomly sample indices
            idx = np.random.randint(0, n, sample_size)
            sample_numeric = numeric_seqs[idx]

            # Compute counts efficiently using numpy
            # For each position, count occurrences of 0, 1, 2, 3
            # Resulting array: (L, 4)
            counts = np.zeros((L, 4))
            for i in range(4):
                counts[:, i] = np.sum(sample_numeric == i, axis=0)

            # Convert counts to background probabilities q
            q = counts / sample_size + 1e-12

            # Compute KL divergence: sum(p * log2(p/q)) along nucleotides (axis 1)
            # kl shape: (L,)
            kl = np.sum(p * np.log2(p / q), axis=1)
            boot[b] = kl

        mean = boot.mean(axis=0)
        low = np.percentile(boot, (100 - ci) / 2, axis=0)
        high = np.percentile(boot, 100 - (100 - ci) / 2, axis=0)

        return mean, low, high


class KLAnalysis:
    def __init__(self, df):
        self.KL = df

    def get_KSI(self):
        """average across Kozak length (6bp before ATG)

        Average across 6 bp (Kozak sequence)
        """
        return sum(self.KL.query("position<0 and position>=-6")["KL_divergence"]) / 6

    def get_total_information(self, min=-1e6, max=0):
        """Area under the curve for position<0

        This removes dilution from averaging.
        Independent of window scaling. Measures total constraint.
        """
        return self.KL.query("position<@max and position>=@min")["KL_divergence"].sum()

    def get_peak_strength(self):
        """max peak strength"""
        return self.KL["KL_divergence"].max()

    def get_peak_position(self):
        """max peak position"""
        index = self.KL["KL_divergence"].argmax()
        if "position" in self.KL.columns:
            return self.KL["position"].iloc[index]
        else:
            return self.KL.index[index]

    def get_signal_concentration(self):
        # peak concentration measures how concentrated is the signal. sharp pike == high C (eukar), broader = lower C (proka)
        Ipeak = self.get_peak_strength()
        Itot = self.get_total_information()
        signal_concentration = Ipeak / Itot
        return signal_concentration

    def get_III(self):
        # Alternative Single-Value Metric that balances global and local structure:
        from math import sqrt

        Ipeak = self.get_peak_strength()
        Itot = self.get_total_information()
        return sqrt(Ipeak * Itot)

    def get_power(self):
        # power metric: emphasize sharpness
        Ipeak = self.get_peak_strength()
        Itot = self.get_total_information()
        P = Ipeak * Ipeak / Itot
        return P

    def compute_W50(self):
        df_up = self.KL.query("position < 0")
        KL_values = df_up["KL_divergence"].values

        I_total = KL_values.sum()
        if I_total == 0:
            return 0

        sorted_KL = np.sort(KL_values)[::-1]
        cumulative = np.cumsum(sorted_KL)

        W50 = np.searchsorted(cumulative, 0.5 * I_total) + 1
        return W50


class ConsensusBuilder:
    def __init__(self, df):

        self.df = df.copy()

    def _iupac(self, bases, notation="iupac"):
        """Encode a set of bases.

        :param bases: a list (e.g. ``['A', 'G']``) or a string of bases.
        :param notation: ``iupac`` returns the single-character ambiguity code
            (e.g. ``R`` for A/G); ``expanded`` returns the bases joined by a
            slash (e.g. ``A/G``), following the original Cavener notation.
        """
        # Normalise to a sorted list so the result is order-independent.
        bases = sorted(bases)

        if notation == "expanded":
            return "/".join(bases)
        elif notation == "iupac":
            from sequana.iuapc import dna_ambiguities_r

            # The lookup dict keys are also normalised because some are not
            # alphabetically sorted (e.g. 'GC' instead of 'CG').
            lookup = {"".join(sorted(k)): v for k, v in dna_ambiguities_r.items()}
            return lookup.get("".join(bases), "N")
        else:
            raise ValueError(f"Unknown notation: {notation}. Use 'iupac' or 'expanded'.")

    def _information_content(self, row):
        """Shannon information (DNA max = 2 bits)."""
        p = row.values
        p = p[p > 0]
        H = -np.sum(p * np.log2(p))
        return 2 - H

    def get_consensus(
        self,
        mode="majority",
        threshold=0.25,
        relative=0.8,
        strong=0.6,
        majority=0.5,
        info_threshold=1.0,
        notation="iupac",
    ):
        """
        Parameters
        ----------
        mode : str
            majority | threshold | relative | information | max_only
        threshold : float
            used in threshold mode
        relative : float
            keep bases >= relative * max_frequency
        strong : float
            uppercase if max >= strong
        info_threshold : float
            uppercase if information >= this value
        notation : str
            how to encode ambiguous positions: ``iupac`` (single ambiguity
            code, e.g. ``R``) or ``expanded`` (slash-joined bases, e.g. ``A/G``)
        """

        consensus = []

        for _, row in self.df.iterrows():

            if mode == "max_only":
                base = row.idxmax()
                consensus.append(base)
                continue

            elif mode == "majority":
                if row.max() > majority:
                    base = row.idxmax()
                    consensus.append(base.upper())
                else:
                    consensus.append("n")

            elif mode == "threshold":
                bases = row[row >= threshold].index.tolist()
                letter = self._iupac(bases, notation=notation)

                if row.max() >= strong:
                    consensus.append(letter.upper())
                else:
                    consensus.append(letter.lower())

            elif mode == "relative":
                m = row.max()
                bases = row[row >= relative * m].index.tolist()
                letter = self._iupac(bases, notation=notation)

                if m >= strong:
                    consensus.append(letter.upper())
                else:
                    consensus.append(letter.lower())

            elif mode == "information":
                info = self._information_content(row)
                m = row.max()

                # select dominant bases (relative strategy)
                bases = row[row >= 0.8 * m].index.tolist()
                letter = self._iupac(bases, notation=notation)

                if info >= info_threshold:
                    consensus.append(letter.upper())
                else:
                    consensus.append(letter.lower())

            else:
                raise ValueError(f"Unknown mode: {mode}")

        return "".join(consensus)

    def get_consensus_cavener(self, notation="iupac"):
        """Consensus following the Cavener (1987) rule.

        For each position:

        - if a single nucleotide exceeds 50% frequency **and** is more than
          twice as frequent as the next one, it is reported as a capital
          letter;
        - otherwise, if the two most frequent nucleotides have a combined
          frequency above 75%, a co-consensus is reported as the uppercase
          encoding of that pair;
        - otherwise the most frequent nucleotide is reported as a lowercase
          letter.

        :param notation: how to encode the co-consensus pair: ``iupac`` (single
            ambiguity code, e.g. ``R``) or ``expanded`` (slash-joined bases,
            e.g. ``A/G``, as in the original Cavener notation).

        ::

            >>> import pandas as pd
            >>> from sequana.kozak import ConsensusBuilder
            >>> df = pd.DataFrame({
            ...     "A": [0.90, 0.45, 0.30, 0.55, 0.10],
            ...     "C": [0.04, 0.40, 0.30, 0.20, 0.10],
            ...     "G": [0.03, 0.10, 0.30, 0.15, 0.70],
            ...     "T": [0.03, 0.05, 0.10, 0.10, 0.10],
            ... })
            >>> cb = ConsensusBuilder(df)
            >>> cb.get_consensus_cavener()
            'AMaAG'
            >>> cb.get_consensus_cavener(notation="expanded")
            'AA/CaAG'

        Position 1 (A=0.90) is a capital because A exceeds 50% and is more than
        twice the next base. Position 2 (A=0.45, C=0.40) is a co-consensus
        (combined 0.85 > 75%). Position 3 (A=C=G=0.30) is lowercase. Position 5
        (G=0.70) is a capital again.

        Reference: Cavener DR (1987), Nucleic Acids Research 15(4):1353-1361.
        """
        consensus = []

        for _, row in self.df.iterrows():
            ranked = list(row.sort_values(ascending=False).items())
            top_base, top = ranked[0]
            second_base, second = ranked[1]

            if top > 0.5 and top > 2 * second:
                consensus.append(top_base.upper())
            elif (top + second) > 0.75:
                consensus.append(self._iupac([top_base, second_base], notation=notation).upper())
            else:
                consensus.append(top_base.lower())

        return "".join(consensus)

    def all_consensus(
        self, threshold=0.25, relative=0.8, strong=0.6, majority=0.5, info_threshold=1.0, notation="iupac"
    ):

        for mode in ["majority", "threshold", "relative", "information", "max_only"]:
            res = self.get_consensus(
                mode=mode,
                threshold=threshold,
                relative=relative,
                strong=strong,
                majority=majority,
                info_threshold=info_threshold,
                notation=notation,
            )
            print(f"{mode}: {res}")

        print(f"cavener: {self.get_consensus_cavener(notation=notation)}")


class KozakAddon(Kozak):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @functools.lru_cache(maxsize=None)  # None means no limit on cache size
    def builddata(self):
        dm, dp = self.get_kmer_counts()

        df = pd.DataFrame({"plus": dp, "minus": dm})
        df = df.sum(axis=1)
        df = df.fillna(0).sort_values(ascending=False)

        self._df = df
        return df

    def _get_df(self):
        if self._df is None:
            self.builddata()
        return self._df

    df = property(_get_df)

    def find_kmers(self, sequence, pattern=r"GCC[AG]CC"):
        """
        >>> k.find_kmers("GCCACC")
        True
        >>> k.find_kmers("AAAAAA")
        False
        """

        kozak_pattern = re.compile(pattern)
        match = kozak_pattern.search(sequence)
        if match:
            return True
        else:
            return False

    @functools.lru_cache(maxsize=None)  # None means no limit on cache size
    def get_all_kmer_counts(self, k=6, reverse=False, quiet=True):
        """Get all kmers from the entire genome"""
        counts = defaultdict(int)
        for sequence in tqdm(self.fasta.sequences, disable=quiet):
            count = self._get_kmer_from_sequence(sequence, k=k, reverse=reverse)

            for x, y in count.items():
                counts[x] += y
        return counts

    def _get_kmer_from_sequence(self, sequence, k=6, reverse=False):
        if reverse:
            sequence = reverse_complement(sequence)
        counts = defaultdict(int)
        for i in range(0, len(sequence) - k + 1):
            seq = sequence[i : i + k]
            if "N" not in seq:
                counts[sequence[i : i + k]] += 1
        return counts

    @functools.lru_cache(maxsize=None)  # None means no limit on cache size
    def get_all_kmer_counts_genes_only(self, k=6, genetic_type="gene", reverse=False):
        counts = defaultdict(int)

        df = self.gff.df.query("genetic_type==@genetic_type and strand=='+'")

        # here we loop through all genes
        sequences = self.fasta.sequences
        names = self.fasta.names

        for index, row in df.iterrows():
            # we find the sequence corresponding to the seqid
            try:
                index_seq = names.index(row.seqid)
                # and the given gene
                sequence = sequences[index_seq][row.start - 1 : row.stop]
                # to extract all kmers inside the gene
                count = self._get_kmer_from_sequence(sequence, k=k, reverse=reverse)

                for x, y in count.items():
                    counts[x] += y
            except ValueError:
                pass
        return counts

    @functools.lru_cache(maxsize=None)  # None means no limit on cache size
    def get_all_kmer_counts_by_chromosome(self, k=6, reverse=False):
        counts = {}
        for chrom, sequence in tqdm(zip(self.fasta.names, self.fasta.sequences)):
            count = self._get_kmer_from_sequence(sequence, k=k, reverse=reverse)
            counts[chrom] = count
        return counts

    """def get_all_kmers_normalised(self, k=6):
        kmers = self.get_all_kmer_counts_by_chromosome(k=k)
        kmers = pd.DataFrame(kmers)
        # kmers = kmers.divide(self.fasta.get_lengths_as_dict())*1000000
        return kmers
    """

    @functools.lru_cache(maxsize=None)  # None means no limit on cache size
    def get_odd_ratio(self, mode="all"):
        odds = []

        counts = self.df
        Ngenes = self.df.sum()
        N_kmers = len(self.df)
        print(f"Number of genes : {Ngenes}")
        print(f"Number of unique kmers : {N_kmers}")

        if mode == "all":
            kmers = self.get_all_kmer_counts()
            genome_size = sum(list(kmers.values()))
        elif mode == "gene":
            kmers = self.get_all_kmer_counts_genes_only()
            genome_size = sum(list(kmers.values()))
        elif mode == "chromosome":
            pass

        odds = []
        for kmer in counts.index:
            count = counts[kmer]

            A = count / Ngenes
            B = kmers[kmer] / genome_size
            odds.append([kmer, count, kmers[kmer], A, B, A / B])

        odds = pd.DataFrame(odds)
        odds.columns = ["kmer", "count", "count_genome", "ratio1", "ratio2", "odds"]
        return odds

    def plot_scatter_odds_ratio_gene_vs_genome(self):
        # odds ratio are computed on the same set of kmer but compared
        # to different random kmer distribution
        odds = self.get_odd_ratio(mode="all")
        odds_gene = self.get_odd_ratio(mode="gene")

        # so we can plot a scatter plot
        # odds['count;[ and odds_gene['count'] are the same
        pylab.scatter(odds["odds"], odds_gene["odds"], c=odds["count"], alpha=0.5, s=4 * odds["count"])
        pylab.xlabel("odds ratio (vs all entire genome)")
        pylab.ylabel("odds ratio (vs all entire coding genes)")
        pylab.colorbar()

    def plot_scatter_odds_ratio_annotated(self, pattern=r"GCC[AG]CC"):
        odds = self.get_odd_ratio(mode="gene")
        found = [10 if self.find_kmers(x, pattern=pattern) else 0 for x in odds["kmer"]]
        odds["found"] = found

        print(odds.query("found==10"))
        pylab.clf()
        pylab.scatter(odds["count"], odds["odds"], c=odds["found"], alpha=0.5, s=4 * odds["odds"])
        pylab.xlabel("Count")
        pylab.ylabel("odds")
        pylab.colorbar()
        return odds

    def _get_kmer_counts_strand(self, strand):
        df = self.get_data()
        chroms = df.chrom.unique()
        kozak = defaultdict(int)

        for chrom in tqdm(chroms):
            for seq in df.query("chrom==@chrom and strand==@strand").kozak_left:
                kozak[seq] += 1

        return kozak

    def get_kmer_counts(self):
        """Get kmer counts for both strands. Returns (minus_strand_counts, plus_strand_counts)"""
        dm = self._get_kmer_counts_strand("-")
        dp = self._get_kmer_counts_strand("+")
        return dm, dp

    def plot_logo_all_kmers(self):
        counts = self.get_all_kmer_counts()

        pos = [defaultdict(int) for x in range(6)]
        for sequence, N in counts.items():
            for i, letter in enumerate(sequence):
                pos[i][letter] += N

        logo_data = pd.DataFrame(pos)
        dd = logo_data.divide(logo_data.sum(axis=1), axis=0)
        self._plot_logo(dd)
        return dd

    def plot_cumulated(self):
        cs = np.cumsum(self.df)
        mid = cs.max() / 2
        pylab.plot(cs.values)
        pylab.axhline(mid)
        pylab.xlabel("Number of unique 6-mers")
        pylab.ylabel("Number of genes")


class KozakWeightScore:
    """Compute Kozak Similarity Score using weight matrix approach.

    Weight-matrix based scoring for translation initiation prediction.
    Implements the algorithm from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0256411

    Flexible left/right flank parameters allow custom window sizes.

    Example::

        from sequana.kozak import KozakWeightScore
        kss = KozakWeightScore(left_flank=10, right_flank=10)
        score = kss.score("CGCCGCCACCATGGCGGCGGAGG")
    """

    def __init__(self, weight_matrix=None, left_flank=10, right_flank=10):
        """Initialize weight-based scorer.

        Parameters
        ----------
        weight_matrix : np.ndarray, optional
            Shape (n_positions, 5) where columns: A, T, G, C, N.
            Default: canonical ATG 23bp matrix.
        left_flank : int, default 10
            Bases upstream of codon.
        right_flank : int, default 10
            Bases downstream of codon.
        """
        self.left_flank = left_flank
        self.right_flank = right_flank
        self.codon_len = 3

        if weight_matrix is None:
            self.weight_matrix = self._default_atg_matrix()
        else:
            self.weight_matrix = np.asarray(weight_matrix, dtype=float)

        if self.weight_matrix.shape[1] != 5:
            raise ValueError("Weight matrix must have 5 columns (A, T, G, C, N)")

        self._validate_params()

    def _validate_params(self):
        """Check flanks match matrix size."""
        required = self.left_flank + self.codon_len + self.right_flank
        if self.weight_matrix.shape[0] < required:
            raise ValueError(
                f"Weight matrix rows={self.weight_matrix.shape[0]}, "
                f"but {required} required for left={self.left_flank}, right={self.right_flank}"
            )

    @staticmethod
    def _default_atg_matrix():
        """Default 23bp ATG matrix."""
        return np.array(
            [
                [0.04210526, 0.0, 0.03157895, 0.05263158, 0.0],
                [0.04210526, 0.05263158, 0.10526316, 0.0625, 0.0],
                [0.03157895, 0.04210526, 0.05263158, 0.07368421, 0.0],
                [0.03157895, 0.01052632, 0.04210526, 0.05263158, 0.0],
                [0.08421053, 0.07368421, 0.18947368, 0.10526316, 0.0],
                [0.04210526, 0.05263158, 0.05263158, 0.08421053, 0.0],
                [0.12631579, 0.0625, 0.12631579, 0.21052632, 0.0],
                [0.83157895, 0.12631579, 0.65263158, 0.16842105, 0.0],
                [0.15789474, 0.06315789, 0.11578947, 0.2, 0.0],
                [0.21052632, 0.09473684, 0.31578947, 0.51578947, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [0.24210526, 0.16666667, 0.53684211, 0.13684211, 0.0],
                [0.15789474, 0.09473684, 0.09473684, 0.24210526, 0.0],
                [0.05263158, 0.08421053, 0.14736842, 0.09473684, 0.0],
                [0.07216495, 0.05263158, 0.10526316, 0.06315789, 0.0],
                [0.0, 0.0, 0.0, 0.05263158, 0.0],
                [0.05263158, 0.05263158, 0.10526316, 0.09473684, 0.0],
                [0.04210526, 0.03157895, 0.05263158, 0.04210526, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [0.04210526, 0.04210526, 0.08421053, 0.07368421, 0.0],
                [0.0625, 0.04210526, 0.09473684, 0.05263158, 0.0],
            ]
        )

    def _encode_seq(self, seq):
        """Encode: A=0, T=1, G=2, C=3, U=1, N=4."""
        mapping = {"A": 0, "T": 1, "G": 2, "C": 3, "U": 1, "N": 4}
        return np.array([mapping.get(b.upper(), 4) for b in seq], dtype=int)

    def score(self, sequence):
        """Compute KSS for sequence.

        Parameters
        ----------
        sequence : str
            Length left_flank + 3 + right_flank. Codon centered.

        Returns
        -------
        float
            Normalized score [0, 1].
        """
        expected_len = self.left_flank + self.codon_len + self.right_flank
        assert len(sequence) == expected_len, (
            f"Length {len(sequence)} != {expected_len} " f"(left={self.left_flank}, right={self.right_flank})"
        )

        # Extract relevant weight rows
        start_row = 10 - self.left_flank
        end_row = start_row + expected_len
        weights = self.weight_matrix[start_row:end_row]

        # Encode and score
        codes = self._encode_seq(sequence)
        score = sum(weights[i, codes[i]] for i in range(len(codes)))

        # Normalize
        max_score = np.sum(weights.max(axis=1))
        return score / max_score if max_score > 0 else 0.0

    def score_batch(self, sequences):
        """Score multiple sequences.

        Parameters
        ----------
        sequences : list of str

        Returns
        -------
        np.ndarray
            Scores for each sequence.
        """
        return np.array([self.score(seq) for seq in sequences])


def kozak_weight_score(sequence, weight_matrix=None, left_flank=10, right_flank=10):
    """Quick-score function (convenience wrapper).

    Parameters
    ----------
    sequence : str
        DNA sequence centered on codon.
    weight_matrix : np.ndarray, optional
        Weight matrix. Uses default ATG if None.
    left_flank : int, default 10
    right_flank : int, default 10

    Returns
    -------
    float
        Normalized KSS [0, 1].
    """
    scorer = KozakWeightScore(weight_matrix=weight_matrix, left_flank=left_flank, right_flank=right_flank)
    return scorer.score(sequence)


def _iupac_match(observed, pattern):
    """True when each observed IUPAC code is a subset of the pattern code.

    Both ``observed`` and ``pattern`` are equal-length IUPAC strings. A
    consensus base may itself be ambiguous (e.g. ``R``); it matches when the
    set of bases it represents is contained in the pattern's set.
    """
    from sequana.iuapc import dna_ambiguities

    sets = {code: set(expansion.strip("[]")) for code, expansion in dna_ambiguities.items()}
    for obs, pat in zip(observed.upper(), pattern):
        if obs not in sets or pat not in sets:
            return False
        if not sets[obs] <= sets[pat]:
            return False
    return True


def classify_kozak_strength(motif, left=6):
    """Classify a Kozak consensus motif as optimal, strong, adequate or weak.

    The tiers follow Meijer and Thomas (2002), Biochem J,
    doi:10.1042/bj20011706. The two key nucleotides are a purine
    (``R`` = A/G) at the -3 position and a ``G`` at the +4 position (the base
    just after the ATG codon). U is treated as T. The 6 upstream positions
    (-6..-1) plus the +4 position are compared, the ATG codon being fixed:

    - ``GCCRCC``-AUG-``G`` -> ``"optimal"`` (full optimal context, both key
      nucleotides present);
    - ``NNNRNN``-AUG-``G`` -> ``"strong"`` (both key nucleotides present);
    - ``NNNRNN``-AUG-(A/C/U) or ``NNN(C/U)NN``-AUG-``G`` -> ``"adequate"``
      (only one of the two key nucleotides present; termed ``"moderate"`` in
      Hernandez et al. 2019, Cell);
    - ``NNN(C/U)NN``-AUG-(A/C/U) -> ``"weak"`` (both key nucleotides absent).

    :param str motif: a Kozak motif. The ATG start codon is located
        automatically (first ``ATG`` substring); if absent, ``left`` is used
        as the number of upstream positions and the motif is assumed to be the
        6 upstream bases immediately followed by the +4 base (no ATG codon).
    :param int left: number of upstream positions when no ATG is found.
    :return: one of ``"optimal"``, ``"strong"``, ``"adequate"`` or ``"weak"``.

    ::

        >>> from sequana.kozak import classify_kozak_strength
        >>> classify_kozak_strength("GCCRCCATGG")
        'optimal'
        >>> classify_kozak_strength("AAARAAATGG")
        'strong'
        >>> classify_kozak_strength("AAACAAATGG")
        'adequate'
        >>> classify_kozak_strength("AAACAAATGA")
        'weak'
    """
    # work in DNA alphabet (U -> T) and upper case
    motif = motif.upper().replace("U", "T")

    idx = motif.find("ATG")
    if idx >= 6:
        upstream = motif[idx - 6 : idx]
        plus4 = motif[idx + 3] if len(motif) > idx + 3 else ""
    else:
        # no ATG codon: motif is the upstream window followed by the +4 base
        upstream = motif[:left][-6:]
        plus4 = motif[left] if len(motif) > left else ""

    if len(upstream) < 6 or not plus4:
        raise ValueError("motif too short: need 6 upstream bases plus the +4 base")

    # the two key nucleotides: purine (R) at -3 (upstream index 3) and G at +4
    purine_m3 = _iupac_match(upstream[3], "R")
    g_p4 = _iupac_match(plus4, "G")

    if purine_m3 and g_p4:
        # both key nucleotides present: optimal if the full GCCRCC context
        # matches, otherwise strong
        return "optimal" if _iupac_match(upstream, "GCCRCC") else "strong"
    if purine_m3 or g_p4:
        return "adequate"
    return "weak"
