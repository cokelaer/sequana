#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import itertools

from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana.lazy import numpy as np

import colorlog

logger = colorlog.getLogger(__name__)

from sequana.macs3 import MACS3Reader


class PandasReader:
    def __init__(self, filename, sep="\t", **kwargs):
        self.df = pd.read_csv(filename, sep=sep, **kwargs)
        # If there is a header, let us strip it from spaces
        try:
            self.df.columns = [x.strip() for x in self.df.columns]
        except:
            pass

        # let us strip strings from spaces
        for x in self.df.columns:
            try:
                self.df[x] = self.df[x].str.strip()
            except:
                pass


class MultiHomer:
    def __init__(self):
        pass


class Homer(PandasReader):
    def __init__(self, filename):
        header = open("annotate_peaks/narrow/1_vs_6_7.txt").readline().strip().split("\t")[1:]
        super(Homer, self).__init__(filename, sep="\t", skiprows=1, header=None)
        self.df.columns = ["ID"] + header

    def pie_annotation(self):
        from collections import Counter

        counts = Counter(self.df.Annotation)
        labels = [x.split()[0] for x in counts.keys()]
        pylab.pie(counts.values(), labels=labels)


class IDR(PandasReader):
    """Reader for the output generated by idr package

    Can read the narrow or broad output transparently.

    Note that signalValue = rep1_signal + rep2_signal

    The score columns contains the scaled IDR value, min(int(log2(-125IDR), 1000).
    THis means that peaks with an IDR of 0 have a score of 1000,
    or peaks with an IDR of 0.05 0.05 have a score of int(-125log2(0.05)) = 540,
    and peaks with IDR of 1.0 have a score of 0.

        IDR, score
        0.0039, 1000
        0.01, 830
        0.05, 540
        0.1, 415
        0.5, 125

    This allows to differentiate those that crosses an IDR or 0.05.
    The final IDR is stored in global_idr

    """

    def __init__(self, filename, threshold=0.05):
        super(IDR, self).__init__(filename, sep="\t", header=None)
        self.threshold = threshold

        narrow_columns = [
            "chrom",
            "start",
            "end",
            "region_name",
            "score",
            "strand",
            "signalValue",
            "pvalue",
            "qvalue",
            "summit",
            "local_idr",
            "global_idr",
            "rep1_chrom_start",
            "rep1_chrom_end",
            "rep1_signal",
            "rep1_summit",
            "rep2_chrom_start",
            "rep2_chrom_end",
            "rep2_signal",
            "rep2_summit",
        ]
        broad_columns = [
            "chrom",
            "start",
            "end",
            "region_name",
            "score",
            "strand",
            "signalValue",
            "pvalue",
            "qvalue",
            "local_idr",
            "global_idr",
            "rep1_chrom_start",
            "rep1_chrom_end",
            "rep1_signal",
            "rep2_chrom_start",
            "rep2_chrom_end",
            "rep2_signal",
        ]
        try:
            self.df.columns = narrow_columns
        except:
            self.df.columns = broad_columns

        # add ranks for rep1/rep2
        self.df["rep1_rank"] = self.df["rep1_signal"].rank(ascending=False)
        self.df["rep2_rank"] = self.df["rep2_signal"].rank(ascending=False)
        self.df["idr"] = 10 ** -self.df["local_idr"]

    def _get_N_significant_peaks(self):
        return len(self.df.query("idr<@self.threshold"))

    N_significant_peaks = property(_get_N_significant_peaks)

    def IDR2score(self, IDR):
        if IDR == 0:
            return 1000
        return min(int(-125 * pylab.log2(IDR)), 1000)

    def score2IDR(self, score):
        return 2 ** (score / -125)

    def plot_ranks(self, filename=None, savefig=False):
        # ranks
        # the *score* columns contains the scaled IDR value, min(int(log2(-125IDR), 1000).
        # e.g. peaks with an IDR of 0 have a score of 1000, idr 0.05 have a score of
        # int(-125log2(0.05)) = 540, and idr 1.0 has a score of 0.
        df1 = self.df.query("score>540")
        df2 = self.df.query("score<=540")
        pylab.clf()
        pylab.plot(df1.rep1_rank, df1.rep2_rank, "ko", alpha=0.5, label="<0.05 IDR")
        pylab.plot(df2.rep1_rank, df2.rep2_rank, "ro", alpha=0.5, label=">=0.05 IDR")
        pylab.xlabel("Peak rank - replicate 1")
        pylab.ylabel("Peak rank - replicate 2")
        N = len(self.df)
        pylab.plot([0, N], [0, N], color="blue", alpha=0.5, ls="--")
        # pylab.xlim([0,1.05])
        # pylab.ylim([0,1.05])
        pylab.legend(loc="lower right")
        if savefig:
            pylab.savefig(filename)

    def plot_scores(self, filename=None, savefig=False):
        # scores
        from pylab import log10

        pylab.clf()
        pylab.plot(
            log10(self.df.query("score>540")["rep1_signal"]),
            log10(self.df.query("score>540")["rep2_signal"]),
            "ko",
            alpha=0.5,
            label="<0.05 IDR",
        )
        pylab.plot(
            log10(self.df.query("score<540")["rep1_signal"]),
            log10(self.df.query("score<540")["rep2_signal"]),
            "ro",
            alpha=0.5,
            label=">=0.05 IDR",
        )
        N = pylab.ylim()[1]
        pylab.plot([0, N], [0, N], color="blue", alpha=0.5, ls="--")
        pylab.xlabel("Rep1 log10 score")
        pylab.ylabel("Rep2 log10 score")
        pylab.legend(loc="lower right")
        if savefig:
            pylab.savefig(filename)

    def plot_rank_vs_idr_score(self, filename=None, savefig=False):
        # rank versus IDR scores
        f, axes = pylab.subplots(2, 1)
        df = self.df
        axes[0].plot(range(len(df)), df.sort_values(by="rep1_rank", ascending=False)["local_idr"], "o")
        axes[0].set_ylabel("log10 IDR for replicate 1")
        axes[0].axvline(len(self.df) - self.N_significant_peaks, color="b", ls="--")
        axes[1].plot(range(len(df)), df.sort_values(by="rep2_rank", ascending=False)["local_idr"], "ro")
        axes[1].set_ylabel("log10 IDR for replicate 2")
        axes[1].axvline(len(self.df) - self.N_significant_peaks, color="b", ls="--")
        if savefig:
            pylab.savefig(filename)

    def plot_idr_vs_peaks(self, filename=None, savefig=False):

        pylab.clf()
        X1 = pylab.linspace(0, self.threshold, 100)
        X2 = pylab.linspace(self.threshold, 1, 100)
        # convert local idr to proba

        df1 = self.df.query("idr<@self.threshold")
        df2 = self.df.query("idr>=@self.threshold")

        pylab.plot([sum(df1["idr"] < x) for x in X1], X1, "-", color="r", lw=2)
        shift = len(df1)

        pylab.plot([shift + sum(df2["idr"] < x) for x in X2], X2, "-", color="k", lw=2)
        pylab.xlabel("Number of significant peaks")
        pylab.ylabel("IDR")
        pylab.axhline(0.05, color="b", ls="--")
        pylab.axvline(self.N_significant_peaks, color="b", ls="--")
        if savefig:
            pylab.savefig(filename)


class FRiP(PandasReader):
    """Reader for output of FRiP computation in chipseq pipeline

    expected format::

        bamfile,count,in_peaks, FRiP,
        1_S1_L001.sorted.bam,3548090,53673,0.015127293839784221
        2_S2_L001.sorted.bam,3868608,58292,0.01506795209026089
        3_S3_L001.sorted.bam,4092990,50219,0.01226951446253228

    """

    def __init__(self, filename, design=None):
        super(FRiP, self).__init__(filename, sep=",")
        if design:
            self.design = ChIPExpDesign(design)
        else:
            self.design = None

    def plot(self):
        """"""
        if self.design:
            self.df["label"] = self.design.df["type"] + "/" + self.design.df["condition"]

        pylab.clf()
        MX = self.df.FRiP.max()
        MY = self.df["in_peaks"].max()
        pylab.plot([0, MX], [0, MY], ls="--", color="b", alpha=0.5)
        for label in self.df["label"].unique():
            self.df.query("label==@label").plot(x="FRiP", y="in_peaks", marker="o", lw=0, label=label, ax=pylab.gca())
        pylab.ylabel("Reads in peaks")
        pylab.xlabel("FRiP")
        pylab.xlim(0, pylab.xlim()[1])
        pylab.ylim(0, pylab.ylim()[1])
        pylab.grid()


class ChIPExpDesign(PandasReader):
    def __init__(self, filename):
        super(ChIPExpDesign, self).__init__(filename, sep=",")
        for x in ["ID", "condition", "type", "replicat"]:
            assert x in self.df.columns, x
        if self.df["ID"].duplicated().sum() > 0:
            duplicated = self.df["ID"][self.df["ID"].duplicated()].values
            raise ValueError("Input design file as duplicated IDs: {}".format(duplicated))
        # must have IP, the immunioprecipated samples
        # others names can be defined by the users and are considered as
        # different types of controls.
        if sum(self.df["type"] == "IP") == 0:
            raise ValueError(
                """No rows has the requested type name 'IP', which
is your immunoprecipated sample. Your design file should look like:

ID, type,    condition, replicat
1,  IP,      EXP1,      1
2,  control, EXP1,      2
3,  control, EXP1,      1


"""
            )

    def _get_conditions(self):
        return self.df["condition"].unique()

    conditions = property(_get_conditions)

    def _get_types(self):
        return self.df["type"].unique()

    types = property(_get_types)

    def get_idr_design(self, strict=True):
        """

        :param bool strict: if true returns IDR for replicates with same
            conditions, otherwise all combos of IP samples.
        """
        # IDR consists in comparing two inputs.
        results = {}
        if strict:
            groups = self.df[self.df["type"] == "IP"].groupby("condition")
            for cond, samples in groups.groups.items():
                if len(samples) == 2:
                    IDs = self.df.loc[groups.groups[cond]].ID.values
                    name = "_vs_".join([str(x) for x in IDs])
                    results[cond] = IDs
                elif len(samples) == 1:
                    print(f"condition {cond} has only 1 replicate. not included in IDR")
                else:
                    print(f"condition {cond} has more than 2 replicates. Not included in IDR")
                    # all combos ?
        return results

    def get_simple_comparisons(self):

        comparisons = {}
        for condition in self.conditions:
            df = self.df.query("condition==@condition")
            # TODO check there is at least one IP and one control
            IPs = df.query("type=='IP'")
            controls = df.query("type!='IP'")

            combos_IPs = []
            for x in range(len(IPs)):
                combos_IPs += list(itertools.combinations(IPs["ID"], x + 1))

            for control in controls["type"].unique():
                for combo in combos_IPs:
                    inputs = list(combo)
                    outputs = list(df.query("type==@control")["ID"].values)
                    name = "{}_vs_{}".format("_".join([str(x) for x in inputs]), "_".join([str(x) for x in outputs]))

                    inputs = [self.df.query("ID==@x").sample_name.values[0].strip() for x in inputs]
                    outputs = [self.df.query("ID==@x").sample_name.values[0].strip() for x in outputs]

                    comparisons[name] = {"inputs": inputs, "controls": outputs}
        return comparisons


class PeakConsensus:
    def __init__(self, f1, f2):
        self.df1 = MACS3Reader(f1).df
        self.df1["category"] = "first"
        self.df2 = MACS3Reader(f2).df
        self.df2["category"] = "second"
        self.df_merged = self.merge()

    def merge(self, overlap=0.2):
        df = pd.concat([self.df1, self.df2]).sort_values(["chr", "start"])
        # if overlap at least one base, we merge the peaks and label them with
        # common information, otherwise we report the original peak

        merged = []
        prev = None
        overlaps = 0
        N1 = 0
        N2 = 0
        N12 = 0
        skip_next = True
        for k, current in df.iterrows():
            if skip_next:
                prev = current
                skip_next = False
                continue

            # if current overlaps the prev start or end, there is overlap
            # or if current included in prev there current and prev overlaps
            if current["start"] <= prev["start"] and current["end"] >= prev["start"]:
                overlap = True
                N12 += 1
            elif current["start"] <= prev["end"] and current["end"] >= prev["end"]:
                overlap = True
                N12 += 1
            elif current["start"] >= prev["start"] and current["end"] <= prev["end"]:
                overlap = True
                N12 += 1
            else:
                overlap = False
                if prev["name"].startswith("1_vs_6_7"):
                    N1 += 1
                elif prev["name"].startswith("2_vs_6_7"):
                    N2 += 1

            if overlap:
                m = min(current["start"], prev["start"])
                M = max(current["end"], prev["end"])
                data = current.copy()
                data["start"] = m
                data["end"] = M
                data["stop"] = M  # FIXME same as end. decided on one value
                data["category"] = "both"
                merged.append(data)
                skip_next = True
            else:
                m = min(current["start"], prev["start"])
                M = max(current["end"], prev["end"])
                merged.append(prev)
                skip_next = False

            prev = current
        df = pd.DataFrame(merged)
        df = df.reset_index(drop=True)
        return df

    def plot_venn(self, title="", labels=[]):
        from sequana.viz.venn import plot_venn

        if labels is None:
            labels = [cond1, cond2]

        plot_venn(
            (
                set(self.df_merged.query("category in ['first', 'both']").index),
                set(self.df_merged.query("category in ['second', 'both']").index),
            ),
            labels=(("first", "second")),
        )

    def to_saf(self, filename):
        df = self.df_merged.reset_index()
        df["strand"] = "+"
        df["GeneID"] = df["index"]
        df[["GeneID", "chr", "start", "end", "strand", "category"]].to_csv(filename, sep="\t", index=None)

    def to_bed(self, filename):
        # output to be used by homer
        # annotatePeaks.pl file.bed file.fa -gid -gff file.gff
        # GeneID seems to have to start with Interval_ ??

        df = self.df_merged.reset_index()
        df["strand"] = "+"
        df["GeneID"] = ["Interval_{}".format(x) for x in df["index"]]
        mapper = {"both": 3, "first": 1, "second": 2}
        df["score"] = [mapper[x] for x in df["category"]]
        df[["chr", "start", "end", "GeneID", "score", "strand"]].to_csv(filename, sep="\t", index=None, header=False)


class PhantomPeaksReader:
    """Manipulate output of PhantomPeaks


    1_S1_L001.sorted.bam	3534910	150,250,275	0.97029159945329,0.968713084636106,0.968468832890883	80	0.9682
    442	1500	0.9669921	1.003412	2.6352	2

        The metrics file is tabulated

        * Filename
        * numReads: effective sequencing depth i.e. total number of mapped reads in the input file
        * estFragLen: comma separated strand cross-correlation peak(s) in decreasing order of correlation. In almost all cases, the top (rst) value in the list represents the predominant fragment length.
        * corr estFragLen: comma separated strand cross-correlation value(s) in decreasing order (col3 follows the same order)
        * phantomPeak: Read length/phantom peak strand shift
        * corr phantomPeak: Correlation value at phantom peak
        * argmin corr: strand shift at which cross-correlation is lowest
        * min corr: minimum value of cross-correlation
        * Normalized strand cross-correlation coecient (NSC) = COL4 / COL8. ;1=no
    enrichment. NSC >1.1 is good
        * Relative strand cross-correlation coecient (RSC) = (COL4 - COL8) / (COL6 -
    COL8); RSC=0 means no signal, <1 low quality and >1 means high enrichment.
    should aim at RSC>0.8
        *  QualityTag: Quality tag based on thresholded RSC (codes: -2:veryLow; -1:Low; 0:Medium; 1:High; 2:veryHigh)

    """

    def __init__(self, filename):
        self.df = pd.read_csv(filename, sep="\t", header=None)
        self.df.columns = [
            "filename",
            "num_reads",
            "estimated_fragment_length",
            "corr",
            "phantom_peak",
            "corr_phantom_peak",
            "argmin_corr",
            "min_corr",
            "NSC",
            "RSC",
            "quality_tag",
        ]

    def read(self, filename):
        df = pd.read_csv(filename, sep="\t", header=None)
        df.columns = [
            "filename",
            "num_reads",
            "estimated_fragment_length",
            "corr",
            "phantom_peak",
            "corr_phantom_peak",
            "argmin_corr",
            "min_corr",
            "NSC",
            "RSC",
            "quality_tag",
        ]
        self.df = self.df.append(df)
        self.df.reset_index(inplace=True, drop=True)

    def plot_RSC(self):
        self.df.RSC.plot(kind="bar")
        pylab.axhline(0.8, lw=2, color="r", ls="--")


class Phantom:
    """

    prints
    - the reference chromosome,
    - the starting position-1
    - position -1 + sequence length
    - N
    - 1000
    - + if flag is 16 (or equivalent) otherwise prints -


    c = Phantom()
    c.chromosomes
    data = c.get_data('NC_002506.1')
    mask = c.remove_anomalies(data)
    data = data[mask]
    c.scc(data)

    # Note that binning is set to 5 by default so this will comput from -500 to
    # 500
    X = range(-100,100)
    Y = [c.cor(x) for x in X]
    plot([x*5 for x in X], Y)


     c = chipseq.Phantom(binning=10, start=-500, stop=500)
    c.read_align("test.align")
    results, df = c.run()
    c.stats(results, df)



    """

    def __init__(self, bamfile=None, binning=5, start=-500, stop=1500):

        self.start = start
        self.stop = stop
        self.binning = binning

        cmd = """
        samtools view -F 0x0204 -o | awk 'BEGIN{FS="\t";OFS="\t"} {if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),"N","1000","-"} else {print $3,($4-1),($4-1+length($10)),"N","1000","+"}}' - > test.align
        """
        # shell(cmd)
        self.df = None

    def read_align(self, readfile):

        self.df = pd.read_csv(readfile, sep="\t", header=None)
        self.df.columns = ["ref", "start", "end", "dummy", "quality", "strand"]
        from pylab import median

        self.read_length = round(median(self.df["end"] - self.df["start"]))
        self.chromosomes = self.df["ref"].unique()

    def run(self):
        ## 10% of the time in self.get_data and 90 in cor()
        if self.df is None:
            print("call read_align() method to read alignement file")
            return
        m = int(self.start / self.binning)
        M = int(self.stop / self.binning)
        results = {}
        # because bins is set to 5, we actually go from m*5 to M*5
        X = range(m, M + 1, 1)
        Xreal = np.array(range(m * self.binning, (M + 1) * self.binning, self.binning))

        for chrom in self.chromosomes:
            # logger.info("Processing {}".format(chrom))
            data = self.get_data(chrom)
            L = len(data)
            self.scc(data)

            # shift correlation
            Y = [self.cor(x) for x in X]
            results[chrom] = {"data_length": L, "Y": np.array(Y), "X": Xreal}

        # weighted average usng orginal length of the chrmosomes
        weights = np.array([results[x]["data_length"] for x in self.chromosomes])
        weights = weights / sum(weights)

        self.results = results
        self.weights = weights
        # now the weighted cross correlation
        df_avc = pd.DataFrame([w * results[x]["Y"] for w, x in zip(weights, self.chromosomes)])
        df_avc = df_avc.T
        df_avc.index = Xreal
        return results, df_avc

    def stats(self, results, df_avc, bw=1):

        stats = {}
        stats["read_fragments"] = len(self.df)
        stats["fragment_length"] = self.read_length

        # average cross correlation across all chromosomes
        print("Read {} fragments".format(stats["read_fragments"]))
        print("ChIP data mean length: {}".format(self.read_length))
        # df_avc.sum(axis=1).plot()
        df = df_avc.sum(axis=1)
        corr_max = df.max()
        shift_max = df.idxmax()
        # note that in phantomPeak, they use the last value as min... not the
        # actual min. Not very important.
        corr_min = df.min()
        shift_min = df.idxmin()
        print("Maximum cross-correlation value: {:.5f}".format(corr_max))
        print("Maximum cross-correlation shift: {}".format(shift_max))
        print("Minimum cross-correlation value: {:.5f}".format(corr_min))
        print("Minimum cross-correlation shift: {}".format(shift_min))
        stats["shift_max"] = int(shift_max)  # to make it json serialisable
        stats["corr_max"] = corr_max

        # original code phantomPeak but always equal to 1 it range max >5 ??
        # default is 500 so sbw=1 whatsoever
        # sbw = 2 * floor(ceil(5/15000) / 2) + 1
        sbw = 1

        # here we could use a rolling mean
        # df.rolling(window=5, center=True).mean()

        # so following runnin mean is useless
        # cc$y <- runmean(cc$y,sbw,alg="fast")
        #

        # again, computation of bw but always equal to 1 ....
        # Compute cross-correlation peak
        #  bw <- ceiling(2/iparams$sep.range[2]) # crosscorr[i] is compared to crosscorr[i+/-bw] to find peaks
        # bw = 1
        # search for local peaks within bandwidth of bw = 1
        peakidx = df.diff(periods=bw) > 0
        peakidx = peakidx.astype(int).diff(periods=bw) == -1

        # the final bw points are NA and filled with False
        peakidx = peakidx.shift(-bw).fillna(False)

        df_peaks = df[peakidx]
        # when searching for max, exclude peaks from the excluded region
        exclusion_range = [10, self.read_length + 10]
        mask = np.logical_or(df_peaks.index < exclusion_range[0], df_peaks.index > exclusion_range[1])
        df_peaks = df_peaks[mask]

        #
        max_peak = df_peaks.max()
        shift_peak = df_peaks.idxmax()

        # now, we select peaks that are at least 90% of main peak and with shift
        # higher than main shift. why higher ?
        mask = np.logical_and(df_peaks > max_peak * 0.9, df_peaks.index >= shift_peak)
        best_df_peaks = df_peaks[mask]
        best = best_df_peaks.sort_values(ascending=False)[0:3]

        values = ",".join(["{:.5f}".format(x) for x in best.values])
        pos = ",".join([str(x) for x in best.index])
        print("Top 3 cross-correlation values: {}".format(values))
        print("Top 3 estimates for fragment length: {}".format(pos))

        # now the real window half size according to phantom peaks, not spp ...
        # min + (max-min)/3
        threshold = (df_peaks.max() - corr_min) / 3 + corr_min
        whs = df[df > threshold].index.max()

        # coming back to real cross correlation, identify peak in window
        # readlength +- 2*binning  !! not symmetry in phantompeak
        # x >= ( chip.data$read.length - round(2*binning) &
        # x <= ( chip.data$read.length + round(1.5*binning)

        binning = self.binning
        ph_min = self.read_length - round(2 * binning)
        ph_max = self.read_length + round(1.5 * binning)
        phantom = df[np.logical_and(df.index >= ph_min, df.index <= ph_max)]
        print("Phantom peak range detection:{}-{}".format(ph_min, ph_max))
        print("Phantom peak location:{}".format(phantom.idxmax()))
        print("Phantom peak Correlation: {:.5f}".format(phantom.max()))
        stats["phantom_corr"] = phantom.max()
        stats["phantom_location"] = int(phantom.idxmax())  # for json

        NSC = df_peaks.max() / phantom.max()
        # error in phatompeaks ?? is encoded as follows but no link with phantom
        # peak...
        # Another difference with phantom peak is that the min in phantom peak
        # is not the min but last value on the RHS so
        #    phantom_coeff = df_peaks.max() /  df.min()
        # is
        #    phantom_coeff = df_peaks.max() /  df.iloc[-1]
        NSC_spp = df_peaks.max() / df.iloc[-1]
        print("Normalized Strand cross-correlation coefficient (NSC): {:.5f} [{:.5f}]".format(NSC, NSC_spp))
        RSC = (df_peaks.max() - df.min()) / (phantom.max() - df.min())
        RSC_spp = (df_peaks.max() - df.iloc[-1]) / (phantom.max() - df.iloc[-1])
        print("Relative Strand cross-correlation Coefficient (RSC): {:.5f} [{:.5f}]".format(RSC, RSC_spp))

        if RSC > 0 and RSC < 0.25:
            tag = -2
        elif RSC >= 0.25 and RSC < 0.5:
            tag = -1
        elif RSC >= 0.5 and RSC < 1:
            tag = 0
        elif RSC >= 1 and RSC < 1.5:
            tag = 1
        elif RSC >= 1.5:
            tag = 2
        print("Phantom Peak Quality Tag: {}".format(tag))

        pylab.clf()
        df.plot()
        ##df_peaks.plot(marker="o", lw=0)
        ylim = pylab.ylim()
        # pylab.axvline(whs, ls='--', color='k', lw=1)
        Y0, Y1 = pylab.ylim()
        pylab.plot([phantom.idxmax(), phantom.idxmax()], [Y0, phantom.max()], ls="--", color="k", lw=1)
        pylab.plot([df.idxmax(), df.idxmax()], [Y0, df.max()], ls="--", color="r", lw=2)
        # pylab.fill_betweenx(ylim, 10,85, color='grey', alpha=0.5)
        pylab.ylim(ylim)
        pylab.ylabel("Cross-correlation")
        pylab.xlabel("strand-shift: {}bp\nNSC={:.5f}, RSC={:.5f}, Qtag={}".format(best.index[0], NSC, RSC, tag))
        pylab.xlim(self.start, self.stop)
        pylab.grid(True, zorder=-20)
        try:
            pylab.tight_layout()
        except:
            pass
        stats["NSC"] = NSC
        stats["RSC"] = RSC
        stats["Qtag"] = tag
        return stats

    def get_data(self, chrname, remove_anomalies=True):

        # Could be done once for all in read_alignment
        df = self.df.query("ref==@chrname")

        # first the fragment position, shifting - strand by fragment length
        data = np.array([x if z == "+" else -y for x, y, z in zip(df["start"], df["end"], df["strand"])])

        # sort by absolute position
        res = pd.DataFrame(data)
        res.columns = ["x"]
        res["abs"] = res["x"].abs()
        res = res.sort_values("abs")
        del res["abs"]

        if remove_anomalies:
            mask = self.remove_anomalies(res)
            res = res[mask]

        return res

    def remove_anomalies(self, data, bin=1, trim_fraction=1e-3, z=5, return_indecies=False):

        zo = z * 3

        x = data["x"]

        from numpy import floor, ceil, sqrt

        # the frequencies, sorted from smaller to larger
        tt = floor(x / bin).value_counts().sort_values()

        # sort and select first 99.9%
        # floor to agree with phantom
        stt = tt.iloc[0 : int(floor(len(tt) * (1 - trim_fraction)))]

        mtc = stt.mean()
        var_base = 0.1
        tcd = sqrt(stt.var() + var_base)

        thr = max(1, ceil(mtc + z * tcd))
        thr_o = max(1, ceil(mtc + zo * tcd))

        # filter tt
        tt = tt[tt > thr]

        # get + and - tags
        tp = tt.index

        pti = set(tp[tp > 0])
        pti2 = set([-x for x in tp[tp < 0]])
        it = pti.intersection(pti2)

        it = sorted(it.union(set(tt[tt > thr_o].index)))
        # print(len(it), sum(it))

        # sit <- c(it,(-1)*it);
        # sit =  [-x for x in it] + list(it)
        sit = list(it)[::-1] + [-x for x in list(it)[::-1]]
        # print(len(sit), sum(sit))

        if bin > 1:
            # From 0 to 5+1 to agree with phantom R code
            sit = sorted([x * 5 + y for y in range(0, 6) for x in sit])

        sit = set(sit)
        return [False if x in sit else True for x in data["x"]]

    def scc(self, data, llim=10):

        tt = (np.sign(data["x"]) * np.floor(data["x"].abs() / self.binning + 0.5)).value_counts()

        mu = tt.mean()
        tt = tt[tt < llim * mu]
        #
        tc = tt.index

        pdata = tt[tc > 0]
        ndata = tt[tc < 0]

        ntv = ndata.values
        ptv = pdata.values

        self.pdata = pdata
        self.ndata = ndata
        r1 = min(-1 * ndata.index.max(), pdata.index.min())
        r2 = max(-1 * ndata.index.min(), pdata.index.max())
        l = r2 - r1 + 1
        self.L = l
        mp = ptv.sum() * self.binning / l
        mn = ntv.sum() * self.binning / l
        ptv = ptv - mp
        ntv = ntv - mn
        self.mp = mp
        self.mn = mn
        ntc = -tt.index[[True if x < 0 else False for x in tt.index]].values
        ptc = tt.index[[True if x > 0 else False for x in tt.index]].values

        self.ntc = ntc
        self.ptc = ptc

        self.ntv = ntv
        self.ptv = ptv
        from numpy import sqrt

        ss = sqrt((sum(ptv ** 2) + (l - len(ptv)) * mp ** 2) * (sum(ntv ** 2) + (l - len(ntv)) * mn ** 2))
        self.ss = ss
        self.nn = pd.DataFrame(self.ntv, index=self.ntc)

        logger.info("mp={}; mn={}, ss={}".format(self.mp, self.mn, self.ss))
        logger.info(
            "ptv sum={}; ptv len ={}, ntv sum={} ntv len {}".format(
                sum(self.ptv), len(self.ptv), sum(self.ntv), len(self.ntv)
            )
        )
        logger.info(
            "ptc sum={}; ptc len ={}, ntc sum={} ntc len {}".format(
                sum(self.ptc), len(self.ptc), sum(self.ntc), len(self.ntc)
            )
        )

    def cor(self, s):
        logger.info("shift {}".format(s))
        p = pd.DataFrame(self.ptv, index=self.ptc + s)
        n = self.nn

        self.X = p[p.index.isin(n.index)]
        self.compX = p[p.index.isin(n.index) == False]
        self.Y = n[n.index.isin(p.index)]
        self.compY = n[n.index.isin(p.index) == False]

        R = self.L - len(self.ptv) - len(self.ntc) + len(self.X)

        # using dataframe
        # !!!!!!!!!!!! multiplication takes into account indexing
        XY = (self.X * self.Y).sum()[0]
        XX = self.compX.sum()[0] * self.mn
        YY = self.compY.sum()[0] * self.mp
        corr = (XY - XX - YY + self.mp * self.mn * R) / self.ss

        return corr
