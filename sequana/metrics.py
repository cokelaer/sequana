def lempel_ziv_complexity(seq):
    """

    Lempel-Ziv complexity is a measure of sequence compressibility—how "random" or "repetitive" a sequence is. In the context of DNA, it captures how diverse the patterns are. A low LZ complexity suggests repeats or low complexity regions, while a high LZ complexity suggests a more random or information-rich sequence.



    """

    seq = seq.upper()
    n = len(seq)
    if n == 0:
        return 0.0
    i = 0
    complexity = 1
    dictionary = set()
    w = seq[0]

    while i < n - 1:
        i += 1
        w += seq[i]
        if w not in dictionary:
            dictionary.add(w)
            complexity += 1
            w = ""
            if i + 1 < n:
                w = seq[i + 1]
                i += 1
    return complexity / len(seq)


def compute_bendability(seq, scale, window=3):
    """

    DNA bendability measures the local flexibility of a DNA sequence—its ability to bend or deform. This is a crucial feature in:


            Brukner et al. (1995)

    """

    seq = seq.upper()
    scores = []
    for i in range(len(seq) - window + 1):
        tri = seq[i : i + window]
        score = scale.get(tri, None)
        if score is not None:
            scores.append(score)
        else:
            scores.append(None)
    return scores


brukner_flexibility = {
    "AAA": 0.069,
    "AAC": 0.055,
    "AAG": 0.059,
    "AAT": 0.063,
    "ACA": 0.072,
    "ACC": 0.054,
    "ACG": 0.051,
    "ACT": 0.058,
    "AGA": 0.070,
    "AGC": 0.056,
    "AGG": 0.058,
    "AGT": 0.064,
    "ATA": 0.066,
    "ATC": 0.053,
    "ATG": 0.060,
    "ATT": 0.065,
    "CAA": 0.073,
    "CAC": 0.057,
    "CAG": 0.061,
    "CAT": 0.067,
    "CCA": 0.070,
    "CCC": 0.052,
    "CCG": 0.050,
    "CCT": 0.057,
    "CGA": 0.069,
    "CGC": 0.055,
    "CGG": 0.056,
    "CGT": 0.062,
    "CTA": 0.065,
    "CTC": 0.051,
    "CTG": 0.059,
    "CTT": 0.066,
    "GAA": 0.068,
    "GAC": 0.054,
    "GAG": 0.058,
    "GAT": 0.064,
    "GCA": 0.071,
    "GCC": 0.053,
    "GCG": 0.050,
    "GCT": 0.057,
    "GGA": 0.067,
    "GGC": 0.055,
    "GGG": 0.057,
    "GGT": 0.063,
    "GTA": 0.064,
    "GTC": 0.052,
    "GTG": 0.059,
    "GTT": 0.065,
    "TAA": 0.067,
    "TAC": 0.053,
    "TAG": 0.060,
    "TAT": 0.066,
    "TCA": 0.069,
    "TCC": 0.054,
    "TCG": 0.052,
    "TCT": 0.058,
    "TGA": 0.068,
    "TGC": 0.056,
    "TGG": 0.059,
    "TGT": 0.064,
    "TTA": 0.065,
    "TTC": 0.051,
    "TTG": 0.060,
    "TTT": 0.067,
}

helix_twist = {
    "AA": 36.0,
    "AC": 34.5,
    "AG": 36.6,
    "AT": 32.2,
    "CA": 35.6,
    "CC": 33.1,
    "CG": 30.0,
    "CT": 35.2,
    "GA": 35.7,
    "GC": 36.9,
    "GG": 32.8,
    "GT": 34.3,
    "TA": 34.0,
    "TC": 35.4,
    "TG": 36.1,
    "TT": 36.0,
}


def compute_helix_twist(seq, scale):
    """

    Helix twist refers to the rotational angle between adjacent base pairs in the DNA double helix, typically measured in degrees. It reflects how tightly the DNA is twisted and is important for:

    - DNA supercoiling
    - Protein-DNA interactions
    - Nucleosome positioning
    - Local helical structure variations (e.g., bends, kinks)

    Typical Values: Canonical B-DNA: ~36° per base pair step
    Varies with sequence: ~32°–40° depending on dinucleotide
    """

    seq = seq.upper()
    scores = []
    for i in range(len(seq) - 1):
        dinuc = seq[i : i + 2]
        score = scale.get(dinuc, None)
        if score is not None:
            scores.append(score)
        else:
            scores.append(None)
    return scores


# Dinucleotide Roll angles (degrees), free-DNA / crystal averages from
# Olson et al. (1998) PNAS 95:11163 Table 2. Roll is the main sequence-dependent
# component of intrinsic (static) curvature; Tilt is small and near-symmetric.
# Palindromic complements share a value (5'-XY-3' and its reverse complement).
dinucleotide_roll = {
    "AA": 0.7,
    "TT": 0.7,
    "AC": 0.7,
    "GT": 0.7,
    "AG": 4.5,
    "CT": 4.5,
    "AT": 1.1,
    "CA": 4.7,
    "TG": 4.7,
    "CC": 3.6,
    "GG": 3.6,
    "CG": 5.4,
    "GA": 1.9,
    "TC": 1.9,
    "GC": 0.3,
    "TA": 3.3,
}


def compute_curvature(roll, twist=None, tilt=None, window=11, bp_per_turn=10.5):
    """Local intrinsic (static) curvature of the DNA helix axis.

    Static curvature arises when successive base-pair-step deflections (mainly
    Roll, with a small Tilt component) add *in phase* along the helix. Each step
    is a small bend whose direction rotates with the accumulated helical twist;
    A-tracts (or any bend source) phased at the helical repeat (~10.5 bp) sum
    constructively -> high curvature, the classic intrinsically-curved-DNA
    signature. Roll-vector / wedge model: Trifonov & Sussman (1980) PNAS 77:3816;
    Bolshoy et al. (1991) PNAS 88:2312; Olson et al. (1998) PNAS 95:11163.

    Parameters
    ----------
    roll : sequence of float
        Per-step Roll angle in degrees (length = n_steps = len(seq) - 1), e.g. a
        DNAshapeR Roll profile or a :data:`dinucleotide_roll` lookup. NaN -> 0.
    twist : sequence of float, optional
        Per-step helical Twist (degrees) used for phasing. If None, a canonical
        constant twist of 360 / `bp_per_turn` deg/step is used.
    tilt : sequence of float, optional
        Per-step Tilt (degrees), the orthogonal bend component. If None, a
        roll-only model is used (tilt = 0).
    window : int
        Sliding window length in steps (default 11 ~ one helical turn). The
        curvature at position i is the magnitude of the vector sum of the step
        bends inside the window centred on i.
    bp_per_turn : float
        Helical repeat used to derive the default constant twist.

    Returns
    -------
    numpy.ndarray
        Per-position curvature magnitude (degrees accumulated over `window`),
        one value per step; positions closer than window//2 to either end are NaN.
    """
    import numpy as np

    roll = np.nan_to_num(np.asarray(roll, dtype=float), nan=0.0)
    n = len(roll)
    if n == 0:
        return np.array([], dtype=float)

    if twist is None:
        phase = np.arange(n) * np.deg2rad(360.0 / bp_per_turn)
    else:
        tw = np.nan_to_num(np.asarray(twist, dtype=float), nan=360.0 / bp_per_turn)
        # phase at step i = cumulative twist up to (but not including) step i
        phase = np.concatenate([[0.0], np.cumsum(np.deg2rad(tw))[:-1]])
        phase = phase[:n]

    # each step contributes a bend vector: Roll along the phase direction,
    # Tilt orthogonal to it (phase + 90 deg).
    vx = roll * np.cos(phase)
    vy = roll * np.sin(phase)
    if tilt is not None:
        tl = np.nan_to_num(np.asarray(tilt, dtype=float), nan=0.0)[:n]
        vx = vx + tl * np.cos(phase + np.pi / 2)
        vy = vy + tl * np.sin(phase + np.pi / 2)

    half = window // 2
    csx = np.concatenate([[0.0], np.cumsum(vx)])
    csy = np.concatenate([[0.0], np.cumsum(vy)])
    out = np.full(n, np.nan)
    for i in range(half, n - half):
        a, b = i - half, i + half + 1
        out[i] = np.hypot(csx[b] - csx[a], csy[b] - csy[a])
    return out


# Nearest-neighbour base-pair opening free energies dG37 (kcal/mol per 5'-XY-3'
# step). These are the *cost to separate* a duplex step and equal minus the
# SantaLucia & Hicks (2004) unified duplex-formation dG. AT-rich steps are
# cheapest to open -> the first to melt under negative superhelical stress.
NN_OPENING_DG = {
    "AA": 1.00,
    "TT": 1.00,
    "AT": 0.88,
    "TA": 0.58,
    "CA": 1.45,
    "TG": 1.45,
    "GT": 1.44,
    "AC": 1.44,
    "CT": 1.28,
    "AG": 1.28,
    "GA": 1.30,
    "TC": 1.30,
    "CG": 2.17,
    "GC": 2.24,
    "GG": 1.84,
    "CC": 1.84,
}


def compute_sidd(
    seq,
    sigma=-0.06,
    temperature=310.15,
    nucleation=10.84,
    superhelical_constant=2350.0,
    bp_per_turn=10.5,
    max_run=1000,
    opening_energy=None,
):
    """Approximate stress-induced duplex destabilization (SIDD) profile.

    Two-state, single-open-run approximation of Benham's SIDD analysis (Benham
    1993 PNAS 90:2999; Bi & Benham 2004; Zhabinskaya & Benham SIST, 2015). In a
    topologically closed domain of length N held at superhelical density `sigma`,
    negative supercoiling stores torsional free energy that can be relieved by
    locally opening (denaturing) the duplex. AT-rich, low-stability runs open
    first. For a run of n base pairs the incremental free energy relative to the
    fully closed state is

        dG(run) = a + sum_i b_i + (K/2)[(alpha + n/h)^2 - alpha^2]

    where
      a         nucleation cost of the two open/closed junctions (kcal/mol),
      b_i       nearest-neighbour opening energy of step i (:data:`NN_OPENING_DG`),
      alpha     = sigma * N / h, the imposed linking difference (turns),
      h         base pairs per helical turn,
      K         = superhelical_constant * RT / N, the residual torsional stiffness.

    Opening n bp relaxes the duplex by ~n/h turns, driving (alpha + n/h) toward 0
    (alpha < 0), which *releases* torsional energy and offsets the base-pair
    opening cost. Where the release exceeds the cost, dG(x) < 0 -> a spontaneously
    destabilized SIDD site. The per-base value returned is the minimum run energy
    over all runs covering that base.

    This is the tractable two-state approximation (one open run), appropriate for
    *ranking* destabilization susceptibility along or across sequences. For
    absolute equilibrium opening probabilities use the full partition-function
    SIDD (WebSIDD / SIST).

    Parameters
    ----------
    seq : str
        DNA sequence (ACGT; other characters use the mean opening energy).
    sigma : float
        Superhelical density (negative for underwinding; -0.06 is physiological).
    temperature : float
        Temperature in kelvin (used only via RT in K).
    nucleation : float
        Junction nucleation free energy `a` (kcal/mol).
    superhelical_constant : float
        Benham's torsional coefficient (RT units); K = this * RT / N.
    bp_per_turn : float
        Helical repeat `h`.
    max_run : int
        Longest open run considered (caps cost at O(N * max_run)).
    opening_energy : dict, optional
        Override the nearest-neighbour opening-energy table.

    Returns
    -------
    numpy.ndarray
        Per-base incremental free energy dG(x) in kcal/mol (length = len(seq)).
        Low / negative values mark easily destabilized (SIDD) sites.
    """
    import numpy as np

    seq = seq.upper()
    table = opening_energy if opening_energy is not None else NN_OPENING_DG
    n = len(seq)
    if n < 2:
        return np.full(n, np.nan)

    mean_b = float(np.mean(list(table.values())))
    b = np.array([table.get(seq[i : i + 2], mean_b) for i in range(n - 1)], dtype=float)
    prefix = np.concatenate([[0.0], np.cumsum(b)])  # prefix[k] = sum of first k steps

    RT = 0.0019872041 * temperature  # kcal/mol
    K = superhelical_constant * RT / n
    alpha = sigma * n / bp_per_turn  # imposed linking difference (turns)

    steps = n - 1
    lengths = np.arange(1, min(max_run, steps) + 1)  # run length in *steps*
    # torsional term as a function of run length (bp opened ~ run length in steps)
    torsion = 0.5 * K * ((alpha + lengths / bp_per_turn) ** 2 - alpha**2)

    gpos = np.full(n, np.inf)
    for i in range(steps):
        L = min(len(lengths), steps - i)
        # dG of runs starting at step i, for every run length 1..L
        run_dg = nucleation + (prefix[i + 1 : i + 1 + L] - prefix[i]) + torsion[:L]
        best = float(np.min(run_dg))
        j = int(np.argmin(run_dg))  # best run spans steps i..i+j
        # assign the deepest well from this start to the bases it covers
        span = slice(i, i + j + 2)  # j+1 steps -> j+2 bases
        cur = gpos[span]
        np.minimum(cur, best, out=cur)
        gpos[span] = cur
    gpos[np.isinf(gpos)] = np.nan
    return gpos
