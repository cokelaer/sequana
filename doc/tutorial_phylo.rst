.. _tutorial_phylo:

Phylogenetic analysis
=====================

**Time:** 10 min | **Requires:** phylogenetic trees (Newick) and alignments

Sequana provides tools for working with phylogenetic trees and sequence
alignments. Parse, manipulate, and analyze evolutionary relationships.


Parsing phylogenetic trees (Newick format)
===========================================

Load a tree from Newick format::

    from sequana.phylo import Tree

    # Standard Newick: (A:1.0,B:1.0)root:0.0;
    tree = Tree.from_newick("(A:1.0,B:0.5)root:0.0;")

    # Get tree statistics
    print(f"Leaves (taxa): {tree.leaves()}")
    print(f"Num nodes: {tree.count_nodes()}")

    # Get root node
    root = tree.root
    print(f"Root: {root.name}")
    print(f"Children: {[c.name for c in root.children]}")

Access nodes and compute distances::

    # All leaf nodes
    leaves = tree.get_leaves()
    for leaf in leaves:
        print(f"  {leaf.name}: branch length {leaf.branch_length}")

    # Distance between two taxa
    dist_a_b = tree.distance("A", "B")
    print(f"Distance A-B: {dist_a_b:.2f}")

    # Distance from node to root
    dist_to_root = tree.distance_to_root(leaves[0])
    print(f"Distance to root: {dist_to_root:.2f}")

Re-root tree and get subtrees::

    # Re-root at taxon A
    new_tree = tree.reroot("A")
    print(f"New root: {new_tree.root.name}")

    # Get subtree containing specific taxa
    subtree = tree.subtree(["A", "B"])
    print(f"Subtree leaves: {subtree.leaves()}")


Multiple sequence alignments
=============================

Parse alignment files (PHYLIP, Stockholm, Nexus, FASTA)::

    from sequana.alignment import Alignment, parse_phylip, parse_stockholm

    # Parse PHYLIP format
    aln = parse_phylip("alignment.phy")
    print(f"Sequences: {len(aln)}")
    print(f"Alignment length: {aln.length()}")
    print(f"Sequence names: {aln.names}")

Access sequences and compute properties::

    # Get individual sequences
    for name in aln.names:
        seq = aln.sequences[name]
        gc_content = seq.count("G") + seq.count("C") / len(seq)
        print(f"  {name}: GC = {gc_content:.1%}")

    # Consensus sequence (most common residue per position)
    consensus = aln.consensus()
    print(f"Consensus: {consensus}")

    # Identify conserved columns
    for pos in range(aln.length()):
        col = aln.get_column(pos)
        if len(set(col.values())) == 1:  # All same residue
            print(f"Conserved at position {pos}: {list(col.values())[0]}")

Get alignment statistics::

    # Residue frequency per sequence
    for name, seq in aln.sequences.items():
        num_gaps = seq.count("-")
        num_n = seq.count("N")
        print(f"  {name}: {num_gaps} gaps, {num_n} ambiguous")

    # Pairwise identity
    seq1, seq2 = aln.names[0], aln.names[1]
    s1, s2 = aln.sequences[seq1], aln.sequences[seq2]
    matches = sum(1 for a, b in zip(s1, s2) if a == b and a != "-")
    identity = matches / len(s1)
    print(f"Identity {seq1}-{seq2}: {identity:.1%}")


Combining trees and alignments
===============================

Link phylogenetic relationships to sequences::

    from sequana.phylo import Tree
    from sequana.alignment import parse_phylip

    tree = Tree.from_newick("(human:0.1,chimp:0.08)primate:0.0;")
    aln = parse_phylip("primates.phy")

    # Verify alignment contains all taxa
    taxa = tree.leaves()
    missing = [t for t in taxa if t not in aln.names]
    if missing:
        print(f"Warning: missing sequences for {missing}")

    # Extract alignment subset for clade
    clade = tree.subtree(["human", "chimp"])
    clade_seqs = {name: aln.sequences[name] for name in clade.leaves()}
    print(f"Clade alignment: {len(clade_seqs)} sequences")


See :ref:`references_genomics` for detailed API documentation.
