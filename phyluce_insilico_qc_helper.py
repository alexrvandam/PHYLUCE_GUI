#!/usr/bin/env python3
"""
Quick QC helper for PHYLUCE in-silico runs.

Features
--------
1. Parse `phyluce_assembly_get_fastas_from_match_counts.log` to get
   UCE loci per taxon and:
   - write a TSV table
   - draw a barplot (taxon vs UCE loci)

2. Parse `phyluce_align_get_align_summary_data.log` to read lines like:
   [Matrix 50%]          1048 alignments
   and:
   - write a TSV table of percent vs loci
   - draw a simple line plot (matrix % vs number of loci)

3. Build an `adjusted_genomes.conf` by keeping only taxa with
   >= --min-loci UCE loci, so you can re-run the in-silico pipeline
   without poorly sequenced taxa.

Typical usage
-------------
python phyluce_insilico_qc_helper.py \
    --fastas-log phyluce_assembly_get_fastas_from_match_counts.log \
    --align-summary-log phyluce_align_get_align_summary_data.log \
    --genomes-conf ../../probe-design/genomes.conf \
    --min-loci 100 \
    --out-prefix insilico_qc

This will create:
- insilico_qc_uce_counts.tsv
- insilico_qc_uce_counts.png
- insilico_qc_matrix_summary.tsv
- insilico_qc_matrix_summary.png
- adjusted_genomes.conf
"""

import argparse
import collections
import math
import os
import re
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt


# ---------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------


def parse_fastas_log(path: str) -> Dict[str, int]:
    """
    Parse phyluce_assembly_get_fastas_from_match_counts.log to get
    UCE loci per taxon.

    We try a couple of regex patterns to be robust to minor format
    differences. You can tweak the patterns if your log looks
    different.
    """
    counts = collections.OrderedDict()

    # Example patterns we try to catch:
    #   genus_species1: 1031 UCE loci
    #   There are 1031 UCE loci for genus_species1
    pat1 = re.compile(r"(\S+):\s+(\d+)\s+UCE", re.IGNORECASE)
    pat2 = re.compile(
        r"There are\s+(\d+)\s+UCE.*?\sfor\s+(\S+)", re.IGNORECASE
    )

    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            m1 = pat1.search(line)
            m2 = pat2.search(line)

            if m1:
                taxon = m1.group(1)
                n = int(m1.group(2))
                counts[taxon] = n
            elif m2:
                n = int(m2.group(1))
                taxon = m2.group(2)
                counts[taxon] = n

    if not counts:
        raise RuntimeError(
            f"Could not find any 'UCE loci' lines in {path}. "
            "Open the log and adjust the regex in parse_fastas_log()."
        )

    return counts


def parse_align_summary_log(path: str) -> List[Tuple[int, int]]:
    """
    Parse phyluce_align_get_align_summary_data.log to extract lines like:

        [Matrix 50%]          1048 alignments

    Returns a list of (percent, n_alignments).
    """
    pattern = re.compile(
        r"\[Matrix\s+(\d+)%\]\s+(\d+)\s+alignments", re.IGNORECASE
    )
    entries: List[Tuple[int, int]] = []

    with open(path, "r") as fh:
        for line in fh:
            m = pattern.search(line)
            if m:
                pct = int(m.group(1))
                n = int(m.group(2))
                entries.append((pct, n))

    if not entries:
        raise RuntimeError(
            f"Did not find any '[Matrix X%] Y alignments' lines in {path}. "
            "Open the log and adjust parse_align_summary_log() if needed."
        )

    # Sort by percent ascending, just to be nice
    entries.sort(key=lambda x: x[0])
    return entries


# ---------------------------------------------------------------------
# Genomes.conf rewriting
# ---------------------------------------------------------------------


def write_adjusted_genomes(
    genomes_conf: str, out_conf: str, good_taxa: List[str]
) -> None:
    """
    Write adjusted_genomes.conf keeping only taxa in `good_taxa`.

    Expects genomes.conf lines like:
        taxon_name:/full/path/to/taxon_name.2bit
    plus blank lines or '#' comments.
    """
    good_set = set(good_taxa)

    with open(genomes_conf, "r") as fin, open(out_conf, "w") as fout:
        for line in fin:
            stripped = line.strip()
            # Pass through blank lines and comments
            if not stripped or stripped.startswith("#"):
                fout.write(line)
                continue

            # Everything before the first ":" is the taxon
            if ":" in stripped:
                taxon = stripped.split(":", 1)[0]
            else:
                # Weird line, just keep it
                fout.write(line)
                continue

            if taxon in good_set:
                fout.write(line)
            else:
                # quietly drop bad taxa
                pass


# ---------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------


def make_barplot(
    data_dict: Dict[str, int],
    title: str,
    xlabel: str,
    ylabel: str,
    out_png: str,
    rotate_xticks: bool = True,
) -> None:
    taxa = list(data_dict.keys())
    values = [data_dict[t] for t in taxa]

    if not taxa:
        return

    # Make figure width scale with number of taxa
    fig_width = max(6, len(taxa) * 0.4)
    fig, ax = plt.subplots(figsize=(fig_width, 4))
    ax.bar(range(len(taxa)), values)
    ax.set_xticks(range(len(taxa)))
    ax.set_xticklabels(
        taxa,
        rotation=90 if rotate_xticks else 0,
        fontsize=6,
    )
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def make_matrix_plot(
    matrix_entries: List[Tuple[int, int]], out_png: str
) -> None:
    """
    Simple line plot: x = matrix percent, y = number of loci/alignments.
    """
    if not matrix_entries:
        return

    percents = [p for (p, _) in matrix_entries]
    loci = [n for (_, n) in matrix_entries]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(percents, loci, marker="o")
    ax.set_xlabel("Matrix completeness (%)")
    ax.set_ylabel("Number of loci (alignments)")
    ax.set_title("Matrix completeness vs number of loci")
    ax.grid(True, linestyle="--", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(
        description="QC helper for PHYLUCE in-silico pipeline "
        "(UCE per taxon + matrix occupancy + adjusted genomes.conf)."
    )

    ap.add_argument(
        "--fastas-log",
        required=True,
        help="Path to phyluce_assembly_get_fastas_from_match_counts.log",
    )
    ap.add_argument(
        "--align-summary-log",
        required=False,
        help="Optional path to phyluce_align_get_align_summary_data.log "
        "for matrix occupancy summary.",
    )
    ap.add_argument(
        "--genomes-conf",
        required=True,
        help="Path to original genomes.conf used for in-silico run.",
    )
    ap.add_argument(
        "--min-loci",
        type=int,
        default=100,
        help="Minimum UCE loci per taxon to keep in adjusted_genomes.conf "
        "(default: 100).",
    )
    ap.add_argument(
        "--out-prefix",
        default="insilico_qc",
        help="Prefix for output files (default: insilico_qc).",
    )

    args = ap.parse_args()

    # 1. UCE counts per taxon (pre-MAFFT)
    print(f"[INFO] Parsing UCE counts per taxon from: {args.fastas_log}")
    counts = parse_fastas_log(args.fastas_log)

    # Write TSV
    tsv_path = f"{args.out_prefix}_uce_counts.tsv"
    with open(tsv_path, "w") as out:
        out.write("taxon\tuce_loci_pre_mafft\n")
        for taxon, n in counts.items():
            out.write(f"{taxon}\t{n}\n")
    print(f"[INFO] Wrote UCE counts table: {tsv_path}")

    # Barplot
    png_bar = f"{args.out_prefix}_uce_counts.png"
    make_barplot(
        counts,
        title="UCE loci per taxon (pre-MAFFT)",
        xlabel="Taxon",
        ylabel="Number of UCE loci",
        out_png=png_bar,
    )
    print(f"[INFO] Wrote barplot: {png_bar}")

    # 2. Adjusted genomes.conf
    good_taxa = [t for t, n in counts.items() if n >= args.min_loci]
    print(
        f"[INFO] Taxa with >= {args.min_loci} UCE loci: "
        f"{len(good_taxa)} / {len(counts)}"
    )
    for t in good_taxa:
        print(f"  - {t}")

    adjusted_conf = "adjusted_genomes.conf"
    write_adjusted_genomes(args.genomes_conf, adjusted_conf, good_taxa)
    print(f"[INFO] Wrote adjusted genomes.conf: {adjusted_conf}")

    # 3. Matrix completeness summary (optional)
    if args.align_summary_log:
        print(
            f"[INFO] Parsing matrix completeness from: {args.align_summary_log}"
        )
        matrix_entries = parse_align_summary_log(args.align_summary_log)

        matrix_tsv = f"{args.out_prefix}_matrix_summary.tsv"
        with open(matrix_tsv, "w") as out:
            out.write("percent\talignments\n")
            for pct, n in matrix_entries:
                out.write(f"{pct}\t{n}\n")
        print(f"[INFO] Wrote matrix summary table: {matrix_tsv}")

        matrix_png = f"{args.out_prefix}_matrix_summary.png"
        make_matrix_plot(matrix_entries, matrix_png)
        print(f"[INFO] Wrote matrix summary plot: {matrix_png}")

        print("\n[INFO] Use this table/plot to decide a target threshold.")
        print("      e.g. if 50% -> ~500 loci and 70% -> ~200 loci,")
        print("      pick the percent that gives you ~400–500 loci.")


if __name__ == "__main__":
    main()

