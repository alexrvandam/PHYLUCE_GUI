#!/usr/bin/env python3
"""
Post-UCE / post-phylogeny helper

- Reads insilico-incomplete.incomplete matrix from PHYLUCE.
- Computes number of loci per taxon.
- Computes locus occupancy (proportion of taxa per locus).
- Produces:
    * Bar chart of loci per taxon.
    * Histogram of locus occupancy.
- Writes an adjusted genomes.conf that keeps only taxa above a minimum loci threshold.

Can be used both:
    - as a CLI script, and
    - via run_post_uce_analysis(...) from the GUI.
"""

import argparse
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import matplotlib
matplotlib.use("Agg")  # headless, safe for GUI
import matplotlib.pyplot as plt


def load_incomplete_matrix(matrix_path: Path) -> Tuple[List[str], List[List[int]]]:
    """
    Load the insilico-incomplete.incomplete file.

    Assumes:
        first row: "locus <taxon1> <taxon2> ..."
        subsequent rows: locusID <0/1/...> <0/1/...> ...

    Returns:
        taxa: list of taxon names (columns)
        matrix: list of rows; each row is a list of ints (0/1) per taxon
    """
    taxa: List[str] = []
    matrix: List[List[int]] = []

    with matrix_path.open() as fh:
        header = fh.readline().strip().split()
        if len(header) < 2:
            raise ValueError(f"Matrix header looks wrong in {matrix_path}")
        taxa = header[1:]  # skip 'locus'

        for line in fh:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) != len(header):
                # Skip malformed lines quietly
                continue
            values = []
            for raw in parts[1:]:
                try:
                    v = int(raw)
                    v = 1 if v != 0 else 0
                except ValueError:
                    # Treat obvious missing codes as 0
                    if raw in {".", "NA", "NaN", "none", "None"}:
                        v = 0
                    else:
                        v = 1
                values.append(v)
            matrix.append(values)

    return taxa, matrix


def compute_counts_and_occupancy(
    taxa: List[str],
    matrix: List[List[int]]
) -> Tuple[Dict[str, int], List[float]]:
    """
    From a 0/1 matrix:

    - count loci per taxon
    - compute locus occupancy (fraction of taxa present)
    """
    taxon_counts: Dict[str, int] = {t: 0 for t in taxa}
    locus_occupancy: List[float] = []

    n_taxa = len(taxa)
    if n_taxa == 0:
        return taxon_counts, locus_occupancy

    for row in matrix:
        if len(row) != n_taxa:
            continue
        # per locus
        present = sum(row)
        locus_occupancy.append(present / float(n_taxa))
        # per taxon
        for t, v in zip(taxa, row):
            if v:
                taxon_counts[t] += 1

    return taxon_counts, locus_occupancy


def plot_taxon_locus_counts(
    taxon_counts: Dict[str, int],
    out_png: Path
) -> None:
    """
    Bar plot: taxa (x) vs loci count (y).
    """
    if not taxon_counts:
        return

    # Sort taxa by loci descending for readability
    taxa_sorted = sorted(taxon_counts.keys(), key=lambda t: taxon_counts[t], reverse=True)
    counts = [taxon_counts[t] for t in taxa_sorted]

    fig_width = max(8, len(taxa_sorted) * 0.3)

    plt.figure(figsize=(fig_width, 6))
    plt.bar(range(len(taxa_sorted)), counts)
    plt.xticks(range(len(taxa_sorted)), taxa_sorted, rotation=90)
    plt.ylabel("Number of loci")
    plt.title("Loci per taxon")
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=300)
    plt.close()


def plot_locus_occupancy(
    locus_occupancy: List[float],
    out_png: Path
) -> None:
    """
    Histogram of per-locus occupancy (0–1).
    """
    if not locus_occupancy:
        return

    plt.figure(figsize=(6, 4))
    plt.hist(locus_occupancy, bins=20)
    plt.xlabel("Proportion of taxa present per locus")
    plt.ylabel("Number of loci")
    plt.title("Locus occupancy")
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=300)
    plt.close()


def write_adjusted_genomes_conf(
    genomes_conf: Path,
    taxon_counts: Dict[str, int],
    min_loci: int,
    out_conf: Path,
) -> None:
    """
    Read the original genomes.conf and write a filtered version
    with only taxa that have >= min_loci loci.
    """
    keep_taxa = {t for t, n in taxon_counts.items() if n >= min_loci}

    lines: List[str] = []
    with genomes_conf.open() as fh:
        for line in fh:
            if not line.strip():
                lines.append(line)
                continue
            if line.lstrip().startswith("#"):
                lines.append(line)
                continue
            parts = line.strip().split()
            if not parts:
                continue
            taxon = parts[0]
            if taxon in keep_taxa:
                lines.append(line)

    out_conf.parent.mkdir(parents=True, exist_ok=True)
    with out_conf.open("w") as out:
        out.write("# Auto-generated adjusted genomes.conf\n")
        out.write(f"# Only taxa with >= {min_loci} loci retained\n\n")
        for line in lines:
            out.write(line)

    return


def run_post_uce_analysis(
    project_root: str,
    min_loci: int = 100,
    logger=None,
    insilico_dir: Optional[str] = None,
    genomes_conf: Optional[str] = None,
    output_dir: Optional[str] = None,
) -> None:
    """
    Main entry point for GUI and CLI.

    project_root: project folder (with probe-design/, probe-design-test/ etc.).
    min_loci: minimum loci per taxon for inclusion in adjusted genomes.conf.
    """
    pr = Path(project_root).expanduser().resolve()

    if insilico_dir is None:
        insilico_dir = pr / "probe-design-test" / "taxon-sets" / "insilico-incomplete"
    else:
        insilico_dir = Path(insilico_dir)

    if genomes_conf is None:
        genomes_conf = pr / "probe-design" / "genomes.conf"
    else:
        genomes_conf = Path(genomes_conf)

    if output_dir is None:
        output_dir = insilico_dir
    else:
        output_dir = Path(output_dir)

    log = logger or DummyLogger()

    matrix_file = insilico_dir / "insilico-incomplete.incomplete"
    if not matrix_file.exists():
        raise FileNotFoundError(f"Matrix file not found: {matrix_file}")

    if not genomes_conf.exists():
        raise FileNotFoundError(f"genomes.conf not found: {genomes_conf}")

    log.info(f"[STEP 4] Loading incomplete matrix: {matrix_file}")
    taxa, matrix = load_incomplete_matrix(matrix_file)
    log.info(f"[STEP 4] Found {len(taxa)} taxa and {len(matrix)} loci in matrix")

    taxon_counts, locus_occupancy = compute_counts_and_occupancy(taxa, matrix)

    # Some quick summary
    if locus_occupancy:
        avg_occ = sum(locus_occupancy) / len(locus_occupancy)
        log.info(f"[STEP 4] Mean locus occupancy: {avg_occ:.3f}")
    else:
        log.info("[STEP 4] No loci in matrix; nothing to plot.")

    # Plots
    bar_png = output_dir / "taxon_locus_counts.png"
    occ_png = output_dir / "locus_occupancy_hist.png"
    log.info(f"[STEP 4] Writing taxon-locus bar plot: {bar_png}")
    plot_taxon_locus_counts(taxon_counts, bar_png)

    log.info(f"[STEP 4] Writing locus-occupancy histogram: {occ_png}")
    plot_locus_occupancy(locus_occupancy, occ_png)

    # Adjusted genomes.conf
    base_name = genomes_conf.stem
    adjusted_conf = output_dir / f"{base_name}.min{min_loci}.conf"
    log.info(f"[STEP 4] Writing adjusted genomes.conf: {adjusted_conf}")
    write_adjusted_genomes_conf(genomes_conf, taxon_counts, min_loci, adjusted_conf)

    log.info("[STEP 4] Post-UCE analysis complete.")


class DummyLogger:
    """Fallback logger if none is provided."""

    def info(self, msg):
        print(msg)

    def error(self, msg):
        print(msg)

    def exception(self, msg):
        print(msg)


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Post-UCE helper: plots + adjusted genomes.conf"
    )
    p.add_argument(
        "--project-root",
        required=True,
        help="Project root (with probe-design/ and probe-design-test/).",
    )
    p.add_argument(
        "--insilico-dir",
        default=None,
        help="Optional explicit path to taxon-sets/insilico-incomplete.",
    )
    p.add_argument(
        "--genomes-conf",
        default=None,
        help="Optional explicit path to genomes.conf.",
    )
    p.add_argument(
        "--output-dir",
        default=None,
        help="Optional explicit output dir for PNGs + adjusted genomes.conf "
             "(default: insilico-incomplete).",
    )
    p.add_argument(
        "--min-loci",
        type=int,
        default=100,
        help="Minimum loci per taxon to include in adjusted genomes.conf (default: 100).",
    )

    args = p.parse_args(argv)

    run_post_uce_analysis(
        project_root=args.project_root,
        min_loci=args.min_loci,
        logger=None,
        insilico_dir=args.insilico_dir,
        genomes_conf=args.genomes_conf,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()

