#!/usr/bin/env python3
"""
phyluce_one_passV2.py

Simple GUI + CLI wrapper around parts of the phyluce Tutorial IV workflow:
- Starting from cleaned draft genomes
- Simulating reads with ART
- Aligning reads to a base genome with stampy + samtools
- Converting BAM -> BED -> merged/filtered BED
- Building multi-merge tables and designing temporary baits with phyluce

This is intentionally opinionated and focuses on automating directory
creation and common parameters, while still exposing key thresholds that
control "overlap" / sharing of loci across taxa.
"""

import argparse
import logging
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional
import shutil
import tempfile

try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    TK_AVAILABLE = True
except Exception:
    TK_AVAILABLE = False




# -------------------- Dataclasses & helpers -------------------- #
# Resolve default Stampy path relative to this script, so that
# a bundled ./stampy-1.0.28/stampy.py works regardless of cwd.
SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_STAMPY_PATH = str(SCRIPT_DIR / "stampy-1.0.28" / "stampy.py")
# Full command to run Stampy via its conda env, with the script path quoted
DEFAULT_STAMPY_CMD = f'conda run -n stampy_py27 python "{DEFAULT_STAMPY_PATH}"'

@dataclass
class PipelineConfig:
    genomes_root: Path
    project_root: Path
    base_taxon: str
    exemplar_taxa: List[str]

    # Core resources
    threads: int = 8

    # ART read simulation
    art_read_len: int = 100
    art_coverage: float = 2.0
    art_insert_mean: int = 200
    art_insert_sd: int = 150

    # Stampy
    stampy_substitution_rate: float = 0.05
    stampy_insert_size: int = 400
    stampy_extra_opts: str = ""   # e.g. '--sensitive --prunebelow=0.01'

    # BED / locus filtering
    bedtools_merge_distance: int = 0  # 0 = only overlapping intervals
    strip_min_length: int = 80
    strip_filter_mask: float = 0.25

    # How strict to be about shared loci across exemplar taxa
    multi_merge_specific_counts: int = 5  # e.g. base+5
    genome_seq_buffer_to: int = 160      # bp surrounding each locus

    # Temporary bait design
    tiling_density: int = 3
    temp_masking: float = 0.25
    probe_prefix: str = "uce-"
    design_name: str = "custom-v1"
    designer_name: str = "user"

    # In-silico contig<->bait step (phyluce_assembly_match_contigs_to_probes)
    # (parameter is exposed but not used yet; easy to add later)
    min_coverage: int = 67

    # External tools (assumed on PATH by default)
    bwa_path: str = "conda run -n phyluce-1.7.3 bwa"
    samtools_path: str = "conda run -n phyluce-1.7.3 samtools"
    bamtools_path: str = "conda run -n phyluce-1.7.3 bamtools"
    bcftools_path: str = "conda run -n phyluce-1.7.3 bcftools"
    art_path: str = "conda run -n phyluce-1.7.3 art_illumina"
    stampy_path: str = DEFAULT_STAMPY_CMD
    bedtools_path: str = "conda run -n phyluce-1.7.3 bedtools"
    fatotwobit_path: str = "conda run -n phyluce-1.7.3 faToTwoBit"

    # If phyluce scripts are not on PATH, set a prefix like "/path/to/phyluce/bin/"
    phyluce_prefix: str = "conda run -n phyluce-1.7.3"

    # Alignment restart behaviour
    # If False (default): reuse existing SAMs and only rebuild BAMs.
    # If True: delete/rebuild SAMs with stampy as well.
    overwrite_existing_sam: bool = False

    # Execution control
    dry_run: bool = False

#    def phyluce_cmd(self, script: str) -> str:
##        """Return the full command name for a phyluce script."""
#        if self.phyluce_prefix:
#            return str(Path(self.phyluce_prefix) / script)
#        return script

    def phyluce_cmd(self, script: str) -> str:
        """Return the full command (with optional prefix) for a phyluce script."""
        if self.phyluce_prefix:
             # Allow prefix with or without trailing space
            return f"{self.phyluce_prefix.strip()} {script}"
        return script


def setup_logger(log_path: Optional[Path] = None) -> logging.Logger:
    logger = logging.getLogger("phyluce_wrapper")
    logger.setLevel(logging.INFO)
    logger.handlers = []

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    ch = logging.StreamHandler(stream=sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if log_path is not None:
        fh = logging.FileHandler(log_path)
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


def run_command(
    cmd: str,
    cwd: Optional[Path],
    cfg: PipelineConfig,
    logger: logging.Logger,
) -> None:
    """Run a shell command with logging, optionally in dry-run mode."""
    logger.info(f"CMD: {cmd}")
    if cfg.dry_run:
        return
    result = subprocess.run(
        cmd,
        shell=True,
        cwd=str(cwd) if cwd is not None else None,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {result.returncode}: {cmd}")


def sanitize_taxon_name(path: Path, logger: logging.Logger) -> str:
    """Sanitize a genome filename into a phyluce-friendly taxon name.

    - lower-case
    - no spaces
    - no weird punctuation (non-alnum -> underscore)
    - reasonably short
    """
    base = path.stem
    sanitized = []
    for ch in base:
        if ch.isalnum():
            sanitized.append(ch.lower())
        elif ch in (" ", ".", "-"):
            sanitized.append("_")
        else:
            sanitized.append("_")
    name = "".join(sanitized)
    if len(name) > 15:
        logger.warning(
            f"Taxon name '{name}' derived from '{base}' is long; "
            f"truncating to 15 characters for phyluce compatibility."
        )
        name = name[:15]
    if any(c.isupper() for c in base) or " " in base:
        logger.info(
            f"Sanitized genome name '{base}' -> '{name}' "
            "(phyluce prefers short, lowercase, no spaces)."
        )
    return name


def discover_genomes(genomes_root: Path, logger: logging.Logger) -> Dict[str, Path]:
    """Find FASTA-like genomes under genomes_root and map to taxon names."""
    exts = {".fa", ".fasta", ".fna"}
    genomes: Dict[str, Path] = {}
    for path in genomes_root.rglob("*"):
        if path.is_file() and path.suffix.lower() in exts:
            taxon = sanitize_taxon_name(path, logger)
            if taxon in genomes:
                logger.warning(
                    f"Multiple genome files map to taxon '{taxon}' "
                    f"({genomes[taxon]} and {path}); keeping the first."
                )
                continue
            genomes[taxon] = path
    if not genomes:
        raise FileNotFoundError(
            f"No genome FASTA files (.fa/.fasta/.fna) found under {genomes_root}"
        )
    logger.info(f"Discovered {len(genomes)} genomes: {', '.join(sorted(genomes))}")
    return genomes


# -------------------- Pipeline step helpers -------------------- #

def setup_project_structure(
    cfg: PipelineConfig,
    genomes: Dict[str, Path],
    logger: logging.Logger,
) -> Dict[str, Path]:
    """Create project directory structure and symlink genomes.

    Returns a mapping from taxon -> genome path within the project.
    """
    cfg.project_root.mkdir(parents=True, exist_ok=True)
    genomes_dir = cfg.project_root / "genomes"
    genomes_dir.mkdir(exist_ok=True)

    project_genomes: Dict[str, Path] = {}

    for taxon, src in genomes.items():
        taxon_dir = genomes_dir / taxon
        taxon_dir.mkdir(exist_ok=True)
        dest = taxon_dir / f"{taxon}.fasta"
        if not dest.exists():
            try:
                os.symlink(src, dest)
                logger.info(f"Symlinked genome {src} -> {dest}")
            except OSError:
                logger.info(f"Copying genome {src} -> {dest}")
                if not cfg.dry_run:
                    dest.write_bytes(src.read_bytes())
        else:
            logger.info(f"Genome already present in project: {dest}")
        project_genomes[taxon] = dest

    return project_genomes


def convert_genomes_to_2bit(
    cfg: PipelineConfig,
    project_genomes: Dict[str, Path],
    logger: logging.Logger,
) -> None:
    """Convert each genome to UCSC .2bit format using faToTwoBit."""
    logger.info("Converting genomes to .2bit format")
    for taxon, fasta in project_genomes.items():
        twobit = fasta.with_suffix(".2bit")
        if twobit.exists():
            logger.info(f"{taxon}: .2bit already exists ({twobit}), skipping")
            continue
        cmd = f"{cfg.fatotwobit_path} {fasta} {twobit}"
        run_command(cmd, cwd=fasta.parent, cfg=cfg, logger=logger)


def simulate_reads(
    cfg: PipelineConfig,
    project_genomes: Dict[str, Path],
    logger: logging.Logger,
) -> None:
    """
    Simulate reads from exemplar genomes using ART, as in Tutorial IV.

    This version is restart-friendly: if the final gzipped read file for a
    taxon already exists, that taxon is skipped so that we don't re-run ART
    or hit gzip "file exists" prompts when resuming after a crash.
    """
    reads_dir = cfg.project_root / "reads"
    reads_dir.mkdir(exist_ok=True)
    logger.info("Simulating reads with ART for exemplar taxa")

    for taxon in cfg.exemplar_taxa:
        # Don't simulate reads from the base taxon
        if taxon == cfg.base_taxon:
            continue

        # ---- NEW: resume-safety check ------------------------------------
        # If the final gzipped file already exists, assume this taxon is done
        final_gz = reads_dir / f"{taxon}-pe{cfg.art_read_len}-reads.fq.gz"
        if final_gz.exists():
            logger.info(
                f"Reads already simulated for {taxon}: {final_gz} exists, skipping."
            )
            continue
        # ------------------------------------------------------------------

        genome = project_genomes[taxon]
        out_prefix = reads_dir / f"{taxon}-pe{cfg.art_read_len}-reads"
        cmd = (
            f"{cfg.art_path} "
            f"--paired "
            f"--in {genome} "
            f"--out {out_prefix} "
            f"--len {cfg.art_read_len} "
            f"--fcov {cfg.art_coverage} "
            f"--mflen {cfg.art_insert_mean} "
            f"--sdev {cfg.art_insert_sd} "
            f"-ir 0.0 -ir2 0.0 -dr 0.0 -dr2 0.0 "
            f"-qs 100 -qs2 100 -na"
        )
        run_command(cmd, cwd=reads_dir, cfg=cfg, logger=logger)

        # Merge R1/R2 and gzip, as in the tutorial
        fq1 = reads_dir / f"{taxon}-pe{cfg.art_read_len}-reads1.fq"
        fq2 = reads_dir / f"{taxon}-pe{cfg.art_read_len}-reads2.fq"
        merged = reads_dir / f"{taxon}-pe{cfg.art_read_len}-reads.fq"

        cmd_merge = f"cat {fq1} {fq2} > {merged}"
        run_command(cmd_merge, cwd=reads_dir, cfg=cfg, logger=logger)
        run_command(f"rm {fq1} {fq2}", cwd=reads_dir, cfg=cfg, logger=logger)
        run_command(f"gzip {merged}", cwd=reads_dir, cfg=cfg, logger=logger)



def prepare_base_genome(cfg: PipelineConfig,
                        project_genomes: Dict[str, str],
                        logger: logging.Logger) -> None:
    """
    Copy base genome into base/ directory and build Stampy index+hash.
    This is a one-time preparation step per base_taxon.
    """
    logger.info("==== Preparing base genome for stampy (one-time) ====")

    if cfg.base_taxon not in project_genomes:
        raise RuntimeError(
            f"Configured base_taxon '{cfg.base_taxon}' not found among discovered genomes: "
            f"{', '.join(sorted(project_genomes.keys()))}"
        )

    genome_path = Path(project_genomes[cfg.base_taxon])
    base_dir = cfg.project_root / "base"
    base_dir.mkdir(parents=True, exist_ok=True)

    base_fa = base_dir / genome_path.name

    if not base_fa.exists():
        logger.info(f"Copying base genome {genome_path} -> {base_fa}")
        base_fa.write_bytes(genome_path.read_bytes())
    else:
        logger.info(f"Base genome already present: {base_fa}")

    # --- NEW: skip Stampy indexing if index/hash already exist ---
    stidx = base_dir / f"{cfg.base_taxon}.stidx"
    sthash = base_dir / f"{cfg.base_taxon}.sthash"

    if stidx.exists() and sthash.exists():
        logger.info(
            f"Stampy index ({stidx.name}) and hash ({sthash.name}) already exist; "
            "skipping base genome preparation."
        )
        return

    if stidx.exists():
        logger.info(
            f"Stampy index {stidx.name} already exists; will NOT rerun 'stampy -G'."
        )
    if sthash.exists():
        logger.info(
            f"Stampy hash {sthash.name} already exists; will NOT rerun 'stampy -H'."
        )

    logger.info("Preparing base genome for stampy")

    # Build genome index (.stidx) if needed
    if not stidx.exists():
        cmd_G = (
            f'{cfg.stampy_path} --species="{cfg.base_taxon}" '
            f'--assembly="{cfg.base_taxon}" '
            f"-G {cfg.base_taxon} {base_fa.name}"
        )
        run_command(cmd_G, cwd=base_dir, cfg=cfg, logger=logger)
    else:
        logger.info("Skipping 'stampy -G' because .stidx is already present.")

    # Build hash (.sthash) if needed
    if not sthash.exists():
        cmd_H = f"{cfg.stampy_path} -g {cfg.base_taxon} -H {cfg.base_taxon}"
        run_command(cmd_H, cwd=base_dir, cfg=cfg, logger=logger)
    else:
        logger.info("Skipping 'stampy -H' because .sthash is already present.")


def align_reads_with_stampy(
    cfg: PipelineConfig,
    logger: logging.Logger,
) -> None:
    """
    Align simulated reads to the base genome using stampy, then convert SAM -> BAM
    with samtools.

    Restart-friendly behaviour:
      - If a SAM file already exists and overwrite_existing_sam == False,
        we *skip stampy* and only rebuild the BAM with samtools view -Sb -o.
      - If overwrite_existing_sam == True and a SAM exists, we delete it and
        rerun stampy.
      - BAMs are always regenerated from SAM, so corrupted BAMs from older
        runs get overwritten automatically.
    """
    reads_dir = cfg.project_root / "reads"
    base_dir = cfg.project_root / "base"
    align_dir = cfg.project_root / "alignments"
    align_dir.mkdir(exist_ok=True)

    logger.info("Aligning reads to base genome with stampy + samtools")
    for taxon in cfg.exemplar_taxa:
        if taxon == cfg.base_taxon:
            continue

        taxon_dir = align_dir / taxon
        taxon_dir.mkdir(exist_ok=True)

        reads = reads_dir / f"{taxon}-pe{cfg.art_read_len}-reads.fq.gz"
        if not reads.exists():
            logger.warning(f"Reads not found for {taxon}: {reads} (skipping taxon)")
            continue

        sam_path = taxon_dir / f"{taxon}-to-{cfg.base_taxon}.sam"
        bam_path = taxon_dir / f"{taxon}-to-{cfg.base_taxon}.bam"

        # Decide whether to (re)run stampy
        if sam_path.exists():
            if cfg.overwrite_existing_sam:
                logger.info(
                    f"SAM exists for {taxon} ({sam_path}) and "
                    f"overwrite_existing_sam=True; deleting and re-running stampy."
                )
                if not cfg.dry_run:
                    try:
                        sam_path.unlink()
                    except OSError as e:
                        logger.warning(f"Could not remove existing SAM {sam_path}: {e}")
                run_stampy = True
            else:
                logger.info(
                    f"SAM already exists for {taxon}: {sam_path}; "
                    f"skipping stampy and only rebuilding BAM with samtools."
                )
                run_stampy = False
        else:
            run_stampy = True

        # Run stampy if needed to create SAM
        if run_stampy:
            logger.info(
                f"Aligning reads for {taxon} against base {cfg.base_taxon} with stampy."
            )
            cmd_stampy = (
                f"{cfg.stampy_path} "
                f"--maxbasequal 93 "
                f"-g {base_dir / cfg.base_taxon} "
                f"-h {base_dir / cfg.base_taxon} "
                f"--substitutionrate={cfg.stampy_substitution_rate} "
                f"-t{cfg.threads} "
                f"--insertsize={cfg.stampy_insert_size} "
                f"{cfg.stampy_extra_opts} "
                f"-M {reads} "
                f"-o {sam_path}"
            )
            run_command(cmd_stampy, cwd=cfg.project_root, cfg=cfg, logger=logger)

        # At this point, we expect the SAM to exist
        if not sam_path.exists():
            raise FileNotFoundError(
                f"Expected SAM file for {taxon} not found after stampy step: {sam_path}"
            )

        # Always rebuild BAM from SAM using samtools view -Sb -o
        logger.info(f"Converting SAM to BAM for {taxon} with samtools.")
        cmd_samtools = (
            f"{cfg.samtools_path} view -Sb {sam_path} -o {bam_path}"
        )
        run_command(cmd_samtools, cwd=cfg.project_root, cfg=cfg, logger=logger)
        logger.info(f"Finished BAM for {taxon}: {bam_path}")


def filter_bams_and_link(
    cfg: PipelineConfig,
    logger: logging.Logger,
) -> None:
    """Filter BAMs to mapped reads and symlink them into alignments/all."""
    align_dir = cfg.project_root / "alignments"
    all_dir = align_dir / "all"
    all_dir.mkdir(exist_ok=True)

    logger.info("Filtering BAMs to mapped reads (-F 4) and linking into alignments/all")

    for taxon in cfg.exemplar_taxa:
        if taxon == cfg.base_taxon:
            continue
        bam = align_dir / taxon / f"{taxon}-to-{cfg.base_taxon}.bam"
        m_bam = align_dir / taxon / f"{taxon}-to-{cfg.base_taxon}-MAPPING.bam"

        cmd = f"{cfg.samtools_path} view -h -F 4 -b {bam} -o {m_bam}"
        run_command(cmd, cwd=cfg.project_root, cfg=cfg, logger=logger)
        run_command(f"rm {bam}", cwd=cfg.project_root, cfg=cfg, logger=logger)

        link = all_dir / f"{taxon}-to-{cfg.base_taxon}-MAPPING.bam"
        if not link.exists():
            try:
                os.symlink(m_bam, link)
            except OSError:
                if not cfg.dry_run:
                    link.write_bytes(m_bam.read_bytes())


def bams_to_beds(
    cfg: PipelineConfig,
    logger: logging.Logger,
) -> None:
    """Convert each BAM in alignments/all to BED12 using bedtools bamtobed."""
    bed_dir = cfg.project_root / "bed"
    bed_dir.mkdir(exist_ok=True)
    all_dir = cfg.project_root / "alignments" / "all"

    logger.info("Converting BAMs to BED with bedtools bamtobed")
    for bam in all_dir.glob("*.bam"):
        bed_path = bed_dir / f"{bam.name}.bed"
        if bed_path.exists():
            logger.info(f"BED already exists, skipping: {bed_path}")
            continue
        cmd = f"{cfg.bedtools_path} bamtobed -i {bam} -bed12 > {bed_path}"
        run_command(cmd, cwd=cfg.project_root, cfg=cfg, logger=logger)


def sort_beds(cfg: PipelineConfig, logger: logging.Logger) -> None:
    """Sort BED files (bamtobed output) by coordinate.

    We *only* sort the primary BEDs that come directly from BAMs ("*.bam.bed").
    Derived files such as *.clean3col.bed, *.strip.bed, or any *.merge.bed
    are intentionally ignored so they do not get re-sorted/re-merged and
    accidentally fed back into later steps.
    """
    bed_dir = cfg.project_root / "bed"
    bed_dir.mkdir(exist_ok=True)

    logger.info("Sorting BED files by coordinate")

    # Only the original bamtobed outputs:
    for bed in sorted(bed_dir.glob("*.bam.bed")):
        sort_path = bed_dir / f"{bed.stem}.sort.bed"
        if sort_path.exists():
            logger.info(f"Sorted BED already exists, skipping: {sort_path}")
            continue

        cmd = f"{cfg.bedtools_path} sort -i {bed} > {sort_path}"
        run_command(cmd, cwd=cfg.project_root, cfg=cfg, logger=logger)


def merge_beds(cfg: PipelineConfig, logger: logging.Logger) -> None:
    """Merge sorted BED intervals for each taxon.

    We only merge the coordinate-sorted BEDs that come directly from BAMs
    ("*.bam.sort.bed"). This prevents recursively merging derived BEDs such as
    *.clean3col.sort.bed or *.strip.sort.bed, which would otherwise lead to
    filenames like "*.sort.merge.sort.merge.bed" and redundant processing.
    """
    bed_dir = cfg.project_root / "bed"
    bed_dir.mkdir(exist_ok=True)

    logger.info("Merging sorted BED intervals with bedtools merge")

    # Only the primary sorted BEDs
    for sort_bed in sorted(bed_dir.glob("*.bam.sort.bed")):
        merge_path = bed_dir / f"{sort_bed.stem}.merge.bed"
        if merge_path.exists():
            logger.info(f"Merged BED already exists, skipping: {merge_path}")
            continue

        cmd = f"{cfg.bedtools_path} merge -i {sort_bed} > {merge_path}"
        run_command(cmd, cwd=cfg.project_root, cfg=cfg, logger=logger)


def strip_masked_loci(cfg: PipelineConfig, logger: logging.Logger) -> Path:
    """Strip heavily masked / short loci from merged BEDs using phyluce.

    For each taxon-specific merged BED ("*.bam.sort.merge.bed") we:
      1. Create a 3-column, whitespace-cleaned BED ("*.clean3col.bed").
      2. Run phyluce_probe_strip_masked_loci_from_set on that cleaned BED.
      3. Collect the resulting "*.strip.bed" files into a simple "beds" config
         file that later phyluce steps can consume.
    """
    bed_dir = cfg.project_root / "bed"
    base_twobit = (
        cfg.project_root / "genomes" / cfg.base_taxon / f"{cfg.base_taxon}.2bit"
    )

    logger.info(
        f"Stripping masked loci with phyluce_probe_strip_masked_loci_from_set "
        f"(filter-mask={cfg.strip_filter_mask}, min-length={cfg.strip_min_length})"
    )

    # Only the primary merged BEDs derived from BAMs
    for merge_bed in sorted(bed_dir.glob("*.bam.sort.merge.bed")):
        # 1) Normalise to a strict 3-column BED so phyluce doesn't choke
        clean_bed = bed_dir / f"{merge_bed.stem}.clean3col.bed"
        if not clean_bed.exists():
            logger.info(f"Creating 3-column BED from {merge_bed} -> {clean_bed}")
            with merge_bed.open() as inp, clean_bed.open("w") as outp:
                for raw in inp:
                    line = raw.strip()
                    if not line:
                        continue
                    parts = line.split()
                    if len(parts) < 3:
                        # skip weird / malformed lines
                        continue
                    chromo, start, end = parts[:3]
                    outp.write(f"{chromo}\t{start}\t{end}\n")
        else:
            logger.info(f"Cleaned 3-column BED already exists, using: {clean_bed}")

        # 2) Run phyluce on the cleaned BED
        strip_path = bed_dir / f"{merge_bed.stem}.strip.bed"
        if strip_path.exists():
            logger.info(f"Stripped BED already exists, skipping: {strip_path}")
        else:
            cmd = (
                f"{cfg.phyluce_cmd('phyluce_probe_strip_masked_loci_from_set')} "
                f"--bed {clean_bed} "
                f"--twobit {base_twobit} "
                f"--output {strip_path} "
                f"--filter-mask {cfg.strip_filter_mask} "
                f"--min-length {cfg.strip_min_length}"
            )
            run_command(cmd, cwd=cfg.project_root, cfg=cfg, logger=logger)

    # 3) Build beds.conf listing all stripped BEDs for downstream phyluce
    conf_path = bed_dir / "beds.conf"
    logger.info(f"Writing beds.conf listing stripped BEDs to: {conf_path}")
    lines = ["[beds]\n"]
    for taxon in cfg.exemplar_taxa:
        if taxon == cfg.base_taxon:
            continue
        fname = f"{taxon}-to-{cfg.base_taxon}-MAPPING.bam.sort.merge.strip.bed"
        lines.append(f"{taxon}:{fname}\n")
    conf_path.write_text("".join(lines))
    return conf_path



def write_bed_conf(
    cfg: PipelineConfig,
    logger: logging.Logger,
) -> Path:
    """Write bed-files.conf listing *.sort.merge.strip.bed per taxon."""
    bed_dir = cfg.project_root / "bed"
    conf_path = bed_dir / "bed-files.conf"
    logger.info(f"Writing BED config: {conf_path}")
    lines = ["[beds]\n"]
    for taxon in cfg.exemplar_taxa:
        if taxon == cfg.base_taxon:
            continue
        fname = f"{taxon}-to-{cfg.base_taxon}-MAPPING.bam.sort.merge.strip.bed"
        lines.append(f"{taxon}:{fname}\n")
    conf_path.write_text("".join(lines))
    return conf_path


def multi_merge_table(
    cfg: PipelineConfig,
    logger: logging.Logger,
) -> Path:
    """Run phyluce_probe_get_multi_merge_table to build a sqlite DB.

    Idempotent: if the DB already exists, reuse it.
    """
    bed_dir = cfg.project_root / "bed"
    db_name = f"{cfg.base_taxon}-multi-merge.sqlite"
    db_path = bed_dir / db_name

    if db_path.exists():
        logger.info(
            "Multi-merge sqlite already exists, skipping rebuild: %s",
            db_path,
        )
        return db_path

    logger.info(f"Building multi-merge table: {db_path}")
    cmd = (
        f"{cfg.phyluce_cmd('phyluce_probe_get_multi_merge_table')} "
        f"--conf bed-files.conf "
        f"--base-taxon {cfg.base_taxon} "
        f"--output {db_name}"
    )
    run_command(cmd, cwd=bed_dir, cfg=cfg, logger=logger)
    return db_path



def query_multi_merge_table(
    cfg: PipelineConfig,
    db_path: Path,
    logger: logging.Logger,
) -> Path:
    """Query multi-merge table for loci shared across N taxa, output BED."""
    bed_dir = cfg.project_root / "bed"
    out_bed_name = f"{cfg.base_taxon}+{cfg.multi_merge_specific_counts}.bed"
    logger.info(
        f"Querying multi-merge table for loci shared by base + "
        f"{cfg.multi_merge_specific_counts} taxa -> {out_bed_name}"
    )
    cmd = (
        f"{cfg.phyluce_cmd('phyluce_probe_query_multi_merge_table')} "
        f"--db {db_path.name} "
        f"--base-taxon {cfg.base_taxon} "
        f"--output {out_bed_name} "
        f"--specific-counts {cfg.multi_merge_specific_counts}"
    )
    run_command(cmd, cwd=bed_dir, cfg=cfg, logger=logger)
    return bed_dir / out_bed_name


def preview_multi_merge_counts(
    cfg: PipelineConfig,
    db_path: Path,
    logger: logging.Logger,
) -> str:
    """Run phyluce_probe_query_multi_merge_table in summary mode and return its stdout.

    This runs the command without --output/--specific-counts so that PHYLUCE
    prints the distribution of loci shared by the base taxon + N taxa. The
    textual summary is also logged for the GUI.
    """
    bed_dir = cfg.project_root / "bed"
    cmd = (
        f"{cfg.phyluce_cmd('phyluce_probe_query_multi_merge_table')} "
        f"--db {db_path.name} "
        f"--base-taxon {cfg.base_taxon}"
    )
    logger.info("Previewing shared-locus counts (phyluce_probe_query_multi_merge_table)...")
    if cfg.dry_run:
        logger.info("Dry-run mode is enabled; not executing preview command.")
        return "DRY-RUN: preview not executed because dry_run=True."

    result = subprocess.run(
        cmd,
        shell=True,
        cwd=str(bed_dir),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if result.returncode != 0:
        logger.error(result.stdout)
        raise RuntimeError(
            f"Preview command failed with exit code {result.returncode}: {cmd}"
        )

    logger.info("Shared-locus summary from PHYLUCE:\n" + result.stdout.strip())
    return result.stdout


def extract_genome_sequences_from_bed(
    cfg: PipelineConfig,
    bed_path: Path,
    logger: logging.Logger,
) -> Path:
    """Extract FASTA sequences from the base genome for bait design."""
    bed_dir = cfg.project_root / "bed"
    base_twobit = (
        cfg.project_root / "genomes" / cfg.base_taxon / f"{cfg.base_taxon}.2bit"
    )
    out_fasta_name = f"{bed_path.stem}.fasta"
    logger.info(
        f"Extracting FASTA sequences from base genome ({base_twobit}) "
        f"for BED {bed_path.name} -> {out_fasta_name}"
    )
    cmd = (
        f"{cfg.phyluce_cmd('phyluce_probe_get_genome_sequences_from_bed')} "
        f"--bed {bed_path.name} "
        f"--twobit {base_twobit} "
        f"--buffer-to {cfg.genome_seq_buffer_to} "
        f"--output {out_fasta_name}"
    )
    run_command(cmd, cwd=bed_dir, cfg=cfg, logger=logger)
    return bed_dir / out_fasta_name

def design_temp_baits(
    cfg: PipelineConfig,
    fasta_path: Path,
    logger: logging.Logger,
) -> Path:
    """Design temporary baits from base genome FASTA using phyluce."""
    probe_dir = cfg.project_root / "probe-design"
    probe_dir.mkdir(exist_ok=True)

    # Use relative path from probe_dir to fasta_path
    rel_fasta = os.path.relpath(fasta_path, probe_dir)
    temp_probes_name = f"{fasta_path.stem}.temp.probes"
    temp_probes_path = probe_dir / temp_probes_name

    # === NEW: resume-friendly behaviour to avoid interactive overwrite ===
    if temp_probes_path.exists():
        logger.info(
            "Temporary probes already exist, re-using and "
            "skipping phyluce_probe_get_tiled_probes: %s",
            temp_probes_path,
        )
        return temp_probes_path
    # =====================================================================

    logger.info(
        "Designing temporary baits from %s -> %s (tiling-density=%s, masking=%s)",
        fasta_path,
        temp_probes_name,
        cfg.tiling_density,
        cfg.temp_masking,
    )

    cmd = (
        f"{cfg.phyluce_cmd('phyluce_probe_get_tiled_probes')} "
        f"--input {rel_fasta} "
        f'--probe-prefix "{cfg.probe_prefix}" '
        f"--design {cfg.design_name} "
        f"--designer {cfg.designer_name} "
        f"--tiling-density {cfg.tiling_density} "
        f"--two-probes "
        f"--overlap middle "
        f"--masking {cfg.temp_masking} "
        f"--remove-gc "
        f"--output {temp_probes_name}"
    )
    run_command(cmd, cwd=probe_dir, cfg=cfg, logger=logger)
    return temp_probes_path


def dupe_screen_temp_baits(cfg, temp_probes, logger):
    """
    Run lastz self-alignment on the temporary probes and remove duplicate hits.

    We avoid the 'space in path' bug by:
      - running phyluce_probe_easy_lastz in a temp directory under $HOME
        (which has no spaces), and
      - copying the resulting .lastz file back into the project probe-design dir.

    Parameters
    ----------
    cfg : Config
        Configuration object with project_root and phyluce_cmd().
    temp_probes : str or Path
        Path to the temporary probes FASTA (e.g. zannascaffolds+14.temp.probes).
    logger : logging.Logger
        Logger for status messages.

    Returns
    -------
    Path
        Path to the duplicate-screened temporary probes file.
    """
    project_root = Path(cfg.project_root)
    probe_dir = project_root / "probe-design"

    temp_probes = Path(temp_probes)
    if not temp_probes.is_absolute():
        temp_probes = probe_dir / temp_probes

    if not temp_probes.exists():
        raise FileNotFoundError(f"Temporary probes file not found: {temp_probes}")

    # Where the final lastz output should live in your project
    lastz_out = probe_dir / f"{temp_probes.name}-TO-SELF-PROBES.lastz"

    # --- 1) Run lastz in a temp directory under $HOME (no spaces) ---
    # This avoids the '/media/localuser/T7 Shield' space issue entirely.
    home_tmp_root = Path.home() / "phyluce_tmp"
    tmp_probe_dir = home_tmp_root / "probe-design"
    tmp_probe_dir.mkdir(parents=True, exist_ok=True)

    tmp_probes = tmp_probe_dir / temp_probes.name
    tmp_lastz = tmp_probe_dir / f"{temp_probes.name}-TO-SELF-PROBES.lastz"

    # If lastz output already exists in the project, skip re-running lastz entirely
    if lastz_out.exists():
        logger.info(f"Lastz self-alignments already exist, skipping phyluce_probe_easy_lastz: {lastz_out}")
    else:
        # Ensure the probes file exists in the temp dir (symlink if possible, otherwise copy)
        if not tmp_probes.exists():
            try:
                if tmp_probes.exists():
                    tmp_probes.unlink()
                tmp_probes.symlink_to(temp_probes)
                logger.info(f"Symlinked probes to temp dir for lastz: {tmp_probes} -> {temp_probes}")
            except OSError:
                shutil.copy2(temp_probes, tmp_probes)
                logger.info(f"Copied probes to temp dir for lastz: {tmp_probes}")

        # Build easy_lastz command using ONLY the filename, not an absolute path
        cmd_lastz = (
            f"{cfg.phyluce_cmd('phyluce_probe_easy_lastz')} "
            f"--target {tmp_probes.name} "
            f"--query {tmp_probes.name} "
            f"--identity 50 --coverage 50 "
            f"--output {tmp_lastz.name}"
        )

        logger.info(
            f"Running lastz self-alignment in temp dir without spaces: {tmp_probe_dir}\n"
            f"  Command: {cmd_lastz}"
        )
        # Run in the temp directory so paths are clean and space-free
        run_command(cmd_lastz, cwd=tmp_probe_dir, cfg=cfg, logger=logger)

        if not tmp_lastz.exists():
            raise FileNotFoundError(
                f"Expected lastz output not found in temp dir: {tmp_lastz}"
            )

        # Copy the lastz output back to the real project probe-design directory
        shutil.copy2(tmp_lastz, lastz_out)
        logger.info(f"Copied lastz output back to project dir: {lastz_out}")

    # --- 2) Run the duplicate-removal step in the REAL probe_dir ---
    # Name of the duplicate-screened probes file – you can tweak the suffix if you like.
    dedup_probes = probe_dir / f"{temp_probes.stem}.no-dupes.probes"

    if dedup_probes.exists():
        logger.info(f"Duplicate-screened temp probes already exist, skipping: {dedup_probes}")
    else:
        cmd_rmdupes = (
            f"{cfg.phyluce_cmd('phyluce_probe_remove_duplicate_hits_from_probes_using_lastz')} "
            f"--fasta {temp_probes.name} "
            f"--lastz {lastz_out.name} "
            f"--probe-prefix uce- "
            # f"--output {dedup_probes.name}"
            # f"> {dedup_probes.name}"
        )
        logger.info(f"Removing duplicate hits from temp probes using lastz: {cmd_rmdupes}")
        run_command(cmd_rmdupes, cwd=probe_dir, cfg=cfg, logger=logger)


    return dedup_probes



# -------------------- High-level pipeline orchestration -------------------- #

def run_full_pipeline(cfg: PipelineConfig, logger: logging.Logger) -> None:
    logger.info("==== PHYLUC-E Tutorial IV pipeline starting ====")
    logger.info(f"Genomes root: {cfg.genomes_root}")
    logger.info(f"Project root: {cfg.project_root}")
    logger.info(f"Base taxon: {cfg.base_taxon}")
    logger.info(f"Exemplar taxa: {', '.join(cfg.exemplar_taxa)}")
    logger.info(f"Threads: {cfg.threads} (stampy), dry_run={cfg.dry_run}")

    genomes = discover_genomes(cfg.genomes_root, logger)
    if cfg.base_taxon not in genomes:
        raise ValueError(
            f"Base taxon '{cfg.base_taxon}' not found among discovered genomes: "
            f"{', '.join(sorted(genomes))}"
        )

    # Restrict to taxa we actually want to process
    missing = [t for t in cfg.exemplar_taxa if t not in genomes]
    if missing:
        raise ValueError(f"Exemplar taxa not found among genomes: {', '.join(missing)}")

    subset_genomes = {t: genomes[t] for t in cfg.exemplar_taxa}
    # Ensure base taxon is included as well
    subset_genomes[cfg.base_taxon] = genomes[cfg.base_taxon]

    project_genomes = setup_project_structure(cfg, subset_genomes, logger)
    convert_genomes_to_2bit(cfg, project_genomes, logger)
    simulate_reads(cfg, project_genomes, logger)
    prepare_base_genome(cfg, project_genomes, logger)
    align_reads_with_stampy(cfg, logger)
    filter_bams_and_link(cfg, logger)
    bams_to_beds(cfg, logger)
    sort_beds(cfg, logger)
    merge_beds(cfg, logger)
    strip_masked_loci(cfg, logger)
    write_bed_conf(cfg, logger)
    db_path = multi_merge_table(cfg, logger)
    bed_path = query_multi_merge_table(cfg, db_path, logger)
    fasta_path = extract_genome_sequences_from_bed(cfg, bed_path, logger)
    temp_probes = design_temp_baits(cfg, fasta_path, logger)
    dupe_screen_temp_baits(cfg, temp_probes, logger)

    logger.info(f"Temp DUPE-SCREENED probes written to: {dupe_temp}")
    logger.info("==== PHYLUC-E Tutorial IV pipeline completed ====")
    return dupe_temp

# -------------------- Tkinter GUI -------------------- #

class TextHandler(logging.Handler):
    """Logging handler that writes to a Tkinter Text widget."""
    def __init__(self, text_widget: "tk.Text"):
        super().__init__()
        self.text_widget = text_widget

    def emit(self, record: logging.LogRecord) -> None:
        msg = self.format(record)
        self.text_widget.after(0, self._append, msg)

    def _append(self, msg: str) -> None:
        self.text_widget.insert(tk.END, msg + "\n")
        self.text_widget.see(tk.END)


def launch_gui() -> None:
    if not TK_AVAILABLE:
        print("Tkinter is not available; run with --no-gui and CLI arguments instead.")
        sys.exit(1)

    root = tk.Tk()
    root.title("PHYLUC-E GUI wrapper (Tutorial IV)")

    # Main frames
    frm_paths = ttk.LabelFrame(root, text="Paths")
    frm_paths.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)
    frm_taxa = ttk.LabelFrame(root, text="Taxa")
    frm_taxa.grid(row=1, column=0, sticky="nsew", padx=5, pady=5)
    frm_params = ttk.LabelFrame(root, text="Key parameters")
    frm_params.grid(row=2, column=0, sticky="nsew", padx=5, pady=5)
    frm_buttons = ttk.Frame(root)
    frm_buttons.grid(row=3, column=0, sticky="ew", padx=5, pady=5)
    frm_log = ttk.LabelFrame(root, text="Log")
    frm_log.grid(row=0, column=1, rowspan=4, sticky="nsew", padx=5, pady=5)

    root.columnconfigure(0, weight=0)
    root.columnconfigure(1, weight=1)
    root.rowconfigure(0, weight=0)
    root.rowconfigure(1, weight=0)
    root.rowconfigure(2, weight=0)
    root.rowconfigure(3, weight=0)

    # Paths frame
    tk.Label(frm_paths, text="Genomes root:").grid(row=0, column=0, sticky="w")
    entry_genomes = tk.Entry(frm_paths, width=40)
    entry_genomes.grid(row=0, column=1, sticky="ew", padx=2)
    def browse_genomes():
        d = filedialog.askdirectory(title="Select genomes root")
        if d:
            entry_genomes.delete(0, tk.END)
            entry_genomes.insert(0, d)
    ttk.Button(frm_paths, text="Browse", command=browse_genomes).grid(row=0, column=2, padx=2)

    tk.Label(frm_paths, text="Project root:").grid(row=1, column=0, sticky="w")
    entry_project = tk.Entry(frm_paths, width=40)
    entry_project.grid(row=1, column=1, sticky="ew", padx=2)
    def browse_project():
        d = filedialog.askdirectory(title="Select project root (will be created if needed)")
        if d:
            entry_project.delete(0, tk.END)
            entry_project.insert(0, d)
    ttk.Button(frm_paths, text="Browse", command=browse_project).grid(row=1, column=2, padx=2)

    frm_paths.columnconfigure(1, weight=1)

    # Taxa frame
    btn_scan = ttk.Button(frm_taxa, text="Scan genomes")
    btn_scan.grid(row=0, column=0, columnspan=3, pady=2)

    tk.Label(frm_taxa, text="Base taxon:").grid(row=1, column=0, sticky="w")
    combo_base = ttk.Combobox(frm_taxa, values=[], state="readonly", width=20)
    combo_base.grid(row=1, column=1, sticky="ew", padx=2)

    tk.Label(frm_taxa, text="Exemplar taxa:").grid(row=2, column=0, sticky="nw")
    list_exemplars = tk.Listbox(frm_taxa, selectmode=tk.MULTIPLE, height=8, exportselection=False)
    list_exemplars.grid(row=2, column=1, columnspan=2, sticky="nsew", padx=2, pady=2)
    frm_taxa.columnconfigure(1, weight=1)
    frm_taxa.rowconfigure(2, weight=1)

    # Parameters frame
    tk.Label(frm_params, text="Threads (stampy):").grid(row=0, column=0, sticky="w")
    entry_threads = tk.Entry(frm_params, width=6)
    entry_threads.insert(0, "8")
    entry_threads.grid(row=0, column=1, sticky="w", padx=2)

    tk.Label(frm_params, text="Min shared taxa (specific-counts):").grid(row=1, column=0, sticky="w")
    entry_mm = tk.Entry(frm_params, width=6)
    entry_mm.insert(0, "5")
    entry_mm.grid(row=1, column=1, sticky="w", padx=2)

    tk.Label(frm_params, text="BED merge distance (-d):").grid(row=2, column=0, sticky="w")
    entry_merge_d = tk.Entry(frm_params, width=6)
    entry_merge_d.insert(0, "0")
    entry_merge_d.grid(row=2, column=1, sticky="w", padx=2)

    tk.Label(frm_params, text="Min coverage (contigs<->probes):").grid(row=3, column=0, sticky="w")
    entry_min_cov = tk.Entry(frm_params, width=6)
    entry_min_cov.insert(0, "67")
    entry_min_cov.grid(row=3, column=1, sticky="w", padx=2)

    tk.Label(frm_params, text="Stampy substitution rate:").grid(row=4, column=0, sticky="w")
    entry_sub_rate = tk.Entry(frm_params, width=6)
    entry_sub_rate.insert(0, "0.05")
    entry_sub_rate.grid(row=4, column=1, sticky="w", padx=2)

    tk.Label(frm_params, text="Stampy insert size:").grid(row=5, column=0, sticky="w")
    entry_ins_size = tk.Entry(frm_params, width=6)
    entry_ins_size.insert(0, "400")
    entry_ins_size.grid(row=5, column=1, sticky="w", padx=2)

    tk.Label(frm_params, text="Stampy extra options:").grid(row=6, column=0, sticky="w")
    entry_extra_opts = tk.Entry(frm_params, width=20)
    entry_extra_opts.insert(0, "")
    entry_extra_opts.grid(row=6, column=1, sticky="w", padx=2)

    dry_var = tk.BooleanVar(value=False)
    chk_dry = ttk.Checkbutton(frm_params, text="Dry-run (log commands only)", variable=dry_var)
    chk_dry.grid(row=7, column=0, columnspan=2, sticky="w", pady=2)

    # Log frame
    txt_log = tk.Text(frm_log, wrap="word", height=25)
    txt_log.grid(row=0, column=0, sticky="nsew")
    scroll_log = ttk.Scrollbar(frm_log, orient="vertical", command=txt_log.yview)
    scroll_log.grid(row=0, column=1, sticky="ns")
    txt_log.configure(yscrollcommand=scroll_log.set)
    frm_log.rowconfigure(0, weight=1)
    frm_log.columnconfigure(0, weight=1)

    # Logger for GUI
    logger = setup_logger(log_path=None)
    # Replace handlers with TextHandler
    logger.handlers = []
    text_handler = TextHandler(txt_log)
    text_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.addHandler(text_handler)

    def do_scan_genomes():
        root_dir = entry_genomes.get().strip()
        if not root_dir:
            messagebox.showerror("Error", "Please select a genomes root directory first.")
            return
        try:
            genomes = discover_genomes(Path(root_dir), logger)
        except Exception as e:
            messagebox.showerror("Error scanning genomes", str(e))
            return
        taxa = sorted(genomes.keys())
        combo_base["values"] = taxa
        if taxa:
            combo_base.set(taxa[0])
        list_exemplars.delete(0, tk.END)
        for t in taxa:
            list_exemplars.insert(tk.END, t)
        # Pre-select all taxa as exemplars
        for i in range(len(taxa)):
            list_exemplars.select_set(i)
        logger.info(f"Scanned genomes: {', '.join(taxa)}")

    btn_scan.configure(command=do_scan_genomes)

    def prompt_specific_counts(counts_text: str, default_value: int) -> int:
        """Small dialog to let the user choose specific-counts based on PHYLUCE summary."""
        # Store chosen value in a mutable container so we can modify it in nested callbacks
        chosen = {"value": default_value}

        top = tk.Toplevel(root)
        top.title("Choose specific-counts (number of taxa sharing each locus)")

        tk.Label(
            top,
            text=(
                "PHYLUCE summary of loci shared by the base taxon + N taxa.\n"
                "Use this to choose the 'specific-counts' threshold.\n"
                "Higher values = loci shared across more exemplar genomes (more conserved, fewer loci)."
            ),
            justify="left",
            wraplength=600,
        ).grid(row=0, column=0, columnspan=2, sticky="w", padx=5, pady=5)

        txt = tk.Text(top, width=80, height=10)
        txt.grid(row=1, column=0, columnspan=2, padx=5, pady=5, sticky="nsew")
        txt.insert("1.0", counts_text.strip() or "(no output captured)")
        txt.configure(state="disabled")

        tk.Label(top, text="specific-counts to use:").grid(row=2, column=0, sticky="e", padx=5, pady=5)
        entry = tk.Entry(top, width=6)
        entry.grid(row=2, column=1, sticky="w", padx=5, pady=5)
        entry.insert(0, str(default_value))

        def on_ok():
            val = entry.get().strip()
            try:
                iv = int(val)
                if iv < 0:
                    raise ValueError
            except ValueError:
                messagebox.showerror(
                    "Invalid input",
                    "specific-counts must be a non-negative integer."
                )
                return
            chosen["value"] = iv
            top.destroy()

        def on_cancel():
            top.destroy()

        btn_frame = ttk.Frame(top)
        btn_frame.grid(row=3, column=0, columnspan=2, pady=5)
        ttk.Button(btn_frame, text="Use this value", command=on_ok).grid(row=0, column=0, padx=5)
        ttk.Button(btn_frame, text="Cancel (keep current)", command=on_cancel).grid(row=0, column=1, padx=5)

        # Make dialog modal
        top.transient(root)
        top.grab_set()
        root.wait_window(top)
        return chosen["value"]

    def run_pipeline_with_specific_counts_prompt(cfg: PipelineConfig) -> None:
        """Run the full pipeline, but pause to let the user choose specific-counts.

        This mirrors run_full_pipeline(), but between building the multi-merge
        table and querying it, we:
          1. Run phyluce_probe_query_multi_merge_table in summary mode
          2. Show the summary to the user
          3. Let them adjust specific-counts interactively
        """
        logger.info("==== PHYLUC-E Tutorial IV pipeline (GUI) starting ====")
        logger.info(f"Genomes root: {cfg.genomes_root}")
        logger.info(f"Project root: {cfg.project_root}")
        logger.info(f"Base taxon: {cfg.base_taxon}")
        logger.info(f"Exemplar taxa: {', '.join(cfg.exemplar_taxa)}")
        logger.info(f"Threads: {cfg.threads} (stampy), dry_run={cfg.dry_run}")

        genomes = discover_genomes(cfg.genomes_root, logger)
        if cfg.base_taxon not in genomes:
            raise ValueError(
                f"Base taxon '{cfg.base_taxon}' not found among discovered genomes: "
                f"{', '.join(sorted(genomes))}"
            )

        # Restrict to taxa we actually want to process
        missing = [t for t in cfg.exemplar_taxa if t not in genomes]
        if missing:
            raise ValueError(f"Exemplar taxa not found among genomes: {', '.join(missing)}")

        subset_genomes = {t: genomes[t] for t in cfg.exemplar_taxa}
        # Ensure base taxon is included as well
        subset_genomes[cfg.base_taxon] = genomes[cfg.base_taxon]

        project_genomes = setup_project_structure(cfg, subset_genomes, logger)
        convert_genomes_to_2bit(cfg, project_genomes, logger)
        simulate_reads(cfg, project_genomes, logger)
        prepare_base_genome(cfg, project_genomes, logger)
        align_reads_with_stampy(cfg, logger)
        filter_bams_and_link(cfg, logger)
        bams_to_beds(cfg, logger)
        sort_beds(cfg, logger)
        merge_beds(cfg, logger)
        strip_masked_loci(cfg, logger)
        write_bed_conf(cfg, logger)
        db_path = multi_merge_table(cfg, logger)

        # --- NEW: preview shared-locus counts and prompt for specific-counts ---
        try:
            summary = preview_multi_merge_counts(cfg, db_path, logger)
        except Exception as e:
            logger.error(f"Failed to preview shared-locus counts: {e}")
            raise

        # Use the current GUI entry as default, falling back on cfg multi-merge
        try:
            default_mm = int(entry_mm.get().strip())
        except ValueError:
            default_mm = cfg.multi_merge_specific_counts

        chosen_mm = prompt_specific_counts(summary, default_mm)
        cfg.multi_merge_specific_counts = chosen_mm
        entry_mm.delete(0, tk.END)
        entry_mm.insert(0, str(chosen_mm))
        logger.info(f"Using specific-counts={chosen_mm} for downstream BED and temp bait design.")

        # Continue as in run_full_pipeline
        bed_path = query_multi_merge_table(cfg, db_path, logger)
        fasta_path = extract_genome_sequences_from_bed(cfg, bed_path, logger)
        temp_probes = design_temp_baits(cfg, fasta_path, logger)
        dupe_screen_temp_baits(cfg, temp_probes, logger)

        logger.info("==== PHYLUC-E Tutorial IV pipeline (GUI) completed ====")

    def run_gui_pipeline():
        genomes_root = entry_genomes.get().strip()
        project_root = entry_project.get().strip()
        base_taxon = combo_base.get().strip()
        if not genomes_root or not project_root or not base_taxon:
            messagebox.showerror(
                "Missing input",
                "Please set genomes root, project root, and select a base taxon."
            )
            return
        try:
            threads = int(entry_threads.get().strip())
            mm_count = int(entry_mm.get().strip())
            merge_d = int(entry_merge_d.get().strip())
            min_cov = int(entry_min_cov.get().strip())
            sub_rate = float(entry_sub_rate.get().strip())
            ins_size = int(entry_ins_size.get().strip())
        except ValueError:
            messagebox.showerror("Invalid input", "Numeric parameters must be valid numbers.")
            return

        extra_opts = entry_extra_opts.get().strip()

        selected_indices = list_exemplars.curselection()
        exemplar_taxa = [list_exemplars.get(i) for i in selected_indices]
        if base_taxon not in exemplar_taxa:
            exemplar_taxa.append(base_taxon)

        cfg = PipelineConfig(
            genomes_root=Path(genomes_root),
            project_root=Path(project_root),
            base_taxon=base_taxon,
            exemplar_taxa=exemplar_taxa,
            threads=threads,
            multi_merge_specific_counts=mm_count,
            bedtools_merge_distance=merge_d,
            min_coverage=min_cov,
            stampy_substitution_rate=sub_rate,
            stampy_insert_size=ins_size,
            stampy_extra_opts=extra_opts,
            dry_run=dry_var.get(),
        )

        try:
            run_pipeline_with_specific_counts_prompt(cfg)
            messagebox.showinfo("Done", "Pipeline completed. Check the log for details.")
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            messagebox.showerror("Pipeline failed", str(e))

    ttk.Button(frm_buttons, text="Run pipeline", command=run_gui_pipeline).grid(
        row=0, column=0, padx=5
    )
    ttk.Button(frm_buttons, text="Quit", command=root.destroy).grid(
        row=0, column=1, padx=5
    )

    root.mainloop()


# -------------------- CLI entrypoint -------------------- #

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="GUI + CLI wrapper around parts of phyluce Tutorial IV."
    )
    p.add_argument(
        "--no-gui",
        action="store_true",
        help="Run in command-line mode only (no Tkinter GUI).",
    )
    p.add_argument(
        "--genomes-root",
        type=Path,
        help="Root directory containing cleaned genomes (FASTA files).",
    )
    p.add_argument(
        "--project-root",
        type=Path,
        help="Project output directory (will be created if needed).",
    )
    p.add_argument(
        "--base-taxon",
        type=str,
        help="Short, lowercase taxon name to use as base genome (e.g., tricas1).",
    )
    p.add_argument(
        "--exemplar-taxa",
        type=str,
        help="Comma-separated list of exemplar taxa (short names).",
    )
    p.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of threads/cores to give stampy.",
    )
    p.add_argument(
        "--multi-merge-specific-counts",
        type=int,
        default=5,
        help="Required number of taxa sharing a locus in multi-merge table.",
    )
    p.add_argument(
        "--bed-merge-distance",
        type=int,
        default=0,
        help="bedtools merge -d distance (0 = only overlapping).",
    )
    p.add_argument(
        "--min-coverage",
        type=int,
        default=67,
        help="Min coverage threshold for contigs<->probes (used later).",
    )
    p.add_argument(
        "--stampy-substitution-rate",
        type=float,
        default=0.05,
        help="Stampy substitution rate (default: 0.05).",
    )
    p.add_argument(
        "--stampy-insert-size",
        type=int,
        default=400,
        help="Stampy insert size (default: 400).",
    )
    p.add_argument(
        "--stampy-extra-opts",
        type=str,
        default="",
        help="Additional options passed verbatim to stampy (e.g. '--sensitive').",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print and log commands without executing them.",
    )
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)

    if not args.no_gui:
        # If Tkinter is available and user did not force CLI, launch the GUI
        if TK_AVAILABLE:
            launch_gui()
            return
        else:
            print("Tkinter is not available, falling back to CLI mode.", file=sys.stderr)

    # CLI mode
    if not (args.genomes_root and args.project_root and args.base_taxon and args.exemplar_taxa):
        print(
            "In CLI mode you must provide --genomes-root, --project-root, "
            "--base-taxon, and --exemplar-taxa.",
            file=sys.stderr,
        )
        sys.exit(1)

    exemplar_taxa = [t.strip() for t in args.exemplar_taxa.split(",") if t.strip()]

    log_path = Path(args.project_root) / "phyluce_one_passV2.log"
    logger = setup_logger(log_path=log_path)

    cfg = PipelineConfig(
        genomes_root=args.genomes_root,
        project_root=args.project_root,
        base_taxon=args.base_taxon,
        exemplar_taxa=exemplar_taxa,
        threads=args.threads,
        multi_merge_specific_counts=args.multi_merge_specific_counts,
        bedtools_merge_distance=args.bed_merge_distance,
        min_coverage=args.min_coverage,
        stampy_substitution_rate=args.stampy_substitution_rate,
        stampy_insert_size=args.stampy_insert_size,
        stampy_extra_opts=args.stampy_extra_opts,
        dry_run=args.dry_run,
    )

    try:
        run_full_pipeline(cfg, logger)
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
