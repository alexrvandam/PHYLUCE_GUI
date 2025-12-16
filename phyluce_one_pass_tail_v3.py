#!/usr/bin/env python3
"""
phyluce_one_pass_tail.py

Tail pipeline for PHYLUCE Tutorial IV:
From duplicate-free temporary probes to final multi-genome master bait list.

Intended to be run AFTER your main phyluce_one_passV2 pipeline has produced a
duplicate-free temporary probe file (e.g., triCas1+5.temp-DUPE-SCREENED.probes).

Usage example:

    python phyluce_one_pass_tail.py \
        --project-root /path/to/uce-coleoptera \
        --genomes-root /path/to/uce-coleoptera/genomes \
        --base-taxon triCas1 \
        --temp-probes bed/triCas1+5.temp-DUPE-SCREENED.probes \
        --taxa agrPla1 anoGla1 denPon1 lepDec1 ontTau1 triCas1 menMol1 \
        --multi-fasta-specific-counts 4 \
        --cores 16

This will create a 'probe-design' directory under project-root containing:
    - genome-lastz/
    - genome-fasta/
    - multifastas.sqlite
    - <base>+N-back-to-X.conf
    - <design-name>-master-probe-list.fasta
    - <design-name>-master-probe-list-DUPE-SCREENED.fasta
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional


def run_command(cmd, cwd=None, logger=None):
    """
    Wrapper around subprocess.run that logs the command and captures stdout/stderr.

    NOTE: We used to special-case phyluce_probe_run_multiple_lastzs_sqlite here
    and proactively delete its --output directory to avoid interactive prompts.
    Step 2 now uses our own wrapper script instead, so we no longer touch the
    LASTZ output directory here.
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    logger.info("Running: %s", " ".join(cmd))
    result = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd is not None else None,
        capture_output=True,
        text=True,
    )

    if result.stdout:
        logger.info(result.stdout)

    if result.returncode != 0:
        logger.error("Command failed with exit code %s", result.returncode)
        logger.error("STDOUT:\n%s", result.stdout)
        logger.error("STDERR:\n%s", result.stderr)
        raise subprocess.CalledProcessError(
            result.returncode,
            cmd,
            output=result.stdout,
            stderr=result.stderr,
        )

    return result

def safe_cwd_for_lastz(path: Path) -> Path:
    """
    Work around LASTZ breaking on absolute paths containing spaces.

    If the given path contains spaces, try a version where spaces are replaced
    by underscores (e.g. '/media/.../T7 Shield' -> '/media/.../T7_Shield').
    If that path exists (i.e. you created a symlink), use it as cwd;
    otherwise fall back to the original path.
    """
    s = str(path)
    if " " not in s:
        return path

    candidate = Path(s.replace(" ", "_"))
    if candidate.exists():
        return candidate

    # No matching symlink – just return the original and hope for the best
    return path


def phyluce_cmd(cfg: "TailConfig", script: str) -> List[str]:
    """
    Build the base command list for calling a PHYLUCE CLI script.

    - cfg.phyluce_prefix is usually something like 'conda run -n phyluce-1.7.3'
    - `script` is the PHYLUCE script name, e.g. 'phyluce_probe_slice_sequence_from_genomes'.

    We return a list, e.g.:
        ['conda', 'run', '-n', 'phyluce-1.7.3', 'phyluce_probe_slice_sequence_from_genomes']
    or just ['phyluce_probe_slice_sequence_from_genomes'] if no prefix is set.
    """
    prefix = getattr(cfg, "phyluce_prefix", "") or ""
    prefix = prefix.strip()

    if prefix:
        return prefix.split() + [script]
    else:
        return [script]


def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


@dataclass
class TailConfig:
    project_root: Path
    genomes_root: Path
    base_taxon: str
    temp_probes: Path

    taxa: List[str] = field(default_factory=list)

    # phyluce / lastz parameters
    lastz_identity: int = 50
    cores: int = 4

    # Optional prefix to run PHYLUCE tools via another env, e.g.
    # "conda run -n phyluce-1.7.3"
    phyluce_prefix: str = "conda run -n phyluce-1.7.3"

    # slice-from-genomes parameters
    slice_buffer: int = 180

    # multi-fasta decisions
    multi_fasta_specific_counts: Optional[int] = None  # if None: auto (n-2)

    # final probe design parameters
    probe_prefix: str = "uce-"
    designer: str = "user"
    design_name: str = "uce-design-v1"
    tiling_density: int = 3
    masking_threshold: float = 0.25
    overlap_mode: str = "middle"

    # logging
    log_level: str = "INFO"

    # SQLite DB name for step 2 (can be set by GUI)
    sqlite_dbname: str = ""

    # computed paths
    probe_design_dir: Optional[Path] = None
    genome_lastz_dir: Optional[Path] = None
    genome_fasta_dir: Optional[Path] = None


def setup_paths(cfg: TailConfig) -> None:
    """Initialize probe-design paths under project_root."""
    cfg.project_root = cfg.project_root.resolve()
    cfg.genomes_root = cfg.genomes_root.resolve()
    cfg.temp_probes = cfg.temp_probes.resolve()

    probe_design_dir = ensure_dir(cfg.project_root / "probe-design")
    cfg.probe_design_dir = probe_design_dir
    cfg.genome_lastz_dir = ensure_dir(probe_design_dir / "genome-lastz")
    cfg.genome_fasta_dir = ensure_dir(probe_design_dir / "genome-fasta")


def ensure_probe_design_layout(cfg: TailConfig) -> None:
    """
    Ensure that the probe-design directory structure exists and that
    the corresponding fields on `cfg` are populated.

    This is a small helper for older code that expects a function
    named `ensure_probe_design_layout`.
    """
    if (
        cfg.probe_design_dir is None
        or cfg.genome_lastz_dir is None
        or cfg.genome_fasta_dir is None
    ):
        setup_paths(cfg)


def align_temp_probes_to_genomes(cfg: TailConfig, logger: logging.Logger) -> Path:
    """
    Step 2: Align temp probes to genomes using a resumable LASTZ wrapper.

    - Uses phyluce_probe_run_multiple_lastzs_sqlite_resumable.py
      (must live next to this tail script).
    - Does NOT delete/rename existing LASTZ or DB results; the wrapper
      is responsible for skipping already-completed taxa when --resume is set.
    """
    ensure_probe_design_layout(cfg)

    if not cfg.taxa:
        raise ValueError("No genomes selected for Step 2; taxa list is empty.")

    assert cfg.probe_design_dir is not None

    # Decide DB name if GUI hasn’t already set it
    if not getattr(cfg, "sqlite_dbname", None):
        n_taxa = len(cfg.taxa)
        cfg.sqlite_dbname = f"{cfg.base_taxon}+{n_taxa}.sqlite"

    db_path = cfg.probe_design_dir / cfg.sqlite_dbname

    logger.info("==== STEP 2: Starting temp probes → master bait set ====")
    logger.info("Project root: %s", cfg.project_root)
    logger.info("Genomes root: %s", cfg.genomes_root)
    logger.info("Base taxon:   %s", cfg.base_taxon)
    logger.info("Temp probes:  %s", cfg.temp_probes)
    logger.info("DB file:      %s", db_path)

    if cfg.genome_lastz_dir is None:
        cfg.genome_lastz_dir = cfg.probe_design_dir / "genome-lastz"
    logger.info("LASTZ outdir: %s", cfg.genome_lastz_dir)

    # Run LASTZ for *all* taxa, including the base taxon.
    # This ensures slice_sequence_from_genomes finds a .lastz.clean file
    # for every entry in genomes.conf (including the base).
    taxa = sorted(cfg.taxa)
    logger.info("Scanned genomes (Step 2): %s", ", ".join(taxa))

    #taxa = sorted(t for t in cfg.taxa if t != cfg.base_taxon)
    #logger.info("Scanned genomes (Step 2): %s", ", ".join(taxa))

    if not cfg.temp_probes.is_file():
        raise FileNotFoundError(f"Temp probes file not found: {cfg.temp_probes}")

    cfg.genome_lastz_dir.mkdir(parents=True, exist_ok=True)

    # Wrapper script lives next to this tail script
    wrapper_script = Path(__file__).resolve().with_name(
        "phyluce_probe_run_multiple_lastzs_sqlite_resumable.py"
    )
    if not wrapper_script.exists():
        logger.error("Resumable LASTZ wrapper not found at %s", wrapper_script)
        raise FileNotFoundError(
            "phyluce_probe_run_multiple_lastzs_sqlite_resumable.py must be in the same "
            "directory as phyluce_one_pass_tail_v2.py"
        )

    # Build the command. We hard-code the env name here to avoid extra moving parts.
    cmd = [
        "conda",
        "run",
        "-n",
        "phyluce-1.7.3",
        str(wrapper_script),
        "--db",
        str(db_path),
        "--output",
        str(cfg.genome_lastz_dir),
        "--probefile",
        str(cfg.temp_probes),
        "--scaffoldlist",
        *taxa,
        "--genome-base-path",
        str(cfg.genomes_root),
        "--identity",
        str(cfg.lastz_identity),
        "--cores",
        str(cfg.cores),
        "--verbosity",
        "INFO",
        "--resume",
        "--log-path",
        str(cfg.probe_design_dir),
    ]

    run_command(cmd, cwd=cfg.probe_design_dir, logger=logger)
    return db_path


def write_scaffold_conf(cfg: TailConfig, logger: logging.Logger) -> Path:
    """
    Write a [scaffolds] config file pointing to <taxon>.2bit files.

    Uses paths *relative* to probe-design so that any higher-level
    directories with spaces (e.g. "/media/localuser/T7 Shield") never
    appear in genomes.conf.
    """
    assert cfg.probe_design_dir is not None

    conf_path = cfg.probe_design_dir / "genomes.conf"

    try:
        rel_genomes_root = os.path.relpath(cfg.genomes_root, cfg.probe_design_dir)
    except ValueError:
        rel_genomes_root = str(cfg.genomes_root)

    lines = ["[scaffolds]\n"]
    for taxon in cfg.taxa:
        twobit_rel = Path(rel_genomes_root) / taxon / f"{taxon}.2bit"
        lines.append(f"{taxon}:{twobit_rel.as_posix()}\n")

    logger.info("Writing scaffold config to %s", conf_path)
    conf_path.write_text("".join(lines))
    return conf_path



def slice_sequence_from_genomes(
    cfg: TailConfig, logger: logging.Logger, conf_path: Path
) -> None:
    """
    Step 2b: Slice probe-matching sequences from genomes.

    Wraps `phyluce_probe_slice_sequence_from_genomes`.

    IMPORTANT:
      - PHYLUCE will *interactively* ask:
          "[WARNING] Output directory exists, REMOVE [y/n]"
        if the --output directory already exists.
      - In a non-interactive GUI this causes an EOFError.
      - Here we pre-empt that by backing up any existing
        `genome-fasta` directory to a *.bak directory before
        calling PHYLUCE, so that it never needs to prompt.
    """

    assert cfg.probe_design_dir is not None
    assert cfg.genome_lastz_dir is not None

    # Ensure we have a path for the genome-fasta directory
    if cfg.genome_fasta_dir is None:
        cfg.genome_fasta_dir = cfg.probe_design_dir / "genome-fasta"

    # ------------------------------------------------------------------
    # Pre-flight: if the genome-fasta output directory already exists,
    # rename it to a backup so PHYLUCE sees a *fresh* non-existent
    # directory and does NOT trigger the interactive prompt.
    # ------------------------------------------------------------------
    outdir = cfg.probe_design_dir / cfg.genome_fasta_dir.name
    if outdir.exists():
        backup_dir = outdir.with_name(outdir.name + ".bak")
        i = 1
        while backup_dir.exists():
            backup_dir = outdir.with_name(outdir.name + f".bak{i}")
            i += 1

        logger.info(
            "Backing up existing genome-fasta directory before slicing: %s -> %s",
            outdir,
            backup_dir,
        )
        outdir.rename(backup_dir)

    # ------------------------------------------------------------------
    # Build the PHYLUCE command
    #   phyluce_probe_slice_sequence_from_genomes
    #      --conf genomes.conf
    #      --lastz genome-lastz
    #      --probes <buffer>
    #      --name-pattern <probes>_v_{}.lastz.clean
    #      --output genome-fasta
    #
    # NOTE:
    #   - We only pass *relative* names (genome-lastz, genome-fasta,
    #     genomes.conf) and run from cwd = probe-design/, so we do not
    #     expose any spaces in "/media/localuser/T7 Shield/...".
    # ------------------------------------------------------------------
    probes_basename = cfg.temp_probes.name
    name_pattern = f"{probes_basename}_v_{{}}.lastz.clean"

    cmd = phyluce_cmd(cfg, "phyluce_probe_slice_sequence_from_genomes") + [
        "--conf",
        str(conf_path.name),
        "--lastz",
        str(cfg.genome_lastz_dir.name),
        "--probes",
        str(cfg.slice_buffer),
        "--name-pattern",
        name_pattern,
        "--output",
        str(cfg.genome_fasta_dir.name),
    ]

    logger.info(
        "Running: %s",
        " ".join(str(c) for c in cmd),
    )

    # Run from within probe-design so relative paths line up
    run_command(cmd, cwd=cfg.probe_design_dir, logger=logger)


def build_multi_fasta_table(cfg: TailConfig, logger: logging.Logger) -> Path:
    """
    Run phyluce_probe_get_multi_fasta_table to build multifastas.sqlite.

    This will:
      - use the genome-fasta directory as input
      - write multifastas.sqlite in the probe_design_dir
      - remove any existing multifastas.sqlite first to avoid the
        interactive [y/n] overwrite prompt that crashes in a
        non-interactive subprocess.
    """
    assert cfg.probe_design_dir is not None
    assert cfg.genome_fasta_dir is not None

    db_path = cfg.probe_design_dir / "multifastas.sqlite"

    # NEW: avoid interactive [y/n] prompt from phyluce if file already exists
    if db_path.exists():
        logger.info(
            "Existing multi-fasta DB found; removing to avoid interactive "
            "overwrite prompt: %s",
            db_path,
        )
        db_path.unlink()

    cmd = phyluce_cmd(cfg, "phyluce_probe_get_multi_fasta_table") + [
        "--fastas",
        str(cfg.genome_fasta_dir.name),  # was: genome-fasta
        "--output",
        str(db_path.name),               # multifastas.sqlite
        "--base-taxon",
        cfg.base_taxon,
    ]

    logger.info("Running: %s", " ".join(cmd))
    run_command(cmd, cwd=cfg.probe_design_dir, logger=logger)
    return db_path



def query_multi_fasta_table(cfg: TailConfig, logger: logging.Logger, db_path: Path) -> Path:
    """
    phyluce_probe_query_multi_fasta_table:
      --db multifastas.sqlite
      --base-taxon <base>
      --output <base>+N-back-to-X.conf
      --specific-counts X
    """
    assert cfg.probe_design_dir is not None

    n_taxa = len(cfg.taxa)
    if cfg.multi_fasta_specific_counts is None:
        specific_counts = max(n_taxa - 2, 1)
    else:
        specific_counts = cfg.multi_fasta_specific_counts

    logger.info("Using specific-counts = %d (number of taxa sharing each locus).", specific_counts)

    conf_name = f"{cfg.base_taxon}+{n_taxa}-back-to-{specific_counts}.conf"
    conf_path = cfg.probe_design_dir / conf_name

    cmd = phyluce_cmd(cfg, "phyluce_probe_query_multi_fasta_table") + [
        "--db",
        str(db_path.name),
        "--base-taxon",
        cfg.base_taxon,
        "--output",
        conf_name,
        "--specific-counts",
        str(specific_counts),
    ]

    run_command(cmd, cwd=cfg.probe_design_dir, logger=logger)
    return conf_path


def design_master_probe_list(cfg: TailConfig, logger: logging.Logger, conf_path: Path) -> Path:
    """
    phyluce_probe_get_tiled_probe_from_multiple_inputs
    """
    assert cfg.probe_design_dir is not None
    assert cfg.genome_fasta_dir is not None

    out_name = f"{cfg.design_name}-master-probe-list.fasta"
    out_path = cfg.probe_design_dir / out_name

    cmd = phyluce_cmd(cfg, "phyluce_probe_get_tiled_probe_from_multiple_inputs") + [
        "--fastas",
        str(cfg.genome_fasta_dir.name),
        "--multi-fasta-output",
        str(conf_path.name),
        "--probe-prefix",
        cfg.probe_prefix,
        "--designer",
        cfg.designer,
        "--design",
        cfg.design_name,
        "--tiling-density",
        str(cfg.tiling_density),
        "--overlap",
        cfg.overlap_mode,
        "--masking",
        str(cfg.masking_threshold),
        "--remove-gc",
        "--two-probes",
        "--output",
        out_name,
    ]

    run_command(cmd, cwd=cfg.probe_design_dir, logger=logger)
    return out_path


def dupe_screen_master_probe_list(
    cfg: TailConfig, logger: logging.Logger, master_probe_path: Path
) -> Path:
    """
    Screen the master probe list against itself with LASTZ to identify and remove
    duplicate probes.

    Wraps:
      - phyluce_probe_easy_lastz
      - phyluce_probe_remove_duplicate_hits_from_probes_using_lastz

    and works around LASTZ's inability to handle spaces in absolute paths by
    running the LASTZ step in a temporary directory whose path has no spaces.
    """
    assert cfg.probe_design_dir is not None

    logger.info("==== STEP 2: Screening master probe list for duplicate probes ====")

    probe_dir = cfg.probe_design_dir

    # Filenames used by PHYLUCE
    master_name = master_probe_path.name
    lastz_filename = master_name.replace(".fasta", "-TO-SELF-PROBES.lastz")
    dupe_screened_name = master_name.replace(".fasta", "-DUPE-SCREENED.fasta")

    # Final locations (these may live on a volume with spaces in the path)
    final_lastz = probe_dir / lastz_filename
    dupe_screened = probe_dir / dupe_screened_name

    # Decide whether we need a space-free temp dir for LASTZ
    # LASTZ (as called inside phyluce_probe_easy_lastz) does *not* like spaces in
    # the absolute paths it receives.
    use_tmp = " " in str(probe_dir.resolve())
    tmp_dir: Optional[Path] = None

    if use_tmp:
        import tempfile

        tmp_dir = Path(tempfile.mkdtemp(prefix="phyluce_lastz_"))
        logger.info(
            "Probe-design directory contains spaces; running LASTZ in temp dir: %s",
            tmp_dir,
        )

        # Copy the master probe FASTA into the temp dir with the same file name
        tmp_master = tmp_dir / master_name
        shutil.copy2(master_probe_path, tmp_master)

        lastz_cwd = tmp_dir
        target_arg = tmp_master.name
        query_arg = tmp_master.name
        output_arg = lastz_filename
    else:
        lastz_cwd = probe_dir
        target_arg = master_name
        query_arg = master_name
        output_arg = lastz_filename

    # ------------------------------------------------------------------
    # 1) Run phyluce_probe_easy_lastz (self-LASTZ of master probe list)
    # ------------------------------------------------------------------
    cmd_lastz = phyluce_cmd(cfg, "phyluce_probe_easy_lastz") + [
        "--target",
        target_arg,
        "--query",
        query_arg,
        "--identity",
        str(cfg.lastz_identity),
        "--coverage",
        "50",  # PHYLUCE default for this step
        "--output",
        output_arg,
    ]

    logger.info(
        "Running phyluce_probe_easy_lastz from cwd=%s (space-safe cwd=%s)",
        probe_dir,
        lastz_cwd,
    )
    run_command(cmd_lastz, cwd=lastz_cwd, logger=logger)

    # Path to the LASTZ output where it was actually written
    lastz_out_path = lastz_cwd / lastz_filename

    # If we used a temp dir, copy LASTZ output back to the (possibly spacey) probe dir
    if use_tmp:
        shutil.copy2(lastz_out_path, final_lastz)
        logger.info(
            "Copied self-LASTZ output back to probe-design dir: %s",
            final_lastz,
        )
        lastz_out_path = final_lastz

        # Clean up temp directory (optional; comment out if you want to inspect it)
        try:
            shutil.rmtree(tmp_dir, ignore_errors=True)  # type: ignore[arg-type]
        except Exception as e:
            logger.warning("Could not remove temporary LASTZ directory %s: %s", tmp_dir, e)

    # ------------------------------------------------------------------
    # 2) Remove duplicate probes using the LASTZ results
    # ------------------------------------------------------------------
    # This step does not shell out to LASTZ, so spaces in the path are fine.
    cmd_dupe = phyluce_cmd(
        cfg,
        "phyluce_probe_remove_duplicate_hits_from_probes_using_lastz",
    ) + [
        "--fasta",
        master_name,
        "--lastz",
        lastz_out_path.name,
        "--probe-prefix",
        cfg.probe_prefix,
    ]

    logger.info(
        "Running phyluce_probe_remove_duplicate_hits_from_probes_using_lastz in %s",
        probe_dir,
    )
    run_command(cmd_dupe, cwd=probe_dir, logger=logger)

    logger.info("Dupe-screened master probe list written to %s", dupe_screened)
    return dupe_screened


def run_tail_pipeline(cfg: TailConfig) -> Path:
    """Run the full tail pipeline and return the final DUPE-SCREENED master probe list path."""
    log_level = getattr(logging, cfg.log_level.upper(), logging.INFO)
    logging.basicConfig(
        level=log_level,
        format="[%(levelname)s] %(message)s",
        stream=sys.stdout,
    )
    logger = logging.getLogger("phyluce_tail")

    logger.info("Starting PHYLUCE Tutorial-IV tail pipeline")
    logger.info("Project root: %s", cfg.project_root)
    logger.info("Genomes root: %s", cfg.genomes_root)
    logger.info("Base taxon:   %s", cfg.base_taxon)
    logger.info("Temp probes:  %s", cfg.temp_probes)

    if cfg.base_taxon not in cfg.taxa:
        logger.warning("Base taxon %s is not in --taxa; adding it.", cfg.base_taxon)
        cfg.taxa.append(cfg.base_taxon)

    setup_paths(cfg)

    # 1) Align temp probes (resumable)
    db_path = align_temp_probes_to_genomes(cfg, logger)

    # 2) Write scaffold config
    conf_path_scaffolds = write_scaffold_conf(cfg, logger)

    # 3) Slice sequences
    slice_sequence_from_genomes(cfg, logger, conf_path_scaffolds)

    # 4) Build multi-fasta table
    multifasta_db = build_multi_fasta_table(cfg, logger)

    # 5) Query multi-fasta table
    multi_conf = query_multi_fasta_table(cfg, logger, multifasta_db)

    # 6) Design master bait set
    master_probe = design_master_probe_list(cfg, logger, multi_conf)

    # 7) Final duplicate screen
    dupe_screened_master = dupe_screen_master_probe_list(cfg, logger, master_probe)

    logger.info("Finished tail pipeline.")
    logger.info("Final DUPE-SCREENED master probe list: %s", dupe_screened_master)
    return dupe_screened_master


def parse_args() -> TailConfig:
    p = argparse.ArgumentParser(
        description="Tail pipeline for PHYLUCE Tutorial IV: from temp probes to final master bait list."
    )

    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--genomes-root", required=False, type=Path)
    p.add_argument("--base-taxon", required=True)
    p.add_argument("--temp-probes", required=True, type=Path)
    p.add_argument("--taxa", nargs="+")
    p.add_argument("--lastz-identity", type=int, default=50)
    p.add_argument("--cores", type=int, default=4)
    p.add_argument("--slice-buffer", type=int, default=180)
    p.add_argument(
        "--multi-fasta-specific-counts",
        type=int,
        default=None,
    )
    p.add_argument("--probe-prefix", default="uce-")
    p.add_argument("--designer", default="faircloth")
    p.add_argument("--design-name", default="uce-design-v1")
    p.add_argument("--tiling-density", type=int, default=3)
    p.add_argument("--masking-threshold", type=float, default=0.25)
    p.add_argument("--overlap-mode", default="middle")
    p.add_argument("--log-level", default="INFO")
    p.add_argument(
        "--phyluce-prefix",
        default="",
        help="Optional prefix to run PHYLUCE tools (e.g. 'conda run -n phyluce-1.7.3').",
    )

    args = p.parse_args()

    project_root = args.project_root
    genomes_root = args.genomes_root or (project_root / "genomes")

    cfg = TailConfig(
        project_root=project_root,
        genomes_root=genomes_root,
        base_taxon=args.base_taxon,
        temp_probes=args.temp_probes,
        taxa=args.taxa or [],
        lastz_identity=args.lastz_identity,
        cores=args.cores,
        slice_buffer=args.slice_buffer,
        multi_fasta_specific_counts=args.multi_fasta_specific_counts,
        probe_prefix=args.probe_prefix,
        designer=args.designer,
        design_name=args.design_name,
        tiling_density=args.tiling_density,
        masking_threshold=args.masking_threshold,
        overlap_mode=args.overlap_mode,
        log_level=args.log_level,
        phyluce_prefix=args.phyluce_prefix,
    )
    return cfg


def main() -> None:
    cfg = parse_args()
    run_tail_pipeline(cfg)


if __name__ == "__main__":
    main()

