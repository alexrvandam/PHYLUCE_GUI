#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Wrapper to run phyluce_probe_run_multiple_lastzs_sqlite_mem.py
one species at a time, each in its own subprocess, all writing
into the same SQLite DB.

This avoids long-lived RAM growth in a single process.
"""

import argparse
import os
import sys
import subprocess


def get_args():
    p = argparse.ArgumentParser(
        description="Per-species driver for phyluce_probe_run_multiple_lastzs_sqlite_mem.py"
    )
    p.add_argument("--db", required=True, help="SQLite DB path (same as original script)")
    p.add_argument("--output", required=True, help="Output directory for LASTZ files")
    p.add_argument(
        "--probefile",
        required=True,
        help="Probe FASTA to align against genomes (same as original)",
    )
    p.add_argument(
        "--chromolist",
        type=str,
        nargs="+",
        default=[],
        help="List of chromosome-based genomes (optional; usually empty in your case)",
    )
    p.add_argument(
        "--scaffoldlist",
        type=str,
        nargs="+",
        default=[],
        help="List of scaffold-based genomes (this is your taxa list)",
    )
    p.add_argument(
        "--append",
        action="store_true",
        default=False,
        help="Append to an existing DB (ignored; we manage append internally)",
    )
    p.add_argument(
        "--no-dir",
        action="store_true",
        default=False,
        help="Pass-through to inner script",
    )
    p.add_argument(
        "--cores",
        type=int,
        default=1,
        help="Cores per LASTZ run (passed through to inner script)",
    )
    p.add_argument(
        "--genome-base-path",
        dest="base_path",
        type=str,
        default=None,
        help="Base path to directory containing genomes (as in original script)",
    )
    p.add_argument(
        "--coverage",
        type=float,
        default=83.0,
        help="LASTZ coverage threshold",
    )
    p.add_argument(
        "--identity",
        type=float,
        default=92.5,
        help="LASTZ identity threshold",
    )
    p.add_argument(
        "--verbosity",
        type=str,
        choices=["INFO", "WARN", "CRITICAL"],
        default="INFO",
        help="Logging level",
    )
    p.add_argument(
        "--log-path",
        dest="log_path",
        type=str,
        default=None,
        help="Log directory (if any)",
    )
    p.add_argument(
        "--inner-script",
        default="phyluce_probe_run_multiple_lastzs_sqlite_mem.py",
        help="Name or path of the per-species script (default: mem-safe version)",
    )
    return p.parse_args()


def main():
    args = get_args()

    # Resolve paths
    db_path = os.path.abspath(args.db)
    output_dir = os.path.abspath(args.output)
    probefile = os.path.abspath(args.probefile)
    base_path = args.base_path
    # NEW: ensure output directory exists once, non-interactively
    os.makedirs(output_dir, exist_ok=True)

    # Inner script: resolve relative to *this* file, unless already absolute
    if os.path.isabs(args.inner_script):
        inner_script = args.inner_script
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        inner_script = os.path.join(script_dir, args.inner_script)
    inner_script = os.path.abspath(inner_script)

    if not os.path.exists(inner_script):
        raise FileNotFoundError(f"Inner script not found: {inner_script}")

    scaffolds = args.scaffoldlist or []
    chromos = args.chromolist or []

    if chromos:
        # For now, just be explicit; can be extended later if needed.
        raise NotImplementedError(
            "chromolist handling is not implemented in per-species wrapper."
        )

    if not scaffolds:
        raise RuntimeError("No taxa in --scaffoldlist; nothing to run.")

    print(f"[WRAPPER] Will process {len(scaffolds)} scaffold taxa one by one.")
    print(f"[WRAPPER] DB:      {db_path}")
    print(f"[WRAPPER] Output:  {output_dir}")
    print(f"[WRAPPER] Probes:  {probefile}")
    print(f"[WRAPPER] Genomes: {base_path}")
    print(f"[WRAPPER] Inner script: {inner_script}")

    # Run each species in its own subprocess
    for i, taxon in enumerate(scaffolds):
        print(f"\n[WRAPPER] === Taxon {i+1}/{len(scaffolds)}: {taxon} ===")

        # Build command line for this taxon
        cmd = [
            sys.executable,
            inner_script,
            "--db",
            db_path,
            "--output",
            output_dir,
            "--probefile",
            probefile,
            "--scaffoldlist",
            taxon,
            "--genome-base-path",
            base_path,
            "--identity",
            str(args.identity),
            "--coverage",
            str(args.coverage),
            "--cores",
            str(args.cores),
            "--verbosity",
            args.verbosity,
        ]

        if args.no_dir:
            cmd.append("--no-dir")

        if args.log_path:
            cmd.extend(["--log-path", args.log_path])

        # For all but the first taxon, append to the existing DB
        if i > 0:
            cmd.append("--append")

        print("[WRAPPER CMD]", " ".join(cmd))

        # Run the inner script as a separate process
        result = subprocess.run(cmd)
        if result.returncode != 0:
            raise RuntimeError(
                f"Inner script failed for taxon {taxon} with exit code {result.returncode}"
            )

    print("\n[WRAPPER] All taxa completed successfully.")


if __name__ == "__main__":
    main()

