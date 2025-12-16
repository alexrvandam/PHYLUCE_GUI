#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Resumable version of phyluce_probe_run_multiple_lastzs_sqlite.

Key changes vs. the stock script:
    - No interactive prompts from phyluce.helpers.CreateFile/CreateDir.
    - New --resume flag to skip taxa that are already completed.
    - Simple per-taxon .done marker files in the LASTZ output directory.
"""

import os
import sys
import shutil
import sqlite3
import argparse

from phyluce.log import setup_logging
from phyluce.helpers import is_dir, is_file, FullPaths
from phyluce.many_lastz import multi_lastz_runner

from rich import print  # noqa: F401  (kept for consistency with phyluce style)


def get_args():
    """Get arguments from CLI (non-interactive)."""
    parser = argparse.ArgumentParser(
        description=(
            "Align a set of probes against genome sequence(s) in scaffolds or "
            "chromosomes (resumable version)."
        )
    )

    parser.add_argument(
        "--db",
        required=True,
        type=str,
        help="SQLite database file in which to store results.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=str,
        help="Directory in which to store the LASTZ files.",
    )
    parser.add_argument(
        "--probefile",
        required=True,
        type=is_file,
        action=FullPaths,
        help="The probe file to align against the sequences.",
    )
    parser.add_argument(
        "--chromolist",
        type=str,
        nargs="+",
        default=[],
        help="List of organisms with genome sequences in chromosomes.",
    )
    parser.add_argument(
        "--scaffoldlist",
        type=str,
        nargs="+",
        default=[],
        help="List of organisms with genome sequences in scaffolds/contigs.",
    )
    parser.add_argument(
        "--append",
        action="store_true",
        default=False,
        help="Insert results to an existing database (kept for compatibility).",
    )
    parser.add_argument(
        "--no-dir",
        action="store_true",
        default=False,
        help="If genome sequences are not in their own abbr. directory.",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help="Number of compute cores to use.",
    )
    parser.add_argument(
        "--genome-base-path",
        dest="base_path",
        type=str,
        default=None,
        help="Base path to a directory containing genomes sequences.",
    )
    parser.add_argument(
        "--coverage",
        type=float,
        default=83.0,
        help="Coverage threshold to search for using LASTZ.",
    )
    parser.add_argument(
        "--identity",
        type=float,
        default=92.5,
        help="Percent identity threshold for LASTZ.",
    )
    parser.add_argument(
        "--verbosity",
        type=str,
        choices=["INFO", "WARN", "CRITICAL"],
        default="INFO",
        help="Logging level to use.",
    )
    parser.add_argument(
        "--log-path",
        type=str,
        default=None,
        help="Path to a directory to hold logs.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        default=False,
        help=(
            "Resume a previous run: skip taxa that already have a .done marker "
            "and/or existing clean LASTZ output."
        ),
    )

    p = parser.parse_args()
    return p


def create_species_lastz_tables(cur, g):
    # Make table + index idempotent (IF NOT EXISTS) to support resume.
    query = f"""CREATE TABLE IF NOT EXISTS {g} (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            score INTEGER NOT NULL,
            name1 TEXT NOT NULL,
            strand1 TEXT NOT NULL,
            zstart1 INTEGER NOT NULL,
            end1 INTEGER NOT NULL,
            length1 INTEGER NOT NULL,
            name2 TEXT NOT NULL,
            strand2 TEXT NOT NULL,
            zstart2 INTEGER NOT NULL,
            end2 INTEGER NOT NULL,
            length2 INTEGER NOT NULL,
            diff TEXT NOT NULL,
            cigar TEXT NOT NULL,
            identity TEXT NOT NULL,
            percent_identity FLOAT NOT NULL,
            continuity TEXT NOT NULL,
            percent_continuity FLOAT NOT NULL,
            coverage TEXT NOT NULL,
            percent_coverage FLOAT NOT NULL)"""
    cur.execute(query)

    idx_query = f"""CREATE INDEX IF NOT EXISTS '{g}_name2_idx' on {g}(name2)"""
    cur.execute(idx_query)


def insert_species_to_lastz_tables(cur, g, input_path):
    with open(input_path, "r") as data:
        rows = [tuple(line.strip().split("\t")) for line in data]

    query = f"""INSERT INTO {g} (score, name1, strand1, zstart1, end1,
            length1, name2, strand2, zstart2, end2, length2, diff, cigar,
            identity, percent_identity, continuity, percent_continuity,
            coverage, percent_coverage) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,
            ?,?,?,?,?,?,?)"""

    cur.executemany(query, rows)


def clean_lastz_data(raw_path):
    """
    Convert a raw LASTZ file to a '.clean' version with '%' removed.

    - If '<raw>.clean' already exists, we simply return that.
    - If 'raw_path' exists, we create '<raw>.clean' and remove 'raw_path'.
    """
    raw_path = os.path.abspath(raw_path)
    clean_path = raw_path + ".clean"

    if os.path.isfile(clean_path) and not os.path.isfile(raw_path):
        # Already cleaned in a previous run
        return clean_path

    if not os.path.isfile(raw_path):
        raise IOError(
            f"Cannot clean LASTZ output: neither '{raw_path}' nor '{clean_path}' exist."
        )

    with open(raw_path, "r") as source, open(clean_path, "w") as target:
        data = source.read()
        target.write(data.replace("%", ""))

    os.remove(raw_path)
    return clean_path


def create_species_table(conn, cur, append):
    """
    Create the 'species' table if it does not exist.

    Any previous prompt-based behavior has been removed to avoid interactivity.
    """
    cur.execute(
        """CREATE TABLE IF NOT EXISTS species (
                name TEXT PRIMARY KEY,
                description TEXT NULL,
                version TEXT NULL
                )
        """
    )
    conn.commit()


def species_done_marker(output_dir, probefile_basename, taxon):
    """
    Return the path of the '.done' marker for a given taxon.
    """
    raw_name = f"{probefile_basename}_v_{taxon}.lastz"
    return os.path.abspath(os.path.join(output_dir, raw_name + ".done"))


def align_against_scaffolds(log, cur, args, path, resume=False):
    probefile_basename = os.path.basename(args.probefile)

    for g in args.scaffoldlist:
        raw_output = os.path.abspath(
            os.path.join(args.output, f"{probefile_basename}_v_{g}.lastz")
        )
        clean_output = raw_output + ".clean"
        done_marker = raw_output + ".done"
        target = path.format(g)

        if resume and os.path.isfile(done_marker):
            log.info("[RESUME] %s already completed (.done present); skipping.", g)
            continue

        if resume and os.path.isfile(clean_output):
            log.info("[RESUME] Found existing clean LASTZ for %s; skipping alignment.", g)
        else:
            if resume and os.path.isfile(raw_output):
                log.info(
                    "[RESUME] Using existing raw LASTZ output for %s; cleaning only.", g
                )
            else:
                log.info("Aligning against %s scaffolds", g)
                multi_lastz_runner(
                    log,
                    raw_output,
                    args.cores,
                    target,
                    args.probefile,
                    True,  # scaffolds=True
                    args.coverage,
                    args.identity,
                )

            log.info("Cleaning the LASTZ output for %s", g)
            clean_output = clean_lastz_data(raw_output)

        if args.db:
            log.info("Creating the %s table (if needed)", g)
            create_species_lastz_tables(cur, g)

            log.info("Inserting data to the %s table", g)
            insert_species_to_lastz_tables(cur, g, clean_output)

            # Avoid duplicate entries in 'species'
            cur.execute("INSERT OR IGNORE INTO species (name) VALUES (?)", (g,))

        # Mark this taxon as done
        with open(done_marker, "w") as fh:
            fh.write("")

        log.info("Finished processing %s scaffolds.", g)


def align_against_genomes(log, cur, args, path, resume=False):
    probefile_basename = os.path.basename(args.probefile)

    for g in args.chromolist:
        raw_output = os.path.abspath(
            os.path.join(args.output, f"{probefile_basename}_v_{g}.lastz")
        )
        clean_output = raw_output + ".clean"
        done_marker = raw_output + ".done"
        target = path.format(g)

        if resume and os.path.isfile(done_marker):
            log.info("[RESUME] %s already completed (.done present); skipping.", g)
            continue

        if resume and os.path.isfile(clean_output):
            log.info("[RESUME] Found existing clean LASTZ for %s; skipping alignment.", g)
        else:
            if resume and os.path.isfile(raw_output):
                log.info(
                    "[RESUME] Using existing raw LASTZ output for %s; cleaning only.", g
                )
            else:
                log.info("Aligning against %s chromosomes", g)
                multi_lastz_runner(
                    log,
                    raw_output,
                    args.cores,
                    target,
                    args.probefile,
                    False,  # scaffolds=False => chromosomes
                    args.coverage,
                    args.identity,
                )

            log.info("Cleaning the LASTZ output for %s", g)
            clean_output = clean_lastz_data(raw_output)

        if args.db:
            log.info("Creating the %s table (if needed)", g)
            create_species_lastz_tables(cur, g)

            log.info("Inserting data to the %s table", g)
            insert_species_to_lastz_tables(cur, g, clean_output)

            cur.execute("INSERT OR IGNORE INTO species (name) VALUES (?)", (g,))

        with open(done_marker, "w") as fh:
            fh.write("")

        log.info("Finished processing %s chromosomes.", g)


def check_for_all_genome_sequences(chromo, scaffold, path):
    genomes = chromo + scaffold
    for g in genomes:
        file_path = path.format(g)
        if not os.path.isfile(file_path):
            raise IOError(f"{file_path} is not a file")


def main():
    args = get_args()

    # Normalize paths
    args.db = os.path.abspath(args.db)
    args.output = os.path.abspath(args.output)
    if args.base_path is None:
        raise ValueError("--genome-base-path must be specified")

    # Setup logging
    log, my_name = setup_logging(args)

    # Ensure output dir exists but DO NOT delete it
    os.makedirs(args.output, exist_ok=True)

    # Connect to DB
    conn = sqlite3.connect(args.db)
    cur = conn.cursor()

    # Determine genome .2bit path template
    if not args.no_dir:
        path = os.path.join("{0}".format(args.base_path), "{0}/{0}.2bit")
    else:
        path = os.path.join("{0}".format(args.base_path), "{0}.2bit")

    log.info("Checking to ensure all genome sequences exist as 2bit files")
    check_for_all_genome_sequences(args.chromolist, args.scaffoldlist, path)

    log.info("Created species table")
    create_species_table(conn, cur, args.append)

    # Align scaffolds + genomes with optional resume
    align_against_scaffolds(log, cur, args, path, resume=args.resume)
    align_against_genomes(log, cur, args, path, resume=args.resume)

    conn.commit()
    cur.close()
    conn.close()

    text = f" Completed {my_name} "
    log.info(text.center(65, "="))


if __name__ == "__main__":
    main()

