#!/usr/bin/env python3
"""
Safely sanitize genome FASTA filenames by removing '.', '_' and '-'
from the *stem* (the part before .fasta), while keeping '.fasta' intact.

Usage (dry run, no changes):
    python sanitize_genome_names.py "/media/localuser/T7 Shield/colpoptera"

Apply changes:
    python sanitize_genome_names.py "/media/localuser/T7 Shield/colpoptera" --apply

This will also write a mapping file:
    genome_name_mapping.tsv

with two columns: old_name    new_name
"""

import sys
from pathlib import Path

def build_new_name(path: Path) -> str:
    """Return the sanitized filename (string) for a given Path."""
    stem = path.stem          # e.g. "A15-CMM9.scaffolds"
    suffix = path.suffix      # e.g. ".fasta"

    # Remove ., _, - from the stem only
    cleaned_stem = stem.replace(".", "").replace("_", "").replace("-", "")

    return cleaned_stem + suffix


def main():
    if len(sys.argv) < 2:
        print("Usage: python sanitize_genome_names.py <directory> [--apply]")
        sys.exit(1)

    root = Path(sys.argv[1]).expanduser().resolve()
    apply_changes = "--apply" in sys.argv[2:]

    if not root.is_dir():
        print(f"Error: {root} is not a directory")
        sys.exit(1)

    print(f"Scanning directory: {root}")
    fasta_files = [p for p in root.iterdir() if p.is_file() and p.suffix.lower() == ".fasta"]

    if not fasta_files:
        print("No .fasta files found. Nothing to do.")
        sys.exit(0)

    mapping = []
    new_names = {}

    for p in sorted(fasta_files):
        new_name = build_new_name(p)
        if new_name == p.name:
            # No change needed
            continue

        mapping.append((p.name, new_name))

        if new_name in new_names and new_names[new_name] != p.name:
            # Collision: two different old names map to the same new name
            print("ERROR: Name collision detected!")
            print(f"  {new_names[new_name]}  ->  {new_name}")
            print(f"  {p.name}  ->  {new_name}")
            print("Aborting. Please resolve these manually.")
            sys.exit(1)
        else:
            new_names[new_name] = p.name

    if not mapping:
        print("All .fasta filenames are already sanitized. Nothing to rename.")
        sys.exit(0)

    # Write mapping file for reference
    mapping_file = root / "genome_name_mapping.tsv"
    with mapping_file.open("w", encoding="utf-8") as f:
        f.write("old_name\tnew_name\n")
        for old, new in mapping:
            f.write(f"{old}\t{new}\n")

    print("\nPlanned renames (old_name -> new_name):")
    for old, new in mapping:
        print(f"  {old}  ->  {new}")

    print(f"\nA mapping file has been written to: {mapping_file}")

    if not apply_changes:
        print("\nDRY RUN ONLY. No files were renamed.")
        print("If this looks correct, re-run with --apply to perform the renames.")
        sys.exit(0)

    # Apply renaming
    print("\nApplying renames...")
    for old, new in mapping:
        src = root / old
        dst = root / new
        print(f"  renaming {src.name} -> {dst.name}")
        src.rename(dst)

    print("\nDone. All listed files have been renamed.")
    print(f"Mapping is saved in: {mapping_file}")


if __name__ == "__main__":
    main()
