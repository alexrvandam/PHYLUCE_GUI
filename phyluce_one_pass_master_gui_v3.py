#!/usr/bin/env python3
"""
phyluce_one_pass_master_gui.py

Master GUI that orchestrates the three PHYLUCE Tutorial IV wrapper scripts:

  Step 1: phyluce_one_passV4.py
    - From cleaned draft genomes → project structure → temp DUPE-SCREENED probes
    - Here we extend Step 1 to show the multi-merge summary and let the user
      interactively choose the "specific-counts" parameter before BED/FASTA/probes.

  Step 2: phyluce_one_pass_tail_v3.py
    - From temp DUPE-SCREENED probes → multi-genome master DUPE-SCREENED bait set

  Step 3: phyluce_one_pass_to_raxml.py
    - In-silico test of bait set → alignments → concatenated matrix → RAxML tree

Place this script in the same directory as:
  - phyluce_one_passV4.py
  - phyluce_one_pass_tail_v3.py
  - phyluce_one_pass_to_raxml.py

Run with:

  python phyluce_one_pass_master_gui.py

"""

import logging
import os
import sys
import subprocess
from pathlib import Path
import shutil
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, simpledialog
# Import your three pipeline modules
import phyluce_one_passV4 as step1
import phyluce_one_pass_tail_v3 as step2
import phyluce_one_pass_to_raxml as step3
import phyluce_post_uce_occupancy as step4


# ---------------------------------------------------------------------
# Logging & stdout redirection helpers
# ---------------------------------------------------------------------


class TextHandler(logging.Handler):
    """Logging handler that writes to a Tkinter Text widget."""
    def __init__(self, text_widget: tk.Text):
        super().__init__()
        self.text_widget = text_widget

    def emit(self, record: logging.LogRecord) -> None:
        msg = self.format(record)
        self.text_widget.after(0, self._append, msg)

    def _append(self, msg: str) -> None:
        self.text_widget.insert(tk.END, msg + "\n")
        self.text_widget.see(tk.END)


class TextRedirector:
    """File-like object that redirects writes to a Tkinter Text widget."""
    def __init__(self, text_widget: tk.Text, prefix: str = ""):
        self.text_widget = text_widget
        self.prefix = prefix

    def write(self, s: str) -> None:
        if not s:
            return
        if s.strip() == "":
            return

        def append():
            self.text_widget.insert(tk.END, self.prefix + s)
            self.text_widget.see(tk.END)

        self.text_widget.after(0, append)

    def flush(self) -> None:
        pass


# ---------------------------------------------------------------------
# Master GUI
# ---------------------------------------------------------------------


def launch_master_gui() -> None:
    root = tk.Tk()
    root.title("PHYLUC-E One-pass Master GUI (Tutorial IV, 3-step pipeline)")

    # Overall layout: notebook on top, log at bottom
    root.rowconfigure(0, weight=1)
    root.rowconfigure(1, weight=0)
    root.columnconfigure(0, weight=1)

    notebook = ttk.Notebook(root)
    notebook.grid(row=0, column=0, sticky="nsew")

    # Shared variables across tabs
    genomes_root_var = tk.StringVar()
    project_root_var = tk.StringVar()
    base_taxon_var = tk.StringVar()

    # Step 1 specific vars
    step1_threads_var = tk.StringVar(value="8")
    step1_mm_var = tk.StringVar(value="5")  # initial guess; user can override after preview
    step1_merge_d_var = tk.StringVar(value="0")
    step1_min_cov_var = tk.StringVar(value="67")
    step1_sub_rate_var = tk.StringVar(value="0.05")
    step1_ins_size_var = tk.StringVar(value="400")
    step1_extra_opts_var = tk.StringVar(value="")
    step1_dry_run_var = tk.BooleanVar(value=False)

    # Step 2 specific vars
    temp_probes_var = tk.StringVar(value="")
    step2_specific_counts_var = tk.StringVar(value="")  # optional override; default n_taxa-2
    step2_design_name_var = tk.StringVar(value="uce-design-v1")
    step2_tiling_density_var = tk.StringVar(value="3")
    step2_masking_var = tk.StringVar(value="0.25")
    step2_lastz_cores_var = tk.StringVar(value="4")
    step2_slice_buffer_var = tk.StringVar(value="180")

    # Step 3 specific vars
    step3_probe_fasta_var = tk.StringVar(value="")
    step3_genome_conf_var = tk.StringVar(value="")
    step3_min_taxon_prop_var = tk.StringVar(value="0.75")
    step3_align_cores_var = tk.StringVar(value="12")
    step3_lastz_cores_var = tk.StringVar(value="16")
    step3_flank_var = tk.StringVar(value="400")      # NEW: flank size (bp)
    step3_raxml_threads_var = tk.StringVar(value="10")

    # Step 4 specific vars
    step4_min_loci_var = tk.StringVar(value="100")
    step4_insilico_dir_var = tk.StringVar(value="")
    step4_genomes_conf_var = tk.StringVar(value="")
    step4_output_dir_var = tk.StringVar(value="")

    # Cache of taxa names (shared)
    taxa_cache = {"step1": [], "step2": []}

    # -----------------------------------------------------------------
    # Log frame (bottom)
    # -----------------------------------------------------------------
    log_frame = ttk.LabelFrame(root, text="Log")
    log_frame.grid(row=1, column=0, sticky="nsew", padx=5, pady=5)
    log_frame.rowconfigure(0, weight=1)
    log_frame.columnconfigure(0, weight=1)

    log_text = tk.Text(log_frame, wrap="word", height=15)
    log_text.grid(row=0, column=0, sticky="nsew")
    log_scroll = ttk.Scrollbar(log_frame, orient="vertical", command=log_text.yview)
    log_scroll.grid(row=0, column=1, sticky="ns")
    log_text.configure(yscrollcommand=log_scroll.set)

    # Set up logger
    logger = logging.getLogger("phyluce_master_gui")
    logger.setLevel(logging.INFO)
    logger.handlers = []
    th = TextHandler(log_text)
    th.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.addHandler(th)

    # -----------------------------------------------------------------
    # Shared helper: browse directory and browse file
    # -----------------------------------------------------------------
    def browse_dir(target_var: tk.StringVar, title: str) -> None:
        d = filedialog.askdirectory(title=title)
        if d:
            target_var.set(d)
            
    def browse_file(target_var: tk.StringVar, title: str, filetypes=None) -> None:
        f = filedialog.askopenfilename(
            title=title,
            filetypes=filetypes or [("All files", "*.*")]
        )
        if f:
            target_var.set(f)

    # -----------------------------------------------------------------
    # Helper: preview multi-merge counts (Step 1)
    # -----------------------------------------------------------------
    def preview_multi_merge_counts(cfg: step1.PipelineConfig, db_path: Path) -> str:
        """Run phyluce_probe_query_multi_merge_table in summary mode and return stdout."""
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

    # -----------------------------------------------------------------
    # STEP 1 TAB: genomes → temp DUPE-SCREENED probes
    # -----------------------------------------------------------------
    tab1 = ttk.Frame(notebook)
    notebook.add(tab1, text="Step 1: temp baits")

    tab1.rowconfigure(1, weight=1)
    tab1.columnconfigure(0, weight=1)

    frm1_paths = ttk.LabelFrame(tab1, text="Paths")
    frm1_paths.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
    frm1_paths.columnconfigure(1, weight=1)

    ttk.Label(frm1_paths, text="Genomes root:").grid(row=0, column=0, sticky="w")
    entry_genomes = ttk.Entry(frm1_paths, textvariable=genomes_root_var, width=40)
    entry_genomes.grid(row=0, column=1, sticky="ew", padx=2)
    ttk.Button(
        frm1_paths,
        text="Browse",
        command=lambda: browse_dir(genomes_root_var, "Select genomes root"),
    ).grid(row=0, column=2, padx=2)

    ttk.Label(frm1_paths, text="Project root:").grid(row=1, column=0, sticky="w")
    entry_project = ttk.Entry(frm1_paths, textvariable=project_root_var, width=40)
    entry_project.grid(row=1, column=1, sticky="ew", padx=2)
    ttk.Button(
        frm1_paths,
        text="Browse",
        command=lambda: browse_dir(project_root_var, "Select or create project root"),
    ).grid(row=1, column=2, padx=2)

    frm1_taxa = ttk.LabelFrame(tab1, text="Taxa (Step 1)")
    frm1_taxa.grid(row=1, column=0, sticky="nsew", padx=5, pady=5)
    frm1_taxa.columnconfigure(1, weight=1)
    frm1_taxa.rowconfigure(2, weight=1)

    btn_scan1 = ttk.Button(frm1_taxa, text="Scan genomes")

    ttk.Label(frm1_taxa, text="Base taxon:").grid(row=1, column=0, sticky="w")
    combo_base1 = ttk.Combobox(
        frm1_taxa, textvariable=base_taxon_var, values=[], state="readonly", width=20
    )
    combo_base1.grid(row=1, column=1, sticky="ew", padx=2)

    ttk.Label(frm1_taxa, text="Exemplar taxa:").grid(row=2, column=0, sticky="nw")
    list_exemplars1 = tk.Listbox(
        frm1_taxa, selectmode=tk.MULTIPLE, height=8, exportselection=False
    )
    list_exemplars1.grid(row=2, column=1, columnspan=2, sticky="nsew", padx=2, pady=2)

    def do_scan_genomes_step1() -> None:
        root_dir = genomes_root_var.get().strip()
        if not root_dir:
            messagebox.showerror("Error", "Please set the genomes root directory first.")
            return
        try:
            genomes = step1.discover_genomes(Path(root_dir), logger)
        except Exception as e:
            messagebox.showerror("Error scanning genomes", str(e))
            return

        taxa = sorted(genomes.keys())
        taxa_cache["step1"] = taxa
        combo_base1["values"] = taxa
        if taxa:
            base_taxon_var.set(taxa[0])

        list_exemplars1.delete(0, tk.END)
        for t in taxa:
            list_exemplars1.insert(tk.END, t)
        for i in range(len(taxa)):
            list_exemplars1.select_set(i)

        logger.info(f"Scanned genomes (Step 1): {', '.join(taxa)}")

    btn_scan1.config(command=do_scan_genomes_step1)
    btn_scan1.grid(row=0, column=0, columnspan=3, pady=2)

    frm1_params = ttk.LabelFrame(tab1, text="Key parameters (Step 1)")
    frm1_params.grid(row=2, column=0, sticky="ew", padx=5, pady=5)
    frm1_params.columnconfigure(1, weight=1)

    ttk.Label(frm1_params, text="Threads (stampy):").grid(row=0, column=0, sticky="w")
    entry_threads = ttk.Entry(frm1_params, width=6, textvariable=step1_threads_var)
    entry_threads.grid(row=0, column=1, sticky="w", padx=2)

    ttk.Label(frm1_params, text="Min shared taxa (specific-counts):").grid(row=1, column=0, sticky="w")
    entry_mm = ttk.Entry(frm1_params, width=6, textvariable=step1_mm_var)
    entry_mm.grid(row=1, column=1, sticky="w", padx=2)

    ttk.Label(frm1_params, text="BED merge distance (-d):").grid(row=2, column=0, sticky="w")
    entry_merge_d = ttk.Entry(frm1_params, width=6, textvariable=step1_merge_d_var)
    entry_merge_d.grid(row=2, column=1, sticky="w", padx=2)

    ttk.Label(frm1_params, text="Min coverage (contigs<->probes, future):").grid(row=3, column=0, sticky="w")
    entry_min_cov = ttk.Entry(frm1_params, width=6, textvariable=step1_min_cov_var)
    entry_min_cov.grid(row=3, column=1, sticky="w", padx=2)

    ttk.Label(frm1_params, text="Stampy substitution rate:").grid(row=4, column=0, sticky="w")
    entry_sub_rate = ttk.Entry(frm1_params, width=6, textvariable=step1_sub_rate_var)
    entry_sub_rate.grid(row=4, column=1, sticky="w", padx=2)

    ttk.Label(frm1_params, text="Stampy insert size:").grid(row=5, column=0, sticky="w")
    entry_ins_size = ttk.Entry(frm1_params, width=6, textvariable=step1_ins_size_var)
    entry_ins_size.grid(row=5, column=1, sticky="w", padx=2)

    ttk.Label(frm1_params, text="Stampy extra options:").grid(row=6, column=0, sticky="w")
    entry_extra_opts = ttk.Entry(frm1_params, width=20, textvariable=step1_extra_opts_var)
    entry_extra_opts.grid(row=6, column=1, sticky="w", padx=2)

    ttk.Checkbutton(
        frm1_params,
        text="Dry-run (log commands only)",
        variable=step1_dry_run_var,
    ).grid(row=7, column=0, columnspan=2, sticky="w", pady=2)

    frm1_buttons = ttk.Frame(tab1)
    frm1_buttons.grid(row=3, column=0, sticky="e", padx=5, pady=5)

    def prompt_specific_counts(counts_text: str, default_value: int) -> int:
        """Small modal dialog to let the user choose specific-counts."""
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

        top.transient(root)
        top.grab_set()
        root.wait_window(top)
        return chosen["value"]

    def run_step1_with_prompt() -> None:
        genomes_root = genomes_root_var.get().strip()
        project_root = project_root_var.get().strip()
        base_taxon = base_taxon_var.get().strip()

        if not genomes_root or not project_root or not base_taxon:
            messagebox.showerror(
                "Missing input",
                "Please set genomes root, project root, and select a base taxon."
            )
            return

        try:
            threads = int(step1_threads_var.get().strip())
            mm_count = int(step1_mm_var.get().strip())
            merge_d = int(step1_merge_d_var.get().strip())
            min_cov = int(step1_min_cov_var.get().strip())
            sub_rate = float(step1_sub_rate_var.get().strip())
            ins_size = int(step1_ins_size_var.get().strip())
        except ValueError:
            messagebox.showerror("Invalid input", "Numeric parameters must be valid numbers.")
            return

        extra_opts = step1_extra_opts_var.get().strip()

        selected_indices = list_exemplars1.curselection()
        if not selected_indices:
            # Use all taxa by default
            selected_indices = range(list_exemplars1.size())
        exemplar_taxa = [list_exemplars1.get(i) for i in selected_indices]
        if base_taxon not in exemplar_taxa:
            exemplar_taxa.append(base_taxon)

        cfg = step1.PipelineConfig(
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
            dry_run=step1_dry_run_var.get(),
        )

        try:
            logger.info("==== STEP 1: Starting genomes → temp DUPE-SCREENED probes (interactive specific-counts) ====")
            # This mirrors step1.run_full_pipeline, but inserts a preview/prompt
            genomes = step1.discover_genomes(cfg.genomes_root, logger)
            if cfg.base_taxon not in genomes:
                raise ValueError(
                    f"Base taxon '{cfg.base_taxon}' not found among discovered genomes: "
                    f"{', '.join(sorted(genomes))}"
                )

            missing = [t for t in cfg.exemplar_taxa if t not in genomes]
            if missing:
                raise ValueError(f"Exemplar taxa not found among genomes: {', '.join(missing)}")

            subset_genomes = {t: genomes[t] for t in cfg.exemplar_taxa}
            subset_genomes[cfg.base_taxon] = genomes[cfg.base_taxon]

            project_genomes = step1.setup_project_structure(cfg, subset_genomes, logger)
            step1.convert_genomes_to_2bit(cfg, project_genomes, logger)
            step1.simulate_reads(cfg, project_genomes, logger)
            step1.prepare_base_genome(cfg, project_genomes, logger)
            step1.align_reads_with_stampy(cfg, logger)
            step1.filter_bams_and_link(cfg, logger)
            step1.bams_to_beds(cfg, logger)
            step1.sort_beds(cfg, logger)
            step1.merge_beds(cfg, logger)
            step1.strip_masked_loci(cfg, logger)
            step1.write_bed_conf(cfg, logger)
            db_path = step1.multi_merge_table(cfg, logger)

            # Preview the distribution of shared loci
            summary = preview_multi_merge_counts(cfg, db_path)

            # Use current GUI field as default
            try:
                default_mm = int(step1_mm_var.get().strip())
            except ValueError:
                default_mm = cfg.multi_merge_specific_counts

            chosen_mm = prompt_specific_counts(summary, default_mm)
            cfg.multi_merge_specific_counts = chosen_mm
            step1_mm_var.set(str(chosen_mm))
            logger.info(f"Using specific-counts={chosen_mm} for downstream BED and temp bait design.")

            bed_path = step1.query_multi_merge_table(cfg, db_path, logger)
            fasta_path = step1.extract_genome_sequences_from_bed(cfg, bed_path, logger)
            temp_probes = step1.design_temp_baits(cfg, fasta_path, logger)
            dupe_temp = step1.dupe_screen_temp_baits(cfg, temp_probes, logger)

            temp_probes_var.set(str(dupe_temp))
            logger.info(f"STEP 1 produced temp DUPE-SCREENED probes: {dupe_temp}")
            messagebox.showinfo("Step 1 complete", "Step 1 completed successfully. Check the log for details.")
        except Exception as e:
            logger.exception("Step 1 failed")
            messagebox.showerror("Step 1 failed", str(e))

    ttk.Button(frm1_buttons, text="Run Step 1", command=run_step1_with_prompt).grid(row=0, column=0, padx=5)
    ttk.Button(frm1_buttons, text="Quit", command=root.destroy).grid(row=0, column=1, padx=5)

    # -----------------------------------------------------------------
    # STEP 2 TAB: temp probes → master DUPE-SCREENED bait set
    # (unchanged logic, still wired into phyluce_one_pass_tail.py)
    # -----------------------------------------------------------------
    tab2 = ttk.Frame(notebook)
    notebook.add(tab2, text="Step 2: screened bait set")

    tab2.rowconfigure(2, weight=1)
    tab2.columnconfigure(0, weight=1)

    frm2_paths = ttk.LabelFrame(tab2, text="Paths (Step 2)")
    frm2_paths.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
    frm2_paths.columnconfigure(1, weight=1)

    ttk.Label(frm2_paths, text="Project root:").grid(row=0, column=0, sticky="w")
    ttk.Entry(frm2_paths, textvariable=project_root_var, width=40).grid(
        row=0, column=1, sticky="ew", padx=2
    )
    ttk.Button(
        frm2_paths,
        text="Browse",
        command=lambda: browse_dir(project_root_var, "Select project root"),
    ).grid(row=0, column=2, padx=2)

    ttk.Label(frm2_paths, text="Base taxon:").grid(row=1, column=0, sticky="w")
    ttk.Entry(frm2_paths, textvariable=base_taxon_var, width=20).grid(
        row=1, column=1, sticky="w", padx=2
    )

    ttk.Label(frm2_paths, text="Temp DUPE-SCREENED probes:").grid(row=2, column=0, sticky="w")
    entry_temp_probes = ttk.Entry(frm2_paths, textvariable=temp_probes_var, width=40)
    entry_temp_probes.grid(row=2, column=1, sticky="ew", padx=2)

    def guess_temp_probes_from_step1() -> None:
        proj = project_root_var.get().strip()
        base = base_taxon_var.get().strip()
        mm = step1_mm_var.get().strip()
        if proj and base and mm:
            guess = Path(proj) / "probe-design" / f"{base}+{mm}.temp-DUPE-SCREENED.probes"
            temp_probes_var.set(str(guess))
            logger.info(f"Guessed temp probes path for Step 2: {guess}")

    ttk.Button(
        frm2_paths,
        text="Guess from Step 1",
        command=guess_temp_probes_from_step1,
    ).grid(row=2, column=2, padx=2)

    def browse_temp_probes() -> None:
        f = filedialog.askopenfilename(
            title="Select temp DUPE-SCREENED probes file",
            filetypes=[("Probes", "*.probes"), ("All files", "*.*")]
        )
        if f:
            temp_probes_var.set(f)

    ttk.Button(
        frm2_paths,
        text="Browse file",
        command=browse_temp_probes,
    ).grid(row=3, column=2, padx=2, sticky="e")

    frm2_taxa = ttk.LabelFrame(tab2, text="Taxa for Step 2 (genomes/ subdirectories)")
    frm2_taxa.grid(row=1, column=0, sticky="nsew", padx=5, pady=5)
    frm2_taxa.columnconfigure(0, weight=1)
    frm2_taxa.rowconfigure(1, weight=1)

    btn_scan2 = ttk.Button(frm2_taxa, text="Scan project/genomes")
    list_taxa2 = tk.Listbox(frm2_taxa, selectmode=tk.MULTIPLE, height=8, exportselection=False)

    def do_scan_genomes_step2() -> None:
        proj = project_root_var.get().strip()
        if not proj:
            messagebox.showerror("Error", "Please set the project root first.")
            return
        genomes_dir = Path(proj) / "genomes"
        if not genomes_dir.is_dir():
            messagebox.showerror("Error", f"Genomes directory not found: {genomes_dir}")
            return
        taxa = sorted([p.name for p in genomes_dir.iterdir() if p.is_dir()])
        if not taxa:
            messagebox.showerror("Error", f"No taxon subdirectories found in {genomes_dir}")
            return
        taxa_cache["step2"] = taxa
        list_taxa2.delete(0, tk.END)
        for t in taxa:
            list_taxa2.insert(tk.END, t)
        for i in range(len(taxa)):
            list_taxa2.select_set(i)
        logger.info(f"Scanned genomes (Step 2): {', '.join(taxa)}")

    btn_scan2.config(command=do_scan_genomes_step2)
    btn_scan2.grid(row=0, column=0, pady=2, sticky="w")
    list_taxa2.grid(row=1, column=0, sticky="nsew", padx=2, pady=2)

    frm2_params = ttk.LabelFrame(tab2, text="Key parameters (Step 2)")
    frm2_params.grid(row=2, column=0, sticky="ew", padx=5, pady=5)
    frm2_params.columnconfigure(1, weight=1)

    ttk.Label(frm2_params, text="specific-counts (multi-fasta):").grid(row=0, column=0, sticky="w")
    ttk.Entry(frm2_params, width=6, textvariable=step2_specific_counts_var).grid(
        row=0, column=1, sticky="w", padx=2
    )
    ttk.Label(frm2_params, text="(Leave blank for default n_taxa-2)").grid(row=0, column=2, sticky="w")

    ttk.Label(frm2_params, text="Design name (master baits):").grid(row=1, column=0, sticky="w")
    ttk.Entry(frm2_params, width=20, textvariable=step2_design_name_var).grid(
        row=1, column=1, sticky="w", padx=2
    )

    ttk.Label(frm2_params, text="Tiling density:").grid(row=2, column=0, sticky="w")
    ttk.Entry(frm2_params, width=6, textvariable=step2_tiling_density_var).grid(
        row=2, column=1, sticky="w", padx=2
    )

    ttk.Label(frm2_params, text="Masking threshold:").grid(row=3, column=0, sticky="w")
    ttk.Entry(frm2_params, width=6, textvariable=step2_masking_var).grid(
        row=3, column=1, sticky="w", padx=2
    )

    ttk.Label(frm2_params, text="LASTZ cores:").grid(row=4, column=0, sticky="w")
    ttk.Entry(frm2_params, width=6, textvariable=step2_lastz_cores_var).grid(
        row=4, column=1, sticky="w", padx=2
    )

    ttk.Label(frm2_params, text="Slice buffer (bp):").grid(row=5, column=0, sticky="w")
    ttk.Entry(frm2_params, width=6, textvariable=step2_slice_buffer_var).grid(
        row=5, column=1, sticky="w", padx=2
    )

    frm2_buttons = ttk.Frame(tab2)
    frm2_buttons.grid(row=3, column=0, sticky="e", padx=5, pady=5)

    def run_step2() -> None:
        proj = project_root_var.get().strip()
        base = base_taxon_var.get().strip()
        temp_probes = temp_probes_var.get().strip()

        if not proj or not base or not temp_probes:
            messagebox.showerror(
                "Missing input",
                "Please set project root, base taxon, and temp DUPE-SCREENED probes file."
            )
            return

        genomes_root = Path(proj) / "genomes"
        taxa_sel_idx = list_taxa2.curselection()
        if not taxa_sel_idx:
            taxa_sel_idx = range(list_taxa2.size())
        taxa = [list_taxa2.get(i) for i in taxa_sel_idx]
        if base not in taxa:
            taxa.append(base)

        try:
            cores = int(step2_lastz_cores_var.get().strip())
            slice_buffer = int(step2_slice_buffer_var.get().strip())
            masking = float(step2_masking_var.get().strip())
            tiling = int(step2_tiling_density_var.get().strip())
        except ValueError:
            messagebox.showerror("Invalid input", "Step 2 numeric parameters must be valid numbers.")
            return

        # -----------------------------------------------------------
        # Step 2: optional specific-counts
        # -----------------------------------------------------------
        if step2_specific_counts_var.get().strip():
            try:
                specific_counts = int(step2_specific_counts_var.get().strip())
            except ValueError:
                messagebox.showerror("Invalid input", "specific-counts must be an integer.")
                return
        else:
            # Let the tail code use its own default behaviour
            specific_counts = None

        # -----------------------------------------------------------
        # Step 2: choose LASTZ output directory (no terminal prompt)
        # -----------------------------------------------------------
        # PHYLUCE's `phyluce_probe_run_multiple_lastzs_sqlite` will stop and
        # ask interactively:
        #   "[WARNING] Output directory exists, REMOVE [y/n]:"
        # if the output directory (usually "genome-lastz") already exists.
        # In a GUI (no interactive terminal), that raises EOFError.
        #
        # Here we:
        #   * use "genome-lastz" by default
        #   * if it already exists, ask the user for a NEW folder name
        #     and ensure it does not already exist
        #   * pass that name down to the tail pipeline so the CLI never
        #     needs to ask.
        proj_path = Path(proj)
        probe_design_dir = proj_path / "probe-design"

        # Default LASTZ output dir name used by the PHYLUCE CLI
        default_lastz_name = "genome-lastz"
        default_lastz_path = probe_design_dir / default_lastz_name
        lastz_output_dir_name = default_lastz_name

        # If the default LASTZ output directory already exists, ask the user
        # for a NEW folder name to use for THIS run.
        if default_lastz_path.exists():
            while True:
                msg = (
                    "The default LASTZ output directory already exists:\n\n"
                    f"{default_lastz_path}\n\n"
                    "PHYLUCE normally stops and asks in the terminal:\n"
                    "  [WARNING] Output directory exists, REMOVE [y/n]\n"
                    "Because this GUI cannot answer that prompt, Step 2 would fail.\n\n"
                    "Please enter a NEW folder name to use for this run.\n"
                    "The folder will be created under:\n"
                    f"  {probe_design_dir}\n\n"
                    "Example: genome-lastz-run2"
                )
                new_name = simpledialog.askstring(
                    "LASTZ output directory already exists",
                    msg,
                    initialvalue="genome-lastz-run2",
                    parent=root,
                )

                # User hit Cancel or closed the dialog
                if new_name is None:
                    messagebox.showinfo(
                        "Step 2 cancelled",
                        "Step 2 was cancelled; the existing LASTZ directory was left unchanged.",
                    )
                    return

                new_name = new_name.strip()
                if not new_name:
                    messagebox.showerror(
                        "Invalid folder name",
                        "Folder name cannot be empty. Please enter a valid name.",
                    )
                    continue

                candidate = probe_design_dir / new_name
                if candidate.exists():
                    messagebox.showerror(
                        "Folder already exists",
                        f"The folder:\n\n{candidate}\n\nalready exists.\n"
                        "Please choose a different name.",
                    )
                    continue

                # We just remember the name; the PHYLUCE CLI will create
                # the folder when it runs with --output <name>.
                lastz_output_dir_name = new_name
                break

        # Full path to the LASTZ output directory for this run
        lastz_output_dir = probe_design_dir / lastz_output_dir_name

        # -----------------------------------------------------------
        # Build TailConfig (note the genome_lastz_dir argument)
        # -----------------------------------------------------------
        cfg2 = step2.TailConfig(
            project_root=proj_path,
            genomes_root=genomes_root,
            base_taxon=base,
            temp_probes=Path(temp_probes),
            taxa=taxa,
            cores=cores,
            slice_buffer=slice_buffer,
            multi_fasta_specific_counts=specific_counts,
            design_name=step2_design_name_var.get().strip(),
            masking_threshold=masking,
            tiling_density=tiling,
            phyluce_prefix="conda run -n phyluce-1.7.3",
            genome_lastz_dir=lastz_output_dir,  # tell the tail code where to write
        )



        old_stdout, old_stderr = sys.stdout, sys.stderr
        sys.stdout = sys.stderr = TextRedirector(log_text, prefix="")
        try:
            logger.info("==== STEP 2: Starting temp probes → master bait set ====")
            dupe_master = step2.run_tail_pipeline(cfg2)
            step3_probe_fasta_var.set(str(dupe_master))
            genome_conf = Path(proj) / "probe-design" / "genomes.conf"
            step3_genome_conf_var.set(str(genome_conf))
            logger.info(f"STEP 2 produced master DUPE-SCREENED bait set: {dupe_master}")
            messagebox.showinfo("Step 2 complete", "Step 2 completed successfully. Check the log for details.")
        except Exception as e:
            logger.exception("Step 2 failed")
            messagebox.showerror("Step 2 failed", str(e))
        finally:
            sys.stdout, sys.stderr = old_stdout, old_stderr

    ttk.Button(frm2_buttons, text="Run Step 2", command=run_step2).grid(row=0, column=0, padx=5)

    # -----------------------------------------------------------------
    # STEP 3 TAB: in-silico test + RAxML tree
    # -----------------------------------------------------------------
    tab3 = ttk.Frame(notebook)
    notebook.add(tab3, text="Step 3: in-silico + RAxML")

    tab3.rowconfigure(2, weight=1)
    tab3.columnconfigure(0, weight=1)

    frm3_paths = ttk.LabelFrame(tab3, text="Paths (Step 3)")
    frm3_paths.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
    frm3_paths.columnconfigure(1, weight=1)

    ttk.Label(frm3_paths, text="Project root:").grid(row=0, column=0, sticky="w")
    ttk.Entry(frm3_paths, textvariable=project_root_var, width=40).grid(
        row=0, column=1, sticky="ew", padx=2
    )
    ttk.Button(
        frm3_paths,
        text="Browse",
        command=lambda: browse_dir(project_root_var, "Select project root"),
    ).grid(row=0, column=2, padx=2)

    ttk.Label(frm3_paths, text="Base genome (outgroup):").grid(row=1, column=0, sticky="w")
    ttk.Entry(frm3_paths, textvariable=base_taxon_var, width=20).grid(
        row=1, column=1, sticky="w", padx=2
    )

    ttk.Label(frm3_paths, text="Master DUPE-SCREENED bait FASTA:").grid(row=2, column=0, sticky="w")
    entry_probe_fasta = ttk.Entry(frm3_paths, textvariable=step3_probe_fasta_var, width=40)
    entry_probe_fasta.grid(row=2, column=1, sticky="ew", padx=2)

    def browse_probe_fasta() -> None:
        f = filedialog.askopenfilename(
            title="Select master DUPE-SCREENED probe FASTA",
            filetypes=[("FASTA", "*.fa *.fasta *.fna"), ("All files", "*.*")]
        )
        if f:
            step3_probe_fasta_var.set(f)

    ttk.Button(frm3_paths, text="Browse file", command=browse_probe_fasta).grid(row=2, column=2, padx=2)

    ttk.Label(frm3_paths, text="Genome 2bit config (scaffolds):").grid(row=3, column=0, sticky="w")
    entry_genome_conf = ttk.Entry(frm3_paths, textvariable=step3_genome_conf_var, width=40)
    entry_genome_conf.grid(row=3, column=1, sticky="ew", padx=2)

    def browse_genome_conf() -> None:
        f = filedialog.askopenfilename(
            title="Select genome config (scaffolds)",
            filetypes=[("Config", "*.conf *.cfg *.ini"), ("All files", "*.*")]
        )
        if f:
            step3_genome_conf_var.set(f)

    ttk.Button(frm3_paths, text="Browse file", command=browse_genome_conf).grid(row=3, column=2, padx=2)

    frm3_params = ttk.LabelFrame(tab3, text="Key parameters (Step 3)")
    frm3_params.grid(row=1, column=0, sticky="ew", padx=5, pady=5)
    frm3_params.columnconfigure(1, weight=1)

    ttk.Label(frm3_params, text="Min taxon proportion (e.g. 0.75):").grid(row=0, column=0, sticky="w")
    ttk.Entry(frm3_params, width=6, textvariable=step3_min_taxon_prop_var).grid(
        row=0, column=1, sticky="w", padx=2
    )

    ttk.Label(frm3_params, text="Align cores (MAFFT/Gblocks):").grid(row=1, column=0, sticky="w")
    ttk.Entry(frm3_params, width=6, textvariable=step3_align_cores_var).grid(
        row=1, column=1, sticky="w", padx=2
    )

    ttk.Label(frm3_params, text="LASTZ cores (in-silico):").grid(row=2, column=0, sticky="w")
    ttk.Entry(frm3_params, width=6, textvariable=step3_lastz_cores_var).grid(
        row=2, column=1, sticky="w", padx=2
    )

    ttk.Label(frm3_params, text="Flank size (bp around UCE):").grid(row=3, column=0, sticky="w")
    ttk.Entry(frm3_params, width=6, textvariable=step3_flank_var).grid(row=3, column=1, sticky="w", padx=2)

    ttk.Label(frm3_params, text="RAxML threads:").grid(row=4, column=0, sticky="w")
    ttk.Entry(frm3_params, width=6, textvariable=step3_raxml_threads_var).grid(
        row=4, column=1, sticky="w", padx=2
    )

    frm3_buttons = ttk.Frame(tab3)
    frm3_buttons.grid(row=3, column=0, sticky="e", padx=5, pady=5)

    def run_step3() -> None:
        proj = project_root_var.get().strip()
        base = base_taxon_var.get().strip()
        probe_fasta_path = step3_probe_fasta_var.get().strip()
        genome_conf_path = step3_genome_conf_var.get().strip()

        if not proj or not base or not probe_fasta_path or not genome_conf_path:
            messagebox.showerror(
                "Missing input",
                "Please set project root, base genome, probe FASTA, and genome config."
            )
            return

        try:
            min_prop = float(step3_min_taxon_prop_var.get().strip() or "0.75")
            align_cores = int(step3_align_cores_var.get().strip() or "12")
            lastz_cores = int(step3_lastz_cores_var.get().strip() or "16")
            flank_val = int(step3_flank_var.get().strip() or "400")
            rax_threads = int(step3_raxml_threads_var.get().strip() or "10")
        except ValueError:
            messagebox.showerror(
                "Invalid input",
                "Step 3 numeric parameters must be valid numbers."
            )
            return

        proj_path = Path(proj).resolve()
        pf_path = Path(probe_fasta_path).resolve()
        gc_path = Path(genome_conf_path).resolve()

        try:
            probe_rel = os.path.relpath(str(pf_path), str(proj_path))
            genome_conf_rel = os.path.relpath(str(gc_path), str(proj_path))
        except ValueError:
            probe_rel = pf_path.name
            genome_conf_rel = gc_path.name

        cfg3 = step3.TailConfig(
            project_root=str(proj_path),
            base_genome=base,
            probe_fasta=probe_rel,
            genome_conf=genome_conf_rel,
            min_taxon_proportion=min_prop,
            align_cores=align_cores,
            lastz_cores=lastz_cores,
            flank=flank_val,           
            raxml_threads=rax_threads,
            phyluce_prefix="conda run -n phyluce-1.7.3",
        )

        old_stdout, old_stderr = sys.stdout, sys.stderr
        sys.stdout = sys.stderr = TextRedirector(log_text, prefix="")
        try:
            logger.info("==== STEP 3: Starting in-silico test + RAxML tree ====")
            step3.run_insilico_pipeline(cfg3)
            messagebox.showinfo(
                "Step 3 complete",
                "Step 3 completed successfully. Final tree should be under the project directory. Check the log for details."
            )
        except Exception as e:
            logger.exception("Step 3 failed")
            messagebox.showerror("Step 3 failed", str(e))
        finally:
            sys.stdout, sys.stderr = old_stdout, old_stderr

    ttk.Button(frm3_buttons, text="Run Step 3", command=run_step3).grid(row=0, column=0, padx=5)


    #######################
    ###### run step 4 #####
    #######################

    def run_step4():
        """Run post-UCE / post-phylogeny analysis (Step 4)."""
        proj = project_root_var.get().strip()
        if not proj:
            messagebox.showerror(
                "Missing project root",
                "Please select or enter a project root (Step 1) first."
            )
            return

        proj_path = Path(proj).expanduser().resolve()

        # Read optional overrides from GUI
        insilico_override = step4_insilico_dir_var.get().strip()
        genomes_override = step4_genomes_conf_var.get().strip()
        output_override = step4_output_dir_var.get().strip()

        # Default insilico dir if none chosen
        if insilico_override:
            insilico_dir = Path(insilico_override).expanduser().resolve()
        else:
            insilico_dir = proj_path / "probe-design-test" / "taxon-sets" / "insilico-incomplete"

        matrix_file = insilico_dir / "insilico-incomplete.incomplete"

        # Default genomes.conf if none chosen
        if genomes_override:
            genomes_conf = Path(genomes_override).expanduser().resolve()
        else:
            genomes_conf = proj_path / "probe-design" / "genomes.conf"

        # Default output_dir if none chosen
        if output_override:
            output_dir = Path(output_override).expanduser().resolve()
        else:
            output_dir = insilico_dir

        # Basic validation
        if not insilico_dir.is_dir():
            messagebox.showerror(
                "Missing insilico directory",
                f"Could not find directory:\n{insilico_dir}"
            )
            return

        if not matrix_file.exists():
            messagebox.showerror(
                "Missing matrix",
                f"Could not find {matrix_file}.\n\n"
                "Make sure Step 3 (in-silico + RAxML) has been run, "
                "or select the correct insilico directory."
            )
            return

        if not genomes_conf.exists():
            messagebox.showerror(
                "Missing genomes.conf",
                f"Could not find genomes.conf at:\n{genomes_conf}"
            )
            return

        # Create output dir if needed
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            messagebox.showerror(
                "Output directory error",
                f"Could not create or access output directory:\n{output_dir}\n\n{e}"
            )
            return

        try:
            min_loci = int(step4_min_loci_var.get().strip() or "100")
        except ValueError:
            messagebox.showerror(
                "Invalid minimum loci",
                "Minimum loci per taxon must be an integer (e.g. 100)."
            )
            return

        old_stdout, old_stderr = sys.stdout, sys.stderr
        sys.stdout = sys.stderr = TextRedirector(log_text, prefix="")
        try:
            logger.info("==== STEP 4: Post-UCE / post-phylogeny analysis ====")
            step4.run_post_uce_analysis(
                project_root=str(proj_path),
                min_loci=min_loci,
                logger=logger,
                insilico_dir=str(insilico_dir),
                genomes_conf=str(genomes_conf),
                output_dir=str(output_dir),
            )
            messagebox.showinfo(
                "Step 4 complete",
                "Post-UCE analysis finished.\n\n"
                f"Check the output directory ({output_dir}) for:\n"
                "  * taxon_locus_counts.png\n"
                "  * locus_occupancy_hist.png\n"
                "  * genomes.min<min_loci>.conf"
            )
        except Exception as e:
            logger.exception("Step 4 failed")
            messagebox.showerror("Step 4 failed", str(e))
        finally:
            sys.stdout, sys.stderr = old_stdout, old_stderr


    # ------------------------------------------------------------------
    # STEP 4 TAB – Post-UCE / Post-phylogeny analysis
    # ------------------------------------------------------------------
    tab4 = ttk.Frame(notebook)
    notebook.add(tab4, text="Step 4: Post-UCE/Tree")

    frm4 = ttk.Frame(tab4, padding=10)
    frm4.pack(fill="both", expand=True)

    frm4_params = ttk.LabelFrame(frm4, text="Post-UCE / matrix analysis", padding=10)
    frm4_params.pack(fill="x", pady=(0, 10))
    frm4_params.columnconfigure(1, weight=1)

    # Row 0: minimum loci
    ttk.Label(
        frm4_params,
        text="Minimum loci per taxon to keep:"
    ).grid(row=0, column=0, sticky="w")
    ttk.Entry(frm4_params, width=6, textvariable=step4_min_loci_var).grid(
        row=0, column=1, sticky="w", padx=4
    )

    # Row 1: insilico-incomplete dir
    ttk.Label(
        frm4_params,
        text="insilico-incomplete directory:"
    ).grid(row=1, column=0, sticky="w")
    ttk.Entry(frm4_params, textvariable=step4_insilico_dir_var).grid(
        row=1, column=1, sticky="we", padx=4
    )
    ttk.Button(
        frm4_params,
        text="Browse...",
        command=lambda: browse_dir(step4_insilico_dir_var, "Select insilico-incomplete directory")
    ).grid(row=1, column=2, sticky="w", padx=4)

    # Row 2: genomes.conf
    ttk.Label(
        frm4_params,
        text="genomes.conf path:"
    ).grid(row=2, column=0, sticky="w")
    ttk.Entry(frm4_params, textvariable=step4_genomes_conf_var).grid(
        row=2, column=1, sticky="we", padx=4
    )
    ttk.Button(
        frm4_params,
        text="Browse...",
        command=lambda: browse_file(
            step4_genomes_conf_var,
            "Select genomes.conf",
            [("Config files", "*.conf"), ("All files", "*.*")]
        ),
    ).grid(row=2, column=2, sticky="w", padx=4)

    # Row 3: output dir
    ttk.Label(
        frm4_params,
        text="Output directory:"
    ).grid(row=3, column=0, sticky="w")
    ttk.Entry(frm4_params, textvariable=step4_output_dir_var).grid(
        row=3, column=1, sticky="we", padx=4
    )
    ttk.Button(
        frm4_params,
        text="Browse...",
        command=lambda: browse_dir(step4_output_dir_var, "Select output directory")
    ).grid(row=3, column=2, sticky="w", padx=4)

    # Row 4: explanatory text
    ttk.Label(
        frm4_params,
        text=(
            "This step uses insilico-incomplete.incomplete to:\n"
            "  >:|X--> Count loci per taxon\n"
            "  >:|X--> Plot loci-per-taxon and locus-occupancy\n"
            "  >:|X--> Write an adjusted genomes.conf with only well-sampled taxa"
        ),
        justify="left"
    ).grid(row=4, column=0, columnspan=3, sticky="w", pady=(6, 0))

    frm4_buttons = ttk.Frame(frm4)
    frm4_buttons.pack(anchor="w")

    ttk.Button(
        frm4_buttons,
        text="Run post-UCE / tree analysis",
        command=run_step4
    ).grid(row=0, column=0, padx=5)



    root.mainloop()


if __name__ == "__main__":
    launch_master_gui()
