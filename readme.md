
## License Statement:  Distributed under a Creative Commons Attribution-NonCommercial 4.0 International License 
The software, including Stampy components, is provided "AS IS," without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose, and non-infringement. In no event shall the authors or copyright holders be liable for any claim, damages, or other liability, whether in an action of contract, tort, or otherwise, arising from, out of, or in connection with the software or the use or other dealings in the software. Furthermore, use of this software is strictly restricted to non-commercial purposes only (as defined by the accompanying
Creative Commons Attribution-NonCommercial 4.0 International License or similar academic license terms); any commercial use requires obtaining a separate license from the original copyright holders.

## About
This is a simple Graphical User Interface (GUI) for the PHYLUCE pipeline that tries to mirror the Tutorial IV from the PHYLUCE Readthedocs : https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-4.html

If you use this code you must also cite PHYLUCE and Stampy along with a few other helpers LASTAL 

The GUI will help students and users of all levels keep organized and make the bar to entry of doing phylogenomics less intimidating. There are also a few good things that I added such as accepting spaces in some (but not all ) file names (spaces in file names should be avoided) it also solve a memory bloom issue at the lastz alignment step by writing the sql database in a step-wise manner. Other than these small additions it is basically PHYLUCE wrapped in a GUI, that students can run on their laptops. It assumes you already have assembled genomes. OK have fun with it and I hope it helps to improve your phylogenomic workflow.

## Summary and Install
## phyluce-one-pass (GUI + wrappers)

A small GUI + a set of wrappers around **PhyLUCE Tutorial IV** (UCE probe design + in-silico test),
with extra safeguards for non-interactive runs (no `[y/n]` prompts), memory bloom avoidance and for paths containing spaces.

This repository also vendors the exact Stampy version used by the pipeline (`stampy-1.0.28/`).

## What this runs

The pipeline is very similar to and wraps major parts of the official PhyLUCE Tutorial IV workflow:
https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-4.html

Key steps include:
- Probe design / duplicate screening (LASTZ)
- In-silico test of bait design against genomes (LASTZ)
- UCE alignment (MAFFT), trimming (Gblocks), concatenation
- Tree inference (RAxML)

## Requirements

- Linux recommended
- Conda (Miniconda/Miniforge/Mambaforge)
- Three conda envs are used:
  - `phyluce-one-pass-gui`  (GUI runtime)
  - `phyluce-1.7.3`         (PHYLUCE + external bioinformatics tools)
  - `stampy_py27`           (Python 2.7 runtime for Stampy)

## Installation

```bash
git clone https://github.com/alexrvandam/PHYLUCE_GUI.git
cd PHYLUCE_GUI

# Create envs
conda env create -f gui_envs/phyluce-one-pass-gui.yml
conda env create -f gui_envs/phyluce-1.7.3.yml
conda env create -f gui_envs/stampy_py27.yml
```

Then simply open a Terminal and run:
```bash
conda activate phyluce-one-pass-gui
python phyluce_one_pass_master_gui_v3.py
```
