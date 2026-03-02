# PVMSim: A MATLAB App for Reproducible Double-Diode PV Parameter Extraction

## Authors
Liomnis Osorio a, b, c*
Laurent Duchêne c
Víctor Tuninetti d
Mailyn Moreno-Espino e, f
Calos Zalazar b
Rodrigo Irarrázaval b
Yoalbys Retirado-Mediaceja h

## Affiliations
a Department of Industrial Processes, Faculty of Engineering, Universidad Católica de Temuco, Temuco, Chile
b Doctoral Program in Engineering, MacroFaculty of Engineering (UFRO–UBB–UTALCA Consortium), Chile
c ArGEnCo Department, MSM team, University of Liège, Liège, Belgium
d Department of Mechanical Engineering, Universidad de La Frontera, Temuco, Chile
e Faculty of Informatics, Universidad Complutense de Madrid, Madrid, Spain
f Institute of Knowledge Technology, Universidad Complutense de Madrid, Madrid, Spain
h Universidad de Moa, Moa, Cuba

## Overview
PVMSim is a MATLAB application for reproducible parameter extraction of the double-diode photovoltaic model from measured current–voltage curves. It provides an interactive App Designer interface and a scriptable command-line entry point for exploratory use and headless batch runs. Users load measured I–V files, select a PV module definition from a configuration library, and execute a staged optimization to estimate model parameters. Run controls include the seed and iteration budget. Each run exports a configuration snapshot, logs, tabular summaries, MATLAB results, and integrity hashes, enabling traceable reruns and consistent comparisons across I–V datasets and run settings.

## Requirements
- MATLAB (tested on MATLAB R2025b)
- No additional toolboxes required

## Installation
### Option 1: Run from source (recommended for development)
1. Clone or download this repository.
2. Add the repository root to the MATLAB path:
   - MATLAB: Home → Set Path → Add with Subfolders (or `addpath(genpath(pwd))`)
3. Open the app:
   - `app/PVMSim.mlapp`

### Option 2: Command-line (headless / batch)
Run the CLI entry point from MATLAB:

run_main("module","config/modules/<MODULE_FILE>.txt", ...
         "iv","data/iv/<IV_FILE>.txt", ...
         "seed",42,"iters",80000);

## Getting started
1. Select a PV module definition from config/modules/.
2. Load a measured I–V file from data/iv/.
3. Set the seed and iteration budget.
4. Run the staged optimization and review exported results.

## Input data format
Measured I–V files are plain text with two columns:
- Voltage [V]
- Current [A]
Example files are provided in data/iv/.

## Outputs and reproducibility artifacts
Each run creates a timestamped directory under outputs/runs/ including:
- run_config.json (configuration snapshot)
- log.txt (execution log)
- run_summary.csv (tabular summary)
- run_results.mat (MATLAB results)
- checksums.sha256 (integrity hashes)

## License
- Source code: MIT License (see `LICENSE`).
- Example data and non-code assets (if provided): CC0 1.0 (see `LICENSE-CC0.txt`).
- Third-party notices: see `THIRD_PARTY_NOTICES.md`.

## Support
Email: pvmsim.matlab@gmail.com
