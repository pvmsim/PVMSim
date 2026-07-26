# PVMSim: A MATLAB App for Reproducible Double-Diode Parameter Extraction from Measured Photovoltaic Current–Voltage Curves

PVMSim is a MATLAB application for reproducible parameter extraction of the double-diode photovoltaic model from measured current–voltage (I–V) curves. It combines an interactive App Designer interface with a scriptable command-line workflow for controlled single runs, batch execution, structured run exports, and installation verification through an automated smoke test.

## Current Release

- **PVMSim software version:** `0.1.0`
- **MATLAB toolbox version:** `0.1.0`
- **Tested MATLAB release:** `R2025b`
- **Tested operating system:** Windows 11
- **Additional MATLAB toolboxes:** Not required
- **User manual version:** `1.1`, documenting PVMSim `0.1.0`

## Authors

Liomnis Osorio <sup>a,b,c*</sup>  
Mailyn Moreno-Espino <sup>d,e*</sup>  
Laurent Duchêne <sup>c</sup>  
Víctor Tuninetti <sup>f</sup>  
Carlos Zalazar <sup>b</sup>  
Rodrigo Irarrázaval <sup>b</sup>  
Yoalbys Retirado-Mediaceja <sup>g,h</sup>  
Marco Rivera <sup>i,j</sup>  

<sup>*</sup> Corresponding authors

## Affiliations

<sup>a</sup> Department of Industrial Processes, Faculty of Engineering, Universidad Católica de Temuco, Temuco, Chile  
<sup>b</sup> Doctoral Program in Engineering, MacroFaculty of Engineering (UFRO–UBB–UTALCA Consortium), Chile  
<sup>c</sup> ArGEnCo Department, MSM Team, University of Liège, Liège, Belgium  
<sup>d</sup> Faculty of Informatics, Universidad Complutense de Madrid, Madrid, Spain  
<sup>e</sup> Institute of Knowledge Technology, Universidad Complutense de Madrid, Madrid, Spain  
<sup>f</sup> Department of Mechanical Engineering, Universidad de La Frontera, Temuco, Chile  
<sup>g</sup> Universidad de Moa Dr. Antonio Núñez Jiménez, Moa, Cuba  
<sup>h</sup> Section of Technical Sciences, Cuban Academy of Sciences, Havana, Cuba  
<sup>i</sup> Energy Conversion and Power Electronics Laboratory, Universidad de Talca, Curicó, Chile  
<sup>j</sup> Power Electronics, Machines and Control Research Institute, Faculty of Engineering, University of Nottingham, Nottingham, United Kingdom

## Overview

PVMSim supports traceable and repeatable estimation of double-diode model parameters from measured I–V data. Users load measured curves, select a photovoltaic module definition, configure the random seed and iteration budget, and execute a staged optimization workflow. The software provides fitted parameters, convergence information, measured-versus-calculated plots, error metrics, batch statistics, structured run artifacts, and SHA-256 integrity hashes.

## Key Features

- Interactive MATLAB App Designer graphical interface
- Scriptable command-line execution
- Controlled single-run and batch workflows
- Explicit random-seed and iteration controls
- Configurable photovoltaic module definitions
- Import of measured I–V datasets
- Export of logs, summaries, MATLAB results, and integrity hashes
- Automated fixed-seed smoke test
- Installable MATLAB toolbox package
- Versioned citation metadata through `CITATION.cff`

## Repository Structure

```text
PVMSim/
├── run_main.m
├── package_toolbox.m
├── toolbox_identifier.txt
├── README.md
├── CITATION.cff
├── LICENSE
├── LICENSE-CC0.txt
├── THIRD_PARTY_NOTICES.md
├── app/
├── config/
├── core/
├── data/
├── docs/
├── figures/
├── outputs/
├── tests/
├── third_party/
└── Release/
```

See [`layout.txt`](layout.txt) for the detailed repository hierarchy.

## Installation

### Option 1: Install the MATLAB toolbox

Install:

```text
Release/PVMSim.mltbx
```

In MATLAB:

1. Double-click `PVMSim.mltbx`, or select **Home > Add-Ons > Install from File**.
2. Complete the installation.
3. Open the MATLAB **APPS** tab.
4. Launch **PVMSim**.

After installation, select:

```text
Tools > Run Smoke Test
```

A successful smoke test confirms that the benchmark workflow, run artifacts, JSON configuration, and checksum entries are working as expected.

### Option 2: Run from source

1. Clone or download the repository.
2. Set the MATLAB Current Folder to the repository root.
3. Open:

```text
app/PVMSim.mlapp
```

The app adds the required project folders during execution.

### Option 3: Command-line execution

Run the default benchmark:

```matlab
out = run_main;
```

Run an explicit single extraction:

```matlab
out = run_main( ...
    "module", "config/modules/STM6-40_36.txt", ...
    "iv", "data/iv/CurvasIV_STM6-40_36.txt", ...
    "seed", 42, ...
    "maxIter", 80000);
```

Run a batch with consecutive seeds:

```matlab
out = run_main( ...
    "module", "config/modules/STM6-40_36.txt", ...
    "iv", "data/iv/CurvasIV_STM6-40_36.txt", ...
    "seed", 42, ...
    "maxIter", 80000, ...
    "batchN", 30, ...
    "consecutiveSeeds", true);
```

## Automated Smoke Test

The smoke test can be launched from the GUI through:

```text
Tools > Run Smoke Test
```

or from the MATLAB command window while the Current Folder is set to the repository root:

```matlab
addpath(fullfile(pwd, "tests"));
smoke_test;
```

The test runs the STM6-40/36 benchmark with seed `42` and verifies the reference objective value:

```text
RMSE_obj = 1.728401820819e-03 A
```

It also verifies the required run artifacts, the validity of `run_config.json`, and the expected entries in `checksums.sha256`.

## Input Files

### PV module definitions

Module files are stored in:

```text
config/modules/
```

They contain MATLAB-style assignments describing the electrical and thermal quantities required by the model.

### Measured I–V curves

Measured curves are stored in:

```text
data/iv/
```

Each file contains voltage–current pairs and may include metadata such as curve identifier, irradiance, temperature, and date.

## Outputs and Reproducibility Artifacts

Command-line executions create timestamped directories under:

```text
outputs/runs/
```

A standard single-run package includes:

```text
run_config.json
log.txt
run_summary.csv
run_results.mat
checksums.sha256
```

Manual exports from the graphical interface are written to:

```text
outputs/exports/
```

Generated run and export contents are excluded from version control. The directory structure is preserved with `.gitkeep` files.

## Documentation

- [`docs/workflow.md`](docs/workflow.md): GUI workflow, installation options, smoke-test procedure, and reproducibility scope
- `docs/user_manual/`: complete user manual, LaTeX source, and associated figures
- [`Release/README.md`](Release/README.md): toolbox installation and verification guide

## Toolbox Packaging

The MATLAB toolbox is generated with:

```text
package_toolbox.m
```

The permanent toolbox identifier is stored in:

```text
toolbox_identifier.txt
```

This identifier must remain unchanged in future releases so that MATLAB recognizes new packages as updates of the same toolbox.

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff).

When the versioned software release is archived in Zenodo, cite the DOI assigned specifically to PVMSim version `0.1.0`. The DOI of the user manual is a separate documentation record and should not replace the software DOI.

## License

PVMSim is distributed under the following terms:

- Source code: [MIT License](LICENSE)
- Example data and applicable non-code assets: [CC0 1.0](LICENSE-CC0.txt)
- Third-party notices: [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)

## Support

For installation or usage questions, contact:

**pvmsim.matlab@gmail.com**
