# PVMSim GUI Workflow

## Purpose
PVMSim provides a graphical workflow for parameter extraction of the double-diode photovoltaic model from measured current-voltage curves. The GUI is intended for interactive inspection, controlled single-run execution, repeated runs with seed variation, and manual export of logs, plots, and summary tables.

## Expected project structure
The GUI is designed around the following repository layout:

```text
PVMSim/
│
├── run_main.m                         % reproducible CLI entry point
├── build_release.m                    % script used to generate the MATLAB toolbox
├── toolbox_identifier.txt             % permanent MATLAB toolbox identifier
├── README.md                          % repository overview and instructions
├── CITATION.cff                       % software citation metadata
├── LICENSE                            % PVMSim MIT license
├── LICENSE-CC0.txt                    % CC0 license for applicable metadata/content
├── THIRD_PARTY_NOTICES.md             % third-party attribution notices
├── .gitignore                         % generated-file exclusion rules
├── layout.txt                         % repository structure reference
│
├── app/
│   ├── PVMSim.mlapp                   % graphical user interface
│   ├── cancel.svg
│   ├── check.svg
│   ├── copy.svg
│   ├── download.svg
│   ├── export.svg
│   ├── file.svg
│   ├── folderopen.svg
│   ├── import.svg
│   ├── openfile.svg
│   ├── printer.svg
│   ├── reload.svg
│   ├── reset.svg
│   ├── run.svg
│   ├── run2.svg
│   ├── save.svg
│   ├── setting.svg
│   ├── upload.svg
│   └── validate.svg
│
├── config/
│   └── modules/
│       ├── STM6-40_36.txt
│       ├── PWP201.txt
│       └── STP6-120_36.txt
│
├── core/
│   ├── RS7MSA.m
│   ├── RS7MSA_model_current.m
│   ├── objective_rmse_obj.m
│   ├── metrics_kpis.m
│   ├── mutate_RS7MSA_phase1.m
│   ├── mutate_RS7MSA_phase2.m
│   ├── clamp_val.m
│   └── utils_rng.m
│
├── data/
│   └── iv/
│       ├── CurvasIV_STM6-40_36.txt
│       ├── CurvasIV_PWP201.txt
│       └── CurvasIV_STP6-120_36.txt
│
├── docs/
│   ├── Documentation.md
│   ├── workflow.md                    % GUI workflow and usage guide        
│   └── <user_manual_files>            % manual PDF, source, and associated figures
│
├── figures/
│   └── <generated_or_manual_figures>
│
├── outputs/
│   ├── exports/
│   │   ├── .gitkeep
│   │   └── <GUI_exported_files>
│   └── runs/
│       ├── .gitkeep
│       └── <datetime_model_curve_seed>/
│           ├── run_config.json
│           ├── log.txt
│           ├── run_results.mat
│           ├── run_summary.csv
│           └── checksums.sha256
│
├── tests/
│   └── smoke_test.m                   % automated fixed-seed regression test
│
├── third_party/
│   └── licenses/
│       └── uiw-icons-MIT.txt
│
└── Release/
    └── PVMSim.mltbx                   % installable MATLAB toolbox, version 0.1.0
```

The recommended convention is:
- module definition files in `config/modules`
- measured I-V curves in `data/iv`
- manual GUI exports in `outputs/exports`
- reproducible CLI run packages in `outputs/runs`
- the installable MATLAB toolbox in `Release/PVMSim.mltbx`
- the user manual and its source files in `Documentation`

## Input files

### 1. PV module definition
The module definition is loaded from a plain-text `.txt` file containing MATLAB-style assignments such as:

```text
manufacturer = 'Example manufacturer';
model = 'Example module';
Ns = 36;
Isc = 3.87;
Voc = 21.24;
Ki = 0.002;
Tref_C = 25;
```

The GUI parses these values and maps them into the internal panel structure used by the optimization core.

Before proceeding, verify that `Isc`, `Ns`, and `Tref` contain physically valid values, since these entries are especially relevant for configuring the parameter-extraction setup correctly.

### 2. Measured I-V curve
Measured curves are imported from plain-text `.txt` files containing voltage-current pairs. The GUI also supports optional metadata headers, for example:

```text
% ID=STM6_001
% G=1000
% T=25
% Date=2026-03-01
0.0,3.80
0.5,3.79
...
```

Multiple curves can be imported in one action. Duplicate curve identifiers are skipped.

## GUI layout
The interface is organized into a persistent left-side project workspace and four main analysis tabs.

### Persistent project workspace
This area contains:
- project root selection
- I-V folder selection
- module folder selection
- datasheet table
- path validation controls
- project save/import actions
- application run log

This workspace is used to define the data context before running the model.

### Main tabs
1. **Measured I-V Data**
2. **Algorithm**
3. **Plots & Metrics**
4. **Batch Runs**

## Recommended step-by-step workflow

## Installation options

PVMSim can be used directly from the repository or installed as a MATLAB toolbox.

### Option 1. Run from the repository
Set the MATLAB Current Folder to the repository root and open:

```text
app/PVMSim.mlapp
```

The command-line workflow is available through:

```matlab
run_main;
```

### Option 2. Install the MATLAB toolbox
Install the packaged toolbox:

```text
Release/PVMSim.mltbx
```

The current toolbox release is PVMSim `0.1.0` and uses the permanent MATLAB toolbox identifier stored in:

```text
toolbox_identifier.txt
```

After installation, launch PVMSim from the MATLAB Apps gallery and run **Tools > Run Smoke Test** to verify the installation.


### Step 1. Set the project folders
Use the folder-selection controls in the left-side project workspace:
- **Data Folder Browse** sets the project root
- **I-V Curve Folder Browse** sets the measured-curve directory
- **PV Module Folder Browse** sets the module-definition directory

When the selected project root follows the expected layout, the app auto-detects:
- `data/iv`
- `config/modules`

The I-V and PV module folders only need to be selected manually when the default locations must be overridden.

### Step 2. Load the PV module definition
Press **Import** in the project workspace and select one module `.txt` file.

The GUI will:
- parse the file
- normalize the imported values
- populate the datasheet table
- report the action in the application log

If parsing fails, the app reports an error and keeps the current state unchanged.

### Step 3. Import measured I-V curves
Press **Import** in the **Measured I-V Data** tab and select one or more measured-curve `.txt` files.

For each file, the GUI will:
- read metadata and numeric samples
- perform minimal format validation
- skip invalid files
- skip duplicate IDs
- store valid curves in the internal session state
- populate the measured-curve table

### Step 4. Select a measured curve
Click a row in the measured-curve table.

The GUI will update:
- the measured I-V plot
- the measured P-V plot
- the key electrical points table

The key points derived from the selected curve include:
- maximum power point
- open-circuit voltage
- short-circuit current
- fill factor

Efficiency is shown only if enough information is available. If no curve is explicitly selected, the GUI uses the first imported curve.

### Step 5. Validate the project state
Press **Validate** before running the algorithm.

Validation checks the configured folders and enables:
- **Run Algorithm**
- **Run Batch**

If validation fails, both execution buttons remain disabled.

This validation step does not replace the internal checks performed again at execution time. The app still verifies the availability of the required measured-curve and datasheet information before starting a run.


### Step 6. Verify the installation with the automated smoke test
Open **Tools > Run Smoke Test** to run the fixed-seed verification routine.

The smoke test:
- runs the STM6-40/36 benchmark with seed `42`
- verifies the reference objective value `RMSE_obj = 1.728401820819e-03 A` within the configured tolerance
- temporarily generates and verifies the required run artifacts
- checks the expected source-code entries in `checksums.sha256`

A successful run displays a green confirmation alert. A failed run stops at the first unsuccessful check, displays an error alert, and records the detailed diagnostic message in the application log.

The same test can be executed from the MATLAB command window while the Current Folder is set to the PVMSim project root:

```matlab
addpath(fullfile(pwd,'tests'));
smoke_test;
```

The smoke test is intended to verify the installation and detect unintended changes to the benchmark workflow. It is not a substitute for the multi-run statistical analysis used in performance comparisons.

### Step 7. Configure the algorithm controls
In the **Algorithm** tab, set:
- **Random seed**
- **Maximum iterations**

The current default values are:
- `Seed = 42`
- `Max Iteration = 80000`

The `Zoom x axes` control only changes the visible range of the convergence plot after the run.

### Step 8. Run a single extraction
Press **Run Algorithm**.

The GUI then:
1. checks that a measured curve is available
2. uses the selected curve, or defaults to the first loaded curve
3. builds the panel structure from the datasheet table
4. reads the seed and iteration budget from the GUI controls
5. calls the shared optimization core `RS7MSA(...)`
6. updates the extracted-parameter table
7. computes and displays KPI metrics
8. updates convergence and comparison plots
9. writes execution messages to the algorithm log

The single-run execution is interactive and designed for inspection rather than headless reproducibility packaging.

## Optimization logic used by the GUI
The GUI uses the same core routine as the command-line mode:
- `RS7MSA.m`

The current implementation follows a staged strategy:
1. random-search exploration
2. simulated-annealing refinement

The GUI passes the following run controls to the solver:
- `Seed`
- `MaxIter`

## Outputs visible in the GUI
After a successful single run, the GUI presents:
- extracted DDM parameters
- I-V and P-V measured-versus-calculated plots
- residual plot
- current and power error metrics
- convergence history
- algorithm log
- smoke-test status and diagnostics in the application log

The KPI table includes:
- RMSE for current and power
- MAE for current and power
- IAE for current and power
- relative error at maximum power point

## Batch workflow
The **Batch Runs** tab is used to quantify run-to-run variability.

### Batch inputs
The batch mode uses:
- the current module definition
- the selected measured curve
- the number of runs
- the current seed
- optional consecutive-seed mode

### Batch execution
When **Run Batch** is pressed, the GUI executes repeated calls to the optimization core and records:
- run index
- seed
- objective RMSE
- runtime
- completion status

The batch summary table reports:
- best RMSE
- worst RMSE
- mean RMSE
- standard deviation
- average runtime per run

Batch execution can be canceled through **Cancel Batch**.

## Save, import, and session management
The project workspace also provides:
- **Save** to store a lightweight project configuration with paths and algorithm settings
- **Import** to reload the project state
- **Reset** to clear the datasheet and configured folders
- **Reset I-V** to clear all imported curves and key-point displays
- **Reset Params** to clear the extracted-parameter view

These actions support session continuity, but they are not equivalent to the standardized CLI run package.

## Export options available in the GUI
The GUI supports manual export of user-facing artifacts, including:
- application run log
- algorithm log
- extracted parameters
- KPI table
- batch summary
- screenshots and plots through MATLAB axes export tools

Manual GUI exports are written to `outputs/exports`. Generated contents of `outputs/exports` and `outputs/runs` are excluded from version control, while `.gitkeep` files preserve the empty directory structure in the repository.

The GUI can also open the workflow guide from the Help menu.

## Current reproducibility scope of the GUI
The GUI supports controlled reruns of the optimization core through explicit seed and iteration controls. This is useful for interactive reproducibility at the solver level.

However, the current GUI behavior differs from the command-line mode in one important respect:
- normal GUI runs are primarily interactive and export-oriented
- the GUI can launch the automated smoke test, which internally exercises the reproducible command-line workflow and verifies its run artifacts
- the CLI remains the main interface for generating complete standardized run packages for user-selected studies

Therefore, the GUI should be interpreted as:
- a validated interactive front end for controlled execution and inspection
- not yet the main source of the full standardized reproducibility artifact package described for the CLI

## Recommended use
Use the GUI when the goal is to:
- inspect imported data
- verify module-to-curve consistency
- tune seed and iteration settings interactively
- analyze convergence and fit quality visually
- assess run-to-run variability in batch mode
- export tables, logs, and figures manually

Use the CLI when the goal is to:
- generate structured run folders automatically
- archive configuration snapshots and logs systematically
- support traceable reruns with a standardized artifact package
- perform headless or scripted studies

## Minimal GUI workflow checklist
1. Select the project root and, if needed, override the default folders.
2. Import the PV module definition.
3. Import one or more measured I-V curves.
4. Select the target curve.
5. Validate the workspace.
6. Run the automated smoke test when verifying a new installation or software release.
7. Set the seed and maximum iterations.
8. Run the algorithm.
9. Inspect parameters, convergence, plots, metrics, and logs.
10. Run batch analysis if variability assessment is needed.
11. Export logs, metrics, tables, and figures.

## Help menu integration
The GUI Help menu points to:

```text
https://github.com/pvmsim/docs/workflow.md
```

To keep the in-app documentation consistent, this file should be stored in the repository under:

```text
docs/workflow.md
```

## Versioned software release
This workflow corresponds to PVMSim `v0.1.0`. The Git tag, GitHub release, GUI version string, `CITATION.cff`, command-line metadata, MATLAB toolbox version, and archived software record should all identify software version `0.1.0`. The App Designer sharing metadata may display `0.1` because that field accepts only a `Major.Minor` format.

The installable package is:

```text
Release/PVMSim.mltbx
```

Its permanent MATLAB toolbox identifier is stored in:

```text
toolbox_identifier.txt
```

The same identifier must be retained for future PVMSim toolbox updates.

The user manual follows its own document-version sequence and identifies explicitly the PVMSim version that it documents. For example, PVMSim User Manual `v1.1` may document PVMSim software `v0.1.0`.

When a version-specific DOI is assigned to the archived software release, use that DOI to cite the exact code version associated with the manuscript. Do not replace the version-specific DOI with only the generic GitHub repository URL or with the DOI assigned to the user manual.

## Release package maintenance
The MATLAB toolbox is generated with `build_release.m`. Before packaging a new release:

1. confirm that `toolbox_identifier.txt` contains the established toolbox identifier;
2. update the software version consistently;
3. run the automated smoke test;
4. regenerate `Release/PVMSim.mltbx`;
5. install the generated toolbox in a clean MATLAB environment;
6. rerun the GUI smoke test and the command-line benchmark.

The generated toolbox should not include local files from `outputs/runs` or `outputs/exports`.
