# PVMSim Toolbox Release

PVMSim is a MATLAB application for reproducible parameter extraction of the double-diode photovoltaic model from measured current–voltage (I–V) curves. It combines an interactive App Designer interface with a scriptable command-line workflow for controlled single runs, batch execution, and reproducibility checks.

## Release Information

- **PVMSim software version:** `0.1.0`
- **MATLAB toolbox version:** `0.1.0`
- **Tested MATLAB release:** `R2025b`
- **Tested operating system:** Windows 11
- **Architecture used for packaging:** `win64`
- **Additional MATLAB toolboxes:** Not required

## Package Contents

This release folder contains:

```text
Release/
├── PVMSim.mltbx
└── README.md
```

`PVMSim.mltbx` is the installable MATLAB toolbox package.

## Installation

1. Open MATLAB.
2. Double-click `PVMSim.mltbx`, or install it through **Home > Add-Ons > Install from File**.
3. Complete the installation dialog.
4. Open the MATLAB **APPS** tab.
5. Launch **PVMSim** from the Apps gallery.

## Installation Verification

After installation:

1. Launch PVMSim.
2. Select **Tools > Run Smoke Test**.
3. Confirm that the smoke test finishes successfully.

The automated smoke test runs the STM6-40/36 benchmark with seed `42` and verifies:

- the expected objective value;
- the required reproducibility artifacts;
- the validity of `run_config.json`;
- the expected checksum entries.

The reference objective value is:

```text
RMSE_obj = 1.728401820819e-03 A
```

A successful execution displays a green confirmation message. A failed execution stops at the first unsuccessful check and displays the corresponding diagnostic message.

## Basic Use

After installation:

1. Select the PVMSim project root.
2. Load a PV module definition from `config/modules`.
3. Import one or more measured I–V curves from `data/iv`.
4. Select the target curve.
5. Validate the configured workspace.
6. Set the random seed and maximum iteration budget.
7. Run the staged parameter-extraction workflow.
8. Review the extracted parameters, convergence history, plots, metrics, and logs.
9. Use the batch workflow when repeated runs are required.

## Source Code and Documentation

The source code, workflow documentation, user manual, licensing information, and release metadata are available in the project repository:

<https://github.com/pvmsim>

The complete user manual is available in the repository documentation and through its archived Zenodo record.

## License

PVMSim source code is distributed under the MIT License. Example data and applicable non-code assets are distributed under CC0 1.0. Third-party resources remain subject to their respective licenses and attribution notices.

## Authors

Liomnis Osorio, Mailyn Moreno-Espino, Laurent Duchêne, Víctor Tuninetti, Carlos Zalazar, Rodrigo Irarrázaval, Yoalbys Retirado-Mediaceja, and Marco Rivera.

## Support

For installation or usage questions, contact:

**pvmsim.matlab@gmail.com**
