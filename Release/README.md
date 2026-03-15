# PVMSim: A MATLAB App for Reproducible Double-Diode Parameter Extraction from Measured Photovoltaic Current–Voltage Curves

PVMSim is a MATLAB application for reproducible parameter extraction of the double-diode photovoltaic model from measured current–voltage curves. It combines an interactive App Designer interface with a scriptable command-line workflow for exploratory analysis and headless batch execution.

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
PVMSim supports traceable and repeatable estimation of double-diode model parameters from measured current–voltage data. Users load measured I–V files, select a photovoltaic module definition from a configuration library, and execute a staged optimization workflow to estimate the model parameters. Run settings include the random seed and iteration budget. Each execution exports a configuration snapshot, logs, tabular summaries, MATLAB result files, and integrity hashes to support auditable reruns and consistent comparisons across datasets and run conditions.

## Key Features
- Interactive MATLAB App Designer graphical interface
- Scriptable command-line execution for batch and headless workflows
- Reproducible staged optimization for double-diode parameter extraction
- Support for measured I–V datasets and configurable PV module definitions
- Export of logs, summaries, result files, and integrity hashes
- Consistent rerun support through stored configuration snapshots

## Requirements
- MATLAB, tested on **MATLAB R2025b**
- No additional MATLAB toolboxes required for the current release

## Installation

### MATLAB users
Install the toolbox package `PVMSim.mltbx` in MATLAB. Once the installation is completed, PVMSim will appear in the **APPS** tab and can be launched from the MATLAB Apps gallery.

### MATLAB release
The current package was created with:

• **MATLAB (full): 25.2.0.2998904 (R2025b)**
• **Platform: Microsoft Windows 11 Pro**
• **Architecture: win64**

### Standalone version
A standalone Windows installer is not included in the current release unless explicitly provided as a separate distribution artifact.

## Getting Started
After installation:

1. Open MATLAB.
2. Go to the **APPS** tab.
3. Launch **PVMSim** from the installed apps list.
4. Load the project resources, PV module definition, and measured I–V dataset.
5. Configure the run settings, including the seed and iteration budget.
6. Execute the staged parameter extraction workflow.
7. Review the exported logs, summaries, and result files for traceable analysis.

## Repository Contents
The repository may include:
- MATLAB source files for the graphical interface and command-line workflow
- configuration files for photovoltaic modules
- example measured I–V datasets
- documentation and user guidance
- packaged toolbox release files when applicable

## Support
For installation or usage issues, contact:  
**pvmsim.matlab@gmail.com**
