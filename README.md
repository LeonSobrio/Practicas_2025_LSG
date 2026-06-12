# TaskAnalysis.m

**Authors:** León Sobrio García, Guiomar Niso Galán  
**Affiliation:** Cajal Institute (CSIC), Madrid, Spain  
**Date:** May 2026  
**License:** MIT

---

## Overview

`TaskAnalysis.m` is a MATLAB script that performs event-locked EEG analysis within the [Brainstorm](https://neuroimage.usc.edu/brainstorm/) framework. It is the third and final module of a three-part automated EEG pipeline, running after `QualityControl.m` and `Preprocessing.m` (authored by Natalia López, available at [https://github.com/nattloopez/Practicas2025_NLC](https://github.com/nattloopez/Practicas2025_NLC)) on the same Brainstorm protocol.

The script is designed to be fully automated and reproducible: all analysis parameters are defined in a single user configuration panel at the top of the file, and no manual interaction with the Brainstorm GUI is required.

---

## Requirements

### Software
- MATLAB (tested on R2023b and later)
- [Brainstorm](https://neuroimage.usc.edu/brainstorm/) (must be installed and on the MATLAB path)

### MATLAB Toolboxes
- Signal Processing Toolbox
- Image Processing Toolbox
- Computer Vision Toolbox
- Statistics and Machine Learning Toolbox

### Pipeline dependencies
This script must be run **after** `QualityControl.m` and `Preprocessing.m`. It expects a Brainstorm protocol containing preprocessed EEG files tagged with the string defined in `ProcessedTag` (default: `'Average_reference'`).

---

## What the script does

For each participant and task specified in the configuration panel, the script performs the following steps in order:

1. **Event renaming** (optional): Maps raw BIDS trigger codes to semantic event names using a user-defined lookup table.
2. **Event combining** (optional): Creates compound events from stimulus-response pairs (e.g., correct vs. incorrect Go trials) using `process_evt_combine`.
3. **Epoching**: Segments the continuous preprocessed recording into epochs around the specified events, using a user-defined time window and baseline correction.
4. **Bad trial exclusion**: Epochs already marked as bad by `Preprocessing.m` are automatically excluded.
5. **Averaging**: Computes ERPs per event type (individual averages) and/or per user-defined condition groups (grouped averages), using arithmetic mean by default.
6. **Head model computation**: Computes a forward model on the cortical surface using a 3-shell sphere (Berg) head model.
7. **Noise covariance estimation**: Estimates the noise covariance matrix from the baseline segment of all accepted epochs.
8. **Inverse solution**: Computes a shared minimum norm (MN) inverse kernel using `process_inverse_2018` and applies it to all averages, yielding distributed cortical source estimates.
9. **Report generation**: Saves one HTML report and one JSON file per (subject, task) pair in the specified reports directory.

---

## Configuration

All parameters are set in the **USER CONFIGURATION PANEL** section at the top of the script (lines 30-209). Key parameters include:

| Parameter | Description |
|---|---|
| `BrainstormDbDir` | Path to the Brainstorm database directory |
| `ReportsDir` | Output directory for HTML and JSON reports |
| `ProtocolName` | Name of the Brainstorm protocol to process |
| `Subs` | List of subject IDs to process (empty = all) |
| `TaskNames` | List of BIDS task names to process |
| `ProcessedTag` | String identifying the final preprocessed file |
| `EventRename` | Mapping from raw trigger codes to semantic names |
| `EventCombine` | Definition of compound stimulus-response events |
| `EpochEvents` | Comma-separated list of events to epoch around |
| `EpochTime` | Epoch time window in seconds (e.g., `[-0.2, 0.5]`) |
| `BaselineWin` | Baseline correction window in seconds |
| `AvgMode` | Averaging mode: `'individual'`, `'grouped'`, or `'both'` |
| `AvgGroups` | Struct defining custom event groups for averaging |
| `headModel_eeg` | EEG head model type (default: 3-shell sphere) |
| `sources_inverseMethod` | Inverse method (default: `'minnorm'`) |
| `sources_inverseMeasure` | Inverse measure (default: `'amplitude'`) |

Do **not** edit the code below the `DO NOT EDIT BELOW THIS LINE` comment.

---

## Outputs

For each (subject, task) pair, the script produces:

- Epochs in the Brainstorm database
- Individual and/or grouped ERP averages in the Brainstorm database
- Head model and noise covariance matrix in the Brainstorm database
- Cortical source estimates (shared minimum norm kernel) in the Brainstorm database
- `Task-<subject>-<condition>-<protocol>.html`: Brainstorm HTML report with snapshots
- `Task-<subject>-<condition>-<protocol>.json`: Structured JSON summary of analysis parameters and results

---

## Usage

1. Run `QualityControl.m` and `Preprocessing.m` first on the same protocol.
2. Open `TaskAnalysis.m` in MATLAB.
3. Edit the USER CONFIGURATION PANEL to match your study paths, subjects, tasks, and analysis parameters.
4. Run the script. Brainstorm will launch in no-GUI mode automatically.

---

## Citation

If you use this script, please cite:

L. Sobrio García and G. Niso Galán, "Cognitive Control and Habits: An EEG Approach," Undergraduate thesis, Universidad Carlos III de Madrid / Instituto Cajal (CSIC), Madrid, Spain, 2026.
