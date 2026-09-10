# ETR transient heat-transfer simulation

## Overview

This package contains the custom Python code used to simulate two-dimensional
transient heat conduction in the electrified thermal reduction (ETR) reactor.
The computational domain includes argon, a GeO2 sample layer and two
carbon-paper heater layers. The experimentally measured carbon-paper
temperature profile, `T_hot(t)`, is imposed on the heater layers as a
time-dependent Dirichlet condition.

The simulation generates the temperature fields and temperature-time profiles
used for the heat-transfer analysis reported in Fig. 1c,d and Supplementary
Fig. 21. The numerical method and model assumptions are described in
Supplementary Note 1, equations S1-S6.

## Package contents

This repository contains:

- `Supplementary code.py`: Python source code for the heat-transfer simulation.
- `raw_Time_T_every_0p1s.csv`: experimentally measured carbon-paper temperature
  profile used as the model input.
- `README.md`: system requirements, installation and execution instructions.

`Supplementary code.py` and `raw_Time_T_every_0p1s.csv` must be placed in the same
working directory when the program is executed.

## System requirements

The code runs on a standard CPU-based desktop computer. It does not require a
GPU, computing cluster or other non-standard hardware.

The code passed an import test and a shortened end-to-end simulation test in
the following environment:

- Operating system: Microsoft Windows 11, 64-bit, build 26200
- Python: 3.7.1
- NumPy: 1.21.6
- pandas: 1.3.5
- Matplotlib: 3.0.2
- openpyxl: 3.1.3

The Arial font is used for plot formatting. If Arial is unavailable,
Matplotlib may substitute another font. Font substitution does not affect the
numerical results.

## Installation

Install Python and the required packages in a virtual environment. For
example:

```text
python -m venv etr_simulation_environment
etr_simulation_environment\Scripts\activate
python -m pip install numpy pandas matplotlib
```

On macOS or Linux, activate the environment with:

```text
source etr_simulation_environment/bin/activate
```

Installation normally takes approximately 5 minutes on a desktop computer
with Python and internet access. This estimate excludes the time required to
install Python itself.

## Input data

The program reads the following file from the current working directory:

```text
raw_Time_T_every_0p1s.csv
```

The first row must contain the following exact column names:

| Column | Description | Unit/format |
| --- | --- | --- |
| `time_s` | Elapsed experimental time | seconds |
| `T_hot_C` | Measured carbon-paper temperature | degrees Celsius |

The time series should cover the full simulation interval. Values between
experimental time points are obtained by linear interpolation.

The supplied profile contains 640 records from 0.1 to 64.0 s at 0.1-s
intervals, with no missing or duplicated time points. Its measured heater
temperature ranges from 1000.0 to 1715.8 degrees Celsius. Because the first
record occurs at 0.1 s, the model holds its first recorded value (1000.0
degrees Celsius) from t = 0 to 0.1 s; the remaining values are linearly
interpolated.

## Running the simulation

Open a terminal in the directory containing the source code and input file,
then run:

```text
python "Supplementary code.py"
```

No command-line arguments are required. Progress information, the numerical
time step, probe positions, selected sample temperatures and total runtime are
printed in the terminal.

## Expected runtime

The model uses a 201 x 201 grid and approximately 7.26 million explicit time
steps for the complete 64-s simulation. A benchmark using the same grid and
numerical update on the tested Windows system gives an expected runtime of
approximately 2 hours. Runtime will vary with processor speed and available
memory.

## Expected output

The program writes the following files to the current working directory:

### Numerical data

- `heater_profile_used_every_0p1s.csv`: heater temperature used by the model at 0.1-s
  intervals.
- `simulated_sample_and_probe_Time_T_every_0p1s.csv`: heater, sample and probe
  temperatures at 0.1-s intervals.
- `probe_above_carbon_paper_Time_T_every_0p1s.csv`: heater and probe
  temperatures at positions 0.5, 1.0 and 1.5 cm above the upper carbon-paper
  surface.

### Plots

- `probe_above_carbon_paper_T_vs_time.png`
- `probe_above_carbon_paper_T_vs_time.pdf`
- `snapshot_1s.png`
- `snapshot_2s.png`
- `snapshot_5s.png`
- `snapshot_6s.png`
- `snapshot_10s.png`
- `snapshot_20s.png`
- `snapshot_30s.png`
- `snapshot_40s.png`
- `snapshot_60s.png`
- `snapshot_64s.png`

The snapshot files show the two-dimensional temperature field at the indicated
simulation times. The probe plot shows temperature histories at 0.5, 1.0 and
1.5 cm above the upper carbon-paper surface.

## Reproducing the reported simulation

To reproduce the reported heat-transfer calculation:

1. Place the supplied experimental `raw_Time_T_every_0p1s.csv` file in the same
   directory as `Supplementary code.py`.
2. Do not change the geometry, material-property, boundary-condition, mesh or
   probe-location constants defined near the beginning of the source code.
3. Run `python "Supplementary code.py"`.
4. Compare the generated CSV files and temperature-field snapshots with the
   corresponding source data and panels in Fig. 1c,d and Supplementary Fig. 21.

If separate experimental heating profiles were used for different manuscript
figures, run the program separately for each supplied profile after copying
the relevant profile to the required filename `raw_Time_T_every_0p1s.csv`.
Keep the outputs from each run in a separately labelled directory.

## Using another temperature profile or reactor configuration

To simulate another measured heating profile, replace
`raw_Time_T_every_0p1s.csv` with a CSV file containing the required `time_s`
and `T_hot_C` columns. The code also supports an Excel input file containing
`Time` and `T` columns if the `INPUT_FILE` constant is changed accordingly;
Excel input additionally requires the `openpyxl` package.

The following named constants near the beginning of `Supplementary code.py`
control the model:

- domain length and grid size;
- carbon-paper width and thickness;
- sample thickness;
- probe distances;
- thermal conductivity, density and heat capacity of each material;
- Fourier number;
- ambient temperature and effective heat-transfer coefficient;
- requested snapshot and CSV output times.

Changes to these constants describe a different simulation and should be
recorded when the code is reused.

## Reproducibility notes

The calculation is deterministic and does not use random numbers or stochastic
sampling. Repeated runs with the same input file, code and dependency versions
are expected to produce the same numerical results, apart from possible minor
floating-point or plot-rendering differences across systems.

## License and access

The source code and input data are publicly available in this repository. No
open-source software license has yet been assigned; reuse therefore requires
permission from the copyright holders unless and until a license is added.
