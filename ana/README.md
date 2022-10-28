# Sim Technique Analysis

This is a python "package" which is not configured to actually be installed,
but is expected to be included as a unit.

Both the jupyter notebook and the condensed non-interatice plotting scripts assume that the dark brem event libraries
are stored in `../dblib` so that they can generate comparisons to a pure-MadGraph sample.
Patching them to use a different path is done within the `bundle` (or similar) function in the non-interactive plotting script
and the cells calling the `read` function in the notebook.

## Setup
You could do this inside of a Python3 virtual environment.
```
python3 -m venv .venv --prompt sim-tech-ana
source .venv/bin/activate
pip install -r requirements.txt
jupyter-notebook
```

## Table of Contents
- `dev.ipynb` : python notebook initially used to derive analysis and plotting functions
- `plot.py` : non-interactive plotting script for simulated ROOT files and cross section comparisons
- `scaling-plot.py` : non-interactive plotting script for looking at outputs of `g4db-scample`
` `dark_brem_lhe.py` : script for extracting kinematics from dark brem LHE files
