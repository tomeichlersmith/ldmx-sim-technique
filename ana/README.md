# Sim Technique Analysis

This is mainly centered on a jupyter notebook whose dependencies are listed in [requirements.txt](requirements.txt).

Both the jupyter notebook and the condensed non-interatice plotting script assume that the dark brem event libraries
are stored in `../dblib` so that they can generate comparisons to a pure-MadGraph sample.
Patching them to use a different path is done within the `bundle` function in the non-interactive plotting script
and the cells calling the `read` function in the notebook.

## Setup
You could do this inside of a Python3 virtual environment.
```
python3 -m venv .venv --prompt sim-tech-ana
source .venv/bin/activate
pip install -r requirements.txt
jupyter-notebook
```
