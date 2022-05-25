# Sim Technique Paper Notes

## Documentation
Samples generated with [sim/sw](sim/sw).
This branch had four necessary developments for the studies related to this paper.
1. Allowing the simulation to run without any hit collections being produced.
2. Auto generation of GDML corresponding to big hunk of a single material.
3. Addition of new dark brem model using NA64's DMG4
4. Simple Ntuplizer of particle kinematics involved in dark brem interaction

The necessary dark brem event libraries from MadGraph were generated with version 4.0 of 
[dark-brem-lib-gen](https://github.com/tomeichlersmith/dark-brem-lib-gen).

### Table of Contents
- ana : analysis scripts and notebooks
- sim : code for running either DMG4 or MGS dark brem simulations
