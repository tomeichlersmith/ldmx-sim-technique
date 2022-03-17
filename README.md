# Sim Technique Paper Notes

## Start-Up
Production container built ~7 months ago, pulling down to UMN for testing.
```
cd ${LDMX_BASE}
singularity build ldmx_local_sim-technique-dmg4.sif docker://tomeichlersmith/ldmx-sw:sim-technique-dmg4
ldmx use local sim-technique-dmg4
```

## Documentation
Samples generated on `sim-technique` branch of ldmx-sw. 
This branch had three necessary developments for the studies related to this paper.
1. Allowing the simulation to run without any hit collections being produced.
2. Auto generation of GDML corresponding to big hunk of a single material.
3. Addition of new dark brem model using NA64's DMG4
4. (in dev) Simple Ntuplizer of particles involved in dark brem interaction

The necessary dark brem event libraries from MadGraph were generated with version 4.0 of 
[dark-brem-lib-gen](https://github.com/tomeichlersmith/dark-brem-lib-gen).

- ana : analysis scripts
- sim : config scripts for simulations with ldmx-sw

## Development
Due to DMG4's requirement of GSL, we need to use v2.0 of the ldmx/dev container.
```
ldmx use dev v2.0
```
We can install DMG4 to `${LDMX_BASE}/dmg4/install` so that CMake can find it and
then we need to symlink its libraries to the ldmx-sw install location so that
they can be found by the linker at runtime.
```
cd ldmx-sw/install
ln -s ${LDMX_BASE}/dmg4/install/lib/* .
```
