# sim/sw
Software to run the dark brem event library model.

## Setup
Derivate of ldmx-sw so source that environment script.
```
source /full/path/to/ldmx-sw/scripts/ldmx-env.sh
```
Requires v2.0 of dev container and a specific base
```
cd ldmx-sim-technique
ldmx clean env
ldmx base .
ldmx use dev v2.0
```

## Reference
DMG4 was downloaded using
```
wget http://mkirsano.web.cern.ch/mkirsano/DMG4.tar.gz
```
and then unpacked
```
tar xzf DMG4.tar.gz
```
