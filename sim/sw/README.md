# sim/sw
Software to run the dark brem event library model.

## Setup
Derivate of ldmx-sw so source that environment script.
```
source /full/path/to/ldmx-sw/scripts/ldmx-env.sh
```
Requires v2.0 of dev container and a specific base, so
I have written a simple env file to source.
```
cd ldmx-sim-technique
ldmx source sim/sw/ldmx-env
```

## Building
After Setup is done, we can build it similar to how ldmx-sw is built.
```
cd sim/sw/
ldmx cmake -B build -S .
cd build
ldmx make install
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
Many of the other subdirectories in this software are taken from
[ldmx-sw](https://github.com/LDMX-Software/ldmx-sw.git) and so (like DMG4)
there are many components of them that are not used here.
