# Dark Brem Simulation within Geant4

The code necessary to replicate our studies is contained within [sw](sw) while this directory holds the higher-level interfaces for running the simulation.

### config.py
This python script is a "configuration" script for [the processing Framework](https://github.com/LDMX-Software/Framework.git) we use to run the simulation. It has many options and requires the user to choose a dark brem simulation method for it to run. The arguments are parsed with Python's argparse which means you can list them with the `-h` or `--help` command line argument.

### run.sh
This bash script is simply there for _roughly_ parallelizing the generation of the different samples. It is highly specific to the computer we ran on; however, it may be useful if you wish to replicate our results. It assumes that four processes can run in parallel and it puts the samples in an output directory named after `git describe --tags` for automatic separation of data which may change with developments.

### env.sh
This software is a derivative of [ldmx-sw](https://github.com/ldmx-software/ldmx-sw) which is compiled and run from within a [container](https://www.docker.com/resources/what-container/) using either [Docker](https://www.docker.com/) or [singularity](https://sylabs.io/singularity/#). In theory, compiling all of the necessary dependencies of this software and this software itself can be done "on bare metal" (i.e. not inside a container); nevertheless, only usage of the container is documented here since that is what was used for the studies presented in the paper.

**After** installing a container runner (either docker or singularity), you can source the provided environment script within a `bash` terminal which will download the necessary container and define a container interaction function named `ldmx` (again, inherited from this software's parent).

```bash
source env.sh
# on first source, will take some time downloading the container
ldmx help
# printout of commands defined
```

### MadGraph/MadEvent
This technique heavily relies upon MG/ME to function.
The "libraries" of events and other MG/ME samples were generated with v4.2 of [dark-brem-lib-gen]().
