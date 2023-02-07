# Batch Dark Brem Library Generation

These notes were taken _after_ completing all of the studies and so
following these steps may not reproduce exactly the same results.

### condor
I used dark-brem-lib-gen to generate all dark brem libraries used
for simulation and for studies of the MG outgoing kinematics and
estimate cross section. Below is an example condor submission script that
would mimic what was done. The actual procedure was copying this script
into different sub-directories and changing the target, lepton, apmass,
and energy parameters to generate samples of interest.
`/full/path/to/shared/location` is the full path to a directory in
a filesystem that is shared across all nodes.

```
executable = /full/path/to/shared/location/dbgen.sh
transfer_executable = yes
should_transfer_files = yes
when_to_transfer_output = ON_EXIT

# terminal and condor output log files
#   this is helpful for debugging purposes but you can delete these lines
#   for slightly better performance
output = logs/$(Cluster)-$(Process).out
error  = $(output)
log    = logs/$(Cluster)-$(Process).log

# "hold" the job if the script exits with a non-zero exit code
#   this is a helpful way to list which jobs failed
#   we also store the failure-status in the hold reason sub code so you
#   can see it using condor_q
on_exit_hold = ExitCode != 0
on_exit_hold_subcode = ExitCode
on_exit_hold_reason = "Program exited with non-zero error status (stored in HoldReasonSubCode)"

arguments = "$(energy) --run $(Process) --target copper --lepton muon --apmass 1.0"
# do 5 jobs per energy point
queue 5 energy from seq 5 5 1000 |
```

where `dbgen.sh` is similar to what is outlined in the dark-brem-lib-gen README,
but we limit it to looking at a single energy at a time.

```bash
#!/bin/bash

set -ex

source /full/path/to/shared/location/dark-brem-lib-gen/env.sh
dbgen use /full/path/to/shared/location/tomeichlersmith_dark-brem-lib-gen_v4.4.sif
mkdir scratch
dbgen work scratch

# only want to do one energy at a time
dbgen run --max_energy $1 --min_energy $1 ${@:2}

# condor doesn't scan subdirectories so we should move the LHE files here
find \
  -type f \
  -name "*.lhe" \
  -exec mv {} . ';'
```

