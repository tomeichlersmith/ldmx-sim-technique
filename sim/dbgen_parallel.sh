#!/bin/bash

source /export/scratch/users/eichl008/ldmx/dark-brem-lib-gen/env.sh
dbgen use /export/scratch/users/eichl008/ldmx/sim-technique/tomeichlersmith_dark-brem-lib-gen_v4.2.sif
dbgen cache /export/scratch/users/eichl008/ldmx/sim-technique/.singularity
dbgen work /export/scratch/users/eichl008/ldmx/sim-technique/.dbgen-scratch
dbgen dest /export/scratch/users/eichl008/ldmx/sim-technique/dblib/scaling/
dbgen run --max_energy $1 --min_energy $1 ${@:2}
