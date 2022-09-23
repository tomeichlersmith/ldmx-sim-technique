#!/usr/bin/parallel --shebang-wrap /bin/bash 

source /export/scratch/users/eichl008/ldmx/sim-technique/sim/env.sh
ldmx fire ${LDMX_BASE}/sim/config.py $@
