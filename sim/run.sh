#!/bin/bash

# run.sh run four simulations within container
#   Usage: ldmx run.sh {output_dir}

if [ -z $1 ]; then
  echo "Need to specify an output directory."
fi

_output_dir=$(cd $1 && pwd -P)
mkdir -p ${_output_dir}
_log=${_output_dir}/fire.log

echo "MG-Scaling Muons" | tee ${_log}
fire ${LDMX_BASE}/sim/config.py \
  ${LDMX_BASE}/dblib/muon_copper_MaxE_100.0_MinE_2.0_RelEStep_0.1_UndecayedAP_mA_1.0_run_3000/ \
  --out_dir ${_output_dir} | tee -a ${_log}

echo "DMG4 Muons" | tee -a ${_log}
fire ${LDMX_BASE}/sim/dmg4.py \
  -m brass --particle mu- --primary_energy 100. --ap_mass 1000 1 \
  --out_dir ${_output_dir} | tee -a ${_log}

echo "MG-Scaling Electrons" | tee -a ${_log}
fire ${LDMX_BASE}/sim/config.py \
  ${LDMX_BASE}/dblib/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000/ \
  --out_dir ${_output_dir} | tee -a ${_log}

echo "DMG4 Electrons" | tee -a ${_log}
fire ${LDMX_BASE}/sim/dmg4.py \
  -m tungsten --particle e- --primary_energy 4. --ap_mass 100 1 \
  --out_dir ${_output_dir} | tee -a ${_log}
