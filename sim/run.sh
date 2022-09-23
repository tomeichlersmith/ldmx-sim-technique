#!/bin/bash
# run.sh 
#   run four simulations within container

__usage__() {
  cat <<HELP
 USAGE:
  ./sim/run.sh [-o OUT_DIR]

 OPTIONS:
  -o    : Base output directory for data files (default: 'data/<git describe --tags>')

 We do two simulations at once (G4DB and DMG4) and then we loop over the combinations
 of thicknesses and incident particles. The multiple jobs are written to a job listing
 file in the output directory and then "submitted" to parallel via the ldmx_parallel.sh
 script.

HELP
}

__main__() {
  local _tag=$(git describe --tags)
  local _output_dir=$(cd data && pwd -P)/${_tag}
  while [ $# -gt 0 ]; do
    case $1 in
      -o)
        if [ -z $2 ]; then
          echo "ERROR: '$1' requires an argument."
          return 1
        fi
        case $1 in
          -o) _output_dir=$2;;
        esac
        shift
        shift
        ;;
      -h|--help|-?)
        __usage__
        return 0
        ;;
    esac
  done
  
  if ! mkdir -p ${_output_dir}; then
    echo "ERROR: Could not create output directory ${_output_dir}"
    return $?
  fi

  local _muon_thin=100
  local _elec_thin=0.35
  local _elec_thick=18
  local _muon_thick=2000
  
  echo "--out_dir ${_output_dir} --depth ${_elec_thick} g4db" \
    "${LDMX_BASE}/dblib/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000/" \
    "&>> ${_output_dir}/g4db_electron_thick.log"
  
  echo "--out_dir ${_output_dir} --depth ${_elec_thick}" \
    "dmg4 -m tungsten "--particle electron --primary_energy 4. --ap_mass 100 1" \
    "&>> ${_output_dir}/dmg4_electron_thick.log"

  echo "--out_dir ${_output_dir} --depth ${_elec_thin} g4db" \
    "${LDMX_BASE}/dblib/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000/" \
    "&>> ${_output_dir}/g4db_electron_thin.log"
  
  echo "--out_dir ${_output_dir} --depth ${_elec_thin}" \
    "dmg4 -m tungsten "--particle electron --primary_energy 4. --ap_mass 100 1" \
    "&>> ${_output_dir}/dmg4_electron_thin.log"

  echo "--out_dir ${_output_dir} --depth ${_muon_thin} g4db" \
    "${LDMX_BASE}/dblib/muon_copper_MaxE_100.0_MinE_2.0_RelEStep_0.1_UndecayedAP_mA_1.0_run_3000/" \
    "&>> ${_output_dir}/g4db_muon_thin.log"
  
  echo "--out_dir ${_output_dir} --depth ${_muon_thin}" \
    "dmg4 -m brass "--particle muon --primary_energy 100. --ap_mass 1000 1" \
    "&>> ${_output_dir}/dmg4_muon_thin.log"

  echo "--out_dir ${_output_dir} --depth ${_muon_thick} g4db" \
    "${LDMX_BASE}/dblib/muon_copper_MaxE_100.0_MinE_2.0_RelEStep_0.1_UndecayedAP_mA_1.0_run_3000/" \
    "&>> ${_output_dir}/g4db_muon_thick.log"
  
  echo "--out_dir ${_output_dir} --depth ${_muon_thick}" \
    "dmg4 -m brass "--particle muon --primary_energy 100. --ap_mass 1000 1" \
    "&>> ${_output_dir}/dmg4_muon_thick.log"

  _elec_thin=1.

  echo "--out_dir ${_output_dir} --depth ${_elec_thin} g4db" \
    "${LDMX_BASE}/dblib/electron_lead_MaxE_100.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000/" \
    "&>> ${_output_dir}/g4db_electron_100GeV_lead_thin.log"
  
  echo "--out_dir ${_output_dir} --depth ${_elec_thin}" \
    "dmg4 -m lead --particle electron --primary_energy 100. --ap_mass 100 1" \
    "&>> ${_output_dir}/dmg4_electron_100GeV_lead_thin.log"

  _elec_thin=0.035

  echo "--out_dir ${_output_dir} --depth ${_elec_thin} g4db" \
    "${LDMX_BASE}/dblib/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000/" \
    "&>> ${_output_dir}/g4db_electron_extra_thin.log"
  
  echo "--out_dir ${_output_dir} --depth ${_elec_thin}" \
    "dmg4 -m tungsten --particle electron --primary_energy 4. --ap_mass 100 1" \
    "&>> ${_output_dir}/dmg4_electron_extra_thin.log"
}

__main__ $@ | ./sim/ldmx_parallel.sh
