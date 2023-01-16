#!/bin/bash
# run.sh 
#   run four simulations within container

# main generator of commands
# Expected Env Variables
#   _output_dir : directory in which to write files
#   _only_xsec  : true if only xsec, false if everything
__main__() {
  if ! mkdir -p ${_output_dir}; then
    echo "ERROR: Could not create output directory ${_output_dir}"
    return $?
  fi

  for method in "fullww" "iww" "hiww"; do
    echo "g4db-xsec-calc -o ${_output_dir}/g4db_el_xsec_${method}.csv" \
      "--energy 0.2 100 --ap-mass 0.1 --method ${method}" \
      "&>> ${_output_dir}/g4db_el_xsec_${method}.log"
    echo "g4db-xsec-calc -o ${_output_dir}/g4db_mu_xsec_${method}.csv" \
      "--muons --energy 2 1000 1 --ap-mass 1 --target 29 63.55 --method ${method}" \
      "&>> ${_output_dir}/g4db_xsec_${method}.log"
  done

  if ${_only_xsec}; then
    return
  fi

  local _muon_thin=100
  local _elec_thin=0.35
  local _elec_thick=18
  local _muon_thick=2000
  
  echo "fire ${LDMX_BASE}/sim/config.py --out_dir ${_output_dir} --depth ${_elec_thick} g4db" \
    "${LDMX_BASE}/dblib/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000/" \
    "&>> ${_output_dir}/g4db_electron_thick.log"
  
  echo "fire ${LDMX_BASE}/sim/config.py --out_dir ${_output_dir} --depth ${_elec_thick}" \
    "dmg4 -m tungsten --particle electron --primary_energy 4. --ap_mass 100 1" \
    "&>> ${_output_dir}/dmg4_electron_thick.log"

  echo "fire ${LDMX_BASE}/sim/config.py --out_dir ${_output_dir} --depth ${_elec_thin} g4db" \
    "${LDMX_BASE}/dblib/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000/" \
    "&>> ${_output_dir}/g4db_electron_thin.log"
  
  echo "fire ${LDMX_BASE}/sim/config.py --out_dir ${_output_dir} --depth ${_elec_thin}" \
    "dmg4 -m tungsten --particle electron --primary_energy 4. --ap_mass 100 1" \
    "&>> ${_output_dir}/dmg4_electron_thin.log"

  echo "fire ${LDMX_BASE}/sim/config.py --out_dir ${_output_dir} --depth ${_muon_thin} g4db" \
    "${LDMX_BASE}/dblib/muon_copper_MaxE_100.0_MinE_2.0_RelEStep_0.1_UndecayedAP_mA_1.0_run_3000/" \
    "&>> ${_output_dir}/g4db_muon_thin.log"
  
  echo "fire ${LDMX_BASE}/sim/config.py --out_dir ${_output_dir} --depth ${_muon_thin}" \
    "dmg4 -m brass --particle muon --primary_energy 100. --ap_mass 1000 1" \
    "&>> ${_output_dir}/dmg4_muon_thin.log"

  echo "fire ${LDMX_BASE}/sim/config.py --out_dir ${_output_dir} --depth ${_muon_thick} g4db" \
    "${LDMX_BASE}/dblib/muon_copper_MaxE_100.0_MinE_2.0_RelEStep_0.1_UndecayedAP_mA_1.0_run_3000/" \
    "&>> ${_output_dir}/g4db_muon_thick.log"
  
  echo "fire ${LDMX_BASE}/sim/config.py --out_dir ${_output_dir} --depth ${_muon_thick}" \
    "dmg4 -m brass --particle muon --primary_energy 100. --ap_mass 1000 1" \
    "&>> ${_output_dir}/dmg4_muon_thick.log"

  _elec_thin=1.

  echo "fire ${LDMX_BASE}/sim/config.py --out_dir ${_output_dir} --depth ${_elec_thin} g4db" \
    "${LDMX_BASE}/dblib/electron_lead_MaxE_100.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000/" \
    "&>> ${_output_dir}/g4db_electron_100GeV_lead_thin.log"
  
  echo "fire ${LDMX_BASE}/sim/config.py --out_dir ${_output_dir} --depth ${_elec_thin}" \
    "dmg4 -m lead --particle electron --primary_energy 100. --ap_mass 100 1" \
    "&>> ${_output_dir}/dmg4_electron_100GeV_lead_thin.log"

  _elec_thin=0.035

  echo "fire ${LDMX_BASE}/sim/config.py --out_dir ${_output_dir} --depth ${_elec_thin} g4db" \
    "${LDMX_BASE}/dblib/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000/" \
    "&>> ${_output_dir}/g4db_electron_extra_thin.log"
  
  echo "fire ${LDMX_BASE}/sim/config.py --out_dir ${_output_dir} --depth ${_elec_thin}" \
    "dmg4 -m tungsten --particle electron --primary_energy 4. --ap_mass 100 1" \
    "&>> ${_output_dir}/dmg4_electron_extra_thin.log"
}

__usage__() {
  cat <<HELP
 USAGE:
  ./sim/run.sh [-o OUT_DIR] [--only-xsec] [--dry-run] [-- PARALLEL_OPTS]

 OPTIONS:
  -o            : Base output directory for data files (default: 'data/<git describe --tags>')
  --only-xsec   : only do xsec calculations
  --dry-run, -n : print commands to terminal instead of piping to parallel
  PARALLEL_OPTS : all arguments after '--' are given to parallel (e.g. use -- -j 4 to limit
                  the number of cores used).

 We do two simulations at once (G4DB and DMG4) and then we loop over the combinations
 of thicknesses and incident particles. The multiple jobs are written to a job listing
 file in the output directory and then "submitted" to parallel via the ldmx_parallel.sh
 script.

HELP
}

_tag=$(git describe --tags)
_output_dir=$(cd data && pwd -P)/${_tag}
_only_xsec=false
_dry_run=false
while [ $# -gt 0 ]; do
  case $1 in
    -o)
      if [ -z $2 ]; then
        echo "ERROR: '$1' requires an argument."
        exit 1
      fi
      case $1 in
        -o) _output_dir=$2;;
      esac
      shift
      shift
      ;;
    --only-xsec)
      _only_xsec=true
      shift
      ;;
    --dry-run|-n)
      _dry_run=true
      shift
      ;;
    -h|--help)
      __usage__
      exit 0
      ;;
    --)
      # the rest of the options are options for parallel
      shift
      break
      ;;
    *)
      echo "ERROR: Unrecognized option '$1'."
      exit 1
      ;;
  esac
done

if ${_dry_run}; then
  __main__
  echo "with parallel options '$@'"
else 
  __main__ | parallel $@ ./sim/ldmx_parallel.sh
fi

