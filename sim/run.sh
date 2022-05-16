#!/bin/bash
# run.sh 
#   run four simulations within container

__usage__() {
  cat <<HELP
 USAGE:
  ldmx run.sh [-m MUON_DEPTH] [-e ELEC_DEPTH] [-o OUT_DIR]

 OPTIONS:
  -m    : Depth in mm to make brass target for muons (default: 2000)
  -e    : Depth in mm to make tungsten target for electrons (default: 18)
  -o    : Base output directory for data files (default: 'data/<git describe --tags>')
HELP
}

__main__() {
  local _muon_thickness=2000
  local _elec_thickness=18
  local _tag=$(git describe --tags)
  local _output_dir=$(cd ${1:-data} && pwd -P)/${_tag}
  while [ $# -gt 0 ]; do
    case $1 in
      -m|-t|-o)
        if [ -z $2 ]; then
          echo "ERROR: '$1' requires an argument."
          return 1
        fi
        case $1 in
          -m) _muon_thickness=$2;;
          -e) _elec_thickness=$2;;
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
  local _mgs_log=${_output_dir}/mgs_fire.log
  local _dmg4_log=${_output_dir}/dmg4_fire.log
  
  echo "Muons $(date)"

  fire ${LDMX_BASE}/sim/config.py \
    ${LDMX_BASE}/dblib/muon_copper_MaxE_100.0_MinE_2.0_RelEStep_0.1_UndecayedAP_mA_1.0_run_3000/ \
    --depth ${_muon_thickness} \
    --out_dir ${_output_dir} &>> ${_mgs_log} &
  
  fire ${LDMX_BASE}/sim/dmg4.py \
    -m brass --particle mu- --primary_energy 100. --ap_mass 1000 1 \
    --depth ${_muon_thickness} \
    --out_dir ${_output_dir} &>> ${_dmg4_log} &

  wait
  echo "Electrons $(date)"
  
  fire ${LDMX_BASE}/sim/config.py \
    ${LDMX_BASE}/dblib/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000/ \
    --depth ${_elec_thickness} \
    --out_dir ${_output_dir} &>> ${_mgs_log} &
  
  fire ${LDMX_BASE}/sim/dmg4.py \
    -m tungsten --particle e- --primary_energy 4. --ap_mass 100 1 \
    --depth ${_elec_thickness} \
    --out_dir ${_output_dir} &>> ${_dmg4_log} &

  wait
  echo "done $(date)"
}

__main__ $@
