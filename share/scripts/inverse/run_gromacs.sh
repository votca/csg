#! /bin/bash
#
# Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

if [[ $1 = "--help" ]]; then
cat <<EOF
${0##*/}, version %version%
This script runs a gromacs simulation or pre-simulation

Usage: ${0##*/} [--pre]

Used external packages: gromacs

EOF
   exit 0
fi

tpr="$(csg_get_property cg.inverse.gromacs.topol)"

mdp="$(csg_get_property cg.inverse.gromacs.mdp)"
[[ -f $mdp ]] || die "${0##*/}: gromacs mdp file '$mdp' not found (make sure it is in cg.inverse.filelist)"

conf="$(csg_get_property cg.inverse.gromacs.conf)"
[[ -f $conf ]] || die "${0##*/}: gromacs initial configuration file '$conf' not found (make sure it is in cg.inverse.filelist)"

confout="$(csg_get_property cg.inverse.gromacs.conf_out)"
mdrun_opts="$(csg_get_property --allow-empty cg.inverse.gromacs.mdrun.opts)"

index="$(csg_get_property cg.inverse.gromacs.index)"
[[ -f $index ]] || die "${0##*/}: grompp index file '$index' not found (make sure it is in cg.inverse.filelist)"

topol_in="$(csg_get_property cg.inverse.gromacs.topol_in)"
[[ -f $topol_in ]] || die "${0##*/}: grompp text topol file '$topol_in' not found"

grompp_opts="$(csg_get_property --allow-empty cg.inverse.gromacs.grompp.opts)"

grompp=( $(csg_get_property cg.inverse.gromacs.grompp.bin) )
[[ -n "$(type -p ${grompp[0]})" ]] || die "${0##*/}: grompp binary '${grompp[0]}' not found"

traj=$(csg_get_property cg.inverse.gromacs.traj)
if [[ $1 != "--pre" ]]; then
  #in a presimulation usually do care about traj and temperature
  check_temp || die "${0##*/}: check of tempertures failed"
  if  [[ $traj == *.xtc ]]; then
    [[ $(get_simulation_setting nstxtcout 0) -eq 0 ]] && die "${0##*/}: trajectory type (cg.inverse.gromacs.traj) is '${traj##*.}', but nstxtcout is 0 in $mdp. Please check the setting again and remove the current step."
  elif [[ $traj == *.trr ]]; then
    [[ $(get_simulation_setting nstxout 0) -eq 0 ]] && die "${0##*/}: trajectory type (cg.inverse.gromacs.traj) is '${traj##*.}', but nstxout is 0 in $mdp. Please check the setting again and remove the current step."
  else
    die "${0##*/}: error trajectory type '${traj##*.}' (ending from '$traj') is not supported"
  fi
fi

checkpoint="$(csg_get_property cg.inverse.gromacs.mdrun.checkpoint)"

#we want to do a presimulation and we are not currently performing it
if [[ $(csg_get_property cg.inverse.gromacs.pre_simulation) = "yes" && $1 != "--pre" ]]; then
  for i in tpr mdp conf confout index topol_in checkpoint; do
    f=${!i}
    [[ $f = */* ]] && die "${0##*/}: presimulation feature only work with local file (without /) in $i variable. Just try to copy $f to the maindir and add it to cg.inverse.filelist."
  done
  if ! is_done "Presimulation"; then #prepare/do presimulation
    [[ -d pre_simulation ]] || critical mkdir pre_simulation
    critical cp "${conf}" ./pre_simulation
    # a bit hacky but we have no better solution yet
    critical cp table*.xvg ./pre_simulation
    cp=0
    for i in mdp topol_in index; do
      f="$(csg_get_property "cg.inverse.gromacs.pre_simulation.$i" "${!i}")" #filter me away
      [[ $f != ${!i} ]] && ((cp++))
      [[ -f ${f} ]] || die "${0##*/}: file '$f' not found (make sure it is in cg.inverse.filelist)"
      critical cp "${f}" "./pre_simulation/${!i}"
    done
    [[ $cp -eq 0 ]] && die "${0##*/}: mdp, topol_in and index of the presimulation are the same as for the main simulation, that does not make sense - at least one has to be different!"
    cd pre_simulation || die "${0##*/}: cd pre_simulation failed"
    msg "Doing pre simulation"
    do_external presimulation gromacs #easy to overwrite
    simulation_finish --no-traj || exit 0 #time out or crash - checkpoint_exist in inverse.sh will also look for a checkpoint in pre_simulation.
    critical cp "${confout}" ../"${conf}"
    cd .. || die "${0##*/}: cd .. failed"
    mark_done "Presimulation" #after 'cd..' as restart file is here
  fi
  msg "Doing main simulation"
fi

#support for older mdp file, cutoff-scheme = Verlet is default for Gromacs 5.0, but does not work with tabulated interactions
#XXX is returned if cutoff-scheme is not in mdp file
if [[ $(critical ${grompp[@]} -h 2>&1) = *"VERSION 5."[01]* && $(get_simulation_setting cutoff-scheme XXX) = XXX ]]; then
  echo "cutoff-scheme = Group" >> $mdp
  msg --color blue --to-stderr "Automatically added 'cutoff-scheme = Group' to $mdp, tabulated interactions only work with Group cutoff-scheme!"
fi
# Change of mdp file of gromacs in order to optimize Electrostatic interaction in CGMD simulation
if [[ $(csg_get_property cg.inverse.permittivity) = true ]]; then
  permittivity=$(awk '{print $1}' permittivity.current | tail -n 1)
  echo "deleting current value of epsilon_re from $mdp and replcing it with true value"
  sed -i '/epsilon_r/d' $mdp
  echo "epsilon_r       =  $permittivity" >> $mdp;
  msg --color blue --to-stderr "Automatically added permittivity (epsilon_r = $permittivity) to $mdp  !"
fi





if [[ ${CSG_MDRUN_STEPS} ]]; then
  #mdrun <4.6 does not have -nsteps options, remove this block whenever we drop support for <4.6
  if [[ $(critical ${grompp[@]} -h 2>&1) = *"VERSION 4."[05]* ]]; then
    nsteps=$(get_simulation_setting nsteps)
    critical sed -i "/^nsteps/s/=.*/=${CSG_MDRUN_STEPS}/" $mdp
    msg --color blue --to-stderr "Replace nsteps (=$nsteps) in '$mdp' to be ${CSG_MDRUN_STEPS}"
  else
    msg --color blue --to-stderr "Appending -nsteps ${CSG_MDRUN_STEPS} to mdrun options"
    mdrun_opts+=" -nsteps $CSG_MDRUN_STEPS"
  fi
fi

#see can run grompp again as checksum of tpr does not appear in the checkpoint
critical ${grompp[@]} -n "${index}" -f "${mdp}" -p "$topol_in" -o "$tpr" -c "${conf}" ${grompp_opts} 2>&1 | gromacs_log "${grompp[@]} -n "${index}" -f "${mdp}" -p "$topol_in" -o "$tpr" -c "${conf}" ${grompp_opts}"
[[ -f $tpr ]] || die "${0##*/}: gromacs tpr file '$tpr' not found after runing grompp"

mdrun="$(csg_get_property cg.inverse.gromacs.mdrun.command)"
#no check for mdrun, because mdrun_mpi could maybe exist only computenodes

if [[ -n $CSGENDING ]]; then
  #seconds left for the run
  wall_h=$(( $CSGENDING - $(get_time) ))
  #convert to hours
  wall_h=$(csg_calc $wall_h / 3600 )
  echo "${0##*/}: Setting $mdrun maxh option to $wall_h (hours)"
  mdrun_opts="-cpi $checkpoint -maxh $wall_h ${mdrun_opts}"
else
  echo "${0##*/}: No walltime defined, so no time limitation given to $mdrun"
fi

critical $mdrun -s "${tpr}" -c "${confout}" -o "${traj%.*}".trr -x "${traj%.*}".xtc ${mdrun_opts} ${CSG_RUNTEST:+-v} 2>&1 | gromacs_log "$mdrun -s "${tpr}" -c "${confout}" -o "${traj%.*}".trr -x "${traj%.*}".xtc ${mdrun_opts}"

[[ -z "$(sed -n '/[nN][aA][nN]/p' ${confout})" ]] || die "${0##*/}: There is a nan in '${confout}', this seems to be wrong."

