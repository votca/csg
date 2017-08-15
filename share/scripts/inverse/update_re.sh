#! /bin/bash
#
# Copyright 2009-2015 The VOTCA Development Team (http://www.votca.org)
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
#if [[ $(csg_get_property cg.inverse.pressure) = true ]]; then
 
  topol="$(csg_get_property --allow-empty cg.inverse.gromacs.g_energy.topol)"
  [[ -z $topol ]] && topol=$(csg_get_property cg.inverse.gromacs.topol)
  [[ -f $topol ]] || die "${0##*/}: Gromacs tpr file '$topol' not found"

  g_energy=( $(csg_get_property cg.inverse.gromacs.g_energy.bin) )
  [[ -n "$(type -p ${g_energy[0]})" ]] || die "${0##*/}: g_energy binary '${g_energy[0]}' not found"


  opts="$(csg_get_property --allow-empty cg.inverse.gromacs.g_energy.opts)"

  begin="$(calc_begin_time)"
  if [[ ${CSG_RUNTEST} ]] && csg_calc "$begin" ">" "0"; then
     msg --color blue --to-stderr "Automatically setting begin time to 0, because CSG_RUNTEST was set"
     begin=0
  fi

  echo "Running ${g_energy[@]}"
  #no critical here so that we can print the error
  output=$(echo Pressure | ${g_energy[@]} -b "${begin}" -s "${topol}" ${opts} 2>&1)
  ret="$?"
  echo "$output" | gromacs_log "${g_energy[@]} -b "${begin}" -s "${topol}" ${opts}"
  [[ $ret -eq 0 ]] || die "${0##*/}: '${g_energy[@]} -b "${begin}" -s "${topol}" ${opts}' failed"
  #the number pattern '-\?[0-9][^[:space:]]*[0-9]' is ugly, but it supports X X.X X.Xe+X Xe-X and so on
  #awk 'print $2' does not work for older version of g_energy as the format varies between
  #^Pressure XXX (bar) and ^Pressure (bar) XXX
  p_now=$(echo "$output" | sed -n 's/^Pressure[^-0-9]*\(\(-\?[0-9][^[:space:]]*[0-9]\|nan\)\)[[:space:]].*$/\1/p' ) || \
  die "${0##*/}: sed grep of Pressure failed"
  [[ -z $p_now ]] && die "${0##*/}: Could not get pressure from simulation"

  [[ $p_now = nan && $(csg_get_property cg.inverse.gromacs.g_energy.pressure.allow_nan) = no ]] && \
  die "${0##*/}: Pressure was nan, check your simulation (this usually means system has blow up -> use pre simulation)"
  echo "${p_now}" > CGMDPressure.new
  
  # Calculating related to Coulombic Interactions
  echo ";Total Coulombic interaction per frame of CGMD simulation" >  FrameCoul.new
  nstxtcout=$(get_simulation_setting nstxtcout)
  mddt=$(get_simulation_setting dt)
  dnst=$(echo $nstxtcout*$mddt | bc)
  nsteps=$(get_simulation_setting nsteps)
  tmp=$(echo $nsteps*$mddt | bc)
  tmp=$(echo $tmp-$begin | bc)
  nframes=$(echo $tmp/$dnst+1 | bc )
  frame=0
  UCol=0
  U2Col=0
  if [[ $(csg_get_property cg.inverse.permittivity) = true ]]; then
      echo $nframes > Coul_pot.new
      while [ "$frame" -lt "$nframes" ]
      do
       Dt=$(echo $frame*$dnst | bc)
       frame_t=$(echo $begin + $Dt | bc)
       echo Coulomb-'('SR')' | ${g_energy[@]} -b "${frame_t}" -e "${frame_t}" -s "${topol}" ${opts} 2>&1 > tmpcol.new
       echo Coul.-recip. | ${g_energy[@]} -b "${frame_t}" -e "${frame_t}" -s "${topol}" ${opts} 2>&1 >> tmpcol.new
       coulSR=$(sed -n 's/^Coulomb (SR)[^-0-9]*\(\(-\?[0-9][^[:space:]]*[0-9]\|nan\)\)[[:space:]].*$/\1/p' tmpcol.new)
       coulRE=$(sed -n 's/^Coul. recip.[^-0-9]*\(\(-\?[0-9][^[:space:]]*[0-9]\|nan\)\)[[:space:]].*$/\1/p' tmpcol.new)
       CoulTF=$(echo $coulSR + $coulRE | bc)
       frame=$(echo $frame + 1 | bc)
       echo $CoulTF >> Coul_pot.new
       UCol=$(echo $UCol + $CoulTF | bc)
       U2Col=$(echo $CoulTF*$CoulTF+$U2Col | bc)
       rm \#*
       done
       UCol=$(echo $UCol/$nframes | bc -l)
       U2Col=$(echo $U2Col/$nframes | bc -l)
       echo $UCol > UCol.new
       echo $U2Col > U2Col.new
   fi
 
#fi
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

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version %version%
This script implements update step of relative entropy method by csg_reupdate program

Usage: ${0##*/}
EOF
   exit 0
fi

sim_prog="$(csg_get_property cg.inverse.program)"

topol=$(csg_get_property --allow-empty cg.inverse.$sim_prog.re.topol)
[[ -z $topol ]] && topol=$(csg_get_property cg.inverse.$sim_prog.topol)
[[ -f $topol ]] || die "${0##*/}: topol file '$topol' not found, possibly you have to add it to cg.inverse.filelist"

traj=$(csg_get_property cg.inverse.$sim_prog.traj)
[[ -f $traj ]] || die "${0##*/}: traj file '$traj' not found"

equi_time="$(csg_get_property cg.inverse.$sim_prog.equi_time)"
if [[ ${CSG_RUNTEST} ]] && csg_calc "$equi_time" ">" "0"; then
  msg --color blue --to-stderr "Automatically setting equi_time to 0, because CSG_RUNTEST was set"
  equi_time=0
fi

first_frame="$(csg_get_property cg.inverse.$sim_prog.first_frame)"
csg_reupdate_opts="$(csg_get_property --allow-empty cg.inverse.re.csg_reupdate.opts)"
if [[ ${CSG_RUNTEST} ]] ; then
  msg --color blue --to-stderr "Automatically adding '--hessian-check no', because CSG_RUNTEST was set"
  csg_reupdate_opts+=" --hessian-check no"
fi

tasks=$(get_number_tasks)
if is_done "re_update"; then
    echo "RE update is already done"
else
  #copy+resample all target dist in $this_dir
    for_all "non-bonded bonded" do_external resample target '$(csg_get_interaction_property inverse.target)' '$(csg_get_interaction_property name).dist.tgt'

    critical csg_reupdate --nt $tasks --top ${topol} --trj $traj --options $CSGXMLFILE --begin $equi_time --first-frame $first_frame ${csg_reupdate_opts}
    mark_done "re_update"
fi
