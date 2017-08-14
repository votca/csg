#! /bin/bash
#
# Copyright 2009-2017 The VOTCA Development Team (http://www.votca.org)
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
This script runs a gromacs simulation or pre-simulation with '-multidir' option (see gromacs and VOTCA manuals)

Usage: ${0##*/} [--pre]

Used external packages: gromacs

EOF
   exit 0
fi

tasks=$(get_number_tasks)

mdrun_opts="$(csg_get_property --allow-empty cg.inverse.gromacs.mdrun.opts)"
mdirs=$(echo ${mdrun_opts} | grep -o "sim" | wc -l)

#checks for impossible (due to same checks in run_gromacs.sh)

[[ -n ${mdirs} ]] || die "${0##*/}: mdrun '-multidir' option does not contain directory names including pattern 'sim' (check cg.inverse.gromacs.mdrun.opts)"
[[ ${mdirs} -gt 1 ]] || die "${0##*/}: mdrun '-multidir' option is set but the number of directories including pattern 'sim' <= 1 (check cg.inverse.gromacs.mdrun.opts)"
[[ ${mdirs} -le ${tasks} ]] || die "${0##*/}: mdrun '-multidir' option provided presumes the number of separate simulations greater than the number of tasks/threads (check cg.inverse.gromacs.mdrun.opts)"

echo -e "\nCreating/checking simulation data for ${mdirs} separate simulations (run in parallel)"

confout="$(csg_get_property cg.inverse.gromacs.conf_out)"
conf="$(csg_get_property cg.inverse.gromacs.conf)"
mdp="$(csg_get_property cg.inverse.gromacs.mdp)"
index="$(csg_get_property cg.inverse.gromacs.index)"
topol_in="$(csg_get_property cg.inverse.gromacs.topol_in)"
traj="$(csg_get_property cg.inverse.gromacs.traj)"
checkpoint="$(csg_get_property cg.inverse.gromacs.mdrun.checkpoint)"
tpr="$(csg_get_property cg.inverse.gromacs.topol)"

grompp_opts="$(csg_get_property --allow-empty cg.inverse.gromacs.grompp.opts)"
grompp="$(csg_get_property cg.inverse.gromacs.grompp.bin)"
[[ -n "$(type -p ${grompp})" ]] && echo "using grompp binary '${grompp}'" || die "${0##*/}: grompp binary '${grompp}' not found (override by cg.inverse.gromacs.grompp.bin)"

if is_done "Simtasks"; then
    echo -e "\nAll Simulation tasks done - skipping\n"
else
    for ((it=0;it<mdirs;it++)); do
	this_dir="sim_${it}"

	if [[ -d $this_dir ]]; then
	    cd $this_dir || die "${0##*/}: cd $this_dir failed"
	    if is_done "Simulation"; then 
		echo -e "\nSimulation task ${it} done - skipping\n"
		cd ../
		continue; 
	    fi
	    cd ../
	else
	    echo -e "\nSimulation task $it created at $(date)\n"
	    mkdir -p $this_dir || die "${0##*/}: mkdir -p $this_dir failed (unfinished simulation exists?)"
	fi

	cd $this_dir || die "${0##*/}: cd $this_dir failed"

	filelist="$(csg_get_property --allow-empty cg.inverse.filelist)"

	if is_done "Filecopy"; then
        #check for needed files 
	    echo "Filecopy already done"
	    for f in $filelist; do

		[[ -f $f ]] || cp_from_main_dir "$f"

		echo Comparing "$(get_main_dir)/$f" "$f"
		[[ -z $(type -p cmp) ]] && echo "${0##*/}: program 'cmp' not found, comparision skipped" && continue

		cmp "$(get_main_dir)/$f" "$f" && echo "Unchanged" || \
		msg --color blue --to-stderr "${0##*/}: WARNING: file '$f' in the main dir was changed since the last execution, this will have no effect on the current iteration, to take effect remove the current iteration ('${this_dir##*/}')"

	    done
	else
        #get needed files 
	    echo "Copying the needed files"
	    [[ -n ${filelist} ]] && cp_from_main_dir "$filelist"

	    mark_done "Filecopy"
	fi
	
	if [[ -f ${conf}.${it} ]]; then
    	   critical cp -f "${conf}.${it}" "$conf"
	fi
	[[ -f ${conf} ]] || die "${0##*/}: gromacs initial configuration file '${conf}' not found (make sure it is in cg.inverse.filelist)"

	if [[ -f ${mdp}.${it} ]]; then
    	   critical cp -f "${mdp}.${it}" "$mdp"
	fi
	[[ -f $mdp ]] || die "${0##*/}: gromacs mdp file '$mdp' not found (make sure it is in cg.inverse.filelist)"

	if [[ -f ${index}.${it} ]]; then
    	   critical cp -f "${index}.${it}" "$index"
	fi
	[[ -f $index ]] || die "${0##*/}: grompp index file '$index' not found (make sure it is in cg.inverse.filelist)"

	if [[ -f ${topol_in}.${it} ]]; then
    	   critical cp -f "${topol_in}.${it}" "$topol_in"
	fi
	[[ -f $topol_in ]] || die "${0##*/}: grompp text topol file '$topol_in' not found (make sure it is in cg.inverse.filelist)"

	if [[ $1 != "--pre" ]]; then
        #in a presimulation usually do care about traj and temperature
	   check_temp || die "${0##*/}: check of tempertures failed"
	   if [[ ${traj} == *.xtc ]]; then
	      [[ $(get_simulation_setting nstxtcout 0) -eq 0 ]] && die "${0##*/}: trajectory type (cg.inverse.gromacs.traj) is '${traj##*.}', but nstxtcout is 0 in $mdp. Please check the setting again and remove the current step."
	   elif [[ ${traj} == *.trr ]]; then
		[[ $(get_simulation_setting nstxout 0) -eq 0 ]] && die "${0##*/}: trajectory type (cg.inverse.gromacs.traj) is '${traj##*.}', but nstxout is 0 in $mdp. Please check the setting again and remove the current step."
	   else
	       die "${0##*/}: error trajectory type '${traj##*.}' (ending from '${traj}') is not supported"
	   fi
	fi

        [[ -f ${checkpoint} ]] && echo "using checkpoint file '${checkpoint}'"

#support for older mdp files (before ver.5.0): 
#'cutoff-scheme = Verlet' is the default for Gromacs 5.0, but does not work with tabulated interactions
#XXX is returned if cutoff-scheme is not found in mdp file
#	if [[ $(critical ${grompp} -h 2>&1) = *"VERSION 5.0"* && $(get_simulation_setting cutoff-scheme XXX) = XXX ]]; then

	if [[ $(critical ${grompp} -h 2>&1) = *"VERSION 5"* && $(get_simulation_setting cutoff-scheme XXX) = XXX ]]; then
	   echo "cutoff-scheme = Group" >> $mdp
	   msg --color blue --to-stderr "Automatically added 'cutoff-scheme = Group' to $mdp, tabulated interactions only work with Group cutoff-scheme!"
	fi

	critical cp ../tab*xvg ./

#see if we can run grompp again as checksum of tpr does not appear in the checkpoint
	critical ${grompp} -n "${index}" -f "${mdp}" -p "$topol_in" -o "$tpr" -c "${conf}" ${grompp_opts} 2>&1 
	[[ -f $tpr ]] || die "${0##*/}: gromacs tpr file '$tpr' not found after runing grompp"

	cd ../
    done

    mdrun="$(csg_get_property cg.inverse.gromacs.mdrun.command)"
    #no check for mdrun, because mdrun_mpi could maybe exist only on compute nodes

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

    if [[ -e ${traj} ]]; then

	echo "Overall (combined) trajectory found - skipping simulation"

    else

	critical $mdrun -s "${tpr}" -c "${confout}" -o "${traj%.*}".trr -x "${traj%.*}".xtc ${mdrun_opts} 2>&1 

        simok="YES"
        for ((it=0;it<mdirs;it++)); do
	   this_dir="sim_${it}"

	   critical cd "${this_dir}"

	   if [[ ${simok} == "YES" && -f ${confout} ]]; then
              [[ -z "$(sed -n '/[nN][aA][nN]/p' ${confout})" ]] || { simok="NO"; msg --color blue --to-stderr "${0##*/}: There is a nan in '${this_dir}/${confout}', this seems to be wrong!"; }
           else
              simok="NO"
              msg --color blue --to-stderr "${0##*/}: Gromacs output file ${this_dir}/${confout} not found!"
           fi

	   if [[ ${simok} == "YES" && -f ${traj} ]]; then 
              mark_done "Simulation" 
           else
              simok="NO" 
              msg --color blue --to-stderr "${0##*/}: Trajectory file '${this_dir}/${traj}' not found!"
           fi

	   critical cd ../

        done

	if [[ ${simok} == "YES" ]]; then 
           mark_done "Simtasks"
        else 
           die "${0##*/}: Multi-task simulation failed (try to restart)."
        fi
    fi
fi

is_done "Simtasks" || die "${0##*/}: Multi-task simulation appears unfinished (try to restart)."

trjcat="$(csg_get_property cg.inverse.gromacs.trjcat)"
[[ -n "$(type -p ${trjcat})" ]] && echo "using trjcat binary '${trjcat}'" || die "${0##*/}: trjcat binary '${trjcat}' not found (override by cg.inverse.gromacs.trjcat)"

if is_done "Trajectory"; then

    echo -e "\nTrajectory assembly done for all - skipping\n"

elif  [[ -f "${traj}" ]]; then

    echo -e "\nFound combined trajectory - checking for completeness\n"

    for ((it=0;it<mdirs;it++)); do
	this_dir="sim_${it}"

	critical cd "${this_dir}"

	if is_done "Trajectory"; then

	  echo "Trajectory in ${this_dir} has been appended previously - skipping\n"

	elif is_done "Simulation"; then

          if [[ -f "${traj}" ]]; then

	     [[ -z "$(sed -n '/[nN][aA][nN]/p' ${confout})" ]] || die "${0##*/}: There is a nan in '${this_dir}'/'${confout}', this seems to be wrong."

	     critical echo "0" | ${trjcat} -f ../"${traj}" "${traj}" -o ../"${traj}" -n "$index" -cat

	     mark_done "Trajectory"

             cleanlist="$(csg_get_property --allow-empty cg.inverse.cleanlist)"
	     if [[ -n ${cleanlist} ]]; then
	        msg "Cleaning up files in ${this_dir}: $cleanlist"
	        rm -f ${cleanlist}
	     fi

	  else
	     die "${0##*/}: No trajectory '${traj}' found in '${this_dir}' (simulation failed?)"
	  fi

	else
	  die "${0##*/}: Simulation not finished in '${this_dir}' (simulation failed?)"
	fi

	cleanlist="$(csg_get_property --allow-empty cg.inverse.cleanlist)"
	if [[ -n ${cleanlist} ]]; then
	   msg "Cleaning up files in ${this_dir}: $cleanlist"
	   rm -f ${cleanlist}
	fi

	critical cd ../

    done

    mark_done "Trajectory"
    mark_done "Simulation"

else

    echo -e "\nCombined trajectory not found - checking for simulation results\n"

    if [[ -f "sim_0/${traj}" ]]; then

	critical cd sim_0

        if [[ -f ${confout} ]]; then
	   [[ -z "$(sed -n '/[nN][aA][nN]/p' ${confout})" ]] || msg "${0##*/}: There is a nan in sim_0/'${confout}', this seems to be wrong."
        else
	   msg "${0##*/}: Found trajectory file 'sim_0/${traj}' but not configuration file 'sim_0/${confout}'."
        fi

        is_done Simulation || mark_done "Simulation"

	critical cp -f "$confout" ../
	critical cp -f "${traj}" ../

	mark_done "Trajectory"

	critical cd ../

	for ((it=1;it<mdirs;it++)); do
	   this_dir="sim_${it}"

           if [[ -f "${this_dir}/${traj}" ]]; then

	      critical cd "${this_dir}"

              if [[ -f ${confout} ]]; then
	         [[ -z "$(sed -n '/[nN][aA][nN]/p' ${confout})" ]] || msg "${0##*/}: There is a nan in '${this_dir}/${confout}', this seems to be wrong."
              else
	         msg "${0##*/}: Found trajectory file '${this_dir}/${traj}' but not configuration file '${this_dir}/${confout}'."
              fi

              critical echo "0" | ${trjcat} -f ../"${traj}" "${traj}" -o ../"${traj}" -n "$index" -cat

	      mark_done "Trajectory"

	      cleanlist="$(csg_get_property --allow-empty cg.inverse.cleanlist)"
	      if [[ -n ${cleanlist} ]]; then
	         msg "Cleaning up files in ${this_dir}: $cleanlist"
	         rm -f ${cleanlist}
	      fi

	      critical cd ../

	   else
	      die "${0##*/}: No trajectory '${traj}' found in '${this_dir}' (simulation failed?)"
	   fi

	done

        #copying critial output files to current step directory to be able to analyse
	critical cp -fv sim_0/*edr sim_0/*cpt sim_0/*.tpr ./

	critical cd sim_0

	cleanlist="$(csg_get_property --allow-empty cg.inverse.cleanlist)"
	if [[ -n ${cleanlist} ]]; then
	   msg "Cleaning up files in sim_0: $cleanlist"
	   rm -f ${cleanlist}
	fi

	critical cd ../

	mark_done "Trajectory"
	mark_done "Simulation"

    else
	die "${0##*/}: No trajectory '${traj}' found in sim_0 (simulation failed?)"
    fi

fi
