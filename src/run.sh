#!/bin/bash
# script to run a grid of simulations

RUNDIR=${PWD}

array_b=( 0.75 ) # 1.00 1.50 ) # 2.00 )
array_ub=( 0.85 ) # 0.77 0.56 ) # 0.29 )

for ((i=0;i<${#array_b[@]};++i)); do
  # echo "(${array_b[i]}, ${array_ub[i]})"
	for r_input in 0.03 0.1; do
		for incl in 0.0 90.0; do

			echo "up to b="${array_b[i]} "and R_in="$r_input "for i="$incl

			# Update the code with the new variables
			sed -i -e "s/b_input=.*/b_input=${array_b[i]} /g" selfsimilar_solutions.f90
			sed -i -e "s/ub=.*/ub=${array_ub[i]} /g" selfsimilar_solutions.f90
			sed -i -e "s/r_inner=.*/r_inner=$r_input /g" line_profile.f90
			sed -i -e "s/incl_deg=.*/incl_deg=$incl /g" line_profile.f90
			sed -i -e "s/str_i=.*/str_i='$incl' /g" line_profile.f90

			# Update the name in the submission script
			sed -i -e "s/#PBS -N .*/#PBS -N photoes_b${array_b[i]}\_r$r_input\_i$incl/g" submit-job-dial
			sed -i -e "s/#PBS -N .*/#PBS -N photoes_b${array_b[i]}\_r$r_input\_i$incl/g" submit-job-alice

			# Actually compile the code
			ifort -g -check all -fpe0 -warn -traceback -debug extended -qopenmp -o selfsimilar_solutions selfsimilar_solutions.f90
			ifort -g -check all -fpe0 -warn -traceback -debug extended -qopenmp -o photoes line_profile.f90

			# Make a new directory if it doesn't exist
			# if [ ! -d $RUNDIR/../data_b${array_b[i]}\_r$r_input ]; then
			# 	mkdir $RUNDIR/../data_b${array_b[i]}\_r$r_input
			# fi
			#
			# cd $RUNDIR/../data_b${array_b[i]}\_r$r_input
			if [ ! -d $RUNDIR/../data_b${array_b[i]}\_r$r_input/incl_$incl ]; then
				mkdir $RUNDIR/../data_b${array_b[i]}\_r$r_input/incl_$incl
			fi

			# Copy files to this new directory
			cp selfsimilar_solutions photoes submit-job* $RUNDIR/../data_b${array_b[i]}\_r$r_input/incl_$incl

			# Submit the job on dial or alice
			#qsub $RUNDIR/../data_b${array_b[i]}\_r$r_input/incl_$incl/submit-job-alice
			#qsub $RUNDIR/../submit-job-dial

			# Remove unuseful files
			rm *genmod* photoes selfsimilar_solutions

		done
	done
done
