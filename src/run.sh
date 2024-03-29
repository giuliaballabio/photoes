#!/bin/bash
# script to run a grid of simulations

RUNDIR=${PWD}

array_b=( 0.75 1.00 1.50 ) #2.00 )
array_ub=( 0.85 0.77 0.56 ) #0.29 )

r_inner=0.1
r_outer=9.5
species=NeII # REMEMBER to change it also in line_profile.f90
mdot=mdot10e-9 # REMEMBER to change it also in line_profile.f90

for ((i=0;i<${#array_b[@]};++i)); do
	for cs in 3.0d5 5.0d5 10.0d5; do
		for incl in 0.0 5.0 10.0 15.0 20.0 25.0 30.0 35.0 40.0 45.0 50.0 55.0 60.0 65.0 70.0 75.0 80.0 85.0 90.0; do

			echo $species" with b="${array_b[i]}", mdot="$mdot " and cs="$cs[j] " for i="$incl

			# Update the code with the new variables
			sed -i -e "s/b_input=.*/b_input=${array_b[i]} /g" selfsimilar_solutions.f90
			sed -i -e "s/ub=.*/ub=${array_ub[i]} /g" selfsimilar_solutions.f90
			sed -i -e "s/cs=.*/cs=$cs /g" selfsimilar_solutions.f90
			# sed -i -e "s/species_flag=.*/species_flag=$species_flag /g" line_profile.f90
			sed -i -e "s/b_input=.*/b_input=${array_b[i]} /g" line_profile.f90
			sed -i -e "s/ub=.*/ub=${array_ub[i]} /g" line_profile.f90
			sed -i -e "s/cs=.*/cs=$cs /g" line_profile.f90
			sed -i -e "s/r_inner=.*/r_inner=$r_inner /g" line_profile.f90
			sed -i -e "s/r_outer=.*/r_outer=$r_outer /g" line_profile.f90
			sed -i -e "s/incl_deg=.*/incl_deg=$incl /g" line_profile.f90
			sed -i -e "s/str_i=.*/str_i='$incl' /g" line_profile.f90

			# Update the name in the submission script
			sed -i -e "s/#PBS -N .*/#PBS -N photoes_$species\_b${array_b[i]}\_cs$cs\_i$incl/g" submit-job-dial
			sed -i -e "s/#PBS -N .*/#PBS -N photoes_$species\_b${array_b[i]}\_cs$cs\_i$incl/g" submit-job-alice

			# Actually compile the code
			ifort -g -check all -fpe0 -fpp -r8 -warn -traceback -pedantic -traceback -debug extended -qopenmp -o selfsimilar_solutions selfsimilar_solutions.f90
			ifort -g -check all -fpe0 -fpp -r8 -warn -traceback -pedantic -traceback -debug extended -qopenmp -o photoes line_profile.f90

			if [ ! -d $RUNDIR/../cs$cs/$species/$mdot/data_b${array_b[i]}\_r$r_inner\_r$r_outer/incl_$incl ]; then
				if [ ! -d $RUNDIR/../cs$cs/$species/$mdot/data_b${array_b[i]}\_r$r_inner\_r$r_outer ]; then
					if [ ! -d $RUNDIR/../cs$cs/$species/$mdot ]; then
						if [ ! -d $RUNDIR/../cs$cs/$species ]; then
							if [ ! -d $RUNDIR/../cs$cs ]; then
								mkdir $RUNDIR/../cs$cs
							fi
							mkdir $RUNDIR/../cs$cs/$species
						fi
						mkdir $RUNDIR/../cs$cs/$species/$mdot
					fi
					mkdir $RUNDIR/../cs$cs/$species/$mdot/data_b${array_b[i]}\_r$r_inner\_r$r_outer
				fi
				mkdir $RUNDIR/../cs$cs/$species/$mdot/data_b${array_b[i]}\_r$r_inner\_r$r_outer/incl_$incl
			fi

			# Copy files to this new directory
			cp selfsimilar_solutions photoes submit-job* $RUNDIR/../cs$cs/$species/$mdot/data_b${array_b[i]}\_r$r_inner\_r$r_outer/incl_$incl

			# Remove unuseful files
			rm *genmod* photoes selfsimilar_solutions

			# Submit the job on dial or alice
			cd $RUNDIR/../cs$cs/$species/$mdot/data_b${array_b[i]}\_r$r_inner\_r$r_outer/incl_$incl
			qsub submit-job-alice

			cd $RUNDIR

		done
	done
done
echo "Simulations submitted succesfully!"
