#!/bin/bash
# script to run a grid of simulations

RUNDIR=${PWD}

r_inner=0.1
r_outer=9.5
species=NeII # REMEMBER to change it also in line_profile.f90

for incl in 0.0; do
	for ng_norm in 0.1 1.0 10.0; do
		echo "hydro: for "$species" with normalization "$ng_norm" and for i="$incl

		# Update the code with the new variables
		sed -i -e "s/ng_norm=.*/ng_norm=$ng_norm /g" line_profile_hydro.f90
		sed -i -e "s/incl_deg=.*/incl_deg=$incl /g" line_profile_hydro.f90
		sed -i -e "s/str_i=.*/str_i='$incl' /g" line_profile_hydro.f90

		# Update the name in the submission script
		sed -i -e "s/#PBS -N .*/#PBS -N photoes_$species\_b${array_b[i]}\_cs$cs\_i$incl/g" submit-job-dial
		sed -i -e "s/#PBS -N .*/#PBS -N photoes_$species\_b${array_b[i]}\_cs$cs\_i$incl/g" submit-job-alice

		# Actually compile the code
		ifort -g -check all -fpe0 -fpp -r8 -warn -traceback -pedantic -traceback -debug extended -qopenmp -o photoes line_profile_hydro.f90

		if [ ! -d $RUNDIR/../data_hydro/$species/ng$ng_norm/incl_$incl ]; then
			if [ ! -d $RUNDIR/../data_hydro/$species/ng$ng_norm ]; then
				if [ ! -d $RUNDIR/../data_hydro/$species ]; then
					mkdir $RUNDIR/../data_hydro/$species
				fi
				mkdir $RUNDIR/../data_hydro/$species/ng$ng_norm
			fi
			mkdir $RUNDIR/../data_hydro/$species/ng$ng_norm/incl_$incl
		fi

		# Copy files to this new directory
		cp photoes submit-job* $RUNDIR/../data_hydro/$species/ng$ng_norm/incl_$incl

		# Remove unuseful files
		rm photoes

		# Submit the job on dial or alice
		cd $RUNDIR/../data_hydro/$species/ng$ng_norm/incl_$incl
		qsub submit-job-alice

		cd $RUNDIR
	done
done
echo "Simulations submitted succesfully!"
