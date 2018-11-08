#!/bin/bash
# script to run a grid of simulations

RUNDIR=${PWD}
file=selfsimilar_solutions.f90

for set in 1; do
	for r_input in 0.03; do

		# b_input=$(grep -i "b_input=" "$file")
		# echo $b_input

		echo "up to" $set "and" $r_input

		# Update the code with the new variables
		sed -i -e "s/case=.*/case ($set /g" selfsimilar_solutions.f90
		#sed -i -e "s/b_input=.*/b_input=$b_input /g" selfsimilar_solutions.f90
		sed -i -e "s/r_input=.*/r_inner=$r_input /g" line_profile.f90

		# Update the name in the submission script
		sed -i -e "s/#PBS -N .*/#PBS -N photoes_b$set\_r$r_input/g" submit-job

		# Actually compile the code
		ifort -g -check all -fpe0 -warn -traceback -debug extended -qopenmp -o selfsimilar_solutions selfsimilar_solutions.f90
		ifort -g -check all -fpe0 -warn -traceback -debug extended -qopenmp -o photoes line_profile.f90

		# Make a new directory if it doesn't existnon
		if [ ! -d $RUNDIR/../data_b$set\_r$r_input ]; then
				mkdir $RUNDIR/../data_b$set\_r$r_input
		fi

		# Copy files to this new directory
		cp selfsimilar_solutions photoes submit-job $RUNDIR/../data_b$set\_r$r_input

		# Remove unuseful files
		rm *genmod* photoes selfsimilar_solutions

		# Submit the job if on dial
		#qsub $RUNDIR/../submit-job

	done
done
