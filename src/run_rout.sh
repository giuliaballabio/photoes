#!/bin/bash
# script to run a grid of simulations

RUNDIR=${PWD}

array_b=( 0.75 1.00 1.50 2.00 )
array_ub=( 0.85 0.77 0.56 0.29 )

cs=10.0d5
string_cs=10
species=NeII # REMEMBER to change it also in line_profile.f90
mdot=mdot10e-10 # REMEMBER to change it also in line_profile.f90

for ((i=0;i<${#array_b[@]};++i)); do
  # echo "(${array_b[i]}, ${array_ub[i]})"
	for r_inner in 0.1; do
		for r_outer in 30.0; do
			for incl in 0.0 20.0 30.0 45.0 60.0 75.0 90.0; do

				echo "up to b="${array_b[i]}", R_out="$r_outer "and R_in="$r_inner "for i="$incl

				# Update the code with the new variables
				sed -i -e "s/b_input=.*/b_input=${array_b[i]} /g" selfsimilar_solutions.f
				sed -i -e "s/ub=.*/ub=${array_ub[i]} /g" selfsimilar_solutions.f
				# sed -i -e "s/cs=.*/cs=$cs /g" selfsimilar_solutions.f
				# sed -i -e "s/species_flag=.*/species_flag=$species_flag /g" line_profile.f90
				sed -i -e "s/b_input=.*/b_input=${array_b[i]} /g" line_profile.f90
				sed -i -e "s/ub=.*/ub=${array_ub[i]} /g" line_profile.f90
				sed -i -e "s/cs=.*/cs=$cs /g" line_profile.f90
				sed -i -e "s/r_inner=.*/r_inner=$r_inner /g" line_profile.f90
				sed -i -e "s/r_outer=.*/r_outer=$r_outer /g" line_profile.f90
				sed -i -e "s/incl_deg=.*/incl_deg=$incl /g" line_profile.f90
				sed -i -e "s/str_i=.*/str_i='$incl' /g" line_profile.f90

				# Update the name in the submission script
				sed -i -e "s/#PBS -N .*/#PBS -N photoes_b${array_b[i]}\_r$r_inner\_r$r_outer\_i$incl/g" submit-job-dial
				sed -i -e "s/#PBS -N .*/#PBS -N photoes_b${array_b[i]}\_r$r_inner\_r$r_outer\_i$incl/g" submit-job-alice

				# Actually compile the code
				gfortran -Wunused-variable -fno-range-check -Wextra -ffpe-trap=invalid,zero,overflow -finit-real=snan -pedantic -fbounds-check -g -o selfsimilar_solutions selfsimilar_solutions.f
				#ifort -g -check all -fpe0 -fpp -r8 -warn -traceback -pedantic -traceback -debug extended -qopenmp -o selfsimilar_solutions selfsimilar_solutions.f
				ifort -g -check all -fpe0 -fpp -r8 -warn -traceback -pedantic -traceback -debug extended -qopenmp -o photoes line_profile.f90
				# gfortran -Wunused-variable -fdefault-real-8 -fdefault-double-8 -Wextra -ffpe-trap=invalid,zero,overflow -pedantic -finit-real=snan -fbounds-check -g -fopenmp -o selfsimilar_solutions selfsimilar_solutions.f90
				# gfortran -Wunused-variable -fdefault-real-8 -fdefault-double-8 -Wextra -ffpe-trap=invalid,zero,overflow -pedantic -finit-real=snan -fbounds-check -g -fopenmp -o photoes line_profile.f90

				if [ ! -d $RUNDIR/../cs$string_cs\kms/$species/$mdot/data_b${array_b[i]}\_r$r_inner\_r$r_outer/incl_$incl ]; then
					if [ ! -d $RUNDIR/../cs$string_cs\kms/$species/$mdot/data_b${array_b[i]}\_r$r_inner\_r$r_outer ]; then
						if [ ! -d $RUNDIR/../cs$string_cs\kms/$species/$mdot ]; then
							if [ ! -d $RUNDIR/../cs$string_cs\kms/$species ]; then
								if [ ! -d $RUNDIR/../cs$string_cs\kms ]; then
									mkdir $RUNDIR/../cs$string_cs\kms
								fi
								mkdir $RUNDIR/../cs$string_cs\kms/$species
							fi
							mkdir $RUNDIR/../cs$string_cs\kms/$species/$mdot
						fi
						mkdir $RUNDIR/../cs$string_cs\kms/$species/$mdot/data_b${array_b[i]}\_r$r_inner\_r$r_outer
					fi
					mkdir $RUNDIR/../cs$string_cs\kms/$species/$mdot/data_b${array_b[i]}\_r$r_inner\_r$r_outer/incl_$incl
				fi

				# Copy files to this new directory
				cp selfsimilar_solutions photoes submit-job* $RUNDIR/../cs$string_cs\kms/$species/$mdot/data_b${array_b[i]}\_r$r_inner\_r$r_outer/incl_$incl

				# Remove unuseful files
				rm *genmod* photoes selfsimilar_solutions

				# Submit the job on dial or alice
				cd $RUNDIR/../cs$string_cs\kms/$species/$mdot/data_b${array_b[i]}\_r$r_inner\_r$r_outer/incl_$incl
				rm photoes_b${array_b[i]}\_r$r_inner\_r$r_outer\_i$incl.*

				qsub submit-job-alice
				#qsub $RUNDIR/../submit-job-dial
				cd $RUNDIR

			done
		done
	done
done
