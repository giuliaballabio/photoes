#!/bin/bash
# script to run a grid of simulations

RUNDIR=${PWD}

array_b=( 0.75 1.00 1.50 2.00 )
array_ub=( 0.85 0.77 0.56 0.29 )

cs=10
species=NeII
mdot=mdot10e-8

for ((i=0;i<${#array_b[@]};++i)); do
	# echo "(${array_b[i]}, ${array_ub[i]})"
	for r_inner in 0.1 1.0; do
		for r_outer in 5.0 9.5; do
			for incl in 0.0 45.0 90.0; do

        			echo "up to b="${array_b[i]}", R_out="$r_outer "and R_in="$r_inner "for i="$incl

        			# Update the code with the new variables
        			sed -i -e "s/b_input = .*/b_input = ${array_b[i]} /g" convolution.py
        			sed -i -e "s/incl_deg = .*/incl_deg = $incl /g" convolution.py
        			sed -i -e "s/r_inner = .*/r_inner = $r_inner /g" convolution.py
        			sed -i -e "s/r_outer = .*/r_outer = $r_outer /g" convolution.py
							sed -i -e "s/cs = .*/cs = $cs /g" convolution.py
        			python convolution.py

        			cd $RUNDIR/../cs$cs\kms/$species/data_b${array_b[i]}\_r$r_inner\_r$r_outer/incl_$incl
							cp photoes_b${array_b[i]}\_r$r_inner\_r$r_outer\_i$incl.o* photoes_b${array_b[i]}\_r$r_inner\_r$r_outer\_i$incl.txt
							# cd $RUNDIR/../data_hydro/$species/incl_$incl
							# cp photoes_hydro_i$incl.o* photoes_hydro_i$incl.txt

        			# cd $RUNDIR/../plot
        			# python plot_observables.py

        			cd $RUNDIR

			done
		done
	done
done
