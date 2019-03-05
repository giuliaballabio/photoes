#!/bin/bash
# script to run a grid of simulations

RUNDIR=${PWD}

array_b=( 0.75 1.00 1.50 2.00 )
array_ub=( 0.85 0.77 0.56 0.29 )

for ((i=0;i<${#array_b[@]};++i)); do
  # echo "(${array_b[i]}, ${array_ub[i]})"
	for r_inner in 0.1; do
		for r_outer in 0.9 1.0 1.5; do
			for incl in 0.0 1.0 5.0 10.0 20.0 27.0 35.0 45.0 50.0 60.0 68.0 75.0 82.0 90.0; do

        echo "up to b="${array_b[i]} ",R_out="$r_outer "and R_in="$r_inner "for i="$incl

        # Update the code with the new variables
        sed -i -e "s/b_input = .*/b_input = ${array_b[i]} /g" convolution.py
        sed -i -e "s/incl_deg = .*/incl_deg = $incl /g" convolution.py
        sed -i -e "s/r_inner = .*/r_inner = $r_inner /g" convolution.py
        sed -i -e "s/r_outer = .*/r_outer = $r_outer /g" convolution.py
        python convolution.py

        cd $RUNDIR/../data_b${array_b[i]}\_r$r_inner\_r$r_outer/incl_$incl
        cp photoes_b${array_b[i]}\_r$r_inner\_r$r_outer\_i$incl.o* photoes_b${array_b[i]}\_r$r_inner\_r$r_outer\_i$incl.txt

        # cd $RUNDIR/../plot
        # python plot_observables.py
        
        cd $RUNDIR

			done
		done
	done
done
