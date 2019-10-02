#!/bin/bash
# script to analyse a bunch of simulations

RUNDIR=${PWD}

array_b=( 0.75 1.00 1.50 ) #2.00 )
array_ub=( 0.85 0.77 0.56 ) #0.29 )

r_inner=0.1
r_outer=9.5
species=NeII # REMEMBER to change it also in convoluption.py
mdot=mdot10e-9 # REMEMBER to change it also in convolution.py

for ((i=0;i<${#array_b[@]};++i)); do
	for  cs in 3.0d5 5.0d5 10.0d5; do
		for incl in 0.0 2.5 5.0 7.5 10.0 12.5 15.0 17.5 20.0 22.5 25.0 27.5 30.0 32.5 35.0 37.5 40.0 42.5 45.0 47.5 50.0 52.5 55.0 57.5 60.0 62.5 65.0 67.5 70.0 72.5 75.0 77.5 80.0 82.5 85.0 87.5 90.0; do

      			echo "up to b="${array_b[i]}", R_out="$r_outer "and R_in="$r_inner "for i="$incl

        		# Update the code with the new variables
        		sed -i -e "s/b_input = .*/b_input = ${array_b[i]} /g" convolution.py
        		sed -i -e "s/incl_deg = .*/incl_deg = $incl /g" convolution.py
        		sed -i -e "s/r_inner = .*/r_inner = $r_inner /g" convolution.py
        		sed -i -e "s/r_outer = .*/r_outer = $r_outer /g" convolution.py
			sed -i -e "s/cs = .*/cs = $cs /g" convolution.py
        		python convolution.py

			cd $RUNDIR/../cs$cs\kms/$species/$mdot/data_b${array_b[i]}\_r$r_inner\_r$r_outer/incl_$incl
			#cp photoes_b${array_b[i]}\_r$r_inner\_r$r_outer\_i$incl.o* photoes_b${array_b[i]}\_r$r_inner\_r$r_outer\_i$incl.txt
			cp photoes_$species\_b${array_b[i]}\_cs$cs\_i$incl.o* photoes_$species\_b${array_b[i]}\_cs$cs\_i$incl.txt
			#cd $RUNDIR/../data_hydro/$species/incl_$incl
			#cp photoes_hydro_i$incl.o* photoes_hydro_i$incl.txt

        		cd $RUNDIR

		done
	done
done
