#!/bin/bash
# script to analyse a bunch of simulations

# ------------- ANALYSIS FOR THE ANALYTICAL SOLUTIONS

RUNDIR=${PWD}

array_b=( 0.75 1.00 1.50 ) #2.00 )
array_ub=( 0.85 0.77 0.56 ) #0.29 )

r_inner=0.1
r_outer=9.5
species=NeII # REMEMBER to change it also in convolution.py
mdot=mdot10e-10 # REMEMBER to change it also in convolution.py

echo "Start analyzing..."
for ((i=0;i<${#array_b[@]};++i)); do
	for  cs in 3.0d5 5.0d5 10.0d5; do
		for incl in 0.0 5.0 10.0 15.0 20.0 25.0 30.0 35.0 40.0 45.0 50.0 55.0 60.0 65.0 70.0 75.0 80.0 85.0 90.0; do

			echo $species" with b="${array_b[i]}", mdot="$mdot " and cs="$cs " for i="$incl

      			# Update the code with the new variables
      			sed -i -e "s/b_input = .*/b_input = ${array_b[i]} /g" convolution.py
      			sed -i -e "s/incl_deg = .*/incl_deg = $incl /g" convolution.py
      			sed -i -e "s/r_inner = .*/r_inner = $r_inner /g" convolution.py
      			sed -i -e "s/r_outer = .*/r_outer = $r_outer /g" convolution.py
			sed -i -e "s/cs = .*/cs = '$cs' /g" convolution.py
      			python convolution.py

			cd $RUNDIR/../cs$cs/$species/$mdot/data_b${array_b[i]}\_r$r_inner\_r$r_outer/incl_$incl
			#cp photoes_b${array_b[i]}\_r$r_inner\_r$r_outer\_i$incl.o* photoes_b${array_b[i]}\_r$r_inner\_r$r_outer\_i$incl.txt
			cp photoes_$species\_b${array_b[i]}\_cs$cs\_i$incl.o* photoes_$species\_b${array_b[i]}\_cs$cs\_i$incl.txt
			#cd $RUNDIR/../data_hydro/$species/incl_$incl
			#cp photoes_hydro_i$incl.o* photoes_hydro_i$incl.txt

      			cd $RUNDIR

		done
	done
done
echo "Have fun with your data!"


# ------------- ANALYSIS FOR THE HYDRO DATA

# RUNDIR=${PWD}
#
# species=NeII # REMEMBER to change it also in convolution.py
#
# echo "Start analyzing..."
# for incl in 0.0 10.0 20.0 45.0 60.0 75.0 90.0; do
#
# 	echo "Hydro: "$species" for i="$incl
#
#       # Update the code with the new variables
#       sed -i -e "s/incl_deg = .*/incl_deg = $incl /g" convolution.py
#       python2 convolution.py
#
# 	cd $RUNDIR/../data_hydro_midplane/$species/incl_$incl
# 	cp photoes_hydro_midplane_i$incl.o* photoes_hydro_i$incl.txt
#
#       cd $RUNDIR
#
# done
# echo "Have fun with your data!"
