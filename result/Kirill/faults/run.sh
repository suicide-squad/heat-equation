#!/bin/bash 
PATH_PROGRAM=$(find ../../../_build -name 'app_Kirill_euler')
echo "$PATH_PROGRAM"
cp -f "$PATH_PROGRAM" .
array=( 8 16 32 64 128 )
for h in ${array[*]}
do
	echo "*********"
	echo $h
	echo "*********"
	sed -i 's/NX=.*/NX='"$h"'/g' setting3.ini
	sed -i 's/NY=.*/NY='"$h"'/g' setting3.ini
	sed -i 's/NZ=.*/NZ='"$h"'/g' setting3.ini
	python3.6 generation3.py

	for (( dt = 3; dt <= 8; dt++ ))
	do
		echo "-------"
		echo $dt
		echo "-------"
		sed -i 's/dt=.*/dt=1e-'"$dt"'/g' setting3.ini
		./app_Kirill_euler
		python3.6 fault.py expected3.txt ../res.txt >> faults_euler.txt
	done
done
rm ../res.txt
