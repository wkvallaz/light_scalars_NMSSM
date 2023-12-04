#!/usr/bin/bash
# user runs as ./duplicateinp.dat old_file_prefix new_file_prefix
#					$1 		$2
for TAG in "THY" "LEP" "LHC" "BKF"
do
	cp $1$TAG\randinp.dat $2$TAG\randinp.dat
done
