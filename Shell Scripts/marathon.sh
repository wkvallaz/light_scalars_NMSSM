#!/usr/bin/bash -e
	#/usr/bin/python3
	#1/usr/bin/expect
	#1/usr/bin/python3 -u
# shell script if need to run several nmssm inp files
# REQUIRES make init AND make TO ALREADY BE RUN BEFORE
#for run in "wH" "wh" "wd" "wR" "wr" "wK"
if [ -f "/home/wolf/NMSSMTools_6.0.0/main/nmhdecay_rand.f" ]
then	
	echo Initiating runs...
else
	echo Initializing...
	make init
	make
	echo Initiating runs...
fi

for run in "$@"
do
	for con in "THY" "LEP" "LHC" "BKF"
	do
		date
		./run calculations/$run$con\randinp.dat
	done
done
date
echo -e Runs complete."\n- --- - --- - --- - --- -"
