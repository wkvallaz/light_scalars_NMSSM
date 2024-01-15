#!/usr/bin/bash -e
	#/usr/bin/python3
	#1/usr/bin/expect
	#1/usr/bin/python3 -u
# shell script if need to run several nmssm inp files
# REQUIRES make init AND make TO ALREADY BE RUN BEFORE
#for run in "wH" "wh" "wd" "wR" "wr" "wK"

# copy of marathon.sh just to run single files an have date be run around them

echo "Initiating run(s)..."
for runcon in "$@"
do
	date
	echo "- --- - ----- - --- - ----- - --- - ----- - --- -"	#note to self: just SUSY con, 2.6s/point
	./run calculations/$runcon\randinp.dat
	echo "- --- - ----- - --- - ----- - --- - ----- - --- -"
done
date
echo -e Marathon complete."\n= === = ===== = === = ===== = === = ===== = === ="
