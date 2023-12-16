#!/usr/bin/bash
	#/usr/bin/python3
	#1/usr/bin/expect
	#1/usr/bin/python3 -u
if [ $(($# % 3))  -ne 0 ]
then # if num inputs not divisble by 4
	echo -e DOES NOT have a number of inputs divisible by 3."\n"
	exit 1
fi
#echo -e DOES have a number of inputs divisible by 4."\n"

num=$(($#/3))
arg=($@)
for i in $(seq 1 $num)
do
	date	
	python3 "outreader.py" ${arg[((3*$i-3))]} ${arg[((3*$i-2))]} ${arg[((3*$i-1))]}
done
echo -e "All runs plotted.\n"
date
echo

#exec python3 ./helloworld.py
#expect "hello world" { send "stoodis" }
