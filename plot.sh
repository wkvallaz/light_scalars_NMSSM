#!/usr/bin/bash
	#/usr/bin/python3
	#1/usr/bin/expect
	#1/usr/bin/python3 -u

for ppsc in "wH wHdir- 1 1" "wh whdir 1 1" "wd wddir 1 1" "wR wRdir- 1 1" "wr wrdir 1 1" "wK wKdir 1 1"
do
	pwd
	ls
	date
	#echo python3 "outreader.py" $ppsc
	python3 "outreader.py" $ppsc	
done
date
echo finito

#exec python3 ./helloworld.py
#expect "hello world" { send "stoodis" }
