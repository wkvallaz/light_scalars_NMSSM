from time import time
import numpy as np
import csv 
import sys

configname="F2Ds1win1_0-10"
writ_name = f"NMSSM_ctau2D_{configname}.txt"
read_dir_name = "/home/wolf/NMSSMTools_6.0.0/calculations/"
save_dir_name = "/mnt/c/Users/Wolfgang/FORESEE-main/FORESEE-main/Models/NMSSM/model/"

NEIGHBOR_SMOOTHING = False

def roundnearest(num,to):
	ret = round(num/to)*to
	return ret
def closeto(a,b): # if a is close enough to b (0.05% window)
	ret = (0.9995*a<b and b<1.0005*a)
	return ret

NMSSM = list()
with open("{}{}".format(read_dir_name, f"NMSSM_big_valid_events_list_{configname}.txt")) as f:
#with open("{}{}".format(read_dir_name, f"NMSSM_big_valid_events_list_lighthiggs0-4.txt")) as f:
	f_reader = csv.reader(f, delimiter=" ")
	for i,r in enumerate(f_reader):
		NMSSM.append(#[0]+
				[float(x) for x in r if x!=""]) #my outreader.py output has the 0s here already
		if i%1e50<1:print(len([float(x) for x in r if x!=""]),
				[float(x) for x in r if x!=""],"\n")
f.close()
#with open("{}{}".format(read_dir_name, f"NMSSM_big_valid_events_list_{configname}part2.txt")) as f:
#with open("{}{}".format(read_dir_name, f"NMSSM_big_valid_events_list_lighthiggs5-7.txt")) as f:
#	f_reader = csv.reader(f, delimiter=" ")
#	for i,r in enumerate(f_reader):
#		NMSSM.append(#[0]+
#				[float(x) for x in r if x!=""]) #my outreader.py output has the 0s here already
#		if i%1e50<1:print(#[0]+
#				[float(x) for x in r if x!=""],"\n")
#f.close()

print(f"This NMSSMTools out.dat has {len(NMSSM)} events.")

NMSSM = sorted(NMSSM, key = lambda x:x[36])

printres=len(NMSSM)

f = open("{}{}".format(save_dir_name, writ_name), "w")
for i,r in enumerate(NMSSM):
	line = str(1.975e-16/r[135])
	for e in r[1:]:
		line+=" "+str(e)
	line+="\n"
	if (i)%printres==0: print(i,"\t",line[:-1])

	f.write(line)
f.close()


