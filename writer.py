from time import time
import numpy as np
import csv 
import sys

configname="F2Ds1"
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
#with open("{}{}".format(read_dir_name, "NMSSM_big_valid_events_list.txt")) as f: #050607
#	f_reader = csv.reader(f, delimiter=" ")
#	for i,r in enumerate(f_reader):
#		NMSSM.append(#[0]+
#				[float(x) for x in r if x!=""]) #my outreader.py output has the 0s here already
#		if i%1e50<1:print(#[0]+
#				[float(x) for x in r if x!=""],"\n")
#f.close()
#with open("{}{}".format(read_dir_name, "NMSSM_big_valid_events_list_v2.txt")) as f: #v1020304
#	f_reader = csv.reader(f, delimiter=" ")
#	for i,r in enumerate(f_reader):
#		NMSSM.append(#[0]+
#				[float(x) for x in r if x!=""]) #my outreader.py output has the 0s here already
#		if i<1:print(#[0]+
#				[float(x) for x in r if x!=""],"\n")
#f.close()
#with open("{}{}".format(read_dir_name, "NMSSM_big_valid_events_list_v3.txt")) as f: #FORESEE2D
#with open("{}{}".format(read_dir_name, "NMSSM_big_valid_events_list_v4.txt")) as f: #FORESEE2Dv2
with open("{}{}".format(read_dir_name, f"NMSSM_big_valid_events_list_{configname}.txt")) as f: #FORESEE2Dv#3
	f_reader = csv.reader(f, delimiter=" ")
	for i,r in enumerate(f_reader):
		NMSSM.append(#[0]+
				[float(x) for x in r if x!=""]) #my outreader.py output has the 0s here already
		if i%1e50<1:print(#[0]+
				[float(x) for x in r if x!=""],"\n")
f.close()

print(f"This NMSSMTools out.dat has {len(NMSSM)} events.")
#END OF SOLID ARCH
cslambda=0
cskappa=0
csalambda=0
csakappa=0
csmueff=0
for i,r in enumerate(NMSSM):
	if i%1e50<1: print(r)
	cslambda+=r[19]/len(NMSSM)
	cskappa+=r[20]/len(NMSSM)
	csalambda+=r[21]/len(NMSSM)
	csakappa+=r[22]/len(NMSSM)
	csmueff+=r[23]/len(NMSSM)
print(cslambda,cskappa,csalambda,csakappa,csmueff)
#sys.exit()
#RESUME SOLID ARCH
outfile = list()
for r in NMSSM:
	if r[36] < 0.100 or 10. < r[36]: continue
	#		mA    tanB ctau=lifetime (m)    Acomp    mH+
	outfile.append([r[36],r[1],1.975e-16/r[135], r[37]**2, r[42]])
print(f"It has {len(outfile)} events with pseudoscalar mass [100 MeV, 10 GeV].")

outfile = sorted(outfile, key = lambda x: x[0])

printres=5000 # 10 low res ie. print all of them, 100000000 hi res, etc.

f = open("{}{}".format(save_dir_name, writ_name), "w")
for i,r in enumerate(outfile):# outcopy is neighboring averages
#	line = "{:.2e} {:.5e}\n".format(m,1.975e-16/ctau)
	line = f"{r[0]} {r[1]} {r[2]} {r[3]} {r[4]}\n"
	if int(10*printres*(i/len(outfile)))%printres==0: print(i,"\t",line[:-1])
	f.write(line)
f.close()


