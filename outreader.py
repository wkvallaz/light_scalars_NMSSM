import numpy as np
import matplotlib.pyplot as plt
from time import time
import glob
import pandas as pd
import csv 
start_time = time()

#out_file_name = input("What file from calculations would you like to scrape?: ")
#BRUTE FORCE LAZINESS
out_file_name = "test2randout.dat"

threshold_lighthiggs = 50 #GeV

out_file_matrix = list()
events_lighthiggs = list()
ctr_lighthiggs = 0
with open("/home/wolf/NMSSMTools_6.0.0/calculations/{}".format(out_file_name)) as f:
	f_reader = csv.reader(f, delimiter=" ")
	
	for indexrow,fullrow in enumerate(f_reader):
		row = [indexrow] # trim out strange spacing
		for val in fullrow:
			if val != "": row.append(float(val))
		out_file_matrix.append(row)

		params = row[1:24]
		shiggs = row[24:30] # s1mass s1comp s2mass s2comp s3mass s3comp
		phiggs = row[30:34] # p1mass p1comp p2mass p2comp
		if (float(shiggs[0]) < threshold_lighthiggs) or (float(phiggs[0]) < threshold_lighthiggs):
			ctr_lighthiggs += 1
			print("Light Higgs in event {}:".format(indexrow))
			print(params)
			print(shiggs)
			print(phiggs)
			print()
#	print(out_file_matrix)		
f.close()

sortedbys1mass = sorted(out_file_matrix, key = lambda x: x[24])
sortedbyp1mass = sorted(out_file_matrix, key = lambda x: x[30])


#print(sorted(out_file_matrix, key = lambda x: x[24]))
#for r in sortedbys1mass:
#	print(r[24])
for r in sortedbyp1mass:
	print(r[30])

print(ctr_lighthiggs)
print(time()-start_time)
