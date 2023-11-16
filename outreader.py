import numpy as np
import matplotlib.pyplot as plt
from time import time
import glob
import pandas as pd
import csv 
import os
 
out_file_name = input("What file from calculations would you like to scrape?: ")
start_time = time()

try:   # IF IT DOESN'T ALREADY EXIST, CREATE A DIR WITH NAME OF FILE (minus .dat) TO SAVE IMGS INTO
	os.mkdir("/home/wolf/NMSSMTools_6.0.0/calculations/{}/".format(out_file_name))
except OSError as error:
	print(error)

threshold_lighthiggs = 50 #GeV

out_file_matrix = list()
events_lighthiggs = list()
ctr_lighthiggs = 0
with open("/home/wolf/NMSSMTools_6.0.0/calculations/{}.dat".format(out_file_name)) as f:
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
#			print("Light Higgs in event {}:".format(indexrow))
#			print(params)
#			print(shiggs)
#			print(phiggs)
#			print()
#	print(out_file_matrix)		
f.close()

plt.figure(1) # TANB VS DIMENSIONLESS PARAMS
plt.scatter([r[19] for r in out_file_matrix], [r[1] for r in out_file_matrix], 
		alpha=0.2, color='b', s=2, label="lambda")
plt.scatter([r[20] for r in out_file_matrix], [r[1] for r in out_file_matrix], 
		alpha=0.2, color='r', s=2, label="kappa")
plt.title("tanB v lambda/kappa : {}".format(out_file_name))
plt.ylabel("tanB")
plt.xlabel("lambda/kappa")
plt.legend()
plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/{}_tanB_v_lambda_kappa.png".format(out_file_name,
												out_file_name))
plt.figure(2) # KAPPA VERSUS LAMBDA
plt.scatter([r[19] for r in out_file_matrix], [r[20] for r in out_file_matrix],
		alpha=0.3, c=[r[1] for r in out_file_matrix], cmap="viridis", s=2)
plt.title("lambda v kappa : {}".format(out_file_name))
plt.xlabel("lambda")
plt.ylabel("kappa")
cbar = plt.colorbar()
cbar.ax.set_ylabel("tanB")
plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/{}_lambda_v_kappa.png".format(out_file_name,out_file_name))

plt.figure(3) # TANB VERSUS DIMENSIONFUL PARAMS
plt.scatter([r[21] for r in out_file_matrix], [r[1] for r in out_file_matrix],
		alpha=0.2, c='b', s=2, label="Alambda")
plt.scatter([r[22] for r in out_file_matrix], [r[1] for r in out_file_matrix],
		alpha=0.2, c='r', s=2, label="Akappa")
plt.scatter([r[23] for r in out_file_matrix], [r[1] for r in out_file_matrix],
		alpha=0.2, c='g', s=2, label="mueff")
plt.title("tanB v mueff/Alambda/Akappa : {}".format(out_file_name))
plt.ylabel("tanB")
plt.xlabel("mueff/Alambda/Akappa")
plt.legend()
plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/{}_tanB_v_mueff_Alambda_Akappa.png".format(out_file_name,out_file_name))

plt.figure(4) # MUEFF VERSUS ALAMBDA/AKAPPA
plt.scatter([r[21] for r in out_file_matrix], [r[23] for r in out_file_matrix],
		alpha=0.2, c='b', s=2, label="Alambda")
plt.scatter([r[22] for r in out_file_matrix], [r[23] for r in out_file_matrix],
		alpha=0.2, c='r', s=2, label="Akappa")
plt.title("mueff v Alambda/Akappa : {}".format(out_file_name))
plt.ylabel("mueff")
plt.xlabel("Alambda/Akappa")
plt.legend()
plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/{}_mueff_v_Alambda_Akappa.png".format(out_file_name,
												out_file_name))

plt.figure(5) # ALAMBDA VS AKAPPA
plt.scatter([r[21] for r in out_file_matrix], [r[22] for r in out_file_matrix],
		alpha=0.2, c='g', s=2)
plt.title("Akappa v Alambda : {}".format(out_file_name))
plt.ylabel("Akappa")
plt.xlabel("Alambda")
plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/{}_Akappa_v_Alambda.png".format(out_file_name,
												out_file_name))

plt.figure(6) # S1MASS VS TANB
plt.scatter([r[24] for r in out_file_matrix], [r[1] for r in out_file_matrix],
		alpha=0.2, c='g', s=2)
plt.title("tanB v S1MASS : {}".format(out_file_name))
plt.ylabel("tanB")
plt.xlabel("Lightest scalar Higgs mass")
plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/{}_tanB_v_s1mass.png".format(out_file_name,
												out_file_name))

plt.figure(7) # P1MASS VS TANB
plt.scatter([r[30] for r in out_file_matrix], [r[1] for r in out_file_matrix],
		alpha=0.2, c='g', s=2)
plt.title("tanB v P1MASS : {}".format(out_file_name))
plt.ylabel("tanB")
plt.xlabel("Lightest pseudoscalar Higgs mass")
plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/{}_tanB_v_p1mass.png".format(out_file_name,
										out_file_name))

plt.figure(8) # AKAPPA VERSUS ALAMBDA/MUEFF
plt.scatter([r[21] for r in out_file_matrix], [r[22] for r in out_file_matrix],
		alpha=0.2, c='b', s=2, label="Alambda")
plt.scatter([r[23] for r in out_file_matrix], [r[22] for r in out_file_matrix],
		alpha=0.2, c='r', s=2, label="mueff")
plt.title("Akappa v Alambda/mueff : {}".format(out_file_name))
plt.ylabel("Akappa")
plt.xlabel("Alambda/mueff")
plt.legend()
plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/{}_Akappa_v_Alambda_mueff.png".format(out_file_name,
												out_file_name))
plt.figure(9) # S3MASS VS TANB
plt.scatter([r[28] for r in out_file_matrix], [r[1] for r in out_file_matrix],
		alpha=0.2, c='g', s=2)
plt.title("tanB v S3MASS : {}".format(out_file_name))
plt.ylabel("tanB")
plt.xlabel("Heaviest scalar Higgs mass")
plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/{}_tanB_v_s3mass.png".format(out_file_name,
												out_file_name))
# IDEA: AKAPPA V ALAMBDA WITH CMAP-MUEFF VALUE
plt.figure(10) # ALAMBDA VS AKAPPA - CMAP MUEFF
plt.scatter([r[21] for r in out_file_matrix], [r[22] for r in out_file_matrix],
		alpha=0.4, c=[r[23] for r in out_file_matrix], cmap='magma', s=2)
plt.title("Akappa v Alambda : {}".format(out_file_name))
plt.ylabel("Akappa")
plt.xlabel("Alambda")
cbar = plt.colorbar()
cbar.ax.set_ylabel("mueff")
plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/{}_Akappa_v_Alambda_1.png".format(out_file_name,
												out_file_name))

plt.figure(11) # ALAMBDA VS AKAPPA - CMAP TANB
plt.scatter([r[21] for r in out_file_matrix], [r[22] for r in out_file_matrix],
		alpha=0.4, c=[r[1] for r in out_file_matrix],cmap='viridis',s=2)
plt.title("Akappa v Alambda : {}".format(out_file_name))
plt.ylabel("Akappa")
plt.xlabel("Alambda")
cbar = plt.colorbar()
cbar.ax.set_ylabel("tanB")
plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/{}_Akappa_v_Alambda_2.png".format(out_file_name,
												out_file_name))












#sortedbys1mass = sorted(out_file_matrix, key = lambda x: x[24])
#sortedbyp1mass = sorted(out_file_matrix, key = lambda x: x[30])

print(ctr_lighthiggs)
print(time()-start_time)
