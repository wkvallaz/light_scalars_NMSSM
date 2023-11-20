import numpy as np
import matplotlib.pyplot as plt
from time import time
import glob
import pandas as pd
import csv 
import os
import sys

############ EFFORTS TO COMBINE THE CONSTRAINED AND UNCONSTRAINED ALL ONTO A SINGLE PLOT
 
save_dir_name = input("Where would you like to store the scrape?: ")
start_time = time()
DIR = "/home/wolf/NMSSMTools_6.0.0/calculations/"
try:   # IF IT DOESN'T ALREADY EXIST, CREATE A DIR WITH NAME OF FILE (minus .dat) TO SAVE IMGS INTO
	os.mkdir("{}{}/".format(DIR, save_dir_name))
except OSError as error:
	print(error)
	overwrite_req = input("Continuing will overwrite existing files, are you sure? (y/ye/yea/yes): ").lower()
	if overwrite_req not in ["y", "ye", "yea", "yes"]: sys.exit("Execution halted.")
file_names = ["wide3randout","wide3con1randout","wide3con3randout","wide3con2randout"]

threshold_lighthiggs = 20 #GeV

def GeneratePlots():#(out_file_name, file_index, SAVEFIGS):#####NEXT STEP, PUT A LOOP AROUND EACH PLOT SO IT CAN BE DELETED ONCE SAVEFIG'D
	
	Label = file_names
	#Color = ['darkgray', 'cyan', 'yellow', 'magenta']		#CYM COLORING
	Color = ['black', 'b', 'g', 'r']			#RGB COLORING
	Alpha = [1,1,1,1]   
	Size = [5,5,2,.5] #CONCENTRIC SIZING
	#Size = [0,6,5,4] #ALPHA-STACK SIZING
	DPI = 480
	#CURRENT CONSIDERATIONS ARE TO PLOT A SINGULAR TIME FROM THE UNCONSTRAINED SET, AND TAG EACH POINT WITH
	# A CERTAIN COLOR BASED ON WHICH CONSTRAINTS IT SURVIVES	
	#color_list = constraints_survived #alias list of colors survived	

	pltctr = 1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) #TANB VS LAMBDA
	out_file_matrix = file_matrices[0]
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[19] for r in out_file_matrix], [r[1] for r in out_file_matrix], 
			alpha=Alpha[fx], color=Color[fx], s=Size[fx], label=Label[fx], marker='.') 
	plt.title("tanB v lambda")
	plt.ylabel("tanB")
	plt.xlabel("lambda")
	plt.xlim(0,.8)
	plt.legend()
	plt.savefig("{}{}/tanB_v_lambda.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()

	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # TANB VS KAPPA
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[20] for r in out_file_matrix], [r[1] for r in out_file_matrix], 
			alpha=Alpha[fx], color=Color[fx], s=Size[fx], label=Label[fx], marker='.')	 
	plt.title("tanB v kappa")
	plt.ylabel("tanB")
	plt.xlabel("kappa")
	plt.xlim(0,.8)
	plt.legend()
	plt.savefig("{}{}/tanB_v_kappa.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()

	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # KAPPA VERSUS LAMBDA
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[19] for r in out_file_matrix], [r[20] for r in out_file_matrix],
			alpha=Alpha[fx], color=Color[fx], s=Size[fx], label=Label[fx], marker='.')
#			alpha=0.3, c=[r[1] for r in out_file_matrix], cmap="viridis", s=2)
	plt.title("lambda v kappa")
	plt.xlabel("lambda")
	plt.ylabel("kappa")
#	cbar = plt.colorbar()
#	cbar.ax.set_ylabel("tanB")
	plt.ylim(0,.8)
	plt.xlim(0,.8)
	plt.legend()
	plt.savefig("{}{}/lambda_v_kappa.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()
	
	pltctr+=1 
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # TANB VERSUS ALAMBDA
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[21] for r in out_file_matrix], [r[1] for r in out_file_matrix],
			alpha=Alpha[fx], color=Color[fx], s=Size[fx], label=Label[fx], marker='.')
	plt.title("tanB v Alambda")
	plt.ylabel("tanB")
	plt.xlabel("Alambda")
	plt.legend()
	plt.savefig("{}{}/tanB_v_Alambda.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()

	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # TANB VERSUS AKAPPA
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[22] for r in out_file_matrix], [r[1] for r in out_file_matrix],
			alpha=Alpha[fx], c=Color[fx], s=Size[fx], label=Label[fx], marker='.')
	plt.title("tanB v Akappa")
	plt.ylabel("tanB")
	plt.xlabel("Akappa")
	plt.legend()
	plt.savefig("{}{}/tanB_v_Akappa.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()

	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # TANB VERSUS MUEFF
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[23] for r in out_file_matrix], [r[1] for r in out_file_matrix],
			alpha=Alpha[fx], c=Color[fx], s=Size[fx], label=Label[fx], marker='.')
	plt.title("tanB v mueff")
	plt.ylabel("tanB")
	plt.xlabel("mueff")
	plt.legend()
	plt.savefig("{}{}/tanB_v_mueff.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()

	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # MUEFF VERSUS ALAMBDA
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[21] for r in out_file_matrix], [r[23] for r in out_file_matrix],
			alpha=Alpha[fx], c=Color[fx], s=Size[fx], label=Label[fx], marker='.')
	
	plt.title("mueff v Alambda")
	plt.ylabel("mueff")
	plt.xlabel("Alambda")
	plt.legend()
	plt.savefig("{}{}/mueff_v_Alambda.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()
	
	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # MUEFF VERSUS AKAPPA
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[22] for r in out_file_matrix], [r[23] for r in out_file_matrix],
			alpha=Alpha[fx], c=Color[fx], s=Size[fx], label=Label[fx], marker='.')

	plt.title("mueff v Akappa")
	plt.ylabel("mueff")
	plt.xlabel("Akappa")
	plt.legend()
	plt.savefig("{}{}/mueff_v_Akappa.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()
	
	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # ALAMBDA VS AKAPPA
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[21] for r in out_file_matrix], [r[22] for r in out_file_matrix],
			alpha=Alpha[fx], c=Color[fx], s=Size[fx], label=Label[fx], marker='.')
	plt.title("Akappa v Alambda")
	plt.ylabel("Akappa")
	plt.xlabel("Alambda")
	plt.legend()
	plt.savefig("{}{}/Akappa_v_Alambda.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()
	
	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # S1MASS VS TANB
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[24] for r in out_file_matrix], [r[1] for r in out_file_matrix],
			alpha=Alpha[fx], c=Color[fx], s=Size[fx], label=Label[fx], marker='.')
	plt.title("tanB v S1MASS")
	plt.ylabel("tanB")
	plt.xlabel("Lightest scalar Higgs mass")
	plt.legend()
	plt.savefig("{}{}/tanB_v_s1mass.png".format(DIR, save_dir_name), dpi=DPI)
	plt.xlim(112.5,132.5)
	plt.savefig("{}{}/tanB_v_s1mass_zoom.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()	

	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # S1MASS vs AKAPPA
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[24] for r in out_file_matrix], [r[22] for r in out_file_matrix],
			alpha=Alpha[fx], c=Color[fx], s=Size[fx], label=Label[fx], marker='.')
	plt.title("Akappa v S1MASS")
	plt.ylabel("Akappa")
	plt.xlabel("Lightest scalar Higgs mass")
	plt.legend()
	plt.savefig("{}{}/Akappa_v_s1mass.png".format(DIR, save_dir_name), dpi=DPI)
	plt.xlim(112.5,132.5)
	plt.savefig("{}{}/Akappa_v_s1mass_zoom.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()

	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # S1MASS vs ALAMBDA
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[24] for r in out_file_matrix], [r[21] for r in out_file_matrix],
			alpha=Alpha[fx], c=Color[fx], s=Size[fx], label=Label[fx], marker='.')
	plt.title("Alambda v S1MASS")
	plt.ylabel("Alambda")
	plt.xlabel("Lightest scalar Higgs mass")
	plt.legend()
	plt.savefig("{}{}/Alambda_v_s1mass.png".format(DIR, save_dir_name), dpi=DPI)
	plt.xlim(112.5,132.5)
	plt.savefig("{}{}/Alambda_v_s1mass_zoom.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()

	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # S2MASS VS TANB
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[26] for r in out_file_matrix], [r[1] for r in out_file_matrix],
			alpha=Alpha[fx], c=Color[fx], s=Size[fx], label=Label[fx], marker='.')
	plt.title("tanB v S2MASS")
	plt.ylabel("tanB")
	plt.xlabel("Second scalar Higgs mass")
	plt.legend()
	plt.savefig("{}{}/tanB_v_s2mass.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()

	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # S2MASS VS AKAPPA
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[26] for r in out_file_matrix], [r[22] for r in out_file_matrix],
			alpha=Alpha[fx], c=Color[fx], s=Size[fx], label=Label[fx], marker='.')
	plt.title("Akappa v S2MASS")
	plt.ylabel("Akappa")
	plt.xlabel("Second scalar Higgs mass")
	plt.legend()
	plt.savefig("{}{}/Akappa_v_s2mass.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()
	
	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # S2MASS VS ALAMBDA
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[26] for r in out_file_matrix], [r[21] for r in out_file_matrix],
			alpha=Alpha[fx], c=Color[fx], s=Size[fx], label=Label[fx], marker='.')
	plt.title("Alambda v S2MASS")
	plt.ylabel("Alambda")
	plt.xlabel("Second scalar Higgs mass")
	plt.legend()
	plt.savefig("{}{}/Alambda_v_s2mass.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()
	
	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # S3MASS VS TANB
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[28] for r in out_file_matrix], [r[1] for r in out_file_matrix],
			alpha=Alpha[fx], c=Color[fx], s=Size[fx], label=Label[fx], marker='.')
	plt.title("tanB v S3MASS")
	plt.ylabel("tanB")
	plt.xlabel("Heaviest scalar Higgs mass")
	plt.legend()
	plt.savefig("{}{}/tanB_v_s3mass.png".format(DIR, save_dir_name), dpi=DPI)
	plt.xlim(0,35000)
	plt.savefig("{}{}/tanB_v_s3mass_zoom.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()
	
	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # S3MASS VS AKAPPA
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[28] for r in out_file_matrix], [r[22] for r in out_file_matrix],
			alpha=Alpha[fx], c=Color[fx], s=Size[fx], label=Label[fx], marker='.')
	plt.title("Akappa v S3MASS")
	plt.ylabel("Akappa")
	plt.xlabel("Heaviest scalar Higgs mass")
	plt.legend()
	plt.savefig("{}{}/Akappa_v_s3mass.png".format(DIR, save_dir_name), dpi=DPI)
	plt.xlim(0,35000)
	plt.savefig("{}{}/Akappa_v_s3mass_zoom.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()
	
	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # S3MASS VS ALAMBDA
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[28] for r in out_file_matrix], [r[21] for r in out_file_matrix],
			alpha=Alpha[fx], c=Color[fx], s=Size[fx], label=Label[fx], marker='.')
	plt.title("Alambda v S3MASS")
	plt.ylabel("Alambda")
	plt.xlabel("Heaviest scalar Higgs mass")
	plt.legend()
	plt.savefig("{}{}/Alambda_v_s3mass.png".format(DIR, save_dir_name), dpi=DPI)
	plt.xlim(0,35000)
	plt.savefig("{}{}/Alambda_v_s3mass_zoom.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()
	
	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # P1MASS VS TANB
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[30] for r in out_file_matrix], [r[1] for r in out_file_matrix],
			alpha=Alpha[fx], c=Color[fx], s=Size[fx], label=Label[fx], marker='.')	
	plt.ylabel("tanB")
	plt.title("tanB v P1MASS")
	plt.xlabel("Lightest pseudoscalar Higgs mass")	
	plt.legend()
	plt.savefig("{}{}/tanB_v_p1mass.png".format(DIR, save_dir_name), dpi=DPI)
	plt.xlim(0,10000)
	plt.savefig("{}{}/tanB_v_p1mass_zoom".format(DIR, save_dir_name), dpi=DPI)
	plt.close()

	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # P1MASS VS ALAMBDA
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[30] for r in out_file_matrix], [r[21] for r in out_file_matrix],
			alpha=Alpha[fx], c=Color[fx], s=Size[fx], label=Label[fx], marker='.')	
	plt.ylabel("Alambda")
	plt.title("Alambda v P1MASS")
	plt.xlabel("Lightest pseudoscalar Higgs mass")	
	plt.legend()
	plt.savefig("{}{}/Alambda_v_p1mass.png".format(DIR, save_dir_name), dpi=DPI)
	plt.xlim(0,10000)
	plt.savefig("{}{}/Alambda_v_p1mass_zoom".format(DIR, save_dir_name), dpi=DPI)
	plt.close()

	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # P1MASS VS AKAPPA
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[30] for r in out_file_matrix], [r[22] for r in out_file_matrix],
			alpha=Alpha[fx], c=Color[fx], s=Size[fx], label=Label[fx], marker='.')	
	plt.ylabel("Akappa")
	plt.title("Akappa v P1MASS")
	plt.xlabel("Lightest pseudoscalar Higgs mass")	
	plt.legend()
	plt.savefig("{}{}/Akappa_v_p1mass.png".format(DIR, save_dir_name), dpi=DPI)
	plt.xlim(0,10000)
	plt.savefig("{}{}/Akappa_v_p1mass_zoom".format(DIR, save_dir_name), dpi=DPI)
	plt.close()

	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # P2MASS VS ALAMBDA
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[32] for r in out_file_matrix], [r[21] for r in out_file_matrix],
			alpha=Alpha[fx], c=Color[fx], s=Size[fx], label=Label[fx], marker='.')	
	plt.ylabel("Alambda")
	plt.title("Alambda v P2MASS")
	plt.xlabel("Heaviest pseudoscalar Higgs mass")
	plt.legend()
	plt.savefig("{}{}/Alambda_v_p2mass.png".format(DIR, save_dir_name), dpi=DPI)
	plt.xlim(0,31000)
	plt.savefig("{}{}/Alambda_v_p2mass_zoom.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()

	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # P2MASS VS AKAPPA
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[32] for r in out_file_matrix], [r[22] for r in out_file_matrix],
			alpha=Alpha[fx], c=Color[fx], s=Size[fx], label=Label[fx], marker='.')	
	plt.ylabel("Akappa")
	plt.title("Akappa v P2MASS")
	plt.xlabel("Heaviest pseudoscalar Higgs mass")
	plt.legend()
	plt.savefig("{}{}/Akappa_v_p2mass.png".format(DIR, save_dir_name), dpi=DPI)
	plt.xlim(0,31000)
	plt.savefig("{}{}/Akappa_v_p2mass_zoom.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()

	pltctr+=1
	print("Plotting #{}".format(pltctr))
	plt.figure(pltctr) # P2MASS VS TANB
	for fx,out_file_matrix in enumerate(file_matrices):
		plt.scatter([r[32] for r in out_file_matrix], [r[1] for r in out_file_matrix],
			alpha=Alpha[fx], c=Color[fx], s=Size[fx], label=Label[fx], marker='.')		
	plt.ylabel("tanB")
	plt.title("tanB v P2MASS")
	plt.xlabel("Heaviest pseudoscalar Higgs mass")
	plt.legend()
	plt.savefig("{}{}/tanB_v_p2mass.png".format(DIR, save_dir_name), dpi=DPI)
	plt.xlim(0,31000)
	plt.savefig("{}{}/tanB_v_p2mass_zoom.png".format(DIR, save_dir_name), dpi=DPI)
	plt.close()
	



	print("Finished plots.")

# FOR MIN,MAX MASS FILTERED BY con2
(mins1m,maxs1m) = (140,100) # declared as such so 1st data point forces min/max to resolve
mh_is1 = [0,0,0,0] # elem corresp to file_index from enumerating file_names 
mh_is2 = [0,0,0,0]
mh_dne = [0,0,0,0]
mh_1n2 = [0,0,0,0]
mh_is3 = [0,0,0,0]

def NearSM(mass): #temp fnc that checks if mass is near sm higgs mass
	return ((mass > 120) and (mass < 130))

file_matrices = [list(),list(),list(),list()] # for storing "out_file_matrix" of each out file
#constraints_survived = list()# for each element of unconstrained set, note which points survive later constraints

for file_index,out_file_name in enumerate(file_names):
	with open("{}{}.dat".format(DIR, out_file_name)) as f:
		f_reader = csv.reader(f, delimiter=" ")
		out_file_matrix = file_matrices[file_index]
		ctr_lighthiggs = 0
		for indexrow,fullrow in enumerate(f_reader):
			row = [indexrow] # trim out strange spacing
			for val in fullrow:
				if val != "": row.append(float(val))
			out_file_matrix.append(row)

			params = row[1:24]
			shiggs = row[24:30] # s1mass s1comp s2mass s2comp s3mass s3comp
			phiggs = row[30:34] # p1mass p1comp p2mass p2comp

			if (float(shiggs[0])<threshold_lighthiggs) or (float(phiggs[0])<threshold_lighthiggs):
				ctr_lighthiggs += 1			
		#		print("Light Higgs in event {}:".format(indexrow))
		#		print("params:\t", params)
		#		print("shiggs:\t", shiggs)
		#		print("phiggs:\t", phiggs)
		#		print()
		
			# tracking which events have NearSM higgs in s1, s2, both, s3, or none
			if (NearSM(shiggs[0]) and NearSM(shiggs[2])): 
				mh_1n2[file_index]+=1
				print("Wait, really?{}:event{}:s1@{}:s2@{}".format(
					out_file_name,indexrow,shiggs[0],shiggs[2]))
			elif (NearSM(shiggs[0])):
				mh_is1[file_index]+=1
			elif (NearSM(shiggs[2])):
				mh_is2[file_index]+=1
				print("{}event#{}:\ts1@{}\s2@{}".format(out_file_name,indexrow,shiggs[0],shiggs[2]))
			elif (NearSM(shiggs[4])):
				mh_is3[file_index]+=1
				print("No shot...{}:event{}:s3@{}".format(
					out_file_name,indexrow,shiggs[4]))
			else: mh_dne[file_index]+=1


			# ALSO, JUST GET THE MAX AND MIN VALUES ALLOWED BY con2 FOR SM-mh
			if ("con2" in out_file_name):
				(mins1m,maxs1m) = (min(mins1m,shiggs[0]),max(maxs1m,shiggs[0]))

	f.close()
#for index,event in enumerate(file_matrices[0]):	#each event in the unconstrained set...
#	event_1 = 0
#	event_2 = 0
#	event_3 = 0
#	if event in file_matrices[1]: event_3 = 200/255
#	if event in file_matrices[2]: event_2 = 200/255
#	if event in file_matrices[3]: event_1 = 200/255
#	constraints_survived.append([event_1,event_2,event_3])
GeneratePlots() 

#sortedbys1mass = sorted(out_file_matrix, key = lambda x: x[24])
#sortedbyp1mass = sorted(out_file_matrix, key = lambda x: x[30])

print("\ncon2 min/max s1m:\t",mins1m,"\t",maxs1m)

print("\nFILE_NAME\tmh_is1\tmh_is2\tmh_dne\tmh_1&2\tmh_is3")
for file_index,out_file_name in enumerate(file_names):
	print("{}\t{}\t{}\t{}\t{}\t{}".format(out_file_name[:-3], 
	mh_is1[file_index],mh_is2[file_index],mh_dne[file_index],mh_1n2[file_index],mh_is3[file_index]))

print("\n# Light Higgs:\t", ctr_lighthiggs)
print("Runtime(s):\t",time()-start_time)
