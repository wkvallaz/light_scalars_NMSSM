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
num_files = len(file_names)
threshold_lighthiggs = 20 #GeV

def SinglePlot(pltctr, xpar, xind, xmin, xmax, ypar, yind, ymin, ymax,  
		Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI): 
# args ^      (   int,  str,  int,  int,  int,  str,  int,  int,  int, 
#		  arr,  arr,    arr,  arr, str,            arr, int):
	#needs to accept a list of args:
	#	pltctr		# plot count for figuure - considering subplots/axs?
	#	xmin, xmax 	# if xmin!=xmax do zoomed plot also
	#	ymin, ymax	#    ^ analogous
	# 	xpar		# x parameter
	#	ypar		# y parameter 
	#	xind		# x par corresponding column number
	#	yind		# y par corresponding column number
	#	Color 		# Color array
	#	Alpha		# Alpha array
	#	Size		# Size array
	#	Label		# Label array
	#	LOC		# Legend location
	#	BBOX_TO_ANCHOR	# Legend offset from location
	#	DPI		# Figure dpi
	
	DIR = "/home/wolf/NMSSMTools_6.0.0/calculations/"+save_dir_name+"/"
	print("Plotting #{} ...".format(pltctr))

	def Col(index,matrix): #array generator statement which pulls column 'index' from  ea row of 'matrix'
		return [r[index] for r in matrix]

	#fig,ax = plt.subplots
	plt.figure(pltctr)
	for fx,out_file_matrix in enumerate(file_matrices): #will still need workaround for con_surv.d idea
		plt.scatter(Col(xind,out_file_matrix), Col(yind,out_file_matrix),
			alpha=Alpha[fx], color=Color[fx], s=Size[fx], label=Label[fx], marker='.')
	plt.title(ypar+" v "+xpar)
	plt.ylabel(ypar)
	plt.xlabel(xpar)
	plt.legend(loc=LOC, bbox_to_anchor=BBOX_TO_ANCHOR, ncols=num_files)
	
	
	if xpar in ["lambda","kappa"] or ypar in ["lambda","kappa"]:		# If L or K is involved,
		if xpar in ["lambda","kappa"]: plt.xlim(0,0.8)			#  let first plot be 
		if ypar in ["lambda","kappa"]: plt.ylim(0,0.8)			#  confined to (0,0.8)
		plt.savefig("{}{}_v_{}.png".format(DIR, ypar, xpar), dpi=DPI)	#  for the L&/or K axi/es
	else: plt.savefig("{}{}_v_{}.png".format(DIR, ypar, xpar), dpi=DPI)	# Otherwise just use DEF.
	
	if xmin==xmax and ymin == ymax: plt.close()				# If didn't specify to zoom
	else:									#  just close the plot,
		if xmin!=xmax: plt.xlim(xmin,xmax)				# Otherwise enforce chosen
		if ymin!=ymax: plt.ylim(ymin,ymax)				#  x/y windows
		plt.savefig("{}{}_v_{}_zoom.png".format(DIR, ypar, xpar), dpi=DPI)
		plt.close()
	return


def GeneratePlots():#(out_file_name, file_index, SAVEFIGS):#####NEXT STEP, PUT A LOOP AROUND EACH PLOT SO IT CAN BE DELETED ONCE SAVEFIG'D
	
	Label = [x[:-7] for x in file_names]
	#Color = ['darkgray', 'cyan', 'yellow', 'magenta']		#CYM COLORING
	Color = ['black', 'b', 'g', 'r']			#RGB COLORING
	Alpha = [1,1,1,1]   
	Size = [5,5,2,.5] #CONCENTRIC SIZING
	#Size = [0,6,5,4] #ALPHA-STACK SIZING
	LOC = "center"
	BBOX_TO_ANCHOR = [0.5,1.1125]
	DPI = 480
	#CURRENT CONSIDERATIONS ARE TO PLOT A SINGULAR TIME FROM THE UNCONSTRAINED SET, AND TAG EACH POINT WITH
	# A CERTAIN COLOR BASED ON WHICH CONSTRAINTS IT SURVIVES	
	#color_list = constraints_survived #alias list of colors survived	

# DIMENSIONLESS
	pltctr = 1
	SinglePlot(pltctr, "lambda", 19, 0, 0,
                        "tanB", 1, 0, 0,	
			Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI)			

	pltctr+=1
	SinglePlot(pltctr, "kappa", 20, 0, 0,
			"tanB", 1, 0, 0,
			Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI)
	pltctr+=1
	SinglePlot(pltctr, "lambda", 19, 0, 0,
			"kappa", 20, 0, 0,
			Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI)
# LESS-FUL
	pltctr+=1
	SinglePlot(pltctr, "Alambda", 21, 0, 0,
			"tanB", 1, 0, 0,
			Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI)

	pltctr+=1
	SinglePlot(pltctr, "Akappa", 22, 0, 0,
			"tanB", 1, 0, 0,
			Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI)

	pltctr+=1
	SinglePlot(pltctr, "mueff", 23, 0, 0,
			"tanB", 1, 0, 0,
			Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI)
# DIMENSIONFUL
	pltctr+=1
	SinglePlot(pltctr, "Alambda", 21, 0, 0,
			"Akappa", 22, 0, 0,
			Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI)

	pltctr+=1
	SinglePlot(pltctr, "Alambda", 21, 0, 0,
			"mueff", 23, 0, 0,
			Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI)

	pltctr+=1
	SinglePlot(pltctr, "mueff", 23, 0, 0,
			"Akappa", 22, 0, 0,
			Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI)
# END OF PARAM-PARAM PLOTS
##	> this is just a copy-paste empty template
#	pltctr+=1
#	SinglePlot(pltctr, "", 0, 0, 0,
#			"", 0, 0, 0,
#			Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI)

	print("Beginning mass plots")	
	for (higgs,hix) in [("s1",24),("s2",26),("s3",28),("p1",30),("p2",32)]:
		print("{}:".format(higgs))
		for (param,pix) in [("tanB",1), ("Alambda",21), ("Akappa",22), ("{}comp".format(higgs), hix+1)]:
			pltctr+=1
			if higgs == "s1": (mmin, mmax) = (112.5,132.5)
			elif higgs == "s2": (mmin, mmax) = (0, 0)
			elif higgs == "s3": (mmin, mmax) = (0, 37500)
			elif higgs == "p1": (mmin, mmax) = (0, 10000)
			elif higgs == "p2": (mmin, mmax) = (0, 35000)
			SinglePlot(pltctr,"{}mass".format(higgs), hix, mmin, mmax,
					param, pix, 0, 0,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI)	
	print("Finished plots.")
	return
# FOR MIN,MAX MASS FILTERED BY con2
(mins1m,maxs1m) = (140,100) # declared as such so 1st data point forces min/max to resolve
mh_is1 = [0,0,0,0] # elem corresp to file_index from enumerating file_names 
mh_is2 = [0,0,0,0]
mh_dne = [0,0,0,0]
mh_1n2 = [0,0,0,0]
mh_is3 = [0,0,0,0]

def NearSM(mass): #temp fnc that checks if mass is near sm higgs mass
	buffer = 3 #buffer in GeV around central SM Higgs mass 
	return (mass > 125-buffer) and (mass < 125+buffer)

file_matrices = [list(),list(),list(),list()] # for storing "out_file_matrix" of each out file
constraints_survived = list()# for each element of unconstrained set, note which points survive later constraints

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
				print("{}:event{}:s1@{}:s2@{}\t*** 1&2".format(
					out_file_name,indexrow,shiggs[0],shiggs[2]))
			elif (NearSM(shiggs[0])):
				mh_is1[file_index]+=1
			elif (NearSM(shiggs[2])):
				mh_is2[file_index]+=1
				print("{}event#{}:\ts1@{}\s2@{}\t********** 2".format(out_file_name,indexrow,shiggs[0],shiggs[2]))
			elif (NearSM(shiggs[4])):
				mh_is3[file_index]+=1
				print("{}:event{}:s3@{}\t*************** 3".format(
					out_file_name,indexrow,shiggs[4]))
			else: mh_dne[file_index]+=1


			# ALSO, JUST GET THE MAX AND MIN VALUES ALLOWED BY con2 FOR SM-mh
			if ("con2" in out_file_name):
				(mins1m,maxs1m) = (min(mins1m,shiggs[0]),max(maxs1m,shiggs[0]))

	f.close()




print("So it begins...")
then=time()
for index,event in enumerate(file_matrices[0]):	#each event in the unconstrained set...
	if ( index  % 1000 == 0):
		print(index, "@",time()-then)
	event_r = 0
	event_g = 0
	event_b = 0
	if event in file_matrices[1]: event_b = 200/255
	if event in file_matrices[2]: event_g = 200/255
	if event in file_matrices[3]: event_r = 200/255
	constraints_survived.append([event_r,event_g,event_b])
print("That part took\t", time()-then)





# GeneratePlots() 






print("Sorting by lightest SCALAR")
for file_index,out_file_matrix in enumerate(file_matrices):
	sortedbys1mass = sorted(out_file_matrix, key = lambda x: x[24])
	print(file_names[file_index])
	for event_index,event in enumerate(sortedbys1mass):
		if event_index < 5:
			print("{:.4f}\t".format(event[1]),end="")
			print("{:.6f}\t".format(event[19]),end="")
			print("{:.6f}\t".format(event[20]),end="")
			print("{:.2f}\t".format(event[21]),end="")
			print("{:.2f}\t".format(event[22]),end="")
			print("{:.2f}\t".format(event[23]),end="")
			print(event[24])
	print()
print("\nSorting by lightest PSEUDOSCALAR")
for file_index,out_file_matrix in enumerate(file_matrices):
	sortedbyp1mass = sorted(out_file_matrix, key = lambda x: x[30])
	print(file_names[file_index])
	for event_index,event in enumerate(sortedbyp1mass):
		if event_index < 5:
			print("{:.4f}\t".format(event[1]),end="")
			print("{:.6f}\t".format(event[19]),end="")
			print("{:.6f}\t".format(event[20]),end="")
			print("{:.2f}\t".format(event[21]),end="")
			print("{:.2f}\t".format(event[22]),end="")
			print("{:.2f}\t".format(event[23]),end="")
			print(event[30])
	print()




print("con2 min/max s1m:\t",mins1m,"\t",maxs1m)

print("\nFILE_NAME\tmh_is1\tmh_is2\tmh_dne\tmh_1&2\tmh_is3")
for file_index,out_file_name in enumerate(file_names):
	print("{}\t{}\t{}\t{}\t{}\t{}".format(out_file_name[:-3], 
	mh_is1[file_index],mh_is2[file_index],mh_dne[file_index],mh_1n2[file_index],mh_is3[file_index]))







print("\n# Light Higgs:\t", ctr_lighthiggs)
print("Runtime(s):\t",time()-start_time)
