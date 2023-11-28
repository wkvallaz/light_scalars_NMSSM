import numpy as np
import matplotlib.pyplot as plt
from time import time
import glob
import pandas as pd
import csv 
import os
import sys

DO_PARAM = 1
DO_MASS = 1
DO_COMP = 1
DO_HEAT = 1 
DO_MISC = 1
file_prefix = "widep"
 
save_dir_name = input("Where would you like to store the scrape?: ")
start_time = time()
DIR = "/home/wolf/NMSSMTools_6.0.0/calculations/"
YES = ["y","ye","yea","yes"]
try:   # IF IT DOESN'T ALREADY EXIST, CREATE A DIR WITH NAME OF FILE (minus .dat) TO SAVE IMGS INTO
	os.mkdir("{}{}/".format(DIR, save_dir_name))
except OSError as error:
	print(error)
	overwrite_req = input("Continuing will overwrite existing files, are you sure? (y/ye/yea/yes): ").lower()
	if overwrite_req not in YES: sys.exit("Execution halted.")
SAVEPLOTS = input("Do you want to save plots?").lower() in YES
CMYK = input("CMYK mode?") in YES

file_names = [	"{}randout".format(file_prefix),"{}con1randout".format(file_prefix),
		"{}con3randout".format(file_prefix),"{}con2randout".format(file_prefix) ]
num_files = len(file_names)
threshold_lighthiggs = 20 #GeV

def SinglePlot(pltctr, xpar, xind, xmin, xmax, ypar, yind, ymin, ymax,  
		Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, folder, subfolder): 
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
	if folder != "": DIR += folder+"/"
	try:
		os.mkdir(DIR)
	except OSError as error:
		pass
		#print(error)
	if subfolder != "":
		DIR += subfolder+"/"
		try:
			os.mkdir(DIR)
		except OSError as error:
			pass

	print("Plotting #{} ...".format(pltctr))

	def Col(index,matrix): #array generator statement which pulls column 'index' from  ea row of 'matrix'
		return [r[index] for r in matrix]

	#fig,ax = plt.subplots
	plt.figure(pltctr)		#   vvv used to be file_matrices, changed to master_list for m.ex.sets
	if CMYK: relevant_matrix = master_list
	else: relevant_matrix = file_matrices		
	for fx,out_file_matrix in enumerate(relevant_matrix): # file_matrices for RGB, master_list for CMYK
		if ypar == "rt n 3 k Ak mueff div lambda":
			if fx==0:
				plt.plot([0, max(Col(xind,file_matrices[0]))], [0, max(Col(xind,file_matrices[0]))],
					color="gray", alpha=0.5, linestyle="--", linewidth=0.15, label = "y = x")
			fn_arr = [ (-3*r[20]*r[22]*r[23]/r[19])**.5 for r in out_file_matrix]
			plt.scatter(Col(xind,out_file_matrix), fn_arr,
				alpha=Alpha[fx], color=Color[fx], s=Size[fx], label=Label[fx], 
				marker=',', linewidths=0)
			
		else:
			plt.scatter(Col(xind,out_file_matrix), Col(yind,out_file_matrix),
				alpha=Alpha[fx], color=Color[fx], s=Size[fx], label=Label[fx], 
				marker=',', linewidths=0)
	plt.title(ypar+" v "+xpar)
	plt.ylabel(ypar)
	plt.xlabel(xpar)
	plt.legend(loc=LOC, bbox_to_anchor=BBOX_TO_ANCHOR, ncols=num_files+1, frameon=False)
	
	
	if xpar in ["lambda","kappa"] or ypar in ["lambda","kappa"]:		# If L or K is involved,
		if xpar in ["lambda","kappa"]: plt.xlim(0,0.8)			#  let first plot be 
		if ypar in ["lambda","kappa"]: plt.ylim(0,0.8)			#  confined to (0,0.8)
		plt.savefig("{}{}_v_{}.png".format(DIR, ypar, xpar), dpi=DPI)	#  for the L&/or K axi/es
	else: plt.savefig("{}{}_v_{}.png".format(DIR, ypar, xpar), dpi=DPI)	# Otherwise just use DEF.
	
	if xmin==xmax and ymin==ymax: plt.close()				# If didn't specify to zoom
	else:									#  just close the plot,
		if xmin!=xmax: plt.xlim(xmin,xmax)				# Otherwise enforce chosen
		if ymin!=ymax: plt.ylim(ymin,ymax)				#  x/y windows
		plt.savefig("{}{}_v_{}_zoom.png".format(DIR, ypar, xpar), dpi=DPI)
		plt.close()
	return

def HeatPlot(pltctr, cpar, cind, cmap_n, xpar, xind, xmin, xmax, ypar, yind, ymin, ymax,  # fxn reworked ..
		Size, DPI, folder, subfolder): 					#  .. to plot a 2D with a heat map
	
	DIR = "/home/wolf/NMSSMTools_6.0.0/calculations/"+save_dir_name+"/"
	if folder != "": DIR += folder+"/"
	try:
		os.mkdir(DIR)
	except OSError as error:
		pass
		#print(error)
	if subfolder != "":
		DIR += subfolder+"/"
		try:
			os.mkdir(DIR)
		except OSError as error:
			pass
	
	print("Plotting #{} ...".format(pltctr))
	
	def Col(index,matrix): #array generator statement which pulls column 'index' from  ea row of 'matrix'
		return [r[index] for r in matrix]

	#fig,ax = plt.subplots
	plt.figure(pltctr)		#   vvv used to be file_matrices, changed to master_list for m.ex.sets
	if CMYK: relevant_matrix = master_list[-1] #the fully surviving stuff
	else: relevant_matrix = file_matrices[-1] #just LHC CON		
	plt.scatter(Col(xind,relevant_matrix), Col(yind,relevant_matrix),
		c=Col(cind,relevant_matrix), cmap=cmap_n, s=Size[-1], marker=',', linewidths=0)
	plt.title(ypar+" v "+xpar+" c "+cpar)
	plt.ylabel(ypar)
	plt.xlabel(xpar)
	plt.colorbar(label=cpar)
	#plt.legend(loc=LOC, bbox_to_anchor=BBOX_TO_ANCHOR, ncols=num_files, frameon=False)
	
	
	if xpar in ["lambda","kappa"] or ypar in ["lambda","kappa"]:		# If L or K is involved,
		if xpar in ["lambda","kappa"]: plt.xlim(0,0.8)			#  let first plot be 
		if ypar in ["lambda","kappa"]: plt.ylim(0,0.8)			#  confined to (0,0.8)
		plt.savefig("{}{}_v_{}_c_{}.png".format(DIR, ypar, xpar, cpar), dpi=DPI)	#  for the L&/or K axi/es
	else: plt.savefig("{}{}_v_{}_c_{}.png".format(DIR, ypar, xpar, cpar), dpi=DPI)	# Otherwise just use DEF.
	
	if xmin==xmax and ymin==ymax: plt.close()				# If didn't specify to zoom
	else:									#  just close the plot,
		if xmin!=xmax: plt.xlim(xmin,xmax)				# Otherwise enforce chosen
		if ymin!=ymax: plt.ylim(ymin,ymax)				#  x/y windows
		plt.savefig("{}{}_v_{}_c_{}_zoom.png".format(DIR, ypar, xpar, cpar), dpi=DPI)
		plt.close()
	return

def GeneratePlots(DO_PARAM, DO_MASS, DO_COMP, DO_HEAT, DO_MISC):
	# w2: trigger each plot type w separate flag
	if not CMYK:
		Label = [x[:-7] for x in file_names]
		#Color = ['darkgray', 'cyan', 'yellow', 'magenta']		#CYM COLORING
		Color = ['black', 'b', 'g', 'r']			#RGB COLORING
		Alpha = [1,1,1,1]   
		Size = [2,2,.8,.1] #CONCENTRIC SIZING
		#Size = [0,6,5,4] #ALPHA-STACK SIZING
		LOC = "center"
		BBOX_TO_ANCHOR = [0.5,1.1]
		DPI = 480
	#above was for just working with looping over file_matrices and double/trip/quad counting points on plt
	#below, have 8 nonintersecting sets in master_list
	else: #this is== if CMYK
		Label = ["0:None", "1:LEP", "2:LHC", "3:Flav", "12", "13", "23", "123"]
		Color = ["lightgray", "cyan", "magenta", "yellow", "blue", "green","red", "black"]
		Alpha = [1,1,1,1, 1,1,1,1]
		Size= [.25,.25,.25,.25,.25,.25,.25,.25]
		LOC = "center"
		BBOX_TO_ANCHOR = [0.5,1.1]
		DPI = 480
	pltctr = 0
	par_list = [ ("lambda",19), ("kappa",20), ("Alambda",21), ("mueff",23), ("Akappa",22), ("tanB",1) ] 
	mass_list = [ ("s1mass",24), ("s2mass",26), ("s3mass",28), ("p1mass",30), ("p2mass",32) ]
	comp_list = [ ("s1comp",25), ("s2comp",27), ("s3comp",29), ("p1comp",31), ("p2comp",33) ]
	heatmap_list = [#"viridis", "plasma", 
			"inferno", "magma", #"cividis",
			"brg", "rainbow","jet","turbo"] # viridis, plasma, cividis read poorly

	if DO_PARAM:
		print("Beginning parameter plots")
		for i,(xpar,xind) in enumerate(par_list): # ALL PARAM VS
			for j,(ypar,yind) in enumerate(par_list): #PARAM
				if j<=i: continue
				pltctr+=1
				SinglePlot(pltctr, xpar, xind, 0, 0, ypar, yind, 0, 0,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "Parameter", "")
	if DO_MASS:
		print("Beginning mass plots") # PLOT ea Higgs' mass against each parameter also its singlet comp
		for h,(h_mass,hix) in enumerate(mass_list):
			print("{}:".format(h_mass))
			for (param,pix) in par_list+comp_list: #c_l[h] does higgs v own comp, jus c_l v all comps
				pltctr+=1
				if "s1" in h_mass: (mmin, mmax) = (110.0, 140.0)#LHC window
				elif "s2" in h_mass: (mmin, mmax) = (110, 140) #wide 1250, spec2 800
				elif "s3" in h_mass: (mmin, mmax) = (0, 7500)#wide 37.5k, spec2 7.5k
				elif "p1" in h_mass: (mmin, mmax) = (0, 2000) #wide 10k, spec2 2k
				elif "p2" in h_mass: (mmin, mmax) = (0, 7500)#wide 32k, spec2 7500
				SinglePlot(pltctr,h_mass, hix, mmin, mmax,
						param, pix, 0, 0,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "Mass",h_mass)	
	if DO_COMP:
		print("Beginning composition plots") # PLOT each Higgs' singlet comp against each parameter
		for (h_comp,cix) in comp_list:
			print("{}:".format(h_comp))
			for (param,pix) in par_list:
				pltctr+=1
				SinglePlot(pltctr,param,pix,0,0,
						h_comp, cix, 0,0,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "Comp", h_comp)

	if DO_HEAT:
		print("Beginning heat map plots") #heatmaps for s1 to look @ tanB region underneath main LHC blob
		for n,(c_par,c_ix) in enumerate(par_list+comp_list):# params as heatmap choice
			if c_par != "tanB": 	# heatmaps for s1mass @ lo tanB region, LHC blob (exclude tanB)
				pltctr+=1
				HeatPlot(pltctr, c_par, c_ix, heatmap_list[n%len(heatmap_list)],
						"s1mass", 24, 110, 140,
						"tanB", 1, 0, 0, Size, DPI, "Heatmap","s1mass")
			if c_par != "s1comp":
				pltctr+=1	# s1comp v s1mass, what drives survival (exclude s1comp)
				HeatPlot(pltctr, c_par, c_ix, heatmap_list[n%len(heatmap_list)],
						"s1mass", 24, 0, 50,
						"s1comp", 25, .96, 1, Size, DPI, "Heatmap", "s1mass")	
	if DO_MISC:
		print("Comparing LO p1mamss")
		pltctr+=1
		SinglePlot(pltctr, "p1mass", 30, 0,0,
				"rt n 3 k Ak mueff div lambda", 0, 0,0,
			Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "","")

		print("s2mass v s1mass")
		pltctr+=1
		SinglePlot(pltctr, "s1mass", 24, 100, 130,
				"s2mass", 26, 0,0,
			Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "Mass","")
# <empty copy paste template > 
#	pltctr+=1
#	SinglePlot(pltctr, "", 0, 0, 0,
#			"", 0, 0, 0,
#			Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "", "")
#	pltctr+=1
#	HeatPlot(pltctr, cpar, cind, cmap, xpar, xind, xmin, xmax,
#						ypar, yind, ymin, ymax, Size, DPI, "", "")
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

for file_index,out_file_name in enumerate(file_names):
	with open("{}{}.dat".format(DIR, out_file_name)) as f:
		f_reader = csv.reader(f, delimiter=" ")
		out_file_matrix = file_matrices[file_index]
		ctr_lighthiggs = 0
		print("Reading in\t",out_file_name)
		for indexrow,fullrow in enumerate(f_reader):
			row = [0] # trim out strange spacing ---> this used to be the event number
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
		
			# CONTINUE if don't want to count NearSM higgs events nor con2lims
			continue 
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


if CMYK:
	print("Splitting into mutually exclusive sets...")

	def Set(List): # fn is List to set of tuples conversion
		return set(map(tuple, List))
	def List(Set): # fn is Set to List conversion
		return list(map(list, Set))
	for i,e in enumerate(file_matrices[3]):	#BE WEARY THIS IS OVERWRITING ORIGINAL INFO
		e[-2]=0				# ON THE PROD SIGMA THRU GGF / con3 has nonzero

	bset0 = Set(file_matrices[0]) 		# base unc. set
	bset1 = Set(file_matrices[1])		# base con1 set
	bset3 = Set(file_matrices[2])		# base con3 set
	bset2 = Set(file_matrices[3])		# base con2 set
	
	set132 = bset1 & bset3 & bset2 	#union all cons
	set13 = (bset1 & bset3).difference(set132)	#survd 1 and 3, but not 2
	set12 = (bset1 & bset2).difference(set132)	#  .   1 and 2, but not 3
	set23 = (bset2 & bset3).difference(set132)	# .    2 and 3, but not 1
	
	set1 = bset1.difference(bset2,bset3)		#survd 1, but not 2 and 3
	set3 = bset3.difference(bset1,bset2)		#survd 3, but not 1 and 2
	set2 = bset2.difference(bset3,bset1)		#survd 2, but not 1 and 3
	set0 = bset0.difference(bset1,bset3,bset2)	#survd no constraints

	master_list = [ List(set0), List(set1), List(set2), List(set3),
		List(set12),List(set13),List(set23),List(set132) ]

# args (DO_PARAM, DO_MASS, DO_COMP, DO_HEAT, DO_MISC)
if SAVEPLOTS: GeneratePlots(DO_PARAM, DO_MASS, DO_COMP, DO_HEAT, DO_MISC)

print("Sorting by lightest SCALAR")
print("tanB\tlambda\t\tkappa\t\tAlambda\tAkappa\t\tmueff\tMass")
for file_index,out_file_matrix in enumerate(file_matrices):
	sortedbys1mass = sorted(out_file_matrix, key = lambda x: x[24])
	print(file_names[file_index])
	for event_index,event in enumerate(sortedbys1mass):
		if event_index < 3:
			print("{:.4f}\t".format(event[1]),end="")
			print("{:.6f}\t".format(event[19]),end="")
			print("{:.6f}\t".format(event[20]),end="")
			print("{:.2f}\t".format(event[21]),end="")
			print("{:.7f}\t".format(event[22]),end="")
			print("{:.2f}\t".format(event[23]),end="")
			print(event[24])
	print()
print("\nSorting by lightest PSEUDOSCALAR")
print("tanB\tlambda\t\tkappa\t\tAlambda\tAkappa\t\tmueff\tMass")
for file_index,out_file_matrix in enumerate(file_matrices):
	sortedbyp1mass = sorted(out_file_matrix, key = lambda x: x[30])
	print(file_names[file_index])
	for event_index,event in enumerate(sortedbyp1mass):
		if event_index < 3:
			print("{:.4f}\t".format(event[1]),end="")
			print("{:.6f}\t".format(event[19]),end="")
			print("{:.6f}\t".format(event[20]),end="")
			print("{:.2f}\t".format(event[21]),end="")
			print("{:.7f}\t".format(event[22]),end="")
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
