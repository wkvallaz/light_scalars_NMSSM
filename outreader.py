import numpy as np
import matplotlib.pyplot as plt
from time import time
import glob
import pandas as pd
import csv 
import os
import sys

# ARGUMENTS FROM COMMAND LINE [,,,,,], execute this file as: python3 outreader.py file_prefix save_dir_name SAVEPLOTS CMYK
# sys.argv = ['outreader.py', file_prefix, save_dir_name, SAVEPLOTS, CMYK]
## currently have CMYK automatically set as True, not needed in arguments
argv = sys.argv

DO_PARAM = 0
DO_MASS = 0
DO_COMP = 0
DO_HEAT = 1 
DO_MISC = 1
#file_prefix = "--"# widep
file_prefix = argv[1]
#file_tags = ["","con1","con3","con2"]
file_tags = ['THY','LEP','LHC','BKF']
file_names = ["{}{}randout".format(file_prefix, tag) for tag in file_tags]

def Time():
	return time() - start_time
											# DIR+OW #
#save_dir_name = input("Where would you like to store the scrape?: ")			##########
save_dir_name = argv[2]
start_time = time()
DIR = "/home/wolf/NMSSMTools_6.0.0/calculations/"
YES = ["y","ye","yea","yes"]
try:   # IF IT DOESN'T ALREADY EXIST, CREATE A DIR WITH NAME OF FILE (minus .dat) TO SAVE IMGS INTO
	os.mkdir("{}{}/".format(DIR, save_dir_name))
except OSError as error:											# USER INP #
	pass  		# possibly very stupidly, for scripting, ALWAYS overwrite				############
#	print(error)
#	overwrite_req = input("Continuing will overwrite existing files, are you sure? (y/ye/yea/yes): ").lower()
#	if overwrite_req not in YES: sys.exit("Execution halted.")			##########

#SAVEPLOTS = input("Do you want to save plots?").lower() in YES
#CMYK = input("CMYK mode?") in YES							##########
SAVEPLOTS = argv[3]												############
#CMYK = argv[4]													###########
CMYK = True

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

##T	print("Plotting #{} ...".format(pltctr), end="\t")

	def Col(index,matrix): #array generator statement which pulls column 'index' from  ea row of 'matrix'
		return [r[index] for r in matrix]

	#fig,ax = plt.subplots
	plt.figure(pltctr)		#   vvv used to be file_matrices, changed to master_list for m.ex.sets
	if CMYK: relevant_matrix = master_list
	else: relevant_matrix = file_matrices		
	for fx,out_file_matrix in enumerate(relevant_matrix): # file_matrices for RGB, master_list for CMYK
		if ypar == "rt n 3 k Ak mueff div lambda" or (xpar == "s1mass" and ypar == "s2mass"):
			if fx==0: # only plot line y = x on the first iter of loop
				plt.plot([0,max(Col(xind,file_matrices[0]))],[0,max(Col(xind,file_matrices[0]))],
					color="gray", alpha=0.5, linestyle="--", linewidth=0.15, label = "y = x")
			# y values described by fn_arr (function array)
			if ypar == "rt n 3 k Ak mueff div lambda":
				fn_arr = [ np.sqrt(np.sqrt( (-3*r[20]*r[22]*r[23]/r[19])**2)) for r in out_file_matrix]
			else:
				fn_arr = Col(yind, out_file_matrix)
			plt.scatter(Col(xind,out_file_matrix), fn_arr,
				alpha=Alpha[fx],color=Color[fx],s=Size[fx],label=Label[fx],marker=',',linewidths=0)
			
		else:
			plt.scatter(Col(xind,out_file_matrix), Col(yind,out_file_matrix),
				alpha=Alpha[fx], color=Color[fx], s=Size[fx], label=Label[fx], 
				marker=',', linewidths=0)
	plt.title(file_prefix+" : "+ypar+" v "+xpar)
	plt.ylabel(ypar)
	plt.xlabel(xpar)
	leg = plt.legend(loc=LOC, bbox_to_anchor=BBOX_TO_ANCHOR, ncols=8, columnspacing=0.7, frameon=False)
	for x in range(len(Label)): leg.legend_handles[x]._sizes = [10]
	
	if xpar in ["lambda","kappa"] or ypar in ["lambda","kappa"]:		# If L or K is involved,
		if xpar in ["lambda"]: plt.xlim(0,1)			#  let first plot be 
		elif ypar in ["lambda"]: plt.ylim(0,1)			#  confined  (,)
		if xpar in ["kappa"]: plt.xlim(0,1)			#   
		elif ypar in ["kappa"]: plt.ylim(0,1)			#  
		plt.savefig("{}{}_v_{}.png".format(DIR, ypar, xpar), dpi=DPI)	#  for the L&/or K axi/es
	else: plt.savefig("{}{}_v_{}.png".format(DIR, ypar, xpar), dpi=DPI)	# Otherwise just use DEF.
	
	if xmin==xmax and ymin==ymax: plt.close()				# If didn't specify to zoom
	else:									#  just close the plot,
		if xmin!=xmax: plt.xlim(xmin,xmax)				# Otherwise enforce chosen
		if ymin!=ymax: plt.ylim(ymin,ymax)				#  x/y windows
		plt.savefig("{}{}_v_{}_zoom.png".format(DIR, ypar, xpar), dpi=DPI)
		plt.close()
##T	print(Time())
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
	
##T	print("Plotting #{} ...".format(pltctr), end="\t")
	
	def Col(index,matrix): #array generator statement which pulls column 'index' from  ea row of 'matrix'
		return [r[index] for r in matrix]

	plt.figure(pltctr)		#   vvv used to be file_matrices, changed to master_list for m.ex.sets
	if CMYK: relevant_matrix = master_list[-1] #the fully surviving stuff
	else: relevant_matrix = file_matrices[-1] #just LHC CON		
	plt.scatter(Col(xind,relevant_matrix), Col(yind,relevant_matrix),
		c=Col(cind,relevant_matrix), cmap=cmap_n, s=Size[-1], marker=',', linewidths=0)
	plt.title(file_prefix+" : "+ypar+" v "+xpar+" c "+cpar)
	plt.ylabel(ypar)
	plt.xlabel(xpar)
	plt.colorbar(label=cpar)
		
	if xpar in ["lambda","kappa"] or ypar in ["lambda","kappa"]:		# If L or K is involved,
		if xpar in ["lambda","kappa"]: plt.xlim(0,1)			#  let first plot be 
		if ypar in ["lambda","kappa"]: plt.ylim(0,1)			#  confined to (0,0.8)
		plt.savefig("{}{}_v_{}_c_{}.png".format(DIR, ypar, xpar, cpar), dpi=DPI)	#  for the L&/or K axi/es
	else: plt.savefig("{}{}_v_{}_c_{}.png".format(DIR, ypar, xpar, cpar), dpi=DPI)	# Otherwise just use DEF.
	
	if xmin==xmax and ymin==ymax: plt.close()				# If didn't specify to zoom
	else:									#  just close the plot,
		if xmin!=xmax: plt.xlim(xmin,xmax)				# Otherwise enforce chosen
		if ymin!=ymax: plt.ylim(ymin,ymax)				#  x/y windows
		plt.savefig("{}{}_v_{}_c_{}_zoom.png".format(DIR, ypar, xpar, cpar), dpi=DPI)
		plt.close()
##T	print(Time())
	return

def GeneratePlots(DO_PARAM, DO_MASS, DO_COMP, DO_HEAT, DO_MISC):
	if not CMYK:
		Label = [file_prefix + x for x in file_names]
		#Color = ['darkgray', 'cyan', 'yellow', 'magenta']		#CYM COLORING
		Color = ['black', 'b', 'g', 'r']			#RGB COLORING
		Alpha = [1,1,1,1]   
		Size = [5,2,.8,.1] #CONCENTRIC SIZING
		#Size = [0,6,5,4] #ALPHA-STACK SIZING
		LOC = "center"
		BBOX_TO_ANCHOR = [0.475,1.105]
		DPI = 480
	else: #this is== if CMYK
		#Label = ["0:None", "1:LEP", "2:LHC", "3:Flav", "12", "13", "23", "123"]
		Label = [ 'T','1','2','3',
			  'T1','T2','T3','12','13','23',
			  'T12','T13','T23','123',
			  'T123' ]
		NC = len(Label)
		Color =[(.8,.8,.8), (.5,1,1), (1,.5,1), (1,1,.5),
			(0,1,1), (1,0,1), (1,1,0),
			(.2,.6,1), (.125,1,.125), (1,.125,.125),
			(0,0,.9), (0,.6,0), (.6,0,0), (.4,.4,.4), 
			(0,0,0)]

		Alpha = [1 for x in range(NC)]
		Size= [.04 for x in range(NC)]#.25 when doing 8 colors
		LOC = "center"
		BBOX_TO_ANCHOR = [0.475,1.105]
		DPI = 480
	pltctr = 0
	par_list = [ ("lambda",19), ("kappa",20), ("Alambda",21), ("mueff",23), ("Akappa",22), ("tanB",1) ] 
	mass_list = [ ("s1mass",24), ("s2mass",26), ("s3mass",28), ("p1mass",30), ("p2mass",32), ("cmass",34) ]
	comp_list = [ ("s1comp",25), ("s2comp",27), ("s3comp",29), ("p1comp",31), ("p2comp",33) ]
	heatmap_list = [#"viridis", "plasma", 
			"inferno", "magma", #"cividis",
			"brg", "rainbow","jet","turbo"] # viridis, plasma, cividis read poorly

	if DO_PARAM:
		print(Time(),"\tBeginning parameter plots")
		for i,(xpar,xind) in enumerate(par_list): # ALL PARAM VS
			for j,(ypar,yind) in enumerate(par_list): #PARAM
				if j<=i: continue
				pltctr+=1
				SinglePlot(pltctr, xpar, xind, 0, 0, ypar, yind, 0, 0,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "Parameter", "")
	if DO_MASS:
		print(Time(),"\tBeginning mass plots") # PLOT ea Higgs' mass against each parameter also its singlet comp
		for h,(h_mass,hix) in enumerate(mass_list):
			print("{}:".format(h_mass))
			for (param,pix) in par_list+comp_list: #c_l[h] does higgs v own comp, jus c_l v all comps
				pltctr+=1
				if "s1" in h_mass: (mmin, mmax) = (110.0, 130.0)#LHC window
				elif "s2" in h_mass: (mmin, mmax) = (0, 12500) #wide 1250, spec2 800
				elif "s3" in h_mass: (mmin, mmax) = (0, 35000)#wide 37.5k, spec2 7.5k
				elif "p1" in h_mass: (mmin, mmax) = (0, 12500) #wide 10k, spec2 2k
				elif "p2" in h_mass: (mmin, mmax) = (0, 35000)#wide 32k, spec2 7500
				elif "c" in h_mass: (mmin, mmax)=(0,35000)
				SinglePlot(pltctr,h_mass, hix, mmin, mmax,
						param, pix, 0, 0,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "Mass",h_mass)	
	

	if DO_COMP:
		print(Time(),"\tBeginning composition plots") # PLOT each Higgs' singlet comp against each parameter
		for (h_comp,cix) in comp_list:
			print("{}:".format(h_comp))
			for (param,pix) in par_list:
				pltctr+=1
				SinglePlot(pltctr,param,pix,0,0,
						h_comp, cix, 0,0,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "Comp", h_comp)


	if DO_HEAT:
		print(Time(),"\tBeginning heat map plots") #heatmaps for s1 to look @ tanB region underneath main LHC blob
		for n,(c_par,c_ix) in enumerate(par_list+comp_list):# params as heatmap choice
			if c_par != "tanB": 	# heatmaps for s1mass @ lo tanB region, LHC blob (exclude tanB)
				pltctr+=1
				HeatPlot(pltctr, c_par, c_ix, heatmap_list[n%len(heatmap_list)],
						"s1mass", 24, 110, 130,
						"tanB", 1, 0, 0, Size, DPI, "Heatmap","s1mass")
			if c_par != "s1comp":
				pltctr+=1	# s1comp v s1mass, what drives survival (exclude s1comp)
				HeatPlot(pltctr, c_par, c_ix, heatmap_list[n%len(heatmap_list)],
						"s1mass", 24, 0, 50,
						"s1comp", 25, .96, 1, Size, DPI, "Heatmap", "s1mass")	
			# possibly kappa v lambda c mueff
			
	if DO_MISC:
		print(Time(),"\tComparing LO p1mamss")
		pltctr+=1
		SinglePlot(pltctr, "p1mass", 30, 0,12500,
				"rt n 3 k Ak mueff div lambda", 0, 0,12500,
			Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "","")

		print(Time(),"\ts2mass v s1mass")
		pltctr+=1
		SinglePlot(pltctr, "s1mass", 24, 110, 130,
				"s2mass", 26, 110,12500,
			Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "Mass","")
	
		print(Time(),"\ts(1,2,3)comp v s1mass")
		pltctr+=1
		fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)
		comp_color_scheme = ['red','green','blue']
#		for fx,out_file_matrix in enumerate(master_list):	# master_list[-1] <--> out_file_matrix
		for color,(comp,cix) in enumerate(comp_list[0:3]):
			ax.scatter( [r[24] for r in master_list[-1]], [r[cix] for r in master_list[-1]],
				alpha=1, color=comp_color_scheme[color], s=Size[0], label=comp, 
				marker=',', linewidths=0)
		plt.title(file_prefix+" : s(1,2,3)comp v s1mass")
		plt.ylabel("s(1,2,3)comp")
		plt.xlabel("s1mass")
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=3, columnspacing=0.7, frameon=False)
		for x in range(len(comp_color_scheme)): leg.legend_handles[x]._sizes = [10]
		plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/Mass/s1mass/s123comp_v_s1mass.png".format(save_dir_name),dpi=DPI)
		plt.xlim(110,130)
		plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/Mass/s1mass/s123comp_v_s1mass_zoom.png".format(save_dir_name),dpi=DPI)
		plt.close()		



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
mh_is1 = [0,0,0,0] # elem corresp to file_index from enumerating file_names 
mh_is2 = [0,0,0,0]
mh_dne = [0,0,0,0]
mh_1n2 = [0,0,0,0]
mh_is3 = [0,0,0,0]

def NearSM(mass): #temp fnc that checks if mass is near sm higgs mass
	buffer = 3 #buffer in GeV around central SM Higgs mass 
	#buffer = 25 # large buffer for testing new 4constyle
	return (mass > 125-buffer) and (mass < 125+buffer)

file_matrices = [list() for file in file_names] # for storing "out_file_matrix" of each out file

for file_index,out_file_name in enumerate(file_names):
	with open("{}{}.dat".format(DIR, out_file_name)) as f:
		f_reader = csv.reader(f, delimiter=" ")
		out_file_matrix = file_matrices[file_index]
		ctr_lighthiggs = 0
		print("{}\tReading in\t{}".format(Time(), out_file_name))
		for indexrow,fullrow in enumerate(f_reader):
			if indexrow%10000==0: print(indexrow)
			row = [0] # trim out strange spacing ---> this used to be the event number
			for val in fullrow:
				if val != "": row.append(float(val))
			out_file_matrix.append(row)
			
			continue # CONTINUING TO IGNORE COUNTING LIGHT/SMLIKE HIGGS EVENTS
			
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
			
			# tracking which events have NearSM higgs in s1, s2, both, s3, or none
			if (NearSM(shiggs[0]) and NearSM(shiggs[2])): 
				mh_1n2[file_index]+=1
#				print("{}:event{}:s1@{}:s2@{}\t*** 1&2".format(
#					out_file_name,indexrow,shiggs[0],shiggs[2]))
			elif (NearSM(shiggs[0])):
				mh_is1[file_index]+=1
			elif (NearSM(shiggs[2])):
				mh_is2[file_index]+=1
#				print("{}event#{}:\ts1@{}\s2@{}\t********** 2".format(out_file_name,indexrow,shiggs[0],shiggs[2]))
			elif (NearSM(shiggs[4])):
				mh_is3[file_index]+=1
#				print("{}:event{}:s3@{}\t*************** 3".format(
#					out_file_name,indexrow,shiggs[4]))
			else: mh_dne[file_index]+=1
			
	f.close()


if CMYK:
	print(Time(),"\tSplitting into mutually exclusive sets...")

	def Set(List): # fn is List to set of tuples conversion
		return set(map(tuple, List))
	def List(Set): # fn is Set to List conversion
		return list(map(list, Set))
	for i,e in enumerate(file_matrices[2]):	#BE WEARY THIS IS OVERWRITING ORIGINAL INFO
		e[-2]=0				# ON THE PROD SIGMA THRU GGF / con3 has nonzero
			# above arg corresponds to LHC FILE, is [3] for [ "","con1","con3","con2"]

	if False: #leaving this here but copying this architecture for the THY/LEP/LHC/BKF idea...
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
	else:	#below just T/L/L/B idea (as noted in the if False line)
		# file_matrices is [ THYoutmat, LEPoutmat, LHCoutmat, BKFoutmat ]
		# know that len of master_list going to be 2^n with n=num_files=len(file_matrices)
		#  then, later, labels in CMYK mode also going to have to be adjusted if
		#                running with 16 COLORS??? ouch (EDIT: actually 15 since no 0con set)
		
		bsT = Set(file_matrices[0]) 			# base sets exactly 1 con. applied
		bs1 = Set(file_matrices[1])
		bs2 = Set(file_matrices[2])
		bs3 = Set(file_matrices[3])
		print("Len. bsT:",len(bsT)/1)
		print("     bs1:",len(bs1)/1)
		print("     bs2:",len(bs2)/1)
		print("     bs3:",len(bs3)/1)
	
		sT123 = bsT & bs1 & bs2 & bs3 			# union all constraints
		print("   sT123:",len(sT123)/1)

		sT12 = (bsT & bs1 & bs2).difference(sT123) 	# surv.d exactly 3 con.s
		sT13 = (bsT & bs1 & bs3).difference(sT123)
		sT23 = (bsT & bs2 & bs3).difference(sT123)
		s123 = (bs1 & bs2 & bs3).difference(sT123)
		print("    sT12:",len(sT12)/1)
		print("    sT13:",len(sT13)/1)
		print("    sT23:",len(sT23)/1)
		print("    s123:",len(s123)/1)

		sT1 = (bsT & bs1).difference(bs2,bs3)		# surv.d exactly 2 con.s
		sT2 = (bsT & bs2).difference(bs1,bs3)
		sT3 = (bsT & bs3).difference(bs1,bs2)
		s12 = (bs1 & bs2).difference(bsT,bs3)
		s13 = (bs1 & bs3).difference(bsT,bs2)
		s23 = (bs2 & bs3).difference(bsT,bs1)
		print("     sT1:",len(sT1)/1)
		print("     sT2:",len(sT2)/1)
		print("     sT3:",len(sT3)/1)
		print("     s12:",len(s12)/1)
		print("     s13:",len(s13)/1)
		print("     s23:",len(s23)/1)

		sT = bsT.difference(bs1,bs2,bs3)		# surv.d exactly 1 con.
		s1 = bs1.difference(bsT,bs2,bs3)
		s2 = bs2.difference(bsT,bs1,bs3)
		s3 = bs3.difference(bsT,bs1,bs2)
		print("      sT:",len(sT)/1)
		print("      s1:",len(s1)/1)
		print("      s2:",len(s2)/1)
		print("      s3:",len(s3)/1)
	# DNE events with 0 con.s since all input events are assumed to have exactly 1 con.

		master_list = [	List(sT), List(s1), List(s2), List(s3),
				List(sT1), List(sT2), List(sT3), List(s12), List(s13), List(s23),
				List(sT12), List(sT13), List(sT23), List(s123),
				List(sT123) ]


# args (DO_PARAM, DO_MASS, DO_COMP, DO_HEAT, DO_MISC)
if SAVEPLOTS: GeneratePlots(DO_PARAM, DO_MASS, DO_COMP, DO_HEAT, DO_MISC)

print("\n# Light Higgs:\t", ctr_lighthiggs)
print("Runtime(s):\t",Time(),"\n#=#=#=#=#=#=#=#=#=#=#=#=#=#=#")

sys.exit()

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

print("\nFILE_NAME\tmh_is1\tmh_is2\tmh_dne\tmh_1&2\tmh_is3")
for file_index,out_file_name in enumerate(file_names):
	print("{}\t{}\t{}\t{}\t{}\t{}".format(out_file_name[:-3], 
	mh_is1[file_index],mh_is2[file_index],mh_dne[file_index],mh_1n2[file_index],mh_is3[file_index]))


