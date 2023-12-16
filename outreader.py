import numpy as np
import cmath
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

DEBUG_MODE = False #enables print statements used for tracking

DO_PARAM = 1
DO_MASS = 1
DO_COMP = 1
DO_HEAT = 1 
DO_MISC = 1
#=-=# file_prefix = "--"# widep
file_prefix = argv[1]
#file_tags = ["","con1","con3","con2"]
file_tags = ['THY','LEP','LHC','BKF']
file_names = ["{}{}randout".format(file_prefix, tag) for tag in file_tags]

def Time(): #current runtime to 10ms place
	return round(time() - start_time,2)
											# DIR+OW #
#save_dir_name = input("Where would you like to store the scrape?: ")			##########
save_dir_name = argv[2]
start_time = time()
DIR = "/home/wolf/NMSSMTools_6.0.0/calculations/"
YES = ["y","ye","yea","yes", 1, "1", True, "true"]
try:   # IF IT DOESN'T ALREADY EXIST, CREATE A DIR WITH NAME OF FILE (minus .dat) TO SAVE IMGS INTO
	os.mkdir("{}{}/".format(DIR, save_dir_name))
except OSError as error:											# USER INP #
	if DEBUG_MODE: print(error)  		# possibly very stupidly, for scripting, ALWAYS overwrite				############
#	print(error)
#	overwrite_req = input("Continuing will overwrite existing files, are you sure? (y/ye/yea/yes): ").lower()
#	if overwrite_req not in YES: sys.exit("Execution halted.")			##########

#SAVEPLOTS = input("Do you want to save plots?").lower() in YES
#CMYK = input("CMYK mode?") in YES							##########
SAVEPLOTS = (argv[3].lower() in YES)												############
if DEBUG_MODE: print("SAVEPLOTS: ", SAVEPLOTS)
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
		if DEBUG_MODE: print(error)
	if subfolder != "":
		DIR += subfolder+"/"
		try:
			os.mkdir(DIR)
		except OSError as error:
			if DEBUG_MODE: print(error)

	def Col(index,matrix): #array generator statement which pulls column 'index' from  ea row of 'matrix'
		return [r[index] for r in matrix]
	if "comp" in ypar:
		y_expon = 2
	else:
		y_expon = 1
	#fig,ax = plt.subplots
	plt.figure(pltctr)		#   vvv used to be file_matrices, changed to master_list for m.ex.sets
	if CMYK: relevant_matrix = master_list
	else: relevant_matrix = file_matrices	
	extralegelem=0	
	for fx,out_file_matrix in enumerate(relevant_matrix): # file_matrices for RGB, master_list for CMYK
		if (ypar == "rt n 3 k Ak mueff div lambda"  or (xpar == "MA") or 
#		  (xpar == "s1mass" and ypar == "s2mass")  or		
#		  (xpar == "s2mass" and ypar == "s3mass")  or
#		  (xpar == "p1mass" and ypar == "p2mass")  or
		  (xpar[-4:] == "mass" and ypar[-4:] == "mass")	):
			if fx==0 and (xpar != "MA"): #only plot line y = x on the first iter of loop & don't
				extralegelem=1	     #				put y=x for the MA stuff
				plt.plot([0,max(Col(xind,file_matrices[0]))],
					 [0,max(Col(xind,file_matrices[0]))],
				 	 	color="black", alpha=1.0, linestyle="solid", linewidth=0.30,
						 label = "y = x")#lw0.15@dpi480
			# y values described by fn_arr (function array)
			if ypar == "rt n 3 k Ak mueff div lambda":	
				fn_arr = [np.real(cmath.sqrt(
	-3*r[20]*r[22]*r[23]/r[19]						)) for r in out_file_matrix]
			else:
				fn_arr = Col(yind, out_file_matrix)
			# x values descd by fnxarr
			if xpar == "MA":
				fnxarr = [np.real(cmath.sqrt(
	2*r[23]*(r[19]*r[21]+r[20]*r[23])/(r[19]*np.sin(2*np.arctan(r[1])))	)) for r in out_file_matrix]
			else:
				fnxarr = Col(xind,out_file_matrix)
			plt.scatter(fnxarr, fn_arr,
				alpha=Alpha[fx],color=Color[fx],s=Size[fx],
				label=Label[fx],marker=',',linewidths=0)
		else:
			plt.scatter(Col(xind,out_file_matrix), [r[yind]**y_expon for r in out_file_matrix],
				alpha=Alpha[fx], color=Color[fx], s=Size[fx],
				label=Label[fx], marker=',', linewidths=0)
	plt.title(file_prefix+" : "+ypar+" v "+xpar)
	plt.ylabel(ypar)
	
	if (len(Label)+extralegelem) <= 4: Ncols = len(Label)+extralegelem
	else: Ncols = np.ceil( (len(Label)+extralegelem)/2 )
	plt.xlabel(xpar)
	leg = plt.legend(loc=LOC, bbox_to_anchor=BBOX_TO_ANCHOR, ncols=Ncols, columnspacing=0.7, frameon=False)
	for x in range(len(Label)+extralegelem): leg.legend_handles[x]._sizes = [10]
	
	if xpar in ["lambda","kappa"] or ypar in ["lambda","kappa"] or (file_prefix == "13032113" and xpar == "Alambda"):		# If L or K is involved,
		if xpar in ["lambda"]: plt.xlim(0,1)			#  let first plot be 
		elif ypar in ["lambda"]: plt.ylim(0,1)			#  confined  (,)
		if xpar in ["kappa"]: plt.xlim(0,1)			#   
		elif ypar in ["kappa"]: plt.ylim(0,1)			#  
		if (file_prefix == "13032113"): 
			if xpar == "Alambda": plt.xlim(-700,300)
			elif ypar == "Alambda": plt.ylim(-700,300)
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
		if DEBUG_MODE: print(error)
	if subfolder != "":
		DIR += subfolder+"/"
		try:
			os.mkdir(DIR)
		except OSError as error:
			if DEBUG_MODE: print(error)
	
	
	def Col(index,matrix): #array generator statement which pulls column 'index' from  ea row of 'matrix'
		return [r[index] for r in matrix]

	plt.figure(pltctr)		#   vvv used to be file_matrices, changed to master_list for m.ex.sets
	if CMYK: relevant_matrix = master_list[-1] #the fully surviving stuff
	else: relevant_matrix = file_matrices[-1] #just LHC CON
	if "comp" in ypar: y_expon = 2
	else: y_expon = 1
	if "comp" in cpar: c_expon = 2
	else: c_expon = 1		
	plt.scatter([r[xind] for r in relevant_matrix], [r[yind]**y_expon for r in relevant_matrix],
	  c=[r[cind]**c_expon for r in relevant_matrix], cmap=cmap_n, s=Size[-1], marker=',', linewidths=0)
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
		DPI = 240
	else: #this is== if CMYK
		#Label = ["0:None", "1:LEP", "2:LHC", "3:Flav", "12", "13", "23", "123"]
		Label = [ 'T','1','2','3',
			  'T1','T2','T3','12','13','23',
			  'T12','T13','T23','123',
			  'T123' ]
		NC = len(Label)
		Color =[(.8,.8,.8), (.5,1,1), (1,.5,1), (1,1,.5),
			(0,.9,.9), (.9,0,.9), (.9,.9,0),
			(.2,.6,1), (.125,1,.125), (1,.125,.125),
			(0,0,.9), (0,.6,0), (.6,0,0), (.4,.4,.4), 
			(0,0,0)]

		Alpha = [1 for x in range(NC)]
		Size= [.1 for x in range(NC)]#.25 when doing 8 colors
		LOC = "center"
		BBOX_TO_ANCHOR = [0.475,1.105]
		DPI = 240 #for very long time operating at 480 DPI, @1730_14dec23 changed to 240
		######## IF DPI @ 480, SIZE OF 0.04 OK. IF DPI @ 240, DOTS DO NOT RENDER @ THAT SIZE. INC TO 0.1
	pltctr = 0
	par_list = [ ("lambda",19), ("kappa",20), ("Alambda",21), ("mueff",23), ("Akappa",22), ("tanB",1) ] 
	mass_list = [ ("s1mass",24), ("s2mass",28), ("s3mass",32), ("p1mass",36), ("p2mass",39), ("cmass",42) ]
		# after edits to nmhdecay_rand.f these comps are matrix elems not true compositions, ned **2
		# comps also used to just be called comp, but specifically called as singlet comp, change filenm
	comp_list = [ ("s1scomp",27), ("s2scomp",31), ("s3scomp",35), ("p1scomp",38), ("p2scomp",41) ]
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
			print("{}".format(h_mass))
			for (param,pix) in par_list+comp_list: #c_l[h] does higgs v own comp, jus c_l v all comps
				pltctr+=1
				if "s1" in h_mass: (mmin, mmax) = (110.0, 130.0)#LHC window
				else: (mmin, mmax) = (0, 1000)
				SinglePlot(pltctr,h_mass, hix, mmin, mmax,
						param, pix, 0, 0,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "Mass",h_mass)
			for yit,(y_higg,yix) in enumerate(mass_list):	 #for mass v mass plots
				if yit <= h: continue
				if "s1" in h_mass: (xmin, xmax) = (110, 130)
				else: (xmin, xmax) = (0, 1000)
				if "s1" in y_higg: (ymin, ymax) = (110, 130) # should never trigger
				else: (ymin, ymax) = (0, 1000)
				SinglePlot(pltctr, h_mass, hix, xmin, xmax,
						y_higg, yix, ymin, ymax,
						Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "Mass", "")
	if DO_COMP:
		print(Time(),"\tBeginning composition plots") # PLOT each Higgs' singlet comp against each parameter
		for (h_comp,cix) in comp_list:
			print("{}".format(h_comp))
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
			if c_par != "s1scomp":
				pltctr+=1	# s1scomp v s1mass, what drives survival (exclude s1comp)
				HeatPlot(pltctr, c_par, c_ix, heatmap_list[n%len(heatmap_list)],
						"s1mass", 24, 0, 50,
						"s1scomp", 27, .96, 1, Size, DPI, "Heatmap", "s1mass")	
			# possibly kappa v lambda c mueff
			
	if DO_MISC:
		print(Time(),"\tComparing LO p1mamss")
		pltctr+=1
		SinglePlot(pltctr, "p1mass", 36, 0,1000,
				"rt n 3 k Ak mueff div lambda", 0, 0,1000,
			Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "","")
	
		print(Time(),"\ts1(u,d,s)comp v s1mass")
		pltctr+=1
		fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)
		comp_color_scheme = [(.9,0,.9),(0,.85,.85),(.9,.9,0)]
#		for fx,out_file_matrix in enumerate(master_list):	# master_list[-1] <--> out_file_matrix
		for color,(comp,cix) in enumerate([("s1ucomp",25),("s1dcomp",26),("s1scomp",27)]):
			ax.scatter( [r[24] for r in master_list[-1]], [r[cix]**2 for r in master_list[-1]],
				alpha=0.7, color=comp_color_scheme[color], s=1, label=comp, 
				marker=',', linewidths=0)

		plt.title(file_prefix+" : s1(u,d,s)comp v s1mass")
		plt.ylabel("s1(u,d,s)comp")
		plt.xlabel("s1mass")
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=3, columnspacing=0.7, frameon=False)
		for x in range(len(comp_color_scheme)): leg.legend_handles[x]._sizes = [10]
		plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/Mass/s1mass/s1udscomp_v_s1mass.png".format(save_dir_name),dpi=DPI)
		plt.xlim(110,130)
		plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/Mass/s1mass/s1udscomp_v_s1mass_zoom.png".format(save_dir_name),dpi=DPI)
		plt.close()		

		print(Time(),"\ts1(h,H,s)comp v s1mass")
		pltctr+=1
		fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)	

		M2A_list = list()
		M2Z_list = list()
		ahH_list = list()
		s1hsmcomp_list = list() #hsmcomps,
		s1Hbsmcomp_list = list() # and Hbsm comps for each event

		for r in master_list[-1]: #for each event in the set...
			 		  # calculate its alpha value (h-H mixing) and the comps
			M2A = 2*r[23]*(r[19]*r[21]+r[20]*r[23])/(r[19]*np.sin(2*np.arctan(r[1])))
			M2Z = 91.187**2
			ahH = (1/2)*np.arctan(  (2*r[1]/(1-r[1]**2))*( (M2A+M2Z)/(M2A-M2Z) )  )	

			M2A_list.append(M2A)
			M2Z_list.append(M2A)
			ahH_list.append(ahH)

			s1hsmcomp_list.append( r[25]*np.cos(ahH) - r[26]*np.sin(ahH) )
			s1Hbsmcomp_list.append( r[25]*np.sin(ahH) + r[26]*np.cos(ahH) )

		for color,comp in enumerate(["s1scomp", "s1(hsm)comp", "s1(Hbsm)comp"]):
			if comp == "s1(hsm)comp":
				fn_arr = [R**2 for R in s1hsmcomp_list]
			elif comp == "s1(Hbsm)comp":
				fn_arr = [R**2 for R in s1Hbsmcomp_list]
			elif comp == "s1scomp":
				fn_arr = [r[27]**2 for r in master_list[-1]]

			ax.scatter( [r[24] for r in master_list[-1]], fn_arr,
				alpha=.7, color=comp_color_scheme[color-1], s=1, label=comp, 
				marker=',', linewidths=0)
		plt.title(file_prefix+" : s1(h,H,s)comp v s1mass")
		plt.ylabel("s1(h,H,s)comp")
		plt.xlabel("s1mass")
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=3, columnspacing=0.7, frameon=False)
		for x in range(3): leg.legend_handles[x]._sizes = [10]
		plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/Mass/s1mass/s1hHscomp_v_s1mass.png".format(save_dir_name),dpi=DPI)
		plt.xlim(110,130)
		plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/Mass/s1mass/s1hHscomp_v_s1mass_zoom.png".format(save_dir_name),dpi=DPI)
		plt.close()		

		print(Time(),"\tp1(A,s)comp v p1mass")
		pltctr+=1
		fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)	
	
		for color,(comp,cix) in enumerate([("p1Acomp",37),("p1scomp",38)]):
			ax.scatter( [r[36] for r in master_list[-1]], [r[cix]**2 for r in master_list[-1]],
				alpha=.7, color=["blue","red"][color], s=1, label=comp, 
				marker=',', linewidths=0)
		plt.title(file_prefix+" : p1(A,s)comp v p1mass")
		plt.ylabel("p1(A,s)comp")
		plt.xlabel("p1mass")
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=2,columnspacing=0.7, frameon=False)
		for x in range(2): leg.legend_handles[x]._sizes = [10]
		plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/Mass/p1mass/p1Ascomp_v_p1mass.png".format(save_dir_name),dpi=DPI)
		plt.xlim(0,1000)
		plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/Mass/p1mass/p1Ascomp_v_p1mass_zoom.png".format(save_dir_name),dpi=DPI)
		plt.close()		
		
		if file_prefix == "13032113":
			print("13032113-specific plots.")
			pltctr+=1
			SinglePlot(pltctr, "MA", 0, 0, 200,
					"lambda", 19, 0, 0,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "", "")
			pltctr+=1
			SinglePlot(pltctr, "MA", 0, 0, 200,
					"kappa", 20, 0, 0,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "","")



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
			if DEBUG_MODE: 
				if indexrow%100000==0: print(indexrow)
			row = [0] # trim out strange spacing ---> this used to be the event number
			for indexelem,val in enumerate(fullrow):
				if val != "": row.append(float(val))
			out_file_matrix.append(row[0:43])

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
		
			#continue # CONTINUE if don't want to count NearSM higgs events nor con2lims
			
			# tracking which events have NearSM higgs in s1, s2, s1&s2, s3, or none
			if (NearSM(shiggs[0]) and NearSM(shiggs[2])): 
				mh_1n2[file_index]+=1
			elif (NearSM(shiggs[0])):
				mh_is1[file_index]+=1
			elif (NearSM(shiggs[2])):
				mh_is2[file_index]+=1
			elif (NearSM(shiggs[4])):
				mh_is3[file_index]+=1
			else: 
				mh_dne[file_index]+=1
			
	f.close()


if CMYK:
	print(Time(),"\tSplitting into mutually exclusive sets...")

	def Set(List): # fn is List to set of tuples conversion
		return set(map(tuple, List))
	def List(Set): # fn is Set to List conversion
		return list(map(list, Set))
	#for i,e in enumerate(file_matrices[2]):	#BE WEARY THIS IS OVERWRITING ORIGINAL INFO
	#	e[-2]=0				# ON THE PROD SIGMA THRU GGF / con2 has nonzero
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
		print("bsT\t",int(len(bsT)),end="\t|")
		print("bs1\t",int(len(bs1)),end="\t|")
		print("bs2\t",int(len(bs2)),end="\t|")
		print("bs3\t",int(len(bs3)))
	
		sT123 = bsT & bs1 & bs2 & bs3 			# union all constraints
		print("sT123\t",int(len(sT123)),end="\t|\t.\t|\t.\t|\t.\n")

		sT12 = (bsT & bs1 & bs2).difference(sT123) 	# surv.d exactly 3 con.s
		sT13 = (bsT & bs1 & bs3).difference(sT123)
		sT23 = (bsT & bs2 & bs3).difference(sT123)
		s123 = (bs1 & bs2 & bs3).difference(sT123)
		print("sT12\t",len(sT12),end="\t|")
		print("sT13\t",len(sT13),end="\t|")
		print("sT23\t",len(sT23),end="\t|")
		print("s123\t",len(s123))

		sT1 = (bsT & bs1).difference(bs2,bs3)		# surv.d exactly 2 con.s
		sT2 = (bsT & bs2).difference(bs1,bs3)
		sT3 = (bsT & bs3).difference(bs1,bs2)
		s12 = (bs1 & bs2).difference(bsT,bs3)
		s13 = (bs1 & bs3).difference(bsT,bs2)
		s23 = (bs2 & bs3).difference(bsT,bs1)
		print("sT1\t",int(len(sT1)),end="\t|")
		print("sT2\t",int(len(sT2)),end="\t|")
		print("sT3\t",int(len(sT3)),end="\t|\t.")
		print("\ns12\t",int(len(s12)),end="\t|")
		print("s13\t",int(len(s13)),end="\t|")
		print("s23\t",int(len(s23)),end="\t|\t.\n")

		sT = bsT.difference(bs1,bs2,bs3)		# surv.d exactly 1 con.
		s1 = bs1.difference(bsT,bs2,bs3)
		s2 = bs2.difference(bsT,bs1,bs3)
		s3 = bs3.difference(bsT,bs1,bs2)
		print("sT\t",int(len(sT)),end="\t|")
		print("s1\t",int(len(s1)),end="\t|")
		print("s2\t",int(len(s2)),end="\t|")
		print("s3\t",int(len(s3)))
	# DNE events with 0 con.s since all input events are assumed to have exactly 1 con.

		master_list = [	List(sT), List(s1), List(s2), List(s3),
				List(sT1), List(sT2), List(sT3), List(s12), List(s13), List(s23),
				List(sT12), List(sT13), List(sT23), List(s123),
				List(sT123) ]


# args (DO_PARAM, DO_MASS, DO_COMP, DO_HEAT, DO_MISC)
if SAVEPLOTS: 
	if DEBUG_MODE: print(Time(),"\tStarting to plot...")
	GeneratePlots(DO_PARAM, DO_MASS, DO_COMP, DO_HEAT, DO_MISC)

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


