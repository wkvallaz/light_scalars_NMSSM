import numpy as np
import cmath
import matplotlib.pyplot as plt
from time import time
import glob
import pandas as pd
import csv 
import os
import sys
import gc

# ARGS FROM COMMAND LINE [,,,,,], execute this file as: python3 outreader.py file_prefix save_dir_name SAVEPLOTS
# sys.argv = ['outreader.py', file_prefix, save_dir_name, SAVEPLOTS]
## currently have CMYK automatically set as True, not needed in arguments
argv = sys.argv

DEBUG_MODE = 0 #enables print statements used for tracking
MASSTRK = 0 #enables tracking masses near LHC and of light s/o
DO_PARAM = 1
DO_MASS = 1
DO_COMP = 0
DO_HEAT = 0 
DO_MISC = 0
DO_REPL = 0

threshold_lighthiggs = 50 # GeV
file_prefix = argv[1] #=-=# file_prefix = "--"# widep
file_tags = ['THY','LEP','LHC','BKF']#file_tags = ["","con1","con3","con2"]
file_names = ["{}{}randout".format(file_prefix, tag) for tag in file_tags]

N_EXTRA = 0 # number of extra seeds for files (x4 for actual num extra files)
for ie in range(N_EXTRA):
	extra_names = ["{}_{}{}randout".format(file_prefix, ie+2, tag) for tag in file_tags]
	file_names = file_names + extra_names

if "108035020" in file_prefix: (KMIN, KMAX, LMIN, LMAX) = (-.015, .015, 0, .1)	#def plot axis window
elif "s" == file_prefix[0]: (KMIN, KMAX, LMIN, LMAX) = (0,.1,0,.5)
elif "PQ" == file_prefix[0:2]: (KMIN, KMAX, LMIN, LMAX) = (0,.1,0,.7)
else: (KMIN, KMAX, LMIN, LMAX) = (0, 1, 0, 1)

save_dir_name = argv[2]
start_time = time()
DIR = "/home/wolf/NMSSMTools_6.0.0/calculations/"
YES = ["y","ye","yea","yes", 1, "1", True, "true"]
try:   
	os.mkdir("{}{}/".format(DIR, save_dir_name))
except OSError as error:							
	if DEBUG_MODE: print(error)  		

SAVEPLOTS = (argv[3].lower() in YES)
if DEBUG_MODE: print("SAVEPLOTS: ", SAVEPLOTS)
CMYK = True

num_files = len(file_names)

def Time(): #current runtime to 10ms place
	return round(time() - start_time,2)

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
			if xpar == "xnopyt": #THIS USED TO BE if xpar == "MA" but it innaccurate for MA value
				fnxarr = [np.real(cmath.sqrt(
	2*r[23]*(r[19]*r[21]+r[20]*r[23])/(r[19]*np.sin(2*np.arctan(r[1])))	)) for r in out_file_matrix]
			else:
				fnxarr = Col(xind,out_file_matrix)
			plt.scatter(fnxarr, fn_arr,
				alpha=Alpha[fx],color=Color[fx],s=Size[fx],
				label=Label[fx],marker=',',linewidths=0)
		elif ypar=="neu3Hcomp":
			plt.scatter([r[xind] for r in out_file_matrix], [r[59]**2+r[60]**2 for r in out_file_matrix], alpha=Alpha[fx], color=Color[fx], s=Size[fx], label=Label[fx], marker=',', linewidths=0)
		elif ypar=="neu4Hcomp":
			plt.scatter([r[xind] for r in out_file_matrix], [r[65]**2+r[66]**2 for r in out_file_matrix], alpha=Alpha[fx], color=Color[fx], s=Size[fx], label=Label[fx], marker=',', linewidths=0)
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
		if xpar in ["lambda"]: plt.xlim(LMIN,LMAX)			#  let first plot be 
		elif ypar in ["lambda"]: plt.ylim(LMIN,LMAX)			#  confined  (,)
		if xpar in ["kappa"]: plt.xlim(KMIN,KMAX)			#   
		elif ypar in ["kappa"]: plt.ylim(KMIN,KMAX)			#  
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
		Size, DPI, folder, subfolder): 				#  .. to plot a 2D with a heat map
	
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
	
	if "13032113" in file_prefix: 
		LEN = len(master_list)
		relevant_matrix = master_list[LEN-6]+master_list[LEN-5]+master_list[LEN-4]+master_list[LEN-3]+master_list[LEN-2]+master_list[LEN-1]

	

	if "comp" in ypar: y_expon = 2
	else: y_expon = 1
	if "comp" in cpar: c_expon = 2
	else: c_expon = 1		
	if "108035020" in file_prefix and "k div l" == xpar:
		plt.scatter([r[20]/r[19] for r in relevant_matrix],[r[yind]**y_expon for r in relevant_matrix],			c=[r[cind]**c_expon for r in relevant_matrix], cmap=cmap_n, s=Size[-1],marker=',',linewidths=0)
	else:
		plt.scatter([r[xind] for r in relevant_matrix], [r[yind]**y_expon for r in relevant_matrix],		 	c=[r[cind]**c_expon for r in relevant_matrix], cmap=cmap_n, s=Size[-1],marker=',',linewidths=0)
	plt.title(file_prefix+" : "+ypar+" v "+xpar+" c "+cpar)
	plt.ylabel(ypar)
	plt.xlabel(xpar)
	plt.colorbar(label=cpar)
		
	if xpar in ["lambda","kappa"] or ypar in ["lambda","kappa"]:		# If L or K is involved,
		if xpar in ["lambda"]: plt.xlim(LMIN,LMAX)
		elif ypar == "lambda": plt.ylim(LMIN,LMAX)			#  let first plot be 
		if xpar in ["kappa"]: plt.xlim(KMIN,KMAX)
		elif ypar == "kappa": plt.ylim(KMIN,KMAX)			#  confined to (0,0.8)
		plt.savefig("{}{}_v_{}_c_{}.png".format(DIR, ypar, xpar, cpar), dpi=DPI)	#  for the L&/or K axi/es
	else: plt.savefig("{}{}_v_{}_c_{}.png".format(DIR, ypar, xpar, cpar), dpi=DPI)	# Otherwise just use DEF.
	
	if xmin==xmax and ymin==ymax: plt.close()				# If didn't specify to zoom
	else:									#  just close the plot,
		if xmin!=xmax: plt.xlim(xmin,xmax)				# Otherwise enforce chosen
		if ymin!=ymax: plt.ylim(ymin,ymax)				#  x/y windows
		plt.savefig("{}{}_v_{}_c_{}_zoom.png".format(DIR, ypar, xpar, cpar), dpi=DPI)
		plt.close()
	return

def GeneratePlots(DO_PARAM, DO_MASS, DO_COMP, DO_HEAT, DO_MISC, DO_REPL):
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
		dot_size = 0.2;	Size= [dot_size for x in range(NC-1)]+[8*dot_size]#.25 when doing 8 colors
		
		if "13032113" in file_prefix: Size = [10*s for s in Size]
		elif "108035020" in file_prefix: Size = [10*s for s in Size] 
		LOC = "center"
		BBOX_TO_ANCHOR = [0.475,1.105]
		DPI = 360 #for very long time operating at 480 DPI, @1730_14dec23 changed to 240
		######## IF DPI @ 480, SIZE OF 0.04 OK. IF DPI @ 240, DOTS DO NOT RENDER @ THAT SIZE. INC TO 0.1
	pltctr = 0
	par_list = [ ("lambda",19), ("kappa",20), ("Alambda",21), ("mueff",23), ("Akappa",22), ("tanB",1) ]
	par_list = [("MA",43)] + par_list
	if "108035020" in file_prefix: par_list = [("AU3",5),("M1",2),("M2",3),("M3",4)] + par_list
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
			print("{}\t{}".format(Time(),h_mass))
			for (param,pix) in par_list+comp_list: #c_l[h] does higgs v own comp,
				pltctr+=1			# , just c_l is h v all comps
				if "s1" in h_mass: (mmin, mmax) = (110.0, 130.0)#LHC window
				else: (mmin, mmax) = (0, 500)
				SinglePlot(pltctr,h_mass, hix, mmin, mmax,
						param, pix, 0, 0,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "Mass",h_mass)
			for yit,(y_higg,yix) in enumerate(mass_list):	 #for mass v mass plots
				if yit <= h: continue
				if "s1" in h_mass: (xmin, xmax) = (110, 130)
				else: (xmin, xmax) = (0, 500)
				if "s1" in y_higg: (ymin, ymax) = (110, 130) # should never trigger
				else: (ymin, ymax) = (0, 500)
				SinglePlot(pltctr, h_mass, hix, xmin, xmax,
						y_higg, yix, ymin, ymax,
						Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "Mass", "")
	if DO_COMP:
		print(Time(),"\tBeginning composition plots") # PLOT each Higgs' singlet comp against each parameter
		for (h_comp,cix) in comp_list:
			print("{}\t{}".format(Time(),h_comp))
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
			if c_par != "kappa" and c_par != "lambda":
				pltctr+=1
				HeatPlot(pltctr, c_par, c_ix, heatmap_list[n%len(heatmap_list)],
						"lambda", 19, 0, 0,
						"kappa", 20, 0, 0, Size, DPI, "Heatmap", "")
			if c_par != "kappa" and c_par != "Akappa":
				pltctr+=1
				HeatPlot(pltctr, c_par, c_ix, heatmap_list[n%len(heatmap_list)],
						"kappa", 20, 0, 0,
						"Akappa", 22, 0, 0, Size, DPI, "Heatmap", "108035020")

	if DO_MISC:
		print(Time(),"\tComparing LO p1mamss")
		pltctr+=1
		SinglePlot(pltctr, "p1mass", 36, 0,1000,
				"rt n 3 k Ak mueff div lambda", 0, 0,1000,
			Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "","")

		ahH_list=list()		  # calculate its alpha value (h-H mixing)
		for r in master_list[-1]: #for each event in the set...
			M2A = 2*r[23]*(r[19]*r[21]+r[20]*r[23])/(r[19]*np.sin(2*np.arctan(r[1]))) #DOESN'T WORK
			M2A = r[43]**2
			M2Z = 91.187**2
			ahH = (1/2)*np.arctan(  (2*r[1]/(1-r[1]**2))*( (M2A+M2Z)/(M2A-M2Z) )  )	
			ahH_list.append(ahH)
		
		for (smass, sindex) in [("s1mass",24),("s2mass",28),("s3mass",32)]:
			print(Time(),"\t{}(u,d,s)comp v {}".format(smass[0:2],smass))
			pltctr+=1
			fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)
			comp_color_scheme = [(.9,0,.9),(0,.85,.85),(.9,.9,0)]
			for color,(comp,cix) in enumerate([("{}ucomp".format(smass[0:2]),sindex+1),("{}dcomp".format(smass[0:2]),sindex+2),("{}scomp".format(smass[0:2]),sindex+3)]):
				ax.scatter( [r[sindex] for r in master_list[-1]], [r[cix]**2 for r in master_list[-1]],
					alpha=0.5, color=comp_color_scheme[color], s=1, label=comp, 
					marker=',', linewidths=0)
			plt.title(file_prefix+" : {}(u,d,s)comp v {}".format(smass[0:2],smass))
			plt.ylabel("{}(u,d,s)comp".format(smass[0:2]))
			plt.xlabel(smass)
			leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=3, columnspacing=0.7, frameon=False)
			for x in range(len(comp_color_scheme)): leg.legend_handles[x]._sizes = [10]
			plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/Mass/{}/{}udscomp_v_{}.png".format(save_dir_name,smass,smass[0:2],smass),dpi=DPI)

			if smass=="s1mass":plt.xlim(110,130)
			elif smass=="s2mass":plt.xlim(100,600)
			elif smass=="s3mass":plt.xlim(0,1000)
			plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/Mass/{}/{}udscomp_v_{}_zoom.png".format(save_dir_name,smass,smass[0:2],smass),dpi=DPI)
			plt.close()		

			print(Time(),"\t{}(h,H,s)comp v {}".format(smass[0:2],smass))
			pltctr+=1
			fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)	
	
			shsmcomp_list = list() #hsmcomps,
			sHbsmcomp_list = list() # and Hbsm comps for each event
			for i,r in enumerate(master_list[-1]):
				shsmcomp_list.append( r[sindex+1]*np.cos(ahH_list[i]) - r[sindex+2]*np.sin(ahH_list[i]) )
				sHbsmcomp_list.append( r[sindex+1]*np.sin(ahH_list[i]) + r[sindex+2]*np.cos(ahH_list[i]) )
	
			
			for color,comp in enumerate(["{}scomp".format(smass[0:2]), "{}(hsm)comp".format(smass[0:2]), "{}(Hbsm)comp".format(smass[0:2])]):
				if comp == "{}(hsm)comp".format(smass[0:2]):
					fn_arr = [R**2 for R in shsmcomp_list]
				elif comp == "{}(Hbsm)comp".format(smass[0:2]):
					fn_arr = [R**2 for R in sHbsmcomp_list]
				elif comp == "{}scomp".format(smass[0:2]):
					fn_arr = [r[sindex+3]**2 for r in master_list[-1]]

				ax.scatter( [r[sindex] for r in master_list[-1]], fn_arr,
					alpha=.7, color=comp_color_scheme[color-1], s=1, label=comp, 
					marker=',', linewidths=0)
			plt.title(file_prefix+" : {}(h,H,s)comp v {}".format(smass[0:2],smass))
			plt.ylabel("{}(h,H,s)comp".format(smass[0:2]))
			plt.xlabel(smass)
			leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=3, columnspacing=0.7, frameon=False)
			for x in range(3): leg.legend_handles[x]._sizes = [10]
			plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/Mass/{}/{}hHscomp_v_{}.png".format(save_dir_name,smass,smass[0:2],smass),dpi=DPI)
			if smass=="s1mass":plt.xlim(110,130)
			elif smass=="s2mass":plt.xlim(100,600)
			elif smass=="s3mass":plt.xlim(0,1000)
			plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/Mass/{}/{}hHscomp_v_{}_zoom.png".format(save_dir_name,smass,smass[0:2],smass),dpi=DPI)
			plt.close()		
# END OF SUDSCOMP LOOP
# DONT BOTHER TO LOOP P1 AND P2 COMPS, ITS JUST TWO PLOTS AND EASIER THIS WAY
		print(Time(),"\tp1(A,s)comp v p1mass")
		pltctr+=1
		fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)	
	
		for color,(comp,cix) in enumerate([("p1Acomp",37),("p1scomp",38)]):
			ax.scatter( [r[36] for r in master_list[-1]], [r[cix]**2 for r in master_list[-1]],
				alpha=.5, color=["blue","red"][color], s=1, label=comp, 
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
	
		print(Time(),"\tp2(A,s)comp v p2mass")
		pltctr+=1
		fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)	
	
		for color,(comp,cix) in enumerate([("p2Acomp",40),("p2scomp",41)]):
			ax.scatter( [r[39] for r in master_list[-1]], [r[cix]**2 for r in master_list[-1]],
				alpha=.5, color=["blue","red"][color], s=1, label=comp, 
				marker=',', linewidths=0)
		plt.title(file_prefix+" : p2(A,s)comp v p2mass")
		plt.ylabel("p2(A,s)comp")
		plt.xlabel("p2mass")
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=2,columnspacing=0.7, frameon=False)
		for x in range(2): leg.legend_handles[x]._sizes = [10]
		plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/Mass/p2mass/p2Ascomp_v_p2mass.png".format(save_dir_name),dpi=DPI)
		plt.xlim(0,1000)
		plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/Mass/p2mass/p2Ascomp_v_p2mass_zoom.png".format(save_dir_name),dpi=DPI)
		plt.close()		

### START OF S123U , S123D, S123S COMPS LOOP, possibly sep flag
		for scalarmass,scalarindex in [("s1mass",24),("s2mass",28),("s3mass",32)]:#x axis, then full plt
			for sindx,(jcomp, cstep) in enumerate([("ucomp",1),("dcomp",2),("scomp",3)]): #fore jcomp
				print(Time(),"\ts123{} v {}".format(jcomp,scalarmass))
				pltctr+=1
				fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)
				comp_color_scheme = [(.9,0,.9),(0,.85,.85),(.9,.9,0)]
				for color,(sname,snum) in enumerate([("s1",24),("s2",28),("s3",32)]):
					ax.scatter( [r[scalarindex] for r in master_list[-1]], # xaxis mass
						[r[snum+cstep]**2 for r in master_list[-1]],  #all the jcomps
						alpha=0.5, color=comp_color_scheme[color], s=1,
						label="{}{}".format(sname,jcomp), marker=',', linewidths=0)
				plt.title(file_prefix+" : s123{} v {}".format(jcomp,scalarmass))
				plt.ylabel("s123{}".format(jcomp))
				plt.xlabel(scalarmass)
				leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=3, columnspacing=0.7, frameon=False)
				for x in range(len(comp_color_scheme)): leg.legend_handles[x]._sizes = [10]
				
				try:
					os.mkdir("/home/wolf/NMSSMTools_6.0.0/calculations/{}/Comp/s123/".format(save_dir_name))
				except OSError as error:
					if DEBUG_MODE: print(error)
				
				plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/Comp/s123/s123{}_v_{}.png".format(save_dir_name,jcomp,scalarmass),dpi=DPI)
	
				if scalarmass=="s1mass":plt.xlim(110,130)
				elif scalarmass=="s2mass":plt.xlim(100,600)
				elif scalarmass=="s3mass":plt.xlim(0,1000)
				plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/Comp/s123/s123{}_v_{}_zoom.png".format(save_dir_name,jcomp,scalarmass),dpi=DPI)
				plt.close()		
# S123jCOMP V SiMASS LOOP

	if DO_REPL: # for the very specific plots in replic trials
		if "13032113" in file_prefix:
			print(Time(),"\tPlotting Figure (3)...")
			# 3A cmass v s1mass,s2mass,s3mass in s1m126
			pltctr+=1
			fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)	
			for matrix in master_list:
				for r in matrix:
					if 124<r[24] and r[24]<128:
						ax.scatter( r[42], r[24], alpha=.9, color="red", 
							    s=9, label="s1mass", marker=',', linewidths=0)
						ax.scatter( r[42], r[28], alpha=.9, color="green", 
							    s=9, label="s2mass", marker=',', linewidths=0)
						ax.scatter( r[42], r[32], alpha=.9, color="blue", 
							    s=9, label="s3mass", marker=',', linewidths=0)
						# i just realized this legend will suck
			plt.title(file_prefix+" : (3a) s123mass v cmass in s1m126")	
			plt.ylabel("s(1,2,3)mass")
			plt.xlabel("cmass")
			plt.xlim(100,300)
			plt.ylim(0,300)
			plt.savefig("{}{}/s123mass_v_cmass_s1m126.png".format(DIR, save_dir_name))
			plt.close()
			# 3C cmass v s1m,s2m,s3m in s2m126
			pltctr+=1
			fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)	
			for matrix in master_list:
				for r in matrix: #s2m126 case acc to table 3
					if 124<r[28] and r[28]<128:
						ax.scatter( r[42], r[24], alpha=.9, color="red", 
							    s=9, label="s1mass", marker=',', linewidths=0)
						ax.scatter( r[42], r[28], alpha=.9, color="green", 
							    s=9, label="s2mass", marker=',', linewidths=0)
						ax.scatter( r[42], r[32], alpha=.9, color="blue", 
							    s=9, label="s3mass", marker=',', linewidths=0)
						# i just realized this legend will suck
			plt.title(file_prefix+" : (3c) s123mass v cmass in s2m126")	
			plt.ylabel("s(1,2,3)mass")
			plt.xlabel("cmass")
			plt.xlim(100,300)
			plt.ylim(0,300)
			plt.savefig("{}{}/s123mass_v_cmass_s2m126.png".format(DIR, save_dir_name))
			plt.close()
			# 3E cmass v s1m,s2m,s3m in s3m126
			# 3B cm v p1p2cm
			pltctr+=1
			fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)	
			for matrix in master_list:
				for r in matrix:
					if 124<r[24] and r[24]<128:
						ax.scatter( r[42], r[36], alpha=.9, color="magenta", 
							    s=9, label="p1mass", marker=',', linewidths=0)
						ax.scatter( r[42], r[39], alpha=.9, color="orange", 
							    s=9, label="p2mass", marker=',', linewidths=0)
						ax.scatter( r[42], r[42], alpha=.9, color="cyan", 
							    s=9, label="cmass", marker=',', linewidths=0)
						# i just realized this legend will suck
			plt.title(file_prefix+" : (3b) p1p2cmass v cmass in s1m126")	
			plt.ylabel("p1p2cmass")
			plt.xlabel("cmass")
			plt.xlim(100,300)
			plt.ylim(0,300)
			plt.savefig("{}{}/p1p2cmass_v_cmass_s1m126.png".format(DIR, save_dir_name))
			plt.close()
			# 3D cm v p1p2cm
			pltctr+=1
			fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)	
			for matrix in master_list:
				for r in matrix:
					if 124<r[28] and r[28]<128:
						ax.scatter( r[42], r[36], alpha=.9, color="magenta", 
							    s=9, label="p1mass", marker=',', linewidths=0)
						ax.scatter( r[42], r[39], alpha=.9, color="orange", 
							    s=9, label="p2mass", marker=',', linewidths=0)
						ax.scatter( r[42], r[42], alpha=.9, color="cyan", 
							    s=9, label="cmass", marker=',', linewidths=0)
						# i just realized this legend will suck
			plt.title(file_prefix+" : (3d) p1p2cmass v cmass in s2m126")	
			plt.ylabel("p1p2cmass")
			plt.xlabel("cmass")
			plt.xlim(100,300)
			plt.ylim(0,300)
			plt.savefig("{}{}/p1p2cmass_v_cmass_s2m126.png".format(DIR, save_dir_name))
			plt.close()
	#		# 3F cm v p1p2cm
	#	#	#	#	#	#	#	#	#	#	#	#	#
			print(Time(),"\tPlotting Figure (4)...")
			# 4A param plots in s1m126			
			pltctr+=1
			fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)	
			matcolors = ["gray" for arr in master_list] 	# other 	## hopefully removing
			matcolors[7] = "pink"				# s12
			matcolors[13] = "lightgreen"				# s123
			matcolors[14] = "black"				# sT123
			matcolors[7]="gray"		# changed order of these two bc they plotted in
			matcolors[12]="pink"		#  a different order so 7 doesnt get covered
			for matind,matrix in enumerate(master_list):
				if matind==0: continue	# don't plot s1
				elif matind==7: matind=12
				elif matind==12: matind==7
				for r in matrix: # s1m126 window acc to table 2
					if 100<r[24] and r[24]<140:	#cull st only pts which on plot are pltd
						ptcolor = "lightgray"
						if 124<r[24] and r[24]<128:
							ptcolor = "pink"
							if matind == 13: ptcolor = "lightgreen"
						if matind == 14: ptcolor = "black"
						# i don't have ggF sigma due to processing
						ax.scatter( r[19], r[24], alpha=1., color=ptcolor, 
							    s=30, marker='.', linewidths=0)
			plt.title(file_prefix+" : (4a) s1mass v lambda in s1m126")
			plt.ylabel("s1mass")
			plt.xlabel("lambda")
			plt.xlim(0,1)
			plt.ylim(100,140)
			plt.savefig("{}{}/s1mass_v_lambda_s1m126.png".format(DIR, save_dir_name))
			plt.close()
			# 4B param plots in s1m126
			pltctr+=1
			fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)	
			for matind,matrix in enumerate(master_list):
				if matind==0: continue	# don't plot s1
				elif matind==7: matind=12
				elif matind==12: matind==7
				for r in matrix: # s1m126 window acc to table 2
					if 100<r[24] and r[24]<140:	#cull st only pts which on plot are pltd
						ptcolor = "lightgray"
						if 124<r[24] and r[24]<128:
							ptcolor = "pink"
							if matind == 13: ptcolor = "lightgreen"
						if matind == 14: ptcolor = "black"
						# i don't have ggF sigma due to processing
						ax.scatter( r[20], r[24], alpha=1., color=ptcolor, 
							    s=30, label="s1mass", marker='.', linewidths=0)
			plt.title(file_prefix+" : (4b) s1mass v kappa in s1m126")
			plt.ylabel("s1mass")
			plt.xlabel("kappa")
			plt.xlim(0,1)
			plt.ylim(100,140)
			plt.savefig("{}{}/s1mass_v_kappa_s1m126.png".format(DIR, save_dir_name))
			plt.close()
			# 4C param plots in s1m126
			pltctr+=1
			fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)
			for matind,matrix in enumerate(master_list):
				if matind==0: continue	# don't plot s1
				elif matind==7: matind=12
				elif matind==12: matind==7
				for r in matrix: # s1m126 window acc to table 2
					if 100<r[24] and r[24]<140:	#cull st only pts which on plot are pltd
						ptcolor = "lightgray"
						if 124<r[24] and r[24]<128:
							ptcolor = "pink"
							if matind == 13: ptcolor = "lightgreen"
						if matind == 14: ptcolor = "black"
						# i don't have ggF sigma due to processing
						ax.scatter( r[1], r[24], alpha=1., color=ptcolor, 
							    s=30, label="s1mass", marker='.', linewidths=0)
			plt.title(file_prefix+" : (4c) s1mass v tanB in s1m126")
			plt.ylabel("s1mass")
			plt.xlabel("tanB")
			plt.xlim(0,10)
			plt.ylim(100,140)
			plt.savefig("{}{}/s1mass_v_tanB_s1m126.png".format(DIR, save_dir_name))
			plt.close()
			# 4D param plots in s1m126
			pltctr+=1
			fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)
			for matind,matrix in enumerate(master_list):
				if matind==0: continue	# don't plot s1
				elif matind==7: matind=12
				elif matind==12: matind==7
				for r in matrix: # s1m126 window acc to table 2
					if 100<r[24] and r[24]<140:	#cull st only pts which on plot are pltd
						ptcolor = "lightgray"
						if 124<r[24] and r[24]<128:
							ptcolor = "pink"
							if matind == 13: ptcolor = "lightgreen"
						if matind == 14: ptcolor = "black"
						# i don't have ggF sigma due to processing
						ax.scatter( r[23], r[24], alpha=1., color=ptcolor, 
							    s=30, label="s1mass", marker='.', linewidths=0)
			plt.title(file_prefix+" : (4d) s1mass v mueff in s1m126")
			plt.ylabel("s1mass")
			plt.xlabel("mueff")
			plt.xlim(0,1000)
			plt.ylim(100,140)
			plt.savefig("{}{}/s1mass_v_mueff_s1m126.png".format(DIR, save_dir_name))
			plt.close()
	#	#	#	#	#	#	#	#	#	#	#	#	#
			print(Time(),"\tPlotting Figure (5)...")
			# 5A param plots in s1m126
			pltctr+=1
			fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)
			for matind,matrix in enumerate(master_list):
				if matind==0: continue	# don't plot s1
				elif matind==7: matind=12
				elif matind==12: matind==7
				for r in matrix: # s1m126 window acc to table 2
					ptcolor = "lightgray"
					if 124<r[24] and r[24]<128:
						ptcolor = "pink"
						if r[36]>r[24]/2 and True: ptcolor = "lightgreen"
						if r[36]>r[24]/2 and False: ptcolor = "red"
						if r[36]<r[24]/2: ptcolor = "magenta"
					if matind == 14: ptcolor = "black"
					# i don't have ggF sigma due to processing
					ax.scatter( r[43], r[19], alpha=1., color=ptcolor, 
						    s=30, label="", marker='.', linewidths=0)
			plt.title(file_prefix+" : (5a) lambda v MA in s1m126")
			plt.ylabel("lambda")
			plt.xlabel("MA")
			plt.xlim(0,200)
			plt.ylim(0,1)
			plt.savefig("{}{}/lambda_v_MA_s1m126.png".format(DIR, save_dir_name))
			plt.close()
			# 5B param plots in s1m126
			pltctr+=1
			fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)
			for matind,matrix in enumerate(master_list):
				if matind==0: continue	# don't plot s1
				elif matind==7: matind=12
				elif matind==12: matind==7
				for r in matrix: # s1m126 window acc to table 2
					ptcolor = "lightgray"
					if 124<r[24] and r[24]<128:
						ptcolor = "pink"
						if r[36]>r[24]/2 and True: ptcolor = "lightgreen"
						if r[36]>r[24]/2 and False: ptcolor = "red"
						if r[36]<r[24]/2: ptcolor = "magenta"
					if matind == 14: ptcolor = "black"
					# i don't have ggF sigma due to processing
					ax.scatter( r[43], r[20], alpha=1., color=ptcolor, 
						    s=30, label="", marker='.', linewidths=0)
			plt.title(file_prefix+" : (5b) kappa v MA in s1m126")
			plt.ylabel("kappa")
			plt.xlabel("MA")
			plt.xlim(0,200)
			plt.ylim(0,1)
			plt.savefig("{}{}/kappa_v_MA_s1m126.png".format(DIR, save_dir_name))
			plt.close()
			# 5C param plots in s1m126
			pltctr+=1
			fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)
			for matind,matrix in enumerate(master_list):
				if matind==0: continue	# don't plot s1
				elif matind==7: matind=12
				elif matind==12: matind==7
				for r in matrix: # s1m126 window acc to table 2
					ptcolor = "lightgray"
					if 124<r[24] and r[24]<128:
						ptcolor = "pink"
						if r[36]>r[24]/2 and True: ptcolor = "lightgreen"
						if r[36]>r[24]/2 and False: ptcolor = "red"
						if r[36]<r[24]/2: ptcolor = "magenta"
					if matind == 14: ptcolor = "black"
					# i don't have ggF sigma due to processing
					ax.scatter( r[20], r[19], alpha=1., color=ptcolor, 
						    s=30, label="", marker='.', linewidths=0)
			plt.title(file_prefix+" : (5c) lambda v kappa in s1m126")
			plt.ylabel("lambda")
			plt.xlabel("kappa")
			plt.xlim(0,1)
			plt.ylim(0,1)
			plt.savefig("{}{}/lambda_v_kappa_s1m126.png".format(DIR, save_dir_name))
			plt.close()
			# 5D param plots in s1m126
			# 5E param plots in s1m126
			pltctr+=1
			fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)
			for matind,matrix in enumerate(master_list):
				if matind==0: continue	# don't plot s1
				elif matind==7: matind=12
				elif matind==12: matind==7
				for r in matrix: # s1m126 window acc to table 2
					if (r[23]<600 and -700<r[21] and r[21]<250):	#cull
						ptcolor = "lightgray"
						if 124<r[24] and r[24]<128:
							ptcolor = "pink"
							if r[36]>r[24]/2 and True: ptcolor = "lightgreen"
							if r[36]>r[24]/2 and False: ptcolor = "red"
							if r[36]<r[24]/2: ptcolor = "magenta"
						if matind == 14: ptcolor = "black"
						# i don't have ggF sigma due to processing
						ax.scatter( r[23], r[21], alpha=1., 
							color=ptcolor, s=30, label="", marker='.', linewidths=0)
			plt.title(file_prefix+" : (5e) Alambda v mueff in s1m126")
			plt.ylabel("Alambda")
			plt.xlabel("mueff")
			plt.xlim(100,600)
			plt.ylim(-700,250)
			plt.savefig("{}{}/Alambda_v_mueff_s1m126.png".format(DIR, save_dir_name))
			plt.close()
			# 5F param plots in s1m126
			pltctr+=1
			fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)
			for matind,matrix in enumerate(master_list):
				if matind==0: continue	# don't plot s1
				elif matind==7: matind=12
				elif matind==12: matind==7
				for r in matrix: # s1m126 window acc to table 2
					if (r[23]<600 and -1200<r[22] and r[22]<250):	#cull
						ptcolor = "lightgray"
						if 124<r[24] and r[24]<128:
							ptcolor = "pink"
							if r[36]>r[24]/2 and True: ptcolor = "lightgreen"
							if r[36]>r[24]/2 and False: ptcolor = "red"
							if r[36]<r[24]/2: ptcolor = "magenta"
						if matind == 14: ptcolor = "black"
						# i don't have ggF sigma due to processing
						ax.scatter( r[23], r[22], alpha=1., 
							color=ptcolor, s=30, label="", marker='.', linewidths=0)
			plt.title(file_prefix+" : (5f) Akappa v mueff in s1m126")
			plt.ylabel("Akappa")
			plt.xlabel("mueff")
			plt.xlim(100,600)
			plt.ylim(-1200,250)
			plt.savefig("{}{}/Akappa_v_mueff_s1m126.png".format(DIR, save_dir_name))
			plt.close()
	#	#	#	#	#	#	#	#	#	#	#	#	#
	#	#	#	END OF 13032113 REPL	#	#	#	#	#	#	#
	#	#	#	#	#	#	#	#	#	#	#	#	#	

	#keeping for templating
	#		pltctr+=1
	#		fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)	
	#		for matrix in master_list:
	#			for r in matrix:
	#				if r[28]>122 and r[28]<128:
	#					ax.scatter( r[42], r[36], alpha=.5, color="magenta", 
	#						    s=1, label="p1mass", marker=',', linewidths=0)
	#					ax.scatter( r[42], r[39], alpha=.5, color="orange", 
	#						    s=1, label="p2mass", marker=',', linewidths=0)
	#					ax.scatter( r[42], r[42], alpha=.5, color="cyan", 
	#						    s=1, label="cmass", marker=',', linewidths=0)
	#		plt.title(file_prefix+" : ")
	#		plt.ylabel("")
	#		plt.xlabel("")
	#		plt.xlim(,)
	#		plt.ylim(,)
	#		plt.savefig("{}{}/_v_.png".format(DIR, save_dir_name))
	#		plt.close()

		elif "108035020" in file_prefix:
			print(Time(),"\t108035020-specific plots.")
			pltctr+=1
			HeatPlot(pltctr, "mneu1", 44, "viridis",
					"k div l", 0, -.2, .2,
					"mueff", 23, 0, 0, Size, DPI, "Heatmap","108035020")
			pltctr+=1
			HeatPlot(pltctr, "mueff", 23, "viridis",
					"kappa", 20, 0, 0,
					"Akappa", 22, -500, 500, Size, DPI, "Heatmap","108035020")
			pltctr+=1		
			SinglePlot(pltctr, "mneu1", 44, 100, 650,
					"neu1scomp", 49, .4, 1,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "", "")
			pltctr+=1		
			SinglePlot(pltctr, "mneu2", 50, 100, 650,
					"neu2Bcomp", 51, .4, 1,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "", "")
			pltctr+=1		
			SinglePlot(pltctr, "mneu3", 56, 500, 2250,
					"neu3Hcomp", 0, .4, 1,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "", "")
			pltctr+=1		
			SinglePlot(pltctr, "mneu4", 62, 500, 2250,
					"neu4Hcomp", 0, .4, 1,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "", "")

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

#        [ dne, in1, in2, in1and2, in3, in1and3, in2and3, in1n2n3], [none, 1, 12, 123], or [none, 1, 12]
ms_trk = [[0 for j in range(2**3)] for i in range(len(file_tags))] #tracking which bins sca ~ LHC higgs are in
ls_trk = [[0 for j in range(1+3)] for i in range(len(file_tags))] #      tracking which bins LIGHT pseudos are in
lp_trk = [[0 for j in range(1+2)] for i in range(len(file_tags))] #       (or for the pseudos)
def NearSM(mass): #temp fnc that checks if mass is near sm higgs mass
	buffer = 3.02083333 #buffer in GeV around central SM Higgs mass 
	#buffer = 25 # large buffer for testing new 4constyle
	return (mass > 125.173611-buffer) and (mass < 125.173611+buffer)

file_matrices = [list() for file in file_tags] # for storing "out_file_matrix" of each out file

for file_index,out_file_name in enumerate(file_names):
	with open("{}{}.dat".format(DIR, out_file_name)) as f:
		f_reader = csv.reader(f, delimiter=" ")
		ctr_lighthiggs = 0
		print("{}\tReading in\t{}".format(Time(), out_file_name))
		for indexrow,fullrow in enumerate(f_reader):
			if DEBUG_MODE: 
				if indexrow%100000==0: print(indexrow)
			row = [0] # trim out strange spacing ---> this used to be the event number
			
			reject_row = False
			if "108035020" in file_prefix: last_element=73	#trunc out after neu5scomp
			else: last_element=43				#                MA

			for indexelem,val in enumerate(fullrow):
				if val != "": row.append(float(val))
				if "108035020" in file_prefix and len(row)==21:
					if abs(row[20])/row[19] > 0.15: #outside of allowed kappa/lambda ratio
						reject_row = True
						break
				if len(row)>last_element: break
			if not reject_row: file_matrices[file_index % len(file_tags)].append(row)
			if not MASSTRK: continue # CONTINUING TO IGNORE COUNTING LIGHT/SMLIKE HIGGS EVENTS
			
			params = row[1:24]+[row[43]]
			shiggs = row[24:32+1:4] # s1mass s2mass s3mass
			phiggs = [row[36]]+[row[39]] # p1mass p2mass
			
			in_trk = 0 #ms_trk index
			in_trk += 1*int(NearSM(shiggs[0]))
			in_trk += 2*int(NearSM(shiggs[1]))
			in_trk += 4*int(NearSM(shiggs[2]))
			ms_trk[file_index][in_trk]+=1
			if in_trk==7: print("All scalars in LHC range! #{}:{}.".format(indexrow,out_file_name))
			in_trk = 0 #ls_trk
			for mass in shiggs: in_trk += 1*int(mass<threshold_lighthiggs)
			ls_trk[file_index][in_trk]+=1
			in_trk = 0 #lp_trk
			for mass in phiggs: in_trk += 1*int(mass<threshold_lighthiggs)
			lp_trk[file_index][in_trk]+=1
	#	print(len(file_matrices[file_index%4]))
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
		master_list = [None for i in range(15)]

		print(Time(),"\tbsT\t",int(len(bsT)),end="\t|")
		print("bs1\t",int(len(bs1)),end="\t|")
		print("bs2\t",int(len(bs2)),end="\t|")
		print("bs3\t",int(len(bs3)))

		bT1 = bsT & bs1
		print(Time(),"\tbT1\t",len(bT1))
		bT2 = bsT & bs2
		print(Time(),"\tbT2\t",len(bT2))
		b23 = bs2 & bs3
		print(Time(),"\tb23\t",len(b23))
		b12 = bs1 & bs2
		print(Time(),"\tb12\t",len(b12))
		b13 = bs1 & bs3
		print(Time(),"\tb13\t",len(b13))
		bT3 = bsT & bs3
		print(Time(),"\tbT3\t",len(bT3))
		sT123 = bT1 & b23
		master_list[14]=List(sT123)
		print(Time(),"\t*sT123*\t",len(sT123))
		sT12m = bT1 & b12
		sT12 = sT12m.difference(sT123)
		master_list[10] = List(sT12)
		del sT12m
		print(Time(),"\tsT12\t",len(sT12))
		sT1m3 = bT1 & bT3
		sT13 = sT1m3.difference(sT123)
		master_list[11] = List(sT13)
		del sT1m3
		print(Time(),"\tsT13\t",len(sT13))
		sTm23 = bT3 & b23
		sT23 = sTm23.difference(sT123)
		del sTm23
		master_list[12] = List(sT23)
		print(Time(),"\tsT23\t",len(sT23))
		sm123 = b12 & b23
		s123 = sm123.difference(sT123)
		master_list[13] = List(s123)
		del sm123
		print(Time(),"\ts123\t",len(s123))
		sT1 = bT1.difference(sT12,sT123,sT13)
		master_list[4] = List(sT1)
		print(Time(),"\tsT1\t",len(sT1))
		s23 = b23.difference(sT23,sT123,s123)
		master_list[9] = List(s23)
		print(Time(),"\ts23\t",len(s23))		
		s12 = b12.difference(sT12,sT123,s123)
		master_list[7] = List(s12)
		del b12
		print(Time(),"\ts12\t",len(s12)) # KILLED after this print
		sT3 = bT3.difference(sT23,sT123,sT13)
		del bT3
		master_list[6] = List(sT3)
		print(Time(),"\tsT3\t",len(sT3))
		sT2 = bT2.difference(sT12,sT123,sT23)
		master_list[5] = List(sT2)
		print(Time(),"\tsT2\t",len(sT2))
		del bT2
		s13 = b13.difference(sT13,sT123,s123)
		master_list[8] = List(s13)
		print(Time(),"\ts13\t",len(s13))
		del b13
		del sT123
		sT = bsT.difference(bT1,sT2,sT23,sT3)
		master_list[0] = List(sT)
		print(Time(),"\tsT\t",len(sT))
		del bsT
		del sT23
		s1 = bs1.difference(bT1,s12,s123,s13)
		master_list[1] = List(s1)
		print(Time(),"\ts1\t",len(s1))
		del bs1
		del bT1
		del s123
		s2 = bs2.difference(b23,sT2,sT12,s12)
		master_list[2] = List(s2)
		print(Time(),"\ts2\t",len(s2))
		del bs2
		del sT2
		del sT12
		del s12
		s3 = bs3.difference(b23,sT3,sT13,s13)
		master_list[3] = List(s3)
		print(Time(),"\ts3\t",len(s3))
		del bs3
		del b23
		del s3
		del sT3
		del sT13
		del s13

# args (DO_PARAM, DO_MASS, DO_COMP, DO_HEAT, DO_MISC)
if SAVEPLOTS: 
	if DEBUG_MODE: print(Time(),"\tStarting to plot...")
	GeneratePlots(DO_PARAM, DO_MASS, DO_COMP, DO_HEAT, DO_MISC, DO_REPL)

if MASSTRK:
	print("\nSorting by lightest SCALAR")
	print("tanB\tlambda\t\tkappa\t\tAlambda\tAkappa\t\tmueff\tMass")
	for file_index,out_file_matrix in enumerate(file_matrices):
		sortedbys1mass = sorted(out_file_matrix, key = lambda x: x[24])
		print(file_names[file_index], "\t{} events".format(len(out_file_matrix)))
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
		sortedbyp1mass = sorted(out_file_matrix, key = lambda x: x[36])
		print(file_names[file_index], "\t{} events".format(len(out_file_matrix)))
		for event_index,event in enumerate(sortedbyp1mass):
			if event_index < 3:
				print("{:.4f}\t".format(event[1]),end="")
				print("{:.6f}\t".format(event[19]),end="")
				print("{:.6f}\t".format(event[20]),end="")
				print("{:.2f}\t".format(event[21]),end="")
				print("{:.7f}\t".format(event[22]),end="")
				print("{:.2f}\t".format(event[23]),end="")
				print(event[36])
		print()
	
	print("\nTracking masses in LHC BOUNDS (122.15, 128.19)")
	for file_index,out_file_name in enumerate(file_names):
		print(out_file_name, "\t{} events".format(len(file_matrices[file_index])))
		print("LHC Higgs in:\tNone\ts1\ts2\ts1&2\ts3\ts1&3\ts2&3\ts1&2&3",end="\n\t\t")
		for ms in ms_trk[file_index]: print(ms,end="\t")
		print("\nLight s in:\tNone\ts1\ts1&2\ts1&2&3",end="\n\t\t")
		for ls in ls_trk[file_index]: print(ls,end="\t")
		print("\nLight p in:\tNone\tp1\tp1&2",end="\n\t\t")
		for lp in lp_trk[file_index]: print(lp,end="\t")
		print("\n")
#if not MASSTRK: #ln656 temp for spam runs to filter output
	ls_ctr = 0
	lp_ctr = 0
	for r in master_list[-1]:
		if r[24]<threshold_lighthiggs: ls_ctr+=1
		if r[36]<threshold_lighthiggs: lp_ctr+=1
	print("Light Mass Threshold:\t{}\n# Light Scalars:\t{}\n# Light Pseudoscalars:\t{}".format(threshold_lighthiggs,ls_ctr,lp_ctr))
print("Runtime(s):\t{}\n#=#=#=#=#=#=#=#=#=#=#=#=#=#=#".format(Time()))
#sys.exit()
