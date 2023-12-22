import numpy as np
import cmath
import matplotlib.pyplot as plt
from time import time
import glob
import pandas as pd
import csv 
import os
import sys

# ARGS FROM COMMAND LINE [,,,,,], execute this file as: python3 outreader.py file_prefix save_dir_name SAVEPLOTS
# sys.argv = ['outreader.py', file_prefix, save_dir_name, SAVEPLOTS]
## currently have CMYK automatically set as True, not needed in arguments
argv = sys.argv

DEBUG_MODE = 0 #enables print statements used for tracking
MASSTRK = 0 #enables tracking masses near LHC and of light s/o
DO_PARAM = 1
DO_MASS = 1
DO_COMP = 1
DO_HEAT = 1 
DO_MISC = 1

threshold_lighthiggs = 50 # GeV
file_prefix = argv[1] #=-=# file_prefix = "--"# widep
file_tags = ['THY','LEP','LHC','BKF']#file_tags = ["","con1","con3","con2"]
file_names = ["{}{}randout".format(file_prefix, tag) for tag in file_tags]

if "108035020" in file_prefix: (KMIN, KMAX, LMIN, LMAX) = (-.015, .015, 0, .1)	#def plot axis window
elif "s" == file_prefix[0]: (KMIN, KMAX, LMIN, LMAX) = (0,.1,0,.5)
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
		Size= [.1 for x in range(NC-1)]+[10*.1]#.25 when doing 8 colors
		if "13032113" in file_prefix: Size = [10*s for s in Size] 
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
		if "13032113" in file_prefix: extra_pars_list=[("MA",43)]
		elif "108035020" in file_prefix: extra_pars_list = [("AU3", 5), ("M1", 2),("M2", 3),("M3", 4)]
		else: extra_pars_list=[]
		for i,(xpar,xind) in enumerate(extra_pars_list+par_list): # ALL PARAM VS
			for j,(ypar,yind) in enumerate(extra_pars_list+par_list): #PARAM
				if j<=i: continue
				pltctr+=1
				SinglePlot(pltctr, xpar, xind, 0, 0, ypar, yind, 0, 0,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "Parameter", "")
	if DO_MASS:
		print(Time(),"\tBeginning mass plots") # PLOT ea Higgs' mass against each parameter also its singlet comp
		for h,(h_mass,hix) in enumerate(mass_list):
			print("{}\t{}".format(Time(),h_mass))
			if "13032113" in file_prefix: extra_pars_list = [("AU3",5), ("MQ3",13)]
			elif "108035020" in file_prefix: extra_pars_list = [("AU3", 5), ("M1", 2),("M2", 3),("M3", 4)]
			else: extra_pars_list = []
			for (param,pix) in par_list+comp_list+extra_pars_list: #c_l[h] does higgs v own comp, jus c_l v all comps
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
			print("{}\t{}".format(Time(),h_comp))
			if "108035020" in file_prefix: extra_pars_list=[("AU3", 5), ("M1",2),("M2",3),("M3",4)]
			else: extra_pars_list = []
			for (param,pix) in par_list+extra_pars_list:
				pltctr+=1
				SinglePlot(pltctr,param,pix,0,0,
						h_comp, cix, 0,0,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "Comp", h_comp)
	if DO_HEAT:
		print(Time(),"\tBeginning heat map plots") #heatmaps for s1 to look @ tanB region underneath main LHC blob
		if "13032113" in file_prefix: extra_pars_list = [("AU3",5), ("MQ3",13)]
		elif "108035020" in file_prefix: extra_pars_list = [("AU3", 5), ("M1", 2),("M2", 3),("M3", 4)]
		else: extra_pars_list = []
		for n,(c_par,c_ix) in enumerate(par_list+comp_list+extra_pars_list):# params as heatmap choice
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

		if "108035020" in file_prefix:
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
ms_trk = [[0 for j in range(2**3)] for i in range(num_files)] #tracking which bins scalars like LHC higgs are in
ls_trk = [[0 for j in range(1+3)] for i in range(num_files)] #      tracking which bins LIGHT pseudos are in
lp_trk = [[0 for j in range(1+2)] for i in range(num_files)] #       (or for the pseudos)
def NearSM(mass): #temp fnc that checks if mass is near sm higgs mass
	buffer = 3.02083333 #buffer in GeV around central SM Higgs mass 
	#buffer = 25 # large buffer for testing new 4constyle
	return (mass > 125.173611-buffer) and (mass < 125.173611+buffer)

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
			out_file_matrix.append(row[0:75])

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
