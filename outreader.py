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

DEBUG_MODE = 0		#enables print statements used for tracking
MASSTRKFILE = 0		#enables tracking masses near LHC and of light s/o
MASSTRKBOUNDS = 0	# At the end, count higgses below threshold_lighthiggs
BENCH_CHECK = 0

DO_PARAM = 0
DO_MASS = 0
DO_COMP = 0
DO_HEAT = 0
DO_COUP = 0
DO_BR = 0
DO_DC = 0
DO_MISC = 1
DO_REPL = 0

NEU_INFO = 1
SHMIX_INFO = 1
threshold_lighthiggs = 10# GeV
MZ = 91.187
M2Z = MZ**2

file_prefix = argv[1] #=-=# file_prefix = "--"# widep
file_tags = ['THY','LEP','LHC','BKF']#file_tags = ["","con1","con3","con2"]
file_names = ["{}{}randout".format(file_prefix, tag) for tag in file_tags]
save_dir_name = argv[2]

N_EXTRA = 0 # number of extra seeds for files (x4 for actual num extra files)
for ie in range(N_EXTRA):
	extra_names = ["{}_{:0>2}{}randout".format(file_prefix, ie+2, tag) for tag in file_tags]
	file_names = file_names + extra_names

if "108035020" in file_prefix: (KMIN, KMAX, LMIN, LMAX) = (-.015, .015, 0, .1)	#def plot axis window
elif "s" == file_prefix[0]: (KMIN, KMAX, LMIN, LMAX) = (0,.1,0,.5)
elif file_prefix[:6] in ["PQp1v4","PQp1v5"]: (KMIN, KMAX, LMIN, LMAX) = (0,.1,0,.5)
elif "PQ" == file_prefix[0:2]: (KMIN, KMAX, LMIN, LMAX) = (0,.1,0,.7)
else: (KMIN, KMAX, LMIN, LMAX) = (0, 1, 0, 1)
(S1MMIN,S1MMAX,P1MMIN,P1MMAX) = (110,130,0,25)
if "_s2sm" in save_dir_name: (S1MMIN,S1MMAX) = (0,0)
elif "_spread" in save_dir_name: (S1MMIN,S1MMAX) = (0,25)


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
	elap = time() - start_time
	return "{:02d}:{:0>{}}".format(int(elap/60), round(elap%60,1), 4)
	
def AlphaMix(r):	#returns the alpha (hsm Hbsm mixing in scalar Higgs sector) from event r in args
	return (1/2) * np.arctan((2*r[1]/(1-r[1]**2))*( (r[43]**2+M2Z)/(r[43]**2-M2Z) )  )

def FunctionArr(par,expon,ind,matrix):
	if par == "rt n 3 k Ak mueff div lambda":
		return [np.real(cmath.sqrt(-3*r[20]*r[22]*r[23]/r[19])) or r in matrix]
	elif par=="k div l":
		return [r[20]/r[19]			for r in matrix]
	elif par=="neu1Hcomp":
		return [r[47]**2 + r[48]**2	for r in matrix]
	elif par=="neu2Hcomp":
		return [r[53]**2 + r[54]**2		for r in matrix]
	elif par=="neu3Hcomp":
		return [r[59]**2+r[60]**2	for r in matrix]
	elif par=="neu4Hcomp":
		return [r[65]**2 + r[66]**2		for r in matrix]
	elif par=="neu5Hcomp":
		return [r[71]**2 + r[72]**2	for r in matrix]
	elif par=="s1(hsm)comp":
		return [(r[25]*np.cos(AlphaMix(r)) - r[26]*np.sin(AlphaMix(r)))**2 for i,r in enumerate(matrix)]
	elif par=="s1(Hbsm)comp":
		return [(r[25]*np.sin(AlphaMix(r)) + r[26]*np.cos(AlphaMix(r)))**2 for i,r in enumerate(matrix)]
	elif par=="s2(hsm)comp":
		return [(r[29]*np.cos(AlphaMix(r)) - r[30]*np.sin(AlphaMix(r)))**2 for i,r in enumerate(matrix)]
	elif par=="s2(Hbsm)comp":
		return [(r[29]*np.sin(AlphaMix(r)) + r[30]*np.cos(AlphaMix(r)))**2 for i,r in enumerate(matrix)]
	elif par=="s3(hsm)comp":
		return [(r[33]*np.cos(AlphaMix(r)) - r[34]*np.sin(AlphaMix(r)))**2 for i,r in enumerate(matrix)]
	elif par=="s3(Hbsm)comp":
		return [(r[33]*np.sin(AlphaMix(r)) + r[34]*np.cos(AlphaMix(r)))**2 for i,r in enumerate(matrix)]
	elif par=="mcha1 - mneu1":
		return [r[74]-r[44]		for r in matrix]
	elif par=="neu1mass div p1mass":
		return [r[44]/r[36]		for r in matrix]
	elif par=="neu1mass div s1mass":
		return [r[44]/r[24]		for r in matrix]
	else:
		return [r[ind]**expon		for r in matrix]
	

def SinglePlot(pltctr, xpar, xind, xmin, xmax, ypar, yind, ymin, ymax,	
		Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, folder, subfolder): 
	#needs to accept a list of args:
	#	pltctr		# plot count for figuure - considering subplots/axs?
	#	xmin, xmax	# if xmin!=xmax do zoomed plot also
	#	ymin, ymax	#	 ^ analogous
	#	xpar		# x parameter
	#	ypar		# y parameter 
	#	xind		# x par corresponding column number
	#	yind		# y par corresponding column number
	#	Color		# Color array
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
	if "comp" in ypar: y_expon = 2
	else: y_expon = 1
	if "comp" in xpar: x_expon = 2
	else: x_expon = 1
	
	plt.figure(pltctr)
	if CMYK: relevant_matrix = master_list
	else: relevant_matrix = file_matrices	
	extralegelem=0	
	for fx,out_file_matrix in enumerate(relevant_matrix): # file_matrices for RGB, master_list for CMYK	
		if fx==0 and ypar == "rt n 3 k Ak mueff div lambda":	# (only plot y=x once for the LOp1m)
			extralegelem=1
			plt.plot([0,max(Col(xind,file_matrices[-1]))],
				 [0,max(Col(xind,file_matrices[-1]))],
					color="black", alpha=1.0, linestyle="solid",
					linewidth=0.20, label = "y = x")#lw0.15@dpi480
		fn_arr = FunctionArr(ypar,y_expon,yind,out_file_matrix)
		fnxarr = FunctionArr(xpar,x_expon,xind,out_file_matrix)

		plt.scatter(fnxarr, fn_arr,
			alpha=Alpha[fx],color=Color[fx],s=Size[fx],
			label=Label[fx],marker=',',linewidths=0)
	
	plt.title(save_dir_name+" : "+ypar+" v "+xpar)
	plt.ylabel(ypar)
	plt.xlabel(xpar)
	if "dw" == ypar[-2:]: plt.yscale("log")
	elif "XI" == ypar[:2]: plt.yscale("log")
	elif ypar in ["neu1mass div p1mass", "neu1mass div s1mass"]: plt.yscale("log")
	else: plt.yscale("linear")
	if "dw" == xpar[-2:]: plt.xscale("log")
	elif "XI" == xpar[:2]: plt.xscale("log")
	elif xpar in ["neu1mass div p1mass", "neu1mass div s1mass"]: plt.xscale("log")
	else: plt.xscale("linear")

	if (len(Label)+extralegelem) <= 4: Ncols = len(Label)+extralegelem
	else: Ncols = np.ceil( (len(Label)+extralegelem)/2 )
	leg = plt.legend(loc=LOC, bbox_to_anchor=BBOX_TO_ANCHOR, ncols=Ncols, columnspacing=0.7, frameon=False)
	for x in range(len(Label)+extralegelem): leg.legend_handles[x]._sizes = [10]
	
	if (xpar in ["lambda","kappa"] or ypar in ["lambda","kappa"] 
	or (file_prefix == "13032113" and xpar == "Alambda") ):		# If L or K is involved,
		if xpar == "lambda":	plt.xlim(LMIN,LMAX)			#  let first plot be 
		elif ypar == "lambda":	plt.ylim(LMIN,LMAX)			#  confined  (,)
		if xpar == "kappa":		plt.xlim(KMIN,KMAX)			#	
		elif ypar == "kappa":	plt.ylim(KMIN,KMAX)			#  
		if (file_prefix == "13032113"): 
			if xpar == "Alambda": plt.xlim(-700,300)
			elif ypar == "Alambda": plt.ylim(-700,300)
		plt.savefig("{}{}_v_{}.png".format(DIR, ypar, xpar), dpi=DPI)	#  for the L&/or K axi/es
	else: plt.savefig("{}{}_v_{}.png".format(DIR, ypar, xpar), dpi=DPI)	# Otherwise just use DEF.
	
	if xmin==xmax and ymin==ymax: plt.close("all")				# If didn't specify to zoom
	else:									#  just close the plot,
		if xmin!=xmax: plt.xlim(xmin,xmax)				# Otherwise enforce chosen
		if ymin!=ymax: plt.ylim(ymin,ymax)				#  x/y windows
		plt.savefig("{}{}_v_{}_zoom.png".format(DIR, ypar, xpar), dpi=DPI)
		plt.close("all")
	return

def HeatPlot(pltctr, cpar, cind, cmap_n, xpar, xind, xmin, xmax, ypar, yind, ymin, ymax,  # fxn reworked ..
		Size, DPI, folder, subfolder):				#  .. to plot a 2D with a heat map
	
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

	plt.figure(pltctr)		#	vvv used to be file_matrices, changed to master_list for m.ex.sets
	if CMYK: relevant_matrix = master_list[-1] #the fully surviving stuff
	else: relevant_matrix = file_matrices[-1] #just LHC CON
	
	if "13032113" in file_prefix: 
		LEN = len(master_list)
		relevant_matrix = master_list[LEN-6]+master_list[LEN-5]+master_list[LEN-4]+master_list[LEN-3]+master_list[LEN-2]+master_list[LEN-1]

	if "comp" in ypar: y_expon = 2
	else: y_expon = 1
	if "comp" in cpar: c_expon = 2
	else: c_expon = 1		
	if "comp" in xpar: x_expon = 2
	else: x_expon = 1
	
	if cpar == "k div l": Norm = "linear"#Norm="log"# log norm for PQp1v2 and v3 but for v4 try linear
	elif "dw" == cpar[-2:]: Norm = "log"
	elif "XI" == cpar[:2]: Norm = "log"
	elif cpar in ["neu1mass div p1mass", "neu1mass div s1mass"]: Norm = "log"
	else: Norm = "linear"

	if xpar == "p1mass": plt.xscale("linear")#"log")
	elif xpar == "rt n 3 k Ak mueff div lambda": plt.xscale("linear")#"log")
	elif "dw" == xpar[-2:]: plt.xscale("log")
	elif "XI" == xpar[:2]: plt.xscale("log")
	elif xpar in ["neu1mass div p1mass", "neu1mass div s1mass"]: plt.xscale("log")
	else: plt.xscale("linear")

	if ypar == "p1mass": plt.yscale("linear")#"log")
	elif ypar == "rt n 3 k Ak mueff div lambda": plt.yscale("linear")#"log")
	elif "dw" == ypar[-2:]: plt.yscale("log")
	elif "XI" == ypar[:2]: plt.yscale("log")
	elif ypar in ["neu1mass div p1mass", "neu1mass div s1mass"]: plt.yscale("log")
	else: plt.yscale("linear")
	

	
	x_arr = FunctionArr(xpar,x_expon,xind,relevant_matrix)
	y_arr = FunctionArr(ypar,y_expon,yind,relevant_matrix)
	c_arr = FunctionArr(cpar,c_expon,cind,relevant_matrix)

	plt.scatter(x_arr,y_arr,c=c_arr,cmap=cmap_n,norm=Norm,s=Size[-1],marker=',',linewidths=0)
	
	plt.title(save_dir_name+" : "+ypar+" v "+xpar+" c "+cpar)
	plt.ylabel(ypar)
	plt.xlabel(xpar)
	plt.colorbar(label=cpar) #norm keyword
	
	if xpar in ["lambda","kappa"] or ypar in ["lambda","kappa"]:		# If L or K is involved,
		if xpar == "lambda": plt.xlim(LMIN,LMAX)
		elif ypar == "lambda": plt.ylim(LMIN,LMAX)			#  let first plot be 
		if xpar == "kappa": plt.xlim(KMIN,KMAX)
		elif ypar == "kappa": plt.ylim(KMIN,KMAX)			#  confined to (0,0.8)
		plt.savefig("{}{}_v_{}_c_{}.png".format(DIR, ypar, xpar, cpar), dpi=DPI)	#  for the L&/or K axi/es
	else: plt.savefig("{}{}_v_{}_c_{}.png".format(DIR, ypar, xpar, cpar), dpi=DPI)	# Otherwise just use DEF.
	
	if xmin==xmax and ymin==ymax: plt.close("all")				# If didn't specify to zoom
	else:									#  just close the plot,
		if xmin!=xmax: plt.xlim(xmin,xmax)				# Otherwise enforce chosen
		if ymin!=ymax: plt.ylim(ymin,ymax)				#  x/y windows
		plt.savefig("{}{}_v_{}_c_{}_zoom.png".format(DIR, ypar, xpar, cpar), dpi=DPI)
		plt.close("all")
	return

def GeneratePlots(DO_PARAM,DO_MASS,DO_COMP,DO_HEAT,DO_COUP,DO_BR,DO_DC,DO_MISC,DO_REPL):
	# MAJOR IDEAS FOR REFACTORING
	# - REINFORCE THE HEATMAP GENERATION BEING ITS OWN UNIQUE LOOP STRUCTURE
	# -- WOULD POSSIBLY REQUIRE MULTIPLE RUNS FOR TWO DIFFERENT TYPES OF HEATMAPS
	# - COULD POSSIBLY ARCHITECTURE S.T. SINGLEPLOT LOOPS SIMILARLY TO HEATPLOT THRU DO_HEAT
	#    BUT WOULD HAVE TO FIND MORE PRECISE/ELEGANT WAY OF GENERATING X,Y PAIRS RATHER THAN ALL PAIRS
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
		col_scheme = [(.9,0,.9),(0,.85,.85),(.9,.9,0)] # for plots of 3 things against each other

		Alpha = [1 for x in range(NC)]
		dot_size = 0.1;	Size= [dot_size for x in range(NC-1)]+[10*dot_size]#.25 when doing 8 colors
		
		if "13032113" in file_prefix: Size = [10*s for s in Size]
		elif "108035020" in file_prefix: Size = [10*s for s in Size] 
		LOC = "center"
		BBOX_TO_ANCHOR = [0.475,1.105]
		DPI = 240 #for very long time operating at 480 DPI, @1730_14dec23 changed to 240
		######## IF DPI @ 480, SIZE OF 0.04 OK. IF DPI @ 240, DOTS DO NOT RENDER @ THAT SIZE. INC TO 0.1
	pltctr = 0
	par_list = [ ("lambda",19), ("kappa",20), ("Alambda",21), ("mueff",23), ("Akappa",22), ("tanB",1) ]
	if "PQp1" in file_prefix: par_list += [("k div l", 0)]
	#par_list = [("MA",43)] + par_list
	elif "108035020" in file_prefix: par_list = [("AU3",5),("M1",2),("M2",3),("M3",4)] + par_list
	elif "PQv8" in file_prefix: par_list = [("AD3",6),("AU3",5)] + par_list
	elif "PQv9" in file_prefix or "PQv10" == file_prefix: par_list = [("AU3",5)]+par_list
	hmass_list = [ ("s1mass",24), ("s2mass",28), ("s3mass",32), ("p1mass",36), ("p2mass",39), ("cmass",42) ]
	if NEU_INFO:
		neumass_list = [("neu1mass",44), ("neu2mass",50), ("neu3mass",56), 
				("neu4mass",62), ("neu5mass",68), ("cha1mass",74) ]
	elif "PQ" in file_prefix: neumass_list = [("neu1mass",44)]
	else: neumass_list = []
	mass_list = neumass_list[:3] + hmass_list + [("cha1mass",74)]
	
	# after edits to nmhdecay_rand.f these comps are matrix elems not true compositions, ned **2
	# comps also used to just be called comp, but specifically called as singlet comp, change filenm
	scomp_list = [ ("s1scomp",27), ("s2scomp",31), ("s3scomp",35), ("p1scomp",38), ("p2scomp",41) ]
	ucomp_list = [ ("s1ucomp",25), ("s2ucomp",29), ("s3ucomp",33) ]
	dcomp_list = [ ("s1dcomp",26), ("s2dcomp",30), ("s3dcomp",34) ]
	Acomp_list = [ ("p1Acomp",37), ("p2Acomp",40) ]
	
	if SHMIX_INFO: 
		shmix_list = [  ("s1(hsm)comp",0), ("s1(Hbsm)comp",0),
				("s2(hsm)comp",0), ("s2(Hbsm)comp",0),
				("s3(hsm)comp",0), ("s3(Hbsm)comp",0) ]
	else:	shmix_list = []
	
	if NEU_INFO:
		neucomp_list =[ ("neu1Bcomp",45), ("neu1Wcomp",46),("neu1Hcomp",0), ("neu1scomp",49),
				("neu2Bcomp",51), ("neu2Wcomp",52),("neu2Hcomp",0), ("neu2scomp",55),
				("neu3Bcomp",57), ("neu3Wcomp",58),("neu3Hcomp",0), ("neu3scomp",61),
				("neu4Bcomp",63), ("neu4Wcomp",64),("neu4Hcomp",0), ("neu4scomp",67),
				("neu5Bcomp",69), ("neu5Wcomp",70),("neu5Hcomp",0), ("neu5scomp",73) ]
	else:	neucomp_list=[]
	
	heatmap_list = ["turbo"]# #"viridis", "plasma",	#"inferno", "magma", "cividis",	#"brg", "rainbow","jet",
	
	if DO_BR or "PQp1v5"==file_prefix:
		br_list = [ 	("br_neu2_s1neu1", 130),("br_neu2_s2neu1", 131), 
				("br_neu2_p1neu1", 132),("br_neu2_zneu1", 133),
				("br_neu3_s1neu1", 134),("br_neu3_s2neu1", 135),
				("br_neu3_p1neu1", 136),("br_neu3_zneu1", 137),
				("br_cha1_wneu1", 138), ("br_cha1_hcneu1", 139),
				("br_s1_neu1neu1", 120), ("br_p1_neu1neu1", 121) ]
	else: br_list = []
	
	if DO_COUP or "PQp1v5" == file_prefix:
		coup_list =[    ("XIs1u", 86), ("XIs2u", 92),("XIs3u",  98), ("XIp1u", 104),("XIp2u", 110),
				("XIs1d", 87), ("XIs2d", 93),("XIs3d",  99), ("XIp1d", 105),("XIp2d", 111),
				("XIs1z", 88), ("XIs2z", 94),("XIs3z", 100), ("XIp1z", 106),("XIp2z", 112),
				("XIs1gl",89), ("XIs2gl",95),("XIs3gl",101), ("XIp1gl",107),("XIp2gl",113),
				("XIs1ga",90), ("XIs2ga",96),("XIs3ga",102), ("XIp1ga",108),("XIp2ga",114),
				("XIs1b", 91), ("XIs2b", 97),("XIs3b", 103), ("XIp1b", 109),("XIp2b", 115) ]
	else: coup_list = []

	if DO_DC or file_prefix == "PQp1v5":
		dc_list = [("s1dw",140),("s2dw",141),("p1dw",142),("neu2dw",143),("neu3dw",144),("cha1dw",145)]
	else: dc_list = []

	
	if DO_PARAM:
		print(Time(),"Beginning parameter plots")
		for i,(xpar,xind) in enumerate(par_list): # ALL PARAM VS
			for j,(ypar,yind) in enumerate(par_list): #PARAM
				if j<=i: continue
				pltctr+=1
				SinglePlot(pltctr, xpar, xind, 0, 0, ypar, yind, 0, 0,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "Parameter", "")
	if DO_MASS:
		print(Time(),"Beginning mass plots") # PLOT ea mass against chosen params and other masses
		for h,(h_mass,hix) in enumerate(mass_list):
			print(Time(),h_mass)
			
			if "s1" in h_mass: (xmin, xmax) = (S1MMIN,S1MMAX)
			elif "p1" in h_mass: (xmin, xmax) = (P1MMIN,P1MMAX)
			else: (xmin, xmax) = (0, 0)

			for (param,pix) in par_list+scomp_list+ucomp_list+dcomp_list+shmix_list+Acomp_list+neucomp_list[:12]:
				(ymin, ymax) = (0,0)
				if (param,pix) not in par_list:
					if h_mass[:3] in ["neu","cha"] and (param,pix) not in neucomp_list:
						continue
					elif h_mass[:2]!=param[:2]: continue

				pltctr+=1
				SinglePlot(pltctr,h_mass, hix, xmin, xmax,
						param, pix, ymin, ymax,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "Mass",h_mass)
			for yit,(y_higg,yix) in enumerate(mass_list):
				if yit <= h: continue
				
				if "s1" in y_higg: (ymin, ymax) = (S1MMIN,S1MMAX)
				elif "p1" in y_higg: (ymin, ymax) = (P1MMIN,P1MMAX)
				else: (ymin, ymax) = (0, 0)
				
				pltctr+=1
				SinglePlot(pltctr, h_mass, hix, xmin, xmax,
						y_higg, yix, ymin, ymax,
						Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "Mass", "")
	if DO_COMP:
		DO_SCOMP = 0	# Plot singlet comps (of scalars and pseudoscalars)
		DO_UCOMP = 0	#	   u-Higgs comps of scalars
		DO_DCOMP = 0	#	   d-Higgs comps .	.
		DO_SHMIX = 0	#	   Hsm/Hbsm mixing of Hu and Hd
		DO_ACOMP = 0	#	   A_MSSM comps of pseudoscalars
		DO_NCOMP = 1	#	   Neutralino comps

		print(Time(),"Beginning composition plots") # PLOT ea Higgs comps vs each parameter
		for (h_comp,cix) in scomp_list +ucomp_list +dcomp_list +shmix_list +Acomp_list +neucomp_list:
			if not DO_SCOMP and (h_comp,cix) in scomp_list: continue
			if not DO_UCOMP and (h_comp,cix) in ucomp_list: continue
			if not DO_DCOMP and (h_comp,cix) in dcomp_list: continue
			if not DO_SHMIX and (h_comp,cix) in shmix_list: continue
			if not DO_ACOMP and (h_comp,cix) in Acomp_list: continue
			if not DO_NCOMP and (h_comp,cix) in neucomp_list[:12]: continue
			print(Time(),h_comp)
			for (param,pix) in par_list:
				pltctr+=1
				SinglePlot(pltctr,param,pix,0,0,
					h_comp, cix, 0,0,
					Label, Color, Alpha, Size,
					LOC, BBOX_TO_ANCHOR, DPI, "Comp", h_comp)
		
	if DO_HEAT:

		print(Time(),"Beginning heat map plots")
		# (C, X, Y)
		DO_SCOMP = (0,0,0)	# Plotting singlet comps (of scalars and pseudoscalars)
		DO_UCOMP = (0,0,0)	#	   u-Higgs comps of scalars
		DO_DCOMP = (0,0,0)	#	   d-Higgs comps .	.
		DO_SHMIX = (0,0,0)	#	   Hsm/Hbsm mixing of Hu and Hd
		DO_ACOMP = (0,0,0)	#	   A_MSSM comps of pseudoscalars
		DO_NCOMP = (1,0,1)	#	   Neutralino comps
		DO_PARS  = (0,1,0)	#	   Core parameters
		DO_MASSC = (0,0,0)	#	   Neu and H masses
		DO_COUPL = (0,0,0)	#	   reduced couplings
		DO_DECAY = (0,0,0)	#	   decay widths

		toggle_list = []		
		toggle_list += par_list 
		toggle_list += mass_list
		toggle_list += scomp_list
		toggle_list += ucomp_list
		toggle_list += dcomp_list
		toggle_list += shmix_list
		toggle_list += Acomp_list
		toggle_list += neucomp_list[:8]
		toggle_list += coup_list
		toggle_list += dc_list


		togg_labels = ["par" for par in par_list]+["mass" for mass in mass_list]+["scomp" for scomp in scomp_list]+["ucomp" for ucomp in ucomp_list]+["dcomp" for dcomp in dcomp_list]+["shmix" for shmix in shmix_list]+["Acomp" for Acomp in Acomp_list]+["neucomp" for neucomp in neucomp_list[:8]]+["XI" for XI in coup_list]+["DC" for DC in dc_list]
	
		for n,(c_par,c_ix) in enumerate(toggle_list):
			if (c_par,c_ix) in par_list and not DO_PARS[0]: continue
			if (c_par,c_ix) in mass_list and not DO_MASSC[0]: continue
			if (c_par,c_ix) in scomp_list and not DO_SCOMP[0]: continue
			if (c_par,c_ix) in ucomp_list and not DO_UCOMP[0]: continue
			if (c_par,c_ix) in dcomp_list and not DO_DCOMP[0]: continue
			if (c_par,c_ix) in shmix_list and not DO_SHMIX[0]: continue
			if (c_par,c_ix) in Acomp_list and not DO_ACOMP[0]: continue
			if (c_par,c_ix) in neucomp_list and not DO_NCOMP[0]: continue
			if (c_par,c_ix) in coup_list and not DO_COUPL[0]: continue
			if (c_par,c_ix) in dc_list and not DO_DECAY[0]: continue

			print(Time(),"Coloring with",c_par)

			for xind,(x_par,x_ix) in enumerate(toggle_list):

				if (x_par,x_ix) in par_list and not DO_PARS[1]: continue
				if (x_par,x_ix) in mass_list and not DO_MASSC[1]: continue
				if (x_par,x_ix) in scomp_list and not DO_SCOMP[1]: continue
				if (x_par,x_ix) in ucomp_list and not DO_UCOMP[1]: continue
				if (x_par,x_ix) in dcomp_list and not DO_DCOMP[1]: continue
				if (x_par,x_ix) in shmix_list and not DO_SHMIX[1]: continue
				if (x_par,x_ix) in Acomp_list and not DO_ACOMP[1]: continue
				if (x_par,x_ix) in neucomp_list and not DO_NCOMP[1]: continue
				if (x_par,x_ix) in coup_list and not DO_COUPL[1]: continue
				if (x_par,x_ix) in dc_list and not DO_DECAY[1]: continue
			
				if x_par == c_par: continue

				print(Time(),"... Analyzing",x_par)

				if x_par == "s1mass": (x_min, x_max) = (S1MMIN,S1MMAX)
				elif x_par == "p1mass": (x_min, x_max) = (P1MMIN,P1MMAX)
				else: (x_min, x_max) = (0, 0)

				for yind,(y_par,y_ix) in enumerate(toggle_list):
					if (y_par,y_ix) in par_list:
						if not DO_PARS[2]: continue
						if (x_par,x_ix) in par_list:
							if par_list.index((y_par,y_ix))<=par_list.index((x_par,x_ix)):
								continue
					if (y_par,y_ix) in mass_list:
						if not DO_MASSC[2]: continue
						if (x_par,x_ix) in mass_list:
							if mass_list.index((y_par,y_ix))<=mass_list.index((x_par,x_ix)):
								continue
					if (y_par,y_ix) in scomp_list:
						if not DO_SCOMP[2]: continue
						if (x_par,x_ix) in scomp_list:
							if scomp_list.index((y_par,y_ix))<=scomp_list.index((x_par,x_ix)):
								continue
					if (y_par,y_ix) in ucomp_list:
						if not DO_UCOMP[2]: continue
						if (x_par,x_ix) in ucomp_list:
							if ucomp_list.index((y_par,y_ix))<=ucomp_list.index((x_par,x_ix)):
								continue
					if (y_par,y_ix) in dcomp_list:
						if not DO_DCOMP[2]: continue
						if (x_par,x_ix) in dcomp_list:
							if dcomp_list.index((y_par,y_ix))<=dcomp_list.index((x_par,x_ix)):
								continue
					if (y_par,y_ix) in shmix_list:
						if not DO_SHMIX[2]: continue
						if (x_par,x_ix) in shmix_list:
							if shmix_list.index((y_par,y_ix))<=shmix_list.index((x_par,x_ix)):
								continue
					if (y_par,y_ix) in Acomp_list:
						if not DO_ACOMP[2]: continue
						if (x_par,x_ix) in Acomp_list:
							if Acomp_list.index((y_par,y_ix))<=Acomp_list.index((x_par,x_ix)):
								continue
					if (y_par,y_ix) in neucomp_list:
						if not DO_NCOMP[2]: continue
						if (x_par,x_ix) in neucomp_list:
							if (neucomp_list.index((y_par,y_ix))
							<=neucomp_list.index((x_par,x_ix))):
								continue
					if (y_par,y_ix) in coup_list:
						if not DO_COUPL[2]: continue
						if (x_par,x_ix) in coup_list:
							if (coup_list.index((y_par,y_ix))
							<=coup_list.index((x_par,x_ix))):
								continue
					if (y_par,y_ix) in dc_list:
						if not DO_DECAY[2]: continue
						if (x_par,x_ix) in dc_list:
							if (dc_list.index((y_par,y_ix))
							<=dc_list.index((x_par,x_ix))):
								continue

					if y_par == x_par or y_par == c_par: continue
					
					if DEBUG_MODE: print(Time(),"... ... against",y_par)
				
					if y_par == "s1mass": (y_min, y_max) = (S1MMIN,S1MMAX)
					elif y_par == "p1mass": (y_min, y_max) = (P1MMIN,P1MMAX)
					else: (y_min, y_max) = (0, 0)

					sub_dir=togg_labels[xind]+" v "+togg_labels[yind]+" c "+togg_labels[n]

					pltctr+=1
					HeatPlot(pltctr, c_par, c_ix, heatmap_list[n%len(heatmap_list)],
						x_par, x_ix, x_min, x_max,
						y_par, y_ix, y_min, y_max, Size, DPI, "Heatmap",sub_dir)
	if DO_BR:
	
		print(Time(),"Beginning BR plots")
		for (br,brix) in br_list:
			print(Time(),"Evaluating {}...".format(br))
			for (par,pix) in par_list:
				pltctr+=1
				SinglePlot(pltctr, par, pix, 0, 0,
					br, brix, 0, 0,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "BR","Parameter")
				
			for i,(mass,mix) in enumerate([("neu1mass",44),("p1mass",36),("s1mass",24),("cha1mass",74)]):
				print(Time(),"... against {}".format(mass))
				pltctr+=1
				SinglePlot(pltctr, mass, mix, 0, 0,
						br, brix, 0, 0,
						Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "BR","Mass")
				for (comp,cix) in neucomp_list[:12]:
					pltctr+=1
					HeatPlot(pltctr, comp, cix,"turbo",
					mass,mix,0,0,
					br,brix,0,0,
					Size, DPI, "Heatmap","BR v mass c comp")
				for j,(mass2,mix2) in enumerate([("neu1mass",44),("p1mass",36),("s1mass",24),("cha1mass",74)]):
					if j<=i: continue
					pltctr+=1
					HeatPlot(pltctr, mass2, mix2,"turbo",
					mass,mix,0,0,
					br,brix,0,0,
					Size, DPI, "Heatmap","BR v mass c mass")

			pltctr+=1
			HeatPlot(pltctr, br, brix,"turbo",
				"neu1mass",44,0,0,
				"cha1mass",74,0,0,
				 Size, DPI, "Heatmap","Mass v mass c BR")
		
	
		print(Time(),"Starting neu2 BRs by channel") ##################################################
		pltctr+=1
		plt.scatter(	[r[50] for r in master_list[-1]], 
				[r[130]+r[131]+r[132]+r[133] for r in master_list[-1]],
				label="br_sum", color="black", alpha=0.5, s=15, marker='.', linewidths=0)
		for i,(br,brix) in enumerate([("br_neu2_s1neu1",130),("br_neu2_s2neu1",131),("br_neu2_p1neu1",132),("br_neu2_zneu1",133)]):
			plt.scatter([r[50] for r in master_list[-1]], [r[brix] for r in master_list[-1]],
				label=br, alpha=1, s=1, marker=',', linewidths=0)
		plt.title(save_dir_name+" : neu2 BR per channel v neu2mass")
		plt.xlabel("neu2mass")
		plt.ylabel("BR per channel")
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=3,columnspacing=0.7,frameon=False)
		for x in range(5): leg.legend_handles[x]._sizes = [10]
		plt.savefig("{}/{}/BR/neu2 BRs.png".format(DIR,save_dir_name),dpi=DPI)
		plt.close()

		print(Time(),"Starting neu3 BRs by channel") ##################################################
		pltctr+=1
		plt.scatter(	[r[56] for r in master_list[-1]],
				[r[134]+r[135]+r[136]+r[137] for r in master_list[-1]],
				label="br_sum", color="black", alpha=0.5, s=15, marker='.', linewidths=0)
		for i,(br,brix) in enumerate([("br_neu3_s1neu1",134),("br_neu3_s2neu1",135),("br_neu3_p1neu1",136),("br_neu3_zneu1",137)]):
			plt.scatter([r[56] for r in master_list[-1]], [r[brix] for r in master_list[-1]],
				label=br, alpha=1, s=1, marker=',', linewidths=0)
		plt.title(save_dir_name+" : neu3 BR per channel v neu3mass")
		plt.xlabel("neu3mass")
		plt.ylabel("BR per channel")
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=3,columnspacing=0.7,frameon=False)
		for x in range(5): leg.legend_handles[x]._sizes = [10]
		plt.savefig("{}/{}/BR/neu3 BRs.png".format(DIR,save_dir_name),dpi=DPI)
		plt.close()

		print(Time(),"Starting cha1 BRs by channel") ##################################################
		pltctr+=1
		plt.scatter([r[74] for r in master_list[-1]], [r[138]+r[139] for r in master_list[-1]],
				label="br_sum", color="black", alpha=0.5, s=15, marker='.', linewidths=0)
		for i,(br,brix) in enumerate([("br_cha1_wneu1",138),("br_cha1_hcneu1",139)]):
			plt.scatter([r[74] for r in master_list[-1]], [r[brix] for r in master_list[-1]],
				label=br, color=(139-brix,-138+brix,0), alpha=1,s=1, marker=',', linewidths=0)
		plt.title(save_dir_name+" : cha1 BR per channel v cha1mass")
		plt.xlabel("cha1mass")
		plt.ylabel("BR per channel")
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=2,columnspacing=0.7,frameon=False)
		for x in range(3): leg.legend_handles[x]._sizes = [10]
		plt.savefig("{}/{}/BR/cha1 BRs.png".format(DIR,save_dir_name),dpi=DPI)
		plt.close()

		print(Time(),"Starting neu2 BR sums") ##################################################
		pltctr+=1
		plt.scatter(	[r[130]+r[132] for r in master_list[-1]],
				[r[131]+r[133] for r in master_list[-1]],
				color="black", alpha=1, s=1, marker=',', linewidths=0)
		plt.title(save_dir_name+" : neu2 BR sums")
		plt.xlabel("br_neu2_(s1/p1)neu1")
		plt.ylabel("br_neu2_(s2/z)neu1")
		plt.savefig("{}/{}/BR/neu2 BR sums.png".format(DIR,save_dir_name),dpi=DPI)
		plt.close()

		print(Time(),"Starting neu3 BR sums") ##################################################
		pltctr+=1
		plt.scatter(	[r[134]+r[136] for r in master_list[-1]],
				[r[135]+r[137] for r in master_list[-1]],
				color="black", alpha=1, s=1, marker=',', linewidths=0)
		plt.title(save_dir_name+" : neu3 BR sums")
		plt.xlabel("br_neu3_(s1/p1)neu1")
		plt.ylabel("br_neu3_(s2/z)neu1")
		plt.savefig("{}/{}/BR/neu3 BR sums.png".format(DIR,save_dir_name),dpi=DPI)
		plt.close()

	if DO_DC:
		print(Time(),"Doing decay width plots") ##################################################
		dc_masses = [	("s1mass",24),("s2mass",28),("p1mass",36),
				("neu2mass",50),("neu3mass",56),("cha1mass",74)	]

		for i,(dw,dwix) in enumerate(dc_list):
			pltctr+=1
			SinglePlot(pltctr, dc_masses[i][0], dc_masses[i][1], 0, 0, dw, dwix, 0, 0,
				    	Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "DC", "Mass")
			for (par,pix) in par_list:
				SinglePlot(pltctr, par, pix, 0, 0, dw, dwix, 0, 0,
					    	Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "DC", "Parameter")


	if DO_COUP:
		print(Time(),"Beginning XI plots")
		for (coup,cix) in coup_list:
			print(Time(),"Evaluating {}...".format(coup))
			for (par,pix) in par_list:		
				pltctr+=1
				SinglePlot(pltctr, par, pix, 0, 0,
						coup, cix, 0, 0,
						Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "XI","Parameter")
				for (mass,mix,mmin,mmax) in [("s1mass",24,0,50),("p1mass",36,0,25)]:#,("neu1mass",44,0,0)]:		
					pltctr+=1
					HeatPlot(pltctr, mass, mix,"turbo",
						par, pix, 0, 0,
						coup, cix, 0, 0, Size, DPI, "Heatmap","XI v par c mass")
					pltctr+=1
					HeatPlot(pltctr, par,pix,"turbo",
						mass,mix, mmin, mmax,
						coup, cix, 0, 0, Size, DPI, "Heatmap","XI v mass c par")

			for (mass,pix) in [("s1mass",24),("p1mass",36)]:#,("neu1mass",44)]:		
				pltctr+=1
				SinglePlot(pltctr, mass, pix, 0, 0,coup,cix, 0, 0,
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "XI","Mass")
			
			if DO_DC: #definitely a more pretty wy to do this, but plot each particles coup against its decay width
				if coup[2:4]=="s1":
					pltctr+=1
					SinglePlot(pltctr, "s1mass", 24, 0, 0,coup,cix, 0, 0,
						Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI,"XI","DC")
				elif coup[2:4]=="s2":
					pltctr+=1
					SinglePlot(pltctr, "s2mass", 28, 0, 0,coup,cix, 0, 0,
						Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI,"XI","DC")
				elif coup[2:4]=="p1":
					pltctr+=1
					SinglePlot(pltctr, "p1mass", 36, 0, 0,coup,cix, 0, 0,
						Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI,"XI","DC")

	

	if DO_MISC:
		print(Time(),"Initiating miscellany.")
		pltctr+=1
		SinglePlot(pltctr, "neu1mass div p1mass", 0, 0, 0, "p1dw", 142, 0, 0,
				Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "DC","")
		pltctr+=1
		SinglePlot(pltctr, "br_p1_neu1neu1", 121, 0, 0, "p1dw", 142, 0, 0,
				Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "DC","")
		pltctr+=1
		SinglePlot(pltctr, "neu1mass div s1mass", 0, 0, 0, "s1dw", 140, 0, 0,
				Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "DC","")
		pltctr+=1
		SinglePlot(pltctr, "br_s1_neu1neu1", 120, 0, 0, "s1dw", 140, 0, 0,
				Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "DC","")
		print(Time(),"p1 heats")
		pltctr+=1
		HeatPlot(pltctr, "br_p1_neu1neu1", 121, "turbo",
			"neu1mass div p1mass", 0,0,0,
			"p1dw", 142, 0, 0, Size, DPI, "Heatmap","")	
		pltctr+=1
		HeatPlot(pltctr, "br_p1_neu1neu1", 121, "turbo",
			"neu1mass div p1mass", 0,0,0,
			"p1mass", 36, 0, 0, Size, DPI, "Heatmap","")	
		pltctr+=1
		HeatPlot(pltctr, "neu1mass div p1mass", 0, "turbo",
			"p1mass", 36, 0, 0,
			"br_p1_neu1neu1", 121,0,0, Size, DPI, "Heatmap","")	
		print(Time(),"s1 heats")
		pltctr+=1
		HeatPlot(pltctr, "br_s1_neu1neu1", 120, "turbo",
			"neu1mass div s1mass", 0,0,0,
			"s1dw", 140, 0, 0, Size, DPI, "Heatmap","")	
		pltctr+=1
		HeatPlot(pltctr, "br_s1_neu1neu1", 120, "turbo",
			"neu1mass div s1mass", 0,0,0,
			"s1mass", 24, 0, 0, Size, DPI, "Heatmap","")	
		pltctr+=1
		HeatPlot(pltctr, "neu1mass div s1mass", 0, "turbo",
			"s1mass", 24, 0, 0,
			"br_s1_neu1neu1", 120,0,0, Size, DPI, "Heatmap","")		
		print(Time(),"Doing kdl histos")
		pltctr+=1
		plt.figure(pltctr)
		plt.hist([r[20]/r[19] for r in master_list[-1]], color="r", rwidth=.9)
		plt.title(save_dir_name+" : distribution of events in k div l")
		plt.xlabel("k div l")
		plt.ylabel("N")
		plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/Parameter/hist k div l.png".format(save_dir_name),dpi=DPI)
		plt.close()

		return # EARLY RETURN TO NOT OVERFLOW
		print(Time(),"Comparing LO p1mass")
		pltctr+=1
		SinglePlot(pltctr, "p1mass", 36, P1MMIN,P1MMAX,
				"rt n 3 k Ak mueff div lambda", 0, P1MMIN,P1MMAX,
				Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "","")
		c_pars_list = par_list + scomp_list + shmix_list
		for (c_par,c_ix) in c_pars_list:
			HeatPlot(pltctr,c_par, c_ix, "turbo",
				"p1mass",36,P1MMIN,P1MMAX,
				"rt n 3 k Ak mueff div lambda", 0, P1MMIN, P1MMAX, Size, DPI, "Heatmap","")	
		for (smass, sindex) in [("s1mass",24),("s2mass",28),("s3mass",32)]:
			print(Time(),"{}(u,d,s)comp v {}".format(smass[0:2],smass))
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
			plt.close("all")		

			print(Time(),"{}(h,H,s)comp v {}".format(smass[0:2],smass))
			pltctr+=1
			fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)	
	
			shsmcomp = list() #hsmcomps,
			sHbsmcomp = list() # and Hbsm comps for each event
			for i,r in enumerate(master_list[-1]):
				shsmcomp.append(r[sindex+1]*np.cos(AlphaMix(r))-r[sindex+2]*np.sin(AlphaMix(r)))
				sHbsmcomp.append(r[sindex+1]*np.sin(AlphaMix(r))+r[sindex+2]*np.cos(AlphaMix(r)))	
			
			for color,comp in enumerate(["{}scomp".format(smass[0:2]), "{}(hsm)comp".format(smass[0:2]), "{}(Hbsm)comp".format(smass[0:2])]):
				if comp == "{}(hsm)comp".format(smass[0:2]):
					fn_arr = [R**2 for R in shsmcomp]
				elif comp == "{}(Hbsm)comp".format(smass[0:2]):
					fn_arr = [R**2 for R in sHbsmcomp]
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
		print(Time(),"p1(A,s)comp v p1mass")
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
	
		print(Time(),"p2(A,s)comp v p2mass")
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
				print(Time(),"s123{} v {}".format(jcomp,scalarmass))
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
			print(Time(),"Plotting Figure (3)...")
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
			print(Time(),"Plotting Figure (4)...")
			# 4A param plots in s1m126			
			pltctr+=1
			fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)	
			matcolors = ["gray" for arr in master_list]		# other		## hopefully removing
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
			print(Time(),"Plotting Figure (5)...")
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
	#							s=1, label="p1mass", marker=',', linewidths=0)
	#					ax.scatter( r[42], r[39], alpha=.5, color="orange", 
	#							s=1, label="p2mass", marker=',', linewidths=0)
	#					ax.scatter( r[42], r[42], alpha=.5, color="cyan", 
	#							s=1, label="cmass", marker=',', linewidths=0)
	#		plt.title(file_prefix+" : ")
	#		plt.ylabel("")
	#		plt.xlabel("")
	#		plt.xlim(,)
	#		plt.ylim(,)
	#		plt.savefig("{}{}/_v_.png".format(DIR, save_dir_name))
	#		plt.close()

		elif "108035020" in file_prefix:
			print(Time(),"108035020-specific plots.")
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

#		 [ dne, in1, in2, in1and2, in3, in1and3, in2and3, in1n2n3], [none, 1, 12, 123], or [none, 1, 12]
ms_trk = [[0 for j in range(2**3)] for i in range(len(file_tags))] #tracking which bins sca ~ LHC higgs are in
ls_trk = [[0 for j in range(1+3)] for i in range(len(file_tags))] #		 tracking which bins LIGHT pseudos are in
lp_trk = [[0 for j in range(1+2)] for i in range(len(file_tags))] #		  (or for the pseudos)
def NearSM(mass): #temp fnc that checks if mass is near sm higgs mass
	buffer = 3.02083333 #buffer in GeV around central SM Higgs mass 
	#buffer = 25 # large buffer for testing new 4constyle
	return (mass > 125.173611-buffer) and (mass < 125.173611+buffer)

file_matrices = [list() for file in file_tags] # for storing "out_file_matrix" of each out file

for file_index,out_file_name in enumerate(file_names):
	with open("{}{}.dat".format(DIR, out_file_name)) as f:
		f_reader = csv.reader(f, delimiter=" ")
		ctr_lighthiggs = 0
		print(Time(),"Reading in: {}...".format(out_file_name))
		for indexrow,fullrow in enumerate(f_reader):
			if DEBUG_MODE: 
				if indexrow%100000==0: print(indexrow)
			row = [0] # trim out strange spacing ---> this used to be the event number

			reject_row = False
			if "108035020" in file_prefix: last_element=74	#trunc out after MCHA(1)
			elif DO_DC or "PQp1v5" == file_prefix[:6]: last_element=1000	#(don't trunc)
			elif DO_BR: last_element = 140
			elif NEU_INFO: last_element=74			#		 MCHA(1)
			elif "PQ" in file_prefix: last_element=44	#		 neu1mass
			else: last_element=43				#		 MA
			for indexelem,val in enumerate(fullrow):
				if val != "": row.append(float(val))
				# THIS WHERE TO ENFORCE CUTS AND FILTERS
				if "108035020" in file_prefix and len(row)==21:
					if abs(row[20])/row[19] > 0.15:
						reject_row = True
						break
				elif "PQp1v3" == file_prefix and len(row)==21:
					if row[20]/row[19] > 2:
						reject_row = True
						break
				elif file_prefix in ["PQp1v4"] and len(row)==21:
					if row[20]/row[19] > 1:
						reject_row = True
						break
				elif "s1scomp_ge_p9" in save_dir_name and len(row)==28:
					if abs(row[27]) < 0.9:
						reject_row = True
						break
				elif file_prefix in ["PQp1v5","PQp1v6"]:
					if len(row)==21 and row[20]/row[19] > 1: #redundant
						reject_row = True		# since doing the generation
						break				#  with this in mind

					if "_s1sm" in save_dir_name and len(row)==25 and not NearSM(row[24]):
						reject_row = True
						break
					elif "_s2sm" in save_dir_name and len(row)==29 and not NearSM(row[28]):
						reject_row = True
						break
					elif ("_bino" in save_dir_name and		#binolike
					(len(row)==46 and row[45]**2<0.05) ):		# keep 5%+ Bino
						reject_row = True
						break
					elif ("_compress" in save_dir_name and (	#compress
					(len(row)==45 and row[44]<150        ) or	# hino bound 150
					(len(row)==46 and row[45]**2>=0.05   ) or	# reject binolike 5%+
					(len(row)==75 and row[74]-row[44]>=10)	) ):	# neu1~cha1 mass
						reject_row = True
						break
					elif ("_spread" in save_dir_name and (		#spread
					(len(row)==46 and row[45]**2>=0.05   )or	# reject binolike 5%+
					(len(row)==75 and row[74]-row[44]<10 ) ) ):	# reject compressed
						reject_row = True	# want to triple check w group,
						break			# need to reject Hino<400 here??

				if len(row)>last_element: break

			if not reject_row and last_element > 128: row[128]=0 # SCUFFED MANUAL OVERWRITE OF GGF H SIGMA	
			if not reject_row: file_matrices[file_index % len(file_tags)].append(row)
			if not MASSTRKFILE: continue # CONTINUING TO IGNORE COUNTING LIGHT/SMLIKE HIGGS EVENTS
			
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
	print(Time(),"Splitting into sets...")
	
	def Set(List): # fn is List to set of tuples conversion
		return set(map(tuple, List))
	def List(Set): # fn is Set to List conversion
		return list(map(list, Set))
	#for i,e in enumerate(file_matrices[2]):	#BE WEARY THIS IS OVERWRITING ORIGINAL INFO
	#	e[-2]=0				# ON THE PROD SIGMA THRU GGF / con2 has nonzero
			# above arg corresponds to LHC FILE, is [3] for [ "","con1","con3","con2"]

	if False: #leaving this here but copying this architecture for the THY/LEP/LHC/BKF idea...
		bset0 = Set(file_matrices[0])		# base unc. set
		bset1 = Set(file_matrices[1])		# base con1 set
		bset3 = Set(file_matrices[2])		# base con3 set
		bset2 = Set(file_matrices[3])		# base con2 set
		
		set132 = bset1 & bset3 & bset2	#union all cons
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
		#				 running with 16 COLORS??? ouch (EDIT: actually 15 since no 0con set)
		print(Time(),"Base sets and intersection")		
		bT = Set(file_matrices[0])				# base sets exactly 1 con. applied
		b1 = Set(file_matrices[1])
		b2 = Set(file_matrices[2])
		b3 = Set(file_matrices[3])
		master_list = [None for i in range(15)]

		print(Time(),"bT  ) {: >6}".format(len(bT)),end=" |")
		print("b1  ) {: >6}".format(len(b1)),end=" |")
		print("b2  ) {: >6}".format(len(b2)),end=" |")
		print("b3  ) {: >6}".format(len(b3)))

		bT1 = bT & b1
		if DEBUG_MODE: print(Time(),"bT1	",len(bT1))
		bT2 = bT & b2
		if DEBUG_MODE: print(Time(),"bT2	",len(bT2))
		b23 = b2 & b3
		if DEBUG_MODE: print(Time(),"b23	",len(b23))
		b12 = b1 & b2
		if DEBUG_MODE:	print(Time(),"b12	 ",len(b12))
		b13 = b1 & b3
		if DEBUG_MODE:	print(Time(),"b13	 ",len(b13))
		bT3 = bT & b3
		if DEBUG_MODE:	print(Time(),"bT3	 ",len(bT3))
		
		print(Time(),"bT1 ) {: >6}".format(len(bT1)),end=" |")
		print("bT2 ) {: >6}".format(len(bT2)),end=" |")
		print("bT3 ) {: >6}".format(len(bT3)),end=" |")
		print("b12 ) {: >6}".format(len(b12)),end=" |")
		print("b13 ) {: >6}".format(len(b13)),end=" |")
		print("b23 ) {: >6}".format(len(b23)))
		
		print(Time(),"Mutually exclusive subsets")
		sT123 = bT1 & b23
		master_list[14]=List(sT123)
		print(Time(),"SURVIVING ALL CONSTRAINTS: sT123) {: >6}".format(len(sT123)))
		sT12m = bT1 & b12
		sT12 = sT12m.difference(sT123)
		master_list[10] = List(sT12)
		del sT12m
		if DEBUG_MODE:	print(Time(),"sT12	 ",len(sT12))
		sT1m3 = bT1 & bT3
		sT13 = sT1m3.difference(sT123)
		master_list[11] = List(sT13)
		del sT1m3
		if DEBUG_MODE:	print(Time(),"sT13	 ",len(sT13))
		sTm23 = bT3 & b23
		sT23 = sTm23.difference(sT123)
		del sTm23
		master_list[12] = List(sT23)
		if DEBUG_MODE:	print(Time(),"sT23	 ",len(sT23))
		sm123 = b12 & b23
		s123 = sm123.difference(sT123)
		master_list[13] = List(s123)
		del sm123
		if DEBUG_MODE:	print(Time(),"s123	 ",len(s123))
		
		print(Time(),"sT12) {: >6}".format(len(sT12)),end=" |")
		print("sT13) {: >6}".format(len(sT13)),end=" |")
		print("sT23) {: >6}".format(len(sT23)),end=" |")
		print("s123) {: >6}".format(len(s123)))
		
		sT1 = bT1.difference(sT12,sT123,sT13)
		master_list[4] = List(sT1)
		if DEBUG_MODE: print(Time(),"sT1	",len(sT1))
		s23 = b23.difference(sT23,sT123,s123)
		master_list[9] = List(s23)
		if DEBUG_MODE: print(Time(),"s23	",len(s23))		
		s12 = b12.difference(sT12,sT123,s123)
		master_list[7] = List(s12)
		del b12
		if DEBUG_MODE: print(Time(),"s12	",len(s12)) # KILLED after this print
		sT3 = bT3.difference(sT23,sT123,sT13)
		del bT3
		master_list[6] = List(sT3)
		if DEBUG_MODE: print(Time(),"sT3	",len(sT3))
		sT2 = bT2.difference(sT12,sT123,sT23)
		master_list[5] = List(sT2)
		if DEBUG_MODE: print(Time(),"sT2	",len(sT2))
		del bT2
		s13 = b13.difference(sT13,sT123,s123)
		master_list[8] = List(s13)
		if DEBUG_MODE: print(Time(),"s13	",len(s13))
		del b13
		del sT123
	
		print(Time(),"sT1 ) {: >6}".format(len(sT1)),end=" |")
		print("sT2 ) {: >6}".format(len(sT2)),end=" |")
		print("sT3 ) {: >6}".format(len(sT3)),end=" |")
		print("s12 ) {: >6}".format(len(s12)),end=" |")
		print("s13 ) {: >6}".format(len(s13)),end=" |")
		print("s23 ) {: >6}".format(len(s23)))
	
		sT = bT.difference(bT1,sT2,sT23,sT3)
		master_list[0] = List(sT)
		if DEBUG_MODE: print(Time(),"sT		",len(sT))
		del bT
		del sT23
		s1 = b1.difference(bT1,s12,s123,s13)
		master_list[1] = List(s1)
		if DEBUG_MODE: print(Time(),"s1		",len(s1))
		del b1
		del bT1
		del s123
		s2 = b2.difference(b23,sT2,sT12,s12)
		master_list[2] = List(s2)
		if DEBUG_MODE: print(Time(),"s2		",len(s2))
		del b2
		del sT2
		del sT12
		del s12
		s3 = b3.difference(b23,sT3,sT13,s13)
		master_list[3] = List(s3)
		if DEBUG_MODE: print(Time(),"s3		",len(s3))
		del b3
		del b23
		print(Time(),"sT  ) {: >6}".format(len(sT)),end=" |")
		print("s1  ) {: >6}".format(len(s1)),end=" |")
		print("s2  ) {: >6}".format(len(s2)),end=" |")
		print("s3  ) {: >6}".format(len(s3)))
		del s3
		del sT3
		del sT13
		del s13
		del s1
		del s2
		del sT
		# # # # ONCE sT123 DETERMINED...	
		ahH=[list() for i in master_list]			# calculate alpha value (h-H mixing)
		for i,dataset in enumerate(master_list):		# for each event in each set...
			for r in dataset:
				M2A = r[43]**2
				ahH[i].append((1/2)*np.arctan((2*r[1]/(1-r[1]**2))*( (M2A+M2Z)/(M2A-M2Z) )	))

if SAVEPLOTS: 
	if DEBUG_MODE: print(Time(),"Starting to plot...")
	GeneratePlots(DO_PARAM,DO_MASS,DO_COMP,DO_HEAT,DO_COUP,DO_BR,DO_DC,DO_MISC,DO_REPL)

if MASSTRKFILE:
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
if MASSTRKBOUNDS:
	for threshold_lighthiggs in [9]:
		ls_ctr = 0
		lp_ctr = 0
		for r in master_list[-1]:
			if r[24]<threshold_lighthiggs:
				ls_ctr+=1
				print("s1mass {: >4} < {: >2} : ".format(round(r[24],1),threshold_lighthiggs), r[1], r[19], r[20], r[21], r[22], r[23])
			if r[36]<threshold_lighthiggs:
				lp_ctr+=1
				print("p1mass {: >4} < {: >2} : ".format(round(r[36],1),threshold_lighthiggs), r[1], r[19], r[20], r[21], r[22], r[23])

		print("Light Mass Threshold:\t{}\n# Light Scalars:\t{}\n# Light Pseudoscalars:\t{}".format(threshold_lighthiggs,ls_ctr,lp_ctr))
	
	(t_lo, t_hi) = (master_list[-1][0][1], master_list[-1][0][1]) 
	(l_lo, l_hi) = (master_list[-1][0][19],master_list[-1][0][19])
	(k_lo, k_hi) = (master_list[-1][0][20],master_list[-1][0][20])
	(Al_lo,Al_hi)= (master_list[-1][0][21],master_list[-1][0][21])
	(Ak_lo,Ak_hi)= (master_list[-1][0][22],master_list[-1][0][22])
	(mu_lo,mu_hi)= (master_list[-1][0][23],master_list[-1][0][23])
	for r in master_list[-1]: 
		(t_lo, t_hi) = (min(t_lo, r[1]), max(t_hi, r[1]))
		(l_lo, l_hi) = (min(l_lo, r[19]),max(l_hi, r[19]))
		(k_lo, k_hi) = (min(k_lo, r[20]),max(k_hi, r[20]))
		(Al_lo,Al_hi)= (min(Al_lo,r[21]),max(Al_hi,r[21]))
		(Ak_lo,Ak_hi)= (min(Ak_lo,r[22]),max(Ak_hi,r[22]))
		(mu_lo,mu_hi)= (min(mu_lo,r[23]),max(mu_hi,r[23]))
	print("Parameter bounds: ( min ~~ max )")
	print("tanB   ({: >8} ~~ {: >8} )".format(t_lo, t_hi))
	print("lambda ({:.2E} ~~ {:.2E} )".format(l_lo, l_hi))	
	print("kappa  ({:.2E} ~~ {:.2E} )".format(k_lo, k_hi))	
	print("Alambda({: >8} ~~ {: >8} )".format(Al_lo,Al_hi)) 
	print("Akappa ({: >8} ~~ {: >8} )".format(Ak_lo,Ak_hi))
	print("mueff  ({: >8} ~~ {: >8} )".format(mu_lo,mu_hi))
#{:0>{}}
if BENCH_CHECK:
	print("BENCHMARK POINTS BY VALUES:     tanB    lambda        kappa  Alambda    Akappa    mueff")
	for i,r in enumerate(master_list[-1]):
#		if r[130]+r[131]+r[132]+r[133] < 1: #brneu2tot
#		if r[134]+r[135]+r[136]+r[137] < 1: #brneu3tot
#		if r[138]+r[139]<1:	#br_cha1_wneu1 + _hcneu1
#		if r[24]<122 and r[36]<15: #lower s1 p1 masses
#		if r[24]<10:			# small s1mass
#			if FunctionArr("neu1Hcomp",0,0,master_list[-1])[i] < 0.5:
		if r[142]<1E-13:	#tiny p1dw
			print("s1mass {: >5.1f} & p1mass {: >4}: {: >8.5f} {: >9} {: >12} {: >8} {: >9} {: >8}".format(round(r[24],1),round(r[36],1),r[1],r[19],r[20],r[21],r[22],r[23]))
print("{}\tFinished.\n#=#=#=#=#=#=#=#=#=#=#=#=#=#=#".format(Time()))
#sys.exit()
