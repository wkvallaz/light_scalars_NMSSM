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
DO_BR = 1
DO_DC = 0
DO_MISC = 0
DO_REPL = 0

NEU_INFO = 1
SHMIX_INFO = 1
threshold_lighthiggs = 10# GeV
MZ = 91.187
M2Z = MZ**2

file_prefix = argv[1] #=-=# file_prefix = "--"# widep
file_tags = ['THY','LEP','LHC','BKF']#file_tags = ["","con1","con3","con2"]
#file_tags = ['ALL'] # for when want to have only one run that already compensates for all overlapping
file_names = ["{}{}randout".format(file_prefix, tag) for tag in file_tags]
save_dir_name = argv[2]


N_EXTRA=3
extra_names = ["lighthiggs_0{}{}randout".format(j,tag)  for j in [4,5,6] for tag in file_tags]
file_names = file_names+extra_names
#N_EXTRA = 3 # number of extra seeds for files (x4 for actual num extra files)
#for ie in range(N_EXTRA):
#	extra_names = ["{}_{:0>2}{}randout".format(file_prefix, ie+2, tag) for tag in file_tags]
#	file_names = file_names + extra_names



(KMIN, KMAX, LMIN, LMAX) = (0, 1, 0, 1)
if "108035020" in file_prefix: (KMIN, KMAX, LMIN, LMAX) = (-.015, .015, 0, .1)	#def plot axis window
elif "s" == file_prefix[0]: (KMIN, KMAX, LMIN, LMAX) = (0,.1,0,.5)
elif file_prefix[:6] in ["PQp1v4","PQp1v5"]: (KMIN, KMAX, LMIN, LMAX) = (0,.1,0,.5)
elif file_prefix[:6] in ["PQp1v6"]: (KMIN, KMAX, LMIN, LMAX) = (0,.05,0,.5)
elif file_prefix[:6] in ["PQp1v8"]: (KMIN, KMAX, LMIN, LMAX) = (0,.025,0,.5)
elif ("PQp1v9" in file_prefix or
	"lighthiggs" in file_prefix): (KMIN, KMAX, LMIN, LMAX) = (0,.015,0,.3)
elif "PQ" == file_prefix[0:2]: (KMIN, KMAX, LMIN, LMAX) = (0,.1,0,.7)

(S1MMIN,S1MMAX,P1MMIN,P1MMAX) = (110,130,0,25)
if "_s2sm" in save_dir_name: (S1MMIN,S1MMAX) = (0,0)
elif ("PQp1v8" in save_dir_name or
	"PQp1v9" in save_dir_name or
	"lighthiggs" in save_dir_name):
	(S1MMIN,S1MMAX,P1MMIN,P1MMAX) = (0.050,11,0.050,11)
elif ("_spread" in save_dir_name or 
	"_H-s" in save_dir_name): (S1MMIN,S1MMAX) = (0,25)

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
	elif par=="s1mass div p1mass":
		return [r[24]/r[36]		for r in matrix]
	else:
		return [r[ind]**expon		for r in matrix]
	

def SinglePlot(pltctr, xpar, xind, xmin, xmax, xscale, ypar, yind, ymin, ymax, yscale,	
		Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, folder, subfolder, **kwargs): 
	#needs to accept a list of args:
	#	pltctr		# int		: plot count for figuure - considering subplots/axs?
	#	xmin, xmax	# float 0	: if xmin!=xmax do zoomed plot also
	#	ymin, ymax	# int	0	:^ analogous
	#	xpar		# str x parameter
	#	ypar		# str	: y parameter 
	#	xind		# int 0	: x par corresponding column number
	#	yind		# int 0	: y par corresponding column number
	#	xscale		# str	: x par choice for plt.xscale(xscale)
	#	yscale		# str	: y par choice for plt.yscale(yscale)
	#	Color		# [str]	: Color array
	#	Alpha		# [int]	: Alpha array
	#	Size		# [int]	: Size array
	#	Label		# [str]	: Label array
	#	LOC		# float : Legend location
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
	
	titletag = ""
	for key,val in kwargs.items():
		if key == "titletag":
			titletag = val
	
	if "dw" == ypar[-2:]: plt.yscale("log")
	elif "XIp1"==ypar[:4]: plt.yscale("symlog")
	elif "XI" == ypar[:2]: plt.yscale("log")
	elif ypar in ["neu1mass div p1mass", "neu1mass div s1mass","s1mass div p1mass"]: plt.yscale("log")
	else: plt.yscale("linear")
	plt.yscale(yscale)

	if "dw" == xpar[-2:]: plt.xscale("log")
	elif "XIp1"==xpar[:4]: plt.xscale("symlog")
	elif "XI" == xpar[:2]: plt.xscale("log")
	elif xpar in ["neu1mass div p1mass", "neu1mass div s1mass","s1mass div p1mass"]: plt.xscale("log")
	else: plt.xscale("linear")
	plt.xscale(xscale)

	if (len(Label)+extralegelem) <= 4: Ncols = len(Label)+extralegelem
	else: Ncols = np.ceil( (len(Label)+extralegelem)/2 )
	leg = plt.legend(loc=LOC, bbox_to_anchor=BBOX_TO_ANCHOR, ncols=Ncols, columnspacing=0.7, frameon=False)
	for x in range(len(Label)+extralegelem): leg.legend_handles[x]._sizes = [10]
	
	if (xpar in ["lambda","kappa"] or ypar in ["lambda","kappa"] 
	or (file_prefix == "13032113" and xpar == "Alambda") ):		# If L or K is involved,
		if xpar == "lambda":	plt.xlim(LMIN,LMAX)			#  let first plot be 
		elif ypar == "lambda":	plt.ylim(LMIN,LMAX)			#  confined  (,)
		if xpar == "kappa":	plt.xlim(KMIN,KMAX)			#	
		elif ypar == "kappa":	plt.ylim(KMIN,KMAX)			#  
		if (file_prefix == "13032113"): 
			if xpar == "Alambda": plt.xlim(-700,300)
			elif ypar == "Alambda": plt.ylim(-700,300)
		plt.savefig("{}{}_v_{}{}.png".format(DIR, ypar, xpar,titletag), dpi=DPI)	#  for the L&/or K axi/es
	else: plt.savefig("{}{}_v_{}{}.png".format(DIR, ypar, xpar, titletag), dpi=DPI)	# Otherwise just use DEF.
	
	if xmin==xmax and ymin==ymax: plt.close("all")				# If didn't specify to zoom
	else:									#  just close the plot,
		if xmin!=xmax: plt.xlim(xmin,xmax)				# Otherwise enforce chosen
		if ymin!=ymax: plt.ylim(ymin,ymax)				#  x/y windows
	#	if "mass" in xpar: 						# Change scale lin/log on crop
	#		plt.xscale("log")
	#	if "mass" in ypar:
	#		plt.yscale("log")
		plt.savefig("{}{}_v_{}{}_zoom.png".format(DIR, ypar, xpar, titletag), dpi=DPI)
		plt.close("all")
	return

def HeatPlot(pltctr, cpar, cind, cmap_n, xpar, xind, xmin, xmax, xscale, ypar, yind, ymin, ymax, yscale,  # fxn reworked ..
		Size, DPI, folder, subfolder, **kwargs):				#  .. to plot a 2D with a heat map
	
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
	elif cpar in ["p1mass","s1mass"]: Norm = "linear"
	elif "dw" == cpar[-2:]: Norm = "log"
	elif "XIp1"==cpar[:4]: Norm = "symlog"
	elif "XI" == cpar[:2]: Norm = "log"
	elif cpar in ["neu1mass div p1mass", "neu1mass div s1mass","s1mass div p1mass"]: Norm = "log"
	else: Norm = "linear"

	if xpar in ["p1mass","s1mass"]: plt.xscale("log")
	elif xpar == "rt n 3 k Ak mueff div lambda": plt.xscale("linear")#"log")
	elif "dw" == xpar[-2:]: plt.xscale("log")
	elif "XIp1"==xpar[:4]: plt.xscale("symlog")
	elif "XI" == xpar[:2]: plt.xscale("log")
	elif xpar in ["neu1mass div p1mass", "neu1mass div s1mass","s1mass div p1mass"]: plt.xscale("log")
	else: plt.xscale("linear")
	plt.xscale(xscale)

	if ypar in ["p1mass", "s1mass"]: plt.yscale("log")
	elif ypar == "rt n 3 k Ak mueff div lambda": plt.yscale("linear")#"log")
	elif "dw" == ypar[-2:]: plt.yscale("log")
	elif "XIp1"==ypar[:4]: plt.yscale("symlog")
	elif "XI" == ypar[:2]: plt.yscale("log")
	elif ypar in ["neu1mass div p1mass", "neu1mass div s1mass","s1mass div p1mass"]: plt.yscale("log")
	else: plt.yscale("linear")
	plt.yscale(yscale)


	x_arr = FunctionArr(xpar,x_expon,xind,relevant_matrix)
	y_arr = FunctionArr(ypar,y_expon,yind,relevant_matrix)
	c_arr = FunctionArr(cpar,c_expon,cind,relevant_matrix)

	plt.scatter(x_arr,y_arr,c=c_arr,cmap=cmap_n,norm=Norm,s=Size[-1],marker=',',linewidths=0)
	
	plt.title(save_dir_name+" : "+ypar+" v "+xpar+" c "+cpar)
	plt.ylabel(ypar)
	plt.xlabel(xpar)
	plt.colorbar(label=cpar) #norm keyword
	
	titletag = ""
	for key,val in kwargs.items():
		if key == "titletag":
			titletag = val

	if xpar in ["lambda","kappa"] or ypar in ["lambda","kappa"]:		# If L or K is involved,
		if xpar == "lambda": plt.xlim(LMIN,LMAX)
		elif ypar == "lambda": plt.ylim(LMIN,LMAX)			#  let first plot be 
		if xpar == "kappa": plt.xlim(KMIN,KMAX)
		elif ypar == "kappa": plt.ylim(KMIN,KMAX)			#  confined to (0,0.8)
		plt.savefig("{}{}_v_{}_c_{}{}.png".format(DIR, ypar, xpar, cpar, titletag), dpi=DPI)	#  for the L&/or K axi/es
	else: plt.savefig("{}{}_v_{}_c_{}{}.png".format(DIR, ypar, xpar, cpar, titletag), dpi=DPI)	# Otherwise just use DEF.
	
	if xmin==xmax and ymin==ymax: plt.close("all")				# If didn't specify to zoom
	else:									#  just close the plot,
		if xmin!=xmax: plt.xlim(xmin,xmax)				# Otherwise enforce chosen
		if ymin!=ymax: plt.ylim(ymin,ymax)				#  x/y windows
	#	if "mass" in xpar: 						# Change scale lin/log on crop
	#		plt.xscale("log")
	#	if "mass" in ypar:
	#		plt.yscale("log")
		plt.savefig("{}{}_v_{}_c_{}{}_zoom.png".format(DIR, ypar, xpar, cpar, titletag), dpi=DPI)
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
		if num_files > 1:
			Label = [ 'T','1','2','3',
			  'T1','T2','T3','12','13','23',
			  'T12','T13','T23','123',
			  'T123' ]
			Color =[(.8,.8,.8), (.5,1,1), (1,.5,1), (1,1,.5),
			(0,.9,.9), (.9,0,.9), (.9,.9,0),
			(.2,.6,1), (.125,1,.125), (1,.125,.125),
			(0,0,.9), (0,.6,0), (.6,0,0), (.4,.4,.4), 
			(0,0,0)]
		else:
			Label = ["Valid"]
			Color = ["k"]
		NC = len(Label)
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
	if "lighthiggs" in file_prefix or "PQp1" in file_prefix: par_list += [("k div l", 0)]
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
	
	if DO_COUP or file_prefix[:6] in ["PQp1v5","PQp1v6","PQp1v7","PQp1v8","PQp1v9","lighthiggs"]:
		coup_list =[    ("XIs1u", 75), ("XIs2u", 81),("XIs3u", 87), ("XIp1u", 93),("XIp2u",  99),
				("XIs1d", 76), ("XIs2d", 82),("XIs3d", 88), ("XIp1d", 94),("XIp2d", 100),
				("XIs1z", 77), ("XIs2z", 83),("XIs3z", 89), ("XIp1z", 95),("XIp2z", 101),
				("XIs1gl",78), ("XIs2gl",84),("XIs3gl",90), ("XIp1gl",96),("XIp2gl",102),
				("XIs1ga",79), ("XIs2ga",85),("XIs3ga",91), ("XIp1ga",97),("XIp2ga",103),
				("XIs1b", 80), ("XIs2b", 86),("XIs3b", 92), ("XIp1b", 98),("XIp2b", 104) ]
	else: coup_list = []

	if DO_BR or file_prefix[:6] in ["PQp1v5","PQp1v6","PQp1v8","PQp1v9"] or "lighthiggs" in file_prefix:
		br_list = [	("br_s1_hadr",107),	("br_p1_hadr",108),
				("br_s1_bb",109),	("br_p1_bb",110),
				("br_s1_cc",111),	("br_p1_cc",112),
				("br_s1_tata",113),	("br_p1_tata",114),
				("br_s1_mumu",115),	("br_p1_mumu",116),
				("br_s1_ee",117),	("br_p1_ee",118),
				("br_s1_gamgam",119),	("br_p1_gamgam",120),
				("br_s1_neu1neu1",121),	("br_p1_neu1neu1",122),

				("br_s2_s1s1",105),	("br_s1_p1p1",106),

				("br_neu2_s1neu1",123),	("br_neu2_s2neu1",124),
				("br_neu2_p1neu1",125),	("br_neu2_zneu1",126),
				("br_neu3_s1neu1",127),	("br_neu3_s2neu1",128),
				("br_neu3_p1neu1",129),	("br_neu3_zneu1",130),
				("br_cha1_wneu1",131),	("br_cha1_hcneu1",132)	]
	else: br_list = []
	s1_brs = []
	for (br,brix) in br_list:
		if br[3:5]=="s1": s1_brs.append((br,brix))
	p1_brs = []
	for (br,brix) in br_list:
		if br[3:5]=="p1": p1_brs.append((br,brix))

	if DO_DC or file_prefix[:6] in ["PQp1v5","PQp1v6"]:
		dc_list = [("s1dw",133),("s2dw",134),("p1dw",135),("neu2dw",136),("neu3dw",137),("cha1dw",138)]
	else: dc_list = []
	
	pdc_list = [	("GamHhadr",	139),	("GamAhadr",	140),
			("GamH2Pi",	141), 	("GamH2PiC",	142),	("GamHPiE",	143),
			("GamHPiEP",	144),	("GamH2KC",	145),	("GamH2K0",	146),
			("GamH2E",	147),	("GamH2EP",	148),	("GamHEEP",	149),
			("GamHss",	150),	("GamHjj",	151),	("GamHchic1p",	152),
			("GamHhadrcc",	153),	("GamHchib12p",	154),	("GamHhadrbb",	155),
			("GamA3Pi",	156),	("GamAPi3PiC",	157),	("GamAEPi3",	158),
			("GamAEPiC",	159),	("GamAEPPi3",	160),	("GamAEPPiC",	161),
			("GamAPiEE",	162),	("GamAPiEEP",	163),	("GamAPiEPEP",	164),
			("GamA3E",	165),	("GamAE2EP",	166),	("GamAEEP2",	167),
			("GamA3EP",	168),	("GamAPiKC",	169),	("GamAPiK0",	170),
			("GamAPiKCK0",	171),	("GamAEKC",	172),	("GamAEK0",	173),
			("GamAEPKC",	174),	("GamAEPK0",	175),	("GamARhogam",	176),
			("GamAss",	177),	("GamAjj",	178),		
			("GamAetac1s",	179),	("GamAhadrcc",	180),	("GamAetab123s",181),
			("GamAhadrbb",	182) ]
	Hpdc_list = [(gam,gix) for (gam,gix) in pdc_list if gam[:4]=="GamH"]
	Apdc_list = [(gam,gix) for (gam,gix) in pdc_list if gam[:4]=="GamA"]
	
	if DO_PARAM:
		print(Time(),"Beginning parameter plots")
		for i,(xpar,xind) in enumerate(par_list): # ALL PARAM VS
			for j,(ypar,yind) in enumerate(par_list): #PARAM
				if j<=i: continue
				pltctr+=1
				SinglePlot(pltctr, xpar, xind, 0, 0, "linear", ypar, yind, 0, 0, "linear",
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
				SinglePlot(pltctr,h_mass, hix, xmin, xmax, "log",
						param, pix, ymin, ymax, "linear",
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "Mass",h_mass)
			for yit,(y_higg,yix) in enumerate(mass_list):
				if yit <= h: continue
				
				if "s1" in y_higg: (ymin, ymax) = (S1MMIN,S1MMAX)
				elif "p1" in y_higg: (ymin, ymax) = (P1MMIN,P1MMAX)
				else: (ymin, ymax) = (0, 0)
				
				pltctr+=1
				SinglePlot(pltctr, h_mass, hix, xmin, xmax, "log",
						y_higg, yix, ymin, ymax, "log",
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
				SinglePlot(pltctr,param,pix,0,0,"linear",
					h_comp, cix, 0,0,"linear",
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

				if x_par == "s1mass":
					(x_min, x_max) = (S1MMIN,S1MMAX)
					x_sc = "log"
				elif x_par == "p1mass":
					(x_min, x_max) = (P1MMIN,P1MMAX)
					x_sc = "log"
				else:
					(x_min, x_max) = (0, 0)
					x_sc = "linear"

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
				
					if y_par == "s1mass":
						(y_min, y_max) = (S1MMIN,S1MMAX)
						y_sc = "log"
					elif y_par == "p1mass": 
						(y_min, y_max) = (P1MMIN,P1MMAX)
						y_sc = "log"
					else: 
						(y_min, y_max) = (0, 0)
						y_sc = "linear"

					sub_dir=togg_labels[xind]+" v "+togg_labels[yind]+" c "+togg_labels[n]

					pltctr+=1
					HeatPlot(pltctr, c_par, c_ix, heatmap_list[n%len(heatmap_list)],
						x_par, x_ix, x_min, x_max, x_sc,
						y_par, y_ix, y_min, y_max, y_sc, Size, DPI, "Heatmap",sub_dir)
	if DO_BR:
	
		print(Time(),"Beginning BR plots")
		for (br,brix) in br_list:
			
			if br[3:5] not in ["s1","p1"]: continue			#temp to just spam s1/p1

			print(Time(),"Evaluating {}...".format(br))
			for (par,pix) in par_list:	# br proportionality to nmssm parameters
				continue # to spam v mass ignore this
				pltctr+=1
				SinglePlot(pltctr, par, pix, 0, 0,"linear",
					br, brix, 0, 0,"linear",
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "BR","Parameter")
				
			for i,(mass,mix) in enumerate(mass_list):
	
				if mass[:2] not in ["s1","p1"]: continue	#temp to just spam s1/p1
			
				if ((mass[:3] in ["neu","cha"] and mass[:4] not in br ) or		# plotting BRs for
				    (mass[:3] not in ["neu","cha"] and mass[:2] not in br )): continue 	#  involved particles


				print(Time(),"... against {}".format(mass))
				if mass[:2]=="s1": (mmin,mmax) = (S1MMIN,S1MMAX)
				elif mass[:2]=="p1": (mmin,mmax) = (P1MMIN,P1MMAX)
				else: (mmin,mmax)=(0,0)

				pltctr+=1
				SinglePlot(pltctr, mass, mix, mmin,mmax,"log",
						br, brix, 0, 0,"linear",
						Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "BR","Mass")
				continue # temp bc don't cary about below
				if mass=="s1mass":
					pltctr+=1
					HeatPlot(pltctr, "s1mass div p1mass", 0, "turbo", mass, mix, mmin, mmax, "log", br, brix,0,0, "linear",
						Size, DPI, "Heatmap","BR v mass c mass")

				for (comp,cix) in neucomp_list[:12]:
					continue # don't currently care about heatmapping these		# TEMP CONTINUE #
					pltctr+=1
					HeatPlot(pltctr, comp, cix,"turbo",
					mass,mix,0,0,"log",
					br,brix,0,0,"linear",
					Size, DPI, "Heatmap","BR v mass c comp")
				for j,(mass2,mix2) in enumerate([("neu1mass",44),("p1mass",36),("s1mass",24),("cha1mass",74)]):
					continue # don't currently care about heatmapping these		# TEMP CONTINUE #
					if j<=i: continue
					pltctr+=1
					HeatPlot(pltctr, mass2, mix2,"turbo",
					mass,mix,0,0,"log",
					br,brix,0,0,"linear",
					Size, DPI, "Heatmap","BR v mass c mass")

			continue #don't care about this heatmap
			pltctr+=1
			HeatPlot(pltctr, br, brix,"turbo",
				"neu1mass",44,0,0,"linear",
				"cha1mass",74,0,0,"linear",
				 Size, DPI, "Heatmap","Mass v mass c BR")
		# addl crops
		pltctr+=1
		SinglePlot(pltctr, "s1mass", 24, 9, 11, "linear", "br_s1_hadr", 107, 0, 0, "linear",
			    	Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "BR","Mass crop",titletag=
"10G lin-lin")
		pltctr+=1
		SinglePlot(pltctr, "s1mass", 24, 0.250, 4, "linear", "br_s1_hadr", 107, 0, 0,"linear",
			    	Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "BR","Mass crop",titletag="2G lin-lin")
		pltctr+=1
		SinglePlot(pltctr, "p1mass", 36, 9, 11, "linear","br_p1_hadr", 108, 0, 0,"linear",
			    	Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "BR","Mass crop",titletag="10G lin-lin")
		pltctr+=1
		SinglePlot(pltctr, "p1mass", 36, 0.400, 5, "linear", "br_p1_hadr", 108, 0, 0,"linear",
			    	Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "BR","Mass crop",titletag="3G lin-lin")
	
		print(Time(),"Starting s1 BRs by channel") ############
		pltctr+=1
		plt.figure(pltctr)
		plt.scatter(	[r[24] for r in master_list[-1]],
				[ sum([r[brix] for (br,brix) in s1_brs]) for r in master_list[-1]],
				label="br sum", color="black", alpha=1, s=10, marker="_")
		for i,(br,brix) in enumerate(s1_brs):
			plt.scatter( [r[24] for r in master_list[-1]], [r[brix] for r in master_list[-1]],
				label=br[3:], alpha=1,s=1,marker=',',linewidths=0)
		plt.title(save_dir_name+" : s1 BR per channel v s1mass")
		plt.xlabel("s1mass")
		plt.ylabel("BR per channel")
		plt.yscale("log")
		plt.ylim(1E-9,1.5)
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=int((len(s1_brs)/2+1)),columnspacing=0.7,frameon=False)
		for x in range(len(s1_brs)+1): leg.legend_handles[x]._sizes = [10]
		plt.savefig("{}/{}/BR/s1 BRs.png".format(DIR,save_dir_name),dpi=DPI)
		plt.yscale("linear")
		plt.ylim(-.02,1.02)
		plt.savefig("{}/{}/BR/s1 BRs lin-y.png".format(DIR,save_dir_name),dpi=DPI)
		plt.yscale("log")
		plt.ylim(1E-9,1.5)
		plt.xlim(S1MMIN,S1MMAX)
		plt.xscale("log")
		plt.savefig("{}/{}/BR/s1 BRs zoom.png".format(DIR,save_dir_name),dpi=DPI)
		plt.ylim(0.005,2)	# copying LSAF plot
		plt.xlim(0.100,10.0)
		plt.savefig("{}/{}/BR/s1 BRs LSAF match.png".format(DIR,save_dir_name),dpi=DPI)
		plt.xlim(S1MMIN,S1MMAX)	# return to not-LSAF xwindow
		plt.yscale("linear")
		plt.ylim(-.02,1.02)
		plt.savefig("{}/{}/BR/s1 BRs zoom lin-y.png".format(DIR,save_dir_name),dpi=DPI)
		plt.xlim(9,11)
		plt.xscale("linear")
		plt.savefig("{}/{}/BR/s1 BRs 10G lin-y.png".format(DIR,save_dir_name),dpi=DPI)
		plt.yscale("log")
		plt.ylim(1E-9,1.5)
		plt.savefig("{}/{}/BR/s1 BRs 10G.png".format(DIR,save_dir_name),dpi=DPI)
		plt.close()
	

		print(Time(),"Starting p1 BRs by channel") ############
		pltctr+=1
		plt.figure(pltctr)
		plt.scatter(	[r[36] for r in master_list[-1]],
			[ sum([r[brix] for (br,brix) in p1_brs]) for r in master_list[-1]],
			label="br sum", color="black", alpha=1, s=10, marker="_")
		for i,(br,brix) in enumerate(p1_brs):
			plt.scatter( [r[36] for r in master_list[-1]], [r[brix] for r in master_list[-1]],
				label=br[3:], alpha=1,s=1,marker=',',linewidths=0)
		plt.title(save_dir_name+" : p1 BR per channel v p1mass")
		plt.xlabel("p1mass")
		plt.ylabel("BR per channel")
		plt.yscale("log")
		plt.ylim(1E-9,1.5)
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=int((len(p1_brs)/2+1)),columnspacing=0.7,frameon=False)
		for x in range(len(p1_brs)+1): leg.legend_handles[x]._sizes = [10]
		plt.savefig("{}/{}/BR/p1 BRs.png".format(DIR,save_dir_name),dpi=DPI)
		plt.yscale("linear")
		plt.ylim(-.02,1.02)
		plt.savefig("{}/{}/BR/p1 BRs lin-y.png".format(DIR,save_dir_name),dpi=DPI)
		plt.yscale("log")
		plt.ylim(1E-9,1.5)
		plt.xlim(P1MMIN,P1MMAX)
		plt.xscale("log")
		plt.savefig("{}/{}/BR/p1 BRs zoom.png".format(DIR,save_dir_name),dpi=DPI)
		plt.ylim(0.005,2)	# copying LSAF plot
		plt.xlim(0.100,10.0)
		plt.xscale("log")
		plt.savefig("{}/{}/BR/p1 BRs LSAF match.png".format(DIR,save_dir_name),dpi=DPI)
		plt.xlim(P1MMIN,P1MMAX)	# return to not-LSAF xwindow
		plt.yscale("linear")
		plt.ylim(-.02,1.02)
		plt.savefig("{}/{}/BR/p1 BRs zoom lin-y.png".format(DIR,save_dir_name),dpi=DPI)
		plt.xlim(9,11)
		plt.xscale("linear")
		plt.savefig("{}/{}/BR/p1 BRs 10G lin-y.png".format(DIR,save_dir_name),dpi=DPI)
		plt.yscale("log")
		plt.ylim(1E-9,1.5)
		plt.savefig("{}/{}/BR/p1 BRs 10G.png".format(DIR,save_dir_name),dpi=DPI)
		plt.close()
		
		######### TO PLOT BR PER CHANNEL OF P1 BUT IN DIFFERENT RANGES OF TANB ##########

		print(Time(),"Starting p1 BRs by channel (tanB slices)")

		for (tanbmin,tanbmax) in [(2,5),(5,10),(10,50)]:
			pltctr+=1
			plt.figure(pltctr)
			plt.scatter(	[r[36] for r in master_list[-1] if (tanbmin<r[1] and r[1]<=tanbmax)],
				[ sum([r[brix] for (br,brix) in p1_brs]) for r in master_list[-1] if (tanbmin<r[1] and r[1]<=tanbmax)] ,
				label="br sum", color="black", alpha=1, s=10, marker="_")
			for i,(br,brix) in enumerate(p1_brs):
				plt.scatter( [r[36] for r in master_list[-1] if (tanbmin<r[1] and r[1]<=tanbmax)], [r[brix] for r in master_list[-1] if (tanbmin<r[1] and r[1]<=tanbmax)],
					label=br[3:], alpha=1,s=1,marker=',',linewidths=0)
			plt.title(save_dir_name+" : p1 BR per channel v p1mass tanB {}-{}".format(tanbmin,tanbmax))
			plt.xlabel("p1mass")
			plt.ylabel("BR per channel")
			plt.yscale("log")
			plt.ylim(1E-9,1.5)
			leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=int((len(p1_brs)/2+1)),columnspacing=0.7,frameon=False)
			for x in range(len(p1_brs)+1): leg.legend_handles[x]._sizes = [10]
			plt.savefig("{}/{}/BR/p1 BRs tanB {}-{}.png".format(DIR,save_dir_name,tanbmin,tanbmax),dpi=DPI)
			plt.yscale("linear")
			plt.ylim(-.02,1.02)
			plt.savefig("{}/{}/BR/p1 BRs tanb {}-{} lin-y.png".format(DIR,save_dir_name,tanbmin,tanbmax),dpi=DPI)
			plt.yscale("log")
			plt.ylim(1E-9,1.5)
			plt.xlim(P1MMIN,P1MMAX)
			plt.xscale("log")
			plt.savefig("{}/{}/BR/p1 BRs tanb {}-{} zoom.png".format(DIR,save_dir_name,tanbmin,tanbmax),dpi=DPI)
			plt.ylim(0.005,2)	# copying LSAF plot
			plt.xlim(0.100,10.0)
			plt.xscale("log")
			plt.savefig("{}/{}/BR/p1 BRs tanb {}-{} LSAF match.png".format(DIR,save_dir_name,tanbmin,tanbmax),dpi=DPI)
			plt.xlim(P1MMIN,P1MMAX)	# return to not-LSAF xwindow
			plt.yscale("linear")
			plt.ylim(-.02,1.02)
			plt.savefig("{}/{}/BR/p1 BRs tanb {}-{} zoom lin-y.png".format(DIR,save_dir_name,tanbmin,tanbmax),dpi=DPI)

		######### ENDING  BR PER CHANNEL OF P1 BUT IN DIFFERENT RANGES OF TANB ##########



	if False: # to keep the text below but not do it
		print(Time(),"Starting neu2 BRs by channel") ##################################################
		pltctr+=1
		plt.figure(pltctr)
		neu2_brs = []
		for (br,brix) in br_list:
			if br[3:7]=="neu2": neu2_brs.append((br,brix))
		plt.scatter(	[r[44] for r in master_list[-1]],
				[ sum([r[brix] for (br,brix) in neu2_brs]) for r in master_list[-1]],
				label="br sum", color="black", alpha=1, s=10, marker="_")
		for i,(br,brix) in enumerate(neu2_brs):
			plt.scatter( [r[44] for r in master_list[-1]], [r[brix] for r in master_list[-1]],
				label=br[3:], alpha=1,s=1,marker=',',linewidths=0)
		plt.title(save_dir_name+" : neu2 BR per channel v neu2mass")
		plt.xlabel("neu2mass")
		plt.ylabel("BR per channel")
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=3,columnspacing=0.7,frameon=False)
		for x in range(5): leg.legend_handles[x]._sizes = [10]
		plt.savefig("{}/{}/BR/neu2 BRs.png".format(DIR,save_dir_name),dpi=DPI)
		plt.close()

		print(Time(),"Starting neu3 BRs by channel") ##################################################
		pltctr+=1
		plt.figure(pltctr)
		neu3_brs = []
		for (br,brix) in br_list:
			if br[3:7]=="neu3": neu3_brs.append((br,brix))
		plt.scatter(	[r[50] for r in master_list[-1]],
				[ sum([r[brix] for (br,brix) in neu3_brs]) for r in master_list[-1]],
				label="br sum", color="black", alpha=1, s=10, marker="_")
		for i,(br,brix) in enumerate(neu3_brs):
			plt.scatter( [r[50] for r in master_list[-1]], [r[brix] for r in master_list[-1]],
				label=br[3:], alpha=1,s=1,marker=',',linewidths=0)
		plt.title(save_dir_name+" : neu3 BR per channel v neu3mass")
		plt.xlabel("neu3mass")
		plt.ylabel("BR per channel")
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=3,columnspacing=0.7,frameon=False)
		for x in range(5): leg.legend_handles[x]._sizes = [10]
		plt.savefig("{}/{}/BR/neu3 BRs.png".format(DIR,save_dir_name),dpi=DPI)
		plt.close()

		print(Time(),"Starting cha1 BRs by channel") ##################################################
		pltctr+=1
		plt.figure(pltctr)
		cha1_brs = []
		for (br,brix) in br_list:
			if br[3:7]=="cha1": cha1_brs.append((br,brix))
		plt.scatter(	[r[74] for r in master_list[-1]],
				[ sum([r[brix] for (br,brix) in cha1_brs]) for r in master_list[-1]],
				label="br sum", color="black", alpha=1, s=10, marker="_")
		for i,(br,brix) in enumerate(cha1_brs):
			plt.scatter( [r[74] for r in master_list[-1]], [r[brix] for r in master_list[-1]],
				label=br[3:], alpha=1,s=1,marker=',',linewidths=0)
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

	if DO_DC and False: # to skip to back of DC =- extremely temp
		print(Time(),"Doing decay width plots") ##################################################
		dc_masses = [	("s1mass",24,S1MMIN,S1MMAX),("s2mass",28,0,0),("p1mass",36,P1MMIN,P1MMAX),
				("neu2mass",50,0,0),("neu3mass",56,0,0),("cha1mass",74,0,0)	]
		for i,(dw,dwix) in enumerate(dc_list):

			#continue
			
			print(Time(),dw)
			pltctr+=1
			SinglePlot(pltctr, dc_masses[i][0], dc_masses[i][1], dc_masses[i][2], dc_masses[i][3], "log", 
					dw, dwix, 0, 0,"log",
				    	Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "DC", "Mass")
			if "s1" in dw or "p1" in dw:

				pltctr+=1
				SinglePlot(pltctr, dc_masses[i][0], dc_masses[i][1], 9, 11, "linear", 
						dw, dwix, 0, 0,"log",
					    	Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "DC", "Mass crop 10G log-lin")
				if "s1" in dw:
					pltctr+=1
					SinglePlot(pltctr, dc_masses[i][0], dc_masses[i][1], 0.250, 4, "linear", 
							dw, dwix, 0, 0,"log",
						    	Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "DC", "Mass crop 2G log-lin")
				if "p1" in dw:
					pltctr+=1
					SinglePlot(pltctr, dc_masses[i][0], dc_masses[i][1], 0.400, 5, "linear", 
							dw, dwix, 0, 0,"log",
						    	Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "DC", "Mass crop 3G log-lin")
			continue # for speed of spam
			for (par,pix) in par_list:
				pltctr+=1
				SinglePlot(pltctr, par, pix, 0, 0, "linear", dw, dwix, 0, 0,"log",
					    	Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "DC", "Parameter")
			for (br,brix) in br_list:
				if ( (dw[:3] in ["neu","cha"] and dw[:4]==br[3:7]) or 
				     (dw[:3] not in ["neu","cha"] and dw[:2]==br[3:5]) ):
					pltctr+=1
					SinglePlot(pltctr,br,brix,0,0,"linear", dw,dwix,0,0, "log", Label, Color,
						Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "DC", "BR")
	
		pltctr+=1
		SinglePlot(pltctr,"s1mass",24,S1MMIN,S1MMAX,"log", "GamHtot",133,0,0,"log",
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "DC","pdw")
		pltctr+=1
		SinglePlot(pltctr,"s1mass",24,S1MMIN,S1MMAX,"log", "GamHhadr",139,1E-24,1E-4,"log",
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "DC","pdw hadr")
	
		print(Time(),"Starting s1 pdws by channel") ############ S1 BY CHANNEL
		for i,(br,brix) in enumerate(s1_brs):	#### INDIVIDUALLY, THEN STACKED
			pltctr+=1
			plt.figure(pltctr)
			plt.grid(visible=True, which="both",axis="x",alpha=.5,lw=1)
			plt.scatter( 	[r[24] for r in master_list[-1]],
					[r[brix]*r[133] for r in master_list[-1]],
					alpha=1,s=1,marker=',',linewidths=0)
			plt.title(save_dir_name+" : s1 {} pdw v s1mass".format(br[6:]))
			plt.xlabel("s1mass")
			plt.ylabel("GamH{}".format(br[6:]))
			plt.yscale("log")
			plt.xscale("log")
			plt.savefig("{}/{}/DC/pdw/s1 {} pdw v s1mass.png".format(DIR,save_dir_name,br[6:]),dpi=DPI)
			plt.xlim(0.1,10)
			plt.ylim(1E-24,1E-4)
			plt.savefig("{}/{}/DC/pdw/s1 {} pdw v s1mass LSAF match.png".format(DIR,save_dir_name,br[6:]),dpi=DPI)
			plt.close()		

		print(Time(),"stacking s1 pdws")
		pltctr+=1				#### STACKED
		plt.figure(pltctr)
		plt.grid(visible=True, which="both",axis="x",alpha=.5,lw=1)
		for i,(br,brix) in enumerate(s1_brs):
			plt.scatter( 	[r[24] for r in master_list[-1]],
					[r[brix]*r[133] for r in master_list[-1]],
					label="GamH{}".format(br[6:]), alpha=1,s=1,marker=',',linewidths=0)
	#	plt.scatter(	[r[24] for r in master_list[-1]],
	#			[ sum([r[brix]*r[133] for (br,brix) in s1_brs]) for r in master_list[-1]],
	#			label="dw sum (calc.d)", color="black", alpha=1, s=4, marker="_")
	#	plt.scatter(	[r[24] for r in master_list[-1]],
	#			[r[133] for r in master_list[-1]],
	#			label="dw tot (NMSSMT)", color="darkgray", alpha=1, s=4, marker="_")
		plt.title(save_dir_name+" : s1 pdw per channel v s1mass")
		plt.xlabel("s1mass")
		plt.ylabel("pdw per channel")
		plt.yscale("log")
		plt.xscale("log")
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=int((len(s1_brs)/2+1)),columnspacing=0.5,frameon=False,fontsize=6)
		for x in range(len(s1_brs)): leg.legend_handles[x]._sizes = [10]
		plt.savefig("{}/{}/DC/s1 pdws.png".format(DIR,save_dir_name),dpi=DPI)
		plt.xlim(S1MMIN,S1MMAX)
		plt.savefig("{}/{}/DC/s1 pdws zoom.png".format(DIR,save_dir_name),dpi=DPI)
		plt.ylim(1E-24,1E-4)
		plt.xlim(0.1,10)
		plt.savefig("{}/{}/DC/s1 pdws LSAF match.png".format(DIR,save_dir_name),dpi=DPI)
		plt.close()

		print(Time(),"Starting p1 pdws by channel") ############ P1 BY CHANNEL
		for i,(br,brix) in enumerate(p1_brs):	#### INDIVIDUALLY, THEN STACKED
			pltctr+=1
			plt.figure(pltctr)
			plt.grid(visible=True, which="both",axis="x",alpha=.5,lw=1)
			plt.scatter( 	[r[36] for r in master_list[-1]],
					[r[brix]*r[135] for r in master_list[-1]],
					alpha=1,s=1,marker=',',linewidths=0)
			plt.title(save_dir_name+" : p1 {} pdw v p1mass".format(br[6:]))
			plt.xlabel("p1mass")
			plt.ylabel("GamA{}".format(br[6:]))
			plt.yscale("log")
			plt.xscale("log")
			plt.savefig("{}/{}/DC/pdw/p1 {} pdw v p1mass.png".format(DIR,save_dir_name,br[6:]),dpi=DPI)
			plt.xlim(0.1,10)
			plt.ylim(1E-24,1E-2)
			plt.savefig("{}/{}/DC/pdw/p1 {} pdw v p1mass LSAF match.png".format(DIR,save_dir_name,br[6:]),dpi=DPI)
			plt.close()
		print(Time(),"stacking p1 pdws")
		pltctr+=1				#### START STACKED
		plt.figure(pltctr)
		plt.grid(visible=True, which="both",axis="x",alpha=.5,lw=1)
		for i,(br,brix) in enumerate(p1_brs):
			plt.scatter( 	[r[36] for r in master_list[-1]],
					[r[brix]*r[135] for r in master_list[-1]],
					label="GamA{}".format(br[6:]), alpha=1,s=1,marker=',',linewidths=0)
	#	plt.scatter(	[r[36] for r in master_list[-1]],
	#			[ sum([r[brix]*r[135] for (br,brix) in p1_brs]) for r in master_list[-1]],
	#			label="dw sum (calc.d)", color="black", alpha=1, s=4, marker="_")
	#	plt.scatter(	[r[36] for r in master_list[-1]],
	#			[r[135] for r in master_list[-1]],
	#			label="dw tot (NMSSMT)", color="darkgray", alpha=1, s=4, marker="_")
		plt.title(save_dir_name+" : p1 pdw per channel v p1mass")
		plt.xlabel("p1mass")
		plt.ylabel("pdw per channel")
		plt.yscale("log")
		plt.xscale("log")
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=int((len(p1_brs)/2+1)),columnspacing=0.6,frameon=False,fontsize=6)
		for x in range(len(p1_brs)): leg.legend_handles[x]._sizes = [10]
		plt.savefig("{}/{}/DC/p1 pdws.png".format(DIR,save_dir_name),dpi=DPI)
		plt.xlim(P1MMIN,P1MMAX)
		plt.savefig("{}/{}/DC/p1 pdws zoom.png".format(DIR,save_dir_name),dpi=DPI)
		plt.ylim(1E-24,1E-2)
		plt.xlim(0.1,10)
		plt.savefig("{}/{}/DC/p1 pdws LSAF match.png".format(DIR,save_dir_name),dpi=DPI)
		plt.close()

		print(Time(),"s1 hadronic pdws")	######## STARTING S1 HADRONIC BREAKDOWN
		pltctr+=1
		SinglePlot(pltctr, "s1mass", 24, S1MMIN,S1MMAX, "log","GamHhadr", 139, 0, 0, "log",
			Label, Color,Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "DC", "pdw hadr")

		for (gamH,gix) in Hpdc_list:		#### INDIV. S1 HADRONIC
			pltctr+=1
			plt.figure(pltctr)
			plt.grid(visible=True, which="both",axis="x",alpha=.5,lw=1)
			plt.scatter(	[r[24] for r in master_list[-1]],
					[r[139] for r in master_list[-1]],
					color="b",alpha=1,s=4,marker='_',label="GamHhadr")
			if gamH != "GamHhadr":
				if 141<=gix and gix<=149:
					gamHfiltered=[max(min( 4-2*r[24],1),0)*r[gix] for r in master_list[-1]]
				elif 150<=gix and gix<=151:
					gamHfiltered=[min(max(-3+2*r[24],0),1)*r[gix] for r in master_list[-1]]
				elif 152<=gix and gix<=155:
					gamHfiltered=[r[gix] for r in master_list[-1]]
				plt.scatter( 	[r[24] for r in master_list[-1]], gamHfiltered,
					color="r",alpha=1,s=3,marker=',',linewidths=0,label=gamH)
				leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=2,columnspacing=0.6,frameon=False,fontsize=6)
				for x in range(2): leg.legend_handles[x]._sizes = [10]

			plt.title(save_dir_name+" : {} v s1mass".format(gamH))
			plt.xlabel("s1mass")
			plt.ylabel(gamH)
			plt.yscale("log")
			plt.xscale("log")
			plt.savefig("{}/{}/DC/pdw hadr/{} v s1mass.png".format(DIR,save_dir_name,gamH),dpi=DPI)
			plt.xlim(S1MMIN,S1MMAX)
			plt.savefig("{}/{}/DC/pdw hadr/{} v s1mass zoom..png".format(DIR,save_dir_name,gamH),dpi=DPI)
			plt.xlim(0.1,10)
			plt.ylim(1E-24,1E-4)
			plt.savefig("{}/{}/DC/pdw hadr/{} v s1mass LSAF match.png".format(DIR,save_dir_name,gamH),dpi=DPI)
			plt.close()
			

		pltctr+=1	####### INDIVIDUALLY SUMMED PARTS TO hadr FOR S1
		plt.figure(pltctr)
		plt.grid(visible=True, which="both",axis="x",alpha=.5,lw=1)
		plt.scatter(	[r[24] for r in master_list[-1]],
				[r[139] for r in master_list[-1]],
				color="b",alpha=.4,s=3,marker='_',label="GamHhadr (NMSSMT)")
		plt.scatter(	[r[24] for r in master_list[-1]],
				[r[107]*r[133] for r in master_list[-1]],
				color="g",alpha=.4,s=3,marker='_',label="GamHhadr (BR*DC)")
		plt.scatter( 	[r[24] for r in master_list[-1]],
				[	max( min(4-2*r[24],1)  ,0)*sum(r[141:149+1]) +
					min( max(-3+2*r[24],0),1)*sum(r[150:151+1]) +
					sum(r[152:155+1])	for r in master_list[-1] ],
				color="r",alpha=.4,s=3,marker=',',linewidths=0,label="GamHhadr (calc.d)")
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=3,columnspacing=0.5,frameon=False,fontsize=6)
		for x in range(3): leg.legend_handles[x]._sizes = [10]
		plt.title(save_dir_name+" : GamHhadr calc.d/NMSSMT v s1mass")
		plt.xlabel("s1mass")
		plt.ylabel("GamHhadr")
		plt.yscale("log")
		plt.xscale("log")
		plt.savefig("{}/{}/DC/GamHhadr calc v given.png".format(DIR,save_dir_name),dpi=DPI)
		plt.xlim(S1MMIN,S1MMAX)
		plt.savefig("{}/{}/DC/GamHhadr calc v given zoom.png".format(DIR,save_dir_name),dpi=DPI)
		plt.xlim(0.1,10)
		plt.ylim(1E-24,1E-4)
		plt.savefig("{}/{}/DC/GamHhadr calc v given LSAF match.png".format(DIR,save_dir_name),dpi=DPI)
		plt.close()



		pltctr+=1
		plt.figure(pltctr)			#### STACKED S1 HADRONIC PARTIAL DECAY WIDTHS
		plt.grid(visible=True, which="both",axis="x",alpha=.5,lw=1)
		for (gamH,gix) in Hpdc_list:
			if (gix >= 141 and gix <= 149):
				plt.scatter(	[r[24] for r in master_list[-1]],
					[ r[gix]*(1- (r[24]-1.5)/abs(r[24]-1.5))/2 for r in master_list[-1]],	#crop at 1.5G
				#	[ r[gix]*min(max((1-(r[24]-1.5)/(.5)),0),1) for r in master_list[-1]],	#blend thru 1.5-2G
					label=gamH[4:], alpha=1,s=4,marker='>',linewidths=0)
			if (gix >= 150 and gix <= 151):
				plt.scatter(	[r[24] for r in master_list[-1]],
					[ r[gix]*(1 + (r[24]-2)/abs(r[24]-2))/2 for r in master_list[-1]],	#crop at 1.5G
				#	[ r[gix]*min(max(((r[24]-1.5)/(.5)),0),1) for r in master_list[-1]],	#blend thru 1.5-2G
					label=gamH[4:], alpha=1,s=4,marker='<',linewidths=0)
		#	if (gix == 152 or gix == 154):
			if (gix >= 152):
				plt.scatter([r[24] for r in master_list[-1]], [r[gix] for r in master_list[-1]],
					label=gamH[4:], alpha=1,s=1,marker=',',linewidths=0)
		plt.scatter([r[24] for r in master_list[-1]], [r[139] for r in master_list[-1]],
				label="GamHhadr", color="black",alpha=.5,s=4,marker="_")
		plt.title(save_dir_name+" : s1 hadronic pdws v s1mass")
		plt.xlabel("s1mass")
		plt.ylabel("s1 hadronic pdws")
		plt.yscale("log")
		plt.xscale("log")
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=8,columnspacing=0.6,frameon=False,fontsize=6)
		for x in range(len(Hpdc_list)): leg.legend_handles[x]._sizes = [10]
		plt.savefig("{}/{}/DC/s1 hadr pdws.png".format(DIR,save_dir_name),dpi=DPI)
		plt.xlim(S1MMIN,S1MMAX)
		plt.savefig("{}/{}/DC/s1 hadr pdws zoom.png".format(DIR,save_dir_name),dpi=DPI)
		plt.xlim(0.250,4)
		plt.xscale("linear")
		plt.savefig("{}/{}/DC/s1 hadr pdws 2G.png".format(DIR,save_dir_name),dpi=DPI)
		plt.close()

		print(Time(),"p1 hadronic pdws")	######## STARTING P1 HADRONIC BREAKDOWN
		pltctr+=1
		SinglePlot(pltctr, "p1mass", 36, P1MMIN,P1MMAX, "log","GamAhadr", 140, 0, 0, "log",
			Label, Color,Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "DC", "pdw hadr")
		for (gamA,gix) in Apdc_list:		#### INDIV. P1 HADRONIC
			pltctr+=1
			plt.figure(pltctr)
			plt.grid(visible=True, which="both",axis="x",alpha=.5,lw=1)
			plt.scatter(	[r[36] for r in master_list[-1]],
					[r[140] for r in master_list[-1]],
					color="b",alpha=1,s=4,marker='_',label="GamAhadr")
			if gamA != "GamAhadr":
				if 156<=gix and gix<=176:
					gamAfiltered=[max(min( 4-r[36],1),0)*r[gix] for r in master_list[-1]]
				elif 177<=gix and gix<=178:
					gamAfiltered=[min(max(-3+r[36],0),1)*r[gix] for r in master_list[-1]]
				elif 179<=gix and gix<=182:
					gamAfiltered=[r[gix] for r in master_list[-1]]
				plt.scatter( 	[r[36] for r in master_list[-1]], gamAfiltered,
						color="r",alpha=1,s=3,marker=',',linewidths=0,label=gamA)
				leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=2,columnspacing=0.6,frameon=False,fontsize=6)
				for x in range(2): leg.legend_handles[x]._sizes = [10]
			plt.title(save_dir_name+" : {} v p1mass".format(gamA))
			plt.xlabel("p1mass")
			plt.ylabel(gamA)
			plt.yscale("log")
			plt.xscale("log")
			plt.savefig("{}/{}/DC/pdw hadr/{} v p1mass.png".format(DIR,save_dir_name,gamA),dpi=DPI)
			plt.xlim(P1MMIN,P1MMAX)
			plt.savefig("{}/{}/DC/pdw hadr/{} v p1mass zoom..png".format(DIR,save_dir_name,gamA),dpi=DPI)
			plt.xlim(0.1,10)
			plt.ylim(1E-24,1E-2)
			plt.savefig("{}/{}/DC/pdw hadr/{} v p1mass LSAF match.png".format(DIR,save_dir_name,gamA),dpi=DPI)
			plt.close()
	
		pltctr+=1	####### INDIVIDUALLY SUMMED PARTS TO hadr FOR P1
		plt.figure(pltctr)
		plt.grid(visible=True, which="both",axis="x",alpha=.5,lw=1)
		plt.scatter(	[r[36] for r in master_list[-1]],
				[r[140] for r in master_list[-1]],
				color="b",alpha=.4,s=3,marker='_',label="GamAhadr (NMSSMT)")
		plt.scatter(	[r[36] for r in master_list[-1]],
				[r[108]*r[135] for r in master_list[-1]],
				color="g",alpha=.4,s=3,marker='_',label="GamAhadr (BR*DW)")
		plt.scatter( 	[r[36] for r in master_list[-1]],
				[	max( min(4-r[36],1)  ,0)*sum(r[156:176+1]) +
					min( max(-3+r[36],0),1)*sum(r[177:178+1]) +
					sum(r[179:182+1])	for r in master_list[-1] ],
				color="r",alpha=.4,s=3,marker=',',linewidths=0,label="GamAhadr (calc.d)")
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=3,columnspacing=0.5,frameon=False,fontsize=6)
		for x in range(3): leg.legend_handles[x]._sizes = [10]
		plt.title(save_dir_name+" : GamAhadr calc.d/NMSSMT v p1mass")
		plt.xlabel("p1mass")
		plt.ylabel("GamAhadr")
		plt.yscale("log")
		plt.xscale("log")
		plt.savefig("{}/{}/DC/GamAhadr calc v given.png".format(DIR,save_dir_name),dpi=DPI)
		plt.xlim(P1MMIN,P1MMAX)
		plt.savefig("{}/{}/DC/GamAhadr calc v given zoom.png".format(DIR,save_dir_name),dpi=DPI)
		plt.xlim(0.1,10)
		plt.ylim(1E-24,1E-2)
		plt.savefig("{}/{}/DC/GamAhadr calc v given LSAF match.png".format(DIR,save_dir_name),dpi=DPI)
		plt.close()
	
		pltctr+=1
		plt.figure(pltctr)			#### STACKED P1 HADRONIC PARTIAL DECAY WIDTHS
		plt.grid(visible=True, which="both",axis="x",alpha=.5,lw=1)
		for (gamA,gix) in Apdc_list:
			if (gix >= 156 and gix <= 176):
				plt.scatter(	[r[36] for r in master_list[-1]],
					[ r[gix]*(1-(r[36]-3.0001)/abs(r[36]-3.0001))/2 for r in master_list[-1]],	#crop at 3G
				#	[ r[gix]*min(max((4-r[36]),0),1) for r in master_list[-1]],	#blend thru 3-4G
					label=gamA[4:], alpha=1,s=4,marker='>',linewidths=0)
			if (gix >= 177 and gix <= 178):
				plt.scatter(	[r[36] for r in master_list[-1]],
					[ r[gix]*(1+(r[36]-4)/abs(r[36]-4))/2 for r in master_list[-1]],	#crop at 4G
				#	[ r[gix]*min(max((r[36]-3),0),1) for r in master_list[-1]],	#blend thru 3-4G
					label=gamA[4:], alpha=1,s=4,marker='<',linewidths=0)
		#	if (gix == 179 or gix == 181):
			if (gix >= 179):
				plt.scatter([r[36] for r in master_list[-1]], [r[gix] for r in master_list[-1]],
					label=gamA[4:], alpha=1,s=1,marker=',',linewidths=0)
		plt.scatter([r[36] for r in master_list[-1]], [r[140] for r in master_list[-1]],
				label="GamAhadr", color="black",alpha=.5,s=4,marker="_")
		plt.title(save_dir_name+" : p1 hadronic pdws v p1mass")
		plt.xlabel("p1mass")
		plt.ylabel("p1 hadronic pdws")
		plt.yscale("log")
		plt.xscale("log")
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=10,columnspacing=0.6,frameon=False,fontsize=6)
		for x in range(len(Apdc_list)): leg.legend_handles[x]._sizes = [10]
		plt.savefig("{}/{}/DC/p1 hadr pdws.png".format(DIR,save_dir_name),dpi=DPI)
		plt.xlim(P1MMIN,P1MMAX)
		plt.savefig("{}/{}/DC/p1 hadr pdws zoom.png".format(DIR,save_dir_name),dpi=DPI)		
		plt.xlim(0.400,5)
		plt.xscale("linear")
		plt.savefig("{}/{}/DC/p1 hadr pdws 3G.png".format(DIR,save_dir_name),dpi=DPI)
		plt.close()		

	if DO_DC and True:
		########### START A(2HDM) CONVERTED DECAY WIDTHS ##############

		print(Time(),"Starting A(2HDM) pdws by channel") ########## A2HDM BY CHANNEL
		for i,(br,brix) in enumerate(p1_brs):	#### INDIVIDUALLY, THEN STACKED
			pltctr+=1
			plt.figure(pltctr)
			plt.grid(visible=True, which="both",axis="x",alpha=.5,lw=1)
			plt.scatter( 	[r[36] for r in master_list[-1]],
					[r[brix]*r[135]/(r[37]**2) for r in master_list[-1]],
					alpha=1,s=1,marker=',',linewidths=0)
			plt.title(save_dir_name+" : A(2HDM) {} pdw v p1mass".format(br[6:]))
			plt.xlabel("p1mass")
			plt.ylabel("GamA(2HDM){}".format(br[6:]))
			plt.yscale("log")
			plt.xscale("log")
			plt.savefig("{}/{}/DC/pdw/A(2HDM) {} pdw v p1mass.png".format(DIR,save_dir_name,br[6:]),dpi=DPI)
			plt.xlim(0.1,10)
			plt.ylim(1E-24,1E-2)
			plt.savefig("{}/{}/DC/pdw/A(2HDM) {} pdw v p1mass LSAF match.png".format(DIR,save_dir_name,br[6:]),dpi=DPI)
			plt.close()

		for (gamA,gix) in Apdc_list:		#### INDIVIDUAL A(2HDM) HADRONICS
			pltctr+=1
			plt.figure(pltctr)
			plt.grid(visible=True, which="both",axis="x",alpha=.5,lw=1)
			plt.scatter(	[r[36] for r in master_list[-1]],
					[r[140]/(r[37]**2) for r in master_list[-1]],
					color="b",alpha=1,s=4,marker='_',label="GamA(2HDM)hadr")
			if gamA != "GamAhadr":
				if 156<=gix and gix<=176:
					gamAfiltered=[max(min( 4-r[36],1),0)*r[gix]/(r[37]**2) for r in master_list[-1]]
				elif 177<=gix and gix<=178:
					gamAfiltered=[min(max(-3+r[36],0),1)*r[gix]/(r[37]**2) for r in master_list[-1]]
				elif 179<=gix and gix<=182:
					gamAfiltered=[r[gix]/(r[37]**2) for r in master_list[-1]]
				plt.scatter( 	[r[36] for r in master_list[-1]], gamAfiltered,
						color="r",alpha=1,s=3,marker=',',linewidths=0,label=gamA)
				leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=2,columnspacing=0.6,frameon=False,fontsize=6)
				for x in range(2): leg.legend_handles[x]._sizes = [10]
			plt.title(save_dir_name+" : GamA(2HDM){} v p1mass".format(gamA[4:]))
			plt.xlabel("p1mass")
			plt.ylabel("GamA(2HDM){}".format(gamA[4:]))
			plt.yscale("log")
			plt.xscale("log")
			plt.savefig("{}/{}/DC/pdw hadr/GamA(2HDM){} v p1mass.png".format(DIR,save_dir_name,gamA[4:]),dpi=DPI)
			plt.xlim(P1MMIN,P1MMAX)
			plt.savefig("{}/{}/DC/pdw hadr/GamA(2HDM){} v p1mass zoom..png".format(DIR,save_dir_name,gamA[4:]),dpi=DPI)
			plt.xlim(0.1,10)
			plt.ylim(1E-24,1E-2)
			plt.savefig("{}/{}/DC/pdw hadr/GamA(2HDM){} v p1mass LSAF match.png".format(DIR,save_dir_name,gamA[4:]),dpi=DPI)
			plt.close()

		###########  END  A(2HDM) CONVERTED DECAY WIDTHS ##############


	if DO_COUP:
		print(Time(),"Beginning XI plots")
		for (coup,cix) in coup_list:
			if coup[2:4] not in ["p1"]: continue # temp skip to spam s1/p1
			print(Time(),"Evaluating {}...".format(coup))
			for (par,pix) in par_list:	
		
				continue # skip the param plots
		
				pltctr+=1
				SinglePlot(pltctr, par, pix, 0, 0,"linear",
						coup, cix, 0, 0,"symlog",
						Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "XI","Parameter")
				for (mass,mix,mmin,mmax) in [("s1mass",24,S1MMIN,S1MMAX),("p1mass",36,P1MMIN,P1MMAX)]:#,("neu1mass",44,0,0)]:		
					pltctr+=1
					HeatPlot(pltctr, mass, mix,"turbo",
						par, pix, 0, 0,"linear",
						coup, cix, 0, 0, "symlog", Size, DPI, "Heatmap","XI v par c mass")
					pltctr+=1
					HeatPlot(pltctr, par,pix,"turbo",
						mass,mix, mmin, mmax, "log",
						coup, cix, 0, 0, "symlog", Size, DPI, "Heatmap","XI v mass c par")

			for (mass,mix,mmin,mmax) in [("s1mass",24,S1MMIN,S1MMAX),("p1mass",36,P1MMIN,P1MMAX)]:#,("neu1mass",44)]:
				if coup[2:4] not in mass: continue
				pltctr+=1
				SinglePlot(pltctr, mass, mix, mmin, mmax, "log", coup,cix, 0, 0, "symlog",
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "XI","Mass")

	if DO_MISC:
		print(Time(),"Initiating miscellany.")
		pltctr+=1
		SinglePlot(pltctr,"s1mass",24,S1MMIN,S1MMAX,"log", "s1scomp",27,0,0,"linear",
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "BR","pdw")
		pltctr+=1
		SinglePlot(pltctr,"p1mass",36,P1MMIN,P1MMAX,"log", "p1Acomp",37,0,0,"linear",
					Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "BR","pdw hadr")
		
							########## s1, p1 dws are in 133 and 135 ########
		print(Time(),"s1 hadronic BRs")		######## STARTING S1 HADRONIC BREAKDOWN	#########
		for (gamH,gix) in Hpdc_list:		#### INDIV. S1 HADRONIC ### CONVERTING TO BRS ###
			pltctr+=1			#################################################
			plt.figure(pltctr)
			plt.grid(visible=True, which="both",axis="x",alpha=.5,lw=1)
			plt.scatter(	[r[24] for r in master_list[-1]],
					[r[139]/r[133] for r in master_list[-1]],
					color="b",alpha=1,s=4,marker='_',label="br_s1_hadr (calc.d)")
			if gamH != "GamHhadr":
				if 141<=gix and gix<=149:
					brHfiltered=[max(min( 4-2*r[24],1),0)*r[gix]/r[133] for r in master_list[-1]]
				elif 150<=gix and gix<=151:
					brHfiltered=[min(max(-3+2*r[24],0),1)*r[gix]/r[133] for r in master_list[-1]]
				elif 152<=gix and gix<=155:
					brHfiltered=[r[gix]/r[133] for r in master_list[-1]]
				plt.scatter( 	[r[24] for r in master_list[-1]], brHfiltered,
					color="r",alpha=1,s=3,marker=',',linewidths=0,label="br_s1_"+gamH[4:])
				leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=2,columnspacing=0.5,frameon=False,fontsize=6)
				for x in range(2): leg.legend_handles[x]._sizes = [10]

			plt.title(save_dir_name+" : br_s1_{} v s1mass".format(gamH[4:]))
			plt.xlabel("s1mass")
			plt.ylabel("br_s1_"+gamH[4:])
			plt.yscale("linear")
			plt.ylim(-0.02,1.02)
			plt.xscale("log")
			plt.savefig("{}/{}/BR/pdw hadr/br_s1_{} v s1mass.png".format(DIR,save_dir_name,gamH[4:]),dpi=DPI)
			plt.xlim(S1MMIN,S1MMAX)
			plt.ylim(1E-9,1.5)
			plt.yscale("log")
			plt.savefig("{}/{}/BR/pdw hadr/br_s1_{} v s1mass zoom.png".format(DIR,save_dir_name,gamH[4:]),dpi=DPI)
			plt.xlim(0.1,10)
			plt.ylim(.005,1.5)
			plt.savefig("{}/{}/BR/pdw hadr/br_s1_{} v s1mass LSAF match.png".format(DIR,save_dir_name,gamH[4:]),dpi=DPI)
			plt.close()
			

		pltctr+=1	####### INDIVIDUALLY SUMMED PARTS TO hadr FOR S1
		plt.figure(pltctr)
		plt.grid(visible=True, which="both",axis="x",alpha=.5,lw=1)
		plt.scatter(	[r[24] for r in master_list[-1]],
				[r[107] for r in master_list[-1]],
				color="b",alpha=1,s=4,marker='_',label="br_s1_hadr (NMSSMT)")
		plt.scatter( 	[r[24] for r in master_list[-1]],
				[	(max( min(4-2*r[24],1)  ,0)*sum(r[141:149+1]) +
					min( max(-3+2*r[24],0),1)*sum(r[150:151+1]) +
					sum(r[152:155+1]))/r[133]	for r in master_list[-1] ],
				color="r",alpha=1,s=3,marker=',',linewidths=0,label="br_s1_hadr (calc.d)")
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=2,columnspacing=0.5,frameon=False,fontsize=6)
		for x in range(2): leg.legend_handles[x]._sizes = [10]
		plt.title(save_dir_name+" : br_s1_hadr calc.d/NMSSMT v s1mass")
		plt.xlabel("s1mass")
		plt.ylabel("br_s1_hadr")
		plt.ylim(-0.02,1.02)
		plt.xscale("log")
		plt.savefig("{}/{}/BR/br_s1_hadr calc v given.png".format(DIR,save_dir_name),dpi=DPI)
		plt.xlim(S1MMIN,S1MMAX)
		plt.yscale("log")
		plt.ylim(1E-9,1.5)
		plt.savefig("{}/{}/BR/br_s1_hadr calc v given zoom.png".format(DIR,save_dir_name),dpi=DPI)
		plt.xlim(0.1,10)
		plt.savefig("{}/{}/BR/br_s1_hadr calc v given LSAF match.png".format(DIR,save_dir_name),dpi=DPI)
		plt.close()



		pltctr+=1
		plt.figure(pltctr)			#### STACKED S1 HADRONIC PARTIAL DECAY WIDTHS
		plt.grid(visible=True, which="both",axis="x",alpha=.5,lw=1)
		for (gamH,gix) in Hpdc_list:
			if (gix >= 141 and gix <= 149):
				plt.scatter(	[r[24] for r in master_list[-1]],
					[ r[gix]*(1- (r[24]-1.5)/abs(r[24]-1.5))/2/r[133] for r in master_list[-1]],	#crop at 1.5G
				#	[ r[gix]*min(max((1-(r[24]-1.5)/(.5)),0),1) for r in master_list[-1]],	#blend thru 1.5-2G
					label="br_s1_"+gamH[4:], alpha=1,s=4,marker='>',linewidths=0)
			if (gix >= 150 and gix <= 151):
				plt.scatter(	[r[24] for r in master_list[-1]],
					[ r[gix]/r[133]*(1 + (r[24]-2)/abs(r[24]-2))/2 for r in master_list[-1]],	#crop at 1.5G
				#	[ r[gix]*min(max(((r[24]-1.5)/(.5)),0),1) for r in master_list[-1]],	#blend thru 1.5-2G
					label="br_s1_"+gamH[4:], alpha=1,s=4,marker='<',linewidths=0)
		#	if (gix == 152 or gix == 154):
			if (gix >= 152):
				plt.scatter([r[24] for r in master_list[-1]], [r[gix]/r[133] for r in master_list[-1]],
					label="br_s1_"+gamH[4:], alpha=1,s=1,marker=',',linewidths=0)
		plt.scatter([r[24] for r in master_list[-1]], [r[139]/r[133] for r in master_list[-1]],
				label="br_s1_hadr", color="black",alpha=.5,s=4,marker="_")
		plt.title(save_dir_name+" : s1 hadronic BRs v s1mass")
		plt.xlabel("s1mass")
		plt.ylabel("s1 hadronic BRs")
		plt.xscale("log")
		plt.ylim(-0.02,1.02)
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=8,columnspacing=0.6,frameon=False,fontsize=6)
		for x in range(len(Hpdc_list)): leg.legend_handles[x]._sizes = [10]
		plt.savefig("{}/{}/BR/s1 hadr BRs.png".format(DIR,save_dir_name),dpi=DPI)
		plt.xlim(S1MMIN,S1MMAX)
		plt.yscale("log")
		plt.ylim(1E-9,1.5)
		plt.savefig("{}/{}/BR/s1 hadr BRs zoom.png".format(DIR,save_dir_name),dpi=DPI)
		plt.xscale("linear")
		plt.savefig("{}/{}/BR/s1 hadr BRs 2G.png".format(DIR,save_dir_name),dpi=DPI)
		plt.close()

		print(Time(),"p1 hadronic BRs")		######## STARTING P1 HADRONIC BREAKDOWN ##########
		for (gamA,gix) in Apdc_list:		#### INDIV. P1 HADRONIC ### SAME AS ABOVE S1 STYLE
			pltctr+=1
			plt.figure(pltctr)
			plt.grid(visible=True, which="both",axis="x",alpha=.5,lw=1)
			plt.scatter(	[r[36] for r in master_list[-1]],
					[r[140]/r[135] for r in master_list[-1]],
					color="b",alpha=1,s=4,marker='_',label="br_p1_hadr (calc.d)")
			if gamA != "GamAhadr":
				if 156<=gix and gix<=176:
					brAfiltered=[max(min( 4-r[36],1),0)*r[gix]/r[135] for r in master_list[-1]]
				elif 177<=gix and gix<=178:
					brAfiltered=[min(max(-3+r[36],0),1)*r[gix]/r[135] for r in master_list[-1]]
				elif 179<=gix and gix<=182:
					brAfiltered=[r[gix]/r[135] for r in master_list[-1]]
				plt.scatter( 	[r[36] for r in master_list[-1]], brAfiltered,
						color="r",alpha=1,s=3,marker=',',linewidths=0,label="br_p1_"+gamA[4:])
				leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=2,columnspacing=0.6,frameon=False,fontsize=6)
				for x in range(2): leg.legend_handles[x]._sizes = [10]
			plt.title(save_dir_name+" : br_p1_{} v p1mass".format(gamA[4:]))
			plt.xlabel("p1mass")
			plt.ylabel("br_p1_"+gamA[4:])
			plt.ylim(-0.02,1.02)
			plt.xscale("log")
			plt.savefig("{}/{}/BR/pdw hadr/br_p1_{} v p1mass.png".format(DIR,save_dir_name,gamA[4:]),dpi=DPI)
			plt.xlim(P1MMIN,P1MMAX)
			plt.yscale("log")
			plt.ylim(1E-9,2)
			plt.savefig("{}/{}/BR/pdw hadr/br_p1_{} v p1mass zoom..png".format(DIR,save_dir_name,gamA[4:]),dpi=DPI)
			plt.xlim(0.1,10)
			plt.ylim(.005,2)
			plt.savefig("{}/{}/BR/pdw hadr/br_p1_{} v p1mass LSAF match.png".format(DIR,save_dir_name,gamA[4:]),dpi=DPI)
			plt.close()
	
		pltctr+=1	####### INDIVIDUALLY SUMMED PARTS TO hadr FOR P1
		plt.figure(pltctr)
		plt.grid(visible=True, which="both",axis="x",alpha=.5,lw=1)
		plt.scatter(	[r[36] for r in master_list[-1]],
				[r[108] for r in master_list[-1]],
				color="b",alpha=1,s=4,marker='_',label="br_p1_hadr (NMSSMT)")
		plt.scatter( 	[r[36] for r in master_list[-1]],
				[	(max( min(4-r[36],1)  ,0)*sum(r[156:176+1]) +
					min( max(-3+r[36],0),1)*sum(r[177:178+1]) +
					sum(r[179:182+1]))/r[135]	for r in master_list[-1] ],
				color="r",alpha=1,s=3,marker=',',linewidths=0,label="br_p1_hadr (calc.d)")
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=2,columnspacing=0.6,frameon=False,fontsize=6)
		for x in range(2): leg.legend_handles[x]._sizes = [10]
		plt.title(save_dir_name+" : br_p1_hadr calc.d/NMSSMT v p1mass")
		plt.xlabel("p1mass")
		plt.ylabel("br_p1_hadr")
		plt.ylim(-0.02,1.02)
		plt.xscale("log")
		plt.savefig("{}/{}/BR/br_p1_hadr calc v given.png".format(DIR,save_dir_name),dpi=DPI)
		plt.xlim(P1MMIN,P1MMAX)
		plt.yscale("log")
		plt.ylim(1E-9,2)
		plt.savefig("{}/{}/BR/br_p1_hadr calc v given zoom.png".format(DIR,save_dir_name),dpi=DPI)
		plt.xlim(0.1,10)
		plt.savefig("{}/{}/BR/br_p1_hadr calc v given LSAF match.png".format(DIR,save_dir_name),dpi=DPI)
		plt.close()
	
		pltctr+=1
		plt.figure(pltctr)			#### STACKED P1 HADRONIC PARTIAL DECAY WIDTHS
		plt.grid(visible=True, which="both",axis="x",alpha=.5,lw=1)
		for (gamA,gix) in Apdc_list:
			if (gix >= 156 and gix <= 176):
				plt.scatter(	[r[36] for r in master_list[-1]],
					[ r[gix]/r[135]*(1-(r[36]-3)/abs(r[36]-3))/2 for r in master_list[-1]],	#crop at 3G
				#	[ r[gix]*min(max((4-r[36]),0),1) for r in master_list[-1]],	#blend thru 3-4G
					label="br_p1_"+gamA[4:], alpha=1,s=4,marker='>',linewidths=0)
			if (gix >= 177 and gix <= 178):
				plt.scatter(	[r[36] for r in master_list[-1]],
					[ r[gix]/r[135]*(1+(r[36]-4)/abs(r[36]-4))/2 for r in master_list[-1]],	#crop at 4G
				#	[ r[gix]*min(max((r[36]-3),0),1) for r in master_list[-1]],	#blend thru 3-4G
					label="br_p1_"+gamA[4:], alpha=1,s=4,marker='<',linewidths=0)
		#	if (gix == 179 or gix == 181):
			if (gix >= 179):
				plt.scatter([r[36] for r in master_list[-1]], [r[gix] for r in master_list[-1]],
					label="br_p1_"+gamA[4:], alpha=1,s=1,marker=',',linewidths=0)
		plt.scatter([r[36] for r in master_list[-1]], [r[108] for r in master_list[-1]],
				label="br_p1_hadr (NMSSMT)", color="black",alpha=.5,s=4,marker="_")
		plt.title(save_dir_name+" : p1 hadronic BRs v p1mass")
		plt.xlabel("p1mass")
		plt.ylabel("p1 hadronic BRs")
		plt.xscale("log")
		plt.ylim(-.02,1.02)
		leg = plt.legend(loc=LOC,bbox_to_anchor=BBOX_TO_ANCHOR,ncols=10,columnspacing=0.6,frameon=False,fontsize=6)
		for x in range(len(Apdc_list)): leg.legend_handles[x]._sizes = [10]
		plt.savefig("{}/{}/BR/p1 hadr BRs.png".format(DIR,save_dir_name),dpi=DPI)
		plt.xlim(P1MMIN,P1MMAX)
		plt.yscale("log")
		plt.ylim(1E-9,2)
		plt.savefig("{}/{}/BR/p1 hadr BRs zoom.png".format(DIR,save_dir_name),dpi=DPI)		
		plt.xlim(0.400,5)
		plt.xscale("linear")
		plt.savefig("{}/{}/BR/p1 hadr BRs 3G.png".format(DIR,save_dir_name),dpi=DPI)
		plt.close()		

	if DO_MISC and False:
		####### SINGLET COMP; TAU, HADR CHANNELS; WANT TO SEE IF TAU MATCHES LSD AND IF SINGLET AFFECTING HADR
		pltctr+=1
		HeatPlot(pltctr, "s1scomp", 27, "turbo",
			"s1mass", 24,.1,10, "log",
			"br_s1_hadr", 107, -0.02, 1.02, "linear", Size, DPI, "Heatmap","DC")
	
		pltctr+=1
		HeatPlot(pltctr, "s1scomp", 27, "turbo",
			"s1mass", 24,.1,10, "log",
			"br_s1_tata", 113, -0.02, 1.02, "linear", Size, DPI, "Heatmap","DC")

	if DO_MISC and False:	# EARLY RETURN TO NOT OVERFLOW
		pltctr+=1
		SinglePlot(pltctr, "neu1mass div p1mass", 0, 0, 0, "log", "p1dw", 135, 0, 0, "log",
				Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "DC","")
		pltctr+=1
		SinglePlot(pltctr, "neu1mass div p1mass", 0, 0, 0, "log", "br_p1_neu1neu1", 122, 0, 0, "log",
				Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "BR","")
		pltctr+=1
		SinglePlot(pltctr, "br_p1_neu1neu1", 122, 0, 0, "log", "p1dw", 135, 0, 0, "log",
				Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "DC","")
		pltctr+=1
		SinglePlot(pltctr, "neu1mass div s1mass", 0, 0, 0, "log", "s1dw", 133, 0, 0, "log",
				Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "DC","")
		pltctr+=1
		SinglePlot(pltctr, "neu1mass div s1mass", 0, 0, 0, "log", "br_s1_neu1neu1", 121, 0, 0, "log",
				Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "BR","")
		pltctr+=1
		SinglePlot(pltctr, "br_s1_neu1neu1", 121, 0, 0, "log", "s1dw", 133, 0, 0, "log",
				Label, Color, Alpha, Size, LOC, BBOX_TO_ANCHOR, DPI, "DC","")
		print(Time(),"p1 heats")
		pltctr+=1
		HeatPlot(pltctr, "br_p1_neu1neu1", 122, "turbo",
			"neu1mass div p1mass", 0,0,0, "log",
			"p1dw", 135, 0, 0, "log", Size, DPI, "Heatmap","")	
		pltctr+=1
		HeatPlot(pltctr, "br_p1_neu1neu1", 122, "turbo",
			"neu1mass div p1mass", 0,0,0, "log",
			"p1mass", 36, 0, 0, "log", Size, DPI, "Heatmap","")	
		pltctr+=1
		HeatPlot(pltctr, "neu1mass div p1mass", 0, "turbo",
			"p1mass", 36, 0, 0, "log",
			"br_p1_neu1neu1", 122,0,0, "log", Size, DPI, "Heatmap","")	
		print(Time(),"s1 heats")
		pltctr+=1
		HeatPlot(pltctr, "br_s1_neu1neu1", 121, "turbo",
			"neu1mass div s1mass", 0,0,0, "log",
			"s1dw", 133, 0, 0, "log", Size, DPI, "Heatmap","")	
		pltctr+=1
		HeatPlot(pltctr, "br_s1_neu1neu1", 121, "turbo",
			"neu1mass div s1mass", 0,0,0, "log",
			"s1mass", 24, 0, 0, "log", Size, DPI, "Heatmap","")	
		pltctr+=1
		HeatPlot(pltctr, "neu1mass div s1mass", 0, "turbo",
			"s1mass", 24, 0, 0, "log",
			"br_s1_neu1neu1", 121,0,0, "log", Size, DPI, "Heatmap","")		
		print(Time(),"Doing kdl histos")
		pltctr+=1
		plt.figure(pltctr)
		plt.hist([r[20]/r[19] for r in master_list[-1]], color="r", rwidth=.9)
		plt.title(save_dir_name+" : distribution of events in k div l")
		plt.xlabel("k div l")
		plt.ylabel("N")
		plt.savefig("/home/wolf/NMSSMTools_6.0.0/calculations/{}/Parameter/hist k div l.png".format(save_dir_name),dpi=DPI)
		plt.close()

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
def NearSM(mass): # fnc that checks if mass is near sm higgs mass
	buffer = 3.02083333 #buffer in GeV around central SM Higgs mass 
	#buffer = 25 # large buffer for testing 4constyle
	return (mass > 125.173611-buffer) and (mass < 125.173611+buffer)

file_matrices = [set() for file in file_tags] # for storing "out_file_matrix" of each out file	#for RAM, set of list not list of lists
print(Time(),save_dir_name)
for file_index,out_file_name in enumerate(file_names):
	with open("{}{}.dat".format(DIR, out_file_name)) as f:
		f_reader = csv.reader(f, delimiter=" ")
		ctr_lighthiggs = 0
		print(Time(),"Reading in: {}...".format(out_file_name))
		for indexrow,fullrow in enumerate(f_reader):
			if DEBUG_MODE: 
				if indexrow%10000==0: print(Time(),indexrow)
			row = [0] # trim out strange spacing ---> this used to be the event number

			reject_row = False
			if "108035020" in file_prefix: last_element=74	#trunc out after MCHA(1)
			elif (DO_DC or
				file_prefix[:6] in ["PQp1v5","PQp1v8"] or 
				"lighthiggs" in save_dir_name): last_element=182	#(don't trunc)
			elif DO_BR: last_element = 149
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
				elif (file_prefix[:6] in ["PQp1v5","PQp1v6","PQp1v8","PQp1v9"] or
					"lighthiggs" in file_prefix):
					if len(row)==21 and row[20]/row[19] > 1: #redundant
						reject_row = True		# since doing the generation
						break				#  with this in mind

					if ("s1mp1m_le_10" in save_dir_name and (
						(len(row)==25 and row[24]>10) or(len(row)==37 and row[36]>10))):
						reject_row = True
						break
					if ("s1m_le_10" in save_dir_name and (
						(len(row)==25 and row[24]>10) )):
						reject_row = True
						break
					if ("p1m_le_10" in save_dir_name and (
						(len(row)==37 and row[36]>10) )):
						reject_row = True
						break
					if ("s1m_near_p1m" in save_dir_name and (
						(len(row)==37 and (row[24]/row[36] < .5 or row[24]/row[36]>2) ) )):
						reject_row = True
						break
					if ("s1mp1m_lt_2neu1m" in save_dir_name and (
						(len(row)==45 and (row[24] >= 2*row[44] or row[36]>=2*row[44]) ) )):
						reject_row = True
						break
					if ("p1m_lt_2neu1m" in save_dir_name and (
						(len(row)==45 and row[36]>=2*row[44] ) )):
						reject_row = True
						break
					if ("s1m_lt_2neu1m" in save_dir_name and (
						(len(row)==45 and row[24]>=2*row[44] ) )):
						reject_row = True
						break
					if ("s1m_lt_2p1m" in save_dir_name and (
						(len(row)==37 and row[24]>=2*row[36] ) )):
						reject_row = True
						break
					if ("p1m_gt_p5s1m" in save_dir_name and (# same as s1m_lt_2p1m
						(len(row)==37 and row[24]>=2*row[36]))):
						reject_row = True
						break
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
					elif ("_H-s" in save_dir_name and (		#Higgsino-Singlino
					(len(row)==24 and row[23]<400        ) or	# reject low mueff
					(len(row)==46 and row[45]**2>=0.05   ) or	# reject Bino >=5%
					(len(row)==75 and row[74]-row[44]<10)  ) ):	# cha1>neu1+10 mass
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

		##	if not reject_row and last_element > 128: row[128]=0 # SCUFFED MANUAL OVERWRITE OF GGF H SIGMA	
		
			# cast row as tuple as cant have list as content of a set
			row = tuple(row)
			if not reject_row: file_matrices[file_index % len(file_tags)].add(row) # if list(list) append, if set(lists) add
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

def Set(List): # fn is List to set of tuples conversion
	return set(map(tuple, List))
def List(Set): # fn is Set to List conversion
	return list(map(list, Set))

if num_files == 1:
	master_list = List(file_matrices[0])
elif CMYK:
	print(Time(),"Splitting into sets...")
	#for i,e in enumerate(file_matrices[2]):	#BE WEARY THIS IS OVERWRITING ORIGINAL INFO
	#	e[-2]=0				# ON THE PROD SIGMA THRU GGF / con2 has nonzero
			# above arg corresponds to LHC FILE, is [3] for [ "","con1","con3","con2"]

	if False: #leaving this here but copying this architecture for the THY/LEP/LHC/BKF idea...
		bset0 = file_matrices[0]		# base unc. set
		bset1 = file_matrices[1]		# base con1 set
		bset3 = file_matrices[2]		# base con3 set
		bset2 = file_matrices[3]		# base con2 set
		
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
		master_list = [set() for i in range(15)]

		if True:	#super specific toggle since with n4e5 this gets killed handling a bunch of values
			#master_list[0].append(Set(file_matrices[0]))	# THY ==> 0/T
			master_list[0] |= (file_matrices[0])
			print(Time(),"0")
			if not MASSTRKFILE:
				file_matrices[0] = None

			master_list[1] |= (file_matrices[1])	# LEP ==> 1
			print(Time(),"1")
			if not MASSTRKFILE:
				file_matrices[1] = None

			master_list[2] |= (file_matrices[2])	# LHC ==> 2
			print(Time(),"2")
			if not MASSTRKFILE:
				file_matrices[2] = None

			master_list[3] |= (file_matrices[3])	# BKF ==> 3
			print(Time(),"3")
			if not MASSTRKFILE:
				file_matrices[3] = None
				del file_matrices

			master_list[4] = master_list[0] & master_list[1]	# T1
			master_list[0] = master_list[0].difference(master_list[4])
			master_list[1] = master_list[1].difference(master_list[4])
			if DEBUG_MODE: print(Time(),"T1")

			master_list[10] = (master_list[4] & master_list[2])	# T12
			master_list[4] = master_list[4].difference(master_list[10])
			master_list[2] = master_list[2].difference(master_list[10])
			if DEBUG_MODE: print(Time(),"T12")

			master_list[14] = (master_list[10] & master_list[3])	# T123
			master_list[10] = master_list[10].difference(master_list[14])
			master_list[3] = master_list[3].difference(master_list[14])
			if DEBUG_MODE: print(Time(),"T123")

			master_list[11] = (master_list[4] & master_list[3])	# T13
			master_list[4] = master_list[4].difference(master_list[11])
			master_list[3] = master_list[3].difference(master_list[11])
			if DEBUG_MODE: print(Time(),"T13")

			master_list[5] = (master_list[0] & master_list[2])	# T2
			master_list[0] = master_list[0].difference(master_list[5])
			master_list[2] = master_list[2].difference(master_list[5])
			if DEBUG_MODE: print(Time(),"T2")

			master_list[12] = (master_list[3] & master_list[5])	# T23
			master_list[3] = master_list[3].difference(master_list[12])
			master_list[5] = master_list[5].difference(master_list[12])
			if DEBUG_MODE: print(Time(),"T23")

			master_list[6] = (master_list[0] & master_list[3])	# T3
			master_list[0] = master_list[0].difference(master_list[6])
			master_list[3] = master_list[3].difference(master_list[6])
			if DEBUG_MODE: print(Time(),"T3")

			master_list[7] = (master_list[1] & master_list[2])	# 12
			master_list[1] = master_list[1].difference(master_list[7])
			master_list[2] = master_list[2].difference(master_list[7])
			if DEBUG_MODE: print(Time(),"12")

			master_list[13] = (master_list[3] & master_list[7])	# 123
			master_list[3] = master_list[3].difference(master_list[13])
			master_list[7] = master_list[7].difference(master_list[13])
			if DEBUG_MODE: print(Time(),"123")

			master_list[8] = (master_list[1] & master_list[3])	# 13
			master_list[1] = master_list[1].difference(master_list[8])
			master_list[3] = master_list[3].difference(master_list[8])
			if DEBUG_MODE: print(Time(),"13")

			master_list[9] = (master_list[2] & master_list[3])	# 23
			master_list[2] = master_list[2].difference(master_list[9])
			master_list[3] = master_list[3].difference(master_list[9])
			if DEBUG_MODE: print(Time(),"23")
			for i in range(len(master_list)):
				master_list[i] = List(master_list[i])
				if DEBUG_MODE: print(Time(),"{}:\t{}".format(i,len(master_list[i])))
			
			print(Time(),
				"|T123) {: >6}".format(len(master_list[14])))
			print("\t| T12) {: >6}\t| T13) {: >6}\t| T23) {: >6}\t| 123) {: >6}"
				.format(len(master_list[10]),len(master_list[11]),len(master_list[12]),len(master_list[13])))
			print("\t|  T1) {: >6}\t|  T2) {: >6}\t|  T3) {: >6}\t|  12) {: >6}\t|  13) {: >6}\t|  23) {: >6}".format(len(master_list[4]),len(master_list[5]),len(master_list[6]),len(master_list[7]),len(master_list[8]),len(master_list[9])))
			print("\t|   T) {: >6}\t|   1) {: >6}\t|   2) {: >6}\t|   3) {: >6}"
				.format(len(master_list[0]),len(master_list[1]),len(master_list[2]),len(master_list[3])))

####################### PRIOR VER OF SCUFFED VER
#			sT123 = (bT & b1)
#			#master_list[4] = List(sT123)
#			master_list[14] += List(sT123)
#			del bT, b1
#			sT123 = sT123 & b2
#			#master_list[10] = List(sT123)
#			master_list[14] += List(sT123)
#			del b2
#			sT123 = sT123 & b3
#			master_list[14] += List(sT123)
#			del b3, file_matrices
		else:	
			print(Time(),"Base sets and intersection")		

			bT = file_matrices[0]				# base sets exactly 1 con. applied
			b1 = file_matrices[1]
			b2 = file_matrices[2]
			b3 = file_matrices[3]

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
		# promote T1 and T12 to st123 for very specific case
		if "passedT" in save_dir_name:
			master_list[14] += master_list[0] + master_list[4]+ master_list[5]+ master_list[6]+ master_list[10]+ master_list[11]+ master_list[12]
			master_list[0] = []
			master_list[4] = []
			master_list[5] = []
			master_list[6] = []
			master_list[10] = []
			master_list[11] = []
			master_list[12] = []
		elif "sT1sT12" in save_dir_name:
			master_list[14] += master_list[4] + master_list[10]
			master_list[4] = []
			master_list[10] = []

# # # # ONCE SETS DETERMINED...
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
	for file_index,out_file_matrix in enumerate(List(file_matrices)):
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
	for file_index,out_file_matrix in enumerate(List(file_matrices)):
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
	for threshold_lighthiggs in [1,10]:
		ls_ctr = 0
		lp_ctr = 0
		for r in master_list[-1]:
			s1str = "                  "
			p1str = "                  "
			if r[24]<threshold_lighthiggs:
				ls_ctr+=1
				s1str = "s1mass {: >4} < {: >2} :".format(round(r[24],1),threshold_lighthiggs)
			if r[36]<threshold_lighthiggs:
				lp_ctr+=1
				p1str = "p1mass {: >4} < {: >2} :".format(round(r[36],1),threshold_lighthiggs) 
			if ":" in s1str or ":" in p1str:
				print(s1str,p1str,r[1], r[19], r[20], r[21], r[22], r[23])

		print("Light Mass Threshold:\t{}\n# Light Scalars:\t{}\n# Light Pseudoscalars:\t{}".format(threshold_lighthiggs,ls_ctr,lp_ctr))
	# # # # # WITHOUT FILTERING POINTS... # # # # #
	print(Time(),"For all valid events:")
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
	# # # # # SELECTING s1 LIGHT SINGLET EVENTS. # # # # #
	print(Time(),"After picking events with ######")
	(t_lo, t_hi) = (master_list[-1][0][1], master_list[-1][0][1]) 
	(l_lo, l_hi) = (master_list[-1][0][19],master_list[-1][0][19])
	(k_lo, k_hi) = (master_list[-1][0][20],master_list[-1][0][20])
	(Al_lo,Al_hi)= (master_list[-1][0][21],master_list[-1][0][21])
	(Ak_lo,Ak_hi)= (master_list[-1][0][22],master_list[-1][0][22])
	(mu_lo,mu_hi)= (master_list[-1][0][23],master_list[-1][0][23])
	for r in master_list[-1]:
		if r[24] > 25 and r[36] > 25: continue	# ignore heavy s1/p1
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
	row = lambda par, i : FunctionArr(par, 0, 0, master_list[-1])[i] #row(par, i)	#for FunctionArr ease
	print("BENCHMARK POINTS BY VALUES:     tanB     lambda        kappa   Alambda     Akappa    mueff")
	for i,r in enumerate(master_list[-1]):
#		if r[130]+r[131]+r[132]+r[133] < 1: #brneu2tot
#		if r[134]+r[135]+r[136]+r[137] < 1: #brneu3tot
#		if r[138]+r[139]<1:	#br_cha1_wneu1 + _hcneu1
#		if r[24]<122 and r[36]<15: #lower s1 p1 masses
#		if r[24]<10:			# small s1mass
#			if FunctionArr("neu1Hcomp",0,0,master_list[-1])[i] < 0.5:
#		if r[151]<1E-13:	#tiny p1dw
#		if r[44]/r[36] > 0.5 and r[122]>0: #neu1mass/p1mass > 1/2 but br_p1_neu1neu > 0 despite disallowed
#		if r[149] < 1E-16 or r[151] < 1E-16:	#either s1/p1dw small (ctau~>1m)
#		if r[24]>10 and r[36]>10:	# show me events where both heavy
#		if r[24]<0.1 or r[36]<0.1:	# show me events if either SUPER light
#		if r[36]<0.135 and r[140]>0:	# show me events below pion mass MA but with nonzero hadr dec
		if 10<r[36]:
			print("s1mass {: >5.1f} & p1mass {: >4}: {: >8.5f} {: <10} {: >12} {: >9} {: >10} {: >8}".format(round(r[24],1),round(r[36],1),r[1],r[19],r[20],r[21],r[22],r[23]))#"\ts1m:",r[24],"\tp1m:",r[36])	#,"\ts1dw:",r[140],"\tp1dw:",r[142])
print("{}\tFinished.\n#=#=#=#=#=#=#=#=#=#=#=#=#=#=#".format(Time()))
#sys.exit()
