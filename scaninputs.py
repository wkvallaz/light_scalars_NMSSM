import numpy as np
import matplotlib.pyplot as plt
import time
import shutil # utilities for copying files

start_time = time.time()

### THIS FILE TO BE WHAT CREATES CUSTOM PREFIXinpSUFFIX.dat files to be ran with ./run PREFIXinpSUFFIX.dat

PREFIX = "TEST"
SUFFIX = "ITER"

# first objective, to copy inp.dat and be able to rewrite inside of the copy
# possibly not necessary, as can just pull info and write new file
# util.copyfile("/home/wolf/NMSSMTools_6.0.0/calculations/inp.dat", "/home/wolf/NMSSMTools_6.0.0/calculations/{}inp{}.dat".format(PREFIX, SUFFIX))

for ITER in range(1,11):

	f = open("/home/wolf/NMSSMTools_6.0.0/calculations/inp.dat", "r")

	lines = list()
	for index,row in enumerate(f): # MODIFICATIONS TO BLOCK EXTPAR
	#	if "# MSUSY" in row:
	#		lines.append("#\t0\t{}\t# MSUSY (If =/= SQRT(2*MQ1+MU1+MD1)/2) modified\n".format("1000d0"))
	#	elif "# M1" in row:
	#		lines.append("#	1	{}		# M1 (If =/= M2/2) modified\n".format("0d0"))
	#	elif "# M2" in row:
	#		lines.append("	2	{}	# M2 modified\n".format("0.248401E+03"))
	#	elif "# M3" in row:
	#		lines.append("#	3	{}		# M3 (If =/= 3*M2) modified\n".format("0d0"))
	#	elif "# AU3" in row:
	#               lines.append("	11	{}		# AU3 modified\n".format("1500d0"))
	#	elif "# AD3" in row:
	#               lines.append("	12	{}		# AD3 modified\n".format("1500d0"))		
	#	elif "# AE3" in row:
	#               lines.append("	13	{}		# AE3 modified\n".format("1500d0"))
	#	elif "# AE2" in row:
	#		lines.append("#	16	{}		# AE2 = AE1 (If =/= AE3) modified\n".format("0d0"))
	#	elif "# ML3" in row:
	#		lines.append("	33	{}		# ML3 modified\n".format("200d0"))
	#	elif "# ML2" in row:
	#		lines.append("#	32	{}		# ML2 = ML1 (If =/= ML3) modified\n".format("0d0"))
	#	elif "# ME3" in row: 
	#		lines.append("	36	{}		# ME3 modified\n".format("200d0"))
	#	elif "# ME2" in row:
	#		lines.append("#	35	{}		# ME2 = ME1 (If =/= ME3) modified\n".format("0d0")) 
	#	elif "# MQ3" in row: 
	#		lines.append("	43	{}		# MQ3 modified\n".format("1000d0"))
	#	elif "# MQ2" in row: 
	#		lines.append("#	42	{}		# MQ2 = MQ1 (If =/= MQ3) modified\n".format("0d0"))
	#	elif "# MU3" in row: 
	#		lines.append("	46	{}		# MU3 modified\n".format("1000d0"))
	#	elif "# MU2" in row: 
	#		lines.append("#	45	{}		# MU2 = MU1 (If =/= MU3) modified\n".format("0d0"))
	#	elif "# MD3" in row: 
	#		lines.append("	49	{}		# MD3 modified\n".format("1000d0"))
	#	elif "# MD2" in row: 
	#		lines.append("#	48	{}		# MD2 = MD1 (If =/= MD3) modified\n".format("0d0"))
		if "# LAMBDA" in row: 
	#		lines.append("	61	{}	# LAMBDA modified\n".format("0.499793E+00"))
			lines.append("\t61\t{}\t#LAMBDA modified\n".format(.3 + .04*ITER))
	#	elif "# KAPPA" in row: 
	#		lines.append("	62	{}		# KAPPA (If =/= 0) modified\n".format("0d0"))
	#	elif "# ALAMBDA" in row: 
	#		lines.append("	63	{}	# ALAMBDA (If XIF+MA are not inputs) modified\n".format("0.259611E+04"))
	#	elif "# AKAPPA" in row: 
	#		lines.append("#	64	{}		# AKAPPA (If KAPPA =/=0 and XIS+MP are not inputs) modified\n".format("0d0"))
	#	elif "# MUEFF" in row:
	#		lines.append("\t65\t{}\t# MUEFF modified\n".format("0.373538E+03"))
	#	elif "# XIF" in row:
	#		lines.append("#\t66\t{}\t\t# XIF in GeV^2 (If ALAMBDA+MA are not inputs) modified\n".format("0d0"))
	#	elif "# XIS" in row: 
	#		lines.append("#	67	{}	 	# XIS in GeV^3 (If AKAPPA+MP are not inputs) modified\n".format("0d0"))
	#	elif "# MUP" in row: 
	#		lines.append("#	68	{}		# MUP (If =/= 0) modified\n".format("0d0"))
	#	elif "# MSP" in row: 
	#		lines.append("#	69	{}		# MSP in GeV^2 (If =/= 0) modified\n".format("0d0"))
	#	elif "# M3H" in row: 
	#		lines.append("#	72	{}		# M3H in GeV^2 (If =/= 0) modified\n".format("0d0"))
	#	elif "# MA" in row: 
	#		lines.append("	124	{}	# MA (If ALAMBDA+XIF are not inputs) modified\n".format("0.259434E+04"))
	#	elif "# MP" in row: 
	#		lines.append("	125	{}	# MP (If AKAPPA+XIS are not inputs) modified\n".format("0.863261E+02"))
		else:
			lines.append(row)

	#for i in lines:
	#	print(i)
	f.close()

	SUFFIX = str(ITER)
	#WRITES to a new file specified by PREFIX, SUFFIX
	f = open("/home/wolf/NMSSMTools_6.0.0/calculations/{}inp{}.dat".format(PREFIX, SUFFIX), "w")
	for i in lines: f.write(i)
	f.close()








print(time.time()-start_time)
