import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import time

base_dat_file_name = input("What file would you like to analyze (spectr)? ")
pref_dat_file_name = input("What file PREFIX would you like (TEST)? ")

program_start_time = time.time()


list_MH_v_Lambda = list()
list_MH_v_Lambda.append(("MH","MA","Lambda"))

for ITER in range(1,11):
	dat_file_name = "TEST"+base_dat_file_name+"{}".format(ITER)+".dat"

	# open specified file: dat_file_name as f
	f = open('/home/wolf/NMSSMTools_6.0.0/calculations/{}'.format(dat_file_name),"r")

	rec_M = False
	rec_H = False    # should I be recording NMHMIX right now?
	rec_A = False    #   .    .  .    .      NMAMIX .      . ?
	rec_coup = False #   .    .  .    . EFFECTIVE  COUPLINGS . .?
	rec_extpar = False

	# want to pull from: BLOCK MASS, BLOCK NMHMIX, BLOCK NMAMIX, BLOCK  EFFECTIVE_COUPLINGS
	BLOCK_M = list()
	BLOCK_H = list()
	BLOCK_A = list()
	BLOCK_COUP = list()
	BLOCK_EXTPAR = list()
	for row in f:
		if not rec_H and not rec_A and not rec_coup and not rec_M and not rec_extpar:
			if "BLOCK NMHMIX" in row: 
				#print(row)
				BLOCK_H.append(row)
				rec_H = True
			elif "BLOCK NMAMIX" in row:
				#print(row)
				BLOCK_A.append(row)
				rec_A = True
			elif "BLOCK  EFFECTIVE_COUPLINGS" in row:
				#print(row)
				BLOCK_COUP.append(row)
				rec_coup = True
			elif "BLOCK MASS" in row:
				BLOCK_M.append(row)
				rec_M = True
			elif "BLOCK EXTPAR" in row:
				BLOCK_EXTPAR.append(row)
				rec_extpar = True
		else:
			if rec_H:
				if row == "# \n": rec_H = False
				else: BLOCK_H.append(row)
			elif rec_A:
				if row == "# \n": rec_A = False
				else: BLOCK_A.append(row)
			elif rec_coup:
				if row == "# \n": rec_coup = False
				else: BLOCK_COUP.append(row)
			elif rec_M:
				if row == "# \n": rec_M = False
				else: BLOCK_M.append(row)
			elif rec_extpar:
				if row == "# \n": rec_extpar = False
				else: BLOCK_EXTPAR.append(row)

	print("\n{}".format(dat_file_name))
	for i in BLOCK_H: 
		#print(i,end="")
		pass
	for j in BLOCK_A: 
		#print(j,end="")
		pass
	for k in BLOCK_COUP: 
		#print(k,end="")
		pass
	for l in BLOCK_M: 
		print(l,end="")
		if "# lightest neutral scalar" in l:
			MH = l[15:-29]
		elif "# lightest pseudoscalar" in l:
			MA = l[15:-27]
	for r in BLOCK_EXTPAR:
		if "# LAMBDA" in r:
			Lambda = float(r[11:-12])

	list_MH_v_Lambda.append((MH,MA, Lambda))
	f.close()
for i in list_MH_v_Lambda:
	print(i)
	




#df = pd.read_csv('/home/wolf/NMSSMTools_6.0.0/calculations/{}'.format(dat_file_name), delimiter=' ')
#
#bigHead = "# NMSSMTools OUTPUT IN SLHA FORMAT"
#
#print(df.head(), "\n")
#print(df.iat[1200,0], "\n")
#print(type(df.iat[1200,0]), "\n")
#block_num = -1
#for i in range(df.size):
#	if "BLOCK" not in df.iat[i,0][:5] and block_num == -1: #only start caring when you get to the first BLOCK
#		continue
#	if "BLOCK" in df.iat[i,0][:5]: #if starting a new block
#		block_num += 1
#		block_header = df.iat[i,0]
#		block_name = "BLOCK_{}".format(block_num)
#		print(block_header, "\n")
#		exec(block_name + " = pd.DataFrame(columns = ['{}'])".format(block_header))
#	else:
#		exec(block_name + ".append('{}')".format(df.iat[i,0]))
#		print(eval(block_name))
#
#
#		print(eval(block_name))


print("\n-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
program_fin_time = time.time()
print("Runtime (s):\t", program_fin_time - program_start_time)
