import random
import string
import sys
import os

seq_len = 20 #can be changed
num_suffixes = 10 #can be changed
runs = 1 #can be changed

for i in range(0, runs):
	letters = string.ascii_lowercase
	sequence = ''.join(random.choice(letters) for i in range(seq_len))
	    
	seq_file = open("text", "w")
	seq_file.write(sequence)
	seq_file.close()

	suffixes = random.sample(range(0, seq_len), num_suffixes)
	
	suff_file = open("suffixes", "w")
	for i in range(0, len(suffixes)):
		suff_file.write(str(suffixes[i]))
		suff_file.write("\n")
		
	suff_file.close()
	
	print("RUNNING GREG")
	run_greg = r'../greg/./ssa text suffixes ../greg/out'
	os.system(run_greg)
	
	print("RUNNING ALGORITHM3")
	run_alg3 = r'../algorithm3/./ssa text suffixes ../algorithm3/out'
	os.system(run_alg3)
	
	print("RUNNING ACCURATE")
	run_acc = r'./acc_ssa text suffixes out'
	os.system(run_acc)
	
	print("RUNNING RK-LCE")
	run_rk_lce = r'../rk-lce/./sa-rk text ../rk-lce/out suffixes'
	os.system(run_rk_lce)
	
	
	print("COMPARING GREG AND SDSL")
	run_diff_ssa = r'diff out_accurate.ssa ../greg/out.ssa'	
	os.system(run_diff_ssa)
	
	run_diff_lcp = r'diff out_accurate.lcp ../greg/out.lcp'	
	os.system(run_diff_lcp)
	
	
	print("COMPARING ALGORITHM3 AND SDSL")
	run_diff_ssa_alg3 = r'diff out_accurate.ssa ../algorithm3/out.ssa'	
	os.system(run_diff_ssa_alg3)
	
	run_diff_lcp_alg3 = r'diff out_accurate.lcp ../algorithm3/out.lcp'	
	os.system(run_diff_lcp_alg3)
	
	
	print("COMPARING GREG AND RK-LCE")
	run_diff_ssa_rk_lce = r'diff ../greg/out.ssa ../rk-lce/out_rk_lce.ssa'	
	os.system(run_diff_ssa_rk_lce)
	
	run_diff_lcp_rk_lce = r'diff ../greg/out.lcp ../rk-lce/out_rk_lce.lcp'	
	os.system(run_diff_lcp_rk_lce)
	
	
	#run_lorraine = r'.././ssa text suffixes ../out'
	#os.system(run_lorraine)
	
	#run_diff_ssa = r'diff out_accurate.ssa ../out.ssa'	
	#os.system(run_diff_ssa)
	
	#run_diff_lcp = r'diff out_accurate.lcp ../out.lcp'	
	#os.system(run_diff_lcp)
