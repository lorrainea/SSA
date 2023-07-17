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
	
	run_greg = r'../greg/./ssa text suffixes ../greg/out'
	os.system(run_greg)
	
	run_acc = r'./acc_ssa text suffixes out'
	os.system(run_acc)
	
	run_diff_ssa = r'diff out_accurate.ssa ../greg/out.ssa'	
	os.system(run_diff_ssa)
	
	run_diff_lcp = r'diff out_accurate.lcp ../greg/out.lcp'	
	os.system(run_diff_lcp)
	
	
	#run_lorraine = r'.././ssa text suffixes ../out'
	#os.system(run_lorraine)
	
	#run_diff_ssa = r'diff out_accurate.ssa ../out.ssa'	
	#os.system(run_diff_ssa)
	
	#run_diff_lcp = r'diff out_accurate.lcp ../out.lcp'	
	#os.system(run_diff_lcp)
