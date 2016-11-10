'''
Motif PFMs comparison tools : Compares two PFM files and returns similar kmers
Input: List of PFM files 1)known kmers(motifs) 2)new kmers
Output: Motifs in 2) that matched with 1)
'''

import os,sys
import re
import regex

def main(argv):
	in_dir= "/projects_rg/babita/TCGA/pipeline/run12/kmers_analysis/test"
	known_pfms_dir = os.path.join(in_dir,"known_pfms")
	new_pfms_dir = os.path.join(in_dir,"new_pfms")

	'''
	for each known pfms compare it with new ones and create table
	pfm_new	matched_pfms	times_matched	 kmer_length
	'''

	known_pfmfiles = [f for f in os.listdir(known_pfms_dir) if os.path.isfile(os.path.join(known_pfms_dir,f)) ]
	new_pfmfiles =  [f for f in os.listdir(new_pfms_dir) if os.path.isfile(os.path.join(new_pfms_dir,f)) ]

	for fname in known_pfmfiles:
		print fname
		with open(fname) as f:
			content = f.readlines()


'''
class Pfm(object):

	#the class initializer
	
	def __init__(self, len_kmer, motif_id, A_freq, C_freq, G_freq, T_freq):
		self.len_kmer = len_kmer
		self.motif_id = motif_id
		self.A_freq = A_freq
		self.C_freq = C_freq
		self.G_freq = G_freq
		self.T_freq = T_freq

	def A_freq(self):

		cnt = 0
		a = []
		
		for line in content:
				cnt += 1
				if line.startswith('letter'):
					kmer_size = re.search('w=\s[0-9]\s+', line)
					kmer_size = kmer_size.group(0)
					kmer_size = int(kmer_size.split()[1])
					
					for i in range(kmer_size):
						a_line = content[cnt + i]
						a_line.rstrip('\n')
						a_line = float(a_line.split()[0])
						a.append(a_line)
		return a


def make_pfms(len_kmer, motif_id, A_freq, C_freq, T_freq):
	pfm_known = Pfm(len_kmer, motif_id, A_freq, C_freq, T_freq)
	pfm_new = Pfm(len_kmer, motif_id, A_freq, C_freq, T_freq)
	return pfm

'''

if __name__ == "__main__":
	main(sys.argv)