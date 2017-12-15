# 786
from itertools import izip
import os

seqs = []
with open(os.path.dirname(os.path.realpath(__file__)) + '/meho.stats') as f:
	for l in f:
		if l[0] == '>':
			seqs.append(l.strip().split()[1])
s = ''
for j in xrange(len(seqs[0])):
	c = sorted(list(set([t[j] for t in seqs[1:] if t[j] != 'N'])))
	if len(c) == 1: 
		a = c[0]		
	elif len(c) == 2:
		l = ['AC', 'AG', 'AT', 'CG', 'CT', 'GT', '-A', '-C', '-G', '-T'] 
		a = 'BDEFHIJKLMOPQRSUVWXYZ'[l.index(c[0] + c[1])]
	else:    
		print c, j
		raise 'ooops'

	if seqs[0][j] == 'x': #snip
		s += a
	elif seqs[0][j] == '-':
		if seqs[1][j] != '-': #ignore this ins
			s += a.lower()
	elif seqs[0][j] != a:
		if a in 'JKLM' and seqs[1][j] != '-': # hack for LJMLMKJ
			s += '='
		else:
			s += a
	else:
		s += '='
print s
