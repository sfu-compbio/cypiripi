#!/usr/bin/env python
# 786

import os, sys, random
import multiprocessing

chr = ''
with open('data/chr22.fa') as f:
	for l in f: 
		if l[0] != '>': chr += l.strip()
S = 0
while chr[S] == 'N': S += 1
E = len(chr) - 1
while chr[E] == 'N': E -= 1

all = []
F = [ 'data/reference.combined.fa', 'data/sedmice.txt', 'data/osmice.txt' ]
for fn in F:
	a = {}
	with open(fn) as f:
		for l in f:
			n = l.split()[0][1:]
			l = f.next()
			a[n] = l.strip()
	all.append(a)

cyp6 = (42122966 - 1, 42132410)
cyp7 = (42140205 - 1, 42144549)
cyp8 = (42149834 - 1, 42155088)


# 31612965 31622410
# 31630204 31634549
# 31639833 31645088

# size of new genome 40298467

# 44.89

A = chr[S:cyp6[0]]
B = chr[cyp6[1]:cyp7[0]]
C = chr[cyp7[1]:cyp8[0]]
D = chr[cyp8[1]:E]

derivs = [	"CYP2D6*2A",	"CYP2D6*2M",	"CYP2D6*4N",	"CYP2D6*10D",	"CYP2D6*13A",
			"CYP2D6*13B",	"CYP2D6*13Ca",	"CYP2D6*13Cb",	"CYP2D6*13D",	"CYP2D6*13E",	
			"CYP2D6*13Fa",	"CYP2D6*13Fb",	"CYP2D6*13G",	"CYP2D6*13H",	"CYP2D6*31",	
			"CYP2D6*35B",	"CYP2D6*36T",	"CYP2D6*36S",	"CYP2D6*41",	"CYP2D6*51",	
			"CYP2D6*56A",	"CYP2D6*57",	"CYP2D6*58",	"CYP2D6*61", 	"CYP2D6*63",	
			"CYP2D6*68A",	"CYP2D6*82",	"CYP2D6*83", 	"CYP2D7" ]


PAT = 'fasta/F'
def gen_pair_fuse(a, ID, has7, typ):
	rep = all[0]["CYP2D6*1"][:3523]
	y = random.choice(all[1].keys())
	z = random.choice(all[2].keys())
	
	Bx = B
	if not typ:
		Bx = chr[cyp6[1]:(cyp7[0] - 5630)] # remove REP7

	if not has7: 
		y = ''

	for i in xrange(0,2):
		fn = PAT + '{0}_{1}'.format(ID, 'MF'[i])
		with open(fn + '.fa', 'w') as fo:
			x = a[i]
				
			print '*' * 80

			print i, ' '.join(x), y, z
			print >>fo, '>chr22' + 'MF'[i] + ' {0} {1} {2}'.format(' '.join(x), y, z)

			p = ''
			if len(x) == 1:
				p = all[0][x[0]]
			else:
				p = ''.join([A[-100:] + all[0][c] + Bx[:100] for c in x])

			print >>fo, A + p + Bx + ('' if y == '' else all[1][y]) + C + all[2][z] + D
		#gen_reads(fn + '.fa', fn + '.fq', 10)

choices = [ 
	# no 7
	[["CYP2D6*13A"],["CYP2D6*13A"]], 
	[["CYP2D6*1","CYP2D6*13A"],["CYP2D6*1","CYP2D6*13A"]], 
	[["CYP2D6*13Ca"],["CYP2D6*13Ca"]], 
	[["CYP2D6*2A","CYP2D6*13D"],["CYP2D6*2A","CYP2D6*13D"]], 
	[["CYP2D6*13Fb"],["CYP2D6*13Fb"]], 
	[["CYP2D6*1", "CYP2D6*1","CYP2D6*13H"],["CYP2D6*1", "CYP2D6*1","CYP2D6*13H"]],
]

cc = 65
for c in choices:
	gen_pair_fuse(c, cc, False, False)
	cc += 1


for k in [ "68A" ]:
	all[0]["CYP2D6*" + k] = all[0]["CYP2D6*1"][:3523] + all[0]["CYP2D6*" + k] 

choices = [[["CYP2D6*36S"], ["CYP2D6*36S"]]]
for c in choices:
	gen_pair_fuse(c, cc, False, True)
	cc += 1

# has 7
choices = [	[["CYP2D6*4N","CYP2D6*4N","CYP2D6*36S"], ["CYP2D6*4N", "CYP2D6*4N", "CYP2D6*36S"]],
			[["CYP2D6*4A","CYP2D6*68A"],["CYP2D6*4A","CYP2D6*68A"]],
			[["CYP2D6*10A","CYP2D6*57"],["CYP2D6*10A","CYP2D6*57"]],
			[["CYP2D6*82"],["CYP2D6*82"]],
			[["CYP2D6*2A"],["CYP2D6*2A"]]
]

for c in choices:
	gen_pair_fuse(c, cc, True, True)
	cc += 1

exit(1)

for d in derivs:
	del(all[0][d])

def gen_reads(fasta, fastq, cov):
	cmd = 'time simNGS/bin/simLibrary -x {2} {0} | simNGS/bin/simNGS -o fastq -p paired simNGS/data/s_3_4x.runfile > {1}'.format(fasta, fastq, cov)
	os.system(cmd)

PAT = 'fasta/F'
	
def gen_pair(a, ID):
	y = random.choice(all[1].keys())
	z = random.choice(all[2].keys())
	for i in xrange(0,2):
		fn = PAT + '{0}_{1}'.format(ID, 'MF'[i])
		with open(fn + '.fa', 'w') as fo:
			x = a[i]
				
			print '*' * 80
			print i, ' '.join(x), y, z
			print >>fo, '>chr22' + 'MF'[i] + ' {0} {1} {2}'.format(' '.join(x), y, z)

			p = ''
			if len(x) == 1:
				p = all[0][x[0]]
			else:
				p = ''.join([A[-100:] + all[0][c] + B[:100] for c in x])

			print >>fo, A + p + B + all[1][y] + C + all[2][z] + D
		#gen_reads(fn + '.fa', fn + '.fq', 10)

	#os.system('cat ' + PAT + '{0}_[MF].fq > ' + PAT + '{0}.fq'.format(ID))
	#os.system('rm -rm ' + PAT + '{0}_[MF].fq'.format(ID))

keys = all[0].keys()
choices = []

# 1. simple 
p = [ 'CYP2D6*1', 'CYP2D6*15', 'CYP2D6*4M', 'CYP2D6*6A', 'CYP2D6*27', 'CYP2D6*40' ]
random.shuffle(keys)
for k in keys:
	if not k in p:
		p.append(k)
	if len(p) == 20:
		break
for pp in p:
	choices.append([[pp], [pp]])
	
# 2. M-F different
random.shuffle(keys)
for i in xrange(0, len(keys), 2):
	#p = [keys[i], keys[i + 1]]
	choices.append([[keys[i]], [keys[i + 1]]])
	if len(choices) == 40:
		break

# 3. tandems
p = [ 	[ 'CYP2D6*35X', 'CYP2D6*35X' ], 
		[ 'CYP2D6*4A', 'CYP2D6*4A' ], 
		[ 'CYP2D6*9X', 'CYP2D6*9X' ], 
		[ 'CYP2D6*10A', 'CYP2D6*10A' ] ]
for j in [2, 3, 4, 5, 8, 13]:
	p.append(['CYP2D6*2X'] * j)
	p.append(['CYP2D6*1'] * j)
p.append(['CYP2D6*17'] * 8)
p.append([]) # CYP2D6*5

while len(p) < 25:
	r = random.randint(2, 8)
	x = []
	while r > 0:
		x.append(random.choice(keys))
		r -= 1
	p.append(x)

for pp in p:
	choices.append([pp, pp])

# def fun(i):
#  	gen_pair(choices[i], '%02d'%(i+80))

# print choices

#pool = multiprocessing.Pool(3)
#pool.map(fun, range(len(choices)))
i = 0
for f in choices:
	gen_pair(f, '%02d'%i)
	i += 1



