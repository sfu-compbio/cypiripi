#!/usr/bin/env python
#786 

import os, sys, errno, argparse, subprocess, fnmatch

def command_line_process():
	parser = argparse.ArgumentParser(description='CYP2D6 genotyping pipeline!')
	parser.add_argument('--out',   '-o', required=True)
	parser.add_argument('--snpdb', '-d', default=os.path.dirname(os.path.realpath(__file__)) + '/M88833.csv')
	parser.add_argument('--cyp6',  '-6', default=os.path.dirname(os.path.realpath(__file__)) + '/CYP2D6.M33388.fa')
	parser.add_argument('--cyp7',  '-7', default=os.path.dirname(os.path.realpath(__file__)) + '/CYP2D7.M33387.fa')
	return parser

def printerr(*objs):
    print >>sys.stderr, objs

class gene:
	def __init__(self, n, s, e, i):
		self.name = n
		self.seq = s
		self.exons = e
		self.introns = i

class allele:
	def __init__(self):
		self.name = ''
		self.seq = ''
		self.c7 = ''
		self.snps = []
		self.inss = []
		self.dels = []
		self.fuss = []

def parseGene(file):
	name = ''
	seq = ''
	exons = []
	introns = []
	offset = (0, 999999)
	with open(file) as fasta:
		for l in fasta:
			if l[0] == '>':
				items = l.strip().split()
				name = items[0][1:]
				exx = items[2].split(',') if len(items) > 2 else []
				offset = items[3].split('..') if len(items) > 3 else (0, 999999)
				offset = map(int, offset)
			else:
				seq += l.strip()

	si, ei = 0, 0
	if exx == []:
		for i in xrange(len(seq)):
			if seq[i].isupper() and (i == 0 or seq[i - 1].islower()):
				ei = i
				introns.append((si, ei))
				si = i
			if seq[i].islower() and (i == 0 or seq[i - 1].isupper()):
				ei = i
				if ei > 0: exons.append((si, ei))
				si = i
		introns.append((si, len(seq)))
		seq = seq.upper()
	else:
		for e in exx:
			s, e = e.split('..')
			s = int(s) - 1
			e = int(e)
			exons.append((s,e))

			ei = s
			if si > 0: introns.append((si,ei))
			si = e	
		introns.insert(0, (offset[0], exons[0][0]))
		introns.append((exons[-1][1], min(len(seq), offset[1])))

	return gene(name, seq, exons, introns)

def shell(command):
	try:
		output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
	except Exception, e:
		output = str(e.output)
	finished = output.split('\n')
	return finished

def align(a, b):
	if (a, b) in align.mem:
		return align.mem[(a, b)]

	with open("aln1.txt", "w") as al:
		al.write(">SEQ1\n" + a)
	with open("aln2.txt", "w") as al:
		al.write(">SEQ2\n" + b)
	out = shell("./align -w 200 aln1.txt aln2.txt")[6:]
	op = ['', '', '']
	for i in xrange(len(out)):
		if out[i][:4] == 'SEQ1':
			op[0] += out[i][7:]
			op[1] += out[i + 1][7:]
			op[2] += out[i + 2][7:]
	out = op
	assert(len(out[0]) == len(out[-1]))
	assert(len(out[1]) == len(out[0]))

	aln = ''
	for i in xrange(len(out[0])):
		if out[1][i] == ':':
			aln += '='
		elif out[0][i] == '-':
			aln += out[2][i].lower()
		elif out[2][i] == '-':
			aln += '-'
		else:
			aln += out[2][i].upper()

	align.mem[(a, b)] = aln
	os.unlink('aln1.txt')
	os.unlink('aln2.txt')
	return aln
align.mem = {}

def revcomp(seq):
	rc = 'ACGTBDEFHIJKLM'
	rv = 'TGCAIHEFDBMLKJ'
	rc += rc.lower()
	rv += rv.lower()
	return ''.join([rv[rc.index(c)] if c in rc else c for c in seq[::-1]])

args = command_line_process().parse_args()
print args

cyp6 = parseGene(args.cyp6)
cyp7 = parseGene(args.cyp7)

# hic sunt dracones!
cyp2d7_fix = shell("python {0}/meho.py".format(os.path.dirname(os.path.realpath(__file__))))
cyp2d7_fix = ' ' * 11236 + revcomp(cyp2d7_fix[0][2:-26])
cyp67diff = []
cleaned_up = 0
jf = 0
for i in xrange(len(cyp6.exons) + len(cyp6.introns)):
	a, b = (cyp6.introns[i / 2] if i % 2 == 0 else cyp6.exons[i / 2], 
			cyp7.introns[i / 2] if i % 2 == 0 else cyp7.exons[i / 2])
	aln = '' #cyp2d7_fix[slice(*b)] #align(a, b)	
	if jf == 0: jf += b[0]
	l = b[0]
	while l < b[1]:
		c = cyp2d7_fix[jf]
		assert(c != ' ')
		aln += c
		jf += 1
		if c not in '-JKLM': l += 1

	a, b = cyp6.seq[slice(*a)], cyp7.seq[slice(*b)]
	
	# print a
	# print b
	# print 

	#aln = list(aln)
	#for j in xrange(len(aln)):
	#	#if aln[j] == 'X' or aln[j] == 'x':
	#	#	aln[j] = 'x' 
	#aln = ''.join(aln)

	cyp67diff.append(aln)

	# print >>sys.stderr, b
	# print >>sys.stderr, aln
	# print >>sys.stderr, a
	# print >>sys.stderr
cyp67diff[0] = '-' * 1532 + cyp67diff[0]
cyp67diff[-1] += '-' * 3523

# df = ''.join(cyp67diff)
# c6, c6x = '', cyp6.seq
# c7, c7x = '', cyp7.seq[11236:]
# for c in df:
# 	if c in '-JKLM':
# 		c6 += c6x[0]
# 		c6x = c6x[1:]
# 		c7 += '-'
# 	elif c.islower():
# 		c7 += c7x[0]
# 		c7x = c7x[1:]
# 		c6 += '-'
# 	else:
# 		c6 += c6x[0]
# 		c6x = c6x[1:]
# 		c7 += c7x[0]
# 		c7x = c7x[1:]
# print len(c6x), len(c7x)
# print c6
# print df
# print c7

cyp2d6_snps_marked = list(cyp6.seq)

cyp6.alleles = {}
with open(args.snpdb) as snps:
	for l in snps:
		l = l.strip().split()
		al = allele()
		al.name = l[0]
		i = 1
		while i < len(l):
			if l[i] == 'SNP':
				pos = 1618 + int(l[i + 1]) 
				fr = l[i + 2]
				to = l[i + 3]
				if pos < 1618: pos += 1
				
				if pos < 0:
					printerr('ERROR ', al.name, "reports invalid SNP %s:%c>%c" % (l[i+1], fr, to))
					i += 4
					continue
				if cyp6.seq[pos] != fr:
					printerr('WARN  ', al.name, "reports SNP %s:%c>%c, while we have %c" % (l[i+1], fr, to, cyp6.seq[pos]))

				cyp2d6_snps_marked[pos] = 'x'

				al.snps.append((pos, to, fr))
				i += 4
			elif l[i] == 'INS':
				pos1 = 1618 + int(l[i + 1])
				pos2 = 1618 + int(l[i + 2])
				if pos1 < 1618: pos1 += 1
				if pos2 < 1618: pos2 += 1
				pos1 = pos2 
				# fix
				inseq = l[i + 3]
				pos2 = pos1 + len(inseq)
				#if pos2 - pos1 != len(inseq):
				#	print 'WARN  ', al.name, "reports INS %s:%d-%d, while we have length %d" % (inseq, pos1, pos2, len(inseq))
				al.inss.append((pos1, pos2, inseq))
				i += 4
			elif l[i] == 'DEL':
				pos1 = 1618 + int(l[i + 1])
				pos2 = 1618 + int(l[i + 2])
				if pos1 < 1618: pos1 += 1
				if pos2 < 1618: pos2 += 1
				delseq = l[i + 3]
				pos2 += 1
				assert(cyp6.seq[pos1:pos2] == delseq)
				cyp2d6_snps_marked[pos1:pos2] = list('x' * len(delseq))
				al.dels.append((pos1, pos2, delseq))
				i += 4
			elif l[i] == 'CYP2D7':
				if l[i + 1][0] == '?':
					printerr('ERROR ', al.name, "reports unsupported CYP2D7 configuration", l[i + 1][1:])
				else:
					al.fuss.append(l[i + 1])
				i += 2
			elif l[i] == 'END':
				break
			else:
				printerr('ERROR ', al.name, "reports unsupported abberation", l[i])
				i += 1

		parts = [e for e in cyp6.exons] + [i for i in cyp6.introns] 
		parts = sorted(parts, key=lambda x: x[0])
		seqparts = [list(cyp6.seq[i[0]:i[1]]) for i in parts]
		edtparts = [list('=' * (i[1] - i[0])) for i in parts]
		c7parts = [list('6' * (i[1] - i[0])) for i in parts]
		assert(cyp6.seq == ''.join([''.join(s) for s in seqparts]))

		comment = ''
		for s in al.snps:
			found = 0
			for pi in xrange(len(parts)):
				if s[0] >= parts[pi][0] and s[0] < parts[pi][1]:
					seqparts[pi][ s[0] - parts[pi][0] ] = s[1]
					edtparts[pi][ s[0] - parts[pi][0] ] = s[1]
					#comment += 'SNP_%d ' % s[0]
					found = 1
					break
			assert(found == 1)

		for d in al.dels:
			found = 0
			for pi in xrange(len(parts)):
				if d[0] >= parts[pi][0] and d[0] < parts[pi][1] and d[1] <= parts[pi][1]:
		 			for i in xrange(d[0], d[1]):
		 				seqparts[pi][ i - parts[pi][0] ] = '-'
		 				edtparts[pi][ i - parts[pi][0] ] = '-'
		 			#	comment += 'DEL_%d ' % i
					found = 1
					break
			assert(found == 1)

		al.inss = sorted(al.inss, key=lambda x: x[0])
		off = 0
		ppi = -1
		for i in al.inss:
			found = 0
			for pi in xrange(len(parts)):
				if i[0] >= parts[pi][0] and i[0] < parts[pi][1] and i[1] < parts[pi][1]:
					if ppi != pi:
						off = 0
						ppi = pi
					seqparts[pi] = seqparts[pi][:off + i[0] - parts[pi][0]] + list(i[2]) + seqparts[pi][off + i[0] - parts[pi][0]:]
					edtparts[pi] = edtparts[pi][:off + i[0] - parts[pi][0]] + list(i[2].lower()) + edtparts[pi][off + i[0] - parts[pi][0]:]

					#for k in xrange(i[0], i[0] + len(i[2])):
					#	comment += 'INS_%d ' % k
					off += len(i[2])
					found = 1
					break
			assert(found == 1)

		for f in al.fuss:
		 	if f[0] == 'i':
		 		i = ord(f[1]) - ord('0')
		 		edtparts[i * 2] = []
		 		seqparts[i * 2] = list(cyp7.seq[cyp7.introns[i][0]:cyp7.introns[i][1]])
		 	elif f[0] == 'I': # specific region
		 		i = ord(f[1]) - ord('0')
		 		s, e = f[3:].split('..')
		 		s, e = int(s) - 1, int(e) - 1

		 		assert(s==213 and e==244)
		 		#assert(e <= cyp7.introns[i][1])
		 		#cc = (11236-1536+2)*' '+''.join(cyp67diff) #+2 here is hack
		 		#print ''.join(seqparts[i * 2][s:e])
		 		#print ''.join(list(cyp7.seq[cyp7.introns[i][0] + s:cyp7.introns[i][0] + e]))
		 		#print ''.join(list(cc[cyp7.introns[i][0] + s:cyp7.introns[i][0] + e]))
		 		edtparts[i * 2] = edtparts[i * 2][:s] + list('-==============================a') + edtparts[i * 2][e:] #HACK!
		 		seqparts[i * 2] = seqparts[i * 2][:s] + list(cyp7.seq[cyp7.introns[i][0] + s:cyp7.introns[i][0] + e]) + seqparts[i * 2][e:]
		 	elif f[0] == 'e':
		 		i = ord(f[1]) - ord('0')
		 		edtparts[i * 2 - 1] = []
		 		seqparts[i * 2 - 1] = list(cyp7.seq[cyp7.exons[i - 1][0]:cyp7.exons[i - 1][1]]) 
		 	elif f[0] == '>':
		 		l = 2 * (ord(f[2]) - ord('0')) - (1 if f[1] == 'e' else 0)
		 		for i in xrange(l, len(parts)):
		 			interval = cyp7.exons[i / 2] if i % 2 == 1 else cyp7.introns[i / 2]
			 		edtparts[i] = []
		 			seqparts[i] = cyp7.seq[interval[0]:interval[1]]
		 	elif f[0] == '<':
		 		l = 2 * (ord(f[2]) - ord('0')) - (1 if f[1] == 'e' else 0)
		 		for i in xrange(0, l + 1):
		 			interval = cyp7.exons[i / 2] if i % 2 == 1 else cyp7.introns[i / 2]
		 			edtparts[i] = []
		 			seqparts[i] = cyp7.seq[interval[0]:interval[1]]
		 	else:
		 		raise 'error!'

		for i in xrange(len(edtparts)):
			if edtparts[i] == []:
				edtparts[i] = cyp67diff[i]
				c7parts[i] = '7' * len(seqparts[i]) 
				
 		al.seq  = ' '.join([''.join([c for c in s if c != '-']) for s in seqparts])
 		al.edit = ' '.join([''.join(s) for s in edtparts])
 		al.c7 = ' '.join([''.join(s) for s in c7parts])
		cyp6.alleles[al.name] = al
#print revcomp(''.join(cyp2d6_snps_marked))

def key(x):
	if x[5] == '7': return 10000
	x = x[7:]
	val = float(ord(x[-1])) / 100.0
	if len(x) == 1 or not x[1].isdigit(): val += int(x[0]) 
	elif len(x) < 3 or not x[2].isdigit(): val += int(x[:2]) 
	else: val += int(x[:3]) 
	return val

def out_allele(n, fa, fo):
	seq, aln = cyp6.alleles[n].seq, cyp6.alleles[n].edit
	print >>fa, ">%s len:%d" % (n, len(seq))
	print >>fo, ">%s len:%d" % (n, len(seq))
	width = 80
	seq = revcomp(seq)
	aln = revcomp(aln)
	print >>fa, ''.join([c for c in seq if c != ' '])
	print >>fo, seq #''.join([c for c in seq if c != ' ']) #seq
	print >>fo, aln #''.join([c for c in aln if c != ' ']) #aln
	print >>fo, revcomp(cyp6.alleles[n].c7)

names = sorted(cyp6.alleles.keys(), key=key)
with open(args.out + '.cyp2d6.fa', 'w') as fa:
	with open(args.out + '.cyp2d6.align', 'w') as fo:
		for n in names:
			if n != 'CYP2D7': out_allele(n, fa, fo)
			
with open(args.out + '.cyp2d7.fa', 'w') as fa:
	with open(args.out + '.cyp2d7.align', 'w') as fo:
		out_allele('CYP2D7', fa, fo)
