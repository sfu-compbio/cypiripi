#!/usr/bin/env python
#786 

import os, sys, errno, argparse, subprocess, fnmatch

def command_line_process():
	parser = argparse.ArgumentParser(description='CYP2D6 genotyping pipeline!')
	parser.add_argument('--fastq', '-f', default='')
	parser.add_argument('--fasta', '-r', default='')
	parser.add_argument('--sam', '-s', default='')
	parser.add_argument('--cov', '-c', default=0)
	parser.add_argument('--err', '-e', default=500)
	parser.add_argument('--norm', '-n', default=0)
	parser.add_argument('--CN', '-C', default=1)
	parser.add_argument('--force', '-!', action='store_true', default=False)
	return parser

def mkdir_p(path):
	try:
		os.makedirs(path)
	except OSError as e:
		if e.errno == errno.EEXIST and os.path.isdir(path):
			pass
		else: raise

def shell(command, control_file='', log=False):
	print >>sys.stderr, '$', command

	if log:
		print "\tExecuting {0}, control file {1}".format(command, control_file)

	if not log or control_file == '':
		p = subprocess.Popen(command, stdout=sys.stdout, stderr=sys.stdout, shell=True)
	else:
		f = open(control_file, 'w')
		p = subprocess.Popen(command, shell=True, stdout=f, stderr=f)

	ret = p.wait()
	if log and control_file != '':
		f.close()
	if ret != 0:
		raise subprocess.CalledProcessError(ret, command)

def is_exec(f):
	return os.path.isfile(f) and os.access(f, os.X_OK)
	
args = command_line_process().parse_args()

def cleanup():
	if not args.force and os.path.isfile(args.sam + '.discard'):
		return args.sam + '.discard'

	fo6 = open(args.sam + '.cleaned.cyp2d6.fq', 'w')
	fo7 = open(args.sam + '.cleaned.cyp2d7.fq', 'w')
	s6 = set()
	s7 = set()
	with open(args.sam) as f:
		for l in f:
			if l[0] == '@': continue
			l = l.strip().split('\t')
			if l[2][:6] == 'CYP2D6':
				fo = fo6
				if l[0] in s6: continue
				s6.add(l[0])
			else:
				fo = fo7
				if l[0] in s7: continue
				s7.add(l[0])
			if fo:
				print >>fo, '@' + l[0]
				print >>fo, l[9]
				print >>fo, '+'
				print >>fo, l[10]
	fo6.close()
	fo7.close()
	if args.force or not os.path.isfile(args.sam + '.cleaned.sam.sort'):
		shell('mrsfast --index  {0}.cyp2d8.fa'.format(args.fasta))
		shell('mrsfast --search {0}.cyp2d8.fa --crop 100 -e 2 --seq {1}.cleaned.cyp2d6.fq -o {1}.cleaned.cyp2d6.sam -u {1}.cleaned.cyp2d6.unmap'.format(args.fasta, args.sam))
		shell('mrfast  --index  {0}.cyp2d8.fa'.format(args.fasta))
		shell('mrfast  --search {0}.cyp2d8.fa --crop 100 -e 5 --seq {1}.cleaned.cyp2d7.fq -o {1}.cleaned.cyp2d7.sam -u {1}.cleaned.cyp2d7.unmap'.format(args.fasta, args.sam))
		shell('cat {0}.cleaned.cyp2d6.sam {0}.cleaned.cyp2d7.sam > {0}.cleaned.sam'.format(args.sam))
		shell('sort -n -k4 {0}.cleaned.sam > {0}.cleaned.sam.sort'.format(args.sam))
	
	fr7s = open(args.sam + '.discard', 'w')
	s6 = set()
	with open(args.sam + '.cleaned.sam.sort') as f:
		X,Y=0,0
		for l in f:
			l = l.strip()
			p = l.split('\t')
			if len(p) < 4:
				continue
			if int(p[3]) < 1490 or int(p[3]) > 5825 - 30:
				if p[0] in s6 and p[0] in s7:
					c = X
					q = 2 + int(args.CN)
					X += 1
				else:
					c = Y
					q = 2
					Y += 1
				if c % q == q - 1:
					print >>fr7s, l
					s6.add(p[0])
	#print X,Y
	fr7s.close()

	shell('rm -rf {0}*cleaned*'.format(args.sam))
	return args.sam + '.discard'

try:
	#if args.force or not os.path.isfile(args.fasta + '.cyp2d6.align'):
	#	shell('fasta/algen.py --out {0}'.format(args.fasta))
	#if args.force or not os.path.isfile(args.fasta + '.combined.fa'):
	#	shell('cat {0}.cyp2d6.fa {0}.cyp2d7.fa > {0}.combined.fa'.format(args.fasta))
	#if args.force or not os.path.isfile(args.fasta + '.combined.align'):
	#	shell('cat {0}.cyp2d6.align {0}.cyp2d7.align > {0}.combined.align'.format(args.fasta))
	if args.force or not os.path.isfile(args.fasta + '.cyp2d6.fa.index'):
		shell('mrsfast --index  {0}.cyp2d6.fa'.format(args.fasta))
	if args.force or not os.path.isfile(args.fasta + '.cyp2d7.fa.index'):
		shell('mrfast  --index  {0}.cyp2d7.fa'.format(args.fasta))
	
	if args.fastq != '':
		args.sam = args.fastq + '.sam'
	if args.force or not os.path.isfile(args.sam):	
		if args.force or not os.path.isfile(args.fastq + '.cyp2D6.sam'):
			shell('mrsfast --search {0}.cyp2d6.fa -e 2 --seqcomp --seq {1} -o {1}.cyp2D6.sam --disable-sam-header --disable-nohits --threads 4'.format(args.fasta, args.fastq))
		if args.force or not os.path.isfile(args.fastq + '.cyp2D7.sam'):
			shell('mrfast  --search {0}.cyp2d7.fa -e 5 --seqcomp --seq {1} -o {1}.cyp2D7.sam -u {1}.cyp2D6.sam.unmap'.format(args.fasta, args.fastq))
		if args.force or not os.path.isfile(args.fastq + '.sam'):
			shell('cat {0}.cyp2D6.sam {0}.cyp2D7.sam > {0}.sam'.format(args.fastq))
		shell('rm -rf {0}.cyp2D*.sam*'.format(args.fastq))
	if args.force or not os.path.isfile(args.sam + '.paired'):	
		if args.force or not os.path.isfile(args.fastq + '.cyp2D6.sam.paired'):
			shell('mrsfast --search {0}.cyp2d6.fa -e 2 --seqcomp --pe --seq {1} -o {1}.cyp2D6.sam.paired --disable-sam-header --disable-nohits --min 250 --max 650 --threads 4'.format(args.fasta, args.fastq))
		if args.force or not os.path.isfile(args.fastq + '.cyp2D7.sam.paired'):
			shell('mrfast  --search {0}.cyp2d7.fa -e 5 --seqcomp --pe --seq {1} -o {1}.cyp2D7.sam.paired -u {1}.cyp2D6.sam.unmap --min 250 --max 650'.format(args.fasta, args.fastq))
		if args.force or not os.path.isfile(args.fastq + '.sam.paired'):
			shell('cat {0}.cyp2D6.sam.paired {0}.cyp2D7.sam.paired > {0}.sam.paired'.format(args.fastq))
		shell('rm -rf {0}.cyp2D*.sam*'.format(args.fastq))
	remove = cleanup()

	if int(args.norm) == 0:
		args.norm = int(float(args.cov) * .4)
	shell('./cypiripi -f {0}.combined.align -s {1} -C {2} -E {3} -T {4}'.format(args.fasta, args.sam, args.cov, args.err, args.norm))
	
except subprocess.CalledProcessError as e:
	print >>sys.stderr, "{0} failed with exit status {2} and message {1}".format(e.cmd, 'N/A', e.returncode)
