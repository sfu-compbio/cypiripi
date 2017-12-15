
# Cypiripi has been discontinued and deprecated in favor of [Aldy](https://github.com/inumanag/aldy). 

Aldy is much more capable than Cypiripi and covers all aspect of the original algorithm.
This repository is purely for archival purposes. No support requests will be handled anymore.

----

I am finally releasing source code for this project--- mainly for
archival purposes. It's licensed under
[CRAPL license](http://matt.might.net/articles/crapl/), which best
describes its development process and the resulting source code quality (please
don't ask me how to compile it or how to use those scripts in `util` folder). 

As for the original simulations used in the paper, I lost original FASTQ and SAM files
(they probably disappeared when somebody asked me to clean up my
multi-TB home directory on the shared cluster, and only later I realized that 
the purge was too successful). Sorry for that... 

Anyways, if you really need this, good luck (really)! I highly recommend anybody who needs to get the job done to use
[Aldy](http://aldy.csail.mit.edu), which basically started as an attempt to "fix" 
Cypiripi (hopeless attempt, as it turned out), and which can 
do anything Cypiripi currently does in a way more reproducible and stable manner.

----

# Cypiripi

Cypiripi is a computational tool for *CYP2D6* genotyping.

## Introduction

Software manual is located at http://sfu-compbio.github.io/cypiripi/.

Algorithm details are published at ISMB 2015. Publication can be accessed at http://www.ncbi.nlm.nih.gov/pubmed/26072492.

## F.A.Q.

### What kind of data does Cypiripi currently support?

Cypiripi is designed to work with Solexa/Illumina whole-genome sequencing technologies (e.g. Genome Analyzer, HiSeq and similar sequencers).

We have tested Cypiripi on HiSeq simulations and real data (CEPH trio).
Cypiripi can work with other technologies as long as the following conditions are met:

1. **Coverage is uniform**

	This is the most important requirement. Uneven coverage will break copy-number detection, and predictions will suffer.

2. **Coverage is greater or equal than 40x**

	You can get some results on lower coverage, but Cypiripi is not designed to handle such cases.

3. **Read lengths are greater or equal than 75bp**

4. **Read lengths are consistent**

	Read lengths should be of the same length within the sample.

	If they are not, you can add `--crop X` parameter to the mrsfast and mrfast, where X is the length of smallest read in the sample.

	**Be warned**: cropping will resize all long reads! This means that if you have thousand 100bp reads and one 30bp read, all reads will be considered as 30bp reads. This should be fine if you have reads of size 100 and 101bp, but if you have large variations between read lengths, cropping will not produce expected results.

**Warning**: Use original FASTQ files of the sequencing experiment if possible. Do not use FASTQ files created from post-processed SAM/BAM files. These FASTQ files usually omit some reads due to the mapping misses and might even contain hard-clipped reads.

### How to run Cypiripi on large samples

Current versions of mrfast (2.6.1.0 as the time of writing) will load the whole FASTQ file into the memory. This is impractical for large FASTQ files, unless you have 1 TB of RAM.

If you happen to have a large FASTQ file, it is recommended to split it into several smaller files and invoke the mapping manually. Here is a quick step-by-step guide on how to do that.

1. *Split the FASTQ files*

		split -d -l 20000000 input.fq split/input.fq

	This will generate files small FASTQ files in the directory `split`. Each small file will have 5,000,000 records (totaling 20,000,000 lines).

2. *Run mrfast and mrsfast for every small file*

	Before starting, make sure that you have indexed the reference. If you do not have reference indexes, generate them by running:

		mrsfast --index  reference.cyp2d6.fa
		mrfast  --index  reference.cyp2d7.fa

	First, generate single-end mapping (assuming that you have reference files provided with Cypiripi in the same directory):

		for i in split/input.fq*; do
			mrsfast --search reference.cyp2d6.fa -e 2 --seqcomp --seq ${i} -o ${i}.cyp2D6.sam --disable-sam-header --disable-nohits --threads 4 ;
			mrfast  --search reference.cyp2d7.fa -e 5 --seqcomp --seq ${i} -o ${i}.cyp2D7.sam -u ${i}.cyp2D6.sam.unmap ;
		done

	Afterwards, generate paired-end mapping:

		for i in split/input_1.fq.*; do
			mrsfast --search reference.cyp2d6.fa -e 2 --seqcomp --pe --seq ${i} -o ${i}.cyp2D6.sam.paired --disable-sam-header --disable-nohits --min 250 --max 650 --threads 4 ;
			mrfast  --search reference.cyp2d7.fa -e 5 --seqcomp --pe --seq ${i} -o ${i}.cyp2D7.sam.paired -u {1}.cyp2D6.sam.unmap --min 250 --max 650 ;
		done

	In this example, it is assumed that the distance between the pairs is in the range from 250 to 650. If your sequencing experiment has different insert sizes, feel free to adjust the `--min` and `--max` parameters as necessary.

	These commands can be ran in parallel (e.g. on cluster or via GNU Parallel).

3. 	*Concatenate the output*

	In this step, just concatenate the small `.sam` files into a large one by invoking:

			cat split/input.fq*.cyp2D6.sam split/input.fq*.cyp2D7.sam > input.sam
			cat split/input.fq*.cyp2D6.sam.paired split/input.fq*.cyp2D7.sam.paired > input.sam.paired

4. *Run the Cypiripi*

	Assuming that you have `input.sam` and `input.sam.paired` ready, you can now invoke Cypiripi by running:

		python ./cypiripi.py --fasta reference --sam input.sam --cov X

	where `X` is the coverage of the data.
