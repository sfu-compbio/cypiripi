/***	786		***/

#include "common.h"
#include "types.h"
#include "reference.h"
#include "sam.h"
#include "ilp.h"

using namespace std;

int normalizeThreshold = 0;
int expectedCoverage = 20;
int errLimit = 0;
string mulFile = "";
string parFile = "";
int optFusion = 0;

int main (int argc, char **argv) {
	setlocale(LC_ALL, NULL);

	string samFile = "../maps/NA12828/rremap";
	string refFile = "/home/inumanag/snjofavac/git/2d6/fasta/CYP2D6.alleles.align";
	string excFile = "";

	int opt;
	while ((opt = getopt(argc, argv, "T:C:E:s:f:x:X:m:p:F")) != EOF)
		switch (opt) {
			case 'T': normalizeThreshold = atoi(optarg); break;
			case 'C': expectedCoverage = atoi(optarg); break;
			case 'E': errLimit = atoi(optarg); break;
			case 's': samFile = optarg; break;
			case 'f': refFile = optarg; break;
			case 'x': excFile = optarg; excFile = "-" + excFile; break;
			case 'X': excFile = optarg; excFile = "+" + excFile; break;
			case 'm': mulFile = optarg; break;
			case 'p': parFile = optarg; break;
			case 'F': optFusion = 1; break;
			default: McAssert(0, "Wrong parameter supplied: %c", opt);
		}


	if (!normalizeThreshold) normalizeThreshold = expectedCoverage / 3 + 1;
	if (!errLimit) errLimit = 3 * expectedCoverage;

	string sy = samFile + ".cyplog";
	freopen(sy.c_str(), "w", stdout);

	excFile = "-" + samFile + ".discard";
	parFile = samFile + ".paired";

	L("The CYP2D Solver (alpha), powered by The Meho Engine v1.0\nParams: -T(normalization) %d, -C(coverage) %d, -E(error) %d\n", normalizeThreshold, expectedCoverage, errLimit);
	ref_t ref;
	alleles_t all;

	zaman_last();
	readReference(refFile, ref, all);
	readSAMFile(samFile, ref, all, excFile); 
	L("\tReading %d seconds\n", zaman_last());
	solve(ref, all, samFile);
	
	return 0;
}