#ifndef TypesH
#define TypesH

#include "common.h"
using namespace std;

#define DO(i) for(int _=0;_<(i);_++)

template<int W>
struct bitset_t {
	vector<uint64_t> bit;

	bitset_t (char p): bit(W/64+1,p?ULLONG_MAX:0) {}
	bitset_t (): bitset_t(0) {}

	bool get (int x) const {
		return (bit[x>>6]>>(x%64)) & 1;
	}

	void set (int x) {
		bit[x>>6] |= (1ULL<<(x%64));
	}

	void unset (int x) {
		bit[x>>6] &= ~(1ULL<<(x%64));
	}

	void intersect (const bitset_t &b) {
		if (b==bitset_t())return;
		DO(W/64+1) bit[_] &= b.bit[_];
	}

	void add (const bitset_t &b) {
		DO(W/64+1) bit[_] |= b.bit[_];
	}

	bool operator< (const bitset_t &b) const {
		DO(W/64+1) if (bit[_] != b.bit[_]) return (bit[_] < b.bit[_]);
		return false;
	}

	bool operator== (const bitset_t &b) const {
		DO(W/64+1) if (bit[_] != b.bit[_]) return false;
		return true;
	}
};
typedef bitset_t<200> bitset;
namespace std {
	template <>
	struct hash<bitset> {
		size_t operator()(const bitset& k) const {
			size_t h = 0;
			for (auto x: k.bit)
				h ^= x;
			return h;
		}
	};
}

struct abbr_t {
	unordered_set<int> reads;
	unordered_set<int> tentative;
	
	bitset alleles;
	unordered_map<int, unordered_set<int>> coverage;
	bool printStatus;
	
	abbr_t (): printStatus(0) {
	}

	void dump(FILE *fo) {
		fprintf(fo, "%d ", reads.size()); for (auto i: reads) fprintf(fo, "%d ", i);
		fprintf(fo, "%d ", tentative.size()); for (auto i: tentative) fprintf(fo, "%d ", i);
		fprintf(fo, "%d ", coverage.size()); 
		for (auto i: coverage) { 
			fprintf(fo, "%d %d ", i.first, i.second.size());
			for (auto j: i.second) fprintf(fo, "%d ", j);
		}
	}
	void load(FILE *fi) {
		int sz, i, j, k;
		fscanf(fi, "%d", &sz); while (sz--) { fscanf(fi, "%d ", &i); reads.insert(i); }
		fscanf(fi, "%d", &sz); while (sz--) { fscanf(fi, "%d ", &i); tentative.insert(i); }
		fscanf(fi, "%d", &sz); while (sz--) { 
			fscanf(fi, "%d %d ", &i, &j);
			while (j--) { fscanf(fi, "%d ", &k); coverage[i].insert(k); }
		}
	}
};

struct pos_t {
	char ref;
	bool has2D7abbr;
	unordered_map<string, abbr_t> abbrs;

	pos_t (): ref(0), has2D7abbr(0) {
	}

	int catab() { return abbrs.size()-1; }
};
typedef vector<pos_t> ref_t;

struct allele_t {
	int id;
	string seq;
	string edit;
	string c7;
	vector<int> offsets; // A[x] -> where x+A[x] maps to the CYP2D6*1
	map<int, string> abbrs; // <A[x], abbr> -> for X map to
};
typedef unordered_map<string, allele_t> alleles_t;

extern int normalizeThreshold;
extern int expectedCoverage;
extern int errLimit;

#endif