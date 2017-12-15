// 786

#include "ilp.h"
using namespace std;

extern unordered_map<int, unordered_set<int>> pairedReads;
extern int optFusion;
unordered_map<int, string> alleleName;

const int exonicBorder = 5600;

inline bool isMultiSNP (const string &s)
{
	return (s.substr(0,3) == "SNP" && string("ACGTN").find(s[4]) == string::npos);
}

int covAtPos(ref_t &ref, int p) 
{
	int cov = 0; 
	for (auto &ab: ref[p].abbrs) 
		if (ab.first.substr(0,3) == "SNP" && !isMultiSNP(ab.first)) 
			cov += ab.second.reads.size();
	return cov;
}

inline int pairedCoverage (int a, unordered_map<int, unordered_set<int>> &s) 
{
	int H = 0; 
	for (auto &x: s[a]) 
		H += in(x, pairedReads[a]);
	return H;				
}

int posThreshold (int p, ref_t& ref, int readLength = 100) 
{
	int t = normalizeThreshold;
	if (revComp(p) < -1620 + readLength)
		t = t * ((9432 - p) / double(readLength));
	return t;
}

inline int inrevComp (int p) 
{
	p = - p - 1619 + 2 + 9430;
	return p + (bool)(p < 0);
}

bool checkShouldDelete (ref_t &reference, string a)
{
	if (a.size() < 9) return true;
	if (a[7] == '4' && isalpha(a[8])) // 4n
		return (reference[inrevComp(1846)].abbrs["SNPCT"].reads.size() < normalizeThreshold);
	if (a[7] == '3' && isalpha(a[8]))
		return (reference[inrevComp(2548)].abbrs["DEL1"].reads.size() < normalizeThreshold);
	if (a[7] == '6' && isalpha(a[8]))
		return (reference[inrevComp(1707)].abbrs["DEL1"].reads.size() < normalizeThreshold);
	if (a == "CYP2D6*20")
		return (reference[inrevComp(1974)].abbrs["INSC"].reads.size() < normalizeThreshold);
	if (a == "CYP2D6*101")
		return (reference[inrevComp(2926)].abbrs["DEL19"].reads.size() < normalizeThreshold);
	return true;
}

vector<pair<vector<string>, int>> extraContraints (ref_t &reference)
{
	vector<pair<vector<string>, int>> v;
	double p;
	if (p = reference[inrevComp(2548)].abbrs["DEL1"].reads.size())
		v.push_back({vector<string>({"CYP2D6*3A", "CYP2D6*3B"}), round(p / expectedCoverage)});
	if (p = reference[inrevComp(1707)].abbrs["DEL1"].reads.size())
		v.push_back({vector<string>({"CYP2D6*6A", "CYP2D6*6B", "CYP2D6*6C", "CYP2D6*6D"}), round(p / expectedCoverage)});
	if (p = reference[inrevComp(1974)].abbrs["INSC"].reads.size())
		v.push_back({vector<string>({"CYP2D6*20"}), round(p / expectedCoverage)});
	return v;
}

bool checkShouldDeleteUniqueSNP (ref_t &reference, string al)
{
	static unordered_map<string, unordered_map<int, string>> abbrs;
	if (!abbrs.size()) {
		abbrs["CYP2D6*101"][2926] = "DEL19";
		abbrs["CYP2D6*3A"][2548] = "DEL1";
	}
	
	for (auto &x: abbrs[al]) {
		int t = posThreshold(inrevComp(x.first), reference);
		if (reference[inrevComp(x.first)].abbrs[x.second].reads.size() > t) {
			L("%s deletion prevented due to %d:%s\n", al.c_str(), x.first, x.second.c_str());
			return false;
		}
	}
	return true;
}

inline vector<string> unwindMultiSNP (const string &snp) 
{
	vector<string> a(2);
	switch (snp[4]) {
		case 'B': a[0] = S("SNP%c%c", snp[3], 'A'), a[1] = S("SNP%c%c", snp[3], 'C'); break; 
		case 'D': a[0] = S("SNP%c%c", snp[3], 'A'), a[1] = S("SNP%c%c", snp[3], 'G'); break; 
		case 'E': a[0] = S("SNP%c%c", snp[3], 'A'), a[1] = S("SNP%c%c", snp[3], 'T'); break; 
		case 'F': a[0] = S("SNP%c%c", snp[3], 'C'), a[1] = S("SNP%c%c", snp[3], 'G'); break; 
		case 'H': a[0] = S("SNP%c%c", snp[3], 'C'), a[1] = S("SNP%c%c", snp[3], 'T'); break;
		case 'I': a[0] = S("SNP%c%c", snp[3], 'G'), a[1] = S("SNP%c%c", snp[3], 'T'); break;
		default: McAssert(0, "Key %c weird, %d in %s", snp[4], snp.c_str());
	}
	return a;
}

void setUpMultiSNPs (ref_t &reference, alleles_t &alleles)
{
	for (auto &pos: reference) for (auto &ab: pos.abbrs) {
		if (isMultiSNP(ab.first)) {
			auto a = unwindMultiSNP(ab.first);
			DO(2) {
				pos.abbrs[a[_]].reads.insert(ab.second.reads.begin(), ab.second.reads.end());
				for (auto &c: ab.second.coverage) 
					pos.abbrs[a[_]].coverage[c.first].insert(c.second.begin(), c.second.end());
				pos.abbrs[a[_]].printStatus = true; 
			}
			ab.second.printStatus = false; 
		}
	}
}


// TODO fuzije!!!
void resolveMultiSNPsFor2D7 (ref_t &reference, alleles_t &alleles)
{
	for (auto &ref: reference) if (ref.has2D7abbr) {
		int p = int(&ref - &(reference[0]));
		if (revComp(p) > exonicBorder) continue;

		for (auto &ab: ref.abbrs) if (isMultiSNP(ab.first)) {
			auto pp = unwindMultiSNP(ab.first);
			
			vector<int> pm(2);
			for (int i = 0; i < 2; i++) {
				auto abb = ref.abbrs[pp[i]];
					
				int max = 0, maxid = -1;
				for (auto &ax: abb.coverage) 
					if (ax.second.size() > max || (ax.second.size() == max && ax.first == alleles["CYP2D7"].id))
						max = ax.second.size(), maxid = ax.first;
				
				pm[i] = (maxid == alleles["CYP2D7"].id) ? max : -1;
			}
			int todel = 0;
			if (pm[0] != -1 && pm[1] != -1) 
				todel = (pm[0] == pm[1]) ? -1 : (pm[0] > pm[1] ? 1 : 0); 
			else if (pm[0] != -1) todel = 1;
			else if (pm[1] != -1) todel = 0;
			else todel = -1;
			
			if (todel != -1) {
				L("D. MULT:%d\t Removed CYP2D7 from %s (from %s), because %d < %d\n", revComp(p), pp[todel].c_str(), ab.first.c_str(), pm[todel], pm[todel?0:1]);
				ref.abbrs[pp[todel]].alleles.unset(alleles["CYP2D7"].id);
				ref.abbrs[pp[todel]].coverage.erase(alleles["CYP2D7"].id);
				ref.abbrs.erase(ab.first);
				break;
			}
		}
	}
}

bool is2D7derivative (const string &s) 
{
	static unordered_set<string> derivatives({
		"CYP2D6*2A",	"CYP2D6*2M",	"CYP2D6*4N",	"CYP2D6*10D",	"CYP2D6*13A",
		"CYP2D6*13B",	"CYP2D6*13Ca",	"CYP2D6*13Cb",	"CYP2D6*13D",	"CYP2D6*13E",	
		"CYP2D6*13Fa",	"CYP2D6*13Fb",	"CYP2D6*13G",	"CYP2D6*13H",	"CYP2D6*31",	
		"CYP2D6*35B",	"CYP2D6*36T",	"CYP2D6*36S",	"CYP2D6*41",	"CYP2D6*51",	
		"CYP2D6*56A",	"CYP2D6*57",	"CYP2D6*58",	"CYP2D6*61", 	"CYP2D6*63",	
		"CYP2D6*68A",	"CYP2D6*82",	"CYP2D6*83",	"CYP2D7"
	});

	return in(s, derivatives);
}

void prepare2D7 (ref_t &reference, alleles_t &alleles) 
{
	if (!optFusion) 
		return;

	for (auto &al: alleles) {
		for (int i = 0; i < al.second.c7.size(); i++) if (al.second.c7[i] == '7') {
			int p = i + al.second.offsets[i];
			
			for (auto &abbr: reference[p].abbrs) if (abbr.second.alleles.get(alleles["CYP2D7"].id))
			{
				for (auto &c: abbr.second.coverage)
					if (is2D7derivative(alleleName[c.first])) {
						for (auto &x: abbr.second.coverage[alleles["CYP2D7"].id]) {
							c.second.insert(x);
							if (in(x, pairedReads[alleles["CYP2D7"].id]))
								pairedReads[c.first].insert(x);
						}
					}
			}
		}
	}
}

void cleanTentatives (ref_t &reference, alleles_t &alleles)
{
	for (auto &pos: reference) 
	{
		int p = &pos - &(reference[0]);
		if (revComp(p) > exonicBorder) continue;
		int t = posThreshold(p, reference);
		// A. get rid of tentatives! per abbr! 
		for (auto &ab: pos.abbrs) if (ab.second.tentative.size()) {
			L("- A. TENT@%d\t", revComp(p));
			if (ab.second.tentative.size() >= ab.second.reads.size()) {
				L("Remove tentative indel %s because %d > %d\n", ab.first.c_str(), ab.second.tentative.size() , ab.second.reads.size());
				for (auto &c: ab.second.coverage)
					c.second.clear();
				ab.second.reads.clear();
			}
			else {
				L("Add tentative indel %s because %d <= %d\n", ab.first.c_str(), ab.second.tentative.size() , ab.second.reads.size());
				ab.second.reads.insert(ab.second.tentative.begin(), ab.second.tentative.end());
			}
		}
	}
}

void normalize (ref_t &reference, alleles_t &alleles, int readLength = 100) 
{
	for (auto &pos: reference) {
		int p = &pos - &(reference[0]);
		if (revComp(p) > exonicBorder) continue;
		int t = posThreshold(p, reference, readLength);

		// B. clean-up of the non-2D7 paired
		// remove per-allele data
		if (pairedReads.size()) for (auto &ab: pos.abbrs) {
			// int H7 = in(alleles["CYP2D7"].id, ab.second.coverage) ? pairedCoverage(alleles["CYP2D7"].id, ab.second.coverage) : 0;
			for (auto &c: ab.second.coverage) if (c.first != alleles["CYP2D7"].id) { 
				int H = pairedCoverage(c.first, ab.second.coverage); 
				if (H <= t / 2) {
					// L("B. PAIR@%d\tRemove %s from %s because h[%d] < t/2, threshold %d\n", revComp(p), ab.first.c_str(), alleleName[c.first].c_str(), H, t);
					c.second.clear();
				}
			}
		}
		
		// 3. remove per-read lo-cov
		for (auto &abb: pos.abbrs) {
			if (abb.first != S("SNP%c%c",pos.ref,pos.ref) && abb.second.reads.size() <= t) {
				// L("C. THRS@%d\tRemove %s because %d < t, threshold %d\n", revComp(p), abb.first.c_str(),  abb.second.reads.size(), t);
				abb.second.reads.clear();
			}
		}
	}
}

void cleanSubsetInstances (ref_t &reference, alleles_t &alleles)
{
	for (auto &ref: reference) {
		int p = int(&ref - &(reference[0]));
		if (revComp(p) > exonicBorder) continue;
		
		for (auto &ab: ref.abbrs) {
			bool is2D7abbr = in(alleles["CYP2D7"].id, ab.second.coverage);
			for (auto &c: ab.second.coverage) {
				int superset = -1;
				for (auto &cx: ab.second.coverage) if (c.first != cx.first && cx.second.size() > c.second.size()) {
					bool isSuperset = 1;
					for(auto &x: c.second) 
						if (!in(x, cx.second)) { isSuperset = 0; break; }
					if (isSuperset)	{ superset = cx.first; break; }
				}

				if (superset != -1 
						&& !(is2D7abbr && is2D7derivative(alleleName[superset]))) // superset != alleles["CYP2D7"].id)) { 
				{ 
					L("KIMCHI ALL %s removed at %d:%s due to subset of %s!!!\n", alleleName[c.first].c_str(),revComp(p),ab.first.c_str(), alleleName[superset].c_str()); 
					ab.second.alleles.unset(c.first);
				}
			}
			for (auto &a: alleles) if (!ab.second.alleles.get(a.second.id))
				ab.second.coverage.erase(a.second.id);
		}	
	}
}

void deleteAllelesFromSet (ref_t &reference, alleles_t &alleles, const unordered_set<int>& toDelete)
{
	for (auto &ref: reference) {
		int p = int(&ref - &(reference[0]));
		if (revComp(p) > exonicBorder) continue;
		int t = posThreshold(p, reference);

		unordered_set<string> toDeleteAbb;
		for (auto &ab: ref.abbrs) {
			for (auto &an: toDelete) 
				ab.second.alleles.unset(an);

			for (auto &c: ab.second.coverage) {
				int H = pairedCoverage(c.first, ab.second.coverage); // 0; for (auto &x: c.second) H += (pairedReads[c.first].find(x) != pairedReads[c.first].end());
				if (c.second.size() <= t || H <= t/2) 
					ab.second.alleles.unset(c.first);
			}
			for (auto &a: alleles) if (!ab.second.alleles.get(a.second.id))
				ab.second.coverage.erase(a.second.id);

			int sum = 0;
			for (auto &al: alleles)
				sum += ab.second.alleles.get(al.second.id);
			if (!sum) toDeleteAbb.insert(ab.first);
		}
		for (auto &td: toDeleteAbb)
			ref.abbrs.erase(td);
	}
	for (auto &an: toDelete)  
		alleles.erase(alleleName[an]);
	//alleles.erase("CYP2D7X");
}

void cleanUnusedAlleles (ref_t &reference, alleles_t &alleles, bool remove2D7 = false) 
{
	DO(2) {
		cleanTentatives(reference, alleles);
		normalize(reference, alleles);
		resolveMultiSNPsFor2D7(reference, alleles);
		prepare2D7(reference, alleles);

		unordered_set<int> toDelete;
		for (auto &al: alleles) {
			if (!optFusion && (al.first != "CYP2D7" && is2D7derivative(al.first))) {
				toDelete.insert(al.second.id);
				continue;
			}

			vector<int> cyp7InvalidAbbrs;
			for (auto &ab: al.second.abbrs) {
				string abr = ab.second;
				if (isMultiSNP(abr))
					continue;

				int p = ab.first + al.second.offsets[ab.first];
				int t = posThreshold(p, reference);
				if ((!reference[p].abbrs[abr].reads.size() || reference[p].abbrs[abr].coverage[al.second.id].size() <= t) 
						&& reference[p].catab() > 0) 
				{
/*!!!*/					if (is2D7derivative(al.first) && al.second.c7[ab.first] == '7') {
						//L("-- GUMICA %d_%d:%s %d %d\n",p,revComp(p),ab.second.c_str(), reference[p].abbrs.size(), reference[p].abbrs[abr].coverage[al.second.id].size());
						cyp7InvalidAbbrs.push_back(ab.first); 
					}
					else {
						L("E. %s removed due to %d (%d) at %d:%s, thr %d, cc = %d\n", al.first.c_str(), reference[p].abbrs[abr].reads.size(), reference[p].abbrs[abr].coverage[al.second.id].size(),
							revComp(p), abr.c_str(), t, covAtPos(reference,p));
						if (checkShouldDeleteUniqueSNP(reference, al.first)) 
							toDelete.insert(al.second.id); 
						break;
					}
 				} 
			}
			for (auto &x: cyp7InvalidAbbrs)
				al.second.abbrs.erase(x);
		}

		deleteAllelesFromSet(reference, alleles, toDelete);
		if (!_) cleanSubsetInstances(reference, alleles);
	}
}

void removeCYP2D7 (ref_t &reference, alleles_t &alleles, int cnv, int expectedCN) 
{
	int c7id = alleles["CYP2D7"].id;

	// Remove multi-SNPs for CYP2D7 only!
	/*if (!optFusion)*/ for (auto &pos: reference) for (auto &abb: pos.abbrs) {
		if (isMultiSNP(abb.first)) {
			auto a = unwindMultiSNP(abb.first);
			
			auto torem = a[0];
			if (pos.abbrs[a[1]].reads.size() < pos.abbrs[torem].reads.size())
				torem = a[1];

			pos.abbrs[torem].alleles.unset(c7id);
			pos.abbrs[torem].coverage.erase(c7id);
		}
	}

	// Find all CYP2D7-unique reads!
	unordered_set<int> cyp7todelete;
	for (auto &ref: reference) {
		int p = int(&ref - &(reference[0]));
		if (revComp(p) > exonicBorder) continue;
		for (auto &ab: ref.abbrs) {
			unordered_set<int> reads;
			for (auto &c: ab.second.coverage) if (c.first != c7id)
				reads.insert(c.second.begin(), c.second.end());
			if (in(c7id, ab.second.coverage)) {
				for (auto &x: ab.second.coverage[c7id]) 
					if (!in(x, reads))
						cyp7todelete.insert(x);
			}
		}
	}

	for (auto &ref: reference) {
		int p = int(&ref - &(reference[0]));
		if (revComp(p) > exonicBorder) continue;

		double currentCoverage = covAtPos(reference, int(&ref - &(reference[0]))) / double(expectedCN);
		currentCoverage = int(round(currentCoverage) + .1);
		
		for (auto &ab: ref.abbrs) {
			unordered_set<int> reads;
			for (auto &c: ab.second.coverage) if (c.first != c7id)
				reads.insert(c.second.begin(), c.second.end());
			
			if (in(c7id, ab.second.coverage)) {
				vector<int> todel, todel_tmp;
				for (auto &x: ab.second.coverage[c7id]) 
					if (!in(x, reads) || in(x, cyp7todelete))
						todel.push_back(x);
					else
						todel_tmp.push_back(x);
				random_shuffle(todel_tmp.begin(), todel_tmp.end());

				int cnt = min((int)(todel.size() + todel_tmp.size()), (int)currentCoverage * cnv);
				while(todel.size() < cnt) {
					todel.push_back(todel_tmp[todel_tmp.size() - 1]);
					todel_tmp.pop_back();
				}
				
				for (int i = 0; i < todel.size(); i++) {
					ab.second.reads.erase(todel[i]);
					for (auto &c: ab.second.coverage)
						c.second.erase(todel[i]);
				}
				L("SAYONARA [cnt=%d,ccv=%.1lf] %d cleaned from %d:%s\n",cnt,currentCoverage,todel.size(),revComp(p),ab.first.c_str());
				ab.second.alleles.unset(c7id);
				ab.second.coverage.erase(c7id);
			}
			
			for (auto &c: ab.second.coverage) {
				for (auto x: cyp7todelete)
					c.second.erase(x);
			}

			reads.clear();
			for (auto &c: ab.second.coverage) 
				reads.insert(c.second.begin(), c.second.end());
			ab.second.reads = reads;
		}
	}


	// Remove subal
	normalize(reference, alleles);
	
	unordered_set<int> toDelete;
	for (auto &al: alleles) {
		for (auto &ab: al.second.abbrs) if (!isMultiSNP(ab.second)) {
			string abr = ab.second;
			int p = ab.first + al.second.offsets[ab.first];
			int t = posThreshold(p, reference);
			if ((!reference[p].abbrs[abr].reads.size() || reference[p].abbrs[abr].coverage[al.second.id].size() <= t) 
				&& reference[p].catab() > 0) 
			{	
				if (checkShouldDelete(reference, al.first)) {
					L("EX. %s removed due to %d (%d) at %d:%s, thr %d\n", al.first.c_str(), reference[p].abbrs[abr].reads.size(), reference[p].abbrs[abr].coverage[al.second.id].size(),
						revComp(p), abr.c_str(), t);
					toDelete.insert(al.second.id); 
				}
				//break;
			}
		}
		toDelete.insert(alleles["CYP2D7"].id);
		if (optFusion && al.first.size() > 9 && al.first.substr(0,9) == "CYP2D6*13")
			toDelete.insert(al.second.id);
	}

	deleteAllelesFromSet(reference, alleles, toDelete);
}


void printAlleles (ref_t &reference, alleles_t &alleles, int columnWidth = 4)
{
	// Estabilish only used SNPs
	for (auto &al: alleles) {
		vector<string> entries;
		L("+ %s:\n", al.first.c_str());
		for (auto &ab: al.second.abbrs) {
			int p = ab.first + al.second.offsets[ab.first];
			auto &pos = reference[p];

			string refa = S("SNP%c%c", pos.ref, pos.ref);
			if (ab.second.substr(0,3) == "INS") refa = "INS0";

			string en = pos.has2D7abbr ? "." : "_";
			if (ab.second == refa)
				en += S("%5d:SNP*=%c   -/%-3d", revComp(p), revComp(pos.ref), pos.abbrs[refa].reads.size());
			else if (ab.second[0] == 'S' && string("ACGTN").find(ab.second[4]) != string::npos)
				en += S("%5d:SNP*%c%c %3d/%-3d", revComp(p), revComp(pos.ref), revComp(ab.second[4]), pos.abbrs[ab.second].reads.size(), pos.abbrs[refa].reads.size());
			else if (ab.second[0] == 'S' && string("ACGTN").find(ab.second[4]) == string::npos) {
				auto pp = unwindMultiSNP(ab.second);
				en += S("%5d:SNP*%c%c %3d/%-3d", revComp(p), revComp(pos.ref), revComp(pp[0][4]), pos.abbrs[pp[0]].reads.size(), pos.abbrs[refa].reads.size());
				entries.push_back(en + "\t");
				en = pos.has2D7abbr ? "." : "_";
				en += S("%5d:SNP*%c%c %3d/%-3d", revComp(p), revComp(pos.ref), revComp(pp[1][4]), pos.abbrs[pp[1]].reads.size(), pos.abbrs[refa].reads.size());
			}
			else if (ab.second[0] == 'I') {
				int sum = 0;
				for (auto &x: reference[p + 1].abbrs)
					sum += x.second.reads.size();
				sum -= pos.abbrs[ab.second].reads.size();
				en += S("%5d:INS*%s %3d/%-3d", revComp(p), revComp(ab.second.substr(3)).c_str(), pos.abbrs[ab.second].reads.size(), sum);
			}
			else
				en += S("%5d:DEL*%s %3d/%-3d", revComp(p), ab.second.substr(3).c_str(), pos.abbrs[ab.second].reads.size(), reference[p-1].abbrs[S("SNP%c%c",reference[p-1].ref,reference[p-1].ref)].reads.size());
			entries.push_back(en + "\t");
		}
		for (int i = entries.size() - 1; i >= 0; i--) {
			if ((entries.size() - 1 - i) % columnWidth == 0) L("   ");
			L("%s\t", entries[i].c_str());
			if ((entries.size() - 1 - i) % columnWidth == (columnWidth - 1)) L("\n");
		}
		if (!entries.size() || (entries.size() - 1) % columnWidth != (columnWidth - 1)) L("\n");
	}
	L("\n%s\n\n", string(80, '*').c_str());
	for (auto &ref: reference) {
		vector<string> entries;
		string refa = S("SNP%c%c", ref.ref, ref.ref);
		int p = &ref - &(reference[0]);

		bool op=0;
		for (auto &ab: ref.abbrs) if (ab.second.printStatus) {op=1;break;}

		for (auto &ab: ref.abbrs) if (ab.second.printStatus||(ab.first==refa&&op))
		{
			string en;
			if (ab.first == refa)
				en = S("%5d:SNP*%c%c", revComp(p), revComp(ref.ref), revComp(ref.ref));
			else if (ab.first[0] == 'S')
				en = S("%5d:SNP*%c%c", revComp(p), revComp(ref.ref), revComp(ab.first[4]));
			else if (ab.first[0] == 'I')
				en = S("%5d:INS*%s", revComp(p), revComp(ab.first.substr(3)).c_str());
			else
				en = S("%5d:DEL*%s", revComp(p), ab.first.substr(3).c_str());


			L("#%s\t\t", en.c_str());
			L("%d => ",ab.second.reads.size()); for (auto &a: alleles) if (ab.second.alleles.get(a.second.id)) L(" %s", a.first.c_str()); L(" )\n");
			for (auto &c: ab.second.coverage) 
				L("%s[%3d(%3d)]  ", alleleName[c.first].c_str(), c.second.size(), pairedCoverage(c.first, ab.second.coverage));
			L("\n");
		}
	}

	// L("\n%s\n\n", string(80, '~').c_str());
	// for (auto &ref: reference) {
	// 	L("%4d_%d: ", revComp(int(&ref - &(reference[0]))), int(&ref - &(reference[0])));
	// 	for (auto &ab: ref.abbrs)
	// 		L("(%s %d) ", ab.first.c_str(), ab.second.reads.size());
	// 	L("\n");
	// }
}

unordered_map<string, IloNumExprArg> setUpModel 
	(ref_t &reference, alleles_t &alleles, IloEnv &env, IloModel &model, int maxAllele = 20, int expectedCoverage = 40, int expectedCN = 0, int max2D7 = 2) 
{
	unordered_map<string, IloNumExprArg> avar;
	
	unordered_set<int> posToAvoidNormalSNPs;
	// move supersnips
	for (auto &pos: reference) for (auto &abb: pos.abbrs) {
		if (isMultiSNP(abb.first)) {
			auto a = unwindMultiSNP(abb.first);
			abb.second.reads.clear();
			abb.second.coverage.clear();
			abb.second.alleles = bitset(0);
			DO(2) {
				abb.second.reads.insert(pos.abbrs[a[_]].reads.begin(), pos.abbrs[a[_]].reads.end());
				for (auto &c: pos.abbrs[a[_]].coverage) 
					abb.second.coverage[c.first].insert(c.second.begin(), c.second.end());	
				abb.second.alleles.add(pos.abbrs[a[_]].alleles);
				posToAvoidNormalSNPs.insert(int(&pos - &(reference[0])));
				//pos.abbrs.erase(a[_]);
			}
			break;
		}
	}
	//if (in(string("CYP2D7X"), alleles))
		
	if (!max2D7)
		posToAvoidNormalSNPs.clear();

	IloExpr sumOfCyps(env);
	vector<IloIntVar> cyp7vars;
	for (auto &al: alleles) {
		auto a = IloIntVar(env, 0, maxAllele, S("A_%s", al.first == "CYP2D7" ? "2D7" : al.first.substr(7).c_str()).c_str());
		avar[al.first] = a; 
		sumOfCyps += a;
	}
	if (expectedCN)
		model.add(sumOfCyps == expectedCN);
	if (in(string("CYP2D7"), avar))
		model.add(avar["CYP2D7"] <= max2D7);

	auto extra = extraContraints(reference);
	for (auto x: extra) if (x.second)
	{
		IloExpr e(env);
		for (auto v: x.first) if (in(v, avar))
			e += avar[v];
		L("EXTRA %d <=", x.second); for (auto v: x.first)  if (in(v, avar)) L(" + %s", v.c_str()); L("\n");
		model.add(e >= x.second);
	}

	IloExpr errSum(env);
	for (auto &ref: reference) if (ref.catab()) 
	{
		int p = int(&ref - &(reference[0]));
		if (revComp(p) > exonicBorder) continue;
		
		double currentCoverage = expectedCN ? covAtPos(reference, p) / double(expectedCN)
											: expectedCoverage; 
		
		for (auto &ab: ref.abbrs) {
			if (in(p, posToAvoidNormalSNPs) && isMultiSNP(ab.first))
				continue;
		// FIX
			if (!max2D7 && isMultiSNP(ab.first))
				continue;

			///	1. cov(pos) = einrevComp(pos) + l * sum(a)
			IloExpr alsum(env);
			for (auto &n: alleles) if (ab.second.alleles.get(n.second.id)) 
				alsum += avar[n.first];

			IloNumVar err(env, -errLimit, errLimit, S("E_%d_%s", revComp(int(&ref - &(reference[0]))), ab.first.c_str()).c_str());
			model.add(ab.second.reads.size() == currentCoverage * alsum + err);
			errSum += IloAbs(err);
			
			///	2. l * a >= uniq(a) + einrevComp(a)
			/*for (auto &al: alleles) if (ab.second.alleles.get(al.second.id)) {
				int c = ab.second.coverage[al.second.id].size();
				IloNumVar einrevComp(env, 0, errLimit, S("EX_%d_", revComp(int(&ref - &(reference[0])))).c_str());
				//model.add(currentCoverage * avar[al.first] <= c + err);

				IloNumVar err2(env, 0, errLimit, S("ERPE_%d_", revComp(int(&ref - &(reference[0])))).c_str());
				int H(0); for (auto &x: ab.second.coverage[al.second.id]) H += (pairedReads[al.second.id].find(x) != pairedReads[al.second.id].end());
				model.add((currentCoverage/2) * avar[al.first] <= H + err2);
				
				errSum +=  err2;
			}*/
		}
	}

	model.add(IloMinimize(env, errSum));
	return avar;
}

void solve (ref_t &reference, alleles_t &alleles, const string &fn)
{
	for (auto &al: alleles) {
		if (al.first != "CYP2D7X" && al.first != "CYP2D6*9X") 
			alleleName[al.second.id] = al.first;
	}
	alleles.erase("CYP2D7X");
	alleles.erase("CYP2D6*9X");

	setUpMultiSNPs(reference, alleles);
	cleanUnusedAlleles(reference, alleles);
	printAlleles(reference, alleles, 3);
	
	int expectedCN = 0;
	IloEnv env;	


	set<string> banned;
	vector<map<string, int>> prevSol;
	double prevSolVal = -1;
	int numSols = 0;
	int cyp7sol = 0;

	for (int cplexPhase = 0; ; cplexPhase++)
	{
		try
		{
		//	if (cplexPhase >= 2) for (auto &a: banned)
		//		if (av.find(a) != av.end()) { model.add(av[a] == 0);  L("pass %d ban %s\t",cplexPhase, a.c_str()); }
		//	L("---\n");
			IloModel model(env);
			auto av = setUpModel(reference, alleles, env, model, 20, expectedCoverage, expectedCN, cplexPhase == 0 ? 2 : 0);

			if (cplexPhase >= 2) {
				for (auto &x: prevSol) {
					IloExpr ex(env); int c = 0;
					for (auto &a: x) if (a.second) {
						ex += av[a.first];
						c += a.second;
					//	E("+%s",a.first.c_str());
					}
					//E("<=%d\t\t",c-1);
					model.add(ex <= c - 1);
				}
			}

			IloCplex cplex(model);
			cplex.setParam(IloCplex::TiLim, 30);
			//cplex.setParam(IloCplex::TreLim, 2048);
			//cplex.setParam(IloCplex::WorkMem, 2048);
			cplex.setOut(env.getNullStream());
			cplex.setWarning(env.getNullStream());
			cplex.exportModel((fn + ".lp").c_str());
			
			cplex.solve();

			if (cplexPhase >= 2 && fabs(cplex.getObjValue() - prevSolVal) > 0.1)
				break;
			prevSolVal = cplex.getObjValue();

			L("CPLEXer: status %d, optimum value %.2lf, memory %'.1fM, results:\n", cplex.getCplexStatus(), cplex.getObjValue(), env.getTotalMemoryUsage() / (1024 * 1024.0));
			E("CPLEX = %03d, OPT = %7.1lf ", cplex.getCplexStatus(), cplex.getObjValue());
			
			map<string, int> sol;
			for (auto &a: av) {
				int alval = 0;
				try { alval = IloRound(cplex.getValue(a.second)); }
				catch (...) { alval = 0; }
				sol[a.first] = alval;
			}
			if (cplexPhase >= 1)
				prevSol.push_back(sol);

			E("[");

			for (auto &a: sol) if (a.second) {
				if (!cplexPhase) { expectedCN += a.second; }
				string n = a.first;
				// if (n == "CYP2D6*9" && a.second > 2)
				// 	n = "CYP2D6*9X";
				// if (n == "CYP2D6*9X" && a.second <= 2)
				// 	n = "CYP2D6*9";

				if (n.substr(0,6) == "CYP2D7")
					n = "2D7";
				else 
					n = n.substr(6);
				if (a.second > 1) 
					E(" %d×%s", a.second, n.c_str());
				else 
					E(" %s", n.c_str());
				L("%11s: value %d\n", n.c_str(), a.second);
			}
			if (cplexPhase && cyp7sol) {
				if (cyp7sol > 1) 
					E(" %d×2D7", cyp7sol);
				else 
					E(" 2D7");
			}
			E(" ] \n");
			if (cplexPhase >= 1) numSols++;

			if (!cplexPhase) {
				cyp7sol = sol["CYP2D7"];
				E("Solutions:\n");
			}
			if (!cplexPhase && cyp7sol) {
				removeCYP2D7(reference, alleles, cyp7sol, expectedCN);
				L("\n\nMEGAWAHA\n");
				printAlleles(reference, alleles, 3);
				expectedCN -= cyp7sol;
				// break;
			}
		}
		catch (IloException &ex) {
			if (cplexPhase < 2)
				E("CPLEX Error; %s\n", ex.getMessage());
			break;
		}
	}

	if (numSols > 1) E("WARNING: Multiple optimal arrangements found\n");
	E("Done!\n\n");

	env.end();
}