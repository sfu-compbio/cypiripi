/// 786

#include "reference.h"
using namespace std;

void testParsing (const string &n, const string &ref, string al, const vector<int> offsets, const map<int, string> &abbrs) 
{
	if (n == "CYP2D6*1") al = ref;
	for (int i = 0; i < al.size(); i++) {
		McAssert(i + offsets[i] < ref.size(), "%s: offset bound consisitency: %d+%d vs %d", n.c_str(), i,offsets[i], ref.size());
		if (isupper(al[i]) && ref[i + offsets[i]] != al[i]) {
			McAssert(abbrs.find(i) != abbrs.end(), "%s: snp %c->%c validity at %d (%d)", n.c_str(), ref[i + offsets[i]], al[i], i, i+offsets[i]);
		}
	}
}

void readReference (const string &path, ref_t &reference, alleles_t &alleles) 
{
	/* Reference format: 
	 * 	>NAME 	CYP2D[0-9]*[A-Za-z0-9]+
	 *	SEQ 	(I0 E1 I1 ... En In), each [ACGT]+
	 *	EDIT	(I0 E1 I1 ... En In), each [=-ACGTacgt]+
	 */

	ifstream fin(path.c_str());
	string l;
	int id = 0;
	while (getline(fin, l)) {
		string name = "", seq = "", edit = "";
		for (int i = 1; i < l.size() && l[i] != ' '; i++)
			name += l[i];
		getline(fin, l);
		seq = ""; for (char c: l) if (c != ' ') seq += c;
		getline(fin, l);
		edit = ""; for (char c: l) if (c != ' ') edit += c;
		getline(fin, l);
		string c7 = ""; for (char c: l) if (c != ' ') c7 += c;
		
		if (name == "CYP2D6*1") {
			reference.resize(seq.size() + 1);
			for (int i = 0; i < seq.size(); i++) 
				reference[i].ref = seq[i];
		}
		//L("# %s:\t", name.c_str());

		vector<int> offsets(seq.size()+1, 0);
		map<int, string> abbrs;
		int refPos = 0, c = 0;
		for (int i = 0; i < edit.size(); ) {
			McAssert(c <= seq.size(), "%s %d %d", name.c_str(), c, seq.size());
			offsets[c] = c ? offsets[c - 1] : 0;
			McAssert(offsets[c] + c == refPos, "%s: %d+%d not equal to %d", name.c_str(), offsets[c], c, refPos);
			if (string("-JKLM").find(edit[i]) != string::npos) { // DEL
				int ds = 0;
				while (i < edit.size() && string("-JKLM").find(edit[i]) != string::npos) {
					//string s = S("SNP%c%c", reference[refPos].ref, reference[refPos].ref);
					//reference[refPos].abbrs[s]
					offsets[c]++, ds++, refPos++, i++;
				}
				// ignore JKLM
				if (edit[i - ds] == '-' && c && i < edit.size()) { // don't count trailing indels!!! 
					abbrs[c] = S("DEL%d", ds);
					reference[refPos].abbrs[abbrs[c]].alleles.set(id);
					reference[refPos].has2D7abbr |= (name.substr(0,6) == "CYP2D7");
					reference[refPos].abbrs[abbrs[c]].printStatus = true;
					//L("%s: %d_%d:%s\t", name.c_str(), c, revComp(refPos-ds-1), S("DEL%d", ds).c_str());
				}
				//L("%d_%d:%s\t", c, revComp(refPos-ds-1), S("DEL%d", ds).c_str());
			}
			if (i >= edit.size())
				break;
			if (islower(edit[i]) /*&& edit[i] != 'x'*/) { // INS
				string ii;
				while (i < edit.size() && islower(edit[i])) { 
					McAssert(c < offsets.size(), "Insertion edit distance consistency");
					offsets[c] = offsets[c - 1] - 1, ii += seq[c], seq[c] = tolower(seq[c]), c++, i++;
				}
				abbrs[c - ii.size()] = S("INS%s", ii.c_str());
				reference[c - ii.size() + offsets[c - ii.size()]].abbrs[abbrs[c - ii.size()]].alleles.set(id);
				reference[c - ii.size() + offsets[c - ii.size()]].has2D7abbr |= (name.substr(0,6) == "CYP2D7");
				reference[c - ii.size() + offsets[c - ii.size()]].abbrs[abbrs[c - ii.size()]].printStatus = true;

				//L("%d_%d:%s\t", c - ii.size(), revComp(refPos), S("INS%s", ii.c_str()).c_str());
			}
			else {
				if (edit[i] != '=') {
					string s = S("SNP%c%c", reference[refPos].ref, edit[i]);
					abbrs[c] = s; //seq[c]);
					reference[refPos].abbrs[s].alleles.set(id);
					reference[refPos].has2D7abbr |= (name.substr(0,6) == "CYP2D7");
					reference[refPos].abbrs[s].printStatus = true;
					//L("%d_%d:%s\t", c, revComp(refPos), S("SNP%c%c", reference[refPos].ref, seq[c]).c_str());
				}
				else {
					string s = S("SNP%c%c", reference[refPos].ref, reference[refPos].ref);
					//abbrs[c] = s; //seq[c]);

					if (name.substr(0,6) != "CYP2D7" || reference[refPos].has2D7abbr) 
						reference[refPos].abbrs[s].alleles.set(id);
					reference[refPos].abbrs[s].printStatus = false;
				}
				c++, refPos++, i++;
			}
		}
		testParsing(name, alleles["CYP2D6*1"].seq, seq, offsets, abbrs);

		alleles[name] = {id, seq, edit, c7, offsets, abbrs};
		if (name.substr(0,6) != "CYP2D7") id++;
		//L("\n");
	}

	for (auto &pos: reference) for (auto &abb: pos.abbrs) {
		if (abb.first.substr(0, 3) == "SNP" && string("ACGTN").find(abb.first[4]) == string::npos) {
			string a1, a2;
			switch (abb.first[4]) {
				case 'B': a1 = S("SNP%c%c", pos.ref, 'A'), a2 = S("SNP%c%c", pos.ref, 'C'); break; 
				case 'D': a1 = S("SNP%c%c", pos.ref, 'A'), a2 = S("SNP%c%c", pos.ref, 'G'); break; 
				case 'E': a1 = S("SNP%c%c", pos.ref, 'A'), a2 = S("SNP%c%c", pos.ref, 'T'); break; 
				case 'F': a1 = S("SNP%c%c", pos.ref, 'C'), a2 = S("SNP%c%c", pos.ref, 'G'); break; 
				case 'H': a1 = S("SNP%c%c", pos.ref, 'C'), a2 = S("SNP%c%c", pos.ref, 'T'); break;
				case 'I': a1 = S("SNP%c%c", pos.ref, 'G'), a2 = S("SNP%c%c", pos.ref, 'T'); break;
				default: McAssert(0, "Key %c weird, %d:%s", abb.first[4], &pos-&(reference[0]), abb.first.c_str());
			}
			pos.abbrs[a1].alleles.add(abb.second.alleles); // make 2D7-snps properly detected
			pos.abbrs[a2].alleles.add(abb.second.alleles); //
		}
	}
}
