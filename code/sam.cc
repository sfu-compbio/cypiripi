/// 786

#include "sam.h"
#include <inttypes.h>
#include <fstream>

using namespace std;

/*template<typename T>
set<T> set_intersect(const set<T> &A, const set<T> &B) {
	set<T> C;
	set_intersection(A.begin(), A.end(), B.begin(), B.end(), std::inserter(C, C.begin()));
	return C;
}*/
unordered_map<string, int> readnames;
unordered_map<int, unordered_set<int>> pairedReads;

//extern double COV = 0;

//multimap<int, string> allAbbrs;
//unordered_set<int> cyp2D7reads;
extern string parFile;

void parseRead (string allele, int start, int rname, const string &cigar, const string &read, ref_t &reference, alleles_t &alleles) 
{
	//if (atoi(rname.c_str() + 5) > 40000)
	//	return;

	start--;
	if (start < 0) // MRFAST SCREW-UP
		return;
	if (start + read.size() > alleles[allele].seq.size() || cigar == "*") 
		return;

	int starto = start;

	int readHasIns = 0; // hack to prevent wrong 
	int readHasDel = 0;
	
	int num = 0, ri = 0;
	string newRead = "";
	for (int ci = 0; ci < cigar.size(); ci++) 
		if (isdigit(cigar[ci]))
			num = 10 * num + (cigar[ci] - '0');
		else {
			if (cigar[ci] != 'D')  //?
				num = min(num, (int)read.size() - ri);
			switch (cigar[ci]) {
				case 'D': 
					newRead += string(num, '-'), start += num;
					readHasDel += num;
					break;
				case 'M': case '=': case 'X': 
					newRead += read.substr(ri, num), ri += num, start += num;
					break;
				case 'I':  
					ri += num;
					readHasIns += num;
					break;
				default:
					McAssert(0, "CIGAR feature not implemented");
			} 
			num = 0;
		}

	start = starto;
	auto &al = alleles[allele];
	auto lb = al.abbrs.lower_bound(start);
	auto ub = al.abbrs.upper_bound(start + newRead.size());

	bitset supportAlleles(1);
	vector<pair<int, string>> added;
	vector<char> tentative;
	for (auto it = lb; it != ub; it++) {
		int pos = it->first;
		string abbr = it->second;
		McAssert(pos + al.offsets[pos] < reference.size(), "<offset> table is screwed up as usual");

		if (abbr.substr(0, 3) == "DEL") {
			int l = atoi(abbr.c_str() + 3);
			int remaining = min(l, (int)newRead.size() - (pos - start));
			if (readHasIns < remaining) {
				supportAlleles.intersect(reference[pos + al.offsets[pos]].abbrs[abbr].alleles);
				added.push_back({pos + al.offsets[pos], abbr});

				bool pc = supportAlleles.get(alleles[allele].id); 	
				if (allele.substr(0,6) == "CYP2D7" && !reference[pos + al.offsets[pos]].has2D7abbr && pc)
					supportAlleles.set(alleles[allele].id);

				tentative.push_back(bool((int)newRead.size() - (pos - start) < 4) || (pos - start < 4));

				//if (allele == "CYP2D6*19")
				//	E("%s [%s] %s %d %d_%d %d\n", read.c_str(), newRead.c_str(), cigar.c_str(), pos, pos-start, newRead.size()-(pos-start), tentative[tentative.size()-1]); 
			}
		}
		else if (abbr.substr(0, 3) == "INS") {
			string ins = abbr.substr(3);
			int remaining = min(ins.size(), newRead.size()-(pos-start));
			if (ins.substr(0, remaining) == newRead.substr(pos - start, remaining)) {
				for (int k = pos-start; k<pos-start+remaining; k++)
					newRead[k] = tolower(newRead[k]);
				if (readHasDel < remaining) {
					bool pc = supportAlleles.get(alleles[allele].id); 
					supportAlleles.intersect(reference[pos + al.offsets[pos]].abbrs[abbr].alleles);
					if (allele.substr(0,6) == "CYP2D7" && !reference[pos + al.offsets[pos]].has2D7abbr && pc)
						supportAlleles.set(alleles[allele].id);


					added.push_back({pos + al.offsets[pos], abbr});
					
					tentative.push_back(bool((int)newRead.size() - (pos - start + (int)ins.size()) < 4) || (pos - start < 4));

				//	if (allele == "CYP2D6*30")
				//		E("%s [%s] %s %d %d %d\n", read.c_str(), newRead.c_str(), cigar.c_str(), pos, newRead.size()-(pos-start+ ins.size()), tentative[tentative.size()-1]); 
				}
			}
		}
	}

	string lost="";

	bool has1235=0;
	for (int i = 0; i < newRead.size(); i++) {
		int p = start + i + al.offsets[start + i];
		if (isupper(newRead[i])) {
			string abbr = S("SNP%c%c", reference[p].ref, newRead[i]);

			//has1235|=(revComp(p)==138);

			bool pc = supportAlleles.get(alleles[allele].id); 
			supportAlleles.intersect(reference[p].abbrs[abbr].alleles);
			if (allele.substr(0,6) == "CYP2D7" && !reference[p].has2D7abbr && pc)
				supportAlleles.set(alleles[allele].id);

			// if (allele.substr(0,6) == "CYP2D7" && has1235 && pc && !supportAlleles.get(alleles[allele].id)) {
			// 	auto &X = reference[p].abbrs[abbr];
			// 	auto Y =  reference[p].abbrs[abbr].alleles;
			// }

			added.push_back({p, abbr});
			tentative.push_back(false);
		}
	}

	bool is = supportAlleles.get(alleles[allele].id);

	// if (allele.substr(0,6) == "CYP2D7" && has1235) {
	// 	E("pos %s %d -- %d %s %s %s\n", allele.c_str(), revComp(start+ al.offsets[start ]), is, newRead.c_str(), cigar.c_str(), lost.c_str());
	// }

	if (is) for (int i = 0; i < added.size(); i++) {
		reference[added[i].first].abbrs[added[i].second].coverage[alleles[allele].id].insert(rname);
		if (tentative[i])
			reference[added[i].first].abbrs[added[i].second].tentative.insert(rname);
		else
			reference[added[i].first].abbrs[added[i].second].reads.insert(rname);
	}

	/*if (is && allele.substr(0,6) == "CYP2D7") {
		auto iq = bitset(0);
		iq.set(alleles[allele].id);
		if (iq == supportAlleles) 
			cyp2D7reads.insert(rname);
	}*/

}

void readPairedSAMFile (string path, alleles_t &alleles, FILE *fd);
void readDumpFile (FILE *fi, ref_t &reference, alleles_t &alleles) ;


void readSAMFile (string path, ref_t &reference, alleles_t &alleles, string excludePath = "")
{
/*	ifstream fin((path + ".cyp6").c_str());
	if (fin.is_open()) {
		E("Found serialization!\n");
		boost::archive::text_iarchive ia(fin);
    	ia >> reference;
		E("Done!\n");
    	return;
	}*/

	// set<pair<int, string>> aa;
	// for (auto &a: alleles)
	// 	for (auto &ab: a.second.abbrs) {
	// 		aa.insert({ab.first + a.second.offsets[ab.first], ab.second});
	// 	}
	// for (auto &ab: aa)
	// 	allAbbrs.insert(ab);
	// //	L("%d -> %s\n", revComp(ab.first), ab.second.c_str());
	// //exit(0);

    FILE *fd = 0;
    /*fd = fopen(string(path + ".dump").c_str(), "r");
    /*if (fd) {
    	readDumpFile(fd, reference, alleles);
    	fclose(fd);
    	return;
    }
    else fd = 0; //fopen(string(path + ".dump").c_str(), "w");
	*/
	const int MAX_BUFFER = 8 * 1024 * 1024;
	char *sam_name = new char[MAX_BUFFER];

	unordered_set<string> exclude;
	/* exclude format: @name <svasta nesto> */
	if (excludePath != "") {
		ifstream fin(excludePath.substr(1).c_str());
		string s;
		while (getline(fin, s)) {
			sscanf(s.c_str(), "%s", sam_name);
			exclude.insert(sam_name + (sam_name[0] == '@' ? 1 : 0));
		}
	}

	FILE *fi = fopen(path.c_str(), "r");
	char* buffer = new char[MAX_BUFFER];
	char *sam_rname = new char[MAX_BUFFER],
		 *sam_cigar = new char[MAX_BUFFER],
		 *sam_read = new char[MAX_BUFFER],
		 *sam_rnext = new char[MAX_BUFFER];
	int sam_flag, sam_pos;
	int32_t sam_pnext, sam_tlen;
	uint8_t sam_mapq;
	size_t cnt = 0;
	fseek(fi, 0, SEEK_END);
	size_t fsz = ftell(fi);
	fseek(fi, 0, SEEK_SET);

	int id = 0;
	int cntex = 0;
	while (fgets(buffer, MAX_BUFFER, fi)) {
		if (buffer[0] == '@') 
			continue;
		if (!isalpha(buffer[0]))
			continue;
		cnt++;

		char *pbuf = buffer, *pc = buffer;
		int icnt = 0;
		string rname, allele, cigar, read;
		for (; *pc != '\n'; pc++) 
			if (*pc == '\t') {
				*pc = 0;
				if (icnt == 0) 
					rname = pbuf;
				else if (icnt == 1) 
					sam_flag = atoi(pbuf);
				else if (icnt == 2)
					allele = pbuf;
				else if (icnt == 3) 
					sam_pos = atoi(pbuf);  
				else if (icnt == 5)
					cigar = pbuf;
				else if (icnt == 9) {
					read = pbuf;
					break;
				}
				pbuf = pc + 1;
				*pc = '\t';
				icnt++;
			}
		
		//L("--%s -- %s %d\n",buffer, rname.c_str(),sam_flag);
		if (excludePath != "" && excludePath[0] == '+' && exclude.find(rname) == exclude.end()) {
			cntex++;
			continue;
		}
		if (excludePath != "" && excludePath[0] == '-' && exclude.find(rname) != exclude.end()) {
			cntex++;
			continue;
		}

		if (sam_flag & 0x10) // !(sam_flag & 0x8) && !(sam_flag & 0x40))
			rname += "/2";
		else
			rname += "/1";
		//L("%s %d\n",rname.c_str(),sam_flag);
		auto ir = readnames.find(rname);
		if (ir == readnames.end()) {
			int p = readnames.size();
			ir = readnames.insert({rname, p}).first;
		}
		
		parseRead(allele, sam_pos, ir->second, cigar, read, reference, alleles);
		//if (cnt % (1<<16) == 0) 
//			E("\r\tSAM loaded %.2lf%% ... ", 100.0*(ftell(fi))/fsz);
	}
	fclose(fi);
	L("SAM reader: total %'lu reads, %'d excluded\n", cnt, cntex);

	delete[] buffer;
	delete[] sam_name;
	delete[] sam_rname;
	delete[] sam_cigar;
	delete[] sam_read;
	delete[] sam_rnext;

	//L("%d reads unique to CYP2D7\n", cyp2D7reads.size());
	// int id7 = alleles["CYP2D7"].id;
	// for (auto &p: reference)
	// 	for (auto &a: p.abbrs) {
	// 		for (auto &c: a.second.coverage) if (c.first != id7)
	// 			for (auto &x: cyp2D7reads)
	// 				c.second.erase(x);
	// 	}

	if (fd) {
		int s = 0;
		for (auto &r: reference) for (auto &a: r.abbrs)
			s++;
		fprintf(fd, "%d ", s);
		for (auto &r: reference) for (auto &a: r.abbrs)
			fprintf(fd, "%d %s ", int(&r - &(reference[0])), a.first.c_str());
		for (auto &r: reference)
			for (auto &a: r.abbrs)
				a.second.dump(fd);
	}
	readPairedSAMFile(parFile, alleles, fd);
	if (fd) fclose(fd);

	double COV = 0;
	int c = 0;
	for (auto &pos: reference) {
		int p = &pos - &(reference[0]);
		if (revComp(p) > 5000) continue;
		
		for (auto &a: pos.abbrs) {
			COV += a.second.reads.size();
		}
		c++;
	}
	COV /= c;


	
	//E("COVERAGE IS %.2lf\n", COV/4);
}

void readPairedSAMFile (string path, alleles_t &alleles, FILE *fd = 0)
{
	if (path == "") return;
	ifstream fin(path.c_str());
	L("%s loaded as PE file\n", path.c_str());

	string prn, l, rn;
	char rnb[1000], all[1000];

	//E("\r\tReading PE information ... ");
	while (getline(fin, l) && l[0] == '@') ;	
	sscanf(l.c_str(), "%s", rnb); prn = rnb;
	int ln = 1, w = 0;
	while (getline(fin, l)) {
		uint32_t f;
		sscanf(l.c_str(), "%s %u %s", rnb, &f, all); 

		if ((ln & 1) && string(rnb) == prn) {
			if (f & 0x10) // (!(f & 0x8) && !(f & 0x40))
				rn = string(rnb) + "/2", prn += "/1";
			else
				rn = string(rnb) + "/1", prn += "/2";
			
			auto x = readnames.find(rn);
			auto y = readnames.find(prn);
			//E("%s %s %d %d\n",rn.c_str(), prn.c_str(),x!=readnames.end(),y!=readnames.end());
			if (x != readnames.end() && y != readnames.end()) {
				pairedReads[alleles[all].id].insert(x->second);
				pairedReads[alleles[all].id].insert(y->second);
				w++;
			}
		}
		prn = rnb;
		ln++;
		//if(ln==20)exit(1);
	}
	L("\tLoaded %d pairs\n", w);
	fin.close();

	if (fd) {
		for (auto i: pairedReads) { 
			fprintf(fd, "%d %d ", i.first, i.second.size());
			for (auto j: i.second) fprintf(fd, "%d ", j);
		}
	}
}

void readDumpFile (FILE *fi, ref_t &reference, alleles_t &alleles) {
	// pairedReads
	//E("\r\tReading dump information ... ");

	char ax[500];
	int sz, i, j, k;

	fscanf(fi, "%d", &sz);
	while (sz--) { fscanf(fi, "%d %s", &i, ax); if (reference[i].abbrs.find(ax) == reference[i].abbrs.end()) reference[i].abbrs[ax] = abbr_t();  }
	for (auto &r: reference)
		for (auto &a: r.abbrs)
			a.second.load(fi);

	fscanf(fi, "%d", &sz); while (sz--) { 
		fscanf(fi, "%d %d ", &i, &j);
		while (j--) { fscanf(fi, "%d ", &k); pairedReads[i].insert(k); }
	}
}