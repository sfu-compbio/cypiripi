#ifndef Common_H
#define Common_H

#include <ilcplex/ilocplex.h>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <cstdio>
#include <clocale>
#include <fstream>
#include <cctype>
#include <iostream>
#include <algorithm>
#include <cstdarg>
#include <inttypes.h>
#include <unordered_set>
#include <unordered_map>
#include <sys/time.h>

#define max(x,y)\
	(((x)>(y))?(x):(y))
#define L(c,...)\
	fprintf(stdout, c, ##__VA_ARGS__)
#define E(c,...)\
	fprintf(stderr, c, ##__VA_ARGS__)
#define McAssert(expr, fmt, ...)\
	{if(!(expr)) { E("%s @ %s:%d: " fmt "\n",#expr,__FILE__,__LINE__,##__VA_ARGS__);abort(); }}

#define foreach(i,c) \
	for (auto i = (c).begin(); i != (c).end(); ++i)

inline char revComp (char c) 
{
	std::string x = "ACGTacgtBDEFHIJKLM";
	std::string y = "TGCAtgcaIHEFDBMLKJ";
	auto p = x.find(c);
	return p == std::string::npos ? c : y[p];
}

inline std::string revComp (const std::string &s) 
{
	std::string r(s.size(), ' ');
	for (int i = s.size() - 1; i >= 0; i--)
		r[s.size() - 1 - i] = revComp(s[i]);
	return r;
}

inline std::string S (const char* fmt, ...) 
{
	char *ptr = 0;
    va_list args;
    va_start(args, fmt);
    vasprintf(&ptr, fmt, args);
    va_end(args);
    std::string s = ptr;
    free(ptr);
    return s;
}

inline char actgnum (char c) 
{
	const static char table[] = { 0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0 };
	return table[c - 'A'];
}

inline int revComp (int c) 
{
	int x = (9430 - c) - 1619 + 2;
	return x - (bool)(x < 0);
}

inline std::string realpath (const std::string &s) 
{
	char buffer[2048];
	return std::string(::realpath(s.c_str(), buffer));
}

inline uint64_t _zaman() 
{
	struct timeval t;
	gettimeofday(&t, 0);
	return (t.tv_sec * 1000000ll + t.tv_usec);
}

inline uint64_t zaman_last() 
{
	static uint64_t z = 0;
	uint64_t _z = _zaman() - z;
	z = _zaman();
	return _z / 1000000;
}

template<class T, typename K>
inline bool in(const T &t, const std::unordered_map<T,K> &c) 
{
	return c.find(t) != c.end();
}

template<class T>
inline bool in(const T &t, const std::unordered_set<T> &c)
{
	return c.find(t) != c.end();
}


#endif // Common_H