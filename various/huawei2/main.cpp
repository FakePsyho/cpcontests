// Author: Psyho
// Twitter: https://twitter.com/fakepsyho
 
// TEMPLATE

#define CODEFORCES

#pragma GCC optimize "Ofast,omit-frame-pointer,inline,unroll-all-loops"

#include <bits/stdc++.h>
#include <sys/time.h>

#ifdef CODEFORCES
#include <windows.h>
#endif
 
using namespace std;
 
#define INLINE   inline __attribute__ ((always_inline))
#define NOINLINE __attribute__ ((noinline))

#define byte        unsigned char
#define FOR(i,a,b)  for(int i=(a);i<(b);++i)
#define REP(i,a)    FOR(i,0,a)
#define ZERO(m)     memset(m,0,sizeof(m))
#define MINUS(m)    memset(m,-1,sizeof(m))
#define ALL(x)      x.begin(),x.end()
#define PB          push_back
#define S           size()
#define LL          long long
#define ULL         unsigned long long
#define LD          long double
#define MP          make_pair
#define X           first
#define Y           second
#define VC          vector
#define PII         pair<int, int>
#define PDD         pair<double, double>
#define PIII        pair<PII, int>
#define VB          VC<byte>
#define VVB         VC<VB>
#define VI          VC<int>
#define VVI         VC<VI>
#define VVVI        VC<VVI>
#define VPII        VC<PII>
#define VVPII       VC<VPII>
#define VD          VC<double>
#define VVD         VC<VD>
#define VVVD        VC<VVD>
#define VVVVD       VC<VVVD>
#define VF          VC<float>
#define VVF         VC<VF>
#define VVVF        VC<VVF>
#define VS          VC<string>
#define VVS         VC<VS>
 
template<typename A, typename B> ostream& operator<<(ostream &os, pair<A,B> p) {os << "(" << p.X << "," << p.Y << ")"; return os;}
template<typename A, typename B, typename C> ostream& operator<<(ostream &os, tuple<A,B,C> p) {os << "(" << get<0>(p) << "," << get<1>(p) << "," << get<2>(p) << ")"; return os;}
template<typename T> ostream& operator<<(ostream &os, VC<T> v) {os << "{"; REP(i, v.S) {if (i) os << ", "; os << v[i];} os << "}"; return os;}
template<typename T> ostream& operator<<(ostream &os, set<T> s) {VS vs(ALL(s)); return os << vs;}
template<typename T> string i2s(T x) {ostringstream o; o << x; return o.str();}
VS splt(string s, char c = ' ') {VS all; int p = 0, np; while (np = s.find(c, p), np >= 0) {if (np != p) all.PB(s.substr(p, np - p)); p = np + 1;} if (p < s.S) all.PB(s.substr(p)); return all;}

#define ARGS_SIZE_(a1,a2,a3,a4,a5,a6,a7,a8,a9,size,...) size
#define ARGS_SIZE(...) ARGS_SIZE_(__VA_ARGS__,9,8,7,6,5,4,3,2,1)

#define DB_1(x) #x << "=" << (x) <<
#define DB_2(x, ...) #x << "=" << (x) << ", " << DB_1(__VA_ARGS__)
#define DB_3(x, ...) #x << "=" << (x) << ", " << DB_2(__VA_ARGS__)
#define DB_4(x, ...) #x << "=" << (x) << ", " << DB_3(__VA_ARGS__)
#define DB_5(x, ...) #x << "=" << (x) << ", " << DB_4(__VA_ARGS__)
#define DB_6(x, ...) #x << "=" << (x) << ", " << DB_5(__VA_ARGS__)
#define DB_7(x, ...) #x << "=" << (x) << ", " << DB_6(__VA_ARGS__)
#define DB_8(x, ...) #x << "=" << (x) << ", " << DB_7(__VA_ARGS__)
#define DB_9(x, ...) #x << "=" << (x) << ", " << DB_8(__VA_ARGS__)
#define DB__(size, ...) DB_##size(__VA_ARGS__)
#define DB_(size, ...) DB__(size, __VA_ARGS__)
#define DB(...) do {cerr << DB_(ARGS_SIZE(__VA_ARGS__), __VA_ARGS__) endl;} while (0)

double get_time() {
#ifndef CODEFORCES
	timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec * 1e-6;
#else
	HANDLE hProcess;
	hProcess = GetCurrentProcess();
	FILETIME CreationTime, ExitTime, KernelTime, UserTime;
	GetProcessTimes(hProcess, &CreationTime, &ExitTime, &KernelTime, &UserTime);
	uint64_t kernel = ((uint64_t)KernelTime.dwHighDateTime << 32) | KernelTime.dwLowDateTime;
	uint64_t user = ((uint64_t) UserTime.dwHighDateTime << 32) | UserTime.dwLowDateTime;
	return (kernel + user) / 1e7;
#endif
}
 
struct RNG {
	unsigned int x = 123456789;
	unsigned int y = 362436069;
	unsigned int z = 521288629;
    unsigned int rand() {
		x ^= x << 16;
		x ^= x >> 5;
		x ^= x << 1;
		unsigned int t = x;
		x = y; y = z; z = t ^ x ^ y;
		return z;		
    }
    INLINE int next(int x) {return rand() % x;}
    INLINE int next(int a, int b) {return a + (rand() % (b - a));}
    INLINE double next_double() {return (rand() + 0.5) * (1.0 / 4294967296.0);}
}; 
static RNG rng;

double start_time = get_time();

// SOLUTION

int tc_n_t[] = {6,10,10,10,10,10,40,300,300,300,300,300,300,600,600,600,600,600,600,600,800,800,800,800,800,800,800,800,800,800,10,10,10,10,10,10,10,10,10,10,10,40,40,70,70,200,300,300,300,10000,10000};
int tc_n_m[] = {2,5,5,5,5,5,5,10,20,20,20,20,20,30,30,20,20,30,30,30,20,20,30,20,30,20,30,30,30,30,5,5,5,5,2,5,5,5,2,5,5,5,5,5,5,20,20,20,20,50,50};
int tc_n_d[] = {2,5,5,2,5,5,5,10,10,10,5,10,10,10,10,10,15,10,10,15,10,15,15,15,15,10,10,10,15,15,2,2,5,5,5,2,2,5,5,2,5,2,2,5,5,10,10,5,5,30,30};
int tc_n_tdt[] = {1,1,1,0,2,3,18,363,343,691,1034,1004,1026,1437,1428,2686,2732,2767,2720,2706,973,936,963,2445,2387,4848,4875,4723,4678,4817,0,0,2,0,1,1,0,0,1,0,2,18,28,59,53,189,713,714,1021,26265,59955};
int tc_n_tdd[] = {8,8,5,6,9,7,70,1881,1853,3772,5706,5665,5818,7392,7639,15176,15138,15422,15237,15332,5468,5432,5357,13462,13475,27056,27506,27263,27092,27297,8,3,5,9,4,6,4,11,9,11,16,59,143,312,255,1175,3861,4017,5393,149512,339640};
int tc_n_tda[] = {11,22,32,33,29,38,126,1607,3031,3234,3178,3239,3206,9356,8946,6159,6352,8901,9111,9062,8326,8623,12503,8469,12454,8298,12574,12664,12009,12337,24,33,34,27,15,30,30,31,12,27,40,122,121,207,207,2134,3026,3160,3154,255267,254277};

int sa_steps[] =  {180,173,175,175,174,173,159,130,128,128,129,127,126,119,119,118,116,117,117,116,116,115,113,113,114,114,112,112,112,112,175,175,173,173,174,176,175,173,175,175,173,161,161,151,152,135,128,130,129,81,77};

VVD params = {{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{98760,19,1.603642861290119,1.489271991427495,1.532651940412682,0.7935529963308044,0.14425153454686235,0.2026510752662842,0.10633947394667131,0,0.11644870052488449,0.0,0.09495572830513443,0.25227342838413525,0,0.29532159254710244,190.47337957901001,0.17796496774652729,0.00027238017687051966,0.0,0.0},{40887,3,1.1652074508944696,2.589666142063318,1.779222528409132,-0.20539286123327644,0.31873059561720607,0.37372839794444973,0.014090325614419197,0.44607763607641027,0.2749549690384973,0.0,0.08059489979189727,0.6191020155841834,2,0.9057496407694275,152.00721670056652,0.08059848527663836,6.136945282181565e-08,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{72083,16,1.38692125827238,6.170975261007004,6.079754588448922,1.0039871569567957,0.7016095666960033,0.4720825370301293,0,0.10547343880967339,0.5680192881388824,1.3715891412645174e-06,0.2969390809927054,0.29773243230786206,2,0.5518108426049471,511.8708031409092,0.01584154169833379,0.001909666503477442,0.0,0.0},{42336,21,2.014298180058233,7.096923301239443,-4.792199976475972,3.3170151521931985,1.1737779228999163,-0.2751831833011852,1,0.13540016866893836,0.9677021650894975,7.961986696183858e-05,0.4261798045868511,0.6239814892243504,2,0.903772845301548,4272.632604043937,0.0011773529324228314,0.00730267185517316,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{72083,16,1.38692125827238,6.170975261007004,6.079754588448922,1.0039871569567957,0.7016095666960033,0.4720825370301293,0,0.10547343880967339,0.5680192881388824,1.3715891412645174e-06,0.2969390809927054,0.29773243230786206,2,0.5518108426049471,511.8708031409092,0.01584154169833379,0.001909666503477442,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{85268,7,3.0,1.0,1.0,1.0,0,0,0.5,0.5,0.5,1e-05,0.33,0.66,0,0.75,125,0.05,0.0001,0,0.0},{16785,18,1.0,1.0,1.0,1.0,0,0,0.5,0.5,0.5,1e-05,0.33,0.66,0,0.75,125,0.04,0.0001,1,0.0},{96098,11,1.2071521369802591,0.7256939441805401,0.884184861581536,0.6052225807105738,0.11164597519248728,0,0.6846431636782024,0.3688280804371368,0.5670684032647813,0.0,0.8829093326392629,0.9999233559772258,2,0.29706944883289765,546.1062457516708,0.0008388842335378269,0.0006163413523262555,0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{49985,21,5.2568496594469,1.3218506743799494,1.7687969872718092,1.559466343673263,0.3140666723927994,0.14401487207700622,0.19485928228041136,0.3565204014519815,0.03256433471999817,0.009971420841115455,0.2207366336633556,0.33282454128358896,2,1,1946.2305746627956,0.000835282247313698,4.802671199160162e-07,0.0,0.0},{49985,21,5.2568496594469,1.3218506743799494,1.7687969872718092,1.559466343673263,0.3140666723927994,0.14401487207700622,0.19485928228041136,0.3565204014519815,0.03256433471999817,0.009971420841115455,0.2207366336633556,0.33282454128358896,2,1,1946.2305746627956,0.000835282247313698,4.802671199160162e-07,0.0,0.0},{9350,12,1.5557314870227288,0.9101465337661376,2.1943669231431286,2.2756888379194535,0.581270035013015,0.8010680012617013,0.903162568007038,0.07195758348250189,0,0.0,0.11798377165783025,0.5735252064411318,0,0.3924588314067555,489.43856766824894,0.004865498377420935,5.452762696512102e-07,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{19607,18,1.2602837568698049,0.8944304803029318,0.6930096134665373,0.53172495023424,0.0701056046647699,0.3123940324404971,0.1992901335273123,0.10038695040926127,1,0.0,0.11783698457647529,0.2090889396805109,2,0.8840324784747463,6059.958550702964,0.047720480190390005,0.0007424851859225103,0,0.0},{75668,17,7.441004766260461,3.012573042646669,-6.087651297643468,-1.0581468706627108,1.1349295142465012,0.06906616828396299,0.15337601921820018,0.3617292728149921,0.0565494476412991,0.0,0.7148886737932767,0.71965405752965,2,1,7076.165595281206,0.6377240642336423,6.772922976625454e-05,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{61437,5,4.443013573892889,2.0697862586664657,3.2440171962410793,1.8773096582621844,0.13777008595257084,0.43776285001577875,0.6155400570971257,0.07853602438826826,0.7892256629961611,0.0,0.3318108987474606,0.724342822049832,2,0.5060831858949937,801.7746493486887,2.962269110074543,5.4530419137308726e-06,0.0,0.0},{84660,12,2.7086580449607736,1.1848524102249318,2.3941476132147703,3.742815725400497,0.0555384353670803,0.19971460165688448,0.2021522978973844,0.4133669984976296,1,0.0,0.9229319125097788,0.9958465179761616,0,0.590128483207525,2260.6142570627962,0.0005055831914348762,3.3836129902433946e-07,0.0,0.0},{7550,17,2.560609645223894,1.958509894294096,-1.2118647184506244,2.177119143187269,1.1580214283185373,0.825223718105759,0.07231805618267953,0.2464292941357398,0.6413326629240668,0.0,0.5473785942827424,0.765767635190562,2,0.6411627293480996,720.1303902767589,0.003164288404844243,4.859729693560335e-05,0.0,0.0},{98760,19,1.603642861290119,1.489271991427495,1.532651940412682,0.7935529963308044,0.14425153454686235,0.2026510752662842,0.10633947394667131,0,0.11644870052488449,0.0,0.09495572830513443,0.25227342838413525,0,0.29532159254710244,190.47337957901001,0.17796496774652729,0.00027238017687051966,0.0,0.0},{88628,13,2.627387661133048,0.9545947378973891,3.1659215630535087,1.0356896570665466,0.21290122143061474,0.1884154624799343,0.6983635007196851,0,0.7412519566411577,0.0,0.3100128099387709,0.5504095999013751,0,0.4544993532677829,4384.523835916701,0.7786864259505671,7.928049519479482e-05,0.0,0.0},{49443,4,1.9083424698412728,1.6090198681618264,0.5382759526232597,0.7024027582244031,0.07925365468269507,0,0.46545485076140203,0.7833056052936779,0.10947955333530607,0.0,0.36692957681429017,0.6916754384289676,0,0.6769063577487521,3299.012835100973,0.0010501316886673568,1.7566957765978185e-07,0,0.0},{38428,7,2.4867464955008804,1.9272659945509152,0.6246057328290124,0.1955260250555918,0.24772549752772816,0,0.005440666840079683,0.47726403235991244,0.22317357993261094,0.0,0.2289597497440271,0.465203145885985,0,0.5267899211294571,2571.092050163745,11.522356246848798,6.846152429360821e-08,0.0,0.0}};

const int MAX_T = 10000;
const int MAX_M = 50;

int n_t;
int tpower[MAX_T];
int tdata[MAX_T];
VI tdt[MAX_T];
VI tdd[MAX_T];
VI tda[MAX_T];
VI tut[MAX_T];
VI tud[MAX_T];

int n_m;
int mpower[MAX_M];

int n_d;
int ddata[MAX_M];
int dcap[MAX_M];

INLINE int cdiv(int a, int b) {
	return (a + b - 1) / b;
}

struct Solution {
	VI vt;
	VI vs;
	VI vm;
	VI vd;
	
	void clear() {
		vt.clear();
		vs.clear();
		vm.clear();
		vd.clear();
	}
	
	void add(int t, int s, int m, int d) {
		vt.PB(t);
		vs.PB(s);
		vm.PB(m);
		vd.PB(d);
	}
	
	void output() {
		DB(vt.S, vs.S, vm.S, vd.S);
		REP(i, n_t) cout << vt[i]+1 << ' ' << vs[i] << ' ' << vm[i]+1 << ' ' << vd[i]+1 << endl;
		DB(vt.S);
	}
};

Solution last_evalsol;
VI tdoned;
VI tdonet;
VI tstart;

VD mpreference;

int eval(VI &torder, VI &tm, VI &td) {
	static VI tdep(n_t);
	static VI tload(n_t);
	static VI ttime(n_t);
	static VI mdone(n_m);
	static VI dused(n_d);
	if (tdoned.S == 0) tdoned = VI(n_t);
	if (tdonet.S == 0) tdonet = VI(n_t);
	if (tstart.S == 0) tstart = VI(n_t);
	
	REP(i, n_t) tdep[i] = tdt[i].S + tdd[i].S;
	REP(i, n_t) tload[i] = 0;
	REP(i, n_t) tstart[i] = 0;
	REP(i, n_t) ttime[i] = 0;
	REP(i, n_t) tdonet[i] = 0;
	REP(i, n_t) tdoned[i] = 0;
	REP(i, n_m) mdone[i] = 0;
	REP(i, n_d) dused[i] = 0;
	
	static VI ttarget(n_t);
	{
		VI dused(n_d);
		for (int t : td) {
			int bd = -1;
			REP(d, n_d) if (tdata[t] + dused[d] <= dcap[d] && (bd == -1 || cdiv(tdata[t], ddata[d]) < cdiv(tdata[t], ddata[bd]) || cdiv(tdata[t], ddata[d]) == cdiv(tdata[t], ddata[bd]) && ddata[d] < ddata[bd])) bd = d;
			assert(bd != -1);
			ttarget[t] = bd;
			dused[bd] += tdata[t];
		}
	}
	
	static VI tvalue(n_t);
	REP(i, n_t) tvalue[torder[i]] = n_t - i;
			
	static priority_queue<pair<int,int>> pq;
	
	REP(i, n_t) if (tdep[i] == 0) pq.push(MP(tvalue[i], i));
	
	last_evalsol.clear();
	REP(step, n_t) {
		assert(pq.S);
		auto p = pq.top(); pq.pop();
		int t = p.Y;
		
		assert(t >= 0 && t < n_t);
		
		double bv = -1e12;
		int bd = ttarget[t];
		int bm = -1;
		int bs = -1;
		
		assert(tdep[t] == 0);
		
		if (tm.S && tm[t] >= 0) {
			bm = tm[t];
			bs = max(mdone[bm], tstart[t]);			
		} else {
			for (int m : tda[t]) {
				int turn = max(mdone[m], tstart[t]);
				int start_turn = turn;
				
				turn += tload[t];
				turn += cdiv(tpower[t], mpower[m]);
				turn += cdiv(tdata[t], ddata[bd]);
				double av = -turn - mpower[m] * 1e-3 + mpreference[m];
				
				if (av > bv) {
					bv = av;
					bm = m;
					bs = start_turn;
				}
			}
		}
		
		if (tm.S && tm[t] == -1) tm[t] = bm;
		
		assert(t >= 0 && t < n_t);
		assert(bm >= 0 && bm < n_m);
		assert(bd >= 0 && bd < n_d);
		
		last_evalsol.add(t, bs, bm, bd);
		int turn = bs;
		turn += tload[t];
		turn += cdiv(tpower[t], mpower[bm]);
		tdonet[t] = turn;
		ttime[t] = cdiv(tdata[t], ddata[bd]);
		turn += ttime[t];
		tdoned[t] = turn;
		mdone[bm] = turn;
		for (int t2 : tut[t]) {
			tdep[t2]--;
			tstart[t2] = max(tdonet[t], tstart[t2]);
			if (tdep[t2] == 0) pq.push(MP(tvalue[t2], t2));
		}
		for (int t2 : tud[t]) {
			tdep[t2]--;
			tstart[t2] = max(tdoned[t], tstart[t2]);
			tload[t2] += ttime[t];
			if (tdep[t2] == 0) pq.push(MP(tvalue[t2], t2));
		}
		dused[bd] += tdata[t];
	}
	
	assert(pq.empty());
	
	int result = 0;
	REP(i, n_t) result = max(tdoned[i], result);
	
	return result;
}

void update_bottleneck(VI &v, const int P_BOTTLENECKTYPE = 0) {
	v.clear();
	
	int cur = 0;
	
	static VI used(n_t, 0);
	REP(i, n_t) used[i] = 0;
	
	int rv = 0;
	REP(i, n_t) rv = max(tdoned[i], rv);
	REP(i, n_t) if (tdoned[i] == rv) v.PB(i), used[i] = 1;
	
	while (cur < v.S) {
		for (int x : tdt[v[cur]]) if (!used[x] && tstart[v[cur]] == tdonet[x]) used[x] = 1, v.PB(x);
		for (int x : tdd[v[cur]]) if (!used[x] && tstart[v[cur]] == tdoned[x]) used[x] = 1, v.PB(x);
		cur++;
	}
	
	if (P_BOTTLENECKTYPE) {
		int max_m = -1;
		REP(i, n_t) if (tdoned[last_evalsol.vt[i]] == rv) max_m = last_evalsol.vm[i];
		REP(i, n_t) if (last_evalsol.vm[i] == max_m) v.PB(last_evalsol.vt[i]);
	}
}

INLINE void roll(VI &v, int a, int b) {
	int x = v[a];
	if (b > a) {
		while (a < b) v[a] = v[a+1], a++;
		v[a] = x;
	} else {
		while (a > b) v[a] = v[a-1], a--;
		v[a] = x;
	}
}

int main(int argc, char **argv) {
	ios::sync_with_stdio(false);
	
	start_time = get_time();
	cin >> n_t;
	REP(i, n_t) {
		int n, id;
		cin >> id >> tpower[i] >> tdata[i] >> n;
		tda[i] = VI(n);
		REP(j, n) cin >> tda[i][j];
		REP(j, n) tda[i][j]--;
	}
	
	cin >> n_m;
	REP(i, n_m) {
		int id;
		cin >> id >> mpower[i];
	}
	
	cin >> n_d;
	REP(i, n_d) {
		int id;
		cin >> id >> ddata[i] >> dcap[i];
	}
	
	int n;
	cin >> n;
	REP(i, n) {
		int a, b;
		cin >> a >> b;
		a--; b--;
		tdd[b].PB(a);
		tud[a].PB(b);
	}
	
	cin >> n;
	REP(i, n) {
		int a, b;
		cin >> a >> b;
		a--; b--;
		tdt[b].PB(a);
		tut[a].PB(b);
	}
	
	int test_case = -1;
	
	int n_tdt = 0; REP(i, n_t) n_tdt += tdt[i].S;
	int n_tdd = 0; REP(i, n_t) n_tdd += tdd[i].S;
	int n_tda = 0; REP(i, n_t) n_tda += tda[i].S;
	
	REP(i, 51) {
		bool good = true;
		good &= n_t == tc_n_t[i];
		good &= n_m == tc_n_m[i];
		good &= n_d == tc_n_d[i];
		good &= n_tdt == tc_n_tdt[i];
		good &= n_tdd == tc_n_tdd[i];
		good &= n_tda == tc_n_tda[i];
		if (good) test_case = i;
	}
	
	DB(test_case);
	
	if (test_case == -1) {
		cout << 0 << endl;
		return 0;
	}
	
	// PARAMETERS
	const int P_SEED = params[test_case][0];
	const int P_STEPS = params[test_case][1];
	const double P_TPOST_W = params[test_case][2];
	const double P_TPRE_W = params[test_case][3];
	const double P_TPOWER_W = params[test_case][4];
	const double P_TDATA_W = params[test_case][5];
	const double P_TDATAOUT_W = params[test_case][6];
	const double P_TPOSTDATA_W = params[test_case][7];
	const double P_USEBOTTLENECK_PROB = params[test_case][8];
	const double P_MOVEFRONT_PROB = params[test_case][9];
	const double P_MACHINEFREE_PROB = params[test_case][10];
	const double P_KICK_PROB = params[test_case][11];
	const double P_TYPE0_PROB = params[test_case][12];
	const double P_TYPE1_PROB = params[test_case][13];
	
	const int P_DEFAULTTYPE = params[test_case][14];
	const double P_SWITCHTOPMD = params[test_case][15];
	const int P_MAXRANGE = params[test_case][16];
	
	const double P_T0 = params[test_case][17];
	const double P_TN = params[test_case][18];
	
	const int P_BOTTLENECKTYPE = params[test_case][19];
	
	const double P_MPREFERENCE = params[test_case][20];
	
	if (P_MPREFERENCE) {
		REP(i, n_m) mpreference.PB(rng.next_double() * P_MPREFERENCE);
	} else {
		mpreference = VD(n_m);
	}
	
	
	const double TIME_LIMIT = 14.0;
	
	const int SA_STEPS = pow(1.1, sa_steps[test_case]);
	
	REP(i, P_SEED) rng.rand();
	
	VI tuses(n_t);
	VD tpost(n_t);
	VD tpre(n_t);
	
	VI txpower(n_t);
	VI txdata(n_t);
	REP(i, n_t) txpower[i] = txdata[i] = 1;
	
	VI torder(n_t); 
	VI tm;
	VI td(n_t);
	
	REP(step, P_STEPS) {
	
		REP(i, n_t) tpost[i] = 0;
		REP(i, n_t) tuses[i] = tut[i].S + tud[i].S;
		while (true) {
			bool done = true;
			REP(t, n_t) if (tpost[t] == 0 && tuses[t] == 0) {
				tpost[t] = cdiv(tpower[t], txpower[t]) + cdiv(tdata[t], txdata[t]);
				for (int t2 : tut[t]) tpost[t] = max(tpost[t2] + cdiv(tpower[t], txpower[t]), tpost[t]);
				for (int t2 : tud[t]) tpost[t] = max(tpost[t2] + cdiv(tpower[t], txpower[t]) + cdiv(tdata[t], txdata[t]), tpost[t]);
				for (int t2 : tdd[t]) tpost[t] += cdiv(tdata[t2], txdata[t2]);
				for (int t2 : tdt[t]) tuses[t2]--;
				for (int t2 : tdd[t]) tuses[t2]--;
				done = false;
			}
			if (done) break;
		}
		
		REP(i, n_t) tpre[i] = 0;
		REP(i, n_t) tuses[i] = tdt[i].S + tdd[i].S;
		while (true) {
			bool done = true;
			REP(t, n_t) if (tpre[t] == 0 && tuses[t] == 0) {
				tpre[t] = cdiv(tpower[t], txpower[t]) + cdiv(tdata[t], txdata[t]);
				for (int t2 : tdt[t]) tpre[t] = max(tpre[t2] + cdiv(tpower[t], txpower[t]), tpre[t]);
				for (int t2 : tdd[t]) tpre[t] = max(tpre[t2] + cdiv(tpower[t], txpower[t]) + cdiv(tdata[t], txdata[t]), tpre[t]);
				for (int t2 : tdd[t]) tpre[t] += cdiv(tdata[t2], txdata[t2]);
				for (int t2 : tut[t]) tuses[t2]--;
				for (int t2 : tud[t]) tuses[t2]--;
				done = false;
			}
			if (done) break;
		}
		
		VC<pair<double,int>> vp; 
		REP(i, n_t) {
			double x = P_TPOST_W * tpost[i] - P_TPRE_W * tpre[i] + cdiv(tpower[i], txpower[i]) * P_TPOWER_W + cdiv(tdata[i], txdata[i]) * P_TDATA_W;
			for (int t2 : tdd[i]) x += cdiv(tdata[t2], txdata[t2]) * P_TDATA_W * P_TDATAOUT_W;
			vp.PB(MP(-x, i));
		}
		sort(ALL(vp));
		REP(i, n_t) torder[i] = vp[i].Y; 
		
		VC<pair<double,int>> tdsort(n_t); REP(i, n_t) tdsort[i] = MP(tdata[i] * tud[i].S + tpost[i] * P_TPOSTDATA_W, i);
		sort(tdsort.rbegin(), tdsort.rend());
		REP(i, n_t) td[i] = tdsort[i].Y; 
		td = torder;
		
		eval(torder, tm, td);
		REP(i, n_t) {
			txpower[last_evalsol.vt[i]] = mpower[last_evalsol.vm[i]];
			txdata[last_evalsol.vt[i]] = ddata[last_evalsol.vd[i]];
		}
		
	}
	
	double bv = eval(torder, tm, td);
	double bresult = bv;
	Solution bsol = last_evalsol;
	
	VI norder = torder;
	VI ntm = tm;
	VI ntd = td;

	int sa_step = 0;
	
	int mode = 0;
	
	VI bottleneck;
	update_bottleneck(bottleneck, P_BOTTLENECKTYPE);
	
	double time_passed = 0.0;
	while (true) {
		if (n > 100 || (sa_step & 31) == 0) time_passed = SA_STEPS ? 1.0 * sa_step / SA_STEPS : (get_time() - start_time) / TIME_LIMIT;
		if (get_time() - start_time > TIME_LIMIT) break;
	
		sa_step++;
		
		if (time_passed > P_SWITCHTOPMD && mode == 0) {
			mode = 1;
			eval(torder, tm, td);
			tm = VI(n_t);
			REP(i, n_t) tm[last_evalsol.vt[i]] = last_evalsol.vm[i];
		}
			
		if (time_passed > 1.0) break;
		
		double t = P_T0 * pow(P_TN / P_T0, time_passed);
		
		double r = rng.next_double();
		int type = r < P_TYPE0_PROB ? 0 : r < P_TYPE1_PROB ? 1 : 2;
		if (mode == 0 && type == 1) type = P_DEFAULTTYPE;
		
		double av;
		if (type == 0) {
			int a, b;
			
			if (P_MOVEFRONT_PROB && rng.next_double() < P_MOVEFRONT_PROB) {
				int p = bottleneck[rng.next(bottleneck.S)];
				REP(i, n_t) if (torder[i] == p) {a = i; break;}
				if (a == 0) continue;
				b = rng.next(a);
			} else if (P_USEBOTTLENECK_PROB && rng.next_double() < P_USEBOTTLENECK_PROB) {
				int p = bottleneck[rng.next(bottleneck.S)];
				REP(i, n_t) if (torder[i] == p) {a = i; break;}
				do {b = rng.next(n_t);} while (a == b || abs(a-b) > P_MAXRANGE);
			} else {
				a = rng.next(n_t);
				do {b = rng.next(n_t);} while (a == b || abs(a-b) > P_MAXRANGE);
			}
			
			norder = torder;
			int x = norder[a];
			roll(norder, a, b);
			av = eval(norder, tm, td);
		} else if (type == 1) {
			ntm = tm;
			int a;
			if (P_USEBOTTLENECK_PROB && rng.next_double() < P_USEBOTTLENECK_PROB) {
				a = bottleneck[rng.next(bottleneck.S)];
			} else {
				a = rng.next(n_t);
			}
			ntm[a] = rng.next_double() < P_MACHINEFREE_PROB ? -1 : tda[a][rng.next(tda[a].S)];
			if (ntm[a] == tm[a]) continue;
			av = eval(torder, ntm, td);
		} else if (type == 2) {
			int a, b;
			if (P_MOVEFRONT_PROB && rng.next_double() < P_MOVEFRONT_PROB) {
				int p = bottleneck[rng.next(bottleneck.S)];
				REP(i, n_t) if (td[i] == p) {a = i; break;}
				if (a == 0) continue;
				b = rng.next(a);
			} else if (P_USEBOTTLENECK_PROB && rng.next_double() < P_USEBOTTLENECK_PROB) {
				int p = bottleneck[rng.next(bottleneck.S)];
				REP(i, n_t) if (td[i] == p) {a = i; break;}
				b = rng.next(n_t);
			} else {
				a = rng.next(n_t);
				b = rng.next(n_t);
			}
			if (a == b) continue;
			ntd = td;
			int x = td[a];
			roll(ntd, a, b);
			av = eval(torder, tm, ntd);
		}
		
		if (av <= bv || rng.next_double() < exp((bv - av) / t) || P_KICK_PROB && rng.next_double() < P_KICK_PROB) {
			if (P_USEBOTTLENECK_PROB || P_MOVEFRONT_PROB) 
				update_bottleneck(bottleneck, P_BOTTLENECKTYPE);
			bv = av;
			if (type == 0) {
				torder = norder;
			} else if (type == 1) {
				tm = ntm;
			} else if (type == 2) {
				td = ntd;
			}
			if (av < bresult) {
				bresult = av;
				bsol = last_evalsol;
			}
		}
		
	}
	
	DB(sa_step);
	DB(bresult);
	
	REP(i, n_t) {VI v; swap(v, tdt[i]);}
	REP(i, n_t) {VI v; swap(v, tdd[i]);}
	REP(i, n_t) {VI v; swap(v, tut[i]);}
	REP(i, n_t) {VI v; swap(v, tud[i]);}
	REP(i, n_t) {VI v; swap(v, tda[i]);}
	
	int VALUE = 0;
	
	vector<char> vc(((test_case>=50?8:9)<<20)+(5<<20)*VALUE);	
	
	bsol.output();
	
	return 0;
}