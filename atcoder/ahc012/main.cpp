// Author: Psyho
// Twitter: https://twitter.com/fakepsyho

#pragma optimize "Ofast,unroll-all-loops,inline"
 
#include <bits/stdc++.h>
#include <sys/time.h>
 
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
 
double get_time() {timeval tv; gettimeofday(&tv, NULL); return tv.tv_sec + tv.tv_usec * 1e-6;}
 
struct RNG {
    unsigned int MT[624];
    int index;
    RNG(int seed = 1) {init(seed);}
    void init(int seed = 1) {MT[0] = seed; FOR(i, 1, 624) MT[i] = (1812433253UL * (MT[i-1] ^ (MT[i-1] >> 30)) + i); index = 0; }
    void generate() {
        const unsigned int MULT[] = {0, 2567483615UL};
        REP(i, 227) {unsigned int y = (MT[i] & 0x8000000UL) + (MT[i+1] & 0x7FFFFFFFUL); MT[i] = MT[i+397] ^ (y >> 1); MT[i] ^= MULT[y&1]; }
        FOR(i, 227, 623) {unsigned int y = (MT[i] & 0x8000000UL) + (MT[i+1] & 0x7FFFFFFFUL); MT[i] = MT[i-227] ^ (y >> 1); MT[i] ^= MULT[y&1]; }
        unsigned int y = (MT[623] & 0x8000000UL) + (MT[0] & 0x7FFFFFFFUL); MT[623] = MT[623-227] ^ (y >> 1); MT[623] ^= MULT[y&1];
    }
    unsigned int rand() { if (index == 0) generate(); unsigned int y = MT[index]; y ^= y >> 11; y ^= y << 7  & 2636928640UL; y ^= y << 15 & 4022730752UL; y ^= y >> 18; index = index == 623 ? 0 : index + 1; return y;}
    INLINE int next() {return rand(); }
    INLINE int next(int x) {return rand() % x; }
    INLINE int next(int a, int b) {return a + (rand() % (b - a)); }
    INLINE double next_double() {return (rand() + 0.5) * (1.0 / 4294967296.0); }
    INLINE double next_double(double a, double b) {return a + next_double() * (b - a); }
};
 
static RNG rng;

double start_time = get_time();

const int MAX_N = 10000;
const double MAX_TIME = 2.96;
const int FAR = 1e8;
const int MAX_COORD=10020;
int XCUTS=9;
int YCUTS=100-XCUTS;
const double TYPE1_PROB=0.65;
const int TYPE1_RNG=4;
const int ACUTS=100;

int n, k;
int a[11];
PII r[MAX_N];

double av, xv, bv;

VI cutsx, cutsy;
VI cutsxpos, cutsypos;

int cur[ACUTS];
int tcur[ACUTS];
int best[ACUTS];

int cnt[MAX_N];
int ygroup[MAX_COORD*2];
int groups[ACUTS];

int last_og_score;
double moment;
int plo[ACUTS],phi[ACUTS];

double eval(bool ygroup_same=false) {
	REP(i, ACUTS) tcur[i] = cur[i];
	sort(tcur,tcur+XCUTS);
	sort(tcur+XCUTS,tcur+ACUTS);
	FOR(i, 1, 11) cnt[i] = 0;
	
	if (!ygroup_same) {
		REP(i, YCUTS+1) {
			int lo = i == 0 ? 0 : MAX_COORD + cutsy[tcur[XCUTS+i-1]]+1;
			int hi = i == YCUTS ? 2*MAX_COORD : MAX_COORD + cutsy[tcur[XCUTS+i]]+1;
			if (lo == plo[i] && hi == phi[i]) continue;
			plo[i] = lo;
			phi[i] = hi;
			FOR(j, lo, hi) ygroup[j] = i;
		}
	}
	
	REP(i, XCUTS+1) {
		REP(j, YCUTS+1) groups[j] = 0;
		int lo = i == 0 ? 0 : cutsxpos[tcur[i-1]];
		int hi = i == XCUTS ? n : cutsxpos[tcur[i]];
		FOR(j, lo, hi) groups[ygroup[r[j].Y+MAX_COORD]]++;
		REP(j, YCUTS+1) cnt[groups[j]]++;
	}
	
	double rv = 0;
	last_og_score = 0;
	FOR(i, 1, 11) {
		int v = min(cnt[i], a[i]);
		rv += v * pow(i / 10.0, moment*2.85);
		// rv += v;
		last_og_score += v;
	}
	return rv;
}

int main(int argc, char **argv) {
	cin >> n >> k;
	DB(n, k);
	FOR(i, 1, 11) cin >> a[i];
	REP(i, n) cin >> r[i].X >> r[i].Y;
	
	int total = 0;
	FOR(i, 1, 11) total += a[i];
	DB(total);
	
	while ((XCUTS+1)*(YCUTS+1)<total) XCUTS++,YCUTS--;
	DB(XCUTS);
	
	sort(r,r+n);
	REP(i, n-1) if (r[i+1].X > r[i].X) cutsx.PB(r[i].X), cutsxpos.PB(i+1);
	// cutsx.PB(r[n-1].X);
	// cutsxpos.PB(n);
	
	VI vy; REP(i, n) vy.PB(r[i].Y); sort(ALL(vy));
	REP(i, n-1) if (vy[i+1] > vy[i]) cutsy.PB(vy[i]);
	// cutsy.PB(vy.back());
	
	REP(i, XCUTS) cur[i] = rng.next(cutsx.S);
	REP(i, YCUTS) cur[XCUTS+i] = rng.next(cutsy.S);
	
	moment = 1;
	av = eval();
	bv = av;
	xv = 0;
	
	int step = 0;
	int acc  = 0;
	bool ygroup_same = true;
	int maxv[2];
	maxv[0] = cutsx.S-1;
	maxv[1] = cutsy.S-1;
	while (true) {
		step++;
		double time_passed = (get_time() - start_time) / MAX_TIME;
		if (time_passed > 1.0) break;\
		moment = 1 - time_passed;
		
		int type = rng.next_double() < 0.65;
		int d = rng.next(2);
		int p = d == 0 ? rng.next(XCUTS) : rng.next(XCUTS, ACUTS);
		// int v = type == 0 ? rng.next(d == 0 ? cutsx.S : cutsy.S) : max(0, min(maxv[d], cur[p] + rng.next(1, TYPE1_RNG+1)*(2*rng.next(2)-1)));
		int v = type == 0 ? rng.next(d == 0 ? cutsx.S : cutsy.S) : max(0, min(maxv[d], cur[p] + (int)(1+rng.next_double()*rng.next_double()*13.75)*(2*rng.next(2)-1)));
		if (d) ygroup_same = false;
		int old = cur[p];
		if (v == old) continue;
		
		cur[p] = v;
		av = eval(ygroup_same);
		ygroup_same = true;
		
		const double T0 = 0.8;
		const double TN = 0.05;
		double temp = T0 * pow(TN / T0, time_passed);
		
		if (last_og_score > xv) {
			REP(i, ACUTS) best[i] = cur[i];
			xv = last_og_score;
			DB(xv, av, step, time_passed);
		}
		if (av >= bv || rng.next_double() < exp((av-bv)/temp)) {
			acc++;
			bv = av;
		} else {
			cur[p] = old;
			if (d) ygroup_same = false;
		}
	}
	DB(acc);
	DB(step);
	DB(xv);
	double exp_score = 1e6 * xv / total;
	DB(exp_score);
	cout << ACUTS << endl;
	REP(i, XCUTS) cout << cutsx[best[i]] << " " << -FAR << " " << cutsx[best[i]]+1 << " " << FAR << endl;
	REP(i, YCUTS) cout << -FAR << " " << cutsy[best[XCUTS+i]] << " " << FAR << " " << cutsy[best[XCUTS+i]]+1 << endl;
	fflush(stdout);
	return 0;
}