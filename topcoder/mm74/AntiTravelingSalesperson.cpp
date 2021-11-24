#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <map>
#include <set>
#include <list>
#include <deque>
#include <queue>

#ifdef LOCAL
#define TIMESCALE 0.5
#else
#define TIMESCALE 1.0
#endif

#define MAXTIME (9.8 * TIMESCALE)

using namespace std;

#define FOR(i,a,b)  for(int i=(int)(a);i<(int)(b);++i)
#define REP(i,a)    FOR(i,0,a)
#define ZERO(m)     memset(m,0,sizeof(m))
#define ALL(x)      x.begin(),x.end()
#define PB          push_back
#define S           size()
#define LL          long long
#define MP          make_pair
#define X           first
#define Y           second
#define VC          vector
#define VD            VC < double >
#define VVD            VC < VD >
#define PII         pair <int, int>
#define VI          VC < int >
#define VVI            VC < VI >
#define VPII        VC < PII >
#define VS          VC<string>
#define DB(a)        cerr << #a << ": " << a << endl;

void print(VI v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]\n";}
void print(VD v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]\n";}
void print(VS v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]\n";}
template<class T> string i2s(T x) {ostringstream o; o << x; return o.str(); }
VS splt(string s, char c = ' ') {VS rv; int p = 0, np; while (np = s.find(c, p), np >= 0) {if (np != p) rv.PB(s.substr(p, np - p)); p = np + 1;} if (p < (int)s.S) rv.PB(s.substr(p)); return rv;}

double getTime() {
    unsigned LL time;
    __asm__ volatile ("rdtsc" : "=A" (time));
#ifdef LOCAL
    return time / 1.66e9; 
#else
    return time / 3.6e9;
#endif
}

double randDouble() {
    return ((double)(rand() + 0.5) / ((double)RAND_MAX + 1));
}

#define MAXN 10000
#define MAXD 1000000000

int N;
int X[MAXN];
int Y[MAXN];
int TX;
int TY;

double startTime = 0;

void save(int p) {
	TX = X[p];
	TY = Y[p];
}

void restore(int p) {
	X[p] = TX;
	Y[p] = TY;
}

void setRandom(int p) {
	X[p] = rand() % (MAXD + 1);
	Y[p] = rand() % (MAXD + 1);
}

void updateRandom(int p, int d) {
	save(p);
	do {
		restore(p);
		X[p] += rand() % (d + d + 1) - d;
		Y[p] += rand() % (d + d + 1) - d;
	} while (X[p] < 0 || X[p] > MAXD || Y[p] < 0 || Y[p] > MAXD);
}

LL dist(int a, int b) {
	return (LL)(X[a] - X[b]) * (X[a] - X[b]) + (LL)(Y[a] - Y[b]) * (Y[a] - Y[b]);
}

int USED[MAXN];
double test() {
	double bv = 0;
	//REP(i, N) {
	int i = 0;
		memset(USED, 0, N * sizeof(int));
		int p = i;
		double av = 0;
		USED[p] = true;
		REP(j, N - 1) {
			int bp = -1;
			LL bd = 1LL << 61;
			REP(k, N) if (!USED[k]) {
				LL ad = dist(p, k);
				if (ad < bd) bd = ad, bp = k;
			}
			p = bp;
			USED[p] = true;
			av += sqrt(bd);
		}
		av += sqrt(dist(p, i));
		bv = max(bv, av);
	//}
	return bv;	
}

class AntiTravelingSalesperson{public: VI placeLocations(int N, VI SX, VI SY) {
	startTime = getTime();
	
	::N = N;
	
	REP(i, N) setRandom(i);
	int steps = 0;	
	double bv = test();
	double xv = bv;
	
	srand(5);
	
	while (getTime() - startTime < MAXTIME) {
		steps++;
		
		double timeLeft = (MAXTIME - (getTime() - startTime)) / MAXTIME;
		//double temp = 500 + 100000 * timeLeft;
		double temp = 0;
		
		int p = rand() % N;
		save(p);
		//rand() % 2 ? updateRandom(p, 1000000) : setRandom(p);
		//setRandom(p);
		updateRandom(p, (int)(pow(10, 7 + timeLeft * 2)));
		double av = test();
		
		if (av > bv || randDouble() < exp((av - bv) / temp)) {
			bv = av;
			if (bv > xv) {
				xv = bv;
				//DB(temp);
				fprintf(stderr, "%d %.0f\n", steps, xv);
			}
		} else {
			restore(p);
		}
	}
	
	VI rv;
	REP(i, N) rv.PB(X[i]), rv.PB(Y[i]);
	return rv;
}};


#ifdef LOCAL
int main() {
	int N; cin >> N;

	int XN; cin >> XN;
	VI X(XN); REP(i, XN) cin >> X[i];

	int YN; cin >> YN;
	VI Y(YN); REP(i, YN) cin >> Y[i];

	AntiTravelingSalesperson algo;
	VI rv = algo.placeLocations(N, X, Y);

	cout << rv.S << endl;
	REP(i, rv.S) cout << rv[i] << endl;
    fflush(stdout);
	
	return 0;
}
#endif