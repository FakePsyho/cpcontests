// Author: Psyho
// Twitter: https://twitter.com/fakepsyho

#pragma GCC optimize("Ofast,fast-math,inline")
 
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

// Parameters
const double MAX_TIME = 9.8;


// Code

double elapsed() {return get_time() - start_time;}

const int MAX_N = 32;
const int MAX_P = 30;
const int GOFFSET = MAX_N;

const int dd[] = {-MAX_N, -1, 1, MAX_N};
int dt[24][4];


int p1d(int r, int c) {return GOFFSET+r*MAX_N+c;}
PII p2d(int p) {return MP((p-GOFFSET)/MAX_N,(p-GOFFSET)%MAX_N);}

int N, C, P;
int N2;
int gv[MAX_N*MAX_N];
int gc[MAX_N*MAX_N];

VVI endp;

const int MAX_COST = MAX_N*MAX_N/2;

int q[MAX_COST][MAX_N*MAX_N];
int qno[MAX_COST];
int qvs[MAX_N*MAX_N];
int qvs_no = 1;
int vprev[MAX_N*MAX_N];
int vdist[MAX_N*MAX_N][MAX_N*MAX_N];

PII vends[MAX_N*MAX_N];
int vends_no = 0;
struct State {
	int vs[MAX_N*MAX_N];
	int used[MAX_N*MAX_N];
	int enabled[MAX_N*MAX_N];
	VPII conn;
	VVI paths;
	int score;
	
	State() { 
		ZERO(vs);
		ZERO(used);
		ZERO(enabled);
	}
	
	VI find_path(int a) {
		assert(!used[a]);
		assert(gv[a]);
		
		qvs_no++;
		
		const int max_cost = (int)((N*N)*0.30);
		q[0][0] = a;
		vprev[a] = -1;
		qno[0] = 1;
		qvs[a] = qvs_no;
		
		int qc = 0;
		int qv = 0;
		
		int *cdt = dt[rng.next(24)];
		
		const int penalty = 1+rng.next(1+N*P);
		
		vends_no = 0;
		while (qv <= qc) {
			REP(i, qno[qv]) {
				int v = q[qv][i];
				if (gc[v] && v != a) {
					vends[vends_no++] = MP(v, qv);
					continue;
				}
				REP(d, 4) {
					int nv = v + cdt[d];
					if (qvs[nv] >= qvs_no) continue;
					qvs[nv] = qvs_no;
					vprev[nv] = v;
					if (!vs[nv] && (gc[nv] == gc[a] || gc[nv] == 0)) {
						q[qv+1][qno[qv+1]] = nv;
						qno[qv+1]++;
						qc = max(qc, qv+1);
					} else if (gc[nv] == 0) {
						int nc = qc + penalty * vs[nv];
						if (nc > max_cost) continue;
						q[nc][qno[nc]] = nv;
						qno[nc]++;
						qc = max(qc, nc);
					}
				}
			}
			qno[qv] = 0;
			qv++;
		}
		
		if (vends_no == 0) return VI();
		
		int b = -1;
		double bv = 0;
		REP(i, vends_no) {
			PII &p = vends[i];
			double av = 1.0 * gv[p.X] * gv[p.X] * gv[p.X] / sqrt(N+p.Y) * sqrt(rng.next_double());
			if (av > bv) {
				b = p.X;
				bv = av;
			}
		}
		
		VI path;
		while (b != -1) {
			path.PB(b);
			if (path.S > 4*N) return VI();
			b = vprev[b];
		}
		reverse(ALL(path));
		return path;
	}
	
	void add_path(int a, int b, VI &path) {
		assert(!used[a]);
		assert(!used[b]);
		conn.PB(MP(a, b));
		enabled[paths.S] = 1;
		used[a] = 1;
		used[b] = 1;
		paths.PB(path);
		for (int v : path) vs[v]++;
	}
	
	void remove_path(int id) {
		for (int v : paths[id]) vs[v]--;
		used[conn[id].X] = 0;
		used[conn[id].Y] = 0;
		conn[id] = conn.back();
		paths[id] = paths.back();
		conn.pop_back();
		paths.pop_back();
		enabled[id] = enabled[paths.S];
	}
	
	void remove_disabled_paths() {
		REP(id, paths.S) if (!enabled[id]) {
			conn[id] = conn.back();
			paths[id] = paths.back();
			conn.pop_back();
			paths.pop_back();
			enabled[id] = enabled[paths.S];
			id--;
		}
	}
	
	void disable_path(int id) {
		assert(enabled[id]);		
		for (int v : paths[id]) vs[v]--;
		used[conn[id].X] = 0;
		used[conn[id].Y] = 0;
		enabled[id] = 0;
	}
	
	void reenable_paths() {
		REP(i, paths.S) if (!enabled[i]) {
			used[conn[i].X] = 1;
			used[conn[i].Y] = 1;
			for (int v : paths[i]) vs[v]++;
			enabled[i] = 1;
		}
	}
	
	double eval() {
		double rv = 0;
		REP(i, conn.S) if (enabled[i]) rv += gv[conn[i].X] * gv[conn[i].Y];
		REP(r, N) {
			int p = p1d(r,0);
			REP(c, N) if (vs[p+c] > 1) rv -= vs[p+c] * (vs[p+c]-1) / 2 * P;
		}
		return rv;
	}
	
	void print_solution() {
		cout << paths.S << endl;
		for (VI &path : paths) {
			cout << path.S << endl;
			for (int &v : path) cout << p2d(v).X << " " << p2d(v).Y << endl;
		}
	}
	
	
};

int main() {
	// Read Input
	cin >> N >> C >> P;
	REP(r, N) REP(c, N) cin >> gv[p1d(r,c)] >> gc[p1d(r,c)];
	start_time = get_time();
	
	cerr << "[DATA] N = " << N << endl;
	cerr << "[DATA] C = " << C << endl;
	cerr << "[DATA] P = " << P << endl;
	
	N2 = (MAX_N)*(N+2);
	memset(qvs, 0x3F, sizeof(qvs));
	REP(r, N) REP(c, N) qvs[p1d(r,c)] = 0;
	
	// REP(r1, N) REP(c1, N) REP(r2, N) REP(c2, N) vdist[p1d(r1,c1)][p1d(r2,c2)] = abs(r1-r2)+abs(c1-c2);
	REP(r1, N) REP(c1, N) REP(r2, N) REP(c2, N) vdist[p1d(r1,c1)][p1d(r2,c2)] = max(abs(r1-r2),abs(c1-c2));
	
	endp = VVI(C+1);
	REP(i, N2) if (gc[i]) endp[gc[i]].PB(i);
	
	REP(i, C+1) DB(endp[i].S);
	
	VI vd = {dd[0], dd[1], dd[2], dd[3]};
	sort(ALL(vd));
	REP(i, 24) {
		REP(j, 4) dt[i][j] = vd[j];
		next_permutation(ALL(vd));
	}
	
	// Find Solution
	double bv = 0;
	double xv = 0;
	
	State xsol;
	State s;
	int step = 0;
	int acc  = 0;
	
	const double t0 = 30.0;
	const double tn = 4.0;
	double t = t0;
	double time_passed = 0;
	while (true) {
		if ((step & 255) == 0) {
			time_passed = elapsed() / MAX_TIME;
			if (time_passed > 1.0) break;
			t = t0 * pow(tn / t0, time_passed);
		}
		step++;
		
		int color = 1+rng.next(C);
		
		
		//remove paths
		
		int remove_no = rng.next(min(3, (int)s.paths.S+1));
		
		while (remove_no--) {
			int id = rng.next(s.paths.S);
			while (!s.enabled[id]) id = rng.next(s.paths.S);
			if (gc[s.conn[id].X] != color) do {id = rng.next(s.paths.S);} while (!s.enabled[id]);
			if (gc[s.conn[id].X] != color) do {id = rng.next(s.paths.S);} while (!s.enabled[id]);
			if (gc[s.conn[id].X] != color) do {id = rng.next(s.paths.S);} while (!s.enabled[id]);
			if (gc[s.conn[id].X] != color) do {id = rng.next(s.paths.S);} while (!s.enabled[id]);
			if (gc[s.conn[id].X] != color) do {id = rng.next(s.paths.S);} while (!s.enabled[id]);
			assert(s.enabled[id]);
			s.disable_path(id);
		}
		
		//add paths
		
		int paths_added = 0;
		
		static VI unused;
		while (true) {
			int a = -1;
			if (paths_added == 3) break;
			if (paths_added >= 1 && rng.next_double() < 0.5) break;
			
			unused.clear();
			for (int e : endp[color]) if (!s.used[e]) unused.PB(e);
			if (unused.S < 2) break;
			REP(i, 4) {
				int c = unused[rng.next(unused.S)];
				if (a == -1 || gv[c] > gv[a]) a = c;
			}
			VI path = s.find_path(a);
			if (path.S == 0) break;
			int b = path.back();
			s.add_path(a, b, path);
			
			paths_added++;
		}
		
		if (paths_added == 0) {
			s.reenable_paths();
			continue;
		}
		
		double av = s.eval();
		if (av >= bv || rng.next_double() < exp((av-bv)/t)) {
			s.remove_disabled_paths();
			if (av > xv) {
				xv = av;
				DB(step, xv, t, elapsed());
				xsol = s;
			}
			bv = av;
			acc++;
		} else {
			REP(i, paths_added) s.remove_path(s.paths.S - 1);
			s.reenable_paths();
		}
	}
	
	cerr << "[DATA] step = " << step << endl;
	cerr << "[DATA] acc = " << acc << endl;
	
	DB(xsol.eval());
	
	xsol.print_solution();
	
	return 0;

}
