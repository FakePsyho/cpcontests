// Author: Psyho
// Twitter: https://twitter.com/fakepsyho

#pragma GCC optimize("Ofast,fast-math,inline")

// Template:
#include <bits/stdc++.h>
#include <sys/time.h>

#undef assert
#define assert(...) ;
 
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
 
struct FastRNG {
	unsigned long x=123456789, y=362436069, z=521288629;

    FastRNG() { }
	unsigned int rand() { unsigned long t; x ^= x << 16; x ^= x >> 5; x ^= x << 1; t = x; x = y; y = z; z = t ^ x ^ y; return z;}
    INLINE int next() {return rand(); }
    INLINE int next(int x) {return rand() % x; }
    INLINE int next(int a, int b) {return a + (rand() % (b - a)); }
    INLINE double next_double() {return (rand() + 0.5) * (1.0 / 4294967296.0); }
    INLINE double next_double(double a, double b) {return a + next_double() * (b - a); }
};

 
static FastRNG rng;
 
double start_time;
double elapsed() {return get_time() - start_time;}


// Parameters:
const double MAX_TIME = 9.0;
const int BS_RUNS = 20;
const int BS_MUL = 1;

#ifndef LOCAL
const bool USE_TL = true;
#else 
const bool USE_TL = false;
#endif

// Code:
const int MAX_N = 32;
const int GO = MAX_N+1;
const int dd[] = {-MAX_N, -1, 1, MAX_N};

const double MINIMUM = -1e9;

int N, V;
double B;
int grid[MAX_N*MAX_N];

int ds[24][4];

int p1d(int r, int c) {return GO+r*MAX_N+c;}
PII p2d(int p) {return MP((p-GO)/MAX_N,(p-GO)%MAX_N);}
int getd(int p1, int p2) {int dv=p2-p1; return dv==-MAX_N?0:dv==-1?1:dv==1?2:3;}

void generate_dice_states() {
	auto roll_r = [](VI &v) -> VI {return VI{v[5], v[4], v[2], v[3], v[0], v[1]};};
	auto roll_l = [](VI &v) -> VI {return VI{v[4], v[5], v[2], v[3], v[1], v[0]};};
	auto roll_u = [](VI &v) -> VI {return VI{v[3], v[2], v[0], v[1], v[4], v[5]};};
	auto roll_d = [](VI &v) -> VI {return VI{v[2], v[3], v[1], v[0], v[4], v[5]};};
	auto rotate = [](VI &v) -> VI {return VI{v[0], v[1], v[5], v[4], v[2], v[3]};};
	
	VVI base_states = {
		{1,2,3,4,5,6},
		{2,1,6,5,4,3},
		{3,4,2,1,5,6},
		{4,3,1,2,5,6},
		{5,6,3,4,2,1},
		{6,5,3,4,1,2},
	};
	
	map<VI, int> state_map;
	VVI states;
	REP(i, 6) {
		VI state = base_states[i];
		REP(j, 4) {
			state_map[state] = i * 4 + j;
			states.PB(state);
			state = rotate(state);
		}
	}
	REP(i, 24) {
		ds[i][0] = state_map[roll_u(states[i])];
		ds[i][1] = state_map[roll_l(states[i])];
		ds[i][2] = state_map[roll_r(states[i])];
		ds[i][3] = state_map[roll_d(states[i])];
	}
	
	
	//checks:
	VI cnt(24);
	REP(i, 24) REP(j, 4) cnt[ds[i][j]]++;
	REP(i, 24) assert(cnt[i] == 4);
	
	REP(i, 24) REP(j, 4) {
		int s = i;
		REP(k, 4) s = ds[s][j];
		assert(s == i);
	}
	
}

int last_best_sides[6];
bool last_bonus_used;
double path_score(VI &path, bool show=false) {
	if (path.S == 0) return 0;
	
	int cell_cnt[6][10];
	ZERO(cell_cnt);
	
	int state = 0;
	cell_cnt[0][abs(grid[path[0]])] += grid[path[0]];
	int lp = path[0];
	
	FOR(i, 1, path.S) {
		int np = path[i];
		int diff = np - lp;
		int dir = diff == -MAX_N ? 0 : diff == -1 ? 1 : diff == 1 ? 2 : 3;
		state = ds[state][dir];
		cell_cnt[(state>>2)][abs(grid[np])] += grid[np];
		lp = np;
	}
	
	if (show) {
		REP(i, 6) {
			VI v(&cell_cnt[i][1], &cell_cnt[i][V+1]);
			DB(i, v);
		}
	}
	
	double rv = 0;
	REP(i, 6) {
		int bv = -1e9;
		FOR(j, 1, V+1) if (cell_cnt[i][j] > bv) {
			 bv = cell_cnt[i][j];
			 last_best_sides[i] = j;
		}
		rv += bv;
	}
	
	int diff = path[0] - path.back();
	last_bonus_used = false;
	if (diff == -MAX_N || diff == -1 || diff == 1 || diff == MAX_N) {
		rv *= B;
		last_bonus_used = true;
	}
	
	return rv;
}

double zbonus = 0;

int scc_vs_no = 1;
int scc_stack[MAX_N*MAX_N];

const int BIGINT = 0x3F3F3F3F;

double BEST_SCORE = 0;

struct BState {
	int id;
	int start;
	int cur;
	int state;
	int vs[MAX_N*MAX_N];
	int sums[6][10];
	
	void init(int p) {
		memset(vs, 0x3F, sizeof(vs));
		REP(r, N) REP(c, N) vs[p1d(r, c)] = 0;
		ZERO(sums);
		start = p;
		cur = p;
		state = 0;
		vs[p] = BIGINT;
		sums[state>>2][abs(grid[p])] += grid[p];
	}
	
	void copy(BState &s) {
		start = s.start;
		cur = s.cur;
		state = s.state;
		memcpy(vs, s.vs, sizeof(int) * (MAX_N)*(N+2));
		memcpy(sums, s.sums, sizeof(sums));
	}
	
	int go(int p0) {
		if (vs[p0] >= scc_vs_no) return 0;
		
		int stack_no = 1;
		scc_stack[0] = p0;
		int rv = 0;
		while (stack_no) {
			int p = scc_stack[--stack_no];
			rv += grid[p] >= max(1, V-5) ? grid[p] : 0;
			if (vs[p-MAX_N] < scc_vs_no) scc_stack[stack_no++] = p-MAX_N, vs[p-MAX_N] = scc_vs_no;
			if (vs[p-1]     < scc_vs_no) scc_stack[stack_no++] = p-1,     vs[p-1]     = scc_vs_no;
			if (vs[p+1]     < scc_vs_no) scc_stack[stack_no++] = p+1,     vs[p+1]     = scc_vs_no;
			if (vs[p+MAX_N] < scc_vs_no) scc_stack[stack_no++] = p+MAX_N, vs[p+MAX_N] = scc_vs_no;
		}
		return rv;
	}
	
	void move(int d) {
		assert(d >= 0 && d <= 3);
		int np = cur + dd[d];
		assert(vs[np] < BIGINT);
		vs[np] = BIGINT;
		cur = np;
		state = ds[state][d];
		sums[state>>2][abs(grid[np])] += grid[np];
	}
	
	double eval() {
		double rv = 0;
		REP(i, 6) {
			int bv0 = -1e9;
			int bv1 = -1e9;
			int bv2 = -1e9;
			FOR(j, max(1, V-4), V+1) {
				int av = sums[i][j];
				if (av > bv0) {
					bv2 = bv1;
					bv1 = bv0;
					bv0 = av;
				} else if (av > bv1) {
					bv2 = bv1;
					bv1 = av;
				} else if (av > bv2) {
					bv2 = av;
				}
			}
			rv += bv0;
		}
		int options = 0;
		scc_vs_no++;
		
		int value_left = 0;
		REP(d, 4) {
			int np = cur + dd[d];
			value_left = max(value_left, go(np));
		}
		
		bool start_connected = false;
		REP(d, 4) start_connected |= vs[start + dd[d]] == scc_vs_no || start + dd[d] == cur;
		
		if ((rv + value_left) * B < BEST_SCORE) return MINIMUM - 1;
		
		rv += value_left * 0.5;
		if (start_connected) rv *= B;
		
		auto p = p2d(cur);
		rv += min(min(p.X, p.Y), min(N-1-p.X, N-1-p.Y)) * zbonus;
		rv += rng.next_double() * 5.0;
		
		return rv;
	}
	
	double true_eval() {
		double rv = 0;
		REP(i, 6) {
			int bv = -1e9;
			FOR(j, 1, V+1) bv = max(bv, sums[i][j]); 
			rv += bv;
		}
		
		auto s2d = p2d(start);
		auto c2d = p2d(cur);
		if (abs(s2d.X - c2d.X) + abs(s2d.Y - c2d.Y) == 1) rv *= B;
		
		return rv;
	}
	
	double eval_move(int d) {
		int side = ds[state][d]>>2;
		int old_cur = cur;
		int np = cur + dd[d];
		int v = grid[np];
		sums[side][abs(v)] += v;
		vs[np] = BIGINT;
		cur = np;
		double rv = eval();
		cur = old_cur;
		vs[np] = 0;
		sums[side][abs(v)] -= v;
		return rv;
	}
	
};

VI last_path;

const int MAX_BEAM_WIDTH = 30000;
BState bs0[MAX_BEAM_WIDTH];
BState bs1[MAX_BEAM_WIDTH];
double nbs_eval[3*MAX_BEAM_WIDTH];
int order[3*MAX_BEAM_WIDTH];
PII moves[MAX_BEAM_WIDTH*30*30];
PII bmoves[3*MAX_BEAM_WIDTH];
double beam_search(VI starts, double time_mul = 1) {
	int BEAM_WIDTH = min(MAX_BEAM_WIDTH, (int)(2400 * time_mul * BS_MUL / BS_RUNS * 30 * 30 * 30 / N / N / N));
	
	zbonus = -(2*V+1);
	
	last_path.clear();
	if (starts.S == 0) return MINIMUM;
	
	
	int bs_no = 0;
	int nbs_no = 0;
	int max_id = 0;
	
	auto bs = bs0;
	auto nbs = bs1;
	
	double bv = 0;
	int bid = -1;
	
	REP(i, starts.S) {
		int p = starts[i];
		bs[i].init(p);
		bs[i].id = i;
		moves[i] = MP(-1, p);
		max_id++;
		bs_no++;
	}
	
	
	int level = 0;
	while (true) {
		level++;
		nbs_no = 0;
		
		REP(i, bs_no) {
			BState &s = bs[i];
			
			double v = s.true_eval();
			if (v > bv) {
				bv = v;
				bid = s.id;
			}
			REP(d, 4) {
				int np = s.cur + dd[d];
				if (s.vs[np] == BIGINT) continue;
				nbs_eval[nbs_no] = bs[i].eval_move(d);
				if (nbs_eval[nbs_no] < MINIMUM) continue;
				bmoves[nbs_no] = MP(i, d);
				nbs_no++;
			}
		}
		
		if (nbs_no == 0) break;
		
		REP(i, nbs_no) order[i] = i;
		sort(order, order+nbs_no, [&](int a, int b) -> bool {return nbs_eval[a] > nbs_eval[b];});
		
		nbs_no = min(BEAM_WIDTH, nbs_no);
		REP(i, nbs_no) {
			PII &m = bmoves[order[i]];
			nbs[i].copy(bs[m.X]);
			nbs[i].move(m.Y);
			nbs[i].id = max_id;
			moves[max_id] = MP(bs[m.X].id, nbs[i].cur);
			max_id++;
		}
		
		swap(bs, nbs);
		bs_no = nbs_no;
	}
	
	last_path.clear();
	while (bid >= 0) {
		last_path.PB(moves[bid].Y);
		bid = moves[bid].X;
	}
	reverse(ALL(last_path));
	
	return bv;
}


int main() {
	//Read data
	cin >> N >> V >> B;
	start_time = get_time();
	REP(r, N) REP(c, N) cin >> grid[p1d(r, c)];
	
	cerr << "[DATA] N = " << N << endl;
	cerr << "[DATA] V = " << V << endl;
	cerr << "[DATA] B = " << B << endl;
	
	//Generate dice states
	generate_dice_states();
	
	//Find solution
	double bv = 0;
	VI bpath;
	
	VI pts;
	REP(r, N) REP(c, N) pts.PB(p1d(r, c));
	random_shuffle(ALL(pts));
	REP(i, BS_RUNS) {
		VI starts(pts.begin() + i * pts.S / BS_RUNS, pts.begin() + (i+1) * pts.S / BS_RUNS);
		
		double av = beam_search(starts);
		DB(i, av, elapsed());
		if (av > bv) {
			bpath = last_path; 
			bv = av;
			BEST_SCORE = av;
		}
		
		if (USE_TL && elapsed() > MAX_TIME) break;
	}
	
	if (!USE_TL || elapsed() < MAX_TIME) {
		double av = beam_search(VI{bpath[0]}, 1.0);
		DB(av);
		if (av > bv) {
			bpath = last_path;
			bv = av;
			BEST_SCORE = av;
		}
	}
	
	DB(bv);
	DB(elapsed());
	
	//Output solution
	DB(path_score(bpath, true));
	cerr << "[DATA] bonus = " << last_bonus_used << endl;
	REP(i, 6) cout << last_best_sides[i < 2 ? i^1 : i] << endl;
	cout << bpath.S << endl;
	for (int p : bpath) cout << p2d(p).X << ' ' << p2d(p).Y << endl;
	
	return 0;	
}