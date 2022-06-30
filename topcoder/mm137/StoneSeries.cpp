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


#ifdef VM
const double TIME_MUL = 0.66;
#else
const double TIME_MUL = 1.0;
#endif

// Parameters:

const bool TRY_LONGEST_PATH = false;

const bool FIXED_TEMP = false;

const double MAX_TIME = 9.9 * TIME_MUL;

const double TIME_GREEDY = 0.08;
const double TIME_PRESERVE_SA = 0.1;

// Code:
const int MOD = 1e9+7;

const int MAX_ONES = 160000 / 8 / 8;
const int MAX_N = 48;
const int MAX_CELLS = 1600;

#define ITYPE short

int N;
ITYPE blocked[MAX_N*MAX_N];

INLINE int pid(int r, int c) {return 2*MAX_N+r*(N+2)+c;}
INLINE string pstr(int p) {return i2s((p-2*MAX_N)/(N+2)) + " " + i2s((p-2*MAX_N)%(N+2));}
INLINE PII ppii(int p) {return MP((p-2*MAX_N)/(N+2), (p-2*MAX_N)%(N+2));}

int dde[12];
int dd[8];
int dds[4];
int ddx[4];

//state data

ITYPE cells[MAX_CELLS][MAX_CELLS];
ITYPE g[MAX_N*MAX_N];
int sum[MAX_N*MAX_N];
ITYPE n_cells[MAX_CELLS];
pair<ITYPE,ITYPE> cell_pos[MAX_N*MAX_N];
bool v_cnt[MAX_CELLS];
ITYPE maxv;

int last_score;

VI cempty;

VI sol;


void set_option(int p, int v) {
	cell_pos[p] = MP(v, n_cells[v]);
	cells[v][n_cells[v]] = p;
	n_cells[v]++;
}

ITYPE g_backup[MAX_N*MAX_N];

PII gremoved[MAX_CELLS];
int added_no = 0;
int gadded[MAX_CELLS];
int removed_no = 0;

void clear_state() {
	memset(&g[2*MAX_N], 0, (N+2)*N*sizeof(ITYPE));
}

void save_state() {
	memcpy(&g_backup[2*MAX_N], &g[2*MAX_N], (N+2)*N*sizeof(ITYPE));
}

void restore_state() {
	memcpy(&g[2*MAX_N], &g_backup[2*MAX_N], (N+2)*N*sizeof(ITYPE));
}

int score_state() {
	int rv = 0;
	FOR(i, 2*MAX_N, 2*MAX_N+(N+2)*N) rv += g[i];
	return rv;
}

INLINE int greedy() {
	const int max_len = N*N;
	
	memset(n_cells, 0, max_len*sizeof(ITYPE));
	memset(v_cnt, 0, max_len*sizeof(bool));
	
	memset(&sum[2*MAX_N], 0, (N+2)*N*sizeof(int));
	
	int score = 0;
	
	for (int p : cempty) {
		if (g[p]) {
			score += g[p];
			v_cnt[g[p]] = true;
			for (int d : dd) sum[d+p] += (1<<16) + g[p];
		} else {
			cell_pos[p].X = -1;
		}
	}
	
	maxv = 1;
	while (v_cnt[maxv+1]) maxv++;
	
	for (int p : cempty) if ((sum[p]>>16) != 1 && g[p] == 0 && (sum[p]&65535) < max_len) set_option(p, max(1, sum[p]&65535));
	
	for (int v = maxv + 1; v >= 1; v--) {
		if (n_cells[v]) {
			int pos = rng.next(n_cells[v]);
			int p = cells[v][pos];
			g[p] = max(1, sum[p]&65535);
			gadded[added_no++] = p;
			
			v_cnt[v] = true;
			while (v_cnt[maxv+1]) maxv++;
			
			score += g[p];
			
			int id = cell_pos[p].Y;
			n_cells[v]--;
			cells[v][id] = cells[v][n_cells[v]];
			cell_pos[cells[v][id]].Y = id;
			cell_pos[p].X = -1;
			for (int d : dd) {
				int np = d + p;
				if (!blocked[np] && g[np] == 0) {
					sum[np] += (1<<16)+g[p];
					
					if (cell_pos[np].X >= 0) {
						int v0 = cell_pos[np].X;
						int id = cell_pos[np].Y;
						cells[v0][id] = cells[v0][--n_cells[v0]];
						cell_pos[cells[v0][id]].Y = id;
					}
					
					int v0 = sum[np]&65535;
					if (sum[np] >= (2<<16) && v0 < max_len) {
						cell_pos[np] = MP(v0, n_cells[v0]);
						cells[v0][n_cells[v0]] = np;
						n_cells[v0]++;
					} else {
						cell_pos[np].X = -1;
					}
				}
			}
			v = maxv + 2;
		} 
	}
	
	for (int p : cempty) if (g[p] > maxv) {
		gremoved[removed_no++] = MP(p, g[p]);
		score -= g[p];
		g[p] = 0;
	}
	
	last_score = score;
	return score;
}

int v_order[MAX_CELLS];
INLINE int notgreedy() {
	const int max_len = N*N;
	
	memset(n_cells, 0, max_len*sizeof(ITYPE));
	memset(v_cnt, 0, max_len*sizeof(bool));
	
	memset(&sum[2*MAX_N], 0, (N+2)*N*sizeof(int));
	
	int score = 0;
	
	for (int p : cempty) {
		if (g[p]) {
			score += g[p];
			v_cnt[g[p]] = true;
			for (int d : dd) sum[d+p] += (1<<16) + g[p];
		} else {
			cell_pos[p].X = -1;
		}
	}
	
	maxv = 1;
	while (v_cnt[maxv+1]) maxv++;
	
	for (int p : cempty) if ((sum[p]>>16) != 1 && g[p] == 0 && (sum[p]&65535) < max_len) set_option(p, max(1, sum[p]&65535));

	REP(i, N*N/2) v_order[i] = i+2;
	
	int bad = -2;
	while (true) {
		int v;
		if (bad == -2) {
			v = maxv+1;
		} else if (bad == -1) {
			v = 1;
		} else {
			swap(v_order[bad], v_order[rng.next(bad, maxv-1)]);
			v = v_order[bad];
		}
		bad++;
		if (n_cells[v]) {
			bad = -2;
			int pos = rng.next(n_cells[v]);
			int p = cells[v][pos];
			g[p] = max(1, sum[p]&65535);
			gadded[added_no++] = p;
			
			v_cnt[v] = true;
			while (v_cnt[maxv+1]) maxv++;
			
			score += g[p];
			
			int id = cell_pos[p].Y;
			n_cells[v]--;
			cells[v][id] = cells[v][n_cells[v]];
			cell_pos[cells[v][id]].Y = id;
			cell_pos[p].X = -1;
			for (int d : dd) {
				int np = d + p;
				if (!blocked[np] && g[np] == 0) {
					sum[np] += (1<<16)+g[p];
					
					if (cell_pos[np].X >= 0) {
						int v0 = cell_pos[np].X;
						int id = cell_pos[np].Y;
						cells[v0][id] = cells[v0][--n_cells[v0]];
						cell_pos[cells[v0][id]].Y = id;
					}
					
					int v0 = sum[np]&65535;
					if (sum[np] >= (2<<16) && v0 < max_len) {
						cell_pos[np] = MP(v0, n_cells[v0]);
						cells[v0][n_cells[v0]] = np;
						n_cells[v0]++;
					} else {
						cell_pos[np].X = -1;
					}
				}
			}
		} else if (bad >= maxv-1) {
			break;
		}
	}
	
	for (int p : cempty) if (g[p] > maxv) {
		gremoved[removed_no++] = MP(p, g[p]);
		score -= g[p];
		g[p] = 0;
	}
	
	last_score = score;
	return score;
}

void remove_tree(int p0) {
	if (g[p0] == 0) return;
	
	static VI stack;
	stack.PB(p0);
	while (stack.S) {
		int p = stack.back();
		stack.pop_back();
		if (g[p] == 0) continue;
		gremoved[removed_no++] = MP(p, g[p]);
		for (int d : dd) {
			int np = p + d;
			if (g[np] > g[p]) {
				stack.PB(np);
			}
		}
		g[p] = 0;
	}
}

void reconstruct_solution() {
	sol.clear();
	static VPII vp;
	vp.clear();
	for (int p : cempty) if (g[p]) vp.PB(MP(g[p], p));
	sort(ALL(vp));
	for (PII &p : vp) sol.PB(p.Y);
}

void reconstruct_state() {
	clear_state();
	for (int p : sol) {
		for (int d : dd) g[p] += g[p+d];
		if (g[p] == 0) g[p] = 1;
	}
}

ITYPE pre_g[MAX_ONES][MAX_N*MAX_N];
int one_state[MAX_N*MAX_N];
int one_value[MAX_N*MAX_N];
int one_pos[MAX_N*MAX_N];
int one_fc[MAX_N*MAX_N];
int one_fc_pos = 1;
int one_evbad[MAX_N*MAX_N];
int one_evbad_pos = 1;
int build_ones(bool diagonal_ok = false) {
	int rv = 0;
	FOR(r, -1, N+1) FOR(c, -1, N+1) {
		int p = pid(r, c);
		one_state[p] = blocked[p];
		one_value[p] = 0;
		one_pos[p] = -1;
		if (blocked[p]) continue;
		for (int d : dd) {
			int np = p + d;
			one_value[p] += !blocked[np];
		}
	}
	
	int step = 0;
	int last_imp = 0;
	one_evbad_pos++;
	while (true) {
		step++;
		if (step - last_imp > N*N*10) break;
		
		int p = cempty[rng.next(cempty.S)];
		if (one_state[p] || one_evbad[p] == one_evbad_pos) continue;
		
		one_fc_pos++;
		int nv = one_value[p];
		assert(nv >= 0 && nv <= 8);
		for (int d : dd) {
			int np = d + p;
			if (one_pos[np] == -1 || one_fc[one_pos[np]] == one_fc_pos) continue;
			nv -= one_value[one_pos[np]];
			one_fc[one_pos[np]] = one_fc_pos;
			if (nv < 0) break;
		}
		
		if (nv >= 0) {
			one_evbad_pos++;
			rv += nv;
			for (int d : dd) {
				int np = d + p;
				if (one_pos[np] == -1) continue;
				int p2 = one_pos[np];
				assert(one_state[p2] == 2);
				for (int d2 : dd) one_pos[p2 + d2] = -1;
				one_pos[p2] = -1;
				one_state[p2] = 0;
			}
			one_state[p] = 2;
			one_pos[p] = p;
			for (int d : dd) {
				int np = d + p;
				if (!blocked[np]) one_pos[np] = p;
			}
			if (nv > 0) last_imp = step;
		} else {
			one_evbad[p] = one_evbad_pos;
		}
		
	}
	
	return rv;
}

int main() {
	// Read data & Init Stuff
	cin >> N;
	start_time = get_time();
	
	dd[0] = -N-3;
	dd[1] = -N-2;
	dd[2] = -N-1;
	dd[3] = -1;
	dd[4] = +1;
	dd[5] = +N+1;
	dd[6] = +N+2;
	dd[7] = +N+3;
	REP(i, 4) dds[i] = dd[1 + 2*i];
	REP(i, 4) ddx[i] = dd[0 + 2*i];
	REP(i, 8) dde[i] = dd[i];
	dd[8] = -2;
	dd[9] = +2;
	dd[10] = (N+2)*-2;
	dd[11] = (N+2)*+2;
	
	int n_blocked = 0;
	REP(i, MAX_N*MAX_N) blocked[i] = 1;
	REP(r, N) REP(c, N) {
		char x; cin >> x;
		int p = pid(r, c);
		blocked[p] = x == '#';
		n_blocked += blocked[p];
	}	  
	cerr << "[DATA] N = " << N << endl;
	cerr << "[DATA] ratio = " << 1.0*(N*N - n_blocked) / (N*N) << endl;
	
	REP(r, N) REP(c, N) {
		int p = pid(r, c);
		if (!blocked[p]) cempty.PB(p);
	}
	
	const int N_ONES = 160000 / N / N;
	// Pre-compute 1-cover
	REP(i, N_ONES) {
		build_ones();
		for (int p : cempty) pre_g[i][p] = one_state[p] == 2;
	}
	DB(N_ONES, elapsed());
	
	//Runs
	int small_runs = N <= 11 ? 50 : N <= 14 ? 20 : N <= 20 ? 8 : N <= 28 ? 2 : 1;
	int bonus_runs = N > 20 ? 0 : 0;
	int n_runs = small_runs + bonus_runs;
	DB(N, n_runs);
	
	int greedy_step = 0;
	int step = 0;
	int acc = 0;
	int xmaxv = 0;
	
	int xv = 0;
		
	int zero_removals = 0;
	int total_removed = 0;
	int total_removals = 0;
		
	const double T0 = N*sqrt(N)*0.8;
	const double TN = min(T0, N*N*N/400.0);
	DB(T0, TN);
	REP(cur_run, n_runs) {
		bool bonus_run = cur_run >= small_runs;
		
		double run_start = elapsed();
		double run_end = bonus_run ? MAX_TIME : MAX_TIME * (cur_run + 1) / (n_runs + bonus_runs * 2);
		
		int bv = 0;
		
		
		if (!bonus_run) {
			double greedy_done = run_start + (run_end - run_start) * TIME_GREEDY;
			int n_ones_lo = cur_run * N_ONES / small_runs;
			int n_ones_hi = (cur_run + 1) * N_ONES / small_runs;
			while (true) {
				double time_passed = elapsed();
				if (time_passed > greedy_done) break;
				
				added_no = 0;
				clear_state();
				int ones_p = n_ones_lo + greedy_step % (n_ones_hi - n_ones_lo);
				for (int p : cempty) g[p] = pre_g[ones_p][p];
				int av = greedy();
				if (av > bv) {
					bv = av;
					save_state();
					if (av > xv) {
						xv = av;
						reconstruct_solution();
					}
				}
				greedy_step++;
			}
			DB(bv, elapsed());
		} else {
			reconstruct_state();
			DB(score_state());
		}
		
		double time_passed = 0;
		
		restore_state();
		double sa_start = elapsed();
		double sa_end = run_end;
		double t = T0;
		
		int zv = bv;
		while (true) {
			step++;
			
			if ((step & 127) == 0) {
				time_passed = elapsed();
				if (time_passed > sa_end) break;
				time_passed = (time_passed - sa_start) / (sa_end - sa_start);
				t = T0 + (TN - T0) * time_passed;
			}
			
			int big = rng.next_double() < .30;
			
			removed_no = 0;
			if (big) {
				int size = rng.next(2, 6);
				int r0 = rng.next(-size + 2, N-1);
				int c0 = rng.next(-size + 2, N-1);
				
				FOR(r, max(0, r0), min(N, r0 + size)) FOR(c, max(0, c0), min(N, c0 + size)) {
					int np = pid(r, c);
					if (g[np] > 1 || time_passed > TIME_PRESERVE_SA) remove_tree(np);
				}
			} else {
				int no = N > 28 ? rng.next(1, 3) : N > 24 ? rng.next(1, 4) : N > 14 ? rng.next(1, 5) : rng.next(1, 6);
				while (no--) {
					int tries = N*N;
					int p = 0;
					int repeats_no = rng.next_double() < time_passed * .5 + .80 ? 2 : 1;
					REP(repeats, repeats_no) {
						int tp; do {tp = cempty[rng.next(cempty.S)]; tries--;} while (g[tp] == 0 && tries > 0);
						if (g[tp] > g[p]) p = tp;
					}
					if (time_passed < TIME_PRESERVE_SA && g[p] <= 1) continue;
					if (g[p] == 0) continue;
					remove_tree(p);
				}
			}
			total_removed += removed_no;
			total_removals++;
			zero_removals += removed_no == 0;
			
			if (removed_no == 0) continue;
			added_no = 0;
			double av = rng.next_double() < .75 ? greedy() : notgreedy();
			if (av >= bv || rng.next_double() < exp((av - bv) / t)) {
				acc++;
				if (last_score > xv) {
					xv = last_score;
					xmaxv = maxv;
					reconstruct_solution();
				}
				bv = av;
				zv = max(zv, bv);
			} else {
				REP(i, added_no) g[gadded[i]] = 0;
				REP(i, removed_no) g[gremoved[i].X] = gremoved[i].Y;
			}
		}
		DB(zv, bv, xv, elapsed());
	}
	DB(cempty.S);
	DB(total_removed / total_removals);
	DB(total_removals);
	DB(zero_removals);
	
	DB(acc);
	DB(sol.S);
	DB(elapsed());
	cerr << "[DATA] gstep = " << greedy_step << endl;
	cerr << "[DATA] step = " << step << endl;
	cerr << "[DATA] maxv = " << xmaxv << endl;
	cout << sol.S << endl;
	for (int p : sol) cout << pstr(p) << endl;
	return 0;
}