// Author: Psyho
// Twitter: https://twitter.com/fakepsyho
// Site: http://psyho.gg/
 
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
template<typename T> ostream& operator<<(ostream &os, VC<T> v) {os << "{"; REP(i, v.size()) {if (i) os << ", "; os << v[i];} os << "}"; return os;}
template<typename T> ostream& operator<<(ostream &os, set<T> s) {VS vs(ALL(s)); return os << vs;}
template<typename T> string i2s(T x) {ostringstream o; o << x; return o.str();}
VS splt(string s, char c = ' ') {VS all; int p = 0, np; while (np = s.find(c, p), np >= 0) {if (np != p) all.PB(s.substr(p, np - p)); p = np + 1;} if (p < s.size()) all.PB(s.substr(p)); return all;}

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
 
bool USE_BLOCKS = false; 
bool FULL_SOLVE = true;
double TIME_LIMIT = 99.6;
 
static RNG rng;

double global_time = get_time();
double start_time = global_time;

const int MAX_N = 50;
const int MAX_C = 6;
const int MAX_D = 10;
const int MAX_S = 1500;

const int MAX_PATH = MAX_D*MAX_S*2;
const int MAX_PATH_LEN = MAX_N + MAX_N + 20;
const int MAX_SIZE = 64 * (MAX_N + 2);
const int MAX_MOD = MAX_PATH * 20;

const int MOD_MULT = 500;

int N, C, D, S;
VVS tc;
VVI marks;

int path_a[MAX_PATH];
int path_b[MAX_PATH];
int path_start[MAX_PATH];
int path_maxlen[MAX_PATH];
int path_dancer[MAX_D+1];
int path_total = 0;

int dist(int a, int b) {
	return abs((a&63)-(b&63)) + abs((a>>6)-(b>>6));
}

int dd[] = {-64,-1,+1,+64};
int ddv[] = {-64,+64};
int ddh[] = {-1,+1};
int dds[] = {+1,+64,-1,-64};


int mod_data[MAX_PATH * MOD_MULT];
int mod_pos[MAX_PATH * MOD_MULT];
int mod_current[MAX_PATH];
int mod_start[MAX_PATH];
int mod_no[MAX_PATH];

int target_color;
bool target_optimal = false;
VVI clayout;


short outside[MAX_SIZE];

struct State {
	short sum[MAX_SIZE];
	short colors[MAX_SIZE];
	int comp_count[MAX_C];
	
	int comp[MAX_SIZE];
	int comp_no = 0;
	short stack[MAX_N*MAX_N];
	short stack2[MAX_N*MAX_N];
	
	short path_data[MAX_PATH];
	int path_end[MAX_PATH];
	
	int path_copy[MAX_PATH_LEN];
	int path_copy_len;
	
	void init() {
		ZERO(sum);
		ZERO(rest_used);
		ZERO(mod_current);
		REP(i, path_total) path_end[i] = path_start[i];
		REP(i, N) {
			colors[i+1] = -1;
			colors[i+1+(N+1)*64] = -1;
			colors[(i+1)*64] = -1;
			colors[(i+1)*64+N+1] = -1;
			outside[i+1] = -1;
			outside[i+1+(N+1)*64] = -1;
			outside[(i+1)*64] = -1;
			outside[(i+1)*64+N+1] = -1;
		}
	}
	
	void change_sum(int p, int v, bool rem=true) {
		if (rem) rem_component(p);
		sum[p] += v;
		if (colors[p] != -1) colors[p] = tc[(p>>6)-1][(p&63)-1][sum[p] % C] - '0';
	}
	
	void calc_colors() {
		FOR(r, 1, N+1) FOR(c, 1, N+1)
			if (colors[r*64+c] != -1) colors[r*64+c] = tc[r-1][c-1][sum[r*64+c] % C] - '0';
	}
	
	void find_all_components() {
		if (FULL_SOLVE) return;
		FOR(r0, 1, N+1) FOR(c0, 1, N+1) {
			int p0 = r0*64+c0;
			if (comp[p0] || colors[p0] == -1) continue;
			comp_no++;
			comp_count[colors[p0]]++;
			int stack_pos = 0;
			stack[stack_pos++] = p0;
			comp[p0] = comp_no;
			while (stack_pos) {
				int p = stack[--stack_pos];
				if (colors[p] == colors[p-64] && comp[p-64] != comp_no) {if (comp[p-64]) rem_component(p-64, comp_no); else stack[stack_pos++] = p-64, comp[p-64] = comp_no;}
				if (colors[p] == colors[p- 1] && comp[p- 1] != comp_no) {if (comp[p- 1]) rem_component(p- 1, comp_no); else stack[stack_pos++] = p- 1, comp[p- 1] = comp_no;}
				if (colors[p] == colors[p+ 1] && comp[p+ 1] != comp_no) {if (comp[p+ 1]) rem_component(p+ 1, comp_no); else stack[stack_pos++] = p+ 1, comp[p+ 1] = comp_no;}
				if (colors[p] == colors[p+64] && comp[p+64] != comp_no) {if (comp[p+64]) rem_component(p+64, comp_no); else stack[stack_pos++] = p+64, comp[p+64] = comp_no;}
			}
		}
	}
	
	void find_components() {
		if (FULL_SOLVE) return;
		FOR(r0, 1, N+1) FOR(c0, 1, N+1) {
			int p0 = r0*64+c0;
			if (comp[p0] || colors[p0] == -1) continue;
			comp_no++;
			comp_count[colors[p0]]++;
			int stack_pos = 0;
			stack[stack_pos++] = p0;
			comp[p0] = comp_no;
			while (stack_pos) {
				int p = stack[--stack_pos];
				if (colors[p] == colors[p-64] && comp[p-64] != comp_no) {if (comp[p-64]) rem_component(p-64, comp_no); else stack[stack_pos++] = p-64, comp[p-64] = comp_no;}
				if (colors[p] == colors[p- 1] && comp[p- 1] != comp_no) {if (comp[p- 1]) rem_component(p- 1, comp_no); else stack[stack_pos++] = p- 1, comp[p- 1] = comp_no;}
				if (colors[p] == colors[p+ 1] && comp[p+ 1] != comp_no) {if (comp[p+ 1]) rem_component(p+ 1, comp_no); else stack[stack_pos++] = p+ 1, comp[p+ 1] = comp_no;}
				if (colors[p] == colors[p+64] && comp[p+64] != comp_no) {if (comp[p+64]) rem_component(p+64, comp_no); else stack[stack_pos++] = p+64, comp[p+64] = comp_no;}
			}
		}
	}
	
	void rem_component(int p0, int new_val=0) {
		if (FULL_SOLVE) return;
		if (comp[p0] == 0 || colors[p0] == -1) return;
		comp_count[colors[p0]]--;
		int comp_val = comp[p0];
		int stack_pos = 0;
		stack2[stack_pos++] = p0;
		comp[p0] = new_val;
		while (stack_pos) {
			int p = stack2[--stack_pos];
			if (comp[p-64] == comp_val) stack2[stack_pos++] = p-64, comp[p-64] = new_val;
			if (comp[p- 1] == comp_val) stack2[stack_pos++] = p- 1, comp[p- 1] = new_val;
			if (comp[p+ 1] == comp_val) stack2[stack_pos++] = p+ 1, comp[p+ 1] = new_val;
			if (comp[p+64] == comp_val) stack2[stack_pos++] = p+64, comp[p+64] = new_val;
		}
	}
	
	void recalc_components() {
		calc_colors();
		ZERO(comp_count);
		ZERO(comp);
		comp_no = 0;
		find_all_components();
	}
	
	int calc_score(double t) {
		if (FULL_SOLVE) {
			int total = 0;
			
			REP(i, clayout.size()) {
				bool bad = false;
				for (int x : clayout[i]) {
					if (colors[x] == target_color) {
						total++;
					} else {
						bad = true;
					}
				}
				if (bad) break;
			}
			
			return (N*N-total)*1;
		} else {
			int rv = 0;
			REP(i, C) rv += comp_count[i] * comp_count[i];
			return rv;
		}
	}
	
	int recalc_score() {
		recalc_components();
		return calc_score(1);
	}
	
	void copy_path(int id) {
		path_copy_len = 0;
		FOR(i, path_start[id], path_end[id]) path_copy[path_copy_len++] = path_data[i];
	}
	
	void remove_path(int id) {
		FOR(i, path_start[id]+1, path_end[id]) change_sum(path_data[i], -1);
		path_end[id] = path_start[id];
	}
	
	int rest_used[MAX_SIZE];
	short rest_colors[MAX_SIZE];
	int rest_pos[2*MAX_PATH_LEN];
	int rest_used_no = 0;
	int rest_pos_no = 0;
	
	int xolors[MAX_SIZE];
	
	void restore_path(int id) {
		rest_used_no++;
		rest_pos_no = 0;
		
		FOR(r0, 1, N+1) FOR(c0, 1, N+1) xolors[r0*64+c0] = colors[r0*64+c0];
		
		FOR(i, path_start[id]+1, path_end[id]) {
			int p = path_data[i];
			if (rest_used[p] != rest_used_no) {
				rest_used[p] = rest_used_no;
				rest_pos[rest_pos_no++] = p;
				rest_colors[p] = colors[p];
			}
			change_sum(p, -1, false);
		}
		path_end[id] = path_start[id];
		
		REP(i, path_copy_len) path_data[path_end[id]++] = path_copy[i];
		FOR(i, path_start[id]+1, path_end[id]) {
			int p = path_data[i];
			if (rest_used[p] != rest_used_no) {
				rest_used[p] = rest_used_no;
				rest_pos[rest_pos_no++] = p;
				rest_colors[p] = colors[p];
			}
			change_sum(p, +1, false);
		}
		
		REP(i, rest_pos_no) if (rest_colors[rest_pos[i]] != colors[rest_pos[i]])
			swap(rest_colors[rest_pos[i]], colors[rest_pos[i]]);
		REP(i, rest_pos_no) if (rest_colors[rest_pos[i]] != colors[rest_pos[i]])
			rem_component(rest_pos[i]);
		REP(i, rest_pos_no) if (rest_colors[rest_pos[i]] != colors[rest_pos[i]])
			swap(rest_colors[rest_pos[i]], colors[rest_pos[i]]);
	}
	
	void add_random_path(int id) {
		assert(path_start[id] == path_end[id]);
		int p = path_a[id];
		path_data[path_end[id]++] = p;
		
		int type = rng.next(5) == 0;
		int par = rng.next(1, 20);
		int par2 = rng.next(10);
		
		int path_left = path_maxlen[id] - 1;
		
		while (true) {
			
			int ad = dist(p, path_b[id]);
			if (ad == 0) break;
			
			int bp = -1;
			double bv = 1e9;
			
			int atype = type == 0 || type == 1 && rng.next(11) < par2;
			REP(d, 4) {
				int np = p + dd[d];
				int dl = dist(np, path_b[id]) + 1;
				if (dl > path_left || outside[np] == -1) continue;
				double av = 0;
				if (atype) {
					av = dl * 1000.0;
					av += rng.next(par) ? min(min(np/64-1,N-np/64),min(np%64-1,N-np%64)) : 0;
				}
				av += rng.next_double();
				if (av < bv) {
					bv = av;
					bp = np;
				}
			}
			assert(bp != -1);
			path_left--;
			path_data[path_end[id]++] = bp;
			p = bp;
		}
		FOR(i, path_start[id]+1, path_end[id]) change_sum(path_data[i], +1);
	}
	
	INLINE void mod_clear(int id) {
		mod_no[id] = 0;
		mod_pos[mod_start[id]]   = mod_start[id];
		mod_pos[mod_start[id]+1] = mod_start[id];
	}
	
	INLINE void mod_new(int id) {
	}
	
	INLINE void mod_add(int id, int v, int p) {
		mod_data[mod_pos[mod_start[id]+mod_no[id]+1]++] = v;
		mod_data[mod_pos[mod_start[id]+mod_no[id]+1]++] = p;
	}
	
	INLINE void mod_done(int id) {
		mod_no[id]++;
		mod_pos[mod_start[id]+mod_no[id]+1] = mod_pos[mod_start[id]+mod_no[id]];
	}
	
	
	int modify_path(int id) {
		assert(path_start[id] != path_end[id]);
		
		if (!mod_current[id]) {
			mod_clear(id);
			
			FOR(i, path_start[id]+2, path_end[id]) {
				int a = path_data[i-2];
				int b = path_data[i];
				if ((a&63)!=(b&63) && (a>>6)!=(b>>6)) {
					mod_new(id);
					mod_add(id, 0, i-1);
					mod_add(id, (a&63)==(path_data[i-1]&63) ? (b&63)+(a&~63) : (a&63)+(b&~63), i-1);
					mod_done(id);
				}
			}
			
			FOR(repeat, 1, 4) {
				if (C > repeat && path_end[id] - path_start[id] + repeat*2 <= path_maxlen[id]) {
					FOR(i, path_start[id], path_end[id]) {
						REP(d, 4) {
							int np = path_data[i] + dd[d];
							if (outside[np] == -1) continue;
							mod_new(id);
							REP(_, repeat) {
								mod_add(id, np, i);
								mod_add(id, path_data[i], i);
							}
							mod_done(id);
						}
					}
				}
			}		
			
			FOR(repeat, 1, 4) {
				if (C > repeat && path_end[id] - path_start[id] + repeat*2 <= path_maxlen[id]) {
					FOR(i, path_start[id]+1, path_end[id]) {
						int dp = path_data[i] - path_data[i-1];
						int *pdd = abs(dp) == 1 ? ddv : ddh;
						REP(d, 2) {
							int np = path_data[i-1] + pdd[d];
							if (outside[np] == -1) continue;
							mod_new(id);
							REP(_, repeat) {
								mod_add(id, np+dp, i);
								mod_add(id, np, i);
							}
							mod_done(id);
						}
					}
				}
			}
			
			FOR(repeat, 1, 4) {
				FOR(i, path_start[id]+1+repeat*2, path_end[id]) {
					if (dist(path_data[i-1-repeat*2], path_data[i]) == 1) {
						mod_new(id);
						REP(j, repeat) {
							mod_add(id, 0, i-j*2-1);
							mod_add(id, 0, i-j*2-2);
						}
						mod_done(id);
					}
				}
			}
			
			mod_current[id] = 1;
		}
		
		if (mod_no[id] == 0) return -1;
		
		int mod_id = mod_start[id] + rng.next(mod_no[id]);
		
		for (int i = mod_pos[mod_id]; i < mod_pos[mod_id+1]; i += 2) {
			int p = mod_data[i];
			int pos = mod_data[i+1];
			
			if (p == 0) {
				change_sum(path_data[pos], -1);
				if (!FULL_SOLVE) {
					FOR(j, pos, path_end[id])
						path_data[j] = path_data[j+1];
					path_end[id]--;
				}
			} else {
				change_sum(p, +1);
				if (!FULL_SOLVE) {
					for (int j = path_end[id]; j > pos; j--) 
						path_data[j] = path_data[j-1];
					path_data[pos] = p;
					path_end[id]++;
				}
			}
		}
		return mod_id;
	}
	
	VS create_answer() {
		VS dancers;
		REP(i, D) {
			string s(S, '-');
			int sp = 0;
			FOR(id, path_dancer[i], path_dancer[i+1]) {
				FOR(j, path_start[id]+1, path_end[id]) {
					int diff = path_data[j] - path_data[j-1];
					if (diff == -64) s[sp] = 'U';
					if (diff == - 1) s[sp] = 'L';
					if (diff == + 1) s[sp] = 'R';
					if (diff == +64) s[sp] = 'D';
					sp++;
				}
				sp += path_maxlen[id] - (path_end[id] - path_start[id]);
			}
			dancers.PB(s);
		}
		
		VS rv(S, string(D, ' '));
		REP(i, S) REP(j, D) rv[i][j] = dancers[j][i];
		return rv;
	}
};

State best;

VS solve(int N, int C, int D, int S, VVS &tile_colors, VVI &marks) {
	::N = N;
	::C = C;
	::D = D;
	::S = S;
	::tc = tile_colors;
	
	DB(N, C, D, S);
	
	DB(get_time() - start_time);
	
	// Parse paths
	assert(marks.size() == D);
	path_start[0] = 0;
	int path_sum = 0;
	REP(i, D) {
		path_dancer[i] = path_total;
		for (int j = 0; j + 3 < marks[i].size(); j += 3) {
			path_a[path_total] = (marks[i][j+0]+1) + (marks[i][j+1]+1)*64;
			path_b[path_total] = (marks[i][j+3]+1) + (marks[i][j+4]+1)*64;
			path_sum += dist(path_a[path_total], path_b[path_total]);
			path_maxlen[path_total] = marks[i][j+5] - marks[i][j+2] + 1;
			path_total++;
			path_start[path_total] = path_start[path_total-1] + path_maxlen[path_total-1];
			
			mod_start[path_total] = mod_start[path_total-1] + path_maxlen[path_total-1] * MOD_MULT;
			mod_pos[mod_start[path_total]] = mod_start[path_total];
		}
	}
	path_dancer[D] = path_total;
	
	
	if (C % 2) {
		target_color = 0;
		target_optimal = true;
	} else {
		REP(i, C) {
			int sum = 0;
			REP(r, N) REP(c, N) {
				REP(j, C) if (tc[r][c][j] - '0' == i) {
					sum += j;
					break;
				}
			}
			sum += path_sum;
			if (sum % 2 == 0) {
				target_color = i;
				target_optimal = true;
			}
		}
	}
	
	REP(r, N) REP(c, N) {
		int d = min(min(r,N-1-r), min(c,N-1-c)) * 2 + (min(r,N-1-r) != min(c,N-1-c));
		while (clayout.size() <= d) clayout.PB(VI());
		clayout[d].PB((r+1)*64+(c+1));
	}
	while (clayout.size() && clayout.back().size() <= 36) {
		clayout[clayout.size() - 2].insert(clayout[clayout.size() - 2].begin(), ALL(clayout.back()));
		clayout.pop_back();
	}
	
	// Hill Climbing
	State s;
	s.init();
	if (USE_BLOCKS) FOR(r, 2, N) FOR(c, 2, N) s.colors[r*64+c] = -1;
	REP(i, path_total) s.add_random_path(i);
	double bv = s.recalc_score();
	bv = s.calc_score(0.0);
	DB(bv);
	VI vv(s.comp_count, s.comp_count + C);
	DB(vv);
	
	double ov = bv;
	double xv = bv;
	
	
	int step = 0;
	double time_passed = 0;
	double cur_time_passed = 0;
	
	VVI acc(300, VI(2));
	
	int block_step = USE_BLOCKS ? 1 : 0;	
	VI steps(300);
	
	int block_best = 1e9;
	
	int last_imp = 0;
	
	REP(pass, 2) {
		while (true) {
			if ((step & 2047) == 0) {
				time_passed = (get_time() - start_time) / TIME_LIMIT;
				if (time_passed >= 1.0) break;
				
				if (pass == 0 && time_passed > 0.2 && xv > 2) {
					xv = 1e9;
					FULL_SOLVE = false;
					USE_BLOCKS = true;
					block_step = 1;
					if (USE_BLOCKS) FOR(r, 2, N) FOR(c, 2, N) s.colors[r*64+c] = -1;
					bv = s.recalc_score();
					bv = s.calc_score(0.0);
					ov = bv;
					xv = bv;
					DB(xv);
					TIME_LIMIT -= get_time() - start_time;
					start_time = get_time();
					break;
				}
				
				if (USE_BLOCKS) {
					if (block_step && (time_passed > 0.015 * block_step || xv == 1)) {
						DB(block_step, bv, xv, block_best);
						block_best = 1e9;
						
						REP(i, N) {
							s.colors[(block_step+1)+(i+1)*64] = 0;
							s.colors[(N-block_step)+(i+1)*64] = 0;
							s.colors[(block_step+1)*64+(i+1)] = 0;
							s.colors[(N-block_step)*64+(i+1)] = 0;
						}
						block_step++;
						if ((N-block_step) - (block_step+1) <= 8) {
							block_step = 0;
							FOR(r, 1, N) FOR(c, 1, N) s.colors[r*64+c] = 0;
						}
						bv = s.recalc_score();
						xv = 1e9;
						REP(i, path_total) mod_current[i] = 0;
					}
				}
				if (block_step) {
					cur_time_passed = (time_passed - 0.015 * (block_step - 1)) / 0.015; 
				} else {
					// const int TIME_STEPS = 10;
					// int time_step = int(time_passed * TIME_STEPS);
					// cur_time_passed = 1.0 * time_step / TIME_STEPS + (1 - 1.0 * time_step / TIME_STEPS) * (time_passed * TIME_STEPS - int(time_passed * TIME_STEPS));
					cur_time_passed = time_passed;
				}
				
				if (FULL_SOLVE && pass == 0 && xv == 1 && step > last_imp + 1000000) {
					DB("restart");
					last_imp = step;
					s.init();
					REP(i, path_total) s.add_random_path(i);
					bv = s.recalc_score();
					bv = s.calc_score(0.0);
					target_color = rng.next(C);
				}
				
			}
			step++;
			steps[block_step]++;
			
			int id; 
			
			int type = rng.next(FULL_SOLVE || block_step ? 4 : 20) == 0;
			int mod_id = -1;
			
			if (type == 1) {
				id = rng.next(path_total);
				s.copy_path(id);
				s.remove_path(id);
				s.add_random_path(id);
			} else if (type == 0) {
				while (mod_id == -1) {
					id = rng.next(path_total);
					if (!FULL_SOLVE) s.copy_path(id);
					mod_id = s.modify_path(id);
				}
			}
			
			s.find_components();
			double av = s.calc_score(time_passed);
			bool kick = time_passed < 0.9 && rng.next(block_step ? 1000 : time_passed * 25000 + 1000) == 0;
			if (FULL_SOLVE && bv == 1) kick = true;
				
			bool accept;
			if (FULL_SOLVE) {
				// accept = av <= bv || kick && av <= bv * 1.025 + 20;
				accept = av <= bv || kick && av <= bv * 1.050 + 50;
			} else {
			// if (av <= bv || kick && av <= bv * 1.05 + 50 || rng.next_double() <= exp((bv - av) / max(1.0, sqrt(bv)) / (0.500 * (1 - time_passed) * (1 - time_passed)))) {
			// if (av <= bv || kick && av <= bv * 1.05 + 50 || block_step == 0 && rng.next_double() <= exp((xv - av) / max(10.0, bv) / (0.050 * (1 - cur_time_passed) * (1 - cur_time_passed)))) {
				accept = av <= bv || kick && av <= bv * 1.025 + 50 || block_step == 0 && rng.next_double() <= exp((bv - av) / max(10.0, pow(bv, 0.8)) / (0.050 * (1 - cur_time_passed) * (1 - cur_time_passed)));
			}
			
			if (accept) {
				if (FULL_SOLVE && mod_id != -1) {
					for (int i = mod_pos[mod_id]; i < mod_pos[mod_id+1]; i += 2) {
						int p = mod_data[i];
						int pos = mod_data[i+1];
						
						if (p == 0) {
							FOR(j, pos, s.path_end[id])
								s.path_data[j] = s.path_data[j+1];
							s.path_end[id]--;
						} else {
							for (int j = s.path_end[id]; j > pos; j--) 
								s.path_data[j] = s.path_data[j-1];
							s.path_data[pos] = p;
							s.path_end[id]++;
						}
					}
				}
				
				if (block_step == 0 && av < xv) {
					last_imp = step;
					DB(step, av, get_time() - start_time);
					xv = av;
					memcpy(best.path_data, s.path_data, 4*path_start[path_total]);
					memcpy(best.path_end,  s.path_end,  4*path_total);
					if (av == 0 || av == 1 && !target_optimal) goto done;
					
					if (av == 1) {
						FOR(r, 1, N) FOR(c, 1, N) {
							int p = r*64+c;
							if (s.colors[p] != target_color) {
								DB(target_color, s.colors[p]);
								REP(i, C) DB(tc[r-1][c-1][i]);
							}
						}
					}
				}
				
				if (block_step && av < block_best) block_best = av;
				bv = av;
				acc[block_step][type]++;
				mod_current[id] = 0;
				
			} else {
				if (!FULL_SOLVE || mod_id == -1) {
					s.restore_path(id);
				} else {
					for (int i = mod_pos[mod_id]; i < mod_pos[mod_id+1]; i += 2) {
						int p = mod_data[i];
						int pos = mod_data[i+1];
						
						if (p == 0) {
							s.change_sum(s.path_data[pos], +1);
						} else {
							s.change_sum(p, -1);
						}
					}
				}
			}
			
		}
	}
	done:
	
	DB(ov, xv, bv, step, S*D, N*N, C, get_time() - global_time);
	DB(acc[1][0], acc[1][1], steps[1]);
	DB(acc[0][0], acc[0][1], steps[0]);
	return best.create_answer();
}


int main() {
    int N, C, D, S;
    cin >> N >> C >> D >> S;

    VVS tile_colors(N, VS(N));
	REP(r, N) REP(c, N) cin >> tile_colors[r][c];
	
    VVI marks(D);
	REP(i, D) {
        int numMarks;
        cin >> numMarks;
        marks[i].resize(3 * numMarks);
		REP(j, 3*numMarks) 
            cin >> marks[i][j];
    }

    VS rv = solve(N, C, D, S, tile_colors, marks);

    cout << rv.size() << endl;
	for (string &s : rv) cout << s << endl;
    cout.flush();

    return 0;
}
