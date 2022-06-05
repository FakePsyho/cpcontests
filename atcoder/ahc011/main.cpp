// Author: Psyho
// Twitter: https://twitter.com/fakepsyho

//[Template]
 
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

//[Stats]
double timer_bs = 0;
double timer_wm = 0;
double timer_heur = 0;
double timer_hc = 0;
int stat_matching_repeat = 0;
int stat_grids = 0;
int stat_bs_runs = 0;
int stat_hc_steps = 0;

//[Params]
const double TIME_LIMIT_MUL = 1;
const double TIME_LIMIT = 2.9 * TIME_LIMIT_MUL;
int WIDTH = 1000;
const int BASE_KICK_STEPS = 100;

const int PRESOLVE_N = 8;

const int BS_EARLYEXIT = 1;

const int ONLY_GRIDS = 0;

const int BS_BENCHMARK = 0;
const int BSB_N = 10;
const int BSB_RUNS = 20;
const int BSB_WIDTH = 1000;

const int HC_BENCHMARK = 0;

const int FORCE_FAIL = 0;
	

//[Code]
int dr[] = {-1,1,0,0};
int dc[] = {0,0,-1,1};
char dn[] = {'U','D','L','R'};

double cur_time_limit;
const int MAX_N = 10;

int N, T;

int ogrid[MAX_N][MAX_N];
VVI grid;

int vs[MAX_N][MAX_N];

const double BAD_CONN_VALUE = -2.0;
const double OUTER_CONN_VALUE = -2.0;

const double EMPTY_PENALTY = 100;

double MAX_SCORE;
bool grid_fully_solved = false;
bool bad_empty[MAX_N*MAX_N];
bool sg_l[MAX_N*MAX_N];
bool sg_r[MAX_N*MAX_N];
bool sg_d[MAX_N*MAX_N];
bool sg_u[MAX_N*MAX_N];
double score_grid(VI &g, VI &cn, bool show = false) {
	double score = 0;
	
	int sg_vs[MAX_N*MAX_N];
	int stack[MAX_N*MAX_N];
	int stack_pos;
	int biggest = 0;
	ZERO(sg_vs);
	REP(p0, N*N) {
		if (sg_vs[p0]) continue;
		
		if (g[p0] == 0) {
			if (show) DB(p0);
			if (bad_empty[p0]) score -= EMPTY_PENALTY;
			continue;
		}
		
		int cc = 0;
		stack[0] = p0;
		stack_pos = 1;
		sg_vs[p0] = 1;
		while (stack_pos) {
			int p = stack[--stack_pos];
			cc++;
			if ((cn[p] & 1) && !sg_vs[p-1]) sg_vs[p-1] = 1, stack[stack_pos++] = p-1;
			if ((cn[p] & 2) && !sg_vs[p-N]) sg_vs[p-N] = 1, stack[stack_pos++] = p-N;
			if ((cn[p] & 4) && !sg_vs[p+1]) sg_vs[p+1] = 1, stack[stack_pos++] = p+1;
			if ((cn[p] & 8) && !sg_vs[p+N]) sg_vs[p+N] = 1, stack[stack_pos++] = p+N;
		}
		
		if (show) DB(cc);
		// biggest = max(biggest, cc);
		score += cc * log(cc+1);
	}
	// score += biggest * log(biggest+1);
	grid_fully_solved = score > MAX_SCORE;
	
	return score;
}

INLINE void update_connection(VI &g, VI &cn, int p) {
	cn[p] = 0;
	if ((g[p] & 1) && sg_l[p] && (g[p-1] & 4)) cn[p] |= 1;
	if ((g[p] & 2) && sg_u[p] && (g[p-N] & 8)) cn[p] |= 2;
	if ((g[p] & 4) && sg_r[p] && (g[p+1] & 1)) cn[p] |= 4;
	if ((g[p] & 8) && sg_d[p] && (g[p+N] & 2)) cn[p] |= 8;
}

INLINE void update_connection_neighbors(VI &g, VI &cn, int p) {
	cn[p] = 0;
	if (sg_l[p]) { if ((g[p] & 1) && (g[p-1] & 4)) cn[p] |= 1, cn[p-1] |= 4; else cn[p-1] &= ~4; }
	if (sg_u[p]) { if ((g[p] & 2) && (g[p-N] & 8)) cn[p] |= 2, cn[p-N] |= 8; else cn[p-N] &= ~8; }
	if (sg_r[p]) { if ((g[p] & 4) && (g[p+1] & 1)) cn[p] |= 4, cn[p+1] |= 1; else cn[p+1] &= ~1; }
	if (sg_d[p]) { if ((g[p] & 8) && (g[p+N] & 2)) cn[p] |= 8, cn[p+N] |= 2; else cn[p+N] &= ~2; }
}

void show_grid(VVI &g) {
	// VS output(N*3, string(N*3));
	VS output(N*3, string(N*3, ' '));
	REP(r, N) REP(c, N) {
		if (g[r][c]) output[r*3+1][c*3+1] = '*';
		if (g[r][c] & 1) output[r*3+1][c*3+0] = '-';
		if (g[r][c] & 2) output[r*3+0][c*3+1] = '|';
		if (g[r][c] & 4) output[r*3+1][c*3+2] = '-';
		if (g[r][c] & 8) output[r*3+2][c*3+1] = '|';
	}
	for (string s : output) cerr << s << endl;
	cerr << endl;
	
}

void init_grid() {
	grid = VVI(N, VI(N));
	REP(r, N) REP(c, N) grid[r][c] = ogrid[r][c];
}


bool solvable(VVI &g, int ev) {
	int erow = -1;
	VI v;
	REP(r, N) REP(c, N) {
		if (g[r][c] != ev) {
			v.PB(g[r][c]);
		} else {
			erow = r;
		}
	}
	int inv = 0;
	REP(i, v.S) REP(j, i) inv += v[j] > v[i];
	if (N % 2 == 0) inv += erow - (ev / N);
	return inv % 2 == 0;
}

void exec_moves(VVI &g, string moves) {
	PII p = {-1, -1}; REP(r, N) REP(c, N) if (g[r][c] == 0) p = MP(r, c);
	for (char &c : moves) {
		if (c == 'U') {
			swap(g[p.X][p.Y], g[p.X-1][p.Y]);
			p.X--;
		} else if (c == 'D') {
			swap(g[p.X][p.Y], g[p.X+1][p.Y]);
			p.X++;
		} else if (c == 'L') {
			swap(g[p.X][p.Y], g[p.X][p.Y-1]);
			p.Y--;
		} else if (c == 'R') {
			swap(g[p.X][p.Y], g[p.X][p.Y+1]);
			p.Y++;
		} else {
			assert(false);
		}
	}
}

string solve_heuristic(VVI grid, VVI bgrid, int n_lines = N) {
	VVI used(N, VI(N));
	
	auto find_cell = [&](VVI &g, int cell, int rt = -1, int ct = -1) -> PII {
		PII bp = MP(-1, -1);
		double bv = 1e9;
		REP(r, N) REP(c, N) if (!used[r][c] && g[r][c] == cell) {
			int a = abs(r-rt);
			int b = abs(c-ct);
			double av = min(a,b)*6 + (max(a,b)-min(a,b))*5 + rng.next_double()*5;
			if (av < bv) {
				bv = av;
				bp = MP(r, c);
			}
		}
		assert(bp.X != -1);
		return bp;
	};
	
	auto move = [](int r0, int c0, int r1, int c1) -> char {
		return r0 < r1 ? 'D' : r0 > r1 ? 'U' : c0 < c1 ? 'R' : 'L';
	};
	
	auto find_move = [&](VVI &g, int xr, int xc) -> string {
		PII spos = find_cell(g, 0);
		ZERO(vs);			
		queue<int> q;
		q.push(spos.X), q.push(spos.Y);
		vs[spos.X][spos.Y] = -1;
		while (!q.empty()) {
			int r = q.front(); q.pop();
			int c = q.front(); q.pop();
			if (xr == r && xc == c) break;
			REP(d, 4) {
				int nr = r + dr[d];
				int nc = c + dc[d];
				if (nr < 0 || nr >= N || nc < 0 || nc >= N || vs[nr][nc] || used[nr][nc]) continue;
				vs[nr][nc] = r * N + c + 1;
				q.push(nr);
				q.push(nc);
			}
		}
		
		if (!vs[xr][xc]) {
			DB(spos);
			DB(xr, xc);
			REP(r, N) {
				REP(c, N) cerr << used[r][c];
				cerr << endl;
			}
		}
		assert(vs[xr][xc]);
		string rv = "";
		int r = xr;
		int c = xc;
		while (vs[r][c] != -1) {
			int nr = (vs[r][c] - 1) / N;
			int nc = (vs[r][c] - 1) % N;
			rv += move(nr, nc, r, c);
			r = nr;
			c = nc;
		}
		
		reverse(ALL(rv));
		return rv;
	};
	
	auto exec_move = [&](VVI &g, string moves) -> void {
		PII p = find_cell(g, 0);
		for (char &c : moves) {
			if (c == 'U') {
				swap(g[p.X][p.Y], g[p.X-1][p.Y]);
				p.X--;
			} else if (c == 'D') {
				swap(g[p.X][p.Y], g[p.X+1][p.Y]);
				p.X++;
			} else if (c == 'L') {
				swap(g[p.X][p.Y], g[p.X][p.Y-1]);
				p.Y--;
			} else if (c == 'R') {
				swap(g[p.X][p.Y], g[p.X][p.Y+1]);
				p.Y++;
			} else {
				assert(false);
			}
		}
	};
	
	auto move_empty = [&](VVI &g, int rd, int cd) -> string {
		string m = find_move(g, rd, cd);
		exec_move(g, m);
		return m;
	};
	
	auto move_cell_line = [&](VVI &g, int r0, int c0, int rd, int cd) -> string {
		assert(r0==rd || c0==cd);
		string rv = "";
		while (r0!=rd || c0!=cd) {
			used[r0][c0] = 1;
			int rt = r0 + (r0 < rd ? 1 : r0 > rd ? -1 : 0);
			int ct = c0 + (c0 < cd ? 1 : c0 > cd ? -1 : 0);
			string m = find_move(grid, rt, ct);
			used[r0][c0] = 0;
			m += move(rt, ct, r0, c0);
			exec_move(g, m);
			r0 = rt;
			c0 = ct;
			rv += m;
		}
		return rv;
	};
	
	int mc_vs[MAX_N][MAX_N][MAX_N][MAX_N];
	ZERO(mc_vs);
	
	int mc_vs_no = 0;
	auto move_cell = [&](VVI &g, int r0, int c0, int rd, int cd) -> string {
		// DB(r0, c0, rd, cd);
		PII ep = find_cell(g, 0);
		int correct_cell = g[r0][c0];
		mc_vs_no += 10;
		if (mc_vs_no > 2e9) {
			ZERO(mc_vs);
			mc_vs_no = 10;
		}
		#define PX pair<PII,PII>
		queue<PX> q;
		q.push(MP(ep, MP(r0,c0)));
		mc_vs[ep.X][ep.Y][r0][c0] = mc_vs_no;
		
		ep.X = -1;
		int steps = 0;
		while (!q.empty()) {
			PX p = q.front(); q.pop();
			steps++;
			if (p.Y.X == rd && p.Y.Y == cd) {
				ep = p.X;
				break;
			}
			REP(d, 4) {
				int nr = p.X.X + dr[d];
				int nc = p.X.Y + dc[d];
				int br = p.Y.X;
				int bc = p.Y.Y;
				if (nr < 0 || nc < 0 || nr >= N || nc >= N || used[nr][nc]) continue;
				if (nr == br && nc == bc) {
					br -= dr[d];
					bc -= dc[d];
				}
				if (mc_vs[nr][nc][br][bc] < mc_vs_no) {
					mc_vs[nr][nc][br][bc] = mc_vs_no + 1 + d;
					q.push(MP(MP(nr,nc),MP(br,bc)));
				}
			}
		}
		assert(ep.X >= 0);
		
		PII bp = MP(rd, cd);
		string rv = "";
		while (true) {
			int v = mc_vs[ep.X][ep.Y][bp.X][bp.Y];
			assert(v >= mc_vs_no);
			int d = v - mc_vs_no - 1;
			if (d == -1) break;
			rv += dn[d];
			ep.X -= dr[d];
			ep.Y -= dc[d];
			if (ep == bp) {
				bp.X += dr[d];
				bp.Y += dc[d];
			}
		}
		reverse(ALL(rv));
		exec_move(g, rv);
		assert(g[rd][cd] == correct_cell);
		return rv;
	};
	
	auto solve3x3 = [&](VVI &bg, VVI &g, int r0, int c0) -> string {
		cerr << "Solve 3x3" << endl;
		
		VVI mg(3, VI(3));
		VVI ug(3, VI(3));
		
		int ev = 0;
		REP(r, 3) REP(c, 3) {
			int cell = bg[r0+r][c0+c];
			if (cell == 0) ev = r * 3 + c;
			REP(i, 3) REP(j, 3) if (!ug[i][j] && g[r0+i][c0+j] == cell) {
				ug[i][j] = 1;
				mg[i][j] = r * 3 + c;
				goto next;
			}
			next: ;
		}
		
		VI goal(9); REP(i, 9) goal[i] = i;
		VI start(9);  REP(r, 3) REP(c, 3) start[r*3+c] = mg[r][c];
		
		int inv = 0;
		REP(i, 9) REP(j, i) if (start[i] != ev && start[j] != ev) inv += start[j] > start[i];
		if (inv & 1) return "X";
		
		map<VI, char> mp;
		queue<VI> q;
		q.push(start);
		
		mp[start] = 'X';
		
		int steps = 0;
		while (!q.empty()) {
			steps++;
			VI v = q.front(); q.pop();
			if (v == goal) break;

			int p = -1;
			REP(i, 9) if (v[i] == ev) p = i;
			assert(p != -1);
			
			if (p >= 3) {
				swap(v[p], v[p-3]);
				if (!mp.count(v)) {
					mp[v] = 'U', q.push(v);
					if (v == goal) break;
				}
				swap(v[p], v[p-3]);
			}
			if (p % 3 > 0) {
				swap(v[p], v[p-1]);
				if (!mp.count(v)) {
					mp[v] = 'L', q.push(v);
					if (v == goal) break;
				}
				swap(v[p], v[p-1]);
			}
			if (p < 6) {
				swap(v[p], v[p+3]);
				if (!mp.count(v)) {
					mp[v] = 'D', q.push(v);
					if (v == goal) break;
				}
				swap(v[p], v[p+3]);
			}
			if (p % 3 < 2) {
				swap(v[p], v[p+1]);
				if (!mp.count(v)) {
					mp[v] = 'R', q.push(v);
					if (v == goal) break;
				}
				swap(v[p], v[p+1]);
			}
		}				
		
		VI v = goal;
		if (!mp.count(v)) {
			DB(mp.S);
			return "X";
		}
		
		string rv = "";
		while (mp[v] != 'X') {
			char d = mp[v];
			rv += d;
			int p = -1; REP(i, 9) if (v[i] == ev) p = i;
			assert(p != -1);
			
			if (d == 'U') {
				swap(v[p], v[p+3]);
			} else if (d == 'D') {
				swap(v[p], v[p-3]);
			} else if (d == 'L') {
				swap(v[p], v[p+1]);
			} else if (d == 'R') {
				swap(v[p], v[p-1]);
			}
		}
		DB(mp.S);
		DB(rv.S);
		
		reverse(ALL(rv));
		return rv;
	};
	
	string sol = "";
	
	int i = 0;
	while (true) {
		if (i*2 >= min(n_lines, N-3)) break;
		
		PII p0, p1;
		
		for (int c=N-1-i; c>=i+2; c--) {
			PII p = find_cell(grid, bgrid[i][c]);
			sol += move_cell(grid, p.X, p.Y, i, c);
			used[i][c] = 1;
		}
		p0 = find_cell(grid, bgrid[i][i]);
		sol += move_cell(grid, p0.X, p0.Y, i, i+1);
		used[i][i+1] = 1;
		p1 = find_cell(grid, bgrid[i][i+1]);
		if (grid[i][i] == 0 && p1.X == i+1 && p1.Y == i) return "";
		if (p1.X == i && p1.Y == i) return "";
		sol += move_cell(grid, p1.X, p1.Y, i+1, i+1);
		used[i+1][i+1] = 1;
		sol += move_empty(grid, i, i);
		exec_move(grid, "RD");
		sol += "RD";
		used[i+1][i+1] = 0;
		used[i][i] = 1;
		
		FOR(r, i+1, N-2-i) {
			PII p = find_cell(grid, bgrid[r][i]);
			sol += move_cell(grid, p.X, p.Y, r, i);
			used[r][i] = 1;
		}
		p0 = find_cell(grid, bgrid[N-1-i][i]);
		sol += move_cell(grid, p0.X, p0.Y, N-2-i, i);
		used[N-2-i][i] = 1;
		p1 = find_cell(grid, bgrid[N-2-i][i]);
		if (grid[N-1-i][i] == 0 && p1.X == N-1-i && p1.Y == i+1) return "";
		if (p1.X == N-1-i && p1.Y == i) return "";
		sol += move_cell(grid, p1.X, p1.Y, N-2-i, i+1);
		used[N-2-i][i+1] = 1;
		sol += move_empty(grid, N-1-i, i);
		exec_move(grid, "UR");
		sol += "UR";
		used[N-2-i][i+1] = 0;
		used[N-1-i][i] = 1;
		
		if (i*2+1 >= min(n_lines, N-3)) break;
		FOR(c, i+1, N-2-i) {
			PII p = find_cell(grid, bgrid[N-1-i][c]);
			sol += move_cell(grid, p.X, p.Y, N-1-i, c);
			used[N-1-i][c] = 1;
		}
		p0 = find_cell(grid, bgrid[N-1-i][N-1-i]);
		sol += move_cell(grid, p0.X, p0.Y, N-1-i, N-2-i);
		used[N-1-i][N-2-i] = 1;
		p1 = find_cell(grid, bgrid[N-1-i][N-2-i]);
		if (grid[N-1-i][N-1-i] == 0 && p1.X == N-2-i && p1.Y == N-1-i) return "";
		if (p1.X == N-1-i && p1.Y == N-1-i) return "";
		sol += move_cell(grid, p1.X, p1.Y, N-2-i, N-2-i);
		used[N-2-i][N-2-i] = 1;
		sol += move_empty(grid, N-1-i, N-1-i);
		exec_move(grid, "LU");
		sol += "LU";
		used[N-2-i][N-2-i] = 0;
		used[N-1-i][N-1-i] = 1;
		
		for (int r = N-2-i; r >= i+3; r--) {
			PII p = find_cell(grid, bgrid[r][N-1-i]);
			sol += move_cell(grid, p.X, p.Y, r, N-1-i);
			used[r][N-1-i] = 1;
		}
		p0 = find_cell(grid, bgrid[i+1][N-1-i]);
		sol += move_cell(grid, p0.X, p0.Y, i+2, N-1-i);
		used[i+2][N-1-i] = 1;
		p1 = find_cell(grid, bgrid[i+2][N-1-i]);
		if (grid[i+1][N-1-i] == 0 && p1.X == i+1 && p1.Y == N-2-i) return "";
		if (p1.X == i+1 && p1.Y == N-1-i) return "";
		sol += move_cell(grid, p1.X, p1.Y, i+2, N-2-i);
		used[i+2][N-2-i] = 1;
		sol += move_empty(grid, i+1, N-1-i);
		exec_move(grid, "DL");
		sol += "DL";
		used[i+2][N-2-i] = 0;
		used[i+1][N-1-i] = 1;
		i++;
	}
	
	if (n_lines >= N) {
		string ss = solve3x3(bgrid, grid, N/2-1, N/2-1);
		if (ss == "X") return "";
		sol += ss;
	}
	
	return sol;
}

VVI matching_greedy(VVI &g, VVI &t, double rscale=1.0) {
	VVI rv(N, VI(N));
	VVI used(N, VI(N));
	int total_distance = 0;
	VPII order;
	REP(r, N) REP(c, N) {
		int d0 = min(r, N-1-r);
		int d1 = min(c, N-1-c);
		order.PB(MP((min(d0,d1)*2-(d0==d1)) * (N == 10 ? 0 : 1000) + r * N + c, r * N + c));
	}
	sort(ALL(order));
	for (PII &p : order) {
		int r = p.Y / N;
		int c = p.Y % N;
		double bv = 1e30;
		PII bp = MP(-1, -1);
		REP(i, N) REP(j, N) if (!used[i][j] && g[i][j] == t[r][c]) {
			double av = abs(r-i) + abs(c-j) + rng.next_double() * rscale;
			// double av = (r-i)*(r-i) + (c-j)*(c-j) + rng.next_double() * rscale;
			if (av < bv) {
				bv = av;
				bp = MP(i, j);
			}
		}
		assert(bp.X != -1);
		rv[bp.X][bp.Y] = r * N + c;
		used[bp.X][bp.Y] = 1;
		
		total_distance += abs(bp.X - r) + abs(bp.Y - c);
	}
	DB(total_distance);
	return rv;
}

struct Edge {
    int from;
	int to; 
	int cap;
	int cost;
};

struct MCMF {
	VVI adj, cost, capacity;
	VC<Edge> edges;

	const int INF = 1e9;

	void shortest_paths(int n, int v0, vector<int>& d, vector<int>& p) {
		d.assign(n, INF);
		d[v0] = 0;
		VB inq(n, false);
		queue<int> q;
		q.push(v0);
		p.assign(n, -1);

		while (!q.empty()) {
			int u = q.front();
			q.pop();
			inq[u] = false;
			for (int v : adj[u]) {
				if (capacity[u][v] > 0 && d[v] > d[u] + cost[u][v]) {
					d[v] = d[u] + cost[u][v];
					p[v] = u;
					if (!inq[v]) {
						inq[v] = true;
						q.push(v);
					}
				}
			}
		}
	}
	
	void add_edge(int from, int to, int cap, int cost) {
		Edge e = {from, to, cap, cost};
		edges.PB(e);
	}

	int run(int N, int K, int s, int t) {
		adj.assign(N, vector<int>());
		cost.assign(N, vector<int>(N, 0));
		capacity.assign(N, vector<int>(N, 0));
		
		random_shuffle(ALL(edges));
		for (Edge e : edges) {
			adj[e.from].push_back(e.to);
			adj[e.to].push_back(e.from);
			cost[e.from][e.to] = e.cost;
			cost[e.to][e.from] = -e.cost;
			capacity[e.from][e.to] = e.cap;
		}

		int flow = 0;
		int cost = 0;
		vector<int> d, p;
		while (flow < K) {
			shortest_paths(N, s, d, p);
			if (d[t] == INF)
				break;

			flow += 1;
			cost += d[t];
			int cur = t;
			while (cur != s) {
				capacity[p[cur]][cur] -= 1;
				capacity[cur][p[cur]] += 1;
				cur = p[cur];
			}
		}

		if (flow < K)
			return -1;
		else
			return cost;
	}
};

int matching_eval(VVI &matching) {
	int rv = 0;
	REP(r1, N) REP(c1, N) {
		int r2 = matching[r1][c1] / N;
		int c2 = matching[r1][c1] % N;
		rv += (r2-r1)*(r2-r1)+(c2-c1)*(c2-c1);
	}
	return rv;
}

int last_mcmf_cost;
VVI matching_mcmf(VVI &g, VVI &t, int randomness) {
	// DB(get_time() - start_time);
	MCMF mcmf;
	int v_s = 0;
	int v_t = 1;
	REP(i, N*N) mcmf.add_edge(v_s, 2 + i, 1, 0);
	REP(i, N*N) mcmf.add_edge(2 + N*N + i, v_t, 1, 0);
	
	bool type = randomness > 3 ? randomness % 2 : 0;
	REP(r1, N) REP(c1, N) REP(r2, N) REP(c2, N) if (t[r1][c1] == g[r2][c2]) {
		int mul = 1;
		int add = 5;
		int d1 = min(r1, N-1-r1);
		int d2 = min(c1, N-1-c1);
		// if (d1 == 0 || d2 == 0) mul *= 2;
		// if (d1 <= 1 || d2 <= 1) mul *= 2;
		// int av = mul * (abs(r1-r2) + abs(c1-c2));
		if (r1 == r2 && c1 == c2) add = 0;
		
		int av;
		if (type == 0) {
			// av = add + mul * (abs((r1-r2)*(r1-r2)*(r1-r2)) + abs((c1-c2)*(c1-c2)*(c1-c2))) + rng.next(max(1, randomness));
			av = add + mul * (abs((r1-r2)*(r1-r2)) + abs((c1-c2)*(c1-c2))) + rng.next(max(1, randomness));
		} else {
			av = add + mul * (abs(r1-r2) + abs(c1-c2));
		}
		mcmf.add_edge(2 + r1*N+c1, 2 + N*N + r2*N+c2, 1, av);
	}
	
	last_mcmf_cost = mcmf.run(2 + 2*N*N, N*N, v_s, v_t);
	
	VVI rv(N, VI(N));
	REP(i, N*N) REP(j, N*N) {
		if (mcmf.capacity[2 + N*N + j][2 + i])
			rv[j/N][j%N] = i;
	}
	
	return rv;
}

double bs_time = 0;
int bs_valid = 0;


#define GT unsigned char
struct BSState {
	GT g[MAX_N*MAX_N];
	int pm;
	int lp;
	int ld;
	int ep;
	int h;
	double v;
	
	BSState() { }
	
	BSState(int _pm, int _lp, int _ld, int _ep, int _h, double _v) {
		pm = _pm;
		lp = _lp;
		ld = _ld;
		ep = _ep;
		h = _h;
		v = _v;
	}
};

const int MAX_WIDTH=20000;
const int HMOD = 1e9+7;
// int bs_ht[HMOD];
// int bs_ht_no = 0;

BSState bs1[MAX_WIDTH*3];
BSState bs2[MAX_WIDTH*3];
BSState *cs;
BSState *ns;
int n_cs = 0;
int n_ns = 0;
string beam_search(VVI &g, int ev, int width, int early_exit=0) {
	VPII bs_pm;
	VVD dist_l1(N*N, VD(N*N));
	VVD dist_l2(N*N, VD(N*N));
	VVD dist_sq(N*N, VD(N*N));
	VVD dist_xx(N*N, VD(N*N));
	VI hv(N*N);
	REP(r1, N) REP(c1, N) REP(r2, N) REP(c2, N) {
		int a = abs(r1 - r2);
		int b = abs(c1 - c2);
		// int e1 = min(min(r1, N-1-r1), min(c1,N-1-c1));
		// int e2 = min(r2, N-1-r2)+ min(c2,N-1-c2);
		// int mul = 1;//e1 == 0 && e2 == 0 ? 4 : 1;
		// int mul = e2 ? 1 : 4;
		int mul = 1;
		int add = N <= 9 ? 0 : 1;// + (e1==0 && e2==0 ? 20 : 0);
		hv[r2*N+c2] = 1+r2*N+c2;
		dist_l1[r1*N+c1][r2*N+c2] = add+mul*(a+b);
		dist_l2[r1*N+c1][r2*N+c2] = add+mul*(sqrt(a*a+b*b));
		dist_sq[r1*N+c1][r2*N+c2] = add+mul*(a*a+b*b);
		// dist_xx[r1*N+c1][r2*N+c2] = add+mul*pow(a*a+b*b);
		// dist_xx[r1*N+c1][r2*N+c2] = add+mul*pow(a*a+b*b, 1.0);
		dist_xx[r1*N+c1][r2*N+c2] = add+mul*pow(a*a+b*b, 0.8);
	}
	
	hv[ev]=ev+1;
	REP(i, N*N) dist_l1[i][ev] = dist_l1[i][i] = 0;
	REP(i, N*N) dist_l2[i][ev] = dist_l2[i][i] = 0;
	REP(i, N*N) dist_sq[i][ev] = dist_sq[i][i] = 0;
	REP(i, N*N) dist_xx[i][ev] = dist_xx[i][i] = 0;
	
	auto state_score_final = [&](VI &s) -> double { double rv = 0; REP(i, N*N) rv += dist_l1[i][s[i]]; return rv; };
	auto state_score_l1 = [&](GT *s) -> double { double rv = 0; REP(i, N*N) rv += dist_l1[i][s[i]]; return rv; };
	auto state_score_l2 = [&](GT *s) -> double { double rv = 0; REP(i, N*N) rv += dist_l2[i][s[i]]; return rv; };
	auto state_score_sq = [&](GT *s) -> double { double rv = 0; REP(i, N*N) rv += dist_sq[i][s[i]];	return rv; };
	auto state_score_xx = [&](GT *s) -> double { double rv = 0;	REP(i, N*N) rv += dist_xx[i][s[i]];	return rv; };
	
	auto state_score = [&](GT *s) -> double {
		double rv = 0;
		REP(i, N*N) if (i != s[i] && (i == s[s[i]]) && i != ev && s[i] != ev) rv += 5;
		// REP(i, N*N) if (i != s[i] && (i == s[s[s[i]]]) && i != ev && s[i] != ev && s[s[i]] != ev) rv += 6;
		return rv + state_score_xx(s);
	};
	
	auto state_score_cell = [&](GT *s, int p) -> double {
		double rv = dist_xx[p][s[p]];
		if (p != s[p] && p != ev && s[p] != ev) {
			rv += p == s[s[p]] ? 10 : 0;
			// rv += p == s[s[s[p]]] && s[s[p]] != ev ? 18 : 0;
		}
		return rv;
	};
	
	VC<LL> mul; REP(i, N*N) mul.PB(rng.next(HMOD));
	auto state_hash = [&](GT *s) -> int {
		LL hash = 0;
		REP(i, N*N)	hash = (hash + hv[s[i]] * mul[i]) % HMOD;
		return hash;
	};	
	
	unordered_set<int> nh;
	
	cs = bs1;
	ns = bs2;
	n_cs = 0;
	n_ns = 0;
	
	BSState s0; 
	REP(r, N) REP(c, N) s0.g[r*N+c] = g[r][c];
	s0.ld = -1;
	REP(i, N*N) if (s0.g[i] == ev) s0.ep = i;
	s0.h = state_hash(s0.g);
	s0.v = state_score(s0.g);
	s0.pm = 0;
	cs[n_cs++] = s0;
	
	bs_pm.PB(MP(-1, -1));
	
	int moves = 0;
	
	int last_imp = 0;
	double bv = 1e30;
	
	VI bs;
	int bmp = -1;
	
	const int dp[] = {-N, +N, -1, +1};
	VVI good_dir(N*N, VI(4, 1));
	REP(i, N*N) {
		if (i / N == 0)   good_dir[i][0] = 0;
		if (i / N == N-1) good_dir[i][1] = 0;
		if (i % N == 0)   good_dir[i][2] = 0;
		if (i % N == N-1) good_dir[i][3] = 0;
	}
	
	VI order; order.PB(0);
	int cur_width = 1;
	while (true) {
		if (get_time() - start_time > cur_time_limit) return "";
		
		if (early_exit && moves >= early_exit)
			return string(early_exit+1, 'X');
		
		double av = cs[order[0]].v;
		if (av < bv) {
			bs.clear();
			REP(i, N*N)	bs.PB(cs[order[0]].g[i]);
			bv = av;
			last_imp = moves;
			if (av <= 1e-3) {
				bmp = cs[order[0]].pm;
				break;
			}
		}
		
		if (last_imp + 20 < moves) break;
		
		REP(i, cur_width) {
			BSState &s = cs[order[i]];
			auto &g = s.g;
			int ep = s.ep;
			assert(ep != -1);
			
			REP(d, 4) {
				if (!good_dir[ep][d] || s.ld == (d ^ 1)) continue;
				int np = ep+dp[d];
				LL h = s.h;
				double v = s.v;
				h -= hv[g[ep]] * mul[ep];
				h -= hv[g[np]] * mul[np];
				v -= state_score_cell(g, ep);
				v -= state_score_cell(g, np);
				swap(g[ep], g[np]);
				h += hv[g[ep]] * mul[ep];
				h += hv[g[np]] * mul[np];
				v += state_score_cell(g, ep);
				v += state_score_cell(g, np);
				h = (h % HMOD + HMOD) % HMOD;
				if (nh.emplace(h).Y) {
					ns[n_ns++] = BSState(bs_pm.S, order[i], d, np, h, v);
					bs_pm.PB(MP(s.pm, d));
				}
				swap(g[ep], g[np]);
			}
		}
		
		moves++;
		
		cur_width = min(n_ns, width);
		order.clear();
		REP(i, n_ns) order.PB(i);
		nth_element(order.begin(), order.begin() + cur_width, order.end(), [&](int a, int b) -> bool {return ns[a].v < ns[b].v;} );
		sort(order.begin(), order.begin() + cur_width, [&](int a, int b) -> bool {return ns[a].v < ns[b].v;} );
		
		REP(i, cur_width) {
			BSState &s = ns[order[i]];
			memcpy(s.g, cs[s.lp].g, sizeof(GT)*N*N);
			swap(s.g[s.ep], s.g[s.ep-dp[s.ld]]);
		}
		
		n_cs = 0;
		swap(n_cs, n_ns);
		swap(cs, ns);
		nh.clear();
	}
	
	
	int score = state_score_final(bs);
	
	// DB(score);
	// REP(i, N*N) {
		// int v = bs[i];
		// if (i == bs[i]) cerr << "-- "; else	cerr << v/N << v%N << " ";
		// if (i%N==N-1) cerr << endl;
	// }
	// DB(ev/N,ev%N);
	
	// DB(moves);
	
	if (score) return "";
	string rv = "";
	while (bmp) {
		rv += dn[bs_pm[bmp].Y];
		bmp = bs_pm[bmp].X;
	}
	reverse(ALL(rv));
	return rv;
}

void bs_benchmark(int _N, int runs, int width, int seed=1) {
	cur_time_limit = 1e9;
	N = _N;
	int good = 0;
	VI moves;
	double all_times = 0;
	double good_times = 0;
	RNG gen_rng(seed);
	VI v_ev;
	REP(r, N) REP(c, N) {
		if (min(r, c) < N/2-2 || max(r, c) > (N+1)/2+1) continue;
		v_ev.PB(r*N+c);
	}
	REP(run, runs) {
		VI v; REP(i, N*N) v.PB(i);
		REP(i, v.S) swap(v[i], v[gen_rng.next(i, v.S)]);
		VVI g(N, VI(N)); REP(r, N) REP(c, N) g[r][c] = v[r*N+c];
		int ev = v_ev[gen_rng.next(v_ev.S)];
		if (!solvable(g, ev)) {
			run--;
			continue;
		}
		double t = get_time();
		string s = beam_search(g, ev, width);
		all_times += get_time() - t;
		if (s.S) {
			good_times += get_time() - t;
			good++;
			moves.PB(s.S);
		}
	}
	
	cerr << endl << "Summary:" << endl;
	DB(runs);
	DB(good);
	sort(ALL(moves));
	DB(moves);
	int avg = 0; for (int x : moves) avg += x; avg /= max(1, (int)moves.S);
	DB(avg);
	all_times /= runs;
	good_times /= max(1, (int)moves.S);
	DB(all_times);
	DB(good_times);
}


int main(int argc, char **argv) {
	//Run benchmarks?	
	if (BS_BENCHMARK) {
		bs_benchmark(BSB_N, BSB_RUNS, BSB_WIDTH, 1);
		exit(0);	
	}
	
	if (HC_BENCHMARK) {
		//TODO: implement
		// hc_benchmark();
		exit(0);
	}
	
	//Read data
	cin >> N >> T;
	REP(r, N) {
		string s;
		cin >> s;
		REP(c, N) ogrid[r][c] = isdigit(s[c]) ? s[c] - '0' : s[c] - 'a' + 10;
	}
	
	cerr << "[DATA] N = " << N << endl;
	// if (N == 9)  WIDTH = WIDTH * 4 / 5;
	if (N == 10) WIDTH = WIDTH * 3 / 4;
	
	MAX_SCORE = (N*N-1) * log(N*N-1+1) - 1e-10;
	int hc_bc[MAX_N*MAX_N] = {0};
	REP(r, N) REP(c, N)	{
		int p = r*N+c;
		bad_empty[p] = min(r, c) < N/2-2 || max(r, c) > (N+1)/2+1;
		sg_l[p] = c > 0;
		sg_u[p] = r > 0;
		sg_r[p] = c < N-1;
		sg_d[p] = r < N-1;
		if (c == 0) hc_bc[p] |= 1;
		if (r == 0) hc_bc[p] |= 2;
		if (c == N-1) hc_bc[p] |= 4;
		if (r == N-1) hc_bc[p] |= 8;
	}
	
	string bsol;
	VC<pair<int, VVI>> all_matchings;
	VVI best_matching;
	double best_bstime = 0;
	
	const int RERUN_WIDTH = N <= 8 ? 3 : 0;
	while (true) {
		cur_time_limit = TIME_LIMIT - RERUN_WIDTH * best_bstime - (RERUN_WIDTH ? 0.0 : 0.0);
		
		double time_passed = get_time() - start_time;
		if (time_passed > cur_time_limit) break;
		
		VI g(N*N); REP(r, N) REP(c, N) g[r*N+c] = ogrid[r][c];
		VI fixed(N*N, 0); 
		if (g[0] && (g[0] & hc_bc[0]) == 0) fixed[0] = 1; 
		if (g[N-1] && (g[N-1] & hc_bc[N-1]) == 0) fixed[N-1] = 1; 
		if (g[N*N-N] && (g[N*N-N] & hc_bc[N*N-N]) == 0) fixed[N*N-N] = 1; 
		if (g[N*N-1] && (g[N*N-1] & hc_bc[N*N-1]) == 0) fixed[N*N-1] = 1; 
		REP(i, N*N) {
			int p = rng.next(i, g.S);
			if (fixed[i] || fixed[p]) continue;
			swap(g[i], g[p]);
		}
		VI cn(N*N); REP(i, N*N) update_connection(g, cn, i);
		double bv = score_grid(g, cn);
		double xv = bv;
		int step = 0;
		int last_acc = 0;
		
		double timer_hc_start = get_time();
		while (true) {
			double time_passed = (get_time() - start_time) / cur_time_limit;
			if (time_passed > 1.0) break;
			
			double t = (time_passed * 25) - floor(time_passed * 25);
			
			int kick_steps = max(1, (int)(BASE_KICK_STEPS * t * t * N * N * N * N));
			
			step++;
			int p0 = rng.next(N*N);
			int p1 = rng.next(N*N);
			if (fixed[p0] || fixed[p1]) continue;
			if (g[p0] == g[p1]) continue;
			if (g[p1] & hc_bc[p0]) continue;
			if (g[p0] & hc_bc[p1]) continue;
			swap(g[p0], g[p1]);
			update_connection_neighbors(g, cn, p0);
			update_connection_neighbors(g, cn, p1);
			double av = score_grid(g, cn);
			if (av >= bv || step > last_acc + kick_steps) {
				if (av > bv || step > last_acc + kick_steps) {
					last_acc = step;
				}
				if (av > xv) {
					xv = av;
					if (grid_fully_solved) break;
				}
				bv = av;
			} else {
				swap(g[p0], g[p1]);
				update_connection_neighbors(g, cn, p0);
				update_connection_neighbors(g, cn, p1);
			}
		}
		timer_hc += get_time() - timer_hc_start;
		stat_hc_steps += step;
		
		if (!grid_fully_solved) break;
		VVI bgrid(N, VI(N)); REP(r, N) REP(c, N) bgrid[r][c] = g[r*N+c];
		// score_grid(g, true);
		// show_grid(bgrid);
		
		stat_grids++;
		
		if (ONLY_GRIDS) continue;
		
		VVI matching;
		set<VVI> matching_vs;
		const int MAX_TRIES = N <= 8 ? 5 : 10;
		int n_tries = 0;
		PII empty_pos; 
		
		again:
		init_grid();
		n_tries++;
		REP(r, N) REP(c, N) if (ogrid[r][c] == 0) empty_pos = {r, c};
		string presolve = "";
		
		if (PRESOLVE_N && N > PRESOLVE_N && bsol.S == 0 && n_tries < MAX_TRIES && get_time() - start_time < cur_time_limit) {
			init_grid();
			double timer_heur_start = get_time();
			presolve = solve_heuristic(grid, bgrid, N - PRESOLVE_N);
			timer_heur += get_time() - timer_heur_start;
			if (presolve.S == 0) goto again;
			init_grid();
			exec_moves(grid, presolve);
			REP(r, N) REP(c, N) if (grid[r][c] == 0) empty_pos = {r, c};
		}
		while (n_tries < MAX_TRIES && get_time() - start_time < cur_time_limit) {
			double timer_wm_start = get_time();
			matching = matching_mcmf(grid, bgrid, n_tries);
			timer_wm += get_time() - timer_wm_start;
			if (solvable(matching, matching[empty_pos.X][empty_pos.Y]))
				break;
			n_tries++;
		}
		if (n_tries >= MAX_TRIES || get_time() - start_time > cur_time_limit) continue;
		
		// all_matchings.PB(MP(matching_eval(matching), matching));
		// continue;
		
		if (matching_vs.count(matching)) {
			stat_matching_repeat++;
			goto again;
		}
		matching_vs.insert(matching);
		
		double timer_bs_start = get_time();
		string sol = beam_search(matching, matching[empty_pos.X][empty_pos.Y], WIDTH, BS_EARLYEXIT ? bsol.S : 0);
		double bstime = get_time() - timer_bs_start;
		timer_bs += bstime;
		stat_bs_runs++;
		
		if (!sol.S) goto again;
		if (presolve.S) sol = presolve + sol;
		
		if (!presolve.S) cerr << "Solve: " << sol.S << " MEval << " << matching_eval(matching) << " (" << last_mcmf_cost << ") Time: " << get_time() - start_time << endl;
		
		if (bsol.S == 0 || sol.S < bsol.S) {
			bsol = sol;
			if (presolve.S) goto again;
			best_matching = matching;
			best_bstime = bstime;
		}
		
	}
	
	// cur_time_limit = 1e9;
	// sort(ALL(all_matchings));
	// DB(all_matchings.S);
	// for (auto &m : all_matchings) {
		// init_grid();
		// PII empty_pos = {-1, -1}; REP(r, N) REP(c, N) if (ogrid[r][c] == 0) empty_pos = {r, c};
		// int meval = m.X;
		// string bseval1000 = beam_search(m.Y, m.Y[empty_pos.X][empty_pos.Y], 1000);
		// string bseval2000 = beam_search(m.Y, m.Y[empty_pos.X][empty_pos.Y], 2000);
		// string bseval5000 = beam_search(m.Y, m.Y[empty_pos.X][empty_pos.Y], 5000);
		// DB(meval, bseval1000.S, bseval2000.S, bseval5000.S);
		// if (bseval1000.S && (bsol.S == 0 || bseval1000.S < bsol.S)) bsol = bseval1000;
		// if (bseval2000.S && (bsol.S == 0 || bseval2000.S < bsol.S)) bsol = bseval2000;
		// if (bseval5000.S && (bsol.S == 0 || bseval5000.S < bsol.S)) bsol = bseval5000;
	// }
	
	
	if (RERUN_WIDTH && best_matching.S) {
		cur_time_limit = 2.95 * TIME_LIMIT_MUL;
		double time_left = cur_time_limit - (get_time() - start_time);
		cerr << "BS Rerun: Time left: " << time_left << " Expected: " << RERUN_WIDTH * best_bstime << endl;
		init_grid();
		int ev = -1; REP(r, N) REP(c, N) if (grid[r][c] == 0) ev = best_matching[r][c];
		string s = beam_search(best_matching, ev, RERUN_WIDTH * WIDTH);
		cerr << "BS Rerun: " << bsol.S << " -> " << s.S << endl;
		if (s.S && s.S < bsol.S) bsol = s;
	}
	
	cerr << "[DATA] time = " << (get_time() - start_time) << endl;
	cerr << "[DATA] timer_hc = " << timer_hc << endl;
	cerr << "[DATA] timer_bs = " << timer_bs << endl;
	cerr << "[DATA] timer_wm = " << timer_wm << endl;
	cerr << "[DATA] timer_heur = " << timer_heur << endl;
	cerr << "[DATA] stat_matching_repeat = " << stat_matching_repeat << endl;
	cerr << "[DATA] stat_grids = " << stat_grids << endl;
	cerr << "[DATA] stat_bs_runs = " << stat_bs_runs << endl;
	cerr << "[DATA] stat_hc_steps = " << stat_hc_steps << endl;
	DB(best_bstime);

	DB(bsol.S);
	
	if (bsol.S || !FORCE_FAIL) {
		cout << bsol << endl;
	} else {
		cout << "UUUUUUUUUUUUUUUUUUUUUUUU" << endl;
	}
		
	return 0;
}
