// Author: Psyho
// Twitter: https://twitter.com/fakepsyho
 
// TEMPLATE

#pragma GCC optimize "Ofast,omit-frame-pointer,inline,unroll-all-loops"

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

// SOLUTION

const int MAX_N = 32;
const int MAX_DIST = MAX_N*MAX_N;
const int MAX_AREA = MAX_N*MAX_N;

const int dd[] = {-1, -MAX_N, 1, MAX_N};
const char dc[] = {'L', 'U', 'R', 'D'};


int N;
int C;
double elfP;
int start_money;
int money;

char grid[MAX_N*MAX_N];
char prev_grid[MAX_N*MAX_N];

VI allp;

bool no_box_mode = false;
bool no_present_mode = false;

INLINE int p1d(int r, int c) {return (r+1)*MAX_N+(c+1);}
INLINE PII p2d(int p) {return MP(p/MAX_N-1,(p-1)&(MAX_N-1));}

bool is_elf(char c) {
	return c == 'e' || c == 'E' || c == 'B';
}

int near(int p, char c) {
	int rv = 0;
	REP(d, 4) rv += grid[p+dd[d]] == c;
	return rv;
}

int bfs_rv[MAX_N*MAX_N];
void full_bfs(VI &starts, bool ignore_elves=false, bool ignore_empty_elves=false) {
	static int q[MAX_N*MAX_N*2*4];
	int q_st = 0;
	int q_en = 0;
	
	memset(bfs_rv, 0x3F, sizeof(int)*MAX_N*(N+2));
	
	for (int x : starts) {
		q[q_en++] = x;
		bfs_rv[x] = 0;
	}
	
	while (q_st < q_en) {
		int x = q[q_st++];
		int v = bfs_rv[x];
		REP(d, 4) {
			int nx = x + dd[d];
			if (nx < 0 || nx >= MAX_DIST || grid[nx] == 'T' || grid[nx] == 'b' || grid[nx] == 'P' || grid[nx] == 'X') continue;
			if (v+1 < bfs_rv[nx]) {
				// if (ignore_elves || (grid[nx] != 'e' && grid[nx] != 'E' && grid[nx] != 'B')) q[q_en++] = nx;
				if ((ignore_empty_elves || grid[nx] != 'e') && (ignore_elves || (grid[nx] != 'E' && grid[nx] != 'B'))) q[q_en++] = nx;
				bfs_rv[nx] = v+1;
			}
		}
	}
}

void full_bfs(string starts, bool ignore_elves=false, bool ignore_empty_elves=false) {
	static VI vi_starts;
	vi_starts.clear();
	FOR(r, -1, N+1) FOR(c, -1, N+1) {
		int p = p1d(r, c);
		for (char c : starts) if (grid[p] == c) {
			vi_starts.PB(p);
			break;
		}
	}
	full_bfs(vi_starts, ignore_elves, ignore_empty_elves);
}

int find_move(int p, int* dist, bool allow_stop=false, bool balance=false) {
	assert(grid[p] == 'e' || grid[p] == 'E' || grid[p] == 'B');
	int x = dist[p];
	int bd = -1;
	int bv = -1;
	REP(d, 4) {
		int np = p + dd[d];
		if (dist[np] == x - 1 && grid[np] != 'e' && grid[np] != 'E' && grid[np] != 'B') {
			if (no_box_mode && grid[np] == 'b') continue;
			if (no_present_mode && grid[np] == 'P') continue;
			int av = 0; 
			if (balance) REP(d2, 4) av += dist[np+dd[d2]] == x - 2;
			if (av > bv) bd = d, bv = av;
		}
	}
	if (bd != -1) return bd;
	if (allow_stop) return -1;
	assert(false);
}

int find_target(int p, int d0, int* dist) {
	p += dd[d0];
	while (dist[p]) {
		int np = p;
		REP(d, 4) {
			if (dist[p+dd[d]] + 1 == dist[p] && !is_elf(grid[p+dd[d]])) {
				np = p+dd[d];
			}
		}
		p = np;
	}
	return p;
}

void make_move(int p, int d) {
	int np = p + dd[d];
	if (grid[p] == 'e') {
		if (grid[np] == '.') {
			grid[np] = grid[p];
		} else if (grid[np] == 'b') {
			grid[np] = 'B';
		} else if (grid[np] == 'P') {
			grid[np] = 'E';
		} else {
			assert(false);
		}
	} else if (grid[p] == 'B' || grid[p] == 'E') {
		if (grid[np] == '.') {
			grid[np] = grid[p];
		} else {
			assert(grid[np] == 'X');
		}
	} else {
		assert(false);
	}
	grid[p] = '.';
}

void show_grid(char *g) {
	cerr << endl;
	REP(r, N) {
		REP(c, N) cerr << g[p1d(r,c)];
		cerr << endl;
	}
}

int ai1mode(char* g, int turn, int start_money, int start_presents) {
	int N2 = N*N;
    double nLose = (N2-N)*(elfP - 1.0/C) - (double)start_money/C + 3,
           cut1 = (start_presents-nLose)/start_presents,
           cut2 = start_presents/((N2-start_money/C)*elfP*C);
	
	if (!(cut1 > 0 && turn < N2-6)) return 0;
	if (cut2 > 0.08) return 1;
	return 2;
}


VI ai2turn(char *g, int turn, int money) {
	if (money < C) return {};
	
	
	auto ai2score = [&]() -> int {
		static int q[MAX_N*MAX_N*2*4];
		static int reachable[MAX_N*MAX_N];
		int q_st = 0;
		int q_en = 0;
		
		memset(reachable, 0, sizeof(reachable));
		
		REP(r, N) REP(c, N) if (g[p1d(r,c)] == '.' && (r == 0 || r == N-1 || c == 0 || c == N-1)) q[q_en++] = p1d(r,c), reachable[p1d(r,c)] = 1;
		
		while (q_st < q_en) {
			int x = q[q_st++];
			REP(d, 4) {
				int nx = x + dd[d];
				if (g[nx] == 'X' || g[nx] == 'T' || g[nx] == 'b' || g[nx] == 'P' || reachable[nx]) continue;
				q[q_en++] = nx;
				reachable[nx] = 1;
			}
		}
		
		int rv = 0;
		for (int p : allp) if (g[p] == 'P') REP(d, 4) rv += reachable[p+dd[d]];
		return rv;
    };
	
	int best = ai2score();
	int bp = -1;
	FOR(r, 1, N-1) FOR(c, 1, N-1) {
		int p = p1d(r,c);
		if (g[p] != '.') continue;
		g[p] = 'b';
		int av = ai2score();
		if (av < best) {
			best = av;
			bp = p;
		}
		g[p] = '.';
	}
	
	if (bp != -1) return {bp};
	return {};
}

VI ai3turn(char* g, int turn, int money) {
	const int MaxDist = 1000;
	VVI Dist=VVI(N,VI(N));
	
	auto ai3score = [&]() -> int {
		static int q[MAX_N*MAX_N*2*4];
		static int ai3_rv[MAX_N*MAX_N];
		int q_st = 0;
		int q_en = 0;
		
		memset(ai3_rv, 0x3F, sizeof(ai3_rv));
		
		REP(r, N) REP(c, N) if (g[p1d(r,c)] == '.' && (r == 0 || r == N-1 || c == 0 || c == N-1)) q[q_en++] = p1d(r,c), ai3_rv[p1d(r,c)] = 0;
		
		while (q_st < q_en) {
			int x = q[q_st++];
			int v = ai3_rv[x];
			REP(d, 4) {
				int nx = x + dd[d];
				if (g[nx] == 'X' || g[nx] == 'T' || g[nx] == 'b') continue;
				if (v+1 < ai3_rv[nx]) {
					if (g[nx] != 'P') q[q_en++] = nx;
					ai3_rv[nx] = v+1;
				}
			}
		}
		
		REP(r,N) REP(c,N) Dist[r][c] = min(ai3_rv[p1d(r,c)], MaxDist);
		
		int rv=0;
		REP(r,N) REP(c,N) if (g[p1d(r,c)] == 'P') rv += Dist[r][c];
		return rv;
	};
	
	if (money < C) return {};
    int bestScore=ai3score();
	int initial=bestScore;
    int best=-1;

    VVI startDist=Dist;

    bool allSaved=true;
	REP(r,N) REP(c,N) if (g[p1d(r,c)] == 'P' && Dist[r][c] != MaxDist) {
		allSaved=false;
		break;
	}    
	
    if (allSaved) return {};
	
	FOR(r, 1, N-1) FOR(c, 1, N-1) if (g[p1d(r,c)] == '.' && startDist[r][c] < MaxDist) {
		g[p1d(r,c)] = 'b';
		int score = ai3score();
		g[p1d(r,c)] = '.';
		
		if (score > bestScore) {
			bestScore=score;
			best=p1d(r,c);
		}
	}
	
	if (best != -1) return {best};
	
	//placeGreedy() is non-deterministic :/
	
	return {-1};
}

VI ai4turn(char* g, int turn, int money) {
	if (money < C) return {};
	
	REP(c, N) REP(r, N) {
		int p = p1d(r,c);
		if (g[p] != 'P') continue;
		bool bad = false;
		for (int dp : dd) bad |= g[p+dp] == 'e';
		if (bad) continue;
		
		static VI dirs = {-1,+MAX_N,+1,-MAX_N};
		int rv = -1;
		for (int dp : dirs) if (g[p+dp]=='.') for (int dp2 : dirs) if (g[p+dp+dp2] == 'e') rv = (rv == -1 || rv == p+dp) ? p+dp : -2;
		if (rv >= 0) return {rv};
	}
	
	FOR(c, 1, N-2) FOR(r, 1, N-2) {
		int p = p1d(r,c);
		if (g[p] != '.') continue;
		
		bool doomed = false;
		for (int dp : dd) doomed |= g[p+dp] == 'e';
		
		int wings = (g[p-MAX_N-1] == 'P') + (g[p+MAX_N+1] == 'P');
		int exposed = 0;
		if (g[p-MAX_N-1] == 'P' && g[p-1+MAX_N] == 'P') {
			exposed = (g[p-1] == '.') + (g[p-MAX_N] == '.') + (g[p+MAX_N] == '.');
		} else if (g[p-1-MAX_N] == 'P' && g[p+1-MAX_N] == 'P') {
			exposed = (g[p-MAX_N] == '.') + (g[p-1] == '.') + (g[p+1] == '.');
		} else if (g[p-1+MAX_N] == 'P' && g[p+1+MAX_N] == 'P') {
			exposed = (g[p+MAX_N] == '.') + (g[p-1] == '.') + (g[p+1] == '.');
		} else if (g[p+1-MAX_N] == 'P' && g[p+1+MAX_N] == 'P') {
			exposed = (g[p+1] == '.') + (g[p-MAX_N] == '.') + (g[p+MAX_N] == '.');
		}
		
		if (doomed && wings == 1 && exposed > 1) return {p};
	}
	return {};
}

VI ai4good_spots() {
	static char tmpg[MAX_AREA];
	REP(i, MAX_AREA) tmpg[i] = is_elf(grid[i]) ? '.' : grid[i];
	if (ai4turn(tmpg, 0, 1000).S) return {};
	
	VI rv;
	for (int p : allp) if (tmpg[p] == '.' || tmpg[p] == 'P') {
		char c = tmpg[p];
		tmpg[p] = 'e';
		if (ai4turn(tmpg, 0, 1000).S == 0) rv.PB(p);
		tmpg[p] = c;
	}
	return rv;
}

VI ai5turn(char* g, int turn, int money) {
	int c = 1 + ((turn * (7919)) % (N-2));
	int r = 1 + ((turn * (50091)) % (N-2));
	VI rv;
	if (money >= C && g[p1d(r,c)] == '.') rv.PB(p1d(r,c));
	return rv;
}

int main() {
	ios_base::sync_with_stdio(false);
	
	cin >> N >> C >> elfP >> start_money;
	
	cerr << "[DATA] N = " << N << endl;
	cerr << "[DATA] C = " << C << endl;
	cerr << "[DATA] P = " << elfP << endl;
	cerr << "[DATA] M = " << start_money << endl;
	cerr << "[DATA] CP = " << C*elfP << endl;
	
	memset(grid, 'X', sizeof(grid));
	
	int og_boxes[MAX_N*MAX_N];
	ZERO(og_boxes);
	
	int AI = -1;
	
	int start_presents = 0;
	
	VI AI_possible = {0,1,1,1,1,1};
	
	static VVI elves_sidestep;
	REP(turn, N*N) {
		int elapsed_time;
		cin >> elapsed_time >> money;
		
		REP(r, N) REP(c, N) cin >> grid[p1d(r,c)];
		if (turn == 0) {
			REP(r, N) REP(c, N) if (grid[p1d(r,c)] != 'T') allp.PB(p1d(r,c));
			REP(i, MAX_AREA) prev_grid[i] = grid[i];
			REP(i, MAX_AREA) if (prev_grid[i] == 'b') prev_grid[i] = '.';
		}
		
		if (turn == 0) {
			for (int p : allp) start_presents += grid[p] == 'P';
			int n_boxes = 0; for (int p : allp) n_boxes += grid[p] == 'b';
			if (n_boxes == 0) {
				// if (ai1mode(prev_grid, turn, start_money, start_presents) != 0) AI_possible[1] = 0;
				AI_possible[2] = 0;
				AI_possible[3] = 0;
			}
			if (n_boxes > 1) {
				AI_possible[2] = 0;
				AI_possible[3] = 0;
				AI_possible[4] = 0;
				AI_possible[5] = 0;
			}
			cerr << "[DATA] AI1MODE = " << ai1mode(prev_grid, turn, start_money, start_presents) << endl;
		}
		
		if (true || AI == -1) {
			VI new_boxes; 
			for (int p : allp) if (grid[p] == 'b' && prev_grid[p] == '.') new_boxes.PB(p);
			REP(i, N) {
				if (grid[p1d(0,i)] == 'e') prev_grid[p1d(0,i)] = 'e';
				if (grid[p1d(N-1,i)] == 'e') prev_grid[p1d(N-1,i)] = 'e';
				if (grid[p1d(i,0)] == 'e') prev_grid[p1d(i,0)] = 'e';
				if (grid[p1d(i,N-1)] == 'e') prev_grid[p1d(i,N-1)] = 'e';
			}
			REP(r, N) REP(c, N) {
				int p = p1d(r, c);
				if (!(prev_grid[p] == grid[p] || prev_grid[p] == '.' && grid[p] == 'b')) {
					DB(r,c,p,prev_grid[p],grid[p]);
				}
			}
			sort(ALL(new_boxes));
			
			// if (AI_possible[1]) {
				// if (ai1mode(prev_grid, turn, start_money, start_presents) == 0 && new_boxes.S) AI_possible[1] = 0;
			// }
			if (AI_possible[2]) {
				VI pred_boxes = ai2turn(prev_grid, turn, money + new_boxes.S * C);
				sort(ALL(pred_boxes));
				if (pred_boxes != new_boxes) AI_possible[2] = 0;
			}
			if (AI_possible[3]) {
				VI pred_boxes = ai3turn(prev_grid, turn, money + new_boxes.S * C);
				sort(ALL(pred_boxes));
				if (pred_boxes != VI(1,-1) && pred_boxes != new_boxes) AI_possible[3] = 0;
			}
			if (AI_possible[4]) {
				VI pred_boxes = ai4turn(prev_grid, turn, money + new_boxes.S * C);
				sort(ALL(pred_boxes));
				if (pred_boxes != new_boxes) AI_possible[4] = 0;
			}
			if (AI_possible[5]) {
				VI pred_boxes = ai5turn(prev_grid, turn, money + new_boxes.S * C);
				sort(ALL(pred_boxes));
				if (pred_boxes != new_boxes) AI_possible[5] = 0;
			}
			
			AI = 0;
			FOR(i, 1, 6) if (AI_possible[i]) AI = AI == 0 ? i : -1;
			// assert(AI != 0);
		}
		
		if (turn == N*N-1) DB(AI_possible);
		
		VPII moves;
		
		static VI elves_present;
		static VI elves_boxes;
		static VI elves;
		elves_present.clear();
		elves_boxes.clear();
		elves.clear();
		
		unordered_set<int> bad_elves;
		for (VI &v : elves_sidestep) bad_elves.insert(v[0]);
		
		bool empty_border = false;
		REP(i, N) {
			empty_border |= grid[p1d(i,0)] == '.';
			empty_border |= grid[p1d(i,N-1)] == '.';
			empty_border |= grid[p1d(0,i)] == '.';
			empty_border |= grid[p1d(N-1,i)] == '.';
		}
		
		bool reachable_presents = false;
		full_bfs("P");
		for (int p : allp) if (grid[p] == 'e' && bfs_rv[p] < MAX_DIST) reachable_presents = true;
		
		no_box_mode = false;
		no_present_mode = false;
		for (int p : allp) {
			bool wait = false;
			bool greedy_wait = false;
			bool delayed_wait = false;
			if (AI_possible[2]) {
				wait = C*elfP < 1.2 && turn < N*N*.75;
				greedy_wait = true;
				delayed_wait = true;
			} else if (AI_possible[3]) {
				wait = C*elfP < 0.9 && turn < N*N*.8;
				greedy_wait = true;
				delayed_wait = true;
			} else if (AI_possible[4]) {
				wait = C*elfP < 1.0 && turn < N*N*.7;
			}
			if (wait && delayed_wait) {
				wait = false;
				no_box_mode = true;
				// no_present_mode = true;
			}
			if (reachable_presents && greedy_wait) wait = false;
			if (bad_elves.count(p)) continue;
			if (grid[p] == 'e') {
				if (!wait) elves.PB(p);
			} else if (grid[p] == 'E') {
				elves_present.PB(p);
			} else if (grid[p] == 'B') {
				elves_boxes.PB(p);
			}
		}
		
		int presents_left = 0; REP(r, N) REP(c, N) if (grid[p1d(r,c)] == 'P') presents_left++;
		
		int alive_zone[MAX_N*MAX_N];
		int alive_box[MAX_N*MAX_N];
		ZERO(alive_zone);
		ZERO(alive_box);
		full_bfs("P");
		for (int p : allp) {
			alive_zone[p] = bfs_rv[p] < MAX_DIST;
			REP(d, 4) alive_box[p] |= bfs_rv[p+dd[d]] < MAX_DIST;
		}
		
		int reachable[MAX_N*MAX_N];
		full_bfs("X");
		REP(i, MAX_N*MAX_N) reachable[i] = bfs_rv[i];
		
		for (int p : allp) {
			if (og_boxes[p]) {
				if (!alive_box[p]) og_boxes[p] = 0;
				bool all_alive = true;
				for (int dp : dd) all_alive &= grid[p+dp] == 'T' || grid[p+dp] == '.' && alive_zone[p+dp];
				if (all_alive) og_boxes[p] = 0;
			}
		}
		
		full_bfs("X");
		if (turn < N*N-3*N) {
			VPII vp;
			for (int p : allp) if (grid[p] == 'P') {
				int md = MAX_DIST;
				REP(d, 4) md = min(md, bfs_rv[p+dd[d]]);
				if (md < MAX_DIST) vp.PB(MP(md, p));
			}
			sort(ALL(vp));
			int n_remove = (int)(vp.S * .3);
			if (AI_possible[4] == 1) n_remove = 0;
			if (AI_possible[3] == 1) n_remove = (int)vp.S - 5;
			if (AI_possible[2] == 1) n_remove = (int)vp.S - 5;
			REP(i, n_remove) grid[vp[i].Y] = 'T';
		}
		
		while (elves_present.S || elves_boxes.S) {
			// sidestep to make room for elf with present
			full_bfs("X", true, false);
			int best = -1;
			
			REP(i, elves_boxes.S) {
				if (og_boxes[elves_boxes[i]]) {
					int unmoved_elves = 0;
					for (int x : dd) {
						int p = elves_boxes[i] + x;
						if (reachable[p] >= MAX_DIST) for (int p2 : elves_present) unmoved_elves += p == p2;
					}
					if (unmoved_elves > 0) {
						if (best == -1 || bfs_rv[elves_boxes[i]] < bfs_rv[elves_boxes[best]]) best = i;
					}
					
				}
			}
			
			if (best != -1) {
				og_boxes[elves_boxes[best]] = 0;
				int bd = -1;
				REP(d, 4) if (grid[elves_boxes[best]+dd[d]] == '.') bd = d;
				REP(d, 4) if (grid[elves_boxes[best]+dd[d]] == '.' && alive_zone[elves_boxes[best]+dd[d]]) bd = d;
				
				if (bd != -1) {
					elves_sidestep.PB({elves_boxes[best]+dd[bd], elves_boxes[best], -1});
					moves.PB(MP(elves_boxes[best], bd));
					make_move(elves_boxes[best], bd);
				}
				
				elves_boxes[best] = elves_boxes.back();
				elves_boxes.pop_back();
				continue;
			}
			
			// return present
			full_bfs("X", true, false);
			REP(i, elves_present.S) {
				
				// if (bfs_rv[elves_present[i]] > MAX_DIST) {
				if (bfs_rv[elves_present[i]] > MAX_DIST || bfs_rv[elves_present[i]] == 1 && presents_left == 0 && turn < N*N-1) {
					elves_present[i] = elves_present.back();
					elves_present.pop_back();
					i--;
					continue;
				}
				
				if (best == -1 || bfs_rv[elves_present[i]] < bfs_rv[elves_present[best]]) best = i;
			}
			if (best != -1) {
				int d = find_move(elves_present[best], bfs_rv, true);
				
				if (d == -1) {
					full_bfs("X", false, false);
					d = find_move(elves_present[best], bfs_rv, true);
				}
				
				if (d != -1) {
					moves.PB(MP(elves_present[best], d));
					make_move(elves_present[best], d);
				}
				
				elves_present[best] = elves_present.back();
				elves_present.pop_back();
				continue;
			}
			
			// return box
			full_bfs("X", true, false);
			REP(i, elves_boxes.S) {
				bool near_dead_elf = false;
				for (int dp : dd) near_dead_elf |= grid[elves_boxes[i]+dp] == 'e' && !alive_zone[elves_boxes[i]+dp];
				
				if (og_boxes[elves_boxes[i]] && !near_dead_elf) {
					elves_boxes[i] = elves_boxes.back();
					elves_boxes.pop_back();
					i--;
					continue;
				}
				
				// sidestep to make room for elf
				// if (og_boxes[elves_boxes[i]] && near_dead_elf) {
				if (og_boxes[elves_boxes[i]] && near_dead_elf && (AI_possible[4] == 0 || bfs_rv[elves_boxes[i]] > MAX_DIST)) {
					int bd = -1;
					REP(d, 4) if (grid[elves_boxes[i]+dd[d]] == '.') bd = d;
					REP(d, 4) if (grid[elves_boxes[i]+dd[d]] == '.' && !alive_zone[elves_boxes[i]+dd[d]]) bd = d;
					
					if (bd != -1) {
						elves_sidestep.PB({elves_boxes[i]+dd[bd], elves_boxes[i], -1});
						moves.PB(MP(elves_boxes[i], bd));
						make_move(elves_boxes[i], bd);
						og_boxes[elves_boxes[i]] = 0;
						elves_boxes[i] = elves_boxes.back();
						elves_boxes.pop_back();
						i--;
						continue;
					}
				}
				
				if (bfs_rv[elves_boxes[i]] > MAX_DIST) {
					elves_boxes[i] = elves_boxes.back();
					elves_boxes.pop_back();
					i--;
					continue;
				}
				
				if (best == -1 || bfs_rv[elves_boxes[i]] < bfs_rv[elves_boxes[best]]) best = i;
			}
			if (best != -1) {
				og_boxes[elves_boxes[best]] = 0;
				
				int d = find_move(elves_boxes[best], bfs_rv, true);
				if (d == -1) {
					full_bfs("X", false, false);
					d = find_move(elves_boxes[best], bfs_rv, true);
				}
				if (d != -1) {
					moves.PB(MP(elves_boxes[best], d));
					make_move(elves_boxes[best], d);
				}
				
				elves_boxes[best] = elves_boxes.back();
				elves_boxes.pop_back();
				continue;
			}
		}
		
		// grab available present
		while (true) {
			int present_pos = -1; for (int p : allp) if (grid[p] == 'P') present_pos = present_pos == -1 ? p : -2;
			
			// if last present then delay
			int elf_near = -1;
			if (present_pos >= 0) REP(d, 4) if (grid[present_pos + dd[d]] == 'e') elf_near = present_pos + dd[d];
			if (elf_near != -1 && turn < N*N - N) {
				REP(i, elves.S) if (elves[i] == elf_near) {
					elves[i] = elves.back();
					elves.pop_back();
					break;
				}
				break;
			}
			
			full_bfs("P", false, true);
			int best = -1;
			REP(i, elves.S) {
				int p = elves[i];
				if (bfs_rv[p] < MAX_DIST && (best == -1 || bfs_rv[p] < bfs_rv[elves[best]])) best = i;
			}
			
			if (best == -1) break;
			int d = find_move(elves[best], bfs_rv, true);
			if (d != -1) {
				moves.PB(MP(elves[best], d));
				make_move(elves[best], d);
			}
			
			elves[best] = elves.back();
			elves.pop_back();
		}
		
		// go to elf with a box (to switch pos)
		int box_dist[MAX_N*MAX_N];
		for (int min_empty = 2; min_empty >= 2; min_empty--) {
			while (true) {
				full_bfs("b", false, true);
				REP(i, MAX_N*MAX_N) box_dist[i] = bfs_rv[i];
				
				VI starts; 
				for (int p : allp) {
					if (!og_boxes[p]) continue;
					if (near(p, '.') >= min_empty) REP(d, 4) if (grid[p+dd[d]] == '.') starts.PB(p+dd[d]);
				}
				full_bfs(starts);
				
				int best = -1;
				REP(i, elves.S) {
					int p = elves[i];
					if (bfs_rv[p] < MAX_DIST && bfs_rv[p] < box_dist[p] + N/2 && (best == -1 || bfs_rv[p] < bfs_rv[elves[best]])) best = i;
				}
				
				if (best == -1) break;
				int d = find_move(elves[best], bfs_rv, true);
				if (d != -1) {
					moves.PB(MP(elves[best], d));
					make_move(elves[best], d);
				}
				
				elves[best] = elves.back();
				elves.pop_back();
			}
		}
		
		// grab box leading to present
		for (int min_empty = (AI==1?4:1); min_empty >= 0; min_empty--) {
			while (true) {
				VI starts; 
				for (int p : allp) if (grid[p] == 'b' && alive_box[p] && near(p, '.') >= min_empty) starts.PB(p);
				full_bfs(starts, false, true);
				// full_bfs("b");
				
				int best = -1;
				REP(i, elves.S) {
					int p = elves[i];
					if (bfs_rv[p] < MAX_DIST && (best == -1 || bfs_rv[p] < bfs_rv[elves[best]])) best = i;
				}
				
				if (best == -1) break;
				int d = find_move(elves[best], bfs_rv, true);
				if (d != -1) {
					moves.PB(MP(elves[best], d));
					make_move(elves[best], d);
				}
				
				if (grid[elves[best]+dd[d]] == 'B') og_boxes[elves[best]+dd[d]] = 1;
				
				elves[best] = elves.back();
				elves.pop_back();
			}
		}
		
		// grab any box
		while (true) {
			full_bfs("b");
			
			int best = -1;
			REP(i, elves.S) {
				int p = elves[i];
				if (bfs_rv[p] < MAX_DIST && (best == -1 || bfs_rv[p] < bfs_rv[elves[best]])) best = i;
			}
			
			if (best == -1) break;
			int d = find_move(elves[best], bfs_rv, true);
			if (d != -1) {
				moves.PB(MP(elves[best], d));
				make_move(elves[best], d);
			}
			
			elves[best] = elves.back();
			elves.pop_back();
		}
		
		// move back to position
		REP(i, elves_sidestep.S) {
			VI &v = elves_sidestep[i];
			if (grid[v[1]] == '.' && v[2] >= 0) {
				int d = 0; while (v[0]+dd[d] != v[1]) d++;
				moves.PB(MP(v[0], d));
				make_move(v[0], d);
				og_boxes[v[1]] = 1;
				
				elves_sidestep[i] = elves_sidestep.back();
				elves_sidestep.pop_back();
				i--;
			} else {
				v[2]++;
				if (v[2] >= (AI==1?7:2)) {
					elves_sidestep[i] = elves_sidestep.back();
					elves_sidestep.pop_back();
					i--;
				}
			}
		}
		
		REP(i, MAX_AREA) if (grid[i] == 'T' && prev_grid[i] == 'P') grid[i] = 'P';
		REP(i, MAX_AREA) prev_grid[i] = grid[i];
		
		if (moves.S == 0) {
			cout << -1 << endl;
			continue;
		}
		
		for (PII &p : moves) cout << p2d(p.X).X << ' ' << p2d(p.X).Y << ' ' << dc[p.Y] << ' ';
		cout << endl;

	}
	return 0;
}
