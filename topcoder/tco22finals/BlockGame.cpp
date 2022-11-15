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
 
double PI = 2*acos(0);

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

// CODE

const int MAX_N = 20;
const int MAX_T = 5;
const int MAX_C = 5;

int N, T, C;
double F, P;
int turn;
int score;
int elapsed_time;

int grid[MAX_N][MAX_N];
int tiles[MAX_T][3][3];
int rc[MAX_N][MAX_C+1];
int cc[MAX_N][MAX_C+1];

const double TIME_LIMIT = 10000;

int counter_add_tile = 0;
int counter_tile_ok = 0;
int counter_score_board = 0;
int counter_filled_lines = 0;


bool tile_ok(int id, int r0, int c0) {
	counter_tile_ok++;
	REP(r, 3) REP(c, 3) {
		int nr = r0 + r;
		int nc = c0 + c;
		if (tiles[id][r][c] && (nr < 0 || nr >= N || nc < 0 || nc >= N || grid[nr][nc])) return false;
	}
	return true;
}

bool qtile_ok(int id, int r0, int c0) {
	counter_tile_ok++;
	REP(r, 3) REP(c, 3) {
		int nr = r0 + r;
		int nc = c0 + c;
		if (tiles[id][r][c] && grid[nr][nc]) return false;
	}
	return true;
}

INLINE void add_cell(int r, int c, int color) {
	grid[r][c] = color;
	rc[r][color]++;
	rc[r][0]--;
	cc[c][color]++;
	cc[c][0]--;
}

INLINE void rem_cell(int r, int c) {
	if (!grid[r][c]) return;
	rc[r][grid[r][c]]--;
	rc[r][0]++;
	cc[c][grid[r][c]]--;
	cc[c][0]++;
	grid[r][c] = 0;
}

INLINE void add_tile(int id, int r0, int c0) {
	counter_add_tile++;
	REP(r, 3) REP(c, 3) if (tiles[id][r][c]) add_cell(r0+r, c0+c, tiles[id][r][c]);
}

INLINE void rem_tile(int id, int r0, int c0) {
	REP(r, 3) REP(c, 3) if (tiles[id][r][c]) rem_cell(r0+r, c0+c);
}

void init_grid() {
	REP(r, N) REP(c, N) {
		rc[r][grid[r][c]]++;
		cc[c][grid[r][c]]++;
	}
}

int last_lines;
VD row_score(21);
VD hole_penalty(21);

VD holes_penalty = {0.0, 3.5, 3.0, 2.0, 2.0, 1.5};

double score_board() {
	counter_score_board++;
	
	double v = 0;
	
	REP(i, N) {
		int mx;
		mx = 0; FOR(j, 1, C+1) mx = max(mx, rc[i][j]); v += row_score[mx]; 
		mx = 0; FOR(j, 1, C+1) mx = max(mx, cc[i][j]); v += row_score[mx]; 
	}
	
	double holes_cnt = 0;
	REP(r, N) REP(c, N) if (!grid[r][c]) {
		if ((c == 0 || grid[r][c-1]) && (c == N-1 || grid[r][c+1])) holes_cnt += hole_penalty[c];
		if ((r == 0 || grid[r-1][c]) && (r == N-1 || grid[r+1][c])) holes_cnt += hole_penalty[r];
	}
	v -= N * holes_penalty[C] * holes_cnt;
	
	double bonus_score = 0;
	
	int lines = 0;
	REP(i, N) {
		if (rc[i][0] == 0) {
			lines++;
			int mx = 0; FOR(j, 1, C+1) mx = max(mx, rc[i][j]);
			bonus_score += mx * mx;
		}
		if (cc[i][0] == 0) {
			lines++;
			int mx = 0; FOR(j, 1, C+1) mx = max(mx, cc[i][j]);
			bonus_score += mx * mx;
		}
	}
	last_lines = lines;
	
	v += bonus_score * N*N*0.2;
	if (lines == 1) v *= 0.1;
	if (lines > 1) v *= lines;
	
	return v;
}

bool filled_lines() {
	counter_filled_lines++;
	REP(i, N) {
		if (rc[i][0] == 0) return true;
		if (cc[i][0] == 0) return true;
	}
	return false;
}

int count_filled_lines() {
	counter_filled_lines++;
	int rv = 0;
	REP(i, N) {
		if (rc[i][0] == 0) rv++;
		if (cc[i][0] == 0) rv++;
	}
	return rv;
}

double edge_penalty(int r, int c) {
	return -(max(abs(r+1 - N*.5 - .5), abs(c+1 - N*.5 - .5))) * 1.4;
}

double tile_penalty(int id) {
	double v = 0;
	REP(r, 3) REP(c, 3) v += tiles[id][r][c] > 0;
	return -(v*v) * N;
}

LL count_states(VVPII &vv, int level) {
	LL rv = 0;
	level = min(level, T);
	if (level == 0) {
		return 1;
	} else if (level == 1) {
		REP(i, T) rv += vv[i].S;
	} else if (level == 2) {
		REP(i, T) REP(j, i) rv += (LL)vv[i].S * vv[j].S * 2;
	} else if (level == 3) {
		REP(i, T) REP(j, i) REP(k, j) rv += (LL)vv[i].S * vv[j].S * vv[k].S * 6;
	} else if (level == 4) {
		REP(i, T) REP(j, i) REP(k, j) REP(l, k) rv += (LL)vv[i].S * vv[j].S * vv[k].S * vv[l].S * 24LL;
	}
	return rv;
}

int main() {
	cin >> N >> T >> C >> F >> P;
	
	cerr << "[DATA] N = " << N << endl;
	cerr << "[DATA] T = " << T << endl;
	cerr << "[DATA] C = " << C << endl;
	cerr << "[DATA] F = " << F << endl;
	cerr << "[DATA] P = " << P << endl;
  
	//read grid
	REP(r, N) REP(c, N) cin >> grid[r][c];
	init_grid();
	
	
	//data
	VI clears(7);
	int dropped_tiles = 0;
	int dropped_cells = 0;
	int total_cells = 0;
	
	int last_t2_turn = 0;
	int last_t3_turn = 0;
	int last_t4_turn = 0;
	
	int xxx1 = 0;
	int xxx2 = 0;
	int xxx3 = 0;
	int xxx4 = 0;
	int xxx5 = 0;
	
	//read tiles
	REP(i, T) {
		string line;
		cin >> line;
		REP(r, 3) REP(c, 3) tiles[i][r][c] = (int)(line[r*3+c]-'0');
	}
	
	REP(i, 21) row_score[i] = pow(i, 3.8 - P * 3.5 + C * 0.05 + max(0.0, 0.6 - 2 * P) + (P > 0.45 ? -0.1 : 0) + (C == 1 ? -0.1 : 0));
	REP(i, N) hole_penalty[i] = 1.0;
	hole_penalty[0] = hole_penalty[N-1] = 2.0;
	hole_penalty[1] = hole_penalty[N-2] = 1.25;
	double move_penalty = pow(N, 2.25);
	
	turn = 0;
	
	VVPII tile_pos(T);
	
	const double MAX_STATES_L3 = 2e6;
	const double MAX_STATES_L4 = 2e6;
	
	while (true) {
		turn++;
		if (turn > 1000) break;
		
		double bv = -1e9;
		int bt = -1, br = -100, bc = -100;
		REP(i, T) {
			double av = -1e5;
			REP(r, 3) REP(c, 3) av += tiles[i][r][c] > 0;
			if (av > bv) {
				bv = av;
				bt = i;
				br = -100;
				bc = -100;
			}
		}
		REP(i, T) {
			tile_pos[i].clear();
			FOR(r, -2, N) FOR(c, -2, N) if (tile_ok(i, r, c)) tile_pos[i].PB(MP(r, c));
		}
		
		VC<LL> cst; FOR(level, 0, 5) cst.PB(count_states(tile_pos, level));
		// DB(turn, cst);
		
		REP(i, T) for (auto &t1 : tile_pos[i]) {
			// if (!tile_ok(i, t1.X, t1.Y)) continue;
			add_tile(i, t1.X, t1.Y);
			double av = score_board();
			xxx1++;
			av -= edge_penalty(t1.X, t1.Y);
			av -= tile_penalty(i);
			if (av > bv) {
				bv = av;
				bt = i;
				br = t1.X;
				bc = t1.Y;
			}
			rem_tile(i, t1.X, t1.Y);
		}
		
		if (T > 1 && elapsed_time < TIME_LIMIT * 0.9) {
			last_t2_turn = turn;
			REP(i, T) for (auto &t1 : tile_pos[i]) {
				// if (!tile_ok(i, t1.X, t1.Y)) continue;
				add_tile(i, t1.X, t1.Y);
				if (!filled_lines()) {
					// REP(j, T) if (i != j) FOR(t2.X, -2, N) FOR(t2.Y, -2, N) {
					REP(j, T) if (i != j) for (auto &t2 : tile_pos[j]) {
						if (!qtile_ok(j, t2.X, t2.Y)) continue;
						// if (abs(t2.X - t1.X) > 2 && abs(t2.Y - t1.Y) > 2) continue;
						add_tile(j, t2.X, t2.Y);
						if (count_filled_lines() >= 2) {
							xxx2++;
							double av = score_board() - move_penalty;
							av -= edge_penalty(t1.X, t1.Y);
							av -= tile_penalty(i);
							if (av > bv) {
								bv = av;
								bt = i;
								br = t1.X;
								bc = t1.Y;
							}
						}
						rem_tile(j, t2.X, t2.Y);
					}
				}
				rem_tile(i, t1.X, t1.Y);
			}
		}
		
		if (T > 2 && elapsed_time < TIME_LIMIT * 0.7 && cst[3] < MAX_STATES_L3) {
			last_t3_turn = turn;
			REP(i, T) for (auto &t1 : tile_pos[i]) {
				// if (!tile_ok(i, t1.X, t1.Y)) continue;
				add_tile(i, t1.X, t1.Y);
				if (!filled_lines()) {
					REP(j, T) if (i != j) for (auto &t2 : tile_pos[j]) {
						// if (abs(t2.X - t1.X) > 2 && abs(t2.Y - t1.Y) > 2) continue;
						if (!qtile_ok(j, t2.X, t2.Y)) continue;
						add_tile(j, t2.X, t2.Y);
						if (!filled_lines()) {
							REP(k, T) if (i != k && j != k) for (auto &t3 : tile_pos[k]) {
								// if (abs(t3.X - t1.X) > 2 && abs(t3.Y - t1.Y) > 2) continue;
								if (!qtile_ok(k, t3.X, t3.Y)) continue;
								add_tile(k, t3.X, t3.Y);
								if (count_filled_lines() >= 2) {
									xxx3++;
									double av = score_board() - move_penalty * 2;
									av -= edge_penalty(t1.X, t1.Y);
									av -= tile_penalty(i);
									if (av > bv) {
										bv = av;
										bt = i;
										br = t1.X;
										bc = t1.Y;
									}
								}
								rem_tile(k, t3.X, t3.Y);
							}
						}
						rem_tile(j, t2.X, t2.Y);
					}
				}
				rem_tile(i, t1.X, t1.Y);
			}
		}
		
		if (T > 3 && elapsed_time < TIME_LIMIT * 0.5 && cst[4] < MAX_STATES_L4) {
			last_t4_turn = turn;
			REP(i, T) for (auto &t1 : tile_pos[i]) {
				// if (!tile_ok(i, t1.X, t1.Y)) continue;
				add_tile(i, t1.X, t1.Y);
				counter_add_tile--;
				if (!filled_lines()) {
					REP(j, T) if (i != j) for (auto &t2 : tile_pos[j]) {
						if (abs(t2.X - t1.X) > 2 && abs(t2.Y - t1.Y) > 2) continue;
						if (!qtile_ok(j, t2.X, t2.Y)) continue;
						add_tile(j, t2.X, t2.Y);
						if (!filled_lines()) {
							REP(k, T) if (i != k && j != k) for (auto &t3 : tile_pos[k]) {
								if (abs(t3.X - t1.X) > 2 && abs(t3.Y - t1.Y) > 2) continue;
								if (!qtile_ok(k, t3.X, t3.Y)) continue;
								add_tile(k, t3.X, t3.Y);
								if (!filled_lines()) {
									REP(l, T) if (i != l && j != l && k != l) for (auto &t4 : tile_pos[l]) {
										if (abs(t4.X - t1.X) > 2 && abs(t4.Y - t1.Y) > 2) continue;
										if (!qtile_ok(l, t4.X, t4.Y)) continue;
										add_tile(l, t4.X, t4.Y);
										if (count_filled_lines() >= 2) {
											xxx4++;
											double av = score_board() - move_penalty * 3;
											av -= edge_penalty(t1.X, t1.Y);
											av -= tile_penalty(i);
											if (av > bv) {
												bv = av;
												bt = i;
												br = t1.X;
												bc = t1.Y;
											}
										}
										rem_tile(l, t4.X, t4.Y);
									}
								}
								rem_tile(k, t3.X, t3.Y);
							}
						}
						rem_tile(j, t2.X, t2.Y);
					}
				}
				rem_tile(i, t1.X, t1.Y);
			}
		}
		
		assert(bt != -1);
		int used_tile = bt;
		if (br == -100) {
			cout << bt << endl;
		} else {
			cout << bt << ' ' << br << ' ' << bc << endl;
		}
		
		REP(r, 3) REP(c, 3) if (tiles[bt][r][c]) total_cells++;
		if (br != -100) {
			add_tile(bt, br, bc);
			
			VI rows, cols;
			
			FOR(r, br, br+3) if (r >= 0 && r < N) {
				int cnt = 0; REP(c, N) cnt += grid[r][c] > 0;
				if (cnt == N) rows.PB(r);
			}
			
			FOR(c, bc, bc+3) if (c >= 0 && c < N) {
				int cnt = 0; REP(r, N) cnt += grid[r][c] > 0;
				if (cnt == N) cols.PB(c);
			}
			
			int total = rows.S + cols.S;
			clears[total]++;
			
			for (int r : rows) {
				VI cnt(C);
				REP(c, N) cnt[grid[r][c]-1]++;
				int mx = 0; REP(i, C) mx = max(mx, cnt[i]);
				score += total * mx * mx;
			}
			
			for (int c : cols) {
				VI cnt(C);
				REP(r, N) cnt[grid[r][c]-1]++;
				int mx = 0; REP(i, C) mx = max(mx, cnt[i]);
				score += total * mx * mx;
			}
			
			for (int r : rows) REP(c, N) rem_cell(r, c);
			for (int c : cols) REP(r, N) rem_cell(r, c);
		} else {
			dropped_tiles++;
			REP(r, 3) REP(c, 3) if (tiles[bt][r][c]) dropped_cells++;
		}
		
		
		if (turn == 1000) {
			cerr << "[DATA] dt = " << dropped_tiles << endl;
			cerr << "[DATA] dc = " << 1.0 * dropped_cells / total_cells << endl;
			cerr << "[DATA] m1 = " << clears[1] << endl;
			cerr << "[DATA] m2 = " << clears[2] << endl;
			cerr << "[DATA] m3 = " << clears[3] << endl;
			cerr << "[DATA] m4 = " << clears[4] << endl;
			cerr << "[DATA] m5 = " << clears[5] << endl;
			cerr << "[DATA] m6 = " << clears[6] << endl;
			cerr << "[DATA] t2 = " << last_t2_turn << endl;
			cerr << "[DATA] t3 = " << last_t3_turn << endl;
			cerr << "[DATA] t4 = " << last_t4_turn << endl;
			cerr.flush();
		#ifdef VM
			this_thread::sleep_for(chrono::milliseconds(1000));
		#endif
		}
		
		string line;
		cin >> line;
		REP(r, 3) REP(c, 3) tiles[used_tile][r][c] = (int)(line[r*3+c]-'0');
		
		cin >> elapsed_time;
	}

	DB(score);
	DB(counter_tile_ok);
	DB(counter_add_tile);
	DB(counter_filled_lines);
	DB(counter_score_board);
	DB(elapsed_time);
	if (xxx1) DB(xxx1);
	if (xxx2) DB(xxx2);
	if (xxx3) DB(xxx3);
	if (xxx4) DB(xxx4);
	if (xxx5) DB(xxx5);

	//terminate
	cout << "-1" << endl;
	cout.flush();
  
	return 0;
}
