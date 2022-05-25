// Author: Psyho
// Twitter: https://twitter.com/fakepsyho

// #define SILENT
 
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
#define VVVVI       VC<VVVI>
#define VPII        VC<PII>
#define VVPII       VC<VPII>
#define VVVPII      VC<VVPII>
#define VVVVPII     VC<VVVPII>
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
#ifdef SILENT
#define DB(...) ;
#else
#define DB(...) do {cerr << DB_(ARGS_SIZE(__VA_ARGS__), __VA_ARGS__) endl;} while (0)
#endif
 
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

double start_time = get_time();

const int MAX_N = 32;
const int MAX_C = 8;

int N;
int C;
int ogrid[MAX_N*MAX_N];
int grid[MAX_N*MAX_N];
int ins[MAX_N*MAX_N];
int D[MAX_C][MAX_C];

void load_data() {
	cin >> N;
	cin >> C;

	REP(r, N) REP(c, N) cin >> ogrid[r*32+c];
	REP(r, C) REP(c, C) cin >> D[r][c];
}

int dr[] = {-2,-2,-1,-1,1,1,2,2};
int dc[] = {-1,1,-2,2,-2,2,-1,1};
int dp[] = {-2*32-1,-2*32+1,-32-2,-32+2,32-2,32+2,2*32-1,2*32+1};

struct Move {
	int p1;
	int p2;
	
	Move() { }
	Move(int _p1, int _p2) {p1=_p1; p2=_p2;}
	
	string to_string() {
		return i2s(p1>>5) + " " + i2s(p1&31) + " " + i2s(p2>>5) + " " + i2s(p2&31);
	}
};

INLINE bool inside(int p) {
	return p >= 0 && ins[p];
}

int vs[MAX_N*MAX_N];
int vs_no = 1;

#define PX  pair<int, PII>
#define VPX VC<PX>
#define PY  pair<PII, PII>
#define VPY VC<PY>
VI    move_side1[MAX_N*MAX_N][MAX_C];
VPII  move_side2[MAX_N*MAX_N][MAX_C];
VPX   move_side3[MAX_N*MAX_N][MAX_C];
VPY   move_side4[MAX_N*MAX_N][MAX_C];

VI    gr_side1[MAX_N*MAX_N][MAX_C];
VPII  gr_side2[MAX_N*MAX_N][MAX_C];

#ifdef VM
const double TIME_LIMIT = 9.7 * 0.7;
#else
const double TIME_LIMIT = 9.7 * 10;
#endif

pair<int, VC<Move>> greedy(bool fast = false) {
	int total_score = 0;
	VC<Move> rv;
	int use_helpers = N < 8 ? rng.next(2) : 1;
	int use_helpers2 = N < 12 ? rng.next(2) : 1;
	double r_bonus = rng.next_double() * 0.2;
	double c_bonus = rng.next_double() * 0.2;
	
	const int N_TRIES = fast ? 4 : max(4, min(100, N*N*N/20));
	while (true) {
		if (get_time() - start_time > TIME_LIMIT) break;
		double bv = -1e9;
		VC<Move> bpath;
		VC<Move> path;
		int bscore = 0;
		int blen = 0;
		int total_calc = 0;
		
		REP(r, N) REP(c, N) {
			int p = r*32+c;
			REP(i, C) {
				gr_side1[p][i].clear();
				gr_side2[p][i].clear();
			}
		}
		
		REP(r, N) REP(c, N) {
			int p = r*32+c;
			if (grid[p] == -1) continue;
			REP(d, 8) {
				int p2 = p + dp[d];
				if (!inside(p2) || grid[p2] == -1) continue; 
				gr_side1[p][grid[p2]].PB(p2);
				REP(d2, 8) {
					int p3 = p2 + dp[d2];
					if (inside(p3) && grid[p3] != -1 && p != p3) gr_side2[p][grid[p3]].PB(MP(p2, p3));
				}
			}
		}
		
		REP(r1, N) REP(c1, N) {
			int p1 = r1*32+c1;
			if (grid[p1] < 0) continue;
			if (total_calc > 1000) {
				if (get_time() - start_time > TIME_LIMIT) {
					r1 = N;
					c1 = N;
					break;
				}
				total_calc = 0;
			}
			REP(d, 8) {
				REP(tries, N_TRIES) {
					int p2 = p1 + dp[d];
					if (!inside(p2) || grid[p2] < 0) continue;
					
					int t1 = grid[p1];
					int t2 = grid[p2];
					
					path.clear();
					path.PB(Move(p1, p2));
					int score = 0;
					int len = 1;
					if (t1 == t2) {
						score = 1;
						vs_no++;
						int np0 = p2;
						vs[p1] = vs_no;
						vs[p2] = vs_no;
						while (true) {
							int bp1 = -1;
							int bp2 = -1;
							int bp3 = -1;
							int bok_moves2 = 0;
							int ok_moves = 0;
							for (int np1 : gr_side1[np0][t1]) {
								if (vs[np1] == vs_no) continue;
								ok_moves++;
								int ok_moves2 = 0;
								for (int np2 : gr_side1[np1][t1]) {
									ok_moves2 += vs[np2] < vs_no;
								}
								if (bp1 < 0 || rng.next(ok_moves) == 0 || ok_moves2 > 0 && bok_moves2 == 0) {
									bp1 = np1;
									bok_moves2 = ok_moves2;
								}
							}
							if (bp1 < 0) {
								if (!use_helpers) break;
								int bcnt = 0;
								REP(nd1, 8) {
									int np1 = np0 + dp[nd1];
									if (!inside(np1) || grid[np1] == -1 || vs[np1] == vs_no) continue;
									int cnt = 0;
									int ap2 = -1;
									for (int np2 : gr_side1[np1][t1]) {
										if (vs[np2] < vs_no) {
											if (rng.next(++cnt) == 0) ap2 = np2;
										}
									}
									if (cnt > 0 && (cnt > bcnt || cnt == bcnt && rng.next(2) == 0)) {
										bp1 = np1;
										bp2 = ap2;
										bcnt = cnt;
									}
								}
							}
							if (bp1 < 0) {
								if (!use_helpers2) break;
								int bcnt = 0;
								REP(nd1, 8) {
									int np1 = np0 + dp[nd1];
									if (!inside(np1) || grid[np1] == -1 || vs[np1] == vs_no) continue;
									int cnt = 0;
									int ap2 = -1;
									int ap3 = -1;
									for (PII &p : gr_side2[np1][t1]) {
										if (vs[p.X] < vs_no && vs[p.Y] < vs_no) {
											if (rng.next(++cnt) == 0) {
												ap2 = p.X;
												ap3 = p.Y;
											}
										}
									}
									if (cnt > 0 && (cnt > bcnt || cnt == bcnt && rng.next(2) == 0)) {
										bp1 = np1;
										bp2 = ap2;
										bp3 = ap3;
										bcnt = cnt;
									}
								}
							}
							if (bp1 < 0) break;
							int np1 = bp1;
							if (bp2 >= 0) {
								int np2 = bp2;
								path.insert(path.begin(), Move(np2, np1));
								if (bp3 >= 0) {
									int np3 = bp3;
									assert(inside(np3) && grid[np3] == t1 && vs[np3] < vs_no);
									path.insert(path.begin(), Move(np3, np2));
									score += D[t1][grid[np2]];
									vs[np3] = vs_no;
								}
								assert(inside(np2) && grid[np2] != -1 && vs[np2] < vs_no);
								score += D[t1][grid[np1]];
								vs[np2] = vs_no;
							}
							vs[np1] = vs_no;
							path.PB(Move(np0, np1));
							np0 = np1;
							len++;
							score += len;
						}
					} else {
						score = D[t1][t2];
					}
					
					total_calc += bpath.S;
					double av = score;
					av += rng.next_double() * 0.1;
					av += r_bonus * r1 / N;
					av += c_bonus * c1 / N;
					if (av > bv) {
						bv = av;
						bpath = path;
						bscore = score;
						blen = len;
					}
					if (bpath.S == 1) break;
				}
			}
		}
		
		if (!bpath.S || bscore < 0) break;
		total_score += bscore;
		assert(bpath.S >= 1);
		for (Move &m : bpath) {
			grid[m.p2] = grid[m.p1];
			grid[m.p1] = -1;
			rv.PB(m);
		}
	}
	return MP(total_score, rv);
}

int sim(VC<Move> &moves, bool restore_grid = false) {
	int score = 0;
	int lp = -1;
	int combo = 0;
	for (Move &m : moves) {
		if (grid[m.p1] == -1 || grid[m.p2] == -1) return 0;
		if (grid[m.p1] == grid[m.p2]) {
			combo = m.p1 != lp ? 1 : combo + 1;
			lp = m.p2;
			score += combo;
		} else {
			combo = 0;
			lp = -1;
			score += D[grid[m.p1]][grid[m.p2]];
		}
		grid[m.p2] = grid[m.p1];
		grid[m.p1] = -1;
	}
	
	if (restore_grid) REP(r, N) REP(c, N) grid[r*32+c] = ogrid[r*32+c];
	return score;
}


struct SAState {
	VVI paths;
	VVI setups;
	int wasted;
	int xscore;
	double bv;
};



VC<Move> sa_initial_paths(double time_limit) {
	REP(r, N) REP(c, N) grid[r*32+c] = ogrid[r*32+c];
	
	VVI paths = VVI(C);
	VVI setups = VVI(MAX_N*N, VI());
	
	VI cache_sscore(C);
	VI cache_wasted(C);
	
	VVI bpaths;
	VVI bsetups;
	
	VI mark_path(MAX_N*N, 0);
	int mark_path_no = 1;
	
	double bv = -1e9;
	int bscore = 0;
	int xscore = 0;
	
	int gr_score = 0;
	VC<Move> gr_moves;
	
	ZERO(vs);
	
	REP(r, N) REP(c, N) {
		int p = r*32+c;
		REP(i, C) {
			move_side1[p][i].clear();
			move_side2[p][i].clear();
			move_side3[p][i].clear();
			move_side4[p][i].clear();
		}
	}
	DB(get_time() - start_time);
	REP(r, N) REP(c, N) {
		int p = r*32+c;
		if (grid[p] == -1) continue;
		REP(d, 8) {
			int p2 = p + dp[d];
			if (!inside(p2) || grid[p2] == -1) continue; 
			move_side1[p][grid[p2]].PB(p2);
			REP(d2, 8) {
				int p3 = p2 + dp[d2];
				if (!inside(p3) || grid[p3] == -1 || p == p3) continue;
				move_side2[p][grid[p3]].PB(MP(p2, p3));
				REP(d3, 8) {
					int p4 = p3 + dp[d3];
					if (!inside(p4) || grid[p4] == -1 || p2 == p4) continue;
					move_side3[p][grid[p4]].PB(MP(p2, MP(p3, p4)));
					REP(d4, 8) {
						int p5 = p4 + dp[d4];
						if (!inside(p5) || grid[p5] == -1 || p3 == p5 || p == p5) continue;
						move_side4[p][grid[p5]].PB(MP(MP(p2, p3), MP(p4, p5)));
					}
				}
			}
		}
	}
	DB(get_time() - start_time);
	
	int step = 0;
	int acc = 0;
	int acc_extend = 0;
	int found_extend = 0;
	int found_wrong = 0;
	
	auto path_score = [](int len) -> int {return len <= 1 ? 0 : len * (len - 1) / 2;};
		
	auto create_moves = [&](VVI &bp, VVI &bs) -> VC<Move> {
		VC<Move> rv;
		REP(color, C) {
			VI &path = bp[color];
			for (int p : path)
				REP(i, bs[p].S)
					rv.PB(Move(bs[p][i], i + 1 == bs[p].S ? p : bs[p][i+1]));
			REP(i, (int)path.S - 1) rv.PB(Move(path[i], path[i+1]));
		}
		return rv;
	};
	
	auto sim_and_update = [&](bool fast, VVI &bp, VVI &bs) -> void {
		VC<Move> moves = create_moves(bp, bs);
		int base_score = sim(moves);
		auto greedy_rv = greedy(fast);
			
		int real_score = base_score + greedy_rv.X;
		
		REP(r, N) REP(c, N) grid[r*32+c] = ogrid[r*32+c];
		
		if (real_score > gr_score) {
			gr_score = real_score;
			gr_moves = moves;
			gr_moves.insert(gr_moves.end(), ALL(greedy_rv.Y));
			DB(step, real_score);
		}
	};
		
	
	int ignore_color1 = -1;
	int ignore_color2 = -1;
	if (C >= 4) {
		VI color_cnts(C);
		REP(r, N) REP(c, N) color_cnts[grid[r*32+c]]++;
		ignore_color1 = 0;
		REP(i, C) if (color_cnts[i] < color_cnts[ignore_color1]) ignore_color1 = i;
		if (C >= 6) {
			ignore_color2 = ignore_color1 == 0;
			REP(i, C) if (i != ignore_color1 && color_cnts[i] < color_cnts[ignore_color2]) ignore_color2 = i;
		}
	}
	
	const int RESTARTS = N <= 10 ? 25 : N <= 25 ? 6 : 4;
	const double RESTART_AT = N <= 10 ? 0.20 : 0.30;
	int cur_restart = 0;
	VC<SAState> sa_states;
	
	bool extend_path = false;
	
	auto save_sa_state = [&]() -> SAState {
		SAState rv;
		rv.paths = bpaths;
		rv.setups = bsetups;
		rv.xscore = xscore;
		rv.bv = bv;
		rv.wasted = 0;
		REP(i, C) rv.wasted += cache_wasted[i];
		return rv;
	};
	
	double timer1 = 0;
	double timer2 = 0;
	double timer3 = 0;
	double timer4 = 0;
	double timer5 = 0;
	
	double global_time_passed;

	while (true) {
		step++;
		
		if (step % 100 == 1) {
			global_time_passed = (get_time() - start_time) / time_limit;
			if (global_time_passed >= 1.0) break;
		}
		
		double time_passed = global_time_passed * ((1.0 - RESTART_AT) + RESTARTS * RESTART_AT) - cur_restart * RESTART_AT;
		if (cur_restart == RESTARTS) time_passed += RESTART_AT;
		
		if (time_passed >= RESTART_AT && cur_restart < RESTARTS) {
			cur_restart++;
			if (RESTART_AT < 1.0) {
				sa_states.PB(save_sa_state());
				xscore = 0;
				REP(i, C) bpaths[i].clear(), bsetups[i].clear();
			}
			
			bv = 0;
			bscore = 0;
			REP(i, C) paths[i].clear();
			REP(r, N) REP(c, N) setups[r*32+c].clear();
			REP(i, C) cache_sscore[i] = 0;
			REP(i, C) cache_wasted[i] = 0;
			time_passed = 0;
			ZERO(vs);
			
			if (cur_restart == RESTARTS && RESTART_AT < 1.0) {
				int best = -1;
				double rbv = 0;
				REP(i, sa_states.S) {
					VC<Move> tm = create_moves(sa_states[i].paths, sa_states[i].setups);
					VI lens; REP(c, C) lens.PB(sa_states[i].paths[c].S);
					DB(i, sa_states[i].xscore, sim(tm, true), lens);
					double rav = sa_states[i].xscore;
					if (rav > rbv) {
						best = i;
						rbv = rav;
					}
				}
				
				if (best != -1) {
					SAState &s = sa_states[best];
					bpaths = s.paths;
					paths = s.paths;
					bsetups = s.setups;
					setups = s.setups;
					bscore = s.xscore;
					xscore = s.xscore;
					time_passed = RESTART_AT;
					int wasted = 0;
					REP(i, C) {
						for (int p : paths[i]) {
							if (!setups[p].S) continue;
							cache_wasted[i] += setups[p].S;
							FOR(j, 1, setups[p].S) cache_sscore[i] += D[i][grid[setups[p][j]]];
							cache_sscore[i] += D[i][grid[p]];
						}
						wasted += cache_wasted[i];
					}
					REP(i, C) for (int p : paths[i]) {
						vs[p] = vs_no;
						for (int x : setups[p]) vs[x] = vs_no;
					}
					bv = bscore * pow(time_passed, 0.04 * wasted / N);
				}
			}
		}
		
		double max_temp = 0.1 + N*N/C/2.7;
		double t = max_temp + (0.1 - max_temp) * time_passed * sqrt(time_passed);
		
		int color = rng.next(C);
		if ((color == ignore_color1 || color == ignore_color2) && time_passed < 0.95) continue;
		
		if (step % 2500 == 0 && (N < 20 || time_passed > 0.5)) {
			sim_and_update(true, bpaths, bsetups);
			vs_no++;
			REP(i, C) for (int p : paths[i]) {
				vs[p] = vs_no;
				for (int x : setups[p]) vs[x] = vs_no;
			}
		}
		
		if (rng.next(2)) reverse(ALL(paths[color]));
		
		if (rng.next(25) == 0) {
			for (int p : paths[color]) {
				if (grid[p] == color) continue;
				int max_ms = setups[p].S;
				for (int x : setups[p]) vs[x] = 0;
				if (setups[p].S) {
					cache_sscore[color] -= D[color][grid[p]];
					FOR(i, 1, setups[p].S) cache_sscore[color] -= D[color][grid[setups[p][i]]];
				}
				
				int cnt = 0;
				VI ns;
				for (int x : move_side1[p][color]) if (vs[x] < vs_no)
					if (rng.next(++cnt) == 0) ns = {x};
				if (ns.S == 0 && max_ms >= 2) {
					for (PII &x : move_side2[p][color]) if (grid[x.X] != color && vs[x.X] < vs_no && vs[x.Y] < vs_no)
					// for (PII &x : move_side2[p][color]) if (vs[x.X] < vs_no && vs[x.Y] < vs_no)
						if (rng.next(++cnt) == 0) ns = {x.Y, x.X};
				}
				if (ns.S == 0 && max_ms >= 3) {
					for (PX &x : move_side3[p][color]) if (grid[x.X] != color && grid[x.Y.X] != color && vs[x.X] < vs_no && vs[x.Y.X] < vs_no && vs[x.Y.Y] < vs_no)
					// for (PX &x : move_side3[p][color]) if (vs[x.X] < vs_no && vs[x.Y.X] < vs_no && vs[x.Y.Y] < vs_no)
						if (rng.next(++cnt) == 0) ns = {x.Y.Y, x.Y.X, x.X};
				}
				if (ns.S == 0) {
					DB(color, max_ms);
					REP(i, setups[p].S) DB(setups[p][i], grid[setups[p][i]]);
				}
				// assert(ns.S);
				// assert(grid[ns[0]] == color);
				cache_wasted[color] += ns.S - setups[p].S;
				for (int x : ns) vs[x] = vs_no;
				setups[p] = ns;
				if (setups[p].S) {
					cache_sscore[color] += D[color][grid[p]];
					FOR(i, 1, setups[p].S) cache_sscore[color] += D[color][grid[setups[p][i]]];
				}
			}
			continue;
		}
		
		
		int use_ms4 = rng.next(max(5, (int)(200 - time_passed * 500))) == 0;
		int use_ms3 = N <= 12 || use_ms4 || rng.next(max(1, (int)(20 - time_passed * 200))) == 0;
		
		int start = -1;
		static VI new_path;
		static VC<pair<int, VI>> new_setups;
		new_path.clear();
		new_setups.clear();
		
		int pos = 0;
		int nspos = 0;
		
		if (paths[color].S && rng.next(100) < (N <= 10 ? 7 : N <= 25 ? 5 : 3)) {
			//TODO: detect wrap around?
			extend_path = true;
			mark_path_no += N*N;
			REP(i, paths[color].S) mark_path[paths[color][i]] = mark_path_no + i;
			pos = rng.next(paths[color].S);
			assert(mark_path[paths[color][pos]] == mark_path_no + pos);
			
			start = paths[color][pos];
			int p = start;
			int end = -1;
			while (true) {
				if (p != start) {
					vs[p] = vs_no;
					new_path.PB(p);
				}
				int bp1 = -1;
				int bp2 = -1;
				int bp3 = -1;
				int bp4 = -1;
				
				int cnt = 0;
				int bok_moves = -1;
				
				if (bp1 < 0) {
					for (int np1 : move_side1[p][color]) {
						if (vs[np1] == vs_no) continue;
						int ok_moves = 0;
						for (int np2 : move_side1[np1][color]) ok_moves |= vs[np2] < vs_no;
						
						if (ok_moves > bok_moves) {
							cnt = 1;
							bok_moves = ok_moves;
							bp1 = np1;
						} else if (ok_moves == bok_moves) {
							if (rng.next(++cnt) == 0) bp1 = np1;
						}
					}
				}
				
				if (bp1 < 0) {
					for (PII &x : move_side2[p][color]) {
						if (vs[x.X] == vs_no || vs[x.Y] == vs_no) continue;
						if (rng.next(++cnt) == 0) bp1 = x.X, bp2 = x.Y;
					}
				}
				
				if (bp1 < 0 && use_ms3) {
					for (PX &x : move_side3[p][color]) {
						if (vs[x.X] == vs_no || vs[x.Y.X] == vs_no || vs[x.Y.Y] == vs_no) continue;
						if (rng.next(++cnt) == 0) bp1 = x.X, bp2 = x.Y.X, bp3 = x.Y.Y;
					}
				}
				
				if (bp1 < 0 && use_ms4) {
					for (PY &x : move_side4[p][color]) {
						if (vs[x.X.X] == vs_no || vs[x.X.Y] == vs_no || vs[x.Y.X] == vs_no || vs[x.Y.Y] == vs_no) continue;
						if (rng.next(++cnt) == 0) bp1 = x.X.X, bp2 = x.X.Y, bp3 = x.Y.X, bp4 = x.Y.Y;
					}
				}
				
				if (bp1 < 0) {
					for (int np1 : move_side1[p][color]) {
						if (mark_path[np1] < mark_path_no) continue;
						int diff = abs(pos - (mark_path[np1] - mark_path_no));
						if (diff <= 1 || diff > new_path.S+1) continue;
						if (rng.next(++cnt) == 0) {
							bp1 = np1;
							end = mark_path[np1] - mark_path_no;
						}
					}
				}
				
				if (bp1 < 0 || end != -1) break;
				
				p = bp1;
				if (bp2 >= 0) {
					static VI ns;
					ns.clear();
					if (bp3 >= 0) {
						if (bp4 >= 0) {
							ns.PB(bp4);
							vs[bp4] = vs_no;
						}
						ns.PB(bp3);
						vs[bp3] = vs_no;
					}
					ns.PB(bp2);
					vs[bp2] = vs_no;
					new_setups.PB(MP(p, ns));
				}
			}
			
			if (end == -1) {
				for (int p : new_path) vs[p] = 0;
				for (auto &p : new_setups) for (int x : p.Y) vs[x] = 0; 
				continue;
			}
			
			found_extend++;
			found_wrong += end < pos;
			if (end < pos) {
				swap(pos, end);
				reverse(ALL(new_path));
			}
			FOR(i, pos+1, end) {
				int p = paths[color][i];
				vs[p] = 0;
				for (int x : setups[p]) vs[x] = 0;
			}
			
			new_path.insert(new_path.begin(), paths[color].begin(), paths[color].begin() + pos+1);
			new_path.insert(new_path.end(), paths[color].begin() + end, paths[color].end());
			REP(i, pos+1) new_setups.PB(MP(paths[color][i], setups[paths[color][i]]));
			FOR(i, end, paths[color].S) new_setups.PB(MP(paths[color][i], setups[paths[color][i]]));
			
			pos = 0;
			nspos = 0;
		} else {
			extend_path = false;
			if (paths[color].S == 0 || rng.next(10) == 0) {
				for (int p : paths[color]) {
					vs[p] = 0;
					for (int x : setups[p]) vs[x] = 0;
				}
				int cnt = 0;
				REP(r, N) REP(c, N) if (grid[r*32+c] == color && vs[r*32+c] < vs_no && rng.next(++cnt) == 0) start = r*32+c;
				if (cnt == 0) {
					for (int p : paths[color]) {
						vs[p] = vs_no;
						for (int x : setups[p]) vs[x] = vs_no;
					}
					continue;
				}
			} else {
				pos = max(rng.next(paths[color].S), max(rng.next(paths[color].S), rng.next(paths[color].S)));
				REP(i, paths[color].S) {
					int p = paths[color][i];
					if (i<pos) {
						new_path.PB(p);
					} else {
						vs[p] = 0;
					}
					if (i<pos+1) {
						if (setups[p].S) {
							new_setups.PB(MP(p, setups[p]));
							nspos = new_setups.S;
						}
					} else {
						for (int x : setups[p]) vs[x] = 0;
					}
				}
				start = paths[color][pos];
			}
			
			if (rng.next(N <= 15 ? 10 : 50) == 0 && paths[color].S) {
				vs[start] = vs_no;
				new_path.PB(start);
				if (false && rng.next(10) == 0) {
					REP(i, pos+1) {
						int p = paths[color][i];
						vs[p] = 0;
						for (int x : setups[p]) vs[x] = 0;
					}
					pos = 0;
					new_path.clear();
					new_setups.clear();
				}
			} else {
				int p = start;
				while (true) {
					vs[p] = vs_no;
					new_path.PB(p);
					int bp1 = -1;
					int bp2 = -1;
					int bp3 = -1;
					int bp4 = -1;
					
					int cnt = 0;
					int bok_moves = -1;
					for (int np1 : move_side1[p][color]) {
						if (vs[np1] == vs_no) continue;
						int ok_moves = 0;
						for (int np2 : move_side1[np1][color]) ok_moves |= vs[np2] < vs_no;
						
						if (ok_moves > bok_moves) {
							cnt = 1;
							bok_moves = ok_moves;
							bp1 = np1;
						} else if (ok_moves == bok_moves) {
							if (rng.next(++cnt) == 0) bp1 = np1;
						}
					}
					
					if (bp1 < 0) {
						for (PII &x : move_side2[p][color]) {
							if (vs[x.X] == vs_no || vs[x.Y] == vs_no) continue;
							if (rng.next(++cnt) == 0) bp1 = x.X, bp2 = x.Y;
						}
					}
					
					if (bp1 < 0 && use_ms3) {
						for (PX &x : move_side3[p][color]) {
							if (vs[x.X] == vs_no || vs[x.Y.X] == vs_no || vs[x.Y.Y] == vs_no) continue;
							if (rng.next(++cnt) == 0) bp1 = x.X, bp2 = x.Y.X, bp3 = x.Y.Y;
						}
					}
					
					if (bp1 < 0 && use_ms4) {
						for (PY &x : move_side4[p][color]) {
							if (vs[x.X.X] == vs_no || vs[x.X.Y] == vs_no || vs[x.Y.X] == vs_no || vs[x.Y.Y] == vs_no) continue;
							if (rng.next(++cnt) == 0) bp1 = x.X.X, bp2 = x.X.Y, bp3 = x.Y.X, bp4 = x.Y.Y;
						}
					}
					
					if (bp1 < 0) break;
					
					p = bp1;
					if (bp2 >= 0) {
						static VI ns;
						ns.clear();
						if (bp3 >= 0) {
							if (bp4 >= 0) {
								ns.PB(bp4);
								vs[bp4] = vs_no;
							}
							ns.PB(bp3);
							vs[bp3] = vs_no;
						}
						ns.PB(bp2);
						vs[bp2] = vs_no;
						new_setups.PB(MP(p, ns));
					}
				}
			}
		}
		
		int sscore = 0;
		int wasted = 0;
		REP(i, C) if (i != color) {
			wasted += cache_wasted[i];
			sscore += cache_sscore[i];
		}
		
		int color_sscore = 0;
		int color_wasted = 0;
		for (auto &p : new_setups) {
			color_wasted += p.Y.S;
			FOR(i, 1, p.Y.S) color_sscore += D[color][grid[p.Y[i]]];
			color_sscore += D[color][grid[p.X]];
		}
		
		sscore += color_sscore;
		wasted += color_wasted;
		
		double av = 0;
		int score = 0;
		REP(i, C) if (i != color) score += path_score(paths[i].S);
		score += path_score(new_path.S);
		score += sscore;
		av = score;
		// av -= wasted * (wasted + 1) / 2 * (1 - time_passed) * (1 - time_passed) * (1 - time_passed) * 10;
		av *= pow(time_passed, 0.04 * wasted / N);
		// av += rng.next_double();
		
		
		if (av >= bv || rng.next_double() < exp((av - bv) / t)) {
			acc++;
			acc_extend += extend_path;
			bv = av;
			bscore = score;
			
			for (auto &p : paths[color]) setups[p].clear();
			paths[color] = new_path;
			for (auto &x : new_setups) setups[x.X] = x.Y;
			
			cache_sscore[color] = color_sscore;
			cache_wasted[color] = color_wasted;
			
			if (score > xscore) {
				VI lens;
				REP(i, paths.S) lens.PB(paths[i].S);
				bpaths = paths;
				bsetups = setups;
				DB(step, score, av, time_passed, t, lens, wasted, mark_path_no);
				xscore = score;
				
			}
		} else {
			FOR(i, pos, new_path.S) vs[new_path[i]] = 0;
			FOR(i, nspos, new_setups.S) for (int x : new_setups[i].Y) vs[x] = 0; 
			
			FOR(i, pos, paths[color].S) {
				int p = paths[color][i];
				vs[p] = vs_no;
				for (int x : setups[p])
					vs[x] = vs_no;
			}
		}
		
	}
	
	cerr << "step=" << step << endl;
	// DB(step);
	DB(acc);
	DB(acc_extend);
	DB(found_extend);
	DB(found_wrong);
	
	sim_and_update(false, bpaths, bsetups);
	
	DB(gr_moves.S);
	return gr_moves;
}

int main() {
	load_data();
	
	REP(r, N) REP(c, N) ins[r*32+c] = 1;
	
	int max_score = 0;
	int max_chain = 0;
	VI color_cnt(C);
	REP(i, C) {
		int cnt = 0;
		REP(r, N) REP(c, N) cnt += ogrid[r*32+c] == i;
		max_score += cnt * (cnt + 1) / 2;
		max_chain = max(max_chain, cnt);
		color_cnt[i] = cnt;
	}
	
	cerr << "[DATA] N = " << N << endl;
	cerr << "[DATA] C = " << C << endl;
	cerr << "[DATA] MS = " << max_score << endl;
	cerr << "[DATA] MC = " << max_chain << endl;
	DB(color_cnt);
	
	VC<Move> bsol = sa_initial_paths(TIME_LIMIT);
	int score = sim(bsol);
	DB(score);
	
	cout << bsol.S << endl;
	REP(i, bsol.S) cout << bsol[i].to_string() << endl;
	cout.flush();
}