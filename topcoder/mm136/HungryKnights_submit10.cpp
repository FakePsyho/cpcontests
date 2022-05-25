// Author: Psyho
// Twitter: https://twitter.com/fakepsyho

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

const double TIME_LIMIT = 9.6;

pair<int, VC<Move>> greedy() {
	int total_score = 0;
	VC<Move> rv;
	int use_helpers = N < 8 ? rng.next(2) : 1;
	int use_helpers2 = N < 12 ? rng.next(2) : 1;
	double move_penalty = rng.next_double() * 2.0;
	double r_bonus = rng.next_double() * 0.2;
	double c_bonus = rng.next_double() * 0.2;
	int sacrifice_color1 = rng.next(-1, C);
	int sacrifice_color2 = -1;
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
				move_side1[p][i].clear();
				move_side2[p][i].clear();
			}
		}
		
		REP(r, N) REP(c, N) {
			int p = r*32+c;
			if (grid[p] == -1) continue;
			REP(d, 8) {
				int p2 = p + dp[d];
				if (!inside(p2) || grid[p2] == -1) continue; 
				move_side1[p][grid[p2]].PB(d);
				REP(d2, 8) {
					int p3 = p2 + dp[d2];
					if (inside(p3) && grid[p3] != -1 && p != p3) move_side2[p][grid[p3]].PB(MP(d, d2));
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
				REP(tries, max(4, min(100, N*N*N/20))) {
					int p2 = p1 + dp[d];
					if (!inside(p2) || grid[p2] < 0) continue;
					
					int t1 = grid[p1];
					int t2 = grid[p2];
					
					path.clear();
					path.PB(Move(p1, p2));
					int score = 0;
					int len = 1;
					int wasted = 0;
					if (t1 == t2) {
						score = 1;
						vs_no++;
						int np0 = p2;
						vs[p1] = vs_no;
						vs[p2] = vs_no;
						while (true) {
							int bd1 = -1;
							int bd2 = -1;
							int bd3 = -1;
							int bok_moves2 = 0;
							int ok_moves = 0;
							for (int nd1 : move_side1[np0][t1]) {
								int np1 = np0 + dp[nd1];
								if (vs[np1] == vs_no) continue;
								ok_moves++;
								int ok_moves2 = 0;
								for (int nd2 : move_side1[np1][t1]) {
									ok_moves2 += vs[np1 + dp[nd2]] < vs_no;
								}
								if (bd1 < 0 || rng.next(ok_moves) == 0 || ok_moves2 > 0 && bok_moves2 == 0) {
									bd1 = nd1;
									bok_moves2 = ok_moves2;
								}
							}
							if (bd1 < 0) {
								if (!use_helpers) break;
								int bcnt = 0;
								REP(nd1, 8) {
									int np1 = np0 + dp[nd1];
									if (!inside(np1) || grid[np1] == -1 || vs[np1] == vs_no) continue;
									int cnt = 0;
									int ad2 = -1;
									for (int nd2 : move_side1[np1][t1]) {
										int np2 = np1 + dp[nd2];
										if (vs[np2] < vs_no) {
											int v = 1;
											v += grid[np1] == sacrifice_color1;
											cnt += v;
											if (rng.next(cnt) < v) ad2 = nd2;
										}
									}
									if (cnt > 0 && (cnt > bcnt || cnt == bcnt && rng.next(2) == 0)) {
										bd1 = nd1;
										bd2 = ad2;
										bcnt = cnt;
									}
								}
							}
							if (bd1 < 0) {
								if (!use_helpers2) break;
								int bcnt = 0;
								REP(nd1, 8) {
									int np1 = np0 + dp[nd1];
									if (!inside(np1) || grid[np1] == -1 || vs[np1] == vs_no) continue;
									int cnt = 0;
									int ad2 = -1;
									int ad3 = -1;
									for (PII &p : move_side2[np1][t1]) {
										if (vs[np1 + dp[p.X]] < vs_no && vs[np1 + dp[p.X] + dp[p.Y]] < vs_no) {
											int v = 1;
											v += grid[np1] == sacrifice_color1;
											v += grid[np1 + dp[p.X]] == sacrifice_color1;
											cnt += v;
											if (rng.next(cnt) < v) {
												ad2 = p.X;
												ad3 = p.Y;
											}
										}
									}
									if (cnt > 0 && (cnt > bcnt || cnt == bcnt && rng.next(2) == 0)) {
										bd1 = nd1;
										bd2 = ad2;
										bd3 = ad3;
										bcnt = cnt;
									}
								}
							}
							if (bd1 < 0) break;
							int np1 = np0 + dp[bd1];
							if (bd2 >= 0) {
								int np2 = np1 + dp[bd2];
								path.insert(path.begin(), Move(np2, np1));
								wasted += 1;//grid[np1] != sacrifice_color;
								if (bd3 >= 0) {
									int np3 = np2 + dp[bd3];
									assert(bd3 >= 0 && bd3 < 8 && inside(np3) && grid[np3] == t1 && vs[np3] < vs_no);
									path.insert(path.begin(), Move(np3, np2));
									score += D[t1][grid[np2]];
									vs[np3] = vs_no;
									wasted += 1;//grid[np2] != sacrifice_color;
								}
								assert(bd2 >= 0 && bd2 < 8 && inside(np2) && grid[np2] != -1 && vs[np2] < vs_no);
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
					av -= wasted * (wasted + 1) / 2 * move_penalty;
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
		// if (bpath.S > 1) {
			// int color = grid[bpath.back().p2];
			// int color_left = 0;
			// REP(r, N) REP(c, N) color_left += grid[r*32+c] == color;
			// DB(blen, bpath.S, color, blen+color_left, bscore);
		// }
	}
	return MP(total_score, rv);
}

VC<Move> sa_initial_paths(double time_limit) {
	REP(r, N) REP(c, N) grid[r*32+c] = ogrid[r*32+c];
	
	VVI paths = VVI(C);
	VVI setups = VVI(MAX_N*N, VI());
	
	VI cache_sscore(C);
	VI cache_wasted(C);
	
	VVI bpaths;
	VVI bsetups;
	
	double bv = -1e9;
	int bscore = 0;
	int xscore = 0;
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
	
	auto path_score = [](int len) -> int {return len <= 1 ? 0 : len * (len - 1) / 2;};
		
	int step = 0;
	while (true) {
		step++;
		double time_passed = (get_time() - start_time) / time_limit;
		if (time_passed >= 1.0) break;
		
		// double max_temp = max(10.0, bscore / 10.0);
		// double t = max_temp * pow(0.1 / max_temp, time_passed);
		// double max_temp = min(200.0, N * 7.0) / 10.0;
		double max_temp = 100.0;
		double t = max_temp + (0.1 - max_temp) * time_passed;
		
		int color = rng.next(C);
		int use_ms3 = rng.next(3)  == 0;
		int use_ms4 = rng.next(10) == 0;
		
		int start = -1;
		static VI new_path;
		static VC<pair<int, VI>> new_setups;
		new_path.clear();
		new_setups.clear();
		
		if (rng.next(2)) reverse(ALL(paths[color]));
		vs_no++;
		
		REP(i, C) if (i != color) for (int p : paths[i]) {
			vs[p] = vs_no;
			for (int x : setups[p]) vs[x] = vs_no;
		}
		
		if (paths[color].S == 0 || rng.next(10) == 0) {
			int cnt = 0;
			REP(r, N) REP(c, N) if (grid[r*32+c] == color && vs[r*32+c] < vs_no && rng.next(++cnt) == 0) start = r*32+c;
			if (cnt == 0) continue;
		} else {
			int pos = rng.next(paths[color].S);
			REP(i, pos+1) {
				int p = paths[color][i];
				if (i<pos) {
					new_path.PB(p);
					vs[p] = vs_no;
				}
				if (setups[p].S) {
					for (int x : setups[p]) vs[x] = vs_no;
					new_setups.PB(MP(p, setups[p]));
				}
			}
			start = paths[color][pos];
		}
		
		if (false && rng.next(50) == 0 && paths[color].S) {
			vs[start] = vs_no;
			new_path.PB(start);
			if (rng.next(10) == 0) {
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
				int bok_moves = 0;
				for (int np1 : move_side1[p][color]) {
					if (vs[np1] == vs_no) continue;
					int ok_moves = 0;
					for (int np2 : move_side1[np1][color]) ok_moves |= vs[np2] < vs_no;
					
					if (ok_moves > bok_moves) {
						cnt = 0;
						bok_moves = ok_moves;
					}
					if (ok_moves == bok_moves) {
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
		
		// if (new_path.S == 1) continue;
		
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
		av *= pow(time_passed, 2.0 * wasted / N / N);
		// av += rng.next_double();
		
		if (av >= bv || rng.next_double() < exp((av - bv) / t)) {
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
				DB(step, score, av, time_passed, t, lens, wasted);
				xscore = score;
			}
		} else {
			// pass
		}
		
	}
	
	DB(step);
	
	VC<Move> rv;
	REP(color, C) {
		VI &path = bpaths[color];
		for (int p : path)
			REP(i, bsetups[p].S)
				rv.PB(Move(bsetups[p][i], i + 1 == bsetups[p].S ? p : bsetups[p][i+1]));
		REP(i, (int)path.S - 1) rv.PB(Move(path[i], path[i+1]));
	}
	
	DB(rv.S);
	return rv;
}

void sim(VC<Move> &moves) {
	int score = 0;
	int lp = -1;
	int combo = 0;
	for (Move &m : moves) {
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
	DB(score);
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
	
	VC<Move> bsol = sa_initial_paths(TIME_LIMIT - 0.1);
	sim(bsol);
	auto greedy_rv = greedy();
	DB(greedy_rv.X);
	DB(greedy_rv.Y.S);
	bsol.insert(bsol.end(), ALL(greedy_rv.Y));
	
	cout << bsol.S << endl;
	REP(i, bsol.S) cout << bsol[i].to_string() << endl;
	cout.flush();
}