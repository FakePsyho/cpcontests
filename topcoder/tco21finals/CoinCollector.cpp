// Author: Psyho
// Twitter: https://twitter.com/fakepsyho

const double TIME_LIMIT = 9.0;

// #define USE_TIMERS

#ifdef USE_TIMERS
#define TIME(x) x
#else
#define TIME(x) ;
#endif
	
 
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


char dirs[] = {'L','U','R','D'};    
int dc[] = {-1,  0,+1, 0};
int dr[] = { 0, -1, 0,+1};
int N;
int D;
int cr;
int cc;

map<string, VI> moves;
VI dmoves[12];
int dmoves_no[12];
VS rmoves;

VS move_names;
map<string, int> dice2id;

VS cd;
VS grid;


const int MAX_STATES = 12376;

map<VI, int> dices2state[7];
VVI state2dices[7];
int state_trans[MAX_STATES][12];
int state_trans_mid[MAX_STATES][12];

int get_dices2state(VS sdice) {
	VI v;
	for (string &s : sdice) v.PB(dice2id[s]);
	sort(ALL(v));
	return dices2state[v.S][v];
}

int get_dices2state(VI v) {
	sort(ALL(v));
	return dices2state[v.S][v];
}


const int SAFETY_STEPS = 24;
float safety[32][32][MAX_STATES];
float safety_mid[32][32][MAX_STATES];
float safety_mid_dice[32][32][MAX_STATES][12];

int dist(int a, int b) {
	return min(abs(a-b),min(abs(a-b+N),abs(a-b-N)));
}

void get_move_pre() {
	for (auto p : moves) 
		REP(i, p.Y.S) rmoves.PB(p.X);
}


string get_rmove() {
	return rmoves[rng.next(rmoves.S)];
}

double last_bv = 0;

int moves_value[] = {0, 1, 0, 60, 0, 0, 105};

double PENALTY_SCALE = 1.0;

PII choose_greedy(VS &sgrid, VS &sdice, int sscore, int sr, int sc, bool sim = false) {
	REP(dice, D) if (sdice[dice] == "STAY" && sgrid[sr][sc] != 'S') {
		last_bv = 1e12;
		return MP(dice, 0);
	}
	
	double bv = -1;
	int bdir = -1;
	int bdice = -1;
	
	
	REP(dice, D) {
		bool repeat = false;
		REP(i, dice) if (sdice[dice] == sdice[i]) repeat = true;
		if (repeat) continue;
		
		int dup = 0;
		REP(i, D) if (dice != i && sdice[dice] == sdice[i]) dup++;
		
		VS ndice = sdice;
		ndice.erase(ndice.begin()+dice);
		int st = get_dices2state(ndice);
		
		REP(dir, 4) {
			int av = 0;
			int lowest_xv = 1e9;
			for (int m : moves[sdice[dice]]) {
				int xv = 0;
				int nr = (sr + dr[dir] * m + N) % N;
				int nc = (sc + dc[dir] * m + N) % N;
				if (sgrid[nr][nc] == '.') {
					xv += sscore;
				} else if (sgrid[nr][nc] == 'S') {
					xv += sscore / 2;
				} else if (sgrid[nr][nc] == 'C') {
					xv += sscore + 100;
				}
				
				int coins = 0;
				int close_coins = 0;
				if (!sim) {
					REP(r, N) if (r != nr) {
						if (sgrid[r][nc] == 'C') {
							if (dist(r, nr) <= 6) close_coins++;
							coins++;
						}
						
					}
					REP(c, N) if (c != nc) {
						if (sgrid[nr][c] == 'C') {
							if (dist(c, nc) <= 6) close_coins++;
							coins++;
						}
					}
				}
				
				int spike_penalty = (xv/2) * (1 - safety_mid[nr][nc][st] * safety_mid[nr][nc][st]);
				
				xv -= min(xv / 2, (int)(spike_penalty));
				xv += min(5, coins) * 6;
				xv += min(1, coins) * 10;
				
				av += xv;
				lowest_xv = min(lowest_xv, xv);
			}
			
			av /= moves[sdice[dice]].S;
			
			av += moves_value[moves[sdice[dice]].S];
			
			if (av > bv) {
				bv = av;
				bdir = dir;
				bdice = dice;
			}
		}
	}
	
	last_bv = bv;
	return MP(bdice, bdir);
}

void apply_move(VS &sgrid, VS &sdice, int &sscore, int &sr, int &sc, int dice, int dir, int value = -1) {
	int m = value == -1 ? moves[sdice[dice]][rng.next(moves[sdice[dice]].S)] : value;
	sr = (sr + dr[dir] * m + N) % N;
	sc = (sc + dc[dir] * m + N) % N;
	
	if (sgrid[sr][sc] == 'C') {
		sgrid[sr][sc] = '.';
		sscore += 100;
	} else if (sgrid[sr][sc] == 'S') {
		sscore /= 2;
	}
	
	sdice[dice] = get_rmove();
}

const int TURN_PENALTY = 0;

VD sim_scores;

void sim_scores_clear() {
	sim_scores.clear();
}

double sim_scores_get_avg_x() {
	double sum = 0;
	REP(i, sim_scores.S) sum += sim_scores[i];
	return sum / sim_scores.S;
}

double sim(VS xgrid, VS xdice, int xscore, int xr, int xc, int max_steps, int sims = 1) {
	double total = 0;
	REP(csim, sims) {
		VS sgrid = xgrid;
		VS sdice = xdice;
		int sscore = xscore;
		int evscore = xscore;
		int sr = xr;
		int sc = xc;
		
		int max_score = sscore;
		
		REP(turn, max_steps) {
			PII move = choose_greedy(sgrid, sdice, sscore, sr, sc, true);
			
			int new_evscore = 0;
			for (int m : moves[sdice[move.X]]) {
				int nr = (sr + m * dr[move.Y] + N) % N;
				int nc = (sc + m * dc[move.Y] + N) % N;
				if (sgrid[nr][nc] == '.') {
					new_evscore += evscore;
				} else if (sgrid[nr][nc] == 'C') {
					new_evscore += evscore + 100;
				} else if (sgrid[nr][nc] == 'S') {
					new_evscore += evscore / 2;
				}
			}
			evscore = new_evscore / moves[sdice[move.X]].S;
			
			apply_move(sgrid, sdice, sscore, sr, sc, move.X, move.Y);
			max_score = max(max_score, evscore - turn * TURN_PENALTY);
			if (evscore + (max_steps - turn) * (100 - TURN_PENALTY) < max_score) break;
		}
		
		sim_scores.PB(max_score);
		total += max_score;
	}
	return total * 1.0 / sims;
}


int main() {
	moves["STAY"] = {0};
	moves["ONE"] = {1};
	moves["TWO"] = {2};
	moves["THREE"] = {3};
	moves["FOUR"] = {4};
	moves["FIVE"] = {5};
	moves["SIX"] = {6};
	moves["EVEN"] = {2,4,6};
	moves["ODD"] = {1,3,5};
	moves["LOW"] = {1,2,3};
	moves["HIGH"] = {4,5,6};
	moves["RANDOM"] = {1,2,3,4,5,6};
	get_move_pre();
	
	cin >> N >> D >> cr >> cc;
	cd = VS(D);
	REP(i, D) cin >> cd[i];
	
	DB(N, D);
	
	grid = VS(N, string(N, 'X'));
	REP(r, N) REP(c, N) {
		string s;
		cin >> s;
		grid[r][c] = s[0];
	}
	
	int coins_left = 0;
	REP(r, N) REP(c, N) coins_left += grid[r][c] == 'C';
	
	for (auto &p : moves) move_names.PB(p.X);
	for (string &s : move_names) {
		int id = dice2id.S;
		dice2id[s] = id;
		dmoves[id] = moves[s];
		dmoves_no[id] = dmoves[id].S;
	}
	
	//gen state2dices & dices2state
	FOR(a, 0, 12) FOR(b, a, 12) FOR(c, b, 12) FOR(d, c, 12) FOR(e, d, 12) FOR(f, e, 12) {
		VI v = {a,b,c,d,e,f};
		int state_no = state2dices[v.S].S;
		state2dices[v.S].PB(v);
		dices2state[v.S][v] = state_no;
	}
	FOR(a, 0, 12) FOR(b, a, 12) FOR(c, b, 12) FOR(d, c, 12) FOR(e, d, 12) {
		VI v = {a,b,c,d,e};
		int state_no = state2dices[v.S].S;
		state2dices[v.S].PB(v);
		dices2state[v.S][v] = state_no;
	}
	FOR(a, 0, 12) FOR(b, a, 12) FOR(c, b, 12) FOR(d, c, 12) {
		VI v = {a,b,c,d};
		int state_no = state2dices[v.S].S;
		state2dices[v.S].PB(v);
		dices2state[v.S][v] = state_no;
	}
	FOR(a, 0, 12) FOR(b, a, 12) FOR(c, b, 12) {
		VI v = {a,b,c};
		int state_no = state2dices[v.S].S;
		state2dices[v.S].PB(v);
		dices2state[v.S][v] = state_no;
	}
	FOR(a, 0, 12) FOR(b, a, 12) {
		VI v = {a,b};
		int state_no = state2dices[v.S].S;
		state2dices[v.S].PB(v);
		dices2state[v.S][v] = state_no;
	}
	FOR(a, 0, 12) {
		VI v = {a};
		int state_no = state2dices[v.S].S;
		state2dices[v.S].PB(v);
		dices2state[v.S][v] = state_no;
	}
	
	//gen state_trans
	MINUS(state_trans);
	REP(state, state2dices[D].S) {
		REP(i, D) {
			VI nv = state2dices[D][state];
			nv.erase(nv.begin()+i);
			state_trans[state][state2dices[D][state][i]] = get_dices2state(nv);
		}
	}
	REP(state, state2dices[D-1].S) {
		REP(i, 12) {
			VI nv = state2dices[D-1][state];
			nv.PB(i);
			state_trans_mid[state][i] = get_dices2state(nv);
		}
	}
	
	DB(get_time() - start_time);
	
	//calc safety
	int states_no = state2dices[D].S;
	DB(states_no);
	int states_mid_no = state2dices[D-1].S;
	DB(states_mid_no);
	
	REP(r, N) REP(c, N) REP(s, states_no) safety[r][c][s] = (grid[r][c] == 'S' ? 0.0 : 1.0);
	
	FOR(step, 1, SAFETY_STEPS) {
		REP(r, N) REP(c, N) REP(s, states_mid_no) {
			safety_mid[r][c][s] = 0;
			if (grid[r][c] == 'S') continue;
			REP(dice, 12) safety_mid[r][c][s] += safety[r][c][state_trans_mid[s][dice]] * dmoves_no[dice];
			safety_mid[r][c][s] /= 25;
		}			
		DB(step, get_time() - start_time);
		
		if (step == SAFETY_STEPS-1) break;
		if (get_time() - start_time > 5.0) break;
			
		REP(r, N) REP(c, N) REP(s, states_mid_no) REP(dice, 12) {
			safety_mid_dice[r][c][s][dice] = 0;
			if (grid[r][c] == 'S') continue;
			REP(dir, 4) {
				float p = 0.0;
				for (int m : dmoves[dice]) {
					int nr = (r + dr[dir] * m + N) % N;
					int nc = (c + dc[dir] * m + N) % N;
					p += safety_mid[nr][nc][s];
				}
				p /= dmoves[dice].S;
				safety_mid_dice[r][c][s][dice] = max(safety_mid_dice[r][c][s][dice], p);
			}
		}
		
		REP(r, N) REP(c, N) REP(s, states_no) {
			safety[r][c][s] = 0;
			if (grid[r][c] == 'S') continue;
			REP(dice, 12) if (state_trans[s][dice] != -1) {
				int st = state_trans[s][dice];
				safety[r][c][s] = max(safety[r][c][s], safety_mid_dice[r][c][st][dice]);
			}
		}
		
		DB(step, get_time() - start_time);
	}
	
	int score = 0;
	int last_time_elapsed = 0;
	
	REP(turn, N*N) {
		
		
		PII rv = choose_greedy(grid, cd, score, cr, cc);
		double bv = last_bv;
		int bdice = rv.X;
		int bdir = rv.Y;
		
		int evscore = 0;
		sim_scores_clear();
		for (int m : moves[cd[bdice]]) {
			VS sgrid = grid;
			VS sdice = cd;
			int sscore = score;
			int sr = cr;
			int sc = cc;
			apply_move(sgrid, sdice, sscore, sr, sc, bdice, bdir, m);
			sim(sgrid, sdice, sscore, sr, sc, min(10 , N*N-turn-1), (N < 15 ? 36 : N < 25 ? 24 : 12) / moves[cd[bdice]].S);
		}
		evscore = sim_scores_get_avg_x();
		
		if (evscore < score) {
			DB(evscore);
			
			evscore = 0;
			sim_scores_clear();
			for (int m : moves[cd[bdice]]) {
				VS sgrid = grid;
				VS sdice = cd;
				int sscore = score;
				int sr = cr;
				int sc = cc;
				apply_move(sgrid, sdice, sscore, sr, sc, bdice, bdir, m);
				sim(sgrid, sdice, sscore, sr, sc, min(75 , N*N-turn-1), (N < 15 ? 36 : N < 25 ? 24 : 18) / moves[cd[bdice]].S);
			}
			evscore = sim_scores_get_avg_x();
				
			if (evscore < score) {
				DB(evscore);
				
				evscore = 0;
				sim_scores_clear();
				for (int m : moves[cd[bdice]]) {
					VS sgrid = grid;
					VS sdice = cd;
					int sscore = score;
					int sr = cr;
					int sc = cc;
					apply_move(sgrid, sdice, sscore, sr, sc, bdice, bdir, m);
					sim(sgrid, sdice, sscore, sr, sc, min(400 , N*N-turn-1), 36 / moves[cd[bdice]].S);
				}
				evscore = sim_scores_get_avg_x();
			}
		}
		
		DB(turn, score, bv, evscore, last_time_elapsed);
		
		if (coins_left == 0 || evscore < score || last_time_elapsed > 9000) {
			cout << "X" << endl;
			exit(0);
		}
		
		
		cout << dirs[bdir] << endl;
		cout << cd[bdice] << endl;
		cout.flush();    

		int dice_value;
		double time_elapsed;
		
		cin >> time_elapsed;
		cin >> dice_value;
		cin >> cd[bdice];
		
		last_time_elapsed = time_elapsed;
		
		int old_score = score;
		cr = (cr + dr[bdir] * dice_value + N) % N;
		cc = (cc + dc[bdir] * dice_value + N) % N;
		if (grid[cr][cc] == 'C') {
			coins_left--;
			grid[cr][cc] = '.';
			score += 100;
		} else if (grid[cr][cc] == 'S') {
			score /= 2;
		}
	}
}
