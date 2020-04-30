// Author: Psyho
// Twitter: https://twitter.com/fakepsyho
// Site: http://psyho.gg/
// Contest: https://www.topcoder.com/challenges/30122730

#define NDEBUG
 
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
#define R           first
#define C           second
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

#define ARGS_SIZE_(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,size,...) size
#define ARGS_SIZE(...) ARGS_SIZE_(__VA_ARGS__,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)

#define DB_1(x) #x << "=" << (x) <<
#define DB_2(x, ...) #x << "=" << (x) << ", " << DB_1(__VA_ARGS__)
#define DB_3(x, ...) #x << "=" << (x) << ", " << DB_2(__VA_ARGS__)
#define DB_4(x, ...) #x << "=" << (x) << ", " << DB_3(__VA_ARGS__)
#define DB_5(x, ...) #x << "=" << (x) << ", " << DB_4(__VA_ARGS__)
#define DB_6(x, ...) #x << "=" << (x) << ", " << DB_5(__VA_ARGS__)
#define DB_7(x, ...) #x << "=" << (x) << ", " << DB_6(__VA_ARGS__)
#define DB_8(x, ...) #x << "=" << (x) << ", " << DB_7(__VA_ARGS__)
#define DB_9(x, ...) #x << "=" << (x) << ", " << DB_8(__VA_ARGS__)
#define DB_10(x, ...) #x << "=" << (x) << ", " << DB_9(__VA_ARGS__)
#define DB_11(x, ...) #x << "=" << (x) << ", " << DB_10(__VA_ARGS__)
#define DB_12(x, ...) #x << "=" << (x) << ", " << DB_11(__VA_ARGS__)
#define DB_13(x, ...) #x << "=" << (x) << ", " << DB_12(__VA_ARGS__)
#define DB_14(x, ...) #x << "=" << (x) << ", " << DB_13(__VA_ARGS__)
#define DB_15(x, ...) #x << "=" << (x) << ", " << DB_14(__VA_ARGS__)
#define DB__(size, ...) DB_##size(__VA_ARGS__)
#define DB_(size, ...) DB__(size, __VA_ARGS__)

#ifdef NDEBUG
#define DB(...) ;
#else
#define DB(...) do {cerr << DB_(ARGS_SIZE(__VA_ARGS__), __VA_ARGS__) endl;} while (0)
#endif 

double get_time() {timeval tv; gettimeofday(&tv, NULL); return tv.tv_sec + tv.tv_usec * 1e-6;}
 
struct RNG {
    unsigned int MT[624];
    int index;
   
    RNG(int seed = 1) {
        init(seed);
    }
   
    void init(int seed = 1) {
        MT[0] = seed;
        FOR(i, 1, 624) MT[i] = (1812433253UL * (MT[i-1] ^ (MT[i-1] >> 30)) + i);
        index = 0;
    }
   
    void generate() {
        const unsigned int MULT[] = {0, 2567483615UL};
        REP(i, 227) {
            unsigned int y = (MT[i] & 0x8000000UL) + (MT[i+1] & 0x7FFFFFFFUL);
            MT[i] = MT[i+397] ^ (y >> 1);
            MT[i] ^= MULT[y&1];
        }
        FOR(i, 227, 623) {
            unsigned int y = (MT[i] & 0x8000000UL) + (MT[i+1] & 0x7FFFFFFFUL);
            MT[i] = MT[i-227] ^ (y >> 1);
            MT[i] ^= MULT[y&1];
        }
        unsigned int y = (MT[623] & 0x8000000UL) + (MT[0] & 0x7FFFFFFFUL);
        MT[623] = MT[623-227] ^ (y >> 1);
        MT[623] ^= MULT[y&1];
    }
   
    unsigned int rand() {
        if (index == 0) {
            generate();
        }
       
        unsigned int y = MT[index];
        y ^= y >> 11;
        y ^= y << 7  & 2636928640UL;
        y ^= y << 15 & 4022730752UL;
        y ^= y >> 18;
        index = index == 623 ? 0 : index + 1;
        return y;
    }
   
    INLINE int next() {
        return rand();
    }
   
    INLINE int next(int x) {
        return rand() % x;
    }
   
    INLINE int next(int a, int b) {
        return a + (rand() % (b - a));
    }
   
    INLINE double next_double() {
        return (rand() + 0.5) * (1.0 / 4294967296.0);
    }
   
    INLINE double next_double(double a, double b) {
        return a + next_double() * (b - a);
    }
};
 
static RNG rng;

double start_time = get_time();


const int MAX_N = 42;

int N, P;

struct Move {
	int r0;
	int c0;
	int d;
	bool right;
	
	Move() { }
	
	Move(int _r0, int _c0, int _d, bool _right) {
		r0 = _r0;
		c0 = _c0;
		d = _d;
		right = _right;
	}
};

int om[MAX_N][MAX_N];
int m[MAX_N][MAX_N];

int mcost[MAX_N];


//Board Operations

int sm[MAX_N][MAX_N];
void save_board() {
	REP(i, N) REP(j, N) sm[i][j] = m[i][j];
}

void restore_board() {
	REP(i, N) REP(j, N) m[i][j] = sm[i][j];
}

int validate_vs[42<<8];
int validate_vs_no = 0;
void validate_board(int r1, int c1) {
	validate_vs_no++;
	REP(r, r1) REP(c, c1) {
		assert(validate_vs[m[r][c]] != validate_vs_no);
		validate_vs[m[r][c]] = validate_vs_no;
	}
	REP(r, r1) REP(c, c1) assert(validate_vs[(r<<8)+c] == validate_vs_no);
}

void crop_board(int r0, int c0, int r1, int c1, bool rescale = true) {
	if (rescale) {
		REP(r, r1-r0) REP(c, c1-c0) m[r][c] = m[r+r0][c+c0] - (r0<<8) - c0;
	} else {
		REP(r, r1-r0) REP(c, c1-c0) m[r][c] = m[r+r0][c+c0];
	}
}

void reset() {
	REP(i, N) REP(j, N) m[i][j] = om[i][j];
}

int mcopy[MAX_N][MAX_N];
void rotate(int r0, int c0, int d, bool right) {
	assert(r0+d<=N && c0+d<=N);
	FOR(r,r0,r0+d) FOR(c,c0,c0+d) mcopy[r][c] = m[r][c];
	if (right) {
		FOR(r,r0,r0+d) FOR(c,c0,c0+d) m[c+r0-c0][-r+r0+c0+d-1] = mcopy[r][c];
	} else {
		FOR(r,r0,r0+d) FOR(c,c0,c0+d) m[r][c] = mcopy[c+r0-c0][-r+r0+c0+d-1];
	}
}

VI apply_move(VI &state, Move &m, int r1, int c1) {
	VI rv = state;
	if (m.right) {
		FOR(r, m.r0, m.r0+m.d) FOR(c, m.c0, m.c0+m.d) rv[(c+m.r0-m.c0)*c1+(-r+m.r0+m.c0+m.d-1)] = state[r*c1+c];
	} else {
		FOR(r, m.r0, m.r0+m.d) FOR(c, m.c0, m.c0+m.d) rv[r*c1+c] = state[(c+m.r0-m.c0)*c1+(-r+m.r0+m.c0+m.d-1)];
	}
	return rv;	
}

void show(int r1=-1, int c1=-1) {
	if (r1 < 0) r1 = N;
	if (c1 < 0) c1 = N;
	REP(r, r1) {
		REP(c, c1) cerr << ((m[r][c]>>8)*N)+(m[r][c]&255)+1 << "\t";
		cerr << endl;
	}
	cerr << endl;
}



// Move Operations

int calc_movescore(VC<Move> &moves) {
	int rv = 0;
	for (Move &m : moves) rv += mcost[m.d];
	return rv;
}

void shift_moves(VC<Move> &moves, int r0, int c0) {
	for (Move &m : moves) m.r0 += r0, m.c0 += c0;
}

void exec_moves(VC<Move> &moves) {
	for (Move &m : moves) rotate(m.r0, m.c0, m.d, m.right);
}



// Scoring Operations

PII score_split = MP(-1, -1);

const int MODE_NORMAL = 0;
const int MODE_RSPLIT = 1;
const int MODE_CSPLIT = 2;
const int MODE_QSPLIT = 3;

const int SCALE_NORMAL = 0;
const int SCALE_MCOST  = 1;

const int SECOND_NONE       = 0;
const int SECOND_NORMAL     = 1;
const int SECOND_RSPLIT_BOT = 2;
const int SECOND_RSPLIT_TOP = 3;

int lu_score[8][64];
int lu_qscore[8][64][64];
int xcost[128];

template <int mode=MODE_NORMAL, int scale=SCALE_NORMAL, int second=SECOND_NONE>
INLINE int score(int x, int r, int c) {
	int rv;
	
	if (mode == MODE_NORMAL) {
		int tr = x>>8;
		int tc = x&255;
		rv = abs(r-tr)+abs(c-tc);
	}
	
	if (mode == MODE_RSPLIT) rv = lu_score[x][r];
	if (mode == MODE_CSPLIT) rv = lu_score[x][c];
	if (mode == MODE_QSPLIT) rv = lu_qscore[x][r][c];

	if (scale == SCALE_MCOST) rv = xcost[rv];
	
	return rv;
}

int calc_score() {
	int rv = 0;
	REP(r, N) REP(c, N) rv += score(m[r][c], r, c);
	return rv;
}


template <int mode=MODE_NORMAL, int scale=SCALE_NORMAL, int second=SECOND_NONE>
int calc_score(VI &state, int r1, int c1) {
	int rv = 0;
	REP(r, r1) REP(c, c1) rv += score<mode, scale, second>(state[r*c1+c], r, c);
	return rv;
}




int prefixsum[MAX_N][MAX_N];
INLINE int score_square(int r0, int c0, int d) {
	return prefixsum[r0+d][c0+d] - prefixsum[r0+d][c0] - prefixsum[r0][c0+d] + prefixsum[r0][c0];
}

template <int mode, int scale, int second>
void calc_prefixsum(VI &state, int r1, int c1) {
	REP(r, r1) REP(c, c1) prefixsum[r+1][c+1] = score<mode, scale, second>(state[r*c1+c], r, c);
	REP(r, r1) REP(c, c1) prefixsum[r+1][c+1] += prefixsum[r][c+1];
	REP(r, r1) REP(c, c1) prefixsum[r+1][c+1] += prefixsum[r+1][c];
}

template <int mode, int scale, int second>
int eval_rotate(VI &state, int r1, int c1, int r0, int c0, int d, bool right) {
	assert(r0+d<=r1 && c0+d<=c1);
	int rotate_score = 0;
	if (right) {
		FOR(r,r0,r0+d) FOR(c,c0,c0+d) rotate_score += score<mode, scale, second>(state[r*c1+c], c+r0-c0, -r+r0+c0+d-1);
	} else {
		FOR(r,r0,r0+d) FOR(c,c0,c0+d) rotate_score += score<mode, scale, second>(state[(c+r0-c0)*c1+(-r+r0+c0+d-1)], r, c);
	}
	return rotate_score - score_square(r0, c0, d);
}


// Solve: Grid Search 
int total_bs = 0;

const int MAX_HASH = 9999889;
const int HASH_MULT = 79;
bitset<MAX_HASH> hashes;

LL hash_mult[42*42];

INLINE int calc_hash(VI &v) {
	// int hash = 0;
	// REP(i, v.S) hash = (hash * HASH_MULT + v[i]) & (MAX_HASH-1);
	// return hash;
	
	int hash = 0;
	REP(i, v.S) hash = (hash * HASH_MULT + v[i]) % MAX_HASH;
	return hash;
	
	// LL hash = 0;
	// REP(i, v.S) hash += v[i] * hash_mult[i];
	// return int(hash % MAX_HASH);
}

template <int mode, int scale=SCALE_NORMAL, int second=SECOND_NONE>


VC<Move> solve_bs(int r1, int c1, int max_states, int max_d=100) {
	hashes.reset();
	
	double bs_start_time = get_time();
	
	const int MAX_MOVES = 2048;
	VC<VC<pair<int,pair<int,Move>>>> qmove(MAX_MOVES);
	VVVI qstate(MAX_MOVES);
	
	VI ostate; REP(r, r1) REP(c, c1) ostate.PB(m[r][c]);
	
	qmove[0].PB(MP(calc_score<mode, scale, second>(ostate, r1, c1), MP(0, Move())));
	qstate[0].PB(ostate);
	
	int states_vis = 0;
	int collisions = 0;
	
	int longest_move = mcost[min(max_d, min(r1, c1))] + 10;
	
	REP(mscore, MAX_MOVES) {
		
		// if (mscore % 10 == 0) DB(mscore, qmove[mscore][0].X, qmove[mscore].S);
		// DB(mscore, qmove[mscore].S, qmove[mscore][0].X);
		int sorted_no = 0;
		int checked = 0;
		REP(i, (int)qmove[mscore].S) {
			if (i >= sorted_no) {
				int prev_sorted_no = sorted_no;
				if (sorted_no + 5*max_states < qmove[mscore].S) {
					nth_element(qmove[mscore].begin() + sorted_no, qmove[mscore].begin() + sorted_no + 5*max_states, qmove[mscore].end(), [](pair<int,pair<int,Move>> &a, pair<int,pair<int,Move>> &b){return a.X < b.X;});
					sorted_no += 5*max_states;
				} else {
					sorted_no = qmove[mscore].S;
				}
				sort(qmove[mscore].begin() + prev_sorted_no, qmove[mscore].begin() + sorted_no, [](pair<int,pair<int,Move>> &a, pair<int,pair<int,Move>> &b){return a.X < b.X;});
			}
			
			if (checked >= max_states) break;
			
			int calc = qmove[mscore][i].X;
			int prev = qmove[mscore][i].Y.X;
			Move &m = qmove[mscore][i].Y.Y;
			
			if (calc <= (second == SECOND_NONE ? 0 : 0xFFFF) || mscore == MAX_MOVES-1) {
				total_bs += mscore;
				DB("", mscore, states_vis, collisions, get_time()-bs_start_time);
				VC<Move> rv;
				int pcost = mscore;
				int ppos = i;
				while (pcost) {
					rv.PB(qmove[pcost][ppos].Y.Y);
					ppos = qmove[pcost][ppos].Y.X;
					pcost -= mcost[rv.back().d];
				}
				reverse(ALL(rv));
				
				return rv;
			}
			
			VI state;
			if (mscore == 0) {
				state = qstate[0][0];
			} else {
				state = apply_move(qstate[mscore-mcost[m.d]][prev], m, r1, c1);
			}
			
			qstate[mscore].PB(state);
			
			states_vis++;
			int h = calc_hash(state);
			if (hashes.test(h)) {
				collisions++;
				continue;
			} 
			
			checked++;
			
			hashes.set(h);
			
			calc_prefixsum<mode, scale, second>(state, r1, c1);
			
			FOR(d, 2, min(max_d+1, min(r1,c1)+1)) {
				REP(r, r1-d+1) REP(c, c1-d+1) {
					qmove[mscore+mcost[d]].PB(MP(prefixsum[r1][c1]+eval_rotate<mode, scale, second>(state, r1, c1, r, c, d, true ), MP(i, Move(r,c,d,true ))));
					qmove[mscore+mcost[d]].PB(MP(prefixsum[r1][c1]+eval_rotate<mode, scale, second>(state, r1, c1, r, c, d, false), MP(i, Move(r,c,d,false))));
				}
			}
		}
	}
	
	return VC<Move>();
}

VC<Move> allmoves;

bool convert_use_map = false;
PII convert_split = MP(-1, -1);
int cmap[MAX_N][MAX_N];

template <int mode, int scale=SCALE_NORMAL, int second=SECOND_NONE>
void xsolve_bs(int r0, int c0, int r1, int c1, int max_states, int max_d=100) {
	DB("", "xsolve_bs", r0, c0, r1, c1, max_states, max_d, mode, scale, convert_use_map, convert_split, score_split);
	
	save_board();
	
	if (convert_use_map) {
		assert(mode == MODE_RSPLIT || mode == MODE_CSPLIT);
		REP(r, N) REP(c, N) m[r][c] = cmap[r][c];
		FOR(r, r0, r1) FOR(c, c0, c1) assert(m[r][c] == 4 || m[r][c] == 5);
	} else {
		if (convert_split.R != -1 && convert_split.C != -1) {
			assert(mode == MODE_QSPLIT);
			REP(r, N) REP(c, N) m[r][c] = 4 + ((m[r][c]>>8)>=convert_split.R) + ((m[r][c]&255)>=convert_split.C)*2;
		} else if (convert_split.R != -1) {
			assert(mode == MODE_RSPLIT);
			REP(r, N) REP(c, N) m[r][c] = 4 + ((m[r][c]>>8)>=convert_split.R);
		} else if (convert_split.C != -1) {
			assert(mode == MODE_CSPLIT);
			REP(r, N) REP(c, N) m[r][c] = 4 + ((m[r][c]&255)>=convert_split.C);
		} else {
			assert(mode == MODE_NORMAL);
		}
	}
	
	crop_board(r0, c0, r1, c1, mode == MODE_NORMAL);
	
	PII copy_score_split = score_split;
	score_split.R -= r0;
	score_split.C -= c0;
	
	
	if (mode == MODE_RSPLIT) {
		FOR(x, 4, 6) REP(r, N) lu_score[x][r] = ((r < score_split.R) != (x == 4)) ? abs(r - score_split.R) + (r >= score_split.R) : 0;
	} else if (mode == MODE_CSPLIT) {
		FOR(x, 4, 6) REP(c, N) lu_score[x][c] = ((c < score_split.C) != (x == 4)) ? abs(c - score_split.C) + (c >= score_split.C) : 0;
 	} else if (mode == MODE_QSPLIT) {
		FOR(x, 4, 8) REP(r, N) REP(c, N) lu_qscore[x][r][c] = 
			(((r < score_split.R) != !(x&1)) ? abs(r - score_split.R) + (r >= score_split.R) : 0) + 
			(((c < score_split.C) != !(x&2)) ? abs(c - score_split.C) + (c >= score_split.C) : 0);
	}
	
	VC<Move> moves = solve_bs<mode, scale, second>(r1-r0, c1-c0, max_states, max_d);
	
	score_split = copy_score_split;
	
	restore_board();
	shift_moves(moves, r0, c0);
	exec_moves(moves);
	for (Move &m : moves) allmoves.PB(m);
}


int WIDTH_SPLIT;
int WIDTH_ENABLE;
int WIDTH_FULL_SMALL;
int WIDTH_FULL_QSPLIT;


void solve_full(int, int, int, int, double);


int msplit[MAX_N][MAX_N];
void enable_split(int r0, int c0, int r1, int c1, int rsplit, int csplit, int no4, int no5) {
	DB("enable_split", r0, c0, r1, c1, rsplit, csplit, no4, no5);
	
	assert(rsplit == -1 || csplit == -1);
	assert(rsplit != -1 || csplit != -1);
	
	assert(convert_split.R == -1 || convert_split.C == -1);
	assert(convert_split.R != -1 || convert_split.C != -1);
	
	if (convert_split.R != -1) {
		REP(r, N) REP(c, N) msplit[r][c] = 4 + ((m[r][c]>>8)>=convert_split.R);
	} else {
		REP(r, N) REP(c, N) msplit[r][c] = 4 + ((m[r][c]&255)>=convert_split.C);
	}
	
	convert_use_map = true;
	FOR(r, r0, r1) FOR(c, c0, c1) cmap[r][c] = 5;
	if (rsplit != -1) {
		assert(no4+no5 == (rsplit-r0)*(c1-c0));
		int no4found = 0, no5found = 0;
		FOR(r, r0, rsplit) FOR(c, c0, c1) {
			if (msplit[r][c] == 4) no4found++; else no5found++;
			cmap[r][c] = 4;
		}
		
		int cmidl = max(c0, (c0+c1)/2-6);
		int cmidr = min(c1, (c0+c1)/2+6);
		
		for (int r = rsplit-1; r >= r0; r--) FOR(c, cmidl, cmidr) {
			if (no4found <= no4 && no5found <= no5) break;
			if (msplit[r][c] == 4 && no4found > no4) {
				cmap[r][c] = 5;
				no4found--;
			} else if (msplit[r][c] == 5 && no5found > no5) {
				cmap[r][c] = 5;
				no5found--;
			}
		}
		FOR(r, rsplit, r1) FOR(c, cmidl, cmidr) {
			if (no4found >= no4 && no5found >= no5) break;
			if (msplit[r][c] == 4 && no4found < no4) {
				cmap[r][c] = 4;
				no4found++;
			} else if (msplit[r][c] == 5 && no5found < no5) {
				cmap[r][c] = 4;
				no5found++;
			} 
		}
		int start_row = r0, end_row = r1;
		while (true) {
			bool ok = true; FOR(c, c0, c1) if (cmap[start_row][c] == 5) ok = false;
			if (!ok) break;
			start_row++;
		}
		while (true) {
			bool ok = true; FOR(c, c0, c1) if (cmap[end_row-1][c] == 4) ok = false;
			if (!ok) break;
			end_row--;
		}
		score_split = MP(rsplit, -1);
		xsolve_bs<MODE_RSPLIT>(start_row, cmidl, end_row, cmidr, WIDTH_ENABLE, 3);
	} else {
		assert(no4+no5 == (csplit-c0)*(r1-r0));
		int no4found = 0, no5found = 0;
		FOR(r, r0, r1) FOR(c, c0, csplit) {
			if (msplit[r][c] == 4) no4found++; else no5found++;
			cmap[r][c] = 4;
		}
		
		int rmidt = max(r0, (r0+r1)/2-6);
		int rmidb = min(r1, (r0+r1)/2+6);
		
		for (int c = csplit-1; c >= c0; c--) FOR(r, rmidt, rmidb) {
			if (no4found <= no4 && no5found <= no5) break;
			if (msplit[r][c] == 4 && no4found > no4) {
				cmap[r][c] = 5;
				no4found--;
			} else if (msplit[r][c] == 5 && no5found > no5) {
				cmap[r][c] = 5;
				no5found--;
			}
		}
		FOR(c, csplit, c1) FOR(r, rmidt, rmidb) {
			if (no4found >= no4 && no5found >= no5) break;
			if (msplit[r][c] == 4 && no4found < no4) {
				cmap[r][c] = 4;
				no4found++;
			} else if (msplit[r][c] == 5 && no5found < no5) {
				cmap[r][c] = 4;
				no5found++;
			} 
		}
		
		int start_col = c0, end_col = c1;
		while (true) {
			bool ok = true; FOR(r, r0, r1) if (cmap[r][start_col] == 5) ok = false;
			if (!ok) break;
			start_col++;
		}
		while (true) {
			bool ok = true; FOR(r, r0, r1) if (cmap[r][end_col-1] == 4) ok = false;
			if (!ok) break;
			end_col--;
		}
		score_split = MP(-1, csplit);
		xsolve_bs<MODE_CSPLIT>(rmidt, start_col, rmidb, end_col, WIDTH_ENABLE, 3);
	}
	convert_use_map = false;
}

void solve_split(int r0, int c0, int r1, int c1, int rsplit, int csplit, double res = 1.0) {
	// DB("solve_split", r0, c0, r1, c1, rsplit, csplit);
	
	int rows = r1 - r0;
	int cols = c1 - c0;
	
	int rmid = (r0+r1)/2;
	int cmid = (c0+c1)/2;
	
	assert(rsplit == -1 || csplit == -1);
	assert(rsplit != -1 || csplit != -1);
	
	// if (!(rows <= 11 && cols <= 11 && abs(rows-cols) <= 2) && 
	
	if (rows <= 11 && cols <= 11 && abs(rows-cols) <= 3) {
		//small ~square
		if (rsplit != -1) {
			score_split = MP(rsplit, -1);
			xsolve_bs<MODE_RSPLIT, SCALE_MCOST>(r0, c0, r1, c1, WIDTH_SPLIT, 3);
		} else {
			score_split = MP(-1, csplit);
			xsolve_bs<MODE_CSPLIT, SCALE_MCOST>(r0, c0, r1, c1, WIDTH_SPLIT, 3);
		}
	} else if (abs(rows-cols) <= 3) {
		//big ~square
		if (rsplit != -1) {
			enable_split(r0, c0, r1, c1, -1, cmid, (cmid-c0)*(rsplit-r0), (cmid-c0)*(r1-rsplit));
			solve_split(r0, c0, r1, cmid, rsplit, -1, res*0.5);
			solve_split(r0, cmid, r1, c1, rsplit, -1, res*0.5);
		} else {
			enable_split(r0, c0, r1, c1, rmid, -1, (rmid-r0)*(csplit-c0), (rmid-r0)*(c1-csplit));
			solve_split(r0, c0, rmid, c1, -1, csplit, res*0.5);
			solve_split(rmid, c0, r1, c1, -1, csplit, res*0.5);
		}
	} else if (rows > cols) {
		//vertical rect
		assert(rsplit != -1);
		int rmove = rsplit-cols/2;
		enable_split(r0, c0, r1, c1, rsplit+(cols&1), -1, (c1-c0)*(rmove-r0), (c1-c0)*(rsplit+(cols&1)-rmove));
		solve_split(r0, c0, rsplit+(cols&1), c1, rmove,  -1, res*0.5);
		solve_split(rsplit+(cols&1), c0, r1, c1, rmove+cols,  -1, res*0.5);
		allmoves.PB(Move(rmove, c0, cols, true));
		allmoves.PB(Move(rmove, c0, cols, true));
		rotate(rmove, c0, cols, true);
		rotate(rmove, c0, cols, true);
	} else {
		//horizontal rect
		assert(csplit != -1);
		int cmove = csplit-rows/2;
		enable_split(r0, c0, r1, c1, -1, csplit+(rows&1), (r1-r0)*(cmove-c0), (r1-r0)*(csplit+(rows&1)-cmove));
		solve_split(r0, c0, r1, csplit+(rows&1), -1, cmove, res*0.5);
		solve_split(r0, csplit+(rows&1), r1, c1, -1, cmove+rows, res*0.5);
		allmoves.PB(Move(r0, cmove, rows, true));
		allmoves.PB(Move(r0, cmove, rows, true));
		rotate(r0, cmove, rows, true);
		rotate(r0, cmove, rows, true);
	}
}

void solve_full(int r0, int c0, int r1, int c1, double res = 1.0) {
	// DB("solve_full", r0, c0, r1, c1);
	
	int rows = r1 - r0;
	int cols = c1 - c0;
	
	int cmid = (c0+c1)/2;
	int rmid = (r0+r1)/2;
	
	assert(abs(cols-rows)<=1 || abs(2*cols-rows)<=3);
	
	if (rows > 1.5 * cols) {
		convert_split = MP(rmid, -1);
		solve_split(r0, c0, r1, c1, rmid, -1, res);
		convert_split = MP(-1, -1);
		
		solve_full(r0, c0, rmid, c1, res * 0.5);
		solve_full(rmid, c0, r1, c1, res * 0.5);
	} else if (rows<=7) {
		xsolve_bs<MODE_NORMAL>(r0, c0, r1, c1, WIDTH_FULL_SMALL, 3);
	} else if (rows<=14) {
		convert_split = MP(rmid, cmid);
		score_split = MP(rmid, cmid);
		xsolve_bs<MODE_QSPLIT>(r0, c0, r1, c1, WIDTH_FULL_QSPLIT, 3);
		convert_split = MP(-1, -1);
		
		solve_full(r0, c0, rmid, cmid, res * 0.25);
		solve_full(r0, cmid, rmid, c1, res * 0.25);
		solve_full(rmid, c0, r1, cmid, res * 0.25);
		solve_full(rmid, cmid, r1, c1, res * 0.25);
	} else {
		convert_split = MP(-1, cmid);
		solve_split(r0, c0, r1, c1, -1, cmid, res);
		convert_split = MP(-1, -1);
		
		solve_full(r0, c0, r1, cmid, res * 0.5);
		solve_full(r0, cmid, r1, c1, res * 0.5);
	}
}




// Benchmarking
void generate_random(int size) {
	
}

class RotatingNumbers {public: vector<string> findSolution(int N, int P, VI grid) {           
	::N = N;
	::P = P;
	
	FOR(i, 2, 41) mcost[i] = (int)(pow(i - 1, 1.5) + 1e-9);
	
	hash_mult[0] = 1;
	FOR(i, 1, 42*42) hash_mult[i] = (hash_mult[i-1] * HASH_MULT) % MAX_HASH;
	
	xcost[0] = 0;
	FOR(i, 1, 128) xcost[i] = int(pow(i, 1.75) + 1e-9);
	
	double WIDTH_SCALE = 1.0;
	
	if (N <= 29) WIDTH_SCALE = 0.8;
	if (N <= 28) WIDTH_SCALE = 0.6;
	if (N <= 27) WIDTH_SCALE = 0.7;
	if (N <= 20) WIDTH_SCALE = 2.0;
	if (N <= 18) WIDTH_SCALE = 3.0;
	if (N <= 16) WIDTH_SCALE = 5.0;
	if (N <= 15) WIDTH_SCALE = 4.0;
	if (N <= 14) WIDTH_SCALE = 2.5;
	if (N <= 13) WIDTH_SCALE = 3.0;
	if (N <= 12) WIDTH_SCALE = 7.0;
	if (N <= 10) WIDTH_SCALE = 20.0;
	if (N <=  9) WIDTH_SCALE = 30.0;
	if (N <=  8) WIDTH_SCALE = 50.0;
	if (N <=  7) WIDTH_SCALE = 12.5;
	if (N <=  6) WIDTH_SCALE = 20.0;
	
	if (N >= 36) WIDTH_SCALE = 0.8;
	if (N >= 38) WIDTH_SCALE = 0.5;
	
	WIDTH_SPLIT = int(80 * WIDTH_SCALE);
	WIDTH_ENABLE = int(100 * WIDTH_SCALE);
	WIDTH_FULL_SMALL = int(500 * WIDTH_SCALE);
	WIDTH_FULL_QSPLIT = int(50 * WIDTH_SCALE);
	
	REP(i, N) REP(j, N) {
		int x = grid[i*N+j];
		x--;
		m[i][j] = ((x/N)<<8)+(x%N);
	}
	REP(i, N) REP(j, N) om[i][j] = m[i][j];
	
	int orig_score = calc_score();
	
	solve_full(0, 0, N, N);
	
	cerr << "N = " << N << endl;
	cerr << "TotalTime = " << (get_time() - start_time) << endl;
	// DB(get_time() - start_time);
	
	VS vs;
	for (Move &m : allmoves) vs.PB(i2s(m.r0) + " " + i2s(m.c0) + " " + i2s(m.d) + (m.right ? " R" : " L"));
	return vs;
}};



int main() {
  RotatingNumbers prog;
  int N;
  int P;
  int num;
  vector<int> grid;

  cin >> N;
  cin >> P;
  for (int i=0; i<N*N; i++)
  {
    cin >> num;
    grid.push_back(num);
  }
  
  vector<string> ret = prog.findSolution(N,P,grid);
  cout << ret.size() << endl;
  for (int i = 0; i < (int)ret.size(); ++i)
      cout << ret[i] << endl;
  cout.flush();
}


