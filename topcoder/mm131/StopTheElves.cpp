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

int N, C, money;
double elfP;
char grid[32][32];
char pgrid[32][32];
int edge_dist[32][32];
int all_dist[32][32][32][32];
int fastest_steal[32][32];
int inside[32][32];

int dr[4] = {-1,1,0,0};
int dc[4] = {0,0,-1,1};

VPII boxes_to_place;

int vs[32][32];
void calc_dist() {
	memset(all_dist, 0x1F, sizeof(all_dist));
	REP(r, N) REP(c, N) edge_dist[r][c] = N*N-1;
	
	REP(r0, N) REP(c0, N) {
		if (grid[r0][c0] == 'T') continue;
		
		ZERO(vs);
		queue<int> q;
		q.push(0); q.push(r0); q.push(c0);
		while (!q.empty()) {
			int dist = q.front(); q.pop();
			int r = q.front(); q.pop();
			int c = q.front(); q.pop();
			all_dist[r0][c0][r][c] = dist;
			REP(d, 4) {
				int nr = r + dr[d];
				int nc = c + dc[d];
				if (nr < 0 || nr >= N || nc < 0 || nc >= N) {
					edge_dist[r0][c0] = min(edge_dist[r0][c0], dist + 1);
					continue;
				}
				if (grid[nr][nc] != 'T' && !vs[nr][nc]) {
					vs[nr][nc] = 1;
					q.push(dist+1);
					q.push(nr);
					q.push(nc);
				}
			}
		}
	}
}


const int SOURCE = 0;
const int SINK = 1;
	
int tflow[32*32*2+32][32*32*2+32];
int cflow[32*32*2+32][32*32*2+32];
VI fcn[32*32*2];
int f_q[32*32*2+32];
int f_prev[32*32*2+32];
int f_vs[32*32*2+32];
int f_vsno = 0;

void f_add_cn(int a, int b, int va, int vb) {
	fcn[a].PB(b);
	fcn[b].PB(a);
	tflow[a][b] = va;
	tflow[b][a] = vb;
	cflow[a][b] = 0;
	cflow[b][a] = 0;
}

void f_rem_cn(int a, int b) {
	assert(cflow[a][b] == 0);
	assert(cflow[b][a] == 0);
	tflow[a][b] = 0;
	tflow[b][a] = 0;
	REP(i, fcn[a].S) if (fcn[a][i] == b) {
		fcn[a].erase(fcn[a].begin()+i);
		break;
	}
	REP(i, fcn[b].S) if (fcn[b][i] == a) {
		fcn[b].erase(fcn[b].begin()+i);
		break;
	}
}

void f_build_graph(VPII rem = VPII()) {
	int remmap[32][32];
	
	REP(i, N*N*2+2) fcn[i].clear();
	
	ZERO(remmap);
	for (PII &p : rem) remmap[p.X][p.Y] = 1;
	REP(r, N) REP(c, N) if (grid[r][c] == 'P' && remmap[r][c] == 0) {
		int id = (r*N+c)+2;
		f_add_cn(SOURCE, id+N*N, N*N, 0);
	}
	
	REP(r, N) REP(c, N) if (grid[r][c] != 'T' && grid[r][c] != 'b') {
		int id = (r*N+c)+2;
		f_add_cn(id, id+N*N, 1, 0);
		if (r+1 < N && grid[r+1][c] != 'T' && grid[r+1][c] != 'b') f_add_cn(id+N*N, id+N, 4, 0);
		if (c+1 < N && grid[r][c+1] != 'T' && grid[r][c+1] != 'b') f_add_cn(id+N*N, id+1, 4, 0);
		if (r > 0 && grid[r-1][c] != 'T' && grid[r-1][c] != 'b') f_add_cn(id+N*N, id-N, 4, 0);
		if (c > 0 && grid[r][c-1] != 'T' && grid[r][c-1] != 'b') f_add_cn(id+N*N, id-1, 4, 0);
		if (r == 0 || c == 0 || r == N-1 || c == N-1) f_add_cn(id+N*N, SINK, 4, 0);
	}
}

bool f_augment() {
	f_vsno++;
	
	int st = 0;
	int en = 0;
	f_q[en++] = SOURCE;
	f_vs[SOURCE] = f_vsno;
	f_prev[SOURCE] = SOURCE;
	
	while (st<en) {
		int a = f_q[st++];
		for (int &b : fcn[a]) {
			if (f_vs[b] != f_vsno && tflow[a][b] - cflow[a][b] > 0) {
				f_q[en++] = b;
				f_vs[b] = f_vsno;
				f_prev[b] = a;
				if (b == SINK) st=en;
			}
		}
	}
	
	if (f_vs[SINK] != f_vsno) return false;

	int b = SINK;
	do {
		int a = f_prev[b];
		cflow[a][b] += 1;
		cflow[b][a] -= 1;
		b = a;
	} while (b != SOURCE);
	
	return true;
}

int f_remove_present(PII &p) {
	int PRESENT_ID = p.X*N+p.Y+2+N*N;
	
	int rv = -cflow[SOURCE][PRESENT_ID];
	while (cflow[SOURCE][PRESENT_ID]) {
		f_vsno++;
		int st = 0;
		int en = 0;
		f_q[en++] = PRESENT_ID;
		f_vs[PRESENT_ID] = f_vsno;
		f_prev[PRESENT_ID] = PRESENT_ID;
		while (st<en) {
			int a = f_q[st++];
			for (int &b: fcn[a]) {
				if (f_vs[b] != f_vsno && cflow[a][b] > 0) {
					f_q[en++] = b;
					f_vs[b] = f_vsno;
					f_prev[b] = a;
				}
			}
		}
		
		assert(f_vs[SINK] == f_vsno);
		
		int b = SINK;
		do {
			int a = f_prev[b];
			cflow[a][b] -= 1;
			cflow[b][a] += 1;
			b = a;
		} while (b != PRESENT_ID);
		
		cflow[SOURCE][PRESENT_ID] -= 1;
		cflow[PRESENT_ID][SOURCE] += 1;
	}
	
	f_rem_cn(SOURCE, PRESENT_ID);
	
	while (true) {
		if (!f_augment()) break;
		rv++;
	}
	
	return rv;
}

int f_add_present(PII &p) {
	int PRESENT_ID = p.X*N+p.Y+2+N*N;
	
	f_add_cn(SOURCE, PRESENT_ID, N*N, 0);
	
	int rv = 0; 
	while (true) {
		if (!f_augment()) break;
		rv++;
	}
	return rv;
}

VPII f_mincut() {
	f_vsno++;
	int st = 0;
	int en = 0;
	f_q[en++] = SOURCE;
	f_vs[SOURCE] = f_vsno;
	f_prev[SOURCE] = SOURCE;
	
	while (st<en) {
		int a = f_q[st++];
		for (int &b : fcn[a]) {
			if (f_vs[b] != f_vsno && tflow[a][b] - cflow[a][b] > 0) {
				f_q[en++] = b;
				f_vs[b] = f_vsno;
				f_prev[b] = a;
			}
		}
	}
	
	VPII points;
	REP(a, N*N*2+2) {
		if (f_vs[a] != f_vsno) continue;
		for (int &b : fcn[a]) if (f_vs[b] != f_vsno) {
			if (b - a == N*N) {
				points.PB(MP((a-2)/N,(a-2)%N));
			}
		}
	}
	
	return points;
}

VPII f_full(VPII rem = VPII()) {
	f_build_graph(rem);
	
	while (true) if (!f_augment()) break;
	
	return f_mincut();
}

bool SAFEBOX = true;
bool USE_REMOVE = true;
bool PLACE_BOXES = true;
bool USE_P3 = false;

int main() {
	cin >> N >> C >> elfP >> money;
	cerr << "[DATA] N = " << N << endl;
	cerr << "[DATA] C = " << C << endl;
	cerr << "[DATA] EP = " << elfP << endl;
	cerr << "[DATA] CP = " << C*elfP << endl;
	cerr << "[DATA] MC = " << money/C << endl;
	REP(r, N) REP(c, N) {char s[5]; scanf("%s", s); pgrid[r][c] = grid[r][c] = s[0];}
	
	int tno = 0;
	REP(r, N) REP(c, N) tno += grid[r][c] == 'T';
	cerr << "[DATA] TP = " << tno*1.0/N/N << endl;
	
	int pno = 0;
	REP(r, N) REP(c, N) pno += grid[r][c] == 'P';
	cerr << "[DATA] PNO = " << pno << endl;
	
	calc_dist();
	
	int p3 = 0;
	int p4 = 0;
	VPII p3_points;
	FOR(r, 1, N-1) FOR(c, 1, N-1) {
		if (grid[r][c] != 'P') continue;
		int t = 0;
		REP(d, 4) t += grid[r+dr[d]][c+dc[d]] == 'T';
		if (t == 3) p3++;
		if (t == 4) p4++;
		
		if (t == 3) REP(d, 4) if (grid[r+dr[d]][c+dc[d]] != 'T') p3_points.PB(MP(r+dr[d],c+dc[d]));
	}
	
	cerr << "[DATA] P3 = " << p3 << endl;
	cerr << "[DATA] P4 = " << p4 << endl;
	
	double start_time = get_time();
	
	double est_elves = (N*N - N) * elfP;
	double est_boxes = (money + N*N - 2*N) * 1.0 / C;
	
	int REMOVE_NO = 0;
	double CP = C*elfP;
	
	if (CP < 0.4) {
		// pass
	} else if (CP < 0.7) {
		REMOVE_NO = N <= 16 ? 0 : 5;
	} else if (CP < 1.0) {
		REMOVE_NO = N <= 12 ? 0 : N <= 20 ? 5 : N <= 23 ? 10 : 15;
	} else if (CP < 1.3) {
		REMOVE_NO = N <= 12 ? 0 : N <= 15 ? 5 : N <= 17 ? 10 : 15;
	} else if (CP < 1.6) {
		REMOVE_NO = N <= 18 ? 15 : 20;
	} else {
		REMOVE_NO = N <= 18 ? 15 : N <= 22 ? 20 : 25;
	}
	
	if (REMOVE_NO) {
		if (CP < 0.77) 
			REMOVE_NO -= 1;
		else if (CP < 0.92) 
			REMOVE_NO += 3;
		else 
			REMOVE_NO += 4;
	}
	
	REMOVE_NO = max(REMOVE_NO, (int)(est_elves - est_boxes + N/2));
	
	DB(REMOVE_NO);
	
	if (!USE_REMOVE) REMOVE_NO = 0;
	
	
	VPII xallp;
	REP(r, N) REP(c, N) if (grid[r][c] == 'P') xallp.PB(MP(r,c));
	REMOVE_NO = min(REMOVE_NO, (int)xallp.S);
	
	DB(xallp.S);
	
	
	VPII brem;
	int xv = 999;
	
	if (REMOVE_NO) {
		REP(tries, 100) {
			if (get_time() - start_time > 3.0) break;
			if (tries && brem.S == xallp.S) break;
			
			VPII allp = xallp;
			VPII rem;
			
			f_full(rem);
			while (true) {
				
				int bv = f_mincut().S * 100000 + 99999;
				int best = -1;
				
				if (rem.S >= REMOVE_NO) 
					break;
				
				VI dists(N*N, 0);
				for (PII &p : allp) dists[edge_dist[p.X][p.Y]]++;
				
				REP(i, allp.S) {
					f_remove_present(allp[i]);
					int av = f_mincut().S * 100000 + (tries == 0 ? edge_dist[allp[i].X][allp[i].Y] : edge_dist[allp[i].X][allp[i].Y] * 100 + rng.next(500));
					f_add_present(allp[i]);
					assert(edge_dist[allp[i].X][allp[i].Y] > 0);
					if (av < bv) {
						bv = av;
						best = i;
					}
				}
				if (best == -1) break;
				
				f_remove_present(allp[best]);
				rem.PB(allp[best]);
				allp[best] = allp.back();
				allp.pop_back();
				
			}
			
			int av = f_full(rem).S;
			DB(tries, av);
			if (av < xv) {
				xv = av;
				brem = rem;
			}
		}
		
		VI removed(xallp.S, 0);
		int total = brem.S;
		REP(i, xallp.S) REP(j, brem.S) if (xallp[i] == brem[j]) removed[i] = 1;
		
		int bv = f_full(brem).S;
		int step = 0;
		while (true) {
			if (get_time() - start_time > 6.0 || total == 0 || total == xallp.S) break;
			step++;
			
			int ra = -1, rb = -1;
			do {ra = rng.next(xallp.S);} while (removed[ra]);
			do {rb = rng.next(xallp.S);} while (!removed[rb]);
			
			f_remove_present(xallp[ra]);
			f_add_present(xallp[rb]);
			
			int av = f_mincut().S;
			if (av <= bv) {
				removed[ra] = 1;
				removed[rb] = 0;
				if (av < bv) {
					DB(av, step);
					brem.clear();
					REP(i, xallp.S) if (removed[i]) brem.PB(xallp[i]);
				}
				bv = av;
			} else {
				f_add_present(xallp[ra]);
				f_remove_present(xallp[rb]);
			}
		}
		DB(step);
		
		
		brem.clear();
		REP(i, xallp.S) if (removed[i]) brem.PB(xallp[i]);
		DB(total);
		DB(bv);
		DB(brem.S);
	}
	
	VPII rem = brem;
	
	
	int pts = f_full(rem).S;
	
	cerr << "[DATA] PTS = " << pts << endl;
	
	REP(i, rem.S) {
		VPII nrem = rem;
		nrem.erase(nrem.begin() + i);
		if (f_full(nrem).S == pts) {
			rem = nrem;
			i--;
		}
	}
	
	REP(turn, N*N) {
		int last_turn = true;
		REP(r, N) REP(c, N) if (grid[r][c] == 'E' && r > 0 && r < N-1 && c > 0 && c < N-1) last_turn = false;
		REP(r, N) REP(c, N) if (grid[r][c] == 'P') last_turn = false;
		if (turn == N*N-1) last_turn = true;
		
		
		VPII presents;
		REP(r, N) REP(c, N) if (grid[r][c] == 'P') presents.PB(MP(r,c));
		memset(fastest_steal, 0x1F, sizeof(fastest_steal));
		REP(r, N) REP(c, N) {
			if (grid[r][c] == 'e') for (PII &p : presents) 
				fastest_steal[r][c] = min(fastest_steal[r][c], all_dist[r][c][p.X][p.Y] + edge_dist[p.X][p.Y]);
			if (grid[r][c] == 'E')
				fastest_steal[r][c] = edge_dist[r][c];
		}
		
		
		int turns_left = N*N - turn;
		boxes_to_place.clear();
		
		REP(i, rem.S) if (grid[rem[i].X][rem[i].Y] != 'P') {
			rem.erase(rem.begin()+i);
			i--;
		}
		
		VPII points = f_full(rem);
		if (!PLACE_BOXES) points.clear();
		
		int money_left = money;
		for (PII &p : points) if (grid[p.X][p.Y] != 'b') money_left -= C;
		
		// if (USE_P3) for (PII &p : p3_points) points.PB(p);
		if (USE_P3 && p3_points.S) points.PB(p3_points[0]);
		
		int total_elfs = 0;
		REP(r, N) REP(c, N) total_elfs += grid[r][c] == 'e';
		
		double surplus_elves = max(0.0, max(0.0, (turns_left - N*5/4) * elfP) - (turns_left) / C);
		
		int presents_left = 0;
		REP(r, N) REP(c, N) presents_left += grid[r][c] == 'P';
		
		// if (money / C + presents_left - surplus_elves - total_elves < N*N && points.S == 0) {
			// USE_P3 = true;
			// PLACE_BOXES = false;
		// }
			
		
		// if (money >= (points.S + total_elfs + (int)surplus_elves) * C && rem.S) {
			// int bv = 0;
			// int best = -1;
			// REP(i, rem.S) {
				// f_add_present(rem[i]);
				// VPII npoints = f_mincut();
				// f_remove_present(rem[i]);
				// int cur_money = money - (total_elfs + (int)surplus_elves) * C;
				// for (PII &p : npoints) {
					// if (grid[p.X][p.Y] != 'b') cur_money -= C;
					// if (grid[p.X][p.Y] == 'e') cur_money -= 100;
				// }
				// int av = cur_money;
				// if (av > bv || av == bv && best != -1 && edge_dist[rem[i].X][rem[i].Y] > edge_dist[rem[best].X][rem[best].Y]) {
					// bv = av;
					// best = i;
				// }
			// }
			
			// if (best != -1) {
				// rem.erase(rem.begin() + best);
				// points = f_full(rem);
			// }
			
			// DB(turn, rem.S, money);
		// }

		while (money_left >= (total_elfs + (int)surplus_elves) * C && rem.S) {
			int bv = C * max(0, (total_elfs + (int)surplus_elves) - 1);
			int best = -1;
			REP(i, rem.S) {
				f_add_present(rem[i]);
				VPII npoints = f_mincut();
				f_remove_present(rem[i]);
				int cur_money = money;
				for (PII &p : npoints) {
					if (grid[p.X][p.Y] != 'b') cur_money -= C;
					if (grid[p.X][p.Y] == 'e') cur_money -= 100;
				}
				int av = cur_money;
				if (av > bv || av == bv && best != -1 && edge_dist[rem[i].X][rem[i].Y] > edge_dist[rem[best].X][rem[best].Y]) {
					bv = av;
					best = i;
				}
			}
			
			if (best == -1) break;
			
			rem.erase(rem.begin() + best);
			money_left = bv;
			points = f_full(rem);
			
			DB(turn, rem.S, money_left);
		}
		
		if (money_left + 50 < (total_elfs + (int)surplus_elves) * C && rem.S < presents.S && turns_left > N && get_time() - start_time < 6.0 + (turn * 3.0 / (N*N)) && turn % 5 == 0) {
			int bv = f_mincut().S * 1000;
			int best = -1;
			REP(i, presents.S) {
				bool exist = false;
				REP(j, rem.S) exist |= presents[i] == rem[j];
				if (exist) continue;
				
				f_remove_present(presents[i]);
				int av = f_mincut().S * 1000;
				f_add_present(presents[i]);
				
				if (av < bv) {
					bv = av;
					best = i;
				}
			}
			
			if (best != -1) {
				cerr << "Removing: " << presents[best] << endl;
				rem.PB(presents[best]);
				points = f_full(rem);
			}
		}
		
		
		if (rem.S && rem.S < presents.S && get_time() - start_time < 6.0 + (turn * 3.0 / (N*N))) {
			
			VI removed(presents.S, 0);
			int total = rem.S;
			REP(i, presents.S) REP(j, rem.S) if (presents[i] == rem[j]) removed[i] = 1;
			
			int bv = points.S * 1000 + total;
			// DB(turn, money / C - (points.S + total_elfs + (int)surplus_elves), bv);
			REP(step, 500) {
				
				int ra = -1, rb = -1;
				
				if (total > 0 && total < presents.S && rng.next(2)) {
					do {ra = rng.next(presents.S);} while (removed[ra]);
					do {rb = rng.next(presents.S);} while (!removed[rb]);
				} else if (total > 0) {
					do {rb = rng.next(presents.S);} while (!removed[rb]);
				} else {
					continue;
				}
				
				if (ra != -1) f_remove_present(presents[ra]);
				if (rb != -1) f_add_present(presents[rb]);
				
				int av = f_mincut().S * 1000 + total - (ra == -1);
				if (av <= bv) {
					if (ra != -1) removed[ra] = 1, total++;
					if (rb != -1) removed[rb] = 0, total--;
					if (av < bv) {
						DB(turn, av, bv, step);
						rem.clear();
						REP(i, presents.S) if (removed[i]) rem.PB(presents[i]);
						points = f_mincut();
					}
					bv = av;
				} else {
					if (ra != -1) f_add_present(presents[ra]);
					if (rb != -1) f_remove_present(presents[rb]);
				}
			}
		}
		
		//mark_inside
		ZERO(inside);
		for (PII &p : points) inside[p.X][p.Y] = 1;
		REP(r, N) REP(c, N) if (grid[r][c] == 'T' || grid[r][c] == 'b') inside[r][c] = 1;
		for (PII &p : presents) {
			bool removed = false;
			for (PII &r : rem) if (p == r) removed = true;
			if (removed) continue;
			queue<int> q;
			q.push(p.X);
			q.push(p.Y);
			while (!q.empty()) {
				int r = q.front(); q.pop();
				int c = q.front(); q.pop();
				if (r < 0 || r >= N || c < 0 || c >= N) continue;
				if (inside[r][c]) continue;
				inside[r][c] = 1;
				REP(d, 4) {
					q.push(r + dr[d]);
					q.push(c + dc[d]);
				}
			}
		}
		
		// for (PII &p : points) if (grid[p.X][p.Y] == '.') {
			// REP(d, 4) REP(d4, 4) {
				// int nr = p.X + dr[d] + dr[d4];
				// int nc = p.Y + dc[d] + dc[d4];
				// if (nr >= 0 && nr < N && nc >= 0 && nc < N && grid[nr][nc] == 'E' && inside[nr][nc] && money >= C && turns_left >= fastest_steal[nr][nc] && grid[p.X][p.Y] == '.') {
					// boxes_to_place.PB(p);
					// grid[p.X][p.Y] = 'b';
					// money -= C;
				// }
				 // nr = p.X + dr[d];
				 // nc = p.Y + dc[d];
				// if (nr >= 0 && nr < N && nc >= 0 && nc < N && grid[nr][nc] == 'e' && money >= C && turns_left >= fastest_steal[nr][nc] && grid[p.X][p.Y] == '.') {
					// boxes_to_place.PB(p);
					// grid[p.X][p.Y] = 'b';
					// money -= C;
				// }
			// }
		// }
		for (PII &p : points) if (grid[p.X][p.Y] == '.') {
			REP(d, 4) {
				int nr = p.X + dr[d];
				int nc = p.Y + dc[d];
				if (nr >= 0 && nr < N && nc >= 0 && nc < N && (grid[nr][nc] == 'e' || grid[nr][nc] == 'E' && inside[nr][nc]) && money >= C && turns_left >= fastest_steal[nr][nc]) {
					boxes_to_place.PB(p);
					grid[p.X][p.Y] = 'b';
					money -= C;
					
					if (SAFEBOX && CP < 1) {
						REP(d2, 4) {
							int nr2 = nr + dr[d2];
							int nc2 = nc + dc[d2];
							if (nr2 >= 0 && nr2 < N && nc2 >= 0 && nc2 < N && grid[nr2][nc2] == 'e' && money >= C && turns_left >= fastest_steal[nr2][nc2]) {
								REP(d3, 4) {
									int nr3 = nr2 + dr[d3];
									int nc3 = nc2 + dc[d3];
									if (nr3 >= 1 && nr3 < N-1 && nc3 >= 1 && nc3 < N-1 && grid[nr3][nc3] == '.') {
										boxes_to_place.PB(MP(nr3,nc3));
										grid[nr3][nc3] = 'b';
										money -= C;
										break;
									}
								}
							}
						}
					}
					break;
				}
			}
		}
		
		if (last_turn) {
			int boxes_left = 0;
			int elves_presents = 0;
			int elves_presents_edge = 0;
			REP(r, N) REP(c, N) {
				if (grid[r][c] == 'b') boxes_left++;
				if (grid[r][c] == 'E') {
					if (r == 0 || r == N-1 || c == 0 || c == N-1) elves_presents_edge++;
					else elves_presents++;
				}
			}
			cerr << "[DATA] BL = " << boxes_left << endl;
			cerr << "[DATA] EWP = " << elves_presents << endl;
			cerr << "[DATA] BWPE = " << elves_presents_edge << endl;
			this_thread::sleep_for(chrono::milliseconds(100));
			fflush(stderr);
		}
		
		if (boxes_to_place.S) {
			for (PII &p : boxes_to_place) printf("%d %d ", p.X, p.Y);
		} else {
			printf("-1");
		}
		printf("\n");
		fflush(stdout);
		
		int elapsedTime;
		scanf("%d%d", &elapsedTime, &money);
		
		REP(r, N) REP(c, N) pgrid[r][c] = grid[r][c];
		REP(r, N) REP(c, N) {char s[5]; scanf("%s", s); grid[r][c] = s[0];}
	}

	return 0;
}
