// Author: Psyho
// Twitter: https://twitter.com/fakepsyho
// Site: http://psyho.gg/

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

#define SQR(a) ((a)*(a))


double max_c = 0;
double max_r = 0;


INLINE bool circ_col(int r1, int x1, int y1, int r2, int x2, int y2) {
	return SQR(x2-x1)+SQR(y2-y1) < SQR(r1+r2);
}

INLINE bool rect_col(int w1, int h1, int x1, int y1, int w2, int h2, int x2, int y2) {
	return x1 < x2+w2 && x1+w1 > x2 &&
	       y1 < y2+h2 && y1+h1 > y2;
}

const double EPS=1e-6;
INLINE bool mix_col(int r1, int x1, int y1, int w2, int h2, int x2, int y2) {
	double x = abs(x1 - (x2+w2/2.0));
	double y = abs(y1 - (y2+h2/2.0));
	if (x+EPS > w2/2.0 + r1) return false;
	if (y+EPS > h2/2.0 + r1) return false;
	if (x+EPS < w2/2.0) return true;
	if (y+EPS < h2/2.0) return true;
	double corner_dist = SQR(x - w2/2.0) + SQR(y - h2/2.0);
	return corner_dist+EPS < SQR(r1);
}

bool circ_ins(int r, int x, int y) {
	return SQR(x) + SQR(y) <= SQR(100-r);
}

bool rect_ins(int w, int h, int x, int y) {
	return SQR(x)+SQR(y) <= 100*100 &&
	       SQR(x)+SQR(y+h) <= 100*100 &&
	       SQR(x+w)+SQR(y) <= 100*100 &&
	       SQR(x+w)+SQR(y+h) <= 100*100;
}

const double INF=1000;


double xtime;

int C,R,N;
double X;

VI   fc;
VPII fr;

double score(int c, int r) {
	return r + c - abs(X*c - (1-X)*r);
}

PDD get_area(VPII &s) {
	double c = 0;
	double r = 0;
	REP(i,C) if (s[i].X < INF) c += fc[i]*fc[i]*PI;
	REP(i,R) if (s[i+C].X < INF) r += fr[i].X*fr[i].Y;
	return MP(c,r);
}

double score(VPII &s) {
	PII p = get_area(s);
	return score(p.X,p.Y);
}

INLINE bool circ_ok(VPII &s, int id, int x, int y) {
	if (!circ_ins(fc[id],x,y)) return false;
	REP(i,C) if (i!=id && s[i].X < INF && circ_col(fc[id],x,y,fc[i],s[i].X,s[i].Y)) return false;
	REP(i,R) if (s[i+C].X < INF && mix_col(fc[id],x,y,fr[i].X,fr[i].Y,s[i+C].X,s[i+C].Y)) return false;
	return true;
}

INLINE bool rect_ok(VPII &s, int id, int x, int y) {
	if (!rect_ins(fr[id].X,fr[id].Y,x,y)) return false;
	REP(i,C) if (s[i].X < INF && mix_col(fc[i],s[i].X,s[i].Y,fr[id].X,fr[id].Y,x,y)) return false;
	REP(i,R) if (i!=id && s[i+C].X < INF && rect_col(fr[id].X,fr[id].Y,x,y,fr[i].X,fr[i].Y,s[i+C].X,s[i+C].Y)) return false;
	return true;
}

const int BLOCK=10;
const int MAX_FIGS=10;
struct State {
	int gno[200/BLOCK+1][200/BLOCK+1];
	int gid[200/BLOCK+1][200/BLOCK+1][MAX_FIGS];
	int bad[256][256];
	
	int check[1000];
	int checkid = 0;
	
	void clear() {
		ZERO(gno);
		ZERO(gid);
		ZERO(check);
		ZERO(bad);
		FOR(ox,-128,128) FOR(oy,-128,128)
			if (SQR(ox)+SQR(oy)>100*100+EPS) bad[ox+128][oy+128] = 1;
	}
	
	State() {
		clear();
	}
	
	INLINE void add_fig(int id, int bx, int by) {
		gid[bx][by][gno[bx][by]++] = id;
	}
	
	INLINE void add_circ(int id, int x, int y) {
		FOR(bx,(100+x-fc[id])/BLOCK,(100+x+fc[id])/BLOCK+1)
			FOR(by,(100+y-fc[id])/BLOCK,(100+y+fc[id])/BLOCK+1) {
				this->add_fig(id, bx, by);
			}
		FOR(ox,-fc[id]+1,fc[id]) FOR(oy,-fc[id]+1,fc[id]) if (SQR(ox)+SQR(oy)+EPS<SQR(fc[id])) bad[ox+128+x][oy+128+y] = 1;
	}
	
	INLINE void add_rect(int id, int x, int y) {
		FOR(bx,(100+x)/BLOCK,(100+x+fr[id].X)/BLOCK+1)
			FOR(by,(100+y)/BLOCK,(100+y+fr[id].Y)/BLOCK+1) {
				this->add_fig(id+C, bx, by);
			}
		FOR(ox,1,fr[id].X) FOR(oy,1,fr[id].Y) bad[ox+128+x][oy+128+y] = 1;
	}
	
	INLINE void rem_circ(int id, int x, int y) {
		//TODO: add
	}
	
	INLINE void rem_rect(int id, int x, int y) {
		//TODO: add
	}
	
	INLINE bool check_circ(VPII &s, int id, int x, int y) {
		if (bad[x+128][y+128]) return false;
		if (bad[x+128-fc[id]+1][y+128]) return false;
		if (bad[x+128+fc[id]-1][y+128]) return false;
		if (bad[x+128][y+128-fc[id]+1]) return false;
		if (bad[x+128][y+128+fc[id]-1]) return false;
		checkid++;
		if (!circ_ins(fc[id],x,y)) return false;
		FOR(bx,(100+x-fc[id])/BLOCK,(100+x+fc[id])/BLOCK+1)
			FOR(by,(100+y-fc[id])/BLOCK,(100+y+fc[id])/BLOCK+1) 
				REP(i, gno[bx][by]) {
					int p = gid[bx][by][i];
					if (check[p] != checkid) {
						if (p<C) {
							if (circ_col(fc[id],x,y,fc[p],s[p].X,s[p].Y)) return false;
						} else {
							if (mix_col(fc[id],x,y,fr[p-C].X,fr[p-C].Y,s[p].X,s[p].Y)) return false;
						}
						check[p] = checkid;
					}
				}
		return true;
	}
	
	INLINE bool check_rect(VPII &s, int id, int x, int y) {
		if (bad[x+128][y+128]) return false;
		if (bad[x+128][y+128+fr[id].Y]) return false;
		if (bad[x+128+fr[id].X][y+128]) return false;
		if (bad[x+128+fr[id].X][y+128+fr[id].Y]) return false;
		checkid++;
		if (!rect_ins(fr[id].X,fr[id].Y,x,y)) return false;
		FOR(bx,(100+x)/BLOCK,(100+x+fr[id].X)/BLOCK+1)
			FOR(by,(100+y)/BLOCK,(100+y+fr[id].Y)/BLOCK+1)
				REP(i, gno[bx][by]) {
					int p = gid[bx][by][i];
					if (check[p] != checkid) {
						if (p<C) {
							if (mix_col(fc[p],s[p].X,s[p].Y,fr[id].X,fr[id].Y,x,y)) return false;
						} else {
							if (rect_col(fr[id].X,fr[id].Y,x,y,fr[p-C].X,fr[p-C].Y,s[p].X,s[p].Y)) return false;
						}
						check[p] = checkid;
					}
				}
		return true;
	}
};

struct Par {
	int type;
	int type2;
	int startr;
	int offx;
	int offy;
	int offx2;
	int offy2;
	int diff;
	int noiserect;
	double sqrat;
	double sqrat2;
};

ostream& operator<<(ostream& os, const Par& p) {
	return os << "t:" << p.type << ",ox:" << p.offx << ",oy:" << p.offy << ",d:" << p.diff << ",r:" << p.sqrat << ",r2:" << p.sqrat2 << ",?:" << p.type2 << ",sr:" << p.startr;
}

double timer1 = 0;
double timer2 = 0;
double timer3 = 0;
double timer4 = 0;
double timer5 = 0;
double timer6 = 0;
double timer7 = 0;
double timer8 = 0;
double timer9 = 0;

double vval[40000];
double xval[40000];


State state;

PII vpos[40000];
PII xpos[40000];
int vposno;
int xposno;

int vsort[40000];
int xsort[40000];

VI solve_knapsack(int t, VI &v) {
	int tot = 0;
	REP(i, v.S) tot += v.S;
	if (tot < t) return VI();
	
	tot = 0;
	VI used(v.S);
	REP(i, 2000) {
		int p = rng.next(v.S);
		int p2 = rng.next(v.S);
		int p3 = rng.next(v.S);
		int p4 = rng.next(v.S);
		if (tot < t && v[p]<v[p2]) p = p2;
		if (tot < t && v[p]<v[p3]) p = p3;
		if (tot < t && v[p]<v[p4]) p = p4;
		
		if (used[p] == 0 && tot < t) {
			used[p] = 1;
			tot += v[p];
		} else if (used[p] && tot > t) {
			used[p] = 0;
			tot -= v[p];
		}
		if (tot == t) {
			VI rv;
			REP(i, v.S) if (used[i]) rv.PB(i);
			return rv;
		}
	}
	
	return VI();
}

bool solve_hrect(VPII &s, int x1, int y1, int x2, int y2) {
	assert(y2-y1>=5 && y2-y1<=25);
	VI ord, v;
	REP(i,R) if (s[i+C].X == INF && fr[i].Y == y2-y1) ord.PB(i), v.PB(fr[i].X);
	VI rv = solve_knapsack(x2-x1, v);
	if (rv.S == 0) return false;
	int sum = 0;
	REP(i, rv.S) s[C+ord[rv[i]]] = MP(x1+sum,y1), sum += v[rv[i]];
	assert(sum == 0 || sum == x2-x1);
	return true;
}

bool solve_hrect(VPII &s, int y1, int y2, int side = 0) {
	int x1 = 0;
	int x2 = 0;
	while (SQR(x1-1)+SQR(y1)<10000-EPS && SQR(x1-1)+SQR(y2)<10000-EPS) x1--;
	while (SQR(x2+1)+SQR(y1)<10000-EPS && SQR(x2+1)+SQR(y2)<10000-EPS) x2++;
	if (side == -1) x1 = 0;
	if (side == +1) x2 = 0;
	if (x2-x1<5) return true;
	return solve_hrect(s, x1, y1, x2, y2);
}

bool solve_top(VPII &s, int y1, int y2, int side = 0) {
	VI sizes; FOR(i,5,26) sizes.PB(i);
	while (true) {
		if (y2-y1 < 5) break;
		REP(i, sizes.S) {
			int p1 = rng.next(i, sizes.S);
			swap(sizes[i], sizes[p1]);
			if (y1+sizes[i] <= y2 && solve_hrect(s, y1, y1+sizes[i], side)) {
				y1 += sizes[i];
				goto niceh;
			}
		}
		break;
		niceh: ;
	}
	return true;
}

bool solve_diag(VPII &s, int x1, int y1, int x2, int y2) {
	const int MN = 22;
	VI bigc;
	REP(i, C) if (fc[i]==15) bigc.PB(i);
	int bigcpos = 0;
	
	for (int x = x1; x <= x2; x += MN)
		for (int y = y1; y <= y2; y += MN) {
			if ((x-x1+y-y1)%(2*MN)==MN) continue;
			if (!circ_ins(15, x, y)) continue;
			if (bigcpos == bigc.S) return true;
			s[bigc[bigcpos++]] = MP(x, y);
		}
	return true;
}

bool solve_vrect(VPII &s, int x1, int y1, int x2, int y2) {
	assert(x2-x1>=5 && x2-x1<=25);
	VI ord, v;
	REP(i,R) if (s[i+C].X == INF && fr[i].X == x2-x1) ord.PB(i), v.PB(fr[i].Y);
	VI rv = solve_knapsack(y2-y1, v);
	if (rv.S == 0) return false;
	int sum = 0;
	REP(i, rv.S) s[C+ord[rv[i]]] = MP(x1,y1+sum), sum += v[rv[i]];
	assert(sum == 0 || sum == y2-y1);
	return true;
}


VPII gen_greedy(Par par, VPII init=VPII()) {
	TIME(double stimer3 = get_time();)
	state.clear();
	TIME(timer3 += get_time() - stimer3;)
	
	VPII s = init.S ? init : VPII(N, MP(INF,INF));
	int usedr = 0;
	int usedc = 0;
	
	TIME(double stimer1 = get_time();)
	
	if (init.S) {
		REP(i,C) if (s[i].X < INF)   state.add_circ(i, s[i].X, s[i].Y);
		REP(i,R) if (s[i+C].X < INF) state.add_rect(i, s[i+C].X, s[i+C].Y);
	}
	
	VI orderc; REP(i,C) if (s[i].X == INF) orderc.PB(i);
	sort(ALL(orderc), [&](int a, int b) -> bool {return fc[a]>fc[b];});
	
	VI orderr; REP(i,R) if (s[i+C].X == INF) orderr.PB(i);
	VPII xorder;
	if (par.type == 0) {
		for (int i : orderr) xorder.PB(MP(-fr[i].X*100-fr[i].Y,i));
	} else if (par.type == 1) {
		for (int i : orderr) xorder.PB(MP(-fr[i].Y*100-fr[i].X,i));
	} else if (par.type == 2) {
		for (int i : orderr) xorder.PB(MP(-fr[i].Y*fr[i].X,i));
	}
	if (par.noiserect) REP(i, xorder.S) xorder[i].X += rng.next(5);
	sort(ALL(xorder));
	orderr.clear(); REP(i, xorder.S) orderr.PB(xorder[i].Y);
	
	TIME(timer1 += get_time() - stimer1;)
	
	double max_total_area = 30000;
	
	double max_circ_area = max_total_area * (1 - X);
	double max_rect_area = max_total_area * X;
	
	TIME(double stimer2 = get_time();)
	vposno=0;
	FOR(x,-95,96) FOR(y,-95,96) {
		vsort[vposno++] = (x+100)*200+y+100;
		vval[(x+100)*200+y+100] = SQR(x+par.offx2)*par.sqrat+SQR(y+par.offy2);
	}
	sort(vsort,vsort+vposno, [&](int a, int b) -> bool {return vval[a]<vval[b];});
	REP(i, vposno) vpos[i] = MP(vsort[i]/200-100,vsort[i]%200-100);
	
	xposno=0;
	FOR(x,-100,95) FOR(y,-100,95) {
		xsort[xposno++] = (x+100)*200+y+100;
		xval[(x+100)*200+y+100] = abs(x-par.offx)*par.sqrat2+abs(y-par.offy);
	}
	
	sort(xsort,xsort+xposno, [&](int a, int b) -> bool {return xval[a]<xval[b];});
	REP(i, xposno) xpos[i] = MP(xsort[i]/200-100,xsort[i]%200-100);
	TIME(timer2 += get_time() - stimer2;)
	
	
	int last_pos_r = 0;
	int last_pos_c = 0;
	int mom = init.S ? 0 : 25;
	while (true) {
		if (usedr == orderr.S && usedc == orderc.S) break;
		
		if (usedr+usedc>=mom) {
			last_pos_c=0;
			last_pos_r=0;
			
			int vposnew=0;
			REP(i, vposno) if (state.bad[vpos[i].X+128][vpos[i].Y+128] == 0)
				vpos[vposnew++] = vpos[i];
			vposno = vposnew;
			
			int xposnew=0;
			REP(i, xposno) if (state.bad[xpos[i].X+128][xpos[i].Y+128] == 0)
				xpos[xposnew++] = xpos[i];
			xposno = xposnew;
			
			mom += mom < 100 ? 25 : 50;
			if (mom >= 500) mom = 100000;
		}
		
		PDD p = get_area(s);
		double c = p.X; 
		double r = p.Y;
		if (X*c < (1-X)*r + par.diff && r >= par.startr && usedc < orderc.S || usedr == orderr.S && usedc < orderc.S) {
			TIME(double stimer4 = get_time();)
			// greedy circle
			while (c + fc[orderc[usedc]]*fc[orderc[usedc]]*PI >= max_circ_area && usedc+1 < orderc.S && usedr < orderr.S && c + fc[orderc[usedc+1]]*fc[orderc[usedc+1]]*PI >= max_circ_area && max_r > max_total_area )
				usedc++;
			
			FOR(i, last_pos_c, vposno) {
				PII &p = vpos[i];
				if (state.check_circ(s, orderc[usedc], p.X, p.Y)) {
					s[orderc[usedc]] = p;
					state.add_circ(orderc[usedc],p.X,p.Y);
					last_pos_c = i;
					goto donec;
				}
			}
			
			while (usedc+1 < orderc.S && fc[orderc[usedc]] == fc[orderc[usedc+1]]) usedc++;
			
			donec: ;
			if (usedc+1 < orderc.S && fc[orderc[usedc]] != fc[orderc[usedc+1]]) last_pos_c = 0;
			usedc++;
			TIME(timer4 += get_time() - stimer4;)
		} else {
			TIME(double stimer5 = get_time();)
			// greedy rect
			while (r + fr[orderr[usedr]].X*fr[orderr[usedr]].Y >= max_rect_area && usedr+1 < orderr.S && usedc < orderc.S && r + fr[orderr[usedr+1]].X*fr[orderr[usedr+1]].Y >= max_rect_area && max_c > max_total_area)
				usedr++;
			
			FOR(i, last_pos_r, xposno) {
				PII &p = xpos[i];
				if (state.check_rect(s, orderr[usedr], p.X, p.Y)) {
					s[C+orderr[usedr]] = p;
					state.add_rect(orderr[usedr], p.X, p.Y);
					last_pos_r = i;
					goto doner;
				}
			}
			while (usedr+1 < orderr.S && fr[orderr[usedr]] == fr[orderr[usedr+1]]) usedr++;
			
			doner: ;
			if (usedr+1 < orderr.S && fr[orderr[usedr]] != fr[orderr[usedr+1]]) last_pos_r = 0;
			usedr++;
			TIME(timer5 += get_time() - stimer5;)
		}
	}
	
	return s;
}

class TastyPizza {public: VS findSolution(int C, int R, double X, VI &circles, VPII &rects) { 
	xtime = get_time();
	
	DB(R,C,X);
	
	::X = X;
	::R = R;
	::C = C;
	N = R+C;
	fc = circles;
	fr = rects;
	
	REP(i, C) max_c += fc[i]*fc[i]*PI;
	REP(i, R) max_r += fr[i].X*fr[i].Y;
	DB(max_c, max_r);
	
	double bv = 0;
	
	int run = 0;
	
	VPII bsol;
	
	VD res;
	VC<Par> pars;
	int mindiff = -1000;
	int maxdiff = 1000;
	
	int minstartr = 0;
	int maxstartr = X > 0.7 ? min((int)max_r, 30000) : 0;
	
	while (true) {
		double passed_time = (get_time() - xtime) / TIME_LIMIT;
		if (passed_time > 1.0) break;
		int mode = 0;//passed_time >= 0.8;
		
		Par par;
		par.type = rng.next(1,3);
		par.type2 = rng.next(2)*2-1;
		par.offx = rng.next(0,50);
		par.offy = rng.next(0,50);
		par.startr = rng.next(minstartr, maxstartr+1);
		// double dist  = rng.next_double(0,40.0);
		// double angle = rng.next_double(0,2*PI);
		// par.offx = int(dist*cos(angle));
		// par.offy = int(dist*sin(angle));
		par.offx2 = par.offx;
		par.offy2 = par.offy;
		par.diff = rng.next(mindiff,maxdiff);
		par.sqrat = pow(4.0, rng.next_double(-1, 1));
		par.sqrat2 = pow(4.0, rng.next_double(-1, 1));
		par.noiserect = mode ? rng.next(2) : 0;
		
		VPII sol;
		if (mode == 0) {
			if (false && X > .1 && X < .5 && rng.next(2)) {
				VPII xsol(N, MP(INF, INF));
				int x1 = rng.next(-128, 128);
				int x2 = rng.next(-128, 128);
				int y1 = rng.next(-128, 128);
				int y2 = rng.next(-128, 128);
				if (x1 > x2) swap(x1, x2);
				if (y1 > y2) swap(y1, y2);
				solve_diag(xsol, x1, y1, x2, y2);
				sol = gen_greedy(par, xsol);
			} else if (X > .4 && rng.next(4)) {
				VPII xsol(N, MP(INF, INF));
				int bot = rng.next(-100, 100);
				int top = rng.next(75, 90);
				if (rng.next(2)) {
					solve_top(xsol, bot, top, -1);
					solve_top(xsol, bot, top, +1);
				} else {
					solve_top(xsol, bot, top);
				}
				sol = gen_greedy(par, xsol);
			} else {
				sol = gen_greedy(par);
			}
		} else {
			VPII init = bsol;
			int x1 = rng.next(-100,100);
			int x2 = rng.next(2)*2000-1000;
			int y1 = rng.next(-100,100);
			int y2 = rng.next(2)*2000-1000;
			if (x1>x2) swap(x1, x2);
			if (y1>y2) swap(y1, y2);
			REP(i, C) if (bsol[i].X   < INF && rect_col(x2-x1,y2-y1,x1,y1,fc[i]*2,fc[i]*2,bsol[i].X-fc[i],bsol[i].Y-fc[i])) init[i]   = MP(INF, INF);
			REP(i, R) if (bsol[i+C].X < INF && rect_col(x2-x1,y2-y1,x1,y1,fr[i].X,fr[i].Y,bsol[i+C].X,bsol[i+C].Y)) 		init[i+C] = MP(INF, INF);
			
			sol = gen_greedy(par, init);
		}
		
		double av = score(sol);
		if (mode == 1) {
			DB(run, bv, av);
		}
		
		if (av > bv) {
			bv = av;
			bsol = sol;
		}
		// DB(run, av);
		
		res.PB(av);
		pars.PB(par);
		run++;
		
		if (run % 50 == 0) {
			VI br; REP(i, res.S) br.PB(i);
			sort(ALL(br), [&](int a, int b) -> bool {return res[a]>res[b];});
			bool onlypos = true;
			mindiff = 10000;
			maxdiff = -10000;
			minstartr = 100000;
			maxstartr = -100000;
			REP(i, 5) {
				mindiff = min(mindiff, pars[br[i]].diff);
				maxdiff = max(maxdiff, pars[br[i]].diff);
				minstartr = min(minstartr, pars[br[i]].startr);
				maxstartr = max(maxstartr, pars[br[i]].startr);
			}
			mindiff -= 50;
			maxdiff += 50;
			if (maxstartr) maxstartr = min((int)max_r, maxstartr+500);
			if (minstartr) minstartr = max(0, minstartr-500);
		}
	}
	DB(run);
	DB(mindiff, maxdiff);
	
	VI br; REP(i, res.S) br.PB(i);
	sort(ALL(br), [&](int a, int b) -> bool {return res[a]>res[b];});
	REP(i, min(10, (int)br.S)) {
		int x = br[i];
		DB(x, res[x], pars[x]);
	}
	
	PDD p = get_area(bsol);
	DB(p.X,p.Y,abs(X*p.X-(1-X)*p.Y));
	
	DB(get_time() - xtime);
	
	if (timer1) DB(timer1);
	if (timer2) DB(timer2);
	if (timer3) DB(timer3);
	if (timer4) DB(timer4);
	if (timer5) DB(timer5);
	if (timer6) DB(timer6);
	if (timer7) DB(timer7);
	if (timer8) DB(timer8);
	if (timer9) DB(timer9);
	
	// cerr << "TotalTime = " << (get_time() - xtime) << endl;
	
	VS rv;
	REP(i, N) rv.PB(bsol[i].X >= INF ? "NA" : i2s(bsol[i].X)+" "+i2s(bsol[i].Y));
	return rv; 
}};

int main() {
  TastyPizza prog;
  int R, C;
  double X;
  cin >> C >> R >> X;
  VPII aRect(R);
  VI aCircles(C);
  REP(i, C) cin >> aCircles[i];
  REP(i, R) cin >> aRect[i].first >> aRect[i].second;
  VS ret = prog.findSolution(C, R, X, aCircles, aRect);
  cout << ret.S << endl;
  REP(i, ret.S) cout << ret[i] << endl;
  cout.flush();
}