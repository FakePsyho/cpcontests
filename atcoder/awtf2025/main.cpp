// Author: Psyho
// Twitter: twitter.com/fakepsyho
// BlueSky: psyho.bsky.social
 
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
#define DATA(x) {cerr << "[DATA] " << #x << " = " << (x) << endl;}
 
double get_time() {timeval tv; gettimeofday(&tv, NULL); return tv.tv_sec + tv.tv_usec * 1e-6;}
double start_time = get_time();
double elapsed() {return get_time() - start_time;}

struct Timer {
    string name;
    double total = 0;
    int count = 0;
    double start_time;
    Timer(string name="") {this->name = name;}
    void reset() {total = 0; count = 0;}
    void start() {start_time = get_time();}
    void stop() {total += get_time() - start_time; count++;}
    double avg() {return total / count;}
};
 
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
    INLINE int next(int x) {return ((LL)rand() * x) >> 32; }
    INLINE int next(int a, int b) {return a + next(b - a); }
    INLINE double next_double() {return (rand() + 0.5) * (1.0 / 4294967296.0); }
    INLINE double next_double(double a, double b) {return a + next_double() * (b - a); }
};
 
static RNG rng;

// SOLUTION

const bool SILENT = true;

const string dirs = "UDLR";

int _;

const int MAXK = 100;

const int N = 30;
int K;

PII src[MAXK];
PII dst[MAXK];

int owallv[N+2][N+2];
int owallh[N+2][N+2];

PII pos[MAXK];
int wallv[N+2][N+2];
int wallh[N+2][N+2];
int cell[N][N];

int wallv_mark[N+2][N+2];
int wallh_mark[N+2][N+2];

void reset() {
    ZERO(cell);
    memcpy(pos, src, sizeof(src));
    REP(i, K) cell[src[i].Y][src[i].X] = 1;
}

void fast_reset() {
    memcpy(pos, src, sizeof(src));
}

int cell_copy[N][N];
PII pos_copy[MAXK];


int order[N][N];
int n_order[N];

void state_save() {
    memcpy(cell_copy, cell, sizeof(cell_copy));
    memcpy(pos_copy, pos, sizeof(pos_copy));
}

void state_load() {
    memcpy(cell, cell_copy, sizeof(cell_copy));
    memcpy(pos, pos_copy, sizeof(pos_copy));
}


int next_wallu[N][N];
int next_walld[N][N];
int next_walll[N][N];
int next_wallr[N][N];

void rebuild_next_wall_col(int c) {
    REP(r, N) next_wallu[r][c] = wallh[r][c] ? r : next_wallu[r-1][c];
    for (int r = N - 1; r >= 0; r--) next_walld[r][c] = wallh[r+1][c] ? r : next_walld[r+1][c];
}

void rebuild_next_wall_row(int r) {
    REP(c, N) next_walll[r][c] = wallv[r][c] ? c : next_walll[r][c-1];
    for (int c = N - 1; c >= 0; c--) next_wallr[r][c] = wallv[r][c+1] ? c : next_wallr[r][c+1];
}


int next_pos[N];

void fmoveu_xfast(int n) {
    ZERO(n_order); 
    REP(i, N) next_pos[i] = -1;
    REP(i, K) order[pos[i].Y][n_order[pos[i].Y]++] = i;
    REP(j, N) REP(k, n_order[j]) {
        int i = order[j][k]; 
        pos[i].Y = max(pos[i].Y - n, max(next_pos[pos[i].X] + 1, next_wallu[pos[i].Y][pos[i].X]));
        next_pos[pos[i].X] = pos[i].Y;
    }
}

void fmoved_xfast(int n) {
    ZERO(n_order);
    REP(i, N) next_pos[i] = N;
    REP(i, K) order[N-1-pos[i].Y][n_order[N-1-pos[i].Y]++] = i;
    REP(j, N) REP(k, n_order[j]) {
        int i = order[j][k];  
        pos[i].Y = min(pos[i].Y + n, min(next_pos[pos[i].X] - 1, next_walld[pos[i].Y][pos[i].X]));
        next_pos[pos[i].X] = pos[i].Y;
    }
}

void fmovel_xfast(int n) {
    ZERO(n_order); 
    REP(i, N) next_pos[i] = -1;
    REP(i, K) order[pos[i].X][n_order[pos[i].X]++] = i;
    REP(j, N) REP(k, n_order[j]) {
        int i = order[j][k];  
        pos[i].X = max(pos[i].X - n, max(next_pos[pos[i].Y] + 1, next_walll[pos[i].Y][pos[i].X]));
        next_pos[pos[i].Y] = pos[i].X;
    }
}

void fmover_xfast(int n) {
    ZERO(n_order); 
    REP(i, N) next_pos[i] = N;
    REP(i, K) order[N-1-pos[i].X][n_order[N-1-pos[i].X]++] = i;
    REP(j, N) REP(k, n_order[j]) {
        int i = order[j][k];  
        pos[i].X = min(pos[i].X + n, min(next_pos[pos[i].Y] - 1, next_wallr[pos[i].Y][pos[i].X]));
        next_pos[pos[i].Y] = pos[i].X;
    }
}


void fmoveu_markwall() {
    ZERO(n_order); 
    REP(i, K) order[pos[i].Y][n_order[pos[i].Y]++] = i;
    REP(j, N) REP(k, n_order[j]) {int i = order[j][k]; wallh_mark[pos[i].Y][pos[i].X]=1; if (!wallh[pos[i].Y][pos[i].X] && !cell[pos[i].Y - 1][pos[i].X]) {
        cell[pos[i].Y][pos[i].X] = 0;
        pos[i].Y--;
        cell[pos[i].Y][pos[i].X] = 1;
    }}
}

void fmoved_markwall() {
    ZERO(n_order);
    REP(i, K) order[N-1-pos[i].Y][n_order[N-1-pos[i].Y]++] = i;
    REP(j, N) REP(k, n_order[j]) {int i = order[j][k]; wallh_mark[pos[i].Y+1][pos[i].X]=1; if (!wallh[pos[i].Y + 1][pos[i].X] && !cell[pos[i].Y + 1][pos[i].X]) {
        cell[pos[i].Y][pos[i].X] = 0;
        pos[i].Y++;
        cell[pos[i].Y][pos[i].X] = 1;
    }}
}

void fmovel_markwall() {
    ZERO(n_order); 
    REP(i, K) order[pos[i].X][n_order[pos[i].X]++] = i;
    REP(j, N) REP(k, n_order[j]) {int i = order[j][k]; wallv_mark[pos[i].Y][pos[i].X]=1; if (!wallv[pos[i].Y][pos[i].X] && !cell[pos[i].Y][pos[i].X - 1]) {
        cell[pos[i].Y][pos[i].X] = 0;
        pos[i].X--;
        cell[pos[i].Y][pos[i].X] = 1;
    }}
}

void fmover_markwall() {
    ZERO(n_order); 
    REP(i, K) order[N-1-pos[i].X][n_order[N-1-pos[i].X]++] = i;
    REP(j, N) REP(k, n_order[j]) {int i = order[j][k]; wallv_mark[pos[i].Y][pos[i].X+1]=1; if (!wallv[pos[i].Y][pos[i].X + 1] && !cell[pos[i].Y][pos[i].X + 1]) {
        cell[pos[i].Y][pos[i].X] = 0;
        pos[i].X++;
        cell[pos[i].Y][pos[i].X] = 1;
    }}
}

int vs[N][N];
char prv[N][N];
int qdata[N*N*2];
int qst, qen;
PII fp_dest;
string find_path(PII p0, PII p1) {
    if (p0 == p1) return "";

    ZERO(vs);
    ZERO(prv);

    vs[p0.Y][p0.X] = 1;
    prv[p0.Y][p0.X] = -1;
    qdata[0] = p0.X;
    qdata[1] = p0.Y;
    qst = 0;
    qen = 2;

    while (qst < qen) {
        PII p = MP(qdata[qst], qdata[qst + 1]);
        qst += 2;
        if (p == p1) break;

        if (!wallh[p.Y][p.X] && vs[p.Y - 1][p.X] == 0 && cell[p.Y - 1][p.X] == 0)  {
            vs[p.Y - 1][p.X] = vs[p.Y][p.X] + 1;
            prv[p.Y - 1][p.X] = 'U';
            qdata[qen++] = p.X;
            qdata[qen++] = p.Y - 1;
        }
        if (!wallh[p.Y + 1][p.X] && vs[p.Y + 1][p.X] == 0 && cell[p.Y + 1][p.X] == 0) {
            vs[p.Y + 1][p.X] = vs[p.Y][p.X] + 1;
            prv[p.Y + 1][p.X] = 'D';
            qdata[qen++] = p.X;
            qdata[qen++] = p.Y + 1;
        }
        if (!wallv[p.Y][p.X] && vs[p.Y][p.X - 1] == 0 && cell[p.Y][p.X - 1] == 0) {
            vs[p.Y][p.X - 1] = vs[p.Y][p.X] + 1;
            prv[p.Y][p.X - 1] = 'L';
            qdata[qen++] = p.X - 1;
            qdata[qen++] = p.Y;
        }
        if (!wallv[p.Y][p.X + 1] && vs[p.Y][p.X + 1] == 0 && cell[p.Y][p.X + 1] == 0) {
            vs[p.Y][p.X + 1] = vs[p.Y][p.X] + 1;
            prv[p.Y][p.X + 1] = 'R'; 
            qdata[qen++] = p.X + 1;
            qdata[qen++] = p.Y;
        }
    }

    // if (vs[p1.Y][p1.X] == 0) return "";

    if (vs[p1.Y][p1.X] == 0) {
        PII bp = MP(0, 0);
        int bv = 1000000;
        REP(r, N) REP(c, N) if (vs[r][c]) {
            int av = (abs(r - p1.Y) + abs(c - p1.X)) * 100 - vs[r][c];
            if (av <= bv) {
                bv = av;
                bp = MP(c, r);
            }
        }
        p1 = bp;
        assert(vs[p1.Y][p1.X] > 0);
    }

    fp_dest = p1;

    string path = "";
    while (p1 != p0) {
        // DB(p1, prv[p1.Y][p1.X]);
        assert(prv[p1.Y][p1.X]);
        path += prv[p1.Y][p1.X];
        if (prv[p1.Y][p1.X] == 'U') p1.Y++;
        else if (prv[p1.Y][p1.X] == 'D') p1.Y--;
        else if (prv[p1.Y][p1.X] == 'L') p1.X++;
        else if (prv[p1.Y][p1.X] == 'R') p1.X--;
        else assert(false);
    }

    reverse(ALL(path));

    return path;
}


int main(int argc, char **argv) {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);

    cin >> _ >> K;
    REP(i, K) {
        cin >> src[i].Y >> src[i].X;
        cin >> dst[i].Y >> dst[i].X;
    }

    REP(r, N) {
        string s; cin >> s;
        REP(c, N) owallv[r][c+1] = (s[c] == '1');
    }
    
    REP(r, N-1) {
        string s; cin >> s;
        REP(c, N) owallh[r+1][c] = (s[c] == '1');
    }
    
    int W = 0;
    REP(r, N) REP(c, N) {
        if (owallv[r][c+1] && (r == 0 || !owallv[r-1][c+1])) W++;
        if (owallh[r+1][c] && (c == 0 || !owallh[r+1][c-1])) W++;
    }

    DATA(W);

    REP(i, N) owallv[i][0] = owallv[i][N] = 1;
    REP(i, N) owallh[0][i] = owallh[N][i] = 1;

    DATA(K);

    int total_dist = 0;
    REP(i, K) {
        total_dist += abs(src[i].X - dst[i].X) + abs(src[i].Y - dst[i].Y);
    }
    DATA(total_dist);

    int optimal = 0;
    int maxu = 0;
    int maxd = 0;
    int maxl = 0;
    int maxr = 0;
    REP(i, K) {
        maxu = max(maxd, dst[i].Y - src[i].Y);
        maxd = max(maxu, src[i].Y - dst[i].Y);
        maxl = max(maxr, dst[i].X - src[i].X);
        maxr = max(maxl, src[i].X - dst[i].X);
    }

    int change = K < 33 ? -2 : -1;

    maxu += change;
    maxd += change;
    maxl += change;
    maxr += change;
    

    optimal = maxu + maxd + maxl + maxr;
    DATA(optimal);

    int bv = 1e9;
    memcpy(wallv, owallv, sizeof(wallv));
    memcpy(wallh, owallh, sizeof(wallh));

    MINUS(wallv_mark);
    MINUS(wallh_mark);

    int step = 0;

    double TIME_SCALE = 1.0;

    double CUTOFF = 1.85278; // opt 1.6-1.9

    const double TIME_LIMIT = CUTOFF * TIME_SCALE;
    const double TIME_LIMIT_BFS = 1.94 * TIME_SCALE;

    const int ttype = K > 55;

    const double t0 = ttype ? 27.46494 : 12.51129; //opt 0.1-100.0 log
    const double tn = ttype ? 0.01022 : 0.01347; //opt 0.0001-0.1 log
    double t = t0;

    const double tempo = ttype ? 2.8584 : 1.15281; // opt 0.25-5.0 log
    const double removed_factor = ttype ? 0.05508 : 0.11375; // opt 0.01-1.0 log

    double time_passed = 1e-9;

    REP(i, N) {
        rebuild_next_wall_col(i);
        rebuild_next_wall_row(i);
    }

    while (true) {
        step++;
        if ((step & 511) == 0) {
            time_passed = elapsed() / TIME_LIMIT;
            if (time_passed > 1.0) break;
            t = t0 * pow(tn / t0, pow(time_passed, tempo));
        }

        int type = rng.next(2);
        // if (time_passed > .9) type = rng.next(4); // good for high K?
        int r, c;
        char dir;
        int id;
        int ppos = -1;

        bool removed = false;

        if (type == 0) {
            r = rng.next(N);
            c = rng.next(N-1);
            if (owallv[r][c+1]) continue;
            wallv[r][c+1] = 1 - wallv[r][c+1];
            removed |= wallv[r][c+1] == 0;
            rebuild_next_wall_row(r);
        } else if (type == 1){
            r = rng.next(N-1);
            c = rng.next(N);
            if (owallh[r+1][c]) continue;
            wallh[r+1][c] = 1 - wallh[r+1][c];
            removed |= wallh[r+1][c] == 0;
            rebuild_next_wall_col(c);
        } 

        fast_reset();
        fmoveu_xfast(maxu / 2);
        fmovel_xfast(maxl / 2);
        fmoved_xfast(maxd / 2);
        fmover_xfast(maxr);
        fmoved_xfast(maxd - maxd / 2);
        fmovel_xfast(maxl - maxl / 2);
        fmoveu_xfast(maxu - maxu / 2);

        int av = 0; 
        REP(i, K) av += abs(pos[i].X - dst[i].X) + abs(pos[i].Y - dst[i].Y);


        // if (av < bv || (removed || rng.next_double() < removed_factor) && (rng.next_double() < exp((bv - av) / t))) {
        if (av < bv || 
            ttype == 1 && (removed || rng.next_double() < removed_factor) && (av < bv + rng.next_double() * t) ||
            ttype == 0 && (removed || rng.next_double() < removed_factor) && (rng.next_double() < exp((bv - av) / t))) {
            if (!SILENT && av < bv) DB(step, av);
            bv = av;
        } else {
            if (type == 0) {
                wallv[r][c+1] = 1 - wallv[r][c+1];
                rebuild_next_wall_row(r);
            } else if (type == 1) {
                wallh[r+1][c] = 1 - wallh[r+1][c];
                rebuild_next_wall_col(c);
            } 
        }

    }

    DATA(bv);

    DATA(step);

    REP(loop, 2) {
        reset();

        ZERO(wallv_mark);
        ZERO(wallh_mark);
        REP(i, maxu / 2) fmoveu_markwall();
        REP(i, maxl / 2) fmovel_markwall();
        REP(i, maxd / 2) fmoved_markwall();
        REP(i, maxr) fmover_markwall();
        FOR(i, maxd / 2, maxd) fmoved_markwall();
        FOR(i, maxl / 2, maxl) fmovel_markwall();
        FOR(i, maxu / 2, maxu) fmoveu_markwall();
        int walls_removed = 0;
        REP(r, N+2) REP(c, N+2) {
            if (!wallv_mark[r][c] && wallv[r][c] && !owallv[r][c]) {
                wallv[r][c] = 0;
                walls_removed++;
            }
            if (!wallh_mark[r][c] && wallh[r][c] && !owallh[r][c]) {
                wallh[r][c] = 0;
                walls_removed++;
            }
        }
    }

    REP(r, N) {
        REP(c, N-1) cout << wallv[r][c+1];
        cout << endl;
    }

    REP(r, N-1) {
        REP(c, N) cout << wallh[r+1][c];
        cout << endl;
    }

    REP(i, K) cout << 0 << ' ';
    cout << endl;

    REP(i, maxu/2) cout << "g 0 U" << endl;
    REP(i, maxl/2) cout << "g 0 L" << endl;
    REP(i, maxd/2) cout << "g 0 D" << endl;
    REP(i, maxr) cout << "g 0 R" << endl;
    FOR(i, maxd/2,maxd) cout << "g 0 D" << endl;
    FOR(i, maxl/2,maxl) cout << "g 0 L" << endl;
    FOR(i, maxu/2,maxu) cout << "g 0 U" << endl;
    
    
    int ex = optimal;
    
    DB(elapsed());

    state_save();
    
    VI bfs_order;
    int best_bfs_value = optimal;

    int last_good = 0;
    REP(i, K) if (pos[i] == dst[i]) bfs_order.PB(i);

    while (last_good < K) {
        VI order(K); REP(i, K) order[i] = i;
        VI order_value(K);
        REP(i, K) {
            order_value[i] = abs(pos[i].X - dst[i].X) + abs(pos[i].Y - dst[i].Y);
        }
        sort(ALL(order), [&](int a, int b) {return order_value[a] > order_value[b];});

        int good = 0;
        for (int i : order) {
            if (pos[i] == dst[i]) {
                good++;
                continue;
            }   
            
            string path = find_path(pos[i], dst[i]);
            if (path.S == 0) continue;
            best_bfs_value += path.S;
            bfs_order.PB(i);
            // for (char c : path) cout << "i " << i << ' ' << c << endl;
            cell[pos[i].Y][pos[i].X] = 0;
            pos[i] = fp_dest;
            cell[pos[i].Y][pos[i].X] = 1;
            good++;
        }

        if (good == last_good) break;
        last_good = good;

        DB(good);
    }

    REP(i, K) if (pos[i] != dst[i]) {
        bfs_order.PB(i);
        best_bfs_value += (abs(pos[i].X - dst[i].X) + abs(pos[i].Y - dst[i].Y)) * 100;
    }

    DB(best_bfs_value);

    VI bad_pos;

    int bfs_step = 0;

    double bad_pos_prob = 0.64389; // opt 0.1-0.9

    while (elapsed() < TIME_LIMIT_BFS) {
        bfs_step++;
        VI new_bfs_order = bfs_order;
        int pos1 = rng.next(K);
        int pos2 = rng.next(K);
        if (bad_pos.S && rng.next_double() < bad_pos_prob) {
            pos1 = bad_pos[rng.next(bad_pos.S)];
            // pos2 = min(pos2, rng.next(K - 1));
            // if (rng.next_double() < .1) pos2 = -1;
        }
        int xxx = new_bfs_order[pos1];
        new_bfs_order.erase(new_bfs_order.begin() + pos1);
        if (pos2 >= 0) new_bfs_order.insert(new_bfs_order.begin() + pos2, xxx);

        int av = optimal;

        state_load();
        REP(xyz, new_bfs_order.S) {
            if (pos2 == -1) {
                int i = xxx;
                string path = find_path(pos[i], dst[i]);
                if (fp_dest != dst[i]) continue;
                av += path.S;
                cell[pos[i].Y][pos[i].X] = 0;
                pos[i] = fp_dest;
                cell[pos[i].Y][pos[i].X] = 1;
                new_bfs_order.insert(new_bfs_order.begin() + xyz, i);
                pos2 = 0;
                continue;
            }
            
            int i = new_bfs_order[xyz];
            string path = find_path(pos[i], dst[i]);
            if (path.S == 0) continue;
            av += path.S;
            cell[pos[i].Y][pos[i].X] = 0;
            pos[i] = fp_dest;
            cell[pos[i].Y][pos[i].X] = 1;
        }

        if (new_bfs_order.S < K) continue;

        REP(i, K) av += (abs(pos[i].X - dst[i].X) + abs(pos[i].Y - dst[i].Y)) * 100;

        if (av <= best_bfs_value) {
            best_bfs_value = av;
            bfs_order = new_bfs_order;
            
            bad_pos.clear();
            REP(i, K) if (pos[bfs_order[i]] != dst[bfs_order[i]]) bad_pos.PB(i);
        }
    }
    DATA(bfs_step);

    DB(best_bfs_value);

    state_load();
    ex = optimal;

    for (int i : bfs_order) {
        string path = find_path(pos[i], dst[i]);
        if (path.S == 0) continue;
        for (char c : path) cout << "i " << i << ' ' << c << endl;
        ex += path.S;
        cell[pos[i].Y][pos[i].X] = 0;
        pos[i] = fp_dest;
        cell[pos[i].Y][pos[i].X] = 1;
    }

    int failed = 0;
    REP(i, K) {
        if (pos[i] != dst[i]) {
            failed = 1;
            ex += (abs(pos[i].X - dst[i].X) + abs(pos[i].Y - dst[i].Y)) * 100;
        }
    }


    DATA(failed);
    DATA(ex);


    double time = elapsed();
    DATA(time);
	
	return 0;
}