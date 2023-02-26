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

double elapsed() {return get_time() - start_time;}

// SOLUTION

const int MAX_N = 200;
const int MAX_W = 4;
const int MAX_K = 10;

int P1D(int r, int c) {return r * MAX_N + c;}
PII P2D(int p) {return PII(p / MAX_N, p % MAX_N);}
int dd[] = {-MAX_N, -1, +1, +MAX_N};
int dr[] = {-1, +1, 0, 0};
int dc[] = {0, 0, -1, +1};

int N, W, K, C;
VPII VW, VH;
int CLEV;

int extracted[MAX_N][MAX_N];
int good[MAX_N][MAX_N];
int dig_lo[MAX_N][MAX_N];
int dig_hi[MAX_N][MAX_N];
int dig_last[MAX_N][MAX_N];

bool dig_done = false;
int dig(int r, int c, int v) {
    if (dig_done) {
        extracted[r][c] = 1;
        return 2;
    }
    cout << r << ' ' << c << ' ' << v << endl;
    int rv; cin >> rv;
    dig_last[r][c] = v;
    dig_lo[r][c] = dig_hi[r][c] + 1;
    dig_hi[r][c] += v;
    dig_lo[r][c] = min(dig_lo[r][c], 5000);
    dig_hi[r][c] = min(dig_hi[r][c], 5000);
    if (rv == 1 || rv == 2) extracted[r][c] = 1, dig_lo[r][c] = max(dig_lo[r][c], 10);
    if (rv == 2) DB(r, c, v);
    dig_done |= rv == 2;
    return rv;
}

int find_best_power(int lo, int hi) {
    int exp = (lo + hi) / 2;
    double bv = 1e9;
    int power = 0;
    const double V_WASTE_MUL[] = {1.2, 1.2, 1.6, 2.3, 3.2, 3.2, 4.3, 7.0};
    FOR(p, 1, 1001) {
        double v_waste = p * 0.5;
        double c_waste = (exp + p - 1) / p * C;
        double av = v_waste * V_WASTE_MUL[CLEV] + c_waste;
        if (av < bv) {
            bv = av;
            power = p;
        }
    }
    return power;
}

int edge_cost(int r1, int c1, int r2, int c2) {
    int d = abs(r1 - r2) + abs(c1 - c2);
    if (d == 0) return 0;
    int v1 = extracted[r1][c1] ? (dig_lo[r1][c1] + dig_hi[r1][c1]) / 2 + 10 : 160 + (int)(dig_hi[r1][c1] * (105 - 4 * CLEV) / (pow(log(dig_hi[r1][c1]+1), 2.0) + 1));
    int v2 = extracted[r2][c2] ? (dig_lo[r2][c2] + dig_hi[r2][c2]) / 2 + 10 : 160 + (int)(dig_hi[r2][c2] * (105 - 4 * CLEV) / (pow(log(dig_hi[r2][c2]+1), 2.0) + 1));
    return (d - 1) * (1.0 * C + (v1 + v2) * .5) + (extracted[r1][c1] ? 0 : C) + (extracted[r2][c2] ? 0 : C);
}

map<PII, int> vmap;
VPII vv;
VPII vs;
VVI ve;

void add_edge(PII p1, PII p2) {
    if (vmap.count(p1) == 0) {
        vmap[p1] = vv.S;
        vv.PB(p1);
        ve.PB(VI());
        vs.PB(MP(-1, -1));
    }

    if (vmap.count(p2) == 0) {
        vmap[p2] = vv.S;
        vv.PB(p2);
        ve.PB(VI());
        vs.PB(MP(-1, -1));
    }

    int v1 = vmap[p1];
    int v2 = vmap[p2];

    bool exist = false;
    for (int v : ve[v1]) exist |= v == v2;
    if (exist) return;

    ve[v1].PB(v2);
    ve[v2].PB(v1);
}

void build_graph(int SCAN_DIST, VPII targets) {
    int SCAN_OFFSET = 0;
    int scan_sum = 1;
    while (scan_sum + SCAN_DIST < 200) scan_sum += SCAN_DIST;
    SCAN_OFFSET = (200 - scan_sum) / 2;
    
    vmap.clear();
    vv.clear();
    vs.clear();
    ve.clear();

    for (int r = SCAN_OFFSET; r < 200; r += SCAN_DIST/2) for (int c = SCAN_OFFSET; c < 200; c += SCAN_DIST/2) {
        if (r%SCAN_DIST != c%SCAN_DIST) continue;
        if (r < N-SCAN_DIST) add_edge(MP(r, c), MP(r+SCAN_DIST, c));
        if (c < N-SCAN_DIST) add_edge(MP(r, c), MP(r, c+SCAN_DIST));
        if (r < N-SCAN_DIST/2 && c > SCAN_DIST/2) add_edge(MP(r, c), MP(r+SCAN_DIST/2, c-SCAN_DIST/2));
        if (r < N-SCAN_DIST/2 && c < 200-SCAN_DIST/2) add_edge(MP(r, c), MP(r+SCAN_DIST/2, c+SCAN_DIST/2));
        for (PII &target : targets)
            if (abs(r-target.X)<SCAN_DIST && abs(c-target.Y)<SCAN_DIST && abs(r-target.X)+abs(c-target.Y)>0) add_edge(MP(r, c), target);

        PII bstart = MP(-1, -1);
        int bv = 1e9;
        FOR(nr, max(0, r-SCAN_DIST+1), min(N, r+SCAN_DIST)) FOR(nc, max(0, c-SCAN_DIST+1), min(N, c+SCAN_DIST)) {
            int dist = abs(r - nr) + abs(c - nc);
            if (dist >= SCAN_DIST) continue;
            bool ok = false;
            ok |= good[nr][nc];
            for (PII &p : VW) ok |= p == MP(nr, nc);
            if (!ok) continue;
            int av = edge_cost(r, c, nr, nc);
            if (av < bv) {
                bv = av;
                bstart = MP(nr, nc);
            }
        }
        vs[vmap[MP(r, c)]] = bstart;
    }
}

int last_find_path_cost = 0;
VPII find_path(VPII targets) {
    static priority_queue<PII, VPII, greater<PII>> pq;
    static int mdist[200*200];
    static int mprev[200*200];

    memset(mdist, 0x3F, sizeof(mdist[0]) * vv.S);
    while (!pq.empty()) pq.pop();

    REP(i, vv.S) if (vs[i].X != -1) {
        int dst = edge_cost(vv[i].X, vv[i].Y, vs[i].X, vs[i].Y);
        pq.push(MP(dst, i));
        mprev[i] = -1;
        mdist[i] = dst;
    }

    PII final_v = MP(-1, -1);
    while (!pq.empty()) {
        PII p = pq.top(); pq.pop();
        int dst = p.X;
        int cur = p.Y;

        if (mdist[cur] < dst) continue;

        for (PII &target : targets) if (cur == vmap[target]) {
            last_find_path_cost = dst;
            final_v = target;
        }
        if (final_v.X != -1) break;

        for (int next : ve[cur]) {
            int ndst = dst + edge_cost(vv[cur].X, vv[cur].Y, vv[next].X, vv[next].Y);
            if (mdist[next] <= ndst) continue;
            mdist[next] = ndst;
            mprev[next] = cur;
            pq.push(MP(ndst, next));
        }
    }

    int cur = vmap[final_v];
    VPII path;
    while (mprev[cur] != -1) {
        path.PB(vv[cur]);
        cur = mprev[cur];
    }
    path.PB(vv[cur]);
    path.PB(vs[cur]);
    reverse(ALL(path));

    return path;
}


int main(int argc, char **argv) {
    ios_base::sync_with_stdio(false);

    cin >> N >> W >> K >> C;
    cerr << "[DATA] N = " << N << endl;
    cerr << "[DATA] W = " << W << endl;
    cerr << "[DATA] K = " << K << endl;
    cerr << "[DATA] C = " << C << endl;

    VW = VPII(W); REP(i, W) cin >> VW[i].X >> VW[i].Y;
    VH = VPII(K); REP(i, K) cin >> VH[i].X >> VH[i].Y;

    CLEV = 0;
    while ((1 << CLEV) < C) CLEV++;
    DB(CLEV);

    DB(VW);
    DB(VH);

    double SCAN_POW = 1.07 + 0.085 * CLEV;

    // Dig houses
    for (PII &p : VH) {
        int power = max(16, C);
        while (!extracted[p.X][p.Y]) {
            dig(p.X, p.Y, power);
            power = min((int)(power * SCAN_POW), 5000 - dig_hi[p.X][p.Y]);
        }
    }        

    while (true) {
        if (dig_done) break;

        VPII houses_left; 
        for (PII &house : VH) if (!good[house.X][house.Y]) houses_left.PB(house);

        build_graph(18, VPII{houses_left});
        while (true) {
            VPII path = find_path(houses_left);

            double bv = -1e9;
            PII target = MP(-1, -1);
            for (PII p : path) if (!extracted[p.X][p.Y]) {
                double av = -dig_hi[p.X][p.Y];
                if (av >= bv) {
                    bv = av;
                    target = p;
                }
            }

            if (target.X == -1) {
                VPII fpath;
                fpath.PB(path[0]);
                for (PII next : path) {
                    while (fpath.back() != next) {
                        PII p = fpath.back();
                        int diff_r = abs(p.X - next.X);
                        int diff_c = abs(p.Y - next.Y);
                        if (diff_c < diff_r) {
                            fpath.PB(MP(p.X + (next.X > p.X ? 1 : -1), p.Y));
                        } else {
                            fpath.PB(MP(p.X, p.Y + (next.Y > p.Y ? 1 : -1)));
                        }
                    }
                }

                PII prev = MP(-1, -1);
                for (PII &cur : fpath) {
                    int power = prev.X == -1 ? 2 * C + 30 : find_best_power(dig_lo[prev.X][prev.Y], dig_hi[prev.X][prev.Y]);
                    int last_v = prev.X == -1 ? 0 : (dig_lo[prev.X][prev.Y] + dig_hi[prev.X][prev.Y]) / 2;
                    int step = 0;
                    while (!extracted[cur.X][cur.Y]) {
                        step++;
                        int cur_power = power;
                        const double LAST_V_MUL[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.025, 1.05, 1.075};
                        if (last_v && last_v > power) {
                            if (step == 1) {
                                cur_power = (int)(last_v * LAST_V_MUL[CLEV]);
                            } else if (step > 1) {
                                cur_power = (int)(power * .45);
                            }
                        }
                        const int MIN_DIG[] = {0, 0, 0, 12, 13, 15, 19, 25};
                        cur_power = max(MIN_DIG[CLEV], cur_power);
                        cur_power = min(cur_power, 5000 - dig_hi[cur.X][cur.Y]);
                        int rv = dig(cur.X, cur.Y, cur_power);
                    }

                    const double ADJUST_RATIO = C <= 2 ? 2.0 : C <= 16 ? 1.5 : C <= 32 ? 1.4 : C <= 64 ? 1.3 : 1.2;
                    FOR(nr, cur.X-1, cur.X+2) FOR(nc, cur.Y-1, cur.Y+2) {
                        if (nr < 0 || nr >= N || nc < 0 || nc >= N || !extracted[nr][nc]) continue;
                        dig_lo[cur.X][cur.Y] = max(dig_lo[cur.X][cur.Y], (int)(dig_lo[nr][nc] / ADJUST_RATIO));
                        dig_hi[cur.X][cur.Y] = min(dig_hi[cur.X][cur.Y], (int)(dig_hi[nr][nc] * 2.1));
                    }

                    good[cur.X][cur.Y] = 1;
                    prev = cur;
                }

                break;
            }

            int SCAN_START[] = {16, 13, 20, 20, 23, 17, 15, 11};
            int p = dig_last[target.X][target.Y] ? min((int)(dig_last[target.X][target.Y] * SCAN_POW), 40*C) : max(16, C)*SCAN_START[CLEV]/16;
            p = min(p, 5000 - dig_hi[target.X][target.Y]);
            dig(target.X, target.Y, p);
        }
    }

    cerr << "[DATA] time = " << elapsed()*1000 << endl;
	return 0;
}
