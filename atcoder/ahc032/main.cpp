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
double start_time = get_time();
double elapsed() {return get_time() - start_time;}
 
struct RNG {
    uint64_t x = 88172645463325252ULL;
    INLINE uint32_t rand() {x ^= (x << 7); return x ^= (x >> 9);}
    INLINE int next() {return rand(); }
    INLINE int next(int x) {return ((LL)rand() * x) >> 32; }
    INLINE int next(int a, int b) {return a + next(b - a); }
    INLINE double next_double() {return (rand() + 0.5) * (1.0 / 4294967296.0); }
    INLINE double next_double(double a, double b) {return a + next_double() * (b - a); }
};
 
static RNG rng;

#if defined(LOCAL)
const double TIME_SCALE = 1.0;
#elif defined(VM)
const double TIME_SCALE = 1.0;
#else
const double TIME_SCALE = 1.0;
#endif

// SOLUTION

const double TIME_LIMIT = 1.9;
const int MOD = 998244353;
const int N = 9;
const int N2 = 81;
const int M = 21;
const int K = 81;

LL a[N2];
LL b[N2];
LL c[N2];
LL s[M][9];

int n_sol = 0;
uint16_t sol[N2][2];
uint16_t xsol[N2][2];

bool overlap[N2][N2];

void clear_sol() {
    n_sol = 0;
}

void add_sol(int m, int p) {
    sol[n_sol][0] = m;
    sol[n_sol][1] = p;
    n_sol++;
    REP(dr, 3) REP(dc, 3) a[p+dr*N+dc] = (a[p+dr*N+dc] + s[m][dr*3+dc]) % MOD;
}

void save_sol() {
    memcpy(xsol, sol, sizeof(sol));
}

LL xv = 0;

void greedy_block14(int r, int c) {
    int pos = r*N+c;
    LL bv = a[pos] + a[pos+1] + a[pos+2] + a[pos+3];
    VI b(4, 20);
    VI d(4, 20);

    REP(v1, M) REP(d1, 4) {
        LL x11 = a[pos] + (d1 == 0 ? s[v1][0] : 0);
        LL x12 = a[pos+1] + (d1 <= 1 ? s[v1][1-d1] : 0);
        LL x13 = a[pos+2] + (d1 <= 2 ? s[v1][2-d1] : 0);
        LL x14 = a[pos+3] + (d1 ? s[v1][3-d1] : 0);
        REP(v2, M) REP(d2, 4) {
            LL x21 = x11 + (d2 == 0 ? s[v2][0] : 0);
            LL x22 = x12 + (d2 <= 1 ? s[v2][1-d2] : 0);
            LL x23 = x13 + (d2 <= 2 ? s[v2][2-d2] : 0);
            LL x24 = x14 + (d2 ? s[v2][3-d2] : 0);
            REP(v3, M) REP(d3, 4) {
                LL x31 = x21 + (d3 == 0 ? s[v3][0] : 0);
                LL x32 = x22 + (d3 <= 1 ? s[v3][1-d3] : 0);
                LL x33 = x23 + (d3 <= 2 ? s[v3][2-d3] : 0);
                LL x34 = x24 + (d3 ? s[v3][3-d3] : 0);
                REP(v4, M) REP(d4, 4) {
                    LL x = x31 + (d4 == 0 ? s[v4][0] : 0);
                    x += x32 + (d4 <= 1 ? s[v4][1-d4] : 0);
                    x += x33 + (d4 <= 2 ? s[v4][2-d4] : 0);
                    x += x34 + (d4 ? s[v4][3-d4] : 0);
                    if (x > bv) {
                        bv = x;
                        b[0] = v1;
                        b[1] = v2;
                        b[2] = v3;
                        b[3] = v4;
                        d[0] = d1;
                        d[1] = d2;
                        d[2] = d3;
                        d[3] = d4;
                    }
                }
            
            }
        }
    }
    add_sol(b[0], pos+d[0]);
    add_sol(b[1], pos+d[1]);
    add_sol(b[2], pos+d[2]);
    add_sol(b[3], pos+d[3]);
}

void greedy_block22(int r, int c) {
    int pos = r*N+c;
    LL bv = 0;
    VI b(4, 20);
    VI d(4, 0);

    REP(v1, M) REP(d1, 4) {
        LL x11 = a[pos] + (d1 == 0 ? s[v1][0] : 0);
        LL x12 = a[pos+1] + (d1 <= 1 ? s[v1][1-d1] : 0);
        LL x13 = a[pos+N] + (d1 == 0 ? s[v1][3] : d1 == 2 ? s[v1][0] : 0);
        LL x14 = a[pos+N+1] + (d1 <= 1 ? s[v1][4-d1] : s[v1][3-d1]);
        REP(v2, M) REP(d2, 3) {
            LL x21 = x11 + (d2 == 0 ? s[v2][0] : 0);
            LL x22 = x12 + (d2 <= 1 ? s[v2][1-d2] : 0);
            LL x23 = x13 + (d2 == 0 ? s[v2][3] : d2 == 2 ? s[v2][0] : 0);
            LL x24 = x14 + (d2 <= 1 ? s[v2][4-d2] : s[v2][3-d2]);
            REP(v3, M) REP(d3, 2) {
                LL x31 = x21 + (d3 == 0 ? s[v3][0] : 0);
                LL x32 = x22 + (d3 <= 1 ? s[v3][1-d3] : 0);
                LL x33 = x23 + (d3 == 0 ? s[v3][3] : d3 == 2 ? s[v3][0] : 0);
                LL x34 = x24 + (d3 <= 1 ? s[v3][4-d3] : s[v3][3-d3]);
                REP(v4, M) REP(d4, 1) {
                    LL x = (x31 + (d4 == 0 ? s[v4][0] : 0)) % MOD;
                    x += (x32 + (d4 <= 1 ? s[v4][1-d4] : 0)) % MOD;
                    x += (x33 + (d4 == 0 ? s[v4][3] : d4 == 2 ? s[v4][0] : 0)) % MOD;
                    x += (x34 + (d4 <= 1 ? s[v4][4-d4] : s[v4][3-d4])) % MOD;
                    if (x > bv) {
                        bv = x;
                        b[0] = v1;
                        b[1] = v2;
                        b[2] = v3;
                        b[3] = v4;
                        d[0] = d1;
                        d[1] = d2;
                        d[2] = d3;
                        d[3] = d4;
                    }
                }
            
            }
        }
    }
    
    REP(i, 4) add_sol(b[i], pos+(d[i]/2)*N+d[i]%2);
}

void greedy_block13(int r, int c, LL randomness = 0) {
    int pos = r*N+c;
    LL bv = a[pos] + a[pos+1] + a[pos+2];
    VI b(3, 20);
    VI d(3, 0);
    REP(i, M) REP(j, M) REP(k, M) REP(d1, 3) REP(d2, 3) REP(d3, 1) {
        LL x = 0;
        x += (a[pos+0] + (d1 == 0 ? s[i][0] : 0) + (d2 == 0 ? s[j][0] : 0) + (d3 == 0 ? s[k][0] : 0)) % MOD;
        x += (a[pos+1] + (d1 <= 1 ? s[i][1-d1] : 0) + (d2 <= 1 ? s[j][1-d2] : 0) + (d3 <= 1 ? s[k][1-d3] : 0)) % MOD;
        x += (a[pos+2] + s[i][2-d1] + s[j][2-d2] + s[k][2-d3]) % MOD;
        if (randomness) x += rng.next(randomness);
        if (x > bv) {
            bv = x;
            b = {i, j, k};
            d = {d1, d2, d3};
        }
    }
    REP(i, 3) add_sol(b[i], pos+d[i]);
}

void greedy_block13_lastcol(int r, int c, LL randomness = 0) {
    int pos = r*N+c;
    LL bv = a[pos] + a[pos+1] + a[pos+2];
    VI b(3, 20);
    VI d(3, 0);
    REP(i, M) REP(j, M) REP(k, M) {
        int d1 = 0, d2 = 0, d3 = 0;
        LL x = 0;
        x += (a[pos+0] + (d1 == 0 ? s[i][0] : 0) + (d2 == 0 ? s[j][0] : 0) + (d3 == 0 ? s[k][0] : 0)) % MOD;
        x += (a[pos+1] + (d1 <= 1 ? s[i][1-d1] : 0) + (d2 <= 1 ? s[j][1-d2] : 0) + (d3 <= 1 ? s[k][1-d3] : 0)) % MOD;
        x += (a[pos+2] + s[i][2-d1] + s[j][2-d2] + s[k][2-d3]) % MOD;
        if (randomness) x += rng.next(randomness);
        if (x > bv) {
            bv = x;
            b = {i, j, k};
            d = {d1, d2, d3};
        }
    }
    REP(i, 3) add_sol(b[i], pos+d[i]);
}

void greedy_block23_lastcol(int r, int c, LL randomness = 0) {
    int pos = r*N+c;
    LL bv = 0;
    VI b(6, 20);
    VI d(6, 0);
    REP(v0, M) REP(v1, v0+1) REP(v2, v1+1) REP(u0, M) REP(u1, u0+1) REP(u2, u1+1) {
        LL x = 0;
        x += (a[pos+0] + s[v0][0] + s[v1][0] + s[v2][0]) % MOD;
        x += (a[pos+1] + s[v0][1] + s[v1][1] + s[v2][1]) % MOD;
        x += (a[pos+2] + s[v0][2] + s[v1][2] + s[v2][2]) % MOD;
        x += (a[pos+N+0] + s[v0][3] + s[v1][3] + s[v2][3] + s[u0][0] + s[u1][0] + s[u2][0]) % MOD;
        x += (a[pos+N+1] + s[v0][4] + s[v1][4] + s[v2][4] + s[u0][1] + s[u1][1] + s[u2][1]) % MOD;
        x += (a[pos+N+2] + s[v0][5] + s[v1][5] + s[v2][5] + s[u0][2] + s[u1][2] + s[u2][2]) % MOD;
        if (randomness) x += rng.next(randomness);
        if (x > bv) {
            bv = x;
            b = {v0, v1, v2, u0, u1, u2};
            d = {0, 0, 0, N, N, N};
        }
    }
    REP(i, 6) add_sol(b[i], pos+d[i]);
}

void greedy_block31_lastrow(int r, int c) {
    assert(r == N-3);
    int pos = r*N+c;

    LL bv = a[pos]+a[pos+N]+a[pos+2*N];
    VI b(3, 20);
    REP(i, M) REP(j, i+1) REP(k, j+1) {
        LL x = 0;
        x += (a[pos] + s[i][0] + s[j][0] + s[k][0]) % MOD;
        x += (a[pos+N] + s[i][3] + s[j][3] + s[k][3]) % MOD;
        x += (a[pos+2*N] + s[i][6] + s[j][6] + s[k][6]) % MOD;
        if (x > bv) {
            bv = x;
            b[0] = i;
            b[1] = j;
            b[2] = k;
        }
    }
    add_sol(b[0], pos);
    add_sol(b[1], pos);
    add_sol(b[2], pos);
}

void greedy_block32_lastrow(int r, int c) {
    assert(r == N-3);
    int pos = r*N+c;

    LL bv = a[pos]+a[pos+N]+a[pos+2*N];
    VI b(6, 20);
    VI d(6, 0);
    REP(v0, M) REP(v1, v0+1) REP(v2, v1+1) REP(u0, M) REP(u1, u0+1) REP(u2, u1+1) {
        LL x = 0;
        x += (a[pos] + s[v0][0] + s[v1][0] + s[v2][0]) % MOD;
        x += (a[pos+N] + s[v0][3] + s[v1][3] + s[v2][3]) % MOD;
        x += (a[pos+2*N] + s[v0][6] + s[v1][6] + s[v2][6]) % MOD;
        x += (a[pos+1] + s[v0][1] + s[v1][1] + s[v2][1] + s[u0][0] + s[u1][0] + s[u2][0]) % MOD;
        x += (a[pos+N+1] + s[v0][4] + s[v1][4] + s[v2][4] + s[u0][3] + s[u1][3] + s[u2][3]) % MOD;
        x += (a[pos+2*N+1] + s[v0][7] + s[v1][7] + s[v2][7] + s[u0][6] + s[u1][6] + s[u2][6]) % MOD;
        if (x > bv) {
            bv = x;
            b = {v0, v1, v2, u0, u1, u2};
            d = {0, 0, 0, 1, 1, 1};
        }
    }
    REP(v0, M) REP(v1, v0+1) REP(v2, v1+1) REP(v3, v2+1) REP(u0, M) REP(u1, u0+1) {
        LL x = 0;
        x += (a[pos] + s[v0][0] + s[v1][0] + s[v2][0] + s[v3][0]) % MOD;
        x += (a[pos+N] + s[v0][3] + s[v1][3] + s[v2][3] + s[v3][3]) % MOD;
        x += (a[pos+2*N] + s[v0][6] + s[v1][6] + s[v2][6] + s[v3][6]) % MOD;
        x += (a[pos+1] + s[v0][1] + s[v1][1] + s[v2][1] + s[v3][1] + s[u0][0] + s[u1][0]) % MOD;
        x += (a[pos+N+1] + s[v0][4] + s[v1][4] + s[v2][4] + s[v3][4] + s[u0][3] + s[u1][3]) % MOD;
        x += (a[pos+2*N+1] + s[v0][7] + s[v1][7] + s[v2][7] + s[v3][7] + s[u0][6] + s[u1][6]) % MOD;
        if (x > bv) {
            bv = x;
            b = {v0, v1, v2, v3, u0, u1};
            d = {0, 0, 0, 0, 1, 1};
        }
    }
    REP(i, 6) add_sol(b[i], pos+d[i]);
}

void greedy_block33_last(int r, int c) {
    assert(r == N-3 && c == N-3);
    int pos = r*N+c;

    LL bv = 0;
    VI b(9, 20);
    int cnt = 0;
    REP(v1, M) REP(v2, v1+1) REP(v3, v2+1) REP(v4, v3+1) REP(v5, v4+1) REP(v6, v5+1) REP(v7, v6+1) REP(v8, v7+1) REP(v9, v8+1) {
        LL x = 0;
        cnt++;
        x += (a[pos] + s[v1][0] + s[v2][0] + s[v3][0] + s[v4][0] + s[v5][0] + s[v6][0] + s[v7][0] + s[v8][0] + s[v9][0]) % MOD;
        x += (a[pos+1] + s[v1][1] + s[v2][1] + s[v3][1] + s[v4][1] + s[v5][1] + s[v6][1] + s[v7][1] + s[v8][1] + s[v9][1]) % MOD;
        x += (a[pos+2] + s[v1][2] + s[v2][2] + s[v3][2] + s[v4][2] + s[v5][2] + s[v6][2] + s[v7][2] + s[v8][2] + s[v9][2]) % MOD;
        x += (a[pos+N] + s[v1][3] + s[v2][3] + s[v3][3] + s[v4][3] + s[v5][3] + s[v6][3] + s[v7][3] + s[v8][3] + s[v9][3]) % MOD;
        x += (a[pos+N+1] + s[v1][4] + s[v2][4] + s[v3][4] + s[v4][4] + s[v5][4] + s[v6][4] + s[v7][4] + s[v8][4] + s[v9][4]) % MOD;
        x += (a[pos+N+2] + s[v1][5] + s[v2][5] + s[v3][5] + s[v4][5] + s[v5][5] + s[v6][5] + s[v7][5] + s[v8][5] + s[v9][5]) % MOD;
        x += (a[pos+2*N] + s[v1][6] + s[v2][6] + s[v3][6] + s[v4][6] + s[v5][6] + s[v6][6] + s[v7][6] + s[v8][6] + s[v9][6]) % MOD;
        x += (a[pos+2*N+1] + s[v1][7] + s[v2][7] + s[v3][7] + s[v4][7] + s[v5][7] + s[v6][7] + s[v7][7] + s[v8][7] + s[v9][7]) % MOD;
        x += (a[pos+2*N+2] + s[v1][8] + s[v2][8] + s[v3][8] + s[v4][8] + s[v5][8] + s[v6][8] + s[v7][8] + s[v8][8] + s[v9][8]) % MOD;
        if (x > bv) {
            bv = x;
            b = {v1, v2, v3, v4, v5, v6, v7, v8, v9};
        }
    }
    REP(i, 9) add_sol(b[i], pos);
}

LL get_randomness(int r, int c, bool randomize) {
    LL randomness = 0;
    if (randomness && r == 0 && c == 0) randomness = 1<<29;
    // if (randomize && rng.next_double() < 0.2) randomness = 1<<29;
    return randomness;
}

void greedy(int type, bool randomize = false) {
    clear_sol();
    memcpy(a, c, sizeof(a));

    if (type < 2) {
        REP(r, 6) {
            greedy_block13(r, 0, get_randomness(r, 0, randomize));
            greedy_block13(r, 3, get_randomness(r, 3, randomize));
        }
    } else {
        for (int r = 0; r < 6; r += 2) for (int c = 0; c < 6; c += 2) greedy_block22(r, c);
    }

    // REP(r, 6) greedy_block13(r, 6, get_randomness(randomize));

    greedy_block23_lastcol(0, 6, false);
    greedy_block23_lastcol(2, 6, false);
    greedy_block23_lastcol(4, 6, false);

    if (type == 0) {
        REP(c, 6) greedy_block31_lastrow(6, c);
    } else if (type >= 1) {
        greedy_block32_lastrow(6, 0);
        greedy_block32_lastrow(6, 2);
        greedy_block32_lastrow(6, 4);
    }
    greedy_block33_last(6, 6);

    LL av = 0; REP(i, N2) av += a[i];
    if (av > xv) {
        xv = av;
        DB(xv, elapsed());
        save_sol();
    }
}

int main(int argc, char **argv) {
    int _; cin >> _ >> _ >> _;

    REP(i, N2) cin >> a[i];
    REP(i, M-1) REP(j, 9) cin >> s[i][j];

    memcpy(c, a, sizeof(a));

    REP(i, K) sol[i][0] = M-1;
    REP(i, K) xsol[i][0] = M-1;
    xv = 0; REP(i, N2) xv += a[i];

    greedy(2, false);
    DB(elapsed());
    while (elapsed() < 1.7) greedy(1, true);
    DB(elapsed());

    int step = 0;
    int acc = 0;

    REP(i, N2) REP(j, N2) {
        int x1 = i / N;
        int y1 = i % N;
        int x2 = j / N;
        int y2 = j % N;
        overlap[i][j] = abs(x1 - x2) <= 2 && abs(y1 - y2) <= 2;
    }

    double time_passed = 0;
    const int MAX_MOVES = 3;
    static VI v_p(MAX_MOVES);
    static VI v_m(MAX_MOVES);
    static VI v_pos(MAX_MOVES);

    clear_sol();
    memcpy(a, c, sizeof(c));
    REP(i, K) {
        int m = sol[i][0];
        int p = sol[i][1];
        add_sol(m, p);
    }

    LL bv = 0;
    REP(i, N2) bv += a[i];
    LL xv = bv;

    while (true) {
        int n_moves = rng.next(1, MAX_MOVES+1);
        step++;

        if ((step & ((1<<16)-1)) == 0) {
            time_passed = elapsed();
            if (time_passed > TIME_LIMIT) break;
        }

        REP(i, n_moves) {
            while (true) {
                v_p[i] = rng.next(K);
                bool ok = true;
                REP(j, i) ok &= v_p[i] != v_p[j];
                if (ok) break;
            }
            v_m[i] = rng.next(M);
            v_pos[i] = rng.next(N-2)*N + rng.next(N-2);
        }

        memcpy(b, a, sizeof(a));
        REP(i, n_moves) {
            int old_m = sol[v_p[i]][0];
            int old_pos = sol[v_p[i]][1];
            b[old_pos      ] -= s[old_m][0]; b[old_pos      ] += b[old_pos      ] < 0 ? MOD : 0;
            b[old_pos+1    ] -= s[old_m][1]; b[old_pos+1    ] += b[old_pos+1    ] < 0 ? MOD : 0;
            b[old_pos+2    ] -= s[old_m][2]; b[old_pos+2    ] += b[old_pos+2    ] < 0 ? MOD : 0;
            b[old_pos+N    ] -= s[old_m][3]; b[old_pos+N    ] += b[old_pos+N    ] < 0 ? MOD : 0;
            b[old_pos+N+1  ] -= s[old_m][4]; b[old_pos+N+1  ] += b[old_pos+N+1  ] < 0 ? MOD : 0;
            b[old_pos+N+2  ] -= s[old_m][5]; b[old_pos+N+2  ] += b[old_pos+N+2  ] < 0 ? MOD : 0;
            b[old_pos+2*N  ] -= s[old_m][6]; b[old_pos+2*N  ] += b[old_pos+2*N  ] < 0 ? MOD : 0;
            b[old_pos+2*N+1] -= s[old_m][7]; b[old_pos+2*N+1] += b[old_pos+2*N+1] < 0 ? MOD : 0;
            b[old_pos+2*N+2] -= s[old_m][8]; b[old_pos+2*N+2] += b[old_pos+2*N+2] < 0 ? MOD : 0;
        }
        REP(i, n_moves) {
            int m = v_m[i];
            int pos = v_pos[i];
            b[pos      ] += s[m][0]; b[pos      ] -= b[pos      ] >= MOD ? MOD : 0;
            b[pos+1    ] += s[m][1]; b[pos+1    ] -= b[pos+1    ] >= MOD ? MOD : 0;
            b[pos+2    ] += s[m][2]; b[pos+2    ] -= b[pos+2    ] >= MOD ? MOD : 0;
            b[pos+N    ] += s[m][3]; b[pos+N    ] -= b[pos+N    ] >= MOD ? MOD : 0;
            b[pos+N+1  ] += s[m][4]; b[pos+N+1  ] -= b[pos+N+1  ] >= MOD ? MOD : 0;
            b[pos+N+2  ] += s[m][5]; b[pos+N+2  ] -= b[pos+N+2  ] >= MOD ? MOD : 0;
            b[pos+2*N  ] += s[m][6]; b[pos+2*N  ] -= b[pos+2*N  ] >= MOD ? MOD : 0;
            b[pos+2*N+1] += s[m][7]; b[pos+2*N+1] -= b[pos+2*N+1] >= MOD ? MOD : 0;
            b[pos+2*N+2] += s[m][8]; b[pos+2*N+2] -= b[pos+2*N+2] >= MOD ? MOD : 0;
        }
        LL av = 0;
        REP(i, N2) av += b[i];

        if (av >= bv) {
            acc++;
            bv = av;
            REP(i, n_moves) sol[v_p[i]][0] = v_m[i], sol[v_p[i]][1] = v_pos[i];

            memcpy(a, b, sizeof(b));
            if (bv > xv) {
                DB(step, bv, acc, elapsed());
                save_sol();
                xv = bv;
            }
        }
    }

    DB(acc);
    DB(step);
    DB(xv);

    int non_zero = 0;
    REP(i, K) non_zero += xsol[i][0] != 20;
    cout << non_zero << endl;
    REP(i, K) if (xsol[i][0] != 20) cout << (int)xsol[i][0] << ' ' << (int)xsol[i][1]/9 << " " << (int)xsol[i][1]%9 << endl;

	return 0;
}
