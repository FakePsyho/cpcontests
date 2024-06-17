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

const double TIME_LIMIT = 1.9;

const int N = 20;
const int N2 = N*N;
PII p2d(int p) {return MP(p/N, p%N);}
int p1d(PII p) {return p.X*N + p.Y;}

char dist[N2][N2];

const int GDIST = 2;
VI gclose[N2];
const int GDIST2 = 1;
VI gclose2[N2];

int grid0[N2];
char grid0wide[N2*2];
int base_score = 0;

VI sol;

int order_size = 0;

double time_passed;

int last_cost;
int last_diff;
int last_n_diff;
template<int OUTPUT> double sim(VI &s) {
    int cost = 0;
    int cur = 0;
    int p = 0;
    if (OUTPUT) {
        REP(i, order_size) {
            PII p0 = p2d(p);
            PII p1 = p2d(s[i]);
            while (p0.X < p1.X) {p0.X++; cout << "D" << endl;}
            while (p0.X > p1.X) {p0.X--; cout << "U" << endl;}
            while (p0.Y < p1.Y) {p0.Y++; cout << "R" << endl;}
            while (p0.Y > p1.Y) {p0.Y--; cout << "L" << endl;}
            cost += dist[p][s[i]] * (100 + cur);
            p = s[i];
            if (grid0[p] > 0) {
                cur += grid0[p];
                cout << "+" << grid0[p] << endl;
            } else {
                if (-grid0[p] > cur) {
                    if (cur) cout << "-" << cur << endl;
                    cur = 0;
                } else {
                    cout << "-" << -grid0[p] << endl;
                    cur += grid0[p];
                }
            }
        }
    } else {
#pragma GCC unroll 32
        REP(i, order_size) {
            cost += dist[p][s[i]] * (100 + cur);
            p = s[i];
            cur += grid0[p];
            cur = max(0, cur);
        }
    }
    cost += base_score;
    last_cost = cost;
    last_diff = cur;
    last_n_diff = 0;
    return 1e9 * base_score / (cost + cur * 100000LL);
}

const double TP_START = 0.2;
double rv(int c, int d, double tp) {
    return 1e9 * base_score / (c + (d * 500 + (d > 0) * 10000) * pow((tp + TP_START) / (1.0 + TP_START), 5.0));
}

double real_rv(int c, int d) {
    return 1e9 * base_score / (c + (d * 10000));
}

int main(int argc, char **argv) {
    int _; cin >> _;
    REP(r, N) REP(c, N) cin >> grid0[r*N+c];
    REP(r, N) REP(c, N) grid0wide[r*N*2+c] = grid0[r*N+c];

    REP(i, N2) REP(j, N2) dist[i][j] = abs(p2d(i).X - p2d(j).X) + abs(p2d(i).Y - p2d(j).Y);

    REP(i, N2) base_score += abs(grid0[i]);

    REP(i, N2) order_size += grid0[i] != 0;

    VI border;
    VI xorder;

    VI spiral1;
    VI spiral2;
    {
        int r = 0, c = -1;
        int dr = 0, dc = 1;
        int len = N;
        int cur = 0;
        REP(i, N2) {
            r += dr;
            c += dc;
            if (grid0[r*N+c]) spiral1.PB(r*N+c);
            if (grid0[c*N+r]) spiral2.PB(c*N+r);
            cur++;
            if (cur == len) {
                if (dr == 0) len--;
                swap(dr, dc);
                dc = -dc;
                cur = 0;
            }
        }
    }

    double spiral1_rv = sim<0>(spiral1);
    double spiral2_rv = sim<0>(spiral2);
    border = spiral1_rv > spiral2_rv ? spiral1 : spiral2;

    VI tacc(4);
    int step = 0;
    int acc = 0;
    int xv = 0;

    double t0 = 5e4;
    double tn = 1.8e4;
    double t = t0;

    sim<0>(border);
    int bcost = last_cost;
    int bdiff = last_diff;

    REP(i, N2) REP(j, N2) if (dist[i][j] <= GDIST && grid0[j] && i != j) gclose[i].PB(j);
    REP(i, N2) REP(j, N2) if (dist[i][j] <= GDIST2 && grid0[j] && i != j) gclose2[i].PB(j);
    REP(i, N2) if (gclose[i].S == 0) gclose[i].PB(i);
    REP(i, N2) if (gclose2[i].S == 0) gclose2[i].PB(i);

    double time_passed_x = 0.0;
    while (true) {
        step++;
        if ((step & 255) == 0) {
            time_passed = elapsed() / TIME_LIMIT;
            time_passed_x = pow(time_passed, 0.1);
            if (time_passed > 1.0) break;
            t = t0 * pow(tn/t0, time_passed);
        }

        int p0 = rng.next(order_size);
        int p1 = rng.next(order_size-1);
        p1 += p1 >= p0;
        int type = min(rng.next(10), 2);

        if (type == 0) {
            if (rng.next_double() < 0.70 * time_passed_x) {
                int v = border[p0];
                int v2 = gclose[v][rng.next(gclose[v].S)];
                REP(i, order_size) if (border[i] == v2) p1 = i;
            }
            swap(border[p0], border[p1]);
        } else if (type == 1) {
            if (rng.next_double() < 0.30 * time_passed_x) {
                int v = border[p0];
                int v2 = gclose[v][rng.next(gclose[v].S)];
                REP(i, order_size) if (border[i] == v2) p1 = i;
            }
            reverse(border.begin() + p0, border.begin() + p1 + 1);
        } else {
            if (rng.next_double() < 0.80 * time_passed_x) {
                int v = border[p0];
                int v2 = gclose[v][rng.next(gclose[v].S)];
                REP(i, order_size) if (border[i] == v2) p1 = i;
            }
            int v = border[p0];
            if (p0 < p1) {
                FOR(i, p0, p1) border[i] = border[i+1];
            } else {
                for (int i = p0; i > p1; i--) border[i] = border[i-1];
            }
            border[p1] = v;
        }

        sim<0>(border);
        double av = rv(last_cost, last_diff, time_passed);
        double bv = rv(bcost, bdiff, time_passed);
        if (av >= bv || rng.next_double() < exp((av - bv) / t)) {
            int real_score = real_rv(last_cost, last_diff);
            if ((int)real_score > xv) {
                xv = real_score;
                xorder = border;
            }
            bcost = last_cost;
            bdiff = last_diff;
            acc++;
            tacc[type]++;
        } else {
            if (type == 0) {
                swap(border[p0], border[p1]);
            } else if (type == 1) {
                reverse(border.begin() + p0, border.begin() + p1 + 1);
            } else {
                int v = border[p1];
                if (p0 < p1) {
                    for (int i = p1; i > p0; i--) border[i] = border[i-1];
                } else {
                    FOR(i, p1, p0) border[i] = border[i+1];
                }
                border[p0] = v;
            } 
        } 
    }

    DB(step, acc);
    DB(tacc);

    cerr << "[DATA] xv = " << xv << endl;
    cerr << "[DATA] acc = " << acc << endl;
    cerr << "[DATA] step = " << step << endl;


    sim<1>(xorder);
    cerr << "[DATA] diff = " << last_diff << endl;

	
	return 0;
}
