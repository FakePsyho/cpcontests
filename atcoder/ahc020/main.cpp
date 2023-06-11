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
    INLINE int next(int x) {return rand() % x; }
    INLINE int next(int a, int b) {return a + (rand() % (b - a)); }
    INLINE double next_double() {return (rand() + 0.5) * (1.0 / 4294967296.0); }
    INLINE double next_double(double a, double b) {return a + next_double() * (b - a); }
};
 
static RNG rng;


// SOLUTION

const double TIME_LIMIT = 1.9;

const int N = 100;
const int MAXM = 300;
const int MAXK = 5000;

int M, K;

PII vpos[N];
VPII ve[N];

int edge_id[N][N];
int edge_cost[MAXM];
PII edges[MAXM];

VPII kpos;

VI kclosest[MAXK];
VI vclosest[N];

int dist(PII a, PII b) {
    int dx = a.X - b.X;
    int dy = a.Y - b.Y;
    return dx * dx + dy * dy;
}

VI mst(VI need) {
    priority_queue<PII, VPII, greater<PII>> pq;

    assert(need[0] == 1);

    need[0] = 0;

    VI rv(M, 0);

    REP(i, ve[0].S) pq.push(MP(ve[0][i].Y, edge_id[0][ve[0][i].X]));

    int mst_cost = 0;

    double t1 = elapsed();

    while (!pq.empty()) {
        PII p = pq.top(); pq.pop();
        int c = p.X;
        int e = p.Y;
        if (need[edges[e].X] == 0 && need[edges[e].Y] == 0) continue;
        int v = need[edges[e].X] == 0 ? edges[e].Y : edges[e].X;
        rv[e] = 1;
        need[v] = 0;
        mst_cost += c;

        for (PII &r : ve[v]) if (need[r.X]) pq.push(MP(r.Y, edge_id[v][r.X]));
    }

    bool ok = true;
    REP(i, N) if (need[i]) ok = false;
    if (!ok) return VI(M, 1);

    return rv;
}

void find_kclosest(int no) {
    REP(i, K) {
        VPII vp;
        REP(j, N) vp.PB(MP(dist(kpos[i], vpos[j]), j));
        sort(ALL(vp));
        kclosest[i] = VI();
        REP(j, no) {
            kclosest[i].PB(vp[j].Y);
            vclosest[vp[j].Y].PB(i);
        }
    }

}

VI kgreedy(VI vs, VI power, bool init=false) {
    #define PX pair<int, PII>
    priority_queue<PX, VC<PX>, greater<PX>> pq;

    VI kvs(K);

    VI EMPTY_PENALTY = VI(N, 1<<30);
    REP(i, M) {
        EMPTY_PENALTY[edges[i].X] = min(EMPTY_PENALTY[edges[i].X], edge_cost[i]);
        EMPTY_PENALTY[edges[i].Y] = min(EMPTY_PENALTY[edges[i].Y], edge_cost[i]);
    }
    if (init) {
        REP(j, N) EMPTY_PENALTY[j] *= .2;
    } else {
        double mul = pow(0.5, rng.next_double() * 5);
        REP(j, N) EMPTY_PENALTY[j] *= mul;
    }

    REP(i, K) {
        int closest = 1<<30;
        int bv = -1;
        for (int v : kclosest[i]) if (vs[v]) {
            int d = dist(kpos[i], vpos[v]) - power[v] * power[v];
            if (power[v] == 0) d += EMPTY_PENALTY[v];
            if (d < closest) {
                closest = d;
                bv = v;
            }
        }
        if (closest <= 0) kvs[i] = 1;
        else pq.push(MP(closest, MP(i, bv)));
    }

    while (!pq.empty()) {
        PX p = pq.top(); pq.pop();
        int k = p.Y.X;
        int v = p.Y.Y;
        int d = p.X;
        if (kvs[k]) continue;
        kvs[k] = 1;
        int d2 = dist(kpos[k], vpos[v]);
        power[v] = max(power[v], (int)sqrt(d2) + 1);
        for (int i : vclosest[v]) if (!kvs[i]) {
            int d3 = dist(kpos[i], vpos[v]);
            pq.push(MP(d3 - power[v] * power[v] + (power[v] ? 0 : EMPTY_PENALTY[v]), MP(i, v)));
        }

    }

    bool ok = true;
    REP(i, K) if (!kvs[i]) ok = false;
    REP(i, N) if (power[i] > 5000) ok = false;
    if (!ok) return VI(N, 5000);

    return power;
}


int main(int argc, char **argv) {
    int _;
    cin >> _ >> M >> K;

    REP(i, N) cin >> vpos[i].X >> vpos[i].Y;

    DB(N, M, K);

    REP(i, M) {
        int a, b, w;
        cin >> a >> b >> w;
        a--;
        b--;
        ve[a].PB(MP(b, w));
        ve[b].PB(MP(a, w));
        edge_id[a][b] = edge_id[b][a] = i;
        edges[i] = MP(a, b);
        edge_cost[i] = w;
    }

    kpos = VPII(K);
    REP(i, K) cin >> kpos[i].X >> kpos[i].Y;


    cerr << "[DATA] M = " << M << endl;
    cerr << "[DATA] K = " << K << endl;

    find_kclosest(10);

    VI need(N, 1);

    VI power = kgreedy(VI(N, 1), VI(N, 0), true);
    LL ncost = 0; for (int x : power) ncost += x * x;
    VI eused = mst(need);
    LL wcost = 0; REP(i, M) if (eused[i]) wcost += edge_cost[i];

    int step1 = 0;
    while (elapsed() < TIME_LIMIT) {
        step1++;
        int a;
        do {a = rng.next(N);} while (power[a] == 0);
        int val = rng.next(2) == 0 ? 0 : rng.next(power[a]);
        if (rng.next(10) == 0) rng.next(power[a] * 2);

        VI npower = power;
        npower[a] = val;
        npower = kgreedy(VI(N, 1), npower);

        LL nncost = 0; for (int x : npower) nncost += x * x;
        VI need(N, 1);
        VI neused = mst(need);
        LL nwcost = 0; REP(i, M) if (neused[i]) nwcost += edge_cost[i];

        REP(xxx, 5) {
            bool update = false;
            FOR(i, 1, N) if (npower[i] == 0 && need[i] == 1) {
                need[i] = 0;
                VI v = mst(need);
                LL cost = 0; REP(i, M) if (v[i]) cost += edge_cost[i];
                if (cost <= nwcost) {
                    update |= cost < nwcost;
                    neused = v;
                    nwcost = cost;
                } else {
                    need[i] = 1;
                }
            }

            FOR(i, 1, N) if (need[i] == 0) {
                need[i] = 1;
                VI v = mst(need);
                LL cost = 0; REP(i, M) if (v[i]) cost += edge_cost[i];
                if (cost <= nwcost) {
                    update |= cost < nwcost;
                    neused = v;
                    nwcost = cost;
                } else {
                    need[i] = 0;
                }
            }

            if (!update) break;
        }

        if (nncost + nwcost <= ncost + wcost) {
            if (nncost + nwcost < ncost + wcost) DB(step1, nncost+nwcost, elapsed());
            power = npower;
            eused = neused;
            ncost = nncost;
            wcost = nwcost;
        }
    }
    DB(step1);

    DB(wcost + ncost);
    DB(1e6 * (1 + 1e8 / (wcost + ncost + 1e7)));
    cerr << "[DATA] wcost = " << wcost << endl;
    cerr << "[DATA] ncost = " << ncost << endl;
    cerr << "[DATA] time = " << elapsed() << endl;

    REP(i, N) cout << power[i] << ' ';
    cout << endl;

    REP(i, M) cout << eused[i] << ' ';
    cout << endl;

	return 0;
}
