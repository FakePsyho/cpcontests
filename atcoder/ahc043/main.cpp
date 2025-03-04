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
#define PSS         pair<short, short>
#define PDD         pair<double, double>
#define PIII        pair<PII, int>
#define VB          VC<byte>
#define VVB         VC<VB>
#define VI          VC<int>
#define VVI         VC<VI>
#define VVVI        VC<VVI>
#define VPII        VC<PII>
#define VVPII       VC<VPII>
#define VPSS        VC<PSS>
#define VVPSS       VC<VPSS>
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

const double PI = 3.14159265358979323846;

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
    INLINE double next_normal() {return sqrt(-2.0 * log(next_double())) * cos(2 * PI * next_double()); }
};
 
static RNG rng;

// SOLUTION

#undef assert
#define assert(x) ;

const bool SILENT = true;
const bool FIXED_SEED = false;

#if defined(LOCAL)
const double TIME_SCALE = 1.0;
#elif defined(VM)
const double TIME_SCALE = 1.3;
#else
const double TIME_SCALE = 1.0;
#endif

const double TIME_MULT = TIME_SCALE;

const double TIME_LIMIT_GREEDY = 0.45 * TIME_MULT;
const double TIME_LIMIT_HC = 2.45 * TIME_MULT;
const double TIME_LIMIT_TOTAL = 2.95 * TIME_MULT;

const int N = 50;
const int NX = 52;
const int MIN_M = 50;
const int MAX_M = 1600;
const int T = 800;

const int COST_STATION = 5000;
const int COST_RAIL = 100;

const int MAX_STATIONS = 800;

Timer timer0;
Timer timer1;

Timer timer_hc_s0;
Timer timer_hc_sall;

const int dp[] = {1, -1, NX, -NX};
const int rail_id_table[] = {-1, -1, -1, 4, -1, 2, 3, -1, -1, 5, 1, -1, 6, -1, -1, -1};

INLINE int rail_id(int p0, int p1, int p2) {
    int diffa = p0 - p1;
    int diffb = p2 - p1;
    int edges = 0;
    edges |= 1 << (diffa < -1 ? 0 : diffa > 1 ? 2 : 2 + diffa);
    edges |= 1 << (diffb < -1 ? 0 : diffb > 1 ? 2 : 2 + diffb);
    int rv = rail_id_table[edges];
    assert(rv != -1);
    return rv;
}

INLINE int p1d(int r, int c) {return r * NX + c;}
INLINE int p1d(PII p) {return p.X * NX + p.Y;}
INLINE PII p2d(int p) {return {p / NX, p % NX};}
INLINE int dist(PII a, PII b) {return abs(a.X - b.X) + abs(a.Y - b.Y);}
INLINE int dist(int r0, int c0, int r1, int c1) {return abs(r0 - r1) + abs(c0 - c1);}
INLINE int dist(int p0, int p1) {return abs(p0 / NX - p1 / NX) + abs(p0 % NX - p1 % NX);}

struct Trip {int r0, c0, r1, c1;};
struct Order {int type, p;};
struct SOrder {int p, d;};

int M, K;
Trip trips[MAX_M];
int trip_income[MAX_M];

VI station_trips0[NX*NX];
VI station_trips1[NX*NX];
bool useful_station[NX*NX];
int total_income[NX*NX];

Order solution[T];
int best_result = 0;

void print_solution() {
    REP(t, T) {
        if (solution[t].type == -1) {
            cout << -1 << endl;
        } else {
            PII p = p2d(solution[t].p);
            cout << solution[t].type << " " << p.X-1 << " " << p.Y-1 << endl;
        }
    }
}

struct TripManager {
    char status[MAX_M];
    int income = 0;

    void reset() {
        memset(status, 0, sizeof(char)*M);
        income = 0;
    }

    void copy(const TripManager &tm) {
        memcpy(status, tm.status, sizeof(char)*M);
        income = tm.income;
    }

    void add_station(int p) {
        for (int id : station_trips0[p]) {
            if (status[id] == 2) income += trip_income[id];
            status[id] |= 1;
        }

        for (int id : station_trips1[p]) {
            if (status[id] == 1) income += trip_income[id];
            status[id] |= 2;
        }
    }

    int possible_income(int p) {
        int rv = 0;
        for (int id : station_trips0[p]) if (status[id] == 2) rv += trip_income[id];
        for (int id : station_trips1[p]) if (status[id] == 1) rv += trip_income[id];
        return rv;
    }
};

int eval_build_order(const VVI &border, const bool save = false) {
    static VC<Order> order_rail;
    static VC<Order> order_station;
    int n_order_rail = 0;
    int n_order_station = 0;
    if (save) {
        order_rail.clear();
        order_station.clear();
        REP(t, T) {
            if (solution[t].type == 0) {
                order_station.PB(solution[t]);
            } else if (solution[t].type > 0) {
                order_rail.PB(solution[t]);
            }
            solution[t] = {-1, -1};
        }
    }

    int overbuild = 0;
    int money = K;
    int turn = 0;
    int last_income = 0;
    for (const VI &bo : border) {
        int overbuild_used = min(bo[0], overbuild);
        overbuild -= overbuild_used;

        REP(_, bo[0] - overbuild_used + bo[2]) {
            while (money < COST_RAIL) money += last_income, turn++;
            money += last_income - COST_RAIL;
            if (save) solution[turn] = order_rail[n_order_rail++];
            turn++;
        }

        overbuild += bo[2];

        while (money < COST_STATION) money += last_income, turn++;
        money += bo[1] - COST_STATION;
        if (save) solution[turn] = order_station[n_order_station++];
        turn++;
        
        last_income = bo[1];
    }

    if (turn > T) return 0;
    return money + (T - turn) * last_income;
}

void optimize_build_order() {
    TripManager tm;
    tm.reset();

    int cur_res = best_result;

    VVI border;
    
    int n_walls = 0;
    REP(t, T) {
        Order o = solution[t];
        if (o.type == -1) continue;
        if (o.type > 0) {
            n_walls++;
            continue;
        }

        tm.add_station(o.p);
        border.PB({n_walls, tm.income, 0});
        n_walls = 0;
    }

    assert(cur_res == eval_build_order(border));

    FOR(i, 2, border.S) {
        while (true) {
            if (elapsed() > TIME_LIMIT_TOTAL) break;
            border[i][2]++;
            int res = eval_build_order(border);
            if (res < cur_res) {
                border[i][2]--;
                break;
            }
            cur_res = res;
        }
    }

    eval_build_order(border, true);
}

static int path[N*N];
struct State {
    TripManager tm;
    Order orders[T];
    int stations[MAX_STATIONS];
    int n_stations;
    char map[NX*NX];
    int turn;
    int money;
    
    short qdist[NX*NX];
    
    void reset() {
        tm.reset();
        n_stations = 0;
        turn = 0;
        money = K;
        MINUS(map); 
        REP(i, NX) map[p1d(i,0)] = map[p1d(i,NX-1)] = map[p1d(0,i)] = map[p1d(NX-1,i)] = 10;
    }
    
    void station_bfs(int p0 = -1) {
        static int qdata[N*N];
        MINUS(qdist);
        int qst = 0, qen = 0;
        if (p0 == -1) {
            REP(i, n_stations) qdata[i] = stations[i], qdist[stations[i]] = 0;
            qen = n_stations;
        } else {
            qdata[qen++] = p0;
            qdist[p0] = 0;
        }

        while (qst < qen) {
            int p = qdata[qst++];
            if (map[p+ 1] == -1 && qdist[p+ 1] == -1) {qdist[p+ 1] = qdist[p] + 1; qdata[qen++] = p+ 1;}
            if (map[p- 1] == -1 && qdist[p- 1] == -1) {qdist[p- 1] = qdist[p] + 1; qdata[qen++] = p- 1;}
            if (map[p+NX] == -1 && qdist[p+NX] == -1) {qdist[p+NX] = qdist[p] + 1; qdata[qen++] = p+NX;}
            if (map[p-NX] == -1 && qdist[p-NX] == -1) {qdist[p-NX] = qdist[p] + 1; qdata[qen++] = p-NX;}
        }
    }

    void build_station(int p) {
        assert(money >= COST_STATION || tm.income > 0);
        while (money < COST_STATION) {
            money += tm.income;
            if (turn >= T) return;
            orders[turn++].type = -1;
        }
        if (turn >= T) return;
        tm.add_station(p);
        money -= COST_STATION;
        money += tm.income;
        stations[n_stations++] = p;
        orders[turn++] = {0, p};
        map[p] = 0;
    }

    void build_rail_path(int n_path) {
        FOR(i, 1, n_path - 1) {
            while (money < COST_RAIL) {
                money += tm.income;
                if (turn >= T) return;
                orders[turn++].type = -1;
            }
            if (turn >= T) return;
            money -= COST_RAIL;
            money += tm.income;
            int id = rail_id(path[i-1], path[i], path[i+1]);
            map[path[i]] = id;
            orders[turn++] = {id, path[i]};
        }
    }

    bool bfs_path(int p0, int offset = 0) {
        static int qdata[N*N];
        static int qprev[NX*NX];
        static int qprev_cnt = 0;
        if (qprev_cnt == 0) ZERO(qprev);

        qprev_cnt += 4;

        int qst = 0, qen = 0;
        qdata[qen++] = p0;

        int tp = -1;

        while (qst < qen) {
            int p = qdata[qst++];
            if (map[p] == 0) {
                tp = p;
                break;
            }
            REP(dd, 4) {
                int d = (dd + offset) & 3;
                int pp = p + dp[d];
                if (map[pp] > 0 || qprev[pp] >= qprev_cnt) continue;
                qdata[qen++] = pp;
                qprev[pp] = qprev_cnt + d;
            }
        }

        if (tp == -1) return false;

        int n_path = 0;
        while (tp != p0) {
            path[n_path++] = tp;
            int d = qprev[tp] - qprev_cnt;
            tp -= dp[d];
        }
        path[n_path++] = p0;

        build_rail_path(n_path);

        return true;
    }

    bool bfs_path_quick(int p0, const int offset = 0) {
        static int qdata[N*N];
        static int qprev[NX*NX];
        static int qprev_cnt = 0;
        if (qprev_cnt == 0) ZERO(qprev);

        qprev_cnt += 4;

        int qst = 0, qen = 0;
        qdata[qen++] = p0;

        int tp = -1;

        if (offset == 0) {
            while (qst < qen) {
                int p = qdata[qst++];
                if (map[p] == 0) {tp = p; break;}
                if (map[p+ 1] <= 0 && qprev[p+ 1] < qprev_cnt) qdata[qen++] = p+ 1, qprev[p+ 1] = qprev_cnt + 0;
                if (map[p- 1] <= 0 && qprev[p- 1] < qprev_cnt) qdata[qen++] = p- 1, qprev[p- 1] = qprev_cnt + 1;
                if (map[p+NX] <= 0 && qprev[p+NX] < qprev_cnt) qdata[qen++] = p+NX, qprev[p+NX] = qprev_cnt + 2;
                if (map[p-NX] <= 0 && qprev[p-NX] < qprev_cnt) qdata[qen++] = p-NX, qprev[p-NX] = qprev_cnt + 3;
            }
        } else if (offset == 1) {
            while (qst < qen) {
                int p = qdata[qst++];
                if (map[p] == 0) {tp = p; break;}
                if (map[p- 1] <= 0 && qprev[p- 1] < qprev_cnt) qdata[qen++] = p- 1, qprev[p- 1] = qprev_cnt + 1;
                if (map[p+NX] <= 0 && qprev[p+NX] < qprev_cnt) qdata[qen++] = p+NX, qprev[p+NX] = qprev_cnt + 2;
                if (map[p-NX] <= 0 && qprev[p-NX] < qprev_cnt) qdata[qen++] = p-NX, qprev[p-NX] = qprev_cnt + 3;
                if (map[p+ 1] <= 0 && qprev[p+ 1] < qprev_cnt) qdata[qen++] = p+ 1, qprev[p+ 1] = qprev_cnt + 0;
            }
        } else if (offset == 2) {
            while (qst < qen) {
                int p = qdata[qst++];
                if (map[p] == 0) {tp = p; break;}
                if (map[p+NX] <= 0 && qprev[p+NX] < qprev_cnt) qdata[qen++] = p+NX, qprev[p+NX] = qprev_cnt + 2;
                if (map[p-NX] <= 0 && qprev[p-NX] < qprev_cnt) qdata[qen++] = p-NX, qprev[p-NX] = qprev_cnt + 3;
                if (map[p+ 1] <= 0 && qprev[p+ 1] < qprev_cnt) qdata[qen++] = p+ 1, qprev[p+ 1] = qprev_cnt + 0;
                if (map[p- 1] <= 0 && qprev[p- 1] < qprev_cnt) qdata[qen++] = p- 1, qprev[p- 1] = qprev_cnt + 1;
            }
        } else if (offset == 3) {
            while (qst < qen) {
                int p = qdata[qst++];
                if (map[p] == 0) {tp = p; break;}
                if (map[p-NX] <= 0 && qprev[p-NX] < qprev_cnt) qdata[qen++] = p-NX, qprev[p-NX] = qprev_cnt + 3;
                if (map[p+ 1] <= 0 && qprev[p+ 1] < qprev_cnt) qdata[qen++] = p+ 1, qprev[p+ 1] = qprev_cnt + 0;
                if (map[p- 1] <= 0 && qprev[p- 1] < qprev_cnt) qdata[qen++] = p- 1, qprev[p- 1] = qprev_cnt + 1;
                if (map[p+NX] <= 0 && qprev[p+NX] < qprev_cnt) qdata[qen++] = p+NX, qprev[p+NX] = qprev_cnt + 2;
            }
        }

        if (tp == -1) return false;

        int n_path = 0;
        while (tp != p0) {
            path[n_path++] = tp;
            int d = qprev[tp] - qprev_cnt;
            tp -= dp[d];
        }
        path[n_path++] = p0;

        build_rail_path(n_path);

        return true;
    }

    void find_rail_path(int p0, int p1 = -1) {
        if (p1 == -1) swap(p0, p1);

        if (p0 != -1) station_bfs(p0);
        assert(qdist[p1] > 0);

        int n_path = 0;
        int p = p1;
        assert(qdist[p] > 0);
        assert(map[p] == -1);
        path[n_path++] = p;
        while (qdist[p] > 0) {
            int bd = -1;
            REP(d, 4) {
                int pp = p + dp[d];
                if (qdist[pp] != qdist[p] - 1) continue;
                bd = d;
                break;
            }

            assert(bd != -1);
            p += dp[bd];
            path[n_path++] = p;
        }

        build_rail_path(n_path);
    }

    bool find_rail_path_naive(int p0, int p1, const int offset) {
        int r0 = p0 / NX, c0 = p0 % NX;
        int r1 = p1 / NX, c1 = p1 % NX;
        int n_up = r1 > r0 ? r1 - r0 : 0;
        int n_down = r0 > r1 ? r0 - r1 : 0;
        int n_left = c1 > c0 ? c1 - c0 : 0;
        int n_right = c0 > c1 ? c0 - c1 : 0;

        int n_path = 0;
        path[n_path++] = p1;
        if (offset == 0) {
            while (n_up--  )  {p1 -= NX; path[n_path++] = p1; if (map[p1] > 0) return false;}
            while (n_down--)  {p1 += NX; path[n_path++] = p1; if (map[p1] > 0) return false;}
            while (n_left--)  {p1--    ; path[n_path++] = p1; if (map[p1] > 0) return false;}
            while (n_right--) {p1++    ; path[n_path++] = p1; if (map[p1] > 0) return false;}
        } else if (offset == 1) {
            while (n_down--)  {p1 += NX; path[n_path++] = p1; if (map[p1] > 0) return false;}
            while (n_left--)  {p1--    ; path[n_path++] = p1; if (map[p1] > 0) return false;}
            while (n_right--) {p1++    ; path[n_path++] = p1; if (map[p1] > 0) return false;}
            while (n_up--  )  {p1 -= NX; path[n_path++] = p1; if (map[p1] > 0) return false;}
        } else if (offset == 2) {
            while (n_left--)  {p1--    ; path[n_path++] = p1; if (map[p1] > 0) return false;}
            while (n_right--) {p1++    ; path[n_path++] = p1; if (map[p1] > 0) return false;}
            while (n_up--  )  {p1 -= NX; path[n_path++] = p1; if (map[p1] > 0) return false;}
            while (n_down--)  {p1 += NX; path[n_path++] = p1; if (map[p1] > 0) return false;}
        } else if (offset == 3) {
            while (n_right--) {p1++    ; path[n_path++] = p1; if (map[p1] > 0) return false;}
            while (n_up--  )  {p1 -= NX; path[n_path++] = p1; if (map[p1] > 0) return false;}
            while (n_down--)  {p1 += NX; path[n_path++] = p1; if (map[p1] > 0) return false;}
            while (n_left--)  {p1--    ; path[n_path++] = p1; if (map[p1] > 0) return false;}
        }
        build_rail_path(n_path);
        return true;
    }

    void save_solution(int max_turn = T) {
        REP(t, T) solution[t] = (t < max_turn) ? orders[t] : Order{-1, 0};
    }
};

int greedy_vs[NX*NX];
int greedy_vs_cnt = 1;

int do_greedy(const bool allow_no_income) {
    double bv = 0;
    int bp0 = -1, bp1 = -1;
    
    static State s;
    s.reset();

    double type_ratio = pow(rng.next_double(), .1);
    double useful_stations_bonus = pow(0.5, rng.next_double(1.0, 4.0));
    double corner_bonus = rng.next_double() * 0.02;

    timer0.start();

    greedy_vs_cnt++;
    TripManager tm;
    FOR(r0, 1, N+1) FOR(c0, 1, N+1) {
        int p0 = p1d(r0, c0);
        if (station_trips0[p0].S + station_trips1[p0].S == 0) continue;
        tm.reset();
        tm.add_station(p0);
        greedy_vs_cnt++;

        for (int id : station_trips0[p0]) {
            FOR(r1, max(1, trips[id].r1 - 2), min(N+1, trips[id].r1 + 3)) FOR(c1, max(1, trips[id].c1 - 2), min(N+1, trips[id].c1 + 3)) {
                if (dist(r1, c1, trips[id].r1, trips[id].c1) > 2) continue;
                int p1 = p1d(r1, c1);
                if (greedy_vs[p1] == greedy_vs_cnt) continue;
                greedy_vs[p1] = greedy_vs_cnt;
                
                int cost_a = COST_STATION * 2 + (dist(r0, c0, r1, c1) - 1) * COST_RAIL;
                if (cost_a > K) continue;
                
                int poss_income = tm.possible_income(p1);
                if (poss_income == 0) continue;

                int cost_b = COST_STATION + dist(r0, c0, r1, c1) * COST_RAIL;
                int cost = cost_a * type_ratio + cost_b * (1 - type_ratio);

                double useful_stations = total_income[p0] + total_income[p1] - 2 * poss_income;
                useful_stations *= .2;

                int corner_dist = min(r0-1, N - r0) + min(c0-1, N - c0) + min(r1-1, N - r1) + min(c1-1, N - c1);
    
                double av = (poss_income + useful_stations * useful_stations_bonus) * rng.next_double(0.2, 4.0) / cost * (1.0 + corner_dist * corner_bonus);
                if (av > bv) {
                    bv = av;
                    bp0 = p0;
                    bp1 = p1;
                }
            }
        }
    }
    timer0.stop();

    assert(bp0 != -1);

    s.find_rail_path_naive(bp0, bp1, 0);
    s.build_station(bp0);
    s.build_station(bp1);

    timer1.start();
    int last_res = 0;
    while (true) {
        double bv = -1e9;
        int bp = -1;

        s.station_bfs();
        FOR(r, 1, N+1) FOR(c, 1, N+1) {
            int p = p1d(r, c);

            if (!useful_station[p]) continue;
            if (s.qdist[p] <= 0 && s.map[p] <= 0) continue;

            int poss_income = s.tm.possible_income(p);
            if (!allow_no_income && poss_income == 0) continue;

            int cost_a = COST_STATION + (s.qdist[p] - 1) * COST_RAIL;
            int cost_b = s.qdist[p] * COST_RAIL;
            if (s.map[p] > 0) cost_a = COST_STATION, cost_b = COST_RAIL;
            double cost = cost_a * type_ratio + cost_b * (1 - type_ratio);

            if (s.turn + s.qdist[p] + 1 >= T) continue;

            if (1.0 * poss_income / cost * 3 < bv) continue;

            double useful_stations = 0;
            for (int id : station_trips0[p]) useful_stations += (s.tm.status[id] == 0) * trip_income[id];
            for (int id : station_trips1[p]) useful_stations += (s.tm.status[id] == 0) * trip_income[id];
            useful_stations *= 1.2;

            double av = (poss_income + useful_stations * useful_stations_bonus) * rng.next_double(0.95, 1.05) / cost;

            if (av > bv) {
                bv = av;
                bp = p;
            }
        }

        if (bp == -1) break;

        if (s.map[bp] == -1) s.find_rail_path(bp);
        s.build_station(bp);

        int res = s.money + (T - s.turn) * s.tm.income;
        if (last_res > res) break;
        last_res = res;

        if (res > best_result) {
            best_result = res;
            if (!SILENT) DB(best_result, type_ratio, useful_stations_bonus, s.turn, elapsed());
            s.save_solution(s.turn);
        }
    }
    timer1.stop();

    return last_res;
}

VD normalize_weights(VD w) {
    VD w_sum = w; FOR(i, 1, w.S) w_sum[i] += w_sum[i-1];
    VD w_norm = w_sum; REP(i, w.S) w_norm[i] /= w_sum.back();
    return w_norm;
}

void do_hc(double time_limit) {
    double hc_start = elapsed();

    VC<SOrder> order; REP(t, T) if (solution[t].type == 0) order.PB({solution[t].p, 0});

    TripManager alltm;
    alltm.reset();
    for (SOrder o : order) alltm.add_station(o.p);

    State s;
    s.reset();

    int hc_step = 0;
    VI acc(5);

    const double t0 = M >= 800 ? 12000 : 8000;
    const double tn = 100.0;
    double t = t0;
    double time_passed = 0;

    VD move_w = {1, 5, 5, 5, 1};
    VD move_w_norm = normalize_weights(move_w);
    DB(move_w_norm);

    int cur_res = best_result;

    while (true) {
        if ((hc_step & 15) == 0) {
            time_passed = (elapsed() - hc_start) / time_limit;
            if (time_passed > 1.0) break;
            t = t0 * pow(tn / t0, pow(time_passed, 2.0));
        }

        hc_step++;

        VC<SOrder> no = order;

        double type_r = rng.next_double();
        int type = 0; 
        while (type < move_w.S && type_r > move_w_norm[type]) type++;

        auto rand_pos = [&](int mx) {return int(pow(rng.next_double(), max(1.0, 2.0 - 2 * time_passed)) * mx);};

        if (type == 0) {
            // remove an order
            int pos = rand_pos(no.S);
            no.erase(no.begin() + pos);
        } else if (type == 1) {
            // insert a new order
            int pos = rand_pos(no.S+1);
            int station = -1;
            int tries = 0;
            if (pos && rng.next_double() < .5) {
                do {
                    if (++tries > 20) break;
                    int turn = rng.next(T);
                    if (solution[turn].type <= 0) continue;
                    station = solution[turn].p;
                } while (!useful_station[station] || M < 400 && alltm.possible_income(station) == 0);
            } else if (pos && rng.next_double() < .25) {
                do {
                    if (++tries > 20) break;
                    PII origin = p2d(no[rng.next(pos)].p);
                    double d = abs(rng.next_normal()) * 10;
                    d = max(d, 1.0);
                    double angle = rng.next_double() * 2 * PI;
                    origin.X = int(origin.X + .5 + d * cos(angle));
                    origin.Y = int(origin.Y + .5 + d * sin(angle));
                    if (origin.X < 1 || origin.X > N || origin.Y < 1 || origin.Y > N) continue;
                    station = p1d(origin);
                } while (!useful_station[station] || M < 400 && alltm.possible_income(station) == 0);
            } else {
                do {
                    if (++tries > 20) break;
                    station = p1d({rng.next(1,N+1), rng.next(1,N+1)});
                } while (!useful_station[station] || M < 400 && alltm.possible_income(station) == 0);
            }
            if (tries >= 20 || station < 0) continue;
            no.insert(no.begin() + pos, {station, rng.next(4)});
        } else if (type == 2) {
            // adjust a single order
            int pos = rand_pos(no.S);
            int station = -1;
            int tries = 0;
            if (rng.next_double() < .1) {
                do {
                    if (++tries > 20) break;
                    station = p1d({rng.next(1,N+1), rng.next(1,N+1)});
                } while (!useful_station[station]);
            } else {
                do {
                    if (++tries > 20) break;
                    station = no[pos].p + rng.next(-1, 2) + rng.next(-1, 2) * NX;
                } while (!useful_station[station] || station == no[pos].p || s.map[station] == 10);
            }
            if (tries >= 20) continue;
            no[pos] = {station, no[pos].d};
        } else if (type == 3) {
            // move order
            int pos0 = rand_pos(no.S);
            int pos1;
            if (rng.next_double() < .7) {
                pos1 = rng.next(no.S - 1);
            } else {
                pos1 = rng.next(max(0, pos0 - 10), min((int)no.S - 1, pos0 + 10));
            }
            pos1 += pos1 >= pos0;
            SOrder tmp = no[pos0];
            no.erase(no.begin() + pos0);
            no.insert(no.begin() + pos1, tmp);
        } else if (type == 4) {
            // change path dir priority
            int pos = rand_pos(no.S);
            int old_d = no[pos].d;
            int new_d = rng.next(3);
            new_d += new_d >= old_d;
            no[pos].d = new_d;
        }
        
        if (no.S < 2) continue;
        if (COST_STATION * 2 + (dist(no[0].p, no[1].p) - 1) * COST_RAIL > K) continue;

        if (no[0].p == no[1].p) continue;

        s.reset();
        s.find_rail_path_naive(no[1].p, no[0].p, no[1].d);
        s.build_station(no[0].p);
        s.build_station(no[1].p);
        if (s.tm.income == 0) continue;

        int best_res = 0;
        int best_turn = 0;

        FOR(i, 2, no.S) {
            if (s.map[no[i].p] == 0) break;
            if (s.map[no[i].p] < 0) {
                bool used_naive = false;
                if (i < 50) {
                    int closest = -1;
                    int closest_dist = 1000;
                    REP(j, i) {
                        int d = dist(no[j].p, no[i].p);
                        if (d < closest_dist) {
                            closest_dist = d;
                            closest = j;
                        } else if (d == closest_dist) {
                            closest = -1;
                        }
                    }

                    if (closest != -1 && s.find_rail_path_naive(no[i].p, no[closest].p, no[i].d)) {
                        used_naive = true;
                    }

                }
                if (!used_naive) {
                    if (!s.bfs_path_quick(no[i].p, no[i].d)) break;
                }
            }
            s.build_station(no[i].p);

            int res = s.money + (T - s.turn) * s.tm.income;
            if (res > best_res) {
                best_res = res;
                best_turn = s.turn;
            } else {
                no.erase(no.begin() + min((int)no.S, i + 2), no.end());
                break;
            }
        }

        if (best_res >= cur_res - rng.next_double() * t) {
            if (!SILENT && best_res > cur_res) DB(hc_step, best_res, elapsed());
            acc[type]++;
            order = no;
            cur_res = best_res;
            if (M < 400) alltm.copy(s.tm);

            if (best_res > best_result) {
                best_result = best_res;
                s.save_solution(best_turn);
            }
        }
    }

    DATA(hc_step);    
    DB(acc);
    int hc_acc = 0; REP(i, acc.S) hc_acc += acc[i];
    DATA(hc_acc);
}

int main(int argc, char **argv) {
    // Read input
    int _; cin >> _ >> M >> K >> _;
    REP(i, M) {
        int r0, c0, r1, c1; cin >> r0 >> c0 >> r1 >> c1;
        trips[i] = {r0+1, c0+1, r1+1, c1+1};
    }

    DATA(M);
    DATA(K);

    // Preprocess station trips & trip income
    REP(i, M) {
        trip_income[i] = dist(trips[i].r0, trips[i].c0, trips[i].r1, trips[i].c1);

        FOR(r, max(1, trips[i].r0 - 2), min(N+1, trips[i].r0 + 3)) {
            FOR(c, max(1, trips[i].c0 - 2), min(N+1, trips[i].c0 + 3)) {
                if (dist(r, c, trips[i].r0, trips[i].c0) <= 2) {
                    station_trips0[p1d(r, c)].PB(i);
                    useful_station[p1d(r, c)] = true;
                }
            }
        }

        FOR(r, max(1, trips[i].r1 - 2), min(N+1, trips[i].r1 + 3)) {
            FOR(c, max(1, trips[i].c1 - 2), min(N+1, trips[i].c1 + 3)) {
                if (dist(r, c, trips[i].r1, trips[i].c1) <= 2) {
                    station_trips1[p1d(r, c)].PB(i);
                    useful_station[p1d(r, c)] = true;
                }
            }
        }
    }

    FOR(r, 1, N+1) FOR(c, 1, N+1) {
        int p = p1d(r, c);
        for (int id : station_trips0[p]) total_income[p] += trip_income[id];
        for (int id : station_trips1[p]) total_income[p] += trip_income[id];
    }

    // Solve

    rng.init(time(0));
    if (FIXED_SEED) rng.init(42);

    int runs = 0;
    LL gsum = 0;
    while (elapsed() < TIME_LIMIT_GREEDY) {
        gsum += do_greedy(runs == 0);
        runs++;
    }
    DATA(runs);

    int gavg = gsum / runs;
    DATA(gavg);

    if (FIXED_SEED) rng.init(42);

    do_hc(TIME_LIMIT_HC + TIME_LIMIT_GREEDY - elapsed());
    DB(best_result);

    optimize_build_order();

    print_solution();

    double time = elapsed();
    DATA(time);

    int stats = 0; REP(t, T) stats += solution[t].type == 0;
    DATA(stats);

	return 0;
}
