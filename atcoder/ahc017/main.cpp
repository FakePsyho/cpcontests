// Author: Psyho
// Twitter: https://twitter.com/fakepsyho

// TEMPLATE

#include <vector>
#include <set>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <map>
#include <algorithm>
#include <cstring>
#include <cmath>
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
#define VL          VC<LL>
#define VVL         VC<VL>
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

INLINE double elapsed() {return get_time() - start_time;}

// SOLUTION

#ifdef VM
const double TIME_SCALE = 0.6;
#else
const double TIME_SCALE = 1.0;
#endif

const double SCALE = 1.0;

const double TIME_LIMIT = 5.9 * TIME_SCALE * SCALE;

const int MAX_DIST = 1e6;

const int MAX_N = 1000;
const int MAX_M = 3000;
const int MAX_D = 30;

const int MAX_SUBSET = MAX_N;
const int DIST_DIV = 1000;

const bool POST_HC = true;
const bool ADJUST_BEAM = true;
const bool MAX_BEAM = false;
const bool MAX_FRUST = false;
const bool STATS = false;

int N, M, D, K;
VVPII ove;
VPII vp;
int fd[MAX_N][MAX_N];
int ed[MAX_M][MAX_M];

short ddist[MAX_D][MAX_SUBSET][MAX_N];
int ddist_sum[MAX_D][MAX_SUBSET];

double avg_dist;

VPII edges;

VI subset;


INLINE int dist2(int a, int b) {return (int)((vp[a].X - vp[b].X) * (vp[a].X - vp[b].X) + (vp[a].Y - vp[b].Y) * (vp[a].Y - vp[b].Y));}
INLINE int dist2e(int a, int b) {return min(min(dist2(edges[a].X, edges[b].X), dist2(edges[a].X, edges[b].Y)), min(dist2(edges[a].Y, edges[b].X), dist2(edges[a].Y, edges[b].Y)));}

void add_edges(VVPII &ve, VI &edge_ids) {
    for (int eid : edge_ids) {
        PII &e = edges[eid];
        ve[e.X].PB(MP(e.Y, fd[e.X][e.Y]));
        ve[e.Y].PB(MP(e.X, fd[e.X][e.Y]));
    }
}

void remove_edges(VVPII &ve, VI &edge_ids) {
    for (int eid : edge_ids) {
        PII &e = edges[eid];
        REP(i, ve[e.X].S) if (ve[e.X][i].X == e.Y) {
            ve[e.X][i] = ve[e.X].back();
            ve[e.X].pop_back();
            break;
        }
        REP(i, ve[e.Y].S) if (ve[e.Y][i].X == e.X) {
            ve[e.Y][i] = ve[e.Y].back();
            ve[e.Y].pop_back();
            break;
        }
    }
}

int call_count = 0;

int vs[MAX_N];
short vl[MAX_N*4][16];
short n_vl[MAX_N*4];
LL calc_all_distances(VVPII &ve, VI subset_ids, int day=-1) {
    call_count += subset_ids.S;

    LL rv = 0;
    for (int subset_id : subset_ids) {
        int v0 = subset[subset_id];
        memset(vs, 0x3F, sizeof(vs[0])*N);
        vl[0][0] = v0;
        n_vl[0] = 1;
        vs[v0] = 0;
        int cur_sum = 0;
        int last = 0;
        int d = 0;
        int conn = 0;
        while (d <= last) {
            REP(i, n_vl[d]) {
                int v = vl[d][i];
                if (d > vs[v]) continue;
                if (day >= 0) ddist[day][subset_id][v] = d;
                cur_sum += d;
                conn++;
                for (PII &p : ve[v]) {
                    if (d + p.Y < vs[p.X]) {
                        vs[p.X] = d + p.Y;
                        vl[d+p.Y][n_vl[d+p.Y]++] = p.X;
                        last = max(last, d+p.Y);
                    }
                }
            }
            d++;
        }
        memset(n_vl, 0, sizeof(n_vl[0])*(last+1));
        cur_sum +=(N - conn) * MAX_DIST;
        ddist_sum[day][subset_id] = cur_sum;
        rv += cur_sum;
        if (conn < N) return (LL)1e14;
    }
    return rv;
}

int skipped = 0;
int max_vl = 0;
int max_dist = 0;


template <bool UPDATE> LL calc_all_distances_rem_edge(VVPII &ve, VI subset_ids, int day, int rem_edge) {
    VI e = {rem_edge};
    LL rv = 0;

    remove_edges(ve, e);
    for (int subset_id : subset_ids) {

        int v0 = subset[subset_id];

        if (abs(ddist[day][subset_id][edges[rem_edge].X] - ddist[day][subset_id][edges[rem_edge].Y]) < fd[edges[rem_edge].X][edges[rem_edge].Y]) {
            skipped++;
            rv += ddist_sum[day][subset_id];
            continue;
        }

        int startv = ddist[day][subset_id][edges[rem_edge].X] > ddist[day][subset_id][edges[rem_edge].Y] ? edges[rem_edge].X : edges[rem_edge].Y;

        int cur_sum = ddist_sum[day][subset_id];

        static unsigned short q[MAX_N];
        static unsigned short bfs_vs[MAX_N];
        static unsigned short bfs_vs_id = 0;

        if (bfs_vs_id == 0) ZERO(bfs_vs);
        bfs_vs_id++;

        int qen = 1;
        q[0] = startv;
        bfs_vs[startv] = bfs_vs_id;

        REP(i, qen) {
            int v = q[i];
            cur_sum -= ddist[day][subset_id][v];
            for (PII &p : ve[v]) {
                if (bfs_vs[p.X] != bfs_vs_id && ddist[day][subset_id][v] + p.Y == ddist[day][subset_id][p.X]) {
                    bfs_vs[p.X] = bfs_vs_id;
                    q[qen++] = p.X;
                }
            }
        }

        int d = MAX_DIST;
        int last = 0;
        REP(i, qen) {
            int v = q[i];
            int bd = MAX_DIST;
            for (PII &p : ve[v]) if (bfs_vs[p.X] != bfs_vs_id) bd = min(ddist[day][subset_id][p.X] + p.Y, bd);
            vs[v] = bd;
            if (bd < MAX_DIST) {
                vl[bd][n_vl[bd]++] = v;
                d = min(bd, d);
                last = max(bd, last);
            }
        }

        if (last < d) {
            cur_sum += qen * MAX_DIST;
            rv += cur_sum;
            continue;
        }

        int start_d = d;
        while (d <= last) {
            REP(i, n_vl[d]) {
                int v = vl[d][i];
                if (d > vs[v]) continue;
                if (UPDATE) {
                    ddist_sum[day][subset_id] += d - ddist[day][subset_id][v];
                    ddist[day][subset_id][v] = d;
                }
                cur_sum += d;
                for (PII &p : ve[v]) {
                    if (bfs_vs[p.X] == bfs_vs_id && d + p.Y < vs[p.X]) {
                        vs[p.X] = d + p.Y;
                        vl[d+p.Y][n_vl[d+p.Y]++] = p.X;
                        last = max(last, d+p.Y);
                    }
                }
            }
            d++;
        }
        memset(&n_vl[start_d], 0, sizeof(int)*(last-start_d+1));
        rv += cur_sum;
    }
    add_edges(ve, e);

    return rv;
}

template <bool UPDATE> LL calc_all_distances_add_edge(VVPII &ve, VI subset_ids, int day, int rem_edge) {
    VI e = {rem_edge};
    LL rv = 0;

    for (int subset_id : subset_ids) {
        int v0 = subset[subset_id];

        if (abs(ddist[day][subset_id][edges[rem_edge].X] - ddist[day][subset_id][edges[rem_edge].Y]) <= fd[edges[rem_edge].X][edges[rem_edge].Y]) {
            skipped++;
            rv += ddist_sum[day][subset_id];
            continue;
        }

        int startv = ddist[day][subset_id][edges[rem_edge].X] > ddist[day][subset_id][edges[rem_edge].Y] ? edges[rem_edge].X : edges[rem_edge].Y;
        int startd = min(ddist[day][subset_id][edges[rem_edge].X], ddist[day][subset_id][edges[rem_edge].Y]) + fd[edges[rem_edge].X][edges[rem_edge].Y];

        int cur_sum = ddist_sum[day][subset_id];

        int d = startd;
        int last = startd;
        vl[startd][n_vl[startd]++] = startv;


        vs[startv] = startd;

        memset(vs, 0x3F, sizeof(vs[0])*N);

        int start_d = d;
        while (d <= last) {
            REP(i, n_vl[d]) {
                int v = vl[d][i];
                if (d > vs[v]) continue;
                cur_sum += d - ddist[day][subset_id][v];
                if (UPDATE) {
                    ddist_sum[day][subset_id] += d - ddist[day][subset_id][v];
                    ddist[day][subset_id][v] = d;
                }
                for (PII &p : ve[v]) {
                    if (d + p.Y < vs[p.X] && d + p.Y < ddist[day][subset_id][p.X]) {
                        vs[p.X] = d + p.Y;
                        vl[d+p.Y][n_vl[d+p.Y]++] = p.X;
                        last = max(last, d+p.Y);
                    }
                }
            }
            d++;
        }

        memset(&n_vl[start_d], 0, sizeof(int)*(last-start_d+1));
        rv += cur_sum;
    }

    return rv;
}

VI find_subset(int no, int steps) {
    const double EDGE_DIST_W = 1.0 / 0.005;

    VI v(no);
    VI used(N);

    int p = 0;
    REP(i, N) {
        if (used[i]) continue;
        v[p++] = i;
        used[i] = 1;
        if (p == no) break;
    }

    VI closest(no);
    double bv = 0; 
    REP(i, no) {
        int d = 1e9;
        REP(j, no) if (i != j) d = min(d, dist2(v[i], v[j]));
        closest[i] = d;
        bv += 1.0 / d;
    }
    DB(bv);
    // DB(closest);

    int acc = 0;
    int last_acc = 0;
    while (steps--) {
        int a = rng.next(no);
        int old = v[a];
        int x = rng.next(N);
        if (used[x]) continue;

        v[a] = x;

        double diff = 0; 
        int d = 1e9;
        REP(j, no) if (j != a) d = min(d, dist2(x, v[j]));
        diff += 1.0 / d;
        diff -= 1.0 / closest[a];
        REP(j, no) if (j != a && (closest[j] == dist2(old, v[j]) || closest[j] > dist2(x, v[j]))) {
            diff -= 1.0 / closest[j];
            if (dist2(x, v[j]) < closest[j]) {
                diff += 1.0 / dist2(x, v[j]);
            } else {
                int d = 1e9;
                REP(k, no) if (j != k) d = min(d, dist2(v[j], v[k]));
                diff += 1.0 / d;
            }
        }
        if (diff <= 0) {
            bv += diff;
            used[old] = 0;
            used[x] = 1;
            REP(i, no) {
                int d = 1e9;
                REP(j, no) if (i != j) d = min(d, dist2(v[i], v[j]));
                closest[i] = d;
            }
            acc++;
            last_acc = steps;
        } else {
            v[a] = old;
        }
    }
    DB(acc, last_acc, bv);

    return v;
}

int calc_beam_size(VPII &beam) {
    int calls_per_edge = 0;
    REP(i, beam.S) calls_per_edge += (i ? (beam[i].X - beam[i-1].X) : beam[i].X) * (i ? beam[i-1].Y : D);
    return beam.back().X + M * beam.back().X + M * calls_per_edge;
}

int main(int argc, char **argv) {
    // Parse Input & Precompute some stuff
    ios_base::sync_with_stdio(false);

	cin >> N >> M >> D >> K;
    ove = VVPII(N);
    vp = VPII(N);

    REP(i, M) {
        int a, b, c;
        cin >> a >> b >> c;
        a--; b--;
        c = (c + DIST_DIV / 2) / DIST_DIV;
        edges.PB(MP(a, b));
        ove[a].PB(MP(b, c));
        ove[b].PB(MP(a, c));
        fd[a][b] = fd[b][a] = c;
    }
    REP(i, N) {
        int a, b;
        cin >> a >> b;
        vp[i] = MP(a, b);
    }

    cerr << "[DATA] N = " << N << endl;
    cerr << "[DATA] M = " << M << endl;
    cerr << "[DATA] D = " << D << endl;
    cerr << "[DATA] K = " << K << endl;
    cerr << "[DATA] Kratio = " << 1.0 * K * D / M << endl;
    cerr << "[DATA] Mratio = " << 1.0 * M / N << endl;
    cerr << "[DATA] calls = " << (LL)D*M*N << endl;

    REP(i, M) REP(j, M) ed[i][j] = dist2e(i, j);

    // Adjust beam width based on expected number of calls to calc_all_distances
    int calls_total = 4000000000LL * sqrt(1.0 * M / N) / N * SCALE;
    VI beam_days = {4, 4, 3, 2, 1};
    VI beam_mul = {1,4,7,10,50};

    VPII cur_beam; REP(i, beam_days.S) cur_beam.PB(MP(beam_mul[i], beam_days[i]));
    DB(cur_beam);
    VVPII all_beams;
    all_beams.PB(cur_beam);

    FOR(i, 1, MAX_SUBSET+1) {
        REP(lv, cur_beam.S) {
            while (cur_beam[lv].X < beam_mul[lv] * i && (lv + 1 == cur_beam.S || cur_beam[lv].X < cur_beam[lv+1].X) && cur_beam[lv].X < min(MAX_SUBSET, N)) {
                cur_beam[lv].X++;
                all_beams.PB(cur_beam);
            }
        }
    }

    DB(all_beams.S);
    VPII beam;
    int beam_pos = -1;
    REP(i, all_beams.S) {
        int calls_required = calc_beam_size(all_beams[i]);
        if (calls_required <= calls_total) {
            beam = all_beams[i];
            beam_pos = i;
        } else {
            break;
        }
    }

    if (MAX_BEAM) beam = all_beams.back(), beam_pos = all_beams.S - 1;
    
    DB(beam);
    DB(elapsed());

    // Find & Reorder Subset
    if (beam[0].X < N / 4) {
        subset = find_subset(beam[0].X, 10000);
        VI used(N); REP(i, subset.S) used[subset[i]] = 1;
        VI non_used; REP(i, N) if (!used[i]) non_used.PB(i);
        REP(i, beam.back().X - beam[0].X) subset.PB(non_used[i]);
    } else {
        subset = VI(beam.back().X); REP(i, subset.S) subset[i] = i;
    }
    DB(subset.S);
    VI ids(subset.S); REP(i, subset.S) ids[i] = i;

    DB(elapsed());
    
    // Sort Edges by Frustration
    double t_frust = elapsed();
    VC<pair<LL, int>> dist;
    VVPII ve = ove;

    int SORT_EDGES_MAX_SUBSETS = M;
    VI frustration_subset_ids(ids.begin(), ids.begin() + min(SORT_EDGES_MAX_SUBSETS, (int)ids.S));
    VI subset_backup = subset;
    if (MAX_FRUST) {
        frustration_subset_ids = VI(N); REP(i, N) frustration_subset_ids[i] = i;
        subset = frustration_subset_ids;
    }
    calc_all_distances(ve, frustration_subset_ids, 0);
    REP(i, M) dist.PB(MP(calc_all_distances_rem_edge<false>(ve, frustration_subset_ids, 0, i), i));
    sort(dist.rbegin(), dist.rend());
    subset = subset_backup;
    cerr << "[DATA] t_frust = " << elapsed() - t_frust << endl;
    DB(elapsed());

    // Calculate initial ddist & days_dist
    VI bsol(M, -1);
    VVI days_edges(D);
    VVL days_dist(D, VL(subset.S, 0));
    REP(i, subset.S) {
        calc_all_distances(ove, {i}, 0);
        days_dist[0][i] = ddist_sum[0][i] + (i ? days_dist[0][i-1] : 0);
        FOR(d, 1, D) {
            days_dist[d][i] = days_dist[0][i];
            ddist_sum[d][i] = ddist_sum[0][i];
            memcpy(ddist[d][i], ddist[0][i], sizeof(ddist[0][0][0]) * N);
        }
    }

    double greedy_start = elapsed();
    VD greedy_timestamps(M+1);
    VI greedy_calls(M+1);

    greedy_timestamps[0] = elapsed();
    DB(elapsed());

    // Greedy
    int panic_mode = 0;
    VI days_left; 
    VVL days_values(D, VL(beam.S));
    REP(i, M) {
        int v = dist[i].Y;

        if (i >= 25 && ADJUST_BEAM && !MAX_BEAM) {
            int calls = greedy_calls[i] - greedy_calls[i-25];
            double time = greedy_timestamps[i] - greedy_timestamps[i-25];

            double time_left = TIME_LIMIT - greedy_timestamps[i];
            LL min_calls_left = (LL)(time_left * calls / time / (M - i));
            LL max_calls_left = (LL)(time_left * max(1.0, time_left * 16 / TIME_LIMIT) * calls / time / (M - i));

            REP(panic_level, 15) {
                if (panic_mode == panic_level && greedy_timestamps[i] > TIME_LIMIT * (1 - 0.2 * pow(0.7, panic_level) * (M - i) / M)) {
                    panic_mode = panic_level + 1;
                    int new_subsets = max(1, (int)subset.S * 2 / 3);
                    while (subset.S > new_subsets) {
                        subset.pop_back();
                        ids.pop_back();
                    }
                    for (int j = beam.S - 1; j >= 0; j--)
                        beam[j].X = min(beam[j].X, (int)subset.S);
                }
            }

            while (calc_beam_size(beam) + (M - i) * subset.S < min_calls_left && beam_pos+1 < all_beams.S) {
                beam_pos++;
                beam = all_beams[beam_pos];
                for (int j = beam.S - 1; j >= 0; j--)
                    beam[j].X = min(beam[j].X, (int)subset.S);
            }

            while (calc_beam_size(beam) + (M - i) * subset.S > max_calls_left && beam_pos > 0) {
                beam_pos--;
                beam = all_beams[beam_pos];
                for (int j = beam.S - 1; j >= 0; j--)
                    beam[j].X = min(beam[j].X, (int)subset.S);
            }
        }

        
        bool empty_day = false;
        days_left.clear();
        REP(j, D) if (days_edges[j].S < K) {
            if (days_edges[j].S == 0) {
                if (empty_day) continue;
                empty_day = true;
            }
            days_left.PB(j);
        }

        double mw = pow(1.0 * (M - 1 - i) / (M - 1), 3.5);
        REP(lv, beam.S) {
            VI subset_ids(ids.begin() + (lv ? beam[lv-1].X : 0), ids.begin() + beam[lv].X);
            VC<pair<double,int>> beam_scores(days_left.S);
            int no = 0;
            for (int day : days_left) {
                remove_edges(ve, days_edges[day]);
                days_values[day][lv] = calc_all_distances_rem_edge<false>(ve, subset_ids, day, v) + (lv ? days_values[day][lv-1] : 0) ;
                beam_scores[no++] = MP((days_values[day][lv] - days_dist[day][beam[lv].X-1]) * (1.0 - 0.30 * mw * pow(1.0 * (K - days_edges[day].S) / K, 12.0) * (days_edges[day].S == 0 ? 1.10 : 1.0)), day);
                add_edges(ve, days_edges[day]);
            }
            sort(ALL(beam_scores));
            days_left.clear();
            REP(j, min(beam[lv].Y, (int)beam_scores.S)) days_left.PB(beam_scores[j].Y);
            if (beam[lv].X == beam.back().X) break;
        }
        int bday = days_left[0];
        remove_edges(ve, days_edges[bday]);
        calc_all_distances_rem_edge<true>(ve, ids, bday, v);
        add_edges(ve, days_edges[bday]);

        days_edges[bday].PB(v);
        REP(j, subset.S) days_dist[bday][j] = ddist_sum[bday][j] + (j ? days_dist[bday][j-1] : 0);
        bsol[v] = bday;

        greedy_calls[i+1] = calc_beam_size(beam) + greedy_calls[i];
        greedy_timestamps[i+1] = elapsed();
    }
    cerr << "[DATA] t_greedy = " << elapsed() << endl;
    DB(call_count);
    DB(skipped);

    DB(beam, subset.S);


    int HC_MAX_SUBSETS = M;
    VI hc_subset_ids(ids.begin(), ids.begin() + min(HC_MAX_SUBSETS, (int)ids.S));

    VL hc_dist(D);

    REP(i, D) hc_dist[i] = days_dist[i][hc_subset_ids.S-1];
    LL xv = 0; for (LL d : hc_dist) xv += d;
    DB(xv);

    int steps = 0;
    int acc = 0;
    while (POST_HC && elapsed() < TIME_LIMIT) {
        int p = rng.next(M);
        int d1 = bsol[p];
        int d2 = rng.next(D);
        if (d1 == d2 || days_edges[d2].S == K) continue;

        remove_edges(ve, days_edges[d1]);
        LL v1 = calc_all_distances_add_edge<false>(ve, hc_subset_ids, d1, p);
        add_edges(ve, days_edges[d1]);

        remove_edges(ve, days_edges[d2]);
        LL v2 = calc_all_distances_rem_edge<false>(ve, hc_subset_ids, d2, p);
        add_edges(ve, days_edges[d2]);
        
        if (v1 + v2 <= hc_dist[d1] + hc_dist[d2]) {
            xv += v1 + v2 - hc_dist[d1] - hc_dist[d2];
            hc_dist[d1] = v1;
            hc_dist[d2] = v2;
            acc++;

            remove_edges(ve, days_edges[d1]);
            calc_all_distances_add_edge<true>(ve, hc_subset_ids, d1, p);
            add_edges(ve, days_edges[d1]);

            remove_edges(ve, days_edges[d2]);
            calc_all_distances_rem_edge<true>(ve, hc_subset_ids, d2, p);
            add_edges(ve, days_edges[d2]);

            days_edges[d1].erase(find(days_edges[d1].begin(), days_edges[d1].end(), p));
            days_edges[d2].PB(p);
            bsol[p] = d2;
        }
        steps++;
    }

    DB(steps, acc);
    DB(xv);

    if (STATS) cerr << "[DATA] max_vl = " << max_vl << endl;
    if (STATS) cerr << "[DATA] max_dist = " << max_dist << endl;
    cerr << "[DATA] time = " << elapsed() << endl;

    for (int x: bsol) cout << " " << x+1;
    cout << endl;


	return 0;
}
