// Author: Psyho
// Twitter: https://twitter.com/fakepsyho

// TODO:




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
#define PIII        tuple<int, int, int>
#define VB          VC<byte>
#define VVB         VC<VB>
#define VI          VC<int>
#define VVI         VC<VI>
#define VVVI        VC<VVI>
#define VPII        VC<PII>
#define VVPII       VC<VPII>
#define VPIII       VC<PIII>
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
#define DATAX(a, b) {cerr << "[DATA] " << a << " = " << (b) << endl;}
 
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

template <class T, class Compare=less<T>> struct PQueue {
    VC<T> data;
    int size;
    Compare comp;

    PQueue(int n = 0) {reset(n);}

    void reset(int n) {data = VC<T>(n+1); size = 0;}
    void clear() {size = 0;}

    void push(const T& x) {
        data[++size] = x;
        for (int c = size, p = c >> 1; c > 1 && comp(data[c], data[p]); c = p, p = c >> 1) swap(data[c], data[p]);
    }

    void add(const T& x) {
        data[++size] = x;
    }

    void replace_top(const T& x) {
        int p = 1;
        for (int c = 2; c <= size; p = c, c = p << 1) {
            c += c < size && comp(data[c + 1], data[c]);
            if (!comp(data[c], x)) break;
            data[p] = data[c];
        }
        data[p] = x;
    }

    T pop() {
        T x = data[1];
        replace_top(data[size--]);
        return x;
    }

    INLINE T top() {return data[1];}
};

// SOLUTION

const double TIME_LIMIT = 2.9;

const int N = 600;
const int T = 600;

const int MAX_L_A = 1200;
const int MAX_L_B = 24;
const int MAX_M = 3 * N;


// input
int M, L_A, L_B;

int A[MAX_L_A];
int B[MAX_L_B];
VI cn[N];
int target_path[N];
PII coords[N];

// calc
int dist[N][N];

struct AUse {int len, a_pos, b_pos;};

int count_unique_edges(VI &path) {
    set<PII> edges;
    FOR(i, 1, path.S) edges.insert({min(path[i-1], path[i]), max(path[i-1], path[i])});
    return edges.S;
}

VVI split_path(VI &path) {
    VVI rv;
    VI current;
    int t_pos = 0;
    for (int p : path) {
        current.PB(p);
        if (p == target_path[t_pos]) {
            t_pos++;
            rv.PB(current);
            current.clear();
        }
    }
    assert(rv.S == T);
    return rv;
}

void print_path_uses(VI &path, VC<AUse> &a_uses) {
    assert(path.S == a_uses.S);
    int abc = 0; REP(i, a_uses.S) abc += a_uses[i].len != -1; DB(abc);
    VI v_pos(N, -1); REP(i, L_A) v_pos[A[i]] = i;
    VI b(L_B, -1);
    int fixed_cnt = 0;
    REP(i, path.S) {
        AUse &use = a_uses[i];
        if (use.len != -1) {
            REP(j, use.len) b[use.b_pos + j] = A[use.a_pos + j];
            cout << "s " << use.len << " " << use.a_pos << " " << use.b_pos << endl;
        }
        bool found = false; REP(j, L_B) found |= b[j] == path[i];
        if (!found) {
            if (b[L_B-1] != -1) b[L_B-1] = path[i];
            cout << "s " << 1 << " " << v_pos[path[i]] << " " << L_B-1 << endl;
            fixed_cnt++;
        }
        cout << "m " << path[i] << endl;
    }
    DB(fixed_cnt);
}

void assert_route(const VI &route) {
    assert(route.S >= L_B);
    REP(i, route.S - 1) {
        bool connected = false;
        for (int v : cn[route[i]]) if (v == route[i+1]) {connected = true; break;}
        assert(connected);
    }
}

VI last_a;
VC<AUse> optimize_usage_and_array(const VI &path, const VI &_a, double time_limit, const VC<AUse> &a_init = {}) {
    VI a = _a;

    double opt_start = elapsed();

    VVI b(path.S, VI(L_B, -1));
    VC<AUse> a_uses(path.S, {-1, -1, -1});

    VVI vnext2b(path.S, VI(N, L_B * 5));
    REP(i, path.S) FOR(j, i, min(int(path.S), i + L_B * 5)) vnext2b[i][path[j]] = min(vnext2b[i][path[j]], j - i);

    if (a_init.S) a_uses = a_init;

    VI path_cost(path.S, 0);

    VVI v_pos(N);
    REP(i, L_A) if (a[i] != -1) v_pos[a[i]].PB(i);

    int mx = 0; REP(i, L_A) if (a[i] != -1) mx = max(mx, a[i]); DB(mx);

    bitset<N> visited; visited.reset(); for (int v : path) visited.set(v, 1);
    VI vlist; REP(i, N) if (visited[i]) vlist.PB(i);

    int short_routes = 0; for (AUse &use : a_uses) if (use.len == 1) short_routes++;

    LL eval_partial_steps = 0;
    LL eval_partial_runs = 0;

    auto eval_full = [&]() -> int {
        int rv = 0;
        bitset<N> active; active.reset();
        REP(i, path.S) {
            if (i) REP(j, L_B) b[i][j] = b[i-1][j];
            auto [len, a_pos, b_pos] = a_uses[i];
            if (len != -1) {
                REP(j, len) if (b[i][b_pos + j] != -1) active.set(b[i][b_pos + j], 0);
                REP(j, len) b[i][b_pos + j] = a[a_pos + j];
                REP(j, L_B) if (b[i][j] != -1) active.set(b[i][j], 1);
            }
            bool found = active.test(path[i]);
            if (!found) {
                if (b[i][L_B-1] != -1) active.set(b[i][L_B-1], 0);
                b[i][L_B-1] = path[i];
                active.set(path[i], 1);
                REP(j, L_B) if (b[i][j] != -1) active.set(b[i][j], 1);
            }
            path_cost[i] = (!found) + (len != -1);
            rv += path_cost[i];
        }
        return rv;
    };

    auto eval_full_tmp = [&]() -> int {
        static VC<char> active(N);
        VI cur_b(L_B, -1);
        int rv = 0;
        REP(i, path.S) {
            auto [len, a_pos, b_pos] = a_uses[i];
            if (len != -1) {
                REP(j, len) if (cur_b[b_pos + j] != -1) active[cur_b[b_pos + j]]--;
                REP(j, len) cur_b[b_pos + j] = a[a_pos + j];
                REP(j, len) if (cur_b[b_pos + j] != -1) active[cur_b[b_pos + j]]++;
            }
            bool found = active[path[i]] > 0;
            if (!found) {
                if (cur_b[L_B-1] != -1) active[cur_b[L_B-1]]--;
                cur_b[L_B-1] = path[i];
                active[path[i]]++;
            }
            rv += (!found) + (len != -1);
        }
        REP(j, L_B) if (cur_b[j] != -1) active[cur_b[j]] = 0;
        return rv;
    };

    auto eval_partial = [&](int pos, bool update = false) -> int {
        static VC<char> active(N); 
        static VI cur_b(L_B);
        if (pos) {
            REP(j, L_B) cur_b[j] = b[pos-1][j];
        } else {
            REP(j, L_B) cur_b[j] = -1;
        }
        REP(j, L_B) if (cur_b[j] != -1) active[cur_b[j]]++;
        int rv = 0;
        FOR(i, pos, path.S) {
            auto [len, a_pos, b_pos] = a_uses[i];
            if (len != -1) {
                REP(j, len) if (cur_b[b_pos + j] != -1) active[cur_b[b_pos + j]]--;
                REP(j, len) cur_b[b_pos + j] = a[a_pos + j];
                REP(j, len) if (cur_b[b_pos + j] != -1) active[cur_b[b_pos + j]]++;
            }
            bool found = active[path[i]] > 0;
            if (!found) {
                if (cur_b[L_B-1] != -1) active[cur_b[L_B-1]]--;
                cur_b[L_B-1] = path[i];
                active[path[i]]++;
            }
            rv -= path_cost[i];
            rv += (!found) + (len != -1);
            if (!update && rv >= 3) break;
            bool eval_merged = true;
            REP(j, L_B) if (cur_b[j] != b[i][j]) {eval_merged = false; break;}
            if (update) {
                path_cost[i] = (!found) + (len != -1);
                REP(j, L_B) b[i][j] = cur_b[j];
            }
            if (eval_merged) {
                eval_partial_steps += i - pos + 1;
                eval_partial_runs++;
                break;
            }
        }
        REP(j, L_B) if (cur_b[j] != -1) active[cur_b[j]] = 0;
        // REP(j, N) assert(active[j] == 0);
        return rv;
    };

    int step = 0;
    VI acc(5);
    int false_acc = 0;
    int better = 0;
    double time_passed = 0;
    int bv = eval_full();
    int xv = bv;
    DB(bv);

    double temp0 = .10;
    double tempn = .001;
    double temp = temp0;
    while (true) {
        step++;
        if ((step & 31) == 0) {
            time_passed = (elapsed() - opt_start) / time_limit;
            temp = temp0 * pow(tempn / temp0, time_passed);
            if (time_passed > 1) break;
        }

        int type = rng.next(3);
        // if (rng.next_double() < 0.5) type = rng.next(2);

        int t;
        AUse old_t, new_t;
        int old_v, new_v;
        int rev0, rev1;
        if (type == 0) {
            t = rng.next(path.S);
            int pos = v_pos[path[t]][rng.next(v_pos[path[t]].S)];
            assert(a[pos] == path[t]);
            int len = rng.next(2, L_B + 1);
            int a_pos = rng.next(max(0, pos - len + 1), min(L_A - len + 1, pos + 1));
            while (vnext2b[t][a[a_pos]] > 2 * L_B) a_pos++, len--;
            while (vnext2b[t][a[a_pos + len - 1]] > 2 * L_B) len--;
            int b_pos = rng.next(L_B - len + 1);

            new_t = {len, a_pos, b_pos};
            old_t = a_uses[t];
        } else if (type == 1) {
            REP(_, 20) {
                t = rng.next(path.S);
                if (a_uses[t].len != -1) break;
            }
            if (a_uses[t].len == -1) continue;
            new_t = {-1, -1, -1};
            old_t = a_uses[t];
        } else if (type == 2) {
            REP(_, 20) {
                t = rng.next(path.S);
                if (a_uses[t].len != -1) break;
            }
            if (a_uses[t].len == -1) continue;
            old_t = a_uses[t];
            auto [old_len, old_a_pos, old_b_pos] = old_t;
            int len, a_pos, b_pos;
            int reps = 0;
            while (true) {
                reps++;
                if (reps > 50) {len = -1; break;}
                int subtype = rng.next(3);
                len = old_len;
                a_pos = old_a_pos;
                b_pos = old_b_pos;
                if (subtype == 0) {
                    b_pos = rng.next(L_B - len + 1);
                } else if (subtype == 1) {
                    a_pos = rng.next(max(0, a_pos - len + 1), min(L_A - len + 1, a_pos + 1));
                } else {
                    len = rng.next(1, min(L_A - a_pos + 1, L_B - b_pos + 1));
                }
                if (old_len == len && old_a_pos == a_pos && old_b_pos == b_pos) continue;
                bool ok = false;
                REP(j, len) ok |= a[a_pos + j] == path[t];
                if (ok) break;
            }
            if (len == -1) continue;
            while (vnext2b[t][a[a_pos]] > 2 * L_B) a_pos++, len--;
            while (vnext2b[t][a[a_pos + len - 1]] > 2 * L_B) len--;
            if (len == old_len && a_pos == old_a_pos && b_pos == old_b_pos) continue;
            new_t = {len, a_pos, b_pos};

        } 

        int av;
        a_uses[t] = new_t;
        av = bv + eval_partial(t);

        if (av <= bv || rng.next_double() < exp((bv - av) / temp)) {
            acc[type]++;
            eval_partial(t, true);

            bv = av;
            if (av < xv) {
                xv = av;
                DB(step, av, time_passed);
            }
        } else {
            a_uses[t] = old_t;
        }
    }
    last_a = a;

    DB(step, acc, false_acc, better, bv);

    DB(eval_partial_steps * 1.0 / eval_partial_runs, eval_partial_runs);

    DATAX("uastep", step);

    REP(i, path.S) if (path_cost[i] && a_uses[i].len == -1) a_uses[i] = {1, v_pos[path[i]][0], L_B-1};
    return a_uses;
}


VI find_path() {
    int edge_cnt[N][N];
    ZERO(edge_cnt);

    int v_cnt[N];

    VI bpath;

    VVPII part_edges(T);


    REP(pass, 2) {
        int current = 0;
        bpath.clear();
        VI last_visited(N, -100);
        REP(i, T) {
            for (auto &e : part_edges[i]) edge_cnt[e.X][e.Y]--, edge_cnt[e.Y][e.X]--, v_cnt[e.Y]--;

            int target = target_path[i];

            static priority_queue<PIII, VPIII, greater<PIII>> pq;
            static VI prev(N, -1);
            static VI bdist(N, 1 << 20);
            REP(i, N) prev[i] = -1;
            REP(i, N) bdist[i] = 1 << 20;

            pq.push({0, bpath.S, current});
            bdist[current] = 0;

            while (prev[target] == -1) {
                int d, plen, v;
                tie(d, plen, v) = pq.top();
                pq.pop();

                if (d > bdist[v]) continue;

                for (int u : cn[v]) {
                    // int nd = d + 1;
                    int nd = d + 100 + (v_cnt[u] == 0) * 200 + (v_cnt[u] == 1) * 50 + (v_cnt[u] == 2) * 25 + (v_cnt[u] == 3) * 10;
                    if (nd < bdist[u]) {
                        bdist[u] = nd;
                        prev[u] = v;
                        pq.push({nd, plen + 1, u});
                    }
                }
            }

            while (!pq.empty()) pq.pop();

            static VI path;
            path.clear();
            for (int v = target; v != current; v = prev[v]) path.PB(v);
            reverse(ALL(path));

            part_edges[i].clear();

            for (int p : path) {
                edge_cnt[current][p]++;
                edge_cnt[p][current]++;
                v_cnt[p]++;
                last_visited[p] = bpath.S;
                bpath.PB(p);
                part_edges[i].PB({current, p});
                current = p;
            }
        }
        int unique_edges = 0;
        REP(i, N) REP(j, i) unique_edges += edge_cnt[i][j] > 0;
        DB(unique_edges, bpath.S, elapsed());
    }

    return bpath;
}

VVI generate_routes(const VI& path = {}) {
    VI vs_cnt(N, 0);
    VVI routes;

    VI target_cnt(N);
    REP(i, T) target_cnt[target_path[i]]++;

    VI viable_v(N);
    VI tmp_path = path.S ? path : find_path();
    for (int p : tmp_path) viable_v[p] = 1;

    // VI viable_v(N, 0);
    // REP(t, T) {
    //     int v0 = t ? target_path[t-1] : 0;
    //     int v1 = target_path[t];
    //     REP(i, N) viable_v[i] |= dist[v0][v1] == dist[v0][i] + dist[i][v1];
    // }

    int viable_cnt = 0;
    REP(i, N) viable_cnt += viable_v[i];
    // DB(viable_cnt);

    int routes_len = 0;

    int step = 0;
    // DB(elapsed());

    const int TOTAL_STEPS = L_B < 18 ? 50000 : 200000;
    const double OFFSET1_PROB = L_B < 11 ? rng.next_double(50.0, 200.0) : L_B < 18 ? rng.next_double(10.0, 50.0) : rng.next_double(0.2, 2.0);
    // const double OFFSET2_PROB = 10.0;

    int start_step_len = L_B < 18 ? +5 : 0;
    int min_step_len = rng.next(4, 16); 
    int eval_power = 6;

    auto merge_paths = [&](int p1, int p2, int rev1, int rev2, int dist) {
        assert(p1 != p2);
        assert(dist <= 1);
        if (rev1) reverse(ALL(routes[p1]));
        if (rev2) reverse(ALL(routes[p2]));
        routes[p1].insert(routes[p1].end(), routes[p2].begin() + (1 - dist), routes[p2].end());
        routes.erase(routes.begin() + p2);
    };

    while (true) {
        step++;
        if (step > TOTAL_STEPS) break;
        int v0 = -1, v1 = -1;
        if (rng.next_double() < 0.1) {
            REP(_, 10) {
                int x;
                x = rng.next(N); if (viable_v[x] && (v0 == -1 || vs_cnt[x] < vs_cnt[v0])) v0 = x;
                x = rng.next(N); if (viable_v[x] && (v1 == -1 || vs_cnt[x] < vs_cnt[v1])) v1 = x;
            }
        } else {
            REP(_, 2) {
                int t = rng.next(T);
                int a = t ? target_path[t-1] : 0;
                int b = target_path[t];
                if (rng.next(2)) swap(a, b);
                if (v0 == -1 || dist[a][b] > dist[v0][v1]) v0 = a, v1 = b;
            }
        }
        if (v0 < 0 || v1 < 0 || v0 == v1 || !viable_v[v0] || !viable_v[v1] || dist[v0][v1] < L_B) continue;
        VI route; route.PB(v0);
        int p = v0;
        static VI route_vs(N);
        route_vs[v0] = 1;
        double max_dist_prob1 = 1.0 - (1.0 * step / TOTAL_STEPS + OFFSET1_PROB) / (1.0 + OFFSET1_PROB);
        double max_dist_prob2 = 0.01;
        while (p != v1) {
            int next = -1;
            int cnt = 0;
            for (int v : cn[p]) {
                if (!viable_v[v]) continue;
                if (route_vs[v]) continue;
                int max_dist = dist[p][v1] - 1;
                max_dist += rng.next_double() < max_dist_prob1;
                // max_dist += rng.next_double() < max_dist_prob2;
                if (dist[v][v1] > max_dist) continue;
                int priority = vs_cnt[v] == 0 ? 20 : vs_cnt[v] == 1 ? 5 : vs_cnt[v] == 2 ? 2 : 1;
                // int priority = vs_cnt[v] == 0 ? 5 : 1;
                // int priority = 1;
                cnt += priority;
                if (rng.next(cnt) < priority) next = v;
            }
            if (next == -1) break;
            p = next;
            route_vs[p]++;
            route.PB(p);
        }
        REP(pass, 2) {
            if (rng.next(2)) while (true) {
                bool found = false;
                int next = 0;
                for (int v : cn[route.back()]) if (viable_v[v] && !route_vs[v] && vs_cnt[v] == 0) {
                    found = true;
                    next = v;
                    break;
                }
                if (!found) break;
                route_vs[next]++;
                route.PB(next);
                if (rng.next(2)) break;
            }
            reverse(ALL(route));
            break;
        }
        for (int v : route) route_vs[v] = 0;
        if (route.S < max(1, int((start_step_len + L_B + step * min_step_len / TOTAL_STEPS)))) continue;
        // if (route.S < L_B) continue;
        // DB(route);

        for (int r : route) vs_cnt[r]++;
        routes.PB(route);
        routes_len += route.S;

        while (true) {
            int missing = 0; REP(i, N) missing += vs_cnt[i] == 0 && viable_v[i];
            if (routes_len <= L_A - missing) break;
            double bv = -1e30;
            int best = -1;
            REP(i, routes.S) {
                double av = 0;
                for (int v : routes[i]) av -= 1.0 / pow(vs_cnt[v] - 0.99, eval_power) * (.1 + target_cnt[v] * target_cnt[v] * target_cnt[v]) ; 
                av /= routes[i].S;
                if (av > bv) {
                    bv = av;
                    best = i;
                }
            }
            assert(best >= 0);
            for (int v : routes[best]) vs_cnt[v]--;
            routes_len -= routes[best].S;
            routes.erase(routes.begin() + best);
        }

    }
    // DB(elapsed());

    VI show_vs; REP(i, N) if (viable_v[i]) show_vs.PB(vs_cnt[i]);
    int verify_len = 0; for (auto &r : routes) verify_len += r.S;
    sort(ALL(show_vs));
    // DB(show_vs, routes.S, routes_len, verify_len);

    REP(i, N) if (viable_v[i] && vs_cnt[i] == 0) routes.PB({i});

    // try and merge paths
    while (true) {
        REP(i, routes.S) REP(j, i) {
            int v0a = routes[i][0];
            int v0b = routes[i].back();
            int v1a = routes[j][0];
            int v1b = routes[j].back();
            bool found = false;
            if      (dist[v0a][v1a] == 0) merge_paths(i, j, 1, 0, 0), found = true;
            else if (dist[v0a][v1b] == 0) merge_paths(i, j, 1, 1, 0), found = true;
            else if (dist[v0b][v1a] == 0) merge_paths(i, j, 0, 0, 0), found = true;
            else if (dist[v0a][v1b] == 0) merge_paths(i, j, 0, 1, 0), found = true;
            if (found) goto next2;
            if      (dist[v0a][v1a] == 1) merge_paths(i, j, 1, 0, 1), found = true;
            else if (dist[v0a][v1b] == 1) merge_paths(i, j, 1, 1, 1), found = true;
            else if (dist[v0b][v1a] == 1) merge_paths(i, j, 0, 0, 1), found = true;
            else if (dist[v0a][v1b] == 1) merge_paths(i, j, 0, 1, 1), found = true;
            if (found) goto next2;
        }
        break;
        next2: ;
    }

    return routes;
}

VVI optimize_routes(const VI &path, double time_limit) {
    int routes_len = 0;
    VPII routes;

    double temp0 = 2.0;
    double tempn = 0.01;
    double temp = temp0;
    int step = 0;

    auto sim = [&]() -> int {
        static VPII vroutes[N]; REP(i, N) vroutes[i].clear();
        REP(i, routes.S) FOR(j, routes[i].X, routes[i].Y) vroutes[path[j]].PB({i, j});

        int pos = 0;
        int rv = 0;
        while (pos < path.S) {
            int best = 0;
            for (auto [r, p] : vroutes[path[pos]]) {
                int dir_up = 0; while (p + dir_up <  routes[r].Y && path[pos + dir_up] == path[p + dir_up]) dir_up++;
                int dir_dw = 0; while (p - dir_dw >= routes[r].X && path[pos + dir_dw] == path[p - dir_dw]) dir_dw++;
                best = max(best, dir_up);
                best = max(best, dir_dw);
            }

            best = max(best, 1);
            best = min(best, L_B);

            pos += best;
            rv++;
        }
        return rv;
    };

    const int T2_RANGE = 2;
    int bv = sim();
    int xv = bv;
    DB(bv);
    time_limit *= 100;
    while (true) {
        step++;
        if ((step & 31) == 0) {
            double time_passed = elapsed() / time_limit;
            temp = temp0 * pow(tempn / temp0, time_passed);
            if (time_passed > 1) break;
        }

        int type = rng.next(3);

        int r, lo, hi;
        int old_lo, old_hi;

        if (type == 0) {
            int len = rng.next(1, 50);
            if (routes_len + len > L_A) continue;
            int start = rng.next(path.S - len + 1);
            lo = start;
            hi = start + len;
            routes.PB({lo, hi});
        } else if (type == 1) {
            if (routes.S == 0) continue;
            r = rng.next(routes.S);
            old_lo = routes[r].X;
            old_hi = routes[r].Y;
            swap(routes[r], routes.back());
            routes.pop_back();
        } else if (type == 2) {
            if (routes.S == 0) continue;
            r = rng.next(routes.S);
            old_lo = routes[r].X;
            old_hi = routes[r].Y;
            lo = rng.next(max(0, old_lo - T2_RANGE), min((int)path.S, old_lo + T2_RANGE + 1));
            hi = rng.next(max(0, old_hi - T2_RANGE), min((int)path.S, old_hi + T2_RANGE + 1));
            if (lo >= hi) continue;
            if (routes_len + hi - lo - (old_hi - old_lo) > L_A) continue;
            routes[r] = {lo, hi};
        }

        int av = sim();
        if (av <= bv || rng.next_double() < exp((bv - av) / temp)) {
            bv = av;
            if (av < xv) {
                DB(step, av, routes_len, temp);
                xv = av;
            }
            if (type == 0) routes_len += hi - lo;
            else if (type == 1) routes_len -= old_hi - old_lo;
            else if (type == 2) routes_len += hi - lo - (old_hi - old_lo);
        } else {
            if (type == 0) routes.pop_back();
            else if (type == 1) routes.PB({old_lo, old_hi});
            else if (type == 2) routes[r] = {old_lo, old_hi};
        }

    }

    DB(bv, step);

    VVI routes_rv(routes.S);
    REP(i, routes.S) routes_rv[i] = VI(path.begin() + routes[i].X, path.begin() + routes[i].Y);
    return routes_rv;
}


pair<VI, VC<AUse>> find_path_routes(const VVI &routes) {
    // complexity is O(path_len * N * L_A * L_B) or O(N * L_A * L_B); could be optimized to O(N * L_A)?

    // for (const VI &route : routes) assert_route(route);

    int n_routes = routes.S;
    VI route_start(n_routes); FOR(i, 1, n_routes) route_start[i] = route_start[i-1] + routes[i-1].S;
    VI route_vertices; for (const VI &route : routes) route_vertices.insert(route_vertices.end(), ALL(route));
    VI route_length(L_A); REP(r, n_routes) REP(i, max((int)routes[r].S - L_B + 1, 1)) route_length[route_start[r]+i] = min(L_B, (int)routes[r].S - i);
    VVI vertex_routes(N); 
    VC<bitset<MAX_L_A>> vertex_routes_set(N);
    REP(r, n_routes) REP(i, max((int)routes[r].S - L_B + 1, 1)) REP(j, route_length[route_start[r]+i]) {
        vertex_routes[routes[r][i+j]].PB(route_start[r]+i);
        vertex_routes_set[routes[r][i+j]].set(route_start[r]+i);
    }

    VVI parent(T+1, VI(L_A, -1));
    int cur_target = 0;

    int cur_dist = 1;
    bitset<MAX_L_A> start_set;
    for (int r : vertex_routes[0]) start_set.set(r);
    VI cur = vertex_routes[0];
    for (int v : cn[0]) for (int r : vertex_routes[v]) if (!start_set.test(r)) {
        start_set.set(r);
        cur.PB(r);
    }

    for (int r : cur) parent[0][r] = -2; // special starting marker

    while (true) {
        while (cur_target < T) {
            bool target_hit = false;
            int target = target_path[cur_target];
            for (int r : cur) if (vertex_routes_set[target].test(r)) {
                target_hit = true;
                break;
            }

            if (!target_hit) break;

            REP(i, cur.S) if (!vertex_routes_set[target].test(cur[i])) {
                cur[i] = cur.back();
                cur.pop_back();
                i--;
            }

            for (int r : cur) parent[cur_target+1][r] = r;
            cur_target++;
        }
        if (cur_target == T) break;
        cur_dist++;

        assert(cur.S > 0);

        static bitset<N> cur_set; 
        static VPII cur_v; cur_v.clear();
        for (int r : cur) REP(i, route_length[r]) {
            int v = route_vertices[r+i];
            if (cur_set.test(v)) continue;
            cur_set.set(v);
            cur_v.PB({v, r});
        }
        for (auto [v, prev] : cur_v) cur_set.set(v, 0);
        assert(cur_v.S > 0);

        static bitset<N> next_set; 
        static VPII next_v; next_v.clear();
        for (auto [v, prev] : cur_v) for (int u : cn[v]) if (!next_set.test(u)) {
            next_set.set(u);
            next_v.PB({u, prev});
        }
        for (auto [v, prev] : next_v) next_set.set(v, 0);
        assert(next_v.S > 0);

        static bitset<MAX_L_A> next_r_set;
        cur.clear();
        for (auto [v, prev] : next_v) {
            for (int r : vertex_routes[v]) {
                if (next_r_set.test(r)) continue;
                if (parent[cur_target][r] != -1) continue;
                parent[cur_target][r] = prev;
                cur.PB(r);
            }
        }
        for (int r : cur) next_r_set.set(r, 0);

        if (cur.S == 0) throw 1;
        assert(cur.S > 0);
        
    }

    DB(cur_dist);


    // reconstruct path

    VPII order;
    int t = T;
    int r = cur[0];
    while (parent[t][r] != -2) {
        order.PB({t, r});
        int prev = parent[t][r];
        if (prev == r) t--;
        r = prev;
    }
    order.PB({t, r});
    reverse(ALL(order));

    VI path;
    VC<pair<int, AUse>> a_uses;

    auto move_along_route = [&](int v0, int v1, int r) {
        int pos0 = -1; REP(i, route_length[r]) if (route_vertices[r+i] == v0) pos0 = i;
        int pos1 = -1; REP(i, route_length[r]) if (route_vertices[r+i] == v1) pos1 = i;
        assert(pos0 != -1 && pos1 != -1);
        while (pos0 < pos1) path.PB(route_vertices[r + ++pos0]);
        while (pos0 > pos1) path.PB(route_vertices[r + --pos0]);
    };

    auto move_route_to_route = [&](int v0, int r0, int r1) {
        int pos0 = -1; REP(i, route_length[r0]) if (route_vertices[r0+i] == v0) pos0 = i;
        assert(pos0 != -1);
        int v1 = -1, v2 = -1;
        REP(i, route_length[r0]) REP(j, route_length[r1]) if (dist[route_vertices[r0+i]][route_vertices[r1+j]] <= 1) {
            v1 = route_vertices[r0+i];
            v2 = route_vertices[r1+j];
            break;
        }
        int pos1 = -1; REP(i, route_length[r0]) if (route_vertices[r0+i] == v1) pos1 = i;
        assert(pos1 != -1);
        while (pos0 < pos1) path.PB(route_vertices[r0 + ++pos0]);
        while (pos0 > pos1) path.PB(route_vertices[r0 + --pos0]);
        a_uses.PB({(int)path.S, {route_length[r1], r1, 0}});
        if (path.back() != v2) path.PB(v2);
    };

    t = 0;
    r = order[0].Y;

    // move from 0 to r
    int v = 0;
    a_uses.PB({0, {route_length[r], r, 0}});
    bool first_move_found = false;
    REP(i, route_length[r]) if (dist[v][route_vertices[r+i]] <= 1) {
        if (dist[v][route_vertices[r+i]] == 1) {
            v = route_vertices[r+i];
            path.PB(v);
        }
        first_move_found = true;
        break;
    }
    assert(first_move_found);

    FOR(i, 1, order.S) {
        auto [nt, nr] = order[i];

        if (nt > t) {
            assert(r == nr);
            move_along_route(v, target_path[t], r);
            v = path.back();
            t = nt;
            continue;
        }

        assert(nt == t);
        assert(nr != r);
        move_route_to_route(v, r, nr);
        v = path.back();
        r = nr;
    }

    VC<AUse> a_uses_v; for (auto [i, use] : a_uses) {
        while (a_uses_v.S < i) a_uses_v.PB({-1, -1, -1});
        a_uses_v.PB(use);
    }
    while (a_uses_v.S < path.S) a_uses_v.PB({-1, -1, -1});

    return MP(path, a_uses_v);
}

int main(int argc, char **argv) {
    cin.tie(0);
    ios::sync_with_stdio(0);

    int _;
    cin >> _ >> M >> _ >> L_A >> L_B;

    // read data
    REP(i, M) {
        int u, v;
        cin >> u >> v;
        cn[u].PB(v);
        cn[v].PB(u);
    }

    REP(i, N) cin >> target_path[i];

    REP(i, N) cin >> coords[i].X >> coords[i].Y;

    DATA(L_A);
    DATA(L_B);
    DATA(M);

    // precalc stuff
    // dist[N][N]
    REP(v0, N) {
        static int bfs[MAX_M];
        int bfs_s = 0, bfs_e = 0;
        REP(i, N) dist[v0][i] = -1;
        dist[v0][v0] = 0;
        bfs[bfs_e++] = v0;
        while (bfs_s < bfs_e) {
            int v = bfs[bfs_s++];
            for (int u : cn[v]) if (dist[v0][u] == -1) {
                dist[v0][u] = dist[v0][v] + 1;
                bfs[bfs_e++] = u;
            }
        }
    }

    // try with routes
    double route_start = elapsed();
    int route_step = 0;
    int best_result = 1e9;
    VVI broutes;
    VI brpath;
    VC<AUse> bruses;
    while (route_step == 0 || elapsed() + (elapsed() - route_start) / route_step < TIME_LIMIT * .85) {
        try {
            VVI routes = generate_routes(elapsed() > TIME_LIMIT * .2 && brpath.S ? brpath : VI());
            auto [rpath, ruses] = find_path_routes(routes);
            int result = 0; for (AUse& use : ruses) result += use.len != -1;
            route_step++;
            if (result < best_result) {
                best_result = result;
                brpath = rpath;
                bruses = ruses;
                broutes = routes;
            }
        } catch (...) { }
    }
    REP(i, L_A) A[i] = 0;
    int pp = 0; for (VI &v : broutes) for (int x : v) A[pp++] = x;
    bruses = optimize_usage_and_array(brpath, VI(A, A+L_A), TIME_LIMIT - elapsed(), bruses);
    REP(i, L_A) A[i] = last_a[i];
    REP(i, L_A) cout << A[i] << " \n"[i == L_A - 1];
    print_path_uses(brpath, bruses);
    DATAX("rstep", route_step);
    DATAX("time", elapsed());

	return 0;
}
