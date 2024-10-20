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
#define VVVB        VC<VC<VB>>
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
#define DATAX(x, name) {cerr << "[DATA] " << name << " = " << (x) << endl;}
 
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

    INLINE T peek() {return data[1];}
};

template <class T> struct DEPQueue {
    VC<T> data;
    int size;

    DEPQueue(int n = 0) {reset(n);}

    void reset(int n) {data = VC<T>(n+1); size = 0;}
    void clear() {size = 0;}

    void push(const T& x) {data[++size] = x; push_up(size);}

    T pop_min() {
        T rv = data[1];
        data[1] = data[size--];
        push_down_min(1);
        return rv;
    }

    T pop_max() {
        if (size == 1) {
            size = 0;
            return data[1];
        }
        int pos = size == 2 || data[2] > data[3] ? 2 : 3;
        T rv = data[pos];
        data[pos] = data[size--];
        push_down_max(pos);
        return rv;
    }

    T replace_min(const T& x) {
        T rv = data[1];
        data[1] = x;
        push_down_min(1);
        return rv;
    }

    T replace_max(const T& x) {
        int pos = size <= 2 ? size : 2 + (data[3] > data[2]);
        T rv = data[pos];
        data[pos] = x;
        push_down_max(pos);
        return rv;
    }

    T top_min() {return data[1];}
    T top_max() {return data[2 + (size > 2 && data[3] > data[2])];}

    bool is_min_level(int i) {return __builtin_clz(i)&1;}

    void push_up(int i) {
        if (i == 1) return;
        int parent = i >> 1;
        if (is_min_level(i)) {
            if (data[i] > data[parent]) {
                swap(data[i], data[parent]);
                push_up_max(parent);
            } else {
                push_up_min(i);
            }
        } else {
            if (data[i] < data[parent]) {
                swap(data[i], data[parent]);
                push_up_min(parent);
            } else {
                push_up_max(i);
            }
        }
    }

    void push_up_min(int i) {
        while (i > 3) {
            int grandparent = i >> 2;
            if (data[i] >= data[grandparent]) break;
            swap(data[i], data[grandparent]);
            i = grandparent;
        }
    }

    void push_up_max(int i) {
        while (i > 3) {
            int grandparent = i >> 2;
            if (data[i] <= data[grandparent]) break;
            swap(data[i], data[grandparent]);
            i = grandparent;
        }
    }

    void push_down_min(int i) {
        while (2 * i <= size) {
            int m = min_descendant(i);
            if (data[m] >= data[i]) break;
            swap(data[m], data[i]);
            if (m < 4 * i) break;
            if (data[m] > data[m >> 1]) swap(data[m], data[m >> 1]);
            i = m;
        }
    }

    void push_down_max(int i) {
        while (2 * i <= size) {
            int m = max_descendant(i);
            if (data[m] <= data[i]) break;
            swap(data[m], data[i]);
            if (m < 4 * i) break;
            if (data[m] < data[m >> 1]) swap(data[m], data[m >> 1]);
            i = m;
        }
    }

    int min_descendant(int i) {
        if (size < 4 * i) {
            if (size < 2 * i) return i;
            return 2*i+(size > 2*i && data[2*i+1] < data[2*i]);
        }
        if (size < 4*i+4) {
            int rv = 2*i+1;
            FOR(j, 4*i, size+1) if (data[j] < data[rv]) rv = j;
            return rv;
        } else {
            int rv = 4*i;
            if (data[4*i+1] < data[rv]) rv = 4*i+1;
            if (data[4*i+2] < data[rv]) rv = 4*i+2;
            if (data[4*i+3] < data[rv]) rv = 4*i+3;
            return rv;
        }
    }

    int max_descendant(int i) {
        if (size < 4 * i) {
            if (size < 2 * i) return i;
            return 2*i+(size > 2*i && data[2*i+1] > data[2*i]);
        }
        if (size < 4*i+4) {
            int rv = 2*i+1;
            FOR(j, 4*i, size+1) if (data[j] > data[rv]) rv = j;
            return rv;
        } else {
            int rv = 4*i;
            if (data[4*i+1] > data[rv]) rv = 4*i+1;
            if (data[4*i+2] > data[rv]) rv = 4*i+2;
            if (data[4*i+3] > data[rv]) rv = 4*i+3;
            return rv;
        }
    }
};

template<typename Move> struct MoveHistory {
    int previous_id;
    Move move;
};

template<typename ValueT, typename Move> struct Candidate {
    int state_id;
    ValueT value;
    Move move;
};

template<typename ValueT, typename Move> struct CandidatesGroup {
    int size;
    int count;
    VC<Candidate<ValueT, Move>> candidates;
    PQueue<pair<ValueT, int>> pqueue;

    CandidatesGroup(int size = 0) : size(size) {
        candidates.resize(size);
        pqueue.reset(size);
        count = 0;
    }

    void clear() {
        pqueue.clear();
        count = 0;
    }

    ValueT min() {
        return pqueue.top().X;
    }

    bool add(const Candidate<ValueT, Move>& cand) {
        if (count < size) {
            candidates[count] = cand;
            pqueue.push(MP(cand.value, count));
            count++;
            return true;
        } else {
            if (cand.value > pqueue.top().X) {
                int pos = pqueue.top().Y;
                candidates[pos] = cand;
                pqueue.replace_top(MP(cand.value, pos));
                return true;
            } else {
                return false;
            }
        }
    }

    bool can_add(const ValueT value) {
        return count < size || value > pqueue.top().X;
    }
};

// TODO: add ability to reuse ids
template<typename ValueT, typename Move> struct XCandidatesGroup {
    int size;
    int count;
    VC<Candidate<ValueT, Move>> candidates;
    DEPQueue<pair<ValueT, int>> depqueue;

    XCandidatesGroup(int size = 0) : size(size) {
        candidates.resize(size);
        depqueue.reset(size);
        count = 0;
    }

    void clear() {depqueue.clear(); count = 0;}

    ValueT min() {return depqueue.top_min().X;}
    ValueT max() {return depqueue.top_max().X;}

    bool add(const Candidate<ValueT, Move>& cand) {
        if (count < size) {
            candidates[count] = cand;
            depqueue.push(MP(cand.value, count));
            count++;
            return true;
        } else {
            if (cand.value > depqueue.top_min().X) {
                int pos = depqueue.top_min().Y;
                candidates[pos] = cand;
                depqueue.replace_min(MP(cand.value, pos)); 
                return true;
            } else {
                return false;
            }
        }
    }

    bool can_add(const ValueT value) {
        return count < size || value > depqueue.top_min().X;
    }

    Candidate<ValueT, Move> pop_max() {
        assert(depqueue.size > 0);
        int pos = depqueue.pop_max().Y;
        return candidates[pos];
    }
};

// SOLUTION

const double GLOBAL_SCALE = 1.0;

const double TIME_LIMIT = 2.9 * GLOBAL_SCALE;

const int MAX_N = 30;
const int MAX_M = MAX_N * MAX_N / 2;
const int MAX_V = 15;

const int dr[] = {0, 1, 0, -1, 0};
const int dc[] = {1, 0, -1, 0, 0};
const char dn[] = {'R', 'D', 'L', 'U', '.'};

int P1D(int r, int c) {return r * MAX_N + c;}
int P1D(PII &p) {return P1D(p.X, p.Y);}
PII P2D(int p) {return MP(p / MAX_N, p % MAX_N);}

int N, M, V;

int src[MAX_N][MAX_N];
int trg[MAX_N][MAX_N];

VPII vsrc;
VPII vtrg;

bool next_combination(VI &a, const int max_value) {
    int pos = a.S - 1;
    while (pos >= 0) {
        a[pos]++;
        if (a[pos] < max_value) return true;
        a[pos] = 0;
        pos--;
    }
    return false;
}

int dist(int r0, int c0, int r1, int c1) {return abs(r0 - r1) + abs(c0 - c1);}
int dist(PII &p0, PII &p1) {return dist(p0.X, p0.Y, p1.X, p1.Y);}
int total_rotation(VI &arm_pos) {int total_rot = 0; REP(i, arm_pos.S) total_rot += arm_pos[i]; return total_rot & 3;}
int total_rotation_short(VI &arm_pos) {int total_rot = 0; REP(i, (int)(arm_pos.S - 1)) total_rot += arm_pos[i]; return total_rot & 3;}

void go_reachable(VVI &reachable, VI &arm, VI &hand, int lev, int r, int c) {
    if (lev == arm.S) {
        for (int h : hand) REP(d, 4) {
            int rr = r + dr[d] * h;
            int cc = c + dc[d] * h;
            if (rr < 0 || rr >= N || cc < 0 || cc >= N) continue;
            reachable[rr][cc] = 1;
        }
        return;
    }
    REP(d, 4) go_reachable(reachable, arm, hand, lev + 1, r + dr[d] * arm[lev], c + dc[d] * arm[lev]);
}


PII calc_arm_pos(int r0, int c0, VI &arm, VI &arm_pos) {
    int r = r0;
    int c = c0;
    int total_rot = 0;
    REP(i, arm.S) {
        total_rot += arm_pos[i];
        total_rot &= 3;
        r += arm[i] * dr[total_rot];
        c += arm[i] * dc[total_rot];
    }

    return MP(r, c);
}

PII calc_arm_pos_short(int r0, int c0, VI &arm, VI &arm_pos) {
    int r = r0;
    int c = c0;
    int total_rot = 0;
    REP(i, arm.S - 1) {
        total_rot += arm_pos[i];
        total_rot &= 3;
        r += arm[i] * dr[total_rot];
        c += arm[i] * dc[total_rot];
    }

    return MP(r, c);
}

Timer timer_bsmoves;
Timer timer_bsapply;
Timer timer_bsactive;

int checks0 = 0;
int checks1 = 0;
int checks2 = 0;

struct SMove {
    char pos[MAX_V];
    char grab[MAX_V];
    int r, c;
};

struct Solution {
    VPII build;
    VC<SMove> vm;
};

Solution last_solution;

void print_solution(VI &arm, Solution &sol) {
    int r = sol.vm[0].r;
    int c = sol.vm[0].c;
    int pos[MAX_V]; ZERO(pos);

    cout << V << endl;
    for (PII &p : sol.build) cout << p.X << " " << p.Y << endl;
    cout << r << ' ' << c << endl;

    // preprocess pos = -1
    REP(i, sol.vm.S) {
        REP(j, V-1) {
            if (sol.vm[i].pos[j] != -1) continue;
            int prev_rot = i ? sol.vm[i-1].pos[j] : 0;
            int next_rot = i + 1 < sol.vm.S ? sol.vm[i+1].pos[j] : -1;

            if (next_rot == -1) {
                sol.vm[i].pos[j] = prev_rot;
            } else if (prev_rot != (next_rot ^ 2)) {
                sol.vm[i].pos[j] = prev_rot;
            } else {
                sol.vm[i].pos[j] = (prev_rot + 3) & 3;
            }
        }
    }

    for (SMove &m : sol.vm) {
        int turns = max(1, dist(r, c, m.r, m.c));
        REP(i, V-1) if (pos[i] == (m.pos[i] ^ 2)) turns = max(turns, 2);

        VS vs(turns, string(2 * V, '.'));

        REP(i, turns) {
            int bd = 4;
            int bv = dist(r, c, m.r, m.c);
            REP(d, 4) {
                int nr = r + dr[d];
                int nc = c + dc[d];
                if (nr < 0 || nr >= N || nc < 0 || nc >= N) continue;
                int v = dist(nr, nc, m.r, m.c);
                if (v < bv) {
                    bv = v;
                    bd = d;
                }
            }
            c += dc[bd];
            r += dr[bd];
            vs[i][0] = dn[bd];
        }
        assert(r == m.r && c == m.c);

        REP(i, V-1) {
            int start_rot = pos[i];
            int end_rot = m.pos[i];

            if (((end_rot + 1) & 3) == start_rot) {
                vs[0][i + 1] = 'L';
            } else if (((end_rot + 3) & 3) == start_rot) {
                vs[0][i + 1] = 'R';
            } else if (((end_rot + 2) & 3) == start_rot) {
                vs[0][i + 1] = 'R';
                vs[1][i + 1] = 'R';
            }
        }

        FOR(i, arm.S, V-1) if (m.grab[i]) vs[turns - 1][V + 1 + i] = 'P';

        for (string &s : vs) cout << s << endl;

        REP(i, V-1) pos[i] = m.pos[i];
    }
}

int zobrist[MAX_N][MAX_N][2];

double cell_value_cached[MAX_N][MAX_N];

struct BState {
    int r, c;
    int M_left;
    bitset<MAX_N*MAX_N> src;
    bitset<MAX_N*MAX_N> trg;
    char pos[MAX_V];
    char hold[MAX_V];
    int state_id;
    int move_id;
    int hash;
    double value;

    void apply_move(SMove &m, BState &bs, VI &arm, VI &hand, VI &bhand) {
        bs.r = m.r;
        bs.c = m.c;
        bs.M_left = M_left;
        bs.src = src;
        bs.trg = trg;
        bs.value = value;
        
        bs.hash = hash;
        bs.hash ^= zobrist[r][c][1];
        bs.hash ^= zobrist[bs.r][bs.c][1];
        REP(i, V-1) bs.hash ^= zobrist[i+1][1+pos[i]][1];
        REP(i, V-1) bs.hash ^= zobrist[i+1][1+m.pos[i]][1];

        REP(i, V-1) bs.pos[i] = m.pos[i];
        REP(i, V-1) bs.hold[i] = hold[i];

        VI arm_pos(m.pos, m.pos + arm.S);
        PII center = calc_arm_pos(m.r, m.c, arm, arm_pos);
        int total_rot = total_rotation(arm_pos);
        REP(i, hand.S) {
            if (!m.grab[i+arm.S]) continue;
            int nr = center.X + dr[(total_rot + m.pos[i+arm.S]) & 3] * hand[i];
            int nc = center.Y + dc[(total_rot + m.pos[i+arm.S]) & 3] * hand[i];
            assert(nr >= 0 && nr < N && nc >= 0 && nc < N);
            if (bs.hold[i+arm.S]) {
                assert(bs.trg[P1D(nr, nc)]);
                bs.trg.reset(P1D(nr, nc));
                bs.hash ^= zobrist[nr][nc][0];
                bs.value += cell_value_cached[nr][nc];
                // bs.value += 1;
            } else {
                assert(bs.src[P1D(nr, nc)]);
                bs.src.reset(P1D(nr, nc));
                bs.hash ^= zobrist[nr][nc][0];
                bs.value += cell_value_cached[nr][nc];
                // bs.value += 1;
            }
            bs.hash ^= zobrist[0][i][1];
            bs.hold[i+arm.S] = 1 - bs.hold[i+arm.S];
            bs.M_left--;
        }

        if (bhand.S) {
            PII center = calc_arm_pos_short(m.r, m.c, arm, arm_pos);
            int total_rot = total_rotation_short(arm_pos);
            REP(i, bhand.S) {
                if (!m.grab[i+arm.S+hand.S]) continue;
                int nr = center.X + dr[(total_rot + m.pos[i+arm.S+hand.S]) & 3] * bhand[i];
                int nc = center.Y + dc[(total_rot + m.pos[i+arm.S+hand.S]) & 3] * bhand[i];
                assert(nr >= 0 && nr < N && nc >= 0 && nc < N);
                if (bs.hold[i+arm.S+hand.S]) {
                    assert(bs.trg[P1D(nr, nc)]);
                    bs.trg.reset(P1D(nr, nc));
                    bs.hash ^= zobrist[nr][nc][0];
                    bs.value += cell_value_cached[nr][nc];
                    // bs.value += 1;
                } else {
                    assert(bs.src[P1D(nr, nc)]);
                    bs.src.reset(P1D(nr, nc));
                    bs.hash ^= zobrist[nr][nc][0];
                    bs.value += cell_value_cached[nr][nc];
                    // bs.value += 1;
                }
                bs.hash ^= zobrist[0][i+hand.S][1];
                bs.hold[i+arm.S+hand.S] = 1 - bs.hold[i+arm.S+hand.S];
                bs.M_left--;
            }

        }
    }
};

const double EPSILON = 1e-6;

const int MAX_STATES = int(900 * 1000 * GLOBAL_SCALE);

int n_move_history = 0;
int n_state = 0;

XCandidatesGroup<double, SMove> cg[400];
MoveHistory<SMove> move_history[MAX_STATES];
BState beam[MAX_STATES];

const int HASH_SIZE = 1 << 25;
bitset<HASH_SIZE> hashes;

const int MAX_CONFLICTS = 1;

int last_best_state;

int solve_bs(VI arm, VI hand, VI bhand, int max_turns, double time_limit, double scale, bool save_output = false) {
    DB(arm.S, hand.S, bhand.S, V);

    n_move_history = 0;
    n_state = 0;

    if (save_output) {
        last_solution.build.clear();
        last_solution.vm.clear();
        REP(i, arm.S) last_solution.build.PB(MP(i, arm[i]));
        REP(i, hand.S) last_solution.build.PB(MP(arm.S, hand[i]));
        REP(i, bhand.S) last_solution.build.PB(MP(arm.S - 1, bhand[i]));
    }

    max_turns = min(max_turns, 390);
    const int BEAM_WIDTH = MAX_STATES / max_turns - 10;
    DB(BEAM_WIDTH);

    REP(i, max_turns+5) cg[i] = XCandidatesGroup<double, SMove>((int)(scale * BEAM_WIDTH * (V == 5 ? 50 : V <= 8 ? 16 : 12) / 10));
    
    int max_hand = 0; for (int h : hand) max_hand = max(max_hand, h);

    VVVI src_cnt(N + 2 * max_hand, VVI(N + 2 * max_hand, VI(max_hand + 1)));
    VVVI trg_cnt(N + 2 * max_hand, VVI(N + 2 * max_hand, VI(max_hand + 1)));

    VVI good_center(N + 2 * max_hand, VI(N + 2 * max_hand));
    REP(r, N) REP(c, N) if (src[P1D(r,c)]) REP(d, 4) for (int h : hand) {
        int nr = r + dr[d] * h;
        int nc = c + dc[d] * h;
        if (nr >= -max_hand && nr < N + max_hand && nc >= -max_hand && nc < N + max_hand) src_cnt[nr+max_hand][nc+max_hand][h] = 1, good_center[nr+max_hand][nc+max_hand] = 1;
    }

    REP(r, N) REP(c, N) if (trg[P1D(r,c)]) REP(d, 4) for (int h : hand) {
        int nr = r + dr[d] * h;
        int nc = c + dc[d] * h;
        if (nr >= -max_hand && nr < N + max_hand && nc >= -max_hand && nc < N + max_hand) trg_cnt[nr+max_hand][nc+max_hand][h] = 1, good_center[nr+max_hand][nc+max_hand] = 1;
    }


    VVI center_calc(N + 2 * max_hand, VI(N + 2 * max_hand));
    VVD center_value(N + 2 * max_hand, VD(N + 2 * max_hand));

    VC<pair<double,PII>> vstart;
    REP(r, N) REP(c, N) vstart.PB(MP(dist(r, c, N/2, N/2), MP(r, c)));
    sort(ALL(vstart));

    int n_beam = 0;
    for (auto &p : vstart) {
        int r = p.Y.X;
        int c = p.Y.Y;
        BState &bs = beam[n_state++];
        bs.r = r;
        bs.c = c;
        bs.M_left = M * 2;
        bs.src.reset(); for (PII &p : vsrc) bs.src[P1D(p)] = 1;
        bs.trg.reset(); for (PII &p : vtrg) bs.trg[P1D(p)] = 1;
        bs.move_id = -1;
        bs.hash = zobrist[r][c][1];
        bs.value = 0;
        REP(r, N) REP(c, N) if (bs.src[P1D(r,c)] || bs.trg[P1D(r,c)]) bs.hash ^= zobrist[r][c][0];
        REP(i, V-1) bs.hash ^= zobrist[i+1][1][1];
        ZERO(bs.pos);
        ZERO(bs.hold);

        n_beam++;
        if (n_beam == BEAM_WIDTH) break;
    }

    int best_state = -1;
    int best_turns = 1<<20;

    int loop = -1;

    hashes.reset();

    double bs_start = elapsed();

    REP(r, N) REP(c, N) {
        int min_dist = 1<<20;
        int cnt = 0;
        REP(d, 4) for (int h : hand) {
            int nr = r + dr[d] * h;
            int nc = c + dc[d] * h;
            if (nr < 0 || nr >= N || nc < 0 || nc >= N) continue;
            cnt += src[nr][nc] || trg[nr][nc];
        }

        double pos_value = 1.0 + 3.0 * (abs(N*.5-.5-r) + abs(N*.5-.5-c)) / (N + N);

        cell_value_cached[r][c] = pos_value / (1.5 + pow(cnt, .3));
        // cell_value_cached[r][c] = cell_value(r, c);
    }

    while (true) {
        loop++;
        if (n_state == MAX_STATES) break;
        if (elapsed() > time_limit) break;

        int level = -1;
        int target_n_beam = 1;
        if (loop) {
            double time_passed = elapsed() - bs_start;
            double avg_per_move = time_passed / n_state;
            double time_left = time_limit - elapsed();
            double expected_moves = time_left / avg_per_move / min(best_turns, max_turns);
            target_n_beam = max(1, min(BEAM_WIDTH, (int)(expected_moves * .2)));
        }
        while (true) {
            level++;
            if (elapsed() > time_limit) break;
            if (n_state == MAX_STATES) break;

            if (loop && level == 0) continue;
            if (level > max_turns) break;
            if (level >= best_turns) break;            

            if (level) {
                timer_bsapply.start();
                n_beam = 0;
                while (true) {
                    if (cg[level].depqueue.size == 0) break;

                    auto c = cg[level].pop_max();
                    beam[c.state_id].apply_move(c.move, beam[n_state], arm, hand, bhand);

                    if (hashes[beam[n_state].hash]) continue;
                    hashes.set(beam[n_state].hash);

                    MoveHistory<SMove> &h = move_history[n_move_history];
                    h.previous_id = beam[c.state_id].move_id;
                    h.move = c.move;
                    beam[n_state].move_id = n_move_history;
                    beam[n_state].state_id = n_state;
                    n_move_history++;
                    n_state++;
                    n_beam++;

                    if (n_beam == target_n_beam) break;
                    if (n_state == MAX_STATES) break;
                }
                timer_bsapply.stop();
            }

            REP(bid, n_beam) {
                int id = n_state - n_beam + bid;
                BState &bs = beam[id];
                if (bs.M_left == 0) {
                    if (level < best_turns) {
                        DB(loop, level, elapsed());
                        best_turns = level;
                        best_state = id;
                    }
                    break;
                }

                static int loop_id = 1;
                if (V > 9) {
                    timer_bsactive.start();
                    loop_id++;

                    REP(r, N) REP(c, N) if (bs.src[P1D(r,c)]) REP(d, 4) for (int h : hand) {
                        int nr = r + dr[d] * h;
                        int nc = c + dc[d] * h;
                        if (nr >= -max_hand && nr < N + max_hand && nc >= -max_hand && nc < N + max_hand) src_cnt[nr+max_hand][nc+max_hand][h] = loop_id;
                    }

                    REP(r, N) REP(c, N) if (bs.trg[P1D(r,c)]) REP(d, 4) for (int h : hand) {
                        int nr = r + dr[d] * h;
                        int nc = c + dc[d] * h;
                        if (nr >= -max_hand && nr < N + max_hand && nc >= -max_hand && nc < N + max_hand) trg_cnt[nr+max_hand][nc+max_hand][h] = loop_id;
                    }

                    timer_bsactive.stop();
                }

                timer_bsmoves.start();

                static VPII vpos; 
                vpos.clear(); 
                REP(r, N) REP(c, N) if (dist(bs.r, bs.c, r, c) <= 1) vpos.PB(MP(r, c));

                static VI bad_hand(hand.S); REP(i, hand.S) bad_hand[i] = bs.pos[arm.S + i] ^ 2;
                static VI good_hand_lo(hand.S); REP(i, hand.S) good_hand_lo[i] = bs.pos[arm.S + i] == -1 ? 0 : bs.pos[arm.S + i] + 3;
                static VI good_hand_hi(hand.S); REP(i, hand.S) good_hand_hi[i] = bs.pos[arm.S + i] == -1 ? 4 : bs.pos[arm.S + i] + 6;
                static VI bad_bhand(bhand.S); REP(i, bhand.S) bad_bhand[i] = bs.pos[arm.S + hand.S + i] ^ 2;
                static VI good_bhand_lo(bhand.S); REP(i, bhand.S) good_bhand_lo[i] = bs.pos[arm.S + hand.S + i] == -1 ? 0 : bs.pos[arm.S + hand.S + i] + 3;
                static VI good_bhand_hi(bhand.S); REP(i, bhand.S) good_bhand_hi[i] = bs.pos[arm.S + hand.S + i] == -1 ? 4 : bs.pos[arm.S + hand.S + i] + 6;

                static VI barm_pos(arm.S, 0); REP(i, arm.S) barm_pos[i] = 0;
                static VI arm_pos(arm.S, 0);
                do {
                    REP(i, arm.S) arm_pos[i] = barm_pos[i] == (bs.pos[i] ^ 2) ? 3 : barm_pos[i];

                    PII offset = calc_arm_pos(0, 0, arm, arm_pos);
                    int total_rot = total_rotation(arm_pos);
                    PII boffset = bhand.S ? calc_arm_pos_short(0, 0, arm, arm_pos) : MP(-1, -1);
                    int btotal_rot = bhand.S ? total_rotation_short(arm_pos) : -1;

                    const int doffset = rng.next(4);


                    for (PII &p : vpos) {
                        PII center = MP(offset.X + p.X, offset.Y + p.Y);
                        if (center.X < -max_hand || center.X >= N + max_hand || center.Y < -max_hand || center.Y >= N + max_hand) continue;
                        PII bcenter;
                        if (bhand.S) {
                            bcenter = MP(boffset.X + p.X, boffset.Y + p.Y);
                            if (bcenter.X < -max_hand || bcenter.X >= N + max_hand || bcenter.Y < -max_hand || bcenter.Y >= N + max_hand) continue;
                        }

                        // if (!good_center[center.X+max_hand][center.Y+max_hand]) continue;

                        if (V > 9) {
                            if (center_calc[center.X+max_hand][center.Y+max_hand] != loop_id) {
                                center_calc[center.X+max_hand][center.Y+max_hand] = loop_id;
                                center_value[center.X+max_hand][center.Y+max_hand] = 0;

                                int cnt = 0;
                                REP(i, hand.S) {
                                    int h = hand[i];
                                    if (bs.hold[arm.S+i]) {
                                        cnt += (trg_cnt[center.X+max_hand][center.Y+max_hand][h] == loop_id);
                                    } else {
                                        cnt += (src_cnt[center.X+max_hand][center.Y+max_hand][h] == loop_id);
                                    }
                                }
                                center_value[center.X+max_hand][center.Y+max_hand] = cnt * cell_value_cached[max(0, min(N-1, center.X))][max(0, min(N-1, center.Y))];
                            }
                            double takos = center_value[center.X+max_hand][center.Y+max_hand];
                            if (takos == 0) continue;

                            double av = bs.value + takos + rng.next_double() * EPSILON;
                            if (!cg[level + 1].can_add(av)) continue;
                        }

                        double real_takos = 0;

                        static int n_used_cell = 0; n_used_cell++;
                        static int used_cell[MAX_N][MAX_N]; if (n_used_cell == 1) ZERO(used_cell);

                        static SMove m;

                        REP(i, hand.S) {
                            m.grab[arm.S+i] = 0;
                            m.pos[arm.S+i] = -1;
                            int h = hand[i];
                            int bd = -1;
                            int br = -1, bc = -1;
                            double bv = 0;
                            if (!bs.hold[arm.S+i]) {
                                if (src_cnt[center.X+max_hand][center.Y+max_hand][h]==loop_id) {
                                    FOR(d, good_hand_lo[i], good_hand_hi[i]) {
                                        int nr = center.X + dr[(d + total_rot) & 3] * h;
                                        int nc = center.Y + dc[(d + total_rot) & 3] * h;
                                        if (nr < 0 || nr >= N || nc < 0 || nc >= N) continue;
                                        if (!bs.src[P1D(nr,nc)]) continue;
                                        double av = cell_value_cached[nr][nc];
                                        if (bd == -1 || av > bv) bd = d, bv = av, br = nr, bc = nc;
                                    }
                                }
                            } else {
                                if (trg_cnt[center.X+max_hand][center.Y+max_hand][h]==loop_id) {
                                    FOR(d, good_hand_lo[i], good_hand_hi[i]) {
                                        int nr = center.X + dr[(d + total_rot) & 3] * h;
                                        int nc = center.Y + dc[(d + total_rot) & 3] * h;
                                        if (nr < 0 || nr >= N || nc < 0 || nc >= N) continue;
                                        if (!bs.trg[P1D(nr,nc)]) continue;
                                        double av = cell_value_cached[nr][nc];
                                        if (bd == -1 || av > bv) bd = d, bv = av, br = nr, bc = nc;
                                    }
                                }
                            }

                            if (bd != -1) {
                                real_takos += bv;
                                m.pos[arm.S+i] = bd & 3;
                                m.grab[arm.S+i] = 1;
                                used_cell[br][bc] = n_used_cell;
                            }
                        }

                        REP(i, bhand.S) {
                            m.grab[arm.S+hand.S+i] = 0;
                            m.pos[arm.S+hand.S+i] = -1;
                            int h = bhand[i];
                            int bd = -1;
                            int br = -1, bc = -1;
                            double bv = 0;
                            if (!bs.hold[arm.S+hand.S+i]) {
                                if (src_cnt[bcenter.X+max_hand][bcenter.Y+max_hand][h]==loop_id) {
                                    FOR(d, good_bhand_lo[i], good_bhand_hi[i]) {
                                        int nr = bcenter.X + dr[(d + btotal_rot) & 3] * h;
                                        int nc = bcenter.Y + dc[(d + btotal_rot) & 3] * h;
                                        if (nr < 0 || nr >= N || nc < 0 || nc >= N) continue;
                                        if (!bs.src[P1D(nr,nc)]) continue;
                                        if (used_cell[nr][nc] == n_used_cell) continue;
                                        double av = cell_value_cached[nr][nc];
                                        if (bd == -1 || av > bv) bd = d, bv = av, br = nr, bc = nc;
                                    }
                                }
                            } else {
                                if (trg_cnt[bcenter.X+max_hand][bcenter.Y+max_hand][h]==loop_id) {
                                    FOR(d, good_bhand_lo[i], good_bhand_hi[i]) {
                                        int nr = bcenter.X + dr[(d + btotal_rot) & 3] * h;
                                        int nc = bcenter.Y + dc[(d + btotal_rot) & 3] * h;
                                        if (nr < 0 || nr >= N || nc < 0 || nc >= N) continue;
                                        if (!bs.trg[P1D(nr,nc)]) continue;
                                        if (used_cell[nr][nc] == n_used_cell) continue;
                                        double av = cell_value_cached[nr][nc];
                                        if (bd == -1 || av > bv) bd = d, bv = av, br = nr, bc = nc;
                                    }
                                }
                            }

                            if (bd != -1) {
                                real_takos += bv;
                                m.pos[arm.S+hand.S+i] = bd & 3;
                                m.grab[arm.S+hand.S+i] = 1;
                                used_cell[br][bc] = n_used_cell;
                            }
                        }

                        // if (real_takos == 0) continue;
                        REP(i, arm.S) m.pos[i] = arm_pos[i];

                        double av2 = bs.value + real_takos + rng.next_double() * EPSILON;

                        m.r = p.X;
                        m.c = p.Y;

                        cg[level+1].add({id, av2, m});
                    }
                } while (next_combination(barm_pos, 3));
                timer_bsmoves.stop();
            }
        }
    }

    DATA(loop);
    DATAX(n_move_history, "moves");

    DB(loop, best_turns, n_move_history, best_state);

    if (best_state == -1) return 1 << 20;
    assert(best_state != -1);
    if (save_output) {
        int move_id = beam[best_state].move_id;
        VC<SMove> moves;
        while (move_id != -1) {
            moves.PB(move_history[move_id].move);
            move_id = move_history[move_id].previous_id;
        }
        reverse(ALL(moves));
        last_solution.vm = moves;
        DB(last_solution.vm.S);
    }

    last_best_state = best_state;

    return best_turns;
}

const int VALUE_EARLY_EXIT = 1 << 21;
const int VALUE_NON_REACHABLE = 1 << 22;

VC<SMove> last_output;
int solve_greedy(int r0, int c0, int r1, int c1, VI arm, VI hand, int total_turns_cutoff, bool save_output = false) {
    if (save_output) {
        last_solution.build.clear();
        last_solution.vm.clear();
        REP(i, arm.S) last_solution.build.PB(MP(i, arm[i]));
        REP(i, hand.S) last_solution.build.PB(MP(arm.S, hand[i]));
    }

    // check if all are reachable;
    VVI reachable(N, VI(N));
    FOR(r, r0, r1+1) FOR(c, c0, c1+1) go_reachable(reachable, arm, hand, 0, r, c);
    REP(r, N) REP(c, N) if (src[r][c] != trg[r][c] && !reachable[r][c]) return VALUE_NON_REACHABLE;

    // sim
    int max_hand = 0; for (int h : hand) max_hand = max(max_hand, h);

    VI cur_arm(arm.S, 0);
    VI cur_hand(hand.S, 0);
    VI cur_hold(hand.S, 0);
    int cur_r = (r0+r1)/2;
    int cur_c = (c0+c1)/2;
    int total_turns = 0;

    VI hand_sizes(max_hand + 1); for (int h : hand) hand_sizes[h]++;

    VVI src_active(N, VI(N)); for (PII &p : vsrc) src_active[p.X][p.Y] = 1;
    VVI trg_active(N, VI(N)); for (PII &p : vtrg) trg_active[p.X][p.Y] = 1;

    VVVI src_cnt(N + 2 * max_hand, VVI(N + 2 * max_hand, VI(max_hand + 1)));
    VVVI trg_cnt(N + 2 * max_hand, VVI(N + 2 * max_hand, VI(max_hand + 1)));

    int M_left = M;

    while (M_left) {
        if (elapsed() > TIME_LIMIT) return 1 << 20;

        VI size_hold(max_hand + 1);
        REP(i, hand.S) if (cur_hold[i]) size_hold[hand[i]]++;

        int takos_hold = 0; for (int h : size_hold) takos_hold += h;

        static int loop_id = 0; 
        loop_id++;

        // update *_cnt
        REP(r, N) REP(c, N) if (src_active[r][c]) REP(d, 4) FOR(dist, 1, max_hand + 1) {
            int nr = r + dr[d] * dist;
            int nc = c + dc[d] * dist;
            if (nr >= -max_hand && nr < N + max_hand && nc >= -max_hand && nc < N + max_hand) src_cnt[nr+max_hand][nc+max_hand][dist] = loop_id;
        }
        REP(r, N) REP(c, N) if (trg_active[r][c]) REP(d, 4) FOR(dist, 1, max_hand + 1) {
            int nr = r + dr[d] * dist;
            int nc = c + dc[d] * dist;
            if (nr >= -max_hand && nr < N + max_hand && nc >= -max_hand && nc < N + max_hand) trg_cnt[nr+max_hand][nc+max_hand][dist] = loop_id;
        }

        // find best move
        VI best_arm_pos0;
        VI best_hand_pos0;
        PII best_pos0;
        int best_turns0;
        VI best_hand_grab0;
        
        double bv = 0;

        int moves1 = 0;
        VI cache_m1(arm.S);

        VI arm_pos0(arm.S);
        do {
            PII p0_base = calc_arm_pos(0, 0, arm, arm_pos0);
            int total_rot0 = total_rotation(arm_pos0);
            int base_turns0 = 1;
            REP(i, arm.S) if ((cur_arm[i] ^ 2) == arm_pos0[i]) base_turns0 = 2;

            FOR(m0_r, r0, r1+1) FOR(m0_c, c0, c1+1) {
                // if (M_left > M / 2 && dist(cur_r, cur_c, m0_r, m0_c) > 2) continue;

                PII center = MP(p0_base.X + m0_r, p0_base.Y + m0_c);
                if (center.X < -max_hand || center.X >= N + max_hand || center.Y < -max_hand || center.Y >= N + max_hand) continue;

                // calc moves
                int turns0 = max(base_turns0, dist(cur_r, cur_c, m0_r, m0_c));

                int takos = 0;
                for (int h : hand) {
                    if (size_hold[h]) {
                        takos += (trg_cnt[center.X+max_hand][center.Y+max_hand][h] == loop_id);
                    } else {
                        takos += (src_cnt[center.X+max_hand][center.Y+max_hand][h] == loop_id);
                    }
                }

                double av = 1.0 * takos / turns0 + rng.next_double() * 1e-9;
                if (takos == 0 || av < bv) continue;

                VI hand_pos0 = cur_hand;
                VI hand_grab0(hand.S);
                VI hold = cur_hold;

                REP(i, hand.S) {
                    int h = hand[i];
                    int bd = -1;
                    if (!size_hold[h]) {
                        if (src_cnt[center.X+max_hand][center.Y+max_hand][h]==loop_id) {
                            REP(d, 4) {
                                int nr = center.X + dr[(d + total_rot0) & 3] * h;
                                int nc = center.Y + dc[(d + total_rot0) & 3] * h;
                                if (nr < 0 || nr >= N || nc < 0 || nc >= N) continue;
                                if (!src_active[nr][nc]) continue;
                                if (bd == -1 || bd == (cur_hand[i] ^ 2)) bd = d;
                            }
                            if (bd == (cur_hand[i] ^ 2)) turns0 = max(turns0, 2);
                            assert(bd != -1);
                            hand_pos0[i] = bd;
                            hand_grab0[i] = 1;
                        }
                    } else {
                        if (trg_cnt[center.X+max_hand][center.Y+max_hand][h]==loop_id) {
                            REP(d, 4) {
                                int nr = center.X + dr[(d + total_rot0) & 3] * h;
                                int nc = center.Y + dc[(d + total_rot0) & 3] * h;
                                if (nr < 0 || nr >= N || nc < 0 || nc >= N) continue;
                                if (!trg_active[nr][nc]) continue;
                                if (bd == -1 || bd == (cur_hand[i] ^ 2)) bd = d;
                            }
                            if (bd == (cur_hand[i] ^ 2)) turns0 = max(turns0, 2);
                            assert(bd != -1);
                            hand_pos0[i] = bd;
                            hand_grab0[i] = 1;
                        }
                    }
                }

                if (av > bv) {
                    bv = av;
                    best_arm_pos0 = arm_pos0;
                    best_hand_pos0 = hand_pos0;
                    best_hand_grab0 = hand_grab0;
                    best_pos0 = MP(m0_r, m0_c);
                    best_turns0 = turns0;
                }
            }
        } while (next_combination(arm_pos0, 4));

        int grab_total = 0; 
        for (int x : best_hand_grab0) grab_total += x;
        if (grab_total == 0) return (1 << 20) + M_left;

        assert(bv > 0);

        total_turns += best_turns0;

        // DB(total_turns, best_turns0, best_pos0, best_arm_pos0, best_hand_pos0, best_hand_grab0);

        // update *_active & M_left
        int total_rot0 = total_rotation(best_arm_pos0);
        REP(i, hand.S) {
            if (best_hand_grab0[i]) {
                PII p0 = calc_arm_pos(best_pos0.X, best_pos0.Y, arm, best_arm_pos0);
                int nr0 = p0.X + dr[(total_rot0 + best_hand_pos0[i]) & 3] * hand[i];
                int nc0 = p0.Y + dc[(total_rot0 + best_hand_pos0[i]) & 3] * hand[i];
                assert(nr0 >= 0 && nr0 < N && nc0 >= 0 && nc0 < N);
                if (cur_hold[i]) {
                    assert(trg_active[nr0][nc0]);
                    trg_active[nr0][nc0] = 0;
                    M_left--;
                } else {
                    assert(src_active[nr0][nc0]);
                    src_active[nr0][nc0] = 0;
                }
                cur_hold[i] = 1 - cur_hold[i];
            }
        }
        if (total_turns > total_turns_cutoff || total_turns == total_turns_cutoff && M_left) return VALUE_EARLY_EXIT;

        // add optional output
        if (save_output) {
            SMove m;
            ZERO(m.grab);
            REP(i, arm.S) m.pos[i] = best_arm_pos0[i];
            REP(i, hand.S) m.pos[arm.S + i] = best_hand_pos0[i];
            REP(i, hand.S) m.grab[arm.S + i] = best_hand_grab0[i];
            m.r = best_pos0.X;
            m.c = best_pos0.Y;
            last_solution.vm.PB(m);
        }

        cur_c = best_pos0.Y;
        cur_r = best_pos0.X;
        cur_arm = best_arm_pos0;
        cur_hand = best_hand_pos0;
    }

    return total_turns;
}

void generate_arm_hand(int n, int v, VI &arm, VI &hand, VI &bhand) {
    arm = {(n - 7) / 2, (n - 1) / 4, (n + 7) / 8, (n + 15) / 16};
    hand = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
    bhand = {};

    if (v <= 12) {
        arm = {(n - 5) / 2, (n - 1) / 4, (n + 7) / 8};
    }

    if (v <= 8) {
        arm = {(n - 5) / 2, (n - 1) / 4};
    }

    if (v == 5) {
        arm = {(n + 3) / 3};
        int add = (n - 10) / 5;
        hand = {add + 1, add + 2, add + 3};
    }

    if (v >= 12) bhand = {1};
    if (v >= 13) bhand = {1, 2};
    if (v >= 15) bhand = {1, 2, 3};

    while (1 + arm.S + hand.S + bhand.S > V) hand.pop_back();
}


int main(int argc, char **argv) {
    cin >> N >> M >> V;


    REP(r, N) {
        string s; cin >> s;
        REP(c, N) src[r][c] = s[c] - '0';
    }
    REP(r, N) {
        string s; cin >> s;
        REP(c, N) trg[r][c] = s[c] - '0';
    }

    REP(r, N) REP(c, N) {
        if (src[r][c] == trg[r][c]) src[r][c] = trg[r][c] = 0;
        if (src[r][c] == 1) vsrc.PB(MP(r, c));
        if (trg[r][c] == 1) vtrg.PB(MP(r, c));
    }
    M = vsrc.S;

    REP(i, N) REP(j, N) REP(k, 2) zobrist[i][j][k] = rng.next() & (HASH_SIZE - 1);

    DATA(N);
    DATA(V);
    DATA(M);
    double Mrat = 1.0 * M / N / N;
    DATA(Mrat);

    VI arm, hand, bhand;
    generate_arm_hand(N, V, arm, hand, bhand);

    int base_r0 = N / 2 - 3;
    int base_c0 = N / 2 - 3;
    int base_r1 = base_r0 + 7;
    int base_c1 = base_c0 + 7;


    int xv = 1000;
    Solution best_sol;
    VI best_arm = arm;
    VI best_hand = hand;
    VI best_bhand = bhand;
    int best_r0 = base_r0;
    int best_c0 = base_c0;
    int best_r1 = base_r1;
    int best_c1 = base_c1;

    int step = 0;

    const double GREEDY_TIME_LIMIT = TIME_LIMIT * .1;

    while (true) {
        if (elapsed() > GREEDY_TIME_LIMIT) break;

        arm = best_arm;
        hand = best_hand;
        base_r0 = best_r0;
        base_c0 = best_c0;
        base_r1 = best_r1;
        base_c1 = best_c1;

        bool changed = false;
        int n_changes = 1;
        if (step == 0) n_changes = 0, changed = true;
        while (n_changes--) {
            int type = rng.next(4);
            if (type == 3 && rng.next_double() < .95) type = rng.next(3);
            if (type == 0) {
                if (V > 8) continue;
                int pos = rng.next(arm.S);
                int old = arm[pos];
                if (rng.next(2)) {
                    arm[pos] = rng.next(1, N);
                } else {
                    arm[pos] += rng.next(-2, 3);
                }
                if (arm[pos] < 1 || arm[pos] >= N || arm[pos] == old) {
                    arm[pos] = old;
                    continue;
                }
            } else if (type == 1) {
                if (V > 8) continue;
                
                int pos = rng.next(hand.S);
                int old = hand[pos];
                hand[pos] = rng.next(1, max((int)hand.S + 1, N / 2));
                bool repeat = false;
                REP(i, hand.S) REP(j, i) if (hand[i] == hand[j]) repeat = true;
                if (repeat || hand[pos] == old)
                    continue;
            } else if (type == 2) {
                if (true || rng.next(2)) {
                    int offset_r = rng.next(-2, 3);
                    int offset_c = rng.next(-2, 3);
                    if (base_r0 + offset_r < 0 || base_r1 + offset_r >= N || base_c0 + offset_c < 0 || base_c1 + offset_c >= N) continue;
                    base_r0 += offset_r;
                    base_c0 += offset_c;
                    base_r1 += offset_r;
                    base_c1 += offset_c;
                } else {
                    int side = rng.next(4);
                    int extend = rng.next(2);
                    if (side == 0) base_r0 += extend ? 1 : -1;
                    if (side == 1) base_r1 += extend ? 1 : -1;
                    if (side == 2) base_c0 += extend ? 1 : -1;
                    if (side == 3) base_c1 += extend ? 1 : -1;
                    if (base_r0 < 0 || base_r1 >= N || base_c0 < 0 || base_c1 >= N || base_r0 > base_r1 || base_c0 > base_c1) {
                        base_r0 = best_r0;
                        base_c0 = best_c0;
                        base_r1 = best_r1;
                        base_c1 = best_c1;
                        continue;
                    }
                }
            } else if (type == 3) {
                if (elapsed() < GREEDY_TIME_LIMIT * .2) continue;
                if (rng.next(2)) {
                    int pos = rng.next(arm.S);
                    arm.erase(arm.begin() + pos);
                    if (true || rng.next(2)) {
                        int lowest = 1;
                        while (true) {
                            bool repeat = false;
                            for (int h : hand) if (h == lowest) repeat = true;
                            if (!repeat) break;
                            lowest++;
                        }
                        hand.PB(lowest);
                    } else {
                        int value = rng.next(1, N);
                        bool repeat = false;
                        REP(i, hand.S) if (hand[i] == value) repeat = true;
                        if (repeat) continue;
                        hand.PB(value);
                    }
                } else {
                    int pos = rng.next(hand.S);
                    hand.erase(hand.begin() + pos);
                    arm.PB(rng.next(1, N));
                }
                if (hand.S == 0 || arm.S == 0) continue;
            }
            changed = true;
        }
        if (!changed) continue;
        sort(ALL(arm));
        sort(ALL(hand));

        int turns = solve_greedy(base_r0, base_c0, base_r1, base_c1, arm, hand, xv, true);
        if (turns <= xv) {
            bool improved = turns < xv;
            xv = turns;
            break;
        }
        step++;
    }
    step--;
    DB(xv, step, elapsed());

    const double TSCALE = 0.05;

    rng = RNG(1);
    double vv = xv;

    generate_arm_hand(max(15, N - 2), V, arm, hand, bhand);
    double bsrv0 = solve_bs(arm, hand, bhand, xv, elapsed() + TIME_LIMIT * TSCALE, TSCALE, true) + last_best_state * 1e-9;
    if (bsrv0 < vv) best_arm = arm, best_hand = hand, best_bhand = bhand, vv = bsrv0, best_sol = last_solution;

    generate_arm_hand(N, V, arm, hand, bhand);
    double bsrv1 = solve_bs(arm, hand, best_bhand, xv, elapsed() + TIME_LIMIT * TSCALE, TSCALE, true) + last_best_state * 1e-9;
    if (bsrv1 < vv) best_arm = arm, best_hand = hand, best_bhand = bhand, vv = bsrv1, best_sol = last_solution;

    generate_arm_hand(min(30, N + 2), V, arm, hand, bhand);
    double bsrv2 = solve_bs(arm, hand, bhand, xv, elapsed() + TIME_LIMIT * TSCALE, TSCALE, true) + last_best_state * 1e-9;
    if (bsrv2 < vv) best_arm = arm, best_hand = hand, best_bhand = bhand, vv = bsrv2, best_sol = last_solution;

    DB(best_arm, best_hand, best_bhand, vv, xv);
    double bsrv = solve_bs(best_arm, best_hand, best_bhand, xv, TIME_LIMIT, 1.0, true);
    if (bsrv < vv) best_sol = last_solution;


    DATAX(timer_bsmoves.total, "tm");
    DATAX(timer_bsapply.total, "tapp");
    DATAX(timer_bsactive.total, "tact");

    double time = elapsed();
    DATA(time);

    print_solution(best_arm, best_sol);

	return 0;
}
