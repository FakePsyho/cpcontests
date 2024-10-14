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

    INLINE T top() {return data[1];}
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

    ValueT last_value() {
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

// SOLUTION

const double TIME_LIMIT = 2.9;

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

Timer timer_active;
Timer timer_reach;
Timer timer_save;
Timer timer_clear;
Timer timer_c2;
Timer timer_cache;
Timer timer_m1;

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

void print_solution(Solution &sol) {
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

        REP(i, V-1) if (m.grab[i]) vs[turns - 1][V + 1 + i] = 'P';

        for (string &s : vs) cout << s << endl;

        REP(i, V-1) pos[i] = m.pos[i];
    }
}

int zobrist[MAX_N][MAX_N][2];

double cell_value(int r, int c) {
    return 1.0 + 2.0 * dist(N/2, N/2, r, c) * 1.0 / (N + N);
}

struct BState {
    int r, c;
    int M_left;
    bitset<MAX_N*MAX_N> src;
    bitset<MAX_N*MAX_N> trg;
    char pos[MAX_V];
    char hold[MAX_V];
    int move_id;
    int hash;
    double value;

    void apply_move(SMove &m, BState &bs, VI &arm, VI &hand) {
        bs.r = m.r;
        bs.c = m.c;
        bs.M_left = M_left;
        bs.src = src;
        bs.trg = trg;
        bs.value = value;
        
        bs.hash = hash;
        bs.hash ^= zobrist[r][c][1];
        bs.hash ^= zobrist[bs.r][bs.c][1];

        REP(i, V-1) bs.pos[i] = m.pos[i];
        REP(i, V-1) bs.hold[i] = hold[i];

        VI arm_pos(m.pos, m.pos + arm.S);
        PII center = calc_arm_pos(m.r, m.c, arm, arm_pos);
        int total_rot = total_rotation(arm_pos);
        FOR(i, arm.S, V-1) {
            if (!m.grab[i]) continue;
            int nr = center.X + dr[(total_rot + m.pos[i]) & 3] * hand[i-arm.S];
            int nc = center.Y + dc[(total_rot + m.pos[i]) & 3] * hand[i-arm.S];
            assert(nr >= 0 && nr < N && nc >= 0 && nc < N);
            if (bs.hold[i]) {
                assert(bs.trg[P1D(nr, nc)]);
                bs.trg.reset(P1D(nr, nc));
                bs.hash ^= zobrist[nr][nc][0];
                bs.value += cell_value(nr, nc);
                // bs.value += 1;
            } else {
                assert(bs.src[P1D(nr, nc)]);
                bs.src.reset(P1D(nr, nc));
                bs.hash ^= zobrist[nr][nc][0];
                bs.value += cell_value(nr, nc);
                // bs.value += 1;
            }
            bs.hash ^= zobrist[0][i][1];
            bs.hold[i] = 1 - bs.hold[i];
            bs.M_left--;
        }
    }
};

const double EPSILON = 1e-6;

const int MAX_TOTAL_MOVES = 300 * 1000;

int n_move_history = 0;
CandidatesGroup<double, SMove> cg[2000];
int n_beam[2000];
MoveHistory<SMove> move_history[MAX_TOTAL_MOVES];
BState beam[MAX_TOTAL_MOVES];

const int HASH_SIZE = 1 << 23;
bitset<HASH_SIZE> hashes;

const int MAX_CONFLICTS = 1;

int solve_bs(int r0, int c0, int r1, int c1, VI arm, VI hand, int width, int max_turns, bool save_output = false) {
    n_move_history = 0;
    if (save_output) {
        last_solution.build.clear();
        last_solution.vm.clear();
        REP(i, arm.S) last_solution.build.PB(MP(i, arm[i]));
        REP(i, hand.S) last_solution.build.PB(MP(arm.S, hand[i]));
    }

    const int MAX_WIDTH = MAX_TOTAL_MOVES / max_turns;

    const int BEAM_WIDTH = min(width, MAX_WIDTH);
    DB(BEAM_WIDTH);

    int max_width = BEAM_WIDTH;

    REP(i, max_turns+10) cg[i] = CandidatesGroup<double, SMove>(BEAM_WIDTH);
    
    int max_hand = 0; for (int h : hand) max_hand = max(max_hand, h);

    VVVI src_cnt(N + 2 * max_hand, VVI(N + 2 * max_hand, VI(max_hand + 1)));
    VVVI trg_cnt(N + 2 * max_hand, VVI(N + 2 * max_hand, VI(max_hand + 1)));

    VVI center_calc(N + 2 * max_hand, VI(N + 2 * max_hand));
    VVD center_value(N + 2 * max_hand, VD(N + 2 * max_hand));

    int level = 0;
    ZERO(n_beam);
    VC<pair<double,PII>> vstart;
    REP(r, N) REP(c, N) vstart.PB(MP(dist(r, c, (r0+r1)/2, (c0+c1)/2), MP(r, c)));
    sort(ALL(vstart));

    for (auto &p : vstart) {
        int r = p.Y.X;
        int c = p.Y.Y;
        BState &bs = beam[n_beam[0]++];
        bs.r = r;
        bs.c = c;
        bs.M_left = M * 2;
        bs.src.reset(); for (PII &p : vsrc) bs.src[P1D(p)] = 1;
        bs.trg.reset(); for (PII &p : vtrg) bs.trg[P1D(p)] = 1;
        bs.move_id = -1;
        bs.hash = zobrist[r][c][1];
        bs.value = 0;
        REP(r, N) REP(c, N) if (bs.src[P1D(r,c)] || bs.trg[P1D(r,c)]) bs.hash ^= zobrist[r][c][0];
        ZERO(bs.pos);
        ZERO(bs.hold);

        if (n_beam[0] == BEAM_WIDTH) break;
    }

    int fallback = 0;

    while (true) {
        hashes.reset();

        if (elapsed() > TIME_LIMIT * 0.60) max_width = max(1, BEAM_WIDTH * 3 / 4);
        if (elapsed() > TIME_LIMIT * 0.70) max_width = max(1, BEAM_WIDTH * 2 / 3);
        if (elapsed() > TIME_LIMIT * 0.80) max_width = max(1, BEAM_WIDTH / 2);
        if (elapsed() > TIME_LIMIT * 0.90) max_width = max(1, BEAM_WIDTH / 4);
        if (elapsed() > TIME_LIMIT * 0.95) max_width = max(1, BEAM_WIDTH / 8);

        if (level) {
            timer_bsapply.start();
            n_beam[level] = 0;
            REP(i, cg[level].count) {
                int id = cg[level].candidates[i].state_id;
                int olevel = id / BEAM_WIDTH;
                int opos = id % BEAM_WIDTH;
                beam[olevel*BEAM_WIDTH+opos].apply_move(cg[level].candidates[i].move, beam[level*BEAM_WIDTH+n_beam[level]], arm, hand);

                if (hashes[beam[level*BEAM_WIDTH+n_beam[level]].hash]) continue;
                hashes.set(beam[level*BEAM_WIDTH+n_beam[level]].hash);

                MoveHistory<SMove> &h = move_history[n_move_history];
                h.previous_id = beam[olevel*BEAM_WIDTH+opos].move_id;
                h.move = cg[level].candidates[i].move;
                beam[level*BEAM_WIDTH+n_beam[level]].move_id = n_move_history;
                n_beam[level]++;
                n_move_history++;
            }
            timer_bsapply.stop();
        }

        static VI order;
        order.clear();
        REP(i, n_beam[level]) order.PB(i);
        sort(ALL(order), [&](int a, int b) {return beam[level*BEAM_WIDTH+a].M_left < beam[level*BEAM_WIDTH+b].M_left;});

        if (n_beam[level]) {
            if (n_beam[level] > max_width) {
                n_beam[level] = max_width;
                order.resize(max_width);
            }
            int lo = beam[level*BEAM_WIDTH+order[0]].M_left;
            int hi = beam[level*BEAM_WIDTH+order.back()].M_left;
            DB(level, n_beam[level], lo, hi, elapsed());
        }

        if (level && n_beam[level] == 0) {
            DB(fallback);
            fallback += 2;
            level = max(0, level - 5);
            FOR(i, level+1, level+10) cg[i].clear();
            continue;
        }

        REP(bid, n_beam[level]) {
            if (bid % 10 == 0 && elapsed() > TIME_LIMIT) return 1<<20;

            BState &bs = beam[level*BEAM_WIDTH+order[bid]];
            if (bs.M_left == 0) {
                cerr << "bs done: " << level << endl;
                if (save_output) {
                    int move_id = bs.move_id;
                    VC<SMove> moves;
                    while (move_id != -1) {
                        // DB(move_id, move_history[move_id].move.r, move_history[move_id].move.c, (int)move_history[move_id].move.pos[0], (int)move_history[move_id].move.pos[1]);
                        moves.PB(move_history[move_id].move);
                        move_id = move_history[move_id].previous_id;
                    }
                    reverse(ALL(moves));
                    last_solution.vm = moves;
                }
                // REP(i, n_move_history) move_history[i].previous_id = 0, move_history[i].move = {};
                return level;
            }

            int id = level*BEAM_WIDTH + order[bid];

            timer_bsactive.start();

            static int loop_id = 0;
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

            timer_bsmoves.start();
            VI arm_pos(arm.S, 0);

            int max_turns = bs.M_left > M / 2 ? 1 : 2;
            max_turns += fallback;
            static VPII vpos; 
            vpos.clear(); 
            REP(r, N) REP(c, N) if (dist(bs.r, bs.c, r, c) <= max_turns) vpos.PB(MP(r, c));

            do {
                PII offset = calc_arm_pos(0, 0, arm, arm_pos);
                int total_rot = total_rotation(arm_pos);
                int base_turns = 1;
                REP(i, arm.S) if ((bs.pos[i] ^ 2) == arm_pos[i]) base_turns = 2;
                if (base_turns > max_turns) continue;

                for (PII &p : vpos) {
                    PII center = MP(offset.X + p.X, offset.Y + p.Y);
                    if (center.X < -max_hand || center.X >= N + max_hand || center.Y < -max_hand || center.Y >= N + max_hand) continue;

                    int turns = max(base_turns, dist(bs.r, bs.c, p.X, p.Y));

                    if (center_calc[center.X+max_hand][center.Y+max_hand] != loop_id) {
                        center_calc[center.X+max_hand][center.Y+max_hand] = loop_id;
                        center_value[center.X+max_hand][center.Y+max_hand] = 0;
                        REP(i, hand.S) {
                            int h = hand[i];
                            if (bs.hold[arm.S+i]) {
                                center_value[center.X+max_hand][center.Y+max_hand] += (trg_cnt[center.X+max_hand][center.Y+max_hand][h] == loop_id) * cell_value(center.X, center.Y);
                            } else {
                                center_value[center.X+max_hand][center.Y+max_hand] += (src_cnt[center.X+max_hand][center.Y+max_hand][h] == loop_id) * cell_value(center.X, center.Y);
                            }
                        }
                    }
                    double takos = center_value[center.X+max_hand][center.Y+max_hand];
                    if (takos == 0) continue;

                    double real_takos = 0;

                    double av = bs.value + takos + rng.next_double() * EPSILON;
                    if (!cg[level + turns].can_add(av)) continue;

                    SMove m;

                    REP(i, V-1) m.grab[i] = 0;
                    REP(i, arm.S) m.pos[i] = arm_pos[i];
                    FOR(i, arm.S, V-1) m.pos[i] = -1;

                    REP(i, hand.S) {
                        int h = hand[i];

                        int bd = -1;
                        if (!bs.hold[arm.S+i]) {
                            if (src_cnt[center.X+max_hand][center.Y+max_hand][h]==loop_id) {
                                REP(d, 4) {
                                    int nr = center.X + dr[(d + total_rot) & 3] * h;
                                    int nc = center.Y + dc[(d + total_rot) & 3] * h;
                                    if (nr < 0 || nr >= N || nc < 0 || nc >= N) continue;
                                    if (!bs.src[P1D(nr,nc)]) continue;
                                    if (bd == -1 || bd == (bs.pos[arm.S + i] ^ 2)) bd = d;
                                }
                                if (bd == (bs.pos[arm.S + i] ^ 2)) turns = max(turns, 2);
                                real_takos += cell_value(center.X + dr[(bd + total_rot) & 3] * h, center.Y + dc[(bd + total_rot) & 3] * h);
                                assert(bd != -1);
                                m.pos[arm.S+i] = bd;
                                m.grab[arm.S+i] = 1;
                            }
                        } else {
                            if (trg_cnt[center.X+max_hand][center.Y+max_hand][h]==loop_id) {
                                REP(d, 4) {
                                    int nr = center.X + dr[(d + total_rot) & 3] * h;
                                    int nc = center.Y + dc[(d + total_rot) & 3] * h;
                                    if (nr < 0 || nr >= N || nc < 0 || nc >= N) continue;
                                    if (!bs.trg[P1D(nr,nc)]) continue;
                                    if (bd == -1 || bd == (bs.pos[arm.S + i] ^ 2)) bd = d;
                                }
                                if (bd == (bs.pos[arm.S + i] ^ 2)) turns = max(turns, 2);
                                real_takos += cell_value(center.X + dr[(bd + total_rot) & 3] * h, center.Y + dc[(bd + total_rot) & 3] * h);
                                assert(bd != -1);
                                m.pos[arm.S+i] = bd;
                                m.grab[arm.S+i] = 1;
                            }
                        }
                    }

                    av = bs.value + real_takos + rng.next_double() * EPSILON;

                    int total_turns = level + turns;

                    m.r = p.X;
                    m.c = p.Y;

                    cg[total_turns].add({id, av, m});
                }
                // DB(level, cg[level+1].count);
            } while (next_combination(arm_pos, 4));
            timer_bsmoves.stop();
        }
        level++;
    }
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
    timer_reach.start();
    VVI reachable(N, VI(N));
    FOR(r, r0, r1+1) FOR(c, c0, c1+1) go_reachable(reachable, arm, hand, 0, r, c);
    timer_reach.stop();
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
        timer_active.start();
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
        timer_active.stop();

        // find best move
        VI best_arm_pos0;
        VI best_hand_pos0;
        PII best_pos0;
        int best_turns0;
        VI best_hand_grab0;
        
        double bv = 0;

        timer_cache.start();
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
        if (src[r][c] == trg[r][c]) continue;
        if (src[r][c] == 1) vsrc.PB(MP(r, c));
        if (trg[r][c] == 1) vtrg.PB(MP(r, c));
    }
    M = vsrc.S;

    // if (V < 10) {
    //     cout << "1\n-1 -1" << endl;
    //     exit(0);
    // }

    REP(i, N) REP(j, N) REP(k, 2) zobrist[i][j][k] = rng.next() & (HASH_SIZE - 1);

    DATA(N);
    DATA(V);
    DATA(M);
    double Mrat = 1.0 * M / N / N;
    DATA(Mrat);

    VI arm = {(N - 1) / 2, (N + 3) / 4, (N + 7) / 8, (N + 15) / 16};
    VI hand = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};

    int base_r0 = N / 2 - 1;
    int base_c0 = N / 2 - 1;
    int base_r1 = base_r0 + 2;
    int base_c1 = base_c0 + 2;

    if (V <= 11) {
        arm = {(N + 1) / 2, (N + 3) / 4, (N + 7) / 8};
        base_r0 = N / 2 - 2;
        base_c0 = N / 2 - 2;
        base_r1 = base_r0 + 2;
        base_c1 = base_c0 + 3;
    }

    if (V <= 8) {
        arm = {(N + 1) / 2, (N + 3) / 4};
        base_r0 = N / 2 - 3;
        base_c0 = N / 2 - 3;
        base_r1 = base_r0 + 7;
        base_c1 = base_c0 + 7;
    }

    arm[0]--;
    arm[1]--;

    while (1 + arm.S + hand.S > V) hand.pop_back();



    int xv = 1 << 20;
    Solution best_sol;
    VI best_arm = arm;
    VI best_hand = hand;
    int best_r0 = base_r0;
    int best_c0 = base_c0;
    int best_r1 = base_r1;
    int best_c1 = base_c1;

    int step = 0;
    int failed = 0;

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
                    // arm[pos] += arm[pos] >= old;
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
        if (turns >= (1<<20) && turns < (1<<21)) failed++;
        if (turns == VALUE_NON_REACHABLE) step = max(step-1, 1);
        if (turns <= xv) {
            bool improved = turns < xv;
            xv = turns;
            best_sol = last_solution;
            best_arm = arm;
            best_hand = hand;
            best_r0 = base_r0;
            best_c0 = base_c0;
            best_r1 = base_r1;
            best_c1 = base_c1;
            if (improved) {
                VI region = {base_r0, base_c0, base_r1, base_c1};
                DB(step, turns, elapsed(), arm, hand, region);
            }
        }
        step++;
    }
    step--;
    failed--;

    DATA(step);
    // DATA(failed);

    // print_solution(best_sol);
    // exit(0);

    // best_arm = arm;
    // best_hand = hand;

    rng = RNG(1);
    int bsrv = solve_bs(best_r0, best_c0, best_r1, best_c1, best_arm, best_hand, step * 10, xv, true);

    int use_greedy = 0;
    if (bsrv < xv) {
        best_sol = last_solution;
        xv = bsrv;
    } else {
        use_greedy = 1;
    }

    // int width = step * 5;
    // while (elapsed() < TIME_LIMIT * .5) {
    //     width *= 2;
    //     bsrv = solve_bs(best_r0, best_c0, best_r1, best_c1, best_arm, best_hand, width, xv, true);
    //     if (bsrv < xv) {
    //         best_sol = last_solution;
    //         xv = bsrv;
    //     }
    // }

    DB(elapsed());

    DATAX(use_greedy, "greedy");
    DATAX(timer_bsmoves.total, "tm");
    DATAX(timer_bsapply.total, "tapp");
    DATAX(timer_bsactive.total, "tact");

    double time = elapsed();
    DATA(time);

    print_solution(best_sol);

	return 0;
}
