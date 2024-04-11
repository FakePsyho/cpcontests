// Author: Psyho
// Twitter: https://twitter.com/fakepsyho

// TEMPLATE

#pragma GCC optimize "Ofast,omit-frame-pointer,inline,fast-math,unroll-all-loops"
#pragma GCC target "avx,avx2,f16c,fma,sse3,ssse3,sse4.1,sse4.2,bmi,bmi2,lzcnt,popcnt"


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

template <class T, auto op, auto def> struct Segtree {
    int n;
    int log;
    int size;
    VC<T> data;

    Segtree(int n = 0) : Segtree(n, def()) {}
    Segtree(int n, T x) {VC<T> v(n, x); _init(&v[0], n);}
    Segtree(VC<T> &v) {_init(&v[0], v.S);}
    Segtree(T* v, int n) {_init(v, n);}

    void _init(T* v, int n) {
        this->n = n;
        log = 32 - __builtin_clz(n);
        size = 1 << log;
        data.resize(2 * size, def());
        REP(i, n) data[size + i] = v[i];
        for (int i = size - 1; i > 0; i--) data[i] = op(data[i<<1], data[(i<<1) + 1]);
    }

    void clear() {REP(i, 2 * size) data[i] = def();}

    void set(int p, T x) {
        data[p += size] = x;
        for (p >>= 1; p > 0; p >>= 1) data[p] = op(data[p<<1], data[(p<<1) + 1]);
    }

    INLINE T get(int p) {return data[p + size];}

    T range(int l, int r) {
        T sml = def, smr = def;
        for (l += size, r += size; l < r; l >>= 1, r >>= 1) {
            if (l & 1) sml = op(sml, data[l++]);
            if (r & 1) smr = op(data[--r], smr);
        }
        return op(sml, smr);
    }

    INLINE T all() {return data[1];}
};

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

// TODO: add version with hashing?
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

    INLINE bool add(const Candidate<ValueT, Move>& cand) {
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

const double TIME_LIMIT = 1.95;
const int MOD = 998244353;
const int N = 9;
const int N2 = 81;
const int M = 20;
const int K = 81;

const int TOTAL_M = M + 1 + (M+1)*M/2 + (M+1)*M*(M+2)/6 + (M+3)*M*(M+1)*(M+2)/24 + (M+4)*(M+3)*M*(M+1)*(M+2)/120 + (M+5)*(M+4)*(M+3)*M*(M+1)*(M+2)/720 + (LL)(M+6)*(M+5)*(M+4)*(M+3)*M*(M+1)*(M+2)/5040;
const int TOTAL_MX = TOTAL_M * 2;

LL a[N2];

int t0[M][9];
int t[TOTAL_MX][9];
char t_movecount[TOTAL_M];
int s_movestart[10];
int s_movecount[10];
char t_moves[TOTAL_MX][8];

int n_sol = 0;
int sol[N2][2] = {-1};

const int N_BUCKETS = 20;
const int N_BUCKETS_USED = 3;
VI ttr_buckets[8][N_BUCKETS][N_BUCKETS][N_BUCKETS];
VI tbl_buckets[8][N_BUCKETS][N_BUCKETS][N_BUCKETS];

void add_sol(int m, int p) {
    sol[n_sol][0] = m;
    sol[n_sol][1] = p;
    n_sol++;
    REP(dr, 3) REP(dc, 3) a[p+dr*N+dc] = (a[p+dr*N+dc] + t0[m][dr*3+dc]) % MOD;
}

LL xv = 0;

const int N_LEVELS = 49;
int lev_order[N_LEVELS];

VI lev_pos[N_LEVELS];
VI les_off[N_LEVELS];
int lev_k[N_LEVELS];

struct Move {
    int id;
};

struct State {
    int move_id;
    int sleft;
    LL av = 0;
    int v[N2];

    INLINE void apply_move(int level, Move &m, int move_id, State &s) {
        s.sleft = sleft;
        s.av = av;
        s.move_id = move_id;
        memcpy(s.v, v, sizeof(v));

        int pos = lev_order[level];
        s.sleft -= t_movecount[m.id];
        s.v[pos+ 0] = (s.v[pos+ 0] + t[m.id][0]) % MOD;
        s.v[pos+ 1] = (s.v[pos+ 1] + t[m.id][1]) % MOD;
        s.v[pos+ 2] = (s.v[pos+ 2] + t[m.id][2]) % MOD;
        s.v[pos+ 9] = (s.v[pos+ 9] + t[m.id][3]) % MOD;
        s.v[pos+10] = (s.v[pos+10] + t[m.id][4]) % MOD;
        s.v[pos+11] = (s.v[pos+11] + t[m.id][5]) % MOD;
        s.v[pos+18] = (s.v[pos+18] + t[m.id][6]) % MOD;
        s.v[pos+19] = (s.v[pos+19] + t[m.id][7]) % MOD;
        s.v[pos+20] = (s.v[pos+20] + t[m.id][8]) % MOD;

        for (int p : lev_pos[level]) s.av += s.v[p];
    }

    template<int N> LL calc_move(int level, int move) {
        LL x = av;

        if (N == 1) {
            int a = v[lev_pos[level][0]] + t[move][0];
            x += a - (a >= MOD) * MOD;
        } else {
            REP(i, lev_pos[level].S) {
                int p = lev_pos[level][i];
                int o = les_off[level][i];
                int a = v[p] + t[move][o];
                x += a - (a >= MOD) * MOD;
            }
        }

        return x;
    }
};

const int MAX_STATES = 22500;
const int MAX_TOTAL_MOVES = MAX_STATES * N_LEVELS;

const int K_DIFF = 5;

int n_move_history = 0;
MoveHistory<Move> move_history[MAX_TOTAL_MOVES];
State all_states[2][MAX_STATES];

int main(int argc, char **argv) {
    int _; cin >> _ >> _ >> _;

    REP(i, N2) cin >> a[i];
    REP(i, M) REP(j, 9) cin >> t0[i][j];

    s_movestart[0] = 0;
    s_movecount[0] = 1;
    t_movecount[0] = 0;

    int t_start = 1;
    int n_tmove = 0;

    auto sort_moves = [&](int len, int moves) {
        VPII vp; REP(i, len) vp.PB(MP(-t[TOTAL_M+i][0], i));
        sort(ALL(vp));

        s_movestart[moves] = t_start;
        s_movecount[moves] = len;
        REP(i, len) {
            REP(j, 9) t[t_start+i][j] = t[TOTAL_M+vp[i].Y][j];
            REP(j, moves) t_moves[t_start+i][j] = t_moves[TOTAL_M+vp[i].Y][j];
            t_movecount[t_start+i] = moves;
        }
        t_start += len;
    };

    n_tmove = 0;
    REP(m1, M) {
        REP(i, 9) t[TOTAL_M+n_tmove][i] = t0[m1][i];
        t_moves[TOTAL_M+n_tmove][0] = m1;
        n_tmove++;
    }
    sort_moves(n_tmove, 1);

    n_tmove = 0;
    REP(m1, M) REP(m2, m1+1) {
        REP(i, 9) t[TOTAL_M+n_tmove][i] = (t0[m1][i] + t0[m2][i]) % MOD;
        t_moves[TOTAL_M+n_tmove][0] = m1;
        t_moves[TOTAL_M+n_tmove][1] = m2;
        n_tmove++;
    }
    sort_moves(n_tmove, 2);

    n_tmove = 0;
    REP(m1, M) REP(m2, m1+1) REP(m3, m2+1) {
        REP(i, 9) t[TOTAL_M+n_tmove][i] = ((LL)t0[m1][i] + t0[m2][i] + t0[m3][i]) % MOD;
        t_moves[TOTAL_M+n_tmove][0] = m1;
        t_moves[TOTAL_M+n_tmove][1] = m2;
        t_moves[TOTAL_M+n_tmove][2] = m3;
        n_tmove++;
    }
    sort_moves(n_tmove, 3);
    
    n_tmove = 0;
    REP(m1, M) REP(m2, m1+1) REP(m3, m2+1) REP(m4, m3+1) {
        REP(i, 9) t[TOTAL_M+n_tmove][i] = ((LL)t0[m1][i] + t0[m2][i] + t0[m3][i] + t0[m4][i]) % MOD;
        t_moves[TOTAL_M+n_tmove][0] = m1;
        t_moves[TOTAL_M+n_tmove][1] = m2;
        t_moves[TOTAL_M+n_tmove][2] = m3;
        t_moves[TOTAL_M+n_tmove][3] = m4;
        n_tmove++;
    }
    sort_moves(n_tmove, 4);

    n_tmove = 0;
    REP(m1, M) REP(m2, m1+1) REP(m3, m2+1) REP(m4, m3+1) REP(m5, m4+1) {
        REP(i, 9) t[TOTAL_M+n_tmove][i] = ((LL)t0[m1][i] + t0[m2][i] + t0[m3][i] + t0[m4][i] + t0[m5][i]) % MOD;
        t_moves[TOTAL_M+n_tmove][0] = m1;
        t_moves[TOTAL_M+n_tmove][1] = m2;
        t_moves[TOTAL_M+n_tmove][2] = m3;
        t_moves[TOTAL_M+n_tmove][3] = m4;
        t_moves[TOTAL_M+n_tmove][4] = m5;
        n_tmove++;
    }
    sort_moves(n_tmove, 5);

    n_tmove = 0;
    REP(m1, M) REP(m2, m1+1) REP(m3, m2+1) REP(m4, m3+1) REP(m5, m4+1) REP(m6, m5+1){
        REP(i, 9) t[TOTAL_M+n_tmove][i] = ((LL)t0[m1][i] + t0[m2][i] + t0[m3][i] + t0[m4][i] + t0[m5][i] + t0[m6][i]) % MOD;
        t_moves[TOTAL_M+n_tmove][0] = m1;
        t_moves[TOTAL_M+n_tmove][1] = m2;
        t_moves[TOTAL_M+n_tmove][2] = m3;
        t_moves[TOTAL_M+n_tmove][3] = m4;
        t_moves[TOTAL_M+n_tmove][4] = m5;
        t_moves[TOTAL_M+n_tmove][5] = m6;
        n_tmove++;
    }
    sort_moves(n_tmove, 6);

    n_tmove = 0;
    REP(m1, M) REP(m2, m1+1) REP(m3, m2+1) REP(m4, m3+1) REP(m5, m4+1) REP(m6, m5+1) REP(m7, m6+1) {
        REP(i, 9) t[TOTAL_M+n_tmove][i] = ((LL)t0[m1][i] + t0[m2][i] + t0[m3][i] + t0[m4][i] + t0[m5][i] + t0[m6][i] + t0[m7][i]) % MOD;
        t_moves[TOTAL_M+n_tmove][0] = m1;
        t_moves[TOTAL_M+n_tmove][1] = m2;
        t_moves[TOTAL_M+n_tmove][2] = m3;
        t_moves[TOTAL_M+n_tmove][3] = m4;
        t_moves[TOTAL_M+n_tmove][4] = m5;
        t_moves[TOTAL_M+n_tmove][5] = m6;
        t_moves[TOTAL_M+n_tmove][6] = m7;
        n_tmove++;
    }
    sort_moves(n_tmove, 7);

    REP(i, s_movestart[6]+s_movecount[6]) {
        {
            int b0 = (LL)t[i][0] * N_BUCKETS / MOD;
            int b1 = (LL)t[i][1] * N_BUCKETS / MOD;
            int b2 = (LL)t[i][2] * N_BUCKETS / MOD;
            ttr_buckets[t_movecount[i]][b0][b1][b2].PB(i);
        }

        {
            int b0 = (LL)t[i][0] * N_BUCKETS / MOD;
            int b1 = (LL)t[i][3] * N_BUCKETS / MOD;
            int b2 = (LL)t[i][6] * N_BUCKETS / MOD;
            tbl_buckets[t_movecount[i]][b0][b1][b2].PB(i);
        }
    }

    FOR(i, s_movestart[7], s_movestart[7]+s_movecount[7]) {
        int b0 = (LL)t[i][0] * N_BUCKETS / MOD;
        int b1 = (LL)t[i][1] * N_BUCKETS / MOD;
        int b2 = (LL)t[i][2] * N_BUCKETS / MOD;
        ttr_buckets[t_movecount[i]][b0][b1][b2].PB(i);
    }

    REP(n, 8) REP(i, N_BUCKETS) REP(j, N_BUCKETS) REP(k, N_BUCKETS) {
        sort(ALL(ttr_buckets[n][i][j][k]));
        sort(ALL(tbl_buckets[n][i][j][k]));
    }

    DB(sizeof(move_history));
    DB(sizeof(all_states));


    VPII vp; REP(r, 7) REP(c, 7) vp.PB(MP(min(r,c) * 100 + r * 10 + c, r*N+c));
    sort(ALL(vp));
    REP(i, vp.S) lev_order[i] = vp[i].Y;

    VVI vlev_pos(81); REP(r, 9) REP(c, 9) vlev_pos[min(r,6)*N+min(c,6)].PB(r*N+c);
    REP(i, N_LEVELS) {
        lev_pos[i] = vlev_pos[lev_order[i]];
        for (int x : lev_pos[i]) {
            int diff = x - lev_order[i];
            les_off[i].PB((diff/N)*3 + (diff%N));
        }
    }

    VD lev_w(N_LEVELS); 
    REP(i, N_LEVELS) {
        int p = lev_order[i];
        double w = 1.0;
        w += (p%9 == 6) * 1.50;
        w += (p/9 == 6) * 1.50;
        w += (p/9 == 6 && p%9 == 6) * 1.50;
        w = min(w, 4.0);
        lev_w[i] = w;
    }
    FOR(i, 1, N_LEVELS) lev_w[i] += lev_w[i-1];
    REP(i, N_LEVELS) lev_k[i] = (int)(lev_w[i]*N2/lev_w[N_LEVELS-1]+.5);

    DB(VI(lev_k, lev_k+N_LEVELS));

    VD groupw(K_DIFF*2+1);
    REP(i, K_DIFF*2+1) {
        int diff = i - K_DIFF;
        groupw[i] = pow(K_DIFF - abs(diff) + 1, 3.0);
    }
    double totalw = accumulate(ALL(groupw), 0.0);
    VI group_sizes(K_DIFF*2+1); REP(i, K_DIFF*2+1) group_sizes[i] = MAX_STATES * groupw[i] / totalw;
    VI group_offset(K_DIFF*2+1); FOR(i, 1, K_DIFF*2+1) group_offset[i] = group_offset[i-1] + group_sizes[i-1];

    DB(group_sizes);
    DB(group_offset);

    State& s0 = all_states[0][0];
    s0.sleft = K;
    s0.av = 0;
    s0.move_id = -1;
    REP(i, N2) s0.v[i] = a[i];

    double timer1 = 0;
    double timer2 = 0;
    double timer3 = 0;
    double timer4 = 0;
    double timer5 = 0;

    int n_states = 1;
    VC<CandidatesGroup<LL, Move>> cg(2*K_DIFF+1);
    REP(i, 2*K_DIFF+1) cg[i] = CandidatesGroup<LL, Move>(group_sizes[i]);
    DB(elapsed());

    VI adds_type(3);

    REP(level, N_LEVELS) {
        timer4 -= elapsed();
        State *states = all_states[level&1];
        State *nstates = all_states[1-(level&1)];

        int n_moves_start = level * MAX_STATES;

        int pos = lev_order[level];

        int p0 = lev_pos[level][0];
        int o0 = les_off[level][0];


        bool right_edge = pos%9 == 6;
        bool bottom_edge = pos/9 == 6;

        REP(i, 2*K_DIFF+1) cg[i].clear();

        if (right_edge && bottom_edge) REP(i, 2*K_DIFF+1) group_sizes[i] = 1, cg[i] = CandidatesGroup<LL, Move>(1);

        timer4 += elapsed();
        int adds = 0;
        auto generate_new_states = [&](int state_id, const int n_stamps, int cut_off) {
            State &s = states[state_id];
            int ns_group = (N2 - s.sleft) - lev_k[level] + n_stamps + K_DIFF;
            if (s.sleft < n_stamps || ns_group < 0 || ns_group > K_DIFF*2) return;

            int val = s.v[p0];
            int l = 0, r = s_movecount[n_stamps];
            while (l < r) {
                int m = (l + r) / 2;
                if (val + t[s_movestart[n_stamps] + m][o0] >= MOD) {
                    l = m + 1;
                } else {
                    r = m;
                }
            }
            int tpos = l;
            cut_off += (tpos == s_movecount[n_stamps]) * MOD;
            tpos %= s_movecount[n_stamps];

            l = 0, r = s_movecount[n_stamps];
            while (l < r) {
                int m = (l + r) / 2;
                int pos = (tpos + m) % s_movecount[n_stamps];
                if ((val + t[s_movestart[n_stamps] + pos][o0]) % MOD >= cut_off) {
                    l = m + 1;
                } else {
                    r = m;
                }
            }
            int tlen = l;
            while (tlen--) {
                int move = s_movestart[n_stamps] + tpos;
                if (++tpos == s_movecount[n_stamps]) tpos = 0, cut_off += MOD;
                int x = val + t[move][0];
                LL av = s.av + x - (x >= MOD) * MOD;
                if (!cg[ns_group].add({state_id, av, {move}})) break;
            }
        };

        auto generate_new_multi_states = [&](int state_id, const int offset1, const int offset2, const int max_moves, auto &buckets) {
            State &s = states[state_id];
            int b0 = (LL)(MOD - s.v[p0+      0] - 1) * N_BUCKETS / MOD;
            int b1 = (LL)(MOD - s.v[p0+offset1] - 1) * N_BUCKETS / MOD;
            int b2 = (LL)(MOD - s.v[p0+offset2] - 1) * N_BUCKETS / MOD;
            static VVI all_buckets(N_BUCKETS_USED*N_BUCKETS_USED*N_BUCKETS_USED, VI(3));
            REP(b0_add, N_BUCKETS_USED) REP(b1_add, N_BUCKETS_USED) REP(b2_add, N_BUCKETS_USED) {
                int set = b0_add + N_BUCKETS_USED * (b1_add + N_BUCKETS_USED * b2_add);
                all_buckets[set][0] = (b0 + N_BUCKETS - N_BUCKETS_USED + 1 + b0_add) % N_BUCKETS;
                all_buckets[set][1] = (b1 + N_BUCKETS - N_BUCKETS_USED + 1 + b1_add) % N_BUCKETS;
                all_buckets[set][2] = (b2 + N_BUCKETS - N_BUCKETS_USED + 1 + b2_add) % N_BUCKETS;
            }
            int ns_group_start = (N2 - s.sleft) - lev_k[level] + K_DIFF;
            int n_stamps_start = max(0, -ns_group_start);
            int n_stamps_end = min(max_moves, K_DIFF*2 - ns_group_start + 1);
            n_stamps_end = min(n_stamps_end, s.sleft+1);
            FOR(n_stamps, n_stamps_start, n_stamps_end) {
                int ns_group = (N2 - s.sleft) - lev_k[level] + n_stamps + K_DIFF;
                REP(sets, N_BUCKETS_USED*N_BUCKETS_USED*N_BUCKETS_USED) {
                    for (int move : buckets[n_stamps][all_buckets[sets][0]][all_buckets[sets][1]][all_buckets[sets][2]]) {
                        LL av = s.calc_move<2>(level, move);
                        cg[ns_group].add({state_id, av, {move}});
                    }
                }
            }
        };

        if (right_edge && bottom_edge) {
            timer3 -= elapsed();
            REP(i, n_states) {
                // if ((i & 15) == 0 && elapsed() > TIME_LIMIT) break;
                generate_new_multi_states(i, 1,  2, 8, ttr_buckets);
            }
            timer3 += elapsed();
            adds_type[2] += adds;
        } else if (right_edge) {
            timer2 -= elapsed();
            REP(i, n_states) generate_new_multi_states(i, 1,  2, 6, ttr_buckets);
            timer2 += elapsed();
            adds_type[1] += adds;
        } else if (bottom_edge) {
            timer2 -= elapsed();
            REP(i, n_states) generate_new_multi_states(i, 9, 18, 6, tbl_buckets);
            timer2 += elapsed();
            adds_type[1] += adds;
        } else {
            timer1 -= elapsed();
            REP(i, n_states) generate_new_states(i, 0, (LL)MOD * 85 / 100);
            REP(i, n_states) generate_new_states(i, 1, (LL)MOD * 90 / 100);
            REP(i, n_states) generate_new_states(i, 2, (LL)MOD * 92 / 100);
            REP(i, n_states) generate_new_states(i, 3, (LL)MOD * 95 / 100);
            timer1 += elapsed();
            adds_type[0] += adds;
        }
        LL bv = 0;


        timer5 -= elapsed();
        n_states = 0;
        for (auto &cgroup : cg) {
            sort(cgroup.candidates.begin(), cgroup.candidates.begin() + cgroup.count, [&](const Candidate<LL, Move> &a, const Candidate<LL, Move> &b) {return a.value < b.value;});
            REP(i, cgroup.count) {
                Candidate<LL, Move> &ns = cgroup.candidates[i];
                MoveHistory<Move> &h = move_history[n_move_history];
                h.previous_id = states[ns.state_id].move_id;
                h.move = ns.move;
                // nstates[n_states] = states[ns.state_id].apply_move(level, h.move, n_move_history);
                states[ns.state_id].apply_move(level, h.move, n_move_history, nstates[n_states]);
                bv = max(bv, nstates[n_states].av);
                n_move_history++;
                n_states++;
            }
        }
        timer5 += elapsed();
        // DB(level, bv);
    }

    DB(adds_type);
    DB(elapsed());

    int best_state = 0;
    REP(i, n_states) if (all_states[N_LEVELS&1][i].av > all_states[N_LEVELS&1][best_state].av) best_state = i;
    State& bs = all_states[N_LEVELS&1][best_state];

    DB(bs.av, bs.sleft);

    REP(i, K) sol[i][0] = -1;
    int move_id = bs.move_id;
    for (int level = N_LEVELS-1; level >= 0; level--) {
        MoveHistory<Move> &m = move_history[move_id];
        REP(i, t_movecount[m.move.id]) add_sol(t_moves[m.move.id][i], lev_order[level]);
        move_id = m.previous_id;
    }
    reverse(sol, sol+n_sol);
    assert(move_id == -1);

    if (timer1) DB(timer1);
    if (timer2) DB(timer2);
    if (timer3) DB(timer3);
    if (timer4) DB(timer4);
    if (timer5) DB(timer5);

    cerr << "[DATA] time = " << elapsed() << endl;
    cerr << "[DATA] sleft = " << bs.sleft << endl;

    int non_zero = 0;
    REP(i, K) non_zero += sol[i][0] != -1;
    cout << non_zero << endl;
    REP(i, K) if (sol[i][0] != -1) cout << (int)sol[i][0] << ' ' << (int)sol[i][1]/9 << " " << (int)sol[i][1]%9 << endl;

    return 0;
}
