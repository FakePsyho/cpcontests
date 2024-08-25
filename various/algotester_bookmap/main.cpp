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
#define PID         pair<int, double>

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
 
static RNG rng(time(0));

// SOLUTION

int sign(double x) {return x < 0 ? -1 : x > 0 ? 1 : 0;}

const double TIME_SCALE = 1.0;

#ifdef LOCAL
const double TIME_LIMIT = 3.0 * TIME_SCALE;
#elif VM
const double TIME_LIMIT = 3.0 * TIME_SCALE;
#else
const double TIME_LIMIT = 3.8;
#endif


const int MAX_T     = 104;
const int MAX_C     = 1000000;
const int MAX_N     = 104;
const int MAX_L     = 100000;
const int MAX_PRICE = 1000000000;

int N, T, C;
int P0[MAX_N];
int L[MAX_N];
double K[MAX_N][MAX_T];
double E[MAX_N][MAX_T];
double A[MAX_N][MAX_N][MAX_T];
double B[MAX_T];

double L_ratio;

struct State {
    int amount[MAX_N][MAX_T];
    double P[MAX_N][MAX_T];
    double P_delta[MAX_N][MAX_T];
    int total_amount[MAX_T];

    void clear() {
        REP(i, N) REP(t, T) amount[i][t] = 0;
        REP(i, N) total_amount[i] = 0;
    }

    void init_last(VI &allowed) {
        REP(i, N) {
            amount[i][T-1] = allowed[i];
            total_amount[i] = allowed[i];
        }
    }

    void init_random(VI &allowed) {
        REP(i, N) {
            int left = allowed[i];
            total_amount[i] = left;
            while (left) {
                int v = max(1, left / T / 10);
                amount[i][rng.next(T)] += v;
                left -= v;
            }
        }
    }

    void output() {
        REP(t, T) {
            REP(i, N) cout << amount[i][t] << ' ';
            cout << endl;
        }
    }

    void copy(State &s, int min_n = -1, int max_n = -1) {
        if (min_n == -1) {
            min_n = 0;
            max_n = N;
        } else if (max_n == -1) {
            max_n = min_n + 1;
        }
        assert(min_n >= 0 && max_n <= N && min_n < max_n);

        FOR(i, min_n, max_n) REP(t, T+1) {
            amount[i][t] = s.amount[i][t];
            P[i][t] = s.P[i][t];
            P_delta[i][t] = s.P_delta[i][t];
        }
        FOR(i, min_n, max_n) total_amount[i] = s.total_amount[i];
    }

    void verify() {
        int total_sum = 0;
        REP(i, N) {
            int sum = 0;
            REP(t, T) {
                assert(amount[i][t] >= 0);
                sum += amount[i][t];
            }
            assert(sum <= L[i]);
            total_sum += sum;
        }
        assert(total_sum == C);
    }
};

template<int PRECISION> double sim(State &s, int min_n = 0, int max_n = 0, int min_t = 0) {
    if (max_n == 0) max_n = N;

    int L_left[MAX_N];
    REP(i, N) s.P[i][0] = P0[i];
    REP(i, N) L_left[i] = L[i];

    double e = exp(1.0);
    REP(t, min_t) FOR(i, min_n, max_n) L_left[i] -= s.amount[i][t];

    FOR(i, min_n, max_n) {
        FOR(t, min_t, T) {  
            double delta = 0;
            if (s.amount[i][t]) delta += s.P[i][t] * (1.0 - 1.0 / (expf(K[i][t] * s.amount[i][t] / (L_left[i] + 1))));
            delta += s.P[i][t] * E[i][t];
            if (PRECISION >= 1) REP(j, i) delta += (A[i][j][t] * s.P_delta[j][t]) / logf(max(e, abs(s.P[i][t] - s.P[j][t])));
            if (PRECISION >= 2) REP(j, t) delta += B[j] * s.amount[i][t - j - 1];
            s.P_delta[i][t] = delta;
            s.P[i][t+1] = s.P[i][t] + delta;
            L_left[i] -= s.amount[i][t];
        }
    }

    double sum_delta = 0;
    REP(i, max_n) sum_delta += abs(s.P[i][T] - P0[i]);

    return sum_delta;
}

double quick_sim_raw(State &s, int i) {
    double left = L[i] + 1;
    double P = P0[i];

    REP(t, T) {
        double delta = 0;
        if (s.amount[i][t]) delta += P * (1.0 - 1.0 / (expf(K[i][t] * s.amount[i][t] / left)));
        delta += P * E[i][t];
        P += delta;
        left -= s.amount[i][t];
    } 
    return P - P0[i];
}

double quick_sim(State &s, int i) {
    return abs(quick_sim_raw(s, i));
}

double exact_sim(State &s, int min_n = 0, int max_n = 0, int min_t = 0) {
    if (max_n == 0) max_n = N;

    int L_left[MAX_N];
    REP(i, N) s.P[i][0] = P0[i];
    REP(i, N) L_left[i] = L[i];

    double e = exp(1.0);
    REP(t, min_t) FOR(i, min_n, max_n) L_left[i] -= s.amount[i][t];

    FOR(i, min_n, max_n) {
        FOR(t, min_t, T) {
            double delta = 0;
            delta += s.P[i][t] * (1.0 - 1.0 / (exp(K[i][t] * s.amount[i][t] / (L_left[i] + 1))));
            delta += s.P[i][t] * E[i][t];
            REP(j, i) delta += (A[i][j][t] * s.P_delta[j][t]) / log(max(e, abs(s.P[i][t] - s.P[j][t])));
            REP(j, t) delta += B[j] * s.amount[i][t - j - 1];
            s.P_delta[i][t] = delta;
            s.P[i][t+1] = s.P[i][t] + delta;
            L_left[i] -= s.amount[i][t];
        }

    }

    double sum_delta = 0;
    REP(i, max_n) sum_delta += abs(s.P[i][T] - P0[i]);

    return sum_delta;
}

int final_score(State &s) {
    double score = exact_sim(s);

    s.verify();

    double p0_sum = 0; REP(i, N) p0_sum += P0[i];
    return (int)max(0.0, (((10 * p0_sum - score) / (10 * p0_sum)) * 1e7));
}

State bs;
double bxv = 1e300;

int last_mstep = 0;
void quick_sa_single(State &s, double time_limit, int row = -1) {
    double quick_sa_start = elapsed();
    int start = row == -1 ? 0 : row;
    int end = row == -1 ? N : row + 1;

    int gstep = 0;
    FOR(i, start, end) {
        double step_start = elapsed();
        double step_end = quick_sa_start + time_limit * (i + 1) / N;
        step_end = max(step_end, step_start + time_limit / N * 0.5);
        if (row >= 0) step_end = quick_sa_start + time_limit;
        int mstep = 0;
        double time_passed = 0;

        if (s.total_amount[i] == 0) continue;
        double bv = quick_sim(s, i);

        while (true) {
            int t0 = rng.next(T);
            if (s.amount[i][t0] == 0) continue;
            int t1 = rng.next(T-1);
            t1 += t1 >= t0;

            mstep++;
            if ((mstep & 31) == 0) {
                time_passed = (elapsed() - step_start) / (step_end - step_start);
                if (time_passed > 1.0) break;
            }

            int v = rng.next(1, s.amount[i][t0] + 1);

            s.amount[i][t0] -= v;
            s.amount[i][t1] += v;
            double av = quick_sim(s, i);
            if (av <= bv) {
                bv = av;
            } else {
                s.amount[i][t0] += v;
                s.amount[i][t1] -= v;
            }
        }
        gstep += mstep;

    }
    last_mstep = gstep;
}

int global_updates = 0;
int global_steps = 0;

void quick_sa_global(State &s, double time_limit, bool early_exit = true) {
    double quick_sa_start = elapsed();
    int mstep = 0;
    int last_acc = 0;
    int last_update = 0;

    VD bv(N); REP(i, N) bv[i] = quick_sim(s, i);

    const int N_REPS = 25;
    while (true) {
        int p0 = rng.next(N);
        int t0 = rng.next(T);
        if (s.amount[p0][t0] == 0) continue;
        int p1 = rng.next(N);
        int t1 = rng.next(T);

        if (p0 == t0 && p1 == t1) continue;
        int v = rng.next(1, s.amount[p0][t0] + 1);

        if (p0 != p1 && s.total_amount[p1] + v > L[p1]) continue;


        mstep++;
        if ((mstep & 31) == 0) {
            double time_passed = (elapsed() - quick_sa_start) / time_limit;
            if (early_exit && mstep > last_acc + 5 * T * N) break;
            if (time_passed > 1.0) break;
        }

        s.amount[p0][t0] -= v;
        s.amount[p1][t1] += v;

        if (p0 == p1) {
            double av = quick_sim(s, p0);
            if (av < bv[p0]) {
                bv[p0] = av;
                last_acc = mstep;
                s.amount[p0][t0] += v;
                s.amount[p1][t1] -= v;
                int bestv = v;
                REP(i, N_REPS) {
                    v = rng.next(1, s.amount[p0][t0] + 1);
                    s.amount[p0][t0] -= v;
                    s.amount[p1][t1] += v;
                    av = quick_sim(s, p0);
                    s.amount[p0][t0] += v;
                    s.amount[p1][t1] -= v;
                    if (av < bv[p0]) {
                        bestv = v;
                        bv[p0] = av;
                    }
                }
                s.amount[p0][t0] -= bestv;
                s.amount[p1][t1] += bestv;
            } else {
                s.amount[p0][t0] += v;
                s.amount[p1][t1] -= v;
            }
        } else {
            double av0 = quick_sim(s, p0);
            double av1 = quick_sim(s, p1);
            if (av0 + av1 < bv[p0] + bv[p1]) {
                bv[p0] = av0;
                bv[p1] = av1;
                last_acc = mstep;
                s.amount[p1][t1] -= v;
                s.amount[p0][t0] += v;
                int bestv = v;
                REP(i, N_REPS) {
                    v = rng.next(1, s.amount[p0][t0] + 1);
                    if (s.total_amount[p1] + v > L[p1]) continue;
                    s.amount[p0][t0] -= v;
                    s.amount[p1][t1] += v;
                    av0 = quick_sim(s, p0);
                    av1 = quick_sim(s, p1);
                    s.amount[p0][t0] += v;
                    s.amount[p1][t1] -= v;
                    if (av0 + av1 < bv[p0] + bv[p1]) {
                        bestv = v;
                        bv[p0] = av0;
                        bv[p1] = av1;
                    }

                }
                s.amount[p0][t0] -= bestv;
                s.amount[p1][t1] += bestv;
                s.total_amount[p0] -= bestv;
                s.total_amount[p1] += bestv;
            } else {
                s.amount[p0][t0] += v;
                s.amount[p1][t1] -= v;
            }
        }
    }

    int zeros = 0; REP(i, N) REP(t, T) if (s.amount[i][t] == 0) zeros++;
    global_steps += mstep;
}

double quick_est_single(int i, int c, double time_limit, int bestt = -1) {
    static State s;
    REP(t, T) s.amount[i][t] = 0;
    if (bestt == -1)
        REP(t, T) s.amount[i][t] = c / T + (t < c % T);
    else
        s.amount[i][bestt] = c;
    s.total_amount[i] = c;
    quick_sa_single(s, time_limit, i);
    return quick_sim(s, i);
}

VI sa_dp(VC<VC<PID>> &dp) {
    int C_left = C;
    VI state(N);

    while (C_left) {
        int non_zero = 0;
        REP(i, N) if (state[i] < L[i]) non_zero++;
        int x = (C_left + non_zero - 1) / non_zero;
        REP(i, N) {
            int v = min(x, L[i] - state[i]);
            v = min(v, C_left);
            state[i] += v;
            C_left -= v;
        }
    }

    auto get_state_value = [&](int i, int c) -> double{
        int p = 0;
        if (c >= dp[i].back().X) {
            return dp[i].back().Y;
        } else {
            while (dp[i][p+1].X < c) p++;
            double ratio = 1.0 * (c - dp[i][p].X) / (dp[i][p+1].X - dp[i][p].X);
            return dp[i][p].Y + ratio * (dp[i][p+1].Y - dp[i][p].Y);
        }
    };

    double bv = 0; REP(i, N) bv += get_state_value(i, state[i]);
    double xv = bv;

    double temp0 = bv;
    double tempn = 0.01;
    double temp = temp0;

    VI bstate = state;

    const int SA_DP_STEPS = 1000000;
    REP(step, SA_DP_STEPS) {
        int a = rng.next(N);
        if (state[a] == 0) continue;
        int b = rng.next(N-1);
        b += b >= a;
        int v = rng.next(1, state[a] + 1);
        if (state[b] + v > L[b]) continue;

        double av = bv;
        av -= get_state_value(a, state[a]);
        av -= get_state_value(b, state[b]);
        state[a] -= v;
        state[b] += v;
        av += get_state_value(a, state[a]);
        av += get_state_value(b, state[b]);

        if ((step & 127) == 0) {
            temp = temp0 * pow(tempn / temp0, 1.0 * step / SA_DP_STEPS);
        }

        if (av < bv || rng.next_double() < exp((bv - av) / temp)) {
            bv = av;
            if (av < xv) {
                xv = av;
                bstate = state;
            }
        } else {
            state[a] += v;
            state[b] -= v;
        }
    }

    return bstate;
}


int main(int argc, char **argv) {
    cin.sync_with_stdio(0);
    cin.tie(0);


    cin >> N >> T >> C;

    REP(i, N) cin >> P0[i] >> L[i];
    REP(i, N) REP(t, T) cin >> K[i][t];
    REP(i, N) REP(t, T) cin >> E[i][t];
    REP(i, N) REP(j, i) REP(t, T) cin >> A[i][j][t];
    REP(t, T) cin >> B[t];

    cerr << "[DATA] N = " << N << endl;
    cerr << "[DATA] T = " << T << endl;
    cerr << "[DATA] C = " << C << endl;
    int L_sum = 0; REP(i, N) L_sum += L[i];
    L_ratio = 1.0 * L_sum / C;
    cerr << "[DATA] L_ratio = " << L_ratio << endl;

    bool EXTRACT = false;
    int EXTRACT_VALUE = N;
    VC<char> mem;
    if (EXTRACT) {
        DB(EXTRACT_VALUE);
        assert(EXTRACT_VALUE <= 200);
        int MEM_VALUE = 10 + (EXTRACT_VALUE);
        int ARRAY_SIZE = (1<<20)*(2*MEM_VALUE);
        mem = VC<char>(ARRAY_SIZE);
        LL x = 0; for (int i = 0; i < ARRAY_SIZE; i += 100) x += mem[i];
        DB(x);
    }

    if (L_ratio > 8.0) {

        VC<VC<PID>> dp(N);
        const double TL = 0.0005 * TIME_LIMIT;
        const int N_STEPS = 10;
        REP(i, N) {
            int l = 0, r = min(C, L[i]);
            REP(steps, N_STEPS) {
                int m = (l + r) / 2;
                double v = quick_est_single(i, m, TL);
                dp[i].PB(MP(m, v));
                if (v < 1e6) {
                    r = m;
                } else {
                    l = m + 1;
                }
            }
            dp[i].PB(MP(0, quick_est_single(i, 0, 0.0)));
            FOR(step, 1, 3) {
                int c = l * step / 3;
                dp[i].PB(MP(c, quick_est_single(i, c, TL)));
            }
            if (l + 5 <= L[i])
                dp[i].PB(MP(l+5, quick_est_single(i, l+5, TL)));
            sort(ALL(dp[i]));
        }
        DB(elapsed());
        VI dp_res = sa_dp(dp);
        DB(elapsed());


        State dps;
        dps.clear();
        dps.init_random(dp_res);
        quick_sa_single(dps, TIME_LIMIT * 0.05);
        DB(final_score(dps));
        quick_sa_global(dps, TIME_LIMIT * 0.10, false);
        DB(final_score(dps));
        DB(elapsed());
        bs.copy(dps);
    } else {

        VI allowed(N);
        int C_left = C;

        for (int i = N-1; i >= 0; i--) {
            int amount = (int)ceil(1.0 * C * L[i] / L_sum);
            amount = min(amount, C_left);
            allowed[i] = amount;
            C_left -= amount;
        }

        int init_count = L_ratio < 10.0 ? 25 : 4;
        double init_single_length = L_ratio < 10.0 ? 0 : TIME_LIMIT * 0.025;
        double init_global_length = L_ratio < 10.0 ? TIME_LIMIT * 0.02 : TIME_LIMIT * 0.175;

        int init_step = 0;

        double est_time = init_count * (init_single_length + init_global_length);
        double init_start = elapsed();
        while (true) {
            if (init_step >= init_count && elapsed() + init_single_length + init_global_length > init_start + est_time) break;
            State s;
            s.clear();
            if (L_ratio < 10 && rng.next(2)) s.init_last(allowed); else s.init_random(allowed);
            if (init_single_length) quick_sa_single(s, init_single_length);
            if (init_global_length) quick_sa_global(s, init_global_length);
            double v0 = 0; REP(i, N) v0 += quick_sim(s, i);
            double v2 = exact_sim(s);
            DB(v0, v2, v2-v0, final_score(s));
            if (v2 < bxv) {
                bxv = v2;
                bs.copy(s);
            }
            init_step++;
        }
        cerr << "[DATA] init_step = " << init_step << endl;
        cerr << "[DATA] g_upd = " << global_updates << endl;
        cerr << "[DATA] g_step = " << global_steps << endl;
    }
    DB(final_score(bs));


    double sa_start = elapsed();

    VD timestep_weights(N);
    REP(i, N) timestep_weights[i] = pow(i+1, .5);
    double total_weight = accumulate(ALL(timestep_weights), 0.0);
    REP(i, N) timestep_weights[i] /= total_weight;
    FOR(i, 1, N) timestep_weights[i] += timestep_weights[i-1];
    timestep_weights.insert(timestep_weights.begin(), 0);

    VD timesteps(N+1);
    REP(i, N+1) timesteps[i] = sa_start + timestep_weights[i] * (TIME_LIMIT - sa_start);

    auto timestep = [&](int i, double p) {return timesteps[i] * (1 - p) + timesteps[i+1] * p;};

    int step = 0;
    int acc = 0;
    REP(i, N) {
        static double prev_bv = 0;

        State ns, s;
        s.copy(bs);

        double sa_start = elapsed();
        double layers_2_time_limit;
        double layers_4_time_limit;
        double layers_6_time_limit;
        double layers_8_time_limit;        
        int max_layers;
        double full_precision_time_limit = timestep(i, .5);
        if (L_ratio < 4.0) {
            layers_2_time_limit = timestep(i, .1);
            layers_4_time_limit = timestep(i, .25);
            layers_6_time_limit = timestep(i, .5);
            layers_8_time_limit = timestep(i, .75);
            max_layers = 8;
        } else {
            layers_2_time_limit = timestep(i, .2);
            layers_4_time_limit = timestep(i, .5);
            layers_6_time_limit = timestep(i, .8);
            layers_8_time_limit = timestep(i, 1.0);
            max_layers = 6;
        }
        int full_precision = 1;
        double hc_time_limit = timestep(i, 1.0);
        int start_step = step;
        double bv = sim<1>(s, max(0, i-max_layers), i+1);
        ns.copy(s, i);
        
        int layers = 1;
        int mstep = 0;
        while (true) {
            mstep++;
            if ((mstep & 127) == 0) {
                double elapsed_time = elapsed();
                if (elapsed_time > hc_time_limit) break;
                if (elapsed_time > full_precision_time_limit && full_precision == 1) full_precision = 2, bv = sim<2>(s, max(0, i-max_layers+1), i+1);
                if (elapsed_time > layers_2_time_limit) layers = 2;
                if (elapsed_time > layers_4_time_limit) layers = 4;
                if (elapsed_time > layers_6_time_limit) layers = 6;
                if (elapsed_time > layers_8_time_limit) layers = 8;
            }

            int p0 = rng.next(max(0, i-layers+1), i+1);
            int t0 = rng.next(T);
            if (s.amount[p0][t0] == 0) continue;
            int p1 = rng.next(max(0, i-layers+1), i+1);
            int t1 = rng.next(T);
            if (p0 == p1 && t0 == t1) continue;
            int v = 1; //rng.next(1, s.amount[p0][t0] + 1);
            if (p0 != p1 && s.total_amount[p1] + v > L[p1]) continue;
            s.amount[p0][t0] -= v;
            s.amount[p1][t1] += v;
            s.total_amount[p0] -= v;
            s.total_amount[p1] += v;

            int mint = min(t0, t1);

            double av = full_precision == 1 ? sim<1>(s, min(p0, p1), i+1, mint) : sim<2>(s, min(p0, p1), i+1, mint);
            step++;

            if (av < bv) {
                bv = av;
                ns.copy(s, max(0, i-layers+1), i+1);
                acc++;
                int bestv = v;
                s.amount[p0][t0] += v;
                s.amount[p1][t1] -= v;
                s.total_amount[p0] += v;
                s.total_amount[p1] -= v;
                REP(tries, 10) {
                    s.copy(ns, max(0, min(p0, p1)), i+1);
                    s.amount[p0][t0] += bestv;
                    s.amount[p1][t1] -= bestv;
                    s.total_amount[p0] += bestv;
                    s.total_amount[p1] -= bestv;
                    v = tries ? rng.next(1, s.amount[p0][t0] + 1) : s.amount[p0][t0];
                    if (p0 != p1 && s.total_amount[p1] + v > L[p1]) continue;
                    s.amount[p0][t0] -= v;
                    s.amount[p1][t1] += v;
                    s.total_amount[p0] -= v;
                    s.total_amount[p1] += v;
                    av = full_precision == 1 ? sim<1>(s, min(p0, p1), i+1, mint) : sim<2>(s, min(p0, p1), i+1, mint);
                    if (av < bv) {
                        bv = av;
                        bestv = v;
                        ns.copy(s, max(0, min(p0, p1)), i+1);
                    }
                }
            } else {
                s.copy(ns, max(0, i-layers+1), i+1);
            }
        }

        int layer_step = step - start_step;
        double inc = bv - prev_bv;
        prev_bv = bv;
        double diff = s.P[i][T] - P0[i];
        DB(i, inc, bv, L[i], s.total_amount[i], diff, layer_step);

        bs.copy(ns, 0, i+1);
        exact_sim(bs, max(0, i-max_layers+1), i+1);
    }


    DB(step);
    cerr << "[DATA] step = " << step << endl;
    cerr << "[DATA] acc = " << acc << endl;
    cerr << "[DATA] time = " << elapsed() << endl;
    
    int score = final_score(bs);
    cerr << "Score = " << score << endl;

    bs.output();

	return 0;
}
