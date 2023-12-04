// Author: Psyho
// Twitter: https://twitter.com/fakepsyho

// TEMPLATE

#pragma GCC optimize "Ofast,omit-frame-pointer,inline,fast-math,unroll-all-loops,tree-loop-vectorize,tree-slp-vectorize"
#pragma GCC target("avx,avx2,f16c,fma,sse3,ssse3,sse4.1,sse4.2,bmi,bmi2,lzcnt,popcnt")

#include <bits/stdc++.h>
#include <sys/time.h>
 
#ifndef LOCAL
#include <windows.h>
#endif

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
 
double get_time() {
#ifdef LOCAL
    unsigned int lo, hi;
    asm volatile ( "rdtsc\n" : "=a" (lo), "=d" (hi) );
    return (double)((ULL)hi << 32 | lo) / 3.6e9;
#else
	HANDLE hProcess;
	hProcess = GetCurrentProcess();
	FILETIME CreationTime, ExitTime, KernelTime, UserTime;
	GetProcessTimes(hProcess, &CreationTime, &ExitTime, &KernelTime, &UserTime);
	uint64_t kernel = ((uint64_t)KernelTime.dwHighDateTime << 32) | KernelTime.dwLowDateTime;
	uint64_t user = ((uint64_t) UserTime.dwHighDateTime << 32) | UserTime.dwLowDateTime;
	return (kernel + user) / 1e7;
#endif
}

double start_time = get_time();
double elapsed() {return get_time() - start_time;}
 
struct RNG {
	unsigned long x=123456789, y=362436069, z=521288629;
    RNG() { }
	unsigned int rand() { unsigned long t; x ^= x << 16; x ^= x >> 5; x ^= x << 1; t = x; x = y; y = z; z = t ^ x ^ y; return z;}
    INLINE int next() {return rand(); }
    INLINE int next(int x) {return rand() % x; }
    INLINE int next(int a, int b) {return a + (rand() % (b - a)); }
    INLINE double next_double() {return (rand() + 0.5) * (1.0 / 4294967296.0); }
    INLINE double next_double(double a, double b) {return a + next_double() * (b - a); }
};
 
static RNG rng;

// SOLUTION

const int MAX_N = 100;
const int MAX_K = 10;
const int MAX_T = 1000;
const int MAX_R = 10;
const int MAX_J = 5000;

int N, K, T, R, J;
float sinr[MAX_T][MAX_K][MAX_R][MAX_N];
float d[MAX_K][MAX_R][MAX_N][MAX_N];

struct Frame {
    float tbs;
    short id;
    short user;
    short start;
    short len;
};

Frame frames[MAX_J];
int frame_id[MAX_T][MAX_N];

VVI good_users;
VVI user_ok;
VPII frame_order;
int cur_frame_order = 0;

void add_frame() {
    int f = frame_order[cur_frame_order++].Y;
    FOR(i, frames[f].start, frames[f].start + frames[f].len) {
        user_ok[i][frames[f].user] = 1;
        good_users[i].PB(frames[f].user);
        frame_id[i][frames[f].user] = f;
    }
}

const double W = 192;

double g_time_passed = 0.0;

int global_dup = 0;

struct State {
    double total_power;
    float power[MAX_T][MAX_K][MAX_R][MAX_N];
    float power_usage[MAX_T][MAX_K];
    float power_usage_single[MAX_T][MAX_K][MAX_R];

    short total_r[MAX_T][MAX_K][MAX_N];
    short total_kr[MAX_T][MAX_N];

    short power_1dlist_len[MAX_T][MAX_K][MAX_R];
    short power_1dlist[MAX_T][MAX_K][MAX_R][MAX_N];
    short power_1dlist_pos[MAX_T][MAX_K][MAX_R][MAX_N];
    short power_2dlist_len[MAX_T][MAX_R];
    PSS   power_2dlist[MAX_T][MAX_R][MAX_N*MAX_K];
    char power_2dlist_pos[MAX_T][MAX_R][MAX_N][MAX_K];

    float old_user_data[MAX_T][MAX_N];
    float user_data[MAX_T][MAX_N];
    float frame_data[MAX_J+1];

    double current_eval;

    char updated_group[MAX_T][MAX_R];
    float cur_data[MAX_T][MAX_K][MAX_R][MAX_N];

    void init() {
        memset(power, 0, T*MAX_K*MAX_R*MAX_N*sizeof(float));
        memset(power_usage, 0, T*MAX_K*sizeof(float));
        memset(power_usage_single, 0, T*MAX_K*MAX_R*sizeof(float));
        memset(total_r, 0, T*MAX_K*MAX_N*sizeof(short));
        memset(total_kr, 0, T*MAX_N*sizeof(short));
        memset(power_1dlist_len, 0, T*MAX_K*MAX_R*sizeof(short));
        memset(power_1dlist_pos, -1, T*MAX_K*MAX_R*MAX_N*sizeof(short));
        memset(power_2dlist_len, 0, T*MAX_R*sizeof(short));
        memset(power_2dlist_pos, -1, T*MAX_R*MAX_N*MAX_K*sizeof(short));
        memset(old_user_data, 0, T*MAX_N*sizeof(float));
        memset(user_data, 0, T*MAX_N*sizeof(float));
        memset(frame_data, 0, (MAX_J+1)*sizeof(float));
        memset(updated_group, 0, T*MAX_R*sizeof(bool));
        memset(cur_data, 0, T*MAX_K*MAX_R*MAX_N*sizeof(float));


        total_power = 0;
        current_eval = 0;
    }

    State() { }

    INLINE void set_power(int t, int k, int r, int n, float v) {
        if (power[t][k][r][n] == v) return;
        // total_power += v - power[t][k][r][n];
        power_usage[t][k] += v - power[t][k][r][n];
        power_usage_single[t][k][r] += v - power[t][k][r][n];
        if (v == 0) {
            total_r[t][k][n]--;
            total_kr[t][n]--;
            int p2d = power_2dlist_pos[t][r][n][k];
            power_2dlist[t][r][p2d] = power_2dlist[t][r][--power_2dlist_len[t][r]];
            power_2dlist_pos[t][r][power_2dlist[t][r][p2d].X][power_2dlist[t][r][p2d].Y] = p2d;
            power_2dlist_pos[t][r][n][k] = -1;
            int p1d = power_1dlist_pos[t][k][r][n];
            power_1dlist[t][k][r][p1d] = power_1dlist[t][k][r][--power_1dlist_len[t][k][r]];
            power_1dlist_pos[t][k][r][power_1dlist[t][k][r][p1d]] = p1d;
            power_1dlist_pos[t][k][r][n] = -1;
        } else if (power[t][k][r][n] == 0) {
            total_r[t][k][n]++;
            total_kr[t][n]++;
            power_2dlist_pos[t][r][n][k] = power_2dlist_len[t][r];
            power_2dlist[t][r][power_2dlist_len[t][r]++] = MP(n, k);
            power_1dlist_pos[t][k][r][n] = power_1dlist_len[t][k][r];
            power_1dlist[t][k][r][power_1dlist_len[t][k][r]++] = n;
        }
        updated_group[t][r] = 1;
        power[t][k][r][n] = v;
    }

    INLINE double eval_frame(int j) {
        return frame_data[j] >= frames[j].tbs ? 1.0 : 0.10 * pow(frame_data[j] / frames[j].tbs, 3.50);        
    }

    void calc_user_data(int t) {
        for (int n : good_users[t]) {
            current_eval -= eval_frame(frame_id[t][n]);
            frame_data[frame_id[t][n]] = max(frame_data[frame_id[t][n]] - user_data[t][n], 0.0f);
            current_eval += eval_frame(frame_id[t][n]);
            old_user_data[t][n] = user_data[t][n];
            user_data[t][n] = 0;
        }

        for (int n0 : good_users[t]) if (total_kr[t][n0]) REP(k0, K) {

            if (total_r[t][k0][n0] == 0) continue;

            float total_s = 1;
            
            REP(r, R) if (power[t][k0][r][n0] > 0) {
                if (updated_group[t][r]) {
                    float q = 1.0;
                    float p = sinr[t][k0][r][n0] * power[t][k0][r][n0];
                    REP(i, power_2dlist_len[t][r]) {
                        auto n = power_2dlist[t][r][i].X;
                        auto k = power_2dlist[t][r][i].Y;
                        if (k == k0) {
                            p /= d[k0][r][n0][n];
                        } else {
                            q += (n0 != n) * sinr[t][k][r][n0] * power[t][k][r][n] * d[k][r][n0][n];
                        }
                    }

                    cur_data[t][k0][r][n0] = p / q;
                }
                total_s *= cur_data[t][k0][r][n0];

            }
            total_s = pow(total_s, 1.0 / total_r[t][k0][n0]);

            user_data[t][n0] += total_r[t][k0][n0] * log2(1.0 + total_s);

        }

        REP(r, R) updated_group[t][r] = 0;

        for (int n : good_users[t]) {
            current_eval -= eval_frame(frame_id[t][n]);
            frame_data[frame_id[t][n]] += user_data[t][n];
            current_eval += eval_frame(frame_id[t][n]);
        }
    }

    void restore_user_data(int t) {
        for (int n : good_users[t]) if (user_data[t][n] != old_user_data[t][n]) {
            frame_data[frame_id[t][n]] = max(frame_data[frame_id[t][n]] + old_user_data[t][n] - user_data[t][n], 0.0f);
            user_data[t][n] = old_user_data[t][n];
        }

    }

    double eval() {
        return current_eval;
    }

    double reeval() {
        current_eval = 0;
        REP(j, J) current_eval += eval_frame(j);
        return current_eval;
    }

    double recalc_user_data() {
        REP(t, T) calc_user_data(t);
        return eval();
    }

    INLINE double calc_score() {
        int total_frames = 0;
        REP(j, J) total_frames += frame_data[j] >= frames[j].tbs;
        return total_frames;
    }


};

State *s;

int main(int argc, char **argv) {
    srand(time(NULL));
    int offset = rand() % 1000;
    REP(i, offset) rng.next();
    
    // Read Input
    string line;
    ios_base::sync_with_stdio(false);
    getline(cin, line); N = atoi(line.c_str());
    getline(cin, line); K = atoi(line.c_str());
    getline(cin, line); T = atoi(line.c_str());
    getline(cin, line); R = atoi(line.c_str());

    REP(t, T) REP(k, K) REP(r, R) {
        getline(cin, line);
        VS parts = splt(line);
        REP(n, N) sinr[t][k][r][n] = atof(parts[n].c_str());
    }

    double d_sum = 0;
    REP(k, K) REP(r, R) REP(n, N) {
        getline(cin, line);
        VS parts = splt(line);
        REP(i, N) {
            double v = -atof(parts[i].c_str());
            d_sum += v;
            d[k][r][n][i] = exp(v);
        }
    }
    double d_avg = d_sum / (K * R * N * N);

    getline(cin, line); J = atoi(line.c_str());
    REP(j, J) {
        getline(cin, line);
        VS parts = splt(line);
        frames[j].id = atoi(parts[0].c_str());
        frames[j].tbs = atoi(parts[1].c_str())*(1+1e-6) / W; 
        frames[j].user = atoi(parts[2].c_str());
        frames[j].start = atoi(parts[3].c_str());
        frames[j].len = atoi(parts[4].c_str());
    }
    sort(frames, frames + J, [&](const Frame &a, const Frame &b) {return a.start < b.start;});

    REP(k, K) REP(r, R) REP(n, N) d[k][r][n][n] = 1;

    s = new State();
    s->init();

    DB(N, K, T, R, J, frames[0].len);

    user_ok = VVI(T, VI(N, 0));
    good_users = VVI(T);
    frame_order = VPII();

    REP(t, T) REP(n, N) frame_id[t][n] = MAX_J;

    // calc obtainable frames
    VI obtainable(J);
    REP(j, J) {
        double data = 0;
        FOR(t, frames[j].start, frames[j].start + frames[j].len) {
            REP(k, K) {
                float best_sinr = 0;
                REP(r, R) best_sinr = max(best_sinr, sinr[t][k][r][frames[j].user]);
                float best_data = 0;
                int br = -1;
                FOR(r, 1, R+1) {
                    float v = W * r * log2(1 + best_sinr * R / r);
                    if (v > best_data) {
                        best_data = v;
                        br = r;
                    }
                    // best_data = max(best_data, W * r * log2(1 + best_sinr * R / r));
                }
                // if (br != R) DB(br);
                data += best_data;
                if (data >= frames[j].tbs) break;
            }
            if (data >= frames[j].tbs) break;
        }
        obtainable[j] = data >= frames[j].tbs;
        // DB(j, frames[j].tbs, data, data > frames[j].tbs);
    }

    REP(j, J) if (obtainable[j]) frame_order.PB(MP(frames[j].tbs, j));
    sort(ALL(frame_order));
    REP(i, frame_order.S) add_frame();
    REP(t, T) sort(ALL(good_users[t]));

    VI used_k(K, 0);
    VI used_n(N, 0);
    VI k_rows(K, 0);
    VVI used_kn(K, VI(N, 0));
    VVI rlist(R);
    VI rpos(R, -1);

    VVVD dmul = VVVD(R, VVD(N, VD(N, 1.0)));
    REP(k, K) REP(r, R) REP(n, N) REP(n2, N) dmul[r][n][n2] *= d[k][r][n][n2];

    REP(t, T) {
        if (d_avg > 0.5) break;
        if (good_users[t].S == 0) continue;

        VVD rpos_value(R, VD(K, 0));
        REP(r, R) {
            REP(k, K) {
                double total_value = 0;

                double avg_sinr = 1.0;
                for (int n : good_users[t]) avg_sinr *= sinr[t][k][r][n];
                avg_sinr = pow(avg_sinr, 1.0 / good_users[t].S);

                for (int n : good_users[t]) {
                    double cur_d = 1.0;
                    for (int n2 : good_users[t]) if (n != n2) cur_d *= d[k][r][n][n2];
                    double cur_dd = 1.0;
                    if (K > R) for (int n2 : good_users[t]) if (n != n2) cur_dd *= dmul[r][n][n2];
                    total_value += log(1 + avg_sinr / pow(cur_d, 1.0) / pow(cur_dd, 0.1));
                }

                rpos_value[r][k] = total_value;
            }
        }

        REP(r, R) rpos[r] = -1;
        REP(k, K) k_rows[k] = 0;
        REP(rstep, R) {
            double bv = -1e9;
            int br = -1;
            int bk = -1;
            REP(r, R) if (rpos[r] == -1) {
                REP(k, K) {
                    double av = rpos_value[r][k] - k_rows[k] * 1000;
                    if (av > bv) {
                        bv = av;
                        br = r;
                        bk = k;
                    }
                }
            }
            rpos[br] = bk;
            k_rows[bk]++;
        }
        REP(step, 500) {
            if (rng.next(2)) {
                int r1 = rng.next(R);
                int r2 = rng.next(R);
                if (r1 == r2) continue;
                double v1 = rpos_value[r1][rpos[r1]] * rpos_value[r2][rpos[r2]];
                double v2 = rpos_value[r1][rpos[r2]] * rpos_value[r2][rpos[r1]];
                if (v2 >= v1) swap(rpos[r1], rpos[r2]);
            } else {
                int r1 = rng.next(R);
                int nk = rng.next(K);
                if (k_rows[rpos[r1]] <= k_rows[nk]) continue;
                double v1 = rpos_value[r1][rpos[r1]];
                double v2 = rpos_value[r1][nk];
                if (v2 >= v1) rpos[r1] = nk, k_rows[nk]++, k_rows[rpos[r1]]--;
            }
        }

        REP(k, K) used_k[k] = 0;
        REP(n, N) used_n[n] = 0;
        REP(r, R) rlist[r].clear();
        REP(k, K) REP(n, N) used_kn[k][n] = 0;
        VC<pair<double,int>> user_list; for (int n : good_users[t]) user_list.PB(MP(frames[frame_id[t][n]].tbs, n));
        // random_shuffle(ALL(user_list));
        sort(ALL(user_list));
        // reverse(ALL(user_list));
        VPII build_order;
        REP(rep, max(1, min(R, min((int)R, 4000 / (int)good_users[t].S / (int)good_users[t].S)))) {
            for (auto &p : user_list) {
                int n = p.Y;
                int f = frame_id[t][n];
                double bv = -1e9;
                int bk = -1;
                int br = -1;
                REP(r, R) {
                    if (rpos[r] == -1) continue;
                    int k = rpos[r];
                    double badd = 1.0;
                    for (int x : rlist[r]) badd /= d[k][r][n][x]; 
                    REP(r2, R) if (r2 != r && rpos[r2] == k) for (int x : rlist[r2]) badd /= d[k][r][n][x];
                    // double av = pow(log(1 + sinr[t][k][r][n] * badd), pow(0.2, rng.next_double())) - rlist[r].S * 10 - used_k[k] * 100 + used_kn[k][n] * 1000;
                    double av = log(1 + sinr[t][k][r][n] * badd) + rng.next_double() - rlist[r].S * 10 - used_k[k] * 100 + used_kn[k][n] * 1000;
                    // double av = rng.next_double();
                    for (int x : rlist[r]) av -= 1e6 * (n == x);
                    if (av > bv) {
                        bv = av;
                        bk = k;
                        br = r;
                    }
                }
                used_k[bk]++;
                used_n[n]++;
                used_kn[bk][n]++;
                rpos[br] = bk;
                build_order.PB(MP(br, rlist[br].S));
                rlist[br].PB(n);
            }
        }
        REP(r, R) if (rpos[r] != -1) {
            int k = rpos[r];
            double available_power = 4.0;
            int total_k = 0; REP(r2, R) total_k += rpos[r] == rpos[r2];
            available_power = min(available_power, 1.0 * R / total_k);
            for (int n : rlist[r]) {
                if (s->frame_data[frame_id[t][n]] >= frames[frame_id[t][n]].tbs) continue;
                s->set_power(t, k, r, n, available_power / rlist[r].S);
                s->calc_user_data(t);
            }
        }
    }

    double bv = s->recalc_user_data();
    double printed_bv = bv;
    double last_print = 0.0;
    DB(bv);

    int acc = 0;
    int step = 0;
    double time_passed = 0;
    const double t0 = (frames[0].len == 1 ? 0.004 : 1e-5);
    const double tn = (frames[0].len == 1 ? 1e-6 : 1e-9);
    double temp = t0;

    const double TIME_LIMIT = (1.920 - elapsed());
    const double sa_start = get_time();

    VVD move_weights = {
        {5.0, 1.8, .5, .4, .4, .25},
        {5.0, 1.8, .5, .4, .4, .25},
    };

    REP(ph, 2) {
        double wsum = 0; for (double w : move_weights[ph]) wsum += w;
        for (double &d : move_weights[ph]) d /= wsum;
        FOR(i, 1, move_weights[ph].S) move_weights[ph][i] += move_weights[ph][i-1];
        move_weights[ph].back() = 1.0;
    }

    int type = 0;
    int attempt = 0;

    static double l1_p[MAX_N];
    static double l2_p[MAX_N];
    static int la[MAX_N];
    static int la_len = 0;

    VI times_left;
    VI frames_left(T);
    REP(t, T) for (int n : good_users[t]) if (s->frame_data[frame_id[t][n]] < frames[frame_id[t][n]].tbs) frames_left[t]++;
    REP(t, T) if (frames_left[t]) times_left.PB(t);

    int phase = 0;

    while (true) {
        int t;
        
        if (frames[0].len == 1 || rng.next_double() < 0.70) {
            if (times_left.S == 0) break;
            t = times_left[rng.next(times_left.S)];
        } else {
            t = rng.next(T);
            if (good_users[t].S == 0) continue;
        }

        int k = rng.next(K);

        int t2, k2, r, r2, cur_n, new_n, cur_n2, new_n2;
        double cur_power, cur_power2, new_power, new_power2;

        if (attempt == 0 || attempt >= 10) {
            double rtype = rng.next_double();
            type = 0; while (rtype > move_weights[phase][type]) type++;
        }
        attempt++;

        const double RR = 0.55;
        const double T0_ZERO = d_avg > 0.5 ? 0.5 : 0.30;

        if (type == 0) {
            r = rng.next(R);

            if (rng.next_double() < 0.4) {
                cur_n = s->power_1dlist_len[t][k][r] ? s->power_1dlist[t][k][r][rng.next(s->power_1dlist_len[t][k][r])] : -1;
                if (cur_n == -1) continue;
            } else {
                cur_n = good_users[t][rng.next(good_users[t].S)];
            }

            
            cur_power = s->power[t][k][r][cur_n];
            if (cur_power == 0 && s->frame_data[frame_id[t][cur_n]] > frames[frame_id[t][cur_n]].tbs && rng.next_double() < 0.85) continue;
            double max_power = max(0.0, min(R - s->power_usage[t][k] + cur_power, 4.0 - s->power_usage_single[t][k][r] + cur_power));
            new_power = rng.next_double() < T0_ZERO ? 0.0 : pow(rng.next_double(), 0.5) * max_power;
            

            if (cur_power == new_power) continue;

        } else if (type == 1) {
            r = rng.next(R);
            if (s->power_1dlist_len[t][k][r] == 0) continue;

            cur_n = s->power_1dlist[t][k][r][rng.next(s->power_1dlist_len[t][k][r])];
            cur_n2 = rng.next_double() < 0.2 ? good_users[t][rng.next(good_users[t].S)] : s->power_1dlist[t][k][r][rng.next(s->power_1dlist_len[t][k][r])];
            if (cur_n == cur_n2) continue;

            cur_power = s->power[t][k][r][cur_n];
            cur_power2 = s->power[t][k][r][cur_n2];
            double power_diff = frames[0].len > 1 && rng.next_double() < 0.8 || rng.next_double() < RR ? cur_power : rng.next_double() * cur_power;
            if (cur_power == 0 && s->frame_data[frame_id[t][cur_n2]] > frames[frame_id[t][cur_n2]].tbs && rng.next_double() < 0.50) continue;

            new_power = cur_power - power_diff;
            new_power2 = cur_power2 + power_diff;
            new_power2 = max(new_power2, 0.0);
        } else if (type == 2) {
            r = rng.next(R);
            r2 = rng.next(R);
            if (r == r2) continue;
            if (s->power_1dlist_len[t][k][r] == 0) continue;

            cur_n = s->power_1dlist[t][k][r][rng.next(s->power_1dlist_len[t][k][r])];
            cur_power = s->power[t][k][r][cur_n];
            cur_power2 = s->power[t][k][r2][cur_n];

            if (cur_power == 0 && s->frame_data[frame_id[t][cur_n]] > frames[frame_id[t][cur_n]].tbs && rng.next_double() < 0.50) continue;
            double power_diff = rng.next_double() < RR ? cur_power : rng.next_double() * cur_power;

            new_power = cur_power - power_diff;
            new_power2 = cur_power2 + power_diff;
            new_power2 = min(new_power2, 4.0 - s->power_usage_single[t][k][r2] + cur_power2);
            new_power2 = max(new_power2, 0.0);

        } else if (type == 3) {
            r = rng.next(R);
            k2 = rng.next(K);
            if (k == k2) continue;
            if (s->power_1dlist_len[t][k][r] == 0) continue;

            cur_n = s->power_1dlist[t][k][r][rng.next(s->power_1dlist_len[t][k][r])];
            cur_power = s->power[t][k][r][cur_n];
            cur_power2 = s->power[t][k2][r][cur_n];

            if (cur_power == 0 && s->frame_data[frame_id[t][cur_n]] > frames[frame_id[t][cur_n]].tbs && rng.next_double() < 0.50) continue;
            double power_diff = rng.next_double() < RR ? cur_power : rng.next_double() * cur_power;

            new_power = cur_power - power_diff;
            new_power2 = cur_power2 + power_diff;
            new_power2 = min(new_power2, 4.0 - s->power_usage_single[t][k2][r] + cur_power2);
            new_power2 = min(new_power2, R - s->power_usage[t][k2] + cur_power2);
            new_power2 = max(new_power2, 0.0);
        } else if (type == 4) {
            r = rng.next(R);
            k2 = rng.next(K);
            if (k == k2) continue;
            if (s->power_usage[t][k] - s->power_usage_single[t][k][r] + s->power_usage_single[t][k2][r] > R) continue;
            if (s->power_usage[t][k2] - s->power_usage_single[t][k2][r] + s->power_usage_single[t][k][r] > R) continue;
        } else if (type == 5) {
            r2 = rng.next(R);
            if (r2 == r) continue;
        }
        attempt = 0;

        step++;

        if ((step & 31) == 0) {
            static const double PHASE0_LENGTH = frames[0].len == 1 ? 0.2 : 0.0;
            time_passed = (get_time() - sa_start) / TIME_LIMIT;
            phase = time_passed < PHASE0_LENGTH ? 0 : 1;
            g_time_passed = time_passed < 0.2 ? 1.0 - time_passed * 5.0 : (time_passed - 0.2) / 0.8;
            temp = t0 * pow(tn / t0, pow(phase == 0 ? 1.0 : (time_passed - PHASE0_LENGTH), 1.50));
            if (time_passed > 1.0) break;
        }

        if (type == 0) {
            s->set_power(t, k, r, cur_n, new_power);
        } else if (type == 1) {
            s->set_power(t, k, r, cur_n, new_power);
            s->set_power(t, k, r, cur_n2, new_power2);
        } else if (type == 2) {
            s->set_power(t, k, r, cur_n, new_power);
            s->set_power(t, k, r2, cur_n, new_power2);
        } else if (type == 3) {
            s->set_power(t, k, r, cur_n, new_power);
            s->set_power(t, k2, r, cur_n, new_power2);
        } else if (type == 4) {
            la_len = 0;
            REP(n, N) if (s->power[t][k][r][n] > 0 || s->power[t][k2][r][n] > 0) {
                la[la_len++] = n;
                l1_p[n] = s->power[t][k][r][n];
                l2_p[n] = s->power[t][k2][r][n];
                s->set_power(t, k, r, n, l2_p[n]);
                s->set_power(t, k2, r, n, l1_p[n]);
            }
        } else if (type == 5) {
            REP(k, K) if (s->power_usage_single[t][k][r] || s->power_usage_single[t][k][r2]) for (int n : good_users[t]) if (s->power[t][k][r][n] > 0 || s->power[t][k][r2][n] > 0) {
                double p1 = s->power[t][k][r][n];
                double p2 = s->power[t][k][r2][n];
                s->set_power(t, k, r, n, p2);
                s->set_power(t, k, r2, n, p1);
            }
        }
        s->calc_user_data(t);

        double av = s->eval();

        if (av >= bv || rng.next_double() < exp((av - bv) / temp)) {
            acc++;

            if (true || av > bv + 0.8) {
                frames_left[t] = 0;
                for (int n : good_users[t]) if (s->frame_data[frame_id[t][n]] < frames[frame_id[t][n]].tbs) frames_left[t]++;
                if (frames_left[t] == 0) {
                    REP(i, times_left.S) if (times_left[i] == t) {
                        times_left[i] = times_left.back();
                        times_left.pop_back();
                        break;
                    }
                }
            }

            if ((int)av > (int)printed_bv && time_passed > last_print + 0.05) {
                DB(step, av);
                printed_bv = av;
                last_print = time_passed;
            }
            bv = av;
        } else {
            if (type == 0) {
                s->set_power(t, k, r, cur_n, cur_power);
            } else if (type == 1) {
                s->set_power(t, k, r, cur_n, cur_power);
                s->set_power(t, k, r, cur_n2, cur_power2);
            } else if (type == 2) {
                s->set_power(t, k, r, cur_n, cur_power);
                s->set_power(t, k, r2, cur_n, cur_power2);
            } else if (type == 3) {
                s->set_power(t, k, r, cur_n, cur_power);
                s->set_power(t, k2, r, cur_n, cur_power2);
            } else if (type == 4) {
                REP(i, la_len) {
                    int n = la[i];
                    s->set_power(t, k, r, n, l1_p[n]);
                    s->set_power(t, k2, r, n, l2_p[n]);
                }
            } else if (type == 5) {
                REP(k, K) if (s->power_usage_single[t][k][r] || s->power_usage_single[t][k][r2]) for (int n : good_users[t]) if (s->power[t][k][r][n] > 0 || s->power[t][k][r2][n] > 0) {
                    double p1 = s->power[t][k][r][n];
                    double p2 = s->power[t][k][r2][n];
                    s->set_power(t, k, r, n, p2);
                    s->set_power(t, k, r2, n, p1);
                }
            }
            s->current_eval = bv;
            s->restore_user_data(t);
        }
    }

    // Write Output 
    char cline[MAX_N * MAX_R * 12];
    REP(t, T) REP(k, K) {
        int cpos = 0;
        cline[cpos] = 0;
        double sumk = 0;
        REP(r, R) {
            double sumr = 0;
            REP(n, N) {
                if (n) cline[cpos++] = ' ';
                double v = max(0.0, min(4.0 - sumr, min(R - sumk, (double)s->power[t][k][r][n])));
                sumr += v;
                sumk += v;
                if (v < 1e-9) {
                    cline[cpos++] = '0';
                } else {
                    cline[cpos++] = '0' + (int)v;
                    cline[cpos++] = '.';
                    int iv = (int)((v - (int)v) * 10000000 + .5);
                    REP(i, 7) {
                        cline[cpos + 6 - i] = '0' + iv % 10;
                        iv /= 10;
                    }
                    cpos += 7;
                }
            }
            cline[cpos++] = '\n';
        }
        cline[cpos] = 0;
        cout << cline;
    }	
	return 0;
}
