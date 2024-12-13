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
 
const double PI = 2 * acos(0.0);
struct RNG {
    uint64_t x;
    RNG(uint64_t seed = 88172645463325252ULL) {x = seed;}
    INLINE uint32_t rand() {x ^= (x << 7); return x ^= (x >> 9);}
    INLINE int next() {return rand(); }
    INLINE int next(int x) {return ((LL)rand() * x) >> 32; }
    INLINE int next(int a, int b) {return a + next(b - a); }
    INLINE double next_double() {return (rand() + 0.5) * (1.0 / 4294967296.0); }
    INLINE double next_double(double a, double b) {return a + next_double() * (b - a); }
    INLINE double next_normal(double mean, double sigma) {
        double u1 = next_double();
        double u2 = next_double();
        double r = sqrt(-2 * log(u1));
        double theta = 2 * PI * u2;
        return mean + sigma * r * cos(theta);
    }
};
 
static RNG rng;

// SOLUTION

const int MIN_N = 30;
const int MAX_N = 100;
const int MIN_SIGMA = 1000;
const int MAX_SIGMA = 10000;
const int MIN_DIM_LO = 10000;
const int MIN_DIM_HI = 50000;
const int MAX_DIM = 100000;

#if defined(VM)
const double MACHINE_SCALE = 1.12;
#elif defined(LOCAL)
const double MACHINE_SCALE = 0.7;
#else
const double MACHINE_SCALE = 1.0;
#endif

const double TIME_LIMIT = 2.9 * MACHINE_SCALE;


int N;
int T;
int sigma;

int ow[MAX_N];
int oh[MAX_N];

VI allw[MAX_N];
VI allh[MAX_N];

struct Rect {int x1, y1, x2, y2;};

int last_max_x;
int last_max_y;

int sim(const VI& rot, const VI& align) {
    static VC<Rect> rects(N);

    int max_x = 0, max_y = 0;
    REP(i, N) {
        int w = ow[i];
        int h = oh[i];
        if (rot[i]) swap(w, h);

        int x1 = align[i] == -1 ? 0 : rects[align[i]].x2;
        int x2 = x1 + w;
        int y1 = 0; REP(j, i) if (max(rects[j].x1, x1) < min(rects[j].x2, x2)) y1 = max(y1, rects[j].y2);
        int y2 = y1 + h;

        rects[i] = {x1, y1, x2, y2};
        max_x = max(max_x, x2);
        max_y = max(max_y, y2);
    }

    last_max_x = max_x;
    last_max_y = max_y;
    return max_x + max_y;
}

int print_sol(const VI& rot, const VI& dir, const VI& align) {
    assert(rot.S == N);
    assert(align.S == N);
    cout << N << "\n";

    REP(i, N) cout << i << ' ' << rot[i] << (dir[i] ? " L " : " U ") << align[i] << "\n";
    cout.flush();

    int est_w, est_h;
    cin >> est_w >> est_h;
    return est_w + est_h;
}

int main(int argc, char **argv) {
    ios::sync_with_stdio(0);
    cin.tie(0);

    int seed = time(NULL);
    rng = RNG(seed);

    cin >> N >> T >> sigma;
    REP(i, N) cin >> ow[i] >> oh[i];

    double Tr = 1.0 * T / N;

    DATA(N);
    DATA(T);
    DATA(Tr);
    DATA(sigma);

    REP(i, N) allw[i].PB(ow[i]);
    REP(i, N) allh[i].PB(oh[i]);

    int sims = 0;

    VI order(N); REP(i, N) order[i] = i;

    VI all_ids; REP(i, N) all_ids.PB(i);
    VI used(N, 0);

    bool orientation = false;
    while (T > N * .45) {
        const int CNT = 9;

        REP(i, 2*CNT+1) {
            REP(_, 5) {
                int pos = rng.next(i, N);
                if (used[all_ids[pos]] <= used[all_ids[i]]) 
                    swap(all_ids[i], all_ids[pos]);
            }
            used[all_ids[i]]++;
        }
        VI ids(all_ids.begin(), all_ids.begin() + 2*CNT+1);
        sort(ALL(ids));

        int hsum = 0; REP(i, CNT) hsum += orientation ? oh[ids[i]] : ow[ids[i]];
        int vsum = 0; FOR(i, CNT, 2*CNT+1) vsum += orientation ? ow[ids[i]] : oh[ids[i]];

        const int SIM_NO = 1000;
        int hmax_sum = 0;
        REP(_, SIM_NO) {
            int hmax = 0;
            FOR(i, CNT, 2*CNT+1) hmax = max(hmax, (int)rng.next_normal(orientation ? oh[ids[i]] : ow[ids[i]], sigma));
            hmax_sum += hmax;
        }
        hsum += hmax_sum / SIM_NO;

        cout << 2*CNT+1 << "\n";
        REP(i, CNT) cout << ids[i] << " " << orientation << " U " << (i ? ids[i-1] : -1) << "\n";
        FOR(i, CNT, 2*CNT+1) cout << ids[i] << " " << orientation << " U " << ids[CNT-1] << "\n";
        cout.flush();

        int w, h; cin >> w >> h;

        int hdiff = w - hsum;
        int vdiff = h - vsum;

        const double convergence_rate = 0.75;

        if (orientation) {
            REP(i, CNT) oh[ids[i]] += round(1.0 * hdiff / CNT * convergence_rate);
            FOR(i, CNT, 2*CNT+1) ow[ids[i]] += round(1.0 * vdiff / (CNT+1) * convergence_rate);
        } else {
            REP(i, CNT) ow[ids[i]] += round(1.0 * hdiff / CNT * convergence_rate);
            FOR(i, CNT, 2*CNT+1) oh[ids[i]] += round(1.0 * vdiff / (CNT+1) * convergence_rate);
        }
        orientation = !orientation;

        T--;
    }

    double probe_time = elapsed();
    DATA(probe_time);

    // REP(i, N) {
    //     ow[i] = 0; for (int x : allw[i]) ow[i] += x; ow[i] = (ow[i] + allw[i].S / 2) / allw[i].S;
    //     oh[i] = 0; for (int x : allh[i]) oh[i] += x; oh[i] = (oh[i] + allh[i].S / 2) / allh[i].S;
    // }

    LL sum = 0; REP(i, N) sum += ow[i] + oh[i];
    DB(sum / 2 / N);
    int U = max(MIN_DIM_LO, min(MIN_DIM_HI, (int)(2 * sum / 2 / N - MAX_DIM)));
    DATA(U);

    LL total_area = 0; REP(i, N) total_area += (LL)ow[i] * oh[i];
    int avg_width = sqrt(total_area);

    REP(i, N) ow[i] = max(max(MIN_DIM_LO+2000,U+1000), min(ow[i], MAX_DIM-3000));
    REP(i, N) oh[i] = max(max(MIN_DIM_LO+2000,U+1000), min(oh[i], MAX_DIM-3000));

    int exp_lines = int(sqrt(N));

    int sum_runs = 0;

    REP(run, T) {
        static VI best_rot(N);
        static VI best_align(N);
        static VI best_dir(N);
        REP(i, N) best_dir[i] = 0;
        REP(i, N) best_rot[i] = rng.next(2);
        int best_bad = 1 << 20;
        int best_width = 0;
        int bv = 1 << 30;

        while (true) {
            const int width = max(MAX_DIM, (int)(avg_width * pow(2, rng.next_double(0, .2))));

            static VI rot(N); 
            static VI align(N); 

            if (rng.next(2)) {
                REP(i, N) rot[i] = rng.next(2);
            } else if (rng.next(2)) {
                REP(i, N) rot[i] = ow[i] > oh[i];
            } else {
                REP(i, N) rot[i] = oh[i] > ow[i];
            }

            int total_width = 0; REP(i, N) total_width += rot[i] ? oh[i] : ow[i];

            int min_lines = 1 + total_width / width;
            static VI line_width(N);
            static VI line_prev(N);
            int max_line = 0;

            const int nextline_max = rng.next(max(0, min_lines-2), min_lines);
            const double nextline_min_margin = MAX_DIM * pow(2, rng.next_double(-4, 1));
            const double nextline_chance = rng.next_double();

            int max_width = 0;
            REP(i, N) {
                int w = ow[i];
                int h = oh[i];
                if (rot[i]) swap(w, h);

                int cur_line = 0;
                while (cur_line < max_line && line_width[cur_line] + w > width) cur_line++;

                int best_line = cur_line;
                while (cur_line < max_line && cur_line < nextline_max) {
                    cur_line++;
                    if (line_width[cur_line] + w + nextline_min_margin > line_width[cur_line-1]) break;
                    if (rng.next_double() < nextline_chance) best_line = cur_line;
                }

                if (best_line == max_line) {
                    line_width[best_line] = 0;
                    line_width[best_line+1] = 0;
                    line_prev[best_line] = -1;
                    max_line++;
                }

                align[i] = line_prev[best_line];
                line_width[best_line] += w;
                line_prev[best_line] = i;
                max_width = max(max_width, line_width[best_line]);
            }

            if (max_line > exp_lines) continue;

            int av = sim(rot, align);
            if (av < bv) {
                bv = av;
                best_rot = rot;
                best_align = align;
                best_width = width;
            }

            sims++;
            if ((sims & 31) == 0 && elapsed() > TIME_LIMIT * (run + .05) / T) break;
        }

        const double t0 = 50000;
        const double tn = 500;

        const double sa_start = elapsed();
        const double sa_end = max(sa_start + 1e-6, TIME_LIMIT * (run + 1) / T);
        double t = t0;

        static VI xbest_align(N);
        static VI xbest_rot(N);
        int xv = 1 << 30;

        int n_lines = exp_lines;
        static VI box_line(N); 
        
        VI line_width(n_lines);
        VI cur_width(n_lines);

        const double sigma_threshold = N < 50 ? 2.0 : 3.0;

        int sa_acc = 0;
        int W_DELTA_PRUNING = MAX_DIM / 5; 
        double time_passed = 0.0;

        const double t0_w = 2.0;
        const double t1_w = 2.0;
        const double t2_w = 1.0;
        const double t3_w = 1.0;
        const double t4_w = 0.0;
        const double allt_w = t0_w + t1_w + t2_w + t3_w + t4_w;

        const double r0 = t0_w / allt_w;
        const double r1 = r0 + t1_w / allt_w;
        const double r2 = r1 + t2_w / allt_w;
        const double r3 = r2 + t3_w / allt_w;

        int max_line = 0;
        REP(i, N) {
            if (best_align[i] == -1) {
                box_line[i] = max_line++;
            } else {
                box_line[i] = box_line[best_align[i]];
            }
        }
        assert(max_line <= n_lines);

        bv = 1 << 30;
        xv = 1 << 30;

        REP(i, N) cur_width[box_line[i]] += best_rot[i] ? oh[i] : ow[i];
        int max_width = 0; REP(i, n_lines) max_width = max(max_width, cur_width[i]);

        bool first_step = true;
        while (true) {
            static VI align(N);
            static VI prev(N);

            sims++;
            if ((sims & 15) == 0) {
                time_passed = (elapsed() - sa_start) / (sa_end - sa_start);
                t = t0 * pow(tn / t0, pow(time_passed, 0.5));
                if (time_passed > 1.0) break;
                W_DELTA_PRUNING = MAX_DIM / (5 + 35 * time_passed);
            }

            int p0, v0, p1, v1;
            double r = rng.next_double();
            int type = r < r0 ? 0 : r < r1 ? 1 : r < r2 ? 2 : r < r3 ? 3 : 4;

            if (first_step) type = 10;
            if (type == 0) {
                p0 = rng.next(N);
                v0 = box_line[p0];
                int new_line = rng.next(min(n_lines, 1 + p0 * n_lines * 8 / 5 / N));
                if (new_line == v0) continue;
                if (cur_width[new_line] + (best_rot[p0] ? oh[p0] : ow[p0]) > max_width + W_DELTA_PRUNING) continue;
                box_line[p0] = new_line;
            } else if (type == 1) {
                const int t1_range = 8;
                p0 = rng.next(N);
                p1 = rng.next(max(0, p0 - t1_range), min(N, p0 + t1_range + 1));
                if (box_line[p0] == box_line[p1]) continue;
                int w0 = best_rot[p0] ? oh[p0] : ow[p0];
                int w1 = best_rot[p1] ? oh[p1] : ow[p1];
                if (cur_width[box_line[p0]] - w0 + w1 > max_width + W_DELTA_PRUNING) continue;
                if (cur_width[box_line[p1]] - w1 + w0 > max_width + W_DELTA_PRUNING) continue;
                swap(box_line[p0], box_line[p1]);
            } else if (type == 2) {
                p0 = rng.next(N);
                if (cur_width[box_line[p0]] + (best_rot[p0] ? ow[p0] - oh[p0] : oh[p0] - ow[p0]) > max_width + W_DELTA_PRUNING) continue;
                best_rot[p0] = 1 - best_rot[p0];
            } else if (type == 3) {
                p0 = rng.next(N);
                v0 = box_line[p0];
                int new_line = rng.next(min(n_lines, 1 + p0 * n_lines * 8 / 5 / N));
                if (new_line == v0) continue;
                if (cur_width[new_line] + (best_rot[p0] ? ow[p0] : oh[p0]) > max_width + W_DELTA_PRUNING) continue;
                box_line[p0] = new_line;
                best_rot[p0] = 1 - best_rot[p0];
            } 

            REP(i, n_lines) prev[i] = -1, line_width[i] = 0;
            int bad = 0;
            REP(i, N) {
                FOR(j, box_line[i]+1, n_lines) bad += line_width[j] && line_width[j] + sigma_threshold * sigma > line_width[box_line[i]];
                line_width[box_line[i]] += best_rot[i] ? oh[i] : ow[i];
                align[i] = prev[box_line[i]];
                prev[box_line[i]] = i;
            }

            int av = sim(best_rot, align);
            av += bad * 50000;

            if (av <= bv || rng.next_double() < exp((bv - av) / t)) {
                first_step = false;

                sa_acc++;
                bv = av;

                best_bad = bad;
                best_align = align;

                if (av < xv) {
                    xbest_align = align;
                    xbest_rot = best_rot;
                    xv = av;
                }

                REP(i, n_lines) cur_width[i] = line_width[i];
                max_width = 0; REP(i, n_lines) max_width = max(max_width, cur_width[i]);
            } else {
                if (type == 0) {
                    box_line[p0] = v0;
                } else if (type == 1) {
                    swap(box_line[p0], box_line[p1]);
                } else if (type == 2) {
                    best_rot[p0] = 1 - best_rot[p0];
                } else if (type == 3) {
                    box_line[p0] = v0;
                    best_rot[p0] = 1 - best_rot[p0];
                } else if (type == 4) {
                    box_line[p0] = v0;
                }
            }
        }

        DB(sa_acc, bv, xv);
        best_align = xbest_align;
        best_rot = xbest_rot;
        bv = sim(best_rot, best_align);

        int rv = print_sol(best_rot, best_dir, best_align);
        static int last_sims = 0;
        int n_sims = sims - last_sims;
        double width_ratio = 1.0*best_width/avg_width;
        DB(run, rv, bv, 1.0*rv/bv, n_sims, width_ratio, sa_acc);
        last_sims = sims;
        sum_runs += bv;
    }

    int avg_run = sum_runs / T;
    DATA(avg_run);
    DATA(sims);

    double time = elapsed();
    DATA(time);
	
	return 0;
}
