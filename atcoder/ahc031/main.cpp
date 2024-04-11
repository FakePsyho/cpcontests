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
    uint64_t x = 88172645463325252ULL;
    INLINE uint32_t rand() {x ^= (x << 7); return x ^= (x >> 9);}
    INLINE int next() {return rand(); }
    INLINE int next(int x) {return ((LL)rand() * x) >> 32; }
    INLINE int next(int a, int b) {return a + next(b - a); }
    INLINE double next_double() {return (rand() + 0.5) * (1.0 / 4294967296.0); }
    INLINE double next_double(double a, double b) {return a + next_double() * (b - a); }
};
 
static RNG rng;

// SOLUTION

// #define USE_TIMERS
const bool SHOW_IMPROVEMENTS = false;

#if defined(LOCAL)
const double TIME_SCALE = 2.0*60;
#elif defined(VM)
const double TIME_SCALE = 1.3;
#else 
const double TIME_SCALE = 1.45;
#endif


#ifdef USE_TIMERS
#define TIMER_START(x) double x##_start = elapsed();
#define TIMER_END(x) x += elapsed() - x##_start;
#else
#define TIMER_START(x)
#define TIMER_END(x)
#endif

double timer1 = 0;
double timer2 = 0;
double timer3 = 0;

const int MAX_N = 50;
const int MAX_D = 50;

const int W = 1000;
int D;
int N;
VVI va;

struct XRect {
    short x, y1, y2;
    XRect() {}
    XRect(int x, int y1, int y2) : x(x), y1(y1), y2(y2) {}
    friend ostream& operator<<(ostream &os, const XRect &r) {os << "(" << r.x << ", " << r.y1 << ", " << r.y2 << ")"; return os; }
};

const bool CALC_H = true;
VI vwidth;

struct XFastBorders {
    LL vlines[W];
    LL xlines;
    int used_vlines[MAX_N];
    int n_used_vlines;

    XFastBorders() { }

    void init() {
        ZERO(vlines);
        n_used_vlines = 0;
    }

    INLINE void clear() {
        REP(i, n_used_vlines) vlines[used_vlines[i]] = 0;
        n_used_vlines = 0;
        if (CALC_H) xlines = 0;
    }

    INLINE void add(int x, int y1, int y2) {
        if (y2 < W) {
            if (vlines[y2] == 0) used_vlines[n_used_vlines++] = y2;
            vlines[y2] |= 1LL << x;
        }
        if (CALC_H) xlines |= 3LL << x;
    }

    INLINE void add(const XRect &r) {
        add(r.x, r.y1, r.y2);
    }

    INLINE void add(const VC<XRect> &vr) {
        REP(i, N) add(vr[i]);
    }

    INLINE void set(const VC<XRect> &vr) {
        clear();
        add(vr);
    }

    int calc_score(XFastBorders& b) {
        int rv = 0;
        REP(i, n_used_vlines) {
            int l = used_vlines[i];
            LL x = vlines[l] & ~b.vlines[l];
            while (x) {
                int p = __builtin_ctzll(x);
                x &= x - 1;
                rv += vwidth[p];
            }
        }
        REP(i, b.n_used_vlines) {
            int l = b.used_vlines[i];
            LL x = b.vlines[l] & ~vlines[l];
            while (x) {
                int p = __builtin_ctzll(x);
                x &= x - 1;
                rv += vwidth[p];
            }
        }
        if (CALC_H) {
            LL x = xlines ^ b.xlines;
            x &= ~(1 + (1LL << (vwidth.S + 1)));
            rv += W * __builtin_popcountll(x);
        }
        return rv;
    }
};


int calc_penalty(const VC<VC<XRect>>& sol, int d) {
    int penalty = 0;
    REP(k, N) {
        int area = vwidth[sol[d][k].x] * (sol[d][k].y2 - sol[d][k].y1);
        penalty += max(0, va[d][k] - area);
    }
    return penalty;
}

int last_penalty = 0;
LL calc_score(const VC<VC<XRect>>& sol, int penalty = -1, LL early_exit = 0) {
    if (penalty < 0) {
        penalty = 0;
        REP(d, D) REP(k, N) {
            int area = vwidth[sol[d][k].x] * (sol[d][k].y2 - sol[d][k].y1);
            penalty += max(0, va[d][k] - area);
        }
    }
    LL score = (LL)penalty * 100 + 1;
    if (early_exit && score > early_exit) return (LL)1e15;

    static XFastBorders fb[2];
    static int first = true;
    if (first) {
        REP(i, 2) fb[i].init();
        first = false;
    }

    REP(d, D) {
        XFastBorders &cur  = fb[d&1];
        XFastBorders &next = fb[1-(d&1)];
        cur.set(sol[d]);
        if (d) score += cur.calc_score(next);
        if (early_exit && score > early_exit) return (LL)1e15;
    }

    last_penalty = penalty;
    return score;
}

void full_solve() {
    pair<PII, PII> sol[MAX_N];
    int used[MAX_N];
    int left[MAX_N];
    int n_left;

    int cur[MAX_N];

    VI vmax(N, 0);
    REP(d, D) REP(i, N) vmax[i] = max(vmax[i], va[d][i]);

    double FS_TIME_LIMIT = 0.05;

    int total_sum = 0; REP(i, N) total_sum += vmax[i];

    DB(total_sum);

    int step = 0;
    while (true) {
        step++;
        if ((step & 255) == 0) {
            double time_passed = elapsed();
            if (time_passed > FS_TIME_LIMIT) break;
        }

        int x_left = W;
        int y_left = W;
        REP(i, N) left[i] = i;
        n_left = N;

        int sum_left = total_sum; 
        while (n_left) {
            if (sum_left > x_left * y_left) break;
            int group_size = rng.next(1, min(4, n_left)+1);
            REP(i, group_size) {
                int j = rng.next(n_left);
                cur[i] = left[j];
                left[j] = left[--n_left];
                sum_left -= vmax[cur[i]];
            }

            int side = rng.next(2);
            int group_sum = 0; REP(i, group_size) group_sum += vmax[cur[i]];

            int w = side ? x_left : y_left;
            int h = (group_sum + w - 1) / w;
            while (true) {
                int w_used = 0;
                REP(i, group_size) {
                    int ww = (vmax[cur[i]] + h - 1) / h;
                    w_used += ww;
                }
                if (w_used <= w) break;
                h++;
            }

            int w_used = 0;
            if (side) {
                REP(i, group_size) {
                    int ww = (vmax[cur[i]] + h - 1) / h;
                    sol[cur[i]] = MP(MP(w_used, w_used+ww), MP(y_left - h, y_left));
                    w_used += ww;
                }
                y_left -= h;
            } else {
                REP(i, group_size) {
                    int ww = (vmax[cur[i]] + h - 1) / h;
                    sol[cur[i]] = MP(MP(x_left - h, x_left), MP(w_used, w_used+ww));
                    w_used += ww;
                }
                x_left -= h;
            }

            if (n_left == 0 && x_left >= 0 && y_left >= 0) {
                DB(step);
                REP(d, D) REP(i, N) cout << sol[i].X.X << ' ' << sol[i].Y.X << ' ' << sol[i].X.Y << ' ' << sol[i].Y.Y << endl;
                exit(0);
            }
        }
    }
    DB(step);
}

void analyze() {
    REP(i, N+1) {
        int sum = 0;
        REP(j, i) {
            int x = 0; 
            REP(d, D) x = max(x, va[d][j]);
            sum += x;
        }
        int rest = 0;
        REP(d, D) {
            int x = 0;
            FOR(j, i, N) x += va[d][j];
            rest = max(x, rest);
        }
        DB(i, sum+rest, sum, rest); 
    }


}


int ab_max_cost;
int ab_pot_cost;
void analyze_bands(VC<VC<XRect>> &sol, VI &vwidth) {
    int bands = vwidth.S;

    VVI bands_used(D, VI(bands+1, 0));
    REP(d, D) REP(i, N) bands_used[d][sol[d][i].x] = 1, bands_used[d][sol[d][i].x+1] = 1;
    int missing_bands = 0;
    REP(d, D) FOR(x, 1, bands) if (!bands_used[d][x]) missing_bands++;
    // DB(missing_bands);

    int max_cost = 0;
    int current_cost = 0;
    int potential_cost = 0;
    REP(x, bands) {
        int h = vwidth[x];
        VVI v(D);
        VI vleft(D, W);
        REP(d, D) {
            VPII vp;
            REP(i, N) if (sol[d][i].x == x) vp.PB(MP(sol[d][i].y1, i));
            sort(ALL(vp), greater<PII>());
                for (PII &p : vp) {
                v[d].PB((va[d][p.Y] + h - 1) / h);
                // assert(v[d].back() <= sol[d][p.Y].y2 - sol[d][p.Y].y1);
                vleft[d] -= v[d].back();
            }
        }

        
        int total = 0;
        int used = 0;
        REP(d, D-1) {
            REP(i, N) {
                if (sol[d][i].x == x && sol[d][i].y2 < W) {
                    bool found = false;
                    REP(j, N) found |= sol[d+1][j].x == x && sol[d+1][j].y2 == sol[d][i].y2;
                    total++;
                    used += found;
                }
            }
            REP(i, N) {
                if (sol[d+1][i].x == x && sol[d+1][i].y2 < W) {
                    bool found = false;
                    REP(j, N) found |= sol[d][j].x == x && sol[d][j].y2 == sol[d+1][i].y2;
                    total++;
                    used += found;
                }
            }
        }

        // REP(d, D) assert(vleft[d] >= 0);

        int potential = 0;
        REP(d, D-1) {
            {
                int csum = 0;
                REP(i, v[d].S) {
                    if (csum) {
                        int l = csum;
                        int r = csum + vleft[d];
                        int c1sum = 0;
                        REP(j, v[d+1].S) {
                            if (c1sum) {
                                int l1 = c1sum;
                                int r1 = c1sum + vleft[d+1];
                                if (max(l, l1) <= min(r, r1)) {
                                    potential++;
                                    break;
                                }
                            }
                            c1sum += v[d+1][j];
                        }   
                    }
                    csum += v[d][i];
                }
            }

            {
                int csum = 0;
                REP(i, v[d+1].S) {
                    if (csum) {
                        int l = csum;
                        int r = csum + vleft[d+1];
                        int c1sum = 0;
                        REP(j, v[d].S) {
                            if (c1sum) {
                                int l1 = c1sum;
                                int r1 = c1sum + vleft[d];
                                if (max(l, l1) <= min(r, r1)) {
                                    potential++;
                                    break;
                                }
                            }
                            c1sum += v[d][j];
                        }   
                    }
                    csum += v[d+1][i];
                }
            }
        }

        // DB(x, h, total, used, potential, (total - used) * h);
        max_cost += total * h;
        current_cost += (total - used) * h;
        potential_cost += (total - potential) * h;
    }
    // DB(max_cost);
    // DB(current_cost);
    // DB(potential_cost);
    ab_max_cost = max_cost;
    ab_pot_cost = potential_cost;
}

const int BH = 1000;
const int BW = 1000;

int main(int argc, char **argv) {
    int _;
	cin >> _ >> D >> N;
    va = VVI(D, VI(N));
    REP(i, D) REP(j, N) cin >> va[i][j];

    cerr << "[DATA] D = " << D << endl;
    cerr << "[DATA] N = " << N << endl;

    int sum = 0; REP(i, D) REP(j, N) sum += va[i][j];
    double ratio = 1.0 * sum / W / W / D;
    cerr << "[DATA] r = " << ratio << endl;

    double m = 0.0;
    REP(i, N) {
        int mx = 0;
        REP(d, D) mx = max(va[d][i], mx);
        m += mx;
    }
    m /= W * W;
    cerr << "[DATA] m = " << m << endl;

    VI vall; REP(d, D) REP(i, N) vall.PB(va[d][i]);
    sort(ALL(vall), greater<int>());
    DB(VI(vall.begin(), vall.begin() + 10));

    // analyze();
    if (m < 1.0) full_solve();

    const double TIME_LIMIT_P1 = 0.5 * TIME_SCALE;
    const double TIME_LIMIT_P1A = 0.2 * TIME_LIMIT_P1;
    const double TIME_LIMIT_P3 = 0.7 * TIME_SCALE;
    const double TIME_LIMIT_P2 = 2.0 * TIME_SCALE;

    int startD = 0;
    int targetD = D;
    REP(i, startD) va.erase(va.begin());
    D = targetD;
    va.erase(va.begin() + D, va.end());

    DB(TIME_LIMIT_P1, TIME_LIMIT_P1A, TIME_LIMIT_P2);

    VI bvwidth = {W};
    VI xvwidth = {W};
    VC<VC<XRect>> bsol = VC<VC<XRect>>(D, VC<XRect>(N));
    LL bv = 1e15;

    VC<VC<XRect>> xsol = bsol;
    LL xv = bv;

    int p1_skip = 0;
    int p1_step = 0;
    // double p1_t0 = 100.0;
    // double p1_tn = 0.01;
    // double p1_t = p1_t0;
    double p1_time_passed = 0;


    VC<VC<XRect>> sol(D, VC<XRect>(N));

    REP(d, D) reverse(va[d].begin(), va[d].end());
    while (true) {
        again: ;
        if ((p1_step & 63) == 0) {
            p1_time_passed = elapsed() / TIME_LIMIT_P1;
            if (p1_time_passed > 1.0) break;
            static bool p1b = false;
            if (!p1b && p1_time_passed > TIME_LIMIT_P1A / TIME_LIMIT_P1) {
                DB(p1_step, bv, elapsed());
                p1b = true;
            }
            // p1_t = p1_t0 * pow(p1_tn / p1_t0, (p1_time_passed - TIME_LIMIT_P1) * (1.0 / (TIME_LIMIT_P1A - TIME_LIMIT_P1)));
        }
        p1_step++;

        if (p1_time_passed < TIME_LIMIT_P1A / TIME_LIMIT_P1 || bvwidth.S <= 2) {
            int pdiff = ratio < 0.98 ? 3 : 2;
            int min_n = max(1, (int)bvwidth.S - pdiff);
            int max_n = min(N-1, (int)bvwidth.S + pdiff) + 1;
            if (bv > 1e12 || p1_time_passed < 0.01) min_n = 1, max_n = N;
            int n = rng.next(min_n, max_n);
            static VI v;
            int reserved = (vall[0] + W - 1) / W + rng.next(10);
            v.clear();
            if (ratio < 0.9) {
                v.PB(0);
                v.PB(W - reserved);
                v.PB(W);
                REP(i, n - 2) v.PB(rng.next(1, W - reserved));
            } else {
                v.PB(0);
                v.PB(W);
                REP(i, n - 1) v.PB(rng.next(1, W));
            }
            sort(ALL(v));
            vwidth.clear();
            REP(i, v.S - 1) vwidth.PB(v[i+1] - v[i]);
        } else {
            static bool first = true;
            vwidth = bvwidth;
            int moves = rng.next(1, 3);
            while (moves--) {
                int type = min(2, rng.next(20));
                if (type == 0 && vwidth.S > 1) {
                    int a = rng.next(vwidth.S);
                    int b = rng.next(vwidth.S - 1);
                    b += b >= a;
                    vwidth[a] += vwidth[b];
                    vwidth.erase(vwidth.begin() + b);
                } else if (type == 1 && vwidth.S < N-1) {
                    int a = rng.next(vwidth.S);
                    if (vwidth[a] == 1) continue;
                    int p = rng.next(1, vwidth[a]);
                    vwidth.PB(p);
                    vwidth[a] -= p;
                } else {
                    int a = rng.next(vwidth.S);
                    int b = rng.next(vwidth.S - 1);
                    b += b >= a;
                    int p = rng.next(1, vwidth[a]);
                    vwidth[a] -= p;
                    vwidth[b] += p;
                }
            }
        }

        int n_parts = vwidth.S;
        bool ok = true;
        REP(i, n_parts) if (vwidth[i] == 0) ok = false;
        if (!ok) continue;

        int biggest = 0; for (int x : vwidth) biggest = max(biggest, x);
        if (biggest <= vall[0] / 1000) {
            p1_skip++;
            continue;
        }

        // if (ratio < 0.98 || p1_time_passed < TIME_LIMIT_P1A / TIME_LIMIT_P1) sort(ALL(vwidth));
        if (ratio < 0.98) sort(ALL(vwidth));

        int cur_penalty = 0;

        static VI v;
        v = VI(n_parts);
        REP(d, D) {
            REP(i, n_parts) v[i] = BW;
            REP(k, N) {
                int bp;
                int bh;
                int bpen = 1e9;

                REP(i, n_parts) if (v[i] * vwidth[i] >= va[d][k]) {
                    bp = i;
                    bh = (va[d][k] + vwidth[i] - 1) / vwidth[i];
                    bpen = 0;
                    break;
                }

                if (bpen) {
                    REP(i, n_parts) {
                        int h = (va[d][k] + vwidth[i] - 1) / vwidth[i];
                        h = min(h, v[i]);
                        int pen = max(0, va[d][k] - vwidth[i] * h);

                        if (pen < bpen) {
                            bp = i;
                            bh = h;
                            bpen = pen;
                            if (pen == 0) break;
                        }
                    }
                }

                cur_penalty += bpen;
                if (bh == 0 || (LL)cur_penalty * 100 > bv) {
                    goto again;
                }

                sol[d][k] = XRect(bp, v[bp] - bh, v[bp]);
                v[bp] -= bh;
            }
        }

        LL av = calc_score(sol, cur_penalty, bv);
        // if (av <= bv || p1_time_passed > TIME_LIMIT_P1A / TIME_LIMIT_P1 && rng.next_double() < exp((bv - av) / p1_t)) {
        if (av <= bv) {
            if (SHOW_IMPROVEMENTS && av < bv) DB(p1_step, av, n_parts, elapsed());
            bv = av;
            bsol = sol;
            bvwidth = vwidth;
            if (av < xv) {
                xv = av;
                xsol = sol;
                xvwidth = vwidth;
            }
        }

    }
    REP(d, D) reverse(va[d].begin(), va[d].end());
    REP(d, D) reverse(xsol[d].begin(), xsol[d].end());

    DB(p1_step, p1_skip, xv, bv);

    bsol = xsol;
    sol = xsol;
    bv = xv;
    vwidth = xvwidth;


    analyze_bands(xsol, vwidth);
    DB(ab_max_cost, ab_pot_cost);

    VC<VC<short>> borders(D); REP(i, D) REP(j, N) borders[i].PB(N-1-j);
    VC<VC<short>> orders = borders;
    VC<VC<short>> badd(D, VC<short>(N));
    VC<VC<short>> add = badd;

    VC<XRect> sol_backup = sol[0];
    VC<short> order_backup = orders[0];
    VC<short> add_backup = add[0];


    VI tc(D); 
    REP(d, D) {
        VI lc(vwidth.S);
        REP(i, N) lc[xsol[d][i].x]++;
        REP(i, vwidth.S) tc[d] += max(0, lc[i] - 1) * vwidth[i];
    }

    VI tp(D); REP(d, D) tp[d] = calc_penalty(xsol, d);

    int p1c_score = 0; REP(d, D) p1c_score += tp[d] * 100 + tc[d] * (d > 0 && d < D-1 ? 2 : 1);

    double p3_t0 = 100.0;
    double p3_tn = .5;
    double p3_t = p3_t0;

    int step3 = 0;
    int step3_acc = 0;
    if (ratio > 0.9) {
        while (true) {
            step3++;
            again3:
            if ((step3 & 63) == 0) {
                double time_passed = elapsed() / (TIME_LIMIT_P3);
                if (time_passed > 1.0) break;
                p3_t = p3_t0 * pow(p3_tn / p3_t0, time_passed);
            }

            int n_parts = vwidth.S;
            int chosen_d = rng.next(D);
            int n_moves = rng.next(1, 2);

            order_backup = orders[chosen_d];
            REP(i, n_moves) {
                int a = rng.next(N);
                int b = rng.next(N-1);
                b += b >= a;
                swap(order_backup[a], order_backup[b]);
            }

            auto &order = order_backup;
            static VI v(n_parts);
            REP(i, n_parts) v[i] = BW;

            static VI lc(n_parts);
            REP(i, n_parts) lc[i] = 0;

            int tp0 = 0;
            int tc0 = 0;

            REP(k, N) {
                int bp;
                int bh;
                int bpen = 1e9;

                REP(i, n_parts) if (v[i] * vwidth[i] >= va[chosen_d][order[k]]) {
                    bp = i;
                    bh = (va[chosen_d][order[k]] + vwidth[i] - 1) / vwidth[i] + add[chosen_d][order[k]];
                    bh = min(bh, v[i]);
                    bpen = max(0, va[chosen_d][order[k]] - vwidth[i] * bh);
                    break;
                }

                if (bpen) {
                    REP(i, n_parts) {
                        int h = (va[chosen_d][order[k]] + vwidth[i] - 1) / vwidth[i] + add[chosen_d][order[k]];
                        h = min(h, v[i]);
                        int pen = max(0, va[chosen_d][order[k]] - vwidth[i] * h);

                        if (pen < bpen) {
                            bp = i;
                            bh = h;
                            bpen = pen;
                        }
                    }
                }

                tp0 += bpen;
                // if (bh == 0 || av + (LL)new_partial_penalty * 100 > bv) {
                if (bh == 0 || tp0 > tp[chosen_d]) {
                    goto again3;
                }

                sol_backup[order[k]] = XRect(bp, v[bp] - bh, v[bp]);
                if (lc[bp]) tc0 += vwidth[bp];
                lc[bp]++;
                v[bp] -= bh;
            }

            int diff = (tp0 - tp[chosen_d]) * 100 + (tc0 - tc[chosen_d]) * (chosen_d > 0 && chosen_d < D-1 ? 2 : 1);

            if (diff <= 0) {
            // if (diff <= 0 || rng.next_double() < exp(-diff / p3_t)) {
                tp[chosen_d] = tp0;
                tc[chosen_d] = tc0;
                orders[chosen_d] = order_backup;
                sol[chosen_d] = sol_backup;
                step3_acc++;
            }

        }
    }

    bsol = sol;
    borders = orders;

    int p1c_score_after = 0; REP(d, D) p1c_score_after += tp[d] * 100 + tc[d] * (d > 0 && d < D-1 ? 2 : 1);
    DB(p1c_score, p1c_score_after, elapsed());
    DB(step3, step3_acc);
    analyze_bands(bsol, vwidth);
    DB(ab_max_cost, ab_pot_cost);

    int p2_step = 0;
    int p2_calc = 0;
    int p2_acc = 0;

    double p2_t0 = 100.0;
    double p2_tn = .5;
    double p2_t = p2_t0;

    VC<VC<XRect>> tsol(D, VC<XRect>(N));

    VI partial_penalty = VI(D);
    REP(d, D) partial_penalty[d] = calc_penalty(bsol, d);
    VI partial_scores = VI(D+1);

    sol = bsol;
    orders = borders;
    add = badd;

    XFastBorders fb[MAX_N];
    REP(d, D) fb[d].init(), fb[d].set(sol[d]);
    REP(d, D) partial_scores[d] = d ? fb[d].calc_score(fb[d-1]) : 0;

    XFastBorders fb_tmp;
    fb_tmp.init();

    double p2_time_passed = 0;

    VVVI allv(D, VVI(vwidth.S));
    REP(d, D) {
        VPII vp; REP(i, N) vp.PB(MP(sol[d][i].y1, i));
        sort(ALL(vp), greater<PII>());
        REP(i, N) allv[d][sol[d][vp[i].Y].x].PB(vp[i].Y);
    }

    bv = 0; REP(d, D) bv += partial_penalty[d] * 100 + partial_scores[d];
    DB(bv);

    VC<VC<VC<short>>> vband(D, VC<VC<short>>(vwidth.S, VC<short>(N)));
    VC<VC<short>> n_vband(D, VC<short>(vwidth.S, 0));
    REP(d, D) REP(i, N) vband[d][sol[d][i].x][n_vband[d][sol[d][i].x]++] = i;

    int cur_day = 0;
    int cur_dir = +1;
    while (true) {
        again2: ;
        p2_step++;
        TIMER_START(timer1);
        if ((p2_step & 63) == 0) {
            p2_time_passed = (elapsed() - TIME_LIMIT_P3) / (TIME_LIMIT_P2 - TIME_LIMIT_P3);
            if (p2_time_passed >= 1.0) break;
            p2_t = p2_t0 * pow(p2_tn / p2_t0, p2_time_passed);
            cur_day += cur_dir;
            if (cur_day == D-1 || cur_day == 0) cur_dir = -cur_dir;
        }

        int n_parts = vwidth.S;
        static VI v(n_parts);

        if (ratio < 0.98 && rng.next_double() < 0.05 * pow(1 - p2_time_passed, 3)) {
            bvwidth = vwidth;

            int x1 = rng.next(vwidth.S);
            int x2 = rng.next(vwidth.S - 1);
            x2 += x2 >= x1;
            if (vwidth[x1] <= 1) continue;
            if (rng.next_double() < 0.2) {
                int p = rng.next(1, vwidth[x1]);
                vwidth[x1] -= p;
                vwidth[x2] += p;
            } else {
                vwidth[x1]--;
                vwidth[x2]++;
            }

            REP(i, n_parts) assert(vwidth[i] > 0);

            int cur_penalty = 0;

            REP(d, D) {
                REP(i, n_parts) v[i] = BW;
                VC<short> &order = orders[d];

                REP(k, N) {
                    int bp;
                    int bh;
                    int bpen = 1e9;

                    REP(i, n_parts) if (v[i] * vwidth[i] >= va[d][order[k]]) {
                        bp = i;
                        bh = (va[d][order[k]] + vwidth[i] - 1) / vwidth[i] + add[d][order[k]];
                        bh = min(bh, v[i]);
                        bpen = max(0, va[d][order[k]] - vwidth[i] * bh);
                        break;
                    }

                    if (bpen) {
                        REP(i, n_parts) {
                            int h = (va[d][order[k]] + vwidth[i] - 1) / vwidth[i] + add[d][order[k]];
                            h = min(h, v[i]);
                            int pen = max(0, va[d][order[k]] - vwidth[i] * h);

                            if (pen < bpen) {
                                bp = i;
                                bh = h;
                                bpen = pen;
                            }
                        }
                    }

                    cur_penalty += bpen;
                    if (bh == 0 || (LL)cur_penalty * 100 > bv) {
                        vwidth = bvwidth;
                        goto again2;
                    }

                    tsol[d][order[k]] = XRect(bp, v[bp] - bh, v[bp]);
                    v[bp] -= bh;
                }
            }

            LL av = calc_score(tsol, cur_penalty, bv);

            if (av <= bv) {
                sol = tsol;
                REP(d, D) fb[d].set(sol[d]);
                REP(d, D) partial_scores[d] = d ? fb[d].calc_score(fb[d-1]) : 0;
                REP(d, D) partial_penalty[d] = calc_penalty(sol, d);

                int score = 1;
                REP(d, D) score += partial_penalty[d] * 100 + partial_scores[d];

                bv = av;
                if (av < xv) {
                    if (SHOW_IMPROVEMENTS) DB(p2_step, av, elapsed());
                    xvwidth = vwidth;
                    xsol = sol;
                    xv = av;
                }
                REP(d, D) REP(i, n_parts) n_vband[d][i] = 0;
                REP(d, D) REP(i, N) vband[d][sol[d][i].x][n_vband[d][sol[d][i].x]++] = i;
            } else {
                vwidth = bvwidth;
            }
            continue;
        }

        int n_swaps = rng.next_double() < 0.5;
        int n_add = rng.next(0, ratio < 0.95 ? 6 : 4);
        if (n_swaps == 0 && n_add == 0) continue;

        // int chosen_d = (p2_step >> 8) % D;
        int chosen_d = cur_day;
        order_backup = orders[chosen_d];
        add_backup = add[chosen_d];

        REP(i, n_swaps) {
            int a = rng.next(N);
            int b = rng.next(N-1);
            b += b >= a;
            swap(orders[chosen_d][a], orders[chosen_d][b]);
        }
        
        REP(i, n_add) {
            int a = rng.next(N);
            int x = sol[chosen_d][a].x;
            int mn = -add[chosen_d][a] - 1;
            int mx = W;
            int cnt = n_vband[chosen_d][x] - 1;
            REP(j, n_vband[chosen_d][x]) mx = min(mx, (int)sol[chosen_d][vband[chosen_d][x][j]].y1);
            if (cnt == 0 || add[chosen_d][a] <= 0 || rng.next_double() < 0.2) {
                int v = rng.next_double() < 0.2 ? mn + 1 : rng.next(mn, mx + 1);
                add[chosen_d][a] += v;     
            } else {
                int c = rng.next(cnt);
                int b = vband[chosen_d][x][c];
                if (a == b) b = vband[chosen_d][x][c+1];

                int v = rng.next(add[chosen_d][a] + 1);
                add[chosen_d][a] -= v;
                add[chosen_d][b] += v;
            }
        }

        TIMER_END(timer1);

        TIMER_START(timer2);

        LL av = bv;
        av -= partial_penalty[chosen_d] * 100;
        av -= partial_scores[chosen_d] + partial_scores[chosen_d+1];

        int new_partial_penalty = 0;

        REP(i, n_parts) v[i] = BW;
        VC<short> &order = orders[chosen_d];

        REP(k, N) {
            int bp;
            int bh;
            int bpen = 1e9;

            REP(i, n_parts) if (v[i] * vwidth[i] >= va[chosen_d][order[k]]) {
                bp = i;
                bh = (va[chosen_d][order[k]] + vwidth[i] - 1) / vwidth[i] + add[chosen_d][order[k]];
                bh = min(bh, v[i]);
                bpen = max(0, va[chosen_d][order[k]] - vwidth[i] * bh);
                break;
            }

            if (bpen) {
                REP(i, n_parts) {
                    int h = (va[chosen_d][order[k]] + vwidth[i] - 1) / vwidth[i] + add[chosen_d][order[k]];
                    h = min(h, v[i]);
                    int pen = max(0, va[chosen_d][order[k]] - vwidth[i] * h);

                    if (pen < bpen) {
                        bp = i;
                        bh = h;
                        bpen = pen;
                    }
                }
            }

            new_partial_penalty += bpen;
            if (bh == 0 || av + (LL)new_partial_penalty * 100 > bv) {
                // sol[chosen_d] = sol_backup;
                orders[chosen_d] = order_backup;
                add[chosen_d] = add_backup;
                TIMER_END(timer2);
                goto again2;
            }

            sol_backup[order[k]] = XRect(bp, v[bp] - bh, v[bp]);
            v[bp] -= bh;
        }
        TIMER_END(timer2);

        TIMER_START(timer3);
        av += new_partial_penalty * 100;
        fb_tmp.set(sol_backup);
        int score_d1 = chosen_d ? fb_tmp.calc_score(fb[chosen_d-1]) : 0;
        int score_d2 = chosen_d < D-1 ? fb_tmp.calc_score(fb[chosen_d+1]) : 0;
        av += score_d1 + score_d2;
        TIMER_END(timer3);

        p2_calc++;
        
        // LL zv = p2_time_passed < 0.2 ? av - score_d2 + partial_scores[chosen_d+1] : av;
        // if (zv <= bv || rng.next_double() < exp((bv - zv) / p2_t)) {
        if (av <= bv || rng.next_double() < exp((bv - av) / p2_t)) {
        // if (av <= bv) {
            p2_acc++;
            bv = av;
            partial_penalty[chosen_d] = new_partial_penalty;
            partial_scores[chosen_d] = score_d1;
            partial_scores[chosen_d+1] = score_d2;
            sol[chosen_d] = sol_backup;
            fb[chosen_d].set(sol_backup);
            // REP(i, n_parts) allv[chosen_d][i].clear();
            // REP(i, N) allv[chosen_d][sol_backup[order[i]].x].PB(order[i]);
            if (av < xv) {
                if (SHOW_IMPROVEMENTS) DB(p2_step, av, elapsed());
                xvwidth = vwidth;
                xsol = sol;
                xv = av;
            }
            REP(i, n_parts) n_vband[chosen_d][i] = 0;
            REP(i, N) vband[chosen_d][sol_backup[i].x][n_vband[chosen_d][sol_backup[i].x]++] = i;
        } else {
            orders[chosen_d] = order_backup;
            add[chosen_d] = add_backup;
        }

    }

    vwidth = xvwidth;

    if (timer1) DB(timer1);
    if (timer2) DB(timer2);
    if (timer3) DB(timer3);

    DB(bv, xv, p2_step, p2_calc, p2_acc);

    DB(vwidth);
    int width_sum = 0; REP(i, vwidth.S) width_sum += vwidth[i];
    DB(width_sum);

    int my_score = calc_score(xsol);
    DB(calc_score(xsol));

    cerr << "[DATA] parts = " << vwidth.S << endl;
    cerr << "[DATA] s_p1 = " << p1_step << endl;
    cerr << "[DATA] s_p2 = " << p2_step << endl;
    cerr << "[DATA] pen = " << last_penalty << endl;

    analyze_bands(xsol, vwidth);
    cerr << "[DATA] p2_mc = " << ab_max_cost << endl;
    cerr << "[DATA] p2_pc = " << ab_pot_cost << endl;

    VI vline; vline.PB(0); REP(i, vwidth.S) vline.PB(vline[i] + vwidth[i]);

    cerr << "[DATA] time = " << elapsed() << endl;
    cerr << "[DATA] S = " << my_score << endl;
    REP(d, D) REP(i, N) {
        bool first = true;
        REP(j, N) if (xsol[d][i].x == xsol[d][j].x && xsol[d][i].y1 > xsol[d][j].y1) first = false;
        cout << vline[xsol[d][i].x] << ' ' << (first ? 0 : xsol[d][i].y1) << ' ' << vline[xsol[d][i].x+1] << ' ' << xsol[d][i].y2 << endl;
    }

	return 0;
}
