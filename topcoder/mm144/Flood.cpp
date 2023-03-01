// Author: Psyho
// Twitter: https://twitter.com/fakepsyho

// TODO
// - bfs avoid flower?
// - hc/sa - move&expand block into neighbors?
// - try multiple start points?
// - alternative init paths based on water edge after X turns?
// - split board into separate regions (is it a good idea?)
// - early exit from sim?
// - test speed on TC machine
// - quickly simulate optimistic score?
// - try A* instead of bfs?
// - add a way to build blocks behind?
 
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
double elapsed() {return get_time() - start_time;}

// SOLUTION

const bool USE_TIMERS = false;

const int MAX_N = 32;
const int MAX_T = 5;
const int MAX_B = 5;

int H;
int W;
int T;
int B;
int X;
VS ogrid;
VPII oB;
VPII oT;

const int dr[] = {0, 1, 0, -1, -1, -1, 1, 1};
const int dc[] = {1, 0, -1, 0, -1, 1, 1, -1};
const char dx[] = "RDLU";

string get_dir(int r, int c) {return r > 0 ? "D" : r < 0 ? "U" : c > 0 ? "R" : "L";}
PII get_move(char c) {return c == 'R' ? MP(0, 1) : c == 'L' ? MP(0, -1) : c == 'D' ? MP(1, 0) : MP(-1, 0);}

double timer1 = 0;
double timer2 = 0;

static PII q[32*32][32*32];
char astar(VS &grid, PII s, PII t, int meta) {
    assert(grid[s.X][s.Y] == 'B');
    assert(grid[t.X][t.Y] != '#');
    if (grid[t.X][t.Y] == 'T') return 0;

    static int first_value = 0;
    static int first[32][32];
    static int qno[32*32];

    if (first_value == 0) ZERO(first);
    first_value += 6;

    int pos0 = abs(s.X - t.X) + abs(s.Y - t.Y);
    q[pos0][0] = s;
    qno[pos0] = 1;
    int pos = pos0;
    int max_pos = pos0;
    while (pos <= max_pos) {
        REP(i, qno[pos]) {
            PII p = q[pos][i];
            if (p == t) goto done;
            int val = pos + 1 - (abs(p.X - t.X) + abs(p.Y - t.Y));
            REP(od, 4) {
                int d = (od + (meta & 3)) & 3;
                int r = p.X + dr[d];
                int c = p.Y + dc[d];
                if (r < 0 || r >= H || c < 0 || c >= W || grid[r][c] == '#' || grid[r][c] == 'T' || first[r][c] > first_value) continue;
                first[r][c] = first[p.X][p.Y] > first_value ? first[p.X][p.Y] : first_value + d + 1;
                int npos = val + abs(r - t.X) + abs(c - t.Y);
                q[npos][qno[npos]++] = MP(r, c);
                max_pos = max(max_pos, npos);
            }
        }
        qno[pos] = 0;
        pos++;
    }
    done: ;
    FOR(i, pos, max_pos+1) qno[i] = 0;

    return first[t.X][t.Y] > first_value ? dx[first[t.X][t.Y] - first_value - 1] : 0;
}

VVS sim_solution;
template <bool GENERATE_SOLUTION> int sim(VVPII &targets, VI &meta) {
    assert(targets.S == oB.S);

    VS grid = ogrid;
    VPII vB = oB;
    VI pos(vB.S);
    if (GENERATE_SOLUTION) sim_solution.clear();

    int turn = 0;

    static PII water_tiles1[MAX_N * MAX_N];
    static PII water_tiles2[MAX_N * MAX_N];
    int n_water_tiles = 0;
    int n_water_tiles_next = 0;
    PII *water_tiles = water_tiles1;
    PII *water_tiles_next = water_tiles2;
    REP(i, oT.S) water_tiles[n_water_tiles++] = oT[i];

    while (true) {
        if (USE_TIMERS) timer1 -= get_time();

        turn++;
        if (GENERATE_SOLUTION) sim_solution.PB(VS());
        REP(i, vB.S) {
            while (pos[i] < targets[i].S && (grid[targets[i][pos[i]].X][targets[i][pos[i]].Y] == 'T' || grid[targets[i][pos[i]].X][targets[i][pos[i]].Y] == '#')) pos[i]++;
            if (pos[i] == targets[i].S) continue;
            if (grid[vB[i].X][vB[i].Y] != 'B') continue;
            
            PII p = targets[i][pos[i]];
            char path = astar(grid, vB[i], p, meta[i]);
            if (path == 0) continue;
            char c = grid[vB[i].X + get_move(path).X][vB[i].Y + get_move(path).Y];
            assert(c == '.' || c == '*' || c == 'B');
            if (c == 'B') continue;
            if (abs(vB[i].X - p.X) + abs(vB[i].Y - p.Y) == 1) {
                if (GENERATE_SOLUTION) sim_solution.back().PB(i2s(vB[i].X) + " " + i2s(vB[i].Y) + " B " + string(1, path));
                grid[p.X][p.Y] = '#';
                pos[i]++;
            } else {
                if (GENERATE_SOLUTION) sim_solution.back().PB(i2s(vB[i].X) + " " + i2s(vB[i].Y) + " M " + string(1, path));
                grid[vB[i].X][vB[i].Y] = '.';
                grid[vB[i].X + get_move(path).X][vB[i].Y + get_move(path).Y] = 'B';
                vB[i].X += get_move(path).X;
                vB[i].Y += get_move(path).Y;
            }
        }
        if (USE_TIMERS) timer1 += get_time();

        if (turn >= X) {
            if (USE_TIMERS) timer2 -= get_time();
            REP(i, n_water_tiles) REP(d, 4) {
                int nr = water_tiles[i].X + dr[d];
                int nc = water_tiles[i].Y + dc[d];
                if (nr < 0 || nr >= H || nc < 0 || nc >= W) continue;
                if (grid[nr][nc] != '.' && grid[nr][nc] != '*' && grid[nr][nc] != 'B') continue;
                grid[nr][nc] = 'T';
                water_tiles_next[n_water_tiles_next++] = MP(nr, nc);
            }
            swap(water_tiles, water_tiles_next);
            n_water_tiles = n_water_tiles_next;
            n_water_tiles_next = 0;

            if (USE_TIMERS) timer2 += get_time();
            if (n_water_tiles == 0) {
                REP(i, B) if (pos[i] < targets[i].S) targets[i].resize(pos[i]);
                break;
            }
        }
    }

    int rv = 0;
    REP(r, H) REP(c, W) {
        rv += (grid[r][c] == '.') * 1;
        rv += (grid[r][c] == '*') * 3;
        rv += (grid[r][c] == 'B') * 6;
    }

    if (GENERATE_SOLUTION) {
        REP(r, H) {
            REP(c, W) cerr << grid[r][c];
            cerr << endl;
        }
    }
    return rv;
}

void output_solution(VVS &sim_solution) {
    cout << sim_solution.S << endl;
    for (VS &vs : sim_solution) {
        cout << vs.S << endl;
        for (string &s : vs) cout << s << endl;
    }
}

int main(int argc, char **argv) {
    // Read data

    cin >> H >> W >> T >> B >> X;
    ogrid = VS(H, string(W, ' '));
    REP(r, H) REP(c, W) cin >> ogrid[r][c];

    cerr << "[DATA] H = " << H << endl;
    cerr << "[DATA] W = " << W << endl;
    cerr << "[DATA] HW = " << H*W << endl;
    cerr << "[DATA] T = " << T << endl;
    cerr << "[DATA] B = " << B << endl;
    cerr << "[DATA] S = " << X << endl;

    REP(r, H) REP(c, W) {
        if (ogrid[r][c] == 'B') oB.PB(MP(r, c));
        if (ogrid[r][c] == 'T') oT.PB(MP(r, c));
    }

    const double TIME_SCALE = 1.0;
    const double TIME_LIMIT_RANDOM = 7.0 * TIME_SCALE;
    const double TIME_LIMIT_HC = 9.5 * TIME_SCALE;

    int bv = 0;
    int xv = 0;
    VVPII btargets = VVPII(B);
    VVPII xtargets = btargets;
    VI bmeta = VI(B, 0);
    VI xmeta = bmeta;
    int step = 0;

    double time_passed = 0;
    double rnd_t0 = 10.0;
    double rnd_tn = 1e-3;
    double rnd_t = rnd_t0;
    while (true) {
        if ((step & 31) == 0) {
            time_passed = elapsed() / TIME_LIMIT_RANDOM;
            rnd_t = rnd_t0 * pow(rnd_tn / rnd_t0, time_passed);
        }
        if (time_passed > 1.0) break;
        step++;        

        VVPII targets = btargets;
        VI meta = bmeta;

        bool mode_single = step > 0 && rng.next_double() < pow(time_passed, 0.25);
        int builder = rng.next(B);

        REP(i, B) {
            if (mode_single && i != builder) continue;
            targets[i].clear();
            meta[i] = rng.rand();
            const int RELOC_DIST = 6;
            int sr = rng.next(max(0, oB[i].X - RELOC_DIST), min(H, oB[i].X + RELOC_DIST + 1));
            int sc = rng.next(max(0, oB[i].Y - RELOC_DIST), min(W, oB[i].Y + RELOC_DIST + 1));
            int r = sr, c = sc;
            int type = rng.next(100);
            if (type < 3) {
                int d0 = rng.next(4);
                REP(dd, 4) {
                    int d = (d0 + dd) % 4;
                    int r = sr + dr[d];
                    int c = sc + dc[d];
                    if (r >= 0 && r < H && c >= 0 && c < W && (ogrid[r][c] == '.' || ogrid[r][c] == '*' || ogrid[r][c] == 'B')) targets[i].PB(MP(r, c));
                }
            } else if (type < 5) {
                int d0 = rng.next(8);
                REP(dd, 8) {
                    int d = (d0 + dd) % 8;
                    int r = sr + dr[d] * 2;
                    int c = sc + dc[d] * 2;
                    if (r >= 0 && r < H && c >= 0 && c < W && (ogrid[r][c] == '.' || ogrid[r][c] == '*' || ogrid[r][c] == 'B')) targets[i].PB(MP(r, c));
                }
            } else if (type < 20) {
                int sd = rng.next(0, 8);
                while (true) {
                    if (ogrid[r][c] == '.' || ogrid[r][c] == '*' || ogrid[r][c] == 'B') targets[i].PB(MP(r, c));
                    r += dr[sd];
                    c += dc[sd];
                    if (r < 0 || r >= H || c < 0 || c >= W) break;
                }
            } else if (type < 50) {
                int sd = rng.next(0, 8);
                int len = rng.next(1, 10);
                int sd2 = ((sd % 4 + (rng.next(2) ? 1 : -1) + 4) % 4) + (sd >= 4 ? 4 : 0);
                while (len--) {
                    if ((len > 0 || sd >= 4) && (ogrid[r][c] == '.' || ogrid[r][c] == '*' || ogrid[r][c] == 'B')) targets[i].PB(MP(r, c));
                    r += dr[sd];
                    c += dc[sd];
                    if (r < 0 || r >= H || c < 0 || c >= W) break;
                }
                if (r >= 0 && r < H && c >= 0 && c < W) {
                    while (true) {
                        if (ogrid[r][c] == '.' || ogrid[r][c] == '*' || ogrid[r][c] == 'B') targets[i].PB(MP(r, c));
                        r += dr[sd2];
                        c += dc[sd2];
                        if (r < 0 || r >= H || c < 0 || c >= W) break;
                    }
                }
            } else {
                int sd = rng.next(0, 8);
                int len = rng.next(1, 10);
                int sd2 = ((sd % 4 + (rng.next(2) ? 1 : -1) + 4) % 4) + (sd >= 4 ? 4 : 0);
                int len2 = rng.next(1, 10);
                int sd3 = ((sd2 % 4 + (rng.next(2) ? 1 : -1) + 4) % 4) + (sd2 >= 4 ? 4 : 0);
                while (len--) {
                    if ((len > 0 || sd >= 4) && (ogrid[r][c] == '.' || ogrid[r][c] == '*' || ogrid[r][c] == 'B')) targets[i].PB(MP(r, c));
                    r += dr[sd];
                    c += dc[sd];
                    if (r < 0 || r >= H || c < 0 || c >= W) break;
                }
                if (r >= 0 && r < H && c >= 0 && c < W) {
                    while (len2--) {
                        if ((len2 > 0 || sd2 >= 4) && (ogrid[r][c] == '.' || ogrid[r][c] == '*' || ogrid[r][c] == 'B')) targets[i].PB(MP(r, c));
                        r += dr[sd2];
                        c += dc[sd2];
                        if (r < 0 || r >= H || c < 0 || c >= W) break;
                    }
                }
                if (r >= 0 && r < H && c >= 0 && c < W) {
                    while (true) {
                        if (ogrid[r][c] == '.' || ogrid[r][c] == '*' || ogrid[r][c] == 'B') targets[i].PB(MP(r, c));
                        r += dr[sd3];
                        c += dc[sd3];
                        if (r < 0 || r >= H || c < 0 || c >= W) break;
                    }
                }
            }
        }
        int av = sim<false>(targets, meta);
        if (av >= bv || rng.next_double() < exp((av - bv) / rnd_t)) {
            if (av > xv) {
                DB(step, elapsed(), av);
                xv = av;
                xtargets = targets;
                xmeta = meta;
            }
            bv = av;
            btargets = targets;
            bmeta = meta;
        }
    }

    int hc_step = 0;
    double sa_t0 = 0.5;
    double sa_tn = 1e-6;
    double sa_t = sa_t0;
    int lowest_accept = 0;
    while (true) {
        if ((hc_step & 31) == 0) {
            time_passed = (elapsed() - TIME_LIMIT_RANDOM) / (TIME_LIMIT_HC - TIME_LIMIT_RANDOM);
            sa_t = sa_t0 * pow(sa_tn / sa_t0, time_passed);
        }
        if (time_passed > 1.0) break;

        hc_step++;

        VVPII targets = btargets;
        VI meta = bmeta;

        int b = rng.next(B);
        int type = rng.next(4);
        if (type == 2 && time_passed < 0.1) continue;
        if (type == 3 && time_passed < 0.2) continue;
        if (type == 0 && targets[b].S > 0) {
            int pos = rng.next(targets[b].S);
            targets[b].erase(targets[b].begin() + pos);
        } else if (type == 1 && targets[b].S) {
            int pos = rng.next(targets[b].S);
            int d = rng.next(4);
            targets[b][pos].X += dr[d];
            targets[b][pos].Y += dc[d];
            if (targets[b][pos].X < 0 || targets[b][pos].X >= H || targets[b][pos].Y < 0 || targets[b][pos].Y >= W || ogrid[targets[b][pos].X][targets[b][pos].Y] == '#') continue;
        } else if (type == 2) {
            int pos = rng.next(targets[b].S + 1);
            int r, c;
            do {
                r = rng.next(H);
                c = rng.next(W);
            } while (ogrid[r][c] == '#');
            targets[b].insert(targets[b].begin() + pos, MP(r, c));
        } else if (type == 3 && targets[b].S > 1) {
            int pos = rng.next(targets[b].S);
            int pos2 = rng.next(targets[b].S);
            swap(targets[b][pos], targets[b][pos2]);
        } else if (type == 4 && targets[b].S) {
            int pos = rng.next(targets[b].S);
            int d = rng.next(4);
            int r0 = targets[b][pos].X;
            int c0 = targets[b][pos].Y;
            targets[b][pos].X += dr[d];
            targets[b][pos].Y += dc[d];
            if (targets[b][pos].X < 0 || targets[b][pos].X >= H || targets[b][pos].Y < 0 || targets[b][pos].Y >= W || ogrid[targets[b][pos].X][targets[b][pos].Y] == '#') continue;
            REP(d2, 4) {
                int r = targets[b][pos].X + dr[d2];
                int c = targets[b][pos].Y + dc[d2];
                if (r == r0 && c == c0) continue;
                if (rng.next(2) && r >= 0 && r < H && c >= 0 && c < W && (ogrid[r][c] == '.' || ogrid[r][c] == '*' || ogrid[r][c] == 'B')) targets[b].insert(targets[b].begin() + pos + 1, MP(r, c));
            }
            targets[b].erase(targets[b].begin() + pos);
        } else if (type == 5 && targets[b].S) {
            int pos = rng.next(targets[b].S);
            int r = rng.next(H);
            int c = rng.next(W);
            targets[b][pos] = MP(r, c);
            if (ogrid[r][c] == '#') continue;
        } else if (type == 6 && B > 1) {
            int b2 = rng.next(B-1);
            b2 += b2 >= b;
            if (targets[b].S == 0) continue;
            if (targets[b2].S == 0) continue;
            int pos1 = rng.next(targets[b].S);
            int pos2 = rng.next(targets[b2].S);
            swap(targets[b][pos1], targets[b2][pos2]);
        }

        int av = sim<false>(targets, meta);
        if (av >= bv || rng.next_double() < exp((av - bv) / sa_t)) {
            if (av > xv) {
                DB(hc_step, elapsed(), av);
                xv = av;
                xtargets = targets;
                xmeta = meta;
            }
            lowest_accept = max(lowest_accept, bv - av);
            bv = av;
            btargets = targets;
            bmeta = xmeta;
        }
    }
    DB(xtargets);
    DB(lowest_accept);

    cerr << "[DATA] t_bfs = " << timer1 << endl;
    cerr << "[DATA] t_water = " << timer2 << endl;
    cerr << "[DATA] rnd_steps = " << step << endl;
    cerr << "[DATA] hc_steps = " << hc_step << endl;

    DB(step, hc_step, xv);

    sim<true>(xtargets, xmeta);
    output_solution(sim_solution);
	return 0;
}
