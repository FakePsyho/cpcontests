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
 
struct RNG {
    uint64_t x;
    RNG(uint64_t seed = 88172645463325252ULL) {x = seed;}
    INLINE uint32_t rand() {x ^= (x << 7); return x ^= (x >> 9);}
    INLINE int next() {return rand(); }
    INLINE int next(int x) {return ((LL)rand() * x) >> 32; }
    INLINE int next(int a, int b) {return a + next(b - a); }
    INLINE double next_double() {return (rand() + 0.5) * (1.0 / 4294967296.0); }
    INLINE double next_double(double a, double b) {return a + next_double() * (b - a); }
};
 
static RNG rng;

// SOLUTION

const int N = 5000;
const int MAX_COORD = 1e5;

const int dx[] = {1, 0, -1, 0};
const int dy[] = {0, -1, 0, 1};

VPII pgood(N);
VPII pbad(N);

const int MAX_GRID = 192;
int GRID = 24;

int xlines[MAX_GRID+1];
int ylines[MAX_GRID+1];

int grid[MAX_GRID][MAX_GRID];
int vline[MAX_GRID][MAX_GRID+1];
int hline[MAX_GRID+1][MAX_GRID];

int blocks_add[512];
int edges_add[512];
int blocks_rem[512];
int edges_rem[512];


int state[MAX_GRID][MAX_GRID];

void add_grid(VPII &vp, int value) {
    for (PII &p : vp) {
        int x = p.X;
        int y = p.Y;
        int gx = lower_bound(xlines, xlines+GRID+1, x) - xlines - 1;
        int gy = lower_bound(ylines, ylines+GRID+1, y) - ylines - 1;
        // DB(x, y, gx, gy, xlines[gx], xlines[gx+1], ylines[gy], ylines[gy+1]);
        assert(gx >= 0 && gx < GRID);
        assert(p.X >= xlines[gx] && p.X <= xlines[gx+1]);
        assert(gy >= 0 && gy < GRID);
        assert(p.Y >= ylines[gy] && p.Y <= ylines[gy+1]);
        grid[gx][gy] += value;
        // vline[gx][gy] += value;
        // vline[gx][gy+1] += value;
        // hline[gx][gy] += value;
        // hline[gx+1][gy] += value;
    }
}

int mini_state[3][3];

int calc_minicomps() {
    int comps = 0;
    VVI vs(3, VI(3, 0));
    REP(i, 3) REP(j, 3) {
        if (mini_state[i][j] == 0) continue;
        if (vs[i][j]) continue;
        comps++;
        vs[i][j] = 1;
        queue<PII> q;
        q.push({i, j});
        while (!q.empty()) {
            PII p = q.front(); q.pop();
            int x = p.X;
            int y = p.Y;
            REP(d, 4) {
                int nx = x + dx[d];
                int ny = y + dy[d];
                if (nx < 0 || nx >= 3 || ny < 0 || ny >= 3) continue;
                if (mini_state[nx][ny] == 0) continue;
                if (vs[nx][ny]) continue;
                vs[nx][ny] = 1;
                q.push({nx, ny});
            }
        }
    }
    return comps;
}

int calc_miniedges() {
    int edges = 0;
    REP(i, 3) REP(j, 3) {
        if (mini_state[i][j] == 0) continue;
        REP(d, 4) {
            int nx = i + dx[d];
            int ny = j + dy[d];
            if (nx < 0 || nx >= 3 || ny < 0 || ny >= 3) continue;
            if (mini_state[nx][ny] == 0) edges++;
        }
    }
    return edges;
}

void gen_moves() {
    REP(i, 512) {
        mini_state[0][0] = (i >> 0) & 1;
        mini_state[0][1] = (i >> 1) & 1;
        mini_state[0][2] = (i >> 2) & 1;
        mini_state[1][0] = (i >> 3) & 1;
        mini_state[1][1] = (i >> 4) & 1;
        mini_state[1][2] = (i >> 5) & 1;
        mini_state[2][0] = (i >> 6) & 1;
        mini_state[2][1] = (i >> 7) & 1;
        mini_state[2][2] = (i >> 8) & 1;

        int orig_comps = calc_minicomps();
        int orig_edges = calc_miniedges();
        if (mini_state[1][1] == 0) {
            mini_state[1][1] = 1;
            blocks_add[i] = calc_minicomps() - orig_comps;
            edges_add[i] = calc_miniedges() - orig_edges;
        } else {
            mini_state[1][1] = 0;
            blocks_rem[i] = calc_minicomps() - orig_comps;
            edges_rem[i] = calc_miniedges() - orig_edges;

            if (mini_state[1][0] && mini_state[1][2] && mini_state[0][1] && mini_state[2][1]) {
                blocks_rem[i] = -100;
                edges_rem[i] = 0;
            }
        }
    }
}

void setup_grid(int n) {
    GRID = n;

    REP(i, GRID+1) xlines[i] = (LL)MAX_COORD * i / GRID;
    REP(i, GRID+1) ylines[i] = (LL)MAX_COORD * i / GRID;

    REP(i, GRID) REP(j, GRID) grid[i][j] = 0;

    add_grid(pgood, +1);
    add_grid(pbad, -1);
}

int main(int argc, char **argv) {
    int _; cin >> _;
    REP(i, N) cin >> pgood[i].X >> pgood[i].Y;
    REP(i, N) cin >> pbad[i].X >> pbad[i].Y;

    gen_moves();

    set<int> unique_x; REP(i, N) unique_x.insert(pgood[i].X), unique_x.insert(pbad[i].X);
    set<int> unique_y; REP(i, N) unique_y.insert(pgood[i].Y), unique_y.insert(pbad[i].Y);
    DB(unique_x.S, unique_y.S);

    int tsteps = 0;

    state[GRID/2][GRID/2] = 1;
    int blocks = 1;
    int comps = 1;
    int edges = 4;
    int bv = grid[GRID/2][GRID/2];
    int xv = bv;
    int best_state[MAX_GRID][MAX_GRID];

    auto upscale_state = [&](int scale) {
        VVI new_state(GRID*scale, VI(GRID*scale, 0));
        REP(i, GRID) REP(j, GRID) {
            if (state[i][j] == 0) continue;
            REP(k, scale) REP(l, scale) new_state[i*scale+k][j*scale+l] = 1;
        }
        REP(i, GRID*scale) REP(j, GRID*scale) state[i][j] = new_state[i][j];

        edges *= scale;
    };

    auto run_sa = [&](const double time_limit, const double temp0, const double tempn, bool save_solution = false) {
        const int MAX_BLOCKS = GRID * GRID;
        const int MAX_EDGES = 4 * GRID - 1;

        double temp = temp0;

        int step = 0;

        int cache_ms[MAX_GRID][MAX_GRID];

        REP(i, GRID) REP(j, GRID) {
            int move_state = 0;
            if (i > 0 && j > 0) move_state |= state[i-1][j-1] << 0;
            if (j > 0) move_state |= state[i][j-1] << 1;
            if (i < GRID-1 && j > 0) move_state |= state[i+1][j-1] << 2;
            if (i > 0) move_state |= state[i-1][j] << 3;
            move_state |= state[i][j] << 4;
            if (i < GRID-1) move_state |= state[i+1][j] << 5;
            if (i > 0 && j < GRID-1) move_state |= state[i-1][j+1] << 6;
            if (j < GRID-1) move_state |= state[i][j+1] << 7;
            if (i < GRID-1 && j < GRID-1) move_state |= state[i+1][j+1] << 8;
            cache_ms[i][j] = move_state;
        }

        const double sa_start = elapsed();

        while (true) {
            step++;
            int x = rng.next(GRID);
            int y = rng.next(GRID);

            if ((step & 255) == 0) {
                double time_passed = (elapsed() - sa_start) / time_limit;
                if (time_passed >= 1) break;
                temp = temp0 * pow(tempn / temp0, time_passed);
            }

            int move_state = cache_ms[x][y];

            int diff_score;
            int diff_blocks;
            int diff_edges;
            int diff_comps;

            if (state[x][y] == 0) {
                if (blocks_add[move_state]) continue;
                diff_score = grid[x][y];
                // diff_blocks = 1;
                diff_edges = edges_add[move_state];
            } else {
                if (blocks_rem[move_state]) continue;
                diff_score = -grid[x][y];
                // diff_blocks = -1;
                diff_edges = edges_rem[move_state];
            }

            // if (blocks + diff_blocks > MAX_BLOCKS) continue;
            if (edges + diff_edges > MAX_EDGES) continue;

            if (diff_score >= 0 || rng.next_double() < exp(diff_score / temp)) {
                // test checkersboard
                if (x > 0 && y > 0) if (state[x][y] != state[x-1][y-1] && state[x][y] == state[x-1][y] && state[x][y] == state[x][y-1]) continue;
                if (x < GRID-1 && y > 0) if (state[x][y] != state[x+1][y-1] && state[x][y] == state[x+1][y] && state[x][y] == state[x][y-1]) continue;
                if (x > 0 && y < GRID-1) if (state[x][y] != state[x-1][y+1] && state[x][y] == state[x-1][y] && state[x][y] == state[x][y+1]) continue;
                if (x < GRID-1 && y < GRID-1) if (state[x][y] != state[x+1][y+1] && state[x][y] == state[x+1][y] && state[x][y] == state[x][y+1]) continue;

                state[x][y] ^= 1;
                bv += diff_score;
                edges += diff_edges;

                // update cache_ms
                if (x > 0 && y > 0) cache_ms[x-1][y-1] ^= 1 << 8;
                if (y > 0) cache_ms[x][y-1] ^= 1 << 7;
                if (x < GRID-1 && y > 0) cache_ms[x+1][y-1] ^= 1 << 6;
                if (x > 0) cache_ms[x-1][y] ^= 1 << 5;
                cache_ms[x][y] ^= 1 << 4;
                if (x < GRID-1) cache_ms[x+1][y] ^= 1 << 3;
                if (x > 0 && y < GRID-1) cache_ms[x-1][y+1] ^= 1 << 2;
                if (y < GRID-1) cache_ms[x][y+1] ^= 1 << 1;
                if (x < GRID-1 && y < GRID-1) cache_ms[x+1][y+1] ^= 1 << 0;

                if (bv > xv && edges <= MAX_EDGES && save_solution) {
                    xv = bv;
                    REP(i, GRID) REP(j, GRID) best_state[i][j] = state[i][j];
                }
                // if (score > best1) {
                //     best1 = score;
                //     DB(step, score, best1);
                // }
                // if (blocks > best) {
                //     best = blocks;
                //     bestGRID = blocks;
                //     DB(step, blocks, best);
                // }
            }
        }

        DB(GRID, elapsed(), step, xv, bv, edges);

        tsteps += step;


    };

    setup_grid(24);
    run_sa(0.8, 100, 0.01);
    upscale_state(2);
    setup_grid(48);
    run_sa(0.7, 20, 0.01);
    upscale_state(2);
    setup_grid(96);
    run_sa(0.2, 1, 0.01);
    upscale_state(2);
    setup_grid(192);
    run_sa(0.2, 1, 0.01, true);

    // DB(step, bv, xv, blocks, edges, comps);
    // DATA(bv);
    // DATA(blocks);
    // DATA(step);

    REP(i, GRID) REP(j, GRID) state[i][j] = best_state[i][j];

    // REP(y, GRID) {
    //     REP(x, GRID) cerr << (state[x][y] ? "#" : ".");
    //     cerr << endl;
    // }

    // draw grid

    VVI pixels(GRID+1, VI(GRID+1, 0));
    REP(i, GRID) REP(j, GRID) if (state[i][j]) REP(k, 2) REP(l, 2) pixels[i+k][j+l] = 1;

    VPII output;
    int ox = -1, oy = -1;
    REP(i, GRID+1) REP(j, GRID+1) if (pixels[i][j]) if (ox == -1) ox = i, oy = j;

    assert(ox != -1);

    int curx = ox, cury = oy;
    int curd = 0;

    do {
        output.PB({xlines[curx], ylines[cury]});
        int turnd = (curd + 3) % 4;
        int turnx = curx + dx[turnd];
        int turny = cury + dy[turnd];
        if (turnx >= 0 && turnx <= GRID && turny >= 0 && turny <= GRID && pixels[turnx][turny]) {
            curd = turnd;
            curx = turnx;
            cury = turny;
            continue;
        }

        int nextx = curx + dx[curd];
        int nexty = cury + dy[curd];
        if (nextx >= 0 && nextx <= GRID && nexty >= 0 && nexty <= GRID && pixels[nextx][nexty]) {
            curx = nextx;
            cury = nexty;
            continue;
        }

        int cornerd = (curd + 1) % 4;
        int cornerx = curx + dx[cornerd];
        int cornery = cury + dy[cornerd];
        assert(cornerx >= 0 && cornerx <= GRID && cornery >= 0 && cornery <= GRID && pixels[cornerx][cornery]);
        curd = cornerd;
        curx = cornerx;
        cury = cornery;
    } while (curx != ox || cury != oy);

    int len = 0;
    REP(i, output.S) {
        PII p0 = output[i];
        PII p1 = output[(i+1) % output.S];
        len += abs(p0.X - p1.X) + abs(p0.Y - p1.Y);
    }
    DATA(len);

    if (len > 4 * MAX_COORD) {
        output.clear();
        output.PB({0, 0});
        output.PB({MAX_COORD, 0});
        output.PB({MAX_COORD, MAX_COORD});
        output.PB({0, MAX_COORD});
    }

    double time = elapsed();
    DATA(time);
    DATA(tsteps);
    int vert = output.S;
    DATA(vert);
    DATA(xv);

    cout << output.S << endl;
    for (PII &p : output) cout << p.X << " " << p.Y << endl;
	
	return 0;
}
