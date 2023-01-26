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
 
struct FastRNG {
	unsigned long x=123456789, y=362436069, z=521288629;
    FastRNG() { }
	unsigned int rand() { unsigned long t; x ^= x << 16; x ^= x >> 5; x ^= x << 1; t = x; x = y; y = z; z = t ^ x ^ y; return z;}
    INLINE int next() {return rand(); }
    INLINE int next(int x) {return rand() % x; }
    INLINE int next(int a, int b) {return a + (rand() % (b - a)); }
    INLINE double next_double() {return (rand() + 0.5) * (1.0 / 4294967296.0); }
    INLINE double next_double(double a, double b) {return a + next_double() * (b - a); }
};

 
static FastRNG rng;
 
double start_time = get_time();

// SOLUTION

const bool AUTO_WIDTH = false;
const double TIME_SCALE = 1.0;
const double TIME_CUTOFF_START = 8.5 * TIME_SCALE;
const double TIME_CUTOFF_END = 9.75 * TIME_SCALE;

const int MAX_N = 30;
const int MAX_C = 4;

int N, C, P;

int N2aligned;

INLINE PII P2D(int p) {return MP(p / N, p % N);}
INLINE int P1D(int r, int c) {return r * N + c;}

// Vertical Beam Search

#define BYTE signed char

BYTE g0[MAX_N*MAX_N];
BYTE g0rot[4][MAX_N*MAX_N];

int column_pen[MAX_C][MAX_N];

const int COL_ADJ = 1;

int zobrist[MAX_N*MAX_N][64];
int tile_zobrist[64];

double hpen[64];

struct VState {
    BYTE g[MAX_N*MAX_N];
    BYTE edges[MAX_N*MAX_N];
    short fu_p[MAX_N*MAX_N];

    int comp[MAX_C];
    int pos[MAX_N];
    int n_bad;
    int n_zeroedge;
    int n_oneedge;
    int tile;
    int pen;
    int moves;
    int xscore;
    int xscore2;
    int hash;
    int rotation;

    void init() {
        ZERO(g);
        ZERO(edges);
        ZERO(fu_p);
        ZERO(comp);
        ZERO(pos);
        n_bad = 0;
        n_zeroedge = 0;
        n_oneedge = 0;
        tile = 0;
        pen = 0;
        moves = 0;
        xscore = 0;
        xscore2 = 0;
        hash = 0;
        rotation = 0;
    }

    void copy(VState &s) {  
        memcpy(g, s.g, sizeof(BYTE)*N2aligned);
        memcpy(edges, s.edges, sizeof(BYTE)*N2aligned);
        memcpy(fu_p, s.fu_p, sizeof(short)*N2aligned);

        memcpy(comp, s.comp, sizeof(comp));
        memcpy(pos, s.pos, sizeof(pos));
        n_bad = s.n_bad;
        n_zeroedge = s.n_zeroedge;
        n_oneedge = s.n_oneedge;
        tile = s.tile;
        pen = s.pen;
        moves = s.moves;
        xscore = s.xscore;
        xscore2 = s.xscore2;
        hash = s.hash;
        rotation = s.rotation;
    }


    INLINE void fu_make(int x) {
        fu_p[x] = -1;
    }

    int fu_find(int x) {
        if (fu_p[x] < 0) return x;
        return fu_p[x] = fu_find(fu_p[x]);
    }

    INLINE bool fu_join(int a, int b) {
        int fa = fu_find(a);
        int fb = fu_find(b);
        if (fa == fb) return false;

        if (fu_p[fa] > fu_p[fb]) swap(fa, fb);

        if (edges[fa] == 1) n_oneedge--;
        if (edges[fb] == 1) n_oneedge--;
        edges[fa] += edges[fb];

        fu_p[fa]--; // faster?
        // fu_p[fa] += fu_p[fb];
        fu_p[fb] = fa;
        return true;
    }

    INLINE void set_edges(int p, int e) {
        edges[p] = e;
        if (e == 0) n_zeroedge++;
        if (e == 1) n_oneedge++;
    }

    INLINE void update_edges(int p, int x) {
        if (edges[p] == 1) {
            n_oneedge--;
            edges[p] = 0;
            n_zeroedge++;
        } else {
            edges[p] += x;
            if (edges[p] == 0) n_zeroedge++;
            if (edges[p] == 1) n_oneedge++;
        }
    }

    void sim(int move) {
        assert(pos[move] < N);

        int gpos = P1D(pos[move], move);

        if (C > 1) {
            if (move > 0   && !fu_p[gpos-1] && g[gpos-1] == 0 && (tile & 8) && (move > 1   && (g[gpos-2] & 2) && (g[gpos-2] & 0xF0) != (tile & 0xF0) || pos[move] > 0 && (g[gpos-N-1] & 4) && (g[gpos-N-1] & 0xF0) != (tile & 0xF0)))
                n_bad++, fu_p[gpos-1] = 1;
            if (move < N-1 && !fu_p[gpos+1] && g[gpos+1] == 0 && (tile & 2) && (move < N-2 && (g[gpos+2] & 8) && (g[gpos+2] & 0xF0) != (tile & 0xF0) || pos[move] > 0 && (g[gpos-N+1] & 4) && (g[gpos-N+1] & 0xF0) != (tile & 0xF0)))
                n_bad++, fu_p[gpos+1] = 1;
            if (pos[move] < N-1 && !fu_p[gpos+N] && (tile & 4) && (move > 0 && (g[gpos+N-1] & 2) && (g[gpos+N-1] & 0xF0) != (tile & 0xF0) || move < N-1 && (g[gpos+N+1] & 8) && (g[gpos+N+1] & 0xF0) != (tile & 0xF0)))
                n_bad++, fu_p[gpos+N] = 1;
        }
        if (fu_p[gpos]) n_bad--, fu_p[gpos] = 0;

        fu_make(gpos);
        set_edges(gpos, ((tile & 1) > 0) + ((tile & 2) > 0) + ((tile & 4) > 0) + ((tile & 8) > 0));

        xscore += move > 0   && (g[gpos-1] & 0xF0) != (tile & 0xF0);
        xscore += move < N-1 && (g[gpos+1] & 0xF0) != (tile & 0xF0);
        xscore += pos[move]  && (g[gpos-N] & 0xF0) != (tile & 0xF0);

        xscore2 += move > 0   && pos[move] && (g[gpos-N-1] & 0xF0) != (tile & 0xF0);
        xscore2 += move < N-1 && pos[move] && (g[gpos-N+1] & 0xF0) != (tile & 0xF0);
        xscore2 += move > 0   && pos[move] < N-1 && (g[gpos+N-1] & 0xF0) != (tile & 0xF0);
        xscore2 += move < N-1 && pos[move] < N-1 && (g[gpos+N+1] & 0xF0) != (tile & 0xF0);

        comp[(tile>>4)-1]++;
        int fpos = gpos;

        if (move == 0) {
            if (tile & 8) pen++, update_edges(fpos, -1);
        } else if (g[gpos-1]) {
            if ((g[gpos-1] & 2) && (tile & 8) && (g[gpos-1] & 0xF0) == (tile & 0xF0)) {
                if (fu_join(gpos-1, gpos)) comp[(tile>>4)-1]--;
                fpos = fu_find(gpos);
                update_edges(fpos, -2);
            } else if ((tile & 8) || (g[gpos-1] & 2)) {
                if (tile & 8) update_edges(fu_find(gpos), -1);
                if (g[gpos-1] & 2) update_edges(fu_find(gpos-1), -1);
                pen++;
            }
        }

        if (move == N-1) {
            if (tile & 2) pen++, update_edges(fpos, -1);
        } else if (g[gpos+1]) {
            if ((g[gpos+1] & 8) && (tile & 2) && (g[gpos+1] & 0xF0) == (tile & 0xF0)) {
                if (fu_join(gpos+1, gpos)) comp[(tile>>4)-1]--;
                fpos = fu_find(gpos);
                update_edges(fpos, -2);
            } else if ((tile & 2) || (g[gpos+1] & 8)) {
                if (tile & 2) update_edges(fpos, -1);
                if (g[gpos+1] & 8) update_edges(fu_find(gpos+1), -1);
                pen++;
            }
        }

        if (pos[move] == 0) {
            if (tile & 1) pen++, update_edges(fpos, -1);
        } else {
            if ((g[gpos-N] & 4) && (tile & 1) && (g[gpos-N] & 0xF0) == (tile & 0xF0)) {
                if (fu_join(gpos-N, gpos)) comp[(tile>>4)-1]--;
                fpos = fu_find(gpos);
                update_edges(fpos, -2);
            } else if ((tile & 1) || (g[gpos-N] & 4)) {
                if (tile & 1) update_edges(fpos, -1);
                if (g[gpos-N] & 4) update_edges(fu_find(gpos-N), -1);
                pen++;
            }
        }

        if (pos[move] == N-1 && (tile & 4)) pen++, update_edges(fpos, -1);

        {
            if (move > 0   && !fu_p[gpos-1] && g[gpos-1] == 0 && (tile & 8) && (move > 1   && (g[gpos-2] & 2) && (g[gpos-2] & 0xF0) == (tile & 0xF0) && fu_find(gpos-2) == fpos || pos[move] > 0 && (g[gpos-N-1] & 4) && (g[gpos-N-1] & 0xF0) == (tile & 0xF0) && fu_find(gpos-N-1) == fpos))
                n_bad++, fu_p[gpos-1] = 1;
            if (move < N-1 && !fu_p[gpos+1] && g[gpos+1] == 0 && (tile & 2) && (move < N-2 && (g[gpos+2] & 8) && (g[gpos+2] & 0xF0) == (tile & 0xF0) && fu_find(gpos+2) == fpos || pos[move] > 0 && (g[gpos-N+1] & 4) && (g[gpos-N+1] & 0xF0) == (tile & 0xF0) && fu_find(gpos-N+1) == fpos))
                n_bad++, fu_p[gpos+1] = 1;
            if (pos[move] < N-1 && !fu_p[gpos+N] && (tile & 4) && (move > 0 && (g[gpos+N-1] & 2) && (g[gpos+N-1] & 0xF0) == (tile & 0xF0) && fu_find(gpos+N-1) == fpos || move < N-1 && (g[gpos+N+1] & 8) && (g[gpos+N+1] & 0xF0) == (tile & 0xF0) && fu_find(gpos+N+1) == fpos))
                n_bad++, fu_p[gpos+N] = 1;
        }



        hash ^= tile_zobrist[tile-16];
        hash ^= zobrist[gpos][tile-16];

        g[gpos] = tile;
        tile = g0rot[rotation][gpos];

        hash ^= tile_zobrist[tile-16];
        pos[move]++;
        moves++;
    }

    double sim_score(int move) {
        assert(pos[move] < N);

        int gpos = P1D(pos[move], move);

        int old_pen = pen;

        static int fus[3];
        int n_fus = 0;
        int t_edge = ((tile & 1) > 0) + ((tile & 2) > 0) + ((tile & 4) > 0) + ((tile & 8) > 0);
        int blocked = 0;
        bool conn_top = false;

        if (move == 0) {
            if (tile & 8) pen++, t_edge--;
        } else if (g[gpos-1]) {
            if ((g[gpos-1] & 2) && (tile & 8) && (g[gpos-1] & 0xF0) == (tile & 0xF0)) {
                fus[n_fus++] = fu_find(gpos-1);
                t_edge -= 2;
            } else if ((tile & 8) || (g[gpos-1] & 2)) {
                if (tile & 8) t_edge--;
                if ((g[gpos-1] & 2) && edges[fu_find(gpos-1)] == 1) blocked++;
                pen++;
            }
        }

        if (move == N-1) {
            if (tile & 2) pen++, t_edge--;
        } else if (g[gpos+1]) {
            if ((g[gpos+1] & 8) && (tile & 2) && (g[gpos+1] & 0xF0) == (tile & 0xF0)) {
                fus[n_fus++] = fu_find(gpos+1);
                t_edge -= 2;
            } else if ((tile & 2) || (g[gpos+1] & 8)) {
                if (tile & 2) t_edge--;
                if ((g[gpos+1] & 8) && edges[fu_find(gpos+1)] == 1) blocked++;
                pen++;
            }
        }

        if (pos[move] == 0) {
            if (tile & 1) pen++, t_edge--;
        } else {
            if ((g[gpos-N] & 4) && (tile & 1) && (g[gpos-N] & 0xF0) == (tile & 0xF0)) {
                conn_top = true;
                fus[n_fus++] = fu_find(gpos-N);
                t_edge -= 2;
            } else if ((tile & 1) || (g[gpos-N] & 4)) {
                if (tile & 1) t_edge--;
                if ((g[gpos-N] & 4) && edges[fu_find(gpos-N)] == 1) blocked++;
                pen++;
            }
        }

        if (pos[move] == N-1 && (tile & 4)) pen++, t_edge--;



        int old_comp = comp[(tile>>4)-1];
        int new_comp = comp[(tile>>4)-1];

        REP(i, n_fus) t_edge += edges[fu_find(fus[i])];
        
        int cur_n_bad = n_bad;

        new_comp++;
        if (n_fus) new_comp--;
        if (n_fus == 2) {
            if (fus[0] != fus[1]) new_comp--;
            else cur_n_bad++;
        } else if (n_fus == 3) {
            if (fus[0] != fus[1] && fus[0] != fus[2] && fus[1] != fus[2]) new_comp -= 2;
            else if (fus[0] != fus[1] || fus[0] != fus[2] || fus[1] != fus[2]) new_comp -= 1, cur_n_bad++;
            else cur_n_bad += 2;
        }
        comp[(tile>>4)-1] = new_comp;

        double rv = 0;
        const int comp_w[] = {0, 0, 4, 8, 13};
        REP(i, C) rv += comp_w[C] * comp[i];
        rv += pen * (20 + P);

        pen = old_pen;
        comp[(tile>>4)-1] = old_comp;

        int fpos = conn_top ? fu_find(gpos-N) : -1;
        if (move > 0   && !fu_p[gpos-1] && g[gpos-1] == 0 && (tile & 8) && (move > 1   && (g[gpos-2] & 2) && (g[gpos-2] & 0xF0) != (tile & 0xF0) || pos[move] > 0 && (g[gpos-N-1] & 4) && (g[gpos-N-1] & 0xF0) != (tile & 0xF0)))
            cur_n_bad++;
        else if (move > 0   && !fu_p[gpos-1] && g[gpos-1] == 0 && (tile & 8) && (move > 1   && (g[gpos-2] & 2) && (g[gpos-2] & 0xF0) == (tile & 0xF0) && fu_find(gpos-2) == fpos || pos[move] > 0 && (g[gpos-N-1] & 4) && (g[gpos-N-1] & 0xF0) == (tile & 0xF0) && fu_find(gpos-N-1) == fpos))
            cur_n_bad++;

        if (move < N-1 && !fu_p[gpos+1] && g[gpos+1] == 0 && (tile & 2) && (move < N-2 && (g[gpos+2] & 8) && (g[gpos+2] & 0xF0) != (tile & 0xF0) || pos[move] > 0 && (g[gpos-N+1] & 4) && (g[gpos-N+1] & 0xF0) != (tile & 0xF0)))
            cur_n_bad++;
        else if (move < N-1 && !fu_p[gpos+1] && g[gpos+1] == 0 && (tile & 2) && (move < N-2 && (g[gpos+2] & 8) && (g[gpos+2] & 0xF0) == (tile & 0xF0) && fu_find(gpos+2) == fpos || pos[move] > 0 && (g[gpos-N+1] & 4) && (g[gpos-N+1] & 0xF0) == (tile & 0xF0) && fu_find(gpos-N+1) == fpos))
            cur_n_bad++;

        if (pos[move] < N-1 && !fu_p[gpos+N] && (tile & 4) && (move > 0 && (g[gpos+N-1] & 2) && (g[gpos+N-1] & 0xF0) != (tile & 0xF0) || move < N-1 && (g[gpos+N+1] & 8) && (g[gpos+N+1] & 0xF0) != (tile & 0xF0)))
            cur_n_bad++;
        else if (pos[move] < N-1 && !fu_p[gpos+N] && (tile & 4) && (move > 0 && (g[gpos+N-1] & 2) && (g[gpos+N-1] & 0xF0) == (tile & 0xF0) && fu_find(gpos+N-1) == fpos || move < N-1 && (g[gpos+N+1] & 8) && (g[gpos+N+1] & 0xF0) == (tile & 0xF0) && fu_find(gpos+N+1) == fpos))
            cur_n_bad++;
            
        if (fu_p[gpos])
            cur_n_bad--;
        rv += 35 * cur_n_bad;

        pos[move]++;
        double h_score = 0;
        REP(i, N-1) h_score += hpen[32 + pos[i] - pos[i+1]];    
        rv += h_score * (C <= 2 ? 2.0 : 1.2);
        pos[move]--;

        rv += n_zeroedge * 17;
        rv += (t_edge == 0) * 17;
        rv += blocked * 17;
        rv += n_oneedge * 3;

        rv += xscore * 3.0;
        rv += (move > 0   && (g[gpos-1] & 0xF0) != (tile & 0xF0)) * 3;
        rv += (move < N-1 && (g[gpos+1] & 0xF0) != (tile & 0xF0)) * 3;
        rv += (pos[move]  && (g[gpos-N] & 0xF0) != (tile & 0xF0)) * 3;

        if (C==1) rv += xscore2 * 1;

        rv += rng.next_double();

        return rv;
    }

    double score() {
        double rv = moves;
        REP(i, C) rv += comp[i] * comp[i];
        rv += pen * P;
        return rv;
    }

    void show_stats() {
        cerr << "Moves: " << moves << " Mismatches: " << pen;
        REP(i, C) cerr << " Col #" << i << ": " << comp[i];
        cerr << endl;
    }
};

INLINE int rotate_tile(int tile) {
    return (tile & 0xF0) | ((tile & 7) << 1) | ((tile & 8) >> 3);
}

void rotate_all() {
    int tmp_g0[MAX_N*MAX_N];
    REP(r, N) REP(c, N) tmp_g0[c*N+N-1-r] = rotate_tile(g0[r*N+c]);
    REP(i, N*N) g0[i] = tmp_g0[i];
}

const int MAX_STATES = 120000;
const int MAX_NEW_STATES = MAX_STATES * MAX_N;
VState beam[2][MAX_STATES];
BYTE sol[2][MAX_STATES][MAX_N*MAX_N];

PII new_states[MAX_NEW_STATES];
double bscores[MAX_NEW_STATES];
int order[MAX_NEW_STATES];

const int MAX_HASHES = 1<<24;
short hashes[MAX_HASHES];

double lv_timestamps[MAX_N*MAX_N];
int lv_widths[MAX_N*MAX_N];

int main(int argc, char **argv) {
    ios_base::sync_with_stdio(false);

    VState s_start;
    s_start.init();
    cin >> N >> C >> P;

    REP(r, N) REP(c, N) {
        int a, b;
        cin >> a >> b;
        g0[P1D(r,c)] = a + b * 16 + 16;
    }

    N2aligned = N*N;
    while ((N2aligned & 15) && N2aligned < N*N) N2aligned++;

    int tile, tile_color;
    cin >> tile >> tile_color;
    s_start.tile = tile + tile_color * 16 + 16;

    cerr << "[DATA] N = " << N << endl;
    cerr << "[DATA] C = " << C << endl;
    cerr << "[DATA] P = " << P << endl;

    REP(i, N*N) REP(j, 64) zobrist[i][j] = rng.next() & (MAX_HASHES - 1);
    REP(i, 64) tile_zobrist[i] = rng.next() & (MAX_HASHES - 1);

    // beam search

    const int orig_width = 10000000 / N / N / sqrt(N) * TIME_SCALE;
    int max_width = orig_width;

    assert(max_width < MAX_STATES);

    DB(max_width);

    if (C == 1) {
        REP(i, 64) hpen[i] = abs(i - 32);
    } else if (C == 2) {
        REP(i, 64) hpen[i] = abs(i - 32) * sqrt(abs(i - 32));
    } else {
        REP(i, 64) hpen[i] = (i - 32) * (i - 32);
    }

    double xv = 1e100;
    VI bsol;
    int brot;
    VState bs;

    VState rot_start[4];
    REP(rot, 4) {
        rot_start[rot] = s_start;
        rot_start[rot].rotation = rot;
        rot_start[rot].hash = rng.next() & (MAX_HASHES - 1);
        REP(i, N*N) g0rot[rot][i] = g0[i];
        rotate_all();
        s_start.tile = rotate_tile(s_start.tile);
    }

    REP(i, 4) beam[0][i] = rot_start[i];
    MINUS(hashes);

    int n_cur_beam = 4;

    int autoscale_start = 0;
    REP(level, N*N) {
        if (AUTO_WIDTH && (level % N == 0 || level >= N*N-N)) DB(get_time()-start_time, level, n_cur_beam);

        lv_widths[level] = n_cur_beam + (level ? lv_widths[level-1] : 0);
        lv_timestamps[level] = get_time() - start_time;

        if (AUTO_WIDTH) {
            if (autoscale_start == 0 && n_cur_beam == orig_width) autoscale_start = level;
            if (autoscale_start && level > autoscale_start) {
                int max_avg = N / 2;
                int n_avg = 1;
                while (n_avg < max_avg && level - n_avg >= autoscale_start) n_avg++;
                double finish_time = (TIME_CUTOFF_START + (TIME_CUTOFF_END - TIME_CUTOFF_START) * (1.0 * level / (N*N)));
                int exp_evals = (int)((finish_time - lv_timestamps[level]) * (lv_widths[level-1] - lv_widths[level-1-n_avg]) / (lv_timestamps[level] - lv_timestamps[level-n_avg]));
                max_width = min(MAX_STATES, max(orig_width / 5, exp_evals / (N*N - level)));
            }
        }

        int n_new_states = 0;

        auto &cur_beam = beam[level & 1];
        auto &next_beam = beam[(level & 1) ^ 1];
        auto &cur_sol = sol[level & 1];
        auto &next_sol = sol[(level & 1) ^ 1];

        REP(i, n_cur_beam) {
            if (hashes[cur_beam[i].hash] == level) continue;
            hashes[cur_beam[i].hash] = level;

            REP(move, N) {
                if (cur_beam[i].pos[move] == N) continue;

                double av = cur_beam[i].sim_score(move);
                bscores[n_new_states] = av;
                new_states[n_new_states++] = MP(i, move);
            }
        }

        REP(i, n_new_states) order[i] = i;
        if (n_new_states > max_width) {
            nth_element(order, order + max_width, order + n_new_states, [](int a, int b) { return bscores[a] < bscores[b]; });
            n_new_states = max_width;
        }

        int level_aligned = level;
        while (level_aligned & 15) level_aligned++;

        REP(i, n_new_states) {
            int id = new_states[order[i]].X;
            int move = new_states[order[i]].Y;
            next_beam[i] = cur_beam[id];
            next_beam[i].sim(move);
            memcpy(next_sol[i], cur_sol[id], sizeof(BYTE) * level);
            next_sol[i][level] = move;
        }

        if (level == N*N-1) {
            REP(i, n_new_states) {
                double av = next_beam[i].score();
                if (av < xv) {
                    xv = av;
                    bsol = VI(next_sol[i], next_sol[i] + N*N);
                    brot = next_beam[i].rotation;
                    bs = next_beam[i];
                }
            }
        }

        n_cur_beam = n_new_states;
    }

    DB(get_time() - start_time);

    cout << bsol.S << endl;
    if (brot == 0) {
        for (int &x : bsol) cout << "U " << x << endl;
    } else if (brot == 1) {
        for (int &x : bsol) cout << "L " << (N-1-x) << endl;
    } else if (brot == 2) {
        for (int &x : bsol) cout << "D " << (N-1-x) << endl;
    } else if (brot == 3) {
        for (int &x : bsol) cout << "R " << x << endl;
    }

	return 0;
}
