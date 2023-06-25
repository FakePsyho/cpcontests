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
 
static RNG rng(3);

// SOLUTION

double timer1 = 0;
double timer2 = 0;
double timer3 = 0;
double timer4 = 0;

// #define USE_TIMERS

#ifdef USE_TIMERS
#define TIMER_START(t) t -= get_time()
#define TIMER_STOP(t)  t += get_time()
#else
#define TIMER_START(t) ;
#define TIMER_STOP(t)  ;
#endif


const double TIME_SCALE = 1.0;

#ifdef VM
const double MACHINE_SCALE = 0.5;
#else
const double MACHINE_SCALE = 1.0;
#endif

const double TIME_LIMIT_P1 = 1.2 * TIME_SCALE * MACHINE_SCALE;
const double TIME_LIMIT_P2 = 8.5 * TIME_SCALE * MACHINE_SCALE;
const double TIME_LIMIT_P3 = 9.8 * TIME_SCALE * MACHINE_SCALE;

const int MAXN = 32;
const int MAXN2 = MAXN * MAXN;

const int dd[] = {-MAXN, -1, 1, MAXN};

int P1D(int r, int c) {return r * MAXN + c + MAXN + 1;}
PII P2D(int p) {return MP(p / MAXN - 1, (p - 1) & (MAXN - 1));}

int N, C, K, P;
char og[MAXN2];
char g[MAXN2];

short dist[MAXN2][MAXN2];
short closest[MAXN2][MAXN2];
short dmove[MAXN2][MAXN2];
short distall[MAXN2][MAXN2];
short closestall[MAXN2][MAXN2];
short dmoveall[MAXN2][MAXN2];

int last_middle = 0;
VI build_distgraph(int p0 = -1) {
  int vs[MAXN2] = {0};

  int bp = -1;
  int bv = 0;
  if (p0 != -1) {
    bp = p0;
  } else {
    REP(r, N) REP(c, N) if (og[P1D(r,c)]) {
      int av = min(r, N-1-r) * min(c, N-1-c);
      // int av = r+c;
      if (av > bv) {
        bv = av;
        bp = P1D(r, c);
      }
    }
  }
  last_middle = bp;

  VI order;
  queue<PII> q;
  map<int, int> m;
  q.push(MP(bp, 0));
  ZERO(vs);

  DB(bp, P2D(bp));

  while (!q.empty()) {
    PII pair = q.front(); q.pop();
    int p = pair.X;
    int dist = pair.Y;
    if (vs[p]) continue;
    vs[p] = 1;
    m[p] = dist;
    order.PB(p);

    REP(d, 4) {
      int np = p + dd[d];
      if (og[np] == 0) continue;
      q.push(MP(np, dist+1));
    }
  }
  reverse(ALL(order));

  FOR(r0, -1, N+1) FOR(c0, -1, N+1) FOR(r1, -1, N+1) FOR(c1, -1, N+1) dist[P1D(r0,c0)][P1D(r1,c1)] = 1<<14;
  FOR(r0, -1, N+1) FOR(c0, -1, N+1) FOR(r1, -1, N+1) FOR(c1, -1, N+1) dmove[P1D(r0,c0)][P1D(r1,c1)] = -2;
  ZERO(vs);
  for (int p0 : order) {
    queue<PII> q;
    q.push(MP(p0, -1));
    int cnt = 0;
    while (!q.empty()) {
      PII pair = q.front(); q.pop();
      int p = pair.X;
      int prev = pair.Y;
      if (dist[p0][p] != (1<<14)) continue;
      if (vs[p]) continue;
      dist[p0][p] = prev != -1 ? dist[p0][prev] + 1 : 0;
      dmove[p0][p] = prev;
      closest[p0][cnt++] = p;
      if (og[p-MAXN]) q.push(MP(p-MAXN, p));
      if (og[p-1   ]) q.push(MP(p-1   , p));
      if (og[p+1   ]) q.push(MP(p+1   , p));
      if (og[p+MAXN]) q.push(MP(p+MAXN, p));
    }
    vs[p0] = 1;
  }

  FOR(r0, -1, N+1) FOR(c0, -1, N+1) FOR(r1, -1, N+1) FOR(c1, -1, N+1) distall[P1D(r0,c0)][P1D(r1,c1)] = 1<<14;
  FOR(r0, -1, N+1) FOR(c0, -1, N+1) FOR(r1, -1, N+1) FOR(c1, -1, N+1) dmoveall[P1D(r0,c0)][P1D(r1,c1)] = -2;
  for (int p0 : order) {
    queue<PII> q;
    q.push(MP(p0, -1));
    int cnt = 0;
    while (!q.empty()) {
      PII pair = q.front(); q.pop();
      int p = pair.X;
      int prev = pair.Y;
      if (distall[p0][p] != (1<<14)) continue;
      if (m[p0] < m[p]) continue;
      distall[p0][p] = prev != -1 ? distall[p0][prev] + 1 : 0;
      dmoveall[p0][p] = prev;
      closestall[p0][cnt++] = p;
      if (og[p-MAXN]) q.push(MP(p-MAXN, p));
      if (og[p-1   ]) q.push(MP(p-1   , p));
      if (og[p+1   ]) q.push(MP(p+1   , p));
      if (og[p+MAXN]) q.push(MP(p+MAXN, p));
    }
  }

  return order;
}


struct Components {
  VI cnt;
  int mx;
  double score;

  int last_comp_change[MAXN2] = {0};
  int last_comp_change_no = 0;

  Components() {
    cnt = VI(N*N, 0);
    mx = 0;
    score = 0;
  }

  void reset() {
    REP(i, mx+1) cnt[i] = 0;
    mx = 0;
    score = 0;
  }

  void add(int x) {
    cnt[x]++;
    if (x > mx) mx = x;
    score += comp_value(x);
  }

  void rem(int x) {
    cnt[x]--;
    while (cnt[mx] == 0 && mx) mx--;
    score -= comp_value(x);
  }

  void xadd(int x) {
    last_comp_change[last_comp_change_no++] = x;
    add(x);
  }

  void xrem(int x) {
    last_comp_change[last_comp_change_no++] = -x;
    rem(x);
  }

  double comp_value(int x) {
    return x == K ? 1 : x > K ? (K-x) : -.1;
  }

  void clear_history() {
    last_comp_change_no = 0;
  }

  void restore() {
    REP(i, last_comp_change_no) {
      int v = last_comp_change[i];
      if (v > 0) {
        rem(v);
      } else {
        add(-v);
      }
    }
  }
};

void count_components(Components &comps) {
  TIMER_START(timer1);
  comps.reset();
  int rv = 0;
  static char cc_vs[MAXN2] = {0};
  static short q[MAXN2];

  ZERO(cc_vs);
  
  REP(r0, N) REP(c0, N) {
    int p0 = P1D(r0, c0);
    if (g[p0] == 0 || cc_vs[p0]) continue;
    int col = g[p0];
    q[0] = p0;
    int qpos = 1;
    int tot = 0;
    while (qpos--) {
      int p = q[qpos];
      if (cc_vs[p]) continue;
      cc_vs[p] = 1;
      tot++;
      if (g[p-MAXN] == col) q[qpos++] = p-MAXN;
      if (g[p-1   ] == col) q[qpos++] = p-1   ;
      if (g[p+1   ] == col) q[qpos++] = p+1   ;
      if (g[p+MAXN] == col) q[qpos++] = p+MAXN;
    }
    comps.add(tot);
  }
  TIMER_STOP(timer1);

}

void dyn_count_components(Components &comps, int p1, int p2) {
  TIMER_START(timer2);
  assert(p1 != p2);
  static int starts[16];
  int starts_cnt = 2;
  starts[0] = p1;
  starts[1] = p2;

  comps.clear_history();

  REP(d, 4) {
    int np1 = p1 + dd[d];
    if (g[np1] == g[p2]) starts[starts_cnt++] = np1;
    int np2 = p2 + dd[d];
    if (g[np2] == g[p1]) starts[starts_cnt++] = np2;
  }

  int removed_cnt = 0;

  static short removed[MAXN2];
  static short q[MAXN2];

  int rv = 0;
  REP(i, starts_cnt) {
    int p0 = starts[i];
    if (g[p0] < 0) continue;
    int col = g[p0];
    q[0] = p0;
    g[p0] = -col;
    int qpos = 1;
    int tot = 0;
    while (qpos--) {
      int p = q[qpos];
      tot++;
      removed[removed_cnt++] = p;
      if (g[p-MAXN] == col) {q[qpos++] = p-MAXN; g[p-MAXN] = -col;}
      if (g[p-1   ] == col) {q[qpos++] = p-1   ; g[p-1   ] = -col;}
      if (g[p+1   ] == col) {q[qpos++] = p+1   ; g[p+1   ] = -col;}
      if (g[p+MAXN] == col) {q[qpos++] = p+MAXN; g[p+MAXN] = -col;}
    }
    comps.xrem(tot);
  }

  swap(g[p1], g[p2]);

  REP(i, removed_cnt) {
    int p0 = removed[i];
    if (g[p0] > 0) continue;
    int col = g[p0];
    q[0] = p0;
    g[p0] = -col;
    int qpos = 1;
    int tot = 0;
    while (qpos--) {
      int p = q[qpos];
      tot++;
      if (g[p-MAXN] == col) {q[qpos++] = p-MAXN; g[p-MAXN] = -col;}
      if (g[p-1   ] == col) {q[qpos++] = p-1   ; g[p-1   ] = -col;}
      if (g[p+1   ] == col) {q[qpos++] = p+1   ; g[p+1   ] = -col;}
      if (g[p+MAXN] == col) {q[qpos++] = p+MAXN; g[p+MAXN] = -col;}
    }
    comps.xadd(tot);
  }

  swap(g[p1], g[p2]);
  TIMER_STOP(timer2);
}

void dyn_count_components(Components &comps, int p1, int p2, int p3) {
  TIMER_START(timer2);
  static int starts[24];
  int starts_cnt = 3;
  starts[0] = p1;
  starts[1] = p2;
  starts[2] = p3;

  assert(g[p1] != g[p2]);
  assert(g[p2] != g[p3]);
  assert(g[p1] != g[p3]);

  comps.clear_history();

  REP(d, 4) {
    int np1 = p1 + dd[d];
    if (g[np1] == g[p3]) starts[starts_cnt++] = np1;
    int np2 = p2 + dd[d];
    if (g[np2] == g[p1]) starts[starts_cnt++] = np2;
    int np3 = p3 + dd[d];
    if (g[np3] == g[p2]) starts[starts_cnt++] = np3;
  }

  int removed_cnt = 0;

  static short removed[MAXN2];
  static short q[MAXN2];

  int rv = 0;
  REP(i, starts_cnt) {
    int p0 = starts[i];
    if (g[p0] < 0) continue;
    int col = g[p0];
    q[0] = p0;
    g[p0] = -col;
    int qpos = 1;
    int tot = 0;
    while (qpos--) {
      int p = q[qpos];
      tot++;
      removed[removed_cnt++] = p;
      if (g[p-MAXN] == col) {q[qpos++] = p-MAXN; g[p-MAXN] = -col;}
      if (g[p-1   ] == col) {q[qpos++] = p-1   ; g[p-1   ] = -col;}
      if (g[p+1   ] == col) {q[qpos++] = p+1   ; g[p+1   ] = -col;}
      if (g[p+MAXN] == col) {q[qpos++] = p+MAXN; g[p+MAXN] = -col;}
    }
    comps.xrem(tot);
  }

  swap(g[p1], g[p2]);
  swap(g[p1], g[p3]);

  REP(i, removed_cnt) {
    int p0 = removed[i];
    if (g[p0] > 0) continue;
    int col = g[p0];
    q[0] = p0;
    g[p0] = -col;
    int qpos = 1;
    int tot = 0;
    while (qpos--) {
      int p = q[qpos];
      tot++;
      if (g[p-MAXN] == col) {q[qpos++] = p-MAXN; g[p-MAXN] = -col;}
      if (g[p-1   ] == col) {q[qpos++] = p-1   ; g[p-1   ] = -col;}
      if (g[p+1   ] == col) {q[qpos++] = p+1   ; g[p+1   ] = -col;}
      if (g[p+MAXN] == col) {q[qpos++] = p+MAXN; g[p+MAXN] = -col;}
    }
    comps.xadd(tot);
  }

  swap(g[p1], g[p3]);
  swap(g[p1], g[p2]);

  TIMER_STOP(timer2);
}

int find_neighbor(int p0) {
  static short updated[MAXN2];
  static short nb[MAXN2];
  static short q[MAXN2];
  int n_updated = 1;
  int n_nb = 0;
  int col = g[p0];

  q[0] = p0;
  g[p0] = -col;
  updated[0] = p0;
  int qpos = 1;

  while (qpos--) {
    int p = q[qpos];
    if (g[p-MAXN] > 0) {if (g[p-MAXN] == col) {q[qpos++] = p-MAXN;} else {nb[n_nb++] = p-MAXN;} g[p-MAXN] = -g[p-MAXN]; updated[n_updated++] = p-MAXN; }
    if (g[p-1   ] > 0) {if (g[p-1   ] == col) {q[qpos++] = p-1   ;} else {nb[n_nb++] = p-1   ;} g[p-1   ] = -g[p-1   ]; updated[n_updated++] = p-1   ; }
    if (g[p+1   ] > 0) {if (g[p+1   ] == col) {q[qpos++] = p+1   ;} else {nb[n_nb++] = p+1   ;} g[p+1   ] = -g[p+1   ]; updated[n_updated++] = p+1   ; }
    if (g[p+MAXN] > 0) {if (g[p+MAXN] == col) {q[qpos++] = p+MAXN;} else {nb[n_nb++] = p+MAXN;} g[p+MAXN] = -g[p+MAXN]; updated[n_updated++] = p+MAXN; }
  }

  REP(i, n_updated) g[updated[i]] = -g[updated[i]];

  int rv = nb[rng.next(n_nb)];
  assert(g[rv] && g[rv] != col);
  return rv;
}

VPII last_sol;

int simulate_swaps(VI &order, const bool gen_sol = false) {
  TIMER_START(timer3);
  last_sol.clear();

  int turns = 0;

  static char cg[MAXN2];
  REP(r, N) REP(c, N) cg[P1D(r,c)] = og[P1D(r,c)];

  REP(j, order.S) {
    int p0 = order[j] % MAXN2;

    int bv = 1<<20;
    int p = -1;
    for (int i=0; ;i++) if (g[p0] == cg[closest[p0][i]]) {
      p = closest[p0][i];
      bv = dist[p0][p];
      break;
    }

    int moves = dist[p0][p];
    int np = dmove[p0][p];
    turns += moves;
    while (np != -1) {
      swap(cg[p], cg[np]);
      if (gen_sol) last_sol.PB(MP(p, np));
      p = np;
      np = dmove[p0][p];
    }
  }

  // for (int p : order) if (g[p] != cg[p]) {
  //   TIMER_STOP(timer3);
  //   return 1<<20;
  // }
  TIMER_STOP(timer3);

  return turns;
}

int simulate_swaps2(VI &order, const bool gen_sol = false) {
  TIMER_START(timer4);
  last_sol.clear();

  int turns = 0;

  static char cg[MAXN2];
  static int used[MAXN2] = {0};
  static int used_no = 1;
  used_no++;

  REP(r, N) REP(c, N) cg[P1D(r,c)] = og[P1D(r,c)];

  REP(j, order.S) {
    int p0 = order[j] % MAXN2;

    int bv = 1<<10;
    int p = -1;
    for (int i=0; i<MAXN2 ;i++) if (g[p0] == cg[closestall[p0][i]] && used[closestall[p0][i]] < used_no) {
      p = closestall[p0][i];
      bv = distall[p0][p];
      break;
    }

    if (p == -1) {
      TIMER_STOP(timer4);
      return 1<<21;
    }

    turns += distall[p0][p];
    int dp = order[j] / MAXN2;
    while (distall[p0][p]) {
      bool found = false;
      REP(d, 4) {
        int np = p + dd[(d+dp)&3];
        if (used[np] == used_no) continue;
        if (distall[p0][np] == distall[p0][p]-1) {
          if (gen_sol) last_sol.PB(MP(p, np));
          swap(cg[p], cg[np]);
          p = np;
          found = true;
          break;
        }
      }
      if (!found) {
        TIMER_STOP(timer4);
        return (1<<20) + (1<<10) - j;
      }
    }
    used[p0] = used_no;

    assert(p == p0);
    assert(g[p0] == cg[p0]);
  }

  for (int p : order) if (g[p%MAXN2] != cg[p%MAXN2]) {
    TIMER_STOP(timer4);
    return 1<<21;
  }
  // for (int p : order) assert(g[p] == cg[p]);
  TIMER_STOP(timer4);

  return turns;
}

int main(int argc, char **argv) {
  cin >> N >> C >> K >> P;
  REP(r, N) REP(c, N) {
    int x; cin >> x;
    og[P1D(r,c)] = x;
  }

  cerr << "[DATA] N = " << N << endl;
  cerr << "[DATA] C = " << C << endl;
  cerr << "[DATA] K = " << K << endl;
  cerr << "[DATA] P = " << P << endl;
  cerr << "[DATA] r = " << 1.0*P/K << endl;

  VI order = build_distgraph();

  int max_components = 0;
  VI colors(C + 1);
  REP(r, N) REP(c, N) colors[og[P1D(r,c)]]++;
  FOR(i, 1, C+1) max_components += colors[i] / K;

  DB(max_components);

  REP(r, N) REP(c, N) g[P1D(r,c)] = og[P1D(r,c)];

  VPII bsol;
  int bscore = -(1<<20);
  int bturns = -1;
  char bg[MAXN2];

  auto update_best = [&](int score, auto sim_function) {
    if (score <= bscore) return;
    int turns = sim_function(order, true);
    bscore = score;
    bsol = last_sol;
    bturns = turns;
    REP(r, N) REP(c, N) bg[P1D(r,c)] = g[P1D(r,c)];
  };

  Components comps;
  
  double time_passed = 0;

  const double p1_t0 = .01;
  const double p1_tn = .00008;

  double p1_t = p1_t0;

  count_components(comps);
  int components = comps.cnt[K];

  double bv = comps.score;

  int p1_step = 0;
  const double P1_STARTED = elapsed();

  while (true) {
    p1_step++;

    if ((p1_step & 1023) == 0) {
      time_passed = (elapsed() - P1_STARTED) / (TIME_LIMIT_P1 - P1_STARTED);
      p1_t = p1_t0 * pow(p1_tn / p1_t0, time_passed);
      if (time_passed > 1 || components == max_components) break;
    }

    int p1, p2;
    do {int r1 = rng.next(N); int c1 = rng.next(N); p1 = P1D(r1, c1);} while (g[p1] == 0);

    int type = rng.next(20) == 0;
    if (type == 0) {
      p2 = find_neighbor(p1);
    } else if (type == 1) {
      do {int r2 = rng.next(N); int c2 = rng.next(N); p2 = P1D(r2, c2);} while (g[p2] == 0 || g[p1] == g[p2]);
    }

    dyn_count_components(comps, p1, p2);
    int new_components = comps.cnt[K];
    double av = comps.score;
    if (av >= bv || rng.next_double() < exp((av - bv) / p1_t)) {
      bv = av;
      swap(g[p1], g[p2]);
      if (new_components > components) {
        int turns = simulate_swaps(order);
        int score = new_components * P - turns;
        update_best(score, simulate_swaps);
        if (K <= 3 || C == 2) break;
      }
      components = new_components;
    } else {
      comps.restore();
    }
  }

  REP(r, N) REP(c, N) g[P1D(r,c)] = bg[P1D(r,c)];

  count_components(comps);
  components = comps.cnt[K];

  DB(p1_step, bscore, components);

  const double p2_t0 = 4;
  const double p2_tn = .6;
  double p2_t = p2_t0;

  bv = bscore;

  const double P2_STARTED = elapsed();
  int p2_step = 0;

  while (true) {
    p2_step++;

    if ((p2_step & 1023) == 0) {
      time_passed = (elapsed() - P2_STARTED) / (TIME_LIMIT_P2 - P2_STARTED);
      // p2_t = p2_t0 * pow(p2_tn / p2_t0, time_passed);
      p2_t = p2_t0 + (p2_tn - p2_t0) * time_passed;
      if (time_passed > 1.0) break;
    }

    int p1, p2, p3;
    do {int r1 = rng.next(N); int c1 = rng.next(N); p1 = P1D(r1, c1);} while (g[p1] == 0);
    
    int type = rng.next(20) == 0;
    if (C >= 3 && rng.next(5) == 0) type = 2;
    if (type == 0) {
      p2 = find_neighbor(p1);
    } else if (type == 1) {
      do {int r2 = rng.next(N); int c2 = rng.next(N); p2 = P1D(r2, c2);} while (g[p2] == 0 || g[p1] == g[p2]);
    } else if (type == 2) {
      p2 = find_neighbor(p1);
      p3 = find_neighbor(p1);
      if (g[p2] == g[p3]) continue;
    }

    if (type == 2) {
      dyn_count_components(comps, p1, p2, p3);
    } else {
      dyn_count_components(comps, p1, p2);
    }
 
    int new_components = comps.cnt[K];

    if (new_components >= components) {
      if (type == 2) {
        swap(g[p1], g[p2]);
        swap(g[p1], g[p3]);
      } else {
        swap(g[p1], g[p2]);
      }
      int turns = simulate_swaps(order);
      int score = new_components * P - turns;
      if (score >= bv || rng.next_double() < exp((score - bv) / p2_t)) {
        update_best(score, simulate_swaps);
        bv = score;
        components = new_components;
      } else {
        if (type == 2) {
          swap(g[p1], g[p3]);
          swap(g[p1], g[p2]);
        } else {
          swap(g[p1], g[p2]);
        }
        comps.restore();
      }
    } else {
      comps.restore();
    }
  }

  REP(r, N) REP(c, N) g[P1D(r,c)] = bg[P1D(r,c)];
  count_components(comps);
  components = comps.cnt[K];
  DB(p2_step, bscore, components);

  const double P3_STARTED = elapsed();

  int best_turns = 1<<24;

  int p3_step = 0;
  int p3_invalid = 0;
  while (true) {
    p3_step++;
    if ((p3_step & 1023) == 0) {
      double time_passed = (elapsed() - P3_STARTED) / (TIME_LIMIT_P3 - P3_STARTED);
      if (time_passed > 1.0) break;
    }

    int type = rng.next(5) == 0;
    int p1, p2;
    int nv, lv;
    p1 = rng.next(order.S);
    if (type <= 0) {
      do {p2 = rng.next(order.S);} while (abs(distall[order[p1]%MAXN2][last_middle] - distall[order[p2]%MAXN2][last_middle]) > (best_turns < (1<<20) ? N : 0));
      swap(order[p1], order[p2]);
    } else if (type == 1) {
      lv = order[p1];
      do {nv = order[p1] % MAXN2 + rng.next(4) * MAXN2;} while (lv == nv);
      order[p1] = nv;
    }

    int turns = simulate_swaps2(order);

    if (turns >= 1<<20) p3_invalid++;
    if (turns <= best_turns) {
      if (turns < best_turns) {
        best_turns = turns;
        int score = components * P - best_turns;
        update_best(score, simulate_swaps2);
      }
    } else {
      if (type == 0) {
        swap(order[p1], order[p2]);
      } else {
        order[p1] = lv;
      }
    }
  }

  DB(p1_step);
  DB(p2_step);
  DB(p3_step);
  DB(p3_invalid);

  DB(elapsed());

  count_components(comps);
  DB(VI(comps.cnt.begin(), comps.cnt.begin() + comps.mx + 1));
  cerr << "[DATA] p1_steps = " << p1_step << endl;
  cerr << "[DATA] p2_steps = " << p2_step << endl;
  cerr << "[DATA] p3_steps = " << p3_step << endl;
  cerr << "[DATA] comps = " << comps.cnt[K] << endl;
  cerr << "[DATA] turns = " << bturns << endl;
  cerr << "[DATA] max_comps = " << max_components << endl;

  if (timer1) DB(timer1);
  if (timer2) DB(timer2);
  if (timer3) DB(timer3);
  if (timer4) DB(timer4);

  cout << bsol.S << endl;
  for (auto &p : bsol) {
    PII p1 = P2D(p.X);
    PII p2 = P2D(p.Y);
    cout << p1.X << ' ' << p1.Y << ' ' << p2.X << ' ' << p2.Y << endl;
  }

	return 0;
}