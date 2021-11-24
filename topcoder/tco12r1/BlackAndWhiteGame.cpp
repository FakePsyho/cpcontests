#define SILENT
#define TIME_LIMIT 5.0
//#define SINGLE_PASS
#define TIMERS
//#define ERROR_CHECKING
#define QUICKFIX

#ifndef LOCAL
#undef TIME_LIMIT
#define TIME_LIMIT 9.8
#endif


#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <sstream>
#include <set>
#include <emmintrin.h>
#include <map>

using namespace std;

#define INLINE   __attribute__ ((always_inline))
#define NOINLINE __attribute__ ((noinline))

#define ALIGNED __attribute__ ((aligned(16)))

#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

#define FOR(i,a,b)  for(int i=(a);i<(b);++i)
#define REP(i,a)    FOR(i,0,a)
#define ZERO(m)     memset(m,0,sizeof(m))
#define ALL(x)      x.begin(),x.end()
#define PB          push_back
#define S           size()
#define LL          long long
#define ULL            unsigned long long
#define LD          long double
#define MP          make_pair
#define X           first
#define Y           second
#define VC          vector
#define PII         pair <int, int>
#define VI          VC < int >
#define VVI         VC < VI >
#define VD          VC < double >
#define VVD         VC < VD >
#define VS          VC < string >
#define DB(a)       cerr << #a << ": " << (a) << endl;

#define SSELOAD(a)     _mm_load_si128((__m128i*)&a)
#define SSESTORE(a, b) _mm_store_si128((__m128i*)&a, b)

void print(VI v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]\n";}
void print(VC < LL > v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]\n";}
void print(VD v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]\n";}
void print(VS v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]\n";}
template<class T> string i2s(T x) {ostringstream o; o << x; return o.str(); }
VS splt(string s, char c = ' ') {VS rv; int p = 0, np; while (np = s.find(c, p), np >= 0) {if (np != p) rv.PB(s.substr(p, np - p)); p = np + 1;} if (p < s.S) rv.PB(s.substr(p)); return rv;}

ULL getTicks() {
    ULL time;
    __asm__ volatile ("rdtsc" : "=A" (time));
    return time;
}

stringstream vMove, bMove;

double convertTicks(ULL time) {
#ifdef LOCAL
    return time / 1.66e9; 
#else
    return time / 3.6e9;
#endif
}

double getTime() {
    return convertTicks(getTicks());
}

struct Timer {
#ifdef TIMERS
    ULL ticks;
    const string name;
    bool running;
    int count;
    int count2;

    Timer(string _name) : name(_name) {
        ticks = 0;
        count = 0;
        count2 = 0;
        running = false;
    }
    
    void start() {
        if (running) return;
        ticks = getTicks() - ticks;
        running = true;
    }
    
    void stop(bool add = false) {
        if (!running) return;
        ticks = getTicks() - ticks;
        running = false;
        count++;
        if (add) count2++;
    }
    
    double elapsed() {
        return convertTicks(running ? getTicks() - ticks : ticks);
    }
    
    void show() {
        if (count == 0) return;
        cerr << name << ": " << elapsed() << " count: " << count;
        if (count2) cerr << " - " << count2;
        cerr << endl;
    }
#else
    Timer(string) { }
    void start() { }
    void stop(bool add = false) { }
    void show() { }
#endif
};

template <class T, class E, int SIZE> struct VCand {
    VC < pair < T, E > > cand;
    int curSize;
    T low;
    int lowPos;
    
    VCand() {
        cand = VC < pair < T, E > >(SIZE);
        curSize = 0;
        low = 1 << 30;
    }    
    
    void add(T v, E x) {
        if (curSize < SIZE) {
            if (v < low) {
                low = v;
                lowPos = curSize;
            }
            cand[curSize++] = MP(v, x);
        } else if (v < low) {
            cand[lowPos] = MP(v, x);
            low = 1 << 30;
            REP(i, SIZE) if (cand[i].X < low) {
                low = cand[i].X;
                lowPos = i;
            }
        }
    }
};

#define MAXN 128
#define NBITS 7

#ifdef SILENT
#define VERB(x) ;
#else
#define VERB(x) x;
#endif

#define BIG (1 << 28)

int TOTAL;
int N;
int tB[108][MAXN] ALIGNED;
short tD[108][MAXN] ALIGNED;
short tX[108][MAXN] ALIGNED;
int VNO = BIG;
int HAND = 0;

int curResult;
int curStep;

LL cnt0 = 0, cnt1 = 0, cnt2 = 0, cnt3 = 0, cnt4 = 0, cnt5 = 0;

typedef int (*N2TAB32)[MAXN];
typedef short (*N2TAB16)[MAXN];

int BONUS[MAXN][MAXN] ALIGNED;

N2TAB32 B = (N2TAB32)&tB[1][8];
N2TAB16 X = (N2TAB16)&tX[1][8];
N2TAB16 D = (N2TAB16)&tD[1][8];

int *B1 = (int*)B;
short *D1 = (short*)D;
short *X1 = (short*)X;

int YSUM[MAXN] ALIGNED;
int XSUM[MAXN] ALIGNED;

PII MOVES[MAXN * MAXN] ALIGNED;
int MOVESNO = 0;

#define MAX_SAVES 10
PII SMOVES[MAX_SAVES][MAXN * MAXN] ALIGNED;
int SMOVESNO[MAX_SAVES];
int SDATA[MAX_SAVES];
int STATESNO = 0;
//int BMOVESNO

#define MAXPMOVES 50
PII PMOVES[MAXN * MAXN] ALIGNED;
int PMOVESPOS[MAXPMOVES];
int PMOVESDATA0[MAXPMOVES];
int PMOVESDATA1[MAXPMOVES];
double PMOVESVAL[MAXPMOVES];
int PMOVESNO = 0;

const int dx[] = {-1,0,1,0};
const int dy[] = {0,-1,0,1};
const int d1[] = {-MAXN,-1,MAXN,1};

int STACK[MAXN * MAXN * 2] ALIGNED;
int STACKPOS = 0;
int STACKSTART = 0;

int TABU[MAXN][MAXN];
int TABUNO = 1000;

int MBIGCOMP = -1, XBIGCOMP = -1;
int MBIGLINE;
int MBIGX1, MBIGX2, MBIGY1, MBIGY2;

Timer timerConn("countMaxConn()");
Timer timerDist("calcDist()");
Timer timerMakeMove("makeMove()");
Timer timerHash("hashTable()");
Timer timerCheck("checkMove()");

Timer timerGreedy[6] = {Timer("greedyMove<1>()"), Timer("greedyMove<2>()"), Timer("greedyMove<3>()"), Timer("greedyMove<4>()"), Timer("greedyMove<5>()"), Timer("greedyMove<6>()")};
Timer timerGreedyX("GreedyXMove()");
Timer timerTeleport("TeleportMove()");
Timer timerSlow("SlowMove()");
Timer timerRandom("RandomMove()");
Timer timerCustom("Custom");
Timer timerGreedyFree("GreedyFreeMove()");
Timer timerPre("BegSim");
Timer timerCorrect("Correct");
Timer timerQuickFix("QuickFix");

INLINE bool inside(int x, int y) {
    return x >= 0 && x < N && y >= 0 && y < N;
}

INLINE bool inside(int p) {
    return p >= 0 && p < (N << NBITS) && (p & (MAXN - 1)) < N;
}

void show() {
    REP(i, N) {
        REP(j, N) cerr << (B[i][j] ? 'X' : ',');
        cerr << endl;
    }
}

void saveMoves(int p) {
    SMOVESNO[p] = MOVESNO;
    REP(i, MOVESNO) {
        SMOVES[p][i] = MOVES[i];
    }    
}

void loadMoves(int p) {
    MOVESNO = SMOVESNO[p];
    REP(i, MOVESNO) {
        MOVES[i] = SMOVES[p][i];
    }
}


int PARENT[4][MAXN * MAXN];
template <int t> int getParent(int p) {
    if (PARENT[t][p] < 0) 
        return p;
    else 
        return PARENT[t][p] = getParent<t>(PARENT[t][p]);
} 

void computeSets() {
    REP(i, N) { 
        REP(j, N) if (X[i][j]) {
            int p = (i << NBITS) + j;
            PARENT[0][p] = -1;
            if (X[i][j-1]) {
                int a = getParent<0>(PARENT[0][p-1]);
                PARENT[0][p] += PARENT[0][a];
                PARENT[0][a] = p;
            }
            if (X[i-1][j]) {
                int a = getParent<0>(PARENT[0][p-MAXN]);
                PARENT[0][p] += PARENT[0][a];
                PARENT[0][a] = p;
            }
        }
    }
}

int countMaxConn(bool updateGlobal) {
    timerConn.start();
    VNO--;
    
    int STARTVNO = VNO;
    
    int bv = 0;
    if (updateGlobal) {
        MBIGLINE = 0;
    }
    
    
    FOR(i, MBIGLINE, N) REP(j, N) if (B[i][j] > STARTVNO) {
        VNO--;
        
        int p = (i << NBITS) + j;
        STACK[STACKPOS++] = p;
        
        int av = 0;
        B1[p] = VNO;
        while (STACKPOS) {
            p = STACK[--STACKPOS];
            
            av++;
            
            int np;
            np = p + d1[0];
            if (B1[np] > VNO) {
                B1[np] = VNO;
                STACK[STACKPOS++] = np;
            }
            np = p + d1[1];
            if (B1[np] > VNO) {
                B1[np] = VNO;
                STACK[STACKPOS++] = np;
            }
            np = p + d1[2];
            if (B1[np] > VNO) {
                B1[np] = VNO;
                STACK[STACKPOS++] = np;
            }
            np = p + d1[3];
            if (B1[np] > VNO) {
                B1[np] = VNO;
                STACK[STACKPOS++] = np;
            }
        }
        
        if (av > bv) {
            bv = av;
            MBIGCOMP = VNO;
            if (updateGlobal) {
                MBIGLINE = av > TOTAL / 4 ? i : 0;
            }
            if (av > TOTAL / 2) goto out;
        }
    }
    out: 
    timerConn.stop();
    return bv;
}

int countMaxConn() {
    timerConn.start();
    VNO--;
    
    int STARTVNO = VNO;
    
    int bv = 0;
    
    FOR(i, MBIGX1, MBIGX2 + 1) FOR(j, MBIGY1, MBIGY2 + 1) if (X[i][j] == 2 && B[i][j] > STARTVNO) {
        VNO--;
        
        int p = (i << NBITS) + j;
        STACK[STACKPOS++] = p;
        
        int av = 0;
        B1[p] = VNO;
        while (STACKPOS) {
            p = STACK[--STACKPOS];
            
            av++;
            
            int np;
            np = p + d1[0];
            if (B1[np] > VNO) {
                B1[np] = VNO;
                STACK[STACKPOS++] = np;
            }
            np = p + d1[1];
            if (B1[np] > VNO) {
                B1[np] = VNO;
                STACK[STACKPOS++] = np;
            }
            np = p + d1[2];
            if (B1[np] > VNO) {
                B1[np] = VNO;
                STACK[STACKPOS++] = np;
            }
            np = p + d1[3];
            if (B1[np] > VNO) {
                B1[np] = VNO;
                STACK[STACKPOS++] = np;
            }
        }
        
        if (av > bv) {
            bv = av;
            if (av > curResult / 2) goto out;
        }
    }
    out: 
    timerConn.stop();
    return bv;
}

int countMaxConn(int sx, int sy) {
    timerConn.start();
    VNO--;
    
    int p = (sx << NBITS) + sy;
    STACK[STACKPOS++] = p;
    
    int av = 0;
    B1[p] = VNO;
    while (STACKPOS) {
        p = STACK[--STACKPOS];
        
        av++;
        
        int np;
        np = p + d1[0];
        if (B1[np] > VNO) {
            B1[np] = VNO;
            STACK[STACKPOS++] = np;
        }
        np = p + d1[1];
        if (B1[np] > VNO) {
            B1[np] = VNO;
            STACK[STACKPOS++] = np;
        }
        np = p + d1[2];
        if (B1[np] > VNO) {
            B1[np] = VNO;
            STACK[STACKPOS++] = np;
        }
        np = p + d1[3];
        if (B1[np] > VNO) {
            B1[np] = VNO;
            STACK[STACKPOS++] = np;
        }
    }
        
    timerConn.stop();
    
    return av;
}

int countMaxConn2() {
    while (STACKPOS) {
        int p = STACK[--STACKPOS];
        
    }    
}

void calcDist() {
    timerDist.start();
    
    MBIGX1 = 1 << 20, MBIGX2 = 0;
    MBIGY1 = 1 << 20, MBIGY2 = 0;
    
    VNO = 2 * BIG - VNO + 1;
    MBIGX1 = -1;
    
    __m128i mBIGCOMP = _mm_set1_epi32(MBIGCOMP);
    FOR(x, MBIGLINE, N) for (int y = 0; y < N; y += 4) {
        __m128i m = SSELOAD(B[x][y]);
        __m128i mcmp = _mm_cmpeq_epi32(m, mBIGCOMP);
        if (_mm_movemask_epi8(mcmp)) {
            MBIGX1 = MBIGX1 == -1 ? x : MBIGX1;
            MBIGX2 = x;
            int p = (x << NBITS) + y;
            if (B[x][y+0] == MBIGCOMP) {
                STACK[STACKPOS++] = p+0;
                B[x][y+0] = VNO + (1 << 30);
                MBIGY1 = min(MBIGY1, y+0);
                MBIGY2 = max(MBIGY2, y+0);
            }
            if (B[x][y+1] == MBIGCOMP) {
                STACK[STACKPOS++] = p+1;
                B[x][y+1] = VNO + (1 << 30);
                MBIGY1 = min(MBIGY1, y+1);
                MBIGY2 = max(MBIGY2, y+1);
            }
            if (B[x][y+2] == MBIGCOMP) {
                STACK[STACKPOS++] = p+2;
                B[x][y+2] = VNO + (1 << 30);
                MBIGY1 = min(MBIGY1, y+2);
                MBIGY2 = max(MBIGY2, y+2);
            }
            if (B[x][y+3] == MBIGCOMP) {
                STACK[STACKPOS++] = p+3;
                B[x][y+3] = VNO + (1 << 30);
                MBIGY1 = min(MBIGY1, y+3);
                MBIGY2 = max(MBIGY2, y+3);
            }
        }
    }
    
    STACK[STACKPOS++] = -1;
    
    int curDist = 0;
    STACKSTART = 0;
    
    FOR(i, MBIGY1 - 1, MBIGY2 + 2) B1[(MBIGX1 - 1) * MAXN + i] += 1 << 30;
    FOR(i, MBIGX1, MBIGX2 + 1) {
        B1[i * MAXN + MBIGY1 - 1] += 1 << 30;
        B1[i * MAXN + MBIGY2 + 1] += 1 << 30;
    }
    FOR(i, MBIGY1 - 1, MBIGY2 + 2) B1[(MBIGX2 + 1) * MAXN + i] += 1 << 30;
    
    while (true) {
        int p = STACK[STACKSTART++];
        
        if (p == -1) {
            if (STACKSTART == STACKPOS) break;
            curDist++;
            STACK[STACKPOS++] = -1;
            continue;
        }
        
        D1[p] = curDist;
        
        int np;
        np = p + d1[0];
        if (B1[np] < (1 << 30)) {
            B1[np] += 1 << 30;
            STACK[STACKPOS++] = np;
        }
        np = p + d1[1];
        if (B1[np] < (1 << 30)) {
            B1[np] += 1 << 30;
            STACK[STACKPOS++] = np;
        }
        np = p + d1[2];
        if (B1[np] < (1 << 30)) {
            B1[np] += 1 << 30;
            STACK[STACKPOS++] = np;
        }
        np = p + d1[3];
        if (B1[np] < (1 << 30)) {
            B1[np] += 1 << 30;
            STACK[STACKPOS++] = np;
        }
        
    }
    STACKPOS = 0;
    
    __m128i mVNO = _mm_set1_epi32(1 << 30);
    FOR(i, MBIGX1 - 1, MBIGX2 + 2) {
        int j = MBIGY1 - 1;
        while ((j & 3) && j < MBIGY2 + 2) B[i][j] -= 1 << 30, j++;
        while (j + 4 < MBIGY2 + 2) *(__m128i*)&B[i][j] = _mm_sub_epi32(*(__m128i*)&B[i][j], mVNO), j += 4;
        while (j < MBIGY2 + 2) B[i][j] -= 1 << 30, j++;
    }
    
    VNO = 2 * BIG - VNO;
    
    __m128i m8 = _mm_set1_epi16(4);
    FOR(x, MBIGX1, MBIGX2 + 1) {
        int p;
        p = x * MAXN + MBIGY1;
        while (p & 7) --p, D1[p] = D1[p + 1] + 1;
        if (p & (MAXN - 1)) {
            __m128i m = SSELOAD(D1[p]);
            while (p & (MAXN - 1)) {
                m = _mm_add_epi16(m, m8);
                p -= 8;
                SSESTORE(D1[p], m);
            }
        }
        
        p = x * MAXN + MBIGY2 + 1;
        while (p & 7) D1[p] = D1[p - 1] + 1, ++p;
        if ((p & (MAXN - 1)) < N) {
            __m128i m = SSELOAD(D1[p - 8]);
            while ((p & (MAXN - 1)) < N) {
                m = _mm_add_epi16(m, m8);
                SSESTORE(D1[p], m);
                p += 8;
            }
        }
    }
    
    __m128i m1 = _mm_set1_epi16(1);
    for (int x = MBIGX1 - 1; x >= 0; x--) {
        int p = x * MAXN;
        while ((p & (MAXN - 1)) < N) {
            __m128i m = SSELOAD(D1[p + MAXN]);
            m = _mm_add_epi16(m, m1);
            SSESTORE(D1[p], m);
            p += 8;
        }
    }
    
    FOR(x, MBIGX2 + 1, N) {
        int p = x * MAXN;
        while ((p & (MAXN - 1)) < N) {
            __m128i m = SSELOAD(D1[p - MAXN]);
            m = _mm_add_epi16(m, m1);
            SSESTORE(D1[p], m);
            p += 8;
        }
    }
        
    timerDist.stop();
}


template <int insert, int size, int dir> void makeMove(int x, int y) {
    timerMakeMove.start();

    if (dir == 0) {
        for (int i = N - 1; i >= size; i--) B[i][y] = B[i - size][y];
        REP(i, size) B[i][y] = insert << 28;
    } else if (dir == 1) {
        REP(i, N - size) B[i][y] = B[i+size][y];
        FOR(i, N - size, N) B[i][y] = insert << 28;
    } else if (dir == 2) {
        for (int i = N - 1; i >= size; i--) B[x][i] = B[x][i-size];
        REP(i, size) B[x][i] = insert << 28;
    } else {
        REP(i, N - size) B[x][i] = B[x][i+size];
        FOR(i, N - size, N) B[x][i] = insert << 28;
    }
    
    timerMakeMove.stop();
}

template <int insert, int size> void makeMove(int x, int y) {
    timerMakeMove.start();

    if (x == -1) {
        for (int i = N - 1; i >= size; i--) B[i][y] = B[i - size][y];
        REP(i, size) B[i][y] = insert << 28;
    } else if (x == N) {
        REP(i, N - size) B[i][y] = B[i+size][y];
        FOR(i, N - size, N) B[i][y] = insert << 28;
    } else if (y == -1) {
        for (int i = N - 1; i >= size; i--) B[x][i] = B[x][i-size];
        REP(i, size) B[x][i] = insert << 28;
    } else {
        REP(i, N - size) B[x][i] = B[x][i+size];
        FOR(i, N - size, N) B[x][i] = insert << 28;
    }
    
    timerMakeMove.stop();
}

template <int insert, int size> INLINE void makeMove(PII p) {
    makeMove<insert, size>(p.X, p.Y);
}

void moveToEdge(int &x, int &y, int d) {
    if (d == 0) {
        x = 0;
    } else if (d == 1) {
        y = 0;
    } else if (d == 2) {
        x = N - 1;
    } else {
        y = N - 1;
    }
}

int countCells() {
    int no = 0;
    REP(x, N) REP(y, N) no += B[x][y] > 0;
    return no;
}

template <int insert, int size> INLINE void makeRevMove(int x, int y) {
    makeMove<insert, size>(x == -1 ? N : x == N ? -1 : x, y == -1 ? N : y == N ? -1 : y);
}

template <int insert, int size> INLINE void makeRevMove(PII p) {
    makeRevMove<insert, size>(p.X, p.Y);
}

INLINE void makeMove(int x, int y) {
    int newHAND = (x == -1 ? B[N-1][y] : x == N ? B[0][y] : y == -1 ? B[x][N-1] : B[x][0]) > 0;
    if (!HAND) {
        makeMove<0, 1>(x, y);
    } else {
        makeMove<1, 1>(x, y);
    }
    HAND = newHAND;
}

INLINE void makeRevMove(int x, int y) {
    int newHAND = (x == N ? B[N-1][y] : x == -1 ? B[0][y] : y == N ? B[x][N-1] : B[x][0]) > 0;
    if (!HAND) {
        makeRevMove<0, 1>(x, y);
    } else {
        makeRevMove<1, 1>(x, y);
    }
    HAND = newHAND;
}

INLINE void makeMove(int i) {
    makeMove(MOVES[i].X, MOVES[i].Y);
}

INLINE void makeRevMove(int i) {
    makeRevMove(MOVES[i].X, MOVES[i].Y);
}

INLINE void makeMoves(int a, int b) {
    FOR(i, a, b) makeMove(i);
}

INLINE void makeRevMoves(int a, int b) {
    for (int i = b; i >= a; i--) makeRevMove(i);
}

void addMove(int x, int y, int d) {
    if (d == 0) {
        x = N;
    } else if (d == 1) {
        y = N;
    } else if (d == 2) {
        x = -1;
    } else {
        y = -1;
    }
    MOVES[MOVESNO++] = MP(x, y);
}

void addPMove(int x, int y, int d) {
    if (d == 0) {
        x = N;
    } else if (d == 1) {
        y = N;
    } else if (d == 2) {
        x = -1;
    } else {
        y = -1;
    }
    
    PMOVES[PMOVESPOS[PMOVESNO + 1]++] = MP(x, y);
}

void resetPMove() {
    PMOVESPOS[PMOVESNO + 1] = PMOVESPOS[PMOVESNO];
}

void nextPMove() {
    PMOVESNO++;
}

void clearPMoves() {
    PMOVESNO = 0;
}

int copyPMove(int p) {
    int rv = MOVESNO;
    FOR(i, PMOVESPOS[p], PMOVESPOS[p+1]) {
        MOVES[MOVESNO++] = PMOVES[i];
    }
    return rv; 
}


template < int size, int dir > void checkMove(int x, int y, int &bx, int &by, int &bv, int &no) {
    timerCheck.start();
    
    int sx = x, sy = y;
    int add = 0;
    if (dir == 0) {
        //FOR(x, max(0, MBIGX1 - size), min(N - size, MBIGX2 + 1)) {
        for (int x = min(N - size, MBIGX2); x >= max(0, MBIGX1 - size); x--) {
            if ((X[x][y] ^ X[x+size][y-1]) == 3 || (X[x][y] ^ X[x+size][y+1]) == 3) {
            sx = x+size; 
            goto ok;
          }
        }
    } else if (dir == 1) {
        //for (int x = min(N - 1, MBIGX2 + size); x >= max(size, MBIGX1); x--) {
        FOR(x, max(size, MBIGX1), min(N, MBIGX2 + 1 + size)) {
            if ((X[x][y] ^ X[x-size][y-1]) == 3 || (X[x][y] ^ X[x-size][y+1]) == 3) {
                sx = x-size; 
                goto ok;
            }
        }
    } else if (dir == 2) {
        //FOR(y, max(0, MBIGY1 - size), min(N - size, MBIGY2 + 1)) {
        for (int y = min(N - size, MBIGY2); y >= max(0, MBIGY1 - size); y--) {
            if ((X[x-1][y+size] ^ X[x][y]) == 3 || (X[x+1][y+size] ^ X[x][y]) == 3) {
                sy = y+size; 
                goto ok;
            }
        }
    } else {
        FOR(y, max(size, MBIGY1), min(N, MBIGY2 + 1 + size)) {
        //for (int y = min(N - 1, MBIGY2 + size); y >= max(size, MBIGY1); y--) {
            if ((X[x-1][y-size] ^ X[x][y]) == 3 || (X[x+1][y-size] ^ X[x][y]) == 3) {
                sy = y-size; 
                goto ok;
            }
        }
    }
    timerCheck.stop();
    return;
    
    ok:
    timerCheck.stop();
    makeMove<0, size, dir>(x, y);
    int av = countMaxConn(sx, sy);
    if (curStep > 100 && (curStep & 4) && av > curResult) av += rand() % 25 <= 0;
    if (av > bv || av == bv && rand() % no == 0) {
        no = av > bv ? 2 : no + 1;
        bv = av;
        bx = x;
        by = y;
    }
    makeMove<0, size, dir ^ 1>(x, y);
}

void checkMoveX(int x, int y, int &bx, int &by, int &bv, int &no, int handValue) {
    makeMove(x, y);
    int av = countMaxConn() + HAND * handValue;
    if (av > bv || av == bv && rand() % no == 0) {
        no = av > bv ? 2 : no + 1;
        bv = av;
        bx = x;
        by = y;
    }
    makeRevMove(x, y);
}

bool calcHVSumDone = false;
void calcHVSum() {
    if (calcHVSumDone) return;
    calcHVSumDone = true;

    memset(YSUM, 0, 4 * N);
    memset(XSUM, 0, 4 * N);
    REP(i, N) REP(j, N) {
        XSUM[i] += B[i][j] > 0;
        YSUM[j] += B[i][j] > 0;
    }
}

template <int size> bool greedyMove() {
    timerGreedy[size - 1].start();

    calcHVSum();
    
    int bx = 0, by = 0;
    int bv = 0;
    int no = 2;
    
    
    FOR(i, max(0, MBIGY1 - 1), min(N, MBIGY2 + 2)) if (YSUM[i]) {
        if (!B[N-1][i] && (size <= 1 || !B[N-2][i]) && (size <= 2 || !B[N-3][i]) && (size <= 3 || !B[N-4][i]) && (size <= 4 || !B[N-5][i]) && (size <= 5 || !B[N-6][i]))
            checkMove<size, 0>(-1, i, bx, by, bv, no);
        if (!B[0][i] && (size <= 1 || !B[1][i]) && (size <= 2 || !B[2][i]) && (size <= 3 || !B[3][i]) && (size <= 4 || !B[4][i]) && (size <= 5 || !B[5][i]))
            checkMove<size, 1>(N, i, bx, by, bv, no);
    }
    
    FOR(i, max(0, MBIGX1 - 1), min(N, MBIGX2 + 2)) if (XSUM[i]) {
        if (!B[i][N-1] && (size <= 1 || !B[i][N-2]) && (size <= 2 || !B[i][N-3]) && (size <= 3 || !B[i][N-4]) && (size <= 4 || !B[i][N-5]) && (size <= 5 || !B[i][N-6]))
            checkMove<size, 2>(i, -1, bx, by, bv, no);
        if (!B[i][0] && (size <= 1 || !B[i][1]) && (size <= 2 || !B[i][2]) && (size <= 3 || !B[i][3]) && (size <= 4 || !B[i][4]) && (size <= 5 || !B[i][5]))
            checkMove<size, 3>(i, N, bx, by, bv, no);
    }

    /*
    FOR(i, max(0, MBIGY1 - 1), min(N, MBIGY2 + 2)) if (YSUM[i]) {
        if (!B[N - size][i]) checkMove<size>(-1, i, bx, by, bv, no);
        if (!B[size - 1][i]) checkMove<size>(N, i, bx, by, bv, no);
    }
    
    FOR(i, max(0, MBIGX1 - 1), min(N, MBIGX2 + 2)) if (XSUM[i]) {
        if (!B[i][N - size]) checkMove<size>(i, -1, bx, by, bv, no);
        if (!B[i][size - 1]) checkMove<size>(i, N, bx, by, bv, no);
    }
    */
    
    if (bv <= curResult) {
        timerGreedy[size - 1].stop();
        return false;
    }
    
    //cerr << size << ' ' << bv - curResult << endl;
    
    VERB(vMove << "Greedy" << size << "Move: " << MOVESNO << ' ' << curResult << ' ' << bv << ' ' << bx << ' ' << by << endl);
    
    REP(i, size)
        MOVES[MOVESNO++] = MP(bx, by);
    //REP(i, size) makeMove(bx, by);
    makeMove<0, size>(bx, by);
    
    timerGreedy[size - 1].stop(true);
    return true;
}

int wallDist(int x, int y) {
    return min(min(x + 1, N - x), min(y + 1, N - y));
}


int ORIGINX1, ORIGINX2;
int ORIGINY1, ORIGINY2;
int originDist(int x, int y) {
    return (x < ORIGINX1 ? ORIGINX1 - x : x > ORIGINX2 ? x - ORIGINX2 : 0) + 
        (y < ORIGINY1 ? ORIGINY1 - y : y > ORIGINY2 ? x - ORIGINY2 : 0);
}

double BSUM[MAXN][MAXN];
void calculateOrigin() {
    ZERO(BSUM);
    REP(i, N + 1) REP(j, N + 1) if (B[i][j] > 0) {
        FOR(x, i - 3, i + 4) FOR(y, j - 3, j + 4) if (x >= 0 && x < N && y >= 0 && y <= N) {
            BSUM[x+1][y+1] += 1.0 / (abs(x - i) + abs(y - j) + 1);
        }
    }
    REP(i, N + 1) REP(j, N + 1) BSUM[i+1][j] += BSUM[i][j];
    REP(i, N + 1) REP(j, N + 1) BSUM[i][j+1] += BSUM[i][j];
    
    DB(BSUM[N][N]);
    
    double bv = 0;
    REP(sx, N + 1) REP(sy, N + 1) if (sx * sy == TOTAL) {
        REP(x, N - sx + 1) REP(y, N - sy + 1) {
            double av = BSUM[x+sx][y+sy] - BSUM[x][y+sy] - BSUM[x+sx][y] + BSUM[x][y];
            if (av > bv) {
                bv = av;
                ORIGINX1 = x;
                ORIGINX2 = x + sx - 1;
                ORIGINY1 = y;
                ORIGINY2 = y + sy - 1;
            }
        }
    }
    DB(ORIGINX1);
    DB(ORIGINX2);
    DB(ORIGINY1);
    DB(ORIGINY2);
}

int failedTeleports = 0;

bool originTeleportMove() {
    timerTeleport.start();
    
    int moveTry = 0;

    int md = 1 << 20;
    REP(i, N) REP(j, N) if (D[i][j] == 0) md = min(md, wallDist(i, j));
    md = min(md, ORIGINX1 + 5);
    md = min(md, N - ORIGINX2 + 4);
    md = min(md, ORIGINY1 + 5);
    md = min(md, N - ORIGINY2 + 4);

teleportMoveStart:
    int no = 2;
    int bv = 0;
    int bx = -1, by = -1;
    REP(i, N) REP(j, N) if (B[i][j] && D[i][j]) {
        int zv = min((int)D[i][j], originDist(i, j)) - wallDist(i, j) - md - 5 - rand() % 3;
        int av = 1000 - wallDist(i, j);
        if (zv >= 0 && (av > bv || av == bv && rand() % no == 0)) {
            no = av > bv ? 2 : no + 1;
            bx = i;
            by = j;
            bv = av;
        }
    }
    
    if (bv <= 0) {
        timerTeleport.stop();
        return false;
    }
    
    bv = 1 << 20;
    int bd = -1;
    REP(d, 4) {
        int nx = bx;
        int ny = by;
        moveToEdge(nx, ny, d);
        int av = abs(nx - bx) + abs(ny - by);
        if (av < bv) {
            bv = av;
            bd = d;
        }
    }
    
    int MOVESSNAP = MOVESNO;
    int cx = bx;
    int cy = by;
    REP(i, bv) {
        cx += dx[bd];
        cy += dy[bd];
        if (B[cx][cy]) {
            if (++moveTry <= 5) goto teleportMoveStart;
            timerTeleport.stop();
            return false;
        }
    }
    int movesNo = bv + 1;
    REP(i, movesNo) {
        addMove(bx, by, bd);
        makeMove(MOVESNO - 1);
    }
    
    md += 5;
    //int md2 = curStep % 4 == 0 || rand() % 8 == 0 ? 1 << 20 : md + 2;
    //int md2 = curStep % 2 == 0 ? 1 << 20 : md + 2;
    int md2 = md + 2;
    //int md2 = 1 << 20;
    
    VC < pair < int, pair < int, int > > > vp;
    
teleportMoveRepeat:
    vp.clear();
    vp.reserve(4 * N);
    
    bx = -2, by = -2;
    bv = 0;
    no = 2;
    
    
    REP(i, N) {
        if (D[0][i] < md && (!B[N-1][i] || D[N-1][i] > md2)) {
            int x = D[0][i];
            if (!x) while(B[-x][i]) x--;
            vp.PB(MP(x * MAXN + (rand() & (MAXN - 1)), MP(-1, i)));
        }
        if (D[N - 1][i] < md && (!B[0][i] || D[0][i] > md2)) {
            int x = D[N - 1][i];
            if (!x) while (B[N - 1 + x][i]) x--;
            vp.PB(MP(x * MAXN + (rand() & (MAXN - 1)), MP(N, i)));
        }
        if (D[i][0] < md && (!B[i][N-1] || D[i][N-1] > md2)) {
            int x = D[i][0];
            if (!x) while (B[i][-x]) x--;
            vp.PB(MP(x * MAXN + (rand() & (MAXN - 1)), MP(i, -1)));
        }
        if (D[i][N - 1] < md && (!B[i][0] || D[i][0] > md2)) {
            int x = D[i][N - 1];
            if (!x) while (B[i][N - 1 + x]) x--;
            vp.PB(MP(x * MAXN + (rand() & (MAXN - 1)), MP(i, N)));
        }
    }
    sort(ALL(vp));
    
    int testsNo = min((int)vp.S, 10);
    REP(i, testsNo) {
        checkMoveX(vp[i].Y.X, vp[i].Y.Y, bx, by, bv, no, 1 << 10);
    }
    
    VERB(vMove << "TeleportMove: " << MOVESNO << ' ' << curResult << ' ' << bv << ' ' << bx << ' ' << by << endl);
    
    if (bx == -2) {
        timerTeleport.stop();
        failedTeleports++;
        REP(i, movesNo) makeRevMove(--MOVESNO);
        return false;
    }
    
    MOVES[MOVESNO] = MP(bx, by);
    makeMove(MOVESNO++);

    if (HAND) {
        //countMaxConn(true);
        //calcDist();
        goto teleportMoveRepeat;
    }
    
    timerTeleport.stop(true);
    
    return true;
}

bool teleportMove() {
    timerTeleport.start();
    
    int moveTry = 0;

    int md = 1 << 20;
    REP(i, N) REP(j, N) if (D[i][j] == 0) md = min(md, wallDist(i, j));

teleportMoveStart:
    int no = 2;
    int bv = 0;
    int bx = -1, by = -1;
    REP(i, N) REP(j, N) if (B[i][j] && D[i][j]) {
        int zv = D[i][j] - wallDist(i, j) - md - 2 - rand() % 3;
        int av = zv;//1000 - wallDist(i, j) + rand() % 3;
        if (zv >= 0 && (av > bv || av == bv && rand() % no == 0)) {
            no = av > bv ? 2 : no + 1;
            bx = i;
            by = j;
            bv = av;
        }
    }
    
    if (bv <= 0) {
        timerTeleport.stop();
        return false;
    }
    
    bv = 1 << 20;
    int bd = -1;
    REP(d, 4) {
        int nx = bx;
        int ny = by;
        moveToEdge(nx, ny, d);
        int av = abs(nx - bx) + abs(ny - by);
        if (av < bv) {
            bv = av;
            bd = d;
        }
    }
    
    int MOVESSNAP = MOVESNO;
    int cx = bx;
    int cy = by;
    REP(i, bv) {
        cx += dx[bd];
        cy += dy[bd];
        if (B[cx][cy]) {
            if (++moveTry <= 5) goto teleportMoveStart;
            timerTeleport.stop();
            return false;
        }
    }
    int movesNo = bv + 1;
    REP(i, movesNo) {
        addMove(bx, by, bd);
        makeMove(MOVESNO - 1);
    }
    
    md += 5;
    //int md2 = curStep % 4 == 0 || rand() % 8 == 0 ? 1 << 20 : md + 2;
    int md2 = curStep % 2 == 0 ? 1 << 20 : md + 2;
    //int md2 = 1 << 20;
    
    VC < pair < int, pair < int, int > > > vp;
    
teleportMoveRepeat:
    vp.clear();
    vp.reserve(4 * N);
    
    bx = -2, by = -2;
    bv = 0;
    no = 2;
    
    
    REP(i, N) {
        if (D[0][i] < md && (!B[N-1][i] || D[N-1][i] > md2)) {
            int x = D[0][i];
            if (!x) while(B[-x][i]) x--;
            vp.PB(MP(x * MAXN + (rand() & (MAXN - 1)), MP(-1, i)));
        }
        if (D[N - 1][i] < md && (!B[0][i] || D[0][i] > md2)) {
            int x = D[N - 1][i];
            if (!x) while (B[N - 1 + x][i]) x--;
            vp.PB(MP(x * MAXN + (rand() & (MAXN - 1)), MP(N, i)));
        }
        if (D[i][0] < md && (!B[i][N-1] || D[i][N-1] > md2)) {
            int x = D[i][0];
            if (!x) while (B[i][-x]) x--;
            vp.PB(MP(x * MAXN + (rand() & (MAXN - 1)), MP(i, -1)));
        }
        if (D[i][N - 1] < md && (!B[i][0] || D[i][0] > md2)) {
            int x = D[i][N - 1];
            if (!x) while (B[i][N - 1 + x]) x--;
            vp.PB(MP(x * MAXN + (rand() & (MAXN - 1)), MP(i, N)));
        }
    }
    
    sort(ALL(vp));
    
    int testsNo = min((int)vp.S, 10);
    REP(i, testsNo) {
        checkMoveX(vp[i].Y.X, vp[i].Y.Y, bx, by, bv, no, 1);
    }
    
    VERB(vMove << "TeleportMove: " << MOVESNO << ' ' << curResult << ' ' << bv << ' ' << bx << ' ' << by << endl);
    
    if (bx == -2) {
        timerTeleport.stop();
        failedTeleports++;
        REP(i, movesNo) makeRevMove(--MOVESNO);
        return false;
    }
    
    MOVES[MOVESNO] = MP(bx, by);
    makeMove(MOVESNO++);

    if (HAND) {
        //countMaxConn(true);
        //calcDist();
        goto teleportMoveRepeat;
    }
    
    timerTeleport.stop(true);
    
    return true;
}


const int SLOWMOVE_RANGE = 12;

int BFSD[MAXN * MAXN];
int BFSQ[3][SLOWMOVE_RANGE * 16];
int BFSQNO[3];
void bfsMove(int x, int y) {
    calcHVSum();
    
    XSUM[x]--;
    YSUM[y]--;
    
    int len = min(N / 2, SLOWMOVE_RANGE * 2);
    
    REP(i, N) REP(j, N) BFSD[(i << NBITS) + j] = 1 << 30;
    
    BFSQNO[0] = BFSQNO[1] = BFSQNO[2];
    
    BFSQ[0][0] = (x << NBITS) + y;
    BFSQNO[0] = 1;
    
    REP(curDist, len) {
        int *Q0 = BFSQ[curDist % 3];
        int *Q1 = BFSQ[(curDist + 1) % 3];
        int *Q2 = BFSQ[(curDist + 2) % 3];
        
        int &QNO0 = BFSQNO[curDist % 3];
        int &QNO1 = BFSQNO[(curDist + 1) % 3];
        int &QNO2 = BFSQNO[(curDist + 2) % 3];
        
        QNO2 = 0;
        
        REP(i, QNO0) {
            int p = Q0[i];
            if (BFSD[p] < curDist) continue;
            
            int cost;
            
            cost = XSUM[p >> NBITS] ? 2 : 1;
            if (BFSD[p - MAXN] > curDist + cost) {
                BFSD[p - MAXN] = curDist + cost;
                if (cost == 2) {
                    Q2[QNO2++] = p - MAXN;
                } else {
                    Q1[QNO1++] = p - MAXN;
                }
            }
            if (BFSD[p + MAXN] > curDist + cost) {
                BFSD[p + MAXN] = curDist + cost;
                if (cost == 2) {
                    Q2[QNO2++] = p + MAXN;
                } else {
                    Q1[QNO1++] = p + MAXN;
                }
            }
            
            cost = YSUM[p & (MAXN - 1)] ? 2 : 1;
            if (BFSD[p - 1] > curDist + cost) {
                BFSD[p - 1] = curDist + cost;
                if (cost == 2) {
                    Q2[QNO2++] = p - 1;
                } else {
                    Q1[QNO1++] = p - 1;
                }
            }
            if (BFSD[p + 1] > curDist + cost) {
                BFSD[p + 1] = curDist + cost;
                if (cost == 2) {
                    Q2[QNO2++] = p + 1;
                } else {
                    Q1[QNO1++] = p + 1;
                }
            }
        }
    }
    
    XSUM[x]++;
    YSUM[y]++;
}

bool createBFSPath(int x, int y, int px, int py) {
    
    XSUM[x]--;
    YSUM[y]--;
    
    int curTry = 0;
    
    const int MAXFRAG = 20;
    int pathlen[MAXFRAG];
    int pathdir[MAXFRAG];
    int pathcost[MAXFRAG];
    
retry:
    resetPMove();
    int p = (px << NBITS) + py;
    int dir = -1;
    int fragno = -1;
    
    bool empty = true;
    while (BFSD[p]) {
        int cost;
        int posDir = -1;
        
        cost = XSUM[p >> NBITS] ? 2 : 1;
        if (BFSD[p - MAXN] + cost == BFSD[p])
            posDir |= 1;
        if (BFSD[p + MAXN] + cost == BFSD[p])
            posDir |= 4;
        
        cost = YSUM[p & (MAXN - 1)] ? 2 : 1;
        if (BFSD[p - 1] + cost == BFSD[p])
            posDir |= 2;
        if (BFSD[p + 1] + cost == BFSD[p])
            posDir |= 8;
            
           if (posDir == 0) {
               cerr << "Error in createBFSPath()" << endl;
               return false;
           }
           
           int ndir;
           if (posDir & (1 << dir)) 
               ndir = dir;
           else {
               do {
                   ndir = rand() & 3;
               } while (!(posDir & (1 << ndir)));
           }
           
           if (dir != ndir) {
               fragno++;
               pathdir[fragno] = ndir;
               pathlen[fragno] = 0;
               pathcost[fragno] = BFSD[p] - BFSD[p + d1[ndir]];
           }
           pathlen[fragno]++;
           
           p += d1[ndir];
           dir = ndir;
       
    }
    fragno++;
    if (fragno == 1) return false;
    
    REP(i, fragno - 1) if (pathcost[i] == 2) 
    
    XSUM[x]++;
    YSUM[y]++;
}

int FXSUM[MAXN];
int FYSUM[MAXN];


int randomdir;
bool slowMove(int x, int y) {
	const int SLOWMOVE_CHECK = 6;

    calcHVSum();
    
    double bv = 0.0;//0;
    int bm = 0;
    int bx = -1, by = -1;
    int no = 2;
    
    int len = min(N / 3, SLOWMOVE_RANGE);
    //int len = min((int)D[x][y], min(N / 3, SLOWMOVE_RANGE));
    
    STACKPOS = 0;
    
	int bcost = D[x][y] * 4 + 3;
    REP(ttd, 4) {
        int d = (ttd + randomdir) & 3;
        int nx = x;
        int ny = y;
        
        int mx = x;
        int my = y;
        
        moveToEdge(mx, my, d);
        mx += dx[d];
        my += dy[d];

        int m = 0;
        
		bool fastUsed = false;
			
        while (m < len) {
            m++;
            nx += dx[d];
            ny += dy[d];
            if (!inside(nx, ny)) break;
			
            mx -= dx[d];
            my -= dy[d];
            if (B[mx][my] || B[nx][ny]) continue;
            
            int tt = dx[d] ? YSUM[y] - FYSUM[y] : XSUM[x] - FXSUM[x];
			tt++;
            
            int cost = (tt == 1 ? m : m * 2 + 2);
			if (cost >= bcost) break;
            //double av = (D[nx][ny] == 1) ? 1.0 + 1.0 / cost : (double)(D[x][y] - D[nx][ny]) / cost;
            double av = (double)(D[x][y] - D[nx][ny]) / ((BONUS[nx][ny] && D[x][y] > 2) ? cost - 2 : cost) + (D[nx][ny] == 1);
            //double av = (double)(D[x][y] - D[nx][ny]) / cost + (D[nx][ny] == 1);
            //if (B[nx][ny] == 0 && (av > bv + 1e-6 || fabs(av - bv) < 1e-6 && rand() % no == 0)) {
            if (B[nx][ny] == 0 && av > bv) {
                REP(dd, 2) {
                    int nd = dd == 0 ? ((d + 3) & 3) : ((d + 1) & 3);
                    int vx = nx + dx[nd];
                    int vy = ny + dy[nd];
                    
                    if (!inside(nx - dx[nd], ny - dy[nd]) || B[vx][vy]) continue;
                    int mvx = vx, mvy = vy;
                    moveToEdge(mvx, mvy, nd ^ 2);
                    if (B[mvx][mvy]) continue;
                    
					if (D[nx][ny] == 1) bcost = cost;
                    bm = m;
                    bv = av;
                    resetPMove();
                    REP(i, m) addPMove(x, y, d);
                    if (tt > 1) {
                        addPMove(nx, ny, nd ^ 2);
                        REP(i, m) addPMove(x, y, d ^ 2);
                        addPMove(nx, ny, nd);
                    }
                    bx = nx;
                    by = ny;
                    PMOVESDATA1[PMOVESNO] = cost;
                }
            }
            
            REP(dd, 2) {
                int nd = dd == 0 ? ((d + 3) & 3) : ((d + 1) & 3);
                int m2 = 0;
                int vx = nx;
                int vy = ny;
                
                int mx2 = nx;
                int my2 = ny;
                moveToEdge(mx2, my2, nd ^ 2);
                mx2 -= dx[nd];
                my2 -= dy[nd];
                
                int tt2 = dx[nd] ? YSUM[ny] - FYSUM[ny] : XSUM[nx] - FXSUM[nx];
                				
                while (m2 + m < len) {
                    m2++;
                    vx += dx[nd];
                    vy += dy[nd];
                    if (!inside(vx, vy)) break;
                    
					//if (m2 + m > SLOWMOVE_CHECK && D[vx][vy] > D[x][y]) break;
					
                    mx2 += dx[nd];
                    my2 += dy[nd];
                    if (B[mx2][my2]) continue;
                    
                    int cost = m * (tt > 1 ? 2 : 1) + m2 * (tt2 > 0 ? 2 : 1);
					if (cost >= bcost) break;
                    //double av = (D[vx][vy] == 1) ? 1.0 + 1.0 / cost : (double)(D[x][y] - D[vx][vy]) / cost;
                    double av = (double)(D[x][y] - D[vx][vy]) / ((BONUS[vx][vy] && D[x][y] > 2) ? cost - 2 : cost) + (D[vx][vy] == 1);
                    //double av = (double)(D[x][y] - D[vx][vy]) / cost + (D[vx][vy] == 1);
                    //if (B[nx][ny] == 0 && (av > bv + 1e-6 || fabs(av - bv) < 1e-6 && rand() % no == 0)) {
                    if (B[vx][vy] == 0 && av > bv) {
						if (D[vx][vy] == 1) bcost = cost;
                        //bm = m;
                        bv = av;
                        resetPMove();
                        if (tt2 > 0) 
                          REP(i, m2) addPMove(nx, ny, nd ^ 2);
                        REP(i, m) addPMove(nx, ny, d);
                        REP(i, m2) addPMove(nx, ny, nd);
                        if (tt > 1) 
                          REP(i, m) addPMove(nx, ny, d ^ 2);
                        bx = vx;
                        by = vy;
                        PMOVESDATA1[PMOVESNO] = cost;
						break;
                    }
					
					
					if (tt2 || fastUsed) continue;
					
					
					int zx = vx;
					int zy = vy;
					
					int mx3 = vx;
					int my3 = vy;
					moveToEdge(mx3, my3, d ^ 2);
					
					int tt3 = dx[d] ? YSUM[vy] - FYSUM[vy] : XSUM[vx] - FXSUM[vx];
					int m3 = 0;
					while (m3 + m2 + m < len) {
						m3++;
						zx += dx[d];
						zy += dy[d];
						
						if (!inside(zx, zy)) break;
						
						//if (m3 + m2 + m > SLOWMOVE_CHECK && D[zx][zy] > D[x][y]) break;
						
						mx3 += dx[d];
						my3 += dy[d];
						if (B[mx3][my3]) continue;
						
						int cost = m * (tt > 1 ? 2 : 1) + m2 + m3 * (tt3 > 0 ? 2 : 1);
						if (cost >= bcost) break;
						double av = (double)(D[x][y] - D[zx][vy]) / ((BONUS[zx][zy] && D[x][y] > 2) ? cost - 2 : cost) + (D[zx][zy] == 1);
						if (B[zx][zy] == 0 && av > bv) {
							bv = av;
							if (D[zx][zy] == 1) bcost = cost;
							resetPMove();
							if (tt3 > 0)
								REP(i, m3) addPMove(vx, vy, d ^ 2);
							REP(i, m) addPMove(x, y, d);
							REP(i, m2) addPMove(nx, ny, nd);
							REP(i, m3) addPMove(vx, vy, d);
							if (tt > 1)
								REP(i, m) addPMove(x, y, d ^ 2);
							PMOVESDATA1[PMOVESNO] = cost;
							break;
						}
					}
                }
            }
			
			fastUsed |= (dx[d] ? XSUM[nx] - FXSUM[nx] : YSUM[ny] - FYSUM[ny]) == 0;
        }
        
    }
    
    if (bx != -1) {
        PMOVESDATA0[PMOVESNO] = D[bx][by] == 1;
        
        //VERB(vMove << "SlowMove: (" << x << ", " << y << ") -> (" << bx << ", " << by << ")" << endl);
        PMOVESVAL[PMOVESNO] = bv;
        return true;
    }
    
    VERB(vMove << "SlowMove Failed" << endl);
    
    return false;
}

bool randomMove(int x, int y, int len) {
    timerRandom.start();

    int a = rand() % 4;
    REP(d, 4) {
        int ad = (a + d) & 3;
        int nx = x, ny = y;
        moveToEdge(nx, ny, ad);
        
        int moved = 0;
        REP(i, len) {
            if (B[nx][ny]) break;
            addMove(x, y, ad);
            makeMove(MOVESNO - 1);
            moved++;
        }
        if (moved) {
            VERB(vMove << "RandomMove: " << x << ' ' << y << ' ' << moved << endl);
            timerRandom.stop(true);
            return true;
        }
    }
    
    VERB(vMove << "RandomMove Failed" << endl);
    
    timerRandom.stop();
    return false;
}

bool greedyFreeMove() {
	timerGreedyFree.start();
	int bx = -2, by = -2;
	int bv = 2;
	int no = 2;
	REP(x, N) if (FXSUM[x] > 1 && !B[x][0] && !B[x][N-1]) {
		int ds = 0;
		int d0 = 0;
		int d1 = 0;
		REP(y, N) if (B[x][y] && D[x][y]) {
			ds += D[x][y];
			d0 += D[x][y-1];
			d1 += D[x][y+1];
		}
		d0 = ds - d0;
		d1 = ds - d1;
		if (d0 > bv || d0 == bv && rand() % no == 0) {
			cnt4++;
			makeMove<0, 1>(x, N);
			if (countMaxConn() == curResult) {
				no = d0 > bv ? 2 : no + 1;
				bx = x;
				by = N;
				bv = d0;
			}
			makeRevMove<0, 1>(x, N);
		}
		if (d1 > bv || d1 == bv && rand() % no == 0) {
			cnt4++;
			makeMove<0, 1>(x, -1);
			if (countMaxConn() == curResult) {
				no = d1 > bv ? 2 : no + 1;
				bx = x;
				by = -1;
				bv = d1;
			}
			makeRevMove<0, 1>(x, -1);
		}
	}
	
	REP(y, N) if (FYSUM[y] > 1 && !B[0][y] && !B[N-1][y]) {
		int ds = 0;
		int d0 = 0;
		int d1 = 0;
		REP(x, N) if (B[x][y] && D[x][y]) {
			ds += D[x][y];
			d0 += D[x-1][y];
			d1 += D[x+1][y];
		}
		d0 = ds - d0;
		d1 = ds - d1;
		if (d0 > bv || d0 == bv && rand() % no == 0) {
			cnt4++;
			makeMove<0, 1>(N, y);
			if (countMaxConn() == curResult) {
				no = d0 > bv ? 2 : no + 1;
				bx = N;
				by = y;
				bv = d0;
			}
			makeRevMove<0, 1>(N, y);
		}
		if (d1 > bv || d1 == bv && rand() % no == 0) {
			cnt4++;
			makeMove<0, 1>(-1, y);
			if (countMaxConn() == curResult) {
				no = d1 > bv ? 2 : no + 1;
				bx = -1;
				by = y;
				bv = d1;
			}
			makeRevMove<0, 1>(-1, y);
		}
	}
	
	if (bx == -2) {
		timerGreedyFree.stop();
		return false;
	}
	
    MOVES[MOVESNO++] = MP(bx, by);
    makeMove<0, 1>(bx, by);
	timerGreedyFree.stop(true);
	return true;
}

int MOVESNOHIST[1 << 16];
class BlackAndWhiteGame {public: VS makeConnected(VS board) {
    double startTime = getTime();

    srand(2);
    Timer gTimer("gTimer");
    gTimer.start();
    
    REP(i, MAX_SAVES) SMOVESNO[i] = 1 << 20;
    
    N = board.S;
    REP(i, N) REP(j, N) B[i][j] = (board[i][j] == 'X') << 28;
    HAND = 0;
    
    TOTAL = countCells();
    calculateOrigin();
    
    cerr << "N: " << N << " TOTAL: " << TOTAL << " FILL: " << (100.0 * TOTAL / N / N) << "%" << endl;
    int BMOVESNO = N * N / 2;
    
    curStep = 0;
    int failures = 0;
    int totalFailedTeleports = 0;
    
    while (true) {
        curStep++;
        
        MOVESNO = 0;
        
        REP(i, N) REP(j, N) B[i][j] = (board[i][j] == 'X') << 28;
    
        bool done = false;
        
        bool slowMovePhase = false;
        //bool optimizePhase = STATESNO && rand() % 2;
        
        double timePassed = (getTime() - startTime) / TIME_LIMIT;
        
        bool optimizePhase = STATESNO && timePassed > 0.5;
        
        int randomlen = 1;
        
        VERB(vMove.str(""));
        
        timerPre.start();
        int curState = STATESNO < 2 ? STATESNO : rand() % STATESNO;
        if (optimizePhase) {
            curState = rand() % STATESNO; 
            loadMoves(curState);
            int movesNo = (int)(MOVESNO * pow(timePassed, 2.0) * 0.9);
            MOVESNO = 0;
            while (MOVESNO < movesNo || HAND)
                makeMove(MOVESNO++);
        }
        timerPre.stop();
        int lastCorrect = MOVESNO - 1;
        
        bool originTeleportMovePhase = curStep & 2;
        failedTeleports = 0;
        randomdir = rand() & 3;
        
        int mstep = 0;
        while (true) {
            calcHVSumDone = false;
            clearPMoves();
            
            if (getTime() - startTime > TIME_LIMIT) {
                cerr << "Out of time!" << endl;
                break;
            }
        
            timerCorrect.start();
            while (lastCorrect >= 0 && lastCorrect + 1 < MOVESNO && (MOVES[lastCorrect + 1].X * MOVES[lastCorrect].X == -N && MOVES[lastCorrect + 1].Y == MOVES[lastCorrect].Y || MOVES[lastCorrect + 1].Y * MOVES[lastCorrect].Y == -N && MOVES[lastCorrect + 1].X == MOVES[lastCorrect].X)) {
                FOR(i, lastCorrect + 2, MOVESNO) {
                    MOVES[i - 2] = MOVES[i];
                }
                MOVESNO -= 2;
                lastCorrect--;
            }
            lastCorrect = MOVESNO - 1;
            timerCorrect.stop();
            
            if (mstep > 100 && MOVESNOHIST[mstep - 100] >= MOVESNO) {
                failures++;
                break;
            }
            MOVESNOHIST[mstep++] = MOVESNO;
    
            curResult = countMaxConn(true);
            
            XBIGCOMP = MBIGCOMP;
            REP(i, N) REP(j, N) X[i][j] = (B[i][j] == MBIGCOMP) ? 2 : (B[i][j] > 0);
            if (curResult == TOTAL) {
                done = true;
                break;
            }

            if (curState == STATESNO && MOVESNO >= N * N / 2 - 1) break;
            if (curState != STATESNO && MOVESNO >= SMOVESNO[curState] - 1) break;
            
#ifdef ERROR_CHECKING            
            {
                int no = 0;
                REP(i, N) REP(j, N) no += B[i][j] > 0;
                if (no != TOTAL) {
                    cerr << "Error: Wrong number of white cells" << endl;
                    break;
                }
            }
#endif

            calcDist();
            
            int md = 1 << 20;
            REP(i, N) REP(j, N) if (B[i][j] && D[i][j]) md = min(md, (int)D[i][j]);
            
            VERB(vMove << MOVESNO << ' ' << countMaxConn() << ' ' << md << ' ' << (getTime() - startTime) << endl);
            
            if (failedTeleports < 5 && originTeleportMovePhase && originTeleportMove()) continue;
            originTeleportMovePhase = false;
            if (failedTeleports < 5 && curResult > TOTAL * 3 / 4 && teleportMove()) continue;
            if (md <= 2 && greedyMove<1>()) continue;
            //if (failedTeleports < 5 && curResult <= TOTAL * 4 / 5 && curResult > TOTAL * 3 / 4 && teleportMove()) continue;
            if (md <= 3 && greedyMove<2>()) continue;
            if (failedTeleports < 5 && curResult <= TOTAL * 3 / 4 && teleportMove()) continue;
            if (md <= 4 && greedyMove<3>()) continue;
            if (md <= 5 && greedyMove<4>()) continue;
            if (md <= 6 && greedyMove<5>()) continue;
            //if (md <= 7 && greedyMove<6>()) continue;
            
            slowMovePhase = true;
            
            timerSlow.start();
            bool slowMoveDone = false;
            int smx = -1, smy = -1;
            
            REP(i, N) REP(j, N) BONUS[i][j] = 0;
            REP(i, N) REP(j, N) if (B[i][j] && D[i][j] == 2) REP(d, 4) if (D[i+dx[d]][j+dy[d]] == 1) BONUS[i+dx[d]][j+dy[d]] = 1;
            VCand <int, int, 15> slowCand;
            int totalFree = 0;
			REP(i, N) FXSUM[i] = 0;
			REP(i, N) FYSUM[i] = 0;
            REP(i, N) REP(j, N) if (B[i][j] && D[i][j]) {
                totalFree++;
				FXSUM[i]++;
				FYSUM[j]++;
                int av = (D[i][j] << 0) + (rand() & 1023);
                slowCand.add(av, (i << NBITS) + j);
            }
			
			//if (greedyFreeMove()) continue;

            
            REP(i, min(totalFree / 6 + 3, slowCand.curSize)) {
                int x = slowCand.cand[i].Y >> NBITS;
                int y = slowCand.cand[i].Y & (MAXN - 1);
                if (slowMove(x, y)) {
                    slowMoveDone = true;
                    nextPMove();
                }
            }
            
            timerSlow.stop(slowMoveDone);
            
            int bestMove = 0;
            if (slowMoveDone) {
                int no = 2;
                FOR(i, 1, PMOVESNO) if (PMOVESVAL[i] > PMOVESVAL[bestMove] || fabs(PMOVESVAL[i] - PMOVESVAL[bestMove]) < 1e-6 && rand() % no == 0) {
                //FOR(i, 1, PMOVESNO) if (PMOVESVAL[i] > PMOVESVAL[bestMove]) {
                    no = PMOVESVAL[i] > PMOVESVAL[bestMove] ? 2 : no + 1;
                    bestMove = i;
                }
            }
            
            if (slowMoveDone) {
                timerCustom.start();
                
                int p = copyPMove(bestMove);
                if (PMOVESDATA0[bestMove]) {
                    makeMoves(p, MOVESNO);
                    int exp = countMaxConn();
                    
                    for (int i = MOVESNO - 1; i >= p; i--) {
                        makeRevMove(i);
                        int size = countMaxConn();
                        if (!HAND && size >= exp) {
                                MOVESNO = i;
                                exp = max(size, exp);
                        }
                    }
                    
                    repeat:
                    for (int i = MOVESNO - 1; i >= p; i--) {
                        if (i == p || MOVES[i] != MOVES[i-1]) {
                            FOR(j, p, MOVESNO) if (j != i) makeMove(j);
                            int size = countMaxConn();
                            if (!HAND && size >= exp) {
                                FOR(j, i, MOVESNO - 1) MOVES[j] = MOVES[j+1];
                                MOVESNO--;
                                makeRevMoves(p, MOVESNO - 1);
                                exp = max(size, exp);
                                goto repeat;
                            }
                            makeRevMoves(i + 1, MOVESNO - 1);
                            makeRevMoves(p, i - 1);
                        }
                    }
                }
                makeMoves(p, MOVESNO);
                timerCustom.stop();
                continue;
            }
            
            
            randomlen++;
            if (smx >= 0 && rand() % 2 ? randomMove(smx, smy, rand() % (randomlen / 5 + 1) + 1) : randomMove(rand() % N, rand() % N, rand() % (randomlen / 5 + 1) + 1)) continue;
        }
        
        totalFailedTeleports += failedTeleports;
        
        if (done && MOVESNO < SMOVESNO[curState]/* || curStep == 1*/) {
#ifdef QUICKFIX
			timerQuickFix.start();
			REP(i, N) REP(j, N) B[i][j] = (board[i][j] == 'X') << 28;
			makeMoves(0, MOVESNO);
			curResult = countMaxConn(true);
			timerQuickFix.stop();
			bool correct = (curResult == TOTAL);
#else		
			bool correct = true;
#endif
			if (correct) {
				if (curState == STATESNO) STATESNO++;
				saveMoves(curState);
				if (MOVESNO < BMOVESNO)
					cerr << "New: " << curStep << ' ' << curState << ' ' << MOVESNO << ' ' << (1.0 * (N * N - MOVESNO) / (N * N)) << ' ' << (getTime() - startTime) << endl;
	BMOVESNO = min(MOVESNO, BMOVESNO);
				VERB(bMove.str(vMove.str()));
			}
		}
        
        if (getTime() - startTime > TIME_LIMIT) break;
        
#ifdef SINGLE_PASS
        if (done && MOVESNO < N * N * 0.4) break;
#endif
    }
    
    DB(curStep);
    
    VS rv;
    int bestState = 0;
    REP(i, STATESNO) if (SMOVESNO[i] < SMOVESNO[bestState]) bestState = i;
    REP(i, SMOVESNO[bestState]) {
        string s = i2s(SMOVES[bestState][i].X) + " " + i2s(SMOVES[bestState][i].Y);
        rv.PB(s);
    }
    
    print(VI(SMOVESNO, SMOVESNO + STATESNO));
    
    gTimer.show();
    timerConn.show();
    timerDist.show();
    timerHash.show();
    timerMakeMove.show();
    timerCheck.show();
    timerGreedy[0].show();
    timerGreedy[1].show();
    timerGreedy[2].show();
    timerGreedy[3].show();
    timerGreedy[4].show();
    timerGreedy[5].show();
    //timerGreedyX.show();
    timerTeleport.show();
    timerSlow.show();
    //timerRandom.show();
    timerCustom.show();
	timerGreedyFree.show();
    //timerPre.show();
    //timerCorrect.show();
	timerQuickFix.show();
    
    if (cnt0) DB(cnt0);
    if (cnt1) DB(cnt1);
    if (cnt2) DB(cnt2);
    if (cnt3) DB(cnt3);
    if (cnt4) DB(cnt4);
    if (cnt5) DB(cnt5);
    
    if (totalFailedTeleports) DB(totalFailedTeleports);
    if (failures) DB(failures);
    
    VERB(cerr << bMove.str());
    
    return rv;
}};



#ifdef
int main() {
	int n;
	cin >> n;
	VS board(n);
	REP(i, n) cin >> board[i];
	
	BlackAndWhiteGame algo;
	VS rv = algo.makeConnected(board);
	
	REP(i, 100000000) n += n; //wait
	n++;
	
	if (n) cout << rv.S << endl;
	REP(i, rv.S) cout << rv[i] << endl;
	
	fflush(stdout);
}
#endif
