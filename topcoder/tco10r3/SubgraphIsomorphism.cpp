#define TOTAL_TIME (10.0 + bonusTime)
#define MAX_ESTIMATE_TIME 1.0

//#define USE_FIXEDORDER
#define MINACCRATIO 2.0
#define PRINT_LEVCOUNTER 0
#define DUMP_GRAPHS 1
#define DEBUG_INFO 1
#define ANALYZE_TEST 0
#define ORDEREST_NO 5

#define ENTRYTEST_TIME 0.01
#define FULLTEST_TIME (EP[N] / N < 15 ? 0.0025 : 0.005)
#define MAX_COMP_LEVEL 31

int EODLev[] = {50, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1};
int EODMaxLev = sizeof (EODLev) / sizeof (int);

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <map>
#include <set>
#include <list>
#include <deque>
#include <queue>
#include <sys/time.h>

using namespace std;

#define FOR(i,a,b)  for(int i=(a), i##0 = (b); i<i##0; ++i)
#define REP(i,a)    FOR(i,0,a)
#define FORE(i,a)   for(__typeof(a.begin()) i=a.begin();i!=a.end();i++)
#define ZERO(m)     memset(m,0,sizeof(m))
#define ALL(x)      x.begin(),x.end()
#define PB          push_back
#define SZ          size()
#define UINT        unsigned int
#define LL          long long
#define LD          long double
#define MP          make_pair
#define PII         pair < int, int >
#define PL          pair < PII, PII >
#define X           first
#define Y           second
#define VC          vector
#define VI          VC<int>
#define VD            VC<double>
#define VVI         VC< VI >
#define VS          VC<string>
#define VPII        VC< PII >
#define DB(a)       cerr << #a << ": " << a << endl;

void print(VI v) {cerr << "[";if (v.SZ) cerr << v[0];FOR(i, 1, v.SZ) cerr << ", " << v[i];cerr << "]\n";}
void print(VD v) {cerr << "[";if (v.SZ) cerr << v[0];FOR(i, 1, v.SZ) cerr << ", " << v[i];cerr << "]\n";}
void print(VS v) {cerr << "[";if (v.SZ) cerr << v[0];FOR(i, 1, v.SZ) cerr << ", " << v[i];cerr << "]\n";}
void print(VVI v) {cerr << "[ ---";if (v.SZ) cerr << " ", print(v[0]);FOR(i, 1, v.SZ) cerr << " ", print(v[i]);    cerr << "--- ]\n";}
void print(PII p) {cerr << "{" << p.X << ", " << p.Y << "}";}
void print(VPII v) {cerr << "[";if (v.SZ) print(v[0]);FOR(i, 1, v.SZ)  cerr << ", ", print(v[i]);cerr << "]\n";}

template<class T> string i2s(T x) {ostringstream o; o << x; return o.str(); }
VS splt(string s, char c = ' ') {VS rv; int p = 0, np; while (np = s.find(c, p), np >= 0) {if (np != p) rv.PB(s.substr(p, np - p)); p = np + 1;} if (p < s.SZ) rv.PB(s.substr(p)); return rv;}

#define RNG_MUL 16807
#define RNG_ADD 0

UINT RNGSeed = 0;
void SetSeed(UINT seed) {
    RNGSeed = seed;
}

UINT Rand() {
    return RNGSeed = RNGSeed * RNG_MUL + RNG_ADD;
}

double getTime() {
    unsigned LL time;
    __asm__ volatile ("rdtsc" : "=A" (time));
#ifdef LOCAL
    return time / 1.66e9; 
#else
    return time / 3.6e9;
#endif
}

#define MAXV 20000
#define MAXHV 4000
#define MAXD 60
#define MAXE (MAXV * MAXD)
#define MAXHE (MAXHV * MAXD)
#define MAXDEG 120

int N;
int E[MAXE];
int EP[MAXV + 1];
int ETris[MAXE];

int L, LX;
int S[MAXHE];
int SP[MAXHV + 1];
int STris[MAXHE];

int T[MAXHE];
int TP[MAXHV + 1];

int TV[MAXHV];
int SV[MAXHV];
int EV2[MAXV];
int SV2[MAXHV];

int order[MAXHV];
int revorder[MAXHV];
int tmporder[MAXHV];
int par[MAXHV];
int vs[MAXV];

char subConn[MAXHV * MAXHV];

int VNotUsed[MAXHV];

int Scopy[MAXHE];
int SPcopy[MAXHV + 1];

int levCounter[MAXHV];
UINT hash[MAXHV];

char fullConn[MAXV * MAXV / 8];

int verbose = 0;

LD totalTime = 0;
LD startTime;
LD bonusTime = 0;
int finished = 0;

LD FindGraphTime = 0;
LD FindGraphShuffleTime = 0;
LD ReorderEdgesTime = 0;
LD ConvertSubGraphTime = 0;
LD FindOrderTime = 0;
LD EstimateOrderTime = 0;
LD ProcessGraphTime = 0;
LD timer1 = 0;
LD timer2 = 0;
LD timer3 = 0;

int shuffleLimit;
int maxLevO;

int counter = 0;

int CUR[MAXHV];

inline int checkConnReal(int v0, int v1) {
    int x = v0 * N + v1;
    return (fullConn[x >> 3] >> (x & 7)) & 1;
}

inline int checkConnSub(int v0, int v1) {
    return subConn[v0 * L + v1];
}

VI ScoreOrder() {
    VI v;
    int LX = 0;
    REP(i, L) LX += VNotUsed[i] == 0;
    REP(i, LX) {
        int x = 0;
        FOR(j, SP[order[i]], SP[order[i] + 1]) {
            x++;
            if (!VNotUsed[S[j]] && revorder[S[j]] < i) x += 100;
        }
        v.PB(x / 100);
    }
    return v;
}

int FODist[MAXHV];
int FOPred[MAXHV];
int FOTmp[MAXHV];
int FOScore[MAXHV];
int FONext[MAXHV];
int tmpNo = 0;

LD xxxTimer = 0;

int FindOrder(int forceFirst = -1, VI scoreCutoff = VI()) {
    LD FindOrderTimer = getTime();

    REP(i, L) vs[i] = VNotUsed[i];
    REP(i, L) FOScore[i] = SP[i + 1] - SP[i];
    REP(i, L) FOR(j, SP[i], SP[i + 1]) if (!VNotUsed[S[j]]) FOScore[i] += 100;
    REP(i, L) FONext[i] = -1;
    if (forceFirst != -1) FOScore[forceFirst] += 1 << 20;
    
    int lastI = 0;
    int better = 0;
    int first = -1;
    REP(i, LX) {
    
        FOR(j, lastI, i) {
            FOR(k, SP[order[j]], SP[order[j] + 1]) if (!VNotUsed[S[k]]) {
                FOScore[S[k]] += 10000;
                if (FOScore[S[k]] < 20000) {
                    int p = S[k] - 1;
                    for (; p >= 0; p--) if (FONext[p] != -1)
                        break;
                    if (p == -1) {
                        FONext[S[k]] = first;
                        first = S[k];
                    } else {
                        FONext[S[k]] = FONext[p];
                        FONext[p] = S[k];
                    }
                }
            }
            int p = order[j] - 1;
            for (; p >= 0; p--) if (FONext[p] != -1)
                break;
            if (p == -1) {
                first = FONext[order[j]];
            } else {
                FONext[p] = FONext[order[j]];
            }
            FONext[order[j]] = -1;
            //cerr << first << " "; print(VI(FONext, FONext + L));
        }
        
        if (scoreCutoff.SZ && !better) {
            FOR(j, lastI, i) {
                int score = 0;
                FOR(k, SP[order[j]], SP[order[j] + 1])
                    if (!VNotUsed[S[k]] && vs[S[k]] && revorder[S[k]] < j) score++;
                if (score > scoreCutoff[j]) {
                    better = 1;
                    break;
                } else if (score < scoreCutoff[j]) {
                    FindOrderTime += getTime() - FindOrderTimer;
                    return 1;
                }
            }
        }
        lastI = i;
    
        int v = -1, bv = -1;
        if (first != -1) {
            int p = first;
            while (p >= 0) {
                if (FOScore[p] > bv) {
                    v = p;
                    bv = FOScore[p];
                }
                p = FONext[p];
            };
        } else {
            v = forceFirst;
            bv = FOScore[v];
            FONext[v] = -2;
            first = v;
        }
        
        if (bv < 20000) {
            queue < pair < int, PII > > q;
            REP(j, L) 
                FODist[j] = 1 << 20;
            REP(j, L) if (vs[j] && VNotUsed[j] == 0) {
                FODist[j] = 0;
                FOR(k, SP[j], SP[j + 1]) if (FODist[S[k]] > 1000) {
                    q.push(MP(S[k], MP(1, j)));
                    FODist[S[k]] = 1000;
                }
            }
            int bestA = -1, bestB = -1, bv = (1 << 20);
            while (!q.empty()) {
                pair < int, PII > p = q.front(); q.pop();
                if (bv < p.Y.X * 10000) break;
                if (p.Y.X >= FODist[p.X]) {
                    int sum = FODist[p.X] + FODist[p.Y.Y];
                    if (sum == 0 || sum * 10000 > bv + 10000) continue;
                    int av = 0;
                    tmpNo++;
                    int pos = p.X;
                    VI v;
                    while (FODist[pos] && FOTmp[pos] != tmpNo) {
                        FOR(j, SP[pos], SP[pos + 1]) if (!vs[S[j]]) av--;
                        v.PB(pos);
                        FOTmp[pos] = tmpNo;
                        pos = FOPred[pos];
                        av += 10000;
                    }
                    pos = p.Y.Y;
                    while (FODist[pos] && FOTmp[pos] != tmpNo) {
                        FOR(j, SP[pos], SP[pos + 1]) if (!vs[S[j]]) av--;
                        v.PB(pos);
                        FOTmp[pos] = tmpNo;
                        pos = FOPred[pos];
                        av += 10000;
                    }
                    if (av < bv) {
                        bv = av;
                        bestA = p.X;
                        bestB = p.Y.Y;
                        if (L > 250) break;
                    }
                } else {
                    FODist[p.X] = p.Y.X;
                    FOPred[p.X] = p.Y.Y;
                    FOR(j, SP[p.X], SP[p.X + 1]) if (VNotUsed[S[j]] == 0 && S[j] != p.Y.Y)
                        q.push(MP(S[j], MP(p.Y.X + 1, p.X)));
                }
            }
            
            if (verbose) cerr << bv << ' ' << bestA << ' ' << bestB << endl;
            
            if (bestA != -1) {
                VI vis;
                tmpNo++;
                int pos = bestA;
                while (FODist[pos] && FOTmp[pos] != tmpNo) {
                    FOTmp[pos] = tmpNo;
                    vis.PB(pos);
                    pos = FOPred[pos];
                }
                reverse(ALL(vis));
                int p = vis.SZ;
                int moveOk = FODist[pos] == 0;
                pos = bestB;
                while (FODist[pos] && FOTmp[pos] != tmpNo) {
                    FOTmp[pos] = tmpNo;
                    vis.PB(pos);
                    pos = FOPred[pos];
                }
                reverse(vis.begin() + p, vis.end());
                
                moveOk |= FODist[pos] == 0;
                if (moveOk) {
                    REP(j, vis.SZ) {
                        order[i] = vis[j];
                        revorder[vis[j]] = i;
                        vs[vis[j]] = 1;
                        i++;
                    }
                    i--;
                    continue;
                }
            }
        }
        
        order[i] = v;
        revorder[v] = i;
        vs[v] = 1;
    }
    
    FindOrderTime += getTime() - FindOrderTimer;
    //print(ScoreOrder());
    
    return 0;
}

void FindOrderAddUnused() {
    VI vnu(VNotUsed, VNotUsed + L);
    int pos = 0;
    REP(i, L) pos += vnu[i] == 0;
    while (pos < L) {
        int best = -1;
        int bv = 0;
        REP(i, L) if (vnu[i]) {
            int av = 0;
            FOR(j, SP[i], SP[i + 1]) if (!vnu[S[j]])
                av = (SP[S[j] + 1] - SP[S[j]]) * 100 + SP[i + 1] - SP[i];
            if (av > bv) {
                bv = av;
                best = i;
            }
        }
        vnu[best] = 0;
        order[pos] = best;
        revorder[best] = pos;
        pos++;
    }
}

void ReorderEdges() {
    LD ReorderEdgesTimer = getTime();

    TP[0] = 0;
    REP(i, L) {
        TP[i + 1] = TP[i] + SP[order[i] + 1] - SP[order[i]];
        FOR(j, SP[order[i]], SP[order[i] + 1])
            T[TP[i] + j - SP[order[i]]] = revorder[S[j]];
    }
    memcpy(S, T, SP[L] * sizeof (int));
    memcpy(SP, TP, (L + 1) * sizeof (int));
    
    ReorderEdgesTime += getTime() - ReorderEdgesTimer;
}


LD earlyExit;
UINT GFComp[MAXV];
UINT GFNeed[MAXHV];

#define GFIntro() \
    if ((counter & ((1 << 9) - 1)) == 0 && earlyExit && getTime() > earlyExit) { \
        REP(i, L) vs[CUR[i]] = 0; \
        REP(i, min(maxLevO, level)) FOR(j, EP[CUR[i]], EP[CUR[i] + 1]) GFComp[E[j]] = 0; \
        throw -1; \
    } \
 \
    if (DEBUG_INFO || ANALYZE_TEST) \
        levCounter[level]++; \
    ++counter; \
 \
    CUR[level] = v; \
    if (level == L - 1) throw 1; \
    int ov = level + 1; \
    vs[v] = ov; \
    int cv = CUR[par[ov]];
    
#define GFIntro2() \
    if (counter > shuffleLimit) { \
        REP(i, L) vs[CUR[i]] = 0; \
        REP(i, min(maxLevO, level)) FOR(j, EP[CUR[i]], EP[CUR[i] + 1]) GFComp[E[j]] = 0; \
        throw -1; \
    } \
 \
    if (DEBUG_INFO || ANALYZE_TEST) \
        levCounter[level]++; \
    ++counter; \
 \
    CUR[level] = v; \
    if (level == L - 1) throw 1; \
    int ov = level + 1; \
    vs[v] = ov; \
    int cv = CUR[par[ov]];
    

void goFindGraph(int v, int level) {
    GFIntro();
	
	if (level < maxLevO) {
		UINT levelBit = hash[level];
		FOR(i, EP[v], EP[v + 1])
			GFComp[E[i]] += levelBit;
	}
    
    FOR(i, EP[cv], EP[cv + 1]) {
        int nv = E[i];
        if (vs[nv] || GFComp[nv] != GFNeed[ov] || SV[ov] > EP[nv + 1] - EP[nv] || ETris[i] < STris[ov]) continue;
		/*
        int no = TP[ov + 1] - TP[ov];
        FOR(j, EP[nv], EP[nv + 1]) if (vs[E[j]]) {
            if (!checkConnSub(ov, vs[E[j]] - 1)) goto next2;
            no--;
        }
        if (no) goto next2;*/
        goFindGraph(nv, level + 1);
        next2: ;
    }
        
	if (level < maxLevO) {
		UINT levelBit = hash[level];
		FOR(i, EP[v], EP[v + 1])
			GFComp[E[i]] -= levelBit;
	}
    
    vs[v] = 0;
}

void goFindGraphComp(int v, int level) {
    GFIntro();
        
    UINT levelBit = 1U << level;
    FOR(i, EP[v], EP[v + 1])
        GFComp[E[i]] += levelBit;
    
    if (level < MAX_COMP_LEVEL) {
        FOR(i, EP[cv], EP[cv + 1]) {
            int nv = E[i];
            if (!vs[nv] && GFComp[nv] == GFNeed[ov] && ETris[i] >= STris[ov] && SV[ov] <= EP[nv + 1] - EP[nv]) 
                goFindGraphComp(nv, level + 1);
        }
    } else {
        FOR(i, EP[cv], EP[cv + 1]) {
            int nv = E[i];
            if (!vs[nv] && GFComp[nv] == GFNeed[ov] && ETris[i] >= STris[ov] && SV[ov] <= EP[nv + 1] - EP[nv]) 
                goFindGraph(nv, level + 1);
        }
    }    
        
    FOR(i, EP[v], EP[v + 1])
        GFComp[E[i]] -= levelBit;
    vs[v] = 0;
}

void goFindGraphShuffle(int v, int level) {
    GFIntro2();
    
	if (level < maxLevO) {
		UINT levelBit = hash[level];
		FOR(i, EP[v], EP[v + 1])
			GFComp[E[i]] += levelBit;
	}
    
    int edgeOrder[MAXDEG];
    int edgeNo = EP[cv + 1] - EP[cv];
    REP(i, edgeNo) edgeOrder[i] = i + EP[cv];
    
    REP(i, edgeNo) {
        swap(edgeOrder[i], edgeOrder[i + Rand() % (edgeNo - i)]);
        int nv = E[edgeOrder[i]];
        if (vs[nv] || GFComp[nv] != GFNeed[ov] || SV[ov] > EP[nv + 1] - EP[nv] || ETris[edgeOrder[i]] < STris[ov]) continue;/*
        int no = TP[ov + 1] - TP[ov];
        FOR(j, EP[nv], EP[nv + 1]) if (vs[E[j]]) {
            if (!checkConnSub(ov, vs[E[j]] - 1)) goto next2;
            no--;
        }
        if (no) goto next2;*/
        goFindGraphShuffle(nv, level + 1);
        next2: ;
    }
        
	if (level < maxLevO) {
		UINT levelBit = hash[level];
		FOR(i, EP[v], EP[v + 1])
			GFComp[E[i]] -= levelBit;
	}
    
    vs[v] = 0;
}

void goFindGraphCompShuffle(int v, int level) {
    GFIntro2();
    
    UINT levelBit = 1U << level;
    FOR(i, EP[v], EP[v + 1])
        GFComp[E[i]] += levelBit;
    
    int edgeOrder[MAXDEG];
    int edgeNo = EP[cv + 1] - EP[cv];
    REP(i, edgeNo) edgeOrder[i] = i + EP[cv];
    if (level < MAX_COMP_LEVEL) {
        REP(i, edgeNo) {
            swap(edgeOrder[i], edgeOrder[i + Rand() % (edgeNo - i)]);
            int nv = E[edgeOrder[i]];
            if (!vs[nv] && GFComp[nv] == GFNeed[ov] && ETris[edgeOrder[i]] >= STris[ov] && SV[ov] <= EP[nv + 1] - EP[nv]) 
                goFindGraphCompShuffle(nv, level + 1);
        }
    } else {
        REP(i, edgeNo) {
            swap(edgeOrder[i], edgeOrder[i + Rand() % (edgeNo - i)]);
            int nv = E[edgeOrder[i]];
            if (!vs[nv] && GFComp[nv] == GFNeed[ov] && ETris[edgeOrder[i]] >= STris[ov] && SV[ov] <= EP[nv + 1] - EP[nv]) 
                goFindGraphShuffle(nv, level + 1);
        }
    }    
        
    FOR(i, EP[v], EP[v + 1])
        GFComp[E[i]] -= levelBit;
    vs[v] = 0;
}

int EODVisited[MAXHV];
int EODPossible[MAXHV];
void EstOrderDiffRec(int v, int level) {
    EODVisited[level]++;
    if (level == EODMaxLev - 1) return;
    CUR[level] = v;
    int ov = level + 1;
    vs[v] = ov;
    int cv = CUR[par[ov]];


    UINT levelBit = 1U << level;
    FOR(i, EP[v], EP[v + 1])
        GFComp[E[i]] += levelBit;
    
    int edgeOrder[MAXDEG];
    int edgeNo = 0;
    FOR(i, EP[cv], EP[cv + 1]) {
        int nv = E[i];
        if (!vs[nv] && GFComp[nv] == GFNeed[ov] && ETris[i] >= STris[ov] && SV[ov] <= EP[nv + 1] - EP[nv])
        //if (!vs[nv] && GFComp[nv] == GFNeed[ov])
            edgeOrder[edgeNo++] = nv;    
    }
    
    EODPossible[level] += edgeNo;
    //EODPossible[level] += EP[cv + 1] - EP[cv];
    REP(i, min(edgeNo, EODLev[ov])) {
        swap(edgeOrder[i], edgeOrder[i + Rand() % (edgeNo - i)]);
        EstOrderDiffRec(edgeOrder[i], level + 1);        
    }

    FOR(i, EP[v], EP[v + 1])
        GFComp[E[i]] -= levelBit;
    vs[v] = 0;    
}

double EstOrderDiff(int no = EODLev[0]) {
    REP(i, EODMaxLev) EODVisited[i] = EODPossible[i] = 0;
    REP(i, no) EstOrderDiffRec(rand() % N, 0);
    //print(VI(EODVisited, EODVisited + EODMaxLev));
    //print(VI(EODPossible, EODPossible + EODMaxLev));
    double deg = N;
    double rv = N;
    REP(i, EODMaxLev) {
        deg *= EODPossible[i] == 0 ? 0 : (double)EODPossible[i] / EODVisited[i];
        rv += deg;
    }
    return rv;
}

void ConvertSubGraph(VI &H) {
    LD ConvertSubGraphTimer = getTime();
    memset(subConn, 0, L * L);
    int pos = 0;
    L = H[pos++];
    SP[0] = 0;
    REP(i, L) {
        int D = H[pos++];
        SP[i + 1] = SP[i] + D;
        REP(j, D) {
            int x = i * L + H[pos]; subConn[x] = 1;
            S[SP[i] + j] = H[pos++];
        }
    }
    ConvertSubGraphTime += getTime() - ConvertSubGraphTimer;
}

void DumpGraph() {
    FILE *f = fopen("graphs.txt", "a");
    fprintf(f, "%d", L); 
    REP(i, L) {
        fprintf(f, " %d", SP[i + 1] - SP[i]);
        FOR(j, SP[i], SP[i + 1])
            fprintf(f, " %d", S[j]);
    }
    fprintf(f, "\n");
    fclose(f);
}

void ProcessGraph() {
    LD ProcessGraphTimer = getTime();
    REP(i, L) levCounter[i] = 0;
    REP(i, L) par[i] = -1;
    REP(i, L) FOR(j, SP[i], SP[i + 1]) if (par[S[j]] == -1) par[S[j]] = i;
    memset(subConn, 0, L * L);
    REP(i, L) FOR(j, SP[i], SP[i + 1]) subConn[i * L + S[j]] = 1;
    
    TP[0] = 0;
    REP(i, L) {
        TP[i + 1] = TP[i];
        FOR(j, SP[i], SP[i + 1]) if (S[j] < i) T[TP[i + 1]++] = S[j];
    }
    
    REP(i, L) STris[i] = 0;
    FOR(i, 1, L) {
        int best = -1;
        int bv = -1;
        FOR(j, TP[i], TP[i + 1]) {
            int av = 0;
            FOR(k, SP[i], SP[i + 1]) if (checkConnSub(T[j], S[k])) av++;
            if (best == -1 || av > bv || av == bv && T[j] < best) {
                bv = av;
                best = T[j];
            }
        }
        par[i] = best;
        STris[i] = bv;
    }
    
    memset(vs, 0, N * sizeof(int));
    memset(GFComp, 0, N * sizeof(int));
    memset(GFNeed, 0, L * sizeof(int));
    memset(TV, 0, L * sizeof(int));
    REP(i, L) FOR(j, TP[i], TP[i + 1]) {
        if (T[j] <= 31) 
            GFNeed[i] += 1U << (T[j] % 32);
        else if (T[j] <= maxLevO - 1) 
            GFNeed[i] += hash[T[j]];
		else
            TV[i]++;
    }
    REP(i, L) 
        SV[i] = SP[i + 1] - SP[i];
		
    ProcessGraphTime += getTime() - ProcessGraphTimer;
}        


int fixedOrder[MAXV];
class SubgraphIsomorphism {public:     

int initialize(VI &  G) {
    startTime = getTime();
    
    if (DUMP_GRAPHS || ANALYZE_TEST) {
        FILE *f = fopen("graphs.txt", "w");
        fclose(f);
    }

    int pos = 0;
    N = G[pos++];
    EP[0] = 0;
    REP(i, N) {
        int D = G[pos++];
        EP[i + 1] = EP[i] + D;
        REP(j, D) {
            int x = i * N + G[pos]; fullConn[x >> 3] += 1 << (x & 7);
            //fullConn[i * N + G[pos]] = 1;
            E[EP[i] + j] = G[pos];
            pos++;
        }
    }
    
    int sumA = 0, sumB = 0;
    REP(i, N) FOR(j, EP[i], EP[i + 1]) {
        sumA += EP[i + 1] - (j + 1);
        FOR(k, j + 1, EP[i + 1]) if (checkConnReal(E[j], E[k])) {
            ETris[j]++;
            ETris[k]++;
        }
        sumB++;
    }
	
	hash[0] = 1;
	FOR(i, 1, MAXHV)
		hash[i] = hash[i - 1] * 3;
		
	   
    cerr << "Tri Prob: " << ((double)sumB / sumA * 100) << "% (" << sumB << " / " << sumA << ")" << endl;
        
#ifdef USE_FIXEDORDER
    VC < pair < int, int > > v;
    REP(i, N) v.PB(MP(EP[i + 1] - EP[i], i));
    sort(ALL(v));
    REP(i, N) fixedOrder[i] = v[i].Y;
#endif
    
    totalTime += getTime() - startTime;
    DB(totalTime);
}

VI query(VI &H) {
    startTime = getTime();
    
    ConvertSubGraph(H);

    SetSeed(L * SP[L]);
    
    REP(i, L) VNotUsed[i] = 0;
    VI degUsed(L);
    REP(i, L) degUsed[i] = SP[i + 1] - SP[i];
    
    LX = L;
    queue < int > VNUQ;
    REP(i, L) if (degUsed[i] == 1) VNUQ.push(i);
    while (!VNUQ.empty() && LX > 1) {
        int v = VNUQ.front(); VNUQ.pop();
        VNotUsed[v] = 1;
        LX--;
        FOR(j, SP[v], SP[v + 1])
            if (--degUsed[S[j]] == 1)
                VNUQ.push(S[j]);
    }
	
    

    
    if (ANALYZE_TEST && ANALYZE_TEST == L) {
        memcpy(Scopy, S, SP[L] * sizeof (int));
        memcpy(SPcopy, SP, (L + 1) * sizeof (int));
    }

    int bv = -1;
    int bestOrder = -1;
    VI vvTri(N);
    REP(i, L) if (!VNotUsed[i]) {
        int av = 0;
        FOR(j, SP[i], SP[i + 1]) FOR(k, j, SP[i + 1]) if (checkConnSub(S[j], S[k])) av++;
        vvTri[i] = av;
        av = av * 1000 + (SP[i + 1] - SP[i]);
        if (av > bv) {
            bv = av;
            bestOrder = i;
        }
    }
	
	maxLevO = max(L, 32);
    
    FindOrder(bestOrder);
    VI bestScore = ScoreOrder();
    
    if (FindOrderTime + EstimateOrderTime < MAX_ESTIMATE_TIME) {
        SetSeed(L * SP[L]);
        VC < pair < VI, int > > Vorder;
        REP(i, L) if (!VNotUsed[i]) {
            VI cutOff = Vorder.SZ ? Vorder[Vorder.SZ - 1].X : VI();
            if (!FindOrder(i, cutOff)) {
                VI score = ScoreOrder();
                int ok = 1;
                REP(j, Vorder.SZ) if (Vorder[j].X == score) {
                    ok = 0;
                    break;
                }
                if (ok) {
                    Vorder.PB(MP(ScoreOrder(), i));
                    sort(Vorder.rbegin(), Vorder.rend());
                    if (Vorder.SZ > ORDEREST_NO) Vorder.pop_back();
                }
            }
        }
        LD EstimateOrderTimer = getTime();
        memcpy(Scopy, S, SP[L] * sizeof (int));
        memcpy(SPcopy, SP, (L + 1) * sizeof (int));
        double bv = 0;
        double obv = 0;
        bestOrder = Vorder[0].Y;
        int origBestOrder = Vorder[0].Y;
        int bestPos = 0;
        REP(i, Vorder.SZ) {
            memcpy(S, Scopy, SP[L] * sizeof (int));
            memcpy(SP, SPcopy, (L + 1) * sizeof (int));
            FindOrder(Vorder[i].Y);
            FindOrderAddUnused();
            ReorderEdges();
            ProcessGraph();
            double av = EstOrderDiff();
            if (i == 0) {
                obv = av;
            }
            if (i == 0 || av < bv) {
                bestOrder = Vorder[i].Y;
                bv = av;
                bestPos = i;
            }
        }
        memcpy(S, Scopy, SP[L] * sizeof (int));
        memcpy(SP, SPcopy, (L + 1) * sizeof (int));
        
        if (bv * MINACCRATIO >= obv) {
            bestOrder = origBestOrder;
            bestPos = 0;
        } else {
            DB(obv);
            DB(bv);
            DB(bestPos);
        }
        
        FindOrder(bestOrder);        
        EstimateOrderTime += getTime() - EstimateOrderTimer;
    } else if (L < 250) {
        VI foVS(VNotUsed, VNotUsed + L);
        if (N < 3000) {
            REP(i, LX) {
                if (bestScore[i] != i) break;
                foVS[order[i]] = 1;
            }
        }
        foVS[bestOrder] = 1;
        REP(i, L) if (foVS[i] == 0 &&  (vvTri[i] || vvTri[bestOrder] == 0)) {
            REP(j, L) order[j] = -1;
            if (!FindOrder(i, bestScore)) {
                VI score = ScoreOrder();
                if (N < 3000) {
                    REP(k, LX) {
                        if (score[k] != k) break;
                        foVS[order[k]] = 1;
                    }
                }
                if (bestScore.SZ == 0 || score > bestScore) {
                    bestScore = score;
                    bestOrder = i;
                }
            }
        }
        REP(i, L) order[i] = -1;
        FindOrder(bestOrder);
    }
    
    if (ANALYZE_TEST && ANALYZE_TEST == L) {
        cerr << "Analyze Test #" << L << " LX: " << LX << " E: " << (((double)SP[L] / 2 - (L - LX)) / LX) << " BestOrder: " << bestOrder << endl << endl;
        int no = 5;
        REP(i, L) if (!VNotUsed[i]) {
            memcpy(S, Scopy, SP[L] * sizeof (int));
            memcpy(SP, SPcopy, (L + 1) * sizeof (int));
            FindOrder(i);
            VI score = ScoreOrder();
            FindOrderAddUnused();
            ReorderEdges();
            ProcessGraph();
            DumpGraph();
            REP(j, L) levCounter[j] = 0;
            int testsFull = 0;
            int resFound = 0;
            LD time = 0;
            REP(j, N) {
                LD timer = getTime();
                try {
                    earlyExit = getTime() + FULLTEST_TIME;
                    goFindGraphComp(j, 0);
                    testsFull++;
                } catch (int e) {
                    if (e != -1) {
                        REP(k, L) vs[CUR[k]] = 0;
                        REP(k, min(32, L)) FOR(l, EP[CUR[k]], EP[CUR[k] + 1]) GFComp[E[l]] = 0;
                        resFound++;
                        testsFull++;
                    }
                }
                time += getTime() - timer;
            }
            cerr << "Graph #" << (no++) << " V: " << i << " Diff: " << EstOrderDiff() << " Time: " << time << " Full: " << testsFull << " Found: " << resFound << endl;
            cerr << "Score: "; print(score);
            int sum = 0; REP(j, L) sum += levCounter[j];
            cerr << "LevCounter: (" << sum << ") [";
            REP(j, L) {
                if (j) cerr << ", ";
                cerr << j << ":" << levCounter[j];
            }
            cerr << endl;            
        }
        cerr << endl;
        return VI();
    }
    
    FindOrderAddUnused();
    ReorderEdges();
    ProcessGraph();
    
    if (DUMP_GRAPHS)
        DumpGraph();
    
    if (DEBUG_INFO)
        cerr << "#" << L << " LX: " << LX << " E: " << (((double)SP[L] / 2 - (L - LX)) / LX) << " ";
    
    int found = 0;
    int actEntry = 0;
    
    int entryUsed = 0;
    int entryUsedFull = 0;
    
    LD lastTime = 0;
    SetSeed(L * SP[L]);
    int newSeed = rand();
    int spCount = 0;
	
	shuffleLimit = 5 * L + 10;
	
#ifdef USE_FIXEDORDER
    VI nodeOrder(fixedOrder, fixedOrder + N);
#else
    VI nodeOrder(N);
    REP(i, N) nodeOrder[i] = i;
    REP(i, N) swap(nodeOrder[i], nodeOrder[i + Rand() % (N - i)]);
#endif
again:
    LD timer;
    try {
        while (1) {
            int fullMode = entryUsed >= N;
			if (fullMode) shuffleLimit += 0, counter = 0;
            int actPos = fullMode ? entryUsedFull + entryUsed % (N - entryUsedFull) : entryUsed;
            
/*            if (false && !fullMode && entryUsed > N / 2 && spCount < 10 && entryUsed == entryUsedFull + 1) {
                spCount++;
                actEntry = nodeOrder[entryUsedFull];
                fullMode = 1;
            } else {*/
                actEntry = nodeOrder[actPos];
                entryUsed++;
            //}
            timer = getTime();
            earlyExit = timer + (fullMode ? FULLTEST_TIME : ENTRYTEST_TIME);
            SetSeed(newSeed); rand(); rand(); rand(); newSeed = rand(); 
            if (EP[actEntry + 1] - EP[actEntry] >= SV[0])
                fullMode ? goFindGraphCompShuffle(actEntry, 0) : goFindGraphComp(actEntry, 0);
            fullMode ?
                FindGraphShuffleTime += getTime() - timer :
                FindGraphTime += getTime() - timer;
            swap(nodeOrder[actPos], nodeOrder[entryUsedFull++]);
            if (totalTime + getTime() - startTime >= TOTAL_TIME) break;
        }
    } catch (int e) {
        lastTime = getTime() - timer;
        if ((getTime() - startTime) + totalTime < TOTAL_TIME) {
            entryUsed > N ? 
                FindGraphShuffleTime += lastTime :
                FindGraphTime += lastTime;
            if (e == -1) goto again;
            found = e;
        }
    }
    
    if (DEBUG_INFO) {
        if (!found) cerr << "NOT ";
        cerr.precision(3);
        cerr << "P: " << actEntry << " A: " << entryUsed << " (" << (entryUsedFull) << ") T: " << (getTime() - startTime) << " LT: " << lastTime << " TT: " << (totalTime + (getTime() - startTime) - bonusTime) << endl;
    }

    if (PRINT_LEVCOUNTER) {
        int sum = 0; REP(i, L) sum += levCounter[i];
        cerr << "LevCounter: (" << sum << ") [";
        REP(i, L) {
            if (i) cerr << ", ";
            cerr << i << ":" << levCounter[i];
        }
        cerr << endl;
    }

    VI rv(L);
    REP(i, L) rv[order[i]] = CUR[i];
    
    totalTime += getTime() - startTime;
    
    if (!found || totalTime >= TOTAL_TIME || L == N) {
        DB(counter);
        DB(FindOrderTime);
        DB(FindGraphTime);
        DB(FindGraphShuffleTime);
        DB(ReorderEdgesTime);
        DB(ConvertSubGraphTime);
        DB(EstimateOrderTime);
        DB(ProcessGraphTime);
        DB(bonusTime);
        DB(totalTime);
        finished = 1;
    }
    
    return rv;
}};


#ifdef LOCAL
VI readArray() {
	int x;
	cin >> x;
	VI v;
	REP(i, x) {
		int a;
		cin >> a;
		v.PB(a);
	}
	return v;
}

int main() {
	VI G = readArray();
	SubgraphIsomorphism algo;
	
	algo.initialize(G);
	
	while (1) {
		VI H = readArray();
		VI res = algo.query(H);
		REP(i, res.SZ) cout << res[i] << "\n";
		fflush(stdout);
		if (finished) break;
	}
	
}

#endif