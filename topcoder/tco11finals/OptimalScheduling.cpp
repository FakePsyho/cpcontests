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

#ifdef TEST
#define cerr cout
#endif

#ifdef LOCAL
#define TIMESCALE 0.5
#else
#define TIMESCALE 1.0
#endif

#define MAXTIME (9.8 * TIMESCALE)
#define PART1 (MAXTIME * 0.9)
#define PART2 (MAXTIME * 0.1)

using namespace std;

#define FOR(i,a,b)  for(int i=(int)(a);i<(int)(b);++i)
#define REP(i,a)    FOR(i,0,a)
#define ZERO(m)     memset(m,0,sizeof(m))
#define ALL(x)      x.begin(),x.end()
#define PB          push_back
#define S           size()
#define LL          long long
#define MP          make_pair
#define X           first
#define Y           second
#define VC          vector
#define VD            VC < double >
#define VVD            VC < VD >
#define PII         pair <int, int>
#define VI          VC < int >
#define VVI            VC < VI >
#define VPII        VC < PII >
#define VS          VC<string>
#define DB(a)        cerr << #a << ": " << a << endl;

void print(VI v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]\n";}
void print(VD v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]\n";}
void print(VS v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]\n";}
template<class T> string i2s(T x) {ostringstream o; o << x; return o.str(); }
VS splt(string s, char c = ' ') {VS rv; int p = 0, np; while (np = s.find(c, p), np >= 0) {if (np != p) rv.PB(s.substr(p, np - p)); p = np + 1;} if (p < (int)s.S) rv.PB(s.substr(p)); return rv;}

double getTime() {
    unsigned LL time;
    __asm__ volatile ("rdtsc" : "=A" (time));
#ifdef LOCAL
    return time / 1.66e9; 
#else
    return time / 3.6e9;
#endif
}

double randDouble() {
    return ((double)(rand() + 0.5) / ((double)RAND_MAX + 1));
}

#define MAXN 64
#define MAXG 100

int N, M;
int NUM[MAXN];
int RANK[MAXN];
int AGE[MAXN];
int W[7];
int G;
int SHOW = 0;

double startTime;
int TG[MAXN];
int CG[MAXN];
int sol[MAXG][6];
int cur[MAXG][6];

void CreateRandom() {
    int gtries = 0;
    int pos = 0;
    
    VI order1, order2, order; 
    //REP(i, N) if (TG[i] == M + 1) order.PB(i);
    //REP(i, N) if (TG[i] == M) order.PB(i);
	
	REP(i, N) REP(j, TG[i]) order.PB(i);
    
    REP(i, order.S) sol[pos / 6][pos % 6] = order[i], pos++;
	bool done = false;
    while (1) {
		done = true;
        REP(i, G) {
            REP(j, 6) REP(k, j) if (sol[i][j] == sol[i][k]) {
                swap(sol[i][j], sol[rand() % G][rand() % 6]);
                done = false;
            }
        }
		if (done) break;
    }
}

int TPos[MAXN][6];
double scorePosition() {
    double SAss = 0;
    REP(i, N) {
        int a = TPos[i][0] + TPos[i][1] + TPos[i][2];
        int b = M - a;
        SAss += abs(a - b);
    }
        
    double SPos = 0;
    REP(i, N) {
        double avg = M / 6.0;
        double x = 0;
        REP(j, 6) x += (TPos[i][j] - avg) * (TPos[i][j] - avg);
        SPos += sqrt(x / 6);
    }
    
    double v = 0;
    v += SAss * W[5];
    v += SPos * W[6];
    if (SHOW) {
        DB(SAss);
        DB(SPos);
    }
    return v;
}

double scorePosition(int p) {
    double SAss = 0;
    {
        int a = TPos[p][0] + TPos[p][1] + TPos[p][2];
        int b = M - a;
        SAss += abs(a - b);
    }
        
    double SPos = 0;
    {
        double avg = M / 6.0;
        double x = 0;
        REP(j, 6) x += (TPos[p][j] - avg) * (TPos[p][j] - avg);
        SPos += sqrt(x / 6);
    }
    
    double v = 0;
    v += SAss * W[5];
    v += SPos * W[6];
    if (SHOW) {
        DB(SAss);
        DB(SPos);
    }
    return v;
}

short TAll[MAXN][MAXN];
short TOpp[MAXN][MAXN];
int TAllNo[MAXN];
int TOppNo[MAXN];

int TUsed[MAXN];
int TUsedNo;
int TCH[12];
int TCHNo = 0;

int TBonus[6];
int TBonusNo;
int BonusT[MAXN];

int GPos[MAXN][16];
int XPos[MAXN][MAXG];
int GPosNo[MAXN];

bool isBonusGame(int a, int b) {
	return BonusT[cur[a][b]] && XPos[cur[a][b]][a] == 2;
}

double scoreOrder() {
    int LP[MAXN];
    double IDT[MAXN];
    REP(i, N) IDT[i] = (double)G / TG[i] - 1;
    REP(i, N) LP[i] = -1;
    
    double av = 0;
    REP(i, G) {
        REP(j, 6) {
            int x = cur[i][j];
            if (LP[x] != -1)
                av += abs(IDT[x] - (i - LP[x] - 1));
            LP[x] = i;
        }
    }
    return av * W[4];
}

bool allDiff(int a, int b) {
	TUsedNo++;
	REP(i, 6)
        TUsed[cur[a][i]] = TUsedNo;
	REP(i, 6)
		if (TUsed[cur[b][i]] == TUsedNo) return false;
	return true;
	
}

void initData() {
    ZERO(TUsed);
    TUsedNo = 1;
    
    ZERO(TAll);
    ZERO(TOpp);
    REP(i, G) {
        REP(j, 3) REP(k, 3) if (j != k) TAll[cur[i][j]][cur[i][k]]++;
        REP(j, 3) REP(k, 3) if (j != k) TAll[cur[i][j+3]][cur[i][k+3]]++;
        REP(j, 3) REP(k, 3) TOpp[cur[i][j]][cur[i][k+3]]++;
        REP(j, 3) REP(k, 3) TOpp[cur[i][j+3]][cur[i][k]]++;
    }
	
	REP(i, N) TOpp[i][i] = TAll[i][i] = -1;
    
    ZERO(TAllNo);
    ZERO(TOppNo);
    REP(i, N) REP(j, N) if (i != j) TAllNo[i] += TAll[i][j] > 0;
    REP(i, N) REP(j, N) if (i != j) TOppNo[i] += TOpp[i][j] > 0;
    
    ZERO(GPosNo);
    REP(i, G) REP(j, 6) {
        int x = cur[i][j];
        XPos[x][i] = GPosNo[x];
        GPos[x][GPosNo[x]++] = i;
    }
    REP(i, N) GPos[i][GPosNo[i]] = 1 << 10;
}

void UpdateRem(int g, int p) {
    if (p < 3) {
        REP(i, 3) if (i != p) {
            --TAll[cur[g][i]][cur[g][p]];
            if (--TAll[cur[g][p]][cur[g][i]] == 0) {
                TAllNo[cur[g][p]]--;
                TAllNo[cur[g][i]]--;
            }
        }
        FOR(i, 3, 6) {
            --TOpp[cur[g][i]][cur[g][p]];
            if (--TOpp[cur[g][p]][cur[g][i]] == 0) {
                TOppNo[cur[g][p]]--;
                TOppNo[cur[g][i]]--;
            }
        }
    } else {
        FOR(i, 3, 6) if (i != p) {
            --TAll[cur[g][i]][cur[g][p]];
            if (--TAll[cur[g][p]][cur[g][i]] == 0) {
                TAllNo[cur[g][p]]--;
                TAllNo[cur[g][i]]--;
            }
        }
        REP(i, 3) {
            --TOpp[cur[g][i]][cur[g][p]];
            if (--TOpp[cur[g][p]][cur[g][i]] == 0) {
                TOppNo[cur[g][p]]--;
                TOppNo[cur[g][i]]--;
            }
        }
    }
}

void UpdateAdd(int g, int p) {
    if (p < 3) {
        REP(i, 3) if (i != p) {
            ++TAll[cur[g][i]][cur[g][p]];
            if (++TAll[cur[g][p]][cur[g][i]] == 1) {
                TAllNo[cur[g][p]]++;
                TAllNo[cur[g][i]]++;
            }
        }
        FOR(i, 3, 6) {
            ++TOpp[cur[g][i]][cur[g][p]];
            if (++TOpp[cur[g][p]][cur[g][i]] == 1) {
                TOppNo[cur[g][p]]++;
                TOppNo[cur[g][i]]++;
            }
        }
    } else {
        FOR(i, 3, 6) if (i != p) {
            ++TAll[cur[g][i]][cur[g][p]];
            if (++TAll[cur[g][p]][cur[g][i]] == 1) {
                TAllNo[cur[g][p]]++;
                TAllNo[cur[g][i]]++;
            }
        }
        REP(i, 3) {
            ++TOpp[cur[g][i]][cur[g][p]];
            if (++TOpp[cur[g][p]][cur[g][i]] == 1) {
                TOppNo[cur[g][p]]++;
                TOppNo[cur[g][i]]++;
            }
        }
    }
}

int URA = 0;
int URO = 0;

void UpdateRemX(int g, int p) {
    int *x = cur[g];
    short *TAllP = TAll[x[p]];
    short *TOppP = TOpp[x[p]];
    if (p < 3) {
		URA -= TAllP[x[0]] == 1;
		URA -= TAllP[x[1]] == 1;
		URA -= TAllP[x[2]] == 1;
		URO -= TOppP[x[3]] == 1;
		URO -= TOppP[x[4]] == 1;
		URO -= TOppP[x[5]] == 1;
    } else {
		URO -= TOppP[x[0]] == 1;
		URO -= TOppP[x[1]] == 1;
		URO -= TOppP[x[2]] == 1;
		URA -= TAllP[x[3]] == 1;
		URA -= TAllP[x[4]] == 1;
		URA -= TAllP[x[5]] == 1;
    }
}

void UpdateAddX(int g, int p) {
    int *x = cur[g];
    short *TAllP = TAll[x[p]];
    short *TOppP = TOpp[x[p]];
    if (p < 3) {
		URA += TAllP[x[0]] == 0;
		URA += TAllP[x[1]] == 0;
		URA += TAllP[x[2]] == 0;
		URO += TOppP[x[3]] == 0;
		URO += TOppP[x[4]] == 0;
		URO += TOppP[x[5]] == 0;
    } else {
		URO += TOppP[x[0]] == 0;
		URO += TOppP[x[1]] == 0;
		URO += TOppP[x[2]] == 0;
		URA += TAllP[x[3]] == 0;
		URA += TAllP[x[4]] == 0;
		URA += TAllP[x[5]] == 0;
    }
}

inline void UpdatePos(int t, int a, int b) {
    int p = XPos[t][a];
    XPos[t][b] = p;
    GPos[t][p] = b;
}

double scoreOrder(int t) {
    double idt = (double)G / TG[t] - 1;
    double v = 0;
    FOR(j, 1, TG[t])
        v += abs(idt - (GPos[t][j] - GPos[t][j-1] - 1));
    return v;
}

inline double scoreOrder(int t, int p) {
    double idt = (double)G / TG[t] - 1;
    double v = 0;
    if (p) v += abs(idt - (GPos[t][p] - GPos[t][p - 1] - 1));
    if (p < TG[t] - 1) v += abs(idt - (GPos[t][p + 1] - GPos[t][p] - 1));
    return v;
}

double scoreDeltaX(int a, int b) {
    double v = 0;
    
    double SAge = 0;
    SAge += abs(AGE[cur[a][0]] + AGE[cur[a][1]] + AGE[cur[a][2]] - AGE[cur[a][3]] - AGE[cur[a][4]] - AGE[cur[a][5]]);
    if (a != b) SAge += abs(AGE[cur[b][0]] + AGE[cur[b][1]] + AGE[cur[b][2]] - AGE[cur[b][3]] - AGE[cur[b][4]] - AGE[cur[b][5]]);
    SAge /= 3;
    
    double SRank = 0;
    SRank += abs(RANK[cur[a][0]] + RANK[cur[a][1]] + RANK[cur[a][2]] - RANK[cur[a][3]] - RANK[cur[a][4]] - RANK[cur[a][5]]);
    if (a != b) SRank += abs(RANK[cur[b][0]] + RANK[cur[b][1]] + RANK[cur[b][2]] - RANK[cur[b][3]] - RANK[cur[b][4]] - RANK[cur[b][5]]);
    SRank /= 3;
    v += SAge * W[0];
    v += SRank * W[1];
    
    return v;
}

int bonusExist() {
    if (TBonusNo > 1) {
        REP(i, TBonusNo) REP(j, i) if (GPos[TBonus[i]][2] == GPos[TBonus[j]][2]) return true;
    } 
	return false;    
}

double scoreAgeAndRank() {
	double v = 0;
	int x = TBonusNo;
	TBonusNo = 0;
	REP(i, G) v += scoreDeltaX(i, i);
	TBonusNo = x;
	return v;
}

void swapTeams(int a, int b) {
	REP(i, G) REP(j, 6) {
		if (cur[i][j] == a) 
			cur[i][j] = b;
		else if (cur[i][j] == b)
			cur[i][j] = a;
	}
	
	REP(i, TG[a]) swap(GPos[a][i], GPos[b][i]);
	REP(i, G) swap(XPos[a][i], XPos[b][i]);
}

int swapTeams() {
	double ov = scoreAgeAndRank();
	double bv = ov;
	
	int t0 = -1, t1 = -1;
	VI perm; REP(i, N) perm.PB(i);
	random_shuffle(ALL(perm));
	REP(p0, N) REP(p1, p0) {
		int i = perm[p0];
		int j = perm[p1];
		if (TG[i] != TG[j]) continue;
		
		swap(AGE[i], AGE[j]);
		swap(RANK[i], RANK[j]);
		double v = scoreAgeAndRank();
		swap(AGE[i], AGE[j]);
		swap(RANK[i], RANK[j]);
		
		if (v < bv - 1e-4) {
			t0 = i;
			t1 = j;
			bv = v;
			goto out;
		}
	}
out: ;
	
	if (t0 != -1) {
		swapTeams(t0, t1);
		return true;
	} else
		return false;
}

double scoreMatch() {
    //initData();

    double SAge = 0;
    REP(i, G)
        SAge += abs(AGE[cur[i][0]] + AGE[cur[i][1]] + AGE[cur[i][2]] - AGE[cur[i][3]] - AGE[cur[i][4]] - AGE[cur[i][5]]);
    SAge /= 3;
    
    double SRank = 0;
    REP(i, G)
        SRank += abs(RANK[cur[i][0]] + RANK[cur[i][1]] + RANK[cur[i][2]] - RANK[cur[i][3]] - RANK[cur[i][4]] - RANK[cur[i][5]]);
    SRank /= 3;
    
    double SAll = 2 * M * N;
    REP(i, N) SAll -= TAllNo[i];
    double SOpp = 3 * M * N;
    REP(i, N) SOpp -= TOppNo[i];
    
    double v = 0;
    v += SAge * W[0];
    v += SRank * W[1];
    v += SAll * W[2];
    v += SOpp * W[3];
    
    if (SHOW) {
        DB(SAge);
        DB(SRank);
        DB(SAll);
        DB(SOpp);
    }
    
    return v + scoreOrder();
}

bool inRange(int tg, int p, int g) {
	int lo = p * G / tg;
	int hi = ((p + 1) * G + tg - 1) / tg;
	return g >= lo && g <= hi;
}

void SAMatch(double totalTime) {
    double time = getTime();
    
    REP(i, G) REP(j, 6) cur[i][j] = sol[i][j];
    
    ZERO(TPos);
    REP(i, G) REP(j, 6) if (!isBonusGame(i, j)) TPos[cur[i][j]][j]++;
	
    initData();
    DB(scoreMatch());
    double bv = scoreMatch(); //scoreRest();
    double av = bv;
    DB(bv);
    
    int steps0 = 0;
    int steps1 = 0;
    int steps2 = 0;
    int steps3 = 0;
    
    double curTime = 0;
    double LTemp = 0;
	double temp = 0;
	
	int xyz0 = 0;
	int xyz1 = 0;
	
    while (1) {
        if ((steps0 & 1023) == 0) {
            curTime = ((getTime() - time) / totalTime);
            if (curTime >= 1.0) break;
            LTemp = 300 + 2500 * (1.0 - curTime);
			temp = 100 + bv / 900 * (1.0 - curTime);
        }
		
        steps0++;
    
        int a = rand() % G;
        int b = a - 5 + rand() % 11;
        if (b < 0 || b >= G) continue;
		
        int c = rand() % 6;
        int d = rand() % 6;
				
        if (a == b && c / 3 == d / 3) c = (c + 3) % 6;
        int t0 = cur[a][c];
        int t1 = cur[b][d];
        if (t0 == t1) continue;
        
        if (a != b && !allDiff(a, b)) continue;
        
		int p0 = XPos[t0][a];
		int p1 = XPos[t1][b];
		
        steps1++;
        
        double valueOrder = 0;
        if (a != b) {
            valueOrder -= scoreOrder(t0, p0) + scoreOrder(t1, p1);
            GPos[t0][p0] = b;
            GPos[t1][p1] = a;
            valueOrder += scoreOrder(t0, p0) + scoreOrder(t1, p1);
            GPos[t0][p0] = a;
            GPos[t1][p1] = b;
            if (valueOrder * W[4] > W[4] + LTemp / 2) continue;
        }
        
        steps2++;
                
        double v = 0;
        
		URO = 0;
		URA = 0;
		//v -= scoreDeltaX(a, b);
		v -= bonusExist() * 50000;
		UpdateRemX(a, c);
		UpdateRemX(b, d);
		TOpp[t0][t1]--;
		TOpp[t1][t0]--;
		swap(cur[a][c], cur[b][d]);
		UpdateAddX(a, c);
		UpdateAddX(b, d);
		TOpp[t0][t1]++;
		TOpp[t1][t0]++;
		GPos[t0][p0] = b;
		GPos[t1][p1] = a;
		//v += scoreDeltaX(a, b);
		v += bonusExist() * 50000;
		GPos[t0][p0] = a;
		GPos[t1][p1] = b;
		v -= 2 * URA * W[2];
		v -= 2 * URO * W[3];
		v += valueOrder * W[4];
		
        //if (v <= 0 || randDouble() < exp((-v) / temp)) {
		if (v <= rand() % (int)LTemp) {
			swap(cur[a][c], cur[b][d]);
			UpdateRem(a, c);
			UpdateRem(b, d);
			swap(cur[a][c], cur[b][d]);
			UpdateAdd(a, c);
			UpdateAdd(b, d);
			if (a != b) {
                UpdatePos(t0, a, b);
                UpdatePos(t1, b, a);
            }
            
            steps3++;
            av += v;
            
            if (av < bv - 1e-9) {
                bv = av;
                REP(i, G) REP(j, 6) sol[i][j] = cur[i][j];
            }
        } else {
            swap(cur[a][c], cur[b][d]);
        }
    }
	
    DB(steps0);
    DB(steps1);
    DB(steps2);
    DB(steps3);
	
	DB(xyz0);
	DB(xyz1);
    
    REP(i, G) REP(j, 6) cur[i][j] = sol[i][j];
    DB(bv);
    initData();
    SHOW = 1; scoreMatch(); scorePosition(); SHOW = 0;
	
	DB(scoreMatch());
	while (swapTeams()); 
	REP(i, G) REP(j, 6) sol[i][j] = cur[i][j];
	initData();
	DB(scoreMatch());
}

void SAPosition(double totalTime) {
    double time = getTime();
	
	if (totalTime < 0) totalTime = 0.1;
    
    REP(i, G) REP(j, 6) cur[i][j] = sol[i][j];
	
    ZERO(TPos);
    REP(i, G) REP(j, 6) if (!isBonusGame(i, j)) TPos[cur[i][j]][j]++;
    
    double bv = scorePosition();
    double av = bv;
    
    double curTime = 0;
    double LTemp = 0;
	double temp = 0;
	
    DB(bv);
    int steps = 0;
    while (1) {
        if ((steps & 63) == 0) {
            curTime = ((getTime() - time) / totalTime);
            if (curTime >= 1.0) break;
            LTemp = (min(W[5], W[6]) + (W[5] + W[6]) * (1.0 - curTime)) / 4;
			temp = 100 + bv / 900 * (1.0 - curTime);
        }
        
        int a = rand() % G;
        int b = rand() % 6;
        int c = (b/3)*3 + rand() % 3;
        int t = rand() % 10;
        
        steps++;
        
		if (t && c == b) continue;
		
		double v = 0;
		
        if (t == 0) {
			REP(i, 6) v -= scorePosition(cur[a][i]);
			REP(i, 6) if (!isBonusGame(a, i)) TPos[cur[a][i]][i]--;
			REP(i, 6) if (!isBonusGame(a, i)) TPos[cur[a][i]][i%3+(i<3?3:0)]++;
			REP(i, 6) v += scorePosition(cur[a][i]);
        } else {
			v -= scorePosition(cur[a][b]);
			v -= scorePosition(cur[a][c]);
			if (!isBonusGame(a, b)) {
				TPos[cur[a][b]][b]--;
				TPos[cur[a][b]][c]++;
			}
			if (!isBonusGame(a, c)) {
				TPos[cur[a][c]][c]--;
				TPos[cur[a][c]][b]++;
			}
			v += scorePosition(cur[a][b]);
			v += scorePosition(cur[a][c]);
        }
		
        if (v <= rand() % (int)LTemp) {
		//if (v <= 0) {
			if (t == 0) {
				swap(cur[a][0], cur[a][3]);
				swap(cur[a][1], cur[a][4]);
				swap(cur[a][2], cur[a][5]);
			} else {
				swap(cur[a][b], cur[a][c]);
			}
			
            av += v;
            if (av < bv) {
                bv = av;
                REP(i, G) REP(j, 6) sol[i][j] = cur[i][j];
            }
        } else {
			if (t == 0) {
				REP(i, 6) if (!isBonusGame(a, i)) TPos[cur[a][i]][i]++;
				REP(i, 6) if (!isBonusGame(a, i)) TPos[cur[a][i]][i%3+(i<3?3:0)]--;
			} else {
				if (!isBonusGame(a, b)) {
					TPos[cur[a][b]][b]++;
					TPos[cur[a][b]][c]--;
				}
				if (!isBonusGame(a, c)) {
					TPos[cur[a][c]][c]++;
					TPos[cur[a][c]][b]--;
				}
			}
			
        }
    }
    DB(bv);
    DB(steps);
}

class OptimalScheduling {public: VS createSchedule(int N, int M, VS Z, VI W, VI F) {
    startTime = getTime();
    
    ::N = N;
    ::M = M;
    REP(i, N) sscanf(Z[i].c_str(), "%d %d %d", &NUM[i], &AGE[i], &RANK[i]);
    REP(i, 7) ::W[i] = W[i];
    G = (N * M + 5) / 6;
    
    srand(1);
    REP(i, N) TG[i] = M;
    print(F);
    TBonusNo = 0;
	ZERO(BonusT);
    REP(i, N) REP(j, F.S) if (NUM[i] == F[j]) {
        TG[i]++;
        TBonus[TBonusNo++] = i;
		BonusT[i] = 1;
    }
    
    cerr << "N: " << N <<  " M: " << M << " G: " << G << endl;
    
    CreateRandom();
    SAMatch(PART1);
    SAPosition((MAXTIME - (getTime() - startTime)));
	
    cerr << "Time: " << (getTime() - startTime) << endl;
	//REP(i, 10) swapTeams();
    //cerr << "Time: " << (getTime() - startTime) << endl;
    
	
    VS rv;
    REP(i, G) REP(j, 6) sol[i][j] = NUM[sol[i][j]];
    REP(i, G)
        rv.PB(i2s(sol[i][0]) + " " + i2s(sol[i][1]) + " " + i2s(sol[i][2]) + " : " + i2s(sol[i][3]) + " " + i2s(sol[i][4]) + " " + i2s(sol[i][5]));
    return rv;    
}};


#ifdef LOCAL
int main() {
	int N, M;
	cin >> N >> M;
	
	int ZNO; cin >> ZNO;
	VS Z(ZNO);
	REP(i, ZNO) {
		string a, b, c;
		cin >> a >> b >> c;
		Z[i] = a + " " + b + " " + c;
	}
	
	int WNO; cin >> WNO;
	VI W(WNO);
	REP(i, WNO) cin >> W[i];
	
	int SNO; cin >> SNO;
	VI SZ(SNO);
	REP(i, SNO) cin >> SZ[i];

	//DB(N); DB(M); DB(ZNO); DB(WNO); DB(SNO);
	
	OptimalScheduling algo;
	VS rv = algo.createSchedule(N, M, Z, W, SZ);

	cout << rv.S << endl;
	REP(i, rv.S) cout << rv[i] << endl;
	//REP(i, rv.S) DB(rv[i]);
    fflush(stdout);

	return 0;
}
#endif
