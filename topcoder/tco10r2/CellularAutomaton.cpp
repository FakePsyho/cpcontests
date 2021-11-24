#ifdef LOCAL
#define SCALE 1.0
#else
#define SCALE 1.0
#endif

#define MAX_TIME 9.9 * SCALE

#define THI  (N < 30 ? 4.0 : 2.75)
#define TLO  0.25 
#define TEMP_TEMPO 0.4

#define ENOUGH_STATES 0.875

#define SHOWSTATES false
#define USE_NOWRAP true
#define REMEMBER_STATES true

#define TEMP_MOD 1.05
#define GOOD_RATIO 0.75
#define FULL_SCORE_RATIO 1.5

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
#define S           size()
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

void print(VI v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]\n";}
void print(VD v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]\n";}
void print(VS v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]\n";}
void print(VVI v) {cerr << "[ ---";if (v.S) cerr << " ", print(v[0]);FOR(i, 1, v.S) cerr << " ", print(v[i]);    cerr << "--- ]\n";}
void print(PII p) {cerr << "{" << p.X << ", " << p.Y << "}";}
void print(VPII v) {cerr << "[";if (v.S) print(v[0]);FOR(i, 1, v.S)  cerr << ", ", print(v[i]);cerr << "]\n";}

template<class T> string i2s(T x) {ostringstream o; o << x; return o.str(); }
VS splt(string s, char c = ' ') {VS rv; int p = 0, np; while (np = s.find(c, p), np >= 0) {if (np != p) rv.PB(s.substr(p, np - p)); p = np + 1;} if (p < s.S) rv.PB(s.substr(p)); return rv;}

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

typedef short int gInt;
typedef gInt* pgInt;
#define OFFSET 112
#define OFFSET2 (2 * OFFSET)
gInt xSt[110 * OFFSET];
gInt xCr[110 * OFFSET];
gInt xNew[110 * OFFSET];
gInt xBest[110 * OFFSET];
int next[10 * 2];
int nextBit;
int N, K, K2;

int SX, SY, SX2, SY2;

LD startTime;
int debugOn = false;

gInt xAct1[110 * OFFSET];
gInt xAct2[110 * OFFSET];

gInt simXArray[2 * 42 * 102 * OFFSET] __attribute__((aligned(16)));
pgInt simX[22];
pgInt simXT[22];

//normal 2-tables sim
int sim() {
    gInt *xAct = xAct1;
    gInt *xActT = xAct2;
    memcpy(xAct, xNew, SY2 * OFFSET * sizeof (gInt));
    REP(turn, K) {
        memcpy(&xAct[1], &xAct[OFFSET + 1 + OFFSET * (SY - 1)], SX2 * sizeof (gInt));
        memcpy(&xAct[OFFSET + 1 + OFFSET * SY], &xAct[OFFSET + 1], SX2 * sizeof (gInt));
        FOR(y, -1, SY + 1) xAct[OFFSET + OFFSET * y] = xAct[OFFSET + 1 + OFFSET * y + SX - 1];
        FOR(y, -1, SY + 1) xAct[OFFSET + OFFSET * y + 1 + SX] = xAct[OFFSET + 1 + OFFSET * y];
        REP(y, SY) {
            gInt *pos = &xAct[OFFSET + 1 + y * OFFSET];
            gInt *xPos = &xActT[OFFSET + 1 + y * OFFSET];
            REP(x, SX) {
                int no = 0;    
                no += pos[-OFFSET - 1];
                no += pos[-OFFSET];
                no += pos[-OFFSET + 1];
                no += pos[-1];
                no += pos[1];
                no += pos[OFFSET - 1];
                no += pos[OFFSET];
                no += pos[OFFSET + 1];
                *xPos = next[no + no + pos[0]];
                pos++;
                xPos++;
            }
        }
        swap(xAct, xActT);
    }
    
    int rv = 0;
    gInt *p = &xAct[OFFSET + 1];
    REP(y, SY) {
        REP(x, SX)
            rv += *p++;
        p += OFFSET - SX;
    }
    return rv;
}

int simLocalCreate(gInt* xGrid) {
    memcpy(simXT[0], xGrid, SY2 * OFFSET * sizeof (gInt));
    FOR(turn, 1, K + 1) {
        REP(ry, SY) REP(rx, SX) {
            int no = 0;
            no += simXT[turn - 1][OFFSET + 1 + OFFSET * ((ry + SY - 1) % SY) + ((rx + SX - 1) % SX)];
            no += simXT[turn - 1][OFFSET + 1 + OFFSET * ((ry + SY - 1) % SY) + (rx)];
            no += simXT[turn - 1][OFFSET + 1 + OFFSET * ((ry + SY - 1) % SY) + ((rx + 1) % SX)];
            no += simXT[turn - 1][OFFSET + 1 + OFFSET * (ry) + ((rx + SX - 1) % SX)];
            no += simXT[turn - 1][OFFSET + 1 + OFFSET * (ry) + ((rx + 1) % SX)];
            no += simXT[turn - 1][OFFSET + 1 + OFFSET * ((ry + 1) % SY) + ((rx + SX - 1) % SX)];
            no += simXT[turn - 1][OFFSET + 1 + OFFSET * ((ry + 1) % SY) + (rx)];
            no += simXT[turn - 1][OFFSET + 1 + OFFSET * ((ry + 1) % SY) + ((rx + 1) % SX)];
            simXT[turn][OFFSET + 1 + OFFSET * (ry) + (rx)] = next[no + no + simXT[turn - 1][OFFSET + 1 + OFFSET * (ry) + (rx)]];
        }
    }
    
    REP(turn, K + 1)
        memcpy(simX[turn], simXT[turn], SY2 * OFFSET * sizeof (gInt));
    
    int rv = 0;
    REP(y, SY) REP(x, SX) rv += simXT[K][OFFSET + 1 + OFFSET * y + x];
    return rv;
}

int MODTABLE[300];
int* const MODY = &MODTABLE[25];
int* const MODX = &MODTABLE[175];
void simLocal(PII pt) {
    simXT[0][OFFSET + 1 + pt.X + OFFSET * pt.Y] = 1 - simXT[0][OFFSET + 1 + pt.X + OFFSET * pt.Y];
    FOR(turn, 1, K + 1) {
        int turnSize = turn + turn + 1;
        gInt* pos = &simXT[turn - 1][OFFSET + 1];
        int* pMY = &MODY[(pt.Y - turn)];
        REP(y, turnSize <? SY) {
            int* pMX = &MODX[(pt.X - turn)];
            REP(x, turnSize <? SX) {
                int no = 0;
                no += pos[pMY[-1] + pMX[-1]];
                no += pos[pMY[-1] + pMX[0]];
                no += pos[pMY[-1] + pMX[1]];
                no += pos[pMY[0] + pMX[-1]];
                no += pos[pMY[0] + pMX[1]];
                no += pos[pMY[1] + pMX[-1]];
                no += pos[pMY[1] + pMX[0]];
                no += pos[pMY[1] + pMX[1]];
                //simXT[turn][OFFSET + 1 + pMY[0] + pMX[0]] = next[no + no + pos[pMY[0] + pMX[0]]];
                simXT[turn][OFFSET + 1 + pMY[0] + pMX[0]] = (nextBit >> (no + no + pos[pMY[0] + pMX[0]])) & 1;
				pMX++;
            }
            pMY++;
        }
    }
}

void simXPush(PII pt) {
    REP(turn, K + 1) {
        int ry = ((pt.Y - turn + (SY << 2)) % SY);
        gInt *pX = &simX[turn][OFFSET + 1 + ry * OFFSET];
        gInt *pXT = &simXT[turn][OFFSET + 1 + ry * OFFSET];
        REP(y, (turn + turn + 1) <? SY) {
            int rx = (pt.X - turn + (SX << 2)) % SX;
            REP(x, (turn + turn + 1) <? SX) {
                pX[rx] = pXT[rx];
                rx = rx == SX - 1 ? 0 : rx + 1;
            }
			if (ry == SY - 1) {
				pX -= OFFSET * (SY - 1);
				pXT -= OFFSET * (SY - 1);
			} else {
				pX += OFFSET;
				pXT += OFFSET;
			}
			ry++;
        }
    }    
}

void simXRevert(PII pt) {
    REP(turn, K + 1) {
        int ry = ((pt.Y - turn + (SY << 2)) % SY);
        gInt *pX = &simX[turn][OFFSET + 1 + ry * OFFSET];
        gInt *pXT = &simXT[turn][OFFSET + 1 + ry * OFFSET];
        REP(y, (turn + turn + 1) <? SY) {
            int rx = (pt.X - turn + (SX << 2)) % SX;
            REP(x, (turn + turn + 1) <? SX) {
                pXT[rx] = pX[rx];
                rx = rx == SX - 1 ? 0 : rx + 1;
            }
			if (ry == SY - 1) {
				pX -= OFFSET * (SY - 1);
				pXT -= OFFSET * (SY - 1);
			} else {
				pX += OFFSET;
				pXT += OFFSET;
			}
			ry++;
        }
    }    
}

int simXFullScore() {
    int rv = 0;
    gInt *p = &simXT[K][OFFSET + 1];
    REP(y, SY) {
        rv += *p;
        int *pi = (int*)(p + 1);
        REP(x, (SX - 1) / 2)
            rv += *pi++;
        if ((SX & 1) == 0) rv += *(p + SX - 1);
        p += OFFSET;
    }
    return (rv & ((1 << 16) - 1)) + (rv >> 16);;
}

inline int simXLocalScore(gInt *xGrid, PII pt) {
    int rv = 0;
    int *pMY = &MODY[pt.Y - K];
    REP(y, (K + K + 1) <? SY) {
        gInt *pos = &xGrid[OFFSET + 1 + *pMY];
        int *pMX = &MODX[pt.X - K];
        REP(x, (K + K + 1) <? SX) {
            rv += pos[*pMX++];
        }
        pMY++;
    }
    return rv;
}

void simLocalNoWrap(PII pt) {
    simXT[0][OFFSET + 1 + pt.X + OFFSET * pt.Y] = 1 - simXT[0][OFFSET + 1 + pt.X + OFFSET * pt.Y];
    FOR(turn, 1, K + 1) {
        int turnSize = turn + turn + 1;
        gInt* pos = &simXT[turn - 1][OFFSET + 1 + (pt.Y - turn) * OFFSET + pt.X - turn];
        gInt* xPos = &simXT[turn][OFFSET + 1 + (pt.Y - turn) * OFFSET + pt.X - turn];
        REP(y, turnSize) {
            REP(x, turnSize) {
                int no = 0;
				no += pos[-OFFSET - 1];
				no += pos[-OFFSET];
				no += pos[-OFFSET + 1];
				no += pos[-1];
				no += pos[1];
				no += pos[OFFSET - 1];
				no += pos[OFFSET];
				no += pos[OFFSET + 1];
                //*xPos = next[no + no + pos[0]];
				*xPos = (nextBit >> (no + no + pos[0])) & 1;
                pos++;
                xPos++;
            }
            pos += OFFSET - turnSize;
            xPos += OFFSET - turnSize;
        }
    }
}

void simXPushNoWrap(PII pt) {
    REP(turn, K + 1) {
        int turnSize = turn + turn + 1;
        gInt* pXT = &simXT[turn][OFFSET + 1 + (pt.Y - turn) * OFFSET + pt.X - turn];
        gInt* pX  = &simX[turn][OFFSET + 1 + (pt.Y - turn) * OFFSET + pt.X - turn];
		REP(y, turnSize) {
			*pX++ = *pXT++;
			REP(x, turn) {
				*pX++ = *pXT++;
				*pX++ = *pXT++;
			}
			pX += OFFSET - turnSize;
			pXT += OFFSET - turnSize;
		}
    }    
}

void simXRevertNoWrap(PII pt) {
    REP(turn, K + 1) {
        int turnSize = turn + turn + 1;
        gInt* pXT = &simX [turn][OFFSET + 1 + (pt.Y - turn) * OFFSET + pt.X - turn];
        gInt* pX  = &simXT[turn][OFFSET + 1 + (pt.Y - turn) * OFFSET + pt.X - turn];
		REP(y, turnSize) {
			*pX++ = *pXT++;
			REP(x, turn) {
				*pX++ = *pXT++;
				*pX++ = *pXT++;
			}
			pX += OFFSET - turnSize;
			pXT += OFFSET - turnSize;
		}
    }    
}

inline int simXLocalScoreNoWrap(gInt *xGrid, PII pt) {
    int rv = 0;
    gInt *pos = &xGrid[OFFSET + 1 + pt.X - K + (pt.Y - K) * OFFSET];
    if (((int)pos) & 1) {
        REP(y, K + K + 1) {
            rv += *pos;
            int *pi = (int*)(pos + 1);
            REP(x, K) {
                rv += *pi++;
            }
            pos += OFFSET;
        }
    } else {
        REP(y, K + K + 1) {
            int *pi = (int*)pos;
            REP(x, K) {
                rv += *pi++;
            }
            rv += pos[K + K];
            pos += OFFSET;
        }
    }
    return (rv & ((1 << 16) - 1)) + (rv >> 16);;
}

int vs[101][101];
int vsx[101][101];
int vsc[101][101];
double moveProb[10000];
int moveNo[10000];
double expValues[2000];
class CellularAutomaton {public: VS configure(VS grid, string rules, int N, int K) {
    startTime = getTime();
    
    srand(1);
    
    ZERO(next);
    REP(i, rules.S) {
        if (rules[i] == '-') {
            next[i + i + 0] = 0;
            next[i + i + 1] = 0;
        } else if (rules[i] == '+') {
            next[i + i + 0] = 1;
            next[i + i + 1] = 1;
        } else if (rules[i] == '=') {
            next[i + i + 0] = 0;
            next[i + i + 1] = 1;
        } else if (rules[i] == 'X') {
            next[i + i + 0] = 1;
            next[i + i + 1] = 0;
        }
    }
	
	nextBit = 0;
	REP(i, 20) nextBit += next[i] << i;
    
    ::N = N;
    ::K = K;
    K2 = K + K + 1;
    
    ZERO(xSt);
    SX = grid[0].S;
    SY = grid.S;
    
    SX2 = SX + 4;
    SY2 = SY + 4;
    
    REP(i, 22) simX[i] = &simXArray[i * (SY + 2) * OFFSET + 7];
    REP(i, 22) simXT[i] = &simXArray[(22 + i) * (SY + 2) * OFFSET + 7];

    FOR(y, -25, SY + 26) MODY[y] = ((y + SY) % SY) * OFFSET;
    FOR(x, -25, SX + 26) MODX[x] = (x + SX) % SX;
    
    REP(y, SY) REP(x, SX) xSt[OFFSET + 1 + OFFSET * y + x] = (grid[y][x] == '1');
    //print(grid);
    memcpy(xCr, xSt, SY2 * OFFSET * sizeof (gInt));
    memcpy(xBest, xCr, SY2 * OFFSET * sizeof (gInt));
    ZERO(xAct1);
    ZERO(xAct2);
    
    debugOn = true;
    int best = simLocalCreate(xCr);
    int cur = best;
    debugOn = false;
    
    int steps = 0, accSteps = 0, dupSteps = 0, realStep = 0;
    int localSearch = 1;
    double sizeRatio = (double)(SX * SY) / ((K + K + 1) * (K + K + 1));
    //if (SX <= K + 1 || SY <= K + 1) localSearch = 0;
    if (sizeRatio < GOOD_RATIO) localSearch = 0;
    int fullScoring = !localSearch || sizeRatio <= FULL_SCORE_RATIO;
    cerr << "Size: " << (SX * SY) << " Search Size: " << ((K + K + 1) * (K + K + 1)) << " Ratio: " << (sizeRatio) << endl;
    
    VPII diff;
    
    REP(i, SX) REP(j, SY) vs[i][j] = -1;
    REP(i, SX) REP(j, SY) vsc[i][j] = -1;

    double tempMod = 1.0;
    int lastAcc = 0;
    
    int spMoveCount = 0;
    
    double curTemp = 0;
    
    int compMoves = 0;
    
    int totalCompMoves = 0;
	int expNo = 0;
    
	ZERO(expValues);
	
	PII lastChange;
	
	DB(getTime() - startTime);
	
    while (true) {
        LD curTime = getTime();
        if (curTime - startTime > MAX_TIME) break;
        
        if (steps % 1000 == 0) {
            double time = (curTime - startTime) / MAX_TIME;
            curTemp = THI + (TLO - THI) * pow(time, TEMP_TEMPO);
			expNo++;
        }
        
        steps++, lastAcc++;
                
        PII change;        
        
        int addMove = (diff.S < N);
        int spMove = 0;
        
        if (addMove) {
                while (true) {
                    int rx, ry;
					rx = rand() % SX;
					ry = rand() % SY;
                    gInt *x = localSearch == 1 ? 
                        &simX[0][OFFSET + 1 + OFFSET * ry + rx] :
                        &xCr[OFFSET + 1 + OFFSET * ry + rx];
                    if (*x != xSt[OFFSET + 1 + OFFSET * ry + rx])
                        continue;
                    change = MP(rx, ry);
                    break;
                }
        } else {
            if (compMoves >= N * ENOUGH_STATES) {
                totalCompMoves += lastAcc;
                int no = 0;
                REP(i, diff.S) if (vs[diff[i].X][diff[i].Y] == realStep) {
                    moveNo[no] = i;
					if (expValues[-vsx[diff[i].X][diff[i].Y]] < expNo)
						expValues[-vsx[diff[i].X][diff[i].Y]] = exp(vsx[diff[i].X][diff[i].Y] / curTemp) + expNo;
                    moveProb[no] = expValues[-vsx[diff[i].X][diff[i].Y]] - expNo;
                    no++;
                }
                
                //print(VD(moveProb, moveProb + no));
                
                FOR(i, 1, no) moveProb[i] += moveProb[i - 1];
                REP(i, no) moveProb[i] /= moveProb[no - 1];
                double r = randDouble();
                REP(i, no) if (r < moveProb[i]) {
                    spMove = moveNo[i] + 1;
                    break;
                }
                spMoveCount++;
            }
            int move = spMove ? spMove - 1 : rand() % diff.S;
            swap(diff[0], diff[move]);
            change = diff[0];
        }
        
        if (localSearch == 0) {
            memcpy(xNew, xCr, SY2 * OFFSET * sizeof (gInt));
            xNew[OFFSET + 1 + OFFSET * change.Y + change.X] = 1 - xNew[OFFSET + 1 + OFFSET * change.Y + change.X];
        }
        
        int needLocal = 0;
        
        int noWrap = (USE_NOWRAP && change.X - K >= 1 && change.X + K <= SX - 2 && change.Y - K >= 1 && change.Y + K <= SY - 2);
        
        int act;
        if (REMEMBER_STATES && vs[change.X][change.Y] == realStep) { 
            act = vsx[change.X][change.Y];
            needLocal = 1;
            dupSteps++;
        } else {
            if (localSearch) {
                noWrap ? simLocalNoWrap(change) : simLocal(change);
                act = fullScoring ? simXFullScore() - cur : noWrap ? 
                    simXLocalScoreNoWrap(simXT[K], change) - simXLocalScoreNoWrap(simX[K], change) :
                    simXLocalScore(simXT[K], change) - simXLocalScore(simX[K], change);
            } else {
                act = sim() - cur;
            }
            vs[change.X][change.Y] = realStep;
            vsx[change.X][change.Y] = act;
        }
		compMoves += vsc[change.X][change.Y] != accSteps;
		vsc[change.X][change.Y] = accSteps;
        
        double temp = curTemp;
        //if (!addMove) temp *= 1.5;
		
		if (act < 0 && expValues[-act] < expNo)
			expValues[-act] = expNo + exp(act / temp);
		
        
        if (spMove || act >= 0 || randDouble() < expValues[-act] - expNo) {
            lastAcc = 0;
            
            if (localSearch == 1 && needLocal == 1) {
                noWrap ? simLocalNoWrap(change) : simLocal(change);
                dupSteps--;
            }
            accSteps++;
			
			lastChange = change;
                
            if (addMove) {    
                diff.PB(change);
            } else {
                swap(diff[0], diff[diff.S - 1]), diff.pop_back();
            }
            
            cur += act;
            if (localSearch == 0)
                memcpy(xCr, xNew, SY2 * OFFSET * sizeof (gInt));
            if (localSearch == 1) {
                noWrap ? simXPushNoWrap(change) : simXPush(change);
            }
            
            if (cur > best) {
                best = cur;
                if (localSearch == 0)
                    memcpy(xBest, xNew, SY2 * OFFSET * sizeof (gInt));
                if (localSearch == 1)
                    memcpy(xBest, simXT[0], SY2 * OFFSET * sizeof (gInt));
                if (SHOWSTATES)
                    cerr << "At: " << steps << " Score: " << best << " (" << ((double)best / (SX * SY)) << ")" << endl;
            }
			
			if ((double)accSteps / (steps - accSteps + 1) < 0.05 * sizeRatio) {
				FOR(x, change.X - K - K, change.X + K + K + 1) {
					int rx = (x + (SX << 2)) % SX;
					int ry = (change.Y - K - K + (SY << 2)) % SY;
					REP(y, K + K + K + K + 1) {
						vs[rx][ry]--;
						ry = ry == SY - 1 ? 0 : ry + 1;
					}
				}
			} else {
				realStep++;
			}			
			compMoves = 0;
            
        } else {
            if (localSearch == 1 && needLocal == 0)
                noWrap ? simXRevert(change) : simXRevert(change);
        }
    }
    
    
    cerr << "Steps: " << (steps - dupSteps) << " : " << accSteps << " : " << dupSteps << " : " << ((double)accSteps / steps) << endl;
    DB(best);
    DB(tempMod);
    DB(spMoveCount);
    DB((double)totalCompMoves / spMoveCount);
    //print(accDist);
    
    VS res;
    REP(y, SY) {
        string s;
        REP(x, SX) s += xBest[OFFSET + 1 + OFFSET * y + x] ? '1' : '0';
        res.PB(s);
    }
    return res;
}};


#ifdef LOCAL
int main() {
    int R, N, K;
	cin >> R;
	VS grid(R);
    for (int i=0; i<R; i++)
        cin >> grid[i];
	string rules;
	cin >> rules >> N >> K;
	CellularAutomaton algo;
    VS res = algo.configure(grid, rules, N, K);
    for (int i=0; i<R; i++)
        cout << res[i] << endl;
    fflush(stdout);
}
#endif