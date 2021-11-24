#define TIME_LIMIT 9.8
#define PART0_LIMIT (TIME_LIMIT * 0.05)
#define PART1_LIMIT (TIME_LIMIT * 0.75)
#define PART2_LIMIT (TIME_LIMIT * 0.95)

#define FINDALL_ESTIMATE 0.01
#define FINDALL_TIMELIMIT 3.0

#define FASTATAN_SIZE 400

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
#include <ext/hash_set>

using namespace __gnu_cxx;
using namespace std;

#define FOR(i,a,b)  for(int i=(a);i<(b);++i)
#define REP(i,a)    FOR(i,0,a)
#define ZERO(m)     memset(m,0,sizeof(m))
#define ALL(x)      x.begin(),x.end()
#define PB          push_back
#define S           size()
#define LL          long long
#define LD          long double
#define MP          make_pair
#define X           first
#define Y           second
#define VC          vector
#define PII         pair <int, int>
#define PDD			pair <double, double>
#define VI          VC<int>
#define VVI			VC < VI >
#define VD			VC < double >
#define VPII        VC < PII >
#define VS          VC<string>
#define DB(a)        cerr << #a << ": " << a << endl;

void print(VI v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]\n";}
void print(VD v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]\n";}
void print(VS v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]\n";}
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

#define MAXN 5000
#define MAXSIZE 20

double randDouble() {
    return ((double)(rand() + 0.5) / ((double)RAND_MAX + 1));
}

int N;
int NOk;
PII P[MAXN];

int biggestPoly = 5;

LL SD, SV;
LL RD, RV;

double PI = 2 * acos(0);

double startTime;
double debugTimer = 0;

int D[MAXN][MAXN];

double FT[FASTATAN_SIZE * FASTATAN_SIZE * 2 + 3 * FASTATAN_SIZE];

bool POk[MAXN];

VVI polys;
set<int> polysHash;

int NOW = 0;

VPII BP;

double timer = 0;
double timer2 = 0;

void printTime() {
	cerr << "Current Time: " << (getTime() - startTime) << endl;
}

int hashPoly(VI v) {
	sort(ALL(v));
	int h = 0;
	REP(i, v.S) h = h * 1337 + v[i];
	h = h * 137;
	h &= ~127;
	h += v.S;
	return h;
}

void addPoly(VI &v) {
	int h = hashPoly(v);
	if (!polysHash.count(h)) {
		polysHash.insert(h);
		polys.PB(v);
	}
}

inline int dist(int x1, int y1, int x2, int y2) {
	return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}

inline int dist(int a, int b) {
	return dist(P[a].X, P[a].Y, P[b].X, P[b].Y);
}

inline PDD center(VI &v) {
	PII sum;
	REP(i, v.S) sum.X += P[v[i]].X, sum.Y += P[v[i]].Y;
	return MP((double)sum.X / v.S, (double)sum.Y / v.S);
}

inline int crossSign(int x0, int y0, int x1, int y1, int x2, int y2) {
	int v = (x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - x1);
	return v == 0 ? 0 : v > 0 ? 1 : -1;
}

inline int crossSign(int p0, int p1, int p2) {
	return crossSign(P[p0].X, P[p0].Y, P[p1].X, P[p1].Y, P[p2].X, P[p2].Y);
}

inline double angle(int x0, int y0, int x1, int y1, int x2, int y2) {
	int v = (x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - x1);
	LL d2 = (LL)(x1 - x0) * (y1 - y0) * (x2 - x1) * (y2 - y1);
	return 1.0 - asin((double)v * v / d2) / (2 * PI);
}

inline double angle(int p0, int p1, int p2) {
	return angle(P[p0].X, P[p0].Y, P[p1].X, P[p1].Y, P[p2].X, P[p2].Y);
}

inline double fastatan(double v) {
	if (v < -FASTATAN_SIZE) 
		v = FASTATAN_SIZE - 1 - 1 / v;
	else if (v > FASTATAN_SIZE)
		v = FASTATAN_SIZE + 1 - 1 / v;
	v *= FASTATAN_SIZE;
	int vi = (int)v;
	double d = v - vi;
	int p = FASTATAN_SIZE * FASTATAN_SIZE + FASTATAN_SIZE + vi;
	return FT[p] * (1 - d) + FT[p] * d;
}

void initfastatan() {
	int o = FASTATAN_SIZE * FASTATAN_SIZE + FASTATAN_SIZE;
	FOR(i, -o, o)
		FT[i + o] = atan(i / (double)FASTATAN_SIZE);
}

inline double fastatan2(double y, double x) {
	if (x > 0) {
		return fastatan(y / x);
	} else if (x < 0) {
		return y >= 0 ? PI + fastatan(y / x) : -PI + fastatan(y / x);
	} else 
		return y > 0 ? PI / 2 : -PI / 2;
}

bool isConvex(VI &v) {
	double sum = 0;
	int sign = 0;
	
	REP(i, v.S) {
		double a0 = atan2(P[v[(i+1)%v.S]].Y - P[v[i]].Y, P[v[(i+1)%v.S]].X - P[v[i]].X);
		double a1 = atan2(P[v[(i+v.S-1)%v.S]].Y - P[v[i]].Y, P[v[(i+v.S-1)%v.S]].X - P[v[i]].X);
		double a = a1 - a0;
		while (a < 0) a += 2 * PI;
		while (a > 2 * PI) a -= 2 * PI;
		if (crossSign(v[i], v[(i+1)%v.S], v[(i+2)%v.S]) <= 0) return false;
		sign += a > PI ? 1 : -1;
		sum += a;
	}
	
	if (abs(sign) != v.S) return false;
	if (abs(abs(sum) - (v.S - 2) * PI) > 1e-6 && abs(abs(sum) - (v.S + 2) * PI) > 1e-6) return false;
	return true;
}

bool isConvexSimple(VI &v) {
	REP(i, v.S) {
		if (crossSign(v[i], v[(i+1)%v.S], v[(i+2)%v.S]) <= 0) return false;
	}
	
	return true;
}

inline bool lineIntersect(int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4) {
	int d = (x2 - x1) * (y4 - y3) - (y2 - y1) * (x4 - x3);
	if (d == 0) return false;
	int p1 = (y1 - y3) * (x4 - x3) - (x1 - x3) * (y4 - y3);
	int p2 = (y1 - y3) * (x2 - x1) - (x1 - x3) * (y2 - y1);
	if (d < 0) d = -d, p1 = -p1, p2 = -p2;
	return (p1 >= 0 && p1 <= d && p2 >= 0 && p2 <= d);
}

bool checkSides(VI &v) {
	int mn = 1 << 30, mx = 0;
	mn = min(mn, D[v[0]][v[v.S - 1]]);
	mx = max(mx, D[v[0]][v[v.S - 1]]);
	REP(i, v.S - 1) {
		mn = min(mn, D[v[i]][v[i + 1]]);
		mx = max(mx, D[v[i]][v[i + 1]]);
	}
	return (mn * 10000LL >= mx * SV);
}

bool checkRadii(VI &v) {
	PII sum;
	REP(i, v.S) sum.X += P[v[i]].X, sum.Y += P[v[i]].Y;
	
	int mn = 1 << 30, mx = 0;
	REP(i, v.S) {
		int d = dist(sum.X, sum.Y, P[v[i]].X * v.S, P[v[i]].Y * v.S);
		mn = min(mn, d);
		mx = max(mx, d);
	}
	return (mn * 10000LL >= mx * RV);
}

double propSides(VI &v) {
	int mn = 1 << 30, mx = 0;
	mn = min(mn, D[v[0]][v[v.S - 1]]);
	mx = max(mx, D[v[0]][v[v.S - 1]]);
	REP(i, v.S - 1) {
		mn = min(mn, D[v[i]][v[i + 1]]);
		mx = max(mx, D[v[i]][v[i + 1]]);
	}
	return (double)(mn * 10000LL) / (mx * SV);
}

double propRadii(VI &v) {
	PII sum;
	REP(i, v.S) sum.X += P[v[i]].X, sum.Y += P[v[i]].Y;
	
	int mn = 1 << 30, mx = 0;
	REP(i, v.S) {
		int d = dist(sum.X, sum.Y, P[v[i]].X * v.S, P[v[i]].Y * v.S);
		mn = min(mn, d);
		mx = max(mx, d);
	}
	return (double)(mn * 10000LL) / (mx * RV);
}

int countRVPts(VI &rv) {
	int pts = 0;
	REP(i, rv.S) if (rv[i] != -1) pts += polys[rv[i]].S;
	return pts;
}

void showRVSummary(VI &rv) {
	double t = getTime();
	set<int> polysFound;
	VI sizeNo(MAXSIZE);
	REP(i, rv.S) {
		if (rv[i] == -1 || polysFound.count(rv[i])) continue;
		polysFound.insert(rv[i]);
		sizeNo[polys[rv[i]].S]++;
	}
	
	while (sizeNo.S && sizeNo[sizeNo.S - 1] == 0) sizeNo.pop_back();
	int score = 0;
	int points = 0;
	REP(i, sizeNo.S) score += sizeNo[i] * i * i, points += sizeNo[i] * i;
	if (sizeNo.S) sizeNo.erase(sizeNo.begin(), sizeNo.begin() + 3);	
	cerr << "Summary: " << score << " (" << (N - points) << ") -> "; print(sizeNo);	
	debugTimer += getTime() - t;
}

void showPolySummary() {
	double t = getTime();

	VI sizeNo(1000);
	REP(i, polys.S) sizeNo[polys[i].S]++;
	
	while (sizeNo.S && sizeNo[sizeNo.S - 1] == 0) sizeNo.pop_back();	
	print(sizeNo);	
	
	debugTimer += getTime() - t;	
}

int BPLType = 0;
int BPLCrossSign = 1;
int BPLSteps = 0;
int BPLMinSize = 5;
int BPLMaxSize = 10;
int BPLFound;
int BPLFoundLimit = 0;
int BPLUseTimer = 0;
int BPLMaxDistOpp = 0;
double BPLAngSum = 0;

double BPLAngle[5000];

VD BPLAng;
hash_set<int> BPLCanVisit;

int BPLStepsSize3 = 0;
int BPLStepsSize4 = 0;
int BPLStepsSize5 = 0;
int BPLStepsFinal = 0;


int buildPolyLinear(VI &v, VI &cur, int pos, int mn, int mx, int force = 0) {
	if (pos == v.S) return 0;
	
	if ((BPLSteps++ & 127) == 0 && BPLUseTimer && getTime() - startTime > PART0_LIMIT) throw 1; 
	//if (cur.S >= 3) BPLStepsSize3++;
	//if (cur.S >= 4) BPLStepsSize4++;
	//if (cur.S >= 5) BPLStepsSize5++;
	
	if (force) {
		BPLCanVisit.clear();
		//BPLAngSum = 0;
		//BPLAng.clear();
		cur.PB(v[pos]);
		return buildPolyLinear(v, cur, pos + 1, mn, mx);
	}
	
	
	int cv = 0;
	//double BPLAngSumCopy = BPLAngSum;
	
	int d = D[cur[cur.S - 1]][v[pos]];
	int newmn = min(mn, d);
	int newmx = max(mx, d);
	if (newmn * 10000LL < newmx * SV) goto out;
	if (BPLCrossSign && cur.S >= 2 && (crossSign(cur[cur.S - 2], cur[cur.S - 1], v[pos]) <= 0 || crossSign(cur[cur.S - 1], v[pos], v[0]) <= 0)) goto out;
	
	/*
	if (BPLMaxDistOpp && cur.S >= 3) {
		REP(i, cur.S) if (D[v[pos]][cur[i]] > BPLMaxDistOpp) goto out;
	}*/
	
	/*
	if (false && cur.S >= 2) {
		BPLAngSum += angle(cur[cur.S - 2], cur[cur.S - 1], v[pos]);
		if (BPLAngSum > (double)(cur.S - 1) / BPLMinSize * 1.1) goto out;
		if (BPLAngSum < (double)(cur.S - 1) / 10) goto out;
	}*/
	
	/*
	if (BPLAngle.S && cur.S >= 3) {
		double pr1 = (double)cur.S / BPLMaxSize;
		if (BPLAngle[pos] < pr1 * 2 / 3) return 0;
		if (cur.S < BPLMinSize) {
			double pr2 = (double)cur.S / (BPLMinSize);
			if (BPLAngle[pos] > 1.0 - (1.0 - pr2) / 2) return 0;
		}
	}*/
	
	
	
	cur.PB(v[pos]);
	if (cur.S < BPLMaxSize) {
		int chash = newmn * 137 + newmx * 277 + pos;
		if (!BPLCanVisit.count(chash)) {
			if (!buildPolyLinear(v, cur, pos + 1, newmn, newmx)) {
				BPLCanVisit.insert(chash);
			}
		}
	}
	
	if (cur.S >= BPLMinSize) {
		d = D[v[pos]][cur[0]];
		newmn = min(newmn, d);
		newmx = max(newmx, d);
		if (newmn * 10000LL >= newmx * SV) {
			cv = 1;
			BPLStepsFinal++;
			if (checkRadii(cur) && (BPLCrossSign == 0 || crossSign(cur[cur.S - 2], cur[cur.S - 1], cur[0]) > 0) && isConvex(cur)) {
				addPoly(cur);
				if (BPLFoundLimit && ++BPLFound >= BPLFoundLimit) throw 1;
			}
		}
	}
	
	cur.pop_back();
	
	out: ;
	
	//BPLAngSum = BPLAngSumCopy;
	//while (BPLAng.S >= cur.S) BPLAng.pop_back();
	cv |= buildPolyLinear(v, cur, pos + 1, mn, mx);
	
	return cv;
}

double FAPtimer = 0;
double FPCtimer0 = 0;
double FPCtimer1 = 0;
double FPCtimer2 = 0;
int FAPMin = 5;
int FAPToFind = 1;

int FAPSteps = 0;
int FAPChecks = 0;
int FAPConvex = 0;

int FAPAt[100];
bool findAnyPoly(VI &v, int tries, int size) {
	if (v.S <= size) return false;
	
	int best[MAXSIZE];
	
	int last = 1;
	best[0] = 0;
    FOR(i, 1, size) {
        double target = (double)i / size;
        int moved = 0;
        while (BPLAngle[last] < target && last < v.S - (size - i)) last++, moved = 1;
        if (moved == 0 || BPLAngle[last] - target < target - BPLAngle[last - 1])
            best[i] = last++;
        else
            best[i] = last - 1;
    }    
	
	v.PB(v[0]);	
	best[size] = v.S - 1;
	
	VI dists(size);
	REP(i, size) dists[i] = D[v[best[i]]][v[best[i+1]]];
	
	int failsLimit = 5;
	int found = 0;
	
	int polyFound = 0;
	REP(abc, tries) {
		FAPSteps++;
		int mn = dists[0], mx = dists[0];
		int pmn = 0, pmx = 0;
		FOR(i, 1, size) {
			if (dists[i] < mn) {
				mn = dists[i];
				pmn = i;
			} else if (dists[i] > mx) {
				mx = dists[i];
				pmx = i;
			}
		}
		
		if (!found && abc * 6 > tries) return false;
		
		int prob = 0;
		if (mn * 10000LL >= mx * SV) {
			FAPChecks++;
			prob = 1;
			found = 1;
			
			int sumx = 0;
			int sumy = 0;
			REP(i, size) {
				sumx += P[v[best[i]]].X;
				sumy += P[v[best[i]]].Y;
			}
			int mn2 = 1 << 30, mx2 = 0;
			
			REP(i, size) {
				int d = dist(sumx, sumy, P[v[best[i]]].X * size, P[v[best[i]]].Y * size);
				if (d < mn2) {
					mn2 = d;
					pmn = i;
				} 
				if (d > mx2) {
					mx2 = d;
					pmx = i;
				}
			}
			
			if (mn2 * 10000LL >= mx2 * RV) {
				FAPConvex++;
				bool crossFailed = false;
				REP(i, size) {
					if (crossSign(v[best[i]], v[best[(i+1)%size]], v[best[(i+2)%size]]) <= 0) {
						crossFailed = true;
						pmn = pmn = (i+1)%size;
						break;
					}
				}
				if (!crossFailed) {
					VI poly(size);
					REP(i, size) poly[i] = v[best[i]];
					addPoly(poly);
					//FAPAt[abc / 10]++;
					if (++polyFound > FAPToFind) return true;
					//return true;
				}
				prob = 2;
			}
		}
		
		int type = rand() % 2;
		int bb;
		int bv;
		int vv;
		int d;
		
		if (prob == 0) {
			bb = type == 0 ? pmn : pmx;
			bv = type == 0 ? mn : mx;
		} else if (prob == 1) {
			if (rand() & 1) {
				bb = rand() % size;
			} else {
				bb = pmn == 0 ? pmx : pmx == 0 ? pmn : (rand() & 1) ? pmn : pmx;
			}
			bv = dists[bb];
		} else {
			bb = pmn;			
			bv = dists[bb];
		}
		
		vv = bb == 0 ? 1 : bb == size - 1 ? size - 1 : bb + (rand() & 1);
		d = bb == vv ? -1 : 1;
		if (type == 1) d = -d;
		int v2 = v[best[bb == vv ? vv + 1 : vv - 1]];
		
		int m0 = best[vv - 1];
		int m1 = best[vv + 1];
		
		while (1) {
			best[vv] += d;
			double dd = D[v2][v[best[vv]]];
			
			if (best[vv] == m0 || best[vv] == m1) {
				if (--failsLimit == 0) return false;
				best[vv] -= d;
				bv = -1;
			}
			
			if ((type == 0 && dd > bv || type == 1 && dd < bv) || bv == -1) {
				dists[vv - 1] = D[v[best[vv]]][v[best[vv - 1]]];
				dists[vv + 0] = D[v[best[vv]]][v[best[vv + 1]]];
				break;
			}
		}
		
		
		
		
	}
	
	return false;
	
}

int FPCSteps = 0;
int FPCType = 0;
double FPCPA[MAXN];

bool FPCVASort(int a, int b) {
	return FPCPA[a] < FPCPA[b];
}

void findAllPolyCentered(int x, int y) {
	pair < int, int > FPCPD[MAXN];
	int FPCVA[MAXN];
	int FPCVANEW[MAXN];

	double t0 = getTime();

	int no = 0;
	REP(i, N) if (POk[i]) {
		FPCPD[no++] = MP(dist(x, y, P[i].X, P[i].Y), i);
		FPCPA[i] = fastatan2((double)(P[i].Y - y), (double)(P[i].X - x)) / (2 * PI);
	}
	
	sort(FPCPD, FPCPD + NOk);
	
	int wallDist = 1 << 20;
	wallDist = min(wallDist, dist(x, y, 0, y));
	wallDist = min(wallDist, dist(x, y, 700, y));
	wallDist = min(wallDist, dist(x, y, x, 0));
	wallDist = min(wallDist, dist(x, y, x, 700));
	
	BPLUseTimer = 1;
	BPLMinSize = 4;
	BPLMaxSize = 8;
	BPLType = 0;
	BPLCrossSign = 1;
	//BPLFoundLimit = NOk > 1000 ? 1 : 10;
	BPLFoundLimit = 5;
	
	FPCtimer0 += getTime() - t0;
	
	int lastp = 0;
	int lastx = -1;
	int lastz = -1;
	int vano = 0;
	REP(i, N) {
		double t1 = getTime();
		
		FPCSteps++;
		if (FPCPD[i].X > wallDist) break;
		
		
		int maxDistance = (int)(FPCPD[i].X * 10000LL / RV) * 140 / 100 - FPCPD[i].X * 40 / 100;
		
		int vanewno = 0;
		for (; lastp < N; lastp++) {
			if (FPCPD[lastp].X > maxDistance) break;
			FPCVANEW[vanewno++] = FPCPD[lastp].Y;
		}
		
		//print(va);
		//print(vanew);
		if (lastx != -1) {
			if (lastz == -1) {
				REP(j, vano) if (FPCVA[j] == lastx) {
					lastz = j;
					break;
				}
			}
			vano--;
			FOR(j, lastz, vano) FPCVA[j] = FPCVA[j + 1];
		}
		lastx = FPCPD[i].Y;
		
		if (vanewno) {
			int ap = (int)vano - 1;
			sort(FPCVANEW, FPCVANEW + vanewno, FPCVASort);
			for (int j = (int)vanewno - 1; j >= 0; j--) {
				while (ap >= 0 && FPCPA[FPCVA[ap]] > FPCPA[FPCVANEW[j]]) {
					FPCVA[ap + j + 1] = FPCVA[ap];
					ap--;
				}
				FPCVA[ap + j + 1] = FPCVANEW[j];
			}
			vano += vanewno;
		}		
		//print(va);
		
		int bad = 0;
		
		if (FPCType == 1) {
			int lo = FAPMin, hi = min(MAXSIZE - 1, biggestPoly + 1);
			//hi = min(hi, vano / 2);
			if (lo > hi) {
				bad = 1;
			} else {
				BPLMinSize = lo + rand() % (hi - lo + 1);
				if (vano <= BPLMinSize * 2) bad = 1;
			}
		} else 
			if (vano <= BPLMinSize) bad = 1;
			
		FPCtimer1 += getTime() - t1;
		
		lastz = -1;
		if (bad) continue;
		
		double t2 = getTime();
		
		VI v(vano);
		REP(j, vano) if (FPCVA[j] == FPCPD[i].Y) {
			lastz = j;
			double a0 = FPCPA[FPCVA[j]];
			REP(k, vano) {
				int z = k + j; 
				if (z >= vano) z -= vano;
				v[k] = FPCVA[z];
				double aa = FPCPA[FPCVA[z]] - a0;
				if (aa < 0) aa += 1;
				BPLAngle[k] = aa;
			}
			break;
		}
		
		//print(VD(BPLAngle, BPLAngle + va.S));
		
		FPCtimer2 += getTime() - t2;
		
		double t = getTime();
		try {
			if (FPCType == 0) {
				BPLFound = 0;
				VI cur;
				int dist = (int)(4 * 3.14 * 3.14 * FPCPD[i].X / BPLMaxSize / BPLMaxSize);
				int minDist = 0;
				int maxDist = dist * 2 * 9;
				
				REP(j, vano - 1) minDist = max(minDist, D[v[j]][v[j+1]]);
				minDist = max(minDist, D[v[v.S - 1]][v[0]]);
				buildPolyLinear(v, cur, 0, maxDist, minDist, 1);
			} else if (FPCType == 1) {
				findAnyPoly(v, BPLMinSize * 50, BPLMinSize);
			}
		} catch (int e) { }
		if (polys.S) biggestPoly = max(biggestPoly, (int)polys[polys.S - 1].S);
		if (FPCType == 1) FAPtimer += getTime() - t;
		if (FPCType == 0 && getTime() - startTime > PART0_LIMIT) break;
		if (FPCType == 1 && getTime() - startTime > PART1_LIMIT) break;
	}
}

//void searchPoly

VI findAllTriangles(VI rv = VI(), int add = 0) {
	if (rv.S == 0) {
		rv = VI(N, -1);
		REP(i, N) if (!POk[i]) rv[i] = -2; 
	}
	
	REP(i, N) if (rv[i] == -1) {
		VC < pair < int, int > > v;
		REP(j, i) if (rv[j] == -1) v.PB(MP(D[i][j], j));
		sort(ALL(v));
		REP(j, v.S) {
			if (getTime() - startTime > (TIME_LIMIT * 0.99)) return rv;
			int mx = (v[j].X * 10000LL / SV + 1);
			FOR(k, j + 1, v.S) {
				if (v[k].X > mx) break;
				int d = D[v[j].Y][v[k].Y];
				if (d <= mx && d * 10000LL >= v[k].X * SV) {
					VI tri; tri.PB(i); tri.PB(v[j].Y); tri.PB(v[k].Y);
					if (checkRadii(tri) && checkSides(tri)) {
						addPoly(tri);
						if (add) {
							rv[i] = rv[v[j].Y] = rv[v[k].Y] = polys.S - 1;
							goto next;
						}
					}
				}
			}
		}
		next: ;
	}
	return rv;
}

VI findAllQuads(VI rv = VI(), int add = 0) {
	if (rv.S == 0) {
		rv = VI(N, -1);
		REP(i, N) if (!POk[i]) rv[i] = -2; 
	}
	
	int tries = 0;
	int tries2 = 0;
	int tries3 = 0;	
	
	REP(i, N) if (rv[i] == -1) {
		VC < pair < int, int > > v;
		REP(j, i) if (rv[j] == -1) v.PB(MP(D[i][j], j));
		sort(ALL(v));
		REP(j, v.S) {
			if (getTime() - startTime > (TIME_LIMIT * 0.95)) goto out;
			int mx = (v[j].X * 10000LL / SV + 1);
			FOR(k, j + 1, v.S) {
				if (v[k].X > mx) break;
				if (D[v[j].Y][v[k].Y] < v[k].X * 11 / 10) continue;
				if (D[v[j].Y][v[k].Y] > v[j].X * 3) continue;
				REP(l, i) if (rv[l] == -1 && l != v[k].Y && l != v[j].Y) {
					if (D[i][l] < v[j].X) continue;
					if (D[i][l] > v[j].X * 3) continue;
					int d1 = D[l][v[k].Y];
					int d2 = D[l][v[j].Y];
					int newmx = max(v[k].X, max(d1, d2));
					int newmn = min(v[j].X, min(d1, d2));
					if (newmn * 10000LL < newmx * SV) continue;
					VI p; p.PB(i); p.PB(v[j].Y); p.PB(l); p.PB(v[k].Y);
					tries2++;
					if (!checkRadii(p)) continue;
					tries3++;
					if (lineIntersect(P[p[0]].X, P[p[0]].Y, P[p[1]].X, P[p[1]].Y, P[p[2]].X, P[p[2]].Y, P[p[3]].X, P[p[3]].Y)) continue;
					if (lineIntersect(P[p[1]].X, P[p[1]].Y, P[p[2]].X, P[p[2]].Y, P[p[3]].X, P[p[3]].Y, P[p[0]].X, P[p[0]].Y)) continue;
					addPoly(p);
					if (add) {
						rv[i] = rv[v[j].Y] = rv[v[k].Y] = rv[l] = polys.S - 1;
						goto next;
					}
				}
			}
		}
		next: ;
	}
	out: ;
	//DB(tries);
	//DB(tries2);
	//DB(tries3);
	return rv;
}

bool findAllPolys() {
	BPLType = 0;
	BPLFoundLimit = 0;
	BPLCrossSign = 1;
	BPLMinSize = 5;
	BPLMaxSize = N;
	BPLUseTimer = 0;
	
	VC < pair < int, int > > vx;
	REP(i, N) if (POk[i]) vx.PB(MP((-P[i].X * 1000 + P[i].Y), i));
	sort(ALL(vx));
	VI order;
	REP(i, vx.S) order.PB(vx[i].Y);
	
	VD times;
	
	double start = getTime();
		
	REP(i, order.S) {
		VC < pair < double, int > > va;
		REP(j, i) va.PB(MP(fastatan2((double)P[order[j]].Y - P[order[i]].Y, (double)P[order[j]].X - P[order[i]].X), order[j]));
		va.PB(MP(-1e9, order[i]));
		sort(ALL(va));
		
		VI v(va.S);
		REP(j, va.S) v[j] = va[j].Y;
		
		VI cur;
		buildPolyLinear(v, cur, 0, 1 << 30, 0, 1);
		
		
		double curTime = getTime() - start;
		times.PB(curTime);
		
		if (curTime > FINDALL_ESTIMATE) {
			if (i < 10) return false;
			
			double diff = pow(curTime / times[i - 10], 0.1);
			double est = curTime * pow(diff, (double)(order.S - i));
			
			double curTime15 = min(1.0, curTime * sqrt(curTime));
			double expEst = FINDALL_TIMELIMIT / curTime15;
			
			if (est > expEst) {
				cerr << "FindAllPolys() Failed at " << i << " time: " << curTime << " est: " << est << " maxEst: " << expEst << endl;				
				return false;
			}
			
		}
				
	}

	DB(BPLSteps);
	DB(BPLStepsSize3);
	DB(BPLStepsSize4);
	DB(BPLStepsSize5);
	DB(BPLStepsFinal);
	
	cerr << "Found All Polys" << endl;
	printTime();
	
	return true;
}

VI getSizeNo() {
	VI sizeNo(MAXSIZE);
	REP(i, polys.S) sizeNo[polys[i].S]++;
	return sizeNo;
}

VI createPolyOrder() {
	VI sizeNo(MAXSIZE);
	REP(i, polys.S) sizeNo[polys[i].S]++;
	for (int i = MAXSIZE - 1; i >= 1; i--) sizeNo[i - 1] += sizeNo[i];
	
	VI order(polys.S);
	REP(i, polys.S)	order[--sizeNo[polys[i].S]] = i;
	return order;
}

VI createPolyOrder(VI sizeNo) {
	for (int i = MAXSIZE - 1; i >= 1; i--) sizeNo[i - 1] += sizeNo[i];
	
	VI order(polys.S);
	REP(i, polys.S)	order[--sizeNo[polys[i].S]] = i;
	return order;
}


int FSFMax;
int FSFSteps;
VI FSFBest;
void findSolFullRec(VI &v, VI &cur, int pos, int pts, int left) {
	if (pts > FSFMax) {
		FSFMax = pts;
		FSFBest = cur;
		cerr << FSFSteps << ' ' << FSFMax << endl;
	}
	
	if (pos == v.S) return;
	
	if (pos && polys[v[pos]].S != polys[v[pos-1]].S) {
		if (FSFMax >= pts + (left / polys[v[pos]].S) * polys[v[pos]].S * polys[v[pos]].S) return;
	}
	
	if (++FSFSteps % 100000000 == 0) cerr << FSFSteps << ' ' << FSFMax << endl;
	
	int ok = 1;
	REP(i, polys[v[pos]].S)
		if (cur[polys[v[pos]][i]] != -1) {
			ok = 0;
			break;
		}
		
	if (ok) {
		REP(i, polys[v[pos]].S) cur[polys[v[pos]][i]] = v[pos];
		findSolFullRec(v, cur, pos + 1, pts + polys[v[pos]].S * polys[v[pos]].S, left - polys[v[pos]].S);
		REP(i, polys[v[pos]].S) cur[polys[v[pos]][i]] = -1;
	}
	
	findSolFullRec(v, cur, pos + 1, pts, left);
}

VI findSolFull() {
	VI order = createPolyOrder();
	
	VI cur(N, -1);

	FSFMax = 0;
	FSFSteps = 0;
	FSFBest = VI(N, -1);
	findSolFullRec(order, cur, 0, 0, N);
	
	return FSFBest;
}

int first = 0;
VI findMIS(VVI &p) {
	int n = p.S;
	if (n == 0) return VI();	
	int m = p[0].S;
	
	VVI conn(n, VI());
	VVI cm(n, VI(n));
	REP(i, n) REP(j, i) {
		hash_set<int> v0;
		REP(k, m) v0.insert(p[i][k]);
		int ok = 0;
		REP(k, m) if (v0.count(p[j][k])) {
			ok = 1;
			break;
		}
		if (ok) {
			conn[i].PB(j);
			conn[j].PB(i);
			cm[i][j] = 1;
			cm[j][i] = 1;
		}
	}
	
	VI d(n);
	REP(i, n) d[i] = conn[i].S;
	
	VI s(n);
	REP(i, n) REP(j, conn[i].S) s[i] += d[conn[i][j]];
	
	VI used(n);
	
	if (first) {
		REP(i, n) {
			REP(j, n) cerr << cm[i][j];
			cerr << endl;
		}		
	}
	while (true) {
		int bv = -1;
		int bb = -1;
		REP(i, n) if (!used[i] && (s[i] > bv || s[i] == bv && d[i] > d[bb])) {
			bv = s[i];
			bb = i;
		}
				
		if (bv <= 0) break;
		if (first) DB(bb);
		
		used[bb] = 1;
		
		/*
		d = VI(n);
		REP(i, n) if (!used[i]) REP(j, conn[i].S) d[i] += !used[conn[i][j]];
		
		s = VI(n);
		REP(i, n) if (!used[i]) REP(j, conn[i].S) s[i] += d[conn[i][j]];
		*/
		
		REP(i, conn[bb].S) {
			int v = conn[bb][i];
			if (used[v]) continue;
			d[v]--;
			s[v] -= d[bb];
			REP(j, conn[v].S) s[conn[v][j]]--;
		}
		s[bb] = 0;
		d[bb] = 0;
	}
	
	first = 0;
	REP(i, n) used[i] = !used[i];
	return used;
}

VI findSolGreedyEx(VI cur = VI()) {
	if (cur.S == 0) cur = VI(N, -1);
	
	VI sizeNo = getSizeNo();
	VI order = createPolyOrder(sizeNo);
	
	int curPos = 0;
	
	for (int s = MAXSIZE; s >= 3; s--) {		
		if (curPos == polys.S || polys[order[curPos]].S != s) continue;
		
		if (false && sizeNo[s] <= 250) {
			VVI np;
			VI pn;
			while (curPos < polys.S && polys[order[curPos]].S == s) {
				int p = order[curPos];
				int ok = 1;
				REP(i, polys[p].S) {
					if (cur[polys[p][i]] != -1) {
						ok = 0;
						break;
					}
				}
				if (ok) {
					pn.PB(p);
					np.PB(polys[p]);
				}
				curPos++;
			}
			VI use = findMIS(np);
			int no = 0;
			REP(i, use.S) if (use[i]) {
				REP(j, s) cur[np[i][j]] = pn[i];
				no++;
			}
			int errors = 0;
			REP(i, pn.S) {
				int ok = 1;
				REP(j, s) {
					if (cur[np[i][j]] != -1) {
						ok = 0;
						break;
					}
				}
				if (ok) {
					errors++;
					REP(j, s) cur[np[i][j]] = pn[i];
				}
			}
			//cerr << s << ' ' << sizeNo[s] << ' ' << no << ' ' << errors << endl;
			continue;
		}
		
		VC < pair < int, int > > vScore;
		REP(i, N) vScore.PB(MP(0, i));
		
		VVI conn(N, VI());
		
		while (curPos < polys.S && polys[order[curPos]].S == s) {
			int p = order[curPos];
			int ok = 1;
			REP(i, polys[p].S) {
				if (cur[polys[p][i]] != -1) {
					ok = 0;
					break;
				}
			}
			if (ok) {
				REP(i, polys[p].S) {
					conn[polys[p][i]].PB(p);
					vScore[polys[p][i]].X++;
				}
			}
			curPos++;
		}
		
		sort(ALL(vScore));

		REP(i, N) {
			int v = vScore[i].Y;
			if (cur[v] != -1) continue;
			REP(j, conn[v].S) {
				int x = conn[v][j];
				REP(k, polys[x].S) {
					if (cur[polys[x][k]] != -1) goto out;
				}
				REP(k, polys[x].S) cur[polys[x][k]] = x;
				out: ;
			}
		}
	}
	
	return cur;
}

VI findSolIterGreedy(VI cur = VI(), double time = 1.0, double tt = 1.0) {
    double start = getTime();
    
    if (cur.S == 0) cur = VI(N, -1);
    
    int scoreStart = 0, scoreEnd = 0;
    REP(i, cur.S) if (cur[i] != -1) scoreStart += polys[cur[i]].S;
    
    VI order = createPolyOrder();
    
    VVI pp(N, VI());
    REP(i, polys.S) {
        int p = order[i];
        REP(j, polys[p].S) pp[polys[p][j]].PB(p);
    }
    
    DB(scoreStart);
    
    int steps = 0;
    int stepsAcc = 0;
    
    VI bestCur = cur;
    VI bestRV;
	
    int scoreBest = scoreStart;
    int scoreAct = scoreStart;
	int left = N;
	int newfree[1000];
	int values[1000];
	int valuesp[1000];
	int used[5000];
	ZERO(used);
	double expV[1000];
	int expComp[1000];
	ZERO(expComp);
	int expCompNo = 0;
	double temp;
    while (getTime() - start < time) {
        int v = rand() % N;
        
        if ((rand() % 20 != 0 && cur[v] != -1) || pp[v].S == 0) continue;
        int pno = pp[v][rand() % pp[v].S];
        steps++;
        
        int diff = polys[pno].S * polys[pno].S;
		
		if (steps % 100000 == 1) {
			temp = (1.0 - (getTime() - start) / time) * tt + 0.5;
			expCompNo++;
        }
		
		int leftDiff = -polys[pno].S;
		int newfreeno = 0;
		int valuesno = 0;
        REP(i, polys[pno].S) {
            int x = polys[pno][i];
            if (cur[x] != -1) {
                int p2 = cur[x];
                diff -= polys[p2].S * polys[p2].S;
                REP(j, polys[p2].S) {
					if (used[polys[p2][j]] != steps) {
						used[polys[p2][j]] = steps;
						values[valuesno] = cur[polys[p2][j]];
						valuesp[valuesno++] = polys[p2][j];
					}
                    cur[polys[p2][j]] = -1;
                    newfree[newfreeno++] = polys[p2][j];
                }
            }
			if (used[x] != steps) {
				used[x] = steps;
				values[valuesno] = cur[x];
				valuesp[valuesno++] = x;
			}
            cur[x] = pno;
        }
        	
		/*
        REP(i, newfreeno) {
            int u = newfree[i];
            if (cur[u] != -1) continue;
			REP(abc, 5) {
				int p = pp[u][rand() % pp[u].S];
                int ok = 1;
                REP(k, polys[p].S) ok &= cur[polys[p][k]] == -1;
                if (ok) {
                    diff += polys[p].S * polys[p].S;
                    REP(k, polys[p].S) {
						if (used[polys[p][k]] != steps) {
							used[polys[p][k]] = steps;
							values[valuesno] = cur[polys[p][k]];
							valuesp[valuesno++] = polys[p][k];
						}
						cur[polys[p][k]] = p;
					}
                    break;
                }
            }
        }*/
		
		if (diff < 0 && expComp[-diff] != expCompNo) {
			expV[-diff] = exp(diff / temp);
			expComp[-diff] = expCompNo;
		}
        
        if (diff >= 0 || randDouble() < expV[-diff]) {
			stepsAcc++;
			scoreAct += diff;
			if (scoreAct > scoreBest) {
				scoreBest = scoreAct;
				bestRV = cur;
				cerr << "At: " << steps << " New: " << scoreBest << endl;
			}
        } else {
			REP(i, valuesno) {
				cur[valuesp[i]] = values[i];
			}
		}
        
        
    }
    
    DB(scoreAct);
    DB(steps);
	DB(stepsAcc);
    
    return bestRV;
}

VS createRV(VI &rv) {
	VS vs;
	
	set<int> polysFound;
	REP(i, rv.S) {
		if (rv[i] == -1 || polysFound.count(rv[i])) continue;
		polysFound.insert(rv[i]);
		string s = "";
		REP(j, polys[rv[i]].S) {
			if (j) s += " ";
			s += i2s(polys[rv[i]][j]);
		}
		vs.PB(s);
	}
	
	return vs;
}

class QualityPolygons {public: VS choosePolygons(VI &points, int sidesDiff, int radiiDiff) {
	startTime = getTime();
	
	//parsing
	N = points.S / 2;
	NOk = N;
	REP(i, N) P[i].X = points[i + i], P[i].Y = points[i + i + 1];
	SD = sidesDiff;
	RD = radiiDiff;	
	
	REP(i, N) POk[i] = true;
	
	SV = (100 - SD) * (100 - SD);
	RV = (100 - RD) * (100 - RD);
	
	REP(i, N) REP(j, i) D[i][j] = D[j][i] = dist(i, j);
	
	initfastatan();
	
	printTime();
	
	if (N <= 500 && findAllPolys()) {
		showPolySummary();
		printTime();
		findAllTriangles();
		printTime();
		findAllQuads();
		printTime();
		showPolySummary();
		VI rv = findSolIterGreedy(VI(), 10.0, 3.0);
		//VI rv = findSolReorderEx(VI(), TIME_LIMIT - (getTime() - startTime));
		showRVSummary(rv);
		return createRV(rv);
	}
	
	FPCType = 0;
	int centersTested = 0;
	
	int NOkSteps = 10;
	int NOkCurStep = 0;
	int maxNOk = 2000;
	int minNOk = 1000;
	double NOkMinProp = 0.5;
	
	int oldPolys;

	VC < pair < int, int > > spPts;
/*
#define SSP_MARGIN 100
#define SSP_STEP 20
#define SSP_NO ((700 - SSP_MARGIN * 2) / SSP_STEP + 1)
#define SSP_WANT 10
#define SSP_MD 50

	int sspno[SSP_NO][SSP_NO];
	
	REP(x, SSP_NO) REP(y, SSP_NO) {
		int px = x * SSP_STEP + SSP_MARGIN;
		int py = y * SSP_STEP + SSP_MARGIN;
		VI v(N);
		REP(i, N) v[i] = dist(px, py, P[i].X, P[i].Y);
		sort(ALL(v));
		int mx = 0;
		int last = 0;
		REP(i, N) {
			int maxDist = (int)(v[i] * 10000LL / RV) * 14 / 10 - v[i] * 4 / 10;
			int wallDist = 1 << 20;
			wallDist = min(wallDist, dist(px, py, 0, py));
			wallDist = min(wallDist, dist(px, py, 700, py));
			wallDist = min(wallDist, dist(px, py, px, 0));
			wallDist = min(wallDist, dist(px, py, px, 700));
			if (v[i] > wallDist) break;
			
			while (last < N && v[last] < maxDist) last++; 
			mx = max(mx, last - i);
		}		
		sspno[x][y] = mx;
	}
	REP(i, SSP_WANT) {
		int bv = 0, bx = -1, by = -1;
		REP(x, SSP_NO) REP(y, SSP_NO) {
			if (sspno[x][y] > bv) {
				bx = x, by = y, bv = sspno[x][y];
			}
		}
		cerr << bx << ' ' << by << ' ' << sspno[bx][by] << endl;
		if (bv == 0) break;
		REP(x, SSP_NO) REP(y, SSP_NO) {
			if (dist(bx, by, x, y) * SSP_STEP * SSP_STEP < SSP_MD * SSP_MD) 
				sspno[x][y] = 0;
		}
		spPts.PB(MP(bx * SSP_STEP + SSP_MARGIN, by * SSP_STEP + SSP_MARGIN));		
	}
	reverse(ALL(spPts));
*/
	
	srand(1);
	unsigned int seedA = rand();
	unsigned int seedB = rand() * rand();
	VI curRV = VI(N, -1);
	while (getTime() - startTime < PART1_LIMIT) {
		if (FPCType == 0 && getTime() - startTime > PART0_LIMIT) {	
			FPCType = 1;
			showPolySummary();
			FPCSteps = 0;
			DB(centersTested);
			oldPolys = polys.S;
			seedA = seedB;
		} else if (FPCType == 1 && spPts.S == 0 && polys.S - oldPolys > 500 && N > 1000 && ((getTime() - startTime) - PART0_LIMIT) / (PART1_LIMIT - PART0_LIMIT) * NOkSteps > NOkCurStep + 1) {
			biggestPoly = 4;
			NOkCurStep++;
			double t = getTime();
			VI TRV = findSolGreedyEx(curRV);
			debugTimer += getTime() - t;
			showRVSummary(TRV);
			VI order = createPolyOrder();
			int MinNOk = max(minNOk + (maxNOk - minNOk) * (NOkSteps - NOkCurStep) / NOkSteps, (int)(N * (1.0 - NOkMinProp * NOkCurStep / NOkSteps)));
			REP(i, order.S) {
				int p = order[i];
				if (POk[polys[p][0]] && TRV[polys[p][0]] == p) {
					if (MinNOk < NOk) {
						REP(j, polys[p].S) {
							POk[polys[p][j]] = 0;
							NOk--;
						}
					} else {
						REP(j, polys[p].S)
							TRV[polys[p][j]] = -1;
					}
				} 		
			}
			curRV = TRV;
		}
			
		centersTested++;
		int x, y;
		if (FPCType == 1 && spPts.S) {
			x = spPts[spPts.S - 1].X;
			y = spPts[spPts.S - 1].Y;
			spPts.pop_back();
		} else {
			x = 50 + rand_r(&seedA) % 600;
			y = 50 + rand_r(&seedA) % 600;
		}
		
		findAllPolyCentered(x, y);
		//DB(polys.S);
		//printTime();
	}
	
	showPolySummary();
	DB((polys.S - oldPolys));
	DB(FAPSteps);
	DB(FAPChecks);
	DB(FAPConvex);
	DB(FAPtimer);
	DB(FPCSteps);
	DB(FPCtimer0);
	DB(FPCtimer1);
	DB(FPCtimer2);
	DB(centersTested);
	
	//print(VI(FAPAt, FAPAt + 50));
	
	printTime();
	
	
	VI rv = findSolGreedyEx(curRV);
	showRVSummary(rv);
	
    if (getSizeNo()[5] > 1000 && N > 1000) {
        REP(i, rv.S) if (POk[i] && rv[i] != -1) {
            POk[i] = 0;
            NOk--;
        }
        startTime += 0.2;
        FAPMin = 4;
        biggestPoly = 4;
		
        while (getTime() - startTime < PART1_LIMIT) {
            int x = 25 + rand_r(&seedA) % 650;
            int y = 25 + rand_r(&seedA) % 650;
            findAllPolyCentered(x, y);
        }
		rv = findSolGreedyEx(rv);
        startTime -= 0.2;        
        showRVSummary(rv);
    }

	rv = findAllQuads(rv, 1);
	printTime();	
	rv = findAllTriangles(rv, 1);
	printTime();
	
	showRVSummary(rv);
	
	findAllTriangles();
	findAllQuads();
	showPolySummary();
	printTime();
	//rv2 = findSolRandGreedyEx(PART2_LIMIT - (getTime() - startTime));	
	
	if (getTime() - startTime < TIME_LIMIT) {
		VI rv2(N, -1);
		//rv2 = findSolReorderEx(rv2, TIME_LIMIT - (getTime() - startTime));
		rv2 = findSolIterGreedy(rv, TIME_LIMIT - (getTime() - startTime));
		if (countRVPts(rv2) > countRVPts(rv)) rv = rv2;
	}

	showRVSummary(rv);
	showPolySummary();
	
	VS vs = createRV(rv);
	
	
	
	//adding bonus points
	REP(i, BP.S) {
		vs.PB("X " + i2s(BP[i].X) + string(" ") + i2s(BP[i].Y));
	}
	
	DB(debugTimer);
	printTime();
	
	return vs;
}};


#ifdef LOCAL
int main() {
	int n;
	scanf("%d", &n);
	VI pts(n);
	int sd, rd;
	REP(i, n) scanf("%d", &pts[i]);
	scanf("%d %d", &sd, &rd);
	
	QualityPolygons algo;
	VS rv = algo.choosePolygons(pts, sd, rd);
	
	printf("%d\n", (int)rv.S);
	REP(i, rv.S) printf("%s\n", rv[i].c_str());
	fflush(stdout);
}
#endif
