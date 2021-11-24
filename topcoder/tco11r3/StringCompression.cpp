#define DEBUG_LEVEL -1
#define SHOW_CANDIDATES -1
#define SHOW_VALID
#define USE_COMPLEX
#define USE_SIMPLE


//#define LEAVE_HARD
#define MAX_JOIN 2500
#define FULL_CAND 2000
#define TIME_LIMIT 22.0
#define MAX_TIME 29.6
#define EARLY_EXIT

#ifndef LOCAL
#undef EARLY_EXIT
#endif

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <ext/hash_set>
#include <xmmintrin.h>

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
#define PID            pair <int, double>
#define PDD            pair <double, double>
#define VI          VC<int>
#define VVI            VC < VI >
#define VD            VC < double >
#define VPII        VC < PII >
#define VS          VC<string>
#define PIVS		pair < int, VS >
#define VVS         VC< VS >
#define DB(a)        cerr << #a << ": " << a << endl;

void print(VI v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]\n";}
void print(VC < LL > v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]\n";}
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

double randDouble() {
    return ((double)(rand() + 0.5) / ((double)RAND_MAX + 1));
}

double startTime;
int statesVis = 0;

int SIZE_LIMIT = 1;
bool USE_JOINING = false;

string gData;
VI gLM;
VS gCorrect;
VI gCorrectOrder;
VVI gCorrectMatches;
double gRealError = 0;
set<string> gAllCorrect;
int gOptimal = 1 << 20;

double gError = 0;
bool gErrorSet = false;

double timerA = 0;
double timerB = 0;
double timerC = 0;
double timerD = 0;
double timerLev[10];

void showError(string s) {
    cerr << "Error: " << s << endl;
}

#define MAXCH 34
#define LP 7

char CVT[256];
string CV(string &s) {
    string rv = s;
    REP(i, s.S) rv[i] = CVT[s[i]];
    return rv;
}

VS CV(VS &vs) {
    VS rv(vs.S);
    REP(i, vs.S) rv[i] = CV(vs[i]);
    return rv;
}

char UNT[256];
string UN(string &s) {
    string rv = s;
    REP(i, s.S) rv[i] = UNT[s[i]];
    return rv;
}

VS UN(VS &vs) {
    VS rv(vs.S);
    REP(i, vs.S) rv[i] = UN(vs[i]);
    return rv;
}

bool initCalled = false;
void init() {
    if (initCalled) return;
    FOR(i, 'a', 'z' + 1) CVT[i] = i - 'a' + LP;
    FOR(i, '2', '8' + 1) CVT[i] = i - '1';
    REP(i, 256) UNT[CVT[i]] = i;
}

bool hasDigit(string &s) {
	REP(i, s.S) if (isdigit(s[i])) return true;
	return false;
}

bool hasDigit(string &s, char c) {
	REP(i, s.S) if (s[i] == c) return true;
	return false;
}

pair < string, int > chooseBestRep(string &data, int size, int off = 0) {
    string bs = "";
    int bv = -1;
    REP(i, size) {
        VI cnt(128);
        for (int j = i + off; j < data.S; j += size)
            cnt[data[j]]++;
        int bc = 0;
        REP(j, 128) if (cnt[j] > cnt[bc])
            bc = j;
        bv += cnt[bc];
        bs += (char)bc;
    }
    return MP(bs, bv);
}

VS chooseBestRep(string &data, VI &lm) {
    int bv = -1;
    VS rv;
    
    VI usedSize(1000);
    
    VI order;
    REP(i, lm.S) order.PB(i);
    int lmSum = 0;
    REP(i, lm.S) lmSum += lm[i];
    
    do {
        int total = 1;
        int sum = lmSum - lm.S + 1;
        
        REP(i, lm.S) {
            if (total * sum >= data.S && !usedSize[sum]) {
                usedSize[sum] = 1;
                REP(add, (data.S < 5000 ? lm[0] - 1 : 1)) {
                    if (total / lm[0] * (lm[0] - add) * sum < data.S) break;
                    pair < string, int > pp = chooseBestRep(data, sum, add);
                    pp.Y += add;
                    if (pp.Y > bv) {
                        bv = pp.Y;
                        rv = VS(lm.S);
                        REP(j, i) rv[order[j]] = string(lm[order[j]], '1' + order[j+1]);
                        FOR(j, i, lm.S) {
                            rv[order[j]] = (j == lm.S - 1) ? pp.X : pp.X.substr(0, lm[order[j]] - 1) + (char)('1' + order[j+1]);
                            pp.X = pp.X.substr(lm[order[j]] - 1);
                        }
                        rv[0] = data.substr(0, add) + string(lm[0] - add, '1' + order[1]);
                    }
                }
            }
            
            total *= lm[order[i]];
            sum -= lm[order[i]] - 1;
        }
    } while (next_permutation(order.begin() + 1, order.end()));
    
    return rv;
}

string fill(string s, int size) {
    return s.S < size ? s + string(size - s.S, 'a') : s;
}

VVS findSizeMatchings(int size, VI lm) {
    int n = lm.S;

    VI order;
    REP(i, n) order.PB(lm[i]);
    sort(order.begin() + 1, order.end());
    
    VVS rv;
    int ok = 0;
    do {
        int usno = ((n * n - n) / 2);
        REP(j, 1 << (usno + usno)) {
            int sz[8];
            int us[10];
            int t = j;
            REP(i, usno) {
                us[i] = (t & 3) + 1;
                t >>= 2;
            }
            sz[0] = order[n - 1];
            sz[1] = order[n - 2] - us[0] + us[0] * sz[0];
            sz[2] = order[n - 3] - us[1] - us[2] + us[1] * sz[0] + us[2] * sz[1];
            sz[3] = order[n - 4] - us[3] - us[4] - us[5] + us[3] * sz[0] + us[4] * sz[1] + us[5] * sz[2];
            if (n == 5)
                sz[4] = order[n - 5] - us[6] - us[7] - us[8] - us[9] + us[6] * sz[0] + us[7] * sz[1] + us[8] * sz[2] + us[9] * sz[3];
            
            int total = n == 4  ? sz[3] : sz[4];
            if (total == size) {
                VS vs;
                int p = 0;
                REP(i, n) {
                    string s = "";
                    REP(j, i) {
                        REP(k, us[p]) s += (char)('1' + (n - 1 - j));
                        p++;
                    }
                    s = fill(s, order[n - 1 - i]);
                    vs.PB(s);
                }
                reverse(ALL(vs));
                rv.PB(vs);
            }
        }
    } while (next_permutation(order.begin() + 1, order.end()));
    
    return rv;
}

int countMatch(string &data, string s, int err) {
    int no = 0;
    for (int i = 0; i + s.S <= data.S; i++) {
        int noerr = 0;
        REP(j, s.S) {
            if (s[j] != data[i + j]) {
                noerr++;
                if (noerr > err || isdigit(s[j])) goto next;
            }
        }
        no++;
        next: ;
    }
    return no;
}

double countMatchValue(string &data, string s) {
    double rv = 0;
    
    unsigned int ssig = 0;
    REP(i, s.S) if (isdigit(s[i])) ssig += 1ul << i;
    
    unsigned int dsig = 0;
    REP(i, s.S - 1) if (isdigit(data[i])) dsig += 1ul << (i + 1);
    
    VI dig(s.S);
    REP(i, s.S) dig[i] = isdigit(s[i]);
    
    unsigned int asig = 1ul << (s.S - 1);
    
	double xError = gErrorSet ? gError + (1.0 - gError) * 0.1 : 0.5;
	//double xError = 0.5;
	VD mv;
	mv.PB(1.0);
	while (mv[mv.S - 1] > 0.01) mv.PB(mv[mv.S - 1] * (1.0 - xError));
	
    for (int i = 0; i + s.S <= data.S; i++) {
        dsig >>= 1;
        if (isdigit(data[i + s.S - 1])) dsig += asig;
        if (dsig != ssig) continue;
    
        int noerr = 0;
        REP(j, s.S) {
            if (s[j] != data[i + j]) {
				noerr++;
                if (noerr >= mv.S) goto next;
            }
        }
        rv += mv[noerr];
        next: ;
    }
    return rv;
}

double countMatchValue(string &data, string &s, VI &pos) {
    double rv = 0;
    
	double xError = gErrorSet ? gError + (1.0 - gError) * 0.1 : 0.5;
	//double xError = 0.5;
	VD mv;
	mv.PB(1.0);
	while (mv.S < s.S) mv.PB(mv[mv.S - 1] * (1.0 - xError));
	
	REP(k, pos.S) {    
		int i = pos[k];
        int noerr = 0;
        REP(j, s.S) {
            if (s[j] != data[i + j])
				noerr++;
        }
        rv += mv[noerr];
    }
    return rv;
}

int FMVC[100000];
VI findMatchNormal(string &data, string &s, int err, int mx = -1) {
    if (mx == -1) mx = data.S;
	int FMVCNo = 0;
    for (int i = 0; i + s.S <= mx; i++) {
        int noerr = 0;
        REP(j, s.S) {
            if (s[j] != data[i + j]) {
                noerr++;
                if (noerr > err || isdigit(data[i + j]) || isdigit(s[j])) goto next;
            }
        }
		FMVC[FMVCNo++] = i;
        next: ;
    }
    return VI(FMVC, FMVC + FMVCNo);
}

VI findMatchSig(string &data, string &s, int err, int mx = -1) {
    if (data.S < s.S) return VI();

	int FMVCNo = 0;
    if (mx == -1) mx = data.S;
    
    unsigned int ssig = 0;
    REP(i, s.S) if (isdigit(s[i])) ssig += 1ul << i;
    
    unsigned int dsig = 0;
    REP(i, s.S - 1) if (isdigit(data[i])) dsig += 1ul << (i + 1);
    
    VI dig(s.S);
    REP(i, s.S) dig[i] = isdigit(s[i]) ? (1 << 20) : 1;
    
    unsigned int asig = 1ul << (s.S - 1);
    
    for (int i = 0; i + s.S <= mx; i++) {
        dsig >>= 1;
        if (isdigit(data[i + s.S - 1])) dsig += asig;
        if (dsig != ssig) continue;
    
        int noerr = 0;
        REP(j, s.S) {
            if (s[j] != data[i + j]) {
                noerr += dig[j];
                if (noerr > err) goto next;
            }
        }
		FMVC[FMVCNo++] = i;
        next: ;
    }
    return VI(FMVC, FMVC + FMVCNo);
}

VI findMatchNoDig(string &data, string &s, int err, int mx = -1) {
    if (mx == -1) mx = data.S;
	int FMVCNo = 0;
    for (int i = 0; i + s.S <= mx; i++) {
        int noerr = 0;
        REP(j, s.S) {
            if (s[j] != data[i + j]) {
                noerr++;
                if (noerr > err) goto next;
            }
        }
		FMVC[FMVCNo++] = i;
        next: ;
    }
    return VI(FMVC, FMVC + FMVCNo);
}

VI findMatch(string &data, string &s, int err, int mx = -1) {
    bool dig = false;
    REP(i, s.S) if (isdigit(s[i])) {
        dig = true;
        break;
    }
    if (!dig) 
        return findMatchNoDig(data, s, err, mx);
    else if (s.S <= 32)
        return findMatchSig(data, s, err, mx);
    else
        return findMatchNormal(data, s, err, mx);
}


int FMSErr[100100][8];
VI findMatchSave(string &data, string &s, VI &size, int err, int pos) {
	int FMVCNo = 0;
    for (int i = 0; i + size[pos] <= data.S; i++) {
        int noerr = 0;
		int p = i;
        REP(j, s.S) {
			if (isdigit(s[j])) {
				int d = s[j] - '1';
				noerr += FMSErr[p][d];
				p += size[d];
			} else {
				noerr += s[j] != data[p];
				p++;
			}
        }
		FMSErr[i][pos] = noerr;
		if (noerr <= err) FMVC[FMVCNo++] = i;
    }
    return VI(FMVC, FMVC + FMVCNo);
}


VI findMatch(string &data, string &s, int err, VI &p) {
	int FMVCNo = 0;
    REP(k, p.S) {
        int i = p[k];
        int noerr = 0;
        REP(j, s.S) {
            if (s[j] != data[i + j]) {
                noerr++;
                if (noerr > err || isdigit(data[i + j]) || isdigit(s[j])) goto next;
            }
        }
        FMVC[FMVCNo++] = i;
        next: ;
    }
    return VI(FMVC, FMVC + FMVCNo);
}


VI findMatchRR(string &data, string &s, int err) {
    VI p = findMatch(data, s, err);
	if (p.S == 0) return p;
	
	VI rv;
	rv.reserve(p.S);
	
	int cur = p[0];
    FOR(i, 1, p.S) {
        if (cur + s.S > p[i]) {
            int no1 = 0;
            int no2 = 0;
            REP(j, s.S) no1 += data[cur + j] == s[j];
            REP(j, s.S) no2 += data[p[i] + j] == s[j];
            if (no1 < no2)
				cur = p[i];
        } else {
			rv.PB(cur);
			cur = p[i];
		}
    }
	rv.PB(cur);
	
    return rv;
}

VI findMatchRRSave(string &data, string &s, VI &size, int err, int pos) {
    VI p = findMatchSave(data, s, size, err, pos);
	if (p.S == 0) return p;
	
	VI rv;
	rv.reserve(p.S);
	
	int cur = p[0];
    FOR(i, 1, p.S) {
        if (cur + s.S > p[i]) {
            int no1 = 0;
            int no2 = 0;
            REP(j, s.S) no1 += data[cur + j] == s[j];
            REP(j, s.S) no2 += data[p[i] + j] == s[j];
            if (no1 < no2)
				cur = p[i];
        } else {
			rv.PB(cur);
			cur = p[i];
		}
    }
	rv.PB(cur);
	
    return rv;
}

int findBestMatch(string &data, string &s) {
    int bv = 1 << 20;
    int rv = -1;
    for (int i = 0; i + s.S <= data.S; i++) {
        int noerr = 0;
        REP(j, s.S) {
            if (s[j] != data[i + j]) {
                if (isdigit(s[j]) || isdigit(data[i + j])) goto next;
                noerr++;
                if (noerr > bv) goto next;
            }
        }
        if (noerr < bv || noerr == bv && rv < 16) {
            bv = noerr;
            rv = i;
        }
        next: ;
    }
    return rv;
}

VI countChar(string &data, VI &p, int off) {
    VI cnt(128);
    REP(i, p.S) if (p[i] + off >= 0 && p[i] + off < data.S) 
        cnt[data[p[i] + off]]++;
    return cnt;
}

int CME[MAX_JOIN + 10];
double CMEDP[MAX_JOIN + 10];
int calcMaxErrors(double errorRate, int size) {
	if (!gErrorSet) return size / 2;

    if (CME[size]) return CME[size];
    
    errorRate = errorRate + (1.0 - errorRate) * 0.1;

	memset(CMEDP, 0, sizeof(double) * (size + 1));
    CMEDP[0] = 1.0;
    REP(i, size) {
		for (int j = i; j >= 0; j--) {
			CMEDP[j+1] += CMEDP[j] * errorRate;
			CMEDP[j] *= (1.0 - errorRate);
		}
		
		if (CME[i+1] == 0) {
			double sum = 0;
			REP(j, i - 2) {
				sum += CMEDP[j];
				if (sum > 0.995) {
					CME[i+1] = max(1, j);
					break;
				}
			}
			if (CME[i+1] == 0) CME[i+1] = max(1, i - 2);
		}
	}
    return CME[size];
}

string improveMatch(string &data, VI &p, int size, int off = 0) {
    string rv = "";
    REP(i, size) {
        VI cnt = countChar(data, p, i + off);
        int bc = 0;
        FOR(j, '2', '8' + 1) if (cnt[j] > cnt[bc]) bc = j;
        FOR(j, 'a', 'z' + 1) if (cnt[j] > cnt[bc]) bc = j;
        rv += (char)bc;
    }
    return rv;
}

string improveMatch(string &data, string s, int maxTries = 10) {
    string r = s;
    REP(i, maxTries) {
        VI p = findMatch(data, r, calcMaxErrors(gError, s.S), min((int)data.S, 20000));
        string nr = improveMatch(data, p, r.S);        
        if (nr == r) break;
        r = nr;
        //cerr << p.S << ' ' << r << endl;
    }
    return r;
}

void calcGlobalError(string &data, string &s, double v = 0.5) {
    if (gErrorSet) return;
    
    s = improveMatch(data, s);
    VI p = findMatch(data, s, (int)(s.S * v + 1e-9));
    double sum = 0;
    REP(i, s.S) {
        VI cnt = countChar(data, p, i);
        int bc = 0;
        REP(j, 128) if (cnt[j] > cnt[bc]) bc = j;
        sum += (double)cnt[bc] / p.S;
    }
    gError = 1.0 - sum / s.S;
    
    if ((p.S == 1 || gError > 0.25) && v < 0.6) {
        gErrorSet = false;
        calcGlobalError(data, s, 0.75);
    } else {
        gErrorSet = true;
        DB(gError);
    }
    
}

VS bestRV;
int bestRes = 0;

VI calcSizes(VS vs) {
    VI rv(vs.S, -1);
    
    while (true) {
        bool done = true;
        REP(i, vs.S) if (rv[i] == -1) {
            int v = 0;
            REP(j, vs[i].S) {
                if (isdigit(vs[i][j])) {
                    int d = vs[i][j] - '1';
                    if (rv[d] < 0) {
                        v = -1;
                        break;
                    }
                    v += rv[d];
                } else
                    v++;
            }
            if (v >= 0)
                rv[i] = v;
            else
                done = false;
        }
        if (done) break;
    }
    
    REP(i, rv.S) if (rv[i] < 0) rv[i] = 0;
    return rv;
}

__attribute__((aligned(16))) short SOXV[10000][32];
__attribute__((aligned(16))) int SOXStack[10000];
int simulateOptimalX(VS &rv, int skipSize = 0) {
    int SOXStackPos = 0;
    int no = 0;
	VI pos(rv.S + 1);
	pos[0] = 0;
	REP(i, rv.S) pos[i+1] = pos[i] + rv[i].S;
		
	int abc = 4 * pos[rv.S];
	//REP(i, abc) ((int*)SOXV)[i] = 0;
    __m128i* m = (__m128i*)SOXV;
    __m128i zero = _mm_setzero_si128();
    for (int i = 0; i < abc; i += 2) {
        m[i]   = zero;
        m[i+1] = zero;
    }
    for (int i = rv[0].S - 1; i >= 0; i--) SOXStack[SOXStackPos++] = isdigit(rv[0][i]) ? -(rv[0][i] - '1') : i;
    while (no < gData.S && SOXStackPos) {
        SOXStackPos--;
        if (SOXStack[SOXStackPos] < 0) {
            int x = -SOXStack[SOXStackPos];
            for (int i = rv[x].S - 1; i >= 0; i--)
                SOXStack[SOXStackPos++] = isdigit(rv[x][i]) ? -(rv[x][i] - '1') : pos[x] + i;
        } else {
            SOXV[SOXStack[SOXStackPos]][gData[no]-'a']++;
            no++;
        } 
    }
    
    int matched = 0;
    FOR(i, skipSize, rv.S) REP(j, rv[i].S) if (!isdigit(rv[i][j])) {
        short *pv = SOXV[pos[i]+j];
		__m128i *m1 = (__m128i*)&pv[0];
		__m128i *m2 = (__m128i*)&pv[8];
		__m128i *m3 = (__m128i*)&pv[16];
		__m128i *m4 = (__m128i*)&pv[24];
		__m128i m5 = _mm_max_epi16(*m1, *m2);
		__m128i m6 = _mm_max_epi16(*m3, *m4);
		m1[0] = _mm_max_epi16(m5, m6);
		int mx = pv[0];
		if (mx < pv[1]) mx = pv[1];
		if (mx < pv[2]) mx = pv[2];
		if (mx < pv[3]) mx = pv[3];
		if (mx < pv[4]) mx = pv[4];
		if (mx < pv[5]) mx = pv[5];
		if (mx < pv[6]) mx = pv[6];
		if (mx < pv[7]) mx = pv[7];
		/*
        for (int k = 1; k < 26; k += 5) {
            if (pv[k+0] > mx) mx = pv[k+0];
            if (pv[k+1] > mx) mx = pv[k+1];
            if (pv[k+2] > mx) mx = pv[k+2];
            if (pv[k+3] > mx) mx = pv[k+3];
            if (pv[k+4] > mx) mx = pv[k+4];
        }
		*/
        matched += mx;
    }
    
    return matched;
}

int simulateOptimal(string &data, VS &rv, bool update = false, int skipSize = 0) {
    int SOXStack[32 * 8];
    int SOXStackPos = 0;
    int no = 0;
    short v[64 * 8][26];
    memset(v, 0, sizeof (short) * 26 * 64 * rv.S);
    
    for (int i = rv[0].S - 1; i >= 0; i--) SOXStack[SOXStackPos++] = isdigit(rv[0][i]) ? -(rv[0][i] - '1') : i;
    while (no < data.S && SOXStackPos) {
        SOXStackPos--;
        if (SOXStack[SOXStackPos] < 0) {
            int x = -SOXStack[SOXStackPos];
            for (int i = rv[x].S - 1; i >= 0; i--)
                SOXStack[SOXStackPos++] = isdigit(rv[x][i]) ? -(rv[x][i] - '1') : (x<<6) + i;
        } else {
            v[SOXStack[SOXStackPos]][data[no]-'a']++;
            no++;
        } 
    }
    //_mm_max_epi16
    int matched = 0;
    FOR(i, skipSize, rv.S) REP(j, rv[i].S) if (!isdigit(rv[i][j])) {
        int bc = 0;
        short *pv = v[(i<<6)+j];
        for (int k = 1; k < 26; k += 5) {
            if (pv[k+0] > pv[bc]) bc = k+0;
            if (pv[k+1] > pv[bc]) bc = k+1;
            if (pv[k+2] > pv[bc]) bc = k+2;
            if (pv[k+3] > pv[bc]) bc = k+3;
            if (pv[k+4] > pv[bc]) bc = k+4;
        }
        matched += pv[bc];
        if (update) rv[i][j] = 'a' + bc;
    }
    
    return matched;
}

void simulateUpdMatch(string &data, VS &rv) {
    int SOXStack[32 * 8];
    int SOXStackPos = 0;
    int no = 0;
    
    gCorrectMatches[0].PB(0);
    for (int i = rv[0].S - 1; i >= 0; i--) SOXStack[SOXStackPos++] = rv[0][i];
    
    while (no < data.S && SOXStackPos) {
        SOXStackPos--;
        if (isdigit(SOXStack[SOXStackPos])) {
            int x = SOXStack[SOXStackPos] - '1';
            gCorrectMatches[x].PB(no);
            for (int i = rv[x].S - 1; i >= 0; i--)
                SOXStack[SOXStackPos++] = rv[x][i];
        } else {
            no++;
        } 
    }
}

int simulate(string &data, VS rv, int id = 0, bool earlyExit = true) {
    char SOXStack[32 * 8];
    int SOXStackPos = 0;
    int no = 0;
    int matched = 0;
    for (int i = rv[id].S - 1; i >= 0; i--) SOXStack[SOXStackPos++] = rv[id][i];

    VI sizes;
    if (id != 0) sizes = calcSizes(rv);
    VI curLev;
    curLev.PB(0);
    while (no < data.S && SOXStackPos) {
        SOXStackPos--;
        if (SOXStack[SOXStackPos] == data[no]) {
            matched += isdigit(data[no]) ? sizes[data[no] - '1'] : 1;
            no++;
        } else if (isdigit(SOXStack[SOXStackPos])) {
            int x = SOXStack[SOXStackPos] - '1';
            curLev.PB(x);
            SOXStack[SOXStackPos++] = '-';
            for (int i = rv[x].S - 1; i >= 0; i--)
                SOXStack[SOXStackPos++] = rv[x][i];
        } else if (SOXStack[SOXStackPos] == '-') {
            curLev.pop_back();
        } else {
            no++;
            if (earlyExit && matched + data.S - no <= bestRes) return 0;
        }
    }
    return matched;
}

void addResult(VS rv, bool ee = true) {
	int newRes = simulateOptimal(gData, rv, true);
    if (newRes > bestRes) {
        bestRV = rv;
        bestRes = newRes;
        cerr << newRes << " - "; print(rv);
#ifdef EARLY_EXIT
        if (ee && bestRes >= gOptimal && (gLM.S > 5 || gOptimal > 0.25 * gData.S))
            throw 1;
#endif
    }
}

VS reorderResult(VS &vs) {
    VS rv;
    VI order(gLM.S);
    order[0] = 0;
    rv.PB(vs[0]);
    
    FOR(i, 1, vs.S) FOR(j, 1, vs.S) if (vs[j].S == gLM[i]) {
        rv.PB(vs[j]);
        order[j] = i;
        vs[j] = "";
        break;
    }
    
    REP(i, rv.S) REP(j, rv[i].S) if (isdigit(rv[i][j])) rv[i][j] = '1' + order[rv[i][j] - '1'];
    return rv;
}

VI findOrder(VS &vs) {
    VI rv;
	VI used(vs.S);
    
    while (true) {
        bool done = true;
        REP(i, vs.S) if (!used[i]) {
			bool ok = true;
            REP(j, vs[i].S) {
                if (isdigit(vs[i][j])) {
                    int d = vs[i][j] - '1';
                    if (!used[d]) {
                        ok = false;
                        break;
                    }
                }
            }
			if (ok) {
				used[i] = true;
				rv.PB(i);
                done = false;
			}
        }
        if (done) break;
    }
	
    return rv;
}

VS joinSingle(VS &vs) {
	VS rv = vs;
	VI sizes = calcSizes(rv);
	int mx = 1;
	FOR(i, 2, rv.S) if (sizes[i] > sizes[mx]) mx = i;
	int no = 0;
	REP(i, rv[0].S) if (rv[0][i] == (char)('1' + mx)) no++;
	if (no == 1) {
		REP(i, rv[0].S) if (rv[0][i] == '1' + mx) {
			rv[0] = rv[0].substr(0, i) + rv[mx] + rv[0].substr(i + 1);
			break;
		}
		rv.erase(rv.begin() + mx);
		REP(i, rv.S) REP(j, rv[i].S) if (isdigit(rv[i][j])) if (rv[i][j] > '1' + mx) rv[i][j]--;
	}
	return rv;
}

VS splitSingle(VS &vs) {
	VS rv = vs;
	print(rv);
	if (rv[0].S > gLM[0]) {
		string s1 = rv[0].substr(0, gLM[0] - 1);
		string s2 = rv[0].substr(gLM[0] - 1);
		rv[0] = s1 + (char)('0' + gLM.S);
		rv.PB(s2);
	}
	print(rv);
	return rv;
}

//PIVS genLevels(string &data, VI &lm, double totalTime, VS &order, bool randomize = false) {
PIVS genLevels(string &data, VI &lm, int stepsNo, VS &order, double tempMul = 0.0, bool randomize = false) {
	double timerStart = getTime();

	VS rv(order.S);

	VI restored(order.S);
	VI sizes = calcSizes(order);
	REP(i, order.S) rv[i] = string(sizes[i], '-');
	
	FOR(lev, 1, lm.S) {
	
		int expanded = -1;
	
		REP(i, order.S) if (!restored[i]) {
			bool ok = true;
			REP(j, order.S) if (!restored[j] && i != j) {
				REP(k, order[j].S)	if (order[j][k] == '1' + i) ok = false;
			}
			if (ok) {
				restored[i] = true;
				rv[i] = order[i];
				if (randomize) random_shuffle(ALL(rv[i]));
				expanded = i;
				break;
			}
		}
		
		if (expanded == -1) {
			showError("expanded == -1");
			print(restored);
			print(rv);
		}
		
		VS brv = rv;
	
        int bv = simulateOptimalX(rv);
		int xv = bv;
		int steps = 0;
		int accSteps = 0;
		double timerStep = getTime();
		//double timePerLevel = totalTime / (gLM.S - 1);
		while (true) {
			//double curTime = (getTime() - timerStep) / timePerLevel;
			double curTime = (double)steps / stepsNo;
			if (curTime > 1.0) break;
		
		
			steps++;
			double temp = tempMul * (1.0 - curTime);
			
			VS vs = rv;
            while (true) {
                int a = (curTime < 0.2 && rand() % 2) ? expanded : rand() % vs.S;
				if (!restored[a]) continue;
				
                int b = rand() % vs[a].S;
                int c = rand() % (vs[a].S - 1);
                if (c >= b) c++;
                
                if (!isdigit(vs[a][b]) && !isdigit(vs[a][c])) continue;
                
                if (b < c) { 
                    vs[a][c] = rv[a][b];
                    FOR(i, b + 1, c + 1) vs[a][i - 1] = rv[a][i];
                } else {
                    vs[a][c] = rv[a][b];
                    FOR(i, c, b) vs[a][i + 1] = rv[a][i];
                }
                break;
            }
			int av = simulateOptimalX(vs);
			if (av >= bv || randDouble() < exp((av - bv) / temp)) {
				if (av > xv) {
					xv = av;
					brv = vs;
				}
				bv = av;
				rv = vs;
				accSteps++;
			}
		}
		
		rv = brv;
		
		if (lev > 1) cerr << " : ";
		cerr << "Step #" << lev << " -> " << xv;
		
	}
	
	cerr << endl;
	
	int res = simulateOptimal(gData, rv, true);
	return MP(res, reorderResult(rv));
}

PIVS genSA(string &data, double totalTime, VS order) {
    double timerStart = getTime();
    
    VS rv = order;
	int skipSize = 0;
	int addValue = 0;
	
	
	if (gLM.S <= 5) {
		rv = joinSingle(rv);
		skipSize = 1;
		REP(i, rv[0].S) addValue += !isdigit(rv[0][i]);
	}
	
    int bv = simulateOptimal(gData, rv);
    int xv = bv;
    
    VS xRV = rv;

    int steps = 0;
    int accSteps = 0;
    int impSteps = 0;
    
    while (true) {
        double curTime = (getTime() - timerStart) / totalTime;
        if (curTime >= 1.0) break;
        
        double temp = 10.0 * (1 - curTime);
		
		VS vs = VS(ALL(rv));
		while (true) {
			int a = rand() % vs.S;
			int b = rand() % vs[a].S;
			int c = rand() % (vs[a].S - 1);
			if (c >= b) c++;
			
			if (!isdigit(vs[a][b]) && !isdigit(vs[a][c])) continue;
			
			if (b < c) { 
				vs[a][c] = rv[a][b];
				FOR(i, b + 1, c + 1) vs[a][i - 1] = rv[a][i];
			} else {
				vs[a][c] = rv[a][b];
				FOR(i, c, b) vs[a][i + 1] = rv[a][i];
			}
			break;
		}
				
        steps++;
        int av = simulateOptimalX(vs, skipSize) + addValue;
        if (av >= bv || randDouble() < exp((av - bv) / temp)) {
            if (av > xv) {
                impSteps++;
                xRV = vs;
                xv = av;
                DB(xv);
            }
            accSteps++;
            bv = av;
            rv = vs;
        }
    }
	
	if (gLM.S <= 5) xRV = splitSingle(xRV);
	
	simulateOptimal(data, xRV, true);
	
    DB(steps);
    
    return MP(xv, reorderResult(xRV));
}

PIVS genRandom(string &data, VI &lm, double totalTime, VVS pred = VVS()) {
    double timerStart = getTime();
    VI order;
    REP(i, lm.S) order.PB(lm[i]);
    
    VS rv;
    REP(i, rv.S) sort(ALL(rv[i]));
    int bv = 0;
    int xv = 0;
    
    VS xRV;

    int steps = 0;
    int accSteps = 0;
    int impSteps = 0;
    
    VI xSteps(lm.S);
	
	int N = lm.S;
    
	int type = 0;
	//int skipSize = lm.S <= 5 ? 1 : 0;
	int skipSize = 0;
	int addValue = 0;
	bool useSA = false;
	
	if (pred.S) DB(pred.S);
	
	//srand(1);
    while (true) {
        double curTime = (getTime() - timerStart) / totalTime;
        if (curTime >= 1.0) break;
        
        VS vs;
        int aaa = 0;
        double temp = 0;
        
		if (pred.S) {
			vs = VS(pred[rand() % pred.S]);
			REP(i, N) random_shuffle(ALL(vs[i]));
		} else {
			random_shuffle(order.begin() + 1, order.end());
			vs = VS(N);
			
			REP(i, N) vs[i].reserve(order[i]);				
			
			bool useDig[10];
			ZERO(useDig);
			REP(i, N - 1) {
				if (i && !useDig[i])
					break;
				REP(j, order[i]) {
					int x = rand() % 8;
					if (x > i && x < N) {
						vs[i] += '1' + x;
						useDig[x] = true;
					} else
						vs[i] += '-';
				}
			}
			
			if (vs[N - 2].S == 0) continue;
			
			VI size(N);
			size[N - 1] = order[N - 1];
			for (int i = N - 2; i >= 0; i--)
				REP(j, order[i]) size[i] += vs[i][j] == '-' ? 1 : size[vs[i][j] - '1'];
			
			if (!(size[0] >= data.S && size[0] < data.S * 5))
				continue;
				
			vs[N - 1] = string(order[N - 1], '-');
		}
            
        steps++;
        int av = simulateOptimalX(vs);
        if (av > bv) {
            if (av > xv) {
                impSteps++;
                xRV = vs;
                xv = av;
                DB(xv);
            }
            accSteps++;
            bv = av;
            rv = vs;
        }
    }
	
	simulateOptimal(data, xRV, true);
    DB(steps);
    
    return MP(xv, reorderResult(xRV));
}

void genSimple() {
	VVS pred;
	if (gLM.S <= 5) pred = findSizeMatchings(gData.S, gLM);
	
	double timeLeft = MAX_TIME - (getTime() - startTime);
	
	VS sol;
	
	if (gLM.S == 4) {
		VC < PIVS > patterns(pred.S);
		REP(i, pred.S) patterns[i].Y = pred[i];
		int stepsNo = 5000;
		double temp = 0;
		if (pred.S <= 3) stepsNo *= 2, temp += 2;
		if (pred.S <= 7) stepsNo *= 2, temp += 2;
		if (pred.S <= 15) stepsNo *= 2, temp += 2;
		if (pred.S > 64) stepsNo /= 2;
		if (pred.S > 128) stepsNo /= 2;
		int curBest = 0;
		
		while (patterns.S > 0) {
			DB(patterns.S);
			REP(i, patterns.S) patterns[i] = genLevels(gData, gLM, stepsNo, patterns[i].Y, temp, true);
			sort(ALL(patterns));
			reverse(ALL(patterns));
			if (patterns[0].X > curBest) {
				sol = patterns[0].Y;
				curBest = patterns[0].X;
			}
			int newSize = patterns.S / 2;
			while (patterns.S > newSize) patterns.pop_back();
			stepsNo *= 2;
			temp += 2;
		}
	} else if (gLM.S == 5) {
		DB(pred.S);
		
		/*
		int curBest = 0;
		REP(i, 200) {
			PIVS pivs = genLevels(gData, gLM, 50, pred[rand() % pred.S], true);
			if (pivs.X > curBest) {
				curBest = pivs.X;
				sol = pivs.Y;
			}
		}
		DB(curBest);
		//sol = genRandom(gData, gLM, timeLeft * 0.2).Y;
		//sol = genLevels(gData, gLM, 100000, gCorrect, true).Y;
		//REP(i, sol.S) random_shuffle(ALL(sol[i]));
		*/
		sol = genRandom(gData, gLM, timeLeft * 0.2, pred).Y;
	} else {
		sol = genRandom(gData, gLM, timeLeft * 0.2).Y;
	}
	sol = genSA(gData, MAX_TIME - (getTime() - startTime), sol).Y;
	addResult(sol);
	
	/*
	VI sizes = calcSizes(sol);
	bool ok = true;
	FOR(i, 1, sizes.S) ok &= sizes[i] < 5000;
	if (ok && gLM.S >= 6)
		sol = genLevels(gData, gLM, timeLeft * 0.1, sol).Y;
	*/
	
}

int FDP[100101][33];
VS finalDP(VS &vs) {
	VI order = findOrder(vs);
	VI sizes = calcSizes(vs);
	
	print(vs);
	print(order);
	print(sizes);
	
	REP(o, order.S - 1) {
		int pos = order[o];
		string &s = vs[pos];
		for (int i = 0; i < gData.S; i++) {
			int matched = 0;
			int p = i;
			REP(j, s.S) {
				if (isdigit(s[j])) {
					int d = s[j] - '1';
					matched += FMSErr[p][d];
					p += sizes[d];
				} else {
					matched += s[j] == gData[p];
					p++;
				}
				if (p >= gData.S) break;
			}
			FMSErr[i][pos] = matched;
		}
	}
	
	int limit = gLM[0];
	memset(FDP, -1, (gData.S + 1) * 33 * sizeof (int));
	FDP[0][0] = 0;
	int N = gData.S;
	int M = order.S;
	
	REP(i, N) {
		REP(j, limit) if (FDP[i][j] >= 0) {
			FDP[i + 1][j + 1] = max(FDP[i + 1][j + 1], FDP[i][j] + 1);
			FOR(k, 1, M) {
				int x = min(N, i + sizes[k]);
				FDP[x][j + 1] = max(FDP[x][j + 1], FDP[i][j] + FMSErr[i][k]);
			}
		}
	}
	
	int bl = 0;
	REP(i, limit + 1) if (FDP[N][i] > FDP[N][bl]) bl = i;
	DB(FDP[N][bl]);
	
	//reconstruct;
	string s = "";
	int cv = FDP[N][bl];
	int cp = N;
	if (FDP[N - 1][bl - 1] + 1 == cv) {
		s += '-';
		cp = N - 1;
	} else {
		REP(i, N) FOR(k, 1, M) if (i + sizes[k] >= N && FDP[i][bl - 1] + FMSErr[i][k] == cv) {
			s += '1' + k;
			cp = i;
			goto out;
		}
		showError("DP reconstruction failed (beg)");
		out: ;
	}
	
	for (int i = bl - 1; i > 0; i--) {
		int cv = FDP[cp][i];
		if (cp == 0) break;
		if (FDP[cp - 1][i - 1] + 1 == cv) {
			s += '-';
			cp -= 1;
			goto next;
		}
		FOR(k, 1, M) if (sizes[k] <= cp && FDP[cp - sizes[k]][i - 1] + FMSErr[cp - sizes[k]][k] == cv) {
			s += '1' + k;
			cp -= sizes[k];
			goto next;
		}
		showError("DP reconstruction failed");
		next: ;
	}
	if (cp != 0)
		showError("DP reconstruction failed : cp != 0");
	
	reverse(ALL(s));
	DB(s);
	
	/*
	int sum = 0;	
	int p = 0;
	REP(i, s.S) {
		if (isdigit(s[i])) {
			sum += FMSErr[p][s[i]-'1'];
			p += sizes[s[i]-'1'];
		} else {
			sum++;
			p++;
		}
		p = min(N, p);
		cerr << i << ' ' << p << ' ' << sum << ' ' << FDP[p][i+1] << endl;
	}
	*/
	
	VS rv(ALL(vs));
	rv[0] = s;
	//DB(calcSizes(rv)[0]);
	//print(vs);
	//print(rv);
	DB(simulateOptimal(gData, vs, false));
	DB(simulateOptimal(gData, rv, true));
	
	return rv;
}

PID scoreCand(string &data, string &r) {
	PID value;
	
	int maxErrors = calcMaxErrors(gError, r.S);
	VI pos;
	while (true) {
		pos = pos.S == 0 ? findMatch(data, r, maxErrors) : findMatch(data, r, maxErrors, pos);
		if (pos.S == 0) break;
		
		int overlap = 0;
		REP(i, pos.S - 1) if (pos[i] + r.S > pos[i + 1]) {
			overlap++;
		}
		
		if (overlap * 100 < pos.S) break;
		
		maxErrors--;
		if (maxErrors < 0) {
			pos = VI();
			break;
		}
	}
	
	//VI pos = findMatchRR(data, r, calcMaxErrors(gError, r.S));
	
	value.X = pos.S;
	if (pos.S == 0) return MP(-1, -1);
	value.Y = 0;
	REP(i, r.S) {
		VI cnt = countChar(data, pos, i);
		int bc = 0;
		REP(j, 128) if (cnt[j] > cnt[bc]) bc = j;
		value.Y += (double)cnt[bc] / pos.S;
	}
	value.Y /= r.S;
	return value;
}

char tData[100010];
int tSeen[100000][33];
int freq[17][128][128];
VI curOrder;

void go(string data, VI lmpos, VI lm, bool inf, VS cur) {
    if (getTime() - startTime > TIME_LIMIT) throw 1;
	
    statesVis++;

    int lev = gLM.S - lm.S;
	if (lev == 0) curOrder.clear();
	
    if (lm.S == 1) {
        if (true || inf) {
            int bv = -1;
            VI sizes = calcSizes(cur);
            int cv = 0;
            if (data.S <= lm[0]) {
                bv = data.S;
                cur[0] = data;
            } else {
                REP(i, lm[0]) {
                    string suffix = data.substr(i);
                    FOR(j, 1, gLM.S) {
                        int av = cv + i + simulate(suffix, cur, j, false);
                        if (av > bv) {
                            bv = av;
                            cur[0] = data.substr(0, i) + (char)('1' + j);
                        }
                    }
                    cv += isdigit(data[i]) ? sizes[data[i] - '1'] : 1;
                }
                if (cv > bv)
                    cur[0] = data.substr(0, lm[0]);
            }
            cur[0] = fill(cur[0], lm[0]);
            if (bv == -1)
                showError("bv == -1");
        } else {
            if (lm[0] != data.S) return;
            cur[0] = data;
        }
        addResult(cur);
        return;
    }
    
    //generate matches
    VC < pair < double, pair < string, int > > > vs;
    set<string> usedSimple;
    
    double timerAStart = getTime();
    
    VI cand;
    if (lm.S > 2 && data.S >= FULL_CAND) {
        REP(i, 17) {
            FOR(j, 'a', 'z' + 1) {
                FOR(k, 'a', 'z' + 1) freq[i][j][k] = 0;
                FOR(k, '0', '9' + 1) freq[i][j][k] = 0;
            }
            FOR(j, '0', '9' + 1) {
                FOR(k, 'a', 'z' + 1) freq[i][j][k] = 0;
                FOR(k, '0', '9' + 1) freq[i][j][k] = 0;
            }
        }
        FOR(i, 1, 16) REP(j, (int)data.S - i) freq[i][data[j]][data[j+i]]++;
        
    
        VC < pair < double, int > > ts;
        for (int i = 0; i + 16 <= min((int)lm.S * 32 - 16, (int)data.S); i++) {
            string s = data.substr(i, 16);
			//if (lev > 0 && !hasDigit(s, '1' + curOrder[lev - 1])) continue;
			if (lev > 0 && !hasDigit(s)) continue;
            s = improveMatch(data, s, 1);
            double av = 0;
            REP(j, 16) FOR(k, j + 1, 16) av += freq[k - j][s[j]][s[k]];
			//cerr << s << ' ' << av << endl;
            ts.PB(MP(av, i));
        }
        sort(ALL(ts));
        reverse(ALL(ts));
        if (ts.S) REP(i, ts.S) {
            if (ts[i].X * 1.25 < ts[0].X) break;
            cand.PB(ts[i].Y);
			if (lev <= SHOW_CANDIDATES) cerr << "Cand: " << data.substr(ts[i].Y, 16) << ' ' << ts[i].X << endl;
        }
    } else {
        for (int i = 0; i + 16 <= min((int)lm.S * 32 - 16, (int)data.S); i++)
            cand.PB(i);
    }
    
    //for (int i = 0; i + 16 <= min((int)lm.S * 32 - 16, (int)data.S); i++) {
    REP(c, cand.S) {
	
        int i = cand[c];
        string s = data.substr(i, 16);
		if (lev && !hasDigit(s)) continue;
        if (usedSimple.count(s)) continue;
        usedSimple.insert(s);
        
		double av = 0;
		
		if (data.S > 30000) {
			string es = improveMatch(data, s, 1);
			if (es != s) {
				s = es;
				if (usedSimple.count(s)) continue;
				usedSimple.insert(s);
			}
			av = countMatchValue(data, s);
		} else {
			VI ps = findMatch(data, s, calcMaxErrors(gError, s.S));
			string es = improveMatch(data, ps, s.S);
			if (es != s) {
				s = es;
				if (usedSimple.count(s)) continue;
				usedSimple.insert(s);
			}
			av = countMatchValue(data, s, ps);
		}      
		
		if (lev <= SHOW_CANDIDATES) cerr << "Val: " << s << ' ' << av << endl;
		vs.PB(MP(av, MP(s, i)));
        
    }
    
    timerA += getTime() - timerAStart;
    sort(ALL(vs));
    reverse(ALL(vs));
    
    if (!gErrorSet) {
		calcGlobalError(data, vs[0].Y.X);
#ifdef LEAVE_HARD
		if (gError > 0.68 && gLM.S > 6) throw 1;
#endif
	}
    
    int sizePossible[33];
    VC < pair < PID, string > > sizeMatch[33];
    pair < PID, string > bestMatch[33];
    ZERO(sizePossible);
    int differentSizes = 0;
    int sizesOk = 0;

    FOR(i, 1, lm.S) {
        if (sizePossible[lm[i]]) continue;
        sizePossible[lm[i]] = 1;
        differentSizes++;
    }
        
    map<string, PID> used;
    
    int digExist[10];
    ZERO(digExist);
    FOR(i, 2, gLM.S + 1) digExist[i] = 1;
    REP(i, lmpos.S) digExist[lmpos[i] + 1] = 0;
    
    double timerBStart = getTime();
    FOR(i, 16, 33) bestMatch[i].X.X = -1;
    
    int maxMatches = 0;
    
    REP(m, vs.S) {
        string s = vs[m].Y.X;
        int p0 = vs[m].Y.Y;
        
		
		VI ps = findMatch(data, s, calcMaxErrors(gError, s.S), min((int)data.S, 20000));

        s = improveMatch(data, ps, s.S);
        int p = findBestMatch(data, s);

        if (lev <= SHOW_CANDIDATES) DB(s);
        REP(addLeft, 17) REP(addRight, 17 - addLeft) {
            if (lev == 0 && getTime() - startTime > TIME_LIMIT) throw 1;
            int size = 16 + addLeft + addRight;
            if (!sizePossible[size]) continue;
            if (p0 - addLeft < 0 || p0 - addLeft + size > data.S || p - addLeft < 0 || p - addLeft + size > data.S) continue;
            if (tSeen[p0 - addLeft][size] == statesVis) continue;
            tSeen[p0 - addLeft][size] = statesVis;
            
            string r = data.substr(p - addLeft, size);
			if (data.S < 500) {
				r = improveMatch(data, r);
			} else {
				r = improveMatch(data, ps, r.S, -addLeft);
			}
            
            PID value;
            if (used.count(r)) {
                value = used[r];
            } else {
				value = scoreCand(data, r);
				used[r] = value;
				maxMatches = max(maxMatches, (int)(value.X + 1e-9));
                if (lev <= SHOW_CANDIDATES) cerr << value.Y << ' ' << r << ' ' << value.X << endl;
            }
            
            int minSize = (int)(bestMatch[size].X.X * 0.95);
            int maxSize = (int)(bestMatch[size].X.X * 1.05 + 1);
            if (lm.S == 2) {
                sizeMatch[size].PB(MP(value, r));
                used[r] = MP(-1, -1);
            } else if (value.X > maxSize ||
                    value.X >= bestMatch[size].X.X && value.Y >= bestMatch[size].X.Y || 
                    data.S < FULL_CAND && value.X > bestMatch[size].X.X || 
                    data.S >= FULL_CAND && value.X >= minSize && value.Y > bestMatch[size].X.Y)
                bestMatch[size] = MP(value, r);
        }
        
        FOR(i, 16, 33) if (bestMatch[i].X.X >= 0) {
            sizeMatch[i].PB(bestMatch[i]);
            if (data.S >= FULL_CAND && sizeMatch[i].S >= SIZE_LIMIT) {
                sizePossible[i] = 0;
                sizesOk++;
            }
            used[bestMatch[i].Y] = MP(-1, -1);
            bestMatch[i].X.X = -1;
        }
        if (sizesOk == differentSizes && lm.S != 2) break;
    }

#ifdef SHOW_VALID
	{
		PID pid = scoreCand(data, gCorrect[gCorrectOrder[lev]]);
		if (lev <= SHOW_CANDIDATES) cerr << "!" << pid.Y << ' ' << gCorrect[gCorrectOrder[lev]] << ' ' << pid.X << endl;
	}
#endif	
	
    timerB += getTime() - timerBStart;
    VC < pair < double, string > > vMatch;
    
    double timerCStart = getTime();
    VI sizesUsed;
    for (int size = 32; size >= 16; size--) if (sizeMatch[size].S) {
        sort(ALL(sizeMatch[size]));
        reverse(ALL(sizeMatch[size]));
        
        int no = sizeMatch[size].S;
        if (lm.S > 2) no = min(no, SIZE_LIMIT);        
        REP(i, no) {
            int matches = (int)sizeMatch[size][i].X.X;
            if (matches < maxMatches * 0.9 - 5) continue;
			
			/*
            REP(j, sizesUsed.S) REP(k, sizeMatch[sizesUsed[j]].S) {
                int matches2 = (int)sizeMatch[sizesUsed[j]][k].X.X;
                if (matches2 == matches && sizeMatch[sizesUsed[j]][k].Y.find(sizeMatch[size][i].Y) != string::npos) {
                    //sizeMatch[size][i].X.Y -= 1;
                    //goto bad;
                }
            }*/
			
            vMatch.PB(MP(sizeMatch[size][i].X.X * 1000 + 500 + sizeMatch[size][i].X.Y - 1e9 * i, sizeMatch[size][i].Y));
            continue;
            
            bad:
            sizeMatch[size].erase(sizeMatch[size].begin() + i);
            i--;
        }
        sizesUsed.PB(size);
    }
    
    sort(ALL(vMatch));
    reverse(ALL(vMatch));
    
    REP(m, vMatch.S) {
        double timerDStart = getTime();
        VI newLM(ALL(lm));
        VI newLMPos(ALL(lmpos));
        string s = vMatch[m].Y;
        if (lev <= DEBUG_LEVEL) {
            REP(i, lev) cerr << "  ";
            double v = vMatch[m].X - floor(vMatch[m].X);
            if (v == 0) v = 1;
            cerr << "Lev: " << lev << (gAllCorrect.count(s) ? "!" : "") << " Found: " << s << " Value: " << v << " ";
        }
        
        int pos = -1;
        FOR(i, 1, lm.S) if (lm[i] == s.S) {
            newLM.erase(newLM.begin() + i);
            newLMPos.erase(newLMPos.begin() + i);
            pos = lmpos[i];
            break;
        }
        
        if (pos == -1) {
            showError("pos == -1");
            continue;
        }
        
        if (!inf || lm.S > 2) {
            bool ok = true;
            int digCount[10];
            ZERO(digCount);
            REP(i, s.S) if (isdigit(s[i])) digCount[s[i] - '0']++;
            REP(i, 10) if (digExist[i] && (digCount[i] < 1 || digCount[i] > 4)) ok = false;
            if (!ok) {
                if (lev <= DEBUG_LEVEL) cerr << "Aborted: Wrong Number of Digits" << endl;
                continue;
            }
        }
        VI p;
        
        p = findMatchRR(data, s, calcMaxErrors(gError, s.S));
        if (p.S == 0) {
            if (lev <= DEBUG_LEVEL) cerr << "Aborted: No Matches" << endl;
            continue;
        }
        
        cur[pos] = s;
        VI sizes = calcSizes(cur);
		
        int dpos = 0;
        int ppos = 0;
        REP(i, data.S) {
            if (ppos < p.S && p[ppos] == i) p[ppos++] = dpos;
            dpos += isdigit(data[i]) ? sizes[data[i] - '1'] : 1;
        }
        
        
        if (USE_JOINING && sizes[pos] < MAX_JOIN) {
			VI np = findMatchRRSave(gData, s, sizes, calcMaxErrors(gError, sizes[pos]), pos);
			bool ok = true;
			if (lev == 1) {
				int no = 0;
				REP(i, s.S) no += isdigit(s[i]);
				if (no == 1) ok = false;
			}
			if (ok) p = np;
        }
        
        dpos = 0;
        ppos = 0;
        int xpos = 0;
        int tDataPos = 0;
		string xData = data;
        REP(i, xData.S) {
            if (ppos < p.S && p[ppos] == dpos) {
				tData[tDataPos++] = '1' + pos;
				xpos = p[ppos] + sizes[pos];                    
				ppos++;
			}
			
			int ndpos = dpos + (isdigit(xData[i]) ? sizes[xData[i] - '1'] : 1);
			if (dpos < xpos && ndpos > xpos || ppos < p.S && dpos < p[ppos] && ndpos > p[ppos]) {
				xData = xData.substr(0, i) + cur[xData[i] - '1'] + xData.substr(i + 1);
				i--;
				continue;
			}
				
            if (dpos >= xpos)
                tData[tDataPos++] = xData[i];
			dpos = ndpos;
        }
        tData[tDataPos++] = 0;
        string newData = string(tData);
		
        if (lev <= DEBUG_LEVEL) {
            int ok = 0;
            int mis = 0;
            int bad = 0;
            int p0 = 0;
            int p1 = 0;
            
            if (gCorrectMatches.S) {
                while (p0 < gCorrectMatches[lev].S && p1 < p.S) {
                    if (gCorrectMatches[lev][p0] == p[p1]) {
                        ok++;
                        p0++;
                        p1++;
                    } else if (gCorrectMatches[lev][p0] < p[p1]) {
                        mis++;
                        p0++;
                    } else {
                        bad++;
                        p1++;
                    }
                }
                mis += gCorrectMatches[lev].S - p0;
                bad += p.S - p1;            
            } else {
                ok += p.S;
            }
            
            cerr << "Matches: " << ok << "-" << mis << "-" << bad  << " Size: " << sizes[pos] << " Left: " << newData.S << " ME: " << calcMaxErrors(gError, s.S) << " Time: " << (getTime() - startTime) << endl;
        }
            
        cur[pos] = s;
        
        double timerLevStart = getTime();
		curOrder.PB(pos);
        timerD += getTime() - timerDStart;
        go(newData, newLMPos, newLM, inf, cur);
		curOrder.pop_back();
        timerLev[lev] += getTime() - timerLevStart;
        
        cur[pos] = "";
    }
}

void test() {
    int s[] = {200, 500, 1000, 2000};
    double e[] = {0.32, 0.64, 1.28, 1.92, 2.56};
    
    REP(si, 4) REP(ei, 5) {
        int size = s[si];
        double error = e[ei];
        
        int sum = 0;
        int mn = 1<<20;
        int mx = 0;
        REP(tries, 10000) {
            VI v(size);
            int te = (int)(size * error);
            REP(i, te) v[rand() % size] = rand() % 26;
            int act = 0;
            REP(i, size) act += v[i] == 0;
            mn = min(mn, act);
            mx = max(mx, act);
            sum += act;
        }
        cerr << "Size: " << size << " Error: " << error << " Avg: " << (sum / 10000.0 / size) << " Min: " << (mn * 1.0 / size) << " Max: " << (mx * 1.0 / size) << endl;
    }

    exit(0);
}

class StringCompression{public: VS compress(string data, VI lm) {
	int total = 0;
	int bad = 0;
	
    startTime = getTime();
    
    init();
    
    //findSizeMatchings(data.S, lm);    
    
    gData = data;
    gLM = lm;
    
    VI lmpos; REP(i, lm.S) lmpos.PB(i);
    
#ifdef USE_COMPLEX
    try { 
        SIZE_LIMIT = 1;
        USE_JOINING = false;
        cerr << "SIZE_LIMIT = " << SIZE_LIMIT << ' ' << "USE_JOINING: " << USE_JOINING << ' ' << "ERROR: " << gError << endl;
        go(data, lmpos, lm, data.S == 100000, VS(lm.S, "")); 
        
        if (gError > 0.1) {
            double tError = gError;
            
            USE_JOINING = true;
            SIZE_LIMIT = 1;
            cerr << "SIZE_LIMIT = " << SIZE_LIMIT << ' ' << "USE_JOINING: " << USE_JOINING << ' ' << "ERROR: " << gError << endl;
            go(data, lmpos, lm, data.S == 100000, VS(lm.S, "")); 
            
			ZERO(CME);
            gError = tError * 0.9;
            USE_JOINING = true;
            SIZE_LIMIT = 1;
            cerr << "SIZE_LIMIT = " << SIZE_LIMIT << ' ' << "USE_JOINING: " << USE_JOINING << ' ' << "ERROR: " << gError << endl;
            go(data, lmpos, lm, data.S == 100000, VS(lm.S, "")); 
            
			/*
			ZERO(CME);
            gError = tError + (1.0 - tError) * 0.075;
            USE_JOINING = true;
            SIZE_LIMIT = 1;
            cerr << "SIZE_LIMIT = " << SIZE_LIMIT << ' ' << "USE_JOINING: " << USE_JOINING << ' ' << "ERROR: " << gError << endl;
            go(data, lmpos, lm, data.S == 100000, VS(lm.S, "")); 
            */
            gError = tError;
        }
        
		ZERO(CME);
        USE_JOINING = false;
        SIZE_LIMIT = 2;
        cerr << "SIZE_LIMIT = " << SIZE_LIMIT << ' ' << "USE_JOINING: " << USE_JOINING << ' ' << "ERROR: " << gError << endl;
        go(data, lmpos, lm, data.S == 100000, VS(lm.S, "")); 
        
        USE_JOINING = false;
        SIZE_LIMIT = 8;
        cerr << "SIZE_LIMIT = " << SIZE_LIMIT << ' ' << "USE_JOINING: " << USE_JOINING << ' ' << "ERROR: " << gError << endl;
        go(data, lmpos, lm, data.S == 100000, VS(lm.S, "")); 
    } catch (int e) { }
    
    
    DB(statesVis);
    DB(timerA);
    DB(timerB);
    DB(timerC);
    DB(timerD);
    print(VD(timerLev, timerLev + lm.S - 1));
#endif
    	
    try {
#ifdef USE_SIMPLE
		int actRes = bestRes;
		VS curRV = bestRV;
        addResult(chooseBestRep(data, lm));
		if (actRes) addResult(finalDP(curRV));
#ifdef EARLY_EXIT
        if (bestRes >= gOptimal && (gLM.S >= 6 || gRealError > 0.65)) throw 1;
#endif
		genSimple();
		if (gData.S <= 20000) addResult(finalDP(bestRV));
#endif
    } catch (int e) { }
    
    cerr << "FinalScore = " << bestRes << endl;
    cerr << "TotalTime = " << (getTime() - startTime) << endl;
	
    return bestRV;
}

void setCorrect(string data, VS correct) {
    init();	

    gData = data;
    gCorrect = correct;
    
    gCorrectMatches = VVI(correct.S);    
    simulateUpdMatch(gData, gCorrect);
    VVI vvi = gCorrectMatches;
	gCorrectOrder = findOrder(gCorrect);
    gCorrectMatches.clear();
    while (vvi.S) {
        int mx = 0;
        FOR(i, 1, vvi.S) if (vvi[i].S > vvi[mx].S) mx = i;
        gCorrectMatches.PB(vvi[mx]);
        vvi.erase(vvi.begin() + mx);
    }
    
    gOptimal = simulateOptimal(gData, gCorrect);
    gRealError = (double)gOptimal / gData.S;
	cerr << "Real Error = " << (1.0 - (double)gOptimal / data.S) << endl;
    DB(simulateOptimal(gData, gCorrect));
    
    VI gMatchesNo(gCorrect.S);
    REP(i, gMatchesNo.S) gMatchesNo[i] = gCorrectMatches[i].S;
    print(gMatchesNo);
    
    VI order; REP(i, correct.S) order.PB(i);
    do {
        FOR(i, 1, order.S) if (gCorrect[order[i]].S != gCorrect[i].S) goto next;
        REP(i, order.S) {
            string s = "";
            REP(j, correct[order[i]].S) {
                char c = correct[order[i]][j];
                s += isdigit(c) ? '1' + order[c - '1'] : c;
            }
            gAllCorrect.insert(s);
        }
        next: ;
    } while (next_permutation(order.begin() + 1, order.end()));
}};


#ifdef
char s[100010];
int main() {
	scanf("%s", s);
	int n;
	scanf("%d", &n);
	int x[10];
	REP(i, n) scanf("%d", &x[i]);
	
	VS cor(n);
	REP(i, n) cin >> cor[i];//scanf("%s", &cor[i]);
	
	StringCompression algo;
	algo.setCorrect(string(s), cor);
	VS rv = algo.compress(string(s), VI(x, x + n));
	cout << rv.S << endl;
	REP(i, rv.S) cout << rv[i] << endl;
	
	fflush(stdout);
}

#endif
