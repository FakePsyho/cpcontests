// Author: Psyho
// Twitter: https://twitter.com/fakepsyho

// TEMPLATE

#pragma GCC optimize "Ofast,omit-frame-pointer,inline,fast-math,unroll-all-loops,tree-loop-vectorize,tree-slp-vectorize"
#pragma GCC option("arch=native","tune=native","no-zero-upper")
#pragma GCC target("avx,avx2,f16c,fma,sse3,ssse3,sse4.1,sse4.2,bmi,bmi2,lzcnt,popcnt")

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

#ifdef STATS
#define INC(x) x++
#else
#define INC(x)
#endif
 
double get_time() {return chrono::duration_cast<chrono::duration<double>>(chrono::high_resolution_clock::now().time_since_epoch()).count();}
double start_time = get_time();
double elapsed() {return get_time() - start_time;}

// PARAMETERS

const int DIFF_STARTS = 20;
 
const double TIME_SCALE = 1.0;

#ifdef LOCAL
const double TIME_LIMIT = 1.90 * TIME_SCALE;
#else
const double TIME_LIMIT = 3.70;
#endif

const LL STEPS_LIMIT = 0;


// SOLUTION

int buffer_pos = 0;
char buffer[50 << 20];

int next_int() {
    int res = 0;
    while (buffer[buffer_pos] < '0' || buffer[buffer_pos] > '9') buffer_pos++;
    while (buffer[buffer_pos] >= '0' && buffer[buffer_pos] <= '9') res = res * 10 + buffer[buffer_pos++] - '0';
    return res;
}

void write_int(int x) {
    if (x == 0) {buffer[buffer_pos++] = '0'; return;}
    int len = 0;
    while (x) {buffer[buffer_pos + len++] = '0' + x % 10; x /= 10;}
    reverse(&buffer[buffer_pos], &buffer[buffer_pos + len]);
    buffer_pos += len;
}

const int MAX_N = 2000;
const int MAX_E = 12000;
const int MAX_M = 50000;
const int MAX_L = 10;
const int MAX_SN = 50;
const int MAX_SE = 250;

// graph data
int N, E, M;

uint8_t labels[MAX_N];
short cl[MAX_E*2];
short p_cl[MAX_N+1];
short cl_lab[MAX_E*2];
short p_cl_lab[MAX_N*MAX_L+1];
uint8_t cn[MAX_N][MAX_N/8+1];

LL quick_check[MAX_N*MAX_L];

int type1[MAX_N];
LL type2[MAX_N];
short type1list[MAX_N*(1<<MAX_L)];
short n_type1list[MAX_L*(1<<MAX_L)];
int p_type1list[MAX_L*(1<<MAX_L)+1];

// subgraph data
short sub_N[MAX_M];
short sub_E[MAX_M];
uint8_t sub_labels[MAX_M][MAX_SN];
uint8_t sub_cl[MAX_M][MAX_SE*2];
short sub_p_cl[MAX_M][MAX_SN+1];
LL sub_cn[MAX_M][MAX_SN];
int sub_type1[MAX_M][MAX_SN];
LL sub_type2[MAX_M][MAX_SN];

// solution data
int found = 0;
int sol_sub[MAX_M];
int sol[MAX_M][MAX_SN];

INLINE void set_cn(int a, int b) {cn[a][b>>3] |= 1 << (b&7);}
INLINE bool get_cn(int a, int b) {return cn[a][b>>3] & (1 << (b&7));}

bool verify() {
    REP(i, found) {
        int o = sol_sub[i];
        REP(j, sub_N[o]) REP(k, sub_N[o]) if ((sub_cn[o][j] & (1LL << k)) && !get_cn(sol[i][j], sol[i][k])) return false;
    }
    return true;
}


uint8_t sub_start[MAX_M][DIFF_STARTS];
LL total_steps = 0;

// go data
int steps_left;
int n_prev_conn[MAX_SN];
int prev_conn[MAX_SN][MAX_SN];
int gorder[MAX_SN];
char vs[MAX_N];

int go_o;
short go_vorder[MAX_SN];
short go_labels[MAX_SN];
LL go_type2[MAX_SN];

template<bool FIRST> void go(int level) {
    if (steps_left <= 0) throw 2;

    if (level == sub_N[go_o]) {
        sol_sub[found] = go_o;
        REP(i, sub_N[go_o]) sol[found][go_vorder[i]] = gorder[i];
        found++;
        throw 1;
    }

    if (FIRST) {
        int p = go_labels[level]*(1<<MAX_L)+(go_type2[level]&1023);
        steps_left -= p_type1list[p+1]-p_type1list[p];
        FOR(i, p_type1list[p], p_type1list[p+1]) {
            auto gv = type1list[i];
            if (vs[gv] || (type2[gv] & go_type2[level]) != go_type2[level]) continue;

            vs[gv] = true;
            gorder[level] = gv;
            go<false>(level+1);
            vs[gv] = false;
        }
    } else {
        int p = gorder[prev_conn[level][0]]*MAX_L+go_labels[level];
        steps_left -= p_cl_lab[p+1] - p_cl_lab[p];
        if ((quick_check[p] & go_type2[level]) != go_type2[level]) return;
        FOR(i, p_cl_lab[p], p_cl_lab[p+1]) {
            auto gv = cl_lab[i];
            if (vs[gv] || (type2[gv] & go_type2[level]) != go_type2[level]) continue;

            FOR(j, 1, n_prev_conn[level]) if (!get_cn(gv, gorder[prev_conn[level][j]])) {goto next;}

            vs[gv] = true;
            gorder[level] = gv;
            go<false>(level+1);
            vs[gv] = false;
            next: ;
        }

    }
}


int main(int argc, char **argv) {
    ios_base::sync_with_stdio(false); cin.tie(0);

    cin.read(buffer, sizeof(buffer));
    DB(elapsed());

    N = next_int();
    E = next_int();

    REP(i, N) labels[i] = next_int()-1;

    PII edge_data[MAX_E];
    REP(i, E) {
        int a = next_int()-1;
        int b = next_int()-1;
        edge_data[i] = MP(a, b);
        p_cl[a]++;
        p_cl[b]++;
        p_cl_lab[a*MAX_L+labels[b]]++;
        p_cl_lab[b*MAX_L+labels[a]]++;
        set_cn(a, b);
        set_cn(b, a);
    }

    REP(i, N-1) p_cl[i+1] += p_cl[i];
    for (int i = N; i > 1; i--) p_cl[i] = p_cl[i-2];
    p_cl[0] = 0;
    p_cl[1] = 0;

    REP(i, N*MAX_L) p_cl_lab[i+1] += p_cl_lab[i];
    for (int i = N*MAX_L; i > 1; i--) p_cl_lab[i] = p_cl_lab[i-2];
    p_cl_lab[0] = 0;
    p_cl_lab[1] = 0;

    REP(i, E) {
        int a = edge_data[i].X;
        int b = edge_data[i].Y;
        cl[p_cl[a+1]++] = b;
        cl[p_cl[b+1]++] = a;
        cl_lab[p_cl_lab[a*MAX_L+labels[b]+1]++] = b;
        cl_lab[p_cl_lab[b*MAX_L+labels[a]+1]++] = a;
    }

    M = next_int();

    REP(i, M) {
        sub_N[i] = next_int();
        sub_E[i] = next_int();

        REP(j, sub_N[i]) sub_labels[i][j] = next_int()-1;

        REP(j, sub_E[i]) {
            int a=next_int()-1;
            int b=next_int()-1;
            edge_data[j] = MP(a, b);
            sub_p_cl[i][a]++;
            sub_p_cl[i][b]++;
            sub_cn[i][a] |= 1LL << b;
            sub_cn[i][b] |= 1LL << a;
        }

        REP(j, sub_N[i]-1) sub_p_cl[i][j+1] += sub_p_cl[i][j];
        for (int j = sub_N[i]; j > 1; j--) sub_p_cl[i][j] = sub_p_cl[i][j-2];
        sub_p_cl[i][0] = 0;
        sub_p_cl[i][1] = 0;

        REP(j, sub_E[i]) {
            int a = edge_data[j].X;
            int b = edge_data[j].Y;
            sub_cl[i][sub_p_cl[i][a+1]++] = b;
            sub_cl[i][sub_p_cl[i][b+1]++] = a;
        }
    }    

    cerr << "[DATA] N = " << N << endl;
    cerr << "[DATA] E = " << E << endl;
    cerr << "[DATA] M = " << M << endl;
    cerr << "[DATA] R = " << 1.0*E/N << endl;

    cerr << "Read time = " << elapsed() << endl;

    REP(i, N) {
        FOR(j, p_cl[i], p_cl[i+1]) type1[i] |= 1 << labels[cl[j]];
        for (int j = type1[i]; j; j = (j-1) & type1[i]) n_type1list[labels[i]*(1<<MAX_L)+j]++;
    }

    REP(i, MAX_L*(1<<MAX_L)) p_type1list[i] = n_type1list[i];

    REP(i, MAX_L*(1<<MAX_L)) p_type1list[i+1] += p_type1list[i];
    for (int i = MAX_L*(1<<MAX_L); i > 1; i--) p_type1list[i] = p_type1list[i-2];
    p_type1list[0] = 0;
    p_type1list[1] = 0;
    
    REP(i, N) for (int j = type1[i]; j; j = (j-1) & type1[i]) type1list[p_type1list[labels[i]*(1<<MAX_L)+j+1]++] = i;
    

    REP(i, N) {
        FOR(j, p_cl[i], p_cl[i+1]) {
            if (type2[i] & (1LL << (labels[cl[j]]+40))) type2[i] |= 1LL << (labels[cl[j]]+50);
            if (type2[i] & (1LL << (labels[cl[j]]+30))) type2[i] |= 1LL << (labels[cl[j]]+40);
            if (type2[i] & (1LL << (labels[cl[j]]+20))) type2[i] |= 1LL << (labels[cl[j]]+30);
            if (type2[i] & (1 << (labels[cl[j]]+10))) type2[i] |= 1 << (labels[cl[j]]+20);
            if (type2[i] & (1 << labels[cl[j]])) type2[i] |= 1 << (labels[cl[j]]+10);
            type2[i] |= 1 << labels[cl[j]];
        }
    }

    REP(i, N) FOR(j, p_cl[i], p_cl[i+1]) quick_check[i*MAX_L+labels[cl[j]]] |= type2[cl[j]];

    REP(i, M) REP(j, sub_N[i])
        FOR(k, sub_p_cl[i][j], sub_p_cl[i][j+1]) sub_type1[i][j] |= 1 << sub_labels[i][sub_cl[i][k]];

    REP(o, M) REP(i, sub_N[o]) {
        auto scl = sub_cl[o];
        auto sp_cl = sub_p_cl[o];
        auto slabels = sub_labels[o];
        auto stype2 = sub_type2[o];
        FOR(j, sp_cl[i], sp_cl[i+1]) {
            if (stype2[i] & (1LL << (slabels[scl[j]]+40))) stype2[i] |= 1LL << (slabels[scl[j]]+50);
            if (stype2[i] & (1LL << (slabels[scl[j]]+30))) stype2[i] |= 1LL << (slabels[scl[j]]+40);
            if (stype2[i] & (1LL << (slabels[scl[j]]+20))) stype2[i] |= 1LL << (slabels[scl[j]]+30);
            if (stype2[i] & (1 << (slabels[scl[j]]+10))) stype2[i] |= 1 << (slabels[scl[j]]+20);
            if (stype2[i] & (1 << slabels[scl[j]])) stype2[i] |= 1 << (slabels[scl[j]]+10);
            stype2[i] |= 1 << slabels[scl[j]];
        }
    }

    cerr << "Precalc time = " << elapsed() << endl;

    VD sub_eratio(M); REP(i, M) sub_eratio[i] = 1.0 * (sub_E[i]+1) / sub_N[i];

    VI status(M);
    const int STATUS_UNSOLVED = 0;
    const int STATUS_SOLVED = 1;
    const int STATUS_IMPOSSIBLE = 2;

    int steps = 8000;

    const int MAX_STEPS = 1e7;

    REP(loop_counter, 1000) {
        if (elapsed() > TIME_LIMIT) break;
        int subs_left = 0; REP(i, M) if (status[i] == STATUS_UNSOLVED) subs_left++;
        if (subs_left == 0) break;

        int parity = (loop_counter / DIFF_STARTS + loop_counter) % 2;

        REP(sub, M) {
            int o = sub;
            if (status[o] != STATUS_UNSOLVED) continue;
            if (elapsed() > TIME_LIMIT) break;

            int vorder[MAX_SN];

            auto scl = sub_cl[o];
            auto sp_cl = sub_p_cl[o];
            auto slabels = sub_labels[o];

            {
                double bv = -1e30;
                int bp = -1;

                if (steps < MAX_STEPS) {
                    int bad[MAX_SN]; 
                    REP(i, sub_N[o]) bad[i] = 0;                    
                    REP(j, min(loop_counter, DIFF_STARTS)) if (sub_start[o][j] != -1) bad[sub_start[o][j]] = 1;

                    REP(i, sub_N[o]) {
                        if (bad[i]) continue;

                        double v = 0;
                        int v1 = n_type1list[slabels[i]*(1<<MAX_L)+sub_type1[o][i]];
                        v += v1 * -1000;
                        v += (sp_cl[i+1] - sp_cl[i]) * 1;
                        if (n_type1list[slabels[i]*(1<<MAX_L)+sub_type1[o][i]] == 0) {bp = i; break;}
                        if (v > bv) bv = v, bp = i;
                    }

                    if (n_type1list[slabels[bp]*(1<<MAX_L)+sub_type1[o][bp]] == 0) {
                        status[o] = STATUS_IMPOSSIBLE;
                        continue;
                    }

                    sub_start[o][loop_counter%DIFF_STARTS] = bp;
                    if (bp == -1) continue;
                } else {
                    bp = loop_counter % sub_N[o];
                }

                short conn[MAX_SN]; memset(conn, 0, sizeof(short) * sub_N[o]);
                short ex_conn[MAX_SN]; memset(ex_conn, 0, sizeof(short) * sub_N[o]);
                short nb_conn[MAX_SN]; REP(i, sub_N[o]) {nb_conn[i]=0; FOR(j, sp_cl[i], sp_cl[i+1]) nb_conn[i] += sp_cl[scl[j]+1] - sp_cl[scl[j]];}
                LL last_conn_bonus[MAX_SN]; memset(last_conn_bonus, 0, sizeof(LL) * sub_N[o]);

                vorder[0] = bp;
                conn[bp] = 3;

                int stack[192];
                int stackpos = 0;

                static int vs[MAX_SN];
                static int vs_cnt = 0;
                if (vs_cnt == 0) ZERO(vs);

                REP(i, sub_N[o] - 1) {
                    FOR(j, sp_cl[vorder[i]], sp_cl[vorder[i]+1]) {conn[scl[j]] |= 1; ex_conn[scl[j]]++; nb_conn[scl[j]] -= sp_cl[vorder[i]+1] - sp_cl[vorder[i]]; last_conn_bonus[scl[j]] += 1LL << i;}
                    double bv = -1e30;
                    int bp = -1;

                    if (parity == 0) {
                        REP(j, sub_N[o]) if (conn[j] == 1) {
                            double v = 0;
                            v += ex_conn[j] * 1e4;
                            v += (nb_conn[j] * 1.0 + 0.1 * last_conn_bonus[j] / (1LL << i)) / n_type1list[slabels[j]*(1<<MAX_L)+sub_type1[o][j]];
                            if (v > bv) bv = v, bp = j;
                        }
                    } else {
                        bool good = false; 
                        REP(j, sub_N[o]) if (conn[j] == 1 && ex_conn[j] > 1) good = true;
                        if (good) {
                            REP(j, sub_N[o]) if (conn[j] == 1) {
                                double v = 0;
                                v += 1e2 * last_conn_bonus[j] * ex_conn[j];
                                v += nb_conn[j] * 0.1 / n_type1list[slabels[j]*(1<<MAX_L)+sub_type1[o][j]];
                                if (v > bv) bv = v, bp = j;
                            }
                        } else {
                            REP(j, sub_N[o]) if (conn[j] == 1) {
                                int loop_len = 100;
                                int chain_len = 0;

                                vs_cnt++;
                                stack[0] = j;
                                stack[1] = 0;
                                stack[2] = -1;
                                vs[j] = vs_cnt;
                                stackpos = 3;
                                while (stackpos) {
                                    int v = stack[stackpos-3];
                                    int d = stack[stackpos-2];
                                    int prev = stack[stackpos-1];
                                    stackpos -= 3;

                                    chain_len = max(chain_len, d);
                                    if (conn[v] == 3) {
                                        if (d > 1) loop_len = min(loop_len, d);
                                        continue;
                                    }

                                    FOR(k, sp_cl[v], sp_cl[v+1]) if (scl[k] != prev) {
                                        if (vs[scl[k]] == vs_cnt) {loop_len=min(loop_len,d+2); continue;}
                                        vs[scl[k]] = vs_cnt;
                                        stack[stackpos+0] = scl[k];
                                        stack[stackpos+1] = d+1;
                                        stack[stackpos+2] = v;
                                        stackpos += 3;
                                    }
                                }

                                double v = 0;
                                v += loop_len * -1e6;
                                v += chain_len * -1e4;
                                if (v > bv) bv = v, bp = j;
                            }
                        }
                    }


                    vorder[i+1] = bp;
                    conn[bp] = 3;
                }
            }
            

            REP(i, sub_N[o]) {
                n_prev_conn[i] = 0;
                REP(j, i) if (sub_cn[o][vorder[i]] & (1LL << vorder[j])) {
                    prev_conn[i][n_prev_conn[i]++] = j;
                }
            }

            go_o = o;
            REP(i, sub_N[o]) {
                go_vorder[i] = vorder[i];
                go_labels[i] = slabels[vorder[i]];
                go_type2[i] = sub_type2[o][vorder[i]];
            }

            try {
                steps_left = steps;
                go<true>(0);
                status[o] = STATUS_IMPOSSIBLE;
            } catch (int e) { 
                if (e == 1) status[o] = STATUS_SOLVED;
                REP(i, sub_N[o]) vs[gorder[i]] = false;
            }
            total_steps += steps - steps_left;
            #ifdef LOCAL
            if (STEPS_LIMIT && total_steps > STEPS_LIMIT) break;
            #endif
        }
        if (loop_counter % 5 == 4) steps *= 4;
        steps = min(steps, MAX_STEPS);
        #ifdef LOCAL
        if (STEPS_LIMIT && total_steps > STEPS_LIMIT) break;
        #endif
    }

#ifdef LOCAL
    assert(verify());
#endif

    int unsolved = 0; REP(i, M) if (status[i] == STATUS_UNSOLVED) unsolved++;
    int solved = 0; REP(i, M) if (status[i] == STATUS_SOLVED) solved++;
    int impossible = 0; REP(i, M) if (status[i] == STATUS_IMPOSSIBLE) impossible++;
    DB(solved, unsolved, impossible);

    cerr << "[DATA] lsteps = " << steps << endl;
    cerr << "[DATA] tsteps = " << total_steps << endl;
    // cerr << "[DATA] slv = " << solved << endl;
    cerr << "[DATA] uslv = " << unsolved << endl;
    // cerr << "[DATA] imp = " << impossible << endl;

    buffer_pos = 0;
    write_int(found);
    buffer[buffer_pos++] = '\n';
    REP(i, found) {
        write_int(sol_sub[i]+1);
        REP(j, sub_N[sol_sub[i]]) {
            buffer[buffer_pos++] = ' ';
            write_int(sol[i][j]+1);
        }
        buffer[buffer_pos++] = '\n';
    }
    buffer[buffer_pos++] = '\0';
    cout << buffer << flush;

    cerr << "Score = " << (int)((found)*1e7/M) << endl;
    cerr << "Time = " << (int)(elapsed()*1000) << endl;

	return 0;
}
