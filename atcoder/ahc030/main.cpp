// Author: Psyho
// Twitter: https://twitter.com/fakepsyho

// TODO

// - investigate speed of hashing
// - full solve during phase 2 (beam search with a few candidates?)
// - use type 2 guess to optimize final guesses!?
// - optimize!
// - simplify/refactor code
// - can we use type 2 during phase 1?
// - calc probabilities for type 2 moves: e.g. if I got x, what's the distribution 

// - try solving via beam search?

// - try silly solution with making bunch of type 2 guesses and HC for finding best solution (minimize MSE?)
// - check if valid solution minimizes MSE
// - mix type 1 & type 2? (how?)
// - smarter shape generation? (randomize 10K solutions and choose ones that differentiate the most?)

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
 
double PI = 2 * acos(0);

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
    INLINE int next(int x) {return ((LL)rand() * x) >> 32; }
    INLINE int next(int a, int b) {return a + next(b - a); }
    INLINE double next_double() {return (rand() + 0.5) * (1.0 / 4294967296.0); }
    INLINE double next_double(double a, double b) {return a + next_double() * (b - a); }
    INLINE double normal() {return sqrt(-2 * log(next_double())) * cos(2 * PI * next_double()); }
};
 

struct FastRNG {
	unsigned long x=123456789, y=362436069, z=521288629;
    FastRNG() { }
	unsigned int rand() { unsigned long t; x ^= x << 16; x ^= x >> 5; x ^= x << 1; t = x; x = y; y = z; z = t ^ x ^ y; return z;}
    INLINE int next() {return rand(); }
    INLINE int next(int x) {return ((LL)rand() * x) >> 32; }
    INLINE int next(int a, int b) {return a + next(b - a); }
    INLINE double next_double() {return (rand() + 0.5) * (1.0 / 4294967296.0); }
    INLINE double next_double(double a, double b) {return a + next_double() * (b - a); }
    INLINE double normal() {return sqrt(-2 * log(next_double())) * cos(2 * PI * next_double()); }
};

static FastRNG rng;

// SOLUTION

double timer1 = 0;
double timer2 = 0;
double timer3 = 0;

const int MIN_N = 10;
const int MAX_N = 20;
const int MAX_M = 20;

const int CONFIG_ALL = false;
const int CONFIG_MUL = 1;
const int MAX_CONFIGS = 2000000 * CONFIG_MUL;


int N, M;
double eps;


uint32_t zobrist[MAX_N*MAX_N];

struct Shape {
    int w;
    int h;
    int len;
    int data[MAX_N*MAX_N];
    int points[MAX_N*MAX_N];
};

Shape shapes[MAX_M];
int guesses[MAX_N*MAX_N];
VI guesses_list;

int n_configs = 1;
uint8_t configs[MAX_CONFIGS][MAX_N*MAX_N];
int config_tiles_missing[MAX_CONFIGS];
int config_stats[MAX_N*MAX_N][MAX_M+4];
VI config_shapes;
int config_tiles_left = 0;

bool config_good[MAX_CONFIGS];


LL hash_test() {
    LL h = 0;
    h ^= N;
    h ^= M << 5;
    h ^= (int)(eps * 100 + 1e-6) << 10;
    REP(i, M) {
        RNG hash_rng(1337 * i);
        VC<LL> v; REP(j, N*N) v.PB((LL)abs(hash_rng.next()));
        REP(j, shapes[i].len) {
            h ^= (LL)v[shapes[i].points[j]] << (12 + 2 * i);
        }
    }
    return h;
}

VI load_gtd() {
    ifstream f("gtd.txt", ios::in);
    LL ht = hash_test();
    while (true) {
        LL h; string a, s;
        f >> h >> a >> s;
        if (h == ht) {
            VS vs = splt(s, ',');
            VI rv; for (string s : vs) rv.PB(atoi(s.c_str()));
            f.close();
            return rv;
        }
    }
    assert(false);
}

INLINE double erf(double z) {
    double t = 1.0 / (1.0 + 0.5 * abs(z));
    double ans = 1 - t * exp(-z*z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
    return z >= 0 ? ans : -ans;
}

INLINE double prob_normal(double mean, double stddev, double x) {
    if (x < 0) x = -99;
    double z = (x - mean) / stddev;
    return 1 - 0.5 * (1 + erf(z / sqrt(2)));
}

INLINE double prob_normal_min(double mean, double stddev, double x) {
    return x > mean ? prob_normal(mean, stddev, x) : prob_normal(mean, stddev, mean + mean - x);
}

INLINE double prob_normal_opt(double mean, double stddev, int x) {
    return x > mean ? prob_normal(mean, stddev, x - 0.5) : 1 - prob_normal(mean, stddev, x + 0.5);
}

INLINE double prob_normal_range(double mean, double stddev, int x) {
    return prob_normal(mean, stddev, x - 0.5) - prob_normal(mean, stddev, x + 0.5);
}

bool config_cmp(int a, int b) {
    REP(p, N*N) if (configs[a][p] != configs[b][p]) return 1;
    return 0;
}

void stats_reset() {
    REP(p, N*N) REP(i, M+1) config_stats[p][i] = 0;
}

void stats_remove(int pos) {
    REP(p, N*N) config_stats[p][configs[pos][p]]--;
}

void stats_add(int pos) {
    REP(p, N*N) config_stats[p][configs[pos][p]]++;
}

void config_move(int src, int dest) {
    memcpy(configs[dest], configs[src], N*N*sizeof(uint8_t));
    config_tiles_missing[dest] = config_tiles_missing[src];
}

void config_update(int p, int v, const bool use_stats = true, int max_possible = 0) {
    int tiles_left = 0;
    REP(i, M) if (find(ALL(config_shapes), i) == config_shapes.end()) tiles_left += shapes[i].len;

    REP(i, n_configs) {
        config_good[i] = v >= configs[i][p] && v <= configs[i][p] + max_possible;
        config_tiles_missing[i] += v - configs[i][p];
        if (config_tiles_missing[i] > tiles_left) config_good[i] = false;
        if (use_stats && !config_good[i]) stats_remove(i);
    }
    int posl = 0;
    int posr = n_configs-1;
    while (posl <= posr) {
        while (posl <= posr && config_good[posl]) posl++;
        while (posl <= posr && !config_good[posr]) posr--;
        if (posl <= posr) {
            config_move(posr, posl);
            posl++;
            posr--;
        }
    }
    n_configs = posl;
}

struct Guess {
    int p;
    int v;
};

const int HASH_BITS = 25;
int hash_table[1<<HASH_BITS];

uint32_t offset_hash[MAX_CONFIGS];

bool config_expand_check(int shape, int offset, int pos) {
    REP(i, shapes[shape].len) configs[pos][shapes[shape].points[i] + offset]++;
    bool rv = false;
    int tiles_missing = 0;
    for (int p : guesses_list) if (guesses[p] < configs[pos][p] || guesses[p] > configs[pos][p] + (M - config_shapes.S - 1)) goto skip;
    for (int p : guesses_list) tiles_missing += guesses[p] - configs[pos][p];
    if (tiles_missing > config_tiles_left - shapes[shape].len) goto skip;
    rv = true;
    skip:
    REP(i, shapes[shape].len) configs[pos][shapes[shape].points[i] + offset]--;
    return rv;
}

void config_expand(int shape, VI offsets, const bool use_stats = true) {
    config_shapes.PB(shape);

    config_tiles_left -= shapes[shape].len;

    static VC<Guess> allg; allg.clear(); REP(p, N*N) if (guesses[p] >= 0) allg.PB({p, guesses[p]});

    int prev_n_configs = n_configs;

    static int hash_offset = 0;
    hash_offset += MAX_CONFIGS;

    static VI offset_guesses(MAX_N*MAX_N, 0);

    REP(i, offsets.S) {
        offset_hash[i] = 0;
        offset_guesses[i] = 0;
        REP(j, shapes[shape].len) {
            offset_hash[i] += zobrist[shapes[shape].points[j] + offsets[i]];
            offset_guesses[i] += guesses[shapes[shape].points[j] + offsets[i]] >= 0;
        }
    }

    REP(i, prev_n_configs) {
        uint32_t og_hash = 0;
        REP(p, N*N) og_hash += configs[i][p] * zobrist[p];

        REP(offset_id, offsets.S) {
            int offset = offsets[offset_id];
            if (config_tiles_missing[i] - offset_guesses[offset_id] > config_tiles_left) continue;

            uint32_t hash = og_hash + offset_hash[offset_id];

            REP(j, shapes[shape].len) {
                configs[i][shapes[shape].points[j] + offset]++;
            }
            
            for (auto &g : allg) {
                if (g.v < configs[i][g.p] || g.v > configs[i][g.p] + M - config_shapes.S) goto next;
            }

            hash &= (1<<HASH_BITS)-1;
            while (hash_table[hash] >= hash_offset) {
                if (config_cmp(hash_table[hash] - hash_offset, i) == 0) goto next;
                hash = (hash + 1) & ((1<<HASH_BITS)-1);
            }
            hash_table[hash] = n_configs + hash_offset;

            config_move(i, n_configs);
            config_tiles_missing[n_configs] = config_tiles_missing[i] - offset_guesses[offset_id];
            if (use_stats) stats_add(n_configs);
            n_configs++;

            next: ;
            REP(j, shapes[shape].len) {
                configs[i][shapes[shape].points[j] + offset]--;
            }
        }
        if (use_stats) stats_remove(i);
    }

    REP(i, prev_n_configs) {
        if (n_configs > prev_n_configs) config_move(n_configs-1, i);
        n_configs--;
    }
}

double pre_eval[200][200];

const int RNG_CONFIGS = 200;
uint16_t samples[MAX_N*MAX_N][RNG_CONFIGS];
double pre_log[RNG_CONFIGS+1];
double pre_pnr_eval[400][400];
double pre_pnr2_eval[400][400];

int main(int argc, char **argv) {
	cin >> N >> M >> eps;

    REP(i, M) {
        Shape &s = shapes[i];
        cin >> s.len;
        REP(j, s.len) {
            int r, c;
            cin >> r >> c;
            s.points[j] = r*N+c;
            s.data[r*N+c] = 1;
            s.w = max(s.w, r+1);
            s.h = max(s.h, c+1);
        }
    }

    double C = 1.0;
    REP(i, M) C *= (N - shapes[i].w + 1) * (N - shapes[i].h + 1);
    double logC = log2(C);

    int tiles = 0; REP(i, M) tiles += shapes[i].len;

    REP(i, N*N) zobrist[i] = rng.next();

    cerr << "[DATA] N = " << N << endl;
    cerr << "[DATA] M = " << M << endl;
    cerr << "[DATA] eps = " << eps << endl;
    cerr << "[DATA] logC = " << logC << endl;
    cerr << "[DATA] cv = " << 1.0 * tiles / N/N << endl;

    int CONFIG_LIMIT = CONFIG_ALL ? MAX_CONFIGS : 250000LL * 20 * 20 / N / N * CONFIG_MUL;
    int totallen = 0; REP(i, M) totallen += shapes[i].len;


    VI gtd;
#ifdef LOCAL
    gtd = load_gtd();
#endif

    int SOL_TYPE = 2;
    if (eps <= 0.12 + 1e-9 && logC < 120) SOL_TYPE = 1;
    if (logC < 20) SOL_TYPE = 0;

    // SOL_TYPE = 1;

    if (SOL_TYPE < 0) {
        cerr << "[DATA] time = " << elapsed()*1000 << endl;
        exit(0);
    }
    
    if (SOL_TYPE == 0) {
        REP(p, N*N) guesses[p] = -1;
        config_tiles_left = 1<<20;
        REP(i, M) {
            VI offsets;
            REP(r, N - shapes[i].w + 1) REP(c, N - shapes[i].h + 1) offsets.PB(r*N+c);
            config_expand(i, offsets, false);
        }
        DB(n_configs);

        VD prob(n_configs, 0);

        const double PROB_CUTOFF = 3;
        const double PROB_EXCLUDED = 15;

        int steps_left = N*N*2;
        
        const int N_CELLS = N*N*25/100;
        const double sigma = sqrt(N_CELLS * eps * (1 - eps));
        REP(a, N*N+1) REP(b, N*N+1) pre_pnr_eval[a][b] = -log2(prob_normal_range(a, sigma, b));
        REP(a, N*N+1) REP(b, N*N+1) pre_pnr2_eval[a][b] = -log2(prob_normal_range((N_CELLS - a) * eps + a * (1 - eps), sigma, b));

        while (true) {
            static VI best_guess;
            best_guess.clear();
            static VI guess; 
            guess.clear();

            static VI left_list;
            left_list.clear();
            double min_prob = 1e10; REP(i, n_configs) min_prob = min(min_prob, prob[i]);
            REP(i, n_configs) if (prob[i] <= min_prob + PROB_CUTOFF) left_list.PB(i);

            double bv = -1e9;
            double xv = -1e9;
            int total_reps = max(1, 500000 / (int)left_list.S);
            int reps = total_reps;
            bool first = true;
            while (reps--) {
                static VI v;
                int type = !first && rng.next(2);
                if (type == 0) {
                    first = false;
                    v.clear();
                    REP(i, N*N) v.PB(i);
                    REP(i, N*N) swap(v[i], v[rng.next(i, N*N)]);
                    while (v.S > N_CELLS) v.pop_back();
                } else {
                    v = guess;
                    int p = rng.next(v.S);
                    while (true) {
                        int x = rng.next(N*N);
                        if (find(ALL(v), x) != v.end()) continue;
                        v[p] = x;
                        break;
                    }
                }
                double av = 0;

                int count[MAX_N*MAX_N+1];
                ZERO(count);
                int minsum = 1e9;
                int maxsum = 0;
                for (int i : left_list) {
                    int sum = 0;
                    for (int p : v) sum += configs[i][p];
                    count[sum]++;
                    minsum = min(minsum, sum);
                    maxsum = max(maxsum, sum);
                }
                FOR(a, minsum, maxsum+1) {
                    double tav = 0;
                    // FOR(b, minsum, maxsum+1) tav += (a - b) * (a - b) * (a - b) * (a - b) * count[b];
                    FOR(b, minsum, maxsum+1) tav += abs(a - b) * (a - b) * (a - b) * count[b];
                    // FOR(b, minsum, maxsum+1) tav += (a - b) * (a - b) * count[b];
                    // FOR(b, minsum, maxsum+1) tav += pre_pnr_eval[a][b] * count[b];
                    av += tav * count[a];
                }
                av /= left_list.S;
                av /= left_list.S;

                if (av >= bv - 1e-9 || type && total_reps > 1000 && rng.next(100) == 0) {
                    bv = av;
                    guess = v;
                    if (av > xv) {
                        xv = av;
                        best_guess = v;
                    }
                }
            }
            
            cout << "q " << N_CELLS;
            for (int p : best_guess) cout << " " << p/N << " " << p%N;
            cout << endl;

            int answer;
            cin >> answer;

            REP(i, n_configs) {
                int sum = 0;
                for (int p : best_guess) sum += configs[i][p];
                prob[i] += pre_pnr2_eval[sum][answer];
                // prob[i] += -log2(prob_normal_range((N_CELLS - sum) * eps + sum * (1 - eps), sigma, answer));
            }

            int best_pos = -1;
            double min1_prob = 1e10;
            double min2_prob = 1e10;

            REP(i, n_configs) {
                if (prob[i] < min1_prob) {
                    min2_prob = min1_prob;
                    min1_prob = prob[i];
                    best_pos = i;
                } else if (prob[i] < min2_prob) {
                    min2_prob = prob[i];
                }
            }

            int n_left = 0; REP(i, n_configs) {
                n_left += prob[i] <= min1_prob + PROB_CUTOFF;
                if (prob[i] > min1_prob + PROB_EXCLUDED) {
                    n_configs--;
                    if (best_pos == n_configs) best_pos = i;
                    prob[i] = prob[n_configs];
                    config_move(n_configs, i);
                    i--;
                }
            }

            DB(min1_prob, min2_prob, best_pos, n_configs, n_left);

            steps_left--;
            bool done = false;
            if (min1_prob + PROB_CUTOFF < min2_prob || steps_left <= 2){
                VI rv; REP(p, N*N) if (configs[best_pos][p]) rv.PB(p);
                cout << "a " << rv.S;
                for (int p : rv) cout << " " << p/N << " " << p%N;
                cout << endl;

                int ans; cin >> ans;
                if (ans) done = true;
                prob[best_pos] = 1e9;
                steps_left--;
            }
            done |= steps_left <= 2;
            done |= min1_prob >= 1e9;
            if (done) {
                cerr << "[DATA] time = " << elapsed()*1000 << endl;
                exit(0);
            } 


        }
    }
    // exit(0);

    VVI t2g;
    VI t2v;

    REP(i, RNG_CONFIGS+1) pre_log[i] = log2(1.0 * i / RNG_CONFIGS);
    FOR(times, 1, 8) {
        if (SOL_TYPE != 1) break;
        const int N_GUESSES = (N*N/4)*times;
        const int N_CELLS = N*N*35/100;
        const int STARTS = N_CELLS;

        double sigma = sqrt(N_CELLS * eps * (1 - eps));

        REP(sum, 200) {
            REP(query, 200) {
                double x = sum * (1 - eps) + (N_CELLS - sum) * eps;
                // pre_eval[sum][query] = -log(prob_normal_range(x, sigma, query));
                pre_eval[sum][query] = -log(prob_normal_range(x, sigma, query));
            }
            int mid = 0;
            REP(query, 200) if (pre_eval[sum][query] < pre_eval[sum][mid]) mid = query;
            
            double prev2 = pre_eval[sum][mid];
            double prev = pre_eval[sum][mid];
            FOR(query, mid+1, 200) {
                if (pre_eval[sum][query] > 1e9) pre_eval[sum][query] = prev + (prev - prev2);
                prev2 = prev;
                prev = pre_eval[sum][query];
            }

            prev2 = pre_eval[sum][mid];
            prev = pre_eval[sum][mid];
            for (int query = mid-1; query >= 0; query--) {
                if (pre_eval[sum][query] > 1e9) pre_eval[sum][query] = prev + (prev - prev2);
                prev2 = prev;
                prev = pre_eval[sum][query];
            }
        }

        // REP(i, RNG_CONFIGS+1) pre_log[i] = i ? (1.0 * i / RNG_CONFIGS) * -log2(1.0 * i / RNG_CONFIGS) : 0.0;

        int acc1 = 0;
        int max_sum = 0;
        REP(i, N*N/4) {
            VI g; 
            // REP(p, N*N) g.PB(p);
            // REP(p, N*N) swap(g[p], g[rng.next(p, N*N)]);
            // while (g.S > N_CELLS) g.pop_back();
            VVI z(N, VI(N, 0));
            REP(j, STARTS) {
                while (true) {
                    int r = rng.next(N);
                    int c = rng.next(N);
                    if (!z[r][c]) {
                        z[r][c] = 1;
                        break;
                    }
                }
            }
            REP(j, N_CELLS-STARTS) {
                VPII v;
                REP(r, N) REP(c, N) if (z[r][c] == 0) {
                    int sum = 0;
                    REP(d, 4) {
                        int nr = r + (d == 2) - (d == 3);
                        int nc = c + (d == 0) - (d == 1);
                        if (nr < 0 || nr >= N || nc < 0 || nc >= N) continue;
                        sum += z[nr][nc];
                    }
                    if (sum) v.PB({r, c});
                }
                PII p = v[rng.next(v.S)];
                z[p.X][p.Y] = 1;
            }
            REP(j, N) REP(k, N) if (z[j][k]) g.PB(j*N+k);

            ZERO(samples);
            REP(i, RNG_CONFIGS) {
                REP(m, M) {
                    int r = rng.next(N - shapes[m].w + 1);
                    int c = rng.next(N - shapes[m].h + 1);
                    REP(j, shapes[m].len) samples[shapes[m].points[j] + r*N+c][i]++;
                }
            }

            int count[MAX_N*MAX_N+1];
            uint16_t sumcfg[RNG_CONFIGS];
            ZERO(sumcfg);
            for (int p : g) REP(j, RNG_CONFIGS) sumcfg[j] += samples[p][j];

            ZERO(count); REP(j, RNG_CONFIGS) count[sumcfg[j]]++;
            // double bv = 0; REP(j, totallen+1) bv += -pre_log[count[j]];
            // double bv = 0; REP(j, totallen) if (count[j]) bv += count[j] * -log2(1.0 * (count[j] + (j ? count[j-1] : 0) + (j>1 ? count[j-2] : 0) + count[j+1] + count[j+2]) / RNG_CONFIGS);
            double bv = -1e9;

            double pre_bv = bv;
            acc1 = 0;
            REP(step, N*N/2) {
                if (logC < 40) break;
                int a = rng.next(g.S);
                int old = g[a];

                int np = rng.next(N*N);
                if (z[np/N][np%N]) continue;

                ZERO(count); REP(j, RNG_CONFIGS) count[sumcfg[j] - samples[old][j] + samples[np][j]]++;
                
                int minsum = 1e9;
                int maxsum = 0;
                REP(j, totallen+1) if (count[j]) minsum = min(minsum, j), maxsum = max(maxsum, j);

                double av = 0;
                FOR(j, minsum, maxsum+1) {
                    double tav = 0;
                    FOR(k, minsum, maxsum+1) tav -= (j - k) * (j - k) * count[k];
                    av += tav * count[j];
                }

                if (av >= bv - 1e-9) {
                    g[a] = np;
                    bv = av;
                    z[old/N][old%N] = 0;
                    z[np/N][np%N] = 1;
                    acc1++;
                    REP(j, RNG_CONFIGS) sumcfg[j] -= samples[old][j], sumcfg[j] += samples[np][j];
                }
            }   

            // DB(pre_bv, bv, acc1, VI(count, count+totallen+1));


            t2g.PB(g);
            cout << "q " << N_CELLS;
            for (int p : g) cout << " " << p/N << " " << p%N;
            cout << endl;
            int x;
            cin >> x;
            t2v.PB(x);
        }
        // DB(acc1);
        // DB(N_GUESSES * 1000);
        // DB(max_sum);
        // DB(elapsed());

        VVI ct2(N*N);
        REP(i, N_GUESSES) for (int p : t2g[i]) ct2[p].PB(i);

        double xv = 1e9;
        double bv = 1e12;
        VI sol(M, 0);
        VI sum(N_GUESSES, 0);
        REP(i, M) REP(j, shapes[i].len) for (int g : ct2[sol[i]+shapes[i].points[j]]) sum[g]++;
        VI bsol = sol;

        VI nsum = sum;

        int steps = 0;
        int last_acc = 0;
        VI acc(2);
        double t0 = 10;
        double tn = 0.05;
        double t = t0;

        double start_sa = elapsed();
        while (true) {
            if ((steps & 255) == 0) {
                double time_passed = (elapsed() - start_sa) / (times * 0.1);
                if (time_passed > 1) break;
                t = t0 * pow(tn / t0, pow(time_passed, 0.9));
            }

            int type = rng.next_double() < 0.01;
            int m = rng.next(M);
            int m2;
            nsum = sum;
            int np = sol[m];
            int oldp = np;
            int oldp2;

            if (type == 0) {
                REP(i, shapes[m].len) for (int g : ct2[np+shapes[m].points[i]]) nsum[g]--;
                if (rng.next_double() < 0.5) {
                    np = rng.next(N - shapes[m].w + 1) * N + rng.next(N - shapes[m].h + 1);
                } else {
                    int range = 2;
                    int r = np / N;
                    int c = np % N;
                    int minr = max(0, r - range);
                    int maxr = min(N - shapes[m].w, r + range);
                    int minc = max(0, c - range);
                    int maxc = min(N - shapes[m].h, c + range);
                    np = rng.next(minr, maxr+1) * N + rng.next(minc, maxc+1);
                }
                REP(i, shapes[m].len) for (int g : ct2[np+shapes[m].points[i]]) nsum[g]++;
                sol[m] = np;
            } else {
                m2 = rng.next(M-1);
                m2 += m2 >= m;
                oldp2 = sol[m2];
                REP(i, shapes[m].len) for (int g : ct2[sol[m]+shapes[m].points[i]]) nsum[g]--;
                REP(i, shapes[m2].len) for (int g : ct2[sol[m2]+shapes[m2].points[i]]) nsum[g]--;

                int r1 = sol[m] / N;
                int c1 = sol[m] % N;
                int r2 = sol[m2] / N;
                int c2 = sol[m2] % N;
                const int range = 2;
                int minr1 = max(0, r2 - range);
                int maxr1 = min(N - shapes[m].w, r2 + range);
                int minc1 = max(0, c2 - range);
                int maxc1 = min(N - shapes[m].h, c2 + range);
                int minr2 = max(0, r1 - range);
                int maxr2 = min(N - shapes[m2].w, r1 + range);
                int minc2 = max(0, c1 - range);
                int maxc2 = min(N - shapes[m2].h, c1 + range);
                double av = 1e9;
                FOR(r, minr1, maxr1+1) FOR(c, minc1, maxc1+1) {
                    int np = r * N + c;
                    REP(i, shapes[m].len) for (int g : ct2[np+shapes[m].points[i]]) nsum[g]++;
                    double v = 0;
                    REP(i, t2g.S) v += pre_eval[nsum[i]][t2v[i]];
                    REP(i, shapes[m].len) for (int g : ct2[np+shapes[m].points[i]]) nsum[g]--;
                    if (v < av) {
                        av = v;
                        sol[m] = np;
                    }
                }
                REP(i, shapes[m].len) for (int g : ct2[sol[m]+shapes[m].points[i]]) nsum[g]++;

                av = 1e9;
                FOR(r, minr2, maxr2+1) FOR(c, minc2, maxc2+1) {
                    int np = r * N + c;
                    REP(i, shapes[m2].len) for (int g : ct2[np+shapes[m2].points[i]]) nsum[g]++;
                    double v = 0;
                    REP(i, t2g.S) v += pre_eval[nsum[i]][t2v[i]];
                    REP(i, shapes[m2].len) for (int g : ct2[np+shapes[m2].points[i]]) nsum[g]--;
                    if (v < av) {
                        av = v;
                        sol[m2] = np;
                    }
                }
                REP(i, shapes[m2].len) for (int g : ct2[sol[m2]+shapes[m2].points[i]]) nsum[g]++;
            }

            double av = 0;
            REP(i, t2g.S) av += pre_eval[nsum[i]][t2v[i]];

            if (av <= bv + 1e-9 || rng.next_double() < exp((bv - av) / t)) {
            // if (av <= bv + 1e-9 + rng.next_double() * t) {
                bv = av;
                sum = nsum;
                acc[type]++;
                if (bv < xv) {
                    xv = bv;
                    last_acc = steps;
                    bsol = sol;
                    // DB(steps, xv, elapsed());
                }
            } else {
                if (type == 0) {
                    sol[m] = oldp;
                } else {
                    sol[m] = oldp;
                    sol[m2] = oldp2;
                }
            }
            steps++;

        }

        double gtdv = 0;
        if (gtd.S) {
            REP(i, t2g.S) {
                int sum = 0;
                for (int g : t2g[i]) sum += gtd[g];
                gtdv += pre_eval[sum][t2v[i]];
            }
        }

        double lo = 1.0;
        REP(i, t2g.S) {
            double mean = sum[i] * (1 - eps) + (N_CELLS - sum[i]) * eps;
            lo = min(lo, prob_normal_range(mean, sigma, t2v[i]));
        }

        DB(steps, last_acc, acc, xv, gtdv, bv, lo, elapsed());

        VI cells(N*N, 0);
        REP(p, N*N) cells[p] = 0;
        REP(i, M) REP(j, shapes[i].len) cells[bsol[i] + shapes[i].points[j]]++;
        VI xrv; REP(p, N*N) if (cells[p] > 0) xrv.PB(p);
        cout << "a " << xrv.S;
        for (int p : xrv) cout << " " << p/N << " " << p%N;
        cout << endl;
        int answer; cin >> answer;
        if (answer || times == 7) {
            cerr << "[DATA] time = " << elapsed()*1000 << endl;
            exit(0);
        }
    }

    

    REP(p, N*N) guesses[p] = -1;

    int step = -1;

    stats_reset();
    stats_add(0);

    VVI voffsets(M);
    config_tiles_left = 0; REP(i, M) config_tiles_left += shapes[i].len;

    while (config_shapes.S < M) {
        VVI counts(M, VI(N*N, 0));
        VI total(M, 0);
        VI mdone(M, 0);

        REP(i, M) voffsets[i].clear();

        for (int shape : config_shapes) mdone[shape] = 1;

        VC<uint8_t> suremin(N*N, M);
        VC<uint8_t> suremax(N*N, 0);
        REP(p, N*N) REP(i, M+1) if (config_stats[p][i]) {
            suremin[p] = min(suremin[p], (uint8_t)i);
            suremax[p] = max(suremax[p], (uint8_t)i);
        }

        REP(i, M) if (!mdone[i]) {
            total[i] = 0;
            REP(p, N*N) counts[i][p] = 0;

            Shape &s = shapes[i];
            REP(r_off, N - s.w + 1) REP(c_off, N - s.h + 1) {
                int off = r_off * N + c_off;
                bool ok = true;
                REP(j, s.len) {
                    int p = s.points[j] + off;
                    if (guesses[p] >= 0 && guesses[p] - suremin[p] == 0) {
                        ok = false;
                        break;
                    }
                }
                if (!ok) continue;

                total[i]++;
                voffsets[i].PB(off);
                REP(j, s.len) {
                    int p = s.points[j] + off;
                    counts[i][p]++;
                }
            }
        }

        double bv = 0;
        int bp = -1;
        REP(p, N*N) {
            if (guesses[p] != -1) continue;
            double av = 0;

            av += rng.next_double() * 1e-6;

            VD dp(M+1, 0);
            REP(i, M+1) dp[i] = 1.0 * config_stats[p][i] / n_configs;
            REP(m, M) if (!mdone[m]) {
                double p = 1.0 * counts[m][p] / total[m];
                for (int i = M-1; i >= 0; i--) {
                    dp[i+1] += dp[i] * p;
                    dp[i] *= 1 - p;
                }
            }

            double ev = 0;
            REP(i, M+1) if (dp[i]) ev += dp[i] * -log(dp[i]);
            av += ev;
            

            const int dr[] = {0, 0, 1, -1};
            const int dc[] = {1, -1, 0, 0};
            int r = p / N;
            int c = p % N;
            REP(d, 4) {
                int nr = r + dr[d];
                int nc = c + dc[d];
                if (nr < 0 || nr >= N || nc < 0 || nc >= N) continue;
                int np = nr * N + nc;
                av += 1e-3 * (guesses[np] == -1);
            }
            int cdist = min(r, N - r - 1);
            int rdist = min(c, N - c - 1);
            av += 1e-4 * (cdist + rdist);

            if (av > bv) {
                bv = av;
                bp = p;
            }
        }


        if (bp == -1) break;

        if (step >= 0) {
            cout << "q 1 " << bp/N << " " << bp%N << endl;
            cin >> guesses[bp];
            guesses_list.PB(bp);

            double possibilities = n_configs;
            REP(i, M) if (!mdone[i]) possibilities *= total[i];
            DB(step, possibilities, bv, total, n_configs, config_shapes, bp/N, bp%N, elapsed());

            int max_possible = 0; REP(i, M) if (!mdone[i]) max_possible += counts[i][bp] > 0;
            config_update(bp, guesses[bp], true, max_possible);
        }
        step++;

        while (config_shapes.S < M) {
            VI used_shapes(M);
            for (int shape : config_shapes) used_shapes[shape] = 1;

            int best_shape = -1;
            double bv = 1e9;
            LL bes = 1LL << 60;
            REP(shape, M) if (!used_shapes[shape]) {
                
                double timer1_start = elapsed();
                int sim = 0;
                REP(i, 100) sim += config_expand_check(shape, voffsets[shape][rng.next(voffsets[shape].S)], rng.next(n_configs));
                timer1 += elapsed() - timer1_start;
                int expected_shapes = (LL)n_configs * total[shape] * sim / 100;

                double av = expected_shapes - 1e-6 * shapes[shape].len;
                if (av < bv) {
                    bv = av;
                    best_shape = shape;
                    bes = expected_shapes;
                }
            }

            if (bes > CONFIG_LIMIT) break;
            if (n_configs + bes > MAX_CONFIGS) break;

            double timer1_start = elapsed();
            int sim = 0;
            REP(i, 1000) sim += config_expand_check(best_shape, voffsets[best_shape][rng.next(voffsets[best_shape].S)], rng.next(n_configs));
            timer1 += elapsed() - timer1_start;
            bes = (LL)n_configs * total[best_shape] * sim / 1000;
            if (bes > CONFIG_LIMIT * 0.9 || n_configs + bes > MAX_CONFIGS * 0.9) break;

            int prev_n_configs = n_configs;
            config_expand(best_shape, voffsets[best_shape], true);
            DB(prev_n_configs, prev_n_configs*total[best_shape], bes, n_configs);
        }

    }

    while (n_configs > 4) {
        // bool early_exit = true;
        // REP(p, N*N) if (config_stats[p][0] > 0 && config_stats[p][0] < n_configs) early_exit = false;
        // if (early_exit) break;

        double bv = -1e9;
        int bp = -1;
        REP(p, N*N) if (guesses[p] == -1) {
            double av = 0;
            REP(i, M) if (config_stats[p][i]) {
                double v = 1.0 * config_stats[p][i] / n_configs;
                av += v * -log(v);
            }

            if (av > bv) {
                bv = av;
                bp = p;
            }
        }

        if (bp == -1) break;

        VI stats = VI(config_stats[bp], config_stats[bp]+M+1);

        cout << "q 1 " << bp/N << " " << bp%N << endl;
        cin >> guesses[bp];
        guesses_list.PB(bp);

        config_update(bp, guesses[bp], elapsed());

        DB(step, n_configs, bv, stats);
        step++;
    }

    cerr << "[DATA] time = " << elapsed()*1000 << endl;

    if (timer1) DB(timer1);
    if (timer2) DB(timer2);
    if (timer3) DB(timer3);

    REP(i, n_configs) {
        VI rv; REP(p,N*N) if (configs[i][p] > 0) rv.PB(p);
        cout << "a " << rv.S;
        for (int p : rv) cout << " " << p/N << " " << p%N;
        cout << endl;
        int answer;
        cin >> answer;
        if (answer) break;
    }

	return 0;
}
