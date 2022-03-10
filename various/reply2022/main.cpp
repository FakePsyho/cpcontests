// Author: Psyho
// Twitter: https://twitter.com/fakepsyho
 
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
#define DB(...) do {cout << DB_(ARGS_SIZE(__VA_ARGS__), __VA_ARGS__) endl;} while (0)
 
double get_time() {timeval tv; gettimeofday(&tv, NULL); return tv.tv_sec + tv.tv_usec * 1e-6;}
 
struct RNG {
    unsigned int MT[624];
    int index;
   
    RNG(int seed = 1) {
        init(seed);
    }
   
    void init(int seed = 1) {
        MT[0] = seed;
        FOR(i, 1, 624) MT[i] = (1812433253UL * (MT[i-1] ^ (MT[i-1] >> 30)) + i);
        index = 0;
    }
   
    void generate() {
        const unsigned int MULT[] = {0, 2567483615UL};
        REP(i, 227) {
            unsigned int y = (MT[i] & 0x8000000UL) + (MT[i+1] & 0x7FFFFFFFUL);
            MT[i] = MT[i+397] ^ (y >> 1);
            MT[i] ^= MULT[y&1];
        }
        FOR(i, 227, 623) {
            unsigned int y = (MT[i] & 0x8000000UL) + (MT[i+1] & 0x7FFFFFFFUL);
            MT[i] = MT[i-227] ^ (y >> 1);
            MT[i] ^= MULT[y&1];
        }
        unsigned int y = (MT[623] & 0x8000000UL) + (MT[0] & 0x7FFFFFFFUL);
        MT[623] = MT[623-227] ^ (y >> 1);
        MT[623] ^= MULT[y&1];
    }
   
    unsigned int rand() {
        if (index == 0) {
            generate();
        }
       
        unsigned int y = MT[index];
        y ^= y >> 11;
        y ^= y << 7  & 2636928640UL;
        y ^= y << 15 & 4022730752UL;
        y ^= y >> 18;
        index = index == 623 ? 0 : index + 1;
        return y;
    }
   
    INLINE int next() {
        return rand();
    }
   
    INLINE int next(int x) {
        return rand() % x;
    }
   
    INLINE int next(int a, int b) {
        return a + (rand() % (b - a));
    }
   
    INLINE double next_double() {
        return (rand() + 0.5) * (1.0 / 4294967296.0);
    }
   
    INLINE double next_double(double a, double b) {
        return a + next_double() * (b - a);
    }
};
 
static RNG rng;

double start_time = get_time();
int best_saved_score = 0;

const string RESULTS_DIR = "results";
string test_case = "";
string run_name = "";
string load_path = "";

double last_save;

int s_start;
int s_max;
int t_total;

struct Demon {
	int s_lost;
	int s_gain;
	int s_turn;
	VI pts;
};
VC<Demon> dm;


struct Test {
	
	void read(string fn) {
		ifstream fs(fn);
		int dn;
		fs >> s_start >> s_max >> t_total >> dn;
		REP(i, dn) {
			Demon d;
			int n;
			fs >>  d.s_lost >> d.s_turn >> d.s_gain >> n;
			d.pts = VI(n);
			REP(j, n) fs >> d.pts[j];
			dm.PB(d);
		}
	}
	
	void precalc() {
		for (Demon &d : dm) {
			FOR(i, 1, d.pts.S) d.pts[i] += d.pts[i-1];
			if (d.pts.S == 0) d.pts.PB(0);
		}
	}
	
	void analyze() {
		DB(test_case);
		DB(s_start, s_max, t_total, dm.S);
		int max_score = 0;
		VI scores;
		for (Demon &d : dm) scores.PB(d.pts[min((int)d.pts.S - 1, t_total)]);
		sort(ALL(scores));
		reverse(ALL(scores));
		REP(i, t_total) max_score += scores[i];
		DB(max_score);
		cout << endl;
	}
};

Test t;

struct State {
	VI sol;
	int bscore;
	
	VI s_inc;
	
	void init() {
		s_inc = VI(t_total);
	}
	
	int sim(VI order) {
		REP(i, t_total) s_inc[i] = 0;
		int score = 0;
		int t = -1;
		int s = s_start;
		int p = 0;
		while (true) {
			int d = order[p];
			t++;
			if (t >= t_total) break;
			s = min(s_max, s + s_inc[t]);
			if (s < dm[d].s_lost) continue;
			s -= dm[d].s_lost;
			if (t + dm[d].s_turn < t_total) s_inc[t + dm[d].s_turn] += dm[d].s_gain;
			score += dm[d].pts[min((int)dm[d].pts.S - 1, t_total - t - 1)];
			p++;
		}
		return score;
	}
	
	int debug_sim(VI order) {
		REP(i, t_total) s_inc[i] = 0;
		int score = 0;
		int t = -1;
		int s = s_start;
		int p = 0;
		
		int wasted_stamina = 0;
		while (true) {
			int d = order[p];
			t++;
			if (t >= t_total) break;
			wasted_stamina += max(0, s + s_inc[t] - s_max);
			s = min(s_max, s + s_inc[t]);
			if (s < dm[d].s_lost) continue;
			s -= dm[d].s_lost;
			if (t + dm[d].s_turn < t_total) s_inc[t + dm[d].s_turn] += dm[d].s_gain;
			score += dm[d].pts[min((int)dm[d].pts.S - 1, t_total - t - 1)];
			p++;
		}
		DB(t, p, score, wasted_stamina);
		return score;
	}
	
	void sa() {
		VI cur;
		if (sol.S == 0) {
			REP(i, dm.S) cur.PB(i);
			random_shuffle(ALL(cur));
			sol = cur;
		} else {
			cur = sol;
		}
		
		int bv = sim(cur);
		
		int xv = bv;
		int step = 0;
		
		int check = 0; //RUNNER//
		
		double t0 = 200; //RUNNER//
		double mt = 60; //RUNNER//
		while (true) {
			int a = rng.next(sol.S);
			int b = rng.next(sol.S);
			if (a >= t_total && b >= t_total) continue;
			step++;
			double t = max(1e-6, (mt - (get_time() - start_time)) / mt * t0);
			VI ncur = cur;
			int x = ncur[a];
			ncur.erase(ncur.begin() + a);
			ncur.insert(ncur.begin() + b, x);
			
			int av = sim(ncur);
			
			if (check && step % check == 0)
				debug_sim(cur);
			
			if (av >= bv || t0 && rng.next_double() < exp((av - bv) / t)) {
				cur = ncur;
				if (av > xv) {
					xv = av;
					sol = cur;
					update_results(*this, xv);
					DB(step, t, xv);
				}
				bv = av;
			} else {
				//pass
			}
		}
	}
	
	
	void show_stats() {
		//TODO: implement
	}
	
	int recalc_score() {
		init();
		bscore = sim(sol);
		return 0;
	}		
	
	void load(string fn) {
		ifstream fs(fn);
		sol = VI();
		while (true) {
			if (fs.eof()) break;
			int x = -1;
			fs >> x;
			if (x < 0) break;
			sol.PB(x);
		}
		// set<int> left;
		// REP(i, dm.S) left.insert(i);
		// REP(i, sol.S) left.erase(sol[i]);
		// for (int x : left) sol.PB(x);
		init();
		debug_sim(sol);
	}
	
	void save(string fn) {
		ofstream fs(fn);
		REP(i, sol.S) fs << sol[i] << endl;
	}
	
	void update_results(State &s, int score) {
		double cur_time = get_time();
		if (score > best_saved_score && cur_time > last_save + 1.0) {
			last_save = cur_time;
			best_saved_score = score;
			s.save(RESULTS_DIR + "/" + test_case + "-" + run_name + "-" + i2s(score) + ".out");
		}
	}
};


State s;

int main(int argc, char **argv) {
	if (argc <= 1) {
		cout << "Need to specify test_case" << endl;
		exit(1);
	}
	
	test_case = argv[1];
	
	FOR(i, 2, argc) {
		string cmd = argv[i];
		assert(cmd.S >= 3);
		assert(cmd[1] == '=');
		
		if (cmd[0] == 'l') {
			load_path = cmd.substr(2);
		} else if (cmd[0] == 'n') {
			run_name = cmd.substr(2);
		} else {
			cout << "Unrecognized cmd-line parameter: " << cmd[0] << endl;
			exit(1);
		}
	}
	
	t.read(test_case + ".txt");
	t.precalc();
	
	int analyze = 0; //RUNNER//
	if (analyze) {
		t.analyze();
		exit(0);
	}
	
	s.init();
	
	if (load_path.S) {
		s.load(load_path);
		// best_saved_score = s.recalc_score();
		s.show_stats();
	}
	
	s.sa();
	
	return 0;
}
