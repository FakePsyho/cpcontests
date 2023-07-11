// Author: Psyho
// Twitter: https://twitter.com/fakepsyho
 
// TEMPLATE

#pragma GCC optimize "Ofast,omit-frame-pointer,inline,unroll-all-loops"

#include <bits/stdc++.h>
#include <sys/time.h>
#include "json.hpp"
 
using namespace std;
using json = nlohmann::json;
 
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
#define VPDD        VC<PDD>
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
 
static RNG rng;

// GLOBAL PARAMETERS

//core
const double total_time = 0.0;
const int n_threads = 0;
const int no_save = 0;
const double target = 0;
const int exit_on_time = 0;
//stats & save
const int sa_stats_upd = 1000;
const int sa_acc_stats_period = 20;
const double sa_best_upd = 10.0;
const double sa_save_upd = 10.0;
//sa
const double sa_t0 = 1e7;
const double sa_tn = 1e3;
const double sa_tmin = 1;
const int sa_tlin = 0;
//problem-specific	
constexpr int test_id = 0;
constexpr int analyze = 0;
constexpr int convert = 0;
constexpr int OG_SCORING = 1;
constexpr double ST2P = 0.0;

// SOLUTION

double start_time = get_time();
atomic<double> time_passed;
atomic<double> sa_t;

int RULES = 0;

const double PI = acos(-1.0);

const LL COLL_PENALTY = 1e9;

const double M_RADIUS = 10.0;
const double EPS = 1e-9;

template<class A, class B> INLINE double dist2(pair<A, A> a, pair<B, B> b) {return (a.X-b.X)*(a.X-b.X)+(a.Y-b.Y)*(a.Y-b.Y);}
template<class A, class B> INLINE double dist(pair<A, A> a, pair<B, B> b) {return sqrt(dist2(a, b));}
INLINE LL scoref(LL x) {return OG_SCORING ? x : max(0LL, x);}

const int PARTITIONS_X = 100;
const int PARTITIONS_Y = 100;
const double MAX_DIST_PROP = 50;

template<class A, class B> INLINE double dist_rect(pair<A, A> r1, pair<A, A> r2, pair<B, B> p) {
    double d = 1e30;
    d = min(d, dist(MP(r1.X, r1.Y), MP(p.X, p.Y)));
    d = min(d, dist(MP(r1.X, r2.Y), MP(p.X, p.Y)));
    d = min(d, dist(MP(r2.X, r1.Y), MP(p.X, p.Y)));
    d = min(d, dist(MP(r2.X, r2.Y), MP(p.X, p.Y)));
    if (p.X >= r1.X && p.X <= r2.X) d = min(abs(p.Y-r1.Y), abs(p.Y-r2.Y));
    if (p.Y >= r1.Y && p.Y <= r2.Y) d = min(abs(p.X-r1.X), abs(p.X-r2.X));
    return d;
}

struct TestCase {
    PII room;
    pair<PII, PII> stage;
    VI musicians;
    VPII a_pos;
    VVD a_pref;
    VPII p_pos;
    VI p_rad;
    int m, n;

    int m_types;
    VVI m_list;
    VI a_list[PARTITIONS_X*PARTITIONS_Y];

    void load(int id) {
        fstream fs(string("tests/test") + i2s(id) + ".json", fstream::in);
        json j = json::parse(fs);
        room = PII(j["room_width"], j["room_height"]);
        stage = MP(PII(j["stage_bottom_left"][0], j["stage_bottom_left"][1]), PII(int(j["stage_bottom_left"][0])+int(j["stage_width"]), int(j["stage_bottom_left"][1])+int(j["stage_height"])));
        musicians = VI(ALL(j["musicians"]));
        for (auto a : j["attendees"]) {
            a_pos.PB(PII(a["x"], a["y"]));
            a_pref.PB(VD(ALL(a["tastes"])));
        }
        for (auto p : j["pillars"]) {
            p_pos.PB(PII(p["center"][0], p["center"][1]));
            p_rad.PB(p["radius"]);
        }

        n = a_pos.S;
        m = musicians.S;

        m_types = 0;
        for (auto m : musicians) m_types = max(m_types, m+1);

        m_list = VVI(m_types);
        REP(i, m) m_list[musicians[i]].PB(i);

        double min_dist = 1e30;
        REP(i, n) min_dist = min(min_dist, dist_rect(MP(stage.X.X + M_RADIUS, stage.X.Y + M_RADIUS), MP(stage.Y.X - M_RADIUS, stage.Y.Y - M_RADIUS), a_pos[i]));
        double max_dist = min_dist * MAX_DIST_PROP;
        DB(min_dist, max_dist);

        REP(px, PARTITIONS_X) REP(py, PARTITIONS_Y) {
            double x1 = stage.X.X + px * (stage.Y.X - stage.X.X) / PARTITIONS_X;
            double x2 = stage.X.X + (px+1) * (stage.Y.X - stage.X.X) / PARTITIONS_X;
            double y1 = stage.X.Y + py * (stage.Y.Y - stage.X.Y) / PARTITIONS_Y;
            double y2 = stage.X.Y + (py+1) * (stage.Y.Y - stage.X.Y) / PARTITIONS_Y;
            x1 = max(x1, stage.X.X + M_RADIUS);
            x2 = min(x2, stage.Y.X - M_RADIUS);
            y1 = max(y1, stage.X.Y + M_RADIUS);
            y2 = min(y2, stage.Y.Y - M_RADIUS);
            static VC<pair<double,int>> vp;
            vp.clear();
            REP(i, n) {
                double x = a_pos[i].X;
                double y = a_pos[i].Y;
                double d = dist_rect(MP(x1, y1), MP(x2, y2), a_pos[i]);
                vp.PB(MP(d, i));
            }
            sort(ALL(vp));
            REP(i, vp.S) {
                if (vp[i].X > max_dist) break;
                // if (vp[i].X > MAX_DIST_PROP * vp[0].X) break;
                a_list[px*PARTITIONS_Y+py].PB(vp[i].Y);
            }
        }

        double a_list_total = 0;
        REP(i, PARTITIONS_X*PARTITIONS_Y) a_list_total += a_list[i].S;
        a_list_total /= PARTITIONS_X*PARTITIONS_Y;

        DB(n, m, m_types, p_rad.S, a_list_total);


        RULES = p_rad.S > 0;
    }
};

bool intersect(PDD m1, PDD m2, PII aud, double radius=5.0) {
    double a = m2.X - m1.X;
    double b = m2.Y - m1.Y;
    double c = aud.X - m1.X;
    double d = aud.Y - m1.Y;

    double dot = a * c + b * d;
    double len_sq = c * c + d * d;
    if (dot < 0 || dot > len_sq) return false;

    double param = dot / len_sq;

    double dx = m2.X - (m1.X + param * c);
    double dy = m2.Y - (m1.Y + param * d);
    return dx * dx + dy * dy <= radius*radius;

}

TestCase t;

int get_partition(PDD &p) {
    int px = int((p.X - t.stage.X.X) * PARTITIONS_X / (t.stage.Y.X - t.stage.X.X));
    int py = int((p.Y - t.stage.X.Y) * PARTITIONS_Y / (t.stage.Y.Y - t.stage.X.Y));
    return px * PARTITIONS_Y + py;
}



struct SA_Move {
    int type;
	int m1;
    int m2;
    PDD npos;
    double score1;
    double score2;
};

struct SA_AccMove {
	SA_Move move;
	double score;
	int pos;
	
	SA_AccMove() { }
	
	SA_AccMove(int _pos, double _score, SA_Move &_move) : pos(_pos), score(_score), move(_move) { }
};

struct State {
    VPDD pos;
    VVI blocks;
    VC<LL> mscores;
    VI colls;
    VD closeness;

    double score = 0;
    int coll = 0;

    string prev_file_name = "";

    void init() {
        pos = VPDD(t.m);
        blocks = VVI(t.m, VI(t.n, 0));
        mscores = VC<LL>(t.m);
        colls = VI(t.m, 0);
        closeness = VD(t.m, 1.0);
    }

    void random_init() {
        init();

        REP(i, t.m) pos[i] = PDD(rng.next_double(t.stage.X.X + M_RADIUS, t.stage.Y.X - M_RADIUS), rng.next_double(t.stage.X.Y + M_RADIUS, t.stage.Y.Y - M_RADIUS));

        update_full();
    }

    void update_full() {
        coll = 0;
        REP(i, t.m) {
            colls[i] = 0;
            REP(j, t.m) if (i != j) colls[i] += dist2(pos[i], pos[j]) < M_RADIUS * M_RADIUS + EPS;
            coll += colls[i];
        }
        coll /= 2;

        REP(i, t.m) {
            mscores[i] = 0;
            REP(j, t.n) {
                blocks[i][j] = 0;
                REP(k, t.m) if (k != i) blocks[i][j] += intersect(pos[i], pos[k], t.a_pos[j]);
                REP(k, t.p_pos.S) blocks[i][j] += intersect(pos[i], t.p_pos[k], t.a_pos[j], t.p_rad[k]);
                if (blocks[i][j] == 0) mscores[i] += (int)ceil(1e6 * t.a_pref[j][t.musicians[i]] / dist2(pos[i], t.a_pos[j]));
            }
        }

        if (RULES == 1) {
            REP(i, t.m) {
                closeness[i] = 1.0;
                for (int j : t.m_list[t.musicians[i]]) if (i != j) {
                    closeness[i] += 1.0 / dist(pos[i], pos[j]);
                }
            }
        }

        score = 0;
        REP(i, t.m) score += scoref(mscores[i]) * closeness[i];
    }

    void update(int m, PDD np) {
        REP(i, t.m) if (m != i) if (dist2(pos[m], pos[i]) < M_RADIUS * M_RADIUS + EPS) {
            colls[m]--;
            colls[i]--;
            coll--;
        }

        assert(colls[m] == 0);

        score -= scoref(mscores[m]) * closeness[m];
        mscores[m] = 0;

        REP(k, t.m) if (k != m) for (int j : t.a_list[get_partition(pos[k])]) if (intersect(pos[k], pos[m], t.a_pos[j])) {
            blocks[k][j]--;
            if (blocks[k][j] == 0) {
                int s = (int)ceil(1e6 * t.a_pref[j][t.musicians[k]] / dist2(pos[k], t.a_pos[j]));
                score -= scoref(mscores[k]) * closeness[k];
                mscores[k] += s;
                score += scoref(mscores[k]) * closeness[k];
            }
        }

        if (RULES == 1) {
            for (int i : t.m_list[t.musicians[m]]) if (i != m) {
                double v = 1.0 / dist(pos[i], pos[m]);
                closeness[i] -= v;
                score -= scoref(mscores[i]) * v;
            }
        }

        pos[m] = np;

        REP(i, t.m) if (m != i) if (dist2(np, pos[i]) < M_RADIUS * M_RADIUS + EPS) {
            colls[m]++;
            colls[i]++;
            coll++;
        }

        if (RULES == 1) {
            closeness[m] = 1.0;
            for (int i : t.m_list[t.musicians[m]]) if (i != m) {
                double v = 1.0 / dist(pos[i], pos[m]);
                closeness[i] += v;
                closeness[m] += v;
                score += scoref(mscores[i]) * v;
            }
        }

        for (int j : t.a_list[get_partition(np)]) if (intersect(np, pos[m], t.a_pos[j])) {
            blocks[m][j] = 0;
            REP(k, t.m) if (k != m) blocks[m][j] += intersect(np, pos[k], t.a_pos[j]);
            REP(k, t.p_pos.S) blocks[m][j] += intersect(np, t.p_pos[k], t.a_pos[j], t.p_rad[k]);
            if (blocks[m][j] == 0) mscores[m] += (int)ceil(1e6 * t.a_pref[j][t.musicians[m]] / dist2(np, t.a_pos[j]));
        }
        score += scoref(mscores[m]) * closeness[m];

        REP(k, t.m) if (k != m) for (int j : t.a_list[get_partition(pos[k])]) if (intersect(pos[k], pos[m], t.a_pos[j])) {
            blocks[k][j]++;
            if (blocks[k][j] == 1) {
                int s = (int)ceil(1e6 * t.a_pref[j][t.musicians[k]] / dist2(pos[k], t.a_pos[j]));
                score -= scoref(mscores[k]) * closeness[k];
                mscores[k] -= s;
                score += scoref(mscores[k]) * closeness[k];
            }
        }
    }

    void update_swap(int m1, int m2, double score1 = -1e9, double score2 = -1e9) {
        assert(t.musicians[m1] != t.musicians[m2]);

        if (RULES == 1) {
            for (int i : t.m_list[t.musicians[m1]]) if (i != m1) {
                double v = 1.0 / dist(pos[i], pos[m1]);
                closeness[i] -= v;
                score -= scoref(mscores[i]) * v;
            }
            for (int i : t.m_list[t.musicians[m2]]) if (i != m2) {
                double v = 1.0 / dist(pos[i], pos[m2]);
                closeness[i] -= v;
                score -= scoref(mscores[i]) * v;
            }
        }

        score -= scoref(mscores[m1]) * closeness[m1];
        score -= scoref(mscores[m2]) * closeness[m2];

        swap(pos[m1], pos[m2]);
        swap(colls[m1], colls[m2]);

        mscores[m1] = 0;
        mscores[m2] = 0;

        if (RULES == 1) {
            closeness[m1] = 1.0;
            closeness[m2] = 1.0;
            for (int i : t.m_list[t.musicians[m1]]) if (i != m1) {
                double v = 1.0 / dist(pos[i], pos[m1]);
                closeness[i] += v;
                closeness[m1] += v;
                score += scoref(mscores[i]) * v;
            }
            for (int i : t.m_list[t.musicians[m2]]) if (i != m2) {
                double v = 1.0 / dist(pos[i], pos[m2]);
                closeness[i] += v;
                closeness[m2] += v;
                score += scoref(mscores[i]) * v;
            }
        }

        if (score1 > -1e9) {
            REP(i, t.n) swap(blocks[m1][i], blocks[m2][i]);
            mscores[m1] = score1;
            mscores[m2] = score2;
        } else {

            REP(i, t.n) {
                if (blocks[m1][i] == 0) mscores[m2] += (int)ceil(1e6 * t.a_pref[i][t.musicians[m2]] / dist2(pos[m2], t.a_pos[i]));
                if (blocks[m2][i] == 0) mscores[m1] += (int)ceil(1e6 * t.a_pref[i][t.musicians[m1]] / dist2(pos[m1], t.a_pos[i]));
                swap(blocks[m1][i], blocks[m2][i]);
            }
        }
        score += scoref(mscores[m1]) * closeness[m1];; 
        score += scoref(mscores[m2]) * closeness[m2];;

    }

    bool no_collision(int m, PDD np) {
        REP(i, t.m) if (m != i) if (dist2(np, pos[i]) < M_RADIUS * M_RADIUS + EPS) return false;
        return true;
    }

    bool outside(PDD np) {
        return np.X - M_RADIUS < t.stage.X.X || np.X + M_RADIUS > t.stage.Y.X || np.Y - M_RADIUS < t.stage.X.Y || np.Y + M_RADIUS > t.stage.Y.Y;
    }

    SA_Move gen_move(int type_seed = -1) {
        int type2 = max(1, min(t.n, t.m) / 5);

        int type = type_seed == -1 ? rng.next(1 + type2) > 0 : (type_seed % (1 + type2)) > 0;
        if (coll > 0) type = 0;

        int m1 = rng.next(t.m);
        int m2 = rng.next(t.m);

        while (coll && colls[m1] == 0) m1 = rng.next(t.m);

        PDD old_p1 = pos[m1];
        PDD old_p2 = pos[m2];
        PDD new_p;

        if (type == 0) {
            int subtype = rng.next(7);
            
            if (subtype == 0) {
                new_p = PDD(rng.next_double(t.stage.X.X + M_RADIUS, t.stage.Y.X - M_RADIUS), rng.next_double(t.stage.X.Y + M_RADIUS, t.stage.Y.Y - M_RADIUS));
            } else if (subtype == 1) {
                double min_x = t.stage.X.X + M_RADIUS;
                double max_x = t.stage.Y.X - M_RADIUS;
                double min_y = t.stage.X.Y + M_RADIUS;
                double max_y = t.stage.Y.Y - M_RADIUS;
                REP(tries, 100) {
                    min_x = max(min_x, old_p1.X - 25);
                    max_x = min(max_x, old_p1.X + 25);
                    min_y = max(min_y, old_p1.Y - 25);
                    max_y = min(max_y, old_p1.Y + 25);
                    new_p = PDD(rng.next_double(min_x, max_x), rng.next_double(min_y, max_y));
                    if (no_collision(m1, new_p)) break;
                }
            } else if (subtype == 2) {
                REP(tries, 100) {
                    new_p = old_p1;
                    double angle = rng.next_double(0, 2 * PI);
                    double d = pow(25.0, rng.next_double()) - 1.0;
                    new_p.X += d * cos(angle);
                    new_p.Y += d * sin(angle);
                    new_p.X = max(new_p.X, t.stage.X.X + M_RADIUS);
                    new_p.X = min(new_p.X, t.stage.Y.X - M_RADIUS);
                    new_p.Y = max(new_p.Y, t.stage.X.Y + M_RADIUS);
                    new_p.Y = min(new_p.Y, t.stage.Y.Y - M_RADIUS);
                    if (no_collision(m1, new_p)) break;
                }
            } else if (subtype == 3) {
                REP(tries, 100) {
                    int side = rng.next(2);
                    if (side == 0) {
                        new_p = PDD(rng.next_double(t.stage.X.X + M_RADIUS, t.stage.Y.X - M_RADIUS), rng.next(2) ? t.stage.X.Y + M_RADIUS : t.stage.Y.Y - M_RADIUS);
                    } else {
                        new_p = PDD(rng.next(2) ? t.stage.X.X + M_RADIUS : t.stage.Y.X - M_RADIUS, rng.next_double(t.stage.X.Y + M_RADIUS, t.stage.Y.Y - M_RADIUS));
                    }
                    if (no_collision(m1, new_p)) break;
                }
            } else if (subtype == 4) {
                REP(tries, 100) {
                    int side = rng.next(2);
                    if (side == 0) {
                        new_p = PDD(rng.next_double(t.stage.X.X + M_RADIUS, t.stage.Y.X - M_RADIUS), rng.next(2) ? t.stage.X.Y + M_RADIUS * 2 : t.stage.Y.Y - M_RADIUS * 2);
                    } else {
                        new_p = PDD(rng.next(2) ? t.stage.X.X + M_RADIUS : t.stage.Y.X - M_RADIUS, rng.next_double(t.stage.X.Y + M_RADIUS * 2, t.stage.Y.Y - M_RADIUS * 2));
                    }
                    if (no_collision(m1, new_p)) break;
                }
                if (outside(new_p)) return gen_move();
            } else if (subtype == 5) {
                REP(tries, 100) {
                    m2 = rng.next(t.m);
                    double angle = rng.next_double(0, 2 * PI);
                    new_p.X = pos[m2].X + (M_RADIUS + 2 * EPS) * cos(angle);
                    new_p.Y = pos[m2].Y + (M_RADIUS + 2 * EPS) * sin(angle);
                    if (!outside(new_p) && no_collision(m1, new_p)) break;
                }
            } else if (subtype == 6) {
                double angle = rng.next_double(0, 2 * PI);
                double min_t = 1e9;

                double Dx = cos(angle);
                double Dy = sin(angle);
                REP(i, t.m) if (i != m1) {
                    double t = Dx * (pos[i].X - pos[m1].X) + Dy * (pos[i].Y - pos[m1].Y);
                    if (t < 0) continue;
                    double Ex = pos[m1].X + t * Dx;
                    double Ey = pos[m1].Y + t * Dy;
                    double LEC2 = dist2(pos[i], PDD(Ex, Ey));
                    if (LEC2 > M_RADIUS * M_RADIUS + EPS) continue;
                    double dt = sqrt(M_RADIUS * M_RADIUS - LEC2);
                    if (t - dt > 0) min_t = min(min_t, t - dt);
                    if (t + dt > 0) min_t = min(min_t, t + dt);
                    // if (coll == 0) assert(t - dt >= 0);
                }
                if (Dx > 0) min_t = min(min_t, (t.stage.Y.X - pos[m1].X - M_RADIUS) / Dx);
                if (Dx < 0) min_t = min(min_t, (t.stage.X.X - pos[m1].X + M_RADIUS) / Dx);
                if (Dy > 0) min_t = min(min_t, (t.stage.Y.Y - pos[m1].Y - M_RADIUS) / Dy);
                if (Dy < 0) min_t = min(min_t, (t.stage.X.Y - pos[m1].Y + M_RADIUS) / Dy);
                // assert(min_t >= 0);
                min_t = max(0.0, min_t - EPS - 1e-6);
                new_p.X = pos[m1].X + Dx * min_t;
                new_p.Y = pos[m1].Y + Dy * min_t;
                if (min_t == 0 || outside(new_p) || !no_collision(m1, new_p)) return gen_move();
            }
        } else {
            if (t.musicians[m1] == t.musicians[m2]) return gen_move();
        }

        new_p.X = max(new_p.X, t.stage.X.X + M_RADIUS);
        new_p.X = min(new_p.X, t.stage.Y.X - M_RADIUS);
        new_p.Y = max(new_p.Y, t.stage.X.Y + M_RADIUS);
        new_p.Y = min(new_p.Y, t.stage.Y.Y - M_RADIUS);
        if (!no_collision(m1, new_p)) return gen_move();

        SA_Move move;
        move.type = type;
        move.m1 = m1;
        move.m2 = m2;
        move.npos = new_p;
        return move;
    }

    double eval_state() {
        update_full();
        return score * 10 - coll * COLL_PENALTY;
    }

    double eval_move(SA_Move &move) {
        double rv;
        if (move.type == 0) {
            PDD opos = pos[move.m1];
            update(move.m1, move.npos);
            rv = score * 10 - coll * COLL_PENALTY;
            update(move.m1, opos);
        } else {
            double oscore1 = mscores[move.m1];
            double oscore2 = mscores[move.m2];
            update_swap(move.m1, move.m2);
            rv = score * 10 - coll * COLL_PENALTY;
            move.score1 = mscores[move.m1];
            move.score2 = mscores[move.m2];
            update_swap(move.m1, move.m2, oscore1, oscore2);
        }
        return rv;
    }

    void apply_move(SA_Move move) {
        if (move.type == 0) {
            update(move.m1, move.npos);
        } else {
            update_swap(move.m1, move.m2, move.score1, move.score2);
        }
    }

    void load(string file_name) {
        fstream fs(file_name, fstream::in);
        json j = json::parse(fs);

        REP(i, t.m) pos[i] = PDD(j["placements"][i]["x"], j["placements"][i]["y"]);
    }

    void save(string file_name) {
        prev_file_name = file_name;

        assert(t.musicians.S == pos.S);

        json j;
        REP(i, t.m) j["placements"][i] = {{"x", pos[i].X}, {"y", pos[i].Y}};
        REP(i, t.m) j["volumes"][i] = mscores[i] >= 0 ? 10.0 : 0.0;

        fstream fs(file_name, fstream::out);
        fs << j << endl;
    }

    void analyze() {
        auto vmus = mscores;
        sort(ALL(vmus));
        DB(vmus);

        VI vaud = VI(t.n);
        REP(i, t.n) REP(j, t.m) if (blocks[j][i] == 0) vaud[i] += (int)ceil(1e6 * t.a_pref[i][t.musicians[j]] / dist2(pos[j], t.a_pos[i]));
        sort(ALL(vaud));
        DB(vaud);
    }
};

struct Solver {
	atomic<int> workers_ready;
	
	mutex stats_lock;
	mutex m_lock;
	condition_variable m_cond;
	
	LL m_step;
	LL m_acc;
	State m_state;
	double best_score;
	double cur_score;
	
	VC<LL> acc_period = VC<LL>(sa_acc_stats_period);
	VC<LL> step_period = VC<LL>(sa_acc_stats_period);
	LL sum_acc = 0;
	LL sum_step = 1;
	int period_pos = 0;
	
	const int w_moves_size = 1<<16;
	const int m_moves_size = 1<<16;
	VC<SA_AccMove> m_moves;
	VC<SA_AccMove> w_moves;
	atomic<int> w_moves_pos;
	atomic<int> m_moves_pos;
	VC<mutex> m_moves_locks;
	
	bool need_save = false;
	
	void update_results(State &s, double cur_time, LL score) {
		static double last_save = 0;
		static LL best_saved_score = 0;
		
		if (no_save) return;

		if (score > best_saved_score && cur_time > last_save + sa_save_upd) {

			last_save = cur_time;
			best_saved_score = score;

            string prev_file_name = s.prev_file_name;
            s.save(string("test") + i2s(test_id) + "_" + i2s((LL)(s.score * 10)) + ".txt");
            if (prev_file_name.S) remove(prev_file_name.c_str());
			need_save = false;
		}
	}
	
	
	//SA Stuff
	void mt_sa_worker(int id) {
		State state(m_state);
		int pos = 1;
		double bv = best_score;
		workers_ready++;
		m_cond.notify_one();
		
		rng = RNG(id * 1337);
		
		int w_step = 0;
		int w_acc = 0;
		
		while (true) {
			int w_pos = w_moves_pos;
			while (pos < w_pos) {
				SA_AccMove &acc_move = w_moves[pos & (w_moves_size - 1)];
				state.apply_move(acc_move.move);
				bv = acc_move.score;
				pos++;
				assert(pos == acc_move.pos);
			}
			
			SA_Move move = state.gen_move(w_moves_pos);
			
			double av = state.eval_move(move);
			// double av = bv + state.eval_move(move);
			
			w_step++;
			if (av >= bv || total_time && sa_t0 && rng.next_double() < exp((av - bv) / sa_t)) {
				w_acc++;
				SA_AccMove acc_move = SA_AccMove(pos, av, move);
				int cur_move_pos = m_moves_pos++;
				cur_move_pos &= m_moves_size - 1;
				{
					unique_lock<mutex> lock(m_moves_locks[cur_move_pos]);
					m_moves[cur_move_pos] = acc_move;
				} 
				m_cond.notify_one();
			}
			
			if (w_step >= sa_stats_upd) {
				int new_period_pos = ((LL)get_time()) % sa_acc_stats_period;
				stats_lock.lock();
				if (sa_acc_stats_period) {
					while (period_pos != new_period_pos) {
						period_pos = (period_pos + 1) % sa_acc_stats_period;
						sum_step -= step_period[period_pos];
						sum_acc -= acc_period[period_pos];
						step_period[period_pos] = 0;
						acc_period[period_pos] = 0;
					}
					
					step_period[period_pos] += w_step;
					acc_period[period_pos] += w_acc;
					sum_step += w_step;
					sum_acc += w_acc;
				}
				m_step += w_step;
				m_acc += w_acc;
				
				stats_lock.unlock();
				w_step = 0;
				w_acc = 0;
			}
		}
	}
	
	void init_sa() {
		m_step = 0;
		m_acc = 0;
		best_score = m_state.eval_state();
		cur_score = best_score;
	}
	
	void mt_sa() {
		init_sa();
		int m_pos = 1;
		w_moves_pos = 0;
		m_moves_pos = 0;
		
		workers_ready = 0;
		time_passed = 0.0;
		sa_t = total_time ? max(sa_tmin, sa_t0) : 0;
		m_moves_locks = VC<mutex>(m_moves_size);
		w_moves = VC<SA_AccMove>(w_moves_size);
		m_moves = VC<SA_AccMove>(m_moves_size);
		
		VC<thread> workers;
		REP(i, n_threads) workers.PB(thread(&Solver::mt_sa_worker, this, i));
		
		{
			unique_lock<mutex> lock(m_lock);
			m_cond.wait(lock, [&]() { return workers_ready == n_threads; } );
		}
		
		SA_AccMove acc_move;
		double cur_time = get_time();
		double last_update = cur_time;
		bool new_update = false;
		
		int start_move_pos = 0;
		
		while (true) {
			cur_time = get_time();
			if (total_time) {
				time_passed = (cur_time - start_time) / total_time;
				sa_t = sa_tlin ? max(sa_tmin, sa_t0 + time_passed * (sa_tn - sa_t0)) : sa_t0 * pow(sa_tn / sa_t0, time_passed);
                if (exit_on_time && time_passed >= 1.0)
                    exit(0);
			}
			
			mutex x_lock;
			bool correct_move = false;
			{
				unique_lock<mutex> lock(x_lock);
				m_cond.wait_for(lock, chrono::milliseconds(100), [&]() { return start_move_pos != m_moves_pos; });
				int cur_move_pos = m_moves_pos;
				while (start_move_pos != cur_move_pos) {
					int i = start_move_pos & (m_moves_size - 1);
					unique_lock<mutex> lock(m_moves_locks[i]);
					SA_AccMove am = m_moves[i];
					if (am.pos == m_pos && (!correct_move || am.score > acc_move.score)) {
						correct_move = true;
						acc_move = am;
					}
					start_move_pos++;
				}
				start_move_pos = m_moves_pos;
			}
			if (!correct_move) continue;
			
			acc_move.pos++;
			w_moves[m_pos & (w_moves_size - 1)] = acc_move;
			m_pos++;
			
			w_moves_pos = m_pos;
			
			m_state.apply_move(acc_move.move);
			cur_score = acc_move.score;
			
			new_update = true;
			if (cur_score > best_score) {
				best_score = cur_score;
				new_update = true;
				need_save = true;
				update_results(m_state, cur_time, best_score);
			}
			
			if (new_update && sa_best_upd && cur_time - last_update > sa_best_upd) {
				cerr << "Step: " << m_step << " Acc: " << m_acc << " Acc@" << sa_acc_stats_period << "s: " << 100.0 * sum_acc / sum_step << "% Pos: " << m_pos << " Score: " << (LL)best_score << " : " << (LL)cur_score << " T: " << sa_t << " Time: " << cur_time - start_time << endl;
				new_update = false;
				last_update = cur_time;
				if (target && best_score >= target)
					exit(0);
			}
		}
	}
	
	void sa() {
		init_sa();
		
		double cur_time = get_time();
		double last_update = cur_time;
		bool new_update = false;
		
		while (true) {
			m_step++;
			
			cur_time = get_time();
			if (total_time) {
				time_passed = (cur_time - start_time) / total_time;
				sa_t = sa_tlin ? max(sa_tmin, sa_t0 + time_passed * (sa_tn - sa_t0)) : sa_t0 * pow(sa_tn / sa_t0, time_passed);
                if (exit_on_time && time_passed >= 1.0)
                    exit(0);
			}
			
			SA_Move move = m_state.gen_move();
			double av = m_state.eval_move(move);
			// double av = cur_score + m_state.eval_move(move);


			
			if (sa_acc_stats_period) {
				int new_period_pos = ((LL)cur_time) % sa_acc_stats_period;
				while (period_pos != new_period_pos) {
					period_pos = (period_pos + 1) % sa_acc_stats_period;
					sum_step -= step_period[period_pos];
					sum_acc -= acc_period[period_pos];
					step_period[period_pos] = 0;
					acc_period[period_pos] = 0;
				}
				
				step_period[period_pos]++;
				sum_step++;
			}
			
			if (av >= cur_score || total_time && sa_t0 && rng.next_double() < exp((av - cur_score) / sa_t)) {
				cur_score = av;
				m_acc++;
				if (sa_acc_stats_period) {
					acc_period[period_pos]++;
					sum_acc++;
				}
				m_state.apply_move(move);
				if (cur_score > best_score) {
					best_score = cur_score;
					new_update = true;
					need_save = true;
					update_results(m_state, cur_time, best_score);
				}
			} else {
				
			}
			
			if (new_update && sa_best_upd && cur_time - last_update > sa_best_upd) {
				cerr << "Step: " << m_step << " Acc: " << m_acc << " Acc@" << sa_acc_stats_period << "s: " << 100.0 * sum_acc / sum_step << " Score: " << (LL)best_score << " T: " << sa_t << " Time: " << cur_time - start_time << endl;
				new_update = false;
				last_update = cur_time;
				if (target && best_score >= target)
					exit(0);
			}
			
			
		}
	}
};

int main(int argc, char **argv) {
    string res_file = "";

    assert(test_id > 0);

    t.load(test_id);

	Solver solver;
    solver.m_state.random_init();
    if (res_file.S) {
        solver.m_state.load(res_file);
        solver.m_state.update_full();
        if (analyze) {
            solver.m_state.analyze();
            return 0;
        }
        if (convert) {
            solver.m_state.save(string("test") + i2s(test_id) + "_" + i2s((LL)(solver.m_state.score * 10)) + ".txt");
            return 0;
        }
    }

    DB(solver.m_state.score, solver.m_state.coll, (LL)solver.m_state.eval_state());

	if (n_threads) {
		solver.mt_sa();
	} else {
		solver.sa();
	}

	return 0;
}
