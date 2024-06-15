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
};
 
static RNG rng;

// SOLUTION

INLINE int dst2(int x1, int y1, int x2, int y2) {int dx = x1 - x2; int dy = y1 - y2; return dx*dx + dy*dy;}
INLINE double dst(int x1, int y1, int x2, int y2) {return sqrt(dst2(x1, y1, x2, y2));}
INLINE int in_range(int x1, int y1, int x2, int y2, int r) {return dst2(x1, y1, x2, y2) <= r*r;}

int test_id;

int map_w;
int map_h;
int start_x;
int start_y;
int base_speed;
int base_power;
int base_range;
int level_speed_coeff;
int level_power_coeff;
int level_range_coeff;
int n_turns;

VVVI attacking_monsters;

struct Monster {
    int x;
    int y;
    int hp;
    int gold;
    int exp;
    int attack;
    int range;
};
VC<Monster> monsters;

int test_type;

VI levels;

VVI close_monsters;

void load_data(int id) {
    fstream fs(string("input/0") + i2s(id/10) + i2s(id%10) + ".json", fstream::in);
    json j = json::parse(fs);
    n_turns = j["num_turns"];
    start_x = j["start_x"];
    start_y = j["start_y"];
    map_w = j["width"];
    map_h = j["height"];
    base_speed = j["hero"]["base_speed"];
    base_power = j["hero"]["base_power"];
    base_range = j["hero"]["base_range"];
    level_speed_coeff = j["hero"]["level_speed_coeff"];
    level_power_coeff = j["hero"]["level_power_coeff"];
    level_range_coeff = j["hero"]["level_range_coeff"];

    if (id <= 25) {
        for (auto m : j["monsters"]) monsters.PB({m["x"], m["y"], m["hp"], m["gold"], m["exp"], 0, 0});
    } else {
        for (auto m : j["monsters"]) monsters.PB({m["x"], m["y"], m["hp"], m["gold"], m["exp"], m["attack"], m["range"]});
    }
}

bool good_monster(Monster &m) {return m.hp < 1e5 && (m.exp > 1 || m.gold > 1);}

struct Move {
    int monster;
    int mx;
    int my;
};
VC<Move> sol;

VI cur_order;
int bv;
VI best_order;
int xv;

void save_solution(int id, string fn="") {
    if (fn.S == 0) fn = string("output/0") + i2s(id/10) + i2s(id%10) + ".json";
    fstream fs(fn, fstream::out);
    json j;
    REP(i, sol.S) {
        auto& move = sol[i];
        if (sol[i].monster == -1) 
            j["moves"][i] = {{"type", "move"}, {"target_x", move.mx}, {"target_y", move.my}};
        else 
            j["moves"][i] = {{"type", "attack"}, {"target_id", move.monster}};
    }
    fs << j << endl;
}

const int MAX_DEPTH = 100;
const int MAX_FATIGUE = 1e6;
const int AVOID = 0;
const int XPATHFIND = 0;
const int XSTOP = 0;

int last_order;
int last_turns;
LL last_fatigue;
int max_depth = 10;
template<int TYPE> int simx(VI &order, bool generate = false) {
    if (generate) sol.clear();

    last_turns = n_turns;

    int x = start_x;
    int y = start_y;
    int cur_turn = 0;
    int level = 0;
    int exp = 0;
    int gold = 0;
    LL fatigue = 0;
    int cur_speed = base_speed;
    int cur_power = base_power;
    int cur_range = base_range;
    static VI monster_alive(monsters.S);

    int xpathfind_used = 0;
    
    if (TYPE) REP(i, monsters.S) monster_alive[i] = 1;

    if (test_id == 36 || test_id == 37) {
        while (y + cur_speed + cur_speed < map_h) {
            if (generate) sol.PB({-1, x, y + cur_speed});
            y += cur_speed;
            cur_turn++;
        }
    }

    REP(i, order.S) {
        last_order = i;
        last_fatigue = fatigue;
        auto &m = monsters[order[i]];
        while (!in_range(x, y, m.x, m.y, cur_range)) {
            if (cur_turn >= n_turns) return gold;

            double dx = m.x - x;
            double dy = m.y - y;
            double d = sqrt(dx*dx + dy*dy);
            dx /= d;
            dy /= d;
            double speed = min(cur_speed * 1.0, dst(x, y, m.x, m.y));
            if (XSTOP && m.range < cur_range) speed = min(speed, dst(x, y, m.x, m.y) - m.range - 1);
            dx *= speed;
            dy *= speed;
            int best_x = -1;
            int best_y = -1;
            FOR(tx, -1, 2) FOR(ty, -1, 2) {
                int nx = x + (int)dx + tx;
                int ny = y + (int)dy + ty;
                if (nx < 0 || nx > map_w || ny < 0 || ny > map_h) continue;
                if (!in_range(nx, ny, x, y, cur_speed)) continue;
                if (XSTOP && m.range < cur_range && in_range(nx, ny, m.x, m.y, m.range)) continue;
                if (best_x == -1 || dst2(nx, ny, m.x, m.y) < dst2(best_x, best_y, m.x, m.y)) {
                    best_x = nx;
                    best_y = ny;
                }
            }

            if (xpathfind_used < XPATHFIND && m.range < cur_range) {
                int aggro_monsters = 0;
                for (int m : attacking_monsters[best_x][best_y]) aggro_monsters += monster_alive[m] * monsters[m].attack;
                if (aggro_monsters) {
                    xpathfind_used++;
                    int orig_d2 = dst2(x, y, m.x, m.y);
                    int d2 = dst2(best_x, best_y, m.x, m.y);
                    LL bv = aggro_monsters * (LL)1e9 + d2;  
                    FOR(tx, -cur_speed, cur_speed + 1) FOR(ty, -cur_speed, cur_speed + 1) {
                        int nx = x + tx;
                        int ny = y + ty;
                        if (nx < 0 || nx > map_w || ny < 0 || ny > map_h) continue;
                        if (!in_range(nx, ny, x, y, cur_speed)) continue;
                        if (dst2(nx, ny, m.x, m.y) < orig_d2) continue;
                        int move_aggro_monsters = 0;
                        for (int m : attacking_monsters[nx][ny]) move_aggro_monsters += monster_alive[m] * monsters[m].attack;
                        LL av = move_aggro_monsters * (LL)1e9 + dst2(nx, ny, m.x, m.y);
                        if (av < bv) {
                            bv = av;
                            best_x = nx;
                            best_y = ny;
                        }
                    }
                }
            }

            assert(best_x != -1);
            int depth = 0;
            const int ddx[] = {0, 1, 0, -1, -1, -1, 1, 1};
            const int ddy[] = {1, 0, -1, 0, -1, 1, -1, 1};
            if (MAX_DEPTH) {
                while (true) {
                    if (!in_range(best_x, best_y, monsters[order[i+depth]].x, monsters[order[i+depth]].y, cur_range)) break;
                    depth++;
                    if (i+depth >= order.S) break;
                    int aggro_monsters = 0;
                    if (TYPE && AVOID == 1) REP(j, depth) if (in_range(best_x, best_y, monsters[order[i+j]].x, monsters[order[i+j]].y, monsters[order[i+j]].range)) aggro_monsters++;
                    if (TYPE && AVOID == 2) for (int m : attacking_monsters[best_x][best_y]) aggro_monsters += monster_alive[m];
                    while (true) {
                        bool update = false;
                        REP(d, 8) {
                            int nx = best_x + ddx[d];
                            int ny = best_y + ddy[d];
                            if (nx < 0 || nx > map_w || ny < 0 || ny > map_h) continue;
                            if (!in_range(nx, ny, x, y, cur_speed)) continue;
                            REP(j, depth) if (!in_range(nx, ny, monsters[order[i+j]].x, monsters[order[i+j]].y, cur_range)) {goto next_d;}
                            if (TYPE && AVOID) {
                                int move_aggro_monsters = 0;
                                if (AVOID == 1) REP(j, depth) if (in_range(nx, ny, monsters[order[i+j]].x, monsters[order[i+j]].y, monsters[order[i+j]].range)) move_aggro_monsters++;
                                if (AVOID == 2) for (int m : attacking_monsters[nx][ny]) move_aggro_monsters += monster_alive[m];
                                if (move_aggro_monsters > aggro_monsters) {goto next_d;}
                            }
                            if (dst2(nx, ny, monsters[order[i+depth]].x, monsters[order[i+depth]].y) < dst2(best_x, best_y, monsters[order[i+depth]].x, monsters[order[i+depth]].y)) {
                                best_x = nx;
                                best_y = ny;
                                update = true;
                            }
                            next_d: ;
                        }
                        if (!update) break;
                        if (in_range(best_x, best_y, monsters[order[i+depth]].x, monsters[order[i+depth]].y, cur_range)) break;
                    }
                    if (depth >= MAX_DEPTH) break;
                }
            }

            if (generate) sol.PB({-1, best_x, best_y});
            x = best_x;
            y = best_y;
            cur_turn++;

            if (TYPE) for (int m : attacking_monsters[x][y]) if (monster_alive[m]) fatigue += monsters[m].attack;
        }

        int turns_to_kill = (m.hp + cur_power - 1) / cur_power;
        REP(_, turns_to_kill) {
            if (cur_turn >= n_turns) return gold;
            if (generate) sol.PB({order[i], -1, -1});
            cur_turn++;

            if (TYPE && _ < turns_to_kill - 1) for (int m : attacking_monsters[x][y]) if (monster_alive[m]) fatigue += monsters[m].attack;
        }



        exp += m.exp;
        while (exp >= levels[level]) level++;

        if (TYPE == 0) {
            gold += m.gold;
        } else {
            gold += (LL)m.gold * 1000 / (1000 + fatigue);
        }

        if (TYPE) {
            monster_alive[order[i]] = 0;
            if (TYPE) for (int m : attacking_monsters[x][y]) if (monster_alive[m]) fatigue += monsters[m].attack;
        }
        
        cur_speed = base_speed * (100 + level * level_speed_coeff) / 100;
        cur_power = base_power * (100 + level * level_power_coeff) / 100;
        cur_range = base_range * (100 + level * level_range_coeff) / 100;

        if (TYPE && fatigue > MAX_FATIGUE) return gold;
    }
    last_turns = cur_turn;
    return gold;
}

int sim(VI &order, bool generate = false) {
    if (test_type == 0) {
        return simx<0>(order, generate);
    } else {
        return simx<1>(order, generate);
    } 
}

const int ANALYZE = 0;
void analyze() {
    int indestructible = 0;
    for (auto &m : monsters) if (m.hp > 1e5) indestructible++;
    DB(indestructible);
}

int main(int argc, char **argv) {
    srand(time(NULL));
    int seed = rand() % 10000;
    REP(i, seed) rng.next();

    assert(argc >= 2);

    test_id = atoi(argv[1]);

    test_type = test_id <= 25 ? 0 : 1;
	load_data(test_id);

    // create levels
    REP(i, 100) levels.PB(1000 + (i+1) * i * 50);
    FOR(i, 1, levels.S) levels[i] += levels[i-1];

    // create attacking monsters
    attacking_monsters = VVVI(map_w + 1, VVI(map_h + 1));
    if (test_type) {
        REP(i, monsters.S) {
            auto &m = monsters[i];
            FOR(x, max(0, m.x - m.range), min(map_w, m.x + m.range) + 1) FOR(y, max(0, m.y - m.range), min(map_h, m.y + m.range) + 1) if (in_range(x, y, m.x, m.y, m.range)) attacking_monsters[x][y].PB(i);
        }
    }

    DB(map_w, map_h, base_speed, base_power, base_range, level_speed_coeff, level_power_coeff, level_range_coeff);
    DB(n_turns);
    DB(monsters.S);

    close_monsters = VVI(monsters.S);
    int n_close_monsters = min(10, (int)monsters.S - 1);
    REP(i, monsters.S) {
        VPII vp;
        REP(j, monsters.S) if (i != j && (monsters[j].gold > 1 || monsters[j].exp > 1) && monsters[j].hp < 1e5) vp.PB({dst2(monsters[i].x, monsters[i].y, monsters[j].x, monsters[j].y), j});
        sort(ALL(vp));
        REP(j, n_close_monsters) close_monsters[i].PB(vp[j].Y);
    }

    if (ANALYZE) {analyze(); return 0;}
    
    bv = -1;

    // optional start from existing solution
    if (argc >= 3) {
        fstream fs(argv[2], fstream::in);
        json j = json::parse(fs);
        int last_monster = -1;
        VI monster_used(monsters.S);
        for (auto move : j["moves"]) {
            if (move["type"] == "attack") {
                int monster = move["target_id"];
                if (last_monster != monster) {
                    last_monster = monster;
                    best_order.PB(monster);
                    monster_used[monster] = 1;
                }
            }
        }
        REP(i, monsters.S) if (!monster_used[i] && (monsters[i].exp > 1 || monsters[i].gold > 1) && monsters[i].hp < 1e5) best_order.PB(i);
        xv = sim(best_order);
        DB(best_order.S);
        DB(xv);
    }

    if (best_order.S == 0) {
        int r_step = 0;
        while (r_step < 100000) {
            static VI order;
            if (order.S == 0) 
                REP(i, monsters.S) if ((monsters[i].exp > 1 || monsters[i].gold > 1) && monsters[i].hp < 1e5) order.PB(i);

            REP(i, order.S) swap(order[i], order[rng.next(i, order.S)]);
            
            int v = sim(order);
            if (v > xv) {
                xv = v;
                best_order = order;
                DB(xv, r_step, elapsed());
            }
            r_step++;
        }
    }

    bv = xv;
    cur_order = best_order;

    if (argc >= 3) xv = 0;

    int sa_step = 0;
    double t0 = xv / 10.0;
    double tn = 1;
    if (test_id == 2 || test_id == 12 || test_id == 13 || test_id == 26 || test_id == 27 || test_id >= 41 && test_id <= 45) tn = 0.1;
    double sa_start = elapsed();
    double t = t0;
    int SA_STEPS = 1000000;
    sim(best_order);
    int best_last_order = last_order;

    int NO_TYPE0 = 0;
    int NO_TYPE1 = 0;
    int NO_TYPE2 = 0;
    VI skip_type(3);
    if (NO_TYPE0) skip_type[0] = 1;
    if (NO_TYPE1) skip_type[1] = 1;
    if (NO_TYPE2) skip_type[2] = 1;
    while (sa_step < SA_STEPS) {
        VI order = cur_order;
        int a = -1, b = -1, p = -1;
        int type = -1;
        do {type = rng.next(3);} while (skip_type[type]);

        if (type == 0) {
            a = rng.next(order.S);
            p = rng.next(order.S);
        } else if (type == 1) {
            p = rng.next(best_last_order+1);
            p = min(p, (int)order.S);
            int m = close_monsters[order[p]][rng.next(close_monsters[order[p]].S)];
            a = -1;
            REP(i, order.S) if (order[i] == m) a = i;
            if (a == -1) continue;
        } else {
            a = rng.next(order.S);
            b = rng.next(order.S);
        }

        if (b == -1) {
            if (a > best_last_order && p > best_last_order) continue;
            int v = order[a];
            order.erase(order.begin() + a);
            order.insert(order.begin() + p, v);
        } else {
            if (a > best_last_order && b > best_last_order) continue;
            swap(order[a], order[b]);
        }

        int av = sim(order);
        double t = t0 * pow(tn/t0, sa_step * 1.0 / SA_STEPS);
        if (av >= bv || rng.next_double() < exp((av - bv) / t)) {
            if (av > xv) {
                xv = av;
                best_order = order;
                DB(xv, sa_step, last_order, last_fatigue, elapsed());
            }
            bv = av;
            best_last_order = last_order;
            cur_order = order;
        }
        sa_step++;
    }

    int final_result = sim(best_order, true);
    save_solution(test_id);
    save_solution(test_id, "all_output/0" + i2s(test_id/10) + i2s(test_id%10) + "-" + i2s(final_result) + ".json");

    cout << test_id << " " << final_result << " " << last_fatigue << " " << last_order << " " << elapsed() << endl;

	return 0;
}
