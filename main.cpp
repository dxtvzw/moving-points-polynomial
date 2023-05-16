#include <bits/stdc++.h>

#include "kosaraju.h"
#include "hamiltonian.h"
#include "matrix.h"

using namespace std;

typedef long long ll;

typedef long double real_t;

const real_t eps = 1e-9;
const real_t inf = 1e18;

bool are_eq(real_t a, real_t b) {
    return abs(a - b) <= eps;
}

real_t get_time() {
    return double(clock()) / CLOCKS_PER_SEC;
}

class CallbackTimer {
public:
    CallbackTimer() {
        start_time = get_time();
    }
    ~CallbackTimer() {
        cout << "\nTime elapsed: " << get_time() - start_time << " s.\n";
    }

private:
    real_t start_time;
};

struct Poly {
    int deg;
    real_t coef;
    real_t operator()(real_t x) const {
        return pow(x, deg) * coef;
    }
    Poly& operator+=(const Poly& ot) {
        if (deg > ot.deg) {
            // do nothing
        } else if (deg < ot.deg) {
            *this = ot;
        } else {
            coef += ot.coef;
        }
        return *this;
    }
    Poly& operator*=(const Poly& ot) {
        deg += ot.deg;
        coef *= ot.coef;
        return *this;
    }
    Poly operator+(const Poly& ot) const {
        Poly tmp = *this;
        return tmp += ot;
    }
    Poly operator*(const Poly& ot) const {
        Poly tmp = *this;
        return tmp *= ot;
    }
    Poly operator*(real_t c) const {
        return {deg, coef * c};
    }
};

struct Edge {
    int to;
    real_t w;
};

const int N = 200 + 10;
vector<pair<int, int>> edges;
vector<Edge> g[N];
priority_queue<real_t, vector<real_t>, greater<>> pq[N];
priority_queue<pair<real_t, int>, vector<pair<real_t, int>>, greater<>> all;

void init(int n) {
    edges.clear();
    for (int i = 1; i <= n; i++) {
        g[i].clear();
        while (!pq[i].empty()) {
            pq[i].pop();
        }
    }
    while (!all.empty()) {
        all.pop();
    }
}

void add_edge(int u, int v, real_t w) {
    g[u].push_back({v, w});
    edges.emplace_back(u, v);
}

real_t get(int u, int v) {
    for (const auto& e : g[u]) {
        if (e.to == v) {
            return e.w;
        }
    }
    return -1;
}

/*
 * Read graph from input
 */
tuple<int, int, int> read_input() {
    int n, m;
    cin >> n >> m;
    for (int i = 1; i <= m; i++) {
        int u, v;
        real_t w;
        cin >> u >> v >> w;
        add_edge(u, v, w);
    }
    int src;
    cin >> src;
    return {n, m, src};
}

/*
 * Read graph without weights from input and initialize weights
 */
tuple<int, int, int> read_simple_input() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<real_t> dist(1, 10);
    int n, m;
    cin >> n >> m;
    for (int i = 1; i <= m; i++) {
        int u, v;
        cin >> u >> v;
        real_t w = dist(gen);
        add_edge(u, v, w);
    }
    int src;
    cin >> src;
    return {n, m, src};
}

/*
 * Read graph without weights from input and initialize weights
 */
tuple<int, int, int> read_input_sqrt_primes() {
    std::random_device rd;
    std::mt19937 gen(rd());
    int n, m;
    cin >> n >> m;
    // cout << "edges with weights as sqrt(w):\n";
    for (int i = 1; i <= m; i++) {
        int u, v, w;
        cin >> u >> v >> w;
        // cout << u << " " << v << " " << w << "\n";
        add_edge(u, v, sqrt(real_t(w)));
    }
    int src;
    cin >> src;
    return {n, m, src};
}

/*
 * Generate non-hamiltonian strongly connected directed graph
 */
tuple<int, int, int> gen_sc_graph() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<real_t> dist(1, 10);
    int n, m, src;
    while (true) {
        n = gen() % 4 + 2;
        m = gen() % (n * (n - 1)) + 1;
        m = min(m, 9);

        vector<pair<int, int>> all_edges;
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                if (i != j) {
                    all_edges.emplace_back(i, j);
                }
            }
        }
        shuffle(all_edges.begin(), all_edges.end(), gen);
        all_edges.erase(all_edges.begin() + m, all_edges.end());

        for (auto [u, v] : all_edges) {
            real_t w = dist(gen);
            add_edge(u, v, w);
        }
        src = 1;
        SCCGraph scc(n);
        for (auto [u, v] : all_edges) {
            scc.add_edge(u, v);
        }
        HamiltonianGraph hmg(n);
        for (auto [u, v] : all_edges) {
            hmg.add_edge(u, v);
        }
        if (scc.is_one_scc() && !hmg.is_hamiltonian()) {
            break;
        } else {
            init(n);
        }
    }
    return {n, m, src};
}

/*
 * Generate hamiltonian directed graph
 */
tuple<int, int, int> gen_hamiltonian_graph() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<real_t> dist(1, 10);
    int n, m, src = 1;
    n = gen() % 10 + 2;
    m = gen() % (n * (n - 1) + 1);
    m = min(m, 20);

    m = max(m, n);

    vector<pair<int, int>> all_edges;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i != j && i % n + 1 != j) {
                all_edges.emplace_back(i, j);
            }
        }
    }
    shuffle(all_edges.begin(), all_edges.end(), gen);
    all_edges.erase(all_edges.begin() + m - n, all_edges.end());
    for (int i = 1; i <= n; i++) {
        all_edges.emplace_back(i, i % n + 1);
    }

    for (auto [u, v] : all_edges) {
        real_t w = dist(gen);
        add_edge(u, v, w);
    }
    return {n, m, src};
}

void print_graph(int n, int m, int src) {
    cout << "---- Graph description ----" << "\n";
    cout << n << " " << m << "\n";
    cout << fixed << setprecision(12);
    for (int i = 1; i <= n; i++) {
        for (const auto& e : g[i]) {
            cout << i << " " << e.to << " " << e.w << "\n";
        }
    }
    cout << src << "\n";
    cout << "---------------------------" << "\n";
    cout << fixed << setprecision(6);
}

/*
 * Rotates cycle until it starts at v
 */
vector<int> fix_cycle(vector<int> cycle, int v) {
    while (cycle[0] != v) {
        rotate(cycle.begin(), cycle.begin() + 1, cycle.end());
    }
    return cycle;
}

bool getbit(int mask, int bit) {
    return (mask >> bit) & 1;
}

vector<int> mask_to_cycle(int mask) {
    int m = edges.size();
    int root = 0;
    map<int, int> ptr;
    for (int i = 0; i < m; i++) {
        if (getbit(mask, i)) {
            auto [u, v] = edges[i];
            root = root == 0 ? u : min(root, u);
            ptr[u] = v;
        }
    }
    vector<int> res = {root};
    for (int v = ptr[root]; v != root; v = ptr[v]) {
        res.push_back(v);
    }
    return res;
}

vector<int> get_cycle_masks(int n) {
    int m = edges.size();
    vector<int> ptr_in(n + 1, 0), ptr_out(n + 1, 0);
    vector<bool> used(n + 1, false);
    vector<int> res;
    for (int mask = 1; mask < (1 << m); mask++) {
        for (int i = 1; i <= n; i++) {
            ptr_in[i] = 0;
            ptr_out[i] = 0;
            used[i] = false;
        }
        bool ok = true;
        for (int i = 0; i < m; i++) {
            if (getbit(mask, i)) {
                auto [u, v] = edges[i];
                if (ptr_in[v] != 0 || ptr_out[u] != 0) {
                    ok = false;
                } else {
                    ptr_in[v] = u;
                    ptr_out[u] = v;
                }
            }
        }
        for (int i = 1; i <= n; i++) {
            if ((ptr_in[i] != 0) != (ptr_out[i] != 0)) {
                ok = false;
            }
        }
        if (!ok) {
            continue;
        } else {
            int cnt = 0;
            for (int i = 1; i <= n; i++) {
                if (!used[i] && ptr_out[i]) {
                    cnt++;
                    for (int j = i; !used[j]; j = ptr_out[j]) {
                        used[j] = true;
                    }
                }
            }
            if (cnt == 1) {
                res.push_back(mask);
            }
        }
    }
    return res;
}

vector<vector<int>> get_cycles(int n) {
    vector<int> cycle_masks = get_cycle_masks(n);
    vector<vector<int>> cycles;
    cycles.reserve(cycle_masks.size());
    for (int mask : cycle_masks) {
        cycles.push_back(mask_to_cycle(mask));
    }
    return cycles;
}

Poly calc_simplex(const vector<vector<int>>& cycles) {
    Poly p = {int(cycles.size()), 1};
    for (int i = 1; i <= cycles.size(); i++) {
        p.coef /= i;
    }
    for (const auto& cycle : cycles) {
        real_t cur = 0;
        for (int i = 0; i < cycle.size(); i++) {
            cur += get(cycle[i], cycle[i + 1 == cycle.size() ? 0 : i + 1]);
        }
        p.coef /= cur;
    }
    return p;
}

Poly calc_simplex(const vector<int>& cycle_masks) {
    vector<vector<int>> cycles;
    cycles.reserve(cycle_masks.size());
    for (int mask : cycle_masks) {
        cycles.push_back(mask_to_cycle(mask));
    }
    return calc_simplex(cycles);
}

int cycle_rank(int n) {
    vector<int> rnk((1 << n), 0);
    for (int mask = 1; mask < (1 << n); mask++) {
        vector<int> pos(n + 1, 0);
        int k = 0;
        for (int i = 1; i <= n; i++) {
            if (getbit(mask, i - 1)) {
                pos[i] = ++k;
            }
        }
        SCCGraph cur(k);
        for (int i = 1; i <= n; i++) {
            if (getbit(mask, i - 1)) {
                for (auto e : g[i]) {
                    if (getbit(mask, e.to - 1)) {
                        cur.add_edge(pos[i], pos[e.to]);
                    }
                }
            }
        }
        vector<vector<int>> comps = cur.get_comps();
        if (comps.size() == __builtin_popcount(mask)) { // graph is acyclic
            rnk[mask] = 0;
        } else if (comps.size() == 1) { // graph is strongly connected
            rnk[mask] = INT_MAX;
            for (int i = 1; i <= n; i++) {
                if (getbit(mask, i - 1)) {
                    rnk[mask] = min(rnk[mask], rnk[mask ^ (1 << (i - 1))]);
                }
            }
            rnk[mask]++;
        } else {
            rnk[mask] = 0;
            for (const auto& comp : comps) {
                int submask = 0;
                for (int i : comp) {
                    submask |= (1 << (i - 1));
                }
                rnk[mask] = max(rnk[mask], rnk[submask]);
            }
        }
    }
    return rnk[(1 << n) - 1];
}

int my_rank(int n) {
    int m = edges.size();
    auto cycle_masks = get_cycle_masks(n);

    int k = cycle_masks.size();
    Matrix<Rational<ll>> a(k, m);
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < m; j++) {
            a[i][j] = getbit(cycle_masks[i], j);
        }
    }

    return rank(a);
}

int min_arc_set(int n) {
    int m = edges.size();
    auto cycle_masks = get_cycle_masks(n);
    int k = cycle_masks.size();
    for (int mask = 0; mask < (1 << m); mask++) {
        bool ok = true;
        for (int i = 0; i < k; i++) {
            ok &= (mask & cycle_masks[i]) != 0;
        }
        if (ok) {
            /*
            cout << "min arc set is:\n";
            for (int i = 0; i < m; i++) {
                if (getbit(mask, i)) {
                    cout << edges[i].first << " " << edges[i].second << "\n";
                }
            }
             */
            return __builtin_popcount(mask);
        }
    }
    return -1;
}

/*
 * All cycles
 */
Poly poly_1(int n, int m, int src) {
    vector<vector<int>> cycles = get_cycles(n);

    for (const auto& cycle : cycles) {
        cout << "cycle: ";
        for (int v : cycle) {
            cout << v << " ";
        }
        cout << "\n";
    }

    Poly p = calc_simplex(cycles);
    real_t sum = 0;
    for (int i = 1; i <= n; i++) {
        for (const auto &e: g[i]) {
            sum += e.w;
        }
    }
    p.coef *= sum * p.deg;
    p.deg--;
    return p;
}

/*
 * All basis of cycles
 */
Poly poly_exact(int n, int m, int src) {
    auto cycle_masks = get_cycle_masks(n);

    auto cycles = get_cycles(n);
    for (int i = 0; i < cycles.size(); i++) {
        const auto& cycle = cycles[i];
        cout << "cycle #" << i + 1 << ": ";
        for (int v : cycle) {
            cout << v << " ";
        }
        cout << "\n";
    }

    int rnk = my_rank(n);
    int arc_rnk = min_arc_set(n);

    cout << "rank of subgroup of flows: " << rnk << "\n";
    // cout << "min arc set: " << arc_rnk << " (useless)\n";

    Poly p = {};
    int k = cycle_masks.size();

    if (k > 20) {
        cout << "Too many cycles in the graph\n";
        exit(0);
    }

    int cnt_basis = 0;

    vector<Poly> all_polies;
    for (int mask = 1; mask < (1 << k); mask++) {
        int t = __builtin_popcount(mask);
        if (t != rnk) {
            continue;
        }
        Matrix<Rational<ll>> a(t, m);
        for (int i = 0, j = 0; i < k; i++) {
            if (getbit(mask, i)) { // using cycle number i
                for (int l = 0; l < m; l++) {
                    if (getbit(cycle_masks[i], l)) {
                        a[j][l] = 1;
                    }
                }
                j++;
            }
        }

        vector<vector<int>> kek(m);
        for (int i = 0; i < k; i++) {
            if (getbit(mask, i)) {
                for (int j = 0; j < m; j++) {
                    if (getbit(cycle_masks[i], j)) {
                        kek[j].push_back(i);
                    }
                }
            }
        }
        if (rank(a) == rnk) {
            cnt_basis++;
            vector<int> cur;
            for (int i = 0; i < k; i++) {
                if (getbit(mask, i)) {
                    cur.push_back(cycle_masks[i]);
                }
            }

            cout << "cycle basis:\n";
            for (int i = 0; i < k; i++) {
                if (getbit(mask, i)) {
                    cout << i + 1 << " ";
                    /*
                    cout << "cycle: ";
                    for (int j : mask_to_cycle(cycle_masks[i])) {
                        cout << j << " ";
                    }
                    cout << "\n";
                     */
                }
            }
            cout << "\n";

            p += calc_simplex(cur);
            all_polies.push_back(calc_simplex(cur));
        }
    }

    cout << "total number of bases: " << cnt_basis << "\n";

    if (cnt_basis == 1) {
        cout << "only 1 basis, useless graph\n";
        // exit(0);
    }

    real_t sum = 0;
    for (int i = 1; i <= n; i++) {
        for (const auto &e: g[i]) {
            sum += e.w;
        }
    }

    p.coef *= sum * p.deg;
    p.deg--;

    // cout << "Comparison of impacts of all bases at P(x):\n";
    for (auto& cur : all_polies) {
        cur.coef *= sum * cur.deg;
        cur.deg--;
        // cout << cur(260.142) << "\n";
    }
    // return all_polies[0] + all_polies[1];
    return p;
}

Poly poly_lb(int n, int m, int src) {
    auto cycle_masks = get_cycle_masks(n);

    int rnk = my_rank(n);

    int k = cycle_masks.size();

    vector<Poly> all_polies;
    for (int mask = 1; mask < (1 << k); mask++) {
        int t = __builtin_popcount(mask);
        if (t != rnk) {
            continue;
        }
        Matrix<Rational<ll>> a(t, m);
        for (int i = 0, j = 0; i < k; i++) {
            if (getbit(mask, i)) { // using cycle number i
                for (int l = 0; l < m; l++) {
                    if (getbit(cycle_masks[i], l)) {
                        a[j][l] = 1;
                    }
                }
                j++;
            }
        }

        vector<vector<int>> kek(m);
        for (int i = 0; i < k; i++) {
            if (getbit(mask, i)) {
                for (int j = 0; j < m; j++) {
                    if (getbit(cycle_masks[i], j)) {
                        kek[j].push_back(i);
                    }
                }
            }
        }
        if (rank(a) == rnk) {
            vector<int> cur;
            for (int i = 0; i < k; i++) {
                if (getbit(mask, i)) {
                    cur.push_back(cycle_masks[i]);
                }
            }
            all_polies.push_back(calc_simplex(cur));
        }
    }

    real_t sum = 0;
    for (int i = 1; i <= n; i++) {
        for (const auto &e: g[i]) {
            sum += e.w;
        }
    }

    for (auto& cur : all_polies) {
        cur.coef *= sum * cur.deg;
        cur.deg--;
    }
    return all_polies[0];
}

Poly poly_ub(int n, int m, int src) {
    auto cycle_masks = get_cycle_masks(n);

    int rnk = my_rank(n);
    Poly p = {};
    int k = cycle_masks.size();

    for (int mask = 1; mask < (1 << k); mask++) {
        int t = __builtin_popcount(mask);
        if (t != rnk) {
            continue;
        }
        Matrix<Rational<ll>> a(t, m);
        for (int i = 0, j = 0; i < k; i++) {
            if (getbit(mask, i)) { // using cycle number i
                for (int l = 0; l < m; l++) {
                    if (getbit(cycle_masks[i], l)) {
                        a[j][l] = 1;
                    }
                }
                j++;
            }
        }

        vector<vector<int>> kek(m);
        for (int i = 0; i < k; i++) {
            if (getbit(mask, i)) {
                for (int j = 0; j < m; j++) {
                    if (getbit(cycle_masks[i], j)) {
                        kek[j].push_back(i);
                    }
                }
            }
        }
        if (rank(a) == rnk) {
            vector<int> cur;
            for (int i = 0; i < k; i++) {
                if (getbit(mask, i)) {
                    cur.push_back(cycle_masks[i]);
                }
            }
            p += calc_simplex(cur);
        }
    }

    real_t sum = 0;
    for (int i = 1; i <= n; i++) {
        for (const auto &e: g[i]) {
            sum += e.w;
        }
    }

    p.coef *= sum * p.deg;
    p.deg--;

    return p;
}

Poly poly_3(int n, int m, int src) {
    Poly p = {};
    p.deg = m - n - 1;
    p.coef = 1;
    return p;
}

void stress_check_poly_deg() {
    while (true) {
        auto [n, m, src] = gen_sc_graph();
        if (my_rank(n) == m - n + 1) {
            cout << "OK" << endl;
            init(n);
        } else {
            cout << "WA\n";
            print_graph(n, m, src);
            break;
        }
    }
    exit(0);
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    freopen("input.txt", "r", stdin);
    cout << fixed << setprecision(6);

    CallbackTimer main_timer;

    // stress_check_poly_deg();

    const int graph_type = 3;
    const bool only_print = false;
    const real_t step_size = 20;
    const int mem_last = 10;
    const real_t max_time = 1.0;

    int n, m, src;

    if (graph_type == 1) {
        auto [n2, m2, src2] = read_input();
        n = n2, m = m2, src = src2;
    } else if (graph_type == 2) {
        auto [n2, m2, src2] = read_simple_input();
        n = n2, m = m2, src = src2;
    } else if (graph_type == 3) {
        auto [n2, m2, src2] = read_input_sqrt_primes();
        n = n2, m = m2, src = src2;
    } else if (graph_type == 4) {
        auto [n2, m2, src2] = gen_sc_graph();
        n = n2, m = m2, src = src2;
    } else {
        auto [n2, m2, src2] = gen_hamiltonian_graph();
        n = n2, m = m2, src = src2;
    }

    print_graph(n, m, src);

    auto p_exact = poly_exact(n, m, src);
    cout << "poly exact: " << p_exact.coef << " * T^" << p_exact.deg << "\n";
    cout << "---------------------------" << "\n";

    auto p_lb = poly_lb(n, m, src);
    cout << "poly lb: " << p_lb.coef << " * T^" << p_lb.deg << "\n";
    cout << "---------------------------" << "\n";

    auto p_ub = poly_ub(n, m, src);
    cout << "poly ub: " << p_ub.coef << " * T^" << p_ub.deg << "\n";
    cout << "---------------------------" << "\n";

    auto p_3 = poly_3(n, m, src);
    // cout << "poly 3: " << p_3.coef << " * T^" << p_3.deg << "\n";
    // cout << "---------------------------" << "\n";

    if (only_print) {
        return 0;
    }

    all.emplace(0, src);
    pq[src].push(0);
    real_t last_time = -inf;

    real_t sum_last = 0;
    queue<real_t> last;

    real_t start_time = get_time();
    while (!all.empty() && get_time() - start_time < max_time) {
        auto [cur_time, v] = all.top();
        all.pop();
        if (pq[v].empty() || !are_eq(cur_time, pq[v].top())) {
            continue;
        }
        if (cur_time - last_time >= step_size) {
            last_time = cur_time;
            int cur = 0;
            for (int i = 1; i <= n; i++) {
                cur += pq[i].size();
            }
            // cout << cur_time << "\t" << cur << "\t" << p_1(cur_time) << "\t" << p_2(cur_time) << "\t" << p_3(cur_time) / cur << "\n";
            real_t x = (cur - p_exact(cur_time)) / p_3(cur_time);
            if (last.size() == mem_last) {
                real_t y = last.front();
                last.pop();
                sum_last -= y;
            }
            if (abs(x) <= 123123123) {
                last.push(x);
                sum_last += x;
            }
            cout << cur_time << "\t" << cur << "\t" << p_exact(cur_time) << "\t" << p_lb(cur_time) << "\t" << p_ub(cur_time) << "\n";
            // cout << cur_time << "\t" << cur << "\t" << p_2(cur_time) << "\t" << x << "\t" << sum_last / (last.empty() ? 1 : last.size()) << "\n";
        }
        while (!pq[v].empty() && are_eq(cur_time, pq[v].top())) {
            pq[v].pop();
        }
        for (auto e : g[v]) {
            pq[e.to].push(cur_time + e.w);
            all.emplace(cur_time + e.w, e.to);
        }
    }
}