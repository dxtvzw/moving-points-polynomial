#include <vector>
#include <algorithm>

#ifndef CPLUSPLUS_HAMILTONIAN_H
#define CPLUSPLUS_HAMILTONIAN_H

class HamiltonianGraph {
    /*
     * this class is to check if a graph is hamiltonian
     * 1-based indexing
     */
    const int n;
    std::vector<std::vector<int>> g;

public:
    HamiltonianGraph(int n) : n(n), g(n + 1) {}

    HamiltonianGraph(const std::vector<std::vector<int>>& g_) : n(g_.size() - 1), g(g_) {
    }

    void add_edge(int u, int v) {
        g[u].push_back(v);
    }

    bool is_hamiltonian() {
        std::vector<int> dp((1 << n), 0);
        std::vector<int> adj(n, 0), adjr(n, 0);
        for (int i = 0; i < n; i++) {
            for (int to : g[i + 1]) {
                to--;
                adj[i] |= (1 << to);
                adjr[to] |= (1 << i);
            }
        }
        dp[1] = 1;
        for (int mask = 2; mask < (1 << n); mask++) {
            for (int i = 1; i < n; i++) {
                if (dp[mask ^ (1 << i)] & adjr[i]) {
                    dp[mask] |= (1 << i);
                }
            }
        }
        return dp[(1 << n) - 1] & adjr[0];
    }
};

#endif //CPLUSPLUS_HAMILTONIAN_H
