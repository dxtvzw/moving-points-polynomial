#include <algorithm>
#include <vector>

#ifndef CPLUSPLUS_KOSARAJU_H
#define CPLUSPLUS_KOSARAJU_H

class SCCGraph {
    /*
     * this class is to check if a graph is strongly connected
     * 1-based indexing
     */
    const int n;
    std::vector<std::vector<int>> g, gr;
    std::vector<bool> used;
    std::vector<int> post;
    std::vector<std::vector<int>> comps;

    void dfs_inv(int v) {
        used[v] = true;
        for (int to : gr[v]) {
            if (!used[to]) {
                dfs_inv(to);
            }
        }
        post.push_back(v);
    }

    void dfs(int v, std::vector<int>& comp) {
        used[v] = true;
        comp.push_back(v);
        for (int to : g[v]) {
            if (!used[to]) {
                dfs(to, comp);
            }
        }
    }

public:
    SCCGraph(int n) : n(n), g(n + 1), gr(n + 1), used(n + 1, false) {}

    SCCGraph(const std::vector<std::vector<int>>& g_) : n(g_.size() - 1), g(g_), gr(g_.size()), used(g_.size(), false) {
        for (int v = 1; v <= n; v++) {
            for (int to : g[v]) {
                gr[to].push_back(v);
            }
        }
    }

    void add_edge(int u, int v) {
        g[u].push_back(v);
        gr[v].push_back(u);
    }

    int Kosaraju() {
        post.clear();
        comps.clear();
        for (int i = 1; i <= n; i++) {
            used[i] = false;
        }
        for (int i = 1; i <= n; i++) {
            if (!used[i]) {
                dfs_inv(i);
            }
        }
        reverse(post.begin(), post.end());
        for (int i = 1; i <= n; i++) {
            used[i] = false;
        }
        for (int v : post) {
            if (!used[v]) {
                std::vector<int> comp;
                dfs(v, comp);
                comps.push_back(comp);
            }
        }
        reverse(comps.begin(), comps.end());
        // edges go from left to right
        return comps.size();
    }

    std::vector<std::vector<int>> get_comps() {
        if (comps.empty()) {
            Kosaraju();
        }
        return comps;
    }

    bool is_one_scc() {
        return Kosaraju() == 1;
    }
};


#endif //CPLUSPLUS_KOSARAJU_H
