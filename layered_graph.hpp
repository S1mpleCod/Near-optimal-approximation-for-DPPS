#pragma once
#include "dpss_common.hpp"
#include <filesystem>
#include <unordered_map>
#include <cctype>
#include <span>
struct LayeredGraph {
    int i = 0;
    vector<int> n;
    // CSR adjacency, one flat array per layer: successors of (j,u) are
    // nbrR[j][offR[j][u] .. offR[j][u+1]), predecessors nbrL[j][offL[j][u] .. offL[j][u+1]).
    vector<vector<int>> offR, nbrR;
    vector<vector<int>> offL, nbrL;

    std::span<const int> outR(int j, int u) const {
        return { nbrR[j].data() + offR[j][u], (size_t)(offR[j][u + 1] - offR[j][u]) };
    }
    std::span<const int> outL(int j, int u) const {
        return { nbrL[j].data() + offL[j][u], (size_t)(offL[j][u + 1] - offL[j][u]) };
    }
    int degR(int j, int u) const { return offR[j][u + 1] - offR[j][u]; }
    int degL(int j, int u) const { return offL[j][u + 1] - offL[j][u]; }
    long long edgeCount() const {
        long long e = 0;
        for (int j = 0; j + 1 < i; ++j) e += offR[j][n[j]];
        return e;
    }

    // --- construction: stage edges, then finalize() into CSR ---
    vector<vector<std::pair<int,int>>> staged;   // per layer j: (u, v) in insertion order
    void beginEdges() { staged.assign(i, {}); }
    void addEdge(int j, int u, int v) { staged[j].push_back({ u, v }); }
    // Stable counting sort: successors of u keep their insertion order, so the CSR
    // reproduces the old per-vertex push_back order exactly.
    void finalize() {
        offR.assign(i, {}); nbrR.assign(i, {});
        for (int j = 0; j < i; ++j) {
            offR[j].assign(n[j] + 1, 0);
            if (j + 1 < i) {
                for (auto& e : staged[j]) ++offR[j][e.first + 1];
                for (int u = 0; u < n[j]; ++u) offR[j][u + 1] += offR[j][u];
                nbrR[j].resize(offR[j][n[j]]);
                vector<int> cur(offR[j].begin(), offR[j].end() - 1);
                for (auto& e : staged[j]) nbrR[j][cur[e.first]++] = e.second;
            }
        }
        staged.clear(); staged.shrink_to_fit();
        buildReverse();
    }
    // predecessors from successors, into existing buffers (grow-only, no per-vertex vectors)
    void buildReverse() {
        auto grow = [](auto& v, size_t sz) { if (v.size() < sz) v.resize(sz); };
        grow(offL, (size_t)i); grow(nbrL, (size_t)i);
        for (int j = 0; j < i; ++j) {
            grow(offL[j], (size_t)n[j] + 1);
            std::fill(offL[j].begin(), offL[j].begin() + n[j] + 1, 0);
        }
        for (int j = 0; j + 1 < i; ++j)
            for (int u = 0; u < n[j]; ++u)
                for (int w : outR(j, u)) ++offL[j + 1][w + 1];
        for (int j = 0; j < i; ++j) {
            for (int u = 0; u < n[j]; ++u) offL[j][u + 1] += offL[j][u];
            grow(nbrL[j], (size_t)offL[j][n[j]]);
        }
        for (int j = 1; j < i; ++j) {           // fill pass, cursor scratch reused across calls
            revCursor_.assign(offL[j].begin(), offL[j].end() - 1);
            for (int u = 0; u < n[j - 1]; ++u)
                for (int w : outR(j - 1, u)) nbrL[j][revCursor_[w]++] = u;
        }
    }
private:
    vector<int> revCursor_;
public:
};

static std::expected<LayeredGraph, std::string>
readTypedEdgeList(std::istream& in) {
    LayeredGraph g;
    if (!(in >> g.i) || g.i < 2)
        return std::unexpected("error: first token must be i >= 2");
    g.n.resize(g.i);
    for (int j = 0; j < g.i; ++j)
        if (!(in >> g.n[j]) || g.n[j] <= 0)
            return std::unexpected(std::format("error: bad layer size for type {}", j));
    g.beginEdges();
    long long j, u, v;
    long long edges = 0;
    while (in >> j >> u >> v) {
        if (j < 0 || j + 1 >= g.i)
            return std::unexpected(std::format(
                "error: edge type {} is not between consecutive meta-path types", j));
        if (u < 0 || u >= g.n[j] || v < 0 || v >= g.n[j + 1])
            return std::unexpected(std::format(
                "error: vertex id out of range on edge ({},{},{})", j, u, v));
        g.addEdge((int)j, (int)u, (int)v);
        ++edges;
    }
    g.finalize();
    std::cerr << "read graph: i=" << g.i << ", |E|=" << edges << ", n=(";
    for (int t = 0; t < g.i; ++t)
        std::cerr << g.n[t] << (t + 1 < g.i ? "," : ")\n");
    return g;
}

static std::expected<LayeredGraph, std::string>
readQueryDir(const std::string& dir) {
    namespace fs = std::filesystem;
    auto trim = [](std::string s) {
        size_t a = s.find_first_not_of(" \t\r\n");
        size_t b = s.find_last_not_of(" \t\r\n");
        return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
    };
    auto idxOf = [](const std::string& name, const char* prefix) -> int {
        std::string p = prefix; if (name.rfind(p, 0) != 0) return -1;
        size_t q = p.size(); long k = 0; bool any = false;
        while (q < name.size() && std::isdigit((unsigned char)name[q])) { k = k * 10 + (name[q] - '0'); ++q; any = true; }
        return any ? (int)k : -1;
    };
    auto endsWith = [](const std::string& s, const char* suf) {
        std::string x = suf; return s.size() >= x.size() && s.compare(s.size() - x.size(), x.size(), x) == 0;
    };
    std::vector<std::pair<int, fs::path>> lfiles, sfiles;
    std::error_code ec;
    for (const auto& e : fs::directory_iterator(dir, ec)) {
        if (ec) break;
        if (!e.is_regular_file()) continue;
        std::string name = e.path().filename().string();
        if (endsWith(name, ".nodes")) { int k = idxOf(name, "layer"); if (k >= 0) lfiles.push_back({k, e.path()}); }
        else if (endsWith(name, ".edges")) { int k = idxOf(name, "step"); if (k >= 0) sfiles.push_back({k, e.path()}); }
    }
    if (ec) return std::unexpected("error: cannot open query directory: " + dir);
    if (lfiles.size() < 2) return std::unexpected("error: need >= 2 layer*_*.nodes files in " + dir);
    std::sort(lfiles.begin(), lfiles.end());
    std::sort(sfiles.begin(), sfiles.end());

    LayeredGraph g;
    g.i = (int)lfiles.size();
    g.n.assign(g.i, 0);
    std::vector<std::unordered_map<std::string, int>> id(g.i);
    for (int j = 0; j < g.i; ++j) {
        std::ifstream f(lfiles[j].second);
        if (!f) return std::unexpected("error: cannot read " + lfiles[j].second.string());
        std::string line;
        while (std::getline(f, line)) {
            std::string s = trim(line);
            if (!s.empty() && !id[j].count(s)) id[j].emplace(std::move(s), (int)id[j].size());
        }
        g.n[j] = (int)id[j].size();
        if (g.n[j] == 0) return std::unexpected(std::format("error: empty layer {} in {}", j, dir));
    }
    g.beginEdges();
    long long edges = 0;
    for (auto& [k, path] : sfiles) {
        if (k < 0 || k + 1 >= g.i) continue;
        auto& a = id[k]; auto& b = id[k + 1];
        std::ifstream f(path);
        if (!f) return std::unexpected("error: cannot read " + path.string());
        std::string line, x, y;
        while (std::getline(f, line)) {
            std::istringstream ss(line);
            if (!(ss >> x >> y)) continue;
            auto ax = a.find(x), by = b.find(y);
            if (ax != a.end() && by != b.end()) { g.addEdge(k, ax->second, by->second); ++edges; continue; }
            auto ay = a.find(y), bx = b.find(x);
            if (ay != a.end() && bx != b.end()) { g.addEdge(k, ay->second, bx->second); ++edges; }
        }
    }
    g.finalize();
    std::cerr << "read query dir: i=" << g.i << ", |E|=" << edges << ", n=(";
    for (int t = 0; t < g.i; ++t) std::cerr << g.n[t] << (t + 1 < g.i ? "," : ")\n");
    return g;
}

static long double totalInstances(const LayeredGraph& g) {
    vector<long double> cur(g.n[0], 1.0L);
    for (int j = 1; j < g.i; ++j) {
        vector<long double> nxt(g.n[j], 0.0L);
        for (int u = 0; u < g.n[j]; ++u) {
            long double s = 0;
            for (int w : g.outL(j, u)) s += cur[w];
            nxt[u] = s;
        }
        cur.swap(nxt);
    }
    long double s = 0;
    for (long double x : cur) s += x;
    return s;
}
