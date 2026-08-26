#pragma once
#include "dpss_common.hpp"
#include <filesystem>
#include <unordered_map>
#include <cctype>
struct LayeredGraph {
    int i = 0;
    vector<int> n;
    vector<vector<vector<int>>> adjR;
    vector<vector<vector<int>>> adjL;

    void buildReverse() {
        adjL.assign(i, {});
        for (int j = 0; j < i; ++j) adjL[j].assign(n[j], {});
        for (int j = 0; j + 1 < i; ++j)
            for (int u = 0; u < n[j]; ++u)
                for (int w : adjR[j][u]) adjL[j + 1][w].push_back(u);
    }
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
    g.adjR.assign(g.i, {});
    for (int j = 0; j < g.i; ++j) g.adjR[j].assign(g.n[j], {});
    long long j, u, v;
    long long edges = 0;
    while (in >> j >> u >> v) {
        if (j < 0 || j + 1 >= g.i)
            return std::unexpected(std::format(
                "error: edge type {} is not between consecutive meta-path types", j));
        if (u < 0 || u >= g.n[j] || v < 0 || v >= g.n[j + 1])
            return std::unexpected(std::format(
                "error: vertex id out of range on edge ({},{},{})", j, u, v));
        g.adjR[j][u].push_back((int)v);
        ++edges;
    }
    g.buildReverse();
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
    g.adjR.assign(g.i, {});
    for (int j = 0; j < g.i; ++j) g.adjR[j].assign(g.n[j], {});
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
            if (ax != a.end() && by != b.end()) { g.adjR[k][ax->second].push_back(by->second); ++edges; continue; }
            auto ay = a.find(y), bx = b.find(x);
            if (ay != a.end() && bx != b.end()) { g.adjR[k][ay->second].push_back(bx->second); ++edges; }
        }
    }
    g.buildReverse();
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
            for (int w : g.adjL[j][u]) s += cur[w];
            nxt[u] = s;
        }
        cur.swap(nxt);
    }
    long double s = 0;
    for (long double x : cur) s += x;
    return s;
}
