#pragma once
#include "layered_graph.hpp"
struct ReduceStats {
    long long Vbefore = 0, Vafter = 0, Ebefore = 0, Eafter = 0;
    long double Fbefore = 0, Fafter = 0;
    long double Dafter = 1;
    int rounds = 0;
    bool emptied = false;
};

static ReduceStats coReduce(const LayeredGraph& g, const vector<double>& M,
                            double rho, double safety,
                            vector<vector<long double>>& L,
                            vector<vector<long double>>& R,
                            vector<vector<char>>& active) {
    const int i = g.i;
    ReduceStats st;

    auto sweep = [&]() {
        for (int j = 0; j < i; ++j) std::fill(L[j].begin(), L[j].end(), 0.0L);
        for (int u = 0; u < g.n[0]; ++u) if (active[0][u]) L[0][u] = 1;
        for (int j = 1; j < i; ++j)
            for (int u = 0; u < g.n[j]; ++u) if (active[j][u]) {
                long double s = 0;
                for (int w : g.outL(j, u)) if (active[j - 1][w]) s += L[j - 1][w];
                L[j][u] = s;
            }
        for (int j = 0; j < i; ++j) std::fill(R[j].begin(), R[j].end(), 0.0L);
        for (int u = 0; u < g.n[i - 1]; ++u) if (active[i - 1][u]) R[i - 1][u] = 1;
        for (int j = i - 2; j >= 0; --j)
            for (int u = 0; u < g.n[j]; ++u) if (active[j][u]) {
                long double s = 0;
                for (int w : g.outR(j, u)) if (active[j + 1][w]) s += R[j + 1][w];
                R[j][u] = s;
            }
    };
    auto countVEF = [&](long long& V, long long& E, long double& F) {
        V = 0; E = 0; F = 0;
        for (int j = 0; j < i; ++j) for (int u = 0; u < g.n[j]; ++u) if (active[j][u]) ++V;
        for (int j = 0; j + 1 < i; ++j) for (int u = 0; u < g.n[j]; ++u) if (active[j][u])
            for (int w : g.outR(j, u)) if (active[j + 1][w]) ++E;
        for (int u = 0; u < g.n[i - 1]; ++u) if (active[i - 1][u]) F += L[i - 1][u];
    };

    sweep();
    countVEF(st.Vbefore, st.Ebefore, st.Fbefore);

    bool changed = true;
    while (changed) {
        changed = false; ++st.rounds;
        sweep();
        for (int j = 0; j < i; ++j) {
            const long double floor = (long double)rho * (long double)M[j]
                                      / ((long double)g.i * (long double)safety);
            const long double thr = floor - 1e-9L * std::max((long double)1.0L, floor);
            for (int u = 0; u < g.n[j]; ++u)
                if (active[j][u] && L[j][u] * R[j][u] < thr) { active[j][u] = 0; changed = true; }
        }
        for (int j = 0; j < i; ++j) {
            bool any = false;
            for (int u = 0; u < g.n[j]; ++u) if (active[j][u]) { any = true; break; }
            if (!any) { st.emptied = true; break; }
        }
        if (st.emptied) break;
    }
    if (st.emptied) { st.Vafter = st.Eafter = 0; st.Fafter = 0; }
    else {
        sweep(); countVEF(st.Vafter, st.Eafter, st.Fafter);
        long double dm = 1.0L;
        for (int j = 0; j < i; ++j)
            for (int u = 0; u < g.n[j]; ++u)
                if (active[j][u]) dm = std::max(dm, L[j][u] * R[j][u]);
        st.Dafter = dm;
    }
    return st;
}

struct TraceRec { size_t step; double rhoHat; bool solved; double cellRho; };

struct PruneAgg {
    size_t grid = 0, solved = 0, eliminated = 0, boxSkipped = 0;
    long long peelDecCalls = 0;
    double sumVf = 0, sumEf = 0, sumFf = 0;
    double maxVf = 0, maxEf = 0, maxFf = 0;
    long long Ebefore = 0, Eafter = 0;
    long double Fbefore = 0, Fafter = 0;
    int rounds_max = 0; long long rounds_sum = 0;
    long long decKeyCalls = 0, bktMoves = 0;
    double sumAlphaFull = 0, sumAlphaCell = 0;
    long double sumFfull = 0, sumFprun = 0;
    long long sumTcap = 0, sumUsed = 0, iterCells = 0;
    bool traceOn = false;
    vector<TraceRec> trace;
};

static void buildReducedInto(const LayeredGraph& g, const vector<vector<char>>& active,
                             LayeredGraph& r, vector<vector<int>>& backmap,
                             vector<vector<int>>& fwd) {
    const int i = g.i;
    auto grow = [](auto& v, size_t sz) { if (v.size() < sz) v.resize(sz); };
    r.i = i;
    grow(r.n, i); grow(fwd, i); grow(backmap, i); grow(r.offR, i); grow(r.nbrR, i);
    for (int j = 0; j < i; ++j) {
        grow(fwd[j], (size_t)g.n[j]);
        std::fill(fwd[j].begin(), fwd[j].begin() + g.n[j], -1);
        backmap[j].clear();
        r.n[j] = 0;
        for (int u = 0; u < g.n[j]; ++u)
            if (active[j][u]) { fwd[j][u] = r.n[j]++; backmap[j].push_back(u); }
    }
    // reduced CSR written straight into r's reusable buffers (counting pass + fill pass)
    for (int j = 0; j < i; ++j) {
        grow(r.offR[j], (size_t)r.n[j] + 1);
        std::fill(r.offR[j].begin(), r.offR[j].begin() + r.n[j] + 1, 0);
        if (j + 1 < i) {
            for (int u = 0; u < g.n[j]; ++u) if (active[j][u]) {
                int d = 0;
                for (int w : g.outR(j, u)) if (active[j + 1][w]) ++d;
                r.offR[j][fwd[j][u] + 1] = d;
            }
            for (int u = 0; u < r.n[j]; ++u) r.offR[j][u + 1] += r.offR[j][u];
            grow(r.nbrR[j], (size_t)r.offR[j][r.n[j]]);
            for (int u = 0; u < g.n[j]; ++u) if (active[j][u]) {
                int* out = r.nbrR[j].data() + r.offR[j][fwd[j][u]];
                for (int w : g.outR(j, u)) if (active[j + 1][w]) *out++ = fwd[j + 1][w];
            }
        }
    }
    r.buildReverse();
}
