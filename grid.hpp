#pragma once
#include "layered_graph.hpp"
static vector<vector<double>> irmGrid(const vector<int>& n, double eta) {
    const int i = (int)n.size();
    const double step = std::log1p(eta);
    vector<double> logn(i), lo(i), hi(i);
    for (int j = 0; j < i; ++j) logn[j] = std::log((double)n[j]);
    const double total = std::ranges::fold_left(logn, 0.0, std::plus{});
    for (int j = 0; j < i; ++j) {
        lo[j] = (1.0 - i) / i * logn[j];
        hi[j] = (total - logn[j]) / i;
    }
    vector<std::pair<long long, long long>> range(i - 1);
    for (int j = 0; j < i - 1; ++j)
        range[j] = { (long long)std::floor(lo[j] / step),
                     (long long)std::ceil (hi[j] / step) };

    vector<vector<double>> grid;
    vector<long long> k(i - 1);
    for (int j = 0; j < i - 1; ++j) k[j] = range[j].first;
    while (true) {
        double s = 0;
        for (int j = 0; j < i - 1; ++j) s += k[j] * step;
        const double xLast = -s;
        if (xLast >= lo[i - 1] - 1e-9 && xLast <= hi[i - 1] + 1e-9) {
            vector<double> m(i);
            for (int j = 0; j < i - 1; ++j) m[j] = std::exp(k[j] * step);
            m[i - 1] = std::exp(xLast);
            grid.push_back(std::move(m));
        }
        int j = 0;
        while (j < i - 1) {
            if (++k[j] <= range[j].second) break;
            k[j] = range[j].first;
            ++j;
        }
        if (j == i - 1) break;
    }
    return grid;
}

static vector<double> roundAndProject(const vector<double>& mStar, double eta) {
    const double step = std::log1p(eta);
    vector<double> mt(mStar.size());
    double prod = 1.0;
    for (size_t j = 0; j < mStar.size(); ++j) {
        mt[j] = std::exp(std::ceil(std::log(mStar[j]) / step) * step);
        prod *= mt[j];
    }
    const double c = std::pow(prod, 1.0 / (double)mStar.size());
    for (double& m : mt) m /= c;
    return mt;
}

struct GridMode {
    int  order    = 2;
    bool deltaBox = true;
};

static vector<long double> computeDelta(const LayeredGraph& g) {
    const int i = g.i;
    vector<vector<long double>> L(i), R(i);
    for (int j = 0; j < i; ++j) { L[j].assign(g.n[j], 0.0L); R[j].assign(g.n[j], 0.0L); }
    for (int u = 0; u < g.n[0]; ++u) L[0][u] = 1.0L;
    for (int j = 1; j < i; ++j)
        for (int u = 0; u < g.n[j]; ++u)
            for (int w : g.adjL[j][u]) L[j][u] += L[j-1][w];
    for (int u = 0; u < g.n[i-1]; ++u) R[i-1][u] = 1.0L;
    for (int j = i - 2; j >= 0; --j)
        for (int u = 0; u < g.n[j]; ++u)
            for (int w : g.adjR[j][u]) R[j][u] += R[j+1][w];
    vector<long double> Delta(i, 1.0L);
    for (int j = 0; j < i; ++j) {
        long double mx = 0.0L;
        for (int u = 0; u < g.n[j]; ++u) mx = std::max(mx, L[j][u] * R[j][u]);
        Delta[j] = std::max<long double>(mx, 1.0L);
    }
    return Delta;
}

static vector<vector<double>> irmGridSecondOrder(
        const vector<int>& n, double eta,
        const vector<double>* boxLo = nullptr,
        const vector<double>* boxHi = nullptr) {
    const int i = (int)n.size();
    const double e = std::exp(1.0);
    const double step = std::sqrt(2.0 * eta / e);
    vector<double> logn(i), lo(i), hi(i);
    for (int j = 0; j < i; ++j) logn[j] = std::log((double)n[j]);
    const double total = std::ranges::fold_left(logn, 0.0, std::plus{});
    for (int j = 0; j < i; ++j) {
        lo[j] = (1.0 - i) / i * logn[j];
        hi[j] = (total - logn[j]) / i;
        if (boxLo) lo[j] = std::max(lo[j], (*boxLo)[j]);
        if (boxHi) hi[j] = std::min(hi[j], (*boxHi)[j]);
    }
    vector<std::pair<long long,long long>> range(i - 1);
    for (int j = 0; j < i - 1; ++j)
        range[j] = { (long long)std::floor(lo[j] / step),
                     (long long)std::ceil (hi[j] / step) };
    vector<vector<double>> grid;
    if (i == 1) { grid.push_back({1.0}); return grid; }
    vector<long long> k(i - 1);
    for (int j = 0; j < i - 1; ++j) k[j] = range[j].first;
    while (true) {
        long long ksum = 0;
        for (int j = 0; j < i - 1; ++j) ksum += k[j];
        const long long kLast = -ksum;
        const double xLast = kLast * step;
        if (xLast >= lo[i-1] - 1e-9 && xLast <= hi[i-1] + 1e-9) {
            vector<double> m(i);
            for (int j = 0; j < i - 1; ++j) m[j] = std::exp(k[j] * step);
            m[i-1] = std::exp(xLast);
            grid.push_back(std::move(m));
        }
        int j = 0;
        while (j < i - 1) {
            if (++k[j] <= range[j].second) break;
            k[j] = range[j].first; ++j;
        }
        if (j == i - 1) break;
    }
    return grid;
}

static void deltaBox(const LayeredGraph& g, double  , double rhoHat,
                     const vector<long double>& Delta,
                     vector<double>& boxLo, vector<double>& boxHi) {
    const int i = g.i;
    boxLo.assign(i, 0.0); boxHi.assign(i, 0.0);
    const double lr = std::log(std::max(rhoHat, 1e-12));
    double sumHi = 0.0;
    for (int j = 0; j < i; ++j) {
        boxHi[j] = std::log((double)Delta[j]) - lr;
        if (boxHi[j] < 0) boxHi[j] = 0;
        sumHi += boxHi[j];
    }
    for (int j = 0; j < i; ++j) boxLo[j] = -sumHi;
}

static vector<vector<double>> buildGrid(const LayeredGraph& g, double eta,
                                        GridMode gm, double rhoHat = -1.0) {
    vector<double> boxLo, boxHi;
    const bool useBox = gm.deltaBox && rhoHat > 0;
    if (useBox) {
        vector<long double> Delta = computeDelta(g);
        deltaBox(g, eta, rhoHat, Delta, boxLo, boxHi);
    }
    if (gm.order == 1) {
        if (!useBox) return irmGrid(g.n, eta);
        const int i = g.i;
        const double step = std::log1p(eta);
        vector<double> logn(i), lo(i), hi(i);
        for (int j = 0; j < i; ++j) logn[j] = std::log((double)g.n[j]);
        const double total = std::ranges::fold_left(logn, 0.0, std::plus{});
        for (int j = 0; j < i; ++j) {
            lo[j] = std::max((1.0 - i)/i*logn[j], boxLo[j]);
            hi[j] = std::min((total - logn[j])/i, boxHi[j]);
        }
        vector<std::pair<long long,long long>> range(i - 1);
        for (int j = 0; j < i - 1; ++j)
            range[j] = { (long long)std::floor(lo[j]/step), (long long)std::ceil(hi[j]/step) };
        vector<vector<double>> grid; vector<long long> k(i - 1);
        if (i == 1) { grid.push_back({1.0}); return grid; }
        for (int j = 0; j < i - 1; ++j) k[j] = range[j].first;
        while (true) {
            double s = 0; for (int j = 0; j < i - 1; ++j) s += k[j]*step;
            const double xLast = -s;
            if (xLast >= lo[i-1]-1e-9 && xLast <= hi[i-1]+1e-9) {
                vector<double> m(i);
                for (int j = 0; j < i - 1; ++j) m[j] = std::exp(k[j]*step);
                m[i-1] = std::exp(xLast); grid.push_back(std::move(m));
            }
            int j = 0; while (j < i - 1) { if (++k[j] <= range[j].second) break; k[j]=range[j].first; ++j; }
            if (j == i - 1) break;
        }
        if (grid.empty()) grid.push_back(vector<double>(i, 1.0));
        return grid;
    }
    auto grid = useBox ? irmGridSecondOrder(g.n, eta, &boxLo, &boxHi)
                       : irmGridSecondOrder(g.n, eta);
    if (grid.empty()) grid.push_back(vector<double>(g.i, 1.0));
    return grid;
}
