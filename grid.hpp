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
    int  order       = 2;
    bool deltaBox    = true;
    bool adaptiveBox = true;   // re-tighten the Delta-box as the incumbent grows (Sec 5.4)
};

// lattice pitch in log-coordinates: first-order ln(1+eta), second-order sqrt(2 eta / e)
static inline double gridStep(double eta, int order) {
    return order == 1 ? std::log1p(eta) : std::sqrt(2.0 * eta / std::exp(1.0));
}

static vector<long double> computeDelta(const LayeredGraph& g) {
    const int i = g.i;
    vector<vector<long double>> L(i), R(i);
    for (int j = 0; j < i; ++j) { L[j].assign(g.n[j], 0.0L); R[j].assign(g.n[j], 0.0L); }
    for (int u = 0; u < g.n[0]; ++u) L[0][u] = 1.0L;
    for (int j = 1; j < i; ++j)
        for (int u = 0; u < g.n[j]; ++u)
            for (int w : g.outL(j, u)) L[j][u] += L[j-1][w];
    for (int u = 0; u < g.n[i-1]; ++u) R[i-1][u] = 1.0L;
    for (int j = i - 2; j >= 0; --j)
        for (int u = 0; u < g.n[j]; ++u)
            for (int w : g.outR(j, u)) R[j][u] += R[j+1][w];
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

// Delta-box at incumbent rhoHat (Lemma deltabox): ln m*_j <= ln Delta_j - ln rhoHat, and
// ln m*_j >= max(-sum_{k!=j} hi_k, -ln Delta)  (product-one, and m*_j >= 1/rho* >= 1/Delta).
// Only the bracket rho* in [rhoHat, Delta] enters, so the box always contains m*.
static void deltaBox(const LayeredGraph& g, double  , double rhoHat,
                     const vector<long double>& Delta,
                     vector<double>& boxLo, vector<double>& boxHi) {
    const int i = g.i;
    boxLo.assign(i, 0.0); boxHi.assign(i, 0.0);
    const double lr = std::log(std::max(rhoHat, 1e-12));
    double sumHi = 0.0, lnDmax = 0.0;
    for (int j = 0; j < i; ++j) {
        boxHi[j] = std::log((double)Delta[j]) - lr;
        if (boxHi[j] < 0) boxHi[j] = 0;
        sumHi += boxHi[j];
        lnDmax = std::max(lnDmax, std::log((double)Delta[j]));
    }
    for (int j = 0; j < i; ++j) boxLo[j] = std::max(-(sumHi - boxHi[j]), -lnDmax);
}

// Adaptive Delta-box (Sec 5.4).  All state is sized once in init(); per cell the test is
// i comparisons on precomputed log-coordinates, and the box is recomputed only when the
// incumbent actually moves.  A one-pitch margin per axis keeps the cell of m* (zero-sum
// rounding moves every coordinate by at most one step, Lemma zero-sum rounding).
struct AdaptiveBox {
    int i = 0;
    double tol = 0.0, lnDmax = 0.0, rhoCur = -1.0;
    vector<double> lnDelta, hi, lo;   // per axis
    vector<double> lnGrid;            // cell-major log-coordinates, i per cell

    void init(const vector<vector<double>>& grid, const vector<long double>& Delta, double step) {
        i = (int)Delta.size();
        tol = step + 1e-9;
        lnDelta.resize(i); hi.resize(i); lo.resize(i);
        lnDmax = 0.0;
        for (int j = 0; j < i; ++j) { lnDelta[j] = std::log((double)Delta[j]); lnDmax = std::max(lnDmax, lnDelta[j]); }
        lnGrid.resize(grid.size() * (size_t)i);
        for (size_t c = 0; c < grid.size(); ++c)
            for (int j = 0; j < i; ++j) lnGrid[c * i + j] = std::log(grid[c][j]);
        rhoCur = -1.0;
    }
    // reorder the cached log-coordinates to match a permutation applied to the grid
    void permute(const vector<int>& perm, vector<double>& scratch) {
        scratch.resize(lnGrid.size());
        for (size_t c = 0; c < perm.size(); ++c)
            std::copy_n(lnGrid.begin() + (size_t)perm[c] * i, i, scratch.begin() + c * i);
        lnGrid.swap(scratch);
    }
    inline void tighten(double rho) {
        if (rho <= rhoCur) return;            // box only shrinks; nothing to redo
        rhoCur = rho;
        const double lr = std::log(std::max(rho, 1e-12));
        double sumHi = 0.0;
        for (int j = 0; j < i; ++j) { hi[j] = std::max(0.0, lnDelta[j] - lr); sumHi += hi[j]; }
        for (int j = 0; j < i; ++j) lo[j] = std::max(-(sumHi - hi[j]), -lnDmax);
    }
    inline bool contains(size_t cell) const {
        const double* x = lnGrid.data() + cell * i;
        for (int j = 0; j < i; ++j)
            if (x[j] > hi[j] + tol || x[j] < lo[j] - tol) return false;
        return true;
    }
    inline bool keep(size_t cell, double rho) { tighten(rho); return contains(cell); }
};

// best-first order: nearest-balanced cells first.  One key array, one index sort, one
// in-place permutation -- no per-comparison logs and no per-cell vector copies.
static void sortGridBestFirst(vector<vector<double>>& grid, vector<int>& perm,
                              vector<double>& key, AdaptiveBox* box = nullptr,
                              vector<double>* scratch = nullptr) {
    const size_t n = grid.size();
    key.resize(n); perm.resize(n);
    for (size_t c = 0; c < n; ++c) {
        double s = 0.0;
        for (double m : grid[c]) { const double x = std::log(m); s += x * x; }
        key[c] = s; perm[c] = (int)c;
    }
    std::sort(perm.begin(), perm.end(), [&](int a, int b) { return key[a] < key[b]; });
    // apply the permutation by swapping vector buffers (moves, no element copies)
    vector<vector<double>> tmp(n);
    for (size_t c = 0; c < n; ++c) tmp[c] = std::move(grid[perm[c]]);
    grid.swap(tmp);
    if (box && scratch) box->permute(perm, *scratch);
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
