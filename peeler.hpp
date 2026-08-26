#pragma once
#include "grid.hpp"
struct PeelResult {
    double rho = -1.0;
    vector<vector<int>> survivors;
};

class Peeler {
public:
    static inline double delta = 0.01;

    explicit Peeler(const LayeredGraph& g) { init(g); }

    void init(const LayeredGraph& g) {
        i_ = g.i; n_ = g.n;
        auto grow = [](auto& v, size_t sz) { if (v.size() < sz) v.resize(sz); };
        grow(hR_, i_); grow(hL_, i_); grow(oR_, i_); grow(oL_, i_); grow(posR_, i_); grow(posL_, i_);
        grow(alive_, i_); grow(L_, i_); grow(R_, i_); grow(dbuf_, i_); grow(touched_, i_);
        off_.assign(i_ + 1, 0);
        // CSR adjacency: one contiguous Half array per layer, offsets oR_/oL_ (size n+1).
        // Built in two counting passes -- no per-vertex heap vectors.
        for (int j = 0; j < i_; ++j) {
            grow(oR_[j], (size_t)n_[j] + 1); grow(oL_[j], (size_t)n_[j] + 1);
            grow(posR_[j], (size_t)n_[j]);   grow(posL_[j], (size_t)n_[j]);
            std::fill(oR_[j].begin(), oR_[j].begin() + n_[j] + 1, 0);
            std::fill(oL_[j].begin(), oL_[j].begin() + n_[j] + 1, 0);
            off_[j + 1] = off_[j] + n_[j];
        }
        for (int j = 0; j + 1 < i_; ++j)
            for (int u = 0; u < n_[j]; ++u) {
                oR_[j][u + 1] += (int)g.adjR[j][u].size();
                for (int v : g.adjR[j][u]) ++oL_[j + 1][v + 1];
            }
        for (int j = 0; j < i_; ++j) {
            for (int u = 0; u < n_[j]; ++u) { oR_[j][u + 1] += oR_[j][u]; oL_[j][u + 1] += oL_[j][u]; }
            grow(hR_[j], (size_t)oR_[j][n_[j]]); grow(hL_[j], (size_t)oL_[j][n_[j]]);
            std::fill(posR_[j].begin(), posR_[j].begin() + n_[j], 0);   // fill cursors
            std::fill(posL_[j].begin(), posL_[j].begin() + n_[j], 0);
        }
        for (int j = 0; j + 1 < i_; ++j)
            for (int u = 0; u < n_[j]; ++u)
                for (int v : g.adjR[j][u]) {
                    const int kR = posR_[j][u]++, kL = posL_[j + 1][v]++;
                    hR_[j][oR_[j][u] + kR]     = { v, kL };
                    hL_[j + 1][oL_[j + 1][v] + kL] = { u, kR };
                }
        grow(hR0_, i_); grow(hL0_, i_);   // pristine copies, restored per cell by rebind()
        for (int j = 0; j < i_; ++j) {
            hR0_[j].assign(hR_[j].begin(), hR_[j].begin() + oR_[j][n_[j]]);
            hL0_[j].assign(hL_[j].begin(), hL_[j].begin() + oL_[j][n_[j]]);
        }
        for (int j = 0; j < i_; ++j) {
            grow(alive_[j], (size_t)n_[j]); grow(L_[j], (size_t)n_[j]);
            grow(R_[j], (size_t)n_[j]);     grow(dbuf_[j], (size_t)n_[j]);
            std::fill(dbuf_[j].begin(), dbuf_[j].begin() + n_[j], 0.0L);
            touched_[j].clear(); touched_[j].reserve(64);
        }
        const int N = off_[i_];
        grow(vbkt_, (size_t)N); grow(vpos_, (size_t)N);
        grow(load_, (size_t)N);
        grow(sizes_, (size_t)i_);
        deleted_.reserve(N);
        grow(kept_, (size_t)i_);
        for (int j = 0; j < i_; ++j) grow(kept_[j], (size_t)n_[j]);
    }

    // Restore the pristine edge order for a new cell on the same graph: compaction only
    // permutes within each vertex's segment, so two flat copies (no allocation) put every
    // cell back to the state init() left, exactly as a full rebuild would.
    void rebind() {
        for (int j = 0; j < i_; ++j) {
            std::copy(hR0_[j].begin(), hR0_[j].begin() + oR_[j][n_[j]], hR_[j].begin());
            std::copy(hL0_[j].begin(), hL0_[j].begin() + oL_[j][n_[j]], hL_[j].begin());
        }
    }

    PeelResult runPasses(const vector<double>& M, int T) {
        std::fill(load_.begin(), load_.begin() + off_[i_], 0.0L);
        double bestRho = -1.0;
        bestDeleted_.clear();
        for (int pass = 0; pass < T; ++pass) {
            PassRes r = onePass(M);
            if (r.rho > bestRho) {
                bestRho = r.rho;
                bestDeleted_.assign(deleted_.begin(), deleted_.begin() + r.prefix);
            }
        }
        PeelResult res; res.rho = bestRho;
        return res;
    }

    PeelResult runPassesPrimal(const vector<double>& M, int patience, int Tcap, int* used = nullptr) {
        std::fill(load_.begin(), load_.begin() + off_[i_], 0.0L);
        double bestRho = -1.0;
        bestDeleted_.clear();
        int stall = 0, t = 0;
        while (t < Tcap) {
            PassRes r = onePass(M);
            ++t;
            const double tol = 0.1 * delta * std::max(1.0, bestRho);
            const bool improved = r.rho > bestRho + tol;
            if (r.rho > bestRho) {
                bestRho = r.rho;
                bestDeleted_.assign(deleted_.begin(), deleted_.begin() + r.prefix);
            }
            stall = improved ? 0 : stall + 1;
            if (stall >= patience) break;
            if (deadline::hit()) break;
        }
        if (used) *used = t;
        PeelResult res; res.rho = bestRho;
        return res;
    }

    PeelResult runPassesAdaptive(const vector<double>& M, double delta, int Tmax,
                                 double rhoIncumbent = -1.0, double eta = 0.0,
                                 int* used = nullptr) {
        std::fill(load_.begin(), load_.begin() + off_[i_], 0.0L);
        const int N = off_[i_];
        double bestRho = -1.0;
        long double bestHm = -1.0L;
        bestDeleted_.clear();
        int t = 0;
        while (t < Tmax) {
            PassRes r = onePass(M); ++t;
            if (r.rho > bestRho) {
                bestRho = r.rho;
                bestDeleted_.assign(deleted_.begin(), deleted_.begin() + r.prefix);
            }
            if (r.hm > bestHm) bestHm = r.hm;
            long double dual = 0.0L;
            for (int v = 0; v < N; ++v) dual = std::max(dual, load_[v] / (long double)t);
            if (rhoIncumbent > 0.0 &&
                (long double)i_ * (1.0L + (long double)eta) * dual <= (long double)rhoIncumbent)
                break;
            if (bestHm > 0.0L && (1.0L - (long double)delta) * dual <= bestHm) break;
            if (deadline::hit()) break;
        }
        if (used) *used = t;
        PeelResult res; res.rho = bestRho;
        return res;
    }

    PeelResult runPassesFW(const vector<double>& M, double delta, int Tmax,
                           double rhoIncumbent = -1.0, double eta = 0.0,
                           int* used = nullptr) {
        const int N = off_[i_];
        mv_.assign(N, 1.0L);
        for (int j = 0; j < i_; ++j)
            for (int u = 0; u < n_[j]; ++u) mv_[off_[j] + u] = (long double)M[j];
        atomsReset();
        auto& x = fwX_; auto& z = fwZ_; auto& sBuf = fwS_;
        x.assign(N, 0.0L); z.assign(N, 0.0L); sBuf.assign(N, 0.0L);
        double bestRho = -1.0; long double bestHm = -1.0L;
        bestDeleted_.clear();
        onePassFW(M, z, sBuf, bestRho, bestHm);
        atomAdd(sBuf, 1.0L); x.assign(sBuf.begin(), sBuf.end());
        int t = 0;
        while (t < Tmax) {
            ++t;
            for (int v = 0; v < N; ++v) z[v] = x[v] / mv_[v];
            std::fill(sBuf.begin(), sBuf.end(), 0.0L);
            onePassFW(M, z, sBuf, bestRho, bestHm);
            int a = 0; long double amax = -1e300L;
            for (size_t k = 0; k < nAtoms_; ++k) {
                long double dp = 0.0L; const auto& A = fwAtoms_[k];
                for (int v = 0; v < N; ++v) dp += z[v] * A[v];
                if (dp > amax) { amax = dp; a = (int)k; }
            }
            long double zsFW = 0.0L;
            for (int v = 0; v < N; ++v) zsFW += z[v] * sBuf[v];
            long double num = amax - zsFW;
            long double den = 0.0L;
            const auto& Aa = fwAtoms_[a];
            for (int v = 0; v < N; ++v) { long double dv = sBuf[v] - Aa[v]; den += dv * dv / mv_[v]; }
            long double gamma = den > 0.0L ? num / den : 0.0L;
            if (gamma < 0.0L) gamma = 0.0L;
            if (gamma > fwAlpha_[a]) gamma = fwAlpha_[a];
            for (int v = 0; v < N; ++v) x[v] += gamma * (sBuf[v] - Aa[v]);
            fwAlpha_[a] -= gamma;
            atomAdd(sBuf, gamma);
            atomsCompact();
            long double U = 0.0L;
            for (int v = 0; v < N; ++v) U = std::max(U, x[v] / mv_[v]);
            if (rhoIncumbent > 0.0 &&
                (long double)i_ * (1.0L + (long double)eta) * U <= (long double)rhoIncumbent) break;
            if (bestHm > 0.0L && (1.0L - (long double)delta) * U <= bestHm) break;
            if (deadline::hit()) break;
        }
        if (used) *used = t;
        PeelResult res; res.rho = bestRho;
        return res;
    }

    PeelResult runPassesFWWarm(const vector<double>& M, double delta, int Tmax,
                               double rhoIncumbent, double eta,
                               const vector<long double>& warm, int* used = nullptr) {
        const int N = off_[i_];
        mv_.assign(N, 1.0L);
        for (int j = 0; j < i_; ++j)
            for (int u = 0; u < n_[j]; ++u) mv_[off_[j] + u] = (long double)M[j];
        auto& x = fwX_; auto& z = fwZ_; auto& sBuf = fwS_;
        x.assign(warm.begin(), warm.end()); z.assign(N, 0.0L); sBuf.assign(N, 0.0L);
        double bestRho = -1.0; long double bestHm = -1.0L;
        bestDeleted_.clear();
        int t = 0;
        while (t < Tmax) {
            ++t;
            for (int v = 0; v < N; ++v) z[v] = x[v] / mv_[v];
            std::fill(sBuf.begin(), sBuf.end(), 0.0L);
            onePassFW(M, z, sBuf, bestRho, bestHm);
            long double num = 0.0L, den = 0.0L;
            for (int v = 0; v < N; ++v) { long double d = sBuf[v] - x[v]; num += z[v] * (-d); den += d * d / mv_[v]; }
            long double gamma = den > 0.0L ? num / den : 0.0L;
            if (gamma < 0.0L) gamma = 0.0L;
            if (gamma > 1.0L) gamma = 1.0L;
            for (int v = 0; v < N; ++v) x[v] += gamma * (sBuf[v] - x[v]);
            long double U = 0.0L;
            for (int v = 0; v < N; ++v) U = std::max(U, x[v] / mv_[v]);
            if (rhoIncumbent > 0.0 &&
                (long double)i_ * (1.0L + (long double)eta) * U <= (long double)rhoIncumbent) break;
            if (bestHm > 0.0L && (1.0L - (long double)delta) * U <= bestHm) break;
            if (num <= 1e-15L) break;
            if (deadline::hit()) break;
        }
        if (used) *used = t;
        PeelResult res; res.rho = bestRho;
        return res;
    }

    static constexpr int kWarmPeel = 3;
    PeelResult runPassesCertPeel(const vector<double>& M, double delta, int Tmax,
                                 double rhoIncumbent = -1.0, double eta = 0.0,
                                 int* used = nullptr) {
        const int N = off_[i_];
        mv_.assign(N, 1.0L);
        for (int j = 0; j < i_; ++j)
            for (int u = 0; u < n_[j]; ++u) mv_[off_[j] + u] = (long double)M[j];
        std::fill(load_.begin(), load_.begin() + N, 0.0L);
        double bestRho = -1.0; long double bestHm = -1.0L;
        bestDeleted_.clear();
        auto finish = [&](int t) { if (used) *used = t; PeelResult res; res.rho = bestRho; return res; };
        atomsReset();
        auto& prev = fwPrev_; auto& b = fwS_;
        prev.assign(N, 0.0L); b.assign(N, 0.0L);
        const int k0 = std::max(1, std::min((int)kWarmPeel, Tmax));
        int t = 0;
        for (; t < k0; ++t) {
            PassRes r = onePass(M);
            if (r.rho > bestRho) { bestRho = r.rho; bestDeleted_.assign(deleted_.begin(), deleted_.begin() + r.prefix); }
            if (r.hm > bestHm) bestHm = r.hm;
            for (int v = 0; v < N; ++v) { b[v] = (load_[v] - prev[v]) * mv_[v]; prev[v] = load_[v]; }
            atomAdd(b, 0.0L);
            long double dual = 0.0L;
            for (int v = 0; v < N; ++v) dual = std::max(dual, load_[v] / (long double)(t + 1));
            if (rhoIncumbent > 0.0 &&
                (long double)i_ * (1.0L + (long double)eta) * dual <= (long double)rhoIncumbent) { fwX_.assign(N, 0.0L); return finish(t + 1); }
            if (bestHm > 0.0L && (1.0L - (long double)delta) * dual <= bestHm) { fwX_.assign(N, 0.0L); return finish(t + 1); }
            if (deadline::hit()) { fwX_.assign(N, 0.0L); return finish(t + 1); }
        }
        for (auto& a : fwAlpha_) a = 1.0L / (long double)nAtoms_;
        auto& x = fwX_; auto& z = fwZ_; auto& sBuf = fwS_;
        x.assign(N, 0.0L); z.assign(N, 0.0L); sBuf.assign(N, 0.0L);
        for (int v = 0; v < N; ++v) x[v] = load_[v] * mv_[v] / (long double)t;
        while (t < Tmax) {
            ++t;
            for (int v = 0; v < N; ++v) z[v] = x[v] / mv_[v];
            std::fill(sBuf.begin(), sBuf.end(), 0.0L);
            onePassFW(M, z, sBuf, bestRho, bestHm);
            int a = 0; long double amax = -1e300L;
            for (size_t k = 0; k < nAtoms_; ++k) {
                long double dp = 0.0L; const auto& A = fwAtoms_[k];
                for (int v = 0; v < N; ++v) dp += z[v] * A[v];
                if (dp > amax) { amax = dp; a = (int)k; }
            }
            long double zsFW = 0.0L;
            for (int v = 0; v < N; ++v) zsFW += z[v] * sBuf[v];
            long double num = amax - zsFW;
            long double den = 0.0L;
            const auto& Aa = fwAtoms_[a];
            for (int v = 0; v < N; ++v) { long double dv = sBuf[v] - Aa[v]; den += dv * dv / mv_[v]; }
            long double gamma = den > 0.0L ? num / den : 0.0L;
            if (gamma < 0.0L) gamma = 0.0L;
            if (gamma > fwAlpha_[a]) gamma = fwAlpha_[a];
            for (int v = 0; v < N; ++v) x[v] += gamma * (sBuf[v] - Aa[v]);
            fwAlpha_[a] -= gamma;
            atomAdd(sBuf, gamma);
            atomsCompact();
            long double U = 0.0L;
            for (int v = 0; v < N; ++v) U = std::max(U, x[v] / mv_[v]);
            if (rhoIncumbent > 0.0 &&
                (long double)i_ * (1.0L + (long double)eta) * U <= (long double)rhoIncumbent) break;
            if (bestHm > 0.0L && (1.0L - (long double)delta) * U <= bestHm) break;
            if (num <= 1e-15L) break;
            if (deadline::hit()) break;
        }
        return finish(t);
    }

    const vector<long double>& lastFWPoint() const { return fwX_; }
    long double dualUnderM(const vector<long double>& x, const vector<double>& M) const {
        long double U = 0.0L;
        for (int j = 0; j < i_; ++j)
            for (int u = 0; u < n_[j]; ++u)
                U = std::max(U, x[off_[j] + u] / (long double)M[j]);
        return U;
    }

    void extractSurvivors(vector<vector<int>>& out) {
        if ((int)out.size() < i_) out.resize(i_);
        for (int j = 0; j < i_; ++j)
            std::fill(kept_[j].begin(), kept_[j].begin() + n_[j], (char)1);
        for (const auto& d : bestDeleted_) kept_[d.first][d.second] = 0;
        for (int j = 0; j < i_; ++j) {
            out[j].clear();
            for (int u = 0; u < n_[j]; ++u) if (kept_[j][u]) out[j].push_back(u);
        }
    }

private:
    struct PassRes { double rho = -1.0; int prefix = 0; long double hm = -1.0L; };

    PassRes onePass(const vector<double>& M) {
        resetState();
        for (int j = 0; j < i_; ++j) sizes_[j] = n_[j];
        long double F = 0;
        for (int u = 0; u < n_[0]; ++u) F += cnt(0, u);
        buildBucketsAug(M, load_);

        deleted_.clear();
        double bestRho = -1.0; int bestPrefix = -1;
        long double bestHm = -1.0L;
        auto density = [&]() -> double {
            long double gm = 1.0L;
            for (int j = 0; j < i_; ++j) { if (!sizes_[j]) return -1.0; gm *= sizes_[j]; }
            return (double)(F / std::pow((double)gm, 1.0 / i_));
        };
        auto hmObj = [&]() -> long double {
            long double w = 0.0L;
            for (int j = 0; j < i_; ++j) { if (!sizes_[j]) return -1.0L; w += (long double)M[j] * sizes_[j]; }
            return w > 0.0L ? F / w : -1.0L;
        };
        if (double r = density(); r > bestRho) { bestRho = r; bestPrefix = 0; }
        if (long double h = hmObj(); h > bestHm) bestHm = h;

        int v;
        while ((v = popMin()) >= 0) {
            const int j = layerOf(v), u = v - off_[j];
            load_[v] += cnt(j, u) / (long double)M[j];
            alive_[j][u] = 0; --sizes_[j];
            deleted_.push_back({ j, u });
            F -= cnt(j, u);
            detachAndPropagateAug(j, u, M, load_);
            L_[j][u] = R_[j][u] = 0;
            if (double r = density(); r > bestRho) { bestRho = r; bestPrefix = (int)deleted_.size(); }
            if (long double h = hmObj(); h > bestHm) bestHm = h;
            if (sizes_[j] == 0) break;
            if (deadline::hit()) break;
        }
        return { bestRho, bestPrefix, bestHm };
    }

public:

private:
    struct Half { int v; int m; };

    int i_;
    vector<int> n_;
    vector<int> off_;
    vector<vector<Half>> hR_, hL_;     // CSR halves per layer
    vector<vector<Half>> hR0_, hL0_;   // pristine order (see rebind)
    vector<vector<int>>  oR_, oL_;     // CSR offsets per layer (n_j + 1 entries)
    vector<vector<int>>  posR_, posL_; // live prefix length of each vertex's segment
    inline int degR(int j, int u) const { return oR_[j][u + 1] - oR_[j][u]; }
    inline int degL(int j, int u) const { return oL_[j][u + 1] - oL_[j][u]; }
    inline Half*       segR(int j, int u)       { return hR_[j].data() + oR_[j][u]; }
    inline Half*       segL(int j, int u)       { return hL_[j].data() + oL_[j][u]; }
    inline Half&       nR(int j, int u, int k)       { return hR_[j][oR_[j][u] + k]; }
    inline const Half& nR(int j, int u, int k) const { return hR_[j][oR_[j][u] + k]; }
    inline Half&       nL(int j, int u, int k)       { return hL_[j][oL_[j][u] + k]; }
    inline const Half& nL(int j, int u, int k) const { return hL_[j][oL_[j][u] + k]; }
    vector<vector<char>> alive_;
    vector<vector<long double>> L_, R_;
    vector<vector<long double>> dbuf_;
    vector<vector<int>> touched_;

    vector<long double> load_;
    vector<int> order_;
    vector<long double> mv_;
    vector<long double> fwX_;
    // Frank-Wolfe scratch, sized once per graph and reused across cells and passes
    vector<long double> fwZ_, fwS_, fwPrev_;
    vector<vector<long double>> fwAtoms_;   // buffer pool; only the first nAtoms_ are live
    vector<long double> fwAlpha_;
    size_t nAtoms_ = 0;

    inline void atomsReset() { nAtoms_ = 0; fwAlpha_.clear(); }
    inline void atomAdd(const vector<long double>& src, long double a) {
        if (nAtoms_ == fwAtoms_.size()) fwAtoms_.emplace_back();
        fwAtoms_[nAtoms_].assign(src.begin(), src.end());   // reuses the buffer's capacity
        fwAlpha_.push_back(a); ++nAtoms_;
    }
    inline void atomsCompact() {   // drop atoms with vanishing weight; buffers stay in the pool
        size_t w = 0;
        for (size_t k = 0; k < nAtoms_; ++k)
            if (fwAlpha_[k] > 1e-15L) { if (w != k) { fwAtoms_[w].swap(fwAtoms_[k]); fwAlpha_[w] = fwAlpha_[k]; } ++w; }
        nAtoms_ = w; fwAlpha_.resize(w);
    }
    vector<int> sizes_;
    vector<std::pair<int,int>> deleted_;
    vector<std::pair<int,int>> bestDeleted_;
    vector<vector<char>> kept_;

    vector<vector<int>> bkt_;
    vector<int> vbkt_, vpos_;
    int cur_ = 0, B_ = 0;
    double log2Min_ = 0;
    double lstep2_ = 0;
public:
    long long nDecCalls_ = 0, nBktMoves_ = 0;
private:

    long double cnt(int j, int u) const { return L_[j][u] * R_[j][u]; }
    int layerOf(int v) const {
        int j = 0;
        while (off_[j + 1] <= v) ++j;
        return j;
    }

    static double log2MantTable(int idx) {
        static double tbl[4096];
        static bool init = false;
        if (!init) {
            for (int k = 0; k < 4096; ++k)
                tbl[k] = std::log2(0.5 + (k + 0.5) / 8192.0);
            init = true;
        }
        return tbl[idx];
    }
    int bucketOf(long double key) const {
        if (key <= 0) return 0;
        int e;
        const double f = (double)frexpl(key, &e);
        const int idx = (int)((f - 0.5) * 8192.0);
        const double lg2 = (double)e + log2MantTable(idx < 0 ? 0 : (idx > 4095 ? 4095 : idx));
        int b = 1 + (int)((lg2 - log2Min_) / lstep2_);
        if (b < 1) b = 1;
        if (b > B_ - 1) b = B_ - 1;
        return b;
    }

    inline void removeFromBucket(int v) {
        const int b = vbkt_[v], p = vpos_[v];
        const int last = bkt_[b].back();
        bkt_[b][p] = last;
        vpos_[last] = p;
        bkt_[b].pop_back();
        vbkt_[v] = vpos_[v] = -1;
    }

    inline void decreaseKey(int v, long double key) {
        ++nDecCalls_;
        const int b = bucketOf(key);
        if (b == vbkt_[v]) return;
        ++nBktMoves_;
        removeFromBucket(v);
        vbkt_[v] = b;
        vpos_[v] = (int)bkt_[b].size();
        bkt_[b].push_back(v);
        if (b < cur_) cur_ = b;
    }

    inline int popMin() {
        while (cur_ < B_ && bkt_[cur_].empty()) ++cur_;
        if (cur_ >= B_) return -1;
        const int v = bkt_[cur_].back();
        bkt_[cur_].pop_back();
        vbkt_[v] = vpos_[v] = -1;
        return v;
    }

    void resetState() {
        for (int j = 0; j < i_; ++j) {
            std::fill(alive_[j].begin(), alive_[j].begin() + n_[j], (char)1);
            for (int u = 0; u < n_[j]; ++u) {
                posR_[j][u] = degR(j, u);
                posL_[j][u] = degL(j, u);
            }
            std::fill(L_[j].begin(), L_[j].begin() + n_[j], 0);
            std::fill(R_[j].begin(), R_[j].begin() + n_[j], 0);
            std::fill(dbuf_[j].begin(), dbuf_[j].begin() + n_[j], 0.0L);
            touched_[j].clear();
        }
        for (int u = 0; u < n_[0]; ++u) L_[0][u] = 1;
        for (int j = 1; j < i_; ++j)
            for (int u = 0; u < n_[j]; ++u)
            {
                const Half* seg = segL(j, u); const int p = posL_[j][u];
                long double acc = 0;
                for (int k = 0; k < p; ++k) acc += L_[j - 1][seg[k].v];
                L_[j][u] += acc;
            }
        for (int u = 0; u < n_[i_ - 1]; ++u) R_[i_ - 1][u] = 1;
        for (int j = i_ - 2; j >= 0; --j)
            for (int u = 0; u < n_[j]; ++u)
            {
                const Half* seg = segR(j, u); const int p = posR_[j][u];
                long double acc = 0;
                for (int k = 0; k < p; ++k) acc += R_[j + 1][seg[k].v];
                R_[j][u] += acc;
            }
    }

    inline void compactL(int jw, int w, int at) {
        int& p = posL_[jw][w];
        --p;
        if (at != p) {
            Half* seg = segL(jw, w);
            std::swap(seg[at], seg[p]);
            const Half a = seg[at]; segR(jw - 1, a.v)[a.m].m = at;
            const Half b = seg[p];  segR(jw - 1, b.v)[b.m].m = p;
        }
    }
    inline void compactR(int jw, int w, int at) {
        int& p = posR_[jw][w];
        --p;
        if (at != p) {
            Half* seg = segR(jw, w);
            std::swap(seg[at], seg[p]);
            const Half a = seg[at]; segL(jw + 1, a.v)[a.m].m = at;
            const Half b = seg[p];  segL(jw + 1, b.v)[b.m].m = p;
        }
    }

    inline long double augKey(int j, int u, const vector<long double>& load) const {
        return (load[off_[j] + u] + cnt(j, u)) / 1.0L;
    }

    void buildBucketsAug(const vector<double>& M, const vector<long double>& load) {
        double maxM = 0;
        for (double m : M) maxM = std::max(maxM, m);
        long double maxKey = 0;
        for (int j = 0; j < i_; ++j)
            for (int u = 0; u < n_[j]; ++u) {
                const long double key = (load[off_[j]+u] + cnt(j,u)) / (long double)M[j];
                if (key > maxKey) maxKey = key;
            }
        log2Min_ = std::log2(1.0 / maxM) - 1e-9;
        lstep2_ = std::log1p(delta) / std::log(2.0);
        const double span = (maxKey > 0 ? (double)std::log2((double)maxKey) : log2Min_) - log2Min_;
        B_ = 2 + (int)(std::max(0.0, span) / lstep2_) + 2;
        if ((int)bkt_.size() < B_) bkt_.resize(B_);
        for (int b = 0; b < B_; ++b) bkt_[b].clear();
        cur_ = 0;
        for (int j = 0; j < i_; ++j)
            for (int u = 0; u < n_[j]; ++u) {
                const int v = off_[j] + u;
                const int b = bucketOf((load[v] + cnt(j,u)) / (long double)M[j]);
                vbkt_[v] = b; vpos_[v] = (int)bkt_[b].size();
                bkt_[b].push_back(v);
            }
    }

    void detachAndPropagateAug(int j, int u, const vector<double>& M,
                               const vector<long double>& load) {
        if (j + 1 < i_) {
            const long double dl = L_[j][u];
            const Half* segU = segR(j, u); const int pu = posR_[j][u];
            for (int k = 0; k < pu; ++k) {
                const Half e = segU[k];
                if (dl > 0) addDelta(j + 1, e.v, dl);
                compactL(j + 1, e.v, e.m);
            }
            for (int jj = j + 1; jj < i_ && !touched_[jj].empty(); ++jj) {
                for (size_t t = 0; t < touched_[jj].size(); ++t) {
                    const int w = touched_[jj][t];
                    const long double d = dbuf_[jj][w];
                    L_[jj][w] -= d;
                    if (alive_[jj][w])
                        decreaseKey(off_[jj] + w,
                                    (load[off_[jj]+w] + cnt(jj,w)) / (long double)M[jj]);
                    if (jj + 1 < i_)
                    {
                        const Half* seg = segR(jj, w); const int p = posR_[jj][w];
                        for (int k = 0; k < p; ++k) addDelta(jj + 1, seg[k].v, d);
                    }
                    dbuf_[jj][w] = 0;
                }
                touched_[jj].clear();
            }
        }
        if (j - 1 >= 0) {
            const long double dr = R_[j][u];
            const Half* segU = segL(j, u); const int pu = posL_[j][u];
            for (int k = 0; k < pu; ++k) {
                const Half e = segU[k];
                if (dr > 0) addDelta(j - 1, e.v, dr);
                compactR(j - 1, e.v, e.m);
            }
            for (int jj = j - 1; jj >= 0 && !touched_[jj].empty(); --jj) {
                for (size_t t = 0; t < touched_[jj].size(); ++t) {
                    const int w = touched_[jj][t];
                    const long double d = dbuf_[jj][w];
                    R_[jj][w] -= d;
                    if (alive_[jj][w])
                        decreaseKey(off_[jj] + w,
                                    (load[off_[jj]+w] + cnt(jj,w)) / (long double)M[jj]);
                    if (jj - 1 >= 0)
                    {
                        const Half* seg = segL(jj, w); const int p = posL_[jj][w];
                        for (int k = 0; k < p; ++k) addDelta(jj - 1, seg[k].v, d);
                    }
                    dbuf_[jj][w] = 0;
                }
                touched_[jj].clear();
            }
        }
    }

    inline void addDelta(int jj, int x, long double d) {
        if (dbuf_[jj][x] == 0) touched_[jj].push_back(x);
        dbuf_[jj][x] += d;
    }

    void detachAndPropagateNoKey(int j, int u) {
        if (j + 1 < i_) {
            const long double dl = L_[j][u];
            const Half* segU = segR(j, u); const int pu = posR_[j][u];
            for (int k = 0; k < pu; ++k) {
                const Half e = segU[k];
                if (dl > 0) addDelta(j + 1, e.v, dl);
                compactL(j + 1, e.v, e.m);
            }
            for (int jj = j + 1; jj < i_ && !touched_[jj].empty(); ++jj) {
                for (size_t t = 0; t < touched_[jj].size(); ++t) {
                    const int w = touched_[jj][t];
                    const long double d = dbuf_[jj][w];
                    L_[jj][w] -= d;
                    if (jj + 1 < i_)
                    {
                        const Half* seg = segR(jj, w); const int p = posR_[jj][w];
                        for (int k = 0; k < p; ++k) addDelta(jj + 1, seg[k].v, d);
                    }
                    dbuf_[jj][w] = 0;
                }
                touched_[jj].clear();
            }
        }
        if (j - 1 >= 0) {
            const long double dr = R_[j][u];
            const Half* segU = segL(j, u); const int pu = posL_[j][u];
            for (int k = 0; k < pu; ++k) {
                const Half e = segU[k];
                if (dr > 0) addDelta(j - 1, e.v, dr);
                compactR(j - 1, e.v, e.m);
            }
            for (int jj = j - 1; jj >= 0 && !touched_[jj].empty(); --jj) {
                for (size_t t = 0; t < touched_[jj].size(); ++t) {
                    const int w = touched_[jj][t];
                    const long double d = dbuf_[jj][w];
                    R_[jj][w] -= d;
                    if (jj - 1 >= 0)
                    {
                        const Half* seg = segL(jj, w); const int p = posL_[jj][w];
                        for (int k = 0; k < p; ++k) addDelta(jj - 1, seg[k].v, d);
                    }
                    dbuf_[jj][w] = 0;
                }
                touched_[jj].clear();
            }
        }
    }

    void onePassFW(const vector<double>& M, const vector<long double>& z,
                   vector<long double>& s, double& bestRho, long double& bestHm) {
        resetState();
        for (int j = 0; j < i_; ++j) sizes_[j] = n_[j];
        const int N = off_[i_];
        long double F = 0;
        for (int u = 0; u < n_[0]; ++u) F += cnt(0, u);
        order_.resize(N);
        for (int v = 0; v < N; ++v) order_[v] = v;
        std::sort(order_.begin(), order_.end(),
                  [&](int a, int b) { return z[a] < z[b]; });
        deleted_.clear();
        auto density = [&]() -> double {
            long double gm = 1.0L;
            for (int j = 0; j < i_; ++j) { if (!sizes_[j]) return -1.0; gm *= sizes_[j]; }
            return (double)(F / std::pow((double)gm, 1.0 / i_));
        };
        auto hmObj = [&]() -> long double {
            long double w = 0.0L;
            for (int j = 0; j < i_; ++j) { if (!sizes_[j]) return -1.0L; w += (long double)M[j] * sizes_[j]; }
            return w > 0.0L ? F / w : -1.0L;
        };
        if (double r = density(); r > bestRho) { bestRho = r; bestDeleted_.clear(); }
        if (long double h = hmObj(); h > bestHm) bestHm = h;
        for (int idx = 0; idx < N; ++idx) {
            const int v = order_[idx];
            const int j = layerOf(v), u = v - off_[j];
            if (!alive_[j][u]) continue;
            const long double c = cnt(j, u);
            s[v] = c;
            F -= c;
            alive_[j][u] = 0; --sizes_[j];
            deleted_.push_back({ j, u });
            detachAndPropagateNoKey(j, u);
            L_[j][u] = R_[j][u] = 0;
            if (double r = density(); r > bestRho) {
                bestRho = r;
                bestDeleted_.assign(deleted_.begin(), deleted_.end());
            }
            if (long double h = hmObj(); h > bestHm) bestHm = h;
            if (deadline::hit()) break;
        }
    }

};
