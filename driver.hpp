#pragma once
#include "peeler.hpp"
#include "config.hpp"
#include "coreduce.hpp"
#include "grid.hpp"
static GridResult gridAlgorithm(const LayeredGraph& g, double eta,
                                size_t maxPeels = 0, int passes = 1,
                                GridMode gm = {}) {
    GridResult best;
    Peeler peeler(g);
    peeler.nDecCalls_ = peeler.nBktMoves_ = 0;
    double rhoHat = -1.0;
    if (gm.deltaBox) {
        vector<double> ones(g.i, 1.0);
        rhoHat = peeler.runPasses(ones, 1).rho;
    }
    auto grid = buildGrid(g, eta, gm, rhoHat);
    if (maxPeels && grid.size() > maxPeels) grid.resize(maxPeels);
    best.gridSize = grid.size();
    const bool adaptive = gm.deltaBox && gm.adaptiveBox && rhoHat > 0;
    vector<long double> Dv = computeDelta(g);
    AdaptiveBox box;
    if (adaptive) {   // near-balanced cells first, so the incumbent tightens the box early
        box.init(grid, Dv, gridStep(eta, gm.order));
        vector<int> perm; vector<double> key, scratch;
        sortGridBestFirst(grid, perm, key, &box, &scratch);
    }
    double Dmax = 1.0, alphaF = 1.0;
    long Utot = 0;
    {
        long double dm = 1.0L; for (long double d : Dv) dm = std::max(dm, d);
        Dmax = (double)dm;
        for (int j = 0; j < g.i; ++j) Utot += g.n[j];
        if (rhoHat > 0) alphaF = std::max(1.0, (double)((long double)g.i * dm / rhoHat));
    }
    double rhoStrong = rhoHat;
    int winCell = -1;
    if (gCertify && passes <= 0) {
        for (size_t k = 0; k < grid.size(); ++k) {
            if (deadline::hit()) break;
            if (adaptive && !box.keep(k, rhoStrong)) continue;
            peeler.rebind();
            double rp = peeler.runPassesPrimal(grid[k], gPatience, kAdaptiveTmax).rho;
            if (rp > rhoStrong) { rhoStrong = rp; winCell = (int)k; }
        }
    }
    long uSum = 0; int uMin = (1 << 30), uMax = 0, capMax = 0;
    if (gCertifyShared && passes <= 0) {
        std::vector<std::vector<long double>> pool;
        long pruned = 0;
        std::vector<int> order(grid.size());
        for (size_t k = 0; k < grid.size(); ++k) order[k] = (int)k;
        if (winCell > 0) std::swap(order[0], order[(size_t)winCell]);
        const int Tcap = uniformPassBudget(g.i, Dmax, Utot, Peeler::delta, gIterConst);
        for (int oi = 0; oi < (int)order.size(); ++oi) {
            if (deadline::hit()) break;
            auto& M = grid[order[oi]];
            const double inc = std::max({best.rho, rhoStrong, rhoHat});
            if (adaptive && !box.keep((size_t)order[oi], inc)) { ++best.boxSkipped; continue; }
            int warmIdx = -1;
            if (!pool.empty()) {
                long double Ub = 1e300L;
                for (size_t p = 0; p < pool.size(); ++p) {
                    long double u = peeler.dualUnderM(pool[p], M);
                    if (u < Ub) { Ub = u; warmIdx = (int)p; }
                }
                if ((long double)g.i * (1.0L + (long double)eta) * Ub
                        <= (long double)inc * (1.0L + (long double)gShareEps)) { ++pruned; continue; }
            }
            peeler.rebind();
            int used = 0;
            PeelResult r = (warmIdx >= 0)
                ? peeler.runPassesFWWarm(M, Peeler::delta, Tcap, inc, eta, pool[warmIdx], &used)
                : peeler.runPassesFW(M, Peeler::delta, Tcap, inc, eta, &used);
            pool.push_back(peeler.lastFWPoint());
            uSum += used; uMin = std::min(uMin, used); uMax = std::max(uMax, used);
            capMax = std::max(capMax, Tcap);
            if (r.rho > best.rho) { best.rho = r.rho; peeler.extractSurvivors(best.survivors); best.M = M; }
        }
        std::fprintf(stderr,
            "[shared-cert] FW-solved %zu / pruned %ld of %zu cells (pool=%zu); passes used min=%d mean=%.1f max=%d\n",
            grid.size() - (size_t)pruned, pruned, grid.size(), pool.size(),
            (uMin==(1<<30)?0:uMin), (grid.size()>(size_t)pruned)?(double)uSum/(grid.size()-pruned):0.0, uMax);
        best.decKeyCalls = peeler.nDecCalls_; best.bktMoves = peeler.nBktMoves_;
        return best;
    }
    for (size_t c = 0; c < grid.size(); ++c) {
        if (deadline::hit()) break;
        if (adaptive && !box.keep(c, std::max({best.rho, rhoStrong, rhoHat}))) { ++best.boxSkipped; continue; }
        const auto& M = grid[c];
        peeler.rebind();   // same graph: restore pristine edge order (two memcpys, no rebuild)
        PeelResult r;
        if (passes > 0) { r = peeler.runPasses(M, passes); }
        else {
            const int Tcap = uniformPassBudget(g.i, Dmax, Utot, Peeler::delta, gIterConst);
            int used = 0;
            if (gCertPeel)
                r = peeler.runPassesCertPeel(M, Peeler::delta, Tcap,
                                             std::max({best.rho, rhoStrong, rhoHat}), eta, &used);
            else if (gCertify)
                r = peeler.runPassesFW(M, Peeler::delta, Tcap,
                                       std::max({best.rho, rhoStrong, rhoHat}), eta, &used);
            else if (gSafe)
                r = peeler.runPassesAdaptive(M, Peeler::delta, Tcap,
                                             std::max(best.rho, rhoHat), eta, &used);
            else
                r = peeler.runPassesPrimal(M, gPatience, Tcap, &used);
            uSum += used; uMin = std::min(uMin, used); uMax = std::max(uMax, used);
            capMax = std::max(capMax, Tcap);
        }
        if (r.rho > best.rho) {
            best.rho = r.rho;
            peeler.extractSurvivors(best.survivors);
            best.M = M;
        }
    }
    if (passes <= 0 && !grid.empty())
        std::fprintf(stderr,
            "[passes used] alpha_f=%.2f patience=%d over %zu cells: "
            "min=%d mean=%.1f max=%d sum=%ld  (theoretical cap max=%d)\n",
            alphaF, gPatience, grid.size(), uMin,
            (double)uSum / grid.size(), uMax, uSum, capMax);
    best.decKeyCalls = peeler.nDecCalls_;
    best.bktMoves = peeler.nBktMoves_;
    return best;
}

static GridResult gridAlgorithmPruned(const LayeredGraph& g, double eta, int passes,
                                      double rhoFixed, PruneAgg& agg, GridMode gm = {},
                                      bool bestFirst = true) {
    GridResult best;
    Peeler peeler(g);
    peeler.nDecCalls_ = peeler.nBktMoves_ = 0;
    double rhoBoot = -1.0;
    if (gm.deltaBox) {
        vector<double> ones(g.i, 1.0);
        rhoBoot = peeler.runPasses(ones, 1).rho;
    }
    auto grid = buildGrid(g, eta, gm, rhoBoot);
    agg.grid = grid.size();
    vector<long double> Dv = computeDelta(g);
    const bool adaptive = gm.deltaBox && gm.adaptiveBox;
    AdaptiveBox box;
    if (adaptive) box.init(grid, Dv, gridStep(eta, gm.order));
    if (rhoFixed < 0 && bestFirst) {
        vector<int> perm; vector<double> key, scratch;
        sortGridBestFirst(grid, perm, key, adaptive ? &box : nullptr, &scratch);
    }
    if (agg.traceOn) agg.trace.reserve(grid.size());
    double rhoHat = (rhoFixed >= 0 ? rhoFixed
                                   : std::max(0.0, rhoBoot));
    const double safety = 1.0;
    double alphaF = 1.0;
    {
        long double Dmax = 1.0L; for (long double d : Dv) Dmax = std::max(Dmax, d);
        if (rhoBoot > 0) alphaF = std::max(1.0, (double)((long double)g.i * Dmax / rhoBoot));
    }

    vector<vector<char>> active(g.i);
    vector<vector<long double>> cL(g.i), cR(g.i);
    for (int j = 0; j < g.i; ++j) {
        active[j].assign(g.n[j], 1);
        cL[j].assign(g.n[j], 0.0L);
        cR[j].assign(g.n[j], 0.0L);
    }
    LayeredGraph rg;
    vector<vector<int>> backmap, fwd;
    vector<vector<int>> survScratch;

    for (size_t s = 0; s < grid.size(); ++s) {
        if (deadline::hit()) break;
        const auto& M = grid[s];
        if (adaptive && rhoHat > 0 && !box.keep(s, rhoHat)) {
            agg.boxSkipped++;
            if (agg.traceOn) agg.trace.push_back({ s, rhoHat, false, -1.0 });
            continue;
        }
        for (int j = 0; j < g.i; ++j)
            std::fill(active[j].begin(), active[j].end(), (char)1);
        ReduceStats st = coReduce(g, M, rhoHat, safety, cL, cR, active);

        agg.Ebefore += st.Ebefore; agg.Fbefore += st.Fbefore;
        double vf = st.Vbefore ? 1.0 - (double)st.Vafter / st.Vbefore : 0.0;
        double ef = st.Ebefore ? 1.0 - (double)st.Eafter / st.Ebefore : 0.0;
        double ff = st.Fbefore > 0 ? 1.0 - (double)(st.Fafter / st.Fbefore) : 0.0;
        agg.sumVf += vf; agg.sumEf += ef; agg.sumFf += ff;
        agg.maxVf = std::max(agg.maxVf, vf); agg.maxEf = std::max(agg.maxEf, ef);
        agg.maxFf = std::max(agg.maxFf, ff);
        agg.rounds_max = std::max(agg.rounds_max, st.rounds); agg.rounds_sum += st.rounds;

        if (st.emptied) {
            agg.eliminated++;
            if (agg.traceOn) agg.trace.push_back({ s, rhoHat, false, -1.0 });
            continue;
        }
        agg.Eafter += st.Eafter; agg.Fafter += st.Fafter;
        agg.solved++;

        buildReducedInto(g, active, rg, backmap, fwd);
        peeler.init(rg);
        PeelResult pr;
        if (passes > 0) { pr = peeler.runPasses(M, passes); }
        else {
            const double den = std::max({rhoHat, rhoBoot, 1e-9});
            const double alphaCell = std::max(1.0,
                (double)((long double)g.i * st.Dafter / den));
            const int Tcap = uniformPassBudget(g.i, (double)st.Dafter, (long)st.Vafter,
                                               Peeler::delta, gIterConst);
            int used = 0;
            if (gCertPeel)
                pr = peeler.runPassesCertPeel(M, Peeler::delta, Tcap,
                                              std::max(best.rho, rhoHat), eta, &used);
            else if (gCertify)
                pr = peeler.runPassesFW(M, Peeler::delta, Tcap,
                                        std::max(best.rho, rhoHat), eta, &used);
            else if (gSafe)
                pr = peeler.runPassesAdaptive(M, Peeler::delta, Tcap,
                                              std::max(best.rho, rhoHat), eta, &used);
            else
                pr = peeler.runPassesPrimal(M, gPatience, Tcap, &used);
            agg.sumAlphaFull += alphaF; agg.sumAlphaCell += alphaCell;
            agg.sumFfull += st.Fbefore; agg.sumFprun += st.Fafter;
            agg.sumTcap += Tcap; agg.sumUsed += used; agg.iterCells++;
        }

        if (pr.rho > best.rho) {
            best.rho = pr.rho; best.M = M;
            peeler.extractSurvivors(survScratch);
            if ((int)best.survivors.size() < g.i) best.survivors.resize(g.i);
            for (int j = 0; j < g.i; ++j) {
                best.survivors[j].clear();
                for (int u : survScratch[j]) best.survivors[j].push_back(backmap[j][u]);
            }
        }
        if (rhoFixed < 0) rhoHat = std::max(rhoHat, pr.rho);
        if (agg.traceOn) agg.trace.push_back({ s, rhoHat, true, pr.rho });
    }
    agg.peelDecCalls = peeler.nDecCalls_;
    agg.decKeyCalls = peeler.nDecCalls_; agg.bktMoves = peeler.nBktMoves_;
    best.gridSize = grid.size();
    return best;
}
