
#include "driver.hpp"

static long long countInstancesExplicit(const LayeredGraph& g,
                                        const vector<vector<char>>& in) {
    vector<long double> cur(g.n[0]);
    for (int u = 0; u < g.n[0]; ++u) cur[u] = in[0][u] ? 1 : 0;
    for (int j = 1; j < g.i; ++j) {
        vector<long double> nxt(g.n[j], 0);
        for (int u = 0; u < g.n[j]; ++u)
            if (in[j][u])
                for (int w : g.outL(j, u)) nxt[u] += cur[w];
        cur.swap(nxt);
    }
    long double s = 0;
    for (auto v : cur) s += v;
    return (long long)std::llround((double)s);
}

static double bruteForceOptimum(const LayeredGraph& g) {
    for (int j = 0; j < g.i; ++j) assert(g.n[j] <= 12);
    double best = -1.0;
    vector<int> masks(g.i);
    std::function<void(int)> rec = [&](int j) {
        if (j == g.i) {
            vector<vector<char>> in(g.i);
            long double gm = 1.0L;
            for (int t = 0; t < g.i; ++t) {
                in[t].assign(g.n[t], 0);
                int c = 0;
                for (int u = 0; u < g.n[t]; ++u)
                    if (masks[t] >> u & 1) { in[t][u] = 1; ++c; }
                gm *= c;
            }
            long long F = countInstancesExplicit(g, in);
            best = std::max(best, F / std::pow((double)gm, 1.0 / g.i));
            return;
        }
        for (int m = 1; m < (1 << g.n[j]); ++m) { masks[j] = m; rec(j + 1); }
    };
    rec(0);
    return best;
}

static void runDemo() {
    {
        LayeredGraph g;
        g.i = 3; g.n = { 3, 3, 3 };
        g.beginEdges();
        auto E = [&](int j, int u, int v) { g.addEdge(j, u, v); };
        E(0,0,0); E(0,1,0); E(0,1,1); E(0,1,2); E(0,2,0); E(0,2,1);
        E(1,0,0); E(1,0,1); E(1,0,2); E(1,1,0); E(1,1,1); E(1,2,2);
        g.finalize();
        const double exact = bruteForceOptimum(g);
        for (double eta : { 0.5, 0.2, 0.05 }) {
            GridResult r = gridAlgorithm(g, eta);
            std::printf("[demo i=3] eta=%-5.2f rho=%.4f  exact=%.4f  "
                        "(%zu peels, guarantee %.3f*rho*)\n",
                        eta, r.rho, exact, r.gridSize,
                        1.0 / (3.0 * (1 + eta) * (1 + eta)));
            assert(r.rho <= exact + 1e-9);
            assert(r.rho >= exact /
                   (3.0 * (1 + eta) * (1 + eta) * (1 + Peeler::delta)) - 1e-9);
        }
    }
    {
        std::mt19937 rng(7);
        LayeredGraph g;
        g.i = 4; g.n = { 14, 8, 5, 10 };
        g.beginEdges();
        const double p[3] = { 0.25, 0.4, 0.3 };
        std::uniform_real_distribution<double> U(0, 1);
        for (int j = 0; j < 3; ++j)
            for (int u = 0; u < g.n[j]; ++u)
                for (int v = 0; v < g.n[j + 1]; ++v)
                    if (U(rng) < p[j]) g.addEdge(j, u, v);
        g.finalize();
        GridResult r = gridAlgorithm(g, 0.5);
        std::printf("[demo i=4] eta=0.50  rho=%.4f  sizes=(", r.rho);
        for (int j = 0; j < 4; ++j)
            std::printf("%zu%s", r.survivors[j].size(), j < 3 ? "," : ")");
        long long profiles = std::ranges::fold_left(g.n, 1LL, std::multiplies{});
        std::printf("  (%zu peels; full profile enumeration = %lld)\n",
                    r.gridSize, profiles);
        vector<vector<char>> in(4);
        long double gm = 1.0L;
        for (int j = 0; j < 4; ++j) {
            in[j].assign(g.n[j], 0);
            for (int u : r.survivors[j]) in[j][u] = 1;
            gm *= r.survivors[j].size();
        }
        const long long F = countInstancesExplicit(g, in);
        const double rho2 = F / std::pow((double)gm, 0.25);
        assert(std::abs(rho2 - r.rho) < 1e-6 * (1 + r.rho));
        std::printf("[demo i=4] survivor-family recount OK: |F|=%lld, "
                    "rho=%.4f\n", F, rho2);
    }
    {
        vector<double> mStar = { 0.3257, 1.0856, 4.3423, 0.6514 };
        auto mp = roundAndProject(mStar, 0.1);
        double prod = 1;
        for (size_t j = 0; j < mp.size(); ++j) {
            prod *= mp[j];
            std::printf("[round&project] m%zu* = %.4f -> %.5f  ratio %.4f\n",
                        j + 1, mStar[j], mp[j], mp[j] / mStar[j]);
            assert(mp[j] / mStar[j] < 1.21 && mp[j] / mStar[j] > 1 / 1.21);
        }
        std::printf("[round&project] product = %.12f\n", prod);
        assert(std::abs(prod - 1.0) < 1e-9);
    }
    {
        LayeredGraph g; g.i = 3; g.n = { 4, 4, 4 };
        g.beginEdges();
        auto E = [&](int j, int u, int v) { g.addEdge(j, u, v); };
        E(0,0,1); E(0,1,3); E(0,3,0); E(0,3,2); E(0,3,3);
        E(1,0,0); E(1,0,2); E(1,0,3); E(1,1,0); E(1,1,3);
        E(1,2,1); E(1,2,2); E(1,3,0); E(1,3,3);
        g.finalize();
        const double exact = bruteForceOptimum(g);
        std::printf("[near-opt i=3] exact rho* = %.4f; ratio by pass count "
                    "(eta=0.3):\n", exact);
        for (int T : {1, 2, 4, 8, 16}) {
            GridResult r = gridAlgorithm(g, 0.3, 0, T);
            std::printf("    T=%-3d  rho=%.4f  ratio=%.4f\n",
                        T, r.rho, r.rho / exact);
            assert(r.rho <= exact + 1e-7);
        }
    }
}

int main(int argc, char** argv) {
    bool demo = false;
    std::string path;
    double eta = 0.2;
    size_t maxPeels = 0;
    int passes = 1;
    int numerics = 0;
    int mode = 0;
    bool quiet = false;
    double timeout = 100.0;
    bool traceOn = false;
    bool bestFirst = true;
    bool showAlloc = false;
    GridMode gm;
    for (int a = 1; a < argc; ++a) {
        std::string s = argv[a];
        if (s == "--demo") demo = true;
        else if (s == "--prune") mode = 1;
        else if (s == "--prune-oracle") mode = 2;
        else if (s == "--quiet") quiet = true;
        else if (s == "--deltabox") gm.deltaBox = true;
        else if (s == "--no-deltabox") gm.deltaBox = false;
        else if (s == "--adaptive-box") gm.adaptiveBox = true;
        else if (s == "--no-adaptive-box") gm.adaptiveBox = false;
        else if (s == "--trace") traceOn = true;
        else if (s == "--alloc") showAlloc = true;
        else if (s == "--order=arbitrary") bestFirst = false;
        else if (s == "--order=bestfirst") bestFirst = true;
        else if (s.rfind("--grid=", 0) == 0) gm.order = std::stoi(s.substr(7));
        else if (s.rfind("--passes=", 0) == 0) passes = std::stoi(s.substr(9));
        else if (s.rfind("--timeout=", 0) == 0) timeout = std::stod(s.substr(10));
        else if (s.rfind("--delta=", 0) == 0) Peeler::delta = std::stod(s.substr(8));
        else if (s.rfind("--iterc=", 0) == 0) gIterConst = std::stod(s.substr(8));
        else if (s.rfind("--tmax=", 0) == 0) kAdaptiveTmax = std::stoi(s.substr(7));
        else if (s.rfind("--patience=", 0) == 0) gPatience = std::stoi(s.substr(11));
        else if (s == "--safe") gSafe = true;
        else if (s == "--certify") gCertify = true;
        else if (s == "--certify-shared") { gCertify = true; gCertifyShared = true; }
        else if (s == "--certify-peel") { gCertify = true; gCertPeel = true; }
        else if (s.rfind("--share-eps=", 0) == 0) { gShareEps = std::stod(s.substr(12)); gCertify = true; gCertifyShared = true; }
        else if (path.empty() && s.find_first_not_of("0123456789.") != std::string::npos)
            path = s;
        else if (numerics++ == 0) eta = std::stod(s);
        else maxPeels = (size_t)std::stoull(s);
    }
    if (demo) { runDemo(); if (path.empty()) return 0; }
    if (path.empty()) {
        std::cerr << "usage: dpps_grid <edge-list-file> [eta] [maxpeels] [--passes=T]\n"
                  << "                 [--grid=1|2] [--deltabox|--no-deltabox]\n"
                  << "                 [--prune | --prune-oracle] [--quiet] [--demo]\n"
                  << "  --passes=T     : T>0 fixed passes (T=1 is the 1/i peel); T<=0 near-optimal,\n"
                  << "                   capped at the uniform budget T_max=ceil(C*i*Delta*ln|U|/delta^2)\n"
                  << "  --iterc=C      : constant C in the uniform failsafe cap (default 0.1; not tight)\n"
                  << "  --patience=k   : near-opt primal stop: k sub-delta passes before halting (default 2)\n"
                  << "  --safe         : near-opt uses the SOUND dual-gap certificate stop\n"
                  << "                   (provably (1-delta); the certificate normally fires in a pass or two)\n"
                  << "  --certify      : near-opt CERTIFIED mode -- min-norm/pairwise Frank-Wolfe engine\n"
                  << "                   (linear-rate dual) + strong-incumbent pre-sweep, so only the\n"
                  << "                   winning cell is certified to (1-delta); far fewer peels than --safe\n"
                  << "  --certify-shared: --certify + density-continuity across cells: reuse each certified\n"
                  << "                   cell's min-norm point to bound & WARM-START the rest (B(f) is\n"
                  << "                   weight-free).  ~50x faster on clustered nets (Douban); exact at eps=0\n"
                  << "  --share-eps=E  : shared-cert margin (default 0=exact); skip cells within (1+E) of the\n"
                  << "                   incumbent -- guarantee (1-delta)/((1+eta)(1+E)), more pruning\n"
                  << "  --grid=2       : second-order sqrt(eta) lattice (default); --grid=1 first-order\n"
                  << "  --deltabox     : Delta-box localization (default on); --no-deltabox to disable\n"
                  << "  --adaptive-box : re-tighten the Delta-box to the growing incumbent and skip cells\n"
                  << "                   that leave it (Sec 5.4; default on); --no-adaptive-box to disable\n"
                  << "  --prune        : co-reduction with an online (growing) incumbent\n"
                  << "  --prune-oracle : measure reduction at the baseline optimum rho*\n"
                  << "  --timeout=S    : algorithm-only cutoff in seconds (default 100; 0 = off)\n"
                  << "  --trace        : emit the online incumbent trajectory (Fig C(c), prune modes)\n"
                  << "  --order=arbitrary|bestfirst : online cell order (default bestfirst)\n"
                  << "  --alloc        : report heap-allocation counts (instance-free evidence)\n";
        return path.empty() && !demo ? 1 : 0;
    }
    std::expected<LayeredGraph, std::string> parsed;
    std::error_code fsec;
    if (std::filesystem::is_directory(path, fsec)) {
        parsed = readQueryDir(path);
    } else {
        std::ifstream f(path);
        if (!f) { std::cerr << "cannot open " << path << "\n"; return 1; }
        parsed = readTypedEdgeList(f);
    }
    if (!parsed) { std::cerr << parsed.error() << "\n"; return 1; }
    LayeredGraph g = std::move(*parsed);

    {
        const long double Ftot = totalInstances(g);
        long double enumSets = 1.0L;
        for (int j = 0; j < g.i; ++j) enumSets *= (long double)g.n[j];
        std::printf("instances |F| = %.0Lf\n", Ftot);
        std::printf("enum_sets (prod n_j) = %.0Lf\n", enumSets);
        std::fflush(stdout);
    }

    using AClock = std::chrono::steady_clock;
    deadline::arm(timeout);
    double algo_s = 0.0;
    unsigned long long algo_allocs = 0, algo_bytes = 0;
    auto timeCall = [&](auto&& fn) {
        const auto a0 = allocstat::allocs(), b0 = allocstat::nbytes();
        auto t0 = AClock::now();
        auto res = fn();
        algo_s += std::chrono::duration<double>(AClock::now() - t0).count();
        algo_allocs += allocstat::allocs() - a0;
        algo_bytes  += allocstat::nbytes() - b0;
        return res;
    };
    auto reportTime = [&]() {
        if (deadline::timedOut)
            std::printf("algo_time = %.6f s  [TIMEOUT: cutoff %gs, best-so-far]\n",
                        algo_s, timeout);
        else
            std::printf("algo_time = %.6f s\n", algo_s);
    };
    auto emitAlloc = [&]() {
        if (showAlloc)
            std::printf("alloc: algorithm %llu allocations (%llu bytes); "
                        "process total %llu allocations\n",
                        algo_allocs, algo_bytes, allocstat::allocs());
    };
    auto emitPeel = [&](long long dec, long long bkt) {
        if (std::getenv("DPPS_STATS"))
            std::fprintf(stderr, "updates: %lld decreaseKey calls, %lld bucket moves\n",
                         dec, bkt);
    };
    auto emitTrace = [&](const PruneAgg& agg) {
        if (!traceOn) return;
        std::fprintf(stderr, "# trajectory (order=%s): step rhoHat solved cellRho\n",
                     bestFirst ? "bestfirst" : "arbitrary");
        for (const TraceRec& t : agg.trace)
            std::fprintf(stderr, "TRACE %zu %.6f %d %.6f\n",
                         t.step, t.rhoHat, t.solved ? 1 : 0, t.cellRho);
    };

    auto printStats = [&](const PruneAgg& agg, const char* tag) {
        const double pct = 100.0;
        std::printf("\n--- co-reduction effectiveness (%s) ---\n", tag);
        std::printf("grid points              : %zu\n", agg.grid);
        std::printf("  solved (peel run)      : %zu\n", agg.solved);
        std::printf("  ELIMINATED (no solve)  : %zu  (%.1f%% of grid)\n",
                    agg.eliminated, pct * (double)agg.eliminated / std::max<size_t>(1, agg.grid));
        std::printf("  SKIPPED (adaptive box) : %zu  (%.1f%% of grid)\n",
                    agg.boxSkipped, pct * (double)agg.boxSkipped / std::max<size_t>(1, agg.grid));
        std::printf("per-point removal (mean / max):\n");
        std::printf("  vertices               : %.1f%% / %.1f%%\n",
                    pct * agg.sumVf / std::max<size_t>(1, agg.grid), pct * agg.maxVf);
        std::printf("  edges                  : %.1f%% / %.1f%%\n",
                    pct * agg.sumEf / std::max<size_t>(1, agg.grid), pct * agg.maxEf);
        std::printf("  instances |F|          : %.1f%% / %.1f%%\n",
                    pct * agg.sumFf / std::max<size_t>(1, agg.grid), pct * agg.maxFf);
        std::printf("work avoided (summed over grid):\n");
        std::printf("  edges touched          : %lld -> %lld  (%.1f%% removed)\n",
                    agg.Ebefore, agg.Eafter,
                    agg.Ebefore ? pct * (1.0 - (double)agg.Eafter / agg.Ebefore) : 0.0);
        std::printf("  instances handled      : %.0Lf -> %.0Lf  (%.1f%% removed)\n",
                    agg.Fbefore, agg.Fafter,
                    agg.Fbefore > 0 ? pct * (double)(1.0 - agg.Fafter / agg.Fbefore) : 0.0);
        std::printf("  peeling decreaseKey    : %lld calls (pruned solves only)\n", agg.peelDecCalls);
        std::printf("  reduction rounds       : %lld total, %d max per point\n",
                    agg.rounds_sum, agg.rounds_max);
        if (agg.iterCells) {
            const double n = (double)agg.iterCells;
            std::printf("iteration count (alpha_f measured after pruning):\n");
            std::printf("  alpha_f (Delta/rho*)   : %.2f -> %.2f  (tight; mean over %lld solves)\n",
                        agg.sumAlphaFull / n, agg.sumAlphaCell / n, agg.iterCells);
            const long double Ff = agg.sumFfull / (long double)n, Fp = agg.sumFprun / (long double)n;
            std::printf("  alpha_f (|F| diameter) : %.4Lg -> %.4Lg  (loose; %.1f%% smaller)\n",
                        Ff, Fp, Ff > 0 ? 100.0 * (double)(1.0L - Fp / Ff) : 0.0);
            std::printf("  passes   T_cap / used  : %.2f / %.2f  (mean)\n",
                        (double)agg.sumTcap / n, (double)agg.sumUsed / n);
        }
    };

    if (mode == 0) {
        GridResult r = timeCall([&]{ return gridAlgorithm(g, eta, maxPeels, passes, gm); });
        std::printf("eta=%g  grid=%zu peels  [grid=%d deltabox=%d adaptive=%d]\nbest density rho = %.6f\n",
                    eta, r.gridSize, gm.order, (int)gm.deltaBox, (int)gm.adaptiveBox, r.rho);
        if (gm.deltaBox && gm.adaptiveBox)
            std::printf("adaptive Delta-box: skipped %zu of %zu cells (solved %zu)\n",
                        r.boxSkipped, r.gridSize, r.gridSize - r.boxSkipped);
        const bool haveFamily = r.survivors.size() == (size_t)g.i;
        std::printf("part sizes: (");
        for (int j = 0; j < g.i; ++j)
            std::printf("%zu%s", haveFamily ? r.survivors[j].size() : (size_t)0,
                        j + 1 < g.i ? "," : ")\n");
        if (!quiet && haveFamily)
            for (int j = 0; j < g.i; ++j) {
                std::printf("type %d:", j);
                for (int u : r.survivors[j]) std::printf(" %d", u);
                std::printf("\n");
            }
        reportTime();
        emitAlloc();
        emitPeel(r.decKeyCalls, r.bktMoves);
        return 0;
    }

    if (mode == 1) {
        PruneAgg agg;
        agg.traceOn = traceOn;
        GridResult r = timeCall([&]{
            return gridAlgorithmPruned(g, eta, passes, -1.0, agg, gm, bestFirst);
        });
        std::printf("eta=%g  grid=%zu peels  [ONLINE co-reduction]\n", eta, agg.grid);
        std::printf("best density rho = %.6f\n", r.rho);
        std::printf("part sizes: (");
        for (int j = 0; j < g.i; ++j)
            std::printf("%zu%s", r.survivors.empty() ? 0 : r.survivors[j].size(),
                        j + 1 < g.i ? "," : ")\n");
        printStats(agg, "online");
        reportTime();
        emitTrace(agg);
        emitAlloc();
        emitPeel(agg.decKeyCalls, agg.bktMoves);
        return 0;
    }

    GridResult base = timeCall([&]{ return gridAlgorithm(g, eta, maxPeels, passes, gm); });
    std::printf("baseline rho* = %.6f  (grid=%zu) [grid=%d deltabox=%d]\n", base.rho, base.gridSize, gm.order, (int)gm.deltaBox);
    PruneAgg agg;
    agg.traceOn = traceOn;
    GridResult r = timeCall([&]{
        return gridAlgorithmPruned(g, eta, passes, base.rho * (1.0 - 1e-7), agg, gm, bestFirst);
    });
    const double rel = std::abs(r.rho - base.rho) / std::max(1.0, base.rho);
    const double guarantee = 1.0 / ((1.0 + eta) * (1.0 + Peeler::delta));
    const bool ok = (r.rho >= base.rho * guarantee - 1e-9) && (rel < Peeler::delta);
    std::printf("pruned-at-rho* recovered rho = %.6f\n", r.rho);
    std::printf("  rel. diff to baseline = %.2e  (peeler tol delta=%.2f) -> %s\n",
                rel, Peeler::delta,
                ok ? "within solver tolerance, optimum preserved" : "OUT OF TOLERANCE!");
    std::printf("  >= guarantee floor %.4f*rho* ? %s\n",
                guarantee, r.rho >= base.rho * guarantee - 1e-9 ? "yes" : "NO");
    printStats(agg, "oracle, threshold = rho*");
    reportTime();
    emitTrace(agg);
    emitAlloc();
    emitPeel(agg.decKeyCalls, agg.bktMoves);
    return 0;
}
