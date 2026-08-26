#pragma once
#include "dpss_common.hpp"
struct GridResult {
    double rho = -1.0;
    vector<vector<int>> survivors;
    vector<double> M;
    size_t gridSize = 0;
    size_t boxSkipped = 0;   // cells dropped by the adaptive Delta-box before solving
    long long decKeyCalls = 0, bktMoves = 0;
};

static int kAdaptiveTmax = 1024;

static int gPatience = 2;
static bool gSafe = false;
static bool gCertify = false;
static bool gCertifyShared = false;
static bool gCertPeel = false;
static double gShareEps = 0.0;

static double gIterConst = 0.1;

static inline int uniformPassBudget(int i, double Delta, long U, double delta, double c) {
    const double lnU = std::log((double)std::max<long>(U, 2));
    const double d = std::max(delta, 1e-9);
    long t = (long)std::ceil(c * (double)i * std::max(Delta, 1.0) * lnU / (d * d));
    return (int)std::min<long>(std::max<long>(t, 1), kAdaptiveTmax);
}
