# dpps_core — Near-Optimal Densest P-Partite Subgraph Search (DPSS)

C++23 implementation of the cover-and-conquer algorithm from
*"From Enumeration to Covering: Near-Optimal Densest P-Partite Subgraph Search
over Large Heterogeneous Information Networks"* 

**Problem (DPSS).** Given a HIN G = (V, E) and a meta-path P = (A_1, …, A_i) of
length i = |P|, find the P-family V = (V_1, …, V_i), V_j ⊆ V(A_j), maximizing the
density ρ(V) = |F(V)| / g(V), where F(V) is the set of P-instances of the induced
subgraph G(V) and g(V) = (∏_j |V_j|)^(1/i).

**Algorithm.** Instead of enumerating iRM-sets (m_j = g(V)/|V_j|, ∏_j m_j = 1),
cover the iRM polytope Π_n with a multiplicative η-net — a lattice on the i−1
free log-coordinates, last coordinate back-solved so every representative
satisfies the product constraint exactly — and run one fixed-weight cell solve
per representative: near-optimal weighted iterative peeling on the supermodular
reformulation. The cell containing the optimum's iRM-set certifies a
(1−δ)/(1+η) approximation. Peeling is instance-free: the per-vertex incidence
factors as P(v, G(V)) = L(v)·R(v), maintained by forward/backward DP sweeps, so
instances are never materialized.

## Files

Single translation unit: `dpps_grid.cpp` includes everything, so the header split
has zero codegen/perf impact.

| File | Contents | Paper section |
|---|---|---|
| `dpss_common.hpp` | includes, allocation counter (`--alloc`), deadline cutoff | — |
| `layered_graph.hpp` | `LayeredGraph`, typed-edge / query-directory readers, L·R instance DP | §6 instance-free solver |
| `grid.hpp` | η-net builders: first/second-order nets, Δ-box localization | §5 covering |
| `peeler.hpp` | `Peeler`: 1/i peel, near-opt iterative peeling, dual-gap certificate (Frank–Wolfe), certified peeling | §6 cell solver |
| `config.hpp` | flags/globals, uniform pass budget T, `GridResult` | — |
| `coreduce.hpp` | co-reduction: per-type incidence-floor vertex pruning + cell drops | §7 pruning |
| `driver.hpp` | `gridAlgorithm` + `gridAlgorithmPruned` (cover-and-conquer loops) | §4 overview |
| `dpps_grid.cpp` | brute-force self-checks, demo, `main()` | — |

(§ numbers refer to the paper: §4 Solution Overview, §5 From Enumeration to
Covering, §6 Near-Optimal Approximation for i>2, §7 Vertex and Cell Pruning.)

## Build & run

```sh
make            # = g++ -O2 -std=c++23 -o dpps_grid dpps_grid.cpp
make test       # build + smoke run on graph_1k.txt (incl. --demo self-checks)
make debug      # -O0 -g with ASan/UBSan

./dpps_grid test_queries/graph_1k.txt 0.2             # eta = 0.2
./dpps_grid test_queries/graph_1k.txt 0.2 --demo      # also run built-in self-checks
```

Output: best density ρ found, part sizes of the best family, grid size (number
of cell solves), surviving vertex ids per layer.

## Input

`argv[1]` is either of:

1. **Typed edge list** (the files in `test_queries/`):
   ```
   i                  # meta-path length (number of layers)
   n_1 n_2 ... n_i    # vertices per layer (ids 0..n_j-1)
   j u v              # edge: vertex u of layer j — vertex v of layer j+1 (j 0-based)
   ```
2. **Materialized query directory**: `layer<k>_<TYPE>.nodes` (one node id per
   line) + `step<k>_<REL>.edges` ("u v" per line, edge between layers k and
   k+1). `main()` dispatches on `is_directory`; the reader is `readQueryDir`
   in `layered_graph.hpp` (node ids mapped in first-appearance order). This is
   the format of the paper's benchmark query pools, grouped by i.

## Key flags

Flag combinations reproduce the algorithm
variants A2–A8 of the experiments (Table "tab:algs"); A1 is the prior
state-of-the-art enumeration baseline.

```
eta                 positional, after the input path (e.g. ./dpps_grid <input> 0.2):
                    covering budget of the eta-net, cell loss factor 1+eta
--passes=T          T>0 fixed peeling passes (T=1 = the prior 1/i peel); T<=0 near-optimal
--grid=1|2          first- or second-order eta-net (second-order = sqrt(2*eta/e) spacing)
--deltabox / --no-deltabox   Delta-box localization of the net (Lemma "deltabox")
--prune | --prune-oracle     co-reduction: incidence-floor vertex pruning, incumbent-driven
--delta=d           solver tolerance delta (solver factor beta = 1-delta)
--certify           certified mode: min-norm / pairwise Frank-Wolfe engine proving
                    (1-delta)/(1+eta) on every cell
--certify-shared    certified + density-continuity cell pruning (--share-eps=)
--certify-peel      certified iterative peeling
--patience=k  --iterc=C  --tmax=  --timeout=  --quiet  --trace  --alloc
```

Run with no arguments for full usage text.
