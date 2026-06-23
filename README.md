# sml-mcmc

Markov-Chain Monte Carlo samplers in pure Standard ML — a random-walk
**Metropolis-Hastings** sampler and a **Gibbs** sampler, plus chain diagnostics
(mean / variance / acceptance rate) — built on top of
[`sml-prng`](https://github.com/sjqtentacles/sml-prng) for *all* randomness and
[`sml-stats`](https://github.com/sjqtentacles/sml-stats) for the summary
statistics. No FFI, no IO, no wall-clock, **no ambient randomness**, and
**deterministic**, byte-identically under both [MLton](http://mlton.org/) and
[Poly/ML](https://www.polyml.org/).

[![CI](https://github.com/sjqtentacles/sml-mcmc/actions/workflows/ci.yml/badge.svg)](https://github.com/sjqtentacles/sml-mcmc/actions/workflows/ci.yml)

## MCMC, not plain Monte Carlo

This library is **MCMC** and is deliberately distinct from the sibling
`sml-montecarlo`:

- **Plain Monte Carlo** (`sml-montecarlo`) draws *independent* samples to
  estimate an integral or expectation directly — e.g. averaging `f(x)` over
  uniform draws to approximate `∫ f`.
- **MCMC** (this library) builds a *Markov chain* whose **stationary
  distribution is the target**. Each sample depends on the previous one; the
  chain itself is what you analyse, after discarding an initial **burn-in**.
  Use MCMC when you can only evaluate a target density up to a normalizing
  constant (you supply `logTarget`, an *unnormalized* log-density) and cannot
  sample it directly.

So `sml-montecarlo` answers "what is this integral?"; `sml-mcmc` answers "give
me (correlated) draws from this hard-to-sample distribution."

## Status

- 41 assertions, green on MLton and Poly/ML.
- Basis-library + vendored `sml-prng` + vendored `sml-stats` only;
  deterministic across compilers.
- Vendors **both** dependencies (Layout B) under `lib/github.com/sjqtentacles/`,
  so the repo builds standalone.

## Purity

No FFI, no IO inside the library, no wall-clock, **no ambient randomness**, and
no threads. Randomness comes *only* from an explicitly threaded `sml-prng`
generator state: every sampler takes a state and returns the chain together
with the successor state, so the same seed always reproduces the **exact same
chain** across runs, machines, and compilers. All Gaussian proposals use a
Box-Muller transform over `sml-prng` uniforms (byte-identical across
compilers). All chain statistics are `real`, so tests compare them through an
explicit tolerance (`Support.approxTol`), never string or structural equality;
the determinism tests assert *bit-for-bit* identical chains across two runs.

## Install

With [`smlpkg`](https://github.com/diku-dk/smlpkg):

```
smlpkg add github.com/sjqtentacles/sml-mcmc
smlpkg sync
```

Include the MLB from your own (it pulls in the vendored `sml-prng` and
`sml-stats`):

```
local
  $(SML_LIB)/basis/basis.mlb
  lib/github.com/sjqtentacles/sml-mcmc/... (via smlpkg)
in
  ...
end
```

This brings `structure Mcmc` (and the vendored `SplitMix64`/`Stats`) into
scope.

## Quick start

```sml
(* Random-walk Metropolis-Hastings on a standard normal N(0,1):
   logTarget is the UNNORMALIZED log-density, here -x^2/2. *)
fun logTarget xs = case xs of [x] => ~0.5 * x * x | _ => raise Fail "1-D"

val seed = SplitMix64.seed 0wx5EED1234

val (chain, _) =
  Mcmc.metropolis
    { logTarget = logTarget, proposalStdev = 2.4, init = [0.0], steps = 50000 }
    seed

val post = Mcmc.dropBurnIn 5000 chain     (* discard burn-in            *)
val m    = Mcmc.mean 0 post               (* ~ 0.0  (coordinate 0 mean) *)
val v    = Mcmc.variance 0 post           (* ~ 1.0                      *)
```

## API (`signature MCMC`)

```sml
type rng                                  (* underlying sml-prng state *)
exception Mcmc of string

(* random-walk Metropolis-Hastings (symmetric Gaussian proposal) *)
val metropolis :
  { logTarget : real list -> real, proposalStdev : real
  , init : real list, steps : int } -> rng -> real list list * rng

(* as metropolis, but also reports how many proposals were accepted *)
val metropolisWithAccepts :
  { logTarget : real list -> real, proposalStdev : real
  , init : real list, steps : int }
  -> rng -> { chain : real list list, accepted : int } * rng

(* Gibbs sampling: one conditional sampler per coordinate *)
val gibbs :
  { conditionals : (real list -> rng -> real * rng) vector
  , init : real list, steps : int } -> rng -> real list list * rng

(* diagnostics (thin wrappers over the vendored sml-stats) *)
val coordinate     : int -> real list list -> real list
val dropBurnIn     : int -> real list list -> real list list
val mean           : int -> real list list -> real
val variance       : int -> real list list -> real
val acceptanceRate : { accepted : int, steps : int } -> real
```

### Conventions

- **Points and chains.** A point in state space is a `real list` (its length is
  the dimension). A chain is the `real list list` of successive points, in
  order, of length `steps + 1` — the initial point is the head.
- **`logTarget`.** The *log* of the (possibly unnormalized) target density.
  Only differences matter, so any additive constant is fine. Metropolis uses a
  symmetric Gaussian random-walk proposal, so the acceptance probability
  reduces to `min(1, exp(logTarget(proposed) − logTarget(current)))`.
- **Acceptance test.** A uniform `u ~ U[0,1)` is drawn and the move is accepted
  iff `ln u < logTarget(proposed) − logTarget(current)` (always true when the
  log-ratio is ≥ 0). Both the proposal and this test consume the threaded
  generator state, in that order.
- **Gibbs conditionals.** `conditionals` is a `vector` of length `d` (the
  dimension). `conditionals` sub `i` has type `real list -> rng -> real * rng`:
  given the **current full point** (all coordinates) and a generator state, it
  returns a fresh draw for *its own* coordinate `i` and the successor state.
  Each sampler reads whatever other coordinates it needs out of the supplied
  point, so the interface is uniform regardless of the dependency structure.
  Each "step" is a full sweep that updates coordinates `0..d-1` in order, and
  later updates within a sweep see the freshly updated earlier coordinates.
- **Diagnostics.** `coordinate i` projects a chain onto its `i`-th coordinate;
  `mean`/`variance` then call `Stats.mean`/`Stats.variance` (sample variance,
  `/(n−1)`) on that series. `dropBurnIn n` discards the first `n` points.
  `acceptanceRate {accepted, steps}` is `accepted / steps` (raises `Mcmc` if
  `steps ≤ 0`).
- **Determinism.** Identical seed ⇒ identical chain, bit-for-bit, on every run
  and both compilers. There is no other source of randomness.

## Build & test

```
make test        # MLton
make test-poly   # Poly/ML
make all-tests   # both
make example     # build + run examples/demo.sml
make clean
```

Both compilers run the same strict-TDD suite (41 assertions):

- **Metropolis on a standard normal** (`logTarget = −x²/2`): after burn-in the
  empirical mean ≈ 0 and variance ≈ 1 within generous tolerances; chain shape
  (`length = steps + 1`, head = init, every point 1-D) is checked exactly; and
  argument validation (`steps < 0`, empty `init`, non-positive `proposalStdev`)
  raises.
- **Determinism** — the same seed reproduces the **exact** same chain
  (bit-for-bit), a different seed produces a different chain, and continuing
  from the returned successor state is itself reproducible.
- **2-D Gibbs on a bivariate normal** (means `(1, −2)`, stddevs `(1, 2)`,
  correlation `0.5`) — recovers the marginal means and variances within
  generous tolerances, is reproducible under a fixed seed, and validates its
  arguments.
- **Diagnostics** — `coordinate`/`mean`/`variance` agree with `sml-stats`
  directly; `dropBurnIn` trims correctly; `acceptanceRate` arithmetic is exact;
  and the empirical acceptance rate of a random-walk sampler lies strictly in
  `(0, 1)` and is reproducible.

All reals are compared through a tolerance (`Support.approxTol`), never
stringified.

## Vendoring

This library depends on **two** sibling libraries, both vendored verbatim
(minus their tests) under `lib/github.com/sjqtentacles/`:

- [`sml-prng`](https://github.com/sjqtentacles/sml-prng) —
  `lib/github.com/sjqtentacles/sml-prng/` (`prng.sig`, `prng.sml`, `prng.mlb`).
  The sole source of randomness.
- [`sml-stats`](https://github.com/sjqtentacles/sml-stats) —
  `lib/github.com/sjqtentacles/sml-stats/` (`stats.sig`, `stats.sml`,
  `stats.mlb`). The **descriptive-statistics subset** used for diagnostics
  (`mean`, `variance`, `quantile`, the `Normal` distribution). The full
  upstream library also vendors `sml-specfun` for its t/F/chi-square tests;
  that part is intentionally not vendored here, since MCMC diagnostics never
  touch it.

**Dependency order matters.** `sml-stats`'s `StatsFn` functor is taken over a
`sml-prng` `RANDOM` generator (and the default `Stats` is instantiated with
`SplitMix64`), so `sml-prng` must load **before** `sml-stats`, which loads
before the MCMC sources:

- In the MLton build, `src/mcmc.mlb` lists `../lib/.../sml-prng/prng.mlb`, then
  `../lib/.../sml-stats/stats.mlb` (which itself references the same `prng.mlb`
  — MLton elaborates a given `.mlb` path only once, so the explicit listing
  just pins the order), then `mcmc.sig`/`mcmc.sml`.
- In the Poly/ML build, the `make test-poly` `use`-chain loads
  `prng.sig`/`prng.sml`, then `stats.sig`/`stats.sml`, then `mcmc.sig`/
  `mcmc.sml`, then the test driver — strict dependency order.

`sml.pkg` records **both** dependencies in its `require` block so `smlpkg sync`
can refresh them.

## Example

`make example` runs a Metropolis sampler on a standard normal and a Gibbs
sampler on a correlated bivariate normal, all from fixed seeds (output is
byte-identical under MLton and Poly/ML):

```
=== sml-mcmc demo =============================================

All randomness is a SEEDED sml-prng state, threaded explicitly.

Metropolis-Hastings: target = standard normal N(0,1)
  proposal stddev = 2.4, steps = 50000, burn-in = 5000
  acceptance rate = 0.4420
  empirical mean      = ~0.0014   (target 0.0000)
  empirical variance  = 1.0065   (target 1.0000)

Gibbs: target = bivariate normal, means (1, -2), stddevs (1, 2), rho 0.5
  sweeps = 50000, burn-in = 5000
  coord 0: mean = 0.9964  (target 1.0000)   variance = 0.9921  (target 1.0000)
  coord 1: mean = ~2.0044  (target -2.0000)   variance = 3.9513  (target 4.0000)

Determinism: re-running MH from the same seed reproduces the chain: yes

===============================================================
```

### Poly/ML note

CI builds Poly/ML 5.9.1 from source rather than using the Ubuntu package
(Poly/ML 5.7.1), whose X86 code generator crashes (`asGenReg raised while
compiling`) on some code. See `.github/workflows/ci.yml`.

## License

MIT — see [LICENSE](LICENSE).
