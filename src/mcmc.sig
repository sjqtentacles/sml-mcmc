(* mcmc.sig

   Markov-Chain Monte Carlo samplers in pure Standard ML: a random-walk
   Metropolis-Hastings sampler and a Gibbs sampler, plus a handful of chain
   diagnostics.

   This is MCMC, *not* plain Monte Carlo. Plain Monte Carlo (e.g. the sibling
   `sml-montecarlo`) draws independent samples to estimate an integral or
   expectation directly. MCMC instead builds a *Markov chain* whose stationary
   distribution is the target: each sample depends on the previous one, and the
   chain is what you analyse (after discarding an initial "burn-in"). Use MCMC
   when you can evaluate a target density only up to a normalizing constant and
   cannot sample it directly.

   PURITY / RANDOMNESS. Every sampler is pure and deterministic. Randomness
   comes *only* from an explicitly threaded `sml-prng` generator state (`rng`):
   a sampler takes a state and returns the produced chain together with the
   successor state, so the same seed always reproduces the exact same chain on
   every run, machine, and compiler (MLton and Poly/ML alike). There is no
   ambient randomness, no clock, no IO, no FFI, and no threads. All Gaussian
   draws use a Box-Muller transform over `sml-prng` uniforms, byte-identical
   across compilers.

   STATE / TYPES. A "point" in the state space is a `real list` (its length is
   the dimension). A "chain" is the `real list list` of successive points, in
   order, of length `steps + 1` (the initial point is included as the head).
   The abstract type `rng` is the underlying `sml-prng` generator state. The
   library is a functor over a `RANDOM` generator (see `McmcFn`); the default
   `Mcmc` structure is instantiated with `SplitMix64`, so `rng =
   SplitMix64.state`. *)

signature MCMC =
sig
  (* Generator state of the underlying `sml-prng` instance. *)
  type rng

  (* Raised on invalid arguments (e.g. steps < 0, empty init, non-positive
     proposal standard deviation, dimension mismatch). *)
  exception Mcmc of string

  (* ---- Metropolis-Hastings (symmetric Gaussian random walk) ----

     `metropolis {logTarget, proposalStdev, init, steps} rng`

       logTarget     log of the (possibly unnormalized) target density at a
                     point; only differences matter, so any additive constant
                     is fine.
       proposalStdev standard deviation of the symmetric Gaussian random-walk
                     proposal, applied independently to each coordinate
                     (must be > 0).
       init          the starting point (must be non-empty).
       steps         number of proposal steps to take (must be >= 0).

     Returns `(chain, rng')` where `chain` has length `steps + 1`, head = init.
     Because the proposal is symmetric, the acceptance probability reduces to
     min(1, exp(logTarget(proposed) - logTarget(current))). *)
  val metropolis :
    { logTarget : real list -> real
    , proposalStdev : real
    , init : real list
    , steps : int } -> rng -> real list list * rng

  (* As `metropolis`, but also reports how many proposals were accepted, so
     callers can compute an acceptance rate. The accepted count is in
     [0, steps]. *)
  val metropolisWithAccepts :
    { logTarget : real list -> real
    , proposalStdev : real
    , init : real list
    , steps : int } -> rng -> { chain : real list list, accepted : int } * rng

  (* ---- Gibbs sampling ----

     Gibbs sampling updates one coordinate at a time by drawing it from its
     full conditional distribution given the current values of all the other
     coordinates. The caller supplies one conditional sampler per coordinate.

     A conditional sampler has type `real list -> rng -> real * rng`: given the
     *current full point* (all coordinates, including the one being updated)
     and a generator state, it returns a fresh draw for *its own* coordinate
     and the successor state. The sampler for coordinate i is responsible for
     reading whatever other coordinates it needs out of the supplied point;
     this keeps the interface uniform regardless of the dependency structure.

     `gibbs {conditionals, init, steps} rng`

       conditionals  a vector of length d (= dimension); `conditionals` sub i
                     samples coordinate i given the rest.
       init          the starting point (length must equal d, non-empty).
       steps         number of full sweeps (each sweep updates all d
                     coordinates in index order 0..d-1); must be >= 0.

     Returns `(chain, rng')` where `chain` has length `steps + 1` (one point
     recorded after each full sweep, plus the initial point as head). *)
  val gibbs :
    { conditionals : (real list -> rng -> real * rng) vector
    , init : real list
    , steps : int } -> rng -> real list list * rng

  (* ---- diagnostics (thin wrappers over the vendored sml-stats) ----

     A chain is a `real list list`; coordinate i of every point forms a scalar
     series, summarized below. `coordinate i chain` extracts that series.
     `dropBurnIn n chain` discards the first n points. *)

  val coordinate  : int -> real list list -> real list
  val dropBurnIn  : int -> real list list -> real list list

  (* Mean / variance (sample, /(n-1)) of coordinate i across the chain. *)
  val mean        : int -> real list list -> real
  val variance    : int -> real list list -> real

  (* Acceptance rate = accepted / steps, in [0, 1]. Raises `Mcmc` if
     steps <= 0. *)
  val acceptanceRate : { accepted : int, steps : int } -> real
end
