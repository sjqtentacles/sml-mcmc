(* stats.sig (vendored subset for sml-mcmc)

   The descriptive-statistics core of `sml-stats`: moments, order statistics,
   and the Normal distribution. This is the subset `sml-mcmc` depends on for
   chain diagnostics (mean / variance / quantiles). The full upstream library
   (regression, t/F/chi-square tests, Binomial/Poisson) additionally vendors
   `sml-specfun`; that part is intentionally NOT vendored here, since MCMC
   diagnostics never touch it.

   Everything is built on Basis-library `real`/`int` and the vendored
   `sml-prng` generators, so it is deterministic and behaves identically under
   MLton and Poly/ML. Sampling is pure: every `sample` takes a generator state
   and returns the draw together with the successor state, so callers thread
   the state exactly as they would with `sml-prng`.

   The abstract type `rng` is the random-generator state. The library is a
   functor over a `RANDOM` generator (see `StatsFn`); the default `Stats`
   structure is instantiated with `SplitMix64`, so `rng = SplitMix64.state`.

   Numerical conventions:
     - `variance`/`stddev` are the *sample* statistics (Bessel-corrected,
       divide by n-1); `variancePop`/`stddevPop` are the population forms
       (divide by n).
     - `quantile q xs` uses the linear-interpolation method (R/NumPy "type 7":
       h = (n-1)*q), so `quantile 0.5` agrees with `median`. *)

signature STATS =
sig
  (* Generator state of the underlying `sml-prng` instance. *)
  type rng

  (* Raised by the descriptive statistics on an empty input. *)
  exception Empty

  (* ---- descriptive statistics (input lists must be non-empty) ---- *)
  val sum         : real list -> real
  val mean        : real list -> real
  val variance    : real list -> real   (* sample,     /(n-1) *)
  val variancePop : real list -> real   (* population, /n     *)
  val stddev      : real list -> real   (* sqrt variance      *)
  val stddevPop   : real list -> real   (* sqrt variancePop   *)
  val median      : real list -> real
  (* `quantile q xs`, q in [0,1], linear interpolation (type 7). *)
  val quantile    : real -> real list -> real
  val minimum     : real list -> real
  val maximum     : real list -> real

  (* ---- distributions ---- *)

  structure Normal :
  sig
    type param = { mu : real, sigma : real }     (* sigma > 0 *)
    val pdf    : param -> real -> real
    val cdf    : param -> real -> real
    (* Box-Muller draw. *)
    val sample : param -> rng -> real * rng
  end
end
