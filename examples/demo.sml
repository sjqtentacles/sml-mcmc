(* demo.sml

   A tour of `sml-mcmc`: run a random-walk Metropolis-Hastings sampler on a
   standard normal target, then a 2-D Gibbs sampler on a correlated bivariate
   normal, both driven entirely by a *seeded* `sml-prng` generator (no ambient
   randomness). Chain diagnostics (mean / variance) come from the vendored
   `sml-stats`. The output is byte-identical across MLton and Poly/ML
   (fixed-decimal formatting, deterministic seeded sampling).

   This is MCMC, not plain Monte Carlo: we build correlated Markov chains whose
   stationary distribution is the target, then summarize them after burn-in.

   Build and run with `make example`. *)

structure M = Mcmc

(* Real formatting that is byte-identical across compilers (fixed decimals;
   always includes a decimal point). *)
fun fmt k x = Real.fmt (StringCvt.FIX (SOME k)) x
fun line s = print (s ^ "\n")

val () = line "=== sml-mcmc demo ============================================="
val () = line ""
val () = line "All randomness is a SEEDED sml-prng state, threaded explicitly."
val () = line ""

(* ---- 1. Metropolis-Hastings on a standard normal N(0,1) ---- *)

(* logTarget is the log of the (unnormalized) density: N(0,1) -> -x^2/2. *)
fun stdNormalLogTarget xs =
  case xs of [x] => ~0.5 * x * x | _ => raise Fail "1-D target"

val mhSeed = SplitMix64.seed 0wx5EED1234
val mhSteps = 50000
val mhBurn = 5000

val ({ chain = mhChain, accepted = mhAccepted }, _) =
  M.metropolisWithAccepts
    { logTarget = stdNormalLogTarget
    , proposalStdev = 2.4
    , init = [0.0]
    , steps = mhSteps } mhSeed

val mhPost = M.dropBurnIn mhBurn mhChain

val () = line "Metropolis-Hastings: target = standard normal N(0,1)"
val () = line ("  proposal stddev = 2.4, steps = " ^ Int.toString mhSteps
               ^ ", burn-in = " ^ Int.toString mhBurn)
val () = line ("  acceptance rate = "
               ^ fmt 4 (M.acceptanceRate { accepted = mhAccepted, steps = mhSteps }))
val () = line ("  empirical mean      = " ^ fmt 4 (M.mean 0 mhPost) ^ "   (target 0.0000)")
val () = line ("  empirical variance  = " ^ fmt 4 (M.variance 0 mhPost) ^ "   (target 1.0000)")
val () = line ""

(* ---- 2. Gibbs sampling on a bivariate normal ---- *)

(* means (1, -2), marginal stddevs (1, 2), correlation 0.5. The full
   conditional of each coordinate of a bivariate normal is univariate normal,
   drawn here with the vendored Stats.Normal.sample (Box-Muller). *)
val m1 = 1.0 and m2 = ~2.0
val s1 = 1.0 and s2 = 2.0
val rho = 0.5
val cond = Math.sqrt (1.0 - rho * rho)

fun cond0 point s =
  let
    val x1 = List.nth (point, 1)
    val mu = m1 + rho * (s1 / s2) * (x1 - m2)
  in Stats.Normal.sample { mu = mu, sigma = s1 * cond } s end

fun cond1 point s =
  let
    val x0 = List.nth (point, 0)
    val mu = m2 + rho * (s2 / s1) * (x0 - m1)
  in Stats.Normal.sample { mu = mu, sigma = s2 * cond } s end

val conditionals = Vector.fromList [cond0, cond1]

val gSeed = SplitMix64.seed 0wx0CC4117
val gSteps = 50000
val gBurn = 5000

val (gChain, _) =
  M.gibbs { conditionals = conditionals, init = [0.0, 0.0], steps = gSteps } gSeed

val gPost = M.dropBurnIn gBurn gChain

val () = line "Gibbs: target = bivariate normal, means (1, -2), stddevs (1, 2), rho 0.5"
val () = line ("  sweeps = " ^ Int.toString gSteps ^ ", burn-in = " ^ Int.toString gBurn)
val () = line ("  coord 0: mean = " ^ fmt 4 (M.mean 0 gPost) ^ "  (target 1.0000)"
               ^ "   variance = " ^ fmt 4 (M.variance 0 gPost) ^ "  (target 1.0000)")
val () = line ("  coord 1: mean = " ^ fmt 4 (M.mean 1 gPost) ^ "  (target -2.0000)"
               ^ "   variance = " ^ fmt 4 (M.variance 1 gPost) ^ "  (target 4.0000)")
val () = line ""

(* ---- 3. Determinism: same seed -> identical chain ---- *)

val (again, _) =
  M.metropolis
    { logTarget = stdNormalLogTarget, proposalStdev = 2.4, init = [0.0], steps = mhSteps }
    (SplitMix64.seed 0wx5EED1234)

fun chainsEqual (xs, ys) =
  ListPair.allEq
    (fn (a, b) => ListPair.allEq (fn (x, y) => Real.== (x, y)) (a, b))
    (xs, ys)

val () = line ("Determinism: re-running MH from the same seed reproduces the chain: "
               ^ (if chainsEqual (mhChain, again) then "yes" else "no"))
val () = line ""
val () = line "==============================================================="
