(* test_metropolis.sml -- random-walk Metropolis-Hastings on a standard normal.

   The target is the standard normal N(0,1) up to a constant, so
   `logTarget x = -x^2/2` (we read the single coordinate out of the point).
   With a fixed seed and a sensible proposal width, a long chain's empirical
   mean should be ~0 and its variance ~1 after we drop a burn-in prefix.

   Tolerances are deliberately generous: these are sampling estimates from a
   correlated chain, not closed-form values. We also sanity-check the chain's
   shape (length = steps + 1, head = init) which IS exact. *)

structure MetropolisTests =
struct
  open Support
  structure M = Mcmc

  val seed = SplitMix64.seed 0wx5EED1234

  fun stdNormalLogTarget xs =
    case xs of
      [x] => ~0.5 * x * x
    | _   => raise Fail "1-D target expects a singleton point"

  val steps = 20000
  val burn  = 2000

  val (chain, _) =
    M.metropolis
      { logTarget = stdNormalLogTarget
      , proposalStdev = 2.4
      , init = [0.0]
      , steps = steps } seed

  fun run () =
    let
      val () = Harness.section "metropolis: chain shape (exact)"
      val () = Harness.checkInt "chain length = steps + 1"
                 (steps + 1, List.length chain)
      val () = Harness.checkBool "head of chain is init"
                 (true, (case chain of (x :: _) :: _ => Real.== (x, 0.0) | _ => false))
      val () = Harness.check "every point is 1-dimensional"
                 (List.all (fn p => List.length p = 1) chain)

      val post = M.dropBurnIn burn chain
      val m = M.mean 0 post
      val v = M.variance 0 post

      val () = Harness.section "metropolis: standard normal moments after burn-in"
      (* generous tolerances for a correlated sampling estimate *)
      val () = checkApproxTol 0.1 "empirical mean ~ 0" (0.0, m)
      val () = checkApproxTol 0.15 "empirical variance ~ 1" (1.0, v)

      val () = Harness.section "metropolis: argument validation"
      val () = Harness.checkRaises "steps < 0 raises"
                 (fn () => M.metropolis
                    { logTarget = stdNormalLogTarget, proposalStdev = 1.0
                    , init = [0.0], steps = ~1 } seed)
      val () = Harness.checkRaises "empty init raises"
                 (fn () => M.metropolis
                    { logTarget = stdNormalLogTarget, proposalStdev = 1.0
                    , init = [], steps = 10 } seed)
      val () = Harness.checkRaises "non-positive proposalStdev raises"
                 (fn () => M.metropolis
                    { logTarget = stdNormalLogTarget, proposalStdev = 0.0
                    , init = [0.0], steps = 10 } seed)

      val () = Harness.section "metropolis: zero steps -> just the init point"
      val (c0, _) =
        M.metropolis
          { logTarget = stdNormalLogTarget, proposalStdev = 1.0
          , init = [3.0], steps = 0 } seed
      val () = Harness.checkInt "zero-step chain has length 1" (1, List.length c0)
    in
      ()
    end
end
