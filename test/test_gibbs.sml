(* test_gibbs.sml -- 2-D Gibbs sampler on a bivariate normal.

   Target: bivariate normal with means (m1, m2) = (1.0, ~2.0), marginal
   standard deviations (s1, s2) = (1.0, 2.0), correlation rho = 0.5. The full
   conditional of each coordinate of a bivariate normal is itself univariate
   normal:

     x0 | x1 ~ N( m1 + rho*(s1/s2)*(x1 - m2),  s1*sqrt(1 - rho^2) )
     x1 | x0 ~ N( m2 + rho*(s2/s1)*(x0 - m1),  s2*sqrt(1 - rho^2) )

   Each conditional sampler reads the *other* coordinate out of the supplied
   current point and draws its own coordinate with the vendored
   `Stats.Normal.sample` (Box-Muller over sml-prng). After burn-in the
   empirical marginal means and variances should match (m1, m2) and
   (s1^2, s2^2) within generous tolerances. *)

structure GibbsTests =
struct
  open Support
  structure M = Mcmc

  val m1 = 1.0 and m2 = ~2.0
  val s1 = 1.0 and s2 = 2.0
  val rho = 0.5
  val cond = Math.sqrt (1.0 - rho * rho)

  (* coordinate 0 | coordinate 1 *)
  fun cond0 point s =
    let
      val x1 = List.nth (point, 1)
      val mu = m1 + rho * (s1 / s2) * (x1 - m2)
      val sigma = s1 * cond
    in
      Stats.Normal.sample { mu = mu, sigma = sigma } s
    end

  (* coordinate 1 | coordinate 0 *)
  fun cond1 point s =
    let
      val x0 = List.nth (point, 0)
      val mu = m2 + rho * (s2 / s1) * (x0 - m1)
      val sigma = s2 * cond
    in
      Stats.Normal.sample { mu = mu, sigma = sigma } s
    end

  val conditionals = Vector.fromList [cond0, cond1]

  val seed = SplitMix64.seed 0wx0CC4117
  val steps = 30000
  val burn  = 3000

  val (chain, _) =
    M.gibbs { conditionals = conditionals, init = [0.0, 0.0], steps = steps } seed

  fun run () =
    let
      val () = Harness.section "gibbs: chain shape (exact)"
      val () = Harness.checkInt "chain length = steps + 1"
                 (steps + 1, List.length chain)
      val () = Harness.check "every point is 2-dimensional"
                 (List.all (fn p => List.length p = 2) chain)

      val post = M.dropBurnIn burn chain
      val mean0 = M.mean 0 post and mean1 = M.mean 1 post
      val var0  = M.variance 0 post and var1 = M.variance 1 post

      val () = Harness.section "gibbs: bivariate-normal marginal means after burn-in"
      val () = checkApproxTol 0.1 "mean of coord 0 ~ 1.0" (m1, mean0)
      val () = checkApproxTol 0.1 "mean of coord 1 ~ ~2.0" (m2, mean1)

      val () = Harness.section "gibbs: bivariate-normal marginal variances after burn-in"
      val () = checkApproxTol 0.15 "var of coord 0 ~ 1.0" (s1 * s1, var0)
      val () = checkApproxTol 0.6  "var of coord 1 ~ 4.0" (s2 * s2, var1)

      val () = Harness.section "gibbs: determinism (same seed -> identical chain)"
      val (again, _) =
        M.gibbs { conditionals = conditionals, init = [0.0, 0.0], steps = 1000 }
          (SplitMix64.seed 0wx0CC4117)
      val (again2, _) =
        M.gibbs { conditionals = conditionals, init = [0.0, 0.0], steps = 1000 }
          (SplitMix64.seed 0wx0CC4117)
      val () = Harness.check "two gibbs runs from the same seed agree exactly"
                 (chainEqual (again, again2))

      val () = Harness.section "gibbs: argument validation"
      val () = Harness.checkRaises "init length <> #conditionals raises"
                 (fn () => M.gibbs
                    { conditionals = conditionals, init = [0.0], steps = 10 } seed)
      val () = Harness.checkRaises "negative steps raises"
                 (fn () => M.gibbs
                    { conditionals = conditionals, init = [0.0, 0.0], steps = ~5 } seed)
    in
      ()
    end
end
