(* test_diagnostics.sml -- chain diagnostics and acceptance rate.

   The diagnostics are thin, faithful wrappers over the vendored `sml-stats`:
   `coordinate i` projects a chain onto its i-th coordinate, `mean`/`variance`
   then call `Stats.mean`/`Stats.variance` on that series, and `dropBurnIn`
   trims a prefix. We pin them against hand-built chains (exact `Stats` values,
   compared to a tight epsilon) and against `Stats` directly.

   The acceptance rate of a random-walk Metropolis sampler on a smooth target
   must lie strictly inside (0, 1): a positive proposal width means some
   proposals are rejected (rate < 1) and, for a reasonable width, many are
   accepted (rate > 0). `acceptanceRate` itself is just accepted/steps. *)

structure DiagnosticsTests =
struct
  open Support
  structure M = Mcmc

  (* A tiny, fully explicit chain: three 2-D points. *)
  val chain = [ [1.0, 10.0], [2.0, 20.0], [3.0, 30.0] ]

  fun logTarget xs = case xs of [x] => ~0.5 * x * x | _ => raise Fail "1-D"

  fun run () =
    let
      val () = Harness.section "diagnostics: coordinate projection"
      val c0 = M.coordinate 0 chain
      val c1 = M.coordinate 1 chain
      val () = checkApprox "coord 0 head" (1.0, List.nth (c0, 0))
      val () = checkApprox "coord 0 last" (3.0, List.nth (c0, 2))
      val () = checkApprox "coord 1 head" (10.0, List.nth (c1, 0))
      val () = Harness.checkInt "coordinate series length = chain length"
                 (3, List.length c0)

      val () = Harness.section "diagnostics: mean/variance agree with sml-stats"
      val () = checkApprox "mean coord 0 = Stats.mean" (Stats.mean c0, M.mean 0 chain)
      val () = checkApprox "variance coord 0 = Stats.variance"
                 (Stats.variance c0, M.variance 0 chain)
      val () = checkApprox "mean coord 0 = 2.0" (2.0, M.mean 0 chain)

      val () = Harness.section "diagnostics: dropBurnIn"
      val trimmed = M.dropBurnIn 1 chain
      val () = Harness.checkInt "dropBurnIn 1 leaves 2 points" (2, List.length trimmed)
      val () = checkApprox "dropBurnIn keeps the tail (coord 0 head = 2.0)"
                 (2.0, List.nth (M.coordinate 0 trimmed, 0))
      val () = Harness.checkInt "dropBurnIn beyond length -> empty"
                 (0, List.length (M.dropBurnIn 99 chain))

      val () = Harness.section "diagnostics: acceptanceRate arithmetic"
      val () = checkApprox "accepted/steps = 0.5" (0.5, M.acceptanceRate { accepted = 5, steps = 10 })
      val () = checkApprox "accepted/steps = 0.0" (0.0, M.acceptanceRate { accepted = 0, steps = 10 })
      val () = checkApprox "accepted/steps = 1.0" (1.0, M.acceptanceRate { accepted = 10, steps = 10 })
      val () = Harness.checkRaises "acceptanceRate with steps <= 0 raises"
                 (fn () => M.acceptanceRate { accepted = 0, steps = 0 })

      val () = Harness.section "diagnostics: empirical acceptance rate in (0,1)"
      val seed = SplitMix64.seed 0wxFACEFEED
      val ({ accepted, ... }, _) =
        M.metropolisWithAccepts
          { logTarget = logTarget, proposalStdev = 2.4, init = [0.0], steps = 5000 }
          seed
      val rate = M.acceptanceRate { accepted = accepted, steps = 5000 }
      val () = Harness.check "acceptance rate > 0" (rate > 0.0)
      val () = Harness.check "acceptance rate < 1" (rate < 1.0)
      val () = Harness.check "accepted count within [0, steps]"
                 (accepted >= 0 andalso accepted <= 5000)

      val () = Harness.section "diagnostics: accepted count is reproducible"
      val ({ accepted = acc2, ... }, _) =
        M.metropolisWithAccepts
          { logTarget = logTarget, proposalStdev = 2.4, init = [0.0], steps = 5000 }
          (SplitMix64.seed 0wxFACEFEED)
      val () = Harness.checkInt "same seed -> same accepted count" (accepted, acc2)
    in
      ()
    end
end
