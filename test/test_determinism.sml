(* test_determinism.sml -- the same seed reproduces the EXACT same chain.

   Determinism is the whole point of threading a seeded `sml-prng` state
   instead of using ambient randomness. Two runs of the same sampler from the
   same seed must produce bit-for-bit identical chains; two runs from different
   seeds must (with overwhelming probability) differ. We also confirm the
   returned successor state is usable and that re-running from that state is
   itself reproducible. *)

structure DeterminismTests =
struct
  open Support
  structure M = Mcmc

  fun logTarget xs = case xs of [x] => ~0.5 * x * x | _ => raise Fail "1-D"

  fun runChain s =
    M.metropolis
      { logTarget = logTarget, proposalStdev = 1.5, init = [0.0], steps = 500 } s

  fun run () =
    let
      val sA = SplitMix64.seed 0wxABCDEF
      val sB = SplitMix64.seed 0wxABCDEF      (* same seed *)
      val sC = SplitMix64.seed 0wx123456      (* different seed *)

      val (chainA, nextA) = runChain sA
      val (chainB, _)     = runChain sB
      val (chainC, _)     = runChain sC

      val () = Harness.section "determinism: same seed -> identical chain"
      val () = Harness.checkInt "same length"
                 (List.length chainA, List.length chainB)
      val () = Harness.check "chains are bit-for-bit identical"
                 (chainEqual (chainA, chainB))

      val () = Harness.section "determinism: different seed -> different chain"
      val () = Harness.check "different seed yields a different chain"
                 (not (chainEqual (chainA, chainC)))

      (* Continuing from the returned successor state is itself reproducible. *)
      val () = Harness.section "determinism: threaded successor state is reproducible"
      val (cont1, _) =
        M.metropolis
          { logTarget = logTarget, proposalStdev = 1.5, init = [0.0], steps = 100 }
          nextA
      val (cont2, _) =
        M.metropolis
          { logTarget = logTarget, proposalStdev = 1.5, init = [0.0], steps = 100 }
          nextA
      val () = Harness.check "two runs from the successor state agree exactly"
                 (chainEqual (cont1, cont2))
      val () = Harness.check "continuation differs from the original chain"
                 (not (chainEqual (cont1, chainA)))
    in
      ()
    end
end
