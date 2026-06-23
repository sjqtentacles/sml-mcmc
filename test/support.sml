(* support.sml -- shared helpers for the sml-mcmc tests.

   MCMC output is floating point, and sampling estimates only converge to the
   target moments *approximately*, so every real comparison goes through an
   explicit tolerance rather than string or structural equality: `Real.toString`
   differs between MLton and Poly/ML, and a finite chain never matches the
   target mean/variance exactly.

   Two tolerance scales are used:
     - `eps` (1e-9): tight, for algebraic identities and the *determinism*
       checks, where two runs with the same seed must agree bit-for-bit (we
       still compare with a tiny epsilon rather than `=` on reals to stay
       portable, but the difference is genuinely 0).
     - caller-supplied generous tolerances (via `approxTol`): for sampling
       estimates of means/variances after burn-in. *)

structure Support =
struct
  val eps = 1E~9

  fun approx (a, b) = Real.abs (a - b) <= eps

  (* approx with a caller-supplied tolerance. *)
  fun approxTol tol (a, b) = Real.abs (a - b) <= tol

  fun checkApprox name (expected, actual) =
    Harness.check name (approx (expected, actual))

  fun checkApproxTol tol name (expected, actual) =
    Harness.check name (approxTol tol (expected, actual))

  (* Exact (bit-for-bit) equality of two chains, used by the determinism
     tests. Reals are compared with `Real.==` so that two identically
     computed chains are required to match exactly. *)
  fun chainEqual (xs, ys) =
    ListPair.allEq
      (fn (a, b) =>
         ListPair.allEq (fn (x, y) => Real.== (x, y)) (a, b))
      (xs, ys)
end
