(* stats.sml (vendored subset for sml-mcmc)

   Implementation of `STATS` as a functor over a `sml-prng` RANDOM generator,
   plus a default `Stats` instantiated with SplitMix64. This is the descriptive
   subset of the upstream `sml-stats` library (the part `sml-mcmc` uses for
   chain diagnostics); see stats.sig for what is and is not vendored here.

   Everything is Basis-library `real`/`int` arithmetic threaded through pure
   helpers, so results are deterministic and identical under MLton and
   Poly/ML. The only non-elementary piece is `erf` via the Abramowitz & Stegun
   7.1.26 rational approximation, used by the normal cdf (|error| < 1.5e-7). *)

functor StatsFn (R : RANDOM) :> STATS where type rng = R.state =
struct
  type rng = R.state

  exception Empty

  val pi = Math.pi

  (* ---- descriptive statistics ---- *)

  fun sum xs = List.foldl (op +) 0.0 xs

  fun count xs = real (List.length xs)

  fun mean xs =
    case xs of
      [] => raise Empty
    | _  => sum xs / count xs

  (* Sum of squared deviations from the mean. *)
  fun ss xs =
    let val m = mean xs
    in List.foldl (fn (x, acc) => acc + (x - m) * (x - m)) 0.0 xs end

  fun variancePop xs =
    case xs of [] => raise Empty | _ => ss xs / count xs

  fun variance xs =
    case xs of
      []  => raise Empty
    | [_] => raise Empty
    | _   => ss xs / (count xs - 1.0)

  fun stddevPop xs = Math.sqrt (variancePop xs)
  fun stddev xs = Math.sqrt (variance xs)

  (* Ascending merge sort on a real list (Basis-only, stable). *)
  fun sorted xs =
    let
      fun merge ([], ys) = ys
        | merge (xs, []) = xs
        | merge (x :: xs, y :: ys) =
            if x <= y then x :: merge (xs, y :: ys)
            else y :: merge (x :: xs, ys)
      fun split [] = ([], [])
        | split [a] = ([a], [])
        | split (a :: b :: rest) =
            let val (l, r) = split rest in (a :: l, b :: r) end
      fun msort [] = []
        | msort [a] = [a]
        | msort ys =
            let val (l, r) = split ys
            in merge (msort l, msort r) end
    in
      msort xs
    end

  fun minimum xs =
    case xs of [] => raise Empty | y :: ys => List.foldl Real.min y ys
  fun maximum xs =
    case xs of [] => raise Empty | y :: ys => List.foldl Real.max y ys

  (* Linear-interpolation quantile (R/NumPy "type 7"). *)
  fun quantile q xs =
    case xs of
      [] => raise Empty
    | _  =>
        let
          val arr = Array.fromList (sorted xs)
          val n = Array.length arr
          val qq = if q < 0.0 then 0.0 else if q > 1.0 then 1.0 else q
          val h = real (n - 1) * qq
          val lo = Real.floor h
          val frac = h - real lo
          val xl = Array.sub (arr, lo)
        in
          if lo + 1 < n
          then xl + frac * (Array.sub (arr, lo + 1) - xl)
          else xl
        end

  fun median xs = quantile 0.5 xs

  (* ---- special functions ---- *)

  (* erf via Abramowitz & Stegun 7.1.26 (|error| < 1.5e-7). *)
  fun erf x =
    let
      val sign = if x < 0.0 then ~1.0 else 1.0
      val ax = Real.abs x
      val t = 1.0 / (1.0 + 0.3275911 * ax)
      val poly =
        ((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t
          - 0.284496736) * t + 0.254829592) * t
      val y = 1.0 - poly * Math.exp (~(ax * ax))
    in
      sign * y
    end

  (* ---- distributions ---- *)

  structure Normal =
  struct
    type param = { mu : real, sigma : real }

    fun pdf { mu, sigma } x =
      let val z = (x - mu) / sigma
      in Math.exp (~0.5 * z * z) / (sigma * Math.sqrt (2.0 * pi)) end

    fun cdf { mu, sigma } x =
      0.5 * (1.0 + erf ((x - mu) / (sigma * Math.sqrt 2.0)))

    (* Box-Muller: two uniforms -> one standard normal. *)
    fun sample { mu, sigma } s =
      let
        val (u1, s1) = R.real01 s
        val (u2, s2) = R.real01 s1
        (* guard ln 0 at the open end of [0,1) *)
        val u1' = if u1 <= 0.0 then 1.0E~300 else u1
        val z = Math.sqrt (~2.0 * Math.ln u1') * Math.cos (2.0 * pi * u2)
      in
        (mu + sigma * z, s2)
      end
  end
end

(* Default instantiation: deterministic SplitMix64 sampling. *)
structure Stats = StatsFn (SplitMix64)
