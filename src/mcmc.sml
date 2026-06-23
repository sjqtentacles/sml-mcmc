(* mcmc.sml

   Implementation of `MCMC` as a functor over a `sml-prng` RANDOM generator,
   plus a default `Mcmc` instantiated with SplitMix64.

   Everything is Basis-library `real`/`int` arithmetic threaded through pure
   helpers over an explicitly passed generator state, so results are
   deterministic and identical under MLton and Poly/ML. There is no ambient
   randomness: the *only* source of entropy is the threaded `R.state`.

   Gaussian proposals use the Box-Muller transform over two `sml-prng`
   uniforms, the same construction the vendored `sml-stats` uses for its normal
   sampler, so draws are byte-identical across compilers given the same seed.

   Diagnostics delegate to the vendored `sml-stats` (`Stats.mean`,
   `Stats.variance`); those are pure list statistics that do not involve the
   generator, so the default SplitMix64-instantiated `Stats` is fine regardless
   of which generator `Mcmc` is built over. *)

functor McmcFn (R : RANDOM) :> MCMC where type rng = R.state =
struct
  type rng = R.state

  exception Mcmc of string

  val pi = Math.pi

  (* ---- Box-Muller standard-normal draw over sml-prng uniforms ----
     Two uniforms in [0,1) -> one N(0,1) deviate. Identical construction to
     sml-stats's Normal.sample, so streams match across libraries/compilers. *)
  fun stdNormal s =
    let
      val (u1, s1) = R.real01 s
      val (u2, s2) = R.real01 s1
      (* guard ln 0 at the open end of [0,1) *)
      val u1' = if u1 <= 0.0 then 1.0E~300 else u1
      val z = Math.sqrt (~2.0 * Math.ln u1') * Math.cos (2.0 * pi * u2)
    in
      (z, s2)
    end

  (* ---- Metropolis-Hastings (symmetric Gaussian random walk) ---- *)

  (* Core routine: also tallies accepted proposals. The chain is accumulated in
     reverse and reversed once at the end, so this is O(steps * d). *)
  fun mhCore { logTarget, proposalStdev, init, steps } s0 =
    let
      val () = if steps < 0 then raise Mcmc "steps must be >= 0" else ()
      val () = if List.null init then raise Mcmc "init must be non-empty" else ()
      val () = if proposalStdev <= 0.0
               then raise Mcmc "proposalStdev must be > 0" else ()

      (* Propose a new point: current + N(0, proposalStdev^2) per coordinate. *)
      fun propose (point, s) =
        let
          fun loop ([], acc, s) = (List.rev acc, s)
            | loop (x :: xs, acc, s) =
                let val (z, s') = stdNormal s
                in loop (xs, (x + proposalStdev * z) :: acc, s') end
        in
          loop (point, [], s)
        end

      fun step (current, curLP, s) =
        let
          val (cand, s1) = propose (current, s)
          val candLP = logTarget cand
          val logAlpha = candLP - curLP
          (* accept with prob min(1, exp(logAlpha)); draw u ~ U[0,1) and accept
             iff ln u < logAlpha (always true when logAlpha >= 0). *)
          val (u, s2) = R.real01 s1
          val u' = if u <= 0.0 then 1.0E~300 else u
          val accept = Math.ln u' < logAlpha
        in
          if accept then (cand, candLP, true, s2)
          else (current, curLP, false, s2)
        end

      fun loop (0, current, _, acc, accepted, s) =
            (List.rev (current :: acc), accepted, s)
        | loop (k, current, curLP, acc, accepted, s) =
            let
              val (next, nextLP, didAccept, s') = step (current, curLP, s)
              val accepted' = if didAccept then accepted + 1 else accepted
            in
              loop (k - 1, next, nextLP, current :: acc, accepted', s')
            end

      val (chain, accepted, s') =
        loop (steps, init, logTarget init, [], 0, s0)
    in
      ({ chain = chain, accepted = accepted }, s')
    end

  fun metropolis args s =
    let val (r, s') = mhCore args s in (#chain r, s') end

  fun metropolisWithAccepts args s = mhCore args s

  (* ---- Gibbs sampling ----
     Each full sweep updates coordinate i (in index order 0..d-1) by drawing it
     from `conditionals[i]` applied to the *current full point*, which already
     reflects updates to coordinates < i within the sweep. *)

  fun gibbs { conditionals, init, steps } s0 =
    let
      val d = Vector.length conditionals
      val () = if steps < 0 then raise Mcmc "steps must be >= 0" else ()
      val () = if List.null init then raise Mcmc "init must be non-empty" else ()
      val () = if List.length init <> d
               then raise Mcmc "init length must equal number of conditionals"
               else ()

      val initArr = Array.fromList init

      (* Update coordinate i in `arr` in place using its conditional sampler.
         The sampler is handed the current full point as a list. *)
      fun updateCoord (arr, i, s) =
        let
          val point = Array.foldr (op ::) [] arr
          val sampler = Vector.sub (conditionals, i)
          val (xi, s') = sampler point s
        in
          Array.update (arr, i, xi); s'
        end

      (* One full sweep over all coordinates; returns the new state. *)
      fun sweep (arr, s) =
        let
          fun loop (i, s) =
            if i >= d then s
            else loop (i + 1, updateCoord (arr, i, s))
        in
          loop (0, s)
        end

      fun snapshot arr = Array.foldr (op ::) [] arr

      fun loop (0, acc, s) = (List.rev acc, s)
        | loop (k, acc, s) =
            let val s' = sweep (initArr, s)
            in loop (k - 1, snapshot initArr :: acc, s') end

      (* head = initial point, then one snapshot after each sweep *)
      val (rest, s') = loop (steps, [], s0)
      val chain = init :: rest
    in
      (chain, s')
    end

  (* ---- diagnostics ---- *)

  fun coordinate i chain = List.map (fn p => List.nth (p, i)) chain

  fun dropBurnIn n chain =
    if n <= 0 then chain
    else if n >= List.length chain then []
    else List.drop (chain, n)

  fun mean i chain = Stats.mean (coordinate i chain)
  fun variance i chain = Stats.variance (coordinate i chain)

  fun acceptanceRate { accepted, steps } =
    if steps <= 0 then raise Mcmc "steps must be > 0 for an acceptance rate"
    else real accepted / real steps
end

(* Default instantiation: deterministic SplitMix64 sampling. *)
structure Mcmc = McmcFn (SplitMix64)
