(* entry.sml -- runs every suite and exits with a status code. *)

fun runAllSuites () =
  ( Harness.reset ()
  ; MetropolisTests.run ()
  ; DeterminismTests.run ()
  ; GibbsTests.run ()
  ; DiagnosticsTests.run ()
  ; Harness.run () )

fun main () =
  OS.Process.exit
    (if runAllSuites () then OS.Process.success else OS.Process.failure)
