# Real ORB5 run: functional test report

Environment: MATLAB R2023b. Source files: ORB5_Plot_Suite_GUI/input and orb5_res.h5 (4,561,078,472 bytes). Simulation files were read only.

The input declares deuterium and electrons, with electrons marked is_kinetic = .false. Transport output exists for deuterium only. The potential is electrostatic, with selected toroidal modes n = 0 and n = 40.

## Fixes verified

- Single-precision HDF5 coordinates caused numeric GUI field construction to fail. Control values now use double precision.
- Single-precision normalized timestamps falsely failed the uniform-sampling check. Normalization now starts from double timestamps.
- Selected toroidal modes are read directly: 2 of 41 stored slots for this run, without changing the retained spectral coefficients.
- Time-slice CWT avoids the full frequency-time-radius allocation. Independent complex-signal tests compare it with the full toolbox transform, including boundary samples.
- Cross-correlation redraw no longer attempts to copy unsupported yyaxis axes; two linked panels preserve the signals.

## Individual real-file checks

Times include plotting, repeated selection for interactive plots, and one PNG export. They are one-run observations rather than formal benchmarks.

| Function / plot | Result | Seconds |
| --- | --- | ---: |
| load_potential | PASS | 1.22 |
| phimax_lfs_Cs | PASS | 1.26 |
| radial_st | PASS | 1.48 |
| radial_potential | PASS | 1.93 |
| temporal_potential | PASS | 1.90 |
| plot_potsc | PASS | 2.13 |
| phi_s_omega | PASS | 2.07 |
| f_Cs | PASS | 0.96 |
| cross_correlation | PASS | 1.95 |
| plot_dTi_theory | PASS | 2.08 |
| hfs | PASS | 0.01 |
| get_theory_overlay | PASS | 0.36 |
| load_transport | PASS | 0.97 |
| radial_flux | PASS | 2.44 |
| temporal_flux | PASS | 2.41 |
| radial_transport | PASS | 1.95 |
| temporal_transport | PASS | 2.01 |
| transport_st | PASS | 0.48 |
| transport_st_RoLT_component | PASS | 0.25 |
| eflux0d | PASS | 0.27 |
| energy | PASS | 0.39 |
| convergence | PASS | 0.41 |
| load_equilibrium | PASS | 0.38 |
| equilibrium | PASS | 0.66 |
| load_calculator | PASS | 0.01 |
| calculator | PASS | 1.52 |

## Coverage and limits

- All 19 registered GUI plots were invoked, plus four loaders, plot_dTi_theory, hfs and the optional radial theory overlay.
- Interactive plots were rendered at an initial and a different selected position/time. PNGs beside this report show the real data.
- Synthetic regression tests cover missing species, an additional ion, separate diagnostic clocks, shape errors, single-precision numeric entry, and timer/window cleanup.
- The analytic theory overlay was checked at time sample 31; its physical assumptions were not redefined.
- The supplied file has no A-parallel output. HFS reconstruction was compared numerically using its scalar-potential coefficients; synthetic EM data cover the EM preparation branch.
- Functional success does not constitute independent validation of all physics formulas. Previously documented experimental assumptions remain.
- Hysteresis remains unavailable at the user's request.

## Final GUI and supplemental verification

- test_real_gui: all 19 main-window actions passed using the supplied run.
- run_orb5_tests: passed, including single-precision coordinate entry, switching to index entry, synthetic EM preparation and GAM/VMD decomposition.
- test_cwt_slices: passed for odd/even signal lengths, multiple complex wavelets, and first/middle/last selected times.
- The final GUI run also passed frequency analysis with explicit real-component display and no implicit complex-value plotting warnings.

## Grouped-module refactor verification — 2026-09-24

The 31 private implementation files have been replaced with five grouped modules:
potential, transport, equilibrium, shared data utilities and interface utilities.
The calculator remains in ORB5_calculator.m. Public filenames forward to the grouped
implementations, preserving existing calls. Original-suite and simulation files
were not edited.

The table above was refreshed after the final cleanup. All 19 registered plots,
four loaders, the dTi comparison, HFS reconstruction and optional theory overlay
passed again. All 19 main GUI actions passed. Synthetic checks passed for species
availability, independent clocks, known moments, EM/GAM, coordinate/index controls,
timer cleanup and malformed inputs.

Numerical verification:

- Prepared potential and transport fields and sampled spectral coefficients match
  the saved pre-refactor baseline (double relative tolerance 1e-12; single 5e-6).
- Sparse CWT matches the full toolbox transform for odd/even lengths, two complex
  wavelets and first/middle/last times. Reusing the kernel with changed signal
  values and sampling frequency also passes.
- Vectorized correlation traces and lag curves match the previous per-time loop
  on nonuniform radial grids, including both boundaries and shortened filters.
- MATLAB Code Analyzer reports no findings in the potential or shared data module.
  Remaining advisory findings concern figure creation inside plot loops, tiny
  legacy alias strings, and the nested cleanup callback that releases the GUI
  busy flag. Repeated callback and cleanup behavior passes the GUI tests.

Representative timings, in seconds:

| Operation | Before | After |
| --- | ---: | ---: |
| Correlation plotting, median of 3 | 0.8241 | 0.5090 |
| Selected-time CWT, median of 3 | 0.1694 | 0.0753 |
| Potential preparation, one observation | 5.6294 | 1.7575 |
| Transport preparation, one observation | 2.0348 | 1.2195 |

Repeated correlation and spectrum operations were approximately 38% and 56%
faster in this check. Loading timings are strongly affected by operating-system
file caching and should not be interpreted as an isolated refactor speedup.
The CWT retains one wavelet response, at most 128 MiB, to accelerate repeated
slice selections. Physics formulas remain unchanged pending scientific review.

Reproduction details and raw evidence: ../check_refactor.m,
../refactor_benchmark.json, ../refactor_baseline.mat and ../refactor_verification.log.

