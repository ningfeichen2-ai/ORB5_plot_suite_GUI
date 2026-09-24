# Review of the original scripts

This review concerns the supplied ORB5_Plot_Suite. Original files were not edited.

## Addressed in the new folder

- ORB5_plotter uses a GUI rather than an edited mode index.
- Missing caches no longer break initial EM detection.
- Prepared data use config.path; changing runs does not silently reuse old MAT caches.
- Species are discovered from diagnostics, and equilibrium availability is distinct
  from kinetic transport availability.
- Transport plots include each available species/quantity, including fast particles
  and additional named ion groups.
- Each transport diagnostic retains its own time axis.
- cmd_interface no longer blocks MATLAB with input() or a while loop.
- Slider, coordinate/index entry and timer controls replace the command prompt;
  plot windows are reused on updates.
- Dispatcher branches are registered in a separate file; energy/J dot E now has a route.
- Related implementations are grouped into potential, transport and equilibrium
  modules; shared data and GUI utilities occupy two further private modules.
  Existing public filenames remain as forwarding entry points.
- Numeric namelist parsing no longer evaluates arbitrary MATLAB expressions and
  handles comments, D exponents, quoted strings and multiple assignments per line.
- Harmonic extraction preserves singleton dimensions and matches scalar/vector requests.
- Radius-frequency plots use the selected time normalization; frequency evolution
  no longer always multiplies by the Cs conversion in EM mode.
- Absent field-energy diagnostics do not prevent potential-amplitude plots.
- Equilibrium species use their own coordinates; the original unexplained fast-density
  x100 display multiplier and hardcoded objective q curve are omitted from the generic view.
- Configurable theory coefficients replace the fixed-size zero C profile.
- Spectral/window failures are reported in the GUI while keeping it open.
- Real-file testing exposed single-precision GUI values; controls now convert
  coordinates to doubles, and time normalization also uses double precision.
- Potential loading reads the selected n=0 and n=40 slots directly instead of
  loading all 41 stored slots in the supplied run.
- Time-slice CWT evaluation avoids the full 19 GB frequency-time-radius array,
  with complex numerical equivalence tested against the full toolbox transform.
- Cross-correlation uses two linked time-trace panels because MATLAB cannot
  copy yyaxis axes when updating plot windows. Per-mode traces are retained.
- Correlation filtering processes time columns together, with numerical checks
  against the original loop on nonuniform grids and boundary positions.
- Repeated slice CWT evaluations reuse a bounded impulse-response cache;
  changed wavelets, signal lengths, input values and sampling rates are tested.
- Dead calculations and repeated mode-map construction are removed; HDF5 index
  traversal and diagnostic-clock collection allocate their collections once.
- Module comments describe array dimensions, time/radius conventions, optional
  diagnostics, cache lifetimes and the formulas retained for later review.

## Scientific issues deliberately not redefined

- The temperature expression uses vp2 + u2f - uf^2, divided by density and tau.
  Whether the moments are density weighted must be checked against ORB5 output definitions.
- Electron mass ratio and diffusivity scaling differ from the ion conventions.
  Original defaults are exposed instead of silently standardized.
- Zonal potential is selected by slot 6, not by a validated m=0 lookup.
- The density contribution contains the run-specific factor 1.5668.
- dTi-theory documentation/labels disagree with the active absolute-dTi and A formulas.
- The analytic radial theory and cross-correlation routines contain experimental
  assumptions; some original multi-mode summaries use only the last processed mode.
- VMD mode order is assumed to identify components and has not been physically verified.
- The original hysteresis implementation is missing. It remains unavailable at your request.

The supplied electrostatic run now passes all 19 registered plots; detailed
results are in tests/real_data/TEST_REPORT.md. EM/GAM helpers are covered with
synthetic data. A real EM result and reference plots would support further
scientific comparison.
