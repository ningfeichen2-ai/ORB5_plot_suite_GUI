# ORB5 Plot Suite GUI

A separate MATLAB version of ORB5_Plot_Suite. The original folder is unchanged.
The original script/function filenames are retained. ORB5_plotter.m now opens a GUI.

## Start

In a fresh MATLAB session:

```matlab
cd('E:/Research/reseach/ORB5_script/ORB5_Plot_Suite_GUI')
ORB5_plotter
```

Choose the simulation folder containing `input` and `orb5_res.h5`.
Do not put both suite versions on the MATLAB path; the current directory and
already loaded classes can otherwise select the old version.

You can also start with explicit settings:

```matlab
config = orb5_config(struct('path','D:/runs/my_run'));
config.avg_window = [140 160];
config.time_windows = [75 85];
config.profile_species = 'deuterium'; % reference mapping s -> r/a
ORB5_plotter(config);
```

**Reload / inspect** discards loaded arrays and reports equilibrium and transport
species separately. Data load on demand and stay in memory for this GUI session.
Changing the run or settings invalidates that session. Old MAT caches are not read
or overwritten; this prevents accidentally plotting a stale cache from another run.
Use Reload after a running simulation changes its files. The potential loader
reads only the input-selected toroidal modes. Time-slice wavelet plots evaluate
the selected time without building a full frequency-time-radius cube.

The GUI has four blocks:

| Block | Content |
| --- | --- |
| Potential | Amplitudes, envelopes, radial harmonics, time traces, poloidal sections, spectra, cross-correlation |
| Transport | Particle/energy fluxes, J dot E, temperature, gradients, diffusivity, window averages |
| Equilibrium | Available density/temperature/gradient profiles and safety factor |
| ORB5 calculator | Existing CGS parameters and conversion coefficients, with configurable Z and mass number |

Plots open in normal MATLAB figure windows for zoom, pan, data tips and saving.
Slice plots get a separate controller with a slider, numeric coordinate entry,
index option, Previous/Next and Play/Pause. Numeric coordinates snap to the nearest
stored sample, and the controller shows the sample actually used. Out-of-range
typed values are rejected. Each species' flux curve uses its own nearest time and
reports it in the legend; diagnostics outside their time coverage are skipped.
Closing a controller deletes its timer and its associated plot windows.

## Species and diagnostic availability

Species names are discovered from HDF5 groups; plotting does not assume that
deuterium, electrons or fast particles exist. Equilibrium profiles alone do not
imply kinetic transport output. Every quantity is checked independently.

`read_flux_data` now returns
`data.species.<name>.diagnostics.<quantity>`, with `values`, `t`, and `s`.
Transport matrices have shape **radius x time**. Potential arrays retain
**time x toroidal mode x radius x poloidal slot**.

A missing optional diagnostic is skipped. Existing unreadable datasets, missing
required time axes and incompatible shapes raise an error instead of being
misrepresented as absent species. The reader accepts the original ORB5 paths and
a species-local pfluxf0_rad time fallback. Other output schemas need a reader
adapter, ideally checked against a representative real file.

Raw fluxes can be plotted for arbitrary discovered species. Derived temperatures
require all moments, matching clocks, and an explicit mass ratio.
Default mass ratios retain the original deuterium and electron assumptions.
For example, **only if appropriate for your normalization**:

```matlab
config.mass_ratios.fast = 1;
config.mass_ratios.helium = 2;
config.chi_scales.fast = (720/2)^2;
```

No mass is guessed for additional species. MATLAB-safe species field names use
`matlab.lang.makeValidName`.

The public `orb5_plot_library` methods retain their names. Common legacy fields
such as `f_avI`, `f_ave`, `Ti`, `Te`, `chi_perpI`, and `RoLTi`
are provided by the prepared transport loader. New extensions should use the
species/diagnostic structure. The direct reader return structure has changed.

## Physics assumptions

The existing implemented moment and normalization formulas are retained.
This is a software refactor, not a validation of their physical definitions.

| Setting | Default / purpose |
| --- | --- |
| `Z`, `mu` | 1, 2; reference calculator parameters |
| `is_em` | Empty: detect A-parallel output; true/false overrides |
| `time_units` | `normalized` or `raw` |
| `ad_hoc` | 1 for reference-profile r/a mapping, 0 for s |
| `profile_species` | deuterium; select another available reference profile as needed |
| `mass_ratios` | deuterium = 1, electrons = 1/1836, as in the original formulas |
| `tau` | 1 in the temperature expression |
| `density_factor` | 1.5668 in the density contribution |
| `chi_scales` | ion default (lx/2)^2; electron default 1, preserving the original difference |
| `zonal_slot` | 6; explicitly validate against the run's harmonic convention |
| `avg_window`, `time_windows` | Averaging and correlation windows in displayed time units |
| `theory_overlay` | false; experimental radial theory overlay is opt-in |
| `theory` | Original theory constants, q polynomial and smoothing window |
| `theory_C` | 0, preserving the effective original dTi-theory coefficient; scalar expands to the actual grid |

Supply these in the configuration passed to ORB5_plotter. Save config writes a MAT
file; restore with `saved=load('my_config.mat'); ORB5_plotter(saved.config)`.

Temperature/flux division results that are infinite are marked NaN. Nonpositive
density/temperature is excluded from logarithmic gradients. Short-grid smoothing
windows are reduced to valid sizes.

The original `plot_dTi_theory` implementation compares absolute dTi even though
some original labels say relative dTi, and its active A expression omits the
documented division by Ti0. Those numerical choices remain unchanged pending your
physics review. Its fixed 162-point zero C array is now a configurable scalar or
profile. The analytic radial theory overlay retains its experimental assumptions
and last-mode behavior. Frequency analysis and VMD separation also still need
physical validation on real runs.

Hysteresis intentionally remains unavailable until you provide its definition.

## Add your own plot

The GUI reads `orb5_plot_registry.m`. An entry declares the block, data source,
selector and a callback. The callback signature is **(data, config, index)**;
the index refers to displayed time for a time selector, or radius for a radius selector.
For `none`, ignore it. Draw ordinary MATLAB figures in the callback.

You can add an entry without modifying the GUI:

```matlab
config.extra_plots = struct( ...
    'id','my_density', 'label','My density plot', ...
    'block','Transport', 'source','transport', 'selector','radius', ...
    'callback',@plot_my_density);
ORB5_plotter(config);
```

Create `plot_my_density.m` in this folder:

```matlab
function plot_my_density(data, config, si)
figure; hold on
names = fieldnames(data.species);
for k = 1:numel(names)
    sp = data.species.(names{k});
    if ~isfield(sp.diagnostics,'f_av'), continue; end
    d = sp.diagnostics.f_av;
    plot(d.t,d.values(si,:),'DisplayName',sp.name);
end
xlabel(data.timelabel_str,'Interpreter','latex');
ylabel('f_av'); legend('show','Interpreter','none');
end
```

Available sources: `potential`, `transport`, `equilibrium`, `calculator`.
Available selectors: `none`, `time`, `radius`.
Do not use `close all`, timers, or modify unrelated figures in callbacks.

File responsibilities:

| File | Implementations kept together |
| --- | --- |
| `private/orb5_potential.m` | Potential reader, preparation, all potential plots, harmonic reconstruction, CWT, VMD and analytic overlay |
| `private/orb5_transport.m` | Species/flux reader, moments, all transport plots, convergence and dTi theory |
| `private/orb5_equilibrium.m` | Equilibrium reader and profile plots |
| `ORB5_calculator.m` | Plasma parameters and unit conversions |
| `private/orb5_interface.m` | Plot dispatch, interactive selectors, timers and figure updates |
| `private/orb5_data.m` | Shared HDF5 metadata, normalization, smoothing and radial differences |
| `ORB5_plotter.m` | Main GUI and its loaded-data session |
| `orb5_plot_library.m` | Existing public methods |
| `orb5_plot_registry.m` | Plot names, blocks and selectors |
| `orb5_config.m`, `get_orb5_val.m` | Configuration and input namelist parsing |

The former 31 private implementation files are consolidated into five modules.
Each module has an action dispatcher followed by related local functions. MATLAB's
editor function list and `%%` sections help navigate the larger files. The small
public files such as `read_flux_data.m` and `cwt_complex.m` remain as compatibility
entry points: MATLAB needs those filenames to preserve calls from your scripts.
Edit the grouped implementation, not the forwarding entry point.

To keep a new potential plot in the same file, add a local function to
`private/orb5_potential.m`, add its id to `render_potential`, then register the id
in `orb5_plot_registry.m` with source `potential` and callback `[]`. The registry
binds that id to the potential module. For transport, add the id to
`plot_transport_module` in `private/orb5_transport.m` and register it with source
`transport` and callback `[]`. External callbacks, as illustrated above, are an
optional alternative when a plot belongs to a separate extension.

Correlation smoothing works across all time columns in one operation. Slice CWT
reuses one impulse-response kernel keyed by sample count and wavelet; it retains
at most 128 MiB and contains no simulation values. Changes of wavelet or sample
count rebuild it. The standalone dTi comparison caches its prepared data separately;
call `plot_dTi_theory(config,r,'Reload',true)` after changing simulation files.

## Verification and remaining validation

Developed and exercised with installed MATLAB **R2023b**. Signal Processing
functions are used for Savitzky-Golay filtering and VMD; Wavelet Toolbox is needed
for CWT plots; original `smooth` calls may require Curve Fitting Toolbox.
Basic equilibrium/raw flux views do not perform wavelet transforms.

```matlab
addpath(fullfile(pwd,'tests'))
run_orb5_tests
```

Tests create temporary synthetic HDF5 files and check namelist parsing, species
absence, an additional arbitrary ion, species-specific clocks, known moment
values, compound spectral storage, radial/temporal/poloidal controls, repeated
window updates, timer cleanup, and malformed shapes. GUI exports are saved in
`tests/gui_preview.png` and `tests/slice_preview.png`.

The supplied 4.56 GB simulation file has now been exercised on MATLAB R2023b:
all 19 registered plots pass, including changes of slice in interactive plots.
The input declares deuterium and non-kinetic electrons; only deuterium transport
diagnostics exist. Electron transport curves are therefore intentionally absent.
See `tests/real_data/TEST_REPORT.md` for results, fixes and verification limits.

To reproduce the real-file checks (input and HDF5 remain read-only):

```matlab
run_real_data_tests
test_real_gui
test_cwt_slices
test_correlation_vectorization
```

`check_refactor('compare')` compares prepared real-run fields and CWT coefficients
against the saved `tests/refactor_baseline.mat` and writes timing observations to
`tests/refactor_benchmark.json`. Use `check_refactor('save')` only when intentionally
replacing that reference before future edits. Timings include local machine and
file-cache effects; they are not portable performance guarantees.

GUI coordinates accept single-precision simulation data by converting control
values to doubles. Time normalization is performed in double precision to avoid
false nonuniform-sampling errors. Different ORB5 output schemas and independent
scientific agreement with reference plots still require separate validation.
