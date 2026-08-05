# mMKSTD simulation and performance evaluation

This repository contains the Fortran implementation and simulation programs used to evaluate the modified Maximum Kernel Spatial-Temporal Distance method (mMKSTD) and compare it with CDS, NAS, POS, SASAM, TOPR, TSBSS, and TSSPR. The Monte Carlo workloads use MPI master-worker execution.

The repository is organized as one integrated codebase. Every evaluation program links the same core method library, shared runtime setup, data preparation routines, sample generator, and output helpers. Generated caches and large intermediate traces are not versioned; the programs recreate required cache data when it is absent. Compact reference results are included for inspection and regression comparison.

## Contents

- [Methods and executables](#methods-and-executables)
- [Repository layout](#repository-layout)
- [Validated environment](#validated-environment)
- [Dependencies](#dependencies)
- [Build instructions](#build-instructions)
- [Run instructions](#run-instructions)
- [Experiment entry points](#experiment-entry-points)
- [Solar monitoring case study](#solar-monitoring-case-study)
- [Reproducibility protocol](#reproducibility-protocol)
- [Parameters and method-specific behavior](#parameters-and-method-specific-behavior)
- [Cache and result data](#cache-and-result-data)
- [Validation](#validation)
- [Known constraints](#known-constraints)
- [License](#license)

## Methods and executables

| App name | Executable | Role |
|---|---|---|
| `mmkstd` | `mmkstd_evaluation.exe` | Primary mMKSTD evaluation |
| `cds` | `cds_evaluation.exe` | CDS comparison method |
| `nas` | `nas_evaluation.exe` | NAS comparison method |
| `pos` | `pos_evaluation.exe` | POS comparison method |
| `sasam` | `sasam_evaluation.exe` | SASAM comparison method |
| `topr` | `topr_evaluation.exe` | TOPR comparison method |
| `tsbss` | `tsbss_evaluation.exe` | TSBSS comparison method |
| `tsspr` | `tsspr_evaluation.exe` | TSSPR comparison method |
| `case-study` | `case_study.exe` | Real solar-data monitoring workflow |

The Fortran module and procedure identifiers `mMSTD_mod`, `mMSTD`, `mMSTDII`, and `Run_mMSTD_*` are the stable internal API for mMKSTD. The repository name, executable, result directories, and user documentation use the mMKSTD name.

## Repository layout

```text
.
|-- CMakeLists.txt
|-- CMakePresets.json
|-- README.md
|-- AdaptiveSamplingFigures/
|   `-- figure pipeline and canonical solar data/results
|-- CaseStudy/
|   |-- CaseStudy.f90
|   |-- case_study_data.py
|   `-- README.md
|-- cache/
|   `-- README.md
|-- reference_results/
|   |-- mMKSTD/
|   |-- CDS/
|   |-- SASAM/
|   |-- TOPR/
|   |-- TSBSS/
|   `-- TSSPR/
|-- scripts/
|   |-- configure-windows.cmd
|   |-- build-windows.cmd
|   |-- run-windows.cmd
|   `-- smoke-test-windows.cmd
`-- src/
    |-- apps/<method>/
    |   |-- main.f90
    |   |-- Experiment_Runner_mod.f90
    |   `-- performance_evaluation.f90
    |-- core/
    |   `-- shared method, sampling, matrix, and data modules
    |-- gurobi/
    |   `-- gurobi_interface.cpp
    `-- support/
        |-- ExperimentSupport_mod.f90
        `-- RuntimeConfig_mod.f90
```

`src/core` contains one canonical copy of each module. Method-specific experiment definitions remain under `src/apps` because their parameters, chart-state handling, and output statistics differ.

## Validated environment

The complete repository was built and its eight executables were started successfully in the following environment on 2026-08-04.

| Component | Validated version |
|---|---|
| Operating system | Microsoft Windows 11 x86-64, build 10.0.22631.6199 |
| Visual Studio | Community 2022, 17.10.4 |
| Fortran compiler | Intel oneAPI `ifx` 2024.1.0 |
| C++ compiler | Intel oneAPI `icx` 2024.1.0 |
| MPI | Intel MPI 2021.12.0 |
| IMSL | Fortran Numerical Library 2022.1.0, x64 DLL interface |
| Gurobi | Optimizer 13.0, Windows x64 |
| CMake | 3.28.3-msvc11 |
| Ninja | 1.11.0 |
| CMake configuration | Debug |

The CMake project requires CMake 3.25 or newer. Other tool versions may work, but they have not been validated for this repository. The current build and helper scripts target Windows because IMSL module/library naming and the validated compiler flags are Windows-specific.

## Dependencies

### Intel oneAPI and Intel MPI

Install the Intel oneAPI compiler and Intel MPI components. Every executable imports the Fortran `mpi` module. Start an x64 command environment in which Visual Studio, oneAPI, and Intel MPI have already been initialized. The following commands must all resolve through `PATH`:

```bat
where ifx
where icx
where mpiexec
```

An Intel oneAPI command prompt is the simplest option. Alternatively, call the Visual Studio x64 environment script followed by the oneAPI environment script from their local installation locations. These installation locations are intentionally not stored in this repository.

### Visual Studio, CMake, and Ninja

Install the Visual Studio C++ build tools required by `icx`, CMake 3.25 or newer, and Ninja. Confirm that the initialized shell provides them:

```bat
where cmake
where ninja
```

The scripts use `cmake` from `PATH` by default. Set `CMAKE_EXE` only when a different CMake executable must be selected for the current shell.

### IMSL Fortran Numerical Library

IMSL supplies random-number and numerical interfaces used throughout the code, including `RNSET_INT`, `RNUND_INT`, `RNNOR_INT`, `RNMVN_INT`, `CHFAC_INT`, and `SVRGP_INT`. Set `IMSL_ROOT` to the architecture-specific x64 DLL-interface directory whose relative layout includes:

```text
%IMSL_ROOT%
|-- include\dll\RNSET_INT.mod
`-- lib\imsl_dll.lib
```

The matching IMSL runtime libraries must also be present below `%IMSL_ROOT%\lib`. IMSL is commercial software and is not included in this repository.

### Gurobi Optimizer

`mMSTDII` solves a set-cover problem through `src/gurobi/gurobi_interface.cpp`. Set `GUROBI_ROOT` to the Gurobi installation directory. `GUROBI_HOME`, which the Gurobi installer commonly defines, is accepted as a fallback. The required relative layout is:

```text
%GUROBI_ROOT%
|-- include\gurobi_c.h
|-- lib\gurobi130.lib
`-- bin\gurobi130.dll
```

A valid Gurobi license must be available at runtime. Gurobi is commercial software and is not included in this repository.

The run and smoke-test wrappers prepend these two runtime directories to `PATH` and explicitly propagate that value to Intel MPI worker processes. This prevents MPI workers from depending on workstation-wide DLL search-path settings.

### Environment variables

Machine-specific installation locations belong in the current shell or the user's local environment, not in tracked files.

| Variable | Required | Meaning |
|---|---|---|
| `IMSL_ROOT` | Yes | Architecture-specific IMSL directory shown above |
| `GUROBI_ROOT` | Yes, unless `GUROBI_HOME` is set | Gurobi installation directory shown above |
| `GUROBI_HOME` | Optional fallback | Standard Gurobi installation variable |
| `CMAKE_EXE` | Optional | CMake command or local executable; defaults to `cmake` from `PATH` |
| `MMKSTD_SEED` | Optional at runtime | Seed used only when no command-line seed is supplied |

All paths committed to this repository are repository-relative or environment-variable based. CMake necessarily records resolved local dependency paths in the generated `build` tree; that tree is ignored by Git and must not be uploaded.

## Build instructions

Open an initialized x64 Visual Studio/oneAPI command environment. The scripts locate the repository root relative to their own directory, so they can be invoked from any working directory.

### Configure the local dependency environment

```bat
set "IMSL_ROOT=<architecture-specific IMSL directory>"
set "GUROBI_ROOT=<Gurobi installation directory>"
```

The values above are local installation settings. Do not replace the placeholders in tracked repository files and do not commit a generated CMake cache.

### Build

```bat
scripts\configure-windows.cmd
scripts\build-windows.cmd
```

The configure script verifies the toolchain and expected dependency files before running CMake. It performs a fresh configuration so a moved or newly cloned repository cannot reuse paths from an older CMake cache.

As an alternative to environment variables, pass the two dependency roots for the current invocation:

```bat
scripts\configure-windows.cmd "%IMSL_ROOT%" "%GUROBI_ROOT%"
scripts\build-windows.cmd
```

### Manual CMake invocation

After initializing the Visual Studio and oneAPI environments:

```bat
cmake --preset windows-debug ^
  -DIMSL_ROOT="%IMSL_ROOT%" ^
  -DGUROBI_ROOT="%GUROBI_ROOT%"
cmake --build --preset windows-debug --parallel
```

The Debug preset enables interface warnings, traceback information, bounds checking, stack checking, and heap allocation for automatic arrays. Executables are written to:

```text
build\windows-debug\bin
```

## Run instructions

The wrapper syntax is:

```text
scripts\run-windows.cmd APP [EXPERIMENT] [MPI_PROCESSES] [RESULTS_ROOT] [CACHE_ROOT] [SEED]
```

Example:

```bat
scripts\run-windows.cmd mmkstd circle-comparison 4 results cache 20260804
```

Equivalent direct invocation:

```bat
mpiexec -n 4 build\windows-debug\bin\mmkstd_evaluation.exe ^
  circle-comparison results cache 20260804
```

The executable argument order is:

```text
EXPERIMENT [RESULTS_ROOT] [CACHE_ROOT] [SEED]
```

At least two MPI processes are required by the master-worker evaluation routines. Wrapper defaults are four processes, `results`, `cache`, and seed `20260804`. The `MMKSTD_SEED` environment variable supplies a seed when the direct executable call omits the final argument; an explicit command-line seed takes precedence.

## Experiment entry points

| App | Accepted experiment names | Default |
|---|---|---|
| `mmkstd` | `bspline`, `kernel-circle`, `circle-comparison`, `circle`, `spatial-comparison`, `st-comparison`, `ellipse`, `crescent`, `qnum` | `circle-comparison` |
| `cds` | `kernel-circle` | `kernel-circle` |
| `nas` | `bspline`, `kernel-circle` | `kernel-circle` |
| `pos` | `bspline`, `kernel-circle` | `kernel-circle` |
| `sasam` | `bspline`, `kernel-circle`, `ellipse`, `crescent`, `exploration`, `circle`, `qnum` | `kernel-circle` |
| `topr` | `bspline`, `kernel-circle` | `kernel-circle` |
| `tsbss` | `bspline` | `bspline` |
| `tsspr` | `bspline`, `kernel-circle` | `kernel-circle` |

Each experiment creates a method/scenario subdirectory below `RESULTS_ROOT`. Configuration files such as `Parameter_Settings.txt` and `Config_*.txt` record the effective method and simulation settings beside the numerical output.

## Solar monitoring case study

`case_study.exe` computes the five monitoring-statistic series used by the
solar-data figures.

Prepare the normalized cache and run with the fixed default seed:

```bat
python CaseStudy\case_study_data.py run
```

Generated cache and result files are written below `cache/case-study` and
`results/case-study`, respectively, and are ignored by Git. The raw MAT file
and five published series are under `AdaptiveSamplingFigures/ExternalData`.
Compare a run with those plotting
inputs using:

```bat
python CaseStudy\case_study_data.py verify
```

See [`CaseStudy/README.md`](CaseStudy/README.md) for data preparation, running,
and reproducibility details.

## Reproducibility protocol

For a controlled replication:

1. Use the validated software versions or record every version substitution.
2. Start from an empty cache, or use a cache generated with exactly the same grid, covariance, B-spline, neighbor, and bandwidth settings.
3. Use an explicit nonnegative seed.
4. Keep the MPI process count fixed between compared runs.
5. Keep the experiment name and all constants in its `Experiment_Runner_mod.f90` unchanged.
6. Archive the generated configuration files together with the result files.
7. Compare run-length summaries and quantiles. Do not require byte-identical output across different MPI process counts, compiler versions, or scheduling decisions.

The explicit seed replaces clock-based initialization while preserving the rank- and task-specific seed formulas. MPI workers request tasks dynamically, so scheduling and floating-point reduction order can still affect exact output. The same seed, process count, executable, dependency versions, and cache configuration provide the strongest supported repeatability.

## Parameters and method-specific behavior

The shared defaults are declared in `src/core/GlobalSettings_mod.f90`:

| Parameter | Default |
|---|---:|
| Grid size | `100 x 100` nodes |
| Spatial neighbors | `5000` |
| Sampled nodes | `10` |
| Monte Carlo replications (`simu`) | `5000` |
| Maximum run length (`Maxarl`) | `2000` |
| Warm-up window (`Time_window`) | `50` |
| Limit-range samples (`IC_TestRuns`) | `10000` |
| In-control evaluation runs (`IC_runs`) | `1000000` |
| Nominal in-control ARL (`IcArl`) | `370` |
| Ground-process standard deviation | `3.0` |
| Observation-noise standard deviation | `1.0` |
| Kernel bandwidth | `1.0` |

Experiment runners override these defaults where required. The generated configuration file is authoritative for a specific run.

Important method behavior:

- mMKSTD limit searches and out-of-control evaluations call `mMSTDII`, which uses the Gurobi-backed set-cover step. Its serial in-control diagnostic and static limit-range estimator call `mMSTD`.
- The mMKSTD `bspline` entry uses the calibrated in-control limit `4.17287964774770`, matching the retained `BsplineNoBack-2/Limit_Result-1.txt` setting.
- NAS runners explicitly set `gp=1.0e-4`, `gone=1.0e-15`, `delta=1.0e-4`, `lambda0=1.0e-5`, `lambda=0.002`, and `k_allowance=0.01`.
- POS uses `POS_mod` and constructs both spatial kernel-neighbor tables required by the method. Before limit estimation, it calibrates the charting mean and standard deviation using a warm-up of `max(Time_window, omega_POS)` steps followed by 2000 in-control observations per MPI rank.
- CDS uses the dense-covariance conditional inference implementation required by its evaluation interface.

## Cache and result data

The repository does not store generated `.dat` caches. A first run creates the required `Prepared_settings*` directories under `CACHE_ROOT` and generates node coordinates, neighbor tables, B-spline bases, covariance structures, and method-specific preparation data. These files can occupy many gigabytes and take substantial time to generate.

Keep a generated cache to accelerate repeated runs with identical settings. Delete the affected cache directory after changing any of the following:

- grid dimensions or sampled-node count;
- spatial or covariance neighbor counts;
- kernel or background bandwidth;
- B-spline knot counts or degrees;
- covariance, ground-process, or noise parameters.

Cache files use formatted Fortran I/O. Use the validated compiler family when sharing caches; regenerating locally is safer when compiler/runtime details differ.

`reference_results` contains 2003 compact text result files totaling approximately 6.61 MiB. Large per-process distance traces and other intermediate outputs are not included because they are regenerated by production experiments and are not required to inspect summary results. New runs default to the ignored `results` directory, so retained reference results are not overwritten.

## Validation

Validation uses the Debug preset and the validated environment listed above.

The following targets compile and link:

```text
mmkstd_evaluation.exe
cds_evaluation.exe
nas_evaluation.exe
pos_evaluation.exe
sasam_evaluation.exe
topr_evaluation.exe
tsbss_evaluation.exe
tsspr_evaluation.exe
case_study.exe
```

Run the startup suite after building:

```bat
scripts\smoke-test-windows.cmd
```

The suite launches every executable with two MPI processes, a fixed seed, temporary result/cache roots, and an intentionally invalid experiment name. It verifies executable loading, MPI initialization/finalization, commercial runtime resolution, directory creation, and command dispatch without starting a production simulation.

Full default Monte Carlo experiments are intentionally not part of the smoke suite because they regenerate multi-gigabyte caches and can require long runtimes. The retained results provide numerical reference material, but a publication-quality replication should rerun the exact scenarios and archive its generated configurations.

## Known constraints

- The convolution anomaly generator uses the shared `kernel_bandwidth`. The NAS and POS OCARL interfaces accept `anomaly_bandwidth`, but that argument is currently not applied independently by the sample generator. Treat it as non-operative unless the model definition is explicitly extended and revalidated.
- Full 100 by 100 scenarios allocate large matrices and automatic arrays. Sufficient RAM, storage, and runtime are required.
- Gurobi must find a valid license even though the repository only redistributes the wrapper source.
- Only the Windows/Intel toolchain in this README has been validated.
- Exact bitwise agreement is not expected after changing the MPI process count, compiler, numerical libraries, or task scheduling.

## Repository hygiene

`.gitignore` excludes build products, generated caches, and new result directories. `.gitattributes` and `.editorconfig` define text encoding and line endings. The GitHub Actions workflow rejects compiled artifacts, files larger than 95 MiB, and workstation-specific source paths. The workflow does not compile the project because IMSL and Gurobi require separately licensed installations.

All documentation and source comments are written in English. Source files use UTF-8 text encoding.

## License

No source-code license is included. In the absence of a license, copyright remains with the authors and no permission to copy, modify, or redistribute the code is granted beyond applicable law. Add an author-approved license before inviting reuse or redistribution.
