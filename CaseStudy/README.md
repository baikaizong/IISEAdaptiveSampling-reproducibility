# Solar monitoring case study

This directory contains the real-data workflow used to produce
the five solar monitoring statistics plotted by `AdaptiveSamplingFigures`.

## Data locations

- Raw frames:
  `AdaptiveSamplingFigures/ExternalData/solar_frames/data.mat`
- Published five-series result set:
  `AdaptiveSamplingFigures/ExternalData/solar_statistics/`
- Generated normalized cache:
  `cache/case-study/norm_data.bin` (ignored by Git)
- Run output:
  `results/case-study/` (ignored by Git)

`case_study_data.py` prepares the normalized monitoring input as follows:

1. convert the `232 x 292 x 300` frames to floating point;
2. take consecutive-frame differences;
3. estimate each pixel's mean and sample standard deviation from the first 100
   difference frames;
4. standardize all 299 difference frames;
5. round to six decimal places;
6. write little-endian float64 values in the same Fortran array order.

The binary cache contains the six-decimal normalized values and occupies
162,043,648 bytes.

## Build

Configure the repository as described in the root README, then build:

```bat
scripts\build-windows.cmd
```

The resulting executable is:

```text
build\windows-debug\bin\case_study.exe
```

It links IMSL for the random-number routines.

## Prepare and run

The Python dependencies are already listed in
`AdaptiveSamplingFigures/requirements.txt`.

Prepare only:

```bat
python CaseStudy\case_study_data.py prepare
```

Prepare if needed and run all five methods:

```bat
python CaseStudy\case_study_data.py run
```

The wrapper defaults to seed `20260804`; pass another integer with `--seed`
when intentionally starting a different reproducible random stream. A fixed
seed makes repeated runs comparable. Different seeds can produce different
valid trajectories for methods that use randomized node allocation.

The executable can also be called directly:

```text
case_study.exe [NORMALIZED_CACHE] [OUTPUT_DIRECTORY] [SEED]
```

Create the output directory before a direct call.

## Result verification

Compare a new run with the exact plotting inputs:

```bat
python CaseStudy\case_study_data.py verify
```

The command reports length and maximum absolute difference for each method and
returns a nonzero exit code if any series differs. Randomized node allocation
can produce different valid trajectories for different seeds. The files in
`AdaptiveSamplingFigures/ExternalData/solar_statistics` remain the canonical
published result set.
