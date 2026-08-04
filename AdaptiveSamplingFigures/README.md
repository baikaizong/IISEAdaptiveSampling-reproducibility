# Adaptive Sampling Figures

Publication-ready figure generation code and final research figures for the Adaptive Sampling project.

The project deliberately keeps three layers together:

1. the plotting source code;
2. the minimum input data required to reproduce the figures; and
3. the validated final PNG and EPS artifacts used by the manuscript.

All reproducible figures are exported as paired 600 dpi PNG and EPS files. PNG files provide convenient preview and document integration, while EPS files preserve vector text, axes, curves, and annotations wherever the source data permit.

## Repository layout

| Path | Purpose |
| --- | --- |
| `src/` | Independent plotting scripts, the batch generator, EPS optimization, and output validation. |
| `ExternalData/` | Minimum numerical inputs required by the plotting scripts. |
| `assets/` | Source artwork that cannot be reconstructed from the numerical inputs alone. |
| `FinalVision/` | Validated final PNG/EPS pairs retained as project artifacts. |
| `FIGURE_MANIFEST.csv` | Mapping from each final figure to its generator, inputs, and vector/raster characteristics. |

## Reproducible environment

The complete figure set was validated with Python 3.12.3. Direct dependencies are pinned to the versions used for the final rendering.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

On Windows PowerShell, activate the environment with:

```powershell
.venv\Scripts\Activate.ps1
```

## Generate and validate figures

Run the following commands from this directory:

```bash
make figures
make verify-strict
```

`make figures` overwrites the corresponding tracked artifacts in `FinalVision/`. It never substitutes synthetic data when a required input file is missing. Scripts that intentionally generate illustrative stochastic fields use fixed random seeds.

The batch workflow uses Matplotlib's non-interactive `Agg` backend, so it can run in GitHub Actions and other headless environments. After rendering, embedded raster layers in large EPS files are compressed while vector text, curves, axes, and layout boundaries remain intact. Each EPS artifact is kept below 5 MB.

`make verify-strict` checks all 25 expected figure stems for:

- paired PNG and EPS output;
- valid EPS headers;
- a minimum 2,000-pixel long edge for PNG previews; and
- the 5 MB EPS size limit.

Raster heatmaps and solar frames remain raster layers when embedded in EPS. This is the correct publication representation for inherently pixel-based data; mathematical labels, curves, axes, and annotations remain vector elements.

## Final artifacts

Every file in `FinalVision/` is intentionally retained. This allows readers to inspect and reuse the publication artifacts without first recreating the Python environment.

When plotting code or input data change:

1. run `make figures`;
2. inspect the affected figures at their intended publication size;
3. run `make verify-strict`; and
4. retain the source change and the corresponding updated artifacts together.

## Data provenance and licensing

The public mapping records source identifiers and repository-relative paths without exposing developer-specific absolute paths. See [`ExternalData/README.md`](ExternalData/README.md) and [`ExternalData/PATH_MAPPING.csv`](ExternalData/PATH_MAPPING.csv).

When this directory is published as part of the main repository, the plotting code follows the license declared at the main repository root. The repository maintainer must separately confirm redistribution rights and citation requirements for research data and source artwork. If the root license does not cover data, add an explicit data license or data-use notice before public release.
