# External figure data

This directory is the single input boundary for data required by the final figure pipeline. `PATH_MAPPING.csv` records public source identifiers and repository-relative destinations; it intentionally excludes developer-specific absolute paths.

## Included datasets

- `solar_statistics/`: five charting-statistics files used by the solar monitoring charts.
- `exploration/`: four losslessly compressed experiment files and one random-baseline file.
- `shape/`: four statistics files used by the ellipse and crescent comparison figures. Files were renamed to provide a consistent interface to the plotting code.
- `q_boxplot/`: twenty `OC_*.txt` files covering sampling budgets from `Q=5` through `Q=100` in increments of five.
- `solar_frames/`: `data.mat`, containing 300 solar image frames with spatial dimensions `232 x 292`.

Legacy inputs referenced only by retired notebooks are outside the minimum `FinalVision` reproduction chain and are not distributed with this project.

Do not add `Limit_*.txt` files from the original `Qnum` directory. The final Q-budget figure reads only files matching `OC_K-*_P-*_Q-*.txt`.

## Data handling rules

1. Plotting scripts must resolve inputs relative to the repository root; never commit workstation-specific absolute paths.
2. Preserve raw values. Lossless gzip compression is allowed for large text matrices and is handled transparently by the plotting code.
3. Record every distributed input in `PATH_MAPPING.csv` and every final output in `FIGURE_MANIFEST.csv`.
4. Do not replace missing experimental inputs with simulated values. A missing required input must produce an explicit error or skip message.

These files are included to reproduce the manuscript figures. Before publishing the repository, confirm data ownership, required citations, and redistribution permissions. Inclusion in this directory does not change the ownership or license of the underlying data.
