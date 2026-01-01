# Simulation Study: Comparison of Ranking Methods for Functional Data

This repository contains the R code to reproduce the simulation study presented in the manuscript regarding the **FP-OWA** method.

## File Description
* `FP-OWA.R`: The main script that performs data generation, noise induction, and parallel comparison of 5 ranking methods.

## Methods Compared
1. **FP-OWA** (Adaptive Functional Piecewise Ordered Weighted Averaging)
2. **WLR** (Weighted Linear Rank)
3. **FPCA** (Functional Principal Component Analysis)
4. **H-Mode** (Hybrid Mode Depth)
5. **RTD** (Random Tukey Depth)

## How to Run
1. Open `FP-OWA.R` in R or RStudio.
2. Ensure all required packages listed in the script header are installed.
3. The script is configured to run a parallel simulation (default 50 runs for demonstration). You can adjust `n_runs` in the `run_simulation` call for full reproduction (e.g., 500 runs).
4. Results are saved as `simulation_results_summary.csv` in the working directory.

