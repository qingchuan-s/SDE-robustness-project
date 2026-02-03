# Application: Danish COVID-19 Data

This folder contains the application of the methods proposed in the paper to
real-world COVID-19 data from Denmark.

The analysis compares parameter estimates obtained from deterministic and
stochastic epidemic models under missing data scenarios.

---

## Data
Daily epidemiological data are obtained from Statens Serum Institut (SSI) and include.

Processed data and analysis scripts are provided to reproduce the results.

Processed dataset includes:
- date
- daily diagnosed cases (New_cases)
- approximated infectious (approx_infectious_cases)	
- cumulative diagnosed cases (Cumulative_cases)	
- daily death reported （New_deaths）	
- cumulative deaths (Cumulative_deaths)


---

## Models

The following models are fitted to the data:
- Deterministic SIR and SEIR models
- Stochastic SIR and SEIR models with system noise

---

## Experiments included

The analysis investigates:
- sensitivity of parameter estimates to missing data
- stability under different observation windows
- robustness of SDE-based estimation compared to ODE-based estimation
---

## How to run

Run the analysis scripts in this folder to reproduce the COVID-19 data analysis.
Results are saved to disk and can be compared with those reported in the paper.

---

## Output

- Estimated epidemic parameters under different data subsets
- Sensitivity analysis results
- Figures illustrating robustness to missing data
