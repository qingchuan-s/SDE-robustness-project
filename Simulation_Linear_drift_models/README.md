# Simulation: Linear Drift Model

This folder contains simulation studies for the one-dimensional linear drift model

The model is used to compare parameter estimation under ODE and SDE formulations
in the presence of model misspecification and perturbations.

---

## Model

The linear drift model has drift
\[
f(x) = -a(x - b),
\]
corresponding to:
- an ODE model with measurement noise, and
- an SDE model (Ornstein–Uhlenbeck process) with system noise.

---

## Experiments included

The simulations investigate the effect of:
- finite sample size
- instantaneous perturbations
- randomly varying long-term mean
- misspecified drift terms

Both correctly specified and misspecified estimation scenarios are considered
(ODE vs SDE fitting).

---


## How to run

Run the main simulation scripts in this folder to generate simulated data
and parameter estimates. Figures are saved by the scripts.

---

## Output

- Parameter estimate distributions
- Bias and variance comparisons
- Figures illustrating robustness under perturbations
