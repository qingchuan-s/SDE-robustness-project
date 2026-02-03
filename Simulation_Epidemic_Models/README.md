# Simulation: Epidemic Models 

This folder contains simulation studies for epidemic models used in the paper

The focus is on comparing deterministic (ODE) and stochastic (SDE) formulations
of epidemic models under model misspecification.

---

## Models

The following models are implemented:
- SIR model (ODE and SDE) for data generating and model fitting 
- SEIR model (ODE and SDE) for data generating

Stochastic models include system noise, while deterministic models include
measurement noise.

---

## Experiments included

Simulation studies investigate:
- correct vs misspecified noise structure
- instantaneous perturbations (e.g. superspreading events)
- fitting simplified SIR models to SEIR-generated data

Parameter estimation is performed using:
- least squares estimation (ODE models)
- pseudo-likelihood methods based on Strang splitting (SDE models)

---

## How to run

Run the simulation scripts to generate epidemic trajectories and estimate
model parameters under different scenarios.

Scripts are self-contained and save figures and summaries automatically.

---

## Output

- Parameter estimate distributions for contact and removal rates
- Comparisons between ODE and SDE robustness
- Figures illustrating bias and variance under perturbations
