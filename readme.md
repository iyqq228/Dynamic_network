# Dynamic Network

Code for dynamic network models with correlated edge frailty, including simulation studies and an application to ant interaction networks. 

## Code overview

| File | Purpose |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| `function_FN.R` | Shared functions for simulating network event data with time-varying covariates, fitting edge frailty models, and computing parameter estimates, confidence intervals, and variance component tests. |  |  |  |  |
| `simulation/MC_simulation.R` | Runs parallel Monte Carlo simulations, summarizes estimation bias, empirical standard deviations, estimated standard errors, and confidence interval coverage, and evaluates test performance. |  |  |  |  |
| `real_data/Real_data.R` | Main empirical analysis script: reads the ant data, prepares event-level data and covariates, fits the model, performs parameter inference, and produces plots. |  |  |  |  |

## Data files

- `time1.xlsx`: Daily interaction time records for ant pairs.
- `cov.xlsx`: Individual ant covariates, including body size, age, and group membership across periods.


