# mediation-power-simulation
Simulation to assess power and coverage of mediation parameters. 

There is a separate script for each data-generating mechanism (e.g., mediation_power_sim_dgm1.R). 

The baron_kenny.R, iorw.R, medtmle_intermedvar.R, medtmle_nointermedvar.R, mpower_eq.R, mpower_intermedvar.R, and mpower_nointermedvar.R files should be saved in the working directory.

Then this file will call those functions and perform a simulation parameterized based on the relevant data-generating mechanism. 

The simulation will calculate power and coverage across 1,000 simulations with a user-specified sample size for the Baron and Kenny, Inverse odds ratio weighting, Targeted minimum loss-based estimation, and the analytic equation. 

