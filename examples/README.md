This directory contains python and R scripts for preparing (advanced) splice_sim runs and 
analyzing the resulting data. Please note that those scripts do not run out of the box but require
you to update file paths and possibly other parameters to match your local execution environment.
Nevertheless, they should give a good starting point for preparing/analysing your own splice_sim simulation 
runs. 

- 3end_experiment
  - Scripts for preparing a 3'end experiment as described in our paper
- decay_experiment
  - Scripts for preparing a decay timecourse (exponential decay) for half-live estimation as described in our paper
- pulse_experiment
  - Scripts for preparing a pulse timecourse with increasing conversion rates
- slamdunk
  - Scripts for comparing a [Slamdunk](https://github.com/t-neumann/slamdunk) run with splice_sim 3'end simulation data 
- data_analysis
  - An R-markdown script that demonstrates some common analysis use-cases