### "Event-Specific Probabilities and Distributions" (ESPD) approach for modelling competing risks.

This repository contains the code and data associated with our paper *Estimating Event-Specific Probabilities and Time-to-Event Distributions From Censored Data For Modelling Competing Events in Discrete Event Simulations: A Tutorial.*

------------------------------------------------------------------------

#### 📖 Abstract

**Background**: Although several strategies for modelling competing events in discrete event simulations (DES) exist, a methodological gap for the event-specific probabilities and distributions (ESPD) approach when dealing with censored data remains. This study defines and illustrates the ESPD strategy.

**Methods**: The ESPD approach assumes that events are generated through a two-step process. First, the type of event is selected according to some (unknown) mixture proportions. Next, times of occurrence are sampled from a corresponding survival distribution. Both of these steps can be modelled based on covariates. Performance was evaluated through a simulation study, considering sample size and levels of censoring. Additionally, a case study was conducted to assess the approach's ability to produce realistic results, and to demonstrate its implementation using both frequentist and Bayesian frameworks in R.

**Results**: The simulation study showed good performance of the ESPD approach, with accuracy decreasing as sample size and censoring levels increased. The average relative absolute error of the event probability (95%-confidence interval) ranged from 0.04 (0.00; 0.10) to 0.23 (0.01; 0.66) for 60% censoring and sample size 50. The approach yielded realistic results in the case study.

**Conclusion**: The ESPD approach can be used to model competing events in DES based on censored data. Further research is warranted to compare the approach to other modelling approaches for DES, and to evaluate its usefulness in estimating cumulative event incidences in a broader context.

------------------------------------------------------------------------

#### 🔍 Repository content

``` bash
├── bayes_model  
│   └── weibull_mix_cov.stan      # full stan code of the Bayesian model      
├── ready_made_function
│   └── ESPD_frequentist.R        # ready-to-use ESPD function
├── simulation_study
│   ├── simulation_study.R        # script to test performance & accuracy of ESPD approach
│   └── .RDS and .csv files       # results of the simulation study
├── step_by_step_frequentist.Rmd  # competing risk implementation & case study results, frequentist way
├── case_study_bayesian.Rmd       # case study results, Bayesian way
├── CompetingEvents_ESPD.Rproj    # R project for this repository
└── README.md
```

#### 🎬 Simulation Study

Performance was evaluated through a simulation study, considering sample size and levels of censoring.

The considered hypothetical scenario included k=2 competing events: recurrence (recur) and death before recurrence (death).

#### 🔧 Case study

The repository contains **two ESPD implementations** in R, both incorporating a case study, using the `melanoma` dataset from the `boot` library.

1.  `step_by_step_frequentist.Rmd` is the notebook providing a thorough frequentist step-by-step implementation of the log-likelihood and directly apply it to the melanoma case study. Note that this notebook goes from 'zero to hero', implementing the log-likelihood from its simplest version to incorporating censoring, followed by adding covariates, and finally to considering two competing risks.

2.  `case_study_bayesian.Rmd` is the notebook of the Bayesian implementation, which uses the Stan model located in the *bayes_model* folder.

#### 💪 Our ready-made function

We provide a general frequentist function that can be used to apply the ESPD approach for modelling two competing events: `ESPD_frequentist.R`, in the *ready_made_function* folder.

It provides all the functionalities one may require, such as allowing for different parametric families of distributions for individual risks.

Please get in touch with any questions or suggestions!
