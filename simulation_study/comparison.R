
# Potential modeling approaches:
# - Event-specific probabilities and distributions (ESPD)
# - Event-specific distributions (ESD)
# - Unimodal distribution and regression (UDR)

# Potential scenario variables:
# - Data generating mechanisms: ESPD, ESD, UDR
# - Number of competing events: 2, 3, and 4
# - Levels of overlap: low (10%), moderate (50%), and high (90%)
# - Sample sizes: 50, 100, 200, 500
# - Levels of censoring: no censoring (0%), low (10%), moderate (30%), and high (60%)

# Potential outcomes:
# - Relative absolute incidence difference: | estimated - observed | / p_obs
# - Relative entropy of time-to-event distributions: Kullback-Leibler divergence
# - Combine into single outcome over multiple competing events by weighing by incidence

# Pseudocode for simulation study:
# 
# FOR each different data generating mechanism
# 
#   FOR each number of competing events
# 
#     FOR each level of overlap
#
#       Simulate a large population according to the generating mechanism, number of events,
#         and level of overlap
#
#       FOR each sample size
#
#         FOR each level of censoring
#
#           Draw a sample according to the sample size from the large population
#
#           Censor the sample according to the level of censoring
#
#           Analyze the sample according to the different modeling approaches
#
#           Simulate a cohort of individuals according to the different approaches
#
#           Assess the internal performance of the approaches by comparing the outcomes in the
#             simulated cohort that of the uncensored sample
#            
#           Assess the external performance of the approaches by comparing the outcomes in the
#             simulated cohort that of the large population
#           
#         END FOR
#
#       END FOR
#
#     END FOR
#
#   END FOR
#
# END FOR

# Questions / discussion points:
# - What would be the aim of the simulation study? Demonstrate that the ESPD approach works
#   well across a range of scenarios? Or to assess its performance relative to the other
#   modeling approach? On the latter we could probably write a separate paper and that may
#   be confusing to the reader... On the other hand, leaving the comparison out may yield
#   reviewer comments that no guidance is provided on when the approach should be used.
# - Accounting for heterogeneity is not considered in the above, but was considered in the
#   initial setup for the simulation study... We will have to consider whether it is
#   sufficient to illustrate multivariable specifications to account for heterogeneity in 
#   the case study, or whether it should be part of the simulation study. And, if so, in
#   what way alternative scenarios should be defined, or will just a single scenario be
#   used, e.g. such as the intial simulation study example?
