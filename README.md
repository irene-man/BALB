# BALB - A Bayesian approximate likelihood-based estimation method

### Brief description
STAN and R codes for inferring multiple-type pathogen interactions from longitudial data. The code was used to estimate the interaction parameters in the multi-state models proposed in "Approximate likelihood-based estimation method of multiple-type
pathogen interactions: an application to longitudinal pneumococcal carriage data", *Journal, year (doi: )*.

### Prerequisites
- rstan
- dplyr
- gtools

### Example
- Number of types: p = 4
- Baseline acquisition rates: lambda_i = exp(-3.5)
- Baseline clearance rates: mu_i = exp(-1.5)
- Interaction parameter in acquisition: log(k) = log(0.5)
- Interaction parameter in clearance: log(h) = 0
- Variance in heterogeneity: z_i ~ Gamma(5,5)

- Input file: "example_input.csv"
- STAN file: "stancode.stan"
- R file: "rcode.R"
- Utility file: "utils.R"
