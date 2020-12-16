library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("utils.R")

# Read in data
data <- read.csv(file="example_input.csv")
head(data)
stateDic <- getStateDicSimulation(data)
baselineRates <- getInitBaselineRatesSimulation(data)
nPersons <- length(unique(data$person_id))
nTypes <- ncol(data) - 2
nParameters <- nTypes * 2 + 2 + 1 + nPersons
processed.data <- list(nRows = nrow(data),
                       nTypes = nTypes,
                       nStates = nrow(stateDic),
                       nPersons = nPersons,
                       stateDic = stateDic,
                       y = data,
                       initLambda = baselineRates$lambda,
                       initMu = baselineRates$mu)

# Fit stan model
mod <- stan_model('stancode.stan')
fit <- sampling(mod,
                data = processed.data,
                chains = 4,
                iter = 1000)

# Save estimates to dataframe
summary <- summary(fit, probs = c(0.025, 0.975))
est <- summary$summary[1:nParameters, "mean"]
ci_l <- summary$summary[1:nParameters, c("2.5%")]
ci_u <- summary$summary[1:nParameters, c("97.5%")]

