library(gtools)
library(dplyr)

getHiddenStates <- function(si, ds) {
  stateDic <- data.frame(si)        # Initialize content of stateDic
  changedIndices <- which(!ds == 0) # Indices of the types changed
  n <- length(changedIndices)       # Number of changed types
  m <- length(ds)                   # NUmber of types
  indices <- permutations(n, n, changedIndices)  # All permutation of changed types
  
  # Loop over all permutations to get intermediate states
  for (i in 1:nrow(indices)) {
    s <- si
    for (j in 1:(ncol(indices) - 1)) {
      s[indices[i, j]] <- s[indices[i, j]] + ds[indices[i, j]]
      stateDic <- rbind(stateDic, s)
    }
  }
  
  # Filter out unique states
  stateDic <- stateDic[2:nrow(stateDic),]
  stateDic <- unique(stateDic)
  
  return(stateDic)
}


getStateDicSimulation <- function(data) {
  # Initialize with observed states
  stateDic <- unique(data[, 3:ncol(data)])
  
  # Loop through all row
  for (i in 1:(nrow(data) - 1)) {
    # Check whether two consecutive rows are from the same person
    if (data[i, 1] == data[i + 1, 1]) {
      si <- data[i, 3:ncol(data)]		# Starting state
      sf <- data[i + 1, 3:ncol(data)]	# Ending state
      ds <- sf - si						# Difference between teh starting and ending states
      
      # Get all possible hidden states between starting and ending states
      if (length(which(!ds == 0)) > 1) {
        hiddenStates <- getHiddenStates(si, ds)
        stateDic <- rbind(stateDic, hiddenStates)
      }
    }
  }
  
  # Filter out unique states
  stateDic <- unique(stateDic)
  
  return(stateDic)
}


getInitBaselineRatesSimulation <- function(data) {
  nTypes <- ncol(data) - 2
  lambda <- rep(0, nTypes)
  mu <- rep(0, nTypes)
  for (i in 1:nTypes) {
    seq <- data[,2+i]
    df <- data %>% mutate(diff =  seq - lag(seq))
    deltaT <- (data %>% 
                 group_by(person_id) %>% 
                 mutate(deltaT = time - lag(time)) %>% 
                 summarise(meanDeltaT = mean(deltaT, na.rm = TRUE)) %>%
                 summarise(overallMeanDeltaT = mean(meanDeltaT, na.rm = TRUE)))[[1]]
    
    lambda[i] <- length(df$diff[df$diff==1]) / (length(seq[seq==0]) * deltaT) 
    mu[i] <- length(df$diff[df$diff==-1]) / (length(seq[seq==1]) * deltaT)
  }
  
  return(data.frame(lambda=lambda, mu=mu))
}