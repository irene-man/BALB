functions {

  // Compute the factorial of n recursively
  int factorial(int n);
  int factorial(int n) { 
    int output;

    if (n == 0)
      output = 1;
    else
      output = n * factorial(n - 1); 

    return output; 
  } 

 // Get all permutations of array [1,...,n]
 int[,] getPermutations(int n) {
    int nPermutations = factorial(n);
    int m1[1,1] = {{1}};
    int m2[2,2] = {{1, 2},
                   {2, 1}};
    int m3[6,3] = {{1, 2, 3},
                   {1, 3, 2},
                   {2, 1, 3},
                   {2, 3, 1},
                   {3, 1, 2},
                   {3, 2, 1}};
    int m4[24,4] = {{1, 2, 3, 4},
                    {1, 2, 4, 3},
                    {1, 3, 2, 4},
                    {1, 3, 4, 2},
                    {1, 4, 2, 3},
                    {1, 4, 3, 2},
                    {2, 1, 3, 4},
                    {2, 1, 4, 3},
                    {2, 3, 1, 4},
                    {2, 3, 4, 1},
                    {2, 4, 1, 3},
                    {2, 4, 3, 1},
                    {3, 1, 2, 4},
                    {3, 1, 4, 2},
                    {3, 2, 1, 4},
                    {3, 2, 4, 1},
                    {3, 4, 1, 2},
                    {3, 4, 2, 1},
                    {4, 1, 2, 3},
                    {4, 1, 3, 2},
                    {4, 2, 1, 3},
                    {4, 2, 3, 1},
                    {4, 3, 1, 2},
                    {4, 3, 2, 1}};
    int output[nPermutations, n];

    if (n == 1)
      output = m1;
    else if (n == 2)
      output = m2;
    else if (n == 3)
      output = m3;
    else if (n == 4)
      output = m4;
    else
      output = output;
    return output;
  }

  // Compute permutations of array a
  int[,] permuteArray(int[] a) {
    int nElements = num_elements(a);
    int nPermutations = factorial(nElements);
    int permutations[nPermutations, nElements] = getPermutations(nElements);
    int output[nPermutations, nElements];

    for (i in 1: nPermutations)
        output[i] = a[permutations[i]];

    return output;
  }

  // Get number of transitions and the types changed from a diff vector
  int[] getNTransitionsAndTypes (row_vector diff) {
    int nTransitions = 0;             // variable to track how many types have changed
    int output[num_elements(diff) + 1] = rep_array(0, num_elements(diff) + 1);

    // Loop through each type to see whether is has changed
    for (i in 1:num_elements(diff)) {
      if (diff[i] != 0) {
        // if changed
        nTransitions += 1;            // increase the number of types changed
        output[nTransitions + 1] = i; // save a type index if it has changed
      }
    }
    output[1] = nTransitions;

    // first element indicates number of transitions, the remainder lists the changed types
    return output;
  }

  // Get all compatible paths of min. length between xStart and xEnd
  matrix[] getAllPaths(row_vector xStart, row_vector xEnd) {  
    row_vector[num_elements(xStart)] diff = xEnd - xStart;
    int nTransitionsAndTypes[num_elements(xStart) + 1] = getNTransitionsAndTypes(diff);
    int nTransitions = nTransitionsAndTypes[1];
    int changedTypes[nTransitions] = nTransitionsAndTypes[2:(nTransitions + 1)];
    int nPaths = factorial(nTransitions);
    int pathLength = nTransitions + 1;
    int allPermutations[nPaths, nTransitions] = permuteArray(changedTypes); 
    matrix[nTransitions + 1, num_elements(xStart)] paths[nPaths];

    for (i in 1:size(allPermutations)) {
      paths[i, 1] = xStart;
      paths[i, nTransitions + 1] = xEnd;
      for (j in 2:nTransitions) {
        paths[i, j] = paths[i, j - 1];
        paths[i, j, allPermutations[i, j - 1]] += diff[allPermutations[i, j - 1]];
      }
    }

    return paths;
  }

  // Get stateID of a state, which is the row number of stateDic with row_vector that equals state
  int getStateID(row_vector state, matrix stateDic) {
    int stateID = 1;
    int numEqual;

    while (stateID <= rows(stateDic)) {
      numEqual = 0;

      // increase numEqual if the j-th element of state equals that of stateDic[stateID, :]
      for (j in 1:num_elements(state)) {
        if (state[j] == stateDic[stateID, j])  
          numEqual += 1;
      }

      // if all elements in state equal those in stateDic[stateID], break while-loop and return stateID
      // if not, increase stateID and re-enter while-loop
      if (numEqual == num_elements(state)) 
        break; 
      
      stateID += 1;
    }

    return stateID;
  }

  // Get the diagonal elements of matrix Q corresponding to pathStateIDs
  row_vector getQDiagElements(int[] pathStateIDs, matrix Q) {
    row_vector[num_elements(pathStateIDs)] diag;

    for (i in 1:num_elements(pathStateIDs))
      diag[i] = Q[pathStateIDs[i], pathStateIDs[i]];
    
    return diag;
  }

  // Get the off-diagonal elements of matrix Q corresponding to pathStateIDs
  row_vector getQOffDiagElements(int[] pathStateIDs, matrix Q) {
    row_vector[num_elements(pathStateIDs) - 1] offDiag;

    for (i in 1:(num_elements(pathStateIDs) - 1))
      offDiag[i] = Q[pathStateIDs[i], pathStateIDs[i + 1]];
    
    return offDiag;
  }

  // compute the probability/likelihood of a path
  real computePathProb (real deltaT, row_vector diag, row_vector offDiag) {
    // Require: n>=2
    int n = num_elements(diag);
    row_vector[n-1] a;
    real part1;
    real part2;
    real prob;              

    a = diag[:(n-1)] - diag[n];    
    prob = pow(-1, n-1) / prod(a);    

    for (i in 1:(n-1)){
      part1 = prod(a[i] - a[1:(i - 1)]);
      part2 = prod(a[i] - a[(i + 1):(n - 1)]);
      prob += exp(a[i] * deltaT) / (a[i] * part1 * part2);
    }
    
    prob *= exp(diag[n] * deltaT);
    prob *= prod(offDiag);

    return prob;
  }

  // Compute the likelihood of y given parameters Q
  real computeCustom_lpdf (matrix y, matrix Q, matrix stateDic) {
    int nTypes = cols(y) - 1;
    real deltaT = y[2, 1] - y[1, 1];
    row_vector[nTypes] xStart = y[1, 2:];
    row_vector[nTypes] xEnd = y[2, 2:];
    row_vector[nTypes] diff = xEnd - xStart;
    int nTransitions = getNTransitionsAndTypes(diff)[1];
    int nPaths = factorial(nTransitions);
    int pathLength = nTransitions + 1;
    matrix[pathLength, nTypes] paths[nPaths] = getAllPaths(xStart, xEnd);
    real prob = 0;
    real logProb = 0;
    int pathStateIDs[nTransitions + 1];
    row_vector[nTransitions + 1] diag;
    row_vector[nTransitions] offDiag;

    if (nTransitions == 0) {
      pathStateIDs[1] = getStateID(paths[1, 1], stateDic);
      diag = getQDiagElements(pathStateIDs, Q);
      prob = exp(diag[1] * deltaT);
    } else {
      for (i in 1:nPaths) {
        for (j in 1:pathLength) {
          pathStateIDs[j] = getStateID(paths[i,j], stateDic);
        }         
        diag = getQDiagElements(pathStateIDs, Q);
        offDiag = getQOffDiagElements(pathStateIDs, Q);
        prob += computePathProb(deltaT, diag, offDiag);
      }
    }

    logProb = log(prob);
    return logProb;
  }

  // Compute acquisition transition rate matrix Q with
  matrix computeAcqMatrix(vector lambda, real logk, matrix stateDic, int nTypes, int nStates) {
    matrix[nStates, nStates] Q = rep_matrix(0, nStates, nStates);
    row_vector[nTypes] xStart;
    row_vector[nTypes] xEnd;
    row_vector[nTypes] diff;
    row_vector[nTypes] absDiff;
    int typeIdx;

    // Loop over all pairs of states 
    for(i in 1:nStates){
      for(j in 1:nStates){
        // Get the pair of states and their difference
        xStart = stateDic[i];
        xEnd = stateDic[j];
        diff = xEnd - xStart;
        for(k in 1: nTypes) 
          absDiff[k] = fabs(diff[k]);

        // Continue if xStart and xEnd are not adjacent/connected
        if (sum(absDiff) != 1) {
          Q[i, j] = 0;
          continue;                 
        } 

        // Get the type at interest
        typeIdx = getNTransitionsAndTypes(diff)[2];

        // Fill the off-diagonal entry with the appropriate rate
        if (diff[typeIdx] == 1){
          // Acquisition
          Q[i, j] = lambda[typeIdx];    
          if (sum(xStart) > 0) // In presence of other types
            Q[i, j] *= exp(logk);
        } 
      }
    }

    // Fill the diagonal entry with the appropriate negative value of the total rate
    for(i in 1:nStates)
        Q[i, i] = -sum(Q[i]);
      
    return Q;
  }

  // Compute clearance transition rate matrix Q with
  matrix computeClrMatrix(vector mu, real logh, matrix stateDic, int nTypes, int nStates) {
    matrix[nStates, nStates] Q = rep_matrix(0, nStates, nStates);
    row_vector[nTypes] xStart;
    row_vector[nTypes] xEnd;
    row_vector[nTypes] diff;
    row_vector[nTypes] absDiff;
    int typeIdx;

    // Loop over all pairs of states 
    for(i in 1:nStates){
      for(j in 1:nStates){
        // Get the pair of states and their difference
        xStart = stateDic[i];
        xEnd = stateDic[j];
        diff = xEnd - xStart;
        for(k in 1: nTypes) 
          absDiff[k] = fabs(diff[k]);

        // Continue if xStart and xEnd are not adjacent/connected
        if (sum(absDiff) != 1) {
          Q[i, j] = 0;
          continue;                 
        } 

        // Get the type at interest
        typeIdx = getNTransitionsAndTypes(diff)[2];

        // Fill the off-diagonal entry with the appropriate rate
        if (diff[typeIdx] == -1){
          // Clearance                    
          Q[i, j] = mu[typeIdx]; 
          if (sum(xEnd) > 0)  // In absence of other types
            Q[i, j] *= exp(logh);
        }
      }
    }

    // Fill the diagonal entry with the appropriate negative value of the total rate
    for(i in 1:nStates)
        Q[i, i] = -sum(Q[i]);
      
    return Q;
  }

}

data {
  int<lower=0> nRows;									// number of rows across all sampled individuals
  int<lower=0> nTypes;									// number of pneumococcus serotypes
  int<lower=0> nStates;									// number of carriage states
  int<lower=0> nPersons;								// number of persons in the data set
  matrix<lower=0,upper=1>[nStates, nTypes] stateDic;	// dictionary indicating the binary representations
  matrix<lower=0>[nRows, nTypes+2] y;					// all sampled time points and sampled states
  vector<lower=0>[nTypes] initLambda;					// initial estimate of type-specific baseline acquisition rates 
  vector<lower=0>[nTypes] initMu;						// initial estimate of type-specific clearance baseline rates
  // method does not work on lambda with identical elements, same for mu
}

parameters {
  vector<lower=0>[nTypes] lambda;						// type-specific baseline acquisition rates
  vector<lower=0>[nTypes] mu; 							// type-specific clearance rates
  real<lower=-3, upper=3> logk;							// log-value of interaction parameter in acquisition
  real<lower=-3, upper=3> logh;							// log-value of interaction parameter in clearance
  real<lower=0> alpha;									// hyper-prior parameter for the frailty
  vector<lower=0>[nPersons] z; 							// frailty parameters (random effects)
}

model {
  // define ancillary variables
  matrix[nStates, nStates] Qacq;        				// acquisition part of the transition rate matrix
  matrix[nStates, nStates] Qclr;        				// clearance part of the transition rate matrix
  matrix[nStates, nStates] Q;							// transition rate matrix
  int personID = 1;										// running index personID

  // specify priors
  for (i in 1:nTypes) {   
    lambda[i] ~ gamma(0.00001, 0.00001/initLambda[i]);
    mu[i] ~ gamma(0.00001, 0.00001/initMu[i]);
  }
  logk ~ uniform(-5, 5);
  logh ~ uniform(-5, 5);
  alpha ~ lognormal(0, 2);
  z ~ gamma(alpha, alpha);
  
  // initialize transition rate matrix
  Qacq = computeAcqMatrix(lambda, logk, stateDic, nTypes, nStates); 
  Qclr = computeClrMatrix(mu, logh, stateDic, nTypes, nStates);
  Q = z[personID] * Qacq + Qclr;
  
  // get loglikelihood for each pairs of consecutive observation time points 
  for (i in 1:(nRows-1)) {
    if (y[i,1] != y[i+1,1]) {
      personID += 1;
      Q = z[personID] * Qacq + Qclr;
	} else {
	  // add the contribution of a single pair of consecutive observation time points
	  y[i:(i+1), 2:] ~ computeCustom_lpdf(Q, stateDic);  
	}
  }
}
