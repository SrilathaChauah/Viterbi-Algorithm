numstates <- c("A1","A2","A3","A4","A5") # Defining the names of the states

##Setting the probabilities of switching states
A1 <- c(0.6,0.1,0.1,0.1,0.1)
A2 <- c(0.1,0.6,0.1,0.1,0.1)
A3 <- c(0.1,0.1,0.6,0.1,0.1)
A4 <- c(0.1,0.1,0.1,0.6,0.1)
A5 <- c(0.1,0.1,0.1,0.1,0.6)

## Creating an 5 x 5 matrix
transitionMatrix <- matrix(c(A1,A2,A3,A4,A5), 5, 5, byrow = TRUE)

## Create an observationMatrix
## observation

##Implementing Viterbi Algorithm
viterbi <- function (transitionMatrix, observationMatrix, observation)
{
  
  states <- row.names(transitionMatrix)
  numobservation <- length(observation)
  numstates <- length(states)
  initial.probs <- rep(0.5, numstates)
  names(initial.probs) <- states
  Z <- array(NA, c(numstates, numobservation))
  row.names(Z) <- row.names(transitionMatrix)
  colnames(Z) <- 1:numobservation
  
  for (s in states) 
  {
    Z[s, 1] = log(initial.probs[s] * observationMatrix[s,observation[1]])
  }
  
  for (t in 2:numobservation) 
  {
    for (s1 in states)
    {
      max1 = NULL
      for (s2 in states) 
      {
        temp = Z[s2, t - 1] + log(transitionMatrix[s2,s1])
        max1 = max(max1, temp)
      }
      Z[s1, t] = log(observationMatrix[s1, observation[t]]) + max1
    }
  }
  
  viterbi = rep(NA, numobservation)
  for (s in states) 
  {
    if (max(Z[, numobservation]) == Z[s, numobservation]) 
    {
      viterbi[numobservation] = s
      break
    }
  }
  for (t in (numobservation - 1):1) 
  {
    for (s in states) 
    {
      if (max(Z[, t] + log(transitionMatrix[, viterbi[t + 1]])) == Z[s, t] + log(transitionMatrix[s,viterbi[t + 1]])) 
      {
        viterbi[t] = s
        break
      }
    }
  }
  return(viterbi)
}