# MATH-5090-Assignment
error_estimation <- function(n1, n2, theta){
  # n1,n2 -- sample size in each arm
  # theta -- hypothesis value
  
  # generate y1 and y2 from binomial distribution
  y1 <- rbinom(n=n1, size=1, prob=theta)
  y2 <- rbinom(n=n2, size=1, prob=theta)
  # compute sample means of y1 and y2
  p1 <- mean(y1)
  p2 <- mean(y2)
  # get the critical value to use in the hypothesis test
  p <- (n1*p1+n2*p2)/(n1+n2)
  z <- (p1-p2)/sqrt(p*(1-p)*(1/n1+1/n2))
  # get type 1,2 error
  alpha <- 1-pnorm(z,mean = 0,sd = 1)
  beta <- pnorm(z - qnorm(0.95) )
  
  return(c(alpha, beta))
}
# for example
error_estimation(n1=30, n2=60, theta=0.5)
error_estimation(n1=30, n2=60, theta=0.7)




evaluate_design <- function(lambda, gamma, n1, n2, theta) {
  
  # Estimate the expected sample size of a design defined by its 
  # decision rule parameters (lambda, gamma) and sample size 
  # parameters (n1, n2), along with its standard error.
  set.seed(28376)
  # Set the number of simulations.
  M <- 10^4
  # Create an empty vector to store simulated NS.
  Ns <- rep(NA, M)
  for (i in 1:M) {
    # Simulate theta from its prior, and then the stage 1 data conditional
    # on this theta.
    y1 <- rbinom(1, n1, theta)
    
    # Get posterior Beta(a1, b1) parameters.
    a1 <- 0.5 + y1
    b1 <- 0.5 + n1 - y1
    
    # Probability of futility.
    fut1 <- pbeta(0.5, a1, b1)
    
    # Threshold to determine progression, based on the decision rule.
    C1 <- 1 - lambda * (n1 / n2)^gamma
    
    # Note the final total sample size and store in the vector Ns.
    if (fut1 < C1) {
      Ns[i] <- n1
    } else {
      Ns[i] <- n2
    }
  }
  
  # Return the estimated expected sample size and its estimated standard error.
  return(c(mean(Ns), sqrt(var(Ns)/M)))
}
# For example,
evaluate_design(lambda=0.8, gamma=0.2, n1=30, n2=60, theta=0.5)


evaluate_design_exact <- function(lambda, gamma, n1, n2) {
  
  # Calculate the expected sample size of a design defined by its 
  # decision rule parameters (lambda, gamma) and sample size 
  # parameters (n1, n2), along with its standard error.
  
  # Threshold to determine progression, based on the decision rule.
  C1 <- 1 - lambda * (n1 / n2)^gamma
  
  # Vector of possible stage 1 outcomes.
  y_1s <- 0:n1
  
  # Vector of corresponding progression decisions.
  stops <- pbeta(0.5, y_1s + 0.5, n1 - y_1s + 0.5) < C1
  
  # For each outcome, calculate its probability.
  y_1_probs <- prob_y1(y_1s, n1)
  
  sum(n1 * stops * y_1_probs + n2 * (!stops) * y_1_probs)
}

# For example,
evaluate_design_exact(0.8, 0.2, 30, 60)




prob_y1 <- function(y1, n1) {
  # Calculate the probability of observing y_1 responses in n_1 trials
  # under a Beta-binomial model with Beta(a = 0.5, b = 0.5) prior.
  
  choose(n1, y1) * beta(y1 + 0.5, n1 - y1 + 0.5) / beta(0.5, 0.5)
}

test_prob_y1 <- function() {
  # Check that the probabilities sum to 1.
  n1 <- 30
  s <- sum(prob_y1(0:n1, n1))
  return(all.equal(s, 1))
}








