#1. Sample 1000 times 10 observation from uniform distribution over [0. a] interval
#   for a given value of parameter 'a'. Use those samples to estimate parameter 'a' 
#   using MLE method. Visualize results on histogram.

a = 10
n_obs = 10
n_samples = 1000
Uniform_simulator = runif(n_obs*n_samples, max=a)
Uniform_matrix = matrix(Uniform_simulator,nrow=n_obs)
alpha = 0.01
lower.limits <- apply(Uniform_matrix, 2, max)
hist(lower.limits)


#2. Calculate the confidence interval for sample number 1 using the formula.

lower_limit <- max(Uniform_matrix[,1])
upper_limit <- max(Uniform_matrix[,1])/(alpha^(1/n_obs))

lower_limit
upper_limit

#3. Calculate the confidence interval for each sample and check for how many samples the range 
#   contains the true confidence interval parameter value a.

Uni_check <- as.data.frame(Uniform_matrix[1:3,])

for(i in 1:n_samples){
  dummy_lower_limit <- max(Uniform_matrix[,i])
  dummy_upper_limit <- max(Uniform_matrix[,i])/(alpha^(1/n_obs))
  Uni_check[1,i] <- dummy_lower_limit
  Uni_check[2,i] <- dummy_upper_limit
  if(a > dummy_lower_limit & a < dummy_upper_limit){
    Uni_check[3,i] <- 1
  }
  else{
    Uni_check[3,i] <- 0
  }
}

sum(Uni_check[3,]) # number of attempts whose interval contains the real parameter a
