//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int<lower=1> D;
  matrix[N, D] x;
  array[N] int<lower=0, upper=1> y;
  real<lower=0> sigma;
}

parameters {
  vector[D] beta;
}

model {
  beta ~ normal(0, sigma);
  y ~ bernoulli_logit(x * beta);
}
