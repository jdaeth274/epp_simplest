
  functions{

  //CODING UP THE simple SI model from the Bao et al 2012 paper 
  
  real[] hiv_SI(real t,
  real[] y,
  real[] params,
  real[] x_r,
  int[] x_i) {
  real dydt[4];
  
  dydt[1] = (params[1] * 5000) - (params[2] * y[2] * y[1] / y[3] ) - (params[3] * y[1] ) - ((params[4] * 2000) * y[1] / y[3] ) + ((params[5] * 2000) * y[1] / y[3] );
  
  dydt[2] = (params[2] * y[2] * y[1] / y[3] ) - ((params[3] + params[6]) * y[2] ) - ((params[4] * 2000) * y[2] / y[3] ) + ((params[5] * 2000) * y[2] / y[3] );
  
  dydt[3] = (params[1] * 5000) + (params[5] * 2000) - (params[4] * 2000) - (params[3] * (y[1] + y[2] ) + params[6] * y[2] );
  
  dydt[4] = (( dydt[2] + y[2]) / ( dydt[3] + y[3] )) - y[4];
  
  return dydt;
  
  
  
  }
  }

data {
  int<lower = 1> n_obs; // Number of days sampled
  int<lower = 1> n_params; // Number of model parameters
  int<lower = 1> n_difeq; // Number of differential equations in the system
  int<lower = 1> n_sample; // Number of hosts sampled at each time point.
  int<lower = 1> n_fake; // This is to generate "predicted"/"unsampled" data
  
  int y[n_obs]; // The binomially distributed data
  real t0; // Initial time point (zero)
  real ts[n_obs]; // Time points that were sampled
  
  real fake_ts[n_fake]; // Time points for "predicted"/"unsampled" data
}

transformed data {
  real x_r[0];
  int x_i[0];
}


parameters {
  real<lower = 0> params[n_params];
  real<lower = 0> Y0; // Initial fraction of hosts susceptible
  }
  
transformed parameters{
  real y_hat[n_obs, n_difeq]; // Output from the ODE solver
  real y0[n_difeq]; // Initial conditions for both S and I
  real y_prevalence[n_obs, 1]; // This is a relic from when I didn't explicitly model prevalence.

y0[1] = 1000000 - Y0;
y0[2] = Y0;
y0[3] = 1000000;   //This could be wrong, if flagged up set to 1 million
y0[4] = Y0 / 1000000; //This is the starting prevalence in the system

y_hat = integrate_ode_rk45(hiv_SI, y0, t0, ts, params, x_r, x_i);

for (t in 1:n_obs){

y_prevalence[t,1] = y_hat[t, 2] / y_hat[t, 3];

}

}

model {
  params ~ uniform(0, 1); //constrained at the lower limit of all the paramaters and the upper limit of the et paramter
  Y0 ~ normal(0, 100 ); //constrained to be positive.
  
  y ~ binomial(n_sample, y_hat[, 4]  ); //y_hat[,4] is the prevalence as a fraction of the total population
  
  }
generated quantities {
  // Generate predicted data over the whole time series:
real fake_I[n_fake, n_difeq];



fake_I = integrate_ode_rk45(hiv_SI, y0, t0, fake_ts, params, x_r, x_i);

}


