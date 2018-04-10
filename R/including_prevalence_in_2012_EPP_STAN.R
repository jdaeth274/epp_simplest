### fitting to the model with prevalence as a dydt #############################################################################

require(dplyr)
require(ggplot2)
require(deSolve)
require(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


Y0<-10          # initial fraction infected
Z0<-(1000000-Y0) # initial number susceptible
N0<-(Z0+Y0)  # Total population number intially
P0<-(Y0/N0)



# Assign transmission and pathogen-induced death rates:

et<- 2000  #25000 ## the birth rate of the population
r<-  0.5 #0.3    ## The average infection risk
mu<- 1/1000 ## The death rate in the general population
aged<- 100 ## The age out rate 
mig<- 100 ## the rate of net migration into the population
aided<- 1/10 ## the AIDS related death rate

# We will use the package deSolve to integrate, which requires certain data structures.
# Store parameters and initial values
# Parameters must be stored in a named list.
params_hiv_epp <- list(et=et,r=r,mu=mu,aged=aged,mig=mig,aided=aided)

# Initial conditions are stored in a vector
inits_hiv_prev <- c(Z0,Y0,N0,P0)

# Create a time series over which to integrate.
# Here we have an epidemic that is observed over t_max number of days (or weeks or etc).
t_min = 0
t_max = 50
times = t_min:t_max

# We must create a function for the system of ODEs.
# See the 'ode' function documentation for further insights.
SIR_prevalence <- function(t, y, params_hiv_epp) {
  with(as.list(c(params_hiv_epp, y)), {
    
    dZ = et - ((r*y[2]*y[1])/y[3])-(mu*y[1]) - ((aged*y[1])/y[3]) + ((mig*y[1])/y[3])
    
    dY = ((r*y[2]*y[1])/y[3]) - ((mu+aided)*y[2]) - ((aged*y[2])/y[3]) + ((mig*y[2])/y[3])
    
    dN = et + mig - aged - (mu*(y[1]+y[2])+(aided*y[2])) 
    
    dP = (((dY + y[2]) / (dN + y[3])) - y[4]) 
    
    
    res <- c(dZ,dY,dN,dP)
    list(res)
  })
}

# Run the integration:
out_hiv_epp_prev <- ode(inits_hiv_prev, times, SIR_prevalence, params_hiv_epp, method="ode45")

# Store the output in a data frame:
out_hiv_epp_prev <- data.frame(out_hiv_epp_prev)
colnames(out_hiv_epp_prev) <- c("time", "Z", "Y", "N","prevalence")
#out_hiv_epp$prevalence<-(out_hiv_epp$Y/out_hiv_epp$N)*100
out_hiv_epp_prev$prev_percent<-out_hiv_epp_prev$prevalence * 100

prev_plot<-ggplot(data = out_hiv_epp_prev, aes(x=time,y=prev_percent))+geom_line(colour="blue")+
  labs(x="Time",y="Prevalence (%)")

total_infected_plot<-ggplot(data = out_hiv_epp_prev,aes(x=time,y=Y))+geom_line(colour="red")+
  labs(x="Time",y="Total Number Infected")

total_pop_plot<-ggplot(data = out_hiv_epp_prev,aes(x=time,y=N))+geom_line(colour="forest green")+
  labs(x="Time",y="Total Population size")

require(gridExtra)
require(grid)

grid::grid.draw(rbind(ggplotGrob(prev_plot), ggplotGrob(total_infected_plot), ggplotGrob(total_pop_plot),  size = "last"))

#########################################################################################################################
## So that's the model coded up in R, now we will take samples from that data set and use them to fit our next model ####
#########################################################################################################################

sample_years_hiv = 50 # number of days sampled throughout the epidemic
sample_n = 100 # number of host individuals sampled per day

# Choose which days the samples were taken. 
# Ideally this would be daily, but we all know that is difficult.
sample_time_hiv = sort(sample(1:t_max, sample_years_hiv, replace=F))

# Extract the "true" fraction of the population that is infected on each of the sampled days:
sample_propinf_hiv = out_hiv_epp_prev[out_hiv_epp$time %in% sample_time_hiv, 5]

## this just samples our prevalence, to get a probability that the sample we take is HIV infected then we need to divide
## by 100

#sample_propinf_hiv<-sample_propinf_hiv/100

# Generate binomially distributed data.
# So, on each day we sample a given number of people (sample_n), and measure how many are infected.
# We expect binomially distributed error in this estimate, hence the random number generation.
sample_y_hiv_prev = rbinom(sample_years_hiv, sample_n, sample_propinf_hiv)
sample_prev_hiv<-(sample_y_hiv_prev/sample_n)*100

## lets have a ggplot of the y (infected) and out sample of Y over time 
sample_df_100<-data.frame(cbind(sample_time_hiv,sample_prev_hiv))
sample_df_100

ggplot(data = out_hiv_epp_prev,aes(x=time,y=prev_percent))+geom_line(colour="midnightblue",size=1.2)+
  geom_point(data=sample_df, aes(x=sample_df$sample_time_hiv,y=sample_df$sample_prev_hiv),colour="red",size=1)
ggplot(data = sample_df,aes(x=sample_time_hiv,y=sample_prev_hiv))+geom_point(colour="red",size=1.5)


### Now running the Stan model on the prevalence

stan_d_hiv_prev = list(n_obs = sample_years_hiv,
                  n_params = length(params_hiv_epp),
                  n_difeq = length(inits_hiv_prev),
                  n_sample = sample_n,
                  n_fake = length(1:t_max),
                  y = sample_y_hiv_prev,
                  t0 = 0,
                  ts = sample_time_hiv,
                  fake_ts = c(1:t_max))

# Which parameters to monitor in the model:
params_monitor_hiv = c("y_hat", "y0", "params", "fake_I")

# Test / debug the model:
test_hiv_100_year_dataprev_uniform_param_200_sample_n = stan("hiv_project/stan_c_scripts/including_prev_stan_2012_EPP.stan",
                data = stan_d_hiv_prev,
                pars = params_monitor_hiv,
                chains = 1, iter = 10)

## Run the model 

mod_hiv_prev_uniform_param_200_sample_n = stan(fit = test_hiv_prev_uniform_param_200_sample_n,
               data = stan_d_hiv_prev,
               pars = params_monitor,
               chains = 3,
               warmup = 500,
               iter = 1500,
               control = list(adapt_delta = 0.85))
# Extract the posterior samples to a structured list:
posts_hiv <- extract(mod_hiv_prev_uniform_param_200_sample_n)

apply(posts_hiv$params, 2, median)


apply(posts_hiv$y0, 2, median)[1:4]

# These should match well. 

#################
# Plot model fit:

# Proportion infected from the synthetic data:

#sample_prop = sample_y / sample_n

# Model predictions across the sampling time period.
# These were generated with the "fake" data and time series.
mod_median = apply(posts_hiv$fake_I[,,2], 2, median)
mod_low = apply(posts_hiv$fake_I[,,2], 2, quantile, probs=c(0.025))
mod_high = apply(posts_hiv$fake_I[,,2], 2, quantile, probs=c(0.975))
mod_time = stan_d_hiv$fake_ts

prev_median<-(apply(posts_hiv$fake_I[,,4],2,median))*100
prev_low<-(apply(posts_hiv$fake_I[,,4],2,quantile,probs=c(0.025)))*100
prev_high<-(apply(posts_hiv$fake_I[,,4],2,quantile,probs=c(0.975)))*100

# Combine into two data frames for plotting
#df_sample = data.frame(sample_prop, sample_time)
df_fit_prevalence = data.frame(prev_median, prev_low, prev_high, mod_time)

# Plot the synthetic data with the model predictions
# Median and 95% Credible Interval

n_100_plot<-ggplot(sample_df_100, aes(x=sample_df$sample_time_hiv, y=sample_df$sample_prev_hiv)) +
  geom_point(col="red", shape = 19, size = 1.5) +
  geom_line(data = df_fit_prevalence, aes(x=mod_time,y=prev_median),colour="midnightblue",size=1)+
  geom_ribbon(data = df_fit_prevalence,aes(x=mod_time,ymin=prev_low,ymax=prev_high),
              colour="midnightblue",alpha=0.2,fill="midnightblue")+
coord_cartesian(ylim = c(0,100),xlim=c(0,50))+labs(x="Time",y="Prevalence (%)", title="N = 100 plot")
plot(n_50_plot)

grid.arrange(rbind(ggplotGrob(n_50_plot),ggplotGrob(n_100_plot),ggplotGrob(n_200_plot), size="last"))
