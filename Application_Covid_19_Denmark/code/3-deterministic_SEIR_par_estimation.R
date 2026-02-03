# in this file, we use UKF to estimate parameters with  
# 1. full dataset
# 2. dataset with the last observation removed one time, until 70 observations removed
# 3. dataset with the first observation removed one time, until 50 observations removed
# 4. dataset with observation removed randomly, the number of removal observations (mn): 1~30, 
#     for each mn, 50 incomplete datasets are generated

# about the random seeds in gosolnp
# we tried 10 different random seeds and select the best one in case 1,
# and use this random seeds for the rest scenarios 
# with a consideration of time costs

# ================================== Read Data =================================#

setwd('~/Application_Covid_19_Denmark/estimates')

df_full <- read_excel("Data Denmark.xlsx")
n_pop <- 5.83e6


start_date <- 225
end_date <- 405
df_infected <-df_full[225:405,]

# ============================================================================= #
#************************    Deterministic SEIR Model    ***********************#
# ============================================================================= #

sim.ode <- function(mu, x0, theta, delta_t, t_0 = 0, max_t) {
  t.vec <- seq(t_0, max_t, delta_t)
  n <- length(t.vec) - 1
  X.vec <- matrix(NA, ncol = n + 1, nrow = length(x0))
  X.vec[, 1] <- x0
  for (i in 2:(n + 1)) {
    X.vec[, i] <- X.vec[, i - 1] + mu(X.vec[, i - 1], theta) * delta_t
    X.vec[X.vec[, i] < 0, i] <- 0
  }
  return(list(state = X.vec, time = t.vec))
}


mu.seir.model <- function(x,theta) {
  return(c(-theta[1]*x[1]*x[3],
           theta[1]*x[1]*x[3]-theta[2]*x[2],
           theta[2]*x[2]-theta[3]*x[3]))
}

                                  
sum.of.squares.seir <- function(theta, time, infected_obs, r0) {
  beta <- abs(theta[1])
  sigma <- abs(theta[2])
  gamma <- abs(theta[3])
  e0 <- 0.5*abs(theta[4])
  i0 <- abs(theta[4])
  s0 <- 1 - i0 - e0 - r0
  
  if (any(c(s0, i0, e0) < 0) || s0 > 1) {
    return(Inf)
  }
  
  sim <- sim.ode(
    mu = mu.seir.model,
    x0 = c(s0, e0, i0),
    theta = c(beta, sigma, gamma),
    delta_t = 0.01,
    t_0 = 0,
    max_t = max(time)
  )
  
  i.sim <- sim$state[3, ]
  sim.time <- sim$time
  indices <- which(round(sim.time, 2) %in% round(time, 2))
  i.sim <- i.sim[indices]
  
  ssr <- sum((i.sim - infected_obs)^2)
  return(log(ssr))
}


# LSE - SIR model fitting  

N <-  n_pop  
r0 <- sum(df_full$New_cases[1:(start_date-8)])/N # approximate r0

n_obs <- length(df_infected$approx_infectious_cases)
t.obs <- seq(0,(n_obs-1),1)
observations <- df_infected$approx_infectious_cases/N

to.minimize.ls.seir <- function(theta) {
  sum.of.squares.seir(theta, time = t.obs, infected_obs = observations, r0 = r0)
}

lower_bounds <- c(0.1, 0.1, 0.1, 0)
upper_bounds <- c(3,   3,   3,   1e-1)

opt <- gosolnp(fun = to.minimize.ls.seir, LB = lower_bounds, UB = upper_bounds, 
               n.restarts = 5, n.sim = 500, rseed = 409,
               control = list(trace = T, tol = 1e-6)) 
est <- opt$pars
#9.984568e-01 8.952011e-01 8.709919e-01 4.669617e-05

#==============================================================================#
# delete tail/head

par_list_seir = matrix(NA, nrow = 50, ncol = 4) # delete head
# par_list_seir = matrix(NA, nrow = 70, ncol = 4) # delete tail

for (i in 1:nrow(par_list_seir)){
  
  # if delete head
  start_date_d <- start_date - 1 + i
  end_date_d <- end_date
  
  # if delete tail
  #start_date_d <- start_date
  #end_date_d <- end_date + 1 - i
  
  df_infected <- df_full[start_date_d:end_date_d,]
  
  #======= Initialize data =======#
  
  N <-  n_pop  
  r0 <- sum(df_full$New_cases[1:(start_date-8)])/N  # approximate r0
  
  n_obs <- length(df_infected$approx_infectious_cases)
  t.obs <- seq(0,(n_obs-1),1)
  observations <- df_infected$approx_infectious_cases/N
  
  to.minimize.ls.seir <- function(theta) {
    sum.of.squares.seir(theta, time = t.obs, infected_obs = observations, r0 = r0)
  }
  
  lower_bounds <- c(0.1, 0.1, 0.1, 0)
  upper_bounds <- c(3,   3,   3,   1e-1)
  
  opt <- gosolnp(fun = to.minimize.ls.seir, LB = lower_bounds, UB = upper_bounds, 
                 n.restarts = 5, n.sim = 500) 
  
  par_list_seir[i,]= abs(opt$pars)
  
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#           
#-------- DELETING POINTS RANDOMLY     
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

N <-  n_pop  
i0 <- df_infected$approx_infectious_cases[1]/N
r0 <- sum(df_full$New_cases[1:(start_date-8)])/N

n_obs0 <- length(df_infected$approx_infectious_cases)
t.obs0 <- seq(0,(n_obs0-1),1)
y0 <- df_infected$approx_infectious_cases/N

for (mn in 1:30) {
  
  n.missing <- mn
  list_seir_missing <- read.csv(paste("list_seir_missing.",n.missing,".csv",sep=''))
  
  par_list_seir_missing = matrix(NA, nrow = 50, ncol = 4)
  
  for (i in 1:nrow(par_list_seir_missing)){
    
    #======= Initialize data =======#
    ind.missing <- as.numeric(list_seir_missing[i,])
    
    y <-  y0[-ind.missing]
    t.obs <- t.obs0[-ind.missing]
    
    #======
    to.minimize.ls.seir <- function(theta) {
      sum.of.squares.seir(theta, time = t.obs, infected_obs = y, r0 = r0)
    }
    
    
    lower_bounds <- c(0.1, 0.1, 0.1, 1e-5)
    upper_bounds <- c(3,   3,  3,   1e-1)
    
    opt <- gosolnp(fun = to.minimize.ls.seir, LB = lower_bounds, UB = upper_bounds, 
                   n.restarts = 5, n.sim = 500)
    
    par_list_seir_missing[i,]= abs(opt$pars)
    
    #write.csv(par_list_seir_missing, paste("par_list_seir_",n.missing,"missing_ode.csv",sep=''),row.names=FALSE)
    
  }
  
}
