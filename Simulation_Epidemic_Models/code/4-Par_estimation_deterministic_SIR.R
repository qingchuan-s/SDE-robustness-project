

# =========================================================================================== #
#**************************        Least Square Method (OLS)        **************************#
# =========================================================================================== #
sim.ode <- function(mu, x0, theta,  delta_t,  t_0=0,  max_t){
  t.vec <- seq(t_0,max_t,delta_t)
  n <- length(t.vec)-1
  X.vec <- matrix(NA,ncol=n+1,nrow=length(x0))
  X.vec[,1] <- x0
  for (i in 2:(n+1)) {
    X.vec[,i] <- X.vec[,i-1] + mu(X.vec[,i-1],theta)*delta_t
    X.vec[X.vec[,i]<0,i] <- 0  
  }
  return.list <- list(X.vec,t.vec)
  names(return.list) <- c("state","time")
  return(return.list)
}
##### drift function #####
mu.sir.model <- function(x,theta) {
  return(c(-theta[1]*x[1]*x[2],
           theta[1]*x[1]*x[2]-theta[2]*x[2]))
}

#   ==================   #
sum.of.squares.sir <- function(theta,time,susceptible, infected) {
  # check whether "time" and "infected" are of same length
  if (length(time)!=length(infected)) {
    stop("Parameters time and infected of different lengths.")
  }
  # simulate s and i from ODE with parameter theta
  sim <-  sim.ode(
    mu=mu.sir.model,
    x0=c(susceptible[1],infected[1]), 
    theta=theta,   
    delta_t=0.01,  
    t_0=0,    
    max_t=max(time)   
  )
  s.sim <- sim$state[1,]
  i.sim <- sim$state[2,]
  indices <- which(round(sim$time,2) %in% round(time,2))
  s.sim <- s.sim[indices]
  i.sim <- i.sim[indices]
  ssr <- sum((s.sim - susceptible)^2) + sum((i.sim - infected)^2)
  return(ssr)
}

#   ==================   #

to.minimize.ls.sir.int <- function(theta,time,susceptible,infected) {
  # simulate numbers for given theta
  sim <-  sim.ode(
    mu=mu.sir.model, 
    x0=c(susceptible[1],infected[1]),
    theta=theta,     
    delta_t=0.01,    
    t_0=0,         
    max_t=max(time) 
  )
  # sum of squared.residuals
  return(sum.of.squares.sir(theta,time,susceptible,infected))
}


#============================================================================================#
max_t <- 40
delta1 <- 0.5
time1 <- seq(0,max_t,delta1)

# Import data
ind <- seq(1,1000,1)
s_indices <- 2 * ind - 1
i_indices <- 2 * ind

# data_ode_sub: data_ode_for_ode_sir_sub, data_ode_for_ode_seir_sub, data_ode_for_ode_sir_sub_no_ptb, data_ode_for_ode_sir_sub_ptb
# data_sde_sub: data_sde_for_ode_sir_sub, data_sde_for_ode_seir_sub, data_sde_for_ode_sir_sub_no_ptb, data_sde_for_ode_sir_sub_ptb

# full data = data_ode_sub, data_sde_sub
data_sub_s <- full_data[, s_indices]
data_sub_i <- full_data[, i_indices]

mat_lse <-  matrix(NA,ncol=2,nrow=N_sample)

for (n_sample in 1:1000) {
  data_s <- data_sub_s[,n_sample]
  data_i <- data_sub_i[,n_sample]
  
  to.minimize.ls.sir <- function(theta) {
    to.minimize.ls.sir.int(theta,time= seq(0,max_t-d*delta1,delta1),
                           susceptible = data_s,
                           infected= data_i)
  }
  ols <- optim(par=c(0.6,0.4),
                   fn=to.minimize.ls.sir,
                   method="BFGS"    #L-BFGS-B"
  )
  
  mat_lse[n_sample,] <- ols$par
}
  
