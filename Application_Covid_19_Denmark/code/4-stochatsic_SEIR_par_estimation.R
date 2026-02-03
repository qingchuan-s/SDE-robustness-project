# in this file, we use UKF to estimate parameters with  
# 1. full dataset
# 2. dataset with the last observation removed one time, until 70 observations removed
# 3. dataset with the first observation removed one time, until 50 observations removed
# 4. dataset with observation removed randomly, the number of removal observations (mn): 1~30, 
#     for each mn, 50 incomplete datasets are generated

# about the random seeds in gosolnp
# we tried 15 different random seeds and select the best one in case 1,
# and use this random seeds for the rest scenarios 
# with a consideration of time costs
# best seeds may differ on different computer, 

# also, we have so many parameters in this case, 
# we fixed the variance of meassurement error 
# and the diffusion of the state: e

# ================================== Read Data =================================#
df_full <- read_excel("Data Denmark.xlsx")
n_pop <- 5.83e6


start_date <- 225
end_date <- 405
df_infected <-df_full[225:405,]

# ============================================================================= #
#*************************    Stochastic SEIR Model    *************************#
# ============================================================================= #

f_drift <- function(x, theta) {
  S <- x[1]; E <- x[2]; I <- x[3]
  beta  <- theta[1]; alpha <- theta[2]; gamma <- theta[3]
  c(
    S - beta  * S * I,
    E + beta  * S * I - alpha * E,
    I + alpha * E     - gamma * I
  )
}

sigma_pts <- function(x, P, lambda) {
  n <- length(x)
  U <- chol((n + lambda) * P)            
  
  # make x a column matrix
  x  <- matrix(x, nrow = n, ncol = 1)
  
  # sigma-points:  x,  x + Ui,  x – Ui   (i = 1 … n)
  Xi_plus  <- sweep(U, 1, x, "+")        
  Xi_minus <- sweep(-U, 1, x, "+")       
  
  cbind(x, Xi_plus, Xi_minus)          
}

ukf_weights <- function(n, alpha = 1e-3, kappa = 0, beta = 2) {
  lambda <- alpha^2 * (n + kappa) - n
  Wm <- c(lambda / (n + lambda), rep(1 / (2 * (n + lambda)), 2 * n))
  Wc <- Wm;  Wc[1] <- Wc[1] + (1 - alpha^2 + beta)
  list(Wm = Wm, Wc = Wc, lambda = lambda)
}


negLogLik_UKF <- function(p, y, s0, e0, sigmas, tau,
                          alpha = 1e-3, kappa = 0, betaW = 2) {
  
  if (any(p <= 0) || any(p >= 10)) return(1e20)
  
  theta <- p
  n_obs <- length(y)
  
  R <- tau^2
  
  ## initial state & covariance
  x_hat <- c(S = s0, E = e0, I = y[1])
  P     <- diag(sigmas^2)          
  
  ## UKF weights
  n_state <- length(x_hat)
  w <- ukf_weights(n_state, alpha, kappa, betaW)
  
  ll <- 0
  for (t in 2:n_obs) {
    
    ## ---------- PREDICT ----------
    Xi <- sigma_pts(x_hat, P, w$lambda)          
    Xi_pred <- apply(Xi, 2, f_drift, theta = theta)  
    
    x_pred <- Xi_pred %*% w$Wm                  
    P_pred <- diag(sigmas^2)                     
    for (k in seq_len(ncol(Xi_pred)))
      P_pred <- P_pred + w$Wc[k] *
      tcrossprod(Xi_pred[, k] - x_pred)
    
    ## ---------- UPDATE ----------
    Yi <- Xi_pred[3, ]            
    y_pred <- sum(w$Wm * Yi)
    Svv <- R
    for (k in seq_along(Yi))
      Svv <- Svv + w$Wc[k] * (Yi[k] - y_pred)^2
    
    Cxy <- numeric(n_state)
    for (k in seq_along(Yi))
      Cxy <- Cxy + w$Wc[k] *
      (Xi_pred[, k] - x_pred) * (Yi[k] - y_pred)
    
    K <- Cxy / Svv                            
    innov <- y[t] - y_pred
    x_hat <- x_pred + K * innov
    P     <- P_pred - tcrossprod(K) * Svv
    
    ## ---------- log-likelihood ----------
    ll <- ll + 0.5 * (log(2 * pi) + log(Svv) + innov^2 / Svv)
    if (!is.finite(ll) || any(x_hat < 0) || any(x_hat > 1))
      return(1e20)
  }
  ll                       
}


N <-  n_pop  
i0 <- df_infected$approx_infectious_cases[1]/N      #approximate i0
r0 <- sum(df_full$New_cases[1:(start_date-8)])/N    #approximate r0
s0 <- 1 - 1.5*i0 - r0                               #approximate s0
e0 <- 0.5*i0                                        #approximate e0
Iobs <- df_infected$approx_infectious_cases/N


obj_fun <- function(theta, y, s0, e0) {
  
  sigmas <- c(theta[5], 1e-5, theta[4])   #sigma1, sigma2, sigma3
  tau    <- 1e-3                          # meas. noise
  
  negLogLik_UKF(theta[1:3],      
                y      = y,
                s0     = s0,
                e0     = e0,
                sigmas = sigmas,
                tau    = tau)
}

lower <- c(0.5, 0.1, 0.5, 1e-5, 1e-5)  
upper <- c(2,   2,   2,   1e-1, 1e-2)

fit <- gosolnp(
  fun  = obj_fun,
  LB   = lower,
  UB   = upper,
  y    = Iobs,
  s0   = s0,
  e0   = e0,
  n.restarts = 5,
  n.sim      = 500,
  rseed = 996,
  control = list(trace = 1, tol = 1e-6)
)

# 9.736867e-01 5.006845e-01 8.390785e-01 4.817349e-04 1.478588e-05

#==============================================================================#
# delete tail/head

par_list_seir = matrix(NA, nrow = 50, ncol = 5) # delete head
# par_list_seir = matrix(NA, nrow = 70, ncol = 5) # delete tail

start_date = 225; end_date = 405

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
  i0 <- df_infected$approx_infectious_cases[1]/N      
  r0 <- sum(df_full$New_cases[1:(start_date-8)])/N    
  s0 <- 1 - 1.5*i0 - r0                             
  e0 <- 0.5*i0                                       
  Iobs <- df_infected$approx_infectious_cases/N
  
  
  obj_fun <- function(theta, y, s0, e0) {
    
    sigmas <- c(theta[5], 1e-5, theta[4])   
    tau    <- 1e-3                         
    
    negLogLik_UKF(theta[1:3],      
                  y      = y,
                  s0     = s0,
                  e0     = e0,
                  sigmas = sigmas,
                  tau    = tau)
  }
  
  lower <- c(0.5, 0.1, 0.5, 1e-5, 1e-5)   
  upper <- c(2,   2,   2,   1e-1, 1e-1)
  
  fit <- gosolnp(
    fun  = obj_fun,
    LB   = lower,
    UB   = upper,
    y    = Iobs,
    s0   = s0,
    e0   = e0,
    n.restarts = 5,
    n.sim      = 500,
    control = list(trace = 1, tol = 1e-6)
  )
  
  par_list_seir[i,]= abs(fit$pars)
  
}

#==============================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# 1. Model: Drift function
f_drift <- function(x, theta) {
  S <- x[1]; E <- x[2]; I <- x[3]
  beta  <- theta[1]; alpha <- theta[2]; gamma <- theta[3]
  c(
    S - beta  * S * I,
    E + beta  * S * I - alpha * E,
    I + alpha * E     - gamma * I
  )
}

## ────────────────────────────────────────────────────────────────
##  2.  Unscented-transform utilities
## ────────────────────────────────────────────────────────────────
sigma_pts <- function(x, P, lambda) {
  n <- length(x)
  U <- chol((n + lambda) * P)            # n × n  “square root” of cov
  
  # make x a column matrix
  x  <- matrix(x, nrow = n, ncol = 1)
  
  # sigma-points:  x,  x + Ui,  x – Ui   (i = 1 … n)
  Xi_plus  <- sweep(U, 1, x, "+")        # add x to each column of  U
  Xi_minus <- sweep(-U, 1, x, "+")       # add x to each column of -U
  
  cbind(x, Xi_plus, Xi_minus)            # n × (2n+1)
}

ukf_weights <- function(n, alpha = 1e-3, kappa = 0, beta = 2) {
  lambda <- alpha^2 * (n + kappa) - n
  Wm <- c(lambda / (n + lambda), rep(1 / (2 * (n + lambda)), 2 * n))
  Wc <- Wm;  Wc[1] <- Wc[1] + (1 - alpha^2 + beta)
  list(Wm = Wm, Wc = Wc, lambda = lambda)
}


negLogLik_UKF <- function(p, y, s0, e0, sigmas, tau, t.obs,
                          alpha = 1e-3, kappa = 0, betaW = 2) {
  
  # p = (beta, gamma); reject impossible values
  if (any(p <= 0) || any(p >= 10)) return(1e20)
  theta  <- p
  
  n_obs <- length(y)
  if (length(t.obs) != n_obs)
    stop("length(t.obs) must match length(y)")
  
  Rmeas  <- tau^2
  
  # initial state at first obs-time
  x_hat <- c(S = s0, E = e0, I = y[1])
  P     <- diag(sigmas^2)
  
  n_state <- length(x_hat)
  w       <- ukf_weights(n_state, alpha, kappa, betaW)
  
  ll <- 0
  for (i in 2:n_obs) {
    
    dt <- t.obs[i] - t.obs[i-1]
    if (dt < 1) stop("t.obs must be strictly increasing integers")
    
    Xi <- sigma_pts(x_hat, P, w$lambda)
    
    for (step in seq_len(dt)) {
      Xi <- apply(Xi, 2, f_drift, theta = theta)
    }
    
    ## recombine
    x_pred <- Xi %*% w$Wm
    
    ## rebuild P_pred (you may also wish to inflate for process noise)
    P_pred <- diag(sigmas^2)
    for (k in seq_len(ncol(Xi))) {
      P_pred <- P_pred + w$Wc[k] * tcrossprod(Xi[,k] - x_pred)
    }
    
    ## update (we observe only I)
    Yi     <- Xi[3, ]
    y_pred <- sum(w$Wm * Yi)
    
    Svv <- Rmeas
    for (k in seq_along(Yi)) {
      Svv <- Svv + w$Wc[k] * (Yi[k] - y_pred)^2
    }
    
    Cxy <- numeric(n_state)
    for (k in seq_along(Yi)) {
      Cxy <- Cxy + w$Wc[k] * (Xi[,k] - x_pred) * (Yi[k] - y_pred)
    }
    
    K     <- Cxy / Svv
    innov <- y[i] - y_pred
    x_hat <- x_pred + K * innov
    P     <- P_pred - tcrossprod(K) * Svv
    
    ## accumulate negative log-likelihood
    ll <- ll + 0.5 * (log(2*pi) + log(Svv) + innov^2 / Svv)
    if (!is.finite(ll) || any(x_hat < 0) || any(x_hat > 1))
      return(1e20)
  }
  
  ll
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# 4. Loop over missing data scenarios
N <-  n_pop  
i0 <- df_infected$approx_infectious_cases[1]/N
r0 <- sum(df_full$New_cases[1:(start_date-8)])/N

n_obs0 <- length(df_infected$approx_infectious_cases)
t.obs0 <- seq(0, n_obs0-2, 1)
y0 <- df_infected$approx_infectious_cases[-1]/ N

for (mn in 1:30) {
  
  n.missing <- mn
  list_seir_missing <- read.csv(paste0("list_seir_missing.", n.missing, ".csv"))
  
  par_list_seir_missing <- matrix(NA, nrow = 50, ncol = 5)
  for (i in 1:nrow(par_list_seir_missing)) { 
    
    ind.missing <- as.numeric(list_seir_missing[i,])
    
    s0 <- 1-1.5*i0 - r0
    e0 <- 0.5*i0
    Iobs <- y0
    
    Iobs_incomp <- Iobs[-ind.missing]
    t.obs_incomp <- t.obs0[-ind.missing]
    
    obj_fun <- function(theta, y, s0, e0) {
      
      sigmas <- c(theta[5], 1e-5, theta[4])   # σS, σE (tiny), σI
      tau    <- 1e-3                          # meas. noise (fixed)
      
      negLogLik_UKF(theta[1:3],      # beta, alpha, gamma
                    y      = Iobs_incomp,
                    s0     = s0,
                    e0     = e0,
                    sigmas = sigmas,
                    tau    = tau,
                    t.obs  = t.obs_incomp)
    }
    
    lower <- c(0.5, 0.1, 0.5, 1e-5, 1e-5)   # beta, α, γ, σI, σS
    upper <- c(2,   1,   2,   1e-1, 1e-2)
    
    fit <- gosolnp(
      fun  = obj_fun,
      LB   = lower,
      UB   = upper,
      y    = Iobs_incomp,
      s0   = s0,
      e0   = e0,
      n.restarts = 5,  n.sim = 500, 
      control = list(trace = 1, tol = 1e-6)
    )
    
    par_list_seir_missing[i, ] <- abs(fit$pars) # par_list_seir_missing[i, ] <- c()
    print(paste0("Finished iteration ", i, " for n.missing = ", n.missing))
    print(abs(fit$pars))
    write.csv(par_list_seir_missing, paste0("par_list_seir_", n.missing, "missing_sde.csv"), row.names = FALSE)
    
  }
  
  write.csv(par_list_seir_missing, paste0("par_list_seir_", n.missing, "missing_sde.csv"), row.names = FALSE)
  
}

