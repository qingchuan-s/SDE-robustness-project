
# =========================================================================================== #
#**************************   Maximum Likelihood Estimation (MLE)   **************************#
# =========================================================================================== #

mat_A<- function(alpha, beta) {
  A <- torch_zeros(c(2, 2))
  A[1, 1] <- (-alpha)
  A[2, 1] <- alpha
  A[2, 2] <- (-beta)
  A
}
mat_sigma <- function(sd1,sd2) {
  S <- torch_zeros(c(2, 2))
  S[1, 1] <- sd1^2
  S[2, 2] <- sd2^2
  S
}

mat_Omega <- function(h, alpha, beta, sigma1, sigma2){
  M <- torch_zeros(c(4, 4))
  M[1:2, 1:2] <- mat_A(alpha, beta)
  M[1:2, 3:4] <- mat_sigma(sigma1, sigma2)
  M[3:4, 3:4] <- -torch_transpose(M[1:2, 1:2], 1, 2)
  exphM <- torch_matrix_exp(h * M)
  F1 <- exphM[1:2, 1:2]
  G1 <- exphM[1:2, 3:4]
  return(G1$matmul(F1$t()))
}

f_h <- function(N, x, h, alpha, beta){
  s <- torch_tensor(x[, 1])
  i <- torch_tensor(x[, 2])
  c1 <- s + i
  c2 <- (log((1-i)/s)) /(c1-1)
  f_h1 <- (c1-1)/(1-exp((alpha  *h +  c2)*(c1-1)))
  f_h2 <- c1 - f_h1
  return(torch_cat(c(f_h1,  f_h2 ))$reshape(c(2, N)))
}


d_f_h <- function(N, x, h, alpha){
  s <- torch_tensor(x[, 1])
  i <- torch_tensor(x[, 2])
  c1 <- s + i
  c2 <- log((1-i)/s)
  exp_term <- exp((c1-1)*alpha*h + c2)
  d_f_h1 <- 1/(1-exp_term) + (c1-1)/(1-exp_term)^2 *(exp_term *(alpha*h-1/s)) # (1,1)
  d_f_h2 <-  1-1/(1-exp_term) - (c1-1)/(1-exp_term)^2 *(exp_term *(alpha*h-1/s)) # (2,1)
  d_f_h3 <- 1/(1-exp_term) + (c1-1)/(1-exp_term)^2 *(exp_term *(alpha*h-1/(1-i))) # (1,2)
  d_f_h4 <- 1-1/(1-exp_term) - (c1-1)/(1-exp_term)^2 *(exp_term *(alpha*h-1/(1-i))) # (2,2)
  return(torch_cat(c(d_f_h1 , d_f_h2 ,d_f_h3 ,d_f_h4 ))$reshape(c(4, N)))
}


### Objective function

Negativ_log_lik_LT <- function(data, par, h){
 
  N <- length(data[,1])
  data_old <- data[,3:4]
  data_new <- data[,1:2]
  a <- par[1]
  b <- par[2]
  s1 <-par[3]
  s2 <-par[4]
  step <- h
  A <- mat_A(a, b)
  Omega <- mat_Omega(h, a, b, s1, s2)
  fh <- f_h(N,data_old, step, a, b)  # dim: 2 * 50
  error_term <-  t(data_new) - (A*step)$matrix_exp()$matmul(fh) # dim: 2 * 50
  negativ_log_likelihood <- (N * torch_log(Omega$det()) /2  +  0.5 * torch_trace(error_term$t()$matmul(Omega$inverse())$matmul(error_term)))
  return(negativ_log_likelihood)
}

Negativ_log_lik_S <- function(data, par, h){
  #par <- nnf_softplus(par_old, beta = 1)
  N <- length(data[,1])
  data_old <- data[,3:4]
  data_new <- data[,1:2]
  a <- par[1]
  b <- par[2]
  s1 <-par[3]
  s2 <-par[4]
  #h <- h
  step <- h/2
  step_inv <- -h/2
  A <- mat_A(a, b)
  Omega <- mat_Omega(h, a, b, s1, s2)
  fh <- f_h(N,data_old, step, a, b)  # dim: 2 * 50 # step = h/2
  mu <-  (A*h)$matrix_exp()$matmul(fh) # dim: 2 * 50
  fh_inv <- f_h(N,data_new, step_inv, a, b) # dim: 2 * 50
  z <-  fh_inv - mu # dim: 2 * 50
  dfh_inv <- d_f_h(N,data_new, step_inv, a) # dim: 4 * 50
  det_dfh_inv <- abs(dfh[1,] * dfh[4,] - dfh[2,] * dfh[3,])
  negativ_log_likelihood <- (N * torch_log(Omega$det()) /2  - sum(torch_log(det_dfh_inv)) +
                               0.5 * torch_trace(z$t()$matmul(Omega$inverse())$matmul(z)))
  return(negativ_log_likelihood)
}


### estimators

estimator <- function(obj, data, h, par_start, num_iterations = 1000) {
  par_start <- torch_tensor(par_start, requires_grad = TRUE)
  optimizer <- optim_rprop(par_start,lr = 1e-2, step_sizes = c(1e-7,50))
  calc_loss <- function() {
    optimizer$zero_grad()
    value <- obj(data,nnf_softplus(par_start, beta = 1), h)
    value$backward()
    value
  }
  for (i in 1:num_iterations) {
    par_old = as.matrix(par_start)
    optimizer$step(calc_loss)
    par_new = as.matrix(par_start)
    if(norm(par_new - par_old) < 10^-6) break
  }
  convergence <- 0
  if(i == num_iterations) convergence <- 1
  list(as.numeric(nnf_softplus(par_start, beta = 1)), convergence, i)
}



#============================================================================================#
max_t <- 40
delta1 <- 0.5
time1 <- seq(0,max_t,delta1)

# Import data
# df_ode: df_ode_for_sde_sir_sub, df_ode_for_sde_seir_sub, df_ode_for_sde_sir_sub_no_ptb, df_ode_for_sde_sir_sub_no_ptb
# df_sde: df_ode_for_sde_sir_sub, df_ode_for_sde_seir_sub, df_ode_for_sde_sir_sub_no_ptb, df_ode_for_sde_sir_sub_no_ptb

# df_sub = df_ode or df_sde
# mat_lse: using the same dataset but to fit a deterministic SIR model
mat_mle <-  matrix(NA,ncol=4,nrow=N_sample) 

for (n_sample in 1:1000) { 
  df_sample <- df_sub[,(4*n_sample-3) :(4*n_sample)]
  
  initial_value <- c(mat_lse[n_sample,1:2], 1e-2, 1e-3)
  est <- estimator(obj = Negativ_log_lik_S, data = df_sample, h=delta1 , par_start = initial_value)[[1]]
  
  mat_mle[n_sample,] <- est[1:4]

}




