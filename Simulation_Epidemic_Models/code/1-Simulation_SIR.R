
###################### Simulate data of SIR model
################### SIR-ode
################### SIR-sde
################### Output: S, I

#----------------------------------------------------------------------------#
#------------------------------   Simulate data  ----------------------------#
#----------------------------------------------------------------------------#

theta <- c(0.5, 0.3, 3e-3, 1e-3) # alpha, beta, sigma1, sigma2
x0 <- c(0.99,0.001) # s0, i0

s_star <- uniroot(function(x) theta[1]/theta[2]*(1-x) + log(x/x0[1]), c(1e-6, 1-1e-6))$root

t0 <- 0
max_t <- 40

delta_t <- 0.01
delta1 <- 0.5
time1 <- seq(0,max_t,delta1)

t.vec <- seq(t0,max_t,delta_t)
indices1 <- which(round(t.vec,2) %in% round(time1,2))

n <- length(t.vec)-1

drift <- function(x,theta){
  mu <- numeric(2)
  mu[1] <-  - theta[1]*x[1]*x[2]
  mu[2] <- theta[1]*x[1]*x[2]- theta[2]*x[2]
  return(mu)
}

diffusion <- function(x,theta,random1, random2){
  sigma <- numeric(2)
  sigma[1] <-  theta[3]*sqrt(delta_t)*random1 
  sigma[2] <- theta[4]*sqrt(delta_t)*random2 
  return(sigma)
}

N_sample <- 2500

data_ode1 <- matrix(NA,ncol=2*N_sample,nrow=length(time1))
data_sde1 <- matrix(NA,ncol=2*N_sample,nrow=length(time1))

df_ode1<- data.frame(matrix(NA,ncol=4*N_sample,nrow=length(time1)-1))
df_sde1 <- data.frame(matrix(NA,ncol=4*N_sample,nrow=length(time1)-1))

list_stop <- c()
for (n_sample in 1:N_sample) {
  
  # vector X which stores the path. Rows: components of X; Columns: time points
  X1_sir.vec <- matrix(NA,ncol=n+1,nrow=length(x0))  #ode
  y1_sir.vec <- matrix(NA,ncol=n+1,nrow=length(x0))  #ode + measurement error
  X2_sir.vec <- matrix(NA,ncol=n+1,nrow=length(x0))  #sde
  
  # first column contains starting value
  
  X1_sir.vec[,1] <- x0
  y1_sir.vec[,1] <- x0
  X2_sir.vec[,1] <- x0
  
  for (i in 2:(n+1)) {
    random1 <- rnorm(1,0,1)
    random2 <- rnorm(1,0,1)
    
    X1_sir.vec[,i] <- X1_sir.vec[,i-1] + drift(X1_sir.vec[,i-1],theta)*delta_t
    y1_sir.vec[1,i] <- X1_sir.vec[1,i] + (theta[3] * sqrt(max_t)) *random1
    y1_sir.vec[2,i] <- X1_sir.vec[2,i] + (theta[4] / sqrt(2*(theta[2]-theta[1]*s_star))) *random2
    
    y1_sir.vec[1,i] <- ifelse(y1_sir.vec[1,i] >=1, 1, y1_sir.vec[1,i])
    y1_sir.vec[2,i] <- ifelse(y1_sir.vec[2,i] <=0, 0, y1_sir.vec[2,i])
    
    if (X2_sir.vec[2,i-1] > 0){
      X2_sir.vec[,i] <- X2_sir.vec[,i-1] + drift(X2_sir.vec[,i-1],theta)*delta_t + 
        diffusion(X2_sir.vec[,i-1],theta,random1, random2)
      
      if(X2_sir.vec[2,i] <= 0){
        X2_sir.vec[1,i] <-  X2_sir.vec[1,i-1]
        X2_sir.vec[2,i] = 0 
        list_stop <- append(list_stop, n_sample)
      }
      
    }else {
      X2_sir.vec[1,i] <-  X2_sir.vec[1,i-1]
      X2_sir.vec[2,i] <- 0
    }
  }
  
  s.sim.ode1 <- y1_sir.vec[1,indices1] 
  i.sim.ode1 <- y1_sir.vec[2,indices1]
  
  s.sim.sde1 <- X2_sir.vec[1,indices1] 
  i.sim.sde1 <- X2_sir.vec[2,indices1] 
  
  
  # =========================================================================================== #
  ################################        Data Frame (ODE)        ###############################
  # =========================================================================================== #
  
  data_ode1[,(2*n_sample-1)] <- s.sim.ode1
  data_ode1[,2*n_sample] <- i.sim.ode1
  
  df_ode1[,(4*n_sample-3)] <- s.sim.ode1[-1] 
  df_ode1[,(4*n_sample-2)] <- i.sim.ode1[-1]
  df_ode1[,(4*n_sample-1)] <- s.sim.ode1[-nrow(data_ode1)] 
  df_ode1[,4*n_sample] <- i.sim.ode1[-nrow(data_ode1)]
  
  # =========================================================================================== #
  ################################        Data Frame (SDE)        ###############################
  # =========================================================================================== #
  
  data_sde1[,(2*n_sample-1)] <- s.sim.sde1
  data_sde1[,2*n_sample] <- i.sim.sde1
  
  
  df_sde1[,(4*n_sample-3)] <- s.sim.sde1[-1] 
  df_sde1[,(4*n_sample-2)] <- i.sim.sde1[-1]
  df_sde1[,(4*n_sample-1)] <- s.sim.sde1[-nrow(data_sde1)] 
  df_sde1[,4*n_sample] <- i.sim.sde1[-nrow(data_sde1)]
  
}

#******************************************************************************#
label_list <- setdiff(seq(1, N_sample, 1), list_stop)
N_sample_new <- 1000
label_cont <- label_list[1:N_sample_new]

#datasets for fitting ODE model
index_ode <- index_sde <- as.vector(rbind(2 * label_cont - 1, 2 * label_cont))

data_ode_sub <- data_ode1[, index_ode]
data_sde_sub <- data_sde1[, index_sde]

#write.csv(data_ode_sub, "data_ode_for_ode_sir_sub.csv", row.names=FALSE)
#write.csv(data_sde_sub, "data_sde_for_ode_sir_sub.csv", row.names=FALSE)

#datasets for fitting SDE model
index_df <- as.vector(sapply(label_cont, function(k) (4 * k - 3):(4 * k)))
df_ode_sub <- df_ode1[, index_df]
df_sde_sub <- df_sde1[, index_df]

#write.csv(df_ode_sub, "df_ode_for_sde_sir_sub.csv", row.names=FALSE)
#write.csv(df_sde_sub, "df_sde_for_sde_sir_sub.csv", row.names=FALSE)