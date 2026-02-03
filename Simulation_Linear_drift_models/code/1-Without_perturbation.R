
library(ggplot2)   
library(reshape2)  
library(printr)
library(pracma)
library(rootSolve)
library(nleqslv)
library(gridExtra)
library(grid)
library(lattice)
library(tidyr)
library(tidyverse)
library(deSolve)
library(patchwork)
library(latex2exp)
theme_set(theme_bw()) 

#==============================================================================#
# 1d without perturbation
# simulation can be done from line 26 to 174
# the number of trajectories $n$ can be changed
# from line 180, plots can be done directly by reading files
# rho basic model a 0_05.csv; b basic model a 0_05.csv
# rho basic model a 0_1.csv; b basic model a 0_1.csv

# a = 0.05 or 0.1
a0 <- 0.1
b0 <- 0
x0 <- 5
sigma0 <- 0.05  # real sigma0 in ODE 
sigma <- sqrt(2 * a0* sigma0^2)

## simulation
dt <- 0.01
n <- 1e4

## "ode-ode","ode-sde","sde-sde","sde-ode"
rho_list <- matrix(NA, nrow = 3*n, ncol = 5)
b_list <- matrix(NA, nrow = 3*n, ncol = 5)

# choices: long T long delta; 
#          short T long delta; 
#          short T short delta
for (c in 1:3) { 
  
  if (c == 1){ # long T long delta
    T0 <- 100
    delta <- 2
    t  <- seq(delta,T0,delta) 
    N <- length(t)
    
    times <- seq(from=0,to=T0, by=dt)
  }
  else if(c == 2){ # short T long delta
    T0 <- 50
    delta <- 1
    t  <- seq(delta,T0,delta)
    N <- length(t)
    
    times <- seq(from=0,to=T0, by=dt)
  }
  else if(c == 3){ # short T short delta
    T0 <- 50
    delta <- 2
    t  <- seq(delta,T0,delta)
    N <- length(t)
    
    times <- seq(from=0,to=T0, by=dt)
  }
  
  for (i in 1:n) {
    
    x_ode <- c(x0)
    y_ou <- c(x0)
    
    for (l in 2:length(times)){
      x_ode[l] <- x_ode[l-1] -a0*(x_ode[l-1]-b0)*dt 
      y_ou[l] <- y_ou[l-1] -a0*(y_ou[l-1]-b0)*dt +  sigma * rnorm(1,0,sqrt(dt))
    }
    
    ## data ode
    yt_ode <- x_ode[(t/dt+1)] + rnorm(N, mean = 0, sd = sigma0)
    yt_ode_lag <- c(x0,yt_ode[-N])
    
    # data sde
    yt_sde <- y_ou[(t/dt+1)]
    yt_sde_lag <- c(x0,yt_sde[-N])
    
    ## ode lse
    Y <- matrix(yt_ode,nrow = N, ncol = 1)
    expX <- matrix(exp(-t),nrow = N, ncol = 1)
    
    fmatrix <- function(x){
      #for a
      #matrix(rep(1,N),nrow = 1) %*%(Y -(x0-x[2])*expX^x[1]-x[2])^2
      
      #for rho
      matrix(rep(1,N),nrow = 1) %*%(Y -(x0-x[2])*x[1]^(t/delta)-x[2])^2
    }
    initial1 <- c(exp(-x0*delta/T0),yt_ode[N])
    sol1 <- optim(initial1 ,fmatrix, method = "Nelder-Mead")
    
    rho_list[(c-1)*n+i,1]<- sol1$par[1] /exp(-a0*delta)
    b_list[(c-1)*n+i,1]<- sol1$par[2]
    
    ## ode mle
    estimates <- function(v) {     # v[1]: rho_hat; v[2]: b_hat
      eq <- numeric(2)
      eq[1] <- mean(yt_ode) + v[1]/N/(1-v[1])*(yt_ode[N]-yt_ode_lag[1])- v[2] # b_hat
      eq[2] <- sum((yt_ode-v[2])*(yt_ode_lag-v[2]))/sum((yt_ode_lag-v[2])^2) - v[1] # rho_hat
      eq 
    }
    
    vstart <- c(0.5,0) # initial values
    sol2 <- nleqslv(vstart, estimates)
    b_hat <- sol2$x[2]
    rho_hat <- sol2$x[1] #-log(sol2$x[1])/delta
    
    rho_list[(c-1)*n+i,2]<- rho_hat /exp(-a0*delta)
    b_list[(c-1)*n+i,2]<- b_hat
    
    ## sde mle
    estimates2 <- function(v) {     # v[1]: rho_hat; v[2]: b_hat
      eq <- numeric(2)
      eq[1] <- mean(yt_sde) + v[1]/N/(1-v[1])*(yt_sde[N]-yt_sde_lag[1])-v[2] # b_hat
      eq[2] <- sum((yt_sde-v[2])*(yt_sde_lag-v[2]))/sum((yt_sde_lag-v[2])^2) - v[1] # rho_hat
      eq 
    }
    
    vstart <- c(0.5,0) # c(exp(-x0*delta/T0),yt_sde[N]) # initial values
    sol3 <- nleqslv(vstart, estimates2)
    b_hat <- sol3$x[2]
    rho_hat <- sol3$x[1] #-log(sol3$x[1])/delta
    
    rho_list[(c-1)*n+i,3]<- rho_hat /exp(-a0*delta)
    b_list[(c-1)*n+i,3]<- b_hat
    
    ## sde ols
    Y2 <- matrix(yt_sde,nrow = N, ncol = 1)
    expX2 <- matrix(exp(-t),nrow = N, ncol = 1)
    
    fmatrix2 <- function(x){
      #for a
      #matrix(rep(1,N),nrow = 1) %*%(Y -(x0-x[2])*expX^x[1]-x[2])^2
      
      #for rho
      matrix(rep(1,N),nrow = 1) %*%(Y2 -(x0-x[2])*x[1]^(t/delta)-x[2])^2
    }
    initial4 <- c(exp(-x0*delta/T0),yt_sde[N])
    sol4 <- optim(initial4 ,fmatrix2, method = "Nelder-Mead")
    
    rho_list[(c-1)*n+i,4]<- sol4$par[1] /exp(-a0*delta)
    b_list[(c-1)*n+i,4]<- sol4$par[2]
    
    rho_list[(c-1)*n+i,5]<- paste("T =", T0, ", N =",T0/delta)
    b_list[(c-1)*n+i,5]<- paste("T =", T0, ", N =",T0/delta)
    
  }
}


df_rho <- data.frame(rho_list)
colnames(df_rho) <- c("ODE-ODE","ODE-SDE","SDE-SDE","SDE-ODE","Data_intensity")
df_rho <- melt(df_rho,id = c("Data_intensity"))
colnames(df_rho)[2] <- "Label"

df_b <- data.frame(b_list)
colnames(df_b) <-  c("ODE-ODE","ODE-SDE","SDE-SDE","SDE-ODE","Data_intensity")
df_b <- melt(df_b,id = c("Data_intensity"))
colnames(df_b)[2] <- "Label"

#write.csv(df_rho,"rho basic model a 0_05.csv")
#write.csv(df_b,"b basic model a 0_05.csv")

#write.csv(df_rho,"rho basic model a 0_1.csv")
#write.csv(df_b,"b basic model a 0_1.csv")

#==============================================================================#
setwd('Example 1/estimates')

# a = 0.1
df_rho <- read.csv("rho basic model a 0_1.csv")
df_b <- read.csv("b basic model a 0_1.csv")
# a = 0.05
df_rho <- read.csv("rho basic model a 0_05.csv")
df_b<- read.csv("b basic model a 0_05.csv")


p1 <- ggplot(data = df_rho,aes(x = as.numeric(value), color=Label))+
  geom_density(fill="grey",alpha=0.2, linewidth = 0.8)+
  geom_vline(aes(xintercept= 1 ),color="black", linetype="dashed", linewidth=0.5)+
  labs(title=" ",x= TeX("$\\hat{\\rho} / \\rho $"), y = " ")+
  facet_wrap(~ Data_intensity) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )  + xlim(c(0.985, 1.015)) # for a = 0.1
  #xlim(c(0.99, 1.01)) # for a = 0.05

p2 <- ggplot(data = df_b,aes(x = as.numeric(value), color=Label))+
  geom_density(fill= "grey",alpha=0.2, linewidth = 0.8)+
  geom_vline(aes(xintercept= 0),color="black", linetype="dashed", linewidth=0.5)+
  labs(title=" ",x=TeX("$\\hat{b} $"), y = " ")+
  facet_wrap(~ Data_intensity) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )  + xlim(c(-0.1, 0.1)) # for a = 0.1
  #xlim(c(-0.2, 0.2)) # for a = 0.05

(p1 / p2)  + plot_layout(guides = "collect")& theme(legend.position = "top") 


df_rho1_mean <- df_rho %>% group_by(Data_intensity, Label) %>% 
  summarise(mean_value = log(mean(as.numeric(value))), 
            .groups = 'drop') #%>% as.data.frame()
df_rho1_mean

df_b1_mean <- df_b %>% group_by(Data_intensity, Label) %>% 
  summarise(mean_value = mean(as.numeric(value)), 
            .groups = 'drop') #%>% as.data.frame()
df_b1_mean
