
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
# 1d with randomly varying long term mean
# simulation can be done from line 26 to 183
# the number of trajectories $n$ can be changed
# from line 185, plots can be done directly by reading files
# rho oscillation a 0_05.csv; b oscillation a 0_05.csv
# rho oscillation a 0_1.csv; b oscillation a 0_1.csv

a0 <- 0.05 # a = 0.05, 0.1
b0 <- 0
x0 <- 5
sigma0 <- 0.05 # real sigma0 in ODE 
sigma <- sqrt(2*a0* sigma0^2)

sigmap_list <- c(0,0.5,1,2)

## observations
T0 <- 100
delta <- 2
t  <- seq(delta,T0,delta) 
N <- length(t)

## simulation
dt <- 0.01
times <- seq(from=0,to=T0, by=dt)

n <-1e4

## "ode-ode","ode-sde","sde-sde","sde-ode"
rho_list <- matrix(NA, nrow = length(sigmap_list )*n, ncol = 5)
b_list <- matrix(NA, nrow = length(sigmap_list )*n, ncol = 5)

for (i in 1:n) {
  
  r_dynamic <-  rnorm(length(times),0,1)
  r_b <- rnorm(length(times),0,1)
  
  x_ode <- c(x0)
  y_ode <- c(x0+ sigma0*r_dynamic[1])#
  
  y_ou <- c(x0)
  
  for (nh in 1:length(sigmap_list)){
    sigmap <- sigmap_list[nh]
    
    for (l in 2:length(times)){
      
      if(nh == 1){
        b = b0
      }else{
        b <- b0 +sigmap * r_b[l]
      } 
      
      x_ode[l] <- x_ode[l-1] -a0*(x_ode[l-1]-b)*dt 
      y_ode[l] <- x_ode[l] + sigma0 * r_dynamic[l]
      y_ou[l] <- y_ou[l-1] -a0*(y_ou[l-1]-b)*dt +  sigma*sqrt(dt)*r_dynamic[l]
      
    }
    
    yt_ode <- y_ode[(t/dt+1)]
    yt_ode_lag <- c(y_ode[1],yt_ode[-N])
    
    yt_sde <- y_ou[(t/dt+1)]
    yt_sde_lag <- c(x0,yt_sde[-N])
    
    ## ode lse
    Y <- matrix(yt_ode,nrow = N, ncol = 1)
    #expX <- matrix(exp(-t),nrow = N, ncol = 1)
    
    fmatrix <- function(x){
      #for a
      #matrix(rep(1,N),nrow = 1) %*%(Y -(x0-x[2])*expX^x[1]-x[2])^2
      
      #for rho
      matrix(rep(1,N),nrow = 1) %*%(Y -(x0-x[2])*x[1]^(t/delta)-x[2])^2
    }
    
    initial1 <- c(exp(-x0*delta/T0),yt_ode[N])
    sol1 <- optim(initial1 ,fmatrix, method = "BFGS")
    
    rho_list[(nh-1)*n+i,1]<- sol1$par[1]
    b_list[(nh-1)*n+i,1]<- sol1$par[2]
    
    ## ode mle
    estimates <- function(v) {     # v[1]: rho_hat; v[2]: b_hat
      eq <- numeric(2)
      eq[1] <- mean(yt_ode) + v[1]/N/(1-v[1])*(yt_ode[N]-yt_ode_lag[1])-v[2] # b_hat
      eq[2] <- sum((yt_ode-v[2])*(yt_ode_lag-v[2]))/sum((yt_ode_lag-v[2])^2) - v[1] # rho_hat
      eq 
    }
    vstart <- c(0.5,0) # initial values
    sol2 <- nleqslv(vstart, estimates)
    b_hat <- sol2$x[2]
    rho_hat <- sol2$x[1] #-log(sol2$x[1])/delta
    
    rho_list[(nh-1)*n+i,2]<- rho_hat 
    b_list[(nh-1)*n+i,2]<- b_hat
    
    ## sde mle
    estimates2 <- function(v) {     # v[1]: rho_hat; v[2]: b_hat
      eq <- numeric(2)
      eq[1] <- mean(yt_sde) + v[1]/N/(1-v[1])*(yt_sde[N]-yt_sde_lag[1])-v[2] # b_hat
      eq[2] <- sum((yt_sde-v[2])*(yt_sde_lag-v[2]))/sum((yt_sde_lag-v[2])^2) - v[1] # rho_hat
      eq 
    }
    vstart <- c(0.5,yt_sde[N]) # initial values
    sol3 <- nleqslv(vstart, estimates2)
    b_hat <- sol3$x[2]
    rho_hat <- sol3$x[1] #-log(sol3$x[1])/delta
    
    rho_list[(nh-1)*n+i,3]<- rho_hat 
    b_list[(nh-1)*n+i,3]<- b_hat
    
    ## sde ols
    Y2 <- matrix(yt_sde,nrow = N, ncol = 1)
    #expX2 <- matrix(exp(-t),nrow = N, ncol = 1)
    
    fmatrix2 <- function(x){
      #for a
      #matrix(rep(1,N),nrow = 1) %*%(Y -(x0-x[2])*expX^x[1]-x[2])^2
      
      #for rho
      matrix(rep(1,N),nrow = 1) %*%(Y2 -(x0-x[2])*x[1]^(t/delta)-x[2])^2
    }
    initial4 <- c(exp(-x0*delta/T0),yt_sde[N])
    sol4 <- optim(initial4 ,fmatrix2, method = "BFGS")
    
    rho_list[(nh-1)*n+i,4]<- sol4$par[1] 
    b_list[(nh-1)*n+i,4]<- sol4$par[2]
    
    rho_list[(nh-1)*n+i,5]<- paste("sigma_p=", sigmap)
    b_list[(nh-1)*n+i,5]<- paste("sigma_p=", sigmap)
  }
}


df_rho1 <- data.frame(rho_list)
colnames(df_rho1) <- c("ODE-ODE","ODE-SDE","SDE-SDE","SDE-ODE","Perturbation")
df_rho1 <- melt(df_rho1, id.vars = "Perturbation", variable.name = 'label')


df_b1 <- data.frame(b_list)
colnames(df_b1) <- c("ODE-ODE","ODE-SDE","SDE-SDE","SDE-ODE","Perturbation")
df_b1 <- melt(df_b1, id.vars = "Perturbation", variable.name = 'label')


#write.csv(df_rho1, "rho oscillation a 0_1 original.csv", row.names=FALSE)
#write.csv(df_b1, "b oscillation a 0_1  original.csv", row.names=FALSE)

df_rho1$index <- rep(1:1e4,16)
df_rho <- reshape(df_rho1, idvar = c("index","label"),
                  timevar = "Perturbation", direction = "wide")
df_rho$diff1 <- df_rho[,4] / df_rho[,3]
df_rho$diff2 <- df_rho[,5] / df_rho[,3]
df_rho$diff3 <- df_rho[,6] / df_rho[,3]
df_rho <- df_rho[,c("label" ,"diff1","diff2","diff3")]
df_rho <- melt(df_rho, id.vars = "label", variable.name = "Perturbation")

#write.csv(df_rho, "rho oscillation a 0_05.csv", row.names=FALSE)
#write.csv(df_rho, "rho oscillation a 0_1.csv", row.names=FALSE)

df_b1$index <- rep(1:1e4,16)
df_b <- reshape(df_b1, idvar = c("index","label"),
                  timevar = "Perturbation", direction = "wide")
df_b$diff1 <- df_b[,4] - df_b[,3]
df_b$diff2 <- df_b[,5] - df_b[,3]
df_b$diff3 <- df_b[,6] - df_b[,3]
df_b <- df_b[,c("label" ,"diff1","diff2","diff3")]
df_b <- melt(df_b, id.vars = "label", variable.name = "Perturbation")

#write.csv(df_b, "b oscillation a 0_05.csv", row.names=FALSE)
#write.csv(df_b, "b oscillation a 0_1.csv", row.names=FALSE)
#=============================================================================#
setwd('Example 1/estimates')

df_rho1 <- read.csv("rho oscillation a 0_05.csv")
df_b1 <- read.csv("b oscillation a 0_05.csv")
df_rho1 <- read.csv("rho oscillation a 0_1.csv")
df_b1 <- read.csv("b oscillation a 0_1.csv")

df_rho1$value <- as.numeric(df_rho1$value)
df_b1$value <- as.numeric(df_b1$value)

df_rho1_plot <- df_rho1
df_rho1_plot$data.model <- ifelse(df_rho1_plot$label %in% c("ODE-ODE","ODE-SDE"),"ODE","SDE")
df_rho1_plot$fitting.model <- ifelse(df_rho1_plot$label %in% c("ODE-SDE","SDE-SDE"),"SDE","ODE")

df_b1_plot <- df_b1
df_b1_plot$data.model <- ifelse(df_b1_plot$label %in% c("ODE-ODE","ODE-SDE"),"ODE","SDE")
df_b1_plot$fitting.model <- ifelse(df_b1_plot$label %in% c("ODE-SDE","SDE-SDE"),"SDE","ODE")


p1 <- ggplot(data = df_rho1_plot ,aes(x = value, color=Perturbation))+
  geom_density(fill="grey",alpha=0.2, linewidth = 0.8)+
  geom_vline(aes(xintercept= 1),color="black", linetype="dashed", size=0.5) +
  scale_color_manual(labels = c(TeX("$\\sigma_p = 0.5$ "), TeX("$\\sigma_p = 1$ "),TeX("$\\sigma_p = 2$ ")),values = c(3,6,4))+
  labs(title=" ",x= TeX("$\\hat{\\rho} / \\hat{\\rho}_0 $"), y = " ")+  theme(legend.position = "none")+
  facet_grid(fitting.model ~ data.model)+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15),
    strip.text.y = element_blank()
  )

p2 <- ggplot(data = df_b1_plot,aes(x = as.numeric(value), color=Perturbation))+
  geom_density(fill= "grey",alpha=0.2, linewidth = 0.8)+
  geom_vline(aes(xintercept= 0),color="black", linetype="dashed", size=0.5)+
  scale_color_manual(labels = c(TeX("$\\sigma_p = 0.5$ "), TeX("$\\sigma_p = 1$ "),TeX("$\\sigma_p = 2$ ")),values = c(3,6,4))+
  labs(title=" ",x=TeX("$\\hat{b} - \\hat{b}_0$"), y = " ")+
  facet_grid(fitting.model ~ data.model)+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )
p1<- p1 + 
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Data-generating model", breaks = NULL, labels = NULL))
p2<- p2+ 
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Fitted model", breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Data-generating model", breaks = NULL, labels = NULL))


(p1 | p2)  + plot_layout(guides = "collect")& theme(legend.position = "top") 


################################################################################
df_rho1_mean <- df_rho1 %>% group_by(Perturbation, label) %>% 
  summarise(mean_value = log(mean(as.numeric(value))/exp(-a0*delta)), 
            .groups = 'drop') #%>% as.data.frame()
df_rho1_mean

df_rho1_var <- df_rho1 %>% group_by(Perturbation, label) %>% 
  summarise(var_value = var(as.numeric(value)), 
            .groups = 'drop') #%>% as.data.frame()
df_rho1_var

df_b1_mean <- df_b1 %>% group_by(Perturbation, label) %>% 
  summarise(mean_value = mean(as.numeric(value)), 
            .groups = 'drop') #%>% as.data.frame()
df_b1_mean

df_b1_var <- df_b1 %>% group_by(Perturbation, label) %>% 
  summarise(mean_value = var(as.numeric(value)), 
            .groups = 'drop') #%>% as.data.frame()
df_b1_var