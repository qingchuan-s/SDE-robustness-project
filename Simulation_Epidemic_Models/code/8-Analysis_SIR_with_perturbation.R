setwd('Simulation_Epidemic_Models/estimates')

mat_ode_lse0 <- read.csv('est_ode_lse_xptb.csv')
mat_sde_lse0 <- read.csv('est_sde_lse_xptb.csv')
mat_ode_mle0 <- read.csv('est_ode_mle_xptb.csv')
mat_sde_mle0 <- read.csv('est_sde_mle_xptb.csv')

mat_ode_lse1 <- read.csv('est_ode_lse_ptb.csv')
mat_sde_lse1 <- read.csv('est_sde_lse_ptb.csv')
mat_ode_mle1 <- read.csv('est_ode_mle_ptb.csv')
mat_sde_mle1 <- read.csv('est_sde_mle_ptb.csv')

mat_ode_lse_p1 <- cbind(mat_ode_lse0[,1],mat_ode_lse1[,1])
mat_ode_lse_p2 <- cbind(mat_ode_lse0[,2],mat_ode_lse1[,2])

mat_sde_lse_p1 <- cbind(mat_sde_lse0[,1],mat_sde_lse1[,1])
mat_sde_lse_p2 <- cbind(mat_sde_lse0[,2],mat_sde_lse1[,2])

mat_ode_mle_p1 <- cbind(mat_ode_mle0[,1],mat_ode_mle1[,1])
mat_ode_mle_p2 <- cbind(mat_ode_mle0[,2],mat_ode_mle1[,2])

mat_sde_mle_p1 <- cbind(mat_sde_mle0[,1],mat_sde_mle1[,1])
mat_sde_mle_p2 <- cbind(mat_sde_mle0[,2],mat_sde_mle1[,2])

mat_ode_lse_p1 <- data.frame(mat_ode_lse_p1)
mat_ode_lse_p2 <- data.frame(mat_ode_lse_p2)


N_observation = 1000

df_p1 <- mat_ode_lse_p1
df_p2 <- mat_ode_lse_p2
df_p3 <- mat_ode_lse_p1[1:N_observation,1:2]/ mat_ode_lse_p2[1:N_observation,1:2]
colnames(df_p1) <- colnames(df_p2) <- colnames(df_p3) <- c('V1','V2')
df_p1 <- rbind(df_p1, 
               mat_ode_mle_p1[1:N_observation,1:2], 
               mat_sde_lse_p1[1:N_observation,1:2], 
               mat_sde_mle_p1[1:N_observation,1:2])
df_p2 <- rbind(df_p2, 
               mat_ode_mle_p2[1:N_observation,1:2], 
               mat_sde_lse_p2[1:N_observation,1:2], 
               mat_sde_mle_p2[1:N_observation,1:2])
df_p3 <- rbind(df_p3, 
               mat_ode_mle_p1[1:N_observation,1:2]/mat_ode_mle_p2[1:N_observation,1:2], 
               mat_sde_lse_p1[1:N_observation,1:2]/mat_sde_lse_p2[1:N_observation,1:2], 
               mat_sde_mle_p1[1:N_observation,1:2]/mat_sde_mle_p2[1:N_observation,1:2])
colnames(df_p1) <- colnames(df_p2)  <- colnames(df_p3) <-  c("no perturbation","with perturbation")
df_p1$`Data-geneating model` <- df_p2$`Data-geneating model` <- df_p3$`Data-geneating model`<- c(rep('deterministic SIR model', 2*N_observation),rep('stochastic SIR model', 2*N_observation))
df_p1$`Fitted model` <- df_p2$`Fitted model` <- df_p3$`Fitted model` <- c(
  rep('deterministic SIR model', N_observation),rep('stochastic SIR model', N_observation),
  rep('deterministic SIR model', N_observation),rep('stochastic SIR model', N_observation))
df_p1 <- melt(df_p1, id.vars = c("Data-geneating model","Fitted model"))
df_p2 <- melt(df_p2, id.vars = c("Data-geneating model","Fitted model"))
df_p3 <- melt(df_p3, id.vars = c("Data-geneating model","Fitted model"))


df_p1_ode_lse <- df_p1 %>% filter(`Data-geneating model`=='deterministic SIR model' & `Fitted model`=='deterministic SIR model')
df_p1_ode_mle <- df_p1 %>% filter(`Data-geneating model`=='deterministic SIR model' & `Fitted model`=='stochastic SIR model')
df_p1_sde_lse <- df_p1 %>% filter(`Data-geneating model`=='stochastic SIR model' & `Fitted model`=='deterministic SIR model')
df_p1_sde_mle <- df_p1 %>% filter(`Data-geneating model`=='stochastic SIR model' & `Fitted model`=='stochastic SIR model')
df_p2_ode_lse <- df_p2 %>% filter(`Data-geneating model`=='deterministic SIR model' & `Fitted model`=='deterministic SIR model')
df_p2_ode_mle <- df_p2 %>% filter(`Data-geneating model`=='deterministic SIR model' & `Fitted model`=='stochastic SIR model')
df_p2_sde_lse <- df_p2 %>% filter(`Data-geneating model`=='stochastic SIR model' & `Fitted model`=='deterministic SIR model')
df_p2_sde_mle <- df_p2 %>% filter(`Data-geneating model`=='stochastic SIR model' & `Fitted model`=='stochastic SIR model')

df_p3_ode_lse <- df_p3 %>% filter(`Data-geneating model`=='deterministic SIR model' & `Fitted model`=='deterministic SIR model')
df_p3_ode_mle <- df_p3 %>% filter(`Data-geneating model`=='deterministic SIR model' & `Fitted model`=='stochastic SIR model')
df_p3_sde_lse <- df_p3 %>% filter(`Data-geneating model`=='stochastic SIR model' & `Fitted model`=='deterministic SIR model')
df_p3_sde_mle <- df_p3 %>% filter(`Data-geneating model`=='stochastic SIR model' & `Fitted model`=='stochastic SIR model')



p1 <- ggplot(data = df_p1_ode_lse)+
  geom_density(aes(x= value, color = variable),fill="grey",alpha = 0.2, linewidth = 0.8)+
  geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.5)+
  scale_color_manual("",labels = c('No perturbation', 'With perturbation'),
                     values = c(4,7))+
  xlim(0.38,0.80)+ ylim(c(0,100))+
  facet_grid( `Data-geneating model` ~ `Fitted model`) + 
  labs(x = NULL, y=NULL)+  
    theme(legend.position = "none",
          axis.text.x = element_blank(), 
          strip.text.y = element_blank())+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )
p2 <- ggplot(data = df_p1_ode_mle)+
  geom_density(aes(x= value, color = variable),fill="grey",alpha = 0.2, linewidth = 0.8)+
  geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.5)+
  scale_color_manual("",labels = c('No perturbation', 'With perturbation'),
                     values = c(4,7))+
  xlim(0.38,0.80)+ ylim(c(0,100))+
  facet_grid( `Data-geneating model` ~ `Fitted model`) +
  labs(x = NULL, y=NULL)+  
  theme(legend.position = "none",
        axis.text = element_blank(), 
        strip.text.y = element_blank())+
  theme(
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )
p3 <- ggplot(data = df_p1_sde_lse)+
  geom_density(aes(x= value, color = variable),fill="grey",alpha = 0.2, linewidth = 0.8)+
  geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.5)+
  scale_color_manual("",labels = c('No perturbation', 'With perturbation'),
                     values = c(4,7))+
  xlim(0.38,0.80)+ ylim(c(0,50))+
  facet_grid( `Data-geneating model` ~ `Fitted model`) +
  labs(x = NULL, y=NULL)+  
  theme(legend.position = "none",
        #axis.text.y= element_blank(), 
        strip.text = element_blank())+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12)
  )
p4 <- ggplot(data = df_p1_sde_mle)+
  geom_density(aes(x= value, color = variable),fill="grey",alpha = 0.2, linewidth = 0.8)+
  geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.5)+
  scale_color_manual("",labels = c('No perturbation', 'With perturbation'),
                     values = c(4,7))+
  xlim(0.38,0.80)+ ylim(c(0,50))+
  facet_grid( `Data-geneating model` ~ `Fitted model`) +
  labs(x = NULL, y=NULL)+  
  theme(legend.position = "none",
        axis.text.y= element_blank(), 
        strip.text = element_blank())+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12)
  )

p5 <- ggplot(data = df_p2_ode_lse)+
  geom_density(aes(x= value, color = variable),fill="grey",alpha = 0.2, linewidth = 0.8)+
  geom_vline(xintercept = 0.3, linetype = 2, linewidth = 0.5)+
  scale_color_manual("",labels = c('No perturbation', 'With perturbation'),
                     values = c(4,7))+
  xlim(0.2,0.5)+ ylim(c(0,100))+
  facet_grid( `Data-geneating model` ~ `Fitted model`) +
  labs(x = NULL, y=NULL)+  
  theme(legend.position = "none",
        axis.text = element_blank(), 
        strip.text.y = element_blank())+
  theme(
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )
p6 <- ggplot(data = df_p2_ode_mle)+
  geom_density(aes(x= value, color = variable),fill="grey",alpha = 0.2, linewidth = 0.8)+
  geom_vline(xintercept = 0.3, linetype = 2, linewidth = 0.5)+
  scale_color_manual("",labels = c('No perturbation', 'With perturbation'),
                     values = c(4,7))+
  xlim(0.2,0.5)+ ylim(c(0,100))+
  facet_grid( `Data-geneating model` ~ `Fitted model`) +
  labs(x = NULL, y=NULL)+  
  theme(legend.position = "none",
        axis.text = element_blank(), 
        strip.text.y = element_blank())+
  theme(
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )
p7 <- ggplot(data = df_p2_sde_lse)+
  geom_density(aes(x= value, color = variable),fill="grey",alpha = 0.2, linewidth = 0.8)+
  geom_vline(xintercept = 0.3, linetype = 2, linewidth = 0.5)+
  scale_color_manual("",labels = c('No perturbation', 'With perturbation'),
                     values = c(4,7))+
  xlim(0.2,0.5)+ ylim(c(0,50))+
  facet_grid( `Data-geneating model` ~ `Fitted model`) +
  labs(x = NULL, y=NULL)+  
  theme(legend.position = "none",
        axis.text.y= element_blank(), 
        strip.text = element_blank())+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12)
  )
p8 <- ggplot(data = df_p2_sde_mle)+
  geom_density(aes(x= value, color = variable),fill="grey",alpha = 0.2, linewidth = 0.8)+
  geom_vline(xintercept = 0.3, linetype = 2, linewidth = 0.5)+
  scale_color_manual("",labels = c('No perturbation', 'With perturbation'),
                     values = c(4,7))+
  xlim(0.2,0.5)+ ylim(c(0,50))+
  facet_grid( `Data-geneating model` ~ `Fitted model`) +
  theme(strip.text.x = element_blank())+
  labs(x = NULL, y=NULL)+  
  theme(legend.position = "none",
        axis.text.y= element_blank(), 
        strip.text = element_blank())+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12)
  )


p9 <- ggplot(data = df_p3_ode_lse)+
  geom_density(aes(x= value, color = variable),fill="grey",alpha = 0.2, linewidth = 0.8)+
  geom_vline(xintercept = 5/3, linetype = 2, linewidth = 0.5)+
  scale_color_manual("",labels = c('No perturbation', 'With perturbation'),
                     values = c(4,7))+
  xlim(1.5,1.8)+ ylim(c(0,100))+
  facet_grid( `Data-geneating model` ~ `Fitted model`) +
  labs(x = NULL, y=NULL)+  
  theme(legend.position = "none",
        axis.text= element_blank(), 
        strip.text.y = element_blank())+
  theme(
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )
p10 <- ggplot(data = df_p3_ode_mle)+
  geom_density(aes(x= value, color = variable),fill="grey",alpha = 0.2, linewidth = 0.8)+
  geom_vline(xintercept = 5/3, linetype = 2, linewidth = 0.5)+
  scale_color_manual("",labels = c('No perturbation', 'With perturbation'),
                     values = c(4,7))+
  xlim(1.5,1.8)+ ylim(c(0,100))+
  facet_grid( `Data-geneating model` ~ `Fitted model`) +
  labs(x = NULL, y=NULL)+  
  theme(legend.position = "none",
        axis.text = element_blank())+
  theme(
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )
p11 <- ggplot(data = df_p3_sde_lse)+
  geom_density(aes(x= value, color = variable),fill="grey",alpha = 0.2, linewidth = 0.8)+
  geom_vline(xintercept = 5/3, linetype = 2, linewidth = 0.5)+
  scale_color_manual("",labels = c('No perturbation', 'With perturbation'),
                     values = c(4,7))+
  xlim(1.5,1.8)+ ylim(c(0,100))+
  facet_grid( `Data-geneating model` ~ `Fitted model`) +
  labs(x = NULL, y=NULL)+  
  theme(legend.position = "none",
        axis.text.y= element_blank(), 
        strip.text = element_blank())+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12)
  )
p12 <- ggplot(data = df_p3_sde_mle)+
  geom_density(aes(x= value, color = variable),fill="grey",alpha = 0.2, linewidth = 0.8)+
  geom_vline(xintercept = 5/3, linetype = 2, linewidth = 0.5)+
  scale_color_manual("",labels = c('No perturbation', 'With perturbation'),
                     values = c(4,7))+
  xlim(1.5,1.8)+ ylim(c(0,100))+
  facet_grid( `Data-geneating model` ~ `Fitted model`) +
  labs(x = NULL, y=NULL)+  
  theme(legend.position = "none",
        axis.text.y= element_blank(), 
        strip.text.x = element_blank())+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 15)
  )

p_lab1 <- 
  ggplot() + 
  annotate(geom = "text", x = 1, y = 1, label = TeX("$\\hat{\\alpha}$ "), angle = 0) +
  coord_cartesian(clip = "off")+
  theme_void()
p_lab2 <- 
  ggplot() + 
  annotate(geom = "text", x = 2, y = 1, label = TeX("$\\hat{\\beta}$"), angle = 0) +
  coord_cartesian(clip = "off")+
  theme_void()

p_lab3 <- 
  ggplot() + 
  annotate(geom = "text", x = 3, y = 1, label = TeX("$\\hat{R}_0$"), angle = 0) +
  coord_cartesian(clip = "off")+
  theme_void()

p_lab4 <- 
  ggplot() + 
  annotate(geom = "text", x = 1, y = 1, label = "Data-generating Model", angle = -90) +
  coord_cartesian(clip = "off")+
  theme_void()

p_lab5 <- 
  ggplot() + 
  annotate(geom = "text", x = -1, y = 0, label = "Fitted model", angle = 0) +
  coord_cartesian(clip = "off")+
  theme_void()

p_lab6 <- 
  ggplot() + 
  annotate(geom = "text", x = 1, y = 1, label = " ", angle = 0) +
  coord_cartesian(clip = "off")+
  theme_void()

design <- "
AAAAAAAAAAAAAAAAAAAAAAAA#
BBBBDDDDFFFFHHHHJJJJLLLLQ
CCCCEEEEGGGGIIIIKKKKMMMMQ
NNNNNNNNOOOOOOOOPPPPPPPP#
"
p_lab5 + p1+p3+ p2+p4+p5+p7 +p6+p8+p9+p11 +p10+p12+p_lab1 + p_lab2 + p_lab3 + p_lab4 + p_lab6 +
  plot_layout(heights = c( 0.05,1,1,0.05),design = design,guides = "collect")& theme(legend.position = "top")



