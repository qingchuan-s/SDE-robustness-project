library(cowplot)

setwd("~/Application_Covid_19_Denmark/estimates")

df.SIR <- read.csv('SIR random obs removed.csv')
df.SEIR <- read.csv('SEIR random obs removed.csv')

df.SIR$model <- "SIR model"
df.SEIR$model <- "SEIR model"

df.full0 <- rbind(df.SIR, df.SEIR)


df.full0$Estimator <- recode(df.full0$Estimator,
                             "LSE" = "ODE",
                             "MLE" = "SDE")
df.full0$model <- factor(df.full0$model, levels = c("SIR model", "SEIR model"))
colnames(df.full0)[5] <- 'Model'
df.full0$par <- factor(df.full0$par, levels = c("alpha", "beta",'R', 'gamma'))
df.full0$par <- recode(df.full0$par, "R" = "R0")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

df.SIR <- read.csv('SIR tail obs removed.csv')
df.SEIR <- read.csv('SEIR tail obs removed.csv')

df.SIR$model <- "SIR model"
df.SEIR$model <- "SEIR model"

df.full1 <- rbind(df.SIR, df.SEIR)
colnames(df.full1) <- c('Estimates', 'Estimator',"day",  "date", "parameter", "model")

df.full1$Estimator <- recode(df.full1$Estimator,
                             "LSE" = "ODE",
                             "MLE" = "SDE")
df.full1$model <- factor(df.full1$model, levels = c("SIR model", "SEIR model"))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

df.SIR <- read.csv('SIR head obs removed.csv')
df.SEIR <- read.csv('SEIR head obs removed.csv')


df.SIR$model <- "SIR model"
df.SEIR$model <- "SEIR model"

df.full2 <- rbind(df.SIR, df.SEIR)#, df.SEIR
colnames(df.full2) <- c('Estimates', 'Estimator',"day",  "date", "parameter", "model")

df.full2$Estimator <- recode(df.full2$Estimator,
                             "LSE" = "ODE",
                             "MLE" = "SDE")
df.full2$model <- factor(df.full2$model, levels = c("SIR model", "SEIR model"))

#+================================================================================+#
require(ggplot2)
require(tidyr)
require(plyr)
library(ggpubr)

data.a = df.full1 
data.a$date = as.Date(data.a$date, format = "%Y-%m-%d")
data.a$date2 = as.POSIXct(data.a$date)

data.b = df.full2 
data.b$date = as.Date(data.b$date, format = "%Y-%m-%d")
data.b$date2 = as.POSIXct(data.b$date)

data.ab = rbind(data.a,data.b)

data.ab$Estimator = as.factor(data.ab$Estimator)
data.ab$parameter[data.ab$parameter == "R"] = "R0"
data.ab$parameter = ordered(data.ab$parameter, levels = c("alpha", "beta", "R0", "gamma"))

data.ab$model <- recode(data.ab$model,
                        "SEIR model" = "SEIR",
                        "SIR model" = "SIR")
data.ab$model = ordered(data.ab$model, levels = c("SIR", "SEIR"))

data.ab1 = subset(data.ab, date == "2020-09-12")
data.ab1$data3 = rep(as.Date("2020-12-01"),nrow(data.ab1))

data.c = df.full0
data.c$model <- recode(data.c$model,
                       "SEIR model" = "SEIR",
                       "SIR model" = "SIR")
data.c$par[data.c$par == "R"] = "R0"
colnames(data.c)[5] = 'Model'
data.c$model = as.factor(data.c$model)
data.c$model = ordered(data.c$model, levels = c("SIR", "SEIR"))
data.c$par = ordered(data.c$par, levels = c("alpha", "beta", "R0", "gamma"))

data.ab1dummy = data.ab1
data.ab1dummy$Estimates = c(rep(0.15,4), rep(1,2),
                            rep(1.2,2), 0, 1.3,
                            rep(1.3,2), rep(1.3,2))

data.cdummy = data.frame(Model = rep("ODE",8),
                         par = c(rep("alpha",2), rep("beta",2),rep("gamma",2),rep("R0",2)),
                         missing_points = rep(1,8),
                         mean = c(0.15, 1.5, 0.15, 1.3, 0, 1.2,1, 1.3))
data.cdummy$par = ordered(data.cdummy$par, levels = c("alpha", "beta", "R0", "gamma"))

data.ab$Estimates[data.ab$parameter == 'R0'] <- ifelse(data.ab$Estimates[data.ab$parameter == 'R0'] > 1.3, NA, data.ab$Estimates[data.ab$parameter == 'R0'])
p1 = ggplot(data = data.ab, aes(x = date, y = Estimates, color = Estimator)) +
  geom_point(size = 1.2) +
  geom_blank(data = data.ab1dummy) +
  geom_point(data = data.ab1, aes(x = data3), size = 3) +
  facet_grid(cols = vars(model), rows = vars(parameter), scales = "free_y") +
  xlab("Date") + ylab("Estimate") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 0, vjust = 0.5 ),
        axis.title = element_text(size = 20),
        legend.position = c(0.15,0.1),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.background = element_rect(fill="lightgray",
                                         size=0.5, linetype="solid",
                                         colour ="darkblue"),
        strip.text = element_text(size = 15)) +
  labs(color = "Model") + scale_color_manual(
    values = c(
      "ODE" = "blue",
      "SDE" = "red"
    )
  )

p2 = ggplot(data = data.c, aes(x = missing_points, y = mean, color = Model)) +
  geom_point(size = 1.5) +
  geom_line() +
  geom_blank(data = data.cdummy) +
  geom_errorbar(aes(ymin = min, ymax = max))  +
  xlab("Number of missing points") + ylab("") + #ylab("Estimate") +
  facet_grid(cols = vars(model), rows = vars(par), scales = "free")+
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 0, vjust = 0.5 ),
        axis.title = element_text(size = 20),
        legend.position = c(0.15,0.1),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.background = element_rect(fill="lightgray",
                                         size=0.5, linetype="solid",
                                         colour ="darkblue"),
        strip.text = element_text(size = 15)) +
  labs(color = "Model") +  scale_color_manual(
    values = c(
      "ODE" = "blue",
      "SDE" = "red"
    )
  )

ggarrange(p1, p2 + rremove("y.text"),
          labels = c("A", "B"),
          ncol = 2, nrow = 1)


ggsave("Covid19removalpoints1.png", width = 15, height = 10)
