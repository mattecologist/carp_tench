##### Carp/ Tench. Redfin data
library (ggplot2)
library (reshape2)

library(tidyverse)
library(lubridate)

# read in flow data
flow <- read.csv("./Data/Euston flows 1945 to 2008-mid Murray.csv")
flow <- mutate(flow, date = dmy(paste(day, month, year)))

ggplot(flow, aes(y = flow, x = date)) +
  geom_line() +
  theme_classic()

flow.yr <- group_by(flow, year) %>%
  dplyr::summarise(flow = sum(flow))  

ggplot(flow.yr, aes(y = flow, x = year)) +
  geom_line() +
  theme_classic()

# read in data from Dave
dat <- read.csv("./Data/Tench.csv")
glimpse(dat)
dat <- select(dat, n.fishers, tench.kg, carp.kg, redfin.kg, year, fish.days)

dat <- left_join(dat, flow.yr) 

reid_df <- mutate(dat, carp.cpu = carp.kg / n.fishers,
                  tench.cpu = tench.kg / n.fishers,
                  redfin.cpu = redfin.kg / n.fishers)%>%
  select(year, carp.cpu, tench.cpu, redfin.cpu, flow)

reid_df <- melt (reid_df, id.vars=c("year", "flow"))

ggplot (reid_df, aes(year, value, group=variable))+
  geom_point(aes(colour=variable)) +
  geom_smooth(method="glm", formula="y ~ poly(x, 3)",aes(colour=variable))+
  ggtitle("Count data")


ggplot (reid_df, aes((flow), (value), group=variable))+
  geom_point(aes(colour=variable)) +
  geom_smooth(method="glm", formula="y ~ poly(x, 1)",aes(colour=variable))+
  facet_wrap(~variable, scales ="free")+
  theme_minimal()

ggplot (reid_df, aes(year, log10(value+1), group=variable))+
  geom_point(aes(colour=variable)) +
  geom_smooth(method="glm", formula="y ~ poly(x, 1)",aes(colour=variable))+
  theme_minimal()

ggplot (reid_df, aes(year, log10(value+1), group=variable))+
  geom_point(aes(colour=variable)) +
  geom_smooth(method="glm", formula="y ~ poly(x, 2)",aes(colour=variable))+
  theme_minimal()+
  theme(axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text = element_text(size=12))+
  ggtitle("Logged population counts w/ 2nd degree polynomial")

ggplot(reid_df, aes(flow, value))+
  geom_boxplot(aes(fill=variable))

## Adapted from the example in the BUGS 
library (R2jags)
library (runjags)
library (mcmcplots)

carp.dat  <- dcast (reid_df, year ~ variable)
carp.dat <- carp.dat[complete.cases(carp.dat),]


# Specify model in BUGS language
sink("ssm.jags")
cat("
    model {
    # Priors and constraints
    logN.est[1] ~ dnorm(20, 0.01)       # Prior for initial population size
    mean.r ~ dnorm(0, 0.001)             # Prior for mean growth rate
    sigma.proc ~ dunif(0, 10)             # Prior for sd of state process
    #sigma2.proc <- sigma.proc^2
    tau.proc <- 1/(sigma.proc*sigma.proc)
    sigma.obs ~ dunif(0, 10)              # Prior for sd of observation process
    #sigma2.obs <- sigma.obs^2
    tau.obs <- 1/(sigma.obs*sigma.obs)
    
    # Likelihood
    # State process
    for (t in 1:(T-1)){
    r[t] ~ dnorm(mean.r, tau.proc)  ## estimate for the rate of increase at time [t]
    logN.est[t+1] <- logN.est[t] + r[t] ## estimates on a log scale
    }
    # Observation process
    for (t in 1:T) {
    y[t] ~ dnorm(logN.est[t], tau.obs)
    }
    
    # Population sizes on real scale
    for (t in 1:T) {
    N.est[t] <- exp(logN.est[t]) ## rescaled estimates from log data -> real scale
    }
    }
    ",fill = TRUE)
sink()


carp <- carp.dat$carp.cpu
year <- carp.dat$year

# Bundle data
jags.data <- list(y = log(carp), T = length(year))

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 10), mean.r = rnorm(1), sigma.obs = runif(1, 0, 10), logN.est = c(rnorm(1, 20, 0.1), rep(NA, (length(year)-1))))}

# Parameters monitored
parameters <- c("r", "mean.r", "sigma.obs", "sigma.proc", "N.est")

# MCMC settings
ni <- 2000
nt <- 6
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
carp.ssm <- jags(jags.data, inits, parameters, "ssm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(carp.ssm, digits = 3)


# Draw figure
fitted <- lower <- upper <- numeric()
year <- carp.dat$year
n.years <- length(year)
for (i in 1:n.years){
  fitted[i] <- mean(carp.ssm$BUGSoutput$sims.list$N.est[,i])
  lower[i] <- quantile(carp.ssm$BUGSoutput$sims.list$N.est[,i], 0.025)
  upper[i] <- quantile(carp.ssm$BUGSoutput$sims.list$N.est[,i], 0.975)}

plot.dat <- as.data.frame(cbind(year, fitted, lower, upper))

carp.plot <- ggplot (plot.dat, aes(year, fitted))+
  geom_line(colour="#440154FF", size=1, aes())+
  geom_ribbon(aes(ymin=lower, ymax=upper), fill= "#21908CFF", alpha=0.5)+
  geom_point(data=data.frame(cbind(plot.dat$year, carp)),  aes(year, carp), colour="#FDE725FF")+
  ylab("Population size")+
  xlab("Year")+
  theme_minimal()+
  theme(axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text = element_text(size=12))+
  ggtitle("Carp population estimate")


mcmcplot(carp.ssm)

### Redfin data

redfin <- carp.dat$redfin.cpu
year <- carp.dat$year

# Bundle data
jags.data <- list(y = log(redfin), T = length(year))

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 10), mean.r = rnorm(1), sigma.obs = runif(1, 0, 10), logN.est = c(rnorm(1, 20, 0.1), rep(NA, (length(year)-1))))}

# Parameters monitored
parameters <- c("r", "mean.r", "sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 2000; nt <- 6; nb <- 1000; nc <- 3

# Call JAGS from R (BRT 3 min)
redfin.ssm <- jags(jags.data, inits, parameters, "ssm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(redfin.ssm, digits = 3)


# Draw figure
fitted <- lower <- upper <- numeric()
year <- carp.dat$year
n.years <- length(year)
for (i in 1:n.years){
  fitted[i] <- mean(redfin.ssm$BUGSoutput$sims.list$N.est[,i])
  lower[i] <- quantile(redfin.ssm$BUGSoutput$sims.list$N.est[,i], 0.025)
  upper[i] <- quantile(redfin.ssm$BUGSoutput$sims.list$N.est[,i], 0.975)}

plot.dat <- as.data.frame(cbind(year, fitted, lower, upper))

redfin.plot <- ggplot (plot.dat, aes(year, fitted))+
  geom_line(colour="#440154FF", size=1, aes())+
  geom_ribbon(aes(ymin=lower, ymax=upper), fill= "#21908CFF", alpha=0.5)+
  geom_point(data=data.frame(cbind(plot.dat$year, redfin)),  aes(year, redfin), colour="#FDE725FF")+
  ylab("Population size")+
  xlab("Year")+
  theme_minimal()+
  theme(axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text = element_text(size=12))+
  ggtitle("Redfin population estimate")


### tench data

tench <- carp.dat$tench.cpu
year <- carp.dat$year


tench <- tench + 0.01

# Bundle data
jags.data <- list(y = log(tench), T = length(year))

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 10), mean.r = rnorm(1), sigma.obs = runif(1, 0, 10), logN.est = c(rnorm(1, 20, 0.1), rep(NA, (length(year)-1))))}

# Parameters monitored
parameters <- c("r", "mean.r", "sigma2.obs", "sigma2.proc", "N.est")

# MCMC settingsq
ni <- 2000; nt <- 6; nb <- 1000; nc <- 3

# Call JAGS from R (BRT 3 min)
tench.ssm <- jags(jags.data, inits, parameters, "ssm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(tench.ssm, digits = 3)


# Draw figure
fitted <- lower <- upper <- numeric()
year <- carp.dat$year
n.years <- length(year)
for (i in 1:n.years){
  fitted[i] <- mean(tench.ssm$BUGSoutput$sims.list$N.est[,i])
  lower[i] <- quantile(tench.ssm$BUGSoutput$sims.list$N.est[,i], 0.025)
  upper[i] <- quantile(tench.ssm$BUGSoutput$sims.list$N.est[,i], 0.975)}

plot.dat <- as.data.frame(cbind(year, fitted, lower, upper))

tench.plot <- ggplot (plot.dat, aes(year, fitted))+
  geom_line(colour="#440154FF", size=1, aes())+
  geom_ribbon(aes(ymin=lower, ymax=upper), fill= "#21908CFF", alpha=0.5)+
  geom_point(data=data.frame(cbind(plot.dat$year, tench)),  aes(year, tench), colour="#FDE725FF")+
  ylab("Population size")+
  xlab("Year")+
  theme_minimal()+
  theme(axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text = element_text(size=12))+
  ggtitle("tench population estimate")

cowplot::plot_grid(carp.plot, redfin.plot, tench.plot, ncol=1)


####################### With covariates

carp.dat  <- dcast (reid_df, year + flow  ~ variable)

carp.dat <- carp.dat[complete.cases(carp.dat),]

# Specify model in BUGS language
sink("ssm2.jags")
cat("
    model {
    # Priors and constraints
    logN.est[1] ~ dnorm(20, 0.01)       # Prior for initial population size
    mean.r ~ dnorm(0, 0.001)             # Prior for mean growth rate
    sigma.proc ~ dunif(0, 10)             # Prior for sd of state process
    sigma2.proc <- sigma.proc^2
    tau.proc <- 1/(sigma.proc*sigma.proc)
    sigma.obs ~ dunif(0, 10)              # Prior for sd of observation process
    sigma2.obs <- sigma.obs^2
    tau.obs <- 1/(sigma.obs*sigma.obs)
    beta ~ dnorm (0, 0.001)  ## Prior for coefficient used on covariate (Z)
    
    # Likelihood
    # State process

    for (t in 1:(T-1)){
    r[t] ~ dnorm(mean.r, tau.proc)  ## estimate for the rate of increase at time [t]
    logN.est[t+1] <- logN.est[t] + r[t] + beta*Z[t] ## estimates on a log scale ## Z[t] is the year before's flow as the logN.est is t+1
    }

    # Observation process
    for (t in 1:T) {
    y[t] ~ dnorm(logN.est[t], tau.obs)
    }
    
    # Population sizes on real scale
    for (t in 1:T) {
    N.est[t] <- exp(logN.est[t]) ## rescaled estimates from log data -> real scale
    }
    }
    ",fill = TRUE)
sink()


carp <- carp.dat$carp.cpu
year <- carp.dat$year
flow <- carp.dat$flow/1E6

# Bundle data
jags.data <- list(y = log(carp), T = length(year), Z = flow)

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 10), mean.r = rnorm(1), sigma.obs = runif(1, 0, 10), beta=rnorm(1), logN.est = c(rnorm(1, 20, 0.1), rep(NA, (length(year)-1))))}

# Parameters monitored
parameters <- c("r", "mean.r", "sigma2.obs", "sigma2.proc", "N.est", "beta")

# MCMC settings
ni <- 20000; nt <- 6; nb <- 1000; nc <- 3

# Call JAGS from R (BRT 3 min)
carp.ssm2 <- jags(jags.data, inits, parameters, "ssm2.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(carp.ssm2, digits = 3)


posterior <- carp.ssm2$BUGSoutput$sims.matrix

posterior <- as.data.frame(posterior)

posterior2 <- posterior[,c("mean.r","deviance","sigma2.obs" ,"sigma2.proc","beta")]

posterior2 <- melt(posterior2)

ggplot (posterior2, aes(value))+
geom_histogram(fill="#21908CFF", alpha=0.6)+
  geom_vline(xintercept=0, colour="#440154FF")+
  theme_minimal()+
  facet_wrap(~variable, scales="free")+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size=12),
        strip.text = element_text(size=12))+
  ggtitle("Posteriors")

# Draw figure
fitted <- lower <- upper <- numeric()
year <- carp.dat$year
n.years <- length(year)
for (i in 1:n.years){
  fitted[i] <- mean(carp.ssm2$BUGSoutput$sims.list$N.est[,i])
  lower[i] <- quantile(carp.ssm2$BUGSoutput$sims.list$N.est[,i], 0.025)
  upper[i] <- quantile(carp.ssm2$BUGSoutput$sims.list$N.est[,i], 0.975)}

plot.dat <- as.data.frame(cbind(year, fitted, lower, upper))

ggplot (plot.dat, aes(year, fitted))+
  geom_line(colour="#440154FF", size=1, aes())+
  geom_ribbon(aes(ymin=lower, ymax=upper), fill= "#21908CFF", alpha=0.5)+
  geom_point(data=data.frame(cbind(plot.dat$year, carp)),  aes(year, carp), colour="#FDE725FF")+
  ylab("Population size")+
  xlab("Year")+
  theme_minimal()+
  theme(axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text = element_text(size=12))+
  ggtitle("Carp population estimate")


####################### With covariates AND Tench

carp.dat  <- dcast (reid_df, year + flow  ~ variable)

carp.dat <- carp.dat[complete.cases(carp.dat),]

# Specify model in BUGS language
sink("ssm2.jags")
cat("
    model {
    # Priors and constraints
    logN.est[1] ~ dnorm(20, 0.01)       # Prior for initial population size
    mean.r ~ dnorm(0, 0.001)             # Prior for mean growth rate
    sigma.proc ~ dunif(0, 10)             # Prior for sd of state process
    sigma2.proc <- sigma.proc^2
    tau.proc <- 1/(sigma.proc*sigma.proc)
    sigma.obs ~ dunif(0, 10)              # Prior for sd of observation process
    sigma2.obs <- sigma.obs^2
    tau.obs <- 1/(sigma.obs*sigma.obs)
    # flow prior
    beta ~ dnorm (0, 0.001)  ## Prior for coefficient used on covariate (Z)
    
    ## Tench priors
    t.logN.est[1] ~ dnorm(20, 0.01)
    t.mean.r ~ dnorm(0, 0.001)
    t.sigma.proc ~ dunif(0, 10)
    t.sigma2.proc <- t.sigma.proc^2
    t.tau.proc <- 1/(t.sigma.proc*t.sigma.proc)
    t.sigma.obs ~ dunif(0, 10)             
    t.sigma2.obs <- t.sigma.obs^2
    t.tau.obs <- 1/(t.sigma.obs*t.sigma.obs)
    t.beta ~ dnorm(0, 0.001)
 
    # Likelihood
    # State process
    
    for (t in 1:(T-1)){
   #TENCH
    t.r[t] ~ dnorm(t.mean.r, t.tau.proc)  ## estimate for the rate of increase at time [t]
    t.logN.est[t+1] <- t.logN.est[t] + t.r[t] + t.beta*Z[t] 
    #CARP
    r[t] ~ dnorm(mean.r, tau.proc)  ## estimate for the rate of increase at time [t]
    logN.est[t+1] <- logN.est[t] + r[t] + beta*Z[t] ## estimates on a log scale ## Z[t] is the year before's flow as the logN.est is t+1
 
  
    }
    
    # Observation process
    for (t in 1:T) {
    #CARP   
    y[t] ~ dnorm(logN.est[t], tau.obs)
    #TENCH
    t.y[t] ~ dnorm(t.logN.est[t], t.tau.obs)
    }
    
    # Population sizes on real scale
    for (t in 1:T) {
    #CARP
    N.est[t] <- exp(logN.est[t]) ## rescaled estimates from log data -> real scale
    #TENCH
    t.N.est[t] <- exp(t.logN.est[t]) ## rescaled estimates from log data -> real scale
    }
    
 
    


    }
    ",fill = TRUE)
sink()


carp <- carp.dat$carp.cpu
tench <- carp.dat$tench.cpu
year <- carp.dat$year
flow <- carp.dat$flow/1E6

# Bundle data
jags.data <- list(y = log(carp), t.y = log(tench+1), T = length(year), Z = flow)

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 10), mean.r = rnorm(1), sigma.obs = runif(1, 0, 10), beta=rnorm(1), t.beta=rnorm(1), logN.est = c(rnorm(1, 20, 0.1),rep(NA, (length(year)-1))), t.sigma.proc = runif(1, 0, 10), t.mean.r = rnorm(1), t.sigma.obs = runif (1, 0, 10), t.logN.est = c(rnorm(1, 20, 0.1), rep(NA, (length(year)-1))))}

# Parameters monitored
parameters <- c("r", "mean.r", "sigma2.obs", "sigma2.proc", "N.est", "beta", "t.beta", "t.r", "t.mean.r", "t.sigma2.obs", "t.sigma.proc", "t.N.est")

# MCMC settings
ni <- 20000; nt <- 6; nb <- 19000; nc <- 3

# Call JAGS from R (BRT 3 min)
carp.ssm2 <- jags(jags.data, inits, parameters, "ssm2.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(carp.ssm2, digits = 3)


posterior <- carp.ssm2$BUGSoutput$sims.matrix

posterior <- as.data.frame(posterior)

posterior2 <- posterior[,c("mean.r","deviance","sigma2.obs" ,"sigma2.proc","beta", "t.beta",
                           "t.mean.r", "t.sigma2.obs", "t.sigma.proc")]

posterior2 <- melt(posterior2)

ggplot (posterior2, aes(value))+
  geom_histogram(fill="#21908CFF", alpha=0.6)+
  geom_vline(xintercept=0, colour="#440154FF")+
  theme_minimal()+
  facet_wrap(~variable, scales="free")+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size=12),
        strip.text = element_text(size=12))+
  ggtitle("Posteriors")

# Draw figure
fitted <- lower <- upper <- numeric()
year <- carp.dat$year
n.years <- length(year)
for (i in 1:n.years){
  fitted[i] <- mean(carp.ssm2$BUGSoutput$sims.list$N.est[,i])
  lower[i] <- quantile(carp.ssm2$BUGSoutput$sims.list$N.est[,i], 0.025)
  upper[i] <- quantile(carp.ssm2$BUGSoutput$sims.list$N.est[,i], 0.975)}

plot.dat <- as.data.frame(cbind(year, fitted, lower, upper))

ggplot (plot.dat, aes(year, fitted))+
  geom_line(colour="#440154FF", size=1, aes())+
  geom_ribbon(aes(ymin=lower, ymax=upper), fill= "#21908CFF", alpha=0.5)+
  geom_point(data=data.frame(cbind(plot.dat$year, carp)),  aes(year, carp), colour="#FDE725FF")+
  ylab("Population size")+
  xlab("Year")+
  theme_minimal()+
  theme(axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text = element_text(size=12))+
  ggtitle("Carp population estimate")



