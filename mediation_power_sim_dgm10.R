#-----------------------------------------------------------------------------------------------------------------------------
#     
#  simulation to compare performance of the Baron and Kenny, Inverse odds ratio weighting, and 
#     Targeted minimum loss-based estimation for direct and indirect effects when there is an intermediate variable 
#
#-----------------------------------------------------------------------------------------------------------------------------
rm(list=ls())
library(dplyr)

  source("./mpower_intermedvar.R")
  source("./medtmle_intermedvar.R")
  source("./baron_kenny.R")
  source("./iorw.R")
  source("./mpower_eq.R")

#-----------------------------------------------------------------------------------------------------------------------------
#     
# define the truth using a superpopulation with data-generating mechanism that fits the data
#
#-----------------------------------------------------------------------------------------------------------------------------

# define DGM 

# relationships between key variables
AonZ <- 0.6
AonM <- 0.1
AMonY <- 0
AonY <- 0.15

ZonM <- 0.2
ZonY <- 0.3

MonY <- 0.3

# set the super population sample size, sample size used to estimate performance, and number of simulations
superN <- 5000000
n <- 100
simN <- 1000

# define super population
w <- rbinom(superN, 1, 0.2)
a <- rbinom(superN, 1, prob=0.4 + 0.5*w)
z <- rbinom(superN,1, prob=0.1 + AonZ*a + 0.2*w)

Ma1    <- rbinom(superN, 1, prob=0.125 + ZonM*z + AonM + 0.2*w)
Ma0    <- rbinom(superN, 1, prob=0.125 + ZonM*z        + 0.2*w)

Ma1ns  <- rbinom(superN, 1, prob=0.125 + ZonM*z + AonM + 0.2*w)
Ma0ns  <- rbinom(superN, 1, prob=0.125 + ZonM*z        + 0.2*w)

Ya1m1 <- rbinom(superN,1, prob= 0.05 + MonY*Ma1ns + AMonY*Ma1ns + ZonY*z + AonY + 0.2*w)
Ya1m0 <- rbinom(superN,1, prob= 0.05 + MonY*Ma0ns + AMonY*Ma0ns + ZonY*z + AonY + 0.2*w)
Ya0m0 <- rbinom(superN,1, prob= 0.05 + MonY*Ma0ns + AMonY*Ma0ns + ZonY*z +        0.2*w)


#estimate the observed probability of M given A, Z, and W
mfit1 <- glm(Ma1 ~ a + z + w, family="binomial")
pMa1z1a1 <- predict(mfit1, newdata = data.frame(cbind(a=1, z=1, w=w)), type="response")
pMa1z0a1 <- predict(mfit1, newdata = data.frame(cbind(a=1, z=0, w=w)), type="response")

mfit0 <- glm(Ma0 ~ a + z + w, family="binomial")
pMa0z1a0 <- predict(mfit0, newdata = data.frame(cbind(a=0, z=1, w=w)), type="response")
pMa0z0a0 <- predict(mfit0, newdata = data.frame(cbind(a=0, z=0, w=w)), type="response")

#estimate the observed probability of Z given A and W
zfit <- glm(z ~ a + w, family="binomial")
pZa1 <- predict(zfit, newdata=data.frame(cbind(a=1, w)), type="response")
pZa0 <- predict(zfit, newdata=data.frame(cbind(a=0, w)), type="response")


# empirical distribution of m under a=0, marginalized over z
gma0 <- pMa0z1a0*pZa0 + pMa0z0a0*(1-pZa0)
# empirical distribution of m under a=1, marginalized over z
gma1 <- pMa1z1a1*pZa1 + pMa1z0a1*(1-pZa1)


dat <- data.frame(cbind(Ya1m1,Ya0m0, Ma0ns, Ma1ns, a, z, w, gma0, gma1))
# M, and Y should be defined based on observed data 
# so if A=1, then you have a mediator drawn from distribution under A=1 // else under A=0 

dat$m <- ifelse(dat$a==1,dat$Ma1ns,dat$Ma0ns)
dat$y <- ifelse(dat$a==1, dat$Ya1m1, dat$Ya0m0)

fulldat <- dplyr::select(dat, c(a,w,m,z,y, gma0, gma1))


# define models
amodel <- "a ~ w"
ammodel <- "a ~ m + z + w"
zmodel <- "z ~ a + w"
mmodel <- "m ~ a + z + w"
ymodel <- "y ~ a + m + z + w"
ymmodel <- "y ~ a + z + w"
qmodel <- "w"

ammodel_noz <- "a ~ m + w"
mmodel_noz <- "m ~ a + w"
ymodel_noz <- "y ~ a + m + w"
ymmodel_noz <- "y ~ a + w"

# calculate truth from superpopulation 
truth <- medtmle_intermedvar(obsdat=fulldat, amodel=amodel, zmodel=zmodel, mmodel=mmodel, ymodel=ymodel, qmodel=qmodel)

# run power simulations
  all_results_list_a <- lapply(1:simN, function(x) mpower(x,superN, n, fulldat=fulldat, truth=truth, 
                                                          amodel=amodel, zmodel=zmodel, mmodel=mmodel, 
                                                          ymodel=ymodel, qmodel=qmodel, ammodel=ammodel, 
                                                          ymmodel=ymmodel, ammodel_noz = ammodel_noz, 
                                                          mmodel_noz = mmodel_noz, ymodel_noz = ymodel_noz, 
                                                          ymmodel_noz = ymmodel_noz))

# combine results
all_results_t <- do.call(rbind, all_results_list_a)
# make sure these variables are numeric
all_results_t$ed <- as.numeric(levels(all_results_t$ed))[all_results_t$ed]
all_results_t$coverage <- as.numeric(levels(all_results_t$coverage))[all_results_t$coverage]

# summarize power and coverage 
power_results <- all_results_t %>% group_by(parameter, method) %>% summarise(power = mean(ed), coverage=mean(coverage))

# record true effect size 
power_results$effect_size <- ifelse(power_results$parameter=="NDE",truth$sde,truth$sie)


