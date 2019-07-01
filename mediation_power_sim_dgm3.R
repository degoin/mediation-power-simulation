#-----------------------------------------------------------------------------------------------------------------------------
#     
#  simulation to compare performance of the Baron and Kenny, Inverse odds ratio weighting, and 
#     Targeted minimum loss-based estimation for direct and indirect effects when there is an intermediate variable 
#
#-----------------------------------------------------------------------------------------------------------------------------
rm(list=ls())
setwd("/Users/degoin/Documents/Research projects/Mediation power/scripts/for github")
library(dplyr)

  source("./mpower_nointermedvar.R")
  source("./medtmle_nointermedvar.R")
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
AonM <- 0.1
AMonY <- 0
AonY <- 0.3
MonY <- 0.2

# set the super population sample size, sample size used to estimate performance, and number of simulations
superN <- 5000000
n <- 100
simN <- 1000

# calculate natural direct and indirect effects using Pearl's mediational formula 

ey_a1m0w0 <- 0.15 + AonY                  
ey_a1m0w1 <- 0.15 + AonY + 0.2
ey_a0m0w0 <- 0.15
ey_a0m0w1 <- 0.35
ey_a1m1w0 <- 0.15 + AonY + MonY + AMonY
ey_a1m1w1 <- 0.15 + AonY + MonY + AMonY + 0.2
ey_a0m1w0 <- 0.15 + MonY
ey_a0m1w1 <- 0.15 + MonY + 0.2

pm0_a0w0 <- 0.875
pm0_a0w1 <- 0.875 - 0.2
pm0_a1w0 <- 0.875 - AonM
pm0_a1w1 <- 0.875 - (AonM + 0.2)
pm1_a0w0 <- 0.125
pm1_a0w1 <- 0.325
pm1_a1w0 <- 0.125 + AonM
pm1_a1w1 <- 0.125 + AonM + 0.2


#NDE truth 
nde <- (ey_a1m0w0 - ey_a0m0w0)*pm0_a0w0*0.8 + (ey_a1m0w1 - ey_a0m0w1)*pm0_a0w1*0.2 + 
  (ey_a1m1w0 - ey_a0m1w0)*pm1_a0w0*0.8 + (ey_a1m1w1 - ey_a0m1w1)*pm1_a0w1*0.2

# NIE truth 
nie <- ey_a1m0w0*(pm0_a1w0 - pm0_a0w0)*0.8  + ey_a1m0w1*(pm0_a1w1 - pm0_a0w1)*0.2 + 
  ey_a1m1w0*(pm1_a1w0 - pm1_a0w0)*0.8 + ey_a1m1w1*(pm1_a1w1 - pm1_a0w1)*0.2

truth = data.frame(nde=nde, nie=nie)

# define super population
w <- rbinom(superN, 1, 0.2)
a <- rbinom(superN, 1, prob=0.4 + 0.5*w)

Ma1    <- rbinom(superN, 1, prob=0.125 + AonM + 0.2*w)
Ma0    <- rbinom(superN, 1, prob=0.125 +       0.2*w)

Ma1ns  <- rbinom(superN, 1, prob=0.125 + AonM + 0.2*w)
Ma0ns  <- rbinom(superN, 1, prob=0.125 +        0.2*w)

Ya1m1 <- rbinom(superN,1, prob= 0.15 + MonY*Ma1ns + AMonY*Ma1ns +  AonY + 0.2*w)
Ya1m0 <- rbinom(superN,1, prob= 0.15 + MonY*Ma0ns + AMonY*Ma0ns +  AonY + 0.2*w)
Ya0m0 <- rbinom(superN,1, prob= 0.15 + MonY*Ma0ns + AMonY*Ma0ns +         0.2*w)


#estimate the observed probability of M given A and W
gma1 <- predict(glm(Ma1 ~ a + w, family="binomial"), newdata = data.frame(cbind(a=1, w=w)), type="response")
gma0 <- predict(glm(Ma0 ~ a + w, family="binomial"), newdata = data.frame(cbind(a=0, w=w)), type="response")


dat <- data.frame(cbind(Ya1m1,Ya0m0, Ma0ns, Ma1ns, a, w, gma0, gma1))


# M, and Y should be defined based on observed data 
# so if A=1, then you have a mediator drawn from distribution under A=1 // else under A=0 

dat$m <- ifelse(dat$a==1,dat$Ma1ns,dat$Ma0ns)
dat$y <- ifelse(dat$a==1, dat$Ya1m1, dat$Ya0m0)

fulldat <- dplyr::select(dat, c(a,w,m,y, gma0, gma1))


# define models
amodel <- "a ~ w"
ammodel <- "a ~ m + w"
mmodel <- "m ~ a + w"
ymodel <- "y ~ a + m + w"
ymmodel <- "y ~ a + w"
qmodel <- "w"


# run power simulations
  all_results_list_a <- lapply(1:simN, function(x) mpower(x,superN, n, fulldat=fulldat, truth=truth, 
                                                          amodel=amodel, mmodel=mmodel, 
                                                          ymodel=ymodel, qmodel=qmodel, ammodel=ammodel, 
                                                          ymmodel=ymmodel))

# combine results
all_results_t <- do.call(rbind, all_results_list_a)
# make sure these variables are numeric
all_results_t$ed <- as.numeric(levels(all_results_t$ed))[all_results_t$ed]
all_results_t$coverage <- as.numeric(levels(all_results_t$coverage))[all_results_t$coverage]

# summarize power and coverage 
power_results <- all_results_t %>% group_by(parameter, method) %>% summarise(power = mean(ed), coverage=mean(coverage))

# record true effect size 
power_results$effect_size <- ifelse(power_results$parameter=="NDE",truth$nde,truth$nie)


