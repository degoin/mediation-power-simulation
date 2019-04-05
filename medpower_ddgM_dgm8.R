rm(list=ls())
#-----------------------------------------------------------------------------------------------------------------------------
#     
# mediation simulation to compare B&K, IORW, and TMLE when there is an intermediate variable 
#
#-----------------------------------------------------------------------------------------------------------------------------

library(dplyr)


amodel <- "a ~ w"
ammodel <- "a ~ m + z + w"
zmodel <- "z ~ a + w"
mmodel <- "m ~ a + z + w"
ymodel <- "y ~ a*m + z + w"
ymmodel <- "y ~ a + z + w"
qmodel <- "w"

ammodel_noz <- "a ~ m + w"
mmodel_noz <- "m ~ a + w"
ymodel_noz <- "y ~ a*m + w"
ymmodel_noz <- "y ~ a + w"


AonZ <- 0.6
AonM <- 0.1
AMonY <- 0.2
AonY <- 0.2

ZonM <- 0.2
ZonY <- 0.2

MonY <- 0.15

superN <- 5000000
n <- 100
simN <- 1000

cluster <- F

if (cluster==T){
source("/home/degoin/Mediation_power/scripts/medtmle_intermedvar.R")
source("/home/degoin/Mediation_power/scripts/medtmle_intermedvar_pt.R")
}
if (cluster==F) {
source("/Users/degoin/Documents/Research projects/Mediation power/scripts/medtmle_intermedvar.R")
source("/Users/degoin/Documents/Research projects/Mediation power/scripts/medtmle_intermedvar_pt.R")
}

# Cluster setup: This piece spins up all of your slaves.
if (cluster==T) {
  library("Rmpi")
  
  
  # Spawn as many slaves as possible
  mpi.spawn.Rslaves()
  
  # In case R exits unexpectedly, have it automatically clean up
  # resources taken up by Rmpi (slaves, memory, etc...)
  .Last <- function(){
    if (is.loaded("mpi_initialize")){
      if (mpi.comm.size(1) > 0){
        print("Please use mpi.close.Rslaves() to close slaves.")
        mpi.close.Rslaves()
      }
      print("Please use mpi.quit() to quit R")
      .Call("mpi_finalize")
    }
  }
  
  # Tell all slaves to return a message identifying themselves
  mpi.remote.exec(paste("I am",mpi.comm.rank(),"of",mpi.comm.size()))
}

 
#-----------------------------------------------------------------------------------------------------------------------------
#     
# define the truth using a superpopulation
#
#-----------------------------------------------------------------------------------------------------------------------------
w <- rbinom(superN, 1, 0.2)
a <- rbinom(superN, 1, prob=.4 + .5*w)
z <- rbinom(superN,1, prob=AonZ*a + .2*w)

Ma1    <- rbinom(superN, 1, prob=ZonM*z + AonM + .2*w )
Ma0    <- rbinom(superN, 1, prob=ZonM*z        + .2*w )

Ma1ns  <- rbinom(superN, 1, prob=ZonM*z + AonM + .2*w )
Ma0ns  <- rbinom(superN, 1, prob=ZonM*z        + .2*w )

Ya1m1 <- rbinom(superN,1, prob=  MonY*Ma1ns + AMonY*Ma1ns + ZonY*z + AonY + .2*w)
Ya1m0 <- rbinom(superN,1, prob=  MonY*Ma0ns + AMonY*Ma0ns + ZonY*z + AonY + .2*w)
Ya0m0 <- rbinom(superN,1, prob=  MonY*Ma0ns + AMonY*Ma0ns + ZonY*z +        .2*w)


#estimate the observed probability of M given A, Z, and W
pMa1z1a1 <- predict(glm(Ma1 ~ a + z + w, family="binomial"), newdata = data.frame(cbind(a=1, z=1, w=w)), type="response")
pMa1z0a1 <- predict(glm(Ma1 ~ a + z + w, family="binomial"), newdata = data.frame(cbind(a=1, z=0, w=w)), type="response")

pMa0z1a0 <- predict(glm(Ma0 ~ a + z + w, family="binomial"), newdata = data.frame(cbind(a=0, z=1, w=w)), type="response")
pMa0z0a0 <- predict(glm(Ma0 ~ a + z + w, family="binomial"), newdata = data.frame(cbind(a=0, z=0, w=w)), type="response")

#estimate the observed probability of Z given A and W
pZa1 <- predict(glm(z ~ a + w, family="binomial"), newdata=data.frame(cbind(a=1, w)), type="response")
pZa0 <- predict(glm(z ~ a + w, family="binomial"), newdata=data.frame(cbind(a=0, w)), type="response")


# empirical distribution of m under a=0, marginalized over z
#gma0 <- (mz1a0*za0) + (mz0a0*(1-za0))
gma0 <- pMa0z1a0*pZa0 + pMa0z0a0*(1-pZa0)
# empirical distribution of m under a=1, marginalized over z
#gma1 <- (mz1a1*za1) + (mz0a1*(1-za1))
gma1 <- pMa1z1a1*pZa1 + pMa1z0a1*(1-pZa1)


dat <- data.frame(cbind(Ya1m1,Ya0m0, Ma0ns, Ma1ns, a, z, w, gma0, gma1))
# M, and Y should be defined based on observed data 
# so if A=1, then you have a mediator drawn from distribution under A=1 // else under A=0 

dat$m <- ifelse(dat$a==1,dat$Ma1ns,dat$Ma0ns)
dat$y <- ifelse(dat$a==1, dat$Ya1m1, dat$Ya0m0)

fulldat <- dplyr::select(dat, c(a,w,m,z,y, gma0, gma1))


truth <- medtmle_intermedvar_pt(obsdat=fulldat, amodel=amodel, zmodel=zmodel, mmodel=mmodel, ymodel=ymodel, qmodel=qmodel)
  


#-----------------------------------------------------------------------------------------------------------------------------
#     
# power function -- this returns the confidence intervals for each method and 
#                     determines whether or not an effect was detected
#
#-----------------------------------------------------------------------------------------------------------------------------

mpower <- function(iteration,superN, n, fulldat, truth, amodel, zmodel, mmodel, ymodel, qmodel, ammodel, ymmodel, 
                   ammodel_noz, mmodel_noz, ymodel_noz, ymmodel_noz) {
  print(iteration)
  
  # sample size n from the superpopulation   
  
  obsdat <- fulldat[sample(x = superN, size = n, replace=F), ]
  
  #-----------------------------------------------------------------------------------------------------------------------------
  # baron & kenny approach
  #-----------------------------------------------------------------------------------------------------------------------------
  
  print("BK omit Z")
  fit1 <- glm(formula=ymodel_noz, data=obsdat)
  fit2 <- glm(formula=mmodel_noz, data=obsdat)
  
  # resample to get NIE variance 
  
  bk_boot_nie <- function(iteration) {
    boot_df <- obsdat[sample(nrow(obsdat), replace=T),]
    
    b_fit1 <- glm(formula=ymodel_noz, data=boot_df)
    b_fit2 <- glm(formula=mmodel_noz, data=boot_df) 
    
    b_nie <- coef(summary(b_fit1))[3]*coef(summary(b_fit2))[2]
    return(b_nie)
  }
  
  nie_b_list <- lapply(1:1000, function(x) bk_boot_nie(x))
  nie_b <- do.call(rbind, nie_b_list)
  
  bk_nde_ci <- data.frame(nde_lb = coef(summary(fit1))[2] - 1.96*coef(summary(fit1))[2,2], 
                          nde_ub = coef(summary(fit1))[2] + 1.96*coef(summary(fit1))[2,2])
  bk_nde_power <- data.frame(ifelse(bk_nde_ci[1]>=0,1,0))
  names(bk_nde_power) <- "ed"
  bk_nde_ci$cov <- bk_nde_ci$nde_lb<=truth$sde & bk_nde_ci$nde_ub>=truth$sde
  
  bk_nie_ci <- data.frame(nie_lb = quantile(nie_b,0.025), nie_ub = quantile(nie_b,0.975))
  bk_nie_power  <- data.frame(ifelse(bk_nie_ci[1]>=0,1,0))
  names(bk_nie_power) <- "ed"
  bk_nie_ci$cov <- bk_nie_ci$nie_lb<=truth$sie & bk_nie_ci$nie_ub>=truth$sie
  
  
  bk_nde <- data.frame(cbind(lb = bk_nde_ci$nde_lb, ub = bk_nde_ci$nde_ub, 
                             ed = bk_nde_power, coverage = as.numeric(bk_nde_ci$cov),
                             parameter = "NDE", method="BK omit Z"))
  bk_nie <- data.frame(cbind(lb = bk_nie_ci$nie_lb, ub = bk_nie_ci$nie_ub, 
                             ed = bk_nie_power, coverage = as.numeric(bk_nie_ci$cov), 
                             parameter = "NIE", method="BK omit Z"))
  
  bk_oz <- rbind(bk_nde, bk_nie)
  
  
  
  # ---------------------------------------------------------------------------------------
  print("BK control for Z")
  fit1 <- glm(formula=ymodel, data=obsdat)
  fit2 <- glm(formula=mmodel, data=obsdat)
  
  # resample to get NIE variance 
  
  bk_boot_nie <- function(iteration) {
    boot_df <- obsdat[sample(nrow(obsdat), replace=T),]
    
    b_fit1 <- glm(formula=ymodel, data=boot_df)
    b_fit2 <- glm(formula=mmodel, data=boot_df) 
    
    b_nie <- coef(summary(b_fit1))[3]*coef(summary(b_fit2))[2]
    return(b_nie)
  }
  
  nie_b_list <- lapply(1:1000, function(x) bk_boot_nie(x))
  nie_b <- do.call(rbind, nie_b_list)
  
  bk_nde_ci <- data.frame(nde_lb = coef(summary(fit1))[2] - 1.96*coef(summary(fit1))[2,2], 
                          nde_ub = coef(summary(fit1))[2] + 1.96*coef(summary(fit1))[2,2])
  bk_nde_power <- data.frame(ifelse(bk_nde_ci[1]>=0,1,0))
  names(bk_nde_power) <- "ed"
  bk_nde_ci$cov <- bk_nde_ci$nde_lb<=truth$sde & bk_nde_ci$nde_ub>=truth$sde
  
  bk_nie_ci <- data.frame(nie_lb = quantile(nie_b,0.025), nie_ub = quantile(nie_b,0.975))
  bk_nie_power  <- data.frame(ifelse(bk_nie_ci[1]>=0,1,0))
  names(bk_nie_power) <- "ed"
  bk_nie_ci$cov <- bk_nie_ci$nie_lb<=truth$sie & bk_nie_ci$nie_ub>=truth$sie
  
  
  bk_nde <- data.frame(cbind(lb = bk_nde_ci$nde_lb, ub = bk_nde_ci$nde_ub, 
                             ed = bk_nde_power, coverage = as.numeric(bk_nde_ci$cov),
                             parameter = "NDE", method="BK control Z"))
  bk_nie <- data.frame(cbind(lb = bk_nie_ci$nie_lb, ub = bk_nie_ci$nie_ub, 
                             ed = bk_nie_power, coverage = as.numeric(bk_nie_ci$cov), 
                             parameter = "NIE", method="BK control Z"))
  
  bk_cz <- rbind(bk_nde, bk_nie)
  
  
  #-----------------------------------------------------------------------------------------------------------------------------
  # inverse odds ratio weighting approach
  #-----------------------------------------------------------------------------------------------------------------------------
  print("IORW omit Z")
  am_fit <- glm(formula = ammodel_noz, family="quasibinomial", data=obsdat)
  obsdat$p <- predict(am_fit, type="response")
  
  obsdat$iow <- (1 - obsdat$p)/obsdat$p 
  obsdat$iow[obsdat$a==0] <- 1
  
  nde_fit <- glm(formula = ymmodel_noz, family="quasibinomial", weights=iow, data=obsdat)
  te_fit <- glm(formula = ymmodel_noz, family="quasibinomial", data=obsdat)
  
  dfa1 <- dfa0 <- obsdat 
  dfa1$a <- 1
  dfa0$a <- 0 
  
  pnde1 <- predict(nde_fit, type="response", newdata=dfa1)
  pnde0 <- predict(nde_fit, type="response", newdata=dfa0)
  
  pte1 <- predict(te_fit, type="response", newdata=dfa1)
  pte0 <- predict(te_fit, type="response", newdata=dfa0)
  
  iorw_nde <- mean(pnde1-pnde0)
  iorw_nie = mean(pte1-pte0) - mean(pnde1-pnde0)
  
  iorw_boot <- function(iteration) {
    boot_df <-  obsdat[sample(nrow(obsdat), replace=T),]
    
    am_fit <- glm(formula=ammodel_noz, family="quasibinomial", data=boot_df)
    p <- predict(am_fit, type="response")
    
    iow <- (1 - p)/p 
    iow[boot_df$a==0] <- 1
    
    nde_fit <- glm(formula=ymmodel_noz, family="quasibinomial", weights=iow, data=boot_df)
    te_fit <- glm(formula=ymmodel_noz, family="quasibinomial", data=boot_df)
    
    boot_df1 <- boot_df0 <- boot_df
    boot_df1$a <- 1
    boot_df0$a <- 0
    pnde1 <- predict(nde_fit, type="response", newdata=boot_df1)
    pnde0 <- predict(nde_fit, type="response", newdata=boot_df0)
    
    pte1 <- predict(te_fit, type="response", newdata=boot_df1)
    pte0 <- predict(te_fit, type="response", newdata=boot_df0)
    
    iorw_b <- c(NDE = mean(pnde1-pnde0), NIE = mean(pte1-pte0) - mean(pnde1-pnde0))
    return(iorw_b)
  }
  
  iorw_b_list <- lapply(1:250, function(x) iorw_boot(x))
  iorw_b_df <- data.frame(do.call(rbind, iorw_b_list))
  
  
  iorw_nde_ci <- data.frame(nde_lb = iorw_nde - 1.96*sqrt(var(iorw_b_df$NDE)), nde_ub = iorw_nde + 1.96*sqrt(var(iorw_b_df$NDE)))
  iorw_nde_power <- data.frame(ifelse(iorw_nde_ci[1]>=0,1, ifelse(iorw_nde_ci[2]<=0,1,0)))
  names(iorw_nde_power) <- "ed"
  iorw_nde_ci$cov <- iorw_nde_ci$nde_lb<= truth$sde & iorw_nde_ci$nde_ub>=truth$sde
  
  iorw_nie_ci <- data.frame(nie_lb = iorw_nie - 1.96*sqrt(var(iorw_b_df$NIE)), nie_ub = iorw_nde + 1.96*sqrt(var(iorw_b_df$NIE)))
  iorw_nie_power <- data.frame(ifelse(iorw_nie_ci[1]>=0,1, ifelse(iorw_nie_ci[2]<=0,1,0)))
  names(iorw_nie_power) <- "ed"
  iorw_nie_ci$cov <- iorw_nie_ci$nie_lb<= truth$sie & iorw_nie_ci$nie_ub>=truth$sie
  
  iorw_nde <- data.frame(lb = iorw_nde_ci$nde_lb, 
                         ub = iorw_nde_ci$nde_ub, 
                         ed = iorw_nde_power, 
                         coverage = as.numeric(iorw_nde_ci$cov), 
                         parameter = "NDE", method="IORW omit Z")
  iorw_nie <- data.frame(lb = iorw_nie_ci$nie_lb, 
                         ub = iorw_nie_ci$nie_ub, 
                         ed = iorw_nie_power, 
                         coverage = as.numeric(iorw_nie_ci$cov), 
                         parameter = "NIE", method="IORW omit Z")
  
  iorw_oz <- rbind(iorw_nde, iorw_nie)
  
  # -----------------------------------------------------------------------
  print("IORW control for Z")
  am_fit <- glm(formula = ammodel, family="quasibinomial", data=obsdat)
  obsdat$p <- predict(am_fit, type="response")
  
  obsdat$iow <- (1 - obsdat$p)/obsdat$p 
  obsdat$iow[obsdat$a==0] <- 1
  
  nde_fit <- glm(formula = ymmodel, family="quasibinomial", weights=iow, data=obsdat)
  te_fit <- glm(formula = ymmodel, family="quasibinomial", data=obsdat)
  
  dfa1 <- dfa0 <- obsdat 
  dfa1$a <- 1
  dfa0$a <- 0 
  
  pnde1 <- predict(nde_fit, type="response", newdata=dfa1)
  pnde0 <- predict(nde_fit, type="response", newdata=dfa0)
  
  pte1 <- predict(te_fit, type="response", newdata=dfa1)
  pte0 <- predict(te_fit, type="response", newdata=dfa0)
  
  iorw_nde <- mean(pnde1-pnde0)
  iorw_nie = mean(pte1-pte0) - mean(pnde1-pnde0)
  
  iorw_boot <- function(iteration) {
    boot_df <-  obsdat[sample(nrow(obsdat), replace=T),]
    
    am_fit <- glm(formula=ammodel, family="quasibinomial", data=boot_df)
    p <- predict(am_fit, type="response")
    
    iow <- (1 - p)/p 
    iow[boot_df$a==0] <- 1
    
    nde_fit <- glm(formula=ymmodel, family="quasibinomial", weights=iow, data=boot_df)
    te_fit <- glm(formula=ymmodel, family="quasibinomial", data=boot_df)
    
    boot_df1 <- boot_df0 <- boot_df
    boot_df1$a <- 1
    boot_df0$a <- 0
    pnde1 <- predict(nde_fit, type="response", newdata=boot_df1)
    pnde0 <- predict(nde_fit, type="response", newdata=boot_df0)
    
    pte1 <- predict(te_fit, type="response", newdata=boot_df1)
    pte0 <- predict(te_fit, type="response", newdata=boot_df0)
    
    iorw_b <- c(NDE = mean(pnde1-pnde0), NIE = mean(pte1-pte0) - mean(pnde1-pnde0))
    return(iorw_b)
  }
  
  iorw_b_list <- lapply(1:250, function(x) iorw_boot(x))
  iorw_b_df <- data.frame(do.call(rbind, iorw_b_list))
  
  
  iorw_nde_ci <- data.frame(nde_lb = iorw_nde - 1.96*sqrt(var(iorw_b_df$NDE)), nde_ub = iorw_nde + 1.96*sqrt(var(iorw_b_df$NDE)))
  iorw_nde_power <- data.frame(ifelse(iorw_nde_ci[1]>=0,1, ifelse(iorw_nde_ci[2]<=0,1,0)))
  names(iorw_nde_power) <- "ed"
  iorw_nde_ci$cov <- iorw_nde_ci$nde_lb<= truth$sde & iorw_nde_ci$nde_ub>=truth$sde
  
  iorw_nie_ci <- data.frame(nie_lb = iorw_nie - 1.96*sqrt(var(iorw_b_df$NIE)), nie_ub = iorw_nde + 1.96*sqrt(var(iorw_b_df$NIE)))
  iorw_nie_power <- data.frame(ifelse(iorw_nie_ci[1]>=0,1, ifelse(iorw_nie_ci[2]<=0,1,0)))
  names(iorw_nie_power) <- "ed"
  iorw_nie_ci$cov <- iorw_nie_ci$nie_lb<= truth$sie & iorw_nie_ci$nie_ub>=truth$sie
  
  iorw_nde <- data.frame(lb = iorw_nde_ci$nde_lb, 
                         ub = iorw_nde_ci$nde_ub, 
                         ed = iorw_nde_power, 
                         coverage = as.numeric(iorw_nde_ci$cov), 
                         parameter = "NDE", method="IORW control Z")
  iorw_nie <- data.frame(lb = iorw_nie_ci$nie_lb, 
                         ub = iorw_nie_ci$nie_ub, 
                         ed = iorw_nie_power, 
                         coverage = as.numeric(iorw_nie_ci$cov), 
                         parameter = "NIE", method="IORW control Z")
  
  iorw_cz <- rbind(iorw_nde, iorw_nie)
  
  #-----------------------------------------------------------------------------------------------------------------------------
  # stochastic TMLE approach
  #-----------------------------------------------------------------------------------------------------------------------------
  print("TMLE")
  sm_results <- medtmle_intermedvar(obsdat=obsdat, amodel=amodel, zmodel=zmodel, mmodel=mmodel, ymodel=ymodel, qmodel=qmodel)
  
  tmle_nde_ci <- data.frame(nde_lb = sm_results$sde_lb, 
                            nde_ub = sm_results$sde_ub, 
                            cov = sm_results$sde_lb<= truth$sde & sm_results$sde_ub>=truth$sde)
  tmle_nde_power <- data.frame(ifelse(tmle_nde_ci[1]>=0,1,0))
  names(tmle_nde_power) <- "ed"
  
  tmle_nie_ci <- data.frame(nie_lb = sm_results$sie_lb, 
                            nie_ub = sm_results$sie_ub, 
                            cov = sm_results$sie_lb<= truth$sie & sm_results$sie_ub>=truth$sie)
  tmle_nie_power <- data.frame(ifelse(tmle_nie_ci[1]>=0,1,0))
  names(tmle_nie_power) <- "ed"
  
  tmle_nde <- data.frame(lb = tmle_nde_ci$nde_lb, 
                         ub = tmle_nde_ci$nde_ub, 
                         ed = tmle_nde_power, 
                         coverage = as.numeric(tmle_nde_ci$cov),
                         parameter = "NDE", method="TMLE")
  
  tmle_nie <- data.frame(lb = tmle_nie_ci$nie_lb, 
                         ub = tmle_nie_ci$nie_ub, 
                         ed = tmle_nie_power, 
                         coverage = as.numeric(tmle_nie_ci$cov),
                         parameter = "NIE", method="TMLE")
  
  tmle <- rbind(tmle_nde, tmle_nie)
  
  
  # combine and compare results from all 3 methods 
  results <- data.frame(rbind(bk_oz, bk_cz, iorw_oz, iorw_cz, tmle))
  
  return(results)
}





# If using the cluster, send all the necessary info to the slaves
if (cluster==T) {  
  library(Rmpi)
  mpi.bcast.Robj2slave(superN)
  mpi.bcast.Robj2slave(n)
  mpi.bcast.Robj2slave(fulldat)
  mpi.bcast.Robj2slave(truth)
  mpi.bcast.Robj2slave(simN)
  mpi.bcast.Robj2slave(AonZ)
  mpi.bcast.Robj2slave(AonM)
  mpi.bcast.Robj2slave(AonY)
  mpi.bcast.Robj2slave(ZonM)
  mpi.bcast.Robj2slave(ZonY)
  mpi.bcast.Robj2slave(MonY)
  mpi.bcast.Robj2slave(AMonY)
  mpi.bcast.Robj2slave(amodel)
  mpi.bcast.Robj2slave(zmodel)
  mpi.bcast.Robj2slave(mmodel)
  mpi.bcast.Robj2slave(ymodel)
  mpi.bcast.Robj2slave(qmodel)
  mpi.bcast.Robj2slave(ammodel)
  mpi.bcast.Robj2slave(ymmodel)
  mpi.bcast.Robj2slave(mmodel_noz)
  mpi.bcast.Robj2slave(ymodel_noz)
  mpi.bcast.Robj2slave(ammodel_noz)
  mpi.bcast.Robj2slave(ymmodel_noz)
  mpi.bcast.Robj2slave(mpower)
  mpi.bcast.Robj2slave(medtmle_intermedvar)
  
  all_results_list_a <- mpi.parLapply(1:simN, function(x) mpower(x,superN, n, fulldat=fulldat, truth=truth, 
                                                                 amodel=amodel, zmodel=zmodel, mmodel=mmodel, 
                                                                 ymodel=ymodel, qmodel=qmodel, ammodel=ammodel, 
                                                                 ymmodel=ymmodel, ammodel_noz = ammodel_noz, 
                                                                 mmodel_noz = mmodel_noz, ymodel_noz = ymodel_noz, 
                                                                 ymmodel_noz = ymmodel_noz))
} 

if (cluster==F) {
  # run power simulation 
  all_results_list_a <- lapply(1:simN, function(x) mpower(x,superN, n, fulldat=fulldat, truth=truth, 
                                                          amodel=amodel, zmodel=zmodel, mmodel=mmodel, 
                                                          ymodel=ymodel, qmodel=qmodel, ammodel=ammodel, 
                                                          ymmodel=ymmodel, ammodel_noz = ammodel_noz, 
                                                          mmodel_noz = mmodel_noz, ymodel_noz = ymodel_noz, 
                                                          ymmodel_noz = ymmodel_noz))
  }

print(all_results_list_a)
# combine results
all_results_t <- do.call(rbind, all_results_list_a)


all_results_a <- all_results_t %>% group_by(parameter, method) %>% summarise(power = mean(ed), coverage=mean(coverage))

# record true effect size 
all_results_a$effect_size <- ifelse(all_results_a$parameter=="NDE",truth$sde,truth$sie)

if (cluster==T){
write.csv(all_results_a, file=paste0("/home/degoin/Mediation_power/results/mpower_ddgM_dgm8_",n,".csv"))
}
if (cluster==F) {
write.csv(all_results_a, file=paste0("/Users/degoin/Documents/Research projects/Mediation power/results/mpower_ddgM_dgm8_",n,".csv"))
}

if (cluster==T) { 
  mpi.close.Rslaves()
  mpi.quit() 
}














