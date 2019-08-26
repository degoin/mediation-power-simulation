#-----------------------------------------------------------------------------------------------------------------------------
#     
# power function -- this returns the confidence intervals for each method, 
#                     determines whether or not an effect was detected, 
#                     and assesses whether the confidence interval covers the truth
#
#-----------------------------------------------------------------------------------------------------------------------------

mpower <- function(iteration,superN, n, fulldat, truth, amodel,mmodel, ymodel, qmodel, ammodel, ymmodel) {
  print(iteration)
  
  # sample of size n from the superpopulation   
  
  obsdat <- fulldat[sample(x = superN, size = n, replace=F), ]
  
  #-----------------------------------------------------------------------------------------------------------------------------
  # baron & kenny approach
  #-----------------------------------------------------------------------------------------------------------------------------
  
  bk <- baron_kenny(obsdat=obsdat, ymodel=ymodel, mmodel=mmodel)
  
  # one-sided test for power to detect a direct effect 
  nde_power <- ifelse(bk$nde_lb[1]>=0,1,0)
  
  # calculate whether confidence interval covers truth 
  nde_cov <- as.numeric(bk$nde_lb<=truth$nde & bk$nde_ub>=truth$nde)
  
  # one-sided test for power to detect an indirect effect 
  nie_power  <- ifelse(bk$nie_lb[1]>=0,1,0)
  
  # calculate whether confidence interval covers truth 
  nie_cov <- as.numeric(bk$nie_lb<=truth$nie & bk$nie_ub>=truth$nie)
  
  
  bk_nde <- data.frame(cbind(lb = bk$nde_lb, ub = bk$nde_ub, 
                             ed = nde_power, coverage = nde_cov,
                             parameter = "NDE", method="Baron & Kenny"))
  bk_nie <- data.frame(cbind(lb = bk$nie_lb, ub = bk$nie_ub, 
                             ed = nie_power, coverage = nie_cov, 
                             parameter = "NIE", method="Baron & Kenny"))
  
  bk <- rbind(bk_nde, bk_nie)
  
  #-----------------------------------------------------------------------------------------------------------------------------
  # baron & kenny approach with interaction
  #-----------------------------------------------------------------------------------------------------------------------------
  
  bk_int <- baron_kenny_int(obsdat=obsdat, ymodel=ymodel, mmodel=mmodel)
  
  # one-sided test for power to detect a direct effect 
  nde_power <- ifelse(bk_int$nde_lb[1]>=0,1,0)
  
  # calculate whether confidence interval covers truth 
  nde_cov <- as.numeric(bk_int$nde_lb<=truth$nde & bk_int$nde_ub>=truth$nde)
  
  # one-sided test for power to detect an indirect effect 
  nie_power  <- ifelse(bk_int$nie_lb[1]>=0,1,0)
  
  # calculate whether confidence interval covers truth 
  nie_cov <- as.numeric(bk_int$nie_lb<=truth$nie & bk_int$nie_ub>=truth$nie)
  
  
  bk_int_nde <- data.frame(cbind(lb = bk_int$nde_lb, ub = bk_int$nde_ub, 
                             ed = nde_power, coverage = nde_cov,
                             parameter = "NDE", method="Baron & Kenny Interaction"))
  bk_int_nie <- data.frame(cbind(lb = bk_int$nie_lb, ub = bk_int$nie_ub, 
                             ed = nie_power, coverage = nie_cov, 
                             parameter = "NIE", method="Baron & Kenny Interaction"))
  
  bk_int <- rbind(bk_int_nde, bk_int_nie)
  

  #-----------------------------------------------------------------------------------------------------------------------------
  # inverse odds ratio weighting approach
  #-----------------------------------------------------------------------------------------------------------------------------
 
  iorw <- iorw(obsdat=obsdat, ammodel=ammodel, ymmodel)
  
  # one-sided test for power to detect a direct effect 
  nde_power <- ifelse(iorw$nde_lb[1]>=0,1,0)
  
  # calculate whether confidence interval covers truth 
  nde_cov <- as.numeric(iorw$nde_lb<=truth$nde & iorw$nde_ub>=truth$nde)
  
  # one-sided test for power to detect an indirect effect 
  nie_power  <- ifelse(iorw$nie_lb[1]>=0,1,0)
  
  # calculate whether confidence interval covers truth 
  nie_cov <- as.numeric(iorw$nie_lb<=truth$nie & iorw$nie_ub>=truth$nie)
  
  
  iorw_nde <- data.frame(cbind(lb = iorw$nde_lb, ub = iorw$nde_ub, 
                             ed = nde_power, coverage = nde_cov,
                             parameter = "NDE", method="IORW"))
  iorw_nie <- data.frame(cbind(lb = iorw$nie_lb, ub = iorw$nie_ub, 
                             ed = nie_power, coverage = nie_cov, 
                             parameter = "NIE", method="IORW"))
  
  iorw <- rbind(iorw_nde, iorw_nie)
  

  #-----------------------------------------------------------------------------------------------------------------------------
  #  TMLE, variance calculated using EIC
  #-----------------------------------------------------------------------------------------------------------------------------

  sm_results <- medtmle_nointermedvar(obsdat=obsdat, amodel=amodel, mmodel=mmodel, ymodel=ymodel, qmodel=qmodel)
  
  # one-sided test for power to detect a direct effect 
  nde_power <- ifelse(sm_results$sde_lb[1]>=0,1,0)
  
  # calculate whether confidence interval covers truth 
  nde_cov <- as.numeric(sm_results$sde_lb<=truth$nde & sm_results$sde_ub>=truth$nde)
  
  # one-sided test for power to detect an indirect effect 
  nie_power  <- ifelse(sm_results$sie_lb[1]>=0,1,0)
  
  # calculate whether confidence interval covers truth 
  nie_cov <- as.numeric(sm_results$sie_lb<=truth$nie & sm_results$sie_ub>=truth$nie)
  
  
  tmle_nde <- data.frame(cbind(lb = sm_results$sde_lb, ub = sm_results$sde_ub, 
                               ed = nde_power, coverage = nde_cov,
                               parameter = "NDE", method="TMLE"))
  tmle_nie <- data.frame(cbind(lb = sm_results$sie_lb, ub = sm_results$sie_ub, 
                               ed = nie_power, coverage = nie_cov, 
                               parameter = "NIE", method="TMLE"))
  
  tmle <- rbind(tmle_nde, tmle_nie)
  

  #-----------------------------------------------------------------------------------------------------------------------------
  #  TMLE, variance calculated using bootstrap
  #-----------------------------------------------------------------------------------------------------------------------------

  # define bootstrap function 
  tmle_boot <- function(iteration) {
    boot_df <-  obsdat[sample(nrow(obsdat), replace=T),]
    b_results <- medtmle_nointermedvar(obsdat=boot_df, amodel=amodel, mmodel=mmodel, ymodel=ymodel, qmodel=qmodel)
    return(b_results)
  }
  
  tmle_ls <- lapply(1:250, function(x) tmle_boot(x))
  tmle_df <- do.call(rbind, tmle_ls)
  
  sm_boot_results <- data.frame(sde=sm_results$sde, 
                                sde_lb = sm_results$sde - 1.96*sqrt(var(tmle_df$sde)), 
                                sde_ub = sm_results$sde + 1.96*sqrt(var(tmle_df$sde)), 
                                sie=sm_results$sie, 
                                sie_lb = sm_results$sie - 1.96*sqrt(var(tmle_df$sie)), 
                                sie_ub = sm_results$sie + 1.96*sqrt(var(tmle_df$sie)))
  
  # one-sided test for power to detect a direct effect 
  nde_power <- ifelse(sm_boot_results$sde_lb[1]>=0,1,0)
  # calculate whether confidence interval covers truth 
  nde_cov <- as.numeric(sm_boot_results$sde_lb<=truth$nde & sm_boot_results$sde_ub>=truth$nde)
  
  # one-sided test for power to detect an indirect effect 
  nie_power  <- ifelse(sm_boot_results$sie_lb[1]>=0,1,0)
  # calculate whether confidence interval covers truth 
  nie_cov <- as.numeric(sm_boot_results$sie_lb<=truth$nie & sm_boot_results$sie_ub>=truth$nie)
  
  
  tmle_nde <- data.frame(cbind(lb = sm_boot_results$sde_lb, ub = sm_boot_results$sde_ub, 
                               ed = nde_power, coverage = nde_cov,
                               parameter = "NDE", method="TMLE bootstrap"))
  tmle_nie <- data.frame(cbind(lb = sm_boot_results$sie_lb, ub = sm_boot_results$sie_ub, 
                               ed = nie_power, coverage = nie_cov, 
                               parameter = "NIE", method="TMLE bootstrap"))
  
  tmle_bs <- rbind(tmle_nde, tmle_nie)
  
  
  #-----------------------------------------------------------------------------------------------------------------------------
  #  Equation 
  #-----------------------------------------------------------------------------------------------------------------------------
  eq <- data.frame(mpower_eq(obsdat, mmodel=mmodel, ymodel=ymodel, MonY=MonY))
  eq$method = "Equation"

  # combine and compare results from all 3 methods 
  results <- data.frame(rbind(bk, bk_int, iorw, tmle, tmle_bs, eq))
  
  return(results)
}
