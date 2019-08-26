#-----------------------------------------------------------------------------------------------------------------------------
# baron & kenny estimator of natural direct and indirect effects allowing for interaction
#-----------------------------------------------------------------------------------------------------------------------------
                       
baron_kenny_int <- function(obsdat, ymodel, mmodel) {
  
  # fit outcome model 
  fit1 <- glm(formula=ymodel, data=obsdat)
  # fit mediator model
  fit2 <- glm(formula=mmodel, data=obsdat)
  
  # bootstrap for variance of NIE  
  
  bk_boot_ne <- function(iteration) {
    boot_df <- obsdat[sample(nrow(obsdat), replace=T),]
    
    b_fit1 <- glm(formula=ymodel, data=boot_df)
    b_fit2 <- glm(formula=mmodel, data=boot_df) 
    
    b_nde <- mean(coef(summary(b_fit1))[2] + (coef(summary(b_fit1))[4]*(coef(summary(b_fit2))[1] + coef(summary(b_fit2))[3]*obsdat$w)))

    b_nie <- (coef(summary(b_fit1))[3]*coef(summary(b_fit2))[2]) + (coef(summary(b_fit1))[4]*coef(summary(b_fit2))[2])
    return(c(b_nde, b_nie))
  }
  
  ne_b_list <- lapply(1:1000, function(x) bk_boot_ne(x))
  ne_b <- do.call(rbind, ne_b_list)
  
  # direct effect results
  bk_nde <- data.frame(nde = mean(coef(summary(fit1))[2] + (coef(summary(fit1))[4]*(coef(summary(fit2))[1] + coef(summary(fit2))[3]*obsdat$w))), 
                       nde_lb = quantile(ne_b[,1],0.025), nde_ub = quantile(ne_b[,1],0.975))
  
  # indirect effect results
  bk_nie <- data.frame(nie = (coef(summary(fit1))[3]*coef(summary(fit2))[2]) + (coef(summary(fit1))[4]*coef(summary(fit2))[2]), 
                       nie_lb = quantile(ne_b[,2],0.025), nie_ub = quantile(ne_b[,2],0.975))
  
  
  bk <- cbind(bk_nde, bk_nie) 
  return(bk) 
}
      
