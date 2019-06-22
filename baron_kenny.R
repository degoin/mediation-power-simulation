
#-----------------------------------------------------------------------------------------------------------------------------
# baron & kenny estimator for natural direct and indirect effects 
#-----------------------------------------------------------------------------------------------------------------------------

baron_kenny <- function(obsdat, ymodel, mmodel) {
  
  # fit outcome model 
  fit1 <- glm(formula=ymodel, data=obsdat)
  # fit mediator model
  fit2 <- glm(formula=mmodel, data=obsdat)

# bootstrap for variance of NIE  

  bk_boot_nie <- function(iteration) {
    boot_df <- obsdat[sample(nrow(obsdat), replace=T),]
    
    b_fit1 <- glm(formula=ymodel, data=boot_df)
    b_fit2 <- glm(formula=mmodel, data=boot_df) 
    
    b_nie <- coef(summary(b_fit1))[3]*coef(summary(b_fit2))[2]
    return(b_nie)
  }
  
nie_b_list <- lapply(1:1000, function(x) bk_boot_nie(x))
nie_b <- do.call(rbind, nie_b_list)

# direct effect results
bk_nde <- data.frame(nde = coef(summary(fit1))[2], 
                      nde_lb = coef(summary(fit1))[2] - 1.96*coef(summary(fit1))[2,2], 
                      nde_ub = coef(summary(fit1))[2] + 1.96*coef(summary(fit1))[2,2])

# indirect effect results
bk_nie <- data.frame(nie = coef(summary(fit1))[3]*coef(summary(fit2))[2], 
                        nie_lb = quantile(nie_b,0.025), nie_ub = quantile(nie_b,0.975))


bk <- cbind(bk_nde, bk_nie) 
return(bk) 
}

