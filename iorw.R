
#-----------------------------------------------------------------------------------------------------------------------------
# inverse odds ratio weighting estimator for natural direct and indirect effects on the risk difference scale 
#-----------------------------------------------------------------------------------------------------------------------------

iorw <- function(obsdat, ammodel, ymmodel) {

  # calculate inverse odds weights
  am_fit <- glm(formula = ammodel, family="quasibinomial", data=obsdat)
  obsdat$p <- predict(am_fit, type="response")
  obsdat$iow <- (1 - obsdat$p)/obsdat$p 
  # assign unexposed weight of 1
  obsdat$iow[obsdat$a==0] <- 1

  # direct effect estimation 
  nde_fit <- glm(formula = ymmodel, family="gaussian", weights=iow, data=obsdat)
  te_fit <- glm(formula = ymmodel, family="gaussian", data=obsdat)

  # bootstrap for inference 
  iorw_boot <- function(iteration) {
    boot_df <-  obsdat[sample(nrow(obsdat), replace=T),]
  
    am_fit <- glm(formula=ammodel, family="quasibinomial", data=boot_df)
    p <- predict(am_fit, type="response")
  
    iow <- (1 - p)/p 
    iow[boot_df$a==0] <- 1
  
    nde_fit <- glm(formula=ymmodel, family="gaussian", weights=iow, data=boot_df)
    te_fit <- glm(formula=ymmodel, family="gaussian", data=boot_df)
    
    iorw_b <- c(NDE = coef(nde_fit)[2], NIE = coef(te_fit)[2] - coef(nde_fit)[2])
    return(iorw_b)
  }

  iorw_b_list <- lapply(1:250, function(x) iorw_boot(x))
  iorw_b_df <- data.frame(do.call(rbind, iorw_b_list))

# direct effect results
iorw_nde <- data.frame(nde = coef(nde_fit)[2], 
                       nde_lb = coef(nde_fit)[2] - 1.96*sqrt(var(iorw_b_df$NDE)), 
                       nde_ub = coef(nde_fit)[2] + 1.96*sqrt(var(iorw_b_df$NDE)))

# indirect effect results
iorw_nie <- data.frame(nie = coef(te_fit)[2] - coef(nde_fit)[2], 
                       nie_lb = (coef(te_fit)[2] - coef(nde_fit)[2]) - 1.96*sqrt(var(iorw_b_df$NIE)), 
                       nie_ub = (coef(te_fit)[2] - coef(nde_fit)[2]) + 1.96*sqrt(var(iorw_b_df$NIE)))


iorw <- cbind(iorw_nde, iorw_nie)
return(iorw)
}

