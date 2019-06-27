#-----------------------------------------------------------------------------------------------------------------------------
#     
#  power function using equation from Vittinghof et al 
#
#-----------------------------------------------------------------------------------------------------------------------------

mpower_eq <- function(dat, mmodel, ymodel, MonY) {
  
  # define vars for power calculation
  # calculate multiple correlation 
  corfit <- glm(formula=mmodel, data=dat)
  mhat <- predict(corfit)
  rho <- cov(mhat, dat$m)/(sqrt(var(mhat))*sqrt(var(dat$m)))
  
  # set type 1 error = 0.975
    z_a <- qnorm(0.975) 
  # set type 2 error = 0.2, which implies power=0.8
    z_g <- qnorm(0.8)

  lmfit <- glm(formula=ymodel, data=dat)
  b2 <- MonY 

  sigma_1 <- var(dat$a)
  sigma_2 <- var(dat$m)

  sigma_e <- var(residuals(lmfit))

  req_ss <- ((z_a + z_g)^2*sigma_e)/((b2^2*sigma_2)*(1-rho^2))

  power <- pnorm(sqrt((n*b2^2*sigma_2*(1-rho^2))/sigma_e) - z_a)

results <-  cbind(lb = NA, ub=NA, ed = power, coverage = 1.05, parameter="NIE")
return(results)
}

