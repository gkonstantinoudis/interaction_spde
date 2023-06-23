#:::::::::::::::::::::::::::::
# Project: spde_interaction
# description: initialize R session
#:::::::::::::::::::::::::::::


# Load the libraries
library(INLA)
library(sf)
library(spdep)
library(fields)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(FNN)

path_data <- "data"
path_output <- "output"
path_rfiles <- "rfiles"


##
## INLA global options and priors
hyper.iid <- list(theta=list(prior="pc.prec", param=c(0.1, 0.01)))
hyper.bym <- list(theta1=list("PCprior", c(1, 0.01)), theta2 = list("PCprior", c(0.5, 0.5)))
control.family <- inla.set.control.family.default()


##
## INLA function
INLA_FUN <- function(form, dataset, Alphamat = NULL, theta.mode = NULL, verb = FALSE, rerun = FALSE){
  inla(form, 
       data = dataset, 
       family = "Poisson",
       control.family = control.family, 
       verbose = verb,
       num.threads = 1:1, 
       # set the priors for the fixed effect
       control.fixed = list(
         mean = list(`as.factor(hol)1` = 0, year_cont = 0, pm10_lag0_3 = 0),
         prec = list(`as.factor(hol)1` = 1, year_cont = 1, pm10_lag0_3 = 1),
         # have a proper prior for the intercept
         prec.intercept = 5, mean.intercept = 0),
       control.compute = list(config = TRUE), 
       control.predictor = list(link = 1, compute = TRUE, A = Alphamat), 
       control.mode = list(theta = theta.mode, restart = TRUE)
  ) -> mod
  
  if(rerun == TRUE){
    mod <- mod %>% inla.rerun()
  }
  return(mod)
}


##
## Stack function; the deaths should change here based on the simulation.
stk.function <- function(L, scenario){
  inla.stack(data = list(deaths_sim = simdata[[scenario]][,L]),
             A = list(A, 1),
             effects = list(list(field = 1:rf$n.spde),
                            data.frame(
                              expected = dat$expected,
                              hol = dat$hol,
                              year_cont = dat$year_cont,
                              datid = dat$datid, 
                              space_id = dat$space_id,
                              temp_id = dat$temp_id,
                              pm10_lag0_3 = dat$pm10_lag0_3)
             ), tag = "mod") %>% return()
}




