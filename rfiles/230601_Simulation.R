

# Created 25.05.2023

# Simulate data for SPDEs for synergy


# ------------------------------------------------------------------------------

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

setwd("C:/Users/gkonstan/OneDrive - Imperial College London/ICRF Imperial/Projects/Interactions simulation/")

# Fit the baseline model. We use this baseline model to get estimates for space, time, 
# airpollution, temperature etc.


checklinear <- FALSE
checktimes <- FALSE
runbaseline <- FALSE

path_data <- "data"
path_output <- "output"

##
## Prepare the data
dat <- readRDS(file.path(path_data, "ItalyData"))
shp <- readRDS(file.path(path_data, "shp"))
W.nb <- spdep::poly2nb(shp)
spdep::nb2INLA("W.adj", W.nb)
theta <- NULL

dat$temp_id <- inla.group(dat$temperature_lag0_3, method = "cut", n = 50)
dat$pm_id <- inla.group(dat$pm10_lag0_3, method = "cut", n = 50)
dat$datid <- as.numeric(as.factor(dat$date))
dat$year_cont <-  dat$year - min(dat$year) + 1 
dat$space_id <- as.numeric(as.factor(dat$SIGLA))

hyper.iid <- list(theta=list(prior="pc.prec", param=c(0.1, 0.01)))
hyper.bym <- list(theta1=list("PCprior", c(1, 0.01)), theta2 = list("PCprior", c(0.5, 0.5)))
control.family <- inla.set.control.family.default()


# Define the functions
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
stk.function <- function(L){
  inla.stack(data = list(deaths_sim = simdata[,L]),
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

if(checklinear == TRUE){
  ##
  ## Fit baseline model, check if linear
  form_baseline = 
    deaths ~ 
    offset(log(expected)) + 
    as.factor(hol) + 
    year_cont + 
    f(datid, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
    f(space_id, model = "bym2", graph = "W.adj", scale.model = TRUE, hyper = hyper.bym, constr = TRUE) +
    f(temp_id, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
    f(pm_id, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) 
  
  INLA_FUN(form = form_baseline, dataset = dat, verb = TRUE)-> mod_baseline
  plot(mod_baseline$summary.random$pm_id$`0.5quant`) # the effect seems linear so I would use just a line
}


# I will store the baseline inla model, as there might be fluctuations in some decimal digits when rerunning, and this
# way i am forcing results to correspond to identical and known figures. 
if(runbaseline == TRUE){
  ##
  ## Fit baseline model based on linearity on PM10
  form_baseline = 
    deaths ~ 
    offset(log(expected)) + 
    as.factor(hol) + 
    year_cont + 
    f(datid, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
    f(space_id, model = "bym2", graph = "W.adj", scale.model = TRUE, hyper = hyper.bym, constr = TRUE) +
    f(temp_id, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
    pm10_lag0_3
  
  INLA_FUN(form = form_baseline, dataset = dat)-> mod_baseline
  saveRDS(mod_baseline, file = file.path(path_output, "baselineINLAmodel"))
}else{
  mod_baseline <- readRDS(file.path(path_output, "baselineINLAmodel"))
}


##
##Code for fitting the main models
##

## Linear interaction model
dat$pol_temp <- scale(dat$pm10_lag0_3*dat$temperature_lag0_3)
sd.scale.store <- sd(dat$pm10_lag0_3*dat$temperature_lag0_3)

form_regint = 
  deaths_sim ~ 
  offset(log(expected)) + 
  as.factor(hol) + 
  year_cont + 
  f(datid, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
  f(space_id, model = "bym2", graph = "W.adj", scale.model = TRUE, hyper = hyper.bym, constr = TRUE) +
  f(temp_id, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
  pm10_lag0_3 +
  pol_temp

##
## SPDE

# Prepare data for interaction
cut_var <- 30
dat$temp_id_int <- inla.group(dat$temperature_lag0_3, method = "cut", n = cut_var)
dat$pm_id_int <- inla.group(dat$pm10_lag0_3, method = "cut", n = cut_var)

dat$interaction_id <- paste0(dat$temp_id_int, dat$pm_id_int)
dat$interaction_id <- as.numeric(as.factor(dat$interaction_id))

dat$pm_id_factor <- dat$pm_id_int %>% as.factor() %>% as.numeric()
dat$temp_id_factor <- dat$temp_id_int %>% as.factor() %>% as.numeric()

# Create a mesh
expand.grid(
  x = seq(from = min(dat$temperature_lag0_3), 
          to =max(dat$temperature_lag0_3), 
          length.out = 50),
  y = seq(from = min(dat$pm10_lag0_3), 
          to =max(dat$pm10_lag0_3), 
          length.out = 50)
) -> mock_shp

## ok lets generate the mesh. 
locs <- unique(dat[, c("temp_id_factor", "pm_id_factor")])
colnames(locs) <- c("x", "y")

# mesh
set.seed(11)
mesh <- inla.mesh.2d(loc.domain = locs, max.edge = c(1.5, 10),
                     offset = c(1, 4)) 

plot(mesh)
points(mock_shp$x, 
       mock_shp$y, col = "blue", cex = 0.4, pch = 19)
points(locs, col = "red", pch = 18, cex = 0.8)

# == projector matrix ==
# projector matrix
A <- inla.spde.make.A(mesh = mesh,
                      loc = dat[, c("temp_id_factor", "pm_id_factor")] %>% as.matrix())

# == model definition ==
# Mattern prior
rf <- inla.spde2.pcmatern(mesh = mesh,
                          prior.range = c(0.1,0.01),
                          prior.sigma = c(1,0.01), 
                          constr = TRUE)


form_spde_int = 
  deaths_sim ~ 
  offset(log(expected)) + 
  as.factor(hol) + 
  year_cont + 
  f(datid, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
  f(space_id, model = "bym2", graph = "W.adj", scale.model = TRUE, hyper = hyper.bym, constr = TRUE) +
  f(temp_id, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
  pm10_lag0_3 + 
  f(field, model = rf)


##
##Code for simulating data
##

## The idea now is to use the baseline model and add the different types of interactions
dat <- left_join(dat, data.frame(id_tmp_ing = mod_baseline$summary.random$temp_id$ID, 
                                                  ranef_temp = mod_baseline$summary.random$temp_id$`0.5quant`), 
                                  by = c("temp_id" = "id_tmp_ing"))

# the linear predictor
lp_baseline <- log(dat$expected) + 
  (mod_baseline$summary.fixed$`0.5quant` %*% (dat %>% 
                                                  ungroup() %>% 
                                                  mutate(int = 1) %>% 
                                                  select(int, hol, year_cont, pm10_lag0_3) %>% 
                                                  t()) %>% t()) + 
  mod_baseline$summary.random$datid$`0.5quant`[dat$datid] + 
  mod_baseline$summary.random$space_id$`0.5quant`[dat$space_id] + 
  dat$ranef_temp

# I could have also taken the linear predictor from the model, but it will have the uncertainty of the estimates
# incorporated. I dont want that, I want the truth to be fixed. 

lp_baseline <- lp_baseline %>% as.numeric()


# I will simulate based on 5 scenarios
# a. no interaction
# b. linear interaction
# c. 3 types bi-dimentional interactions based on this paper: https://doi.org/10.1016/j.envint.2023.107825


##
##
# DEFINE THE TRUTH FOR INTERACTIONS

# For the truth I will define a 50 times 50 matrix based on the values of pm and temperature.
N_truth <- 50

expand.grid(
  pm10 = seq(from = min(dat$pm10_lag0_3), 
             to = max(dat$pm10_lag0_3), 
             length.out = N_truth), 
  temperature = seq(from = min(dat$temperature_lag0_3), 
                    to = max(dat$temperature_lag0_3), 
                    length.out = N_truth)
) -> gridint
gridint$ID <- paste0(gridint$pm10, gridint$temperature) %>% as.factor() %>% as.numeric()

FNN::get.knnx(gridint %>% select(pm10, temperature), 
              dat %>% ungroup() %>% select(pm10_lag0_3, temperature_lag0_3), 
              k = 1)$nn.index -> knn_calc

dat$ID_truth <- gridint$ID[as.numeric(knn_calc)]

# Scenario 1
gridint$rr_truth_sc1 <- 0
dat <- left_join(dat, gridint %>% select(ID, rr_truth_sc1), by = c("ID_truth" = "ID"))

# Scenario 2

# Scenario 3

# Scenario 4

# Scenario 5


int_list <- 
  list(
    noint = dat$rr_truth_sc1, 
    linearint = dat$rr_truth_sc2, 
    bidiint_1 = dat$rr_truth_sc3, 
    bidiint_2 = dat$rr_truth_sc4, 
    bidiint_3 = dat$rr_truth_sc5
  )


K <- 1 # here is a loop of the scenarios
lptosim <- lp_baseline + int_list[[K]]
set.seed(11)
lapply(1:300, function(i) rpois(n = lptosim %>% length(), lambda = lptosim %>% exp()) %>% return()) %>% 
  do.call(cbind,.) -> simdata




# check times
if(checktimes == TRUE){
  
  ##
  ## linear model
  t_0 <- Sys.time()
  INLA_FUN(form = form_regint, 
           dataset = dat, 
           verb = FALSE,
           theta.mode = c(-1.19, 5.45, -2.37, 5.92), 
           rerun = TRUE) -> mod_regint
  t_1 <- Sys.time()
  t_regint <- t_1 - t_0
  
  # TIME WITHOUT RERUN: ~1.05 minutes
  # TIME WITH RERUN: ~1.96 minutes
  # ~1.96*300*5 = 2940 minutes or 49 hours or 2.04 days
  
  print("Linear interaction:")
  print(t_regint)
  
  ##
  ## SPDE model
  stk <- stk.function(L = 1)
  
  t_0 <- Sys.time()
  INLA_FUN(form = form_spde_int, 
           dataset = inla.stack.data(stk), 
           Alphamat = inla.stack.A(stk), 
           verb = FALSE,
           theta.mode = c(-1.19, 5.45, -2.37, 5.92, 4.35, -2.66), 
           rerun = TRUE) -> mod_spdeint
  t_1 <- Sys.time()
  t_spdeint <- t_1 - t_0
  
  # TIME WITHOUT RERUN: ~3.58 minutes
  # TIME WITH RERUN: ~6.86 minutes
  # ~6.86*300*5 = 10290 minutes or 171.5 hours or 7.14 days
  
  print("SPDE interaction:")
  print(t_spdeint)
  
}


##
## Store
N_samples_posterior <- 1000 
dat$deaths_sim <- simdata[,1]
INLA_FUN(form = form_regint, 
         dataset = dat, 
         verb = FALSE,
         theta.mode = c(-1.19, 5.45, -2.37, 5.92), 
         rerun = TRUE) -> mod_regint

summary(mod_regint)
points(mod_regint$summary.random$temp_id$mean, col = "red")
plot(mod_baseline$summary.random$temp_id$mean)
summary(mod_baseline)
# here i take the joint sample as I need the intercept
inla.posterior.sample(n=1000, result = mod_regint) -> sam

lapply(sam, function(X){
  int <- X$latent[rownames(X$latent)[startsWith(rownames(X$latent), "(Intercept)")],]
  slop <- X$latent[rownames(X$latent)[startsWith(rownames(X$latent), "pol_temp")],]
  smr_pred <- exp(gridint$pm10*gridint$temperature*slop/sd.scale.store + int)
  
  list(int = int, slop = slop, smr_pred = smr_pred) %>% return()
}) -> res.samples

res.samples[[1]]$slop
sapply(res.samples, function(X) X$slop) %>% hist()


##
## Simulation metrics
# 1. MSE, 2. bias, 3. Coverage probability. 



