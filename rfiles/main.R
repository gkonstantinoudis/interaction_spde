#:::::::::::::::::::::::::::::
# Project: spde_interaction
# description: Main file
#:::::::::::::::::::::::::::::


controls = list(
  cleandata = FALSE, 
  checklinear = FALSE, 
  prepare_data = TRUE, # This should be always true
  baseline_effects = FALSE,
  sim_data = TRUE,
  checktimes = FALSE,
  runmodel = TRUE, 
  simmetrics = TRUE
)

# set-up
setwd("C:/Users/gkonstan/OneDrive - Imperial College London/ICRF Imperial/Projects/Interactions simulation/")
source("rfiles/da_setup.R")


#######
# Block 1: Data wrangling ----
#######

if(controls$cleandata){ 
  ##
  ## Need to do this
}


#######
# Block 2: Prepare and load the data for INLA ----
#######

if(controls$prepare_data){ 
  
  if(!file.exists(file.path(path_output, "datfin"))){
    source(file.path(path_rfiles, "pr_021_dataINLA.R"))
    prep_dat_INLA()
    dat <- readRDS(file.path(path_output, "datfin"))
  }else{
    dat <- readRDS(file.path(path_output, "datfin"))
  }
  
}



#######
# Block 3: Check linearity of air-pollution in the actual data to decide how to include it in the simulation----
#######

if(controls$checklinear){ 
  
  form_rw2_pm = 
    deaths ~ 
    offset(log(expected)) + 
    as.factor(hol) + 
    year_cont + 
    f(datid, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
    f(space_id, model = "bym2", graph = file.path(path_output, "W.adj"), scale.model = TRUE, hyper = hyper.bym, constr = TRUE) +
    f(temp_id, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
    f(pm_id, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) 
  
  INLA_FUN(form = form_rw2_pm, dataset = dat, verb = TRUE)-> mod__rw2_pm

  ggplot(mod__rw2_pm$summary.random$pm_id %>% as.data.frame()) + 
    geom_line(aes(x=ID, y=exp(`0.5quant`))) + 
    geom_line(aes(x=ID, y=exp(`0.975quant`))) +
    geom_line(aes(x=ID, y=exp(`0.025quant`))) + 
    theme_bw() # the effect seems linear so I will use just a line
  
}


#######
# Block 4: Run baseline model to get the truth of the effects ----
#######

if(controls$baseline_effects == TRUE){
  
  if(!file.exists(file.path(path_output, "baselineINLAmodel"))){
  ##
  ## Fit baseline model based on linearity on PM10
  form_baseline = 
    deaths ~ 
    offset(log(expected)) + 
    as.factor(hol) + 
    year_cont + 
    f(datid, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
    f(space_id, model = "bym2", graph = file.path(path_output, "W.adj"), scale.model = TRUE, hyper = hyper.bym, constr = TRUE) +
    f(temp_id, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
    pm10_lag0_3
  
  INLA_FUN(form = form_baseline, dataset = dat)-> mod_baseline
  saveRDS(mod_baseline, file = file.path(path_output, "baselineINLAmodel"))
  
  dat <- left_join(dat, data.frame(id_tmp_ing = mod_baseline$summary.random$temp_id$ID, 
                                   ranef_temp = mod_baseline$summary.random$temp_id$`0.5quant`), 
                   by = c("temp_id" = "id_tmp_ing"))
  
  saveRDS(dat, file.path(path_output, "datfin"))
  
  }else{
    
  mod_baseline <- readRDS(file.path(path_output, "baselineINLAmodel"))
  
  }
  
}


#######
# Block 5: Simulate data ----
#######

if(controls$sim_data){ 
  
  if(!file.exists(file.path(path_output, "baselineINLAmodel"))){
    source(file.path(path_rfiles, "pr_051_SimData.R"))
    Simulate_data_function()
  }else{
    simdata <- readRDS(file.path(path_output, "DataSimulation"))
  }
  
}


#######
# Block 6: Check how much time the simulation will take ----
#######

if(controls$checktimes){ 
  source(file.path(path_rfiles, "pr_061_ModelSettings.R"))
  ##
  ## linear model
  dat$deaths_sim <- simdata$sc1[,1]
  t_0 <- Sys.time()
  INLA_FUN(form = form_regint, 
           dataset = dat, 
           verb = FALSE,
           theta.mode = NULL, 
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
  stk <- stk.function(L = 1, scenario = "sc1")
  
  t_0 <- Sys.time()
  INLA_FUN(form = form_spde_int, 
           dataset = inla.stack.data(stk), 
           Alphamat = inla.stack.A(stk), 
           verb = FALSE,
           theta.mode = NULL, 
           rerun = TRUE) -> mod_spdeint
  t_1 <- Sys.time()
  t_spdeint <- t_1 - t_0
  
  # TIME WITHOUT RERUN: ~3.58 minutes
  # TIME WITH RERUN: ~6.86 minutes
  # ~6.86*300*5 = 10290 minutes or 171.5 hours or 7.14 days
  
  print("SPDE interaction:")
  print(t_spdeint)
  
  theta_list <-
    list(
      theta_regint = mod_regint$mode$theta, 
      theta_spdeint = mod_spdeint$mode$theta
    )
  
  saveRDS(theta_list, file = file.path(path_output, "theta_list"))
  
}


#######
# Block 7: Run the simulation ----
#######

if(controls$runmodel){ 
  
  if(!file.exists(file.path(path_output, "theta_list"))){
    theta_list <- NULL
  }else{
    theta_list <- readRDS(file.path(path_output, "theta_list"))
  }
  
  source(file.path(path_rfiles, "pr_061_ModelSettings.R"))
  source(file.path(path_rfiles, "pr_071_ModelFit.R"))
  
  sc <- paste0("sc", 1:length(simdata))

  LinearInteractionPar <- function(K) LinearInteraction(scenario = sc[j], 
                                                        L=K, 
                                                        N_samples_posterior = 300, 
                                                        theta_mode = theta_list[[1]])
  # t_0 <- Sys.time()
  # tmp_lint <- LinearInteractionPar(K = 1)
  # t_1 <- Sys.time()
  # t_1 - t_0 # 5.260115 mins
  
  SPDEInteractionPar <- function(K) SPDEInteraction(scenario = sc[j], 
                                                    L=K, 
                                                    N_samples_posterior = 300, 
                                                    theta_mode = theta_list[[2]])
  # t_0 <- Sys.time()
  # tmp_spde <- SPDEInteractionPar(K = 1)
  # t_1 <- Sys.time()
  # t_1 - t_0 # 33 minutes
  
  # no_cores = 20
  # M <- 300 # this is the number of simulated datasets to run the models on.
  
  no_cores = 5
  M <- 5 # this is the number of simulated datasets to run the models on.
  
  t_0 <- Sys.time()
  for(j in 1:length(simdata)){
    
    print(sc[[j]])

    print("Launch model on parallel")
    cl <- makeCluster(no_cores) 
    
    clusterExport(cl, c("LinearInteraction", "SPDEInteraction", "sc", "dat", "simdata",
                        "stk.function", "A", "rf", "form_spde_int", "form_regint", "j", "theta_list", 
                        "path_output", "INLA_FUN", "hyper.iid", "hyper.bym", "control.family"))
    
    clusterEvalQ(cl, {
      library(dplyr)
      library(INLA)
    })
    
    print("Running the linear model")
    result_linear <- parLapply(cl, 1:M, LinearInteractionPar)  
    print("Running the SPDE model")
    result_SPDE <- parLapply(cl, 1:M, SPDEInteractionPar)  
    
    stopCluster(cl) 
    gc()
  }
  t_1 <- Sys.time()
  print(t_1 - t_0)
  
  # no_cores = 5
  # M <- 5, and 5 scenarios took ~4h
  
  saveRDS(result_linear, file = file.path(path_output, "result_linear"))
  saveRDS(result_SPDE, file = file.path(path_output, "result_SPDE"))

}


#######
# Block 8: Extract results and calculate simulation metrics ----
#######


if(controls$simmetrics){ 
  
  result_linear <- readRDS(file.path(path_output, "result_linear"))
  result_SPDE <- readRDS(file.path(path_output, "result_SPDE"))
  gridint <- readRDS(file.path(path_output, "gridint"))
  
  result_linear[[1]][[1]]$smr_pred %>% length()
  result_SPDE[[1]][[1]]$smr_pred %>% length()
}



