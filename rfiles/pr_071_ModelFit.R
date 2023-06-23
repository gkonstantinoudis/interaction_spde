#:::::::::::::::::::::::::::::
# Project: spde_interaction
# description: Run models on simulated data
#:::::::::::::::::::::::::::::



# Linear interaction model
LinearInteraction <- function(scenario, L, N_samples_posterior = 1000, theta_mode){
  inla.setOption(inla.timeout=1000) # timeout error ~15min
  
  tryCatch({
    
  gridint <- readRDS(file.path(path_output, "gridint"))
  sd.scale.store <- sd(dat$pm10_lag0_3*dat$temperature_lag0_3)
  
  dat$deaths_sim <- simdata[[scenario]][,L]
  
  INLA_FUN(form = form_regint, 
           dataset = dat, 
           verb = FALSE,
           theta.mode = theta_mode, 
           rerun = TRUE) -> mod
  
  # here i take the joint sample as I need the intercept too
  inla.posterior.sample(n=N_samples_posterior, result = mod) -> sam
  
  lapply(sam, function(X){
    int <- X$latent[rownames(X$latent)[startsWith(rownames(X$latent), "(Intercept)")],]
    slop <- X$latent[rownames(X$latent)[startsWith(rownames(X$latent), "pol_temp")],]
    smr_pred <- exp(gridint$pm10*gridint$temperature*slop/sd.scale.store + int)
    rr_pred <- smr_pred/exp(int)
    list(int = int, slop = slop, smr_pred = data.frame(pm10 = gridint$pm10, 
                                                       temperature = gridint$temperature, 
                                                       ID = gridint$ID,
                                                       smr_pred = smr_pred, 
                                                       rr_pred = rr_pred)
         ) %>% return()
  }) -> ret
  ret %>% return()
  }, error = function(e){
    if (grepl("reached elapsed time limit|reached CPU time limit", e$message)){
      -999
    }else{
      # error not related to timeout
      -999
    }
  }
  )
}

# SPDE interaction model
SPDEInteraction <- function(scenario, L, N_samples_posterior = 1000, theta_mode){
  inla.setOption(inla.timeout=1000) # timeout error ~15min
  gridint <- readRDS(file.path(path_output, "gridint"))
  
  tryCatch({
  stk <- stk.function(L = L, scenario = scenario)
  
  INLA_FUN(form = form_spde_int, 
           dataset = inla.stack.data(stk), 
           Alphamat = inla.stack.A(stk), 
           verb = FALSE,
           theta.mode = theta_mode, 
           rerun = TRUE) -> mod
  
  # here i take the joint sample as I need the intercept too
  inla.posterior.sample(n=N_samples_posterior, result = mod) -> sam
  pgrid0 <- inla.mesh.projector(mesh, loc = coordinates(cbind(gridint$pm10, gridint$temperature)))
  
  lapply(sam, function(X){
    int <- X$latent[rownames(X$latent)[startsWith(rownames(X$latent), "(Intercept)")],]
    field <- X$latent[rownames(X$latent)[startsWith(rownames(X$latent), "field")],]
    smr_pred <- exp(int + field)
    rr_pred <- exp((int + field)/int)
    prd0.smr <- inla.mesh.project(pgrid0, smr_pred)
    prd0.rr <- inla.mesh.project(pgrid0, rr_pred)
    prd0.logfield <- inla.mesh.project(pgrid0, field)
    
    list(int = int, field = field, smr_pred = prd0.smr, logfield = prd0.logfield, rr = prd0.rr) %>% return()
  }) -> ret
  ret %>% return()
  }, error = function(e){
    if (grepl("reached elapsed time limit|reached CPU time limit", e$message)){
      -999
    }else{
      # error not related to timeout
      -999
    }
  }
  )
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
