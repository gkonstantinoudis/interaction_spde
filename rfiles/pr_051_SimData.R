#:::::::::::::::::::::::::::::
# Project: spde_interaction
# description: Simulate data
#:::::::::::::::::::::::::::::

Simulate_data_function <- function(){
  
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
  dat$rr_truth_sc2 <- dat$rr_truth_sc3 <- dat$rr_truth_sc4 <- dat$rr_truth_sc5 <- 0
  
  lptosim <- list(
    sc1 = lp_baseline + dat$rr_truth_sc1, 
    sc2 = lp_baseline + dat$rr_truth_sc2,
    sc3 = lp_baseline + dat$rr_truth_sc3, 
    sc4 = lp_baseline + dat$rr_truth_sc4,
    sc5 = lp_baseline + dat$rr_truth_sc5
  )
  
  t_0 <- Sys.time()
  set.seed(11)
  lapply(
    lptosim, 
    function(Y){
      lapply(1:300, 
             function(i) rpois(n = Y %>% length(), 
                               lambda = Y %>% exp()) %>% return()) %>% 
        do.call(cbind, .) %>% 
        return()
    }
  ) -> simdata
  t_1 <- Sys.time()
  t_1 - t_0 # ~ 30 seconds
  
  saveRDS(simdata, file.path(path_output, "DataSimulation"))
  saveRDS(gridint, file.path(path_output, "gridint"))
  
}



#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
