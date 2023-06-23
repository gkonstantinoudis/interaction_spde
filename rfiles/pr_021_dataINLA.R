#:::::::::::::::::::::::::::::
# Project: spde_interaction
# description: Prepare data for INLA
#:::::::::::::::::::::::::::::

prep_dat_INLA <- function(){
  
  dat <- readRDS(file.path(path_output, "ItalyData"))
  shp <- readRDS(file.path(path_data, "shp"))
  W.nb <- spdep::poly2nb(shp)
  spdep::nb2INLA(file.path(path_output, "W.adj"), W.nb)
  theta <- NULL
  
  dat$temp_id <- inla.group(dat$temperature_lag0_3, method = "cut", n = 50)
  dat$pm_id <- inla.group(dat$pm10_lag0_3, method = "cut", n = 50)
  dat$datid <- as.numeric(as.factor(dat$date))
  dat$year_cont <-  dat$year - min(dat$year) + 1 
  dat$space_id <- as.numeric(as.factor(dat$SIGLA))
  
  # scale pollutants for interaction
  dat$pol_temp <- scale(dat$pm10_lag0_3*dat$temperature_lag0_3)
  sd.scale.store <- sd(dat$pm10_lag0_3*dat$temperature_lag0_3)
  
  # Prepare temperature and pm for spde interaction
  cut_var <- 30
  dat$temp_id_int <- inla.group(dat$temperature_lag0_3, method = "cut", n = cut_var)
  dat$pm_id_int <- inla.group(dat$pm10_lag0_3, method = "cut", n = cut_var)
  
  dat$interaction_id <- paste0(dat$temp_id_int, dat$pm_id_int)
  dat$interaction_id <- as.numeric(as.factor(dat$interaction_id))
  
  dat$pm_id_factor <- dat$pm_id_int %>% as.factor() %>% as.numeric()
  dat$temp_id_factor <- dat$temp_id_int %>% as.factor() %>% as.numeric()
  
  saveRDS(dat, file = file.path(path_output, "datfin"))
  
}
