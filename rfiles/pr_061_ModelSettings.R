#:::::::::::::::::::::::::::::
# Project: spde_interaction
# description: SPDE and linear model settings
#:::::::::::::::::::::::::::::

##
## Linear model settings
form_regint = 
  deaths_sim ~ 
  offset(log(expected)) + 
  as.factor(hol) + 
  year_cont + 
  f(datid, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
  f(space_id, model = "bym2", graph = file.path(path_output, "W.adj"), scale.model = TRUE, hyper = hyper.bym, constr = TRUE) +
  f(temp_id, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
  pm10_lag0_3 +
  pol_temp


##
## SPDE settings
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

# plot(mesh)
# points(mock_shp$x, 
#        mock_shp$y, col = "blue", cex = 0.4, pch = 19)
# points(locs, col = "red", pch = 18, cex = 0.8)

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
  f(space_id, model = "bym2", graph = file.path(path_output, "W.adj"), scale.model = TRUE, hyper = hyper.bym, constr = TRUE) +
  f(temp_id, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
  pm10_lag0_3 + 
  f(field, model = rf)



#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
