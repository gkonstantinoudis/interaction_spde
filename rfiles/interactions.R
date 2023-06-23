
##
## Interactions

library(INLA)
library(sf)
library(spdep)
library(fields)
library(dplyr)
library(ggplot2)

setwd("/Users/gkonstan/Desktop/italy/")

dat <- readRDS("ItalyData")
shp <- readRDS("shp")
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

form_nointeraction = 
  deaths ~ 
  offset(log(expected)) + 
  as.factor(hol) + 
  year_cont + 
  f(datid, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
  f(space_id, model = "bym2", graph = "W.adj", scale.model = TRUE, hyper = hyper.bym, constr = TRUE) +
  f(temp_id, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
  f(pm_id, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) 

inla(form_nointeraction, 
     data = dat, 
     family = "Poisson",
     control.family = control.family, 
     verbose = TRUE,
     num.threads = 2, 
     control.compute = list(config = TRUE), 
     control.predictor = list(link=1), 
     control.mode = list(theta=theta, restart=TRUE)
     ) -> mod1

summary(mod1)

plot(mod1$summary.random$temp_id$mean)
plot(mod1$summary.random$pm_id$mean)
plot(mod1$summary.random$datid$mean)

# looks fine

# Now I will start defining interactions

dat$pol_temp <- scale(dat$pm10_lag0_3*dat$temperature_lag0_3)

form_interaction_1 = 
  deaths ~ 
  offset(log(expected)) + 
  as.factor(hol) + 
  year_cont + 
  f(datid, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
  f(space_id, model = "bym2", graph = "W.adj", scale.model = TRUE, hyper = hyper.bym, constr = TRUE) +
  f(temp_id, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
  pm10_lag0_3 +
  pol_temp

inla(form_interaction_1, 
     data = dat, 
     family = "Poisson",
     control.family = control.family, 
     verbose = TRUE,
     num.threads = 2, 
     control.compute = list(config = TRUE), 
     control.predictor = list(link=1), 
     control.mode = list(theta=theta, restart=TRUE)
) -> mod2

summary(mod2)
plot(mod2$summary.random$temp_id$mean)
plot(mod2$summary.random$datid$mean)

# it seems that the interaction is significant

# Interaction 2
# Since the effect of air-pollution is linear we can think of it being linear in each temperature level
dat$temp_id2 <- dat$temp_id
dat$pm10_lag0_3

form_interaction_2 = 
  deaths ~ 
  offset(log(expected)) + 
  as.factor(hol) + 
  year_cont + 
  f(datid, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
  f(space_id, model = "bym2", graph = "W.adj", scale.model = TRUE, hyper = hyper.bym, constr = TRUE) +
  f(temp_id, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
  pm10_lag0_3 + 
  f(temp_id2, pm10_lag0_3, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE)

inla(form_interaction_2, 
     data = dat, 
     family = "Poisson",
     control.family = control.family, 
     verbose = TRUE,
     num.threads = 2, 
     control.compute = list(config = TRUE), 
     control.predictor = list(link=1), 
     control.mode = list(theta=theta, restart=TRUE)
) -> mod3

summary(mod3)
plot(mod$summary.random$temp_id$mean)
plot(mod4$summary.random$datid$mean)
plot(mod4$summary.random$temp_id2$mean)


# Interaction 3: iid
cut_var <- 15
dat$temp_id_int <- inla.group(dat$temperature_lag0_3, method = "quantile", n = cut_var)
dat$pm_id_int <- inla.group(dat$pm10_lag0_3, method = "quantile", n = cut_var)

dat$interaction_id <- paste0(dat$temp_id_int, dat$pm_id_int)
dat$interaction_id <- as.numeric(as.factor(dat$interaction_id))
max(dat$interaction_id) # this is not cut_var^2 as there is no data in all pixels but this is fine for the iid

form_interaction_3 = 
  deaths ~ 
  offset(log(expected)) + 
  as.factor(hol) + 
  year_cont + 
  f(datid, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
  f(space_id, model = "bym2", graph = "W.adj", scale.model = TRUE, hyper = hyper.bym, constr = TRUE) +
  f(temp_id, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
  pm10_lag0_3 + 
  f(interaction_id, model = "iid", hyper = hyper.iid, constr = TRUE)

inla(form_interaction_3, 
     data = dat, 
     family = "Poisson",
     control.family = control.family, 
     verbose = TRUE,
     num.threads = 2, 
     control.compute = list(config = TRUE), 
     control.predictor = list(link=1), 
     control.mode = list(theta=theta, restart=TRUE)
) -> mod4

summary(mod4)
plot(mod4$summary.random$temp_id$mean)
plot(mod4$summary.random$datid$mean)
mod4$summary.random$interaction_id$mean

# plot interaction 
int_df <- dat[,c("interaction_id", "temp_id_int", "pm_id_int")]
int_df <- int_df[!duplicated(int_df$interaction_id),]
plot(int_df[,c("temp_id_int", "pm_id_int")])
int_df <- left_join(int_df, mod4$summary.random$interaction_id[,c("ID", "mean", "sd")], 
                    by = c("interaction_id" = "ID"))
        
ggplot() + geom_point(data=int_df, aes(x=temp_id_int %>% as.factor() %>% as.numeric(), 
                                       y=pm_id_int %>% as.factor() %>% as.numeric(), 
                                       col=mean), size = 4) + scale_colour_viridis_c()

ggplot() + geom_point(data=int_df, 
                      aes(x=temp_id_int %>% as.factor() %>% as.numeric(), 
                          y=pm_id_int %>% as.factor() %>% as.numeric(), 
                          col=sd), size = 4) + scale_colour_viridis_c()

# Interaction 4: SPDE
# what we can also have is an spde in the interaction.
dat$pm_id_factor <- dat$pm_id_int %>% as.factor() %>% as.numeric()
dat$temp_id_factor <- dat$temp_id_int %>% as.factor() %>% as.numeric()

expand.grid(
  x = unique(dat$pm_id_factor), 
  y = unique(dat$temp_id_factor)
) -> mock_shp

plot(mock_shp$x, 
     mock_shp$y)

## ok lets generate the mesh. We can have two meshes, one to predict everywhere and the second 
# to be limited with the points available. 
locs <- unique(dat[, c("temp_id_factor", "pm_id_factor")])
colnames(locs) <- c("x", "y")
 
# mesh <- inla.mesh.2d(loc = locs, offset = c(1, 5), max.edge = 5)
# plot(mesh)
# points(mock_shp$x, 
#        mock_shp$y, col = "red")


# First I will try to predict nationwide

set.seed(2)
# mesh <- inla.mesh.2d(loc.domain = locs, max.edge = c(1, 15), n = 4,
#                      offset = c(0.5, 0.3)) 
mesh <- inla.mesh.2d(loc.domain = locs, max.edge = c(1.5, 10),
                     offset = c(2, 0.1)) 

mesh <- inla.mesh.2d(loc.domain = locs, max.edge = c(1.5, 10),
                     offset = c(3, 0.1)) 

# mesh <- inla.mesh.2d(loc.domain = locs, max.edge = c(1, 15),
#              offset = c(0.5, 0.3)) 
plot(mesh)
points(mock_shp$x, 
       mock_shp$y, col = "blue", cex = 0.5, pch = 19)
points(locs, col = "red", pch = 18)


# == projector matrix ==
# projector matrix
A <- inla.spde.make.A(mesh = mesh,
                            loc = dat[, c("temp_id_factor", "pm_id_factor")] %>% as.matrix())

# == model definition ==
# spde model
rf <- inla.spde2.pcmatern(mesh = mesh,
                          prior.range = c(5,0.5),
                          prior.sigma = c(1,0.01), 
                          constr = TRUE)


form_interaction_4 = 
  deaths ~ 
  1 +
  offset(log(expected)) + 
  as.factor(hol) + 
  year_cont + 
  f(datid, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
  f(space_id, model = "bym2", graph = "W.adj", scale.model = TRUE, hyper = hyper.bym, constr = TRUE) +
  f(temp_id, model = "rw2", hyper = hyper.iid, constr = TRUE, scale.model = TRUE) + 
  pm10_lag0_3 + 
  f(field, model = rf)

# stack 
stk <- inla.stack(data = list(deaths = dat$deaths),
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
                                 ), tag = "mod")


inla(form_interaction_4, 
     data = inla.stack.data(stk), 
     family = "Poisson",
     control.family = control.family, 
     verbose = TRUE,
     num.threads = 2, 
     control.fixed = list(prec.intercept = 5, mean.intercept = 0),
     control.compute = list(config = TRUE), 
     control.predictor = list(link=1, compute = TRUE, A = inla.stack.A(stk)), 
     control.mode = list(theta=theta, restart=TRUE)
) -> mod5

mod1$summary.fixed
mod5$summary.fixed

summary(mod4)
summary(mod5)

ggplot(data = mod5$summary.random$temp_id) + 
  geom_line(aes(x=ID, y=mean)) + 
  geom_line(aes(x=ID, y=`0.025quant`)) + 
  geom_line(aes(x=ID, y=`0.975quant`))

ggplot(data = mod5$summary.random$space_id[1:101,]) + 
  geom_line(aes(x=ID, y=mean)) + 
  geom_line(aes(x=ID, y=`0.025quant`)) + 
  geom_line(aes(x=ID, y=`0.975quant`))

ggplot(data = mod5$summary.random$datid) + 
  geom_line(aes(x=ID, y=mean)) + 
  geom_line(aes(x=ID, y=`0.025quant`)) + 
  geom_line(aes(x=ID, y=`0.975quant`))

ggplot(data = mod5$summary.random$field) + 
  geom_line(aes(x=ID, y=mean)) + 
  geom_line(aes(x=ID, y=`0.025quant`)) + 
  geom_line(aes(x=ID, y=`0.975quant`))

# the variance is very high probably due to the intercept

local.plot.field = function(field, mesh, xlim=c(0,cut_var), ylim=c(0,cut_var), ...){
  stopifnot(length(field) == mesh$n)
  # - error when using the wrong mesh
  proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(300, 300))
  # - Can project from the mesh onto a 300x300 plotting grid 
  field.proj = inla.mesh.project(proj, field)
  # - Do the projection
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
             xlim = xlim, ylim = ylim, col = plasma(101), ...)  
}

local.plot.field(mod5$summary.random[['field']][['mean']], mesh)
local.plot.field(mod5$summary.random[['field']][['sd']], mesh)



