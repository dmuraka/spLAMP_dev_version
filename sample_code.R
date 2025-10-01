
#install.packages("sf")
library(sf)

#### Data at sample sites
data(meuse)
coords  <- meuse[,c("x","y")]
y       <- log( meuse$zinc )
x       <- data.frame(dist= meuse[,"dist"],ffreq2=ifelse(meuse$ffreq==2,1,0),
                      ffreq3=ifelse(meuse$ffreq==3,1,0))

#### Data at prediction sites
data(meuse.grid)
coords0 <-meuse.grid[,c("x","y")]
x0      <-data.frame(dist= meuse.grid$dist,ffreq2=ifelse(meuse.grid$ffreq==2,1,0),
                    ffreq3=ifelse(meuse.grid$ffreq==3,1,0))

#### Optimization of the number of spatial resolutions/bandwidths
mod_hv  <- lamp_hv(y=y, x=x, coords=coords, kernel="exp",
                   rf=TRUE) # If TRUE, random forest is additionally trained 
                            # to capture non-linear patterns and higher order interactions

mod_hv$sse_hv               # Holdout validation error (sum of squares error)
mod_hv$param$bands          # Kernel bandwidths used in the optimal model

#### Local aggregate multiscale process (LAMP) model
mod     <- lamp(y=y, x=x, coords=coords, x0=x0, coords0=coords0, mod_hv=mod_hv)
mod$beta                    # Coefficients and their standard errors
pred    <- mod$pred         # Predictive means at sample sites (first 10 elements)
pred0   <- mod$pred0        # Predictive means at prediction sites (first 10 elements)
plot(y,mod$pred);abline(0,1)# Check data fitting

#### Monte Carlo simulation (uncertainty modeling)
mod_sim <- lamp_sim(mod, nsim=200)
pred0_ms<- mod_sim$pred0    # Predictive mean and standard deviation 
pred0_qt<- mod_sim$pred0_qt # Predictive quantiles (2.5%, 50%, 97.5%)

##### Prediction results
pdat    <- data.frame(coords0,
                     pred    = pred0,           # Predictive mean
                     pred_sd = pred_ms$pred_se, # Predictive standard deviation 
                     pred_l95= pred_qt$q0.025,  # Predictive 2.5 percentile
                     pred_u95= pred_qt$q0.975,  # Predictive 97.5 percentile
                     pred_len95= pred_qt$q0.975-pred_qt$q0.025) # Width of the 95 % credible interval

# Spatial plot
coordinates(pdat)<-c("x","y")
pdat_sf <-st_as_sf(pdat)
plot(pdat_sf[,c("pred")]   , nbreaks=20, pch=15, axes=TRUE, key.pos = 1)
plot(pdat_sf[,c("pred_sd")], nbreaks=20, pch=15, axes=TRUE, key.pos = 1)
plot(pdat_sf[,c("pred_len95")], nbreaks=20, pch=15, axes=TRUE, key.pos = 1)
plot(pdat_sf[,c("pred_l95","pred_u95")], nbreaks=20, pch=15, axes=TRUE, key.pos = 1)



