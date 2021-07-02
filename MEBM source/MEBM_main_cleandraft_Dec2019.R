#-------------------------------------------------------------------------------#
#                Moist Energy Balance Model (MEBM) -- Initialize                # 
#                                                                               #
# Following the solution of Gerard Roe, Univ. Washington (see Roe et al., 2015) #
# Hydroclimate parameterizations after Held, 2001; from Siler et al., (in Rev.) #
#      Coded and modified in R by Tyler Kukla, Stanford Univ. June, 2018        #
#       Modified Jan 2019 by DEI, note all constants have been brough to top    #
#-------------------------------------------------------------------------------#

# req'd for hydro solve
library(tibble)
library(pracma)
# req'd for main ODE 
library(bvpSolve)
# req'd for plotting
library(ggplot2)
library(ggpubr)

setwd("~/Desktop/Planetary MEBM") # Dan's WD

rm(list=ls())                                    # clear global environment
#dev.off()                                        # Clear Figure

# bring in the solution script a
source('MEBM_solve.R')  

#... generate a list for output 
df_mod <- list()
df_mod.low <- list()
df_mod.high <- list()

# Constants
modCO2 <- 360                     # [ppmv] modern CO2 calibrated by G.Roe and M.Brandon (FIXED)
#... MEBM grid and boundary conditions 
nNodes <- 1000                     # number of nodes in the MEBM grid
xin <- seq(-0.99,0.99,length=nNodes)  # x-values (0-1 is equator-pole); (-1 to 1 version [pole to pole] is not yet working)
LBC <- RBC <- 0                   # zero-flux boundary condition
# Diffusion coefficient
D_HF10 <- 1.06e6                  # [m2 s-1] diffusivity in Hwang and Frierson 2010
D_Sil <- 1.16e6 #1.16                   # [m2 s-1] diffusivity in Siler et al., 2018
D_mod <- 1.16e6 # Dan's calibration, note Q0 is also raised to 340.5
myD <- D_mod                 # select the diffusion coefficient you want to use
Din <- myD * p0/g                 # Here, F = -2pi*Re*cos(phi)*D*dh/dy. For DHF10 F = -2pi*Re*cos(phi)*D*(ps/g)*dh/dy    
#A0 <- 207.3                      # [W m-2] OLR constant - NO CO2 DEPENDENCY
B0 <- 3.65 #2.09 is original                        # [W m-2 degC-1] OLR response NOTE UNITS
# Ice line and albedo
alf_flag <- "step"                # prescribe the albedo formulation:
# OPTIONS: "step" = discrete step from no ice to ice albedo
#          "D18" = hyperbolic (smoothed) transition from no ice to ice albedo
#                  from Dortmans et al., 2018 (Climate of the Past)
#          "prescribed" = use calculated albedo from landfrac data
# omega <- 2.73                     # only needed if alf_flag="D18"; this sets smoothness of no ice to ice transition (higher value=smoother(shallower))
# ---------------------------------------
#... define smoothing bandwidth 
default_bdwth <- 1            # baseline kernel smooth bandwidth (smoothing is unnoticable to the naked eye)
#anom_timeslices <- c(14, 25, 26, 29, 30, 38, 50, 51, 57, 60, 67, 68, 69, 73, 83,86, 91, 95, 97, 106, 108)  # times (iterations) where bandwidth should be modified
#anom_bdwth <- c(50, 250, 100, 50, 50, 10, 200, 150, 10, 200, 300, 200, 100, 10, 10, 10, 10, 10, 10, 10, 100)       # the value we modify to when xin=-.99 to.99 and nNodes=2500
bdwth <- rep(default_bdwth, 110)
bdwth[anom_timeslices] <- anom_bdwth
# ---------------------------------------
ice_Threshold <- -10 #-10, original              # [degC] threshold for ice formation
landAlf <- 0.38       #0.38           # [] average albedo of ice-free land surface
oceanAlf <- 0.13      #0.13             # [] average albedo of ocean
iceAlf <- 0.7                     # [] average albedo of ice 
rel_hum <- 0.8                    # [] relative humidity

# *************************************************************** #
# WRITE SIMULATION LOG SHEET ------------------------------------ #
# source('LogFile_writer.R')
# writeLogFile(userName="T Kukla and DE Ibarra", userNotes="Simulation stepping through the Phanerozoic in 5 myr time slices ONLY CHANGING pCO2. Paleogeography is held CONSTANT AT the youngest timeslice (modern). Model is forced with phanerozoic pco2 data. Changing solar constant with time. No milankovitch forcing. Incorporates a new method for smoothing continental land fraction (ksmooth function). Had to update iterations 67-69 to higher bandwidths. Made no attempt to lower bandwidths that were set in run001. Solved using the bvpcol function.")
# *************************************************************** #

landFrac <- rep(0.7,nNodes) #myLandFrac[[1]]$LandFrac_inLatBelt_forAlbedo
#lat_orig <- myLandFrac[[1]]$Lat


#... set greenhouse forcing 
pco2 <- 300                # [ppmv] climatological CO2 for the simulation 

#... Physical parameters input to MEBM
# Radiation parameterization with dimmer sun
t = 4567 # Modern
Fs_mult = 1/(1+(2/5)*(1-t/4567))
Q0 <- 340.5*Fs_mult #338.5 # 338.5original (should be more like 340.5 Dan thinks)                       # [W m-2]  solar constant n.b. 4Q0 = 1354, which is low...
A0 <- 207.3-5.35*log(pco2/modCO2) # [W m-2] OLR constant depending on pCO2 [Myhre et al., 1998, North and Kim, 2017]

# Bring in land fraction data and collapse hemispheres
#modland <- read.csv("modernland_collapse.csv")
landFrac_interp <- approx(x=landFrac, n=nNodes)[[2]]    # the interpolated fraction of land per lat
landFrac_int <- ksmooth(x=seq(1,length(xin), by=1), y=landFrac_interp, bandwidth = default_bdwth)$y

alf_bylat <- (landFrac_int * landAlf) + ((1-landFrac_int) * oceanAlf)  # [] latitudinal albedo weighted by continent
alf_noice <- (landFrac_int * landAlf) + ((1-landFrac_int) * oceanAlf)  # [] latitudinal albedo weighted by continent
afun_df <- as.data.frame(list(x = xin, alf_noice = alf_noice))
afun <- approxfun(afun_df, rule=2)

# output dataframe
df_mod <- MEBM_solve()


# plot some  example results
par(mfrow=c(2,3))
plot(xin,df_mod$temp_C,ylab="Temp",type="l")
plot(xin,df_mod$latentFlux_W,ylab="Latent Heat Flux",type="l")
plot(xin,df_mod$MSE_J_kg,ylab="Moist Static Energy",type="l")
plot(xin,df_mod$Precip_m_yr,ylab="Precip",type="l")
plot(xin,df_mod$EminusP_m_yr,ylab="E-P",type="l")
plot(xin,df_mod$albedo,ylab="Albedo",type="l")



