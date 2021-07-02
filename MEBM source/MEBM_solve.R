#-------------------------------------------------------------------------------#
#                Moist Energy Balance Model (MEBM) -- Solver                    # 
#                                                                               #
# Following the solution of Gerard Roe, Univ. Washington (see Roe et al., 2015) #
# Hydroclimate parameterizations after Held, 2001; from Siler et al., (in Rev.) #
#      Coded and modified in R by Tyler Kukla, Stanford Univ. June, 2018        #
#      Modified Dec 2018 by D. Ibarra                                           #
#-------------------------------------------------------------------------------#

source('MEBM_constants.R')    # physical constants
source('MEBM_ODEfun.R')       # MEBM ODE to get mse
source('MEBM_hydrofun.R')     # solve hydrologic fluxes

MEBM_solve <- function(){
#---------------------------------CALCULATE THE MODEL SOLUTION------------------------------------------#
  #... run the MEBM
        # yini ; yend -- have structure =c(y1=temp, y2=flux)
        # x -- grid
        # func -- bvp function
        # guess -- initial guess for temperature profile
        # atol -- calibrated based on Matlab output of similar model and computation time
        # method -- "lsoda" is fastest... if specifying a jacobian matrix this might not be the fastest 
  print(
    system.time(
      # sol <- bvpshoot(yini = c(NA, LBC), yend = c(NA, RBC), x = xin, func = MEBMfun,guess=c(temp_guess),
      #                 verbose=F, atol=1e-4, maxiter = 2000, method="lsoda",maxsteps=10^5)
      # sol <- bvptwp(yini = c(NA, LBC), yend = c(NA, RBC), x = xin, func = MEBMfun,atol=1e-14,
      #               nmax=500) #yguess=c(temp_guess),
      sol <- bvpcol(yini = c(NA, LBC), yend = c(NA, RBC), x = xin, func = MEBMfun, 
                    atol=3, nmax=1000000, xguess=c(-.99, 0.99),
                    yguess= matrix(nrow = 2, ncol = 2, data=c(5,5,0,0), byrow=T)) #yguess=c(temp_guess),
      #                 verbose=F, atol=1e-14) #, maxiter = 1000)#, method="lsoda")#,maxstep=10^5)
      
    )
  )
  #diagnostics(sol)
  
  #... solve the hydrologic cycle 
  myTemp <- sol[,2] 
  hydro <- hydrofun(temp=myTemp, x=xin)
  
  #... re-calculate MEBM variables
  # Solve the Top of Atmosphere (TOA) fluxes
  # Source  (ASR)
  if(alf_flag == "D18"){  # calculate albedo using the hyperbolic formulation of Dortmans et al., 2018 (Clim of the Past)
    # omega value determines how smooth the ice - land transition 
    #... (larger = shallower; smaller = steeper)
    alf <- 0.5 * ( (alf_noice + iceAlf) + (alf_noice - iceAlf) * tanh((myTemp)/omega) )
    }
  
  if(alf_flag == "step"){ # if the temperature is below the ice-threshold, declare it ice
    # Classic ice-threshold value is -10 degC (goes back to calibrations of Budyko)
    ice_dex <- which(myTemp <= ice_Threshold)    # find the ice line
    alf <- alf_noice ; alf[ice_dex] <- iceAlf    # calculate albedo
  }
  Src <- Q0*(1-0.241*(3*xin^2-1))*(1-alf)   # [W m-2] follows the 2nd order Legendre polynomial approximation
  # Sink
  Snk <- (A0+B0*myTemp)                     # [W m-2] OLR 
  
  # MSE divergence 
  divF  <- 1/(2*pi*Re^2)*gradient(hydro$mseFlux_W,xin)
  divF[1:2] <- divF[3]
  divF[(length(xin)-1):length(xin)] <- divF[length(xin)-2]
  
  
  # entropy production (integrate divergence)
  sdot = trapz(xin, divF)    
  
  
  #... generate an output dataframe 
  myLat <- rad2deg(asin(xin))
  mydf <- as_tibble(cbind(myLat, hydro, alf, Src, Snk, divF, sdot))
  colnames(mydf) <- c("Lat", "temp_C", "MSE_J_kg", "Precip_m_yr", "Evap_m_yr", "EminusP_m_yr", "mseFlux_W", 
                      "latentFlux_W", "sensibleFlux_W", "latentEddyFlux_W", "latentHadleyFlux_W", "LHflux_Div", "albedo", 
                      "ASR_W_m2", "OLR_W_m2", "MSEfluxDiv", "entropyProd")
  
  return(mydf)
}