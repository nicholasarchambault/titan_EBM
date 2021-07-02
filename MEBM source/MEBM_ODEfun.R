#-------------------------------------------------------------------------------#
#                Moist Energy Balance Model (MEBM) -- Solver                    # 
#                                                                               #
# Following the solution of Gerard Roe, Univ. Washington (see Roe et al., 2015) #
#      Coded and modified in R by Tyler Kukla, Stanford Univ. June, 2018        #
#                     Modified by Daniel Ibarra (Stanford)                      #
#-------------------------------------------------------------------------------#

# function to solve for the moist static energy profile 
MEBMfun<- function(x,y,pars) {
  ## NOTE: This function allows the bvp solver to step through space to reach steady state
  ## Function solves 2 related ODEs, denoted by y[1], dy1 and y[2], dy2 where: 
           # y[1] ; dy1 = temperature (degC) ; dtemp/dx
           # y[2] ; dy2 = Moist static energy (MSE) flux ; dFlux/dx
  
  # OLR parameters and diffusion coefficient
  A <- A0       
  B <- B0
  D <- Din
  
  # saturation vapor pressure
  e_sat = e0*exp(a*y[1]/(b+y[1]))
  
  # useful factor in equations [ denominator is cp + q*dq/dx ]
  fac1 = 1/(cp + (eps*rel_hum*Lv*e_sat/p0)*(a*b/(b+y[1])^2))
  
  # set albedo 
  if(alf_flag == "D18"){  # calculate albedo using the hyperbolic formulation of Dortmans et al., 2018 (Clim of the Past)
    # omega value determines how smooth the ice - land transition 
    #... (larger = shallower; smaller = steeper)
    alf <- 0.5 * ( (afun(x) + iceAlf) + (afun(x) - iceAlf) * tanh((y[1])/omega) )
  }
  
  if(alf_flag == "step"){ # if the temperature is below the ice-threshold, declare it ice
                          # Classic ice-threshold value is -10 degC (goes back to calibrations of Budyko)
    if(y[1] <= ice_Threshold){ alf <- iceAlf
    } else{alf <- afun(x)}
  }
  
  # Solve the Top Of Atmosphere (TOA) fluxes
  S <- (1-0.241*(3*x^2-1))               # insolation function. Legendre Polynomial approximation 
                                           # Change this term for Milankovitch experiments 
  
  SrcSnk <- Q0*S*(1-alf) - (A+B*y[1])    # [W m-2] Net forcing: ASR - OLR in Donohoe language
  
  # Calculate dTemp/dx
    #... At x = 1, there is a singularity in the ode, which L'hopital's rule can be used to overcome
  if(x==1){
   dy1 <- 1/(4*pi*D)*fac1*2*pi*Re^2*SrcSnk
  } else{ 
   dy1 <- -1/(2*pi*D)*1/(1-x^2)*fac1*y[2] 
  }
  
  # Calculate dFlux/dx
  dy2 <- 2*pi*Re^2*SrcSnk
  
  # output 
  return(list(c(dy1, dy2)))
}






