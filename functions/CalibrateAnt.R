CalibrateAnt <- function(ano_exp, ano_obs, alpha=NULL, beta=NULL, gamma=NULL) {
  m_mod <- Mean1Dim(ano_exp,2) # Ensemble mean of model
  m_obs <- Mean1Dim(ano_obs,2) # Ensemble mean of obserations

  trend_mod <- Trend(m_mod,posTR=2)
  trend_obs <- Trend(m_obs,posTR=2)

  m_mod_dt <- trend_mod$detrend
  m_obs_dt <- trend_obs$detrend

  diff_mod <- InsertDim(m_mod - m_mod_dt,2,dim(ano_exp)[2])

  # Ensemble inflation of a forecast 
  m_mod_conf_dt <- InsertDim(m_mod_dt,2,dim(ano_exp)[2])
  m_mod_conf <- InsertDim(m_mod,2,dim(ano_exp)[2])

  ediff <- ano_exp-m_mod_conf # Deviation of members from the ensemble mean

  if (is.null(alpha)){
    rho <- s2dverification::Corr(m_mod_dt,m_obs_dt, posloop = 1, poscor = 2) # Ensemble mean correlation
    rhot <- rho[,,2,,] # Select only central correlation value
    sr <- apply(m_obs_dt,FUN=sd,c(1,3,4)) # Standard deviation of observations
    dr <- InsertDim(sr,1,dim(ano_exp)[1])
    sr_conf <- dr[,1,,] # Conformed to number of model systems
    sem <- apply(m_mod_dt,FUN=sd,c(3,4)) # Standard deviation of ensemble mean
    se <- apply(ediff,FUN=sd,c(4,5)) # Standard deviation of the enemble spread
    arho <- abs(rhot) # Absolute value of correlation

    # Calibration parameters
    alpha <- arho*sr_conf/sem
    beta <- sqrt(1-rhot^2)*sr_conf/se
    gamma <- trend_obs$trend[,2,,]/trend_mod$trend[,2,,]
  } else {
    alpha <- InsertDim(alpha,1,1)
    beta <- InsertDim(beta,1,1)
    gamma <- InsertDim(gamma,1,1)
  }

  # Conform calibration parameters to experiments 
  alpha_conf <- InsertDim(InsertDim(InsertDim(alpha,1,dim(ano_exp)[2]),2,dim(ano_exp)[3]),1,1)
  beta_conf <- InsertDim(InsertDim(InsertDim(beta,1,dim(ano_exp)[2]),2,dim(ano_exp)[3]),1,1)
  gamma_conf <- InsertDim(InsertDim(InsertDim(gamma,1,dim(ano_exp)[2]),2,dim(ano_exp)[3]),1,1)

  # Calibration forecasts
  mod_cal <- alpha_conf*m_mod_conf_dt+beta_conf*ediff+gamma_conf*diff_mod
  cal <- list()
  cal$mod <- mod_cal
  cal$alpha <- alpha
  cal$beta <- beta
  cal$gamma <- gamma
  cal$rho <- rhot
  cal$trd_m <- trend_mod$trend[,2,,]
  cal$trd_o <- trend_obs$trend[,2,,]
  return(cal)
} 

