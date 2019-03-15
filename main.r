# ----------------------------------------------------------------
# Code and sample script presented in the publication:
#
# "Towards reliable extreme weather and climate event attribution"
# by Bellprat, O., Guemas, V., Doblas-Reyes, F. and Donat, M.,.G.
# 
# Nature Communications, NCOMMS-18-00404A 
#
# Sample data can be retrieved from the following source:
#
# Climate Model Simulations from the EUCLEIA Project: 
# https://catalogue.ceda.ac.uk/uuid/99b29b4bfeae470599fb96243e90cde3
#
# Observational Datasets:
# CRU 2m Temperature: https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_3.23/
# NOAA 2m Temperature: https://www.ncdc.noaa.gov/ghcnm/
#
# Additional information on sample data:
# Gridding: Linear interpolation to regular grid 100x50
# Masking: Applying ocean mask determined from CRU3.23
# Averaging: Seasonal averages from June to August
# Period: 1960 - 2010
# Ensemble members: 15 as from original data
# Anomalies: Annual climatology removed for all timeseries
#
# Date: 10 March 2019
# Author: Omar Bellprat (omar.bellprat@gmail.com)
# ----------------------------------------------------------------

rm(list=ls())
gc()

# LOAD LIBRARIES
library(SpecsVerification)
library(s2dverification)
library(RColorBrewer)
library(MASS)
library(abind)

source('functions/CalibrateAnt.R')
source('functions/CalibrateNat.R')

# LOAD SAMPLE DATA
data <- get(load('./data/sample_data.RData'))

# REMOVE CLIMATOLOGIES
cens <- Clim(data$mod,data$obs)
data$mod_ano <- Ano(data$mod,cens$clim_exp)
data$obs_ano <- Ano(data$obs,cens$clim_obs)

# RETAIN ANT - NAT SIGNAL IN CLIMATOLOGIES
data$clim_sig <- cens$clim_exp[1,,,] - cens$clim_exp[2,,,]

# DETERMINE PARAMETERS OF DATA
nobs <- dim(data$obs)[1]
nmod <- dim(data$mod)[1]
nmemb <- dim(data$mod)[2]
nyears <- dim(data$mod)[3]
nlat <- dim(data$mod)[4]
nlon <- dim(data$mod)[5]

# CALIBRATE ENSEMBLE SIMULATIONS
data$mod_cal<-array(NA,dim=c(nobs,nmod,nmemb,nyears,nlat,nlon)) # Allocate array for calibrated data

for (i in 1:nobs){ # Calibrate for each obs
  print('Calibrating hindcasts ...')

  # Select observation and simulation
  obs_tmp <- InsertDim(InsertDim(data$obs_ano[i,,,,],1,1),2,1)
  ant_tmp <- InsertDim(data$mod_ano[1,,,,],1,1)
  nat_tmp <- InsertDim(data$mod_ano[2,,,,],1,1)

  # Calibrate Antropogenic and Natural ensembles
  cal_ant <- CalibrateAnt(ant_tmp,obs_tmp)
  cal_nat <- CalibrateNat(nat_tmp,obs_tmp,alpha=cal_ant$alpha,beta=cal_ant$beta) 
  # Use same calibration parameters for the ensemble spread

  # Add again climatologies
  data$mod_cal[i,1,,,,] <- drop(cal_ant$mod) + InsertDim(cens$clim_exp[1,,,],2,nyears)
  data$mod_cal[i,2,,,,] <- drop(cal_nat$mod) + InsertDim(cens$clim_exp[2,,,],2,nyears)
}

# ILLUSTRATE CALIBRATION EFFECT
gridcell <- c(22,10)
timeframe <- seq(nyears-9,nyears)

source('functions/illustrate_calibration.r')
illustrate_calibration(data,gridcell,timeframe)









