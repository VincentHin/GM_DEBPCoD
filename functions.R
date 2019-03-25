if (!exists("par.defaults")) par.defaults <- par()
par.defaults$cin <- par.defaults$cra <- par.defaults$csi <- par.defaults$cxy <- par.defaults$din <-  par.defaults$page <- NULL

# Load some packages
library('stringr'); library("RColorBrewer"); library('PSPManalysis'); library('readr'); library('tibble'); library('gridExtra')
library('reshape2');  library('ggplot2'); library('plyr'); library('tidyr'); library('Hmisc'); library('ggExtra')
library('scales')

## : -> GGOPTS ##
ggopts <- theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.key.height = unit(9, "pt"),
        legend.position = 'bottom',
        legend.box = "horizontal",
        strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(size = 9, hjust = 0.5))

## : -> PL.OPTS ##
pl.opts <- list(
  lwidth = 3,
  cex.text = 0.65,
  cex.form = 0.65,
  cex.labs = 0.75,
  cex.axs = 0.65)

#### : -> DEFINE NAMES & LABELS ##
names.lh.R0 <- c('par1', 'par2', 'age', 'reserves', 'totalmort', 'IDminloglifespan', 
                 'IDagefirstrepro', 'IDagefirstweaning', 'IDreprofemale', 'IDR0estimated', 'IDearlydeaths',
                 'meanIBI', 'sdIBI', 'minIBI', 'maxIBI', 'obsIBI',
                 'meanIWI', 'sdIWI', 'minIWI', 'maxIWI', 'obsIWI')
names.lh <- c('time', 'resource', 'number', 'age_yrs',
              'reserves', 'totalmort',
              'IDlength', 'IDbones', 'IDweight', 'IDweightM', 'IDfatratio',
              'IDingest', 'IDringest', 'IDmingest', 'IDmaint', 'IDgrowth', 'IDpregcosts', 'IDlactcosts', 'IDnetenergy',
              'IDmortality', 'IDbackground', 'IDstarvation', 'IDminloglifespan',
              'IDpregnant', 'IDpregstart', 'IDpregclock', 'IDlactating', 'IDcalved', 'IDweaned',
              'IDagefirstrepro', 'IDagefirstweaning', 'IDreprofemale', 'IDR0estimated', 'IDearlydeaths', 'pregthreshold',
              'number_calf', 'age_yrs_calf',
              'reserves_calf', 'totalmort_calf', 'IDlength_calf', 'IDbones_calf', 'IDweight_calf', 'IDweightM_calf', 'IDfatratio_calf',
              'IDingest_calf', 'IDringest_calf', 'IDmingest_calf', 'IDmaint_calf', 'IDgrowth_calf', 'IDpregcosts_calf', 'IDlactcosts_calf', 'IDnetenergy_calf',
              'IDmortality_calf', 'IDbackground_calf', 'IDstarvation_calf', 'IDminloglifespan_calf')
names.lh.state <- c('survival', 'age_days', 'reserves', 'totalmort',
                    'IDpregstart', 'IDpregclock', 'IDweaned', 'IDcalved', 'IDpregnant', 'IDlactating',
                    'IDlength','IDbones','IDweight','IDweightM','IDfatratio',
                    'IDingest', 'IDringest', 'IDmingest', 'IDmaint', 'IDgrowth', 'IDpregcosts', 'IDlactcosts', 'IDnetenergy',
                    'IDmortality','IDbackground','IDstarvation',
                    'IDagefirstrepro', 'IDagefirstweaning', 'IDreprofemale', 'IDR0estimated', 'IDearlydeaths', 'IDminloglifespan')

## Set plotting labels
var_labels <- c(IDreprofemale = 'Lifetime \nreproductive output',
                IDagefirstrepro_zero = 'Age at first \nreproduction',
                PercEarlyDeath = 'Percentage pre-weaned \ncalf deaths',
                age = 'Age at death',
                reserves = 'Reserve density \nat death',
                Reproductive = 'Proportion \nreproducing')

labels_dist_resamp <- c(Summer_0 = 'No seasonality',
                        Summer_0.25 = 'Seasonality = 0.25\nSummer disturbance',
                        Winter_0.25 = 'Seasonality = 0.25\nWinter disturbance')
labels_2dist <- c('0' = 'No disturbance',
                  '30' = '30 days disturbance')
labels_3dist <- c('0' = 'No disturbance',
                  '20' = '20 days disturbance',
                  '30' =  '30 days disturbance')
labels_dists <- c('0' = 'No disturbance',
                  '20' = '20 days summer\ndisturbance',
                  '30' =  '30 days summer\ndisturbance')
labels_distw <- c('0' = 'No disturbance',
                  '20' = '20 days winter\ndisturbance',
                  '30' =  '30 days winter\ndisturbance')
labels_resamp <- c('0' = 'No seasonality',
                   '0.25' = 'Seasonality = 0.25')

## Read life history output from a bunch of folders and put everything in a list.
import.lh.output <- function(folder = NULL, parameter = NULL, outnms = names.lh, statenms = names.lh.state){
  if(is.null(folder)) return(cat("Specify folder from which to grab output"))
  if(is.null(parameter)) return(cat("Specify parameter from which to grab output"))
  nm <- file.path(folder, parameter)
  
  out <- list()
  out$out <- read.table(paste0(nm,'.out'), col.names = outnms)
  out$esf <- read.table(paste0(nm, '.esf'), skip=2, col.names = statenms)
  out$rep <- readLines(paste0(nm, '.rep'))
  out$R0 <- read.table(paste0(nm,"_R0.dat"), col.names = names.lh.R0)
  
  return(out)
}
##

## Transposes columns IDagefirstrepro, IDagefirstweaning, IDagefirstreceptive, IBI and IWI from days to years
## Calculate some additional columns that exclude zero. Calculates PercWeaned from IDearlydeaths.
days.to.years <- function(input.df = NULL){
  if(is.null(input.df)) return('Please specify input dataframe')
  
  input.df$IDagefirstrepro <- input.df$IDagefirstrepro / 365
  input.df$IDagefirstrepro_zero <- input.df$IDagefirstrepro
  input.df$IDagefirstrepro_zero[input.df$IDagefirstrepro_zero == 0] <- NA
  input.df$Reproductive <- as.numeric(input.df$IDagefirstrepro > 0)
  if(length(input.df$age) > 0) input.df$age_yrs <- input.df$age / 365
  if(!is.null(input.df$IDearlydeaths)) input.df$PercWeaned <- (1 - input.df$IDearlydeaths / (2 * input.df$IDreprofemale + input.df$IDearlydeaths))
  
  if(length(input.df$IDagefirstweaning) > 0){
    input.df$IDagefirstweaning <- input.df$IDagefirstweaning / 365
    input.df$IDagefirstweaning_zero <- input.df$IDagefirstweaning
    input.df$IDagefirstweaning_zero[input.df$IDagefirstweaning_zero == 0] <- NA
  }
  
  if(length(input.df$IDagefirstreceptive) > 0){
    input.df$IDagefirstreceptive <- input.df$IDagefirstreceptive / 365
    input.df$IDagefirstreceptive_zero <- input.df$IDagefirstreceptive
    input.df$IDagefirstreceptive_zero[input.df$IDagefirstreceptive_zero == 0] <- NA
  }
  
  idx.zero <- match(c('meanIBI', 'sdIBI', 'minIBI', 'maxIBI', 'meanIWI', 'sdIWI', 'minIWI', 'maxIWI'), names(input.df))
  if(!is.na(idx.zero)[1]){
    input.df[,idx.zero] <- input.df[,idx.zero] / 365
    input.df$meanIBI_zero <- input.df$meanIBI
    input.df$minIBI_zero <- input.df$minIBI
    input.df$maxIBI_zero <- input.df$maxIBI
    input.df$meanIWI_zero <- input.df$meanIWI
    input.df$minIWI_zero <- input.df$minIWI
    input.df$maxIWI_zero <- input.df$maxIWI
    
    input.df$meanIBI_zero[input.df$meanIBI_zero == 0] <- NA
    input.df$minIBI_zero[input.df$minIBI_zero == 0] <- NA
    input.df$maxIBI_zero[input.df$maxIBI_zero == 0] <- NA
    input.df$meanIWI_zero[input.df$meanIWI_zero == 0] <- NA
    input.df$minIWI_zero[input.df$minIWI_zero == 0] <- NA
    input.df$maxIWI_zero[input.df$maxIWI_zero == 0] <- NA
  }
  
  return(input.df)
}
##

## Set reproductive status of a dataframe of life history output
set.status <- function(input.df = NULL){
  if(is.null(input.df)) return('Please specify an input dataframe with input.df')
  
  if(is.null(input.df$IDwaiting)){
    cat('IDwaiting column not found, using IDpregclock column\n\n')
    input.df$IDwaiting <- 1
    input.df$IDwaiting[input.df$IDpregclock == 1E9] <- 0
  }
  
  input.df$status <- 'resting'
  input.df$status[input.df$IDwaiting == 1 & input.df$IDlactating == 0 & input.df$IDpregnant == 0] <- 'waiting'
  input.df$status[input.df$IDlactating == 0 & input.df$IDpregnant == 1] <- 'pregnant'
  input.df$status[input.df$IDlactating == 1 & input.df$IDpregnant == 0] <- 'lactating'
  input.df$status[input.df$IDlactating == 1 & input.df$IDpregnant == 1] <- 'preglact'
  
  return(input.df)
}
##

## Set reproductive status of a dataframe of life history output or individual state output from pso file
set.numstatus <- function(input.df = NULL){
  if(is.null(input.df)) return('Please specify an input dataframe with input.df')
  
  if(is.null(input.df$IDwaiting)){
    cat('IDwaiting column not found, using IDpregclock column\n\n')
    input.df$IDwaiting <- 1
    input.df$IDwaiting[input.df$IDpregclock == 1E9] <- 0
  }
  
  ##
  input.df$numstatus <- 0 # Resting
  input.df$numstatus[input.df$IDwaiting == 1 & input.df$IDlactating == 0 & input.df$IDpregnant == 0] <- 1 # Waiting
  input.df$numstatus[input.df$IDlactating == 0 & input.df$IDpregnant == 1] <- 2 # Pregnant
  input.df$numstatus[input.df$IDwaiting == 0 & input.df$IDlactating == 1 & input.df$IDpregnant == 0] <- 3 # Lactating
  input.df$numstatus[input.df$IDwaiting == 1 & input.df$IDlactating == 1 & input.df$IDpregnant == 0] <- 4 # Lactating & Waiting
  input.df$numstatus[input.df$IDlactating == 1 & input.df$IDpregnant == 1] <- 5 # Pregnant & Lactating
  
  return(input.df)
}
##

## Set reproductive status of a dataframe of life history output
set.fieldstatus <- function(input.df = NULL){
  if(is.null(input.df)) return('Please specify an input dataframe with input.df')
  
  if(is.null(input.df$IDwaiting)){
    cat('IDwaiting column not found, using IDpregclock column\n\n')
    input.df$IDwaiting <- 1
    input.df$IDwaiting[input.df$IDpregclock == 1E9] <- 0
  }
  
  input.df$fieldstatus <- 'nonlactating'
  input.df$fieldstatus[input.df$IDlactating == 1] <- 'lactating'
  input.df$fieldstatus[input.df$age_yrs < (1223/365)] <- 'calve'
  
  return(input.df)
}
##

## Set ageclass of a dataframe of life history output
set.ageclass <- function(input.df = NULL){
  if(is.null(input.df)) return('Please specify an input dataframe with input.df')
  
  # Female is initiated as juveniles (age at weaning)
  input.df$ageclass <- 'scenescent'
  input.df$ageclass[input.df$age_yrs < 25] <- 'mature'
  input.df$ageclass[input.df$age_yrs < 15] <- 'maturegrowing'
  input.df$ageclass[input.df$age_yrs < 8] <- 'juvenile'
  input.df$ageclass[input.df$age_yrs < (1223/365)] <- 'calve'
  
  return(input.df)
}
##

## Plot lifehistory output ##
plot.lifehistory <- function(outlist = NULL, xlim = NULL, ylim1 = c(-50,150), ylim2 = c(0,480), ylim3 = c(0,1200), ylim4 = NULL, export = FALSE, pdf.name = 'plot_lifehistory.pdf', calf = FALSE){
  if(!is.list(outlist)) return(cat('Specify outlist as a list'))
  if(export == TRUE) pdf(file = pdf.name, width = 8, height = 12)
  par(mfcol = c(4,1), mar=c(4,5,2,5))
  nms <- names(outlist$out)
  
  
  if(is.null(xlim)) xl <- range(outlist$out$age_yrs) else xl <- xlim
  outname <- strsplit(outlist$rep[5], ":")[[1]][2]
  outname <- paste0(outname, ".out")
  
  if(calf == TRUE) sw <- (which(names(plot.df) == 'age_yrs_calf') - which(names(plot.df) == 'age_yrs')) else sw <- 0
  
  plot.df <- single$out
  
  ## ------- PLOT ENERGETICS ------- ##
  if(is.null(ylim1)) yl1  <- c(0,1.05*max(outlist$out[,which(nms == 'IDingest')+sw])) else yl1 <- ylim1
  plot('', xlim=xl, ylim=yl1, xlab='', ylab='')
  mtext(side=3, cex=1, line=0.5, paste('Energetics plot of', outname))
  mtext(side=1, line=2.5, 'Age (years)')
  mtext(side=2, line=2.5, 'Energetic rate (MJ/day)')
  lines(plot.df$age_yrs, plot.df[,which(nms == 'IDingest')+sw], col='red',lwd=2)
  lines(plot.df$age_yrs, plot.df[,which(nms == 'IDringest')+sw], col='darkgreen',lwd=2)
  lines(plot.df$age_yrs, plot.df[,which(nms == 'IDmingest')+sw], col='purple',lwd=2)
  lines(plot.df$age_yrs, plot.df[,which(nms == 'IDmaint')+sw], col='black',lwd=2)
  lines(plot.df$age_yrs, plot.df[,which(nms == 'IDgrowth')+sw], col='blue',lwd=2)
  lines(plot.df$age_yrs, plot.df[,which(nms == 'IDnetenergy')+sw], col='brown', lty=2,lwd=1)
  abline(h = 0, lwd=0.5)
  if(calf == TRUE){
    legend("topleft", legend = c("Total Assimliation", "Resource Assimilation", "Milk Assimilation", "Maintenance", "Growth costs", "Net energy"), 
           lwd=c(rep(2,5),1), lty=c(rep(1,5),2), col=c('red', 'darkgreen', 'purple', 'black','blue', 'brown'), bty='n')
  } else {
    with(plot.df, lines(age_yrs, IDpregcosts, col='orange', lwd=2))
    with(plot.df, lines(age_yrs, IDlactcosts, col='magenta', lwd=2))
    legend("topleft", legend = c("Total Assimilation", "Resource Assimilation", "Milk Assimilation", "Maintenance", "Growth costs", "Pregnancy costs", "Lactation costs", "Net energy"), 
           lwd=c(rep(2,7),1), lty=c(rep(1,7),2), col=c('red', 'darkgreen', 'purple', 'black','blue','orange','magenta','brown'), bty='n')
  }
  
  
  ## ------- PLOT SURVIVAL/LENGTH ------- ##
  if(is.null(ylim2)) yl2 <- c(0,1.05*max(plot.df[,which(nms == 'length')+sw])) else yl2 <- ylim2
  plot('', xlim=xl, ylim=c(0,1), xlab='', ylab='')
  mtext(side=1, line=2.5, 'Age (years)')
  mtext(side=2, line=2.5, 'Survival and F/W ratio')
  mtext(side=3, cex=1, line=0.5, paste('Survival/Length plot of', outname))
  abline(h = 0.15, lty = 2); abline(h = 0.3, lty = 2)
  lines(plot.df$age_yrs, plot.df[,which(nms == 'IDfatratio')+sw], col='brown', lwd=2)
  lines(plot.df$age_yrs, exp(-plot.df[,which(nms == 'totalmort')+sw]), col='red', lwd=2)
  lines(plot.df$age_yrs, exp(-plot.df[,which(nms == 'IDminloglifespan')+sw]), col='red', lwd=2, lty=2)
  legend("topright", legend = c("Survival", "Survival threshold", "Fat ratio", "Length"), lty=c(1,2,1,1), lwd=2, col=c('red', 'red','brown','black'), bty='n', horiz = FALSE)
  
  par(new = TRUE)
  plot('', xlim=xl, ylim=yl2, xlab = '', ylab='', axes=FALSE)
  mtext(side=4, line=2.5, 'Length (cm)'); axis(side = 4)
  lines(plot.df$age_yrs, plot.df[,which(nms == 'IDlength')+sw], col='black', lwd=2)
  
  
  ## ------- PLOT STRUCTURE/RESERVES ------- ##
  if(is.null(ylim3)) yl3 <- c(0,1.05*max(plot.df[,which(nms == 'weight')+sw])) else yl3 <- ylim3
  plot('', xlim=xl, ylim=yl3, xlab='', ylab='')
  mtext(side=1, line=2.5, 'Age (years)')
  mtext(side=2, line=2.5, 'Mass (kg)')
  mtext(side=3, cex=1, line=0.5, paste('Structure/Reserves plot of', outname))
  lines(plot.df$age_yrs, plot.df[,which(nms == 'IDweight')+sw], col='blue', lwd=2)
  lines(plot.df$age_yrs, plot.df[,which(nms == 'IDbones')+sw], col='red', lwd=2)
  lines(plot.df$age_yrs, plot.df[,which(nms == 'reserves')+sw], col='black', lwd=2)
  if(calf){
    legend("topleft", legend = c("Weight", "Structure", "Reserves"),
           lwd=c(2,2,2), lty=c(1,1,1), col=c('blue','red','black'), bty='n')
  } else {
    if(!is.null(plot.df$IDbones_calf)) lines(plot.df$age_yrs, plot.df$IDbones_calf, col='red', lwd=1, lty=2)
    if(!is.null(plot.df$reserves_calf)) lines(plot.df$age_yrs, plot.df$reserves_calf, col='black', lwd=1, lty=2)
    legend("topleft", legend = c("Weight", "Structure", "Reserves", "Structure (calf)", "Reserves (calf)"),
           lwd=c(2,2,2,1,1), lty=c(1,1,1,2,2), col=c('blue','red','black', 'red', 'black'), bty='n')
  }
  if(!is.null(plot.df$IDagefirstrepro)) legend('topright', legend = c(paste("Age at first repro = ", round(tail(plot.df$IDagefirstrepro / 365,1), 1)), 
                                                                      paste("# of female calves = ", tail(plot.df$IDreprofemale,1)),
                                                                      paste("Estimated R0 = ", round(tail(plot.df$IDR0estimated,1),2)),
                                                                      paste("EarlyDeaths = ", tail(plot.df$IDearlydeaths,1))),
                                               lty = 0, bty = 'n')
  
  ## ------- PLOT RESOURCE/PREGNANCY THRESHOLD ------- ##
  if(is.null(ylim4)) yl4 <- range(plot.df[,which(nms == 'resource')]) else yl4 <- ylim4
  yl4pl <- yl4 + c(0,4)
  plot('', xlim=xl, ylim=yl4pl, xlab='', ylab='', axes=FALSE)
  box();axis(1);axis(2,at=pretty(yl4,2), labels=pretty(yl4,2), las=2)
  mtext(side=1, line=2.5, 'Age (years)')
  mtext(side=2, line=2.5, 'Resource density')
  mtext(side=3, cex=1, line=0.5, paste('Resource plot of', outname))
  lines(plot.df$age_yrs, plot.df[,which(nms == 'resource')], col='darkgreen', lwd=1)
  if(!calf){
    par(new = TRUE)
    yl5 <- range(plot.df[,which(nms == 'pregthreshold')])
    yl5 <- yl5 + c(-15,0)
    plot('', xlim=xl, ylim=yl5, xlab = '', ylab='', axes=FALSE)
    mtext(side=4, line=2.5, 'Pregnancy Threshold'); axis(side = 4)
    
    with(plot.df, lines(age_yrs, pregthreshold, lwd=1))
    abline(h = 0, lty = 2)
    legend("topleft", legend = c("Resource", "Pregnancy Threshold"),
           lwd=1, col=c('darkgreen','black'),
           bty='n', horiz = FALSE)
  } else {
    legend("topleft", legend = "Resource",
           lwd=1, col='darkgreen',
           bty='n', horiz = FALSE)
  }
  
if(export == TRUE) dev.off()
par(par.defaults)
}
##