## Specify your working directory
getwd()
# setwd("")

## Source functions into separate environment
detach(mmf)
mmf <- new.env()
source('functions.R', local = mmf)
attach(mmf); rm(mmf)

## Import output files:
single <- import.lh.output(folder = '.', 
                           parameter = 'single', 
                           outnms = names.lh, 
                           statenms = names.lh.state)
tail(single$out)

## Plot lifehistory of female
plot.lifehistory(outlist = single, 
                 xlim = NULL, 
                 ylim1 = c(-50,150), 
                 ylim2 = c(0,480), 
                 ylim3 = c(0,1200), 
                 ylim4 = NULL, 
                 export = FALSE, 
                 pdf.name = 'plot_lifehistory.pdf', 
                 calf = FALSE)

source('Figure_1.R')