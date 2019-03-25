## Uses list 'single' created in GM_DEBPCoD.R
plot.df <- days.to.years(single$out)

## Add fatratio thresholds to data.frame
plot.df$rho <- 0.3 * plot.df$IDweight
plot.df$rhos <- 0.15 * plot.df$IDweight
plot.df$rho_calf <- 0.3 * plot.df$IDweight_calf
plot.df$rhos_calf <- 0.15 * plot.df$IDweight_calf

## Set reproductive status of the female
plot.df <- set.numstatus(plot.df)
plot.df <- set.status(plot.df)

## Select plotting variables for female and calf
idx.mum <- match(c('resource', 'age_yrs', 'reserves', 'IDfatratio', 'rho', 'rhos', 'status', 'numstatus',
                   'IDingest', 'IDringest', 'IDmingest', 'IDmaint', 'IDgrowth', 'IDpregcosts', 'IDlactcosts', 'IDnetenergy'), names(plot.df))
idx.calf <- match(c('resource', 'age_yrs', 'reserves_calf', 'IDfatratio_calf', 'rho_calf', 'rhos_calf', 'status', 'numstatus',
                    'IDingest_calf', 'IDringest_calf', 'IDmingest_calf', 'IDmaint_calf', 'IDgrowth_calf', 'IDpregcosts_calf', 'IDlactcosts_calf', 'IDnetenergy_calf'), names(plot.df))

## Select and transform dataframe for plotting
LH.ener.mum <- plot.df[,idx.mum]
LH.ener.calf <- plot.df[,idx.calf]
names(LH.ener.calf) <- names(LH.ener.mum)
LH.ener.mum$calfmum <- 'mum'; LH.ener.calf$calfmum <- 'calf'
LH.ener <- rbind(LH.ener.mum, LH.ener.calf);
head(LH.ener); tail(LH.ener)

# Get index of negative rates
idx.min <- which(names(LH.ener) %in% c("IDmaint", "IDgrowth", "IDpregcosts", "IDlactcosts"))
# Set set rate to negative values
LH.ener[,idx.min] <- -1*LH.ener[,idx.min]

melt.df <- melt(subset(LH.ener, select = c('age_yrs','calfmum','IDringest','IDmingest','IDmaint','IDgrowth','IDpregcosts','IDlactcosts','IDnetenergy')), 
                id.vars = c('age_yrs','calfmum'))
head(melt.df); levels(melt.df$variable)

## Reorder levels of factor 'calfmum'
melt.df$calfmum <- as.factor(melt.df$calfmum)
melt.df$calfmum <- factor(melt.df$calfmum, levels(melt.df$calfmum)[c(2,1)])

## Reorder levels of factor 'variable'
melt.df$variable <- factor(melt.df$variable,levels(melt.df$variable)[c(1,2,6,5,3,4,7)])

## Create dataframe for panels indicators 
Ener.text <- data.frame(x = c(3.0), y = c(130,50), letter = c('a','b'), calfmum = c('mum', 'calf'))

## Clean some intermediates
rm(LH.ener.mum, LH.ener.calf, idx.mum, idx.calf, idx.min, LH.ener)

## Create the plot
gg.Energetics_calfmum <-
  ggplot() +
  geom_area(data = subset(melt.df, variable != "IDnetenergy"), aes(x = age_yrs, y = value, fill = variable), alpha=.6) +
  geom_line(data = subset(melt.df, variable == "IDnetenergy"), aes(x = age_yrs, y = value, col = variable), size = 0.25) +
  geom_hline(yintercept = 0, size = .2) +
  facet_wrap(~ calfmum, scales = 'free_y', nrow = 2) +
  ## Set the legend layout
  scale_color_manual(values = c('IDnetenergy' = 'black'), name = element_blank(), labels = 'Net energy') +
  scale_fill_manual(values = c('IDringest' = '#4daf4a', 'IDmingest' = '#984ea3', 'IDgrowth' = '#377eb8', 'IDmaint' = '#e41a1c', 'IDpregcosts' = '#ff7f00', 'IDlactcosts' = '#c16b99'),
                    name = element_blank(), breaks = c("IDringest", "IDmingest", "IDgrowth", "IDmaint", "IDpregcosts", "IDlactcosts"),
                    labels = c("Resource assimilation", "Milk assimilation", "Structural growth costs", "Field metabolic rate", "Fetal development", "Lactation costs")) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  geom_text(data = Ener.text,
            mapping = aes(x = x, y = y, label = letter),
            hjust = 'left', vjust = 'top', size = 4) +
  scale_y_continuous(name = "Energetic rate (MJ / day)", sec.axis = sec_axis(~ . * 1, name = "")) +
  ggopts + theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_blank(), legend.box = "vertical") +
  labs(x = "Age of the female (years)", y = "") +
  coord_cartesian(xlim = c(3.5,31))
##

## Save figure as jpeg or pdf
# ggsave(filename = 'Figure_1.pdf',
ggsave(filename = 'Figure_1.jpeg',
       plot = gg.Energetics_calfmum, 
       width = 4,
       height = 6.5, 
       dpi = 600)
##