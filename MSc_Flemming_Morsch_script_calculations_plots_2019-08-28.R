library("plotly")
library(stringr)
library(dplyr)
library(ggplot2)
library(readxl)
library(RColorBrewer)
library(extrafont) # all fonts
library(lubridate)
library(broom)
#library(xlsx)
library(purrr)
library(tidyr)
library(caret)
library(ggthemes)
library(gridExtra)
library(xtable)
library(stargazer)
library(infer)
library(imager)


#####################
############################################################################
# set working directory
setwd("D:/OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/R/publication")



# Import and sort the data for the XRPD plots ###########################################################################################################################

df <- read.csv("D:/OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/soluplus_Data/data-2019-06-21.csv")
#View(df)
df_all <- df
df <- df %>% group_by(theoret_conc_drug, time) %>% summarise( drug_conc=mean(real_conc_drug), 
                                                              agent_conc=mean(theoret_conc_agent), polymer_conc=mean(real_conc_polymer),
                                                              cryst=mean(crystallinity), sd=sd(crystallinity))
df2 <- df[1:15,]
df3 <- df[c(6:10,16:25),]


df_blank <- read.csv("D:/OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/soluplus_Data/data-2019-06-28_blanks.csv")
#View(df)

df_blank2 <- df_blank %>% group_by(theoret_conc_drug, time) %>% summarise( drug_conc=mean(real_conc_drug), 
                                                                           agent_conc=mean(theoret_conc_agent), polymer_conc=mean(real_conc_polymer),
                                                                           cryst=mean(crystallinity), sd=sd(crystallinity))


df_dp <- rbind(df2)


df_20_30_40 <- rbind(df3)
df_20_30_40$label <- c('a', 'a', 'a', 'a', 'a', 'b','b','b','b','b','c','c','c','c','c')

# XRPD raw data

xrpd_sample <- read_excel('D://OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/soluplus_Data/database_soluplus.xlsx', sheet = '2')
xrpd_blank <- read_excel('D://OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/soluplus_Data/database_soluplus_blanks.xlsx', sheet = '2')
xrpd_ind <- read_excel('D://OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/soluplus_Data/database_soluplus.xlsx', sheet = 'xrpd_IND')
solu_stab <- read_excel('D://OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/soluplus_Data/database_soluplus.xlsx', sheet = 'stability')


solu_stab <- solu_stab %>% mutate(g_15_IND_dp_r3_20min_before = g_15_IND_dp_r3_20min_before + 10000, 
                                  g_15_IND_dp_r3_20min_after = g_15_IND_dp_r3_20min_after +8000, 
                                  g_20_IND_dp_r3_20min_before = g_20_IND_dp_r3_20min_before +6000, 
                                  g_20_IND_dp_r3_20min_after = g_20_IND_dp_r3_20min_after +4000,
                                  g_25_IND_dp_r3_20min_before = g_25_IND_dp_r3_20min_before + 2000,
                                  g_20_IND_30_r3_20min_before = g_20_IND_30_r3_20min_before+ 14000, 
                                  g_20_IND_30_r3_20min_after = g_20_IND_30_r3_20min_after + 12000)
names(solu_stab) <- c('position', 'SOL 25I/15G before', 'SOL 25I/15G after', 'SOL 25I/20G before', 'SOL 25I/20G after',
                      'SOL 25I/25G before', 'SOL 25I/25G after','SOL 38I/20G before', 'SOL 38I/20G after')
solu_stab_long <- gather(solu_stab, sample, value, -position)

ggplot(solu_stab_long, aes(x=position, y=value, colour=sample))+
  #geom_line(aes(y=xrpd$g_0_solu_IND_0.25_20min_r1+5000, colour="10% IND 0 min"))+
  geom_line()+
  
  #scale_color_manual("" , values =c( "brown", "#FF8000","#E41A1C", "#00FFFF", "#FF00FF"),
  #                   labels=c(  'SOL 25I/0G','SOL 25I/15G', 'SOL 25I/20G', 'SOL 25I/25G', 'Indomethacin' ))+
  #my_theme2()+
  my_theme2()+
  scale_x_continuous(expression(paste("Position 2", Theta," [°]",sep="")))+
  scale_y_continuous("Intensity", breaks = c())+
  theme(legend.position = 'top')+   # c(0.2,0.85)
  #labs( tag= "D)")
  
  ggsave("2019-08-27_solu_IND_dp_xrpd_stab.png", dpi=500, width=10, height=7)








# sepearte the smaples by height --> better visualisation

xrpd_dp <- cbind(xrpd_sample[,c(1,37,49,60)], xrpd_blank[,2], xrpd_ind[,2]) 
#names(xrpd_dp) <- c("position"     , "g_15_IND_dp_r3_20min"       "g_20_IND_dp_r3_20min"       "g_25_IND_dp_r2_20min"       "a_g_0_solu_IND_0.25_20min_r1"
#                    "IND" )
xrpd_dp <- xrpd_dp %>% mutate(g_15_IND_dp_r3_20min = g_15_IND_dp_r3_20min+3000, g_20_IND_dp_r3_20min = g_20_IND_dp_r3_20min+2000, 
                              g_0_solu_IND_0.25_20min_r1 = g_0_solu_IND_0.25_20min_r1+5000, IND = IND+9000)


xrpd_dp_long <- gather(xrpd_dp, sample, value, -position)

# plot the diffractograms for drug/polymer ratio = constant, glycerol is changed: 15-25% w/w

x_dp <- ggplot(xrpd_dp_long, aes(x=position, y=value, colour=sample))+
  #geom_line(aes(y=xrpd$g_0_solu_IND_0.25_20min_r1+5000, colour="10% IND 0 min"))+
  geom_line()+
  
  scale_color_manual("" , values =c( "brown", "#FF8000","#E41A1C", "#00FFFF", "#FF00FF"),
                     labels=c(  'SOL 25I/0G','SOL 25I/15G', 'SOL 25I/20G', 'SOL 25I/25G', 'Indomethacin' ))+
  #my_theme2()+
  my_theme2()+
  scale_x_continuous("Position [°2Theta]")+
  scale_y_continuous("Intensity", breaks = c())+
  theme(legend.position = c(0.2,0.85))+
  labs( tag= "D)")

ggsave("2019-07-12_solu_IND_dp_xrpd.png", dpi=500, width=10, height=7)


# XRPD quantified - drug/polymer ratio = constant, glycerol is changed: 15-25% w/w

quant_dp <- ggplot(df_dp, aes(x=time, y=cryst, col=factor(agent_conc)))+
  geom_hline(aes(yintercept = 0), linetype='dotted')+
  geom_point(size=4, alpha=.7 )+
  geom_errorbar(aes(ymin=cryst-sd, ymax=cryst+sd), width=0.2)+
  #facet_wrap(~factor(agent_conc))+
  scale_x_continuous("Time [min]",    breaks = c(0,5,10,15,20))+
  scale_y_continuous("Crystallinity [%]", limits = c(-5,105))+
  scale_color_manual("" , values =c( "#FF8000" ,"#E41A1C", "#00FFFF"),
                     labels=c(  'SOL 25I/15G', 'SOL 25I/20G', 'SOL 25I/25G' ))+
  #my_theme2()+
  my_theme2()+
  
  theme(legend.position = c(0.75, 0.75))+
  labs( tag= "B)")

ggsave("2019-06-09_report11_solu_IND_dp.png", dpi=500, width=9, height=7)

# plot the blanks for the appendix
ggplot(df_blank2, aes(x=time, y=cryst, col=factor(theoret_conc_drug)))+
  geom_hline(aes(yintercept = 0), linetype='dotted')+
  geom_point(size=4, alpha=.7 )+
  geom_errorbar(aes(ymin=cryst-sd, ymax=cryst+sd), width=0.002)+
  #facet_wrap(~factor(agent_conc))+
  scale_x_continuous("Time [min]",    breaks = c(0,5,10,15,20))+
  scale_y_continuous("Crystallinity [%]", limits = c(-5,105))+
  scale_color_manual("" , values =c( 'brown', 'purple', 'grey'),
                     labels=c(  'SOL 25I/0G', 'SOL 38I/0G', 'SOL 51I/0G' ))+
  #my_theme2()+
  my_theme2()+
  
  theme(legend.position = c(0.75, 0.75))+
  labs( tag= "B)")

ggsave("2019-06-09_report11_solu_IND_blank.png", dpi=500, width=9, height=7)
######################################################################################################################################################################

# XRPD diffractograms - plot with increasing drug content 20 - 40 %

# sepearte the smaples by height --> better visualisation

xrpd_20_40 <- cbind(xrpd_sample[,c(1,23,13,49, 35,59)], xrpd_blank[,c(4,7,10)], xrpd_ind[,2]) 
names(xrpd_20_40) <- c( "position" , "f_g_20_IND_40_r1_20min",    "e_g_20_IND_30_r1_20min",   "d_g_20_IND_dp_r1_20min" ,"g_g_15_IND_dp_r1_20min", 
                        "h_g_25_IND_dp_r1_20min","a_g_0_solu_IND_0.25_20min_r1", "b_g_0_solu_IND_0.38_20min_r1" ,"c_g_0_solu_IND_0.51_20min_r1" ,"i_IND" )

xrpd_20_40 <- xrpd_20_40 %>% mutate(a_g_0_solu_IND_0.25_20min_r1 = a_g_0_solu_IND_0.25_20min_r1+17000, b_g_0_solu_IND_0.38_20min_r1 = b_g_0_solu_IND_0.38_20min_r1+12000,
                                    c_g_0_solu_IND_0.51_20min_r1 = c_g_0_solu_IND_0.51_20min_r1+7000, d_g_20_IND_dp_r1_20min = d_g_20_IND_dp_r1_20min+6000,
                                    e_g_20_IND_30_r1_20min = e_g_20_IND_30_r1_20min+4000, f_g_20_IND_40_r1_20min = f_g_20_IND_40_r1_20min+2500, 
                                    g_g_15_IND_dp_r1_20min = g_g_15_IND_dp_r1_20min+1500, h_g_25_IND_dp_r1_20min=h_g_25_IND_dp_r1_20min, i_IND = i_IND+22000)


xrpd_20_40_long <- gather(xrpd_20_40, sample, value, -position)

# plot the diffractograms for drug/polymer ratio = constant, glycerol is changed: 15-25% w/w

x_20_40 <- ggplot(xrpd_20_40_long, aes(x=position, y=value, colour=sample))+
  #geom_line(aes(y=xrpd$g_0_solu_IND_0.25_20min_r1+5000, colour="10% IND 0 min"))+
  geom_line()+
  
  scale_color_manual("" ,  labels= c('SOL 25I/0G','SOL 38I/0G', 'SOL 51I/0G',
                                     'SOL 25I/20G', 'SOL 38I/20G', 'SOL 51I/20G', 
                                     'SOL 25I/15G', 'SOL 25I/25G', 'Indomethacin'),
                     values = c("brown", "purple", "grey", "#E41A1C", "#377EB8", "#4DAF4A", "#FF8000", "#00FFFF", "#FF00FF"))+   
  #my_theme2()+
  my_theme2()+
  
  scale_x_continuous("Position [°2Theta]")+
  scale_y_continuous("Intensity", breaks = c(), limits = c(0,50000))+
  theme(legend.position = c(0.15,0.80), legend.background = element_rect(colour = "transparent", fill = "transparent"))+
  labs( tag= "C)")

ggsave("2019-07-12_solu_IND_20-40_xrpd.png", dpi=500, width=10, height=7)




# XRPD quantification - plot with increasing drug content 20 - 40 %  

quant_20_40 <- ggplot(df_20_30_40, aes(x=time, y=cryst, col=label))+
  geom_hline(aes(yintercept = 0), linetype='dotted')+
  geom_point(size=4, alpha=.7)+
  geom_errorbar(aes(ymin=cryst-sd, ymax=cryst+sd), width=0.2)+
  #facet_wrap(~factor(theoret_conc_drug))+
  scale_x_continuous("Time [min]",    breaks = c(0,5,10,15,20))+
  scale_y_continuous("Crystallinity [%]", limits = c(-5,105), breaks = seq(0,100,25))+
  scale_color_manual("" ,  labels= c('SOL 25I/20G', 'SOL 38I/20G', 'SOL 51I/20G', 
                                     'SOL 25I/0G', 'SOL 38I/0G', 'SOL 51I/0G'),
                     values = c("#E41A1C", "#377EB8", "#4DAF4A", "brown", "purple", "grey"), # it should be drug polymer ratio
  )+
  
  #my_theme2()+
  my_theme2()+
  
  theme(legend.position = c(0.75, 0.75))+
  labs( tag= "A)")
#theme(panel.grid.major.y = element_line(colour="grey", size=0.5)) +
ggsave("2019-06-09_report11_solu_IND_20-40.png", dpi=500, width=9, height=7)


# design a 2 x 2 plot with tags of the XRPD results 

xrpd_plot <- grid.arrange(quant_20_40,  quant_dp,  nrow=2)
ggsave("solu_IND_xrpd_plot_2019-07-13.png", xrpd_plot, dpi=300, width = 9, height = 14 )





######################################################################################################################################################################

# Solubility Curves 

sol_solu <- read_excel("D:/OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/soluplus_Data/Solubility_Solu_IND_2019-07-10.xlsx", 
                       sheet = 'solu')
sol_solu_25 <- read_excel("D:/OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/soluplus_Data/Solubility_Solu_IND_2019-07-10.xlsx", 
                          sheet = 'solu_25')
sol_solu_29 <- read_excel("D:/OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/soluplus_Data/Solubility_Solu_IND_2019-07-10.xlsx", 
                          sheet = 'solu_29')
sol_solu_34 <- read_excel("D:/OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/soluplus_Data/Solubility_Solu_IND_2019-07-10.xlsx", 
                          sheet = 'solu_34')


# define sorting function to get the mean of the datapoints 
# not needed right now
#get_sol_points <- function(df){
#  conc_90 <- df[df[2] <= 0.91 & df[2] >= 0.89,]
#  conc_85 <- df[df[2] <= 0.86 & df[2] >= 0.84,]
#  conc_80 <- df[df[2] <= 0.81 & df[2] >= 0.79,]
#  conc_75 <- df[df[2] <= 0.76 & df[2] >= 0.74,]
#  conc_70 <- df[df[2] <= 0.71 & df[2] >= 0.69,]
#  result <- rbind(conc_90, conc_85, conc_80, conc_75, conc_70)
#  return(result)
#}
## Experimental data as mean of n=2

t_solu <- c(160.11,
            156.375,
            153.195,
            148.8,
            139.46,
            133.285
)
f_solu <- c(1,
            0.897545415,
            0.849173711,
            0.799207315,
            0.74999944,
            0.699130279
)
t_solu_25 <- c(160.11,
               157.1,
               151.775,
               141.765,
               130.835,
               116.015
)
f_solu_25 <- c(1,
               0.899829609,
               0.848840468,
               0.801599401,
               0.751129141,
               0.6987999
)
t_solu_29 <- c(160.11,
               155.835,
               153.49,
               148.555,
               137.245,
               125.375
)
f_solu_29 <- c(1,
               0.900025062,
               0.85060144,
               0.800399897,
               0.750565224,
               0.699327006
)
t_solu_34 <- c(160.11, 155.155,148.865,145.005,131.38,128.945)
f_solu_34 <- c(1, 0.89778739,0.847667584,0.798911708,0.749501991,0.699948878)

# Plot the data
# Soluplus - Indomethacin 

sol1 <-ggplot()+
  geom_line( aes(x=sol_solu$Temp, y=sol_solu$Xdrug))+
  geom_line( aes(x=sol_solu$Temp, y=sol_solu$Xupper), linetype='dotted')+
  geom_line( aes(x=sol_solu$Temp, y=sol_solu$Xlower), linetype='dotted')+
  geom_point(aes(x=t_solu, y=f_solu), colour='grey', size=3 )+
  scale_y_reverse(expression(Chi[Drug]), limits=c(1,0), breaks = seq(1,0,-0.1), expand = c(0,0))+
  scale_x_continuous("Temperature [°C]", limits = c(25,165), breaks = seq(25,160,25), expand = c(0,0))+
  my_theme2()+
  labs(tag="A)")
ggsave("solu_IND_sol.png", dpi=500, width=10, height=5)



# Soluplus - Indomethacin - Glycerol

sol2 <- ggplot()+
  geom_line( aes(x=sol_solu_25$Temp, y=sol_solu_25$Xdrug))+
  geom_line( aes(x=sol_solu_25$Temp, y=sol_solu_25$Xupper), linetype='dotted')+
  geom_line( aes(x=sol_solu_25$Temp, y=sol_solu_25$Xlower), linetype='dotted')+
  geom_point(aes(x=t_solu_25, y=f_solu_25), colour= "#E41A1C", size=3 )+
  scale_y_reverse(expression(Chi[Drug]), limits=c(1,0), breaks = seq(1,0,-0.1), expand = c(0,0))+
  scale_x_continuous("Temperature [°C]", limits = c(25,165), breaks = seq(25,160,25), expand = c(0,0))+
  my_theme2()+
  labs( tag= "B)")
ggsave("solu_IND_sol_25.png", dpi=500, width=10, height=5)




sol3 <- ggplot()+
  geom_line( aes(x=sol_solu_29$Temp, y=sol_solu_29$Xdrug))+
  geom_line( aes(x=sol_solu_29$Temp, y=sol_solu_29$Xupper), linetype='dotted')+
  geom_line( aes(x=sol_solu_29$Temp, y=sol_solu_29$Xlower), linetype='dotted')+
  geom_point(aes(x=t_solu_29, y=f_solu_29), colour= "#377EB8", size=3 )+
  scale_y_reverse(expression(Chi[Drug]), limits=c(1,0), breaks = seq(1,0,-0.1), expand = c(0,0))+
  scale_x_continuous("Temperature [°C]", limits = c(25,165), breaks = seq(25,160,25), expand = c(0,0))+
  my_theme2()+
  labs(tag="C)")
ggsave("solu_IND_sol_29.png", dpi=500, width=10, height=5)



sol4 <- ggplot()+
  geom_line( aes(x=sol_solu_34$Temp, y=sol_solu_34$Xdrug))+
  geom_line( aes(x=sol_solu_34$Temp, y=sol_solu_34$Xupper), linetype='dotted')+
  geom_line( aes(x=sol_solu_34$Temp, y=sol_solu_34$Xlower), linetype='dotted')+
  geom_point(aes(x=t_solu_34, y=f_solu_34), colour= "#4DAF4A", size=3 )+
  scale_y_reverse(expression(Chi[Drug]), limits=c(1,0), breaks = seq(1,0,-0.1), expand = c(0,0))+
  scale_x_continuous("Temperature [°C]", limits = c(25,165), breaks = seq(25,160,25), expand = c(0,0))+
  my_theme2()+
  labs(tag="D)")
ggsave("solu_IND_sol_34.png", dpi=500, width=10, height=5)


sol_all <- grid.arrange(sol1, sol2, sol3, sol4, nrow=2)
ggsave("solu_IND_sol_plot_2019-07-10.png", sol_all, dpi=500, width = 10, height = 10 )

################################################################################################################################

# Solubility parameter plot - note all references you use! 

# this plot is intended to clarify why 20% glycerol give the best solubility and therefore the fastest amorphization ! 
# use Hansen solubility parameter for glycerol, soluplus and Indomethacin and the delta thereof (abs value)

sol_param <- data.frame(gly_frac = seq(0,0.5, by=0.05), delta_sol_param = c(3.83,
                                                                            3.1145,
                                                                            2.399,
                                                                            1.6835,
                                                                            0.968,
                                                                            0.2525,
                                                                            0.463,
                                                                            1.1785,
                                                                            1.894,
                                                                            2.6095,
                                                                            3.325))
#sol_param$gly_frac <- sol_param$gly_frac

sol_chi <- rbind(sol_solu[1,], sol_solu_25[1,], sol_solu_29[1,], sol_solu_34[1,]) %>% cbind.data.frame(gly_frac = c(0, 0.253164557, 0.289855072, 0.338983051))

# plot the solubility parameters and the experimental Chi values 

ggplot()+
  geom_point(aes(y=sol_param$delta_sol_param, x=sol_param$gly_frac), size=3)+
  geom_point(aes(y=sol_chi$Xdrug[1]*10, x=sol_chi$gly_frac[1]), colour='grey', size=3)+
  geom_errorbar(aes( x= sol_chi$gly_frac[1], ymin=sol_chi$Xlower[1]*10, ymax=sol_chi$Xupper[1]*10), colour='grey', width=0.01)+
  
  geom_point(aes(y=sol_chi$Xdrug[2]*10, x=sol_chi$gly_frac[2]), colour="#E41A1C", size=3)+
  geom_errorbar(aes( x= sol_chi$gly_frac[2], ymin=sol_chi$Xlower[2]*10, ymax=sol_chi$Xupper[2]*10), colour="#E41A1C", width=0.01)+
  
  geom_point(aes(y=sol_chi$Xdrug[3]*10, x=sol_chi$gly_frac[3]), colour="#377EB8", size=3)+
  geom_errorbar(aes( x= sol_chi$gly_frac[3], ymin=sol_chi$Xlower[3]*10, ymax=sol_chi$Xupper[3]*10), colour="#377EB8", width=0.01)+
  
  geom_point(aes(y=sol_chi$Xdrug[4]*10, x=sol_chi$gly_frac[4]), colour="#4DAF4A", size=3)+
  geom_errorbar(aes( x= sol_chi$gly_frac[4], ymin=sol_chi$Xlower[4]*10, ymax=sol_chi$Xupper[4]*10), colour="#4DAF4A", width=0.01)+
  
  scale_y_continuous(name =  expression(Delta[delta[t]]),breaks = seq(0,4, by=1) ,sec.axis = sec_axis( trans = ~./10 , breaks = seq(0.4,0.1, by=-0.1) ,
                                                                                                       labels = sort(c('0.4','0.3', '0.2', '0.1' )), name = expression(Chi[Drug])))+
  #scale_y_continuous(name =  expression(Delta[delta[t]]),limits = c( 0.2525, 3.83), breaks = seq(3.83, 0.2525, by=-1) 
  #                   sec.axis = sec_axis( trans = ~./10 , breaks = seq(0.4,0.1, by=-0.1) ,
  #                                                                              name = expression(Chi[Drug])))+
  #scale_y_continuous(name =  expression(Delta[delta[t]]) )+
  scale_x_continuous('Glycerol fraction')+
  my_theme2()

ggsave("solu_IND_sol_param_plot_2019-07-13.png", dpi=500, width = 10, height = 7 )



###########################################################################################################################

# FOT studies

fot <- read_excel("D://OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/soluplus_Data/FOT_results_2019-07-11.xlsx")
fot_blank <- fot[,c(1,7,8,9)]
fot_blank$time <- fot_blank$time/60
fot_blank_long <- gather(fot_blank, sample, value, -time)
fot2 <- fot[c(1,3,5,6,4,2)]
names(fot2) <- c('time', 'Compact 1', 'Compact 2', 'Compact 3', 'Compact 4', 'Compact 5')


fot3 <- gather(fot2, sample, value, -time)
fot3$time <- fot3$time / 60

# The data is plotted using a polynomial fit to smooth the curves

ggplot(fot3, aes(x=time, y=value, colour=sample))+
  #geom_line(size=1)+
  geom_smooth(method = 'lm', formula = y ~ poly(x, 6), se = F)+
  #geom_hline(yintercept = 100, linetype='dotted')+
  my_theme2()+
  theme(legend.position = c(0.75,0.15))+
  scale_color_manual("" ,  labels= c( 'SOL 25I/20G', 'SOL 38I/20G',
                                      'SOL 51I/20G', 'SOL 25I/25G','SOL 25I/15G'),
                     values = c("#E41A1C", "#377EB8", "#4DAF4A",'#00FFFF', '#FF8000' ))+ # it should be drug polymer ratio "#FF8000","#E41A1C", "#00FFFF")+
  scale_x_continuous(breaks = seq(0,20, by=1), limits = c(0,20))+
  labs(y='Temperature [°C]', x='Time [min]')


ggsave("FOT_solu_IND_plot_2019-07-12.png",  dpi=500, width = 8, height = 7 )

# FOT blank plot ##########################################################################################
ggplot(fot_blank_long, aes(x=time, y=value, colour=sample))+
  #geom_line(size=1)+
  geom_smooth(method = 'lm', formula = y ~ poly(x, 6), se = F)+
  #geom_hline(yintercept = 100, linetype='dotted')+
  my_theme2()+
  theme(legend.position = c(0.75,0.15))+
  scale_color_manual("" ,  labels= c( 'SOL 25I/0G', 'SOL 38I/0G',
                                      'SOL 51I/0G'),
                     values = c('brown', 'purple', 'grey' ))+ # it should be drug polymer ratio "#FF8000","#E41A1C", "#00FFFF")+
  scale_x_continuous('Time [min]',breaks = seq(0,20, by=1), limits = c(0,20))+
  scale_y_continuous('Temperature [°C]',breaks = seq(25,125,25), limits = c(25,125))
#labs(y='Temperature [°C]', x='Time [min]')


ggsave("FOT_blank_solu_IND_plot_2019-07-12.png",  dpi=500, width = 8, height = 7 )

# test if there is a significant difference between the max temperature of the different samples - THIS DOES NOT MAKE SENSE --> ONLY n=1

fot3_max <- fot3 %>% group_by(sample) %>% summarise(max = max(value))
anova_fot <- aov(max ~ sample, data = fot3_max) %>% tidy()
anova_fot

#############################################################################################################################################

# print the formulation properties in Latex format 

Tg_data <- read_excel("D:/OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/soluplus_Data/database_soluplus.xlsx", 
                      sheet = 'Tg_publication')

print(xtable(Tg_data, type='latex'))



################################################################################################################################################

# Supporting information

calib_comp2 <- read.csv('D://OneDrive - Vrije Universiteit Brussel/GitHub2/xrpd_app2/cal_curve.csv')
#calib_model_final <- lm(calib_comp2, formula = AUC_weight~0+real_conc_drug)
#model_summary <- summary(calib_model_final)
calib_model_final_predict <- lm(calib_comp2, formula = real_conc_drug~0+AUC_weight)
model_summary <- summary(calib_model_final_predict)
#calib_model_final <- lm(calib_comp2, formula = real_conc_drug~0+AUC_weight)
model_auc <- lm(calib_comp2, formula = AUC_weight~0+real_conc_drug)

# predict - also gives Cooks SD

pred <- augment(model_auc, newdata = calib_comp2[,1:27])

plot_cooks <- function(data){
  plot(data$cooksd, pch="*", cex=2, main="Influential observations by Cooks distance", ylab="Cook's distance")  # plot cook's distance
  abline(h = 4*mean(data$cooksd, na.rm=T), col="red") # add cutoff line
  
}

plot(pred$cooksd, pch="*", cex=2, main="Influential observations by Cooks distance", ylab="Cook's distance")  # plot cook's distance
abline(h = 4*mean(pred$cooksd, na.rm=T), col="red") # add cutoff line
#text(x=1:length(pred$cooksd)+1, y=pred$cooksd, labels=ifelse(pred$cooksd>4*mean(pred$cooksd, na.rm=T),names(pred$cooksd),""), col="red")  # add labels


pred2 <- pred

leave <- 1
while(length(leave) != 0){
  print(paste('mean cooks SD: ', mean(pred2$cooksd)))
  leave <- which(pred2$cooksd > 4*mean(pred2$cooksd) )
  
  keep <- which(pred2$cooksd < 4*mean(pred2$cooksd) )
  
  pred2 <- pred2[keep,]
  print(paste('index of sample to leave out: ',leave))
  print(paste('# samples to leave out: ', length(leave)))
  mod <- lm(pred2, formula = AUC_weight~0+real_conc_drug)
  pred2 <- augment(mod, newdata = pred2[,1:27])
}

# find the sample idx --> then predict the rest
idx <- match(pred2$sample, pred$sample)


# design new model based on only the low cooks SD observations --> then predict the other samples! 
model <- lm(pred2, formula = AUC_weight~0+real_conc_drug)
summary(model)
# predict outliers
outliers <- augment(model, newdata = pred[-idx,1:27]) 
outliers$AUC_weight <- outliers$.fitted


pred_final <- rbind.data.frame(pred2, outliers)

model_final <- lm(pred_final, formula = real_conc_drug~0+AUC_weight)
model_final_summary <- summary(model_final)

model_final_summary$coefficients[2] / model_final_summary$coefficients[1] * 3
model_final_summary$coefficients[2] / model_final_summary$coefficients[1] * 10

xtable(model_final_summary)
stargazer(model_final)



# leave one out CV
train(real_conc_drug~0+AUC_weight, method = "lm", data = pred_final[1:27], trControl = trainControl(method = "LOOCV"))

##################################

# plot new cal_curve
ggplot()+
  geom_point( aes(y=pred_final$real_conc_drug, x=pred_final$AUC_weight),colour='#1f77b4', alpha=.6)+
  #geom_point( aes(y=new$AUC_weight, x=new$.fitted),colour='red', alpha=.8, size=2)+
  #geom_errorbar(aes(ymin=mean-sd, 
  #ymax=mean+sd), width=0.2)+
  geom_abline(slope = model_final$coefficients[1], intercept = 0, colour='#ff7f0e', size=1, alpha=0.7)+
  
  scale_y_continuous("Crystalline Indomethacin [%]", breaks = seq(0,50,5), expand = c(0,0), limits = c(0,52) )+
  scale_x_continuous(expression(paste("AUC/weight ", "[intensity * °2",Theta," / mg]")), expand = c(0,0), limits = c(0,500))+
  #annotate("text", x = 100, y = 50, label = paste('R^2 =', round(model_summary$r.squared, digits = 3), sep = " ")+
  #annotate("text", x = 4, y = 25, label = "Some text")+
  my_theme2()

ggsave("xrpd_cal_curve_2019-07-30.png",  dpi=500, width = 8, height = 7 )

###########################
# plot cal curve after 1 oulier removal

ggplot()+
  geom_point( aes(y=calib_comp2$real_conc_drug, x=calib_comp2$AUC_weight),colour='#1f77b4', alpha=.6)+
  #geom_point( aes(y=new$AUC_weight, x=new$.fitted),colour='red', alpha=.8, size=2)+
  #geom_errorbar(aes(ymin=mean-sd, 
  #ymax=mean+sd), width=0.2)+
  geom_abline(slope = calib_model_final_predict$coefficients[1], intercept = 0, colour='#ff7f0e', size=1, alpha=0.7)+
  
  scale_y_continuous("Crystalline Indomethacin [%]", breaks = seq(0,50,5), expand = c(0,0), limits = c(0,52) )+
  scale_x_continuous(expression(paste("AUC/weight ", "[counts*°2",Theta,"/mg]")), expand = c(0,0), limits = c(0,500))+
  #annotate("text", x = 100, y = 50, label = paste('R^2 =', round(model_summary$r.squared, digits = 3), sep = " ")+
  #annotate("text", x = 4, y = 25, label = "Some text")+
  my_theme2()

ggsave("xrpd_cal_curve_2019-07-14.png",  dpi=500, width = 8, height = 7 )

xtable(model_summary)
stargazer(calib_model_final_predict)

###################################################################################################################

# plot for thesis with crystalline versus amorphous IND

xrpd_ind_long <- gather(xrpd_ind, sample, value, -position)

ggplot(xrpd_ind_long, aes(x=position, y=value, colour=sample))+
  #geom_line(aes(y=xrpd$g_0_solu_IND_0.25_20min_r1+5000, colour="10% IND 0 min"))+
  geom_line()+
  
  scale_color_manual("" ,  labels= c('IND crystalline','IND amorphous'), values = c('red', 'green'))+
  #                                   '25% IND + 20% glycerol', '38% IND + 20% glycerol', '51% IND + 20% glycerol', 
  #                                   '25% IND + 15% glycerol', '25% IND + 25% glycerol', 'Indomethacin'),
  #                   values = c("brown", "purple", "grey", "#E41A1C", "#377EB8", "#4DAF4A", "#FF8000", "#00FFFF", "#FF00FF"))+   
  my_theme2()+
  #scale_x_continuous("Position [°2Theta]")+
  scale_x_continuous(expression(paste("Position 2", Theta," [°]",sep="")), breaks = seq(5,35,5))+
  scale_y_continuous("Intensity", breaks = c())+ #, breaks = c(), limits = c(0,50000))+
  theme(legend.position = c(0.15,0.80), legend.background = element_rect(colour = "transparent", fill = "transparent"))+
  #labs( tag= "C)")
  
  ggsave("2019-07-24_IND_xrpd.png", dpi=500, width=10, height=7)


#######################################################################################################################
# tablet weights - wiht glycerol day 0 and day 21 --> check for weight increase

tab_w <- read_excel('D://OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/soluplus_Data/solu_tablet_weight_2019-08-07.xlsx', sheet = '2019-08-20')
tab_w <- na.omit(tab_w)
w_change <- tab_w %>% group_by(compact) %>% summarise(mean = mean(weight_change)*100, sd = sd(weight_change)*100)

ggplot(w_change, aes(y=mean, x=factor(compact), fill=factor(compact)))+
  geom_col()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.5)+
  my_theme2()+
  scale_y_continuous('Weight change after storage [%]', limits = c(0,5), expand = c(0,0))+
  scale_x_discrete('', labels= c( 
    'SOL 25I/20G (n=13) ', 'SOL 38I/20G (n=14)',     
    'SOL 51I/20G (n=12)', 'SOL 25I/25G (n=10)','SOL 25I/15G (n=9)', 
    'SOL 25I/0G (n=13)' , 'SOL 38I/0G (n=14)', 'SOL 51I/0G (n=10)'))+
  scale_fill_manual(labels= c( 
    'Compact 1 (n=13) ', 'Compact 2 (n=14)',
    'Compact 3 (n=12)', 'Compact 4 (n=10)','Compact 5 (n=9)', 
    'Compact 6 (n=13)', 'Compact 7 (n=14)', 'Compact 8 (n=10)'),
    values = c("#E41A1C", 
               "#377EB8", "#4DAF4A",'#00FFFF','#FF8000', 'brown', 'purple', 'grey' ))+
  theme(legend.position = 'None')+
  theme(axis.text.x = element_text(angle = -30))
ggsave("2019-08-20_water_content_solu.png", dpi=500, width=10, height=7)






# we thus see a tendency -- only 3 replicates -- ANOVA is not very usefull here





# ANOVA of compacts to determine significant difference in water absorption
tab_w$weight_change <- tab_w$weight_change*100
anova_analysis <- aov(weight_change ~ compact, data = tab_w) %>% tidy()
anova_analysis

print(xtable(anova_analysis, type='latex'))
# post-hoc testing Bonferoni correction 
k <- 8
alpha <- 0.05
K <- (k*(k-1))/2 
alpha_asterics <- alpha / K
alpha_asterics
# get pairwise t-tests
t_test_results <- pairwise.t.test(tab_w$weight_change, tab_w$compact, p.adjust.method = 'none')
t_test_results
# Check which p-values are significantly different
t_test_results$p.value < alpha_asterics

#print for latex 
print(t_test_results)

print(xtable(tidy(t_test_results), type='latex'))
#stargazer(t_test_results)
# Tidy the result
tidy(t_test_results)
show(tidy(t_test_results))
spread(tidy(t_test_results), group2, p.value))

# use statistical simulation to see if there's actually a difference between the means 
# set the seed for the simulation
set.seed(42)
# Generate bootstrap distribution of means
w_change_boot <- tab_w %>%
  filter(compact == 'Compact 1') %>%
  # Specify the variable of interest
  specify(response = weight_change) %>%  
  # Generate 15000 bootstrap samples
  generate(reps = 10000, type = "bootstrap") %>% 
  # Calculate the mean of each bootstrap sample
  calculate(stat = "mean")

# View its structure
str(w_change_boot)

# Plot the rent_med_ci statistic
ggplot(w_change_boot, aes(x=stat)) +
  # Make it a histogram with a binwidth of 50
  geom_histogram()  

mean(w_change_boot$stat)


# same as function

bootstrap_dist <- function(data, cond){
  set.seed(42)
  # Generate bootstrap distribution of means
  w_change_boot <- data %>%
    filter(compact == cond) %>%
    # Specify the variable of interest
    specify(response = weight_change) %>%  
    # Generate 15000 bootstrap samples
    generate(reps = 10000, type = "bootstrap") %>% 
    # Calculate the mean of each bootstrap sample
    calculate(stat = "mean")
  
  # View its structure
  #str(w_change_boot)
  
  # Plot the rent_med_ci statistic
  print(ggplot(w_change_boot, aes(x=stat)) +
          # Make it a histogram with a binwidth of 50
          geom_histogram())  
  
  # Calculate bootstrap CI as lower and upper quantiles
  conf_info <- w_change_boot %>%
    summarize(
      l = quantile(stat, 0.025),
      u = quantile(stat, 0.975),
      mean = mean(stat)
    ) 
  return(list(w_change_boot, conf_info))
}

# plot Compact 1 vs Compact 6 to see the difference
C1 <- bootstrap_dist(tab_w, 'Compact 1')
C6 <- bootstrap_dist(tab_w, 'Compact 6')
c_both <- rbind(C1[[1]], C6[[1]]) %>% mutate(cond=c(rep('Compact 1',times=10000 ), rep('Compact 6', times=10000)))

ggplot(c_both, aes(x=stat, fill=cond))+
  geom_histogram(bins = 100)

# define a function to do that ! ###########################################

compare_bootstraps <- function(data, cond1, cond2){
  C1 <- bootstrap_dist(data, cond1)
  C2 <- bootstrap_dist(data, cond2)
  c_both <- rbind(C1[[1]], C2[[1]]) %>% mutate(cond=c(rep(cond1,times=10000 ), rep(cond2, times=10000)))
  
  ggplot(c_both, aes(x=stat, fill=cond))+
    geom_histogram(bins = 100)
  
}

# Calculate bootstrap CI as lower and upper quantiles
w_change_boot %>%
  summarize(
    l = quantile(stat, 0.025),
    u = quantile(stat, 0.975)
  ) 

# hypothesis test - done for between the different compacts #######################################################################
# Calculate observed difference in means between all the different compacts 
diff_mean <- tab_w %>%
  filter(compact == 'Compact 1'| compact== 'Compact 2') %>%
  # Group by treatment group
  group_by(compact) %>%
  # Calculate mean change for each group
  summarize(mean_change = mean(weight_change)) %>% 
  # Pull out the value
  pull() %>%
  # Calculate difference
  diff()

# See the result
diff_mean



diff_mean_ht <- tab_w %>%
  filter(compact == 'Compact 1'| compact== 'Compact 2') %>%
  specify(weight_change ~ compact) %>% 
  hypothesize(null = "independence") %>%  
  generate(reps = 10000, type = "permute") %>% 
  calculate(stat = "diff in means", order = c("Compact 1", "Compact 2"))

diff_mean_ht %>%
  # Filter for simulated test statistics greater than observed
  filter(stat >= diff_mean[1]) %>%
  # Calculate p-value
  summarize(p_val = n() / 10000)

# wrap it in a function #############################################################

simulate_effect <- function(data, n_replicates , cond1, cond2) {
  diff_mean = data %>%
    filter(compact == cond1| compact== cond2) %>%
    # Group by treatment group
    group_by(compact) %>%
    # Calculate mean change for each group
    summarize(mean_change = mean(weight_change)) %>% 
    # Pull out the value
    pull() %>%
    # Calculate difference
    diff()
  
  diff_mean_ht = data %>%
    filter(compact == cond1| compact== cond2) %>%
    specify(weight_change ~ compact) %>% 
    hypothesize(null = "independence") %>%  
    generate(reps = n_replicates, type = "permute") %>% 
    calculate(stat = "diff in means", order = c(cond1, cond2))
  
  return(diff_mean_ht %>%
           # Filter for simulated test statistics greater than observed
           filter(stat >= diff_mean) %>%
           # Calculate p-value
           summarize(p_val = n() / n_replicates))}





#####################
############################################################################
# set working directory
setwd("D:/OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/R/publication")

Tg_data <- read_excel('D://OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/python/Tg_model_comparison_import.xlsx', sheet = 'sheet2')

Tg_data %>% select(drug, polymer, agent, tg)

Tg_data2 <- Tg_data %>% group_by(sample,drug, agent ) %>% summarise( polymer = mean(polymer), tg = mean(value_tg), tg_sd = sd(value_tg)) 


#########################################################

# initial experiments with PVP K 12 

# load database

pvp_data <- read_excel('D://OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/Data/XRPD/database_results_thesis_2019-08-07.xlsx')

pvp_data2 <- na.omit(pvp_data)

pvp_data2 <- pvp_data2 %>% filter(theoret_conc_drug == 30, theoret_conc_agent == 13.8) # %>% group_by(time) %>%  summarise(mean = mean(AUC_weight), sd =sd(AUC_weight))

write.csv(pvp_data2, 'pvp_database_thesis_2019-08-09.csv')

pvp_data3 <- pvp_data2 %>% group_by(time) %>%  summarise(mean = mean(AUC_weight), sd =sd(AUC_weight))


ggplot(pvp_data2, aes(x=time, y=AUC_weight, colour = factor(real_conc_drug)))+
  geom_point()+
  geom_line()+
  #geom_hline(aes(yintercept = 0), linetype='dotted')+
  #geom_errorbar(aes(ymin=pvp_data3$mean-pvp_data3$sd, ymax = pvp_data3$mean + pvp_data3$sd))+
  my_theme2()+
  theme(legend.position = 'None')+
  #scale_x_continuous("Time [min]",    breaks = c(0,5,10,15,20))+
  #scale_y_continuous("Crystallinity [%]", limits = c(-5,105), breaks = seq(0,100,25))+
  ggsave("2019-08-10_PVP_30_IND_deviation.png", dpi=500, width=10, height=7)



###############################################################################################

# AUC script 

AUC <- function(x, y, method=c("trapezoid", "step", "spline"), na.rm = FALSE) {
  
  # calculates Area unter the curve
  # example:
  #   AUC( x=c(1,2,3,5), y=c(0,1,1,2))
  #   AUC( x=c(2,3,4,5), y=c(0,1,1,2))
  
  if(na.rm) {
    idx <- na.omit(cbind(x,y))
    x <- x[idx]
    y <- y[idx]
  }
  
  if (length(x) != length(y))
    stop("length x must equal length y")
  
  idx <- order(x)
  x <- x[idx]
  y <- y[idx]
  
  switch( match.arg( arg=method, choices=c("trapezoid","step","spline") )
          , "trapezoid" = { a <- sum((apply( cbind(y[-length(y)], y[-1]), 1, mean))*(x[-1] - x[-length(x)])) }
          , "step" = { a <- sum( y[-length(y)] * (x[-1] - x[-length(x)])) }
          , "spline" = { a <- integrate(splinefun(x, y, method="natural"), lower=min(x), upper=max(x))$value }
  )
  return(a)
}
# XRPD analysis of PVP K 12 data 

pvp <- read_excel('D:/OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/R/publication/pvp_database_thesis_2019-08-09.xlsx', sheet = 1)
pvp_xray <- read_excel('D:/OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/R/publication/pvp_database_thesis_2019-08-09.xlsx', sheet = 2)

# empty result vector
result2 <- c()

base <- baselinefit(data.frame(pvp_xray[1], pvp_xray[35]), gam = 10, maxwdth = 0,  scl.factor = 1.2, tau = 2.5) # if the sample is fully amorphous --> maxwidth = 0 (best result)



# plot the data to see if it is fitted correctly ! 
ggplot(pvp_xray, aes(x=position))+
  geom_line(aes(y=base$y, colour="raw diffraction"))+
  geom_line(aes(y=base$baseline$basisl , colour="baseline"))+
  geom_line(aes(y=base$pmg$fn, colour="fit"))+
  geom_line(aes(y=base$spl$reg,colour="spline"))+
  geom_line(aes(y=base$baseline$peaks, colour='peaks'))+
  #geom_line(aes(y=diff_test$pks[[20]]$fit, colour="peaks"))
  #geom_line(aes(y=corrected, colour="corrected"))+
  labs(y="Intensity", x="Position [°2Theta]")

# if data = good --> add it to the vector
area <- AUC(y= base$baseline$peaks, x= base$x) -  126.3691  # AUC of fully amorphous Indomethacin
area
result2 <- c(result2, area) 


# as for loop
#result <- c()
#
#for (i in 2:ncol(pvp_xray)) {
#  
#  base <- baselinefit(data.frame(pvp_xray[1], pvp_xray[i]), gam = 10, maxwdth = 5,  scl.factor = 1.2, tau = 2.5)
#  area <- AUC(y= base$baseline$peaks, x= base$x) -  126.3691  # AUC of fully amorphous Indomethacin
#  result <- c(result, area)}
#

#pvp$AUC2 <- result


#write.csv(preds, 'pvp_data_thesis_analysed_2019-08-09.csv')


View(pvp)

pvp$AUC_weight <- pvp$AUC2 / pvp$weight_sample_mg



preds <- (augment(model_final, newdata = pvp)) %>% mutate(crystallinity = .fitted/real_conc_drug*100)

View(preds)

# plot the data for 30 % IND with PVP and 20 % glycerol --> different replicates to show the deviation !! 

ggplot(preds, aes(x=time, y=crystallinity, colour = factor(real_conc_drug)))+
  geom_point(size=3)+
  geom_line()+
  geom_hline(aes(yintercept = 0), linetype='dotted')+
  #geom_errorbar(aes(ymin=pvp_data3$mean-pvp_data3$sd, ymax = pvp_data3$mean + pvp_data3$sd))+
  my_theme2()+
  theme(legend.position = 'None')+
  scale_x_continuous("Time [min]",    breaks = c(0,5,10,15,20))+
  scale_y_continuous("Crystallinity [%]", limits = c(-5,105), breaks = seq(0,100,25))+
  ggsave("2019-08-10_PVP_30_IND_deviation.png", dpi=500, width=10, height=7)


# plot the RH humidity per day #######################################


rh <- data.frame(day = c(seq(as.Date('2019-01-28'),as.Date('2019-02-08'),by = 1), seq(as.Date('2019-02-11'),as.Date('2019-02-13'),by = 1)), 
                 RH = c(50,40,40,40,45,45, 45,45,45,45,45,50, 47, 45, 55))



ggplot(rh, aes(x=day, y= RH))+
  geom_line(colour  = 'red', linetype = 'dashed')+
  my_theme2()+
  theme(legend.position = 'None')+
  theme(axis.text.x = element_text(angle = -30))+
  scale_x_date('Day', breaks = '1 day', date_labels = "%d %b")+
  scale_y_continuous('Relative humidity [%]')+
  theme(axis.text.x = element_text(angle = -30))
ggsave("2019-08-10_PRH_plot.png", dpi=500, width=10, height=7)


# Soluplus study - Tg predictions ##########################################################################################################

# load the Tg functions here ###########################


# convert Â°C to Kelvin and vice versa

convert_C_K <- function(kelvin=NULL, celcius=NULL) {
  
  if (is.null(kelvin)) {
    temperature = 273.15 + celcius
    
  }
  
  else if (is.null(celcius)) {
    temperature = kelvin - 273.15}
  return(temperature)
}

# Gordon-Taylot equation ###############################################################################################################################
G_T <- function(Tg1, Tg2, Tg3, w1, w2, w3, rho1, rho2, rho3){
  k1 = (Tg1*rho1)/(Tg2*rho2)
  #print(k1)
  k2 = (Tg1*rho1)/(Tg3*rho3)
  #print(k2)
  Tg = (w1*Tg1 + k1*w2*Tg2+k2*w3*Tg3) / (w1+k1*w2+k2*w3) 
  return(Tg)
}

# load data with properties of the compacts

solu_tg <- read_excel('D://OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/soluplus_Data/database_soluplus.xlsx', sheet = 'Tg_publication')

# set the density values for the 3 compounds --> mention them in the methods section   
Tg1 <- convert_C_K(celcius = as.numeric(solu_tg[10,7]))
Tg2 <- convert_C_K(celcius = as.numeric(solu_tg[11,7]))
Tg3 <- convert_C_K(celcius = as.numeric(solu_tg[9,7]))
w1 <- 0.39 
w2 <- 0.20
w3 <- 0.40 
rho1 <- 1.08
rho2 <- 1.26
rho3 <- 1.31


#View(solu_tg)  
solu_tg <- solu_tg %>% mutate(IND = IND/100, soluplus = soluplus/100, glycerol = glycerol / 100, mg_st = mg_st/100)  
solu_tg <- solu_tg %>% mutate(Tg_pred = convert_C_K(kelvin = G_T(Tg1, Tg2, Tg3, soluplus, glycerol, IND, 1.08, 1.26, 1.31)),
                              Tg_g_s = convert_C_K(kelvin = G_T(Tg1, Tg2, Tg3, 1- `glycerol fraction in polymer`, `glycerol fraction in polymer`, 0, 1.08, 1.26, 1.31))) 


#solu_tg%Tgmix <- convert_C_K(kelvin = G_T(Tg1, Tg2, Tg3, solu_tg$soluplus, solu_tg$glycerol, solu_tg$IND, 1.08, 1.26, 1.14))   



# Plot all PVP results in one plot with horizontal ranges - time needed to get the sample 100% amorphous
# Then make a table with all samples + properties and a column called amorphous yes/no

# get count summary of the different conditions, where samples got amorphous
pvp_data %>% filter(AUC_weight == 0, polymer=='PVP12', theoret_conc_agent !=0)  %>% group_by(theoret_conc_drug, theoret_conc_agent) %>% count()

# calculate statistics for the plot

pvp_data_am <- pvp_data %>% filter(AUC_weight == 0, polymer=='PVP12', theoret_conc_agent !=0)  %>% 
  group_by(theoret_conc_drug, theoret_conc_agent) %>% summarise(time_am=mean(time), time_sd=sd(time), max_time = max(time), min_time = min(time))

pvp_data_am$compact <- c('Compact 1 (n=3)', 'Compact 2 (n=3)','Compact 3 (n=7)','Compact 4 (n=1)','Compact 6 (n=6)')

ggplot(pvp_data_am, aes(x=time_am, y=compact, colour=compact))+
  geom_point(size=3)+
  geom_errorbarh(aes(xmin=min_time, xmax=max_time ),height=0.2, size=1)+
  #my_theme2()+
  my_theme2()+
  scale_y_discrete('')+
  scale_x_continuous('Time [min]', breaks = seq(8,20,2))+
  theme(legend.position = 'None')
ggsave("2019-08-16_PVP_am_plot.png", dpi=500, width=10, height=7)


# plot melting poiint curve and mDSC with Tg of IND 

melt_IND <- read_excel('D://OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/Data/XRPD/database_results_thesis_2019-08-07.xlsx', 
                       col_types = "numeric", sheet = 'IND_dsc_melt')
mDSC_IND <- read_excel('D://OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/Data/XRPD/database_results_thesis_2019-08-07.xlsx',
                       col_types = "numeric", sheet = 'IND_mDSC')


mDSC_IND2 <- na.omit(mDSC_IND[,1:4]) 
melt_IND2 <- na.omit(melt_IND)


melt_IND2_clean <- melt_IND2 %>% gather(cond, value, -temp)
mDSC_IND2_clean <- mDSC_IND2 %>% gather(cond, value, -temp)

# bind data together 
dsc <- rbind(melt_IND2_clean, mDSC_IND2_clean) 


ggplot(mDSC_IND2_clean, aes(x=temp, y=value, col = cond ))+
  geom_line()+
  my_theme2()+
  scale_y_continuous('')


# predic the theoretical PVP Tg values

pvp_tg <- read_excel('D://OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/Data/XRPD/database_results_thesis_2019-08-07.xlsx', sheet = 'pvp_tg')
Tg1 <- convert_C_K(celcius = as.numeric(pvp_tg[8,8]))
Tg2 <- convert_C_K(celcius = as.numeric(pvp_tg[10,8]))
Tg3 <- convert_C_K(celcius = as.numeric(pvp_tg[9,8]))



pvp_tg <- pvp_tg %>% mutate(IND = IND/100, PVP = PVP/100, glycerol = glycerol / 100, mg_st = mg_st/100)  
pvp_tg <- pvp_tg %>% mutate(Tg_pred = convert_C_K(kelvin = G_T(Tg1, Tg2, Tg3, PVP, glycerol, IND, 1.2, 1.26, 1.31)),
                            Tg_g_s = convert_C_K(kelvin = G_T(Tg1, Tg2, Tg3, 1- `glycerol fraction in polymer`, `glycerol fraction in polymer`, 0, 1.2, 1.26, 1.31)))


# plot for pure glycerol Tg - with DMA and DSC

gly_tg_dsc <- read_excel('D://OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/Data/XRPD/database_results_thesis_2019-08-07.xlsx', 
                         sheet = 'dsc_glycerol')
gly_tg_dsc <- gly_tg_dsc[1700:nrow(gly_tg_dsc),]
gly_tg_dma <- read_excel('D://OneDrive - Københavns Universitet/Pharmaceutical_Sciences/Master Thesis/Data/XRPD/database_results_thesis_2019-08-07.xlsx', 
                         sheet = 'dma_glycerol')

gly_tg <- merge(gly_tg_dsc[,c(2,3)], gly_tg_dma[,c(2,6)], by = 'temp', all = T) %>% gather(cond, value, -temp)

gly_tg <- na.omit(gly_tg)

# plot DSC and DMA data in one plot 
tan <- max(gly_tg_dma$Tan_Delta)
temp_tan <- gly_tg_dma[gly_tg_dma$Tan_Delta == tan,2]

tg <- gly_tg_dsc[gly_tg_dsc$Time == 8.83,][1,]


ggplot(gly_tg, aes(x=temp, y=value, col = cond))+
  geom_line()+
  annotate('text', x=-78, y=-.08, label = expression(paste(T[g], ' = -84.1 °C')), colour = 'red' )+
  annotate('text', x = -80, y=0.23, label= expression(paste(T[g], ' = -80.2 °C')), colour = 'green' )+
  geom_point(aes(y=as.numeric(tan), x= as.numeric(temp_tan)), shape = '|', size=6, colour = 'green')+
  geom_point(aes(x=-84.1, y=-0.082), shape = '|', size = 6, colour = 'red')+
  #geom_point(aes(y=as.numeric(tg$Heat_Flow_Normalized), x= as.numeric(tg$temp), shape = '|', size=4, colour = 'red'))+
  
  my_theme2()+
  scale_x_continuous('Temperature [°C]',limits = c(-110,-40), breaks = seq(-120, -40, 10))+
  scale_y_continuous(expression(paste('Heat flow (W/g)')), limits = c( -0.25, 0.25), sec.axis = sec_axis(trans = ~., 
                                                                                                         name = expression(paste('Tan ', Delta))))+
  theme(legend.position = c(0.15,0.80), legend.background = element_rect(colour = "transparent", fill = "transparent"))+
  scale_color_manual('', labels= c('Glycerol (DSC)','Glycerol (DMA)'), values = c('red', 'green', 'black'))
ggsave("2019-08-23_pure_glycerol.png", dpi=500, width=10, height=7)



# plot DsC data for glycerol
ggplot(gly_tg_dsc[,c(2,3)], aes(x=temp, y=Heat_Flow_Normalized))+
  geom_line()+
  annotate('text', x=-70, y=0, label = expression(paste(T[g], ' = -84.13 °C')), colour = 'red' )+
  scale_x_continuous(limits = c(-110,-40), breaks = seq(-120, -40, 10))+
  scale_y_continuous(limits = c( -0.5, 0.5))+
  my_theme2()




# plot DMA data for glycerol
ggplot(gly_tg_dma, aes(x=temp, y=Tan_Delta))+
  geom_line()+
  geom_point(aes(y=as.numeric(tan), x= as.numeric(temp_tan)), shape = 'x', size=4, colour = 'red')+
  annotate('text', x = -60, y=0.2, label= expression(paste(T[g], ' = -79.9 °C')), colour = 'red' )+
  my_theme2()















