##Summary plots for the example cases in the main paper. Figures.....

#Load Packages
library(ggplot2)
library(BayesFactor)
library(dplyr)
library(tidyverse)
library(cowplot)
library(conover.test)

##Set working directory

#SW Asia (Fertile Crescent, Case ID 21)======================================
##Internal Innovation. 30 year time-steps full sequence of growth
sumpop<- read.csv("data/Percapita2/KDESUmmary.csv") #30 year time-steps
sumswasia<-subset(sumpop, region_id2=="SW Asia (21)" & calBP > 8500)

###Load case study processed radiocarbon KDEs at 30. year time-steps
pc30cp<- read.csv("data/Percapita/LevantPerCap2.csv") #30 year time-steps

#pc30cp<- read.csv("data/Percapita200/LevantPerCap200.csv"). ##200 year time-steps

#Subset to the period of the evolution of food production (Neolithic)
pc30cp2<-subset(pc30cp, calBP<11500 & calBP>8500)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(pc30cp2$MKDE-min(pc30cp2$MKDE))/(max((pc30cp2$MKDE)-min(pc30cp2$MKDE)))
StKDE_hi<-(pc30cp2$hi-min(pc30cp2$MKDE))/(max((pc30cp2$MKDE)-min(pc30cp2$MKDE)))
StKDE_lo<-(pc30cp2$lo-min(pc30cp2$MKDE))/(max((pc30cp2$MKDE)-min(pc30cp2$MKDE)))

##Add the standardized KDE to the Neolithic dataframe

pc30cp3<-cbind(StKDE, StKDE_hi, StKDE_lo, pc30cp2)

##Write food production file
#write.table(pc30cp3, file = "data/Percapita2/LevantPerCap2.csv", sep = ",", col.names=NA)

##Collect number of potential equlibria and plot

df2 <- pc30cp3 %>% select(StKDE, PerCap)

# 1. Exact zeros ---------------------------------------------------------
zero_exact <- df2 %>%
  filter(PerCap == 0) %>%
  select(StKDE, PerCap)

# 2. Interpolated zeros between sign changes -----------------------------
# Take consecutive pairs
interp <- lapply(1:(nrow(df2)-1), function(i) {
  y1 <- df2$PerCap[i]
  y2 <- df2$PerCap[i+1]
  
  # Skip if both sides have same sign (no crossing)
  if (y1 == 0 && y2 == 0) return(NULL)
  if (y1 == 0) return(data.frame(MKDE=df2$StKDE[i], PerCap=0))
  if (y2 == 0) return(data.frame(MKDE=df2$StKDE[i+1], PerCap=0))
  if (sign(y1) == sign(y2)) return(NULL)
  
  # Linear interpolation
  x1 <- df2$StKDE[i]
  x2 <- df2$StKDE[i+1]
  
  x_zero <- x1 - y1 * (x2 - x1) / (y2 - y1)
  
  data.frame(MKDE = x_zero, PerCap = 0)
})

zero_interp <- bind_rows(interp)

# 3. Combine and clean ---------------------------------------------------
zero_points <- bind_rows(zero_exact, zero_interp) %>%
  distinct() %>%
  arrange(MKDE)

anchor <- data.frame(MKDE = 0, PerCap = 0.3)

# build line segments: each segment needs two rows
segments <- lapply(seq_len(nrow(zero_points)), function(i) {
  data.frame(
    MKDE   = c(anchor$MKDE, zero_points$MKDE[i]),
    PerCap = c(anchor$PerCap, zero_points$PerCap[i]),
    seg_id = i   # ggplot grouping
  )
}) %>% bind_rows()

segments2<-data.frame(segments)

###Plots
Cpc2 <- ggplot(pc30cp3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),linewidth=1)+
  theme_bw() +
  scale_x_continuous(breaks=c(0,.25, .5, .75, 1), limits=c(0,1.0))+
  scale_y_continuous(limits=c(-.15,.3), breaks=c(-.1,0,.1,.2,.3))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=22, face = "bold"),
        legend.title = element_text(size = 18),   # bigger title
        legend.text  = element_text(size = 18),   # bigger labels
        legend.key.size = unit(2.0, "lines"),     # bigger symbol boxes
        legend.key.width = unit(1.5, "cm"),       # wider legend keys
        legend.key.height = unit(1.2, "cm"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A. SW Asia KDE Per Capita Growth vs. Density")+
 # geom_vline(xintercept = 0.20)+
  geom_hline(yintercept=0)+
    # New line segments from (0, 0.5) to each zero
    geom_line(data = segments, aes(x = MKDE, y = PerCap, group = seg_id), 
              color = "darkgray", linetype="dashed", linewidth = 1,alpha = 0.9) + 
    geom_point( data=zero_points, aes(x=MKDE, y=PerCap), color="darkgray",size=4)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
Cpc2

###Load crop data from Fuller 2025=========================================
cropdat<-read.csv(file="data/Fuller2025Crops2.csv", header=T)
faunadat<-read.csv(file="data/Fuller2025Fauna.csv", header=T)

#Calculated the median of minimum crops per time period
Meanseq <- cropdat%>% group_by(PhaseTime) %>%
  summarize(Avg = median(PerCropsMin))
Meanseq

###Optional graph of minimum crops over time by site
crops <- ggplot(cropdat,aes(x=(calBP*-1), y=(PerCropsMin))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(), size=3) +
 geom_line(data=Meanseq, aes(x=(PhaseTime), y=Avg), color="red", linewidth=2)+
  geom_point(data=Meanseq, aes(x=(PhaseTime), y=Avg), color="red", shape=22, size=3.5)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(limit=c(11500,8500))+
  # scale_y_continuous(limits=c(-.3,0.38))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="Percent domesticated plants", title = "C. Domesticated Plants vs. Time")+
  geom_smooth(se=FALSE)
  #facet_wrap(~region)
crops

#Plot the median and distributions of minimum crops by time period
peq <- ggplot(cropdat, aes(factor(PhaseTime), (PerCropsMin)))+
  geom_violin(fill="slategrey")+
  # stat_summary(fun=median, geom="point", size=4, color="black")+
  stat_boxplot(geom ='errorbar')+
  geom_jitter(aes(color=factor(Region)),size=3.5, width = 0.05)+
  geom_point(data = Meanseq, mapping = aes(x = factor(PhaseTime), y = Avg), size=4.5, shape=22)+
  geom_line(data = Meanseq, mapping = aes(x = factor(PhaseTime), y = Avg), size=1.2)+
  #scale_fill_gradient(low ="#F8766D", high = "#619CFF" ) +
  #geom_line(data = Means, mapping = aes(x = factor(Innovation), y = Avg,group=0), color="black", size=1.2)+
  #scale_fill_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Internal innovation (1=internal)", y="Percent of generations with positive growth", 
       title = "Generations of Positive Per Capita Growth vs. Innovation")
 # annotate("text", x =2, y = .5, label = "ANOVA BF=82.07", size = 6)
#annotate("text", x =1.55, y = .7, label = "W=-0.64", size = 6)+
# annotate("text", x =2.55, y = .7, label = "W=1.2", size = 6)
peq

##Make sure time period is read as a factor variable for ANOVA
cropdat$PhaseTime3 <- factor(cropdat$PhaseTime)

##Run ANOVA using the BayesFactor package
bfexr <- anovaBF(PerCropsMin ~PhaseTime3, data = cropdat)
summary(bfexr)
#print(bfexr)

###Run Kruskal.Wallis test and Conover
kruskal.test(PerCropsMin~PhaseTime3, data = cropdat)
conover.test(cropdat$PerCropsMin, cropdat$PhaseTime3, kw=TRUE, method="bonferroni")

###Plot the Mean KDE and change in minimum crops per time period on the same graph.
Cpt <- ggplot(pc30cp3,aes(x=(calBP), y=StKDE)) +
  geom_ribbon(aes(ymin = StKDE_lo, ymax = StKDE_hi), fill = "grey70", alpha=0.7) +
  geom_point(aes(color=PeriodID), size=3.5) +
  geom_path(aes(y=StKDE),size=1)+
  geom_line(data=Meanseq, aes(x=(PhaseTime), y=Avg/100), color="red", linewidth=2)+
  geom_point(data=Meanseq, aes(x=(PhaseTime), y=Avg/100), color="red", shape=22, size=3.5)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(limits=c(11500,8500))+
  scale_y_continuous(breaks=c(0, .25, .5, .75, 1),
    #limits = c(0, 1)0,
    name = "Standardized KDE density",
    sec.axis = sec_axis(trans = ~ ., name = "Median percent crops" )) +
  #scale_y_continuous(limits=c(0,1))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=22, face = "bold"),
        legend.title = element_text(size = 18),   # bigger title
        legend.text  = element_text(size = 18),   # bigger labels
        legend.key.size = unit(2.0, "lines"),     # bigger symbol boxes
        legend.key.width = unit(1.5, "cm"),       # wider legend keys
        legend.key.height = unit(1.2, "cm"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "A. SW Asia Standardized Density vs. Time")
#geom_hline(yintercept = 0.20)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
Cpt

###Plot 30-year time steps for SW Asia. Population dynamics plots and crops overlaid

#SW Asia paired plot===============================
FigSWAsia<-plot_grid(Cpt, Cpc2, ncol=1, align="hv", axis = "rl")
FigSWAsia

pdf("data/Figs/ManuscriptExamples/SWAsia.pdf", width=20.55, height=14)
FigSWAsia
dev.off()


###Central Texas====================================
sumct<-subset(sumpop, region_id2=="Central Texas (1)" & calBP > 500)
ct30pc<- read.csv("data/Percapita/CTexPerCap.csv")

#200 year time-steps. 2500 year sequence
#ct30pc<- read.csv("data/Percapita200/CTexPerCap200.csv")

ct30pc2<-subset(ct30pc, calBP<4000 & calBP>500)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(ct30pc2$MKDE-min(ct30pc2$MKDE))/(max((ct30pc2$MKDE)-min(ct30pc2$MKDE)))
#StKDE_hit<-(ct30pc2$hi-min(ct30pc2$MKDE))/(max((ct30pc2$MKDE)-min(ct30pc2$MKDE)))
#StKDE_lot<-(ct30pc2$lo-min(ct30pc2$MKDE))/(max((ct30pc2$MKDE)-min(ct30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe

ct30pc3<-cbind(StKDE, ct30pc2)

##Write food production file
#write.table(ct30pc3, file = "data/Percapita2/CTexPerCap.csv", sep = ",", col.names=NA)

##Collect number of potential equlibria and plot

df2 <- ct30pc3 %>% select(StKDE, PerCap)

# 1. Exact zeros ---------------------------------------------------------
zero_exact <- df2 %>%
  filter(PerCap == 0) %>%
  select(StKDE, PerCap)

# 2. Interpolated zeros between sign changes -----------------------------
# Take consecutive pairs
interp <- lapply(1:(nrow(df2)-1), function(i) {
  y1 <- df2$PerCap[i]
  y2 <- df2$PerCap[i+1]
  
  # Skip if both sides have same sign (no crossing)
  if (y1 == 0 && y2 == 0) return(NULL)
  if (y1 == 0) return(data.frame(MKDE=df2$StKDE[i], PerCap=0))
  if (y2 == 0) return(data.frame(MKDE=df2$StKDE[i+1], PerCap=0))
  if (sign(y1) == sign(y2)) return(NULL)
  
  # Linear interpolation
  x1 <- df2$StKDE[i]
  x2 <- df2$StKDE[i+1]
  
  x_zero <- x1 - y1 * (x2 - x1) / (y2 - y1)
  
  data.frame(MKDE = x_zero, PerCap = 0)
})

zero_interp <- bind_rows(interp)

# 3. Combine and clean ---------------------------------------------------
zero_points <- bind_rows(zero_exact, zero_interp) %>%
  distinct() %>%
  arrange(MKDE)

anchor <- data.frame(MKDE = 0, PerCap = 0.3)

# build line segments: each segment needs two rows
segments <- lapply(seq_len(nrow(zero_points)), function(i) {
  data.frame(
    MKDE   = c(anchor$MKDE, zero_points$MKDE[i]),
    PerCap = c(anchor$PerCap, zero_points$PerCap[i]),
    seg_id = i   # ggplot grouping
  )
}) %>% bind_rows()

segments2<-data.frame(segments)

###Plot KDE vs. per capita growth
pcctex <- ggplot(ct30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_continuous(breaks=c(0,.25, .5, .75, 1), limits=c(0,1.0))+
  scale_y_continuous(limits=c(-.15,.3), breaks=c(-.1,0,.1,.2,.3))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=22, face = "bold"),
        legend.title = element_text(size = 18),   # bigger title
        legend.text  = element_text(size = 18),   # bigger labels
        legend.key.size = unit(2.0, "lines"),     # bigger symbol boxes
        legend.key.width = unit(1.5, "cm"),       # wider legend keys
        legend.key.height = unit(1.2, "cm"))+
  # geom_vline(xintercept = 0.20)+
  geom_hline(yintercept=0)+
  # New line segments from (0, 0.5) to each zero
  geom_line(data = segments, aes(x = MKDE, y = PerCap, group = seg_id), 
            color = "darkgray", linetype="dashed", linewidth = 1,alpha = 0.9) + 
  geom_point( data=zero_points, aes(x=MKDE, y=PerCap), color="darkgray",size=4)+
  labs(x = "Mean KDE (density)", y="KDE per capita growth", title = "B. Central Texas Per Capita Growth vs. Density")+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcctex


##Read in Texas Fire-cracked rock data
###Load Burned rock midden data from Freeman et al. 2024=========================================
cropdat<-read.csv(file="data/SurfaceMeans2.csv", header=T)
cropdat<-subset(cropdat, PhaseTime >500)

###Plot of KDE over time against the Burned Rock Midden Index=================
ctexcpt <- ggplot(ct30pc3,aes(x=(calBP), y=(StKDE))) +
  geom_ribbon(data=sumct, aes(ymin = StKDE_lo, ymax = StKDE_hi), fill = "grey70", alpha=.7) +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  geom_line(data=cropdat, aes(x=(PhaseTime), y=BRMIndex), color="red", linewidth=2)+
  geom_point(data=cropdat, aes(x=(PhaseTime), y=BRMIndex), color="red", shape=22, size=3.5)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  scale_y_continuous(breaks=c(0, .25, .5, .75, 1),
    #limits = c(0, 1), 
    name = "Standardized KDE density",
    sec.axis = sec_axis(trans = ~ ., name = "Burned rock midden index" )) +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=22, face = "bold"),
        legend.title = element_text(size = 18),   # bigger title
        legend.text  = element_text(size = 18),   # bigger labels
        legend.key.size = unit(2.0, "lines"),     # bigger symbol boxes
        legend.key.width = unit(1.5, "cm"),       # wider legend keys
        legend.key.height = unit(1.2, "cm"))+
  labs(x = "Years cal BP", y="Standardized KDE Density", title = "B. Central Texas Standardized Density vs. Time")
  #geom_hline(yintercept=0.25)
ctexcpt

#paired plot===============================
Figctex<-plot_grid(ctexcpt, pcctex, ncol=1, align="hv", axis = "rl")
Figctex

pdf("data/Figs/ManuscriptExamples/Exctex.pdf", width=20.55, height=14)
Figctex
dev.off()

###Sonoran Desert=========================================
sumsd<-subset(sumpop, region_id2=="Sonoran Desert (5)" & calBP > 400)
###Plot mean KDE against the per capita growth rate in the Sonoran desert
son30pc<- read.csv("data/Percapita/SonoranPerCap.csv")

#son30pc<- read.csv("data/Percapita200/SonoranPerCap200.csv")

son30pc2<-subset(son30pc, calBP<4900 & calBP>400)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(son30pc2$MKDE-min(son30pc2$MKDE))/(max((son30pc2$MKDE)-min(son30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
son30pc3<-cbind(StKDE, son30pc2)

##Write food production file
#write.table(son30pc3, file = "data/Percapita2/SonoranPerCap.csv", sep = ",", col.names=NA)

##Collect number of potential equlibria and plot

df2 <- son30pc3 %>% select(StKDE, PerCap)

# 1. Exact zeros ---------------------------------------------------------
zero_exact <- df2 %>%
  filter(PerCap == 0) %>%
  select(StKDE, PerCap)

# 2. Interpolated zeros between sign changes -----------------------------
# Take consecutive pairs
interp <- lapply(1:(nrow(df2)-1), function(i) {
  y1 <- df2$PerCap[i]
  y2 <- df2$PerCap[i+1]
  
  # Skip if both sides have same sign (no crossing)
  if (y1 == 0 && y2 == 0) return(NULL)
  if (y1 == 0) return(data.frame(MKDE=df2$StKDE[i], PerCap=0))
  if (y2 == 0) return(data.frame(MKDE=df2$StKDE[i+1], PerCap=0))
  if (sign(y1) == sign(y2)) return(NULL)
  
  # Linear interpolation
  x1 <- df2$StKDE[i]
  x2 <- df2$StKDE[i+1]
  
  x_zero <- x1 - y1 * (x2 - x1) / (y2 - y1)
  
  data.frame(MKDE = x_zero, PerCap = 0)
})

zero_interp <- bind_rows(interp)

# 3. Combine and clean ---------------------------------------------------
zero_points <- bind_rows(zero_exact, zero_interp) %>%
  distinct() %>%
  arrange(MKDE)

anchor <- data.frame(MKDE = 0, PerCap = 0.3)

# build line segments: each segment needs two rows
segments <- lapply(seq_len(nrow(zero_points)), function(i) {
  data.frame(
    MKDE   = c(anchor$MKDE, zero_points$MKDE[i]),
    PerCap = c(anchor$PerCap, zero_points$PerCap[i]),
    seg_id = i   # ggplot grouping
  )
}) %>% bind_rows()

segments2<-data.frame(segments)


pcSon <- ggplot(son30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_y_continuous(limits=c(-.2,.3), breaks=c(-.1,0,.1,.2,.3))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=22, face = "bold"),
        legend.title = element_text(size = 18),   # bigger title
        legend.text  = element_text(size = 18),   # bigger labels
        legend.key.size = unit(2.0, "lines"),     # bigger symbol boxes
        legend.key.width = unit(1.5, "cm"),       # wider legend keys
        legend.key.height = unit(1.2, "cm"))+
  # geom_vline(xintercept = 0.20)+
  geom_hline(yintercept=0)+
  # New line segments from (0, 0.5) to each zero
  geom_line(data = segments, aes(x = MKDE, y = PerCap, group = seg_id), 
            color = "darkgray", linetype="dashed", linewidth = 1,alpha = 0.9) + 
  geom_point( data=zero_points, aes(x=MKDE, y=PerCap), color="darkgray",size=4)+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "C. Sonoran KDE Per Capita Growth vs. Density")
#geom_vline(xintercept = 0.20)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcSon


#Plot KDE vs time including data from Lesure et al. 2021

#3First, data frame based on table 2---12-16 or greater row variety maize. Also this is based on the supporting data
calBP3<-c(4000, 3000, 2000)
PerHigh<-c(0, .29, .74)

maizedat<-data.frame(cbind(calBP3, PerHigh))

Soncpt <- ggplot(son30pc3,aes(x=(calBP), y=(StKDE))) +
  geom_ribbon(data=sumsd, aes(ymin = StKDE_lo, ymax = StKDE_hi), fill = "grey70", alpha=0.7) +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  geom_line(data=maizedat, aes(x=(calBP3), y=PerHigh), color="red", linewidth=2)+
  geom_point(data=maizedat, aes(x=(calBP3), y=PerHigh), color="red", shape=22, size=3.5)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  scale_y_continuous(breaks=c(0, .25, .5, .75, 1),
    #limits = c(0, 1), 
    name = "Standardized KDE density",
    sec.axis = sec_axis(trans = ~ ., name = "Percent 12-16 row maize" )) +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=22, face = "bold"),
        legend.title = element_text(size = 18),   # bigger title
        legend.text  = element_text(size = 18),   # bigger labels
        legend.key.size = unit(2.0, "lines"),     # bigger symbol boxes
        legend.key.width = unit(1.5, "cm"),       # wider legend keys
        legend.key.height = unit(1.2, "cm"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "C. Sonoran Standardized Density vs. Time")
#  geom_hline(yintercept = 0.20)
#geom_vline(xintercept = 2550)+
#geom_vline(xintercept = 2250)+
#geom_vline(xintercept = 830)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
Soncpt

Figsonora<-plot_grid(Soncpt, pcSon, ncol=1, align="hv", axis = "rl")
Figsonora

pdf("data/Figs/ManuscriptExamples/Sonora.pdf", width=17.55, height=15)
Figsonora
dev.off()



####Southern Britain. Attempt to recreate a domesticated plant KDE index from
#the data provided by Bevan et al. 2017.

##Load
box<- read.csv("data/Bevan2017dates.csv")
box2<- subset(box, Latitude>51 & Latitude<55 & Longitude>-5 & Longitude<2)

#write.table(box2, file = "data/SouthBritdates.csv", sep = ",", col.names=NA)
#box2<-read.csv(file="data/SouthBritdates.csv", header=T)

box2 <- box2[!is.na(box2$SiteName), ]

##3Run KDEs for domesticated species===============
boxdomest<-box2 %>%
  filter(Material %in% c("nutshell", "grain", "fruit", "seed"))

allplant <- calibrate(x = boxdomest$CRA,  errors = boxdomest$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = boxdomest$SiteID, ages = boxdomest$CRA, h = 100)
#write.table(boxbins, file = "data/NorFranceboxbins.csv", sep = ",", col.names=NA)

####Run SPD
plants <- spd(allplant, bins=boxbins, runm=200, timeRange=c(7000,4200))
plot(plants, runm=200, xlim=c(7000,4200), type="simple")


##KDE
####make KDEs
US.randates = sampleDates(allplant, bins=boxbins, nsim=1000,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(7000,4200),bw=25, normalised = FALSE)
plot(D.ckde,type='multiline')


##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)
#I then convert NAs at the end due to KDE bandwidth to zeros to easily cbind the data frames below
Check2 <- replace(Check, is.na(Check), 0)
#Calculate the mean KDE of each time step of the 200 simulations
MKDE<-rowMeans(Check2)
#calculate the 5th and 95th percentile of each time step of the 200 KDE simulations
lo<-apply( Check2, #select columns
           1, # row-wise calcs
           quantile, probs=0.05) # give `quantile`
hi <-apply( Check2, #select columns
            1, # row-wise calcs
            quantile, probs=0.95) # give `quantile`

calBP<-plants$grid$calBP
PrDens<-plants$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<- subset(dd, MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/SouthBritAllPlants.csv", sep = ",", col.names=NA)

###3Calculate KDE for grain plants
boxdomest<-box2 %>%
  filter(Material %in% c("nutshell", "grain", "fruit", "seed"))

#boxgrain<-box2 %>%  filter(Material =="grain")

boxgrain<-box2 %>% 
  filter(Material %in% c("grain", "seed"))

grains <- calibrate(x = boxgrain$CRA,  errors = boxgrain$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = boxgrain$SiteID, ages = boxgrain$CRA, h = 100)
#write.table(boxbins, file = "data/NorFranceboxbins.csv", sep = ",", col.names=NA)

####Run SPD
grainplants <- spd(grains, bins=boxbins, runm=200, timeRange=c(7000,4200))
plot(grainplants, runm=200, xlim=c(7000,4200), type="simple")


##KDE
####make KDEs
US.randates = sampleDates(grains, bins=boxbins, nsim=1000,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(7000,4200),bw=25, normalised = FALSE)
plot(D.ckde,type='multiline')


##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)
#I then convert NAs at the end due to KDE bandwidth to zeros to easily cbind the data frames below
Check2 <- replace(Check, is.na(Check), 0)
#Calculate the mean KDE of each time step of the 200 simulations
MKDE<-rowMeans(Check2)
#calculate the 5th and 95th percentile of each time step of the 200 KDE simulations
lo<-apply( Check2, #select columns
           1, # row-wise calcs
           quantile, probs=0.05) # give `quantile`
hi <-apply( Check2, #select columns
            1, # row-wise calcs
            quantile, probs=0.95) # give `quantile`

calBP<-grainplants$grid$calBP
PrDens<-grainplants$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
dd2<-cbind(calBP,PrDens, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
#dd2<- subset(dd, MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/SouthBritGrains.csv", sep = ",", col.names=NA)

##Read in KDES for all plants and grains

plant1<-read.csv("data/KDEs/SouthBritAllPlants.csv", sep = ",")
grain1<- read.csv("data/KDEs/SouthBritGrains.csv", sep = ",")

grain1<-subset(grain1, calBP > 4499)

index<-grain1$MKDE/(grain1$MKDE+plant1$MKDE)

plant1<-cbind(index, plant1)


plot(grain1$calBP~plant1$index)


###Plot mean KDE against the per capita growth rate in the North
sumbr<-subset(sumpop, region_id2=="Southern Britain (27)" & calBP > 4500)

gr30pc<- read.csv("data/Percapita/SouthBritPerCap.csv")
#gr30pc<- read.csv("data/Percapita200/SouthBritPerCap200.csv")

gr30pc2<-subset(gr30pc, calBP<7000 & calBP>4500)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(gr30pc2$MKDE-min(gr30pc2$MKDE))/(max((gr30pc2$MKDE)-min(gr30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
gr30pc3<-cbind(StKDE,gr30pc2)

##Write food production file
#write.table(gr30pc3, file = "data/Percapita2/SouthBritPerCap.csv", sep = ",", col.names=NA)


Grcpt <- ggplot(gr30pc3,aes(x=(calBP), y=(StKDE))) +
  geom_ribbon(data=sumbr, aes(ymin = StKDE_lo, ymax = StKDE_hi), fill = "grey70", alpha=0.7) +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  geom_line(data=plant1, aes(x=(calBP), y=index), color="red", linewidth=2)+
  #geom_point(data=maizedat, aes(x=(calBP), y=PerHigh), color="red", shape=22, size=3.5)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  scale_y_continuous(breaks=c(0, .25, .5, .75, 1),
    #limits = c(0, 1), 
    name = "Standardized KDE density",
    sec.axis = sec_axis(trans = ~ ., name = "Grain and seed index" )) +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=22, face = "bold"),
        legend.title = element_text(size = 18),   # bigger title
        legend.text  = element_text(size = 18),   # bigger labels
        legend.key.size = unit(2.0, "lines"),     # bigger symbol boxes
        legend.key.width = unit(1.5, "cm"),       # wider legend keys
        legend.key.height = unit(1.2, "cm"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "D. S. Britain Standardized Density vs. Time")
  #geom_hline(yintercept = 0.20)
#geom_vline(xintercept = 2550)+
#geom_vline(xintercept = 2250)+
#geom_vline(xintercept = 830)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
Grcpt

####Graph density vs per capita growth'

##Collect number of potential equlibria and plot

df2 <- gr30pc3 %>% select(StKDE, PerCap)

# 1. Exact zeros ---------------------------------------------------------
zero_exact <- df2 %>%
  filter(PerCap == 0) %>%
  select(StKDE, PerCap)

# 2. Interpolated zeros between sign changes -----------------------------
# Take consecutive pairs
interp <- lapply(1:(nrow(df2)-1), function(i) {
  y1 <- df2$PerCap[i]
  y2 <- df2$PerCap[i+1]
  
  # Skip if both sides have same sign (no crossing)
  if (y1 == 0 && y2 == 0) return(NULL)
  if (y1 == 0) return(data.frame(MKDE=df2$StKDE[i], PerCap=0))
  if (y2 == 0) return(data.frame(MKDE=df2$StKDE[i+1], PerCap=0))
  if (sign(y1) == sign(y2)) return(NULL)
  
  # Linear interpolation
  x1 <- df2$StKDE[i]
  x2 <- df2$StKDE[i+1]
  
  x_zero <- x1 - y1 * (x2 - x1) / (y2 - y1)
  
  data.frame(MKDE = x_zero, PerCap = 0)
})

zero_interp <- bind_rows(interp)

# 3. Combine and clean ---------------------------------------------------
zero_points <- bind_rows(zero_exact, zero_interp) %>%
  distinct() %>%
  arrange(MKDE)

anchor <- data.frame(MKDE = 0, PerCap = 0.3)

# build line segments: each segment needs two rows
segments <- lapply(seq_len(nrow(zero_points)), function(i) {
  data.frame(
    MKDE   = c(anchor$MKDE, zero_points$MKDE[i]),
    PerCap = c(anchor$PerCap, zero_points$PerCap[i]),
    seg_id = i   # ggplot grouping
  )
}) %>% bind_rows()

segments2<-data.frame(segments)

pcGr <- ggplot(gr30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_y_continuous(limits=c(-.1,.3), breaks=c(-.1,0,.1,.2,.3))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=22, face = "bold"),
        legend.title = element_text(size = 18),   # bigger title
        legend.text  = element_text(size = 18),   # bigger labels
        legend.key.size = unit(2.0, "lines"),     # bigger symbol boxes
        legend.key.width = unit(1.5, "cm"),       # wider legend keys
        legend.key.height = unit(1.2, "cm"))+
  # geom_vline(xintercept = 0.20)+
  geom_hline(yintercept=0)+
  # New line segments from (0, 0.5) to each zero
  geom_line(data = segments, aes(x = MKDE, y = PerCap, group = seg_id), 
            color = "darkgray", linetype="dashed", linewidth = 1,alpha = 0.9) + 
  geom_point( data=zero_points, aes(x=MKDE, y=PerCap), color="darkgray",size=4)+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "D. S. Britain KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept = 0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcGr



#paired plots
Figsbrit<-plot_grid(Grcpt, pcGr, ncol=1, align="hv", axis = "rl")
Figsbrit

pdf("data/Figs/ManuscriptExamples/ExSBritain.pdf", width=20.55, height=14)
Figsbrit
dev.off()


#Figure 4. Examples of change over time.

#paired plots
Fig4<-plot_grid(Cpt, ctexcpt, Soncpt, Grcpt, ncol=1, align="hv", axis = "rl")
Fig4

pdf("data/Figs/ManuscriptExamples/TimePlot2.pdf", width=16.55, height=20.55)
Fig4
dev.off()

#Figure 5 Examples of density per capita growth plots.

Fig5<-plot_grid(Cpc2, pcctex, pcSon, pcGr, ncol=1, align="hv", axis = "rl")
Fig5

pdf("data/Figs/ManuscriptExamples/PerCapPlot.pdf", width=16.55, height=20.55)
Fig5
dev.off()

###Supporting information AMDUR Mid-continent US Figure vs. KDEs b
#3Figure S1


###First, run the pieswise linear models in ADMUR. Second,
#Save the best fitting model. Thirs, compare with our methodology. The 
#ADMUR package is very computationally time consuming. We llok at two cases here.

library(ADMUR)
library(ggplot2)
library(dplyr)
library(sf)
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)

#####Mid continent of North America, Case ID 3==============================================
##===================================================
#Load radiocarbon ages from P3k data set.
SPD<-read.csv(file="data/RawP3Kc14.csv", header=T)
#subset the data to the specified region
box<- subset(SPD, Latitude>35 & Latitude<41 & Longitude>-93 & Longitude< -87.5)
#boxy<- subset(SPD, Latitude>35 & Latitude<41 & Longitude>-87.51 & Longitude< -80.5)

#remove NA's from the siteID column
box <- box[!is.na(box$SiteID), ]



### 1. Input raw radiocarbon data for ADMUR
age   <- box$Age     # uncalibrated BP
sd <- box$Error   # 1σ errors
site<-box$SiteID #sites
#phase<-box$SiteID #phases
data<-data.frame(age, sd, site)

#Calabrate and create SPD array
CalArray <- makeCalArray(intcal20, calrange = c(5000, 500))
PD <- phaseCalibrator(data, CalArray, remove.external = TRUE)
SPD <- summedPhaseCalibrator( data=data, calcurve=intcal20, calrange=c(5000, 500) )
plotPD(SPD)

set.seed(12345)
CPLparsToHinges(pars=runif(11), years=5000:500)

library(DEoptimR)
best <- JDEoptim(lower = rep(0,5), 
                 upper = rep(1,5), 
                 fn = objectiveFunction, 
                 PDarray = PD, 
                 type = 'CPL', 
                 NP = 100,
                 trace = TRUE)

CPL <- CPLparsToHinges(pars=best$par, years=5000:500)
SPD <- summedPhaseCalibrator( data=data, calcurve=intcal20, calrange=c(5000,500))
plotPD(SPD)
lines(CPL$year, CPL$pdf, lwd=2, col='firebrick')
legend(x=3000, y=max(CPL$pdf), cex=0.7, lwd=2, col='firebrick', bty='n', legend='best fitted CPL')
text(x=CPL$year, y=CPL$pdf, pos=3, labels=c('H1','H2','H3','H4'))

##Write food production file
write.table(SPD, file = "data/Percapita/CPL_SPDMidCont.csv", sep = ",", col.names=NA)
write.table(CPL, file = "data/Percapita/CPL_FitMidCont.csv", sep = ",", col.names=NA)


###load data for mid continent 
sumpop<- read.csv("data/Percapita2/KDESUmmary.csv") #30 year time-steps
summidus<-subset(sumpop, region_id2=="Mid-Continent US (3)" & calBP < 5001 & calBP >499)

###Load case study processed radiocarbon KDEs at 30. year time-steps
us30pc<- read.csv("data/Percapita/USDomestPerCap.csv")
#Initial Domestication: cir. 4200
us30pc2<-subset(us30pc, calBP < 5001 &  calBP>499)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(us30pc2$MKDE-min(us30pc2$MKDE))/(max((us30pc2$MKDE)-min(us30pc2$MKDE)))
StKDE_hi<-(us30pc2$hi-min(us30pc2$MKDE))/(max((us30pc2$MKDE)-min(us30pc2$MKDE)))
StKDE_lo<-(us30pc2$lo-min(us30pc2$MKDE))/(max((us30pc2$MKDE)-min(us30pc2$MKDE)))

##Add the standardized KDE to the Neolithic dataframe

us30pc3<-cbind(StKDE, StKDE_hi, StKDE_lo, us30pc2)


midspd<-read.csv("data/Percapita/CPL_SPDMidCont.csv")
midcpl<-read.csv(file = "data/Percapita/CPL_FitMidCont.csv")

mid1 <- ggplot(us30pc3,aes(x=(calBP), y=(StKDE))) +
  geom_ribbon(aes(ymin = (StKDE_lo), ymax = (StKDE_hi)), fill = "grey70", alpha=0.7) +
  geom_point(aes(y=(StKDE), color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=22, face = "bold"),
        legend.title = element_text(size = 18),   # bigger title
        legend.text  = element_text(size = 18),   # bigger labels
        legend.key.size = unit(2.0, "lines"),     # bigger symbol boxes
        legend.key.width = unit(1.5, "cm"),       # wider legend keys
        legend.key.height = unit(1.2, "cm"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "A. Mid-Continent US KDE Density vs. Time")+
#geom_hline(yintercept = 0.20)
  geom_vline(xintercept = 3350)+
geom_vline(xintercept = 3000)+
geom_vline(xintercept = 2250)+
geom_vline(xintercept = 1600)+
geom_vline(xintercept = 1050)+
annotate("text", x =3400, y = 1, label = "DT1", size = 8)+
  annotate("text", x =2900, y = 1, label = "DT2", size = 8)+
  annotate("text", x =2200, y = 1, label = "DT3", size = 8)+
  annotate("text", x =1500, y = 1, label = "DT4", size = 8)+
  annotate("text", x =1000, y = 1, label = "DT5", size = 8)
mid1


mid2 <- ggplot(midspd,aes(x=(calBP), y=(SPD))) +
  geom_point(aes(y=(SPD), color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  geom_point(data=midcpl, aes(x=(calBP), y=spd), color="firebrick", shape=20, size=5.5)+
  geom_line(data=midcpl, aes(x=(calBP), y=spd), color="firebrick", linewidth=1.5)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
 # scale_y_continuous(name = "Standardized KDE density",
  #                   sec.axis = sec_axis(transform= ~ ., name = "Summed probability distribution" )) +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=22, face = "bold"),
        legend.title = element_text(size = 18),   # bigger title
        legend.text  = element_text(size = 18),   # bigger labels
        legend.key.size = unit(2.0, "lines"),     # bigger symbol boxes
        legend.key.width = unit(1.5, "cm"),       # wider legend keys
        legend.key.height = unit(1.2, "cm"))+
  labs(x = "Years cal BP", y="Summed probability distribution", title = "B. Mid-Continent US CPL, SPD vs. Time")+
geom_vline(xintercept = 3350)+
  geom_vline(xintercept = 3000)+
  geom_vline(xintercept = 2250)+
  geom_vline(xintercept = 1600)+
  geom_vline(xintercept = 1050)+
annotate("text", x =3400, y = .0008, label = "DT1", size = 8)+
annotate("text", x =2900, y = .0008, label = "DT2", size = 8)+
  annotate("text", x =2200, y = .0008, label = "DT3", size = 8)+
annotate("text", x =1500, y = .0008, label = "DT4", size = 8)+
annotate("text", x =1000, y = .0008, label = "DT5", size = 8)
mid2


SIPlot1<-plot_grid(mid1, mid2, ncol=1, align="hv", axis = "rl")
SIPlot1

pdf("data/Figs/ManuscriptExamples/SIPlot1.pdf", width=16.55, height=20.55)
SIPlot1
dev.off()


###Sonoran Desert ADMUR Comparison=========================
#Sonoran Desert: ====================================================
SPD<-read.csv(file="data/RawP3Kc14.csv", header=T)
box<- subset(SPD, Latitude>23 & Latitude<35 & Longitude>-115 & Longitude< -108)
#write.table(boxsd, file = "NERDv4_0/SonoranDesert.csv", sep = ",", col.names=NA)
#boxsd<-read.csv(file="data/SonoranDesert.csv", header=T)
###MAP
#remove NA's from the siteID column
box <- box[!is.na(box$SiteID), ]


### 1. Input raw radiocarbon data for ADMUR
age   <- box$Age     # uncalibrated BP
sd <- box$Error   # 1σ errors
site<-box$SiteID #sites
#phase<-box$SiteID #phases
data<-data.frame(age, sd, site)

#Calabrate and create SPD array
CalArray <- makeCalArray(intcal20, calrange = c(4900, 400))
PD <- phaseCalibrator(data, CalArray, remove.external = TRUE)
SPD <- summedPhaseCalibrator( data=data, calcurve=intcal20, calrange=c(4900, 400) )
plotPD(SPD)

set.seed(12345)
CPLparsToHinges(pars=runif(11), years=4900:400)

library(DEoptimR)
best <- JDEoptim(lower = rep(0,5), 
                 upper = rep(1,5), 
                 fn = objectiveFunction, 
                 PDarray = PD, 
                 type = 'CPL', 
                 NP = 100,
                 trace = TRUE)

CPL <- CPLparsToHinges(pars=best$par, years=4900:400)
SPD <- summedPhaseCalibrator( data=data, calcurve=intcal20, calrange=c(4900,400))
plotPD(SPD)
lines(CPL$year, CPL$pdf, lwd=2, col='firebrick')
legend(x=3000, y=max(CPL$pdf), cex=0.7, lwd=2, col='firebrick', bty='n', legend='best fitted CPL')
text(x=CPL$year, y=CPL$pdf, pos=3, labels=c('H1','H2','H3','H4'))


##Write food production file
write.table(SPD, file = "data/Percapita/CPL_SPDSonora.csv", sep = ",", col.names=NA)
write.table(CPL, file = "data/Percapita/CPL_FitSonora.csv", sep = ",", col.names=NA)


###Sonoran desert comparison of methods
sumpop<- read.csv("data/Percapita2/KDESUmmary.csv") #30 year time-steps
summidus<-subset(sumpop, region_id2=="Sonoran Desert (5)" & calBP < 5001 & calBP >499)


midspd<-read.csv("data/Percapita/CPL_SPDSonora.csv")
midcpl<-read.csv(file = "data/Percapita/CPL_FitSonora.csv")

mid1 <- ggplot(summidus,aes(x=(calBP), y=(StKDE))) +
  geom_ribbon(aes(ymin = (StKDE_lo), ymax = (StKDE_hi)), fill = "grey70", alpha=0.7) +
  geom_point(aes(y=(StKDE), color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=22, face = "bold"),
        legend.title = element_text(size = 18),   # bigger title
        legend.text  = element_text(size = 18),   # bigger labels
        legend.key.size = unit(2.0, "lines"),     # bigger symbol boxes
        legend.key.width = unit(1.5, "cm"),       # wider legend keys
        legend.key.height = unit(1.2, "cm"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "A. Sonoran Desert KDE Density vs. Time")+
  #geom_hline(yintercept = 0.20)
  geom_vline(xintercept = 4300)+
  geom_vline(xintercept = 3600)+
  geom_vline(xintercept = 3200)+
  geom_vline(xintercept = 1950)+
  geom_vline(xintercept = 1400)+
  geom_vline(xintercept = 1050)+
  annotate("text", x =4200, y = 1, label = "DT1", size = 8)+
  annotate("text", x =3700, y = 1, label = "DT2", size = 8)+
  annotate("text", x =3200, y = 1, label = "DT3", size = 8)+
  annotate("text", x =1950, y = 1, label = "DT4", size = 8)+
  annotate("text", x =1500, y = 1, label = "DT5", size = 8)+
annotate("text", x =1050, y = 1, label = "DT6", size = 8)
mid1


mid2 <- ggplot(midspd,aes(x=(calBP), y=(SPD))) +
  geom_point(aes(y=(SPD), color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  geom_point(data=midcpl, aes(x=(calBP), y=spd), color="firebrick", shape=20, size=5.5)+
  geom_line(data=midcpl, aes(x=(calBP), y=spd), color="firebrick", linewidth=1.5)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(name = "Standardized KDE density",
  #                   sec.axis = sec_axis(transform= ~ ., name = "Summed probability distribution" )) +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=22, face = "bold"),
        legend.title = element_text(size = 18),   # bigger title
        legend.text  = element_text(size = 18),   # bigger labels
        legend.key.size = unit(2.0, "lines"),     # bigger symbol boxes
        legend.key.width = unit(1.5, "cm"),       # wider legend keys
        legend.key.height = unit(1.2, "cm"))+
  labs(x = "Years cal BP", y="Summed probability distribution", title = "B. Sonoran Desert CPL, SPD vs. Time")+
  geom_vline(xintercept = 4300)+
  geom_vline(xintercept = 3600)+
  geom_vline(xintercept = 3200)+
  geom_vline(xintercept = 1950)+
  geom_vline(xintercept = 1400)+
  geom_vline(xintercept = 1050)+
  annotate("text", x =4200, y = .001, label = "DT1", size = 8)+
  annotate("text", x =3700, y = .001, label = "DT2", size = 8)+
  annotate("text", x =3200, y = .001, label = "DT3", size = 8)+
  annotate("text", x =1950, y = .001, label = "DT4", size = 8)+
  annotate("text", x =1500, y = .001, label = "DT5", size = 8)+
  annotate("text", x =1050, y = .001, label = "DT6", size = 8)
mid2

SIPlot2<-plot_grid(mid1, mid2, ncol=1, align="hv", axis = "rl")
SIPlot2

pdf("data/Figs/ManuscriptExamples/SIPlot2.pdf", width=16.55, height=20.55)
SIPlot2
dev.off()

