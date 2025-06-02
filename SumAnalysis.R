#### Demographic Transitions and Cycles Summary Analysis for Freeman and Robinson 2025.
#Set working directory to soure file location
#setwd(.......)

#Load Packages
library(ggplot2)
library(BayesFactor)
library(dplyr)
library(tidyverse)
library(cowplot)
library(conover.test)
##Load data
cyc<- read.csv("data/expansion.csv")

###Demographic Transitions===============================================================

###Increments of Increase in carring capacity per demographic transition

##Calculate the median by innovation group
MeansInc <- cyc%>% group_by(Innovation) %>%
  summarize(Avg = median(Kincrease))
MeansInc

#Plot the median and distributions of increase in carring capacity per demographic transtion
pcar6 <- ggplot(cyc, aes(factor(Innovation), (Kincrease)))+
  geom_violin(fill="slategrey")+
  # stat_summary(fun=median, geom="point", size=4, color="black")+
  stat_boxplot(geom ='errorbar')+
  geom_jitter(aes(color=factor(Continent)),size=3.5, width = 0.05)+
  geom_point(data = MeansInc, mapping = aes(x = (Innovation), y = Avg), size=4.5, shape=22)+
  geom_line(data = MeansInc, mapping = aes(x = Innovation, y = Avg), size=1.2)+
  #scale_fill_gradient(low ="#F8766D", high = "#619CFF" ) +
  #geom_line(data = Means, mapping = aes(x = factor(Innovation), y = Avg,group=0), color="black", size=1.2)+
  #scale_fill_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  # facet_wrap( ~ factor(CompID1))+
  #scale_y_continuous(limits=c(-15.5, -4.5))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Internal innovation (1= most internal)", y="K increase per demographic transition", 
       title = "A. K increase per Demographic Transitions vs. Innovation")+
  annotate("text", x =2, y = .75, label = "ANOVA BF=190.97", size = 6)
#annotate("text", x =1.55, y = .7, label = "W=-0.64", size = 6)+
# annotate("text", x =2.55, y = .7, label = "W=1.2", size = 6)
pcar6

##Make sure that innovation is read as a factor variable for ANOVA
cyc$Innovation3 <- factor(cyc$Innovation)

##Run ANOVA using the BayesFactor package
bfexr <- anovaBF(Kincrease ~Innovation3, data = cyc)
summary(bfexr)
#print(bfexr)

###Run Kruskal.Wallis test and Conover
kruskal.test(Kincrease~Innovation, data = cyc)
conover.test(cyc$Kincrease, cyc$Innovation, kw=TRUE, method="bonferroni")


###Cycle Rate analysis

#Calculate median cycle rate
Meanscr <- cyc%>% group_by(Innovation) %>%
  summarize(Avg = median(CycleRate*100))
Meanscr

#Plot median and distributions of cycles

pcr <- ggplot(cyc, aes(factor(Innovation), (CycleRate*100)))+
  geom_violin(fill="slategrey")+
  # stat_summary(fun=median, geom="point", size=4, color="black")+
  stat_boxplot(geom ='errorbar')+
  geom_jitter(aes(color=factor(Continent)),size=3.5, width = 0.05)+
  geom_point(data = Meanscr, mapping = aes(x = (Innovation), y = Avg), size=4.5, shape=22)+
  geom_line(data = Meanscr, mapping = aes(x = Innovation, y = Avg), size=1.2)+
  #scale_fill_gradient(low ="#F8766D", high = "#619CFF" ) +
  #geom_line(data = Means, mapping = aes(x = factor(Innovation), y = Avg,group=0), color="black", size=1.2)+
  #scale_fill_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  # facet_wrap( ~ factor(CompID1))+
  #scale_y_continuous(limits=c(-15.5, -4.5))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Internal innovation (1= most internal)", y="Rate of cycles", title = "B. Rate of Cycles vs. Innovation",
       fill="Continent", color="Continent")+
  annotate("text", x =2, y = .17, label = "ANOVA BF=184.03", size = 6)
#annotate("text", x =1.55, y = .15, label = "W=2.25**", size = 6)+
# annotate("text", x =2.55, y = .15, label = "W=1.86", size = 6)
pcr

##BayesFactor ANOVAAnova 
bfcr <- anovaBF(CycleRate ~Innovation3, data = cyc)
summary(bfcr)


##Run Kruskal.Wallis test and Conover
kruskal.test(CycleRate~Innovation, data = cyc)
conover.test(cyc$CycleRate, cyc$Innovation, kw=TRUE, method="bonferroni")
#conover.test(dCIII$DeltaN, dCIII$CompID3, kw=TRUE, method="by")

##Plot Figure 4 from the main text:

Fig4<-plot_grid(pcar6, pcr, ncol=2, align="hv", axis = "rl")
Fig4

pdf("data/Figs/Figure4.pdf", width=16.55, height=14)
Fig4
dev.off()


##Additional Analysis: Number of Demographic Transitions by Innovation
#Calculate median number of DTs by innovation group

Means5 <- cyc%>% group_by(Innovation) %>%
  summarize(Avg = median(ExpandK))
Means5

#Plot medians and distributions

pcar5 <- ggplot(cyc, aes(factor(Innovation), (ExpandK)))+
  geom_violin(fill="slategrey")+
  # stat_summary(fun=median, geom="point", size=4, color="black")+
  stat_boxplot(geom ='errorbar')+
  geom_jitter(aes(color=factor(Continent)),size=3.5, width = 0.05)+
  geom_point(data = Means5, mapping = aes(x = (Innovation), y = Avg), size=4.5, shape=22)+
  geom_line(data = Means5, mapping = aes(x = Innovation, y = Avg), size=1.2)+
  #scale_fill_gradient(low ="#F8766D", high = "#619CFF" ) +
  #geom_line(data = Means, mapping = aes(x = factor(Innovation), y = Avg,group=0), color="black", size=1.2)+
  #scale_fill_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  # facet_wrap( ~ factor(CompID1))+
  #scale_y_continuous(limits=c(-15.5, -4.5))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Internal innovation (1= most internal)", y="Frequency of demographic transitions", 
       title = "A. Frequency of Demographic Transitions vs. Innovation")+
  annotate("text", x =2, y = 8.75, label = "ANOVA BF=544.36", size = 6)
#annotate("text", x =1.55, y = 7.6, label = "W=-0.28", size = 6)+
#annotate("text", x =2.55, y = 7.6, label = "W=3.06**", size = 6)
pcar5

#BayesFactor ANOVA
bfK <- anovaBF(ExpandK ~Innovation3, data = cyc)
summary(bfK)

###Run kurskal and pairwise con0ver test
kruskal.test(ExpandK~Innovation, data = cyc)
conover.test(cyc$ExpandK, cyc$Innovation, kw=TRUE, method="bonferroni")

###################Materials and Methods Summary Figures=============================

##Figure 40
contexpand <- ggplot(cyc, aes(factor(Continent), (Kincrease)))+
  geom_violin(fill="lightgray")+
  geom_jitter(aes(color=factor(Innovation)), size=3, width = 0.1)+
  stat_summary(fun=median, geom="point", size=4, color="black")+
  stat_boxplot(geom ='errorbar')+
  # geom_point(data = Means3, mapping = aes(x = (Innovation), y = Avg, fill=factor(CaseID)), size=4.5, shape=22)+
  # geom_line(data = Means3, mapping = aes(x = Innovation, y = Avg, color=factor(CaseID)), size=1.2)+
  #scale_fill_gradient(low ="#F8766D", high = "#619CFF" ) +
  #geom_line(data = Means, mapping = aes(x = factor(Innovation), y = Avg,group=0), color="black", size=1.2)+
  #scale_fill_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  # facet_wrap( ~ factor(Innovation))+
  #scale_y_continuous(limits=c(-15.5, -4.5))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"), legend.position = c(.2, .75))+
  labs(x = "Continent", y="K increase per demographic transition", 
       title = "K Increase Per Demographic Transition vs. Innovation", color="Innovation (1=Internal)")+
  annotate("text", x =4, y = .75, label = "ANOVA BF=6.74", size = 6)
contexpand

cyc$Continent2 <- factor(cyc$Continent)

bfK <- anovaBF(Kincrease ~Continent2, data = cyc)
summary(bfK)

pdf("data/Figs/Figure40.pdf", width=14.55, height=12)
contexpand
dev.off()


#Figure 41
contexpand2 <- ggplot(cyc, aes(factor(Continent), (ExpandK)))+
  geom_violin(fill="lightgray")+
  geom_jitter(aes(color=factor(Innovation)), size=3, width = 0.05)+
  stat_summary(fun=median, geom="point", size=4, color="black")+
  stat_boxplot(geom ='errorbar')+
  # geom_point(data = Means3, mapping = aes(x = (Innovation), y = Avg, fill=factor(CaseID)), size=4.5, shape=22)+
  # geom_line(data = Means3, mapping = aes(x = Innovation, y = Avg, color=factor(CaseID)), size=1.2)+
  #scale_fill_gradient(low ="#F8766D", high = "#619CFF" ) +
  #geom_line(data = Means, mapping = aes(x = factor(Innovation), y = Avg,group=0), color="black", size=1.2)+
  #scale_fill_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  # facet_wrap( ~ factor(Innovation))+
  #scale_y_continuous(limits=c(-15.5, -4.5))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"), legend.position = c(.15, .85))+
  labs(x = "Continent", y="Frequency of demographic transitions", 
       title = "Frequency of Demographic Transitions vs. Innovation", color="Innovation (1=Internal)")+
  annotate("text", x =4, y = .75, label = "ANOVA BF=21.17", size = 6)
contexpand2

bfK <- anovaBF(ExpandK ~Continent2, data = cyc)
summary(bfK)

kruskal.test(ExpandK~Continent, data = cyc)
conover.test(cyc$ExpandK, cyc$Continent, kw=TRUE, method="bonferroni")

pdf("data/Figs/Figure41.pdf", width=14.55, height=12)
contexpand2
dev.off()

##Figure 42
Means4 <- cyc%>% group_by(Continent) %>%
  summarize(Avg = median(CycleRate))
Means4

contcycle <- ggplot(cyc, aes(factor(Continent), (CycleRate)))+
  geom_violin(fill="lightgray")+
  geom_jitter(aes(color=factor(Innovation)), size=3, width = 0.05)+
  stat_summary(fun=median, geom="point", size=4, color="black")+
  stat_boxplot(geom ='errorbar')+
  # geom_point(data = Means3, mapping = aes(x = (Innovation), y = Avg, fill=factor(CaseID)), size=4.5, shape=22)+
  # geom_line(data = Means3, mapping = aes(x = Innovation, y = Avg, color=factor(CaseID)), size=1.2)+
  #scale_fill_gradient(low ="#F8766D", high = "#619CFF" ) +
  #geom_line(data = Means, mapping = aes(x = factor(Innovation), y = Avg,group=0), color="black", size=1.2)+
  #scale_fill_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  # facet_wrap( ~ factor(Innovation))+
  #scale_y_continuous(limits=c(-15.5, -4.5))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"), legend.position = c(.3, .15))+
  labs(x = "Continent", y="Rate of cycling", 
       title = "Rate of Cycling vs. Innovation", color="Innovation (1=Internal)")+
  annotate("text", x =4, y = .002, label = "ANOVA BF=.37", size = 6)
contcycle

bfK <- anovaBF(CycleRate ~Continent2, data = cyc)
summary(bfK)


pdf("data/Figs/Figure42.pdf", width=14.55, height=12)
contcycle
dev.off()


###Additional analysis: Time and increase in carrying capacity per demographic trnastion

##Run a regression of time on increase in carrying capacity per demographic transition
fit<-glm(Kincrease~Time, data=cyc)
summary(fit)

#Pull the residuals and combine the residuals with the cyc dataframe
residuals<-residuals(fit)
cyc2<-cbind(cyc,residuals)

##Graph the residuals-Greater than 0 indicates larger increases than expected after controlling for time

resid <- ggplot(cyc, aes(factor(Innovation), (residuals)))+
  geom_violin(fill="lightgray")+
  geom_jitter(aes(color=factor(Innovation)), size=3, width = 0.1)+
  stat_summary(fun=median, geom="point", size=4, color="black")+
  stat_boxplot(geom ='errorbar')+
  # geom_point(data = Means3, mapping = aes(x = (Innovation), y = Avg, fill=factor(CaseID)), size=4.5, shape=22)+
  # geom_line(data = Means3, mapping = aes(x = Innovation, y = Avg, color=factor(CaseID)), size=1.2)+
  #scale_fill_gradient(low ="#F8766D", high = "#619CFF" ) +
  #geom_line(data = Means, mapping = aes(x = factor(Innovation), y = Avg,group=0), color="black", size=1.2)+
  #scale_fill_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  # facet_wrap( ~ factor(Innovation))+
  #scale_y_continuous(limits=c(-15.5, -4.5))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"), legend.position = c(.2, .75))+
  labs(x = "Innovation", y="Residual K increase per demographic transition", 
       title = "Residual increase in K per Demographic Transition vs. Innovation", color="Innovation (1=Internal)")+
  annotate("text", x =2, y = .75, label = "ANOVA BF=10.41", size = 6)
resid

##Make sure that innovation is read as a factor variable for ANOVA
cyc2$Innovation3 <- factor(cyc2$Innovation)

##Run ANOVA using the BayesFactor package
bf <- anovaBF(residuals~Innovation3, data = cyc2)
summary(bf)

pdf("data/Figs/Figure43.pdf", width=12.55, height=10)
resid
dev.off()

##A regression that includes innovation shows the same result.
fit2<-glm(Kincrease~Time+Innovation2, data=cyc)
summary(fit2)


#Figure 44--domesticated animals
Means1 <- cyc%>% group_by(Innovation,DomAnimal) %>%
  summarize(Avg = median(Kincrease))
Means1

pcar2 <- ggplot(cyc, aes(factor(Innovation), (Kincrease)))+
  geom_violin(fill="slategrey")+
  # stat_summary(fun=median, geom="point", size=4, color="black")+
  stat_boxplot(geom ='errorbar')+
  geom_jitter(aes(color=factor(DomAnimal)),size=2, width = 0.1)+
  geom_point(data = Means1, mapping = aes(x = (Innovation), y = Avg, fill=factor(DomAnimal)), size=4.5, shape=22)+
  geom_line(data = Means1, mapping = aes(x = Innovation, y = Avg, color=factor(DomAnimal)), size=1.2)+
  #scale_fill_gradient(low ="#F8766D", high = "#619CFF" ) +
  #geom_line(data = Means, mapping = aes(x = factor(Innovation), y = Avg,group=0), color="black", size=1.2)+
  #scale_fill_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  # facet_wrap( ~ factor(CompID1))+
  #scale_y_continuous(limits=c(-15.5, -4.5))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"), legend.position = c(.25, .85))+
  labs(x = "Internal innovation (1= most internal)", y="K increase per demographic transition", title = "A.Increase in Carrying Capacity vs. Internal Innovation",
       fill="Domesticated Animals", color="Domesticated Animals")
pcar2


Means2 <- cyc%>% group_by(Innovation,DomAnimal) %>%
  summarize(Avg = median(CycleRate))
Means2

pcar3 <- ggplot(cyc, aes(factor(Innovation), (CycleRate)))+
  geom_violin(fill="slategrey")+
  # stat_summary(fun=median, geom="point", size=4, color="black")+
  stat_boxplot(geom ='errorbar')+
  geom_jitter(aes(color=factor(DomAnimal)),size=2, width = 0.1)+
  geom_point(data = Means2, mapping = aes(x = (Innovation), y = Avg, fill=factor(DomAnimal)), size=4.5, shape=22)+
  geom_line(data = Means2, mapping = aes(x = Innovation, y = Avg, color=factor(DomAnimal)), size=1.2)+
  #scale_fill_gradient(low ="#F8766D", high = "#619CFF" ) +
  #geom_line(data = Means, mapping = aes(x = factor(Innovation), y = Avg,group=0), color="black", size=1.2)+
  #scale_fill_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  # facet_wrap( ~ factor(CompID1))+
  #scale_y_continuous(limits=c(-15.5, -4.5))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"), legend.position = c(.85, .85))+
  labs(x = "Internal innovation (1= most internal)", y="Rate of cycles", title = "B. Rate of Cycles vs. Internal Innovation",
       fill="Domesticated Animals", color="Domesticated Animals")
pcar3

Fig44<-plot_grid(pcar2, pcar3, ncol=2, align="hv", axis = "rl")
Fig44

pdf("data/Figs/Figure44.pdf", width=20.55, height=14)
Fig44
dev.off()



