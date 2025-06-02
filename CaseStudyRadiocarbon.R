
#Script for analyzing the radiocarbon data for 33 archaeological regions associated with
#Freeman and Robinson 2024: Innovations in food production and patterns of long-term population expansion

##Load Packages==================================
library(elevatr)
library(terra)
library(ggplot2)
library(dplyr)
library(zoo)
library(rcarbon)
library(cowplot)
library(tidyverse)
library(sf)
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)

#Case study ID 15 The Central Andes=================================
####Map of central Andes with Elevation===============================================
# Define bounding box for Central Andes
central_andes_bbox <- st_as_sf(st_sfc(
  st_polygon(list(rbind(
    c(-78, -10),
    c(-65, -10),
    c(-65, -25),
    c(-78, -25),
    c(-78, -10)
  ))),
  crs = 4326
))

# Download DEM data using elevatr (z=6 gives lower resolution; z=10 is higher)
dem <- get_elev_raster(locations = central_andes_bbox, z = 6, clip = "locations")

# Save to file if needed
writeRaster(dem, "central_andes_dem.tif", overwrite = TRUE)

# Convert to data frame for ggplot
dem_df <- as.data.frame(dem, xy = TRUE, na.rm = TRUE)
names(dem_df)[3] <- "elevation"


##Add radiocarbon ages
box2<- read.csv("data/SAmericaDomestRCD.csv")

# Plot the DEM
ggplot() +
  geom_raster(data = dem_df, aes(x = x, y = y, fill = elevation)) +
  scale_fill_viridis_c(option = "C", name = "Elevation (m)") +
  coord_fixed() +
  theme_minimal() +
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Country)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)+
  labs(
    title = "Elevation Map of the Central Andes",
    x = "Longitude", y = "Latitude"
  )


###Now extract elevation to radiocarbon ages
# Convert raster to terra format if it's not already

dem_terra <- rast(dem)  # converts RasterLayer to SpatRaster
class(dem_terra)

# Convert box2 to sf and then to SpatVector
box2_sf <- st_as_sf(box2, coords = c("Longitude", "Latitude"), crs = 4326)
box2_vect <- vect(box2_sf)
class(box2_vect)

# Extract elevation values
elev_values <- terra::extract(dem_terra, box2_vect)

# Combine with original box2 data
box2_with_elev <- bind_cols(box2, elevation = elev_values[,2])  # second column is elevation

# Filter
box2_high <- box2_with_elev %>% filter(elevation > 2000)

# Plot the DEM
ggplot() +
  geom_raster(data = dem_df, aes(x = x, y = y, fill = elevation)) +
  scale_fill_viridis_c(option = "C", name = "Elevation (m)") +
  coord_fixed() +
  theme_minimal() +
  geom_point(data=box2_high, aes(Longitude, Latitude, color=factor(Country)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)+
  labs(
    title = "Elevation Map of the Central Andes",
    x = "Longitude", y = "Latitude"
  )

####

###Calibrate the radiocarbon ages==============================
cptcal <- calibrate(x = box2_high$Age,  errors = box2_high$Error, calCurves = "shcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2_high$SiteName, ages = box2_high$Age, h = 100)

#spd.cenand <- spd(cptcal, bins=NA, runm=150, timeRange=c(8200,200))
spd.cenand <- spd(cptcal, bins=boxbins, runm=150, timeRange=c(8200,200))
plot(spd.cenand, runm=150, xlim=c(8200,200), type="simple")

calBP<-spd.cenand$grid$calBP
PrDens<-spd.cenand$grid$PrDens

##KDE=============================================
##KDE
####make KDEs
US.randates = sampleDates(cptcal, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')
#D.ckde$timeRange

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

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
#dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd, file = "data/KDEs/SAmericaKDE50bin2.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/SAmericaKDE50bin2.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........

# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
##Add in the 30 year bin dates
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/SAmericaSumbin2.csv", sep = ",", col.names=NA)

###
###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("data/Sumbin/SAmericaSumbin2.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

##Add ID variable for culture history periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                          breaks=c(240, 500, 750, 1450, 2450, 3450, 8170),
                          labels=c('Inka', 'LIP', 'Tiwanaku', 'Late Formative', 'Early Formative', 'Archaic'))


write.table(pcgrowth, file = "data/Percapita/SAmericaPerCap2.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North

###Plot mean KDE against the per capita growth rate in the North
sa30pc<- read.csv("data/Percapita/SAmericaPerCap2.csv")

sa30pc2<-subset(sa30pc, calBP<6500 & calBP>1500)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(sa30pc2$MKDE-min(sa30pc2$MKDE))/(max((sa30pc2$MKDE)-min(sa30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
sa30pc3<-cbind(StKDE,sa30pc2)

pccenand <- ggplot(sa30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Mean KDE (density)", y="KDE per capita growth", title = "D. Andes KDE Per Capita Growth vs. Density")
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pccenand

cenandcpt <- ggplot(sa30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years BP", y="KDE per capita growth", title = "D. Central Andes Density vs. Time")
#geom_vline(xintercept = 800)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
cenandcpt

##Plot side-by-side per capita and time plots
Figca<-plot_grid(pccenand, cenandcpt, ncol=2, align="hv", axis = "rl")
Figca

pdf("data/Figs/ExCentralAndes1.pdf", width=20.55, height=14)
Figca
dev.off()

###SW Asia (Fertile Crescent) Case ID 21==============================================
#====================================================================================
box<- read.csv("data/nerd.csv")
##subset the data
box<- subset(box, Country %in% c("JO", "IL", "LB", "SY", "IQ", "TR") & Latitude<39 & Longitude >33)
#Due to errors, write the table and manually check the data
#write.table(box, file = "data/Levantdates.csv", sep = ",", col.names=NA)

##Reload the subseted data
#box<- read.csv("data/Levantdates.csv")

# Load world map data as sf object
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define the list of Middle Eastern countries
middle_east_countries <- c(
  "Bahrain", "Cyprus", "Egypt", "Iran", "Iraq", "Israel", "Jordan",
  "Kuwait", "Lebanon", "Oman", "Palestine", "Qatar", "Saudi Arabia",
  "Syria", "Turkey", "United Arab Emirates", "Yemen", "Armenia", "N. Cyprus", "Georgia", 	
  "Azerbaijan"
)

# Filter the data for the Middle East
middle_east <- world[world$name %in% middle_east_countries, ]

# Clean and convert coordinates
box$Longitude <- as.numeric(as.character(box$Longitude))
box$Latitude <- as.numeric(as.character(box$Latitude))

# Remove rows with missing coordinates
box <- box[!is.na(box$Longitude) & !is.na(box$Latitude), ]


# Create a spatial points object from your dataframe using the same CRS
box_sf <- st_as_sf(box, coords = c("Longitude", "Latitude"), crs = st_crs(middle_east))

# Plot
ggplot(data = middle_east) +
  geom_sf(fill = NA, color = "black") +
  geom_sf(data = box_sf, aes(color = factor(Country)), alpha = 0.5, size = 2) +
  theme_minimal() +
  labs(
    title = "Map of the Middle East",
    subtitle = "Countries in the Middle East",
    x = "Longitude",
    y = "Latitude"
  )

###Calibrate the radiocarbon ages
cptcal <- calibrate(x = box$Age,  errors = box$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box$SiteID, ages = box$Age, h = 100)

####Run analysis for component 3 logistic 3500 to 150
spd.fc <- spd(cptcal, bins=boxbins, runm=200, timeRange=c(15000,2000))
plot(spd.fc, runm=150, xlim=c(15000,2000), type="simple")

PrDens<-spd.fc$grid$PrDens
calBP<-spd.fc$grid$calBP

##KDE=============================================
##KDE
####make KDEs
US.randates = sampleDates(cptcal, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(15000,2000),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')
D.ckde$timeRange

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

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
#dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd, file = "data/KDEs/LevantKDE50bin2.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/LevantKDE50bin2.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)


### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
##Add in the 30 year bin dates
calBP<-c(14970,14940,14910,14880,14850,14820,14790,14760,14730,14700,14670,14640,14610,14580,14550,14520,14490,14460,14430,14400,14370,14340,14310,14280,14250,14220,14190,14160,14130,14100,14070,14040,14010,13980,13950,13920,13890,13860,13830,
         13800,13770,13740,13710,13680,13650,13620,13590,13560,13530,13500,13470,13440,13410,13380,13350,13320,13290,13260,13230,13200,13170,13140,13110,13080,13050,13020,12990,12960,12930,12900,12870,12840,12810,12780,12750,12720,12690,12660,12630,12600,12570,12540,12510,12480,12450,12420,12390,12360,12330,12300,12270,
         12240,12210,12180, 12150,12120,12090,12060,12030,12000,11970,11940,11910,11880,11850,11820,11790,11760,11730,11700,11670,11640,11610,11580,11550,11520,11490,11460,11430,11400,11370,11340,11310,11280,11250,11220,11190,11160,11130,11100,11070,11040,11010,10980,10950,10920,10890,10860,10830,10800,10770,10740,10710,10680,10650,10620,10590,10560,10530,10500,10470,10440,10410,10380,10350,10320,10290,10260,10230,10200,10170,10140,
         10110,10080,10050,10020,9990,9960,9930,9900,9870,9840,9810,9780,9750,9720,9690,9660,9630,9600,9570,9540,9510,9480,9450,9420,9390,9360,9330,9300,9270,9240,9210,9180,9150,9120,9090,9060,9030,9000,8970,8940,8910,8880,8850,8820,
         8790,8760,8730,8700,8670,8640,8610,8580,8550,8520,8490,8460,8430,8400,8370,8340,8310,8280,8250,8220,8190,8160,8130,8100,8070,8040,8010,7980,7950,7920,7890,7860,7830,7800,7770,7740,7710,7680,7650,7620,7590,7560,7530,7500,7470,7440,7410,7380,7350,7320,7290,7260,7230,
         7200,7170,7140,7110,7080,7050,7020,6990,6960,6930,6900,6870,6840,6810,6780,6750,6720,6690,6660,6630,6600,6570,6540,6510,6480,6450,6420,6390,6360,
         6330,6300,6270,6240,6210,6180,6150,6120,6090,6060,6030,6000,5970,5940,5910,5880,5850,5820,5790,5760,5730,5700,5670,5640,5610,5580,5550,5520,5490,5460,5430,5400,5370,5340,5310,5280,5250,5220,5190,5160,5130,5100,5070,5040,5010,4980,4950,4920,4890,4860,4830,4800,4770,4740,4710,4680,4650,4620,4590,4560,4530,4500,
         4470,4440,4410,4380,4350,4320,4290,4260,4230,4200,4170,4140,4110,4080,4050,4020,3990,3960,3930,3900,3870,
         3840,3810,3780,3750,3720,3690,3660,3630,3600, 3570,3540,3510,3480,3450,3420,3390,3360,3330,3300,3270,3240,3210,3180,3150,3120,3090,3060,3030,3000,2970,2940,2910,2880,2850,2820,2790,2760,2730,2700,2670,2640,2610,2580,2550,2520,2490,2460,2430,
         2400,2370,2340,2310,2280,2250,2220,2190,2160,2130,2100,2070,2040,2010)

sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/LevantSumbin2.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read.csv("data/Sumbin/LevantSumbin2.csv") 
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

##Add ID variable for culture history periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(2000, 3720,5250, 6450, 8000, 9500, 10100, 10500, 11500, 13590, 14970),
                         labels=c( 'Iron Age','Bronze Age','Calcolithic','PN', 'Late PPNB', 'Middle PPNB', 'Early PPNB', 'PPNA','Late Natufian','Early Natufian'))

write.table(pcgrowth, file = "data/Percapita/LevantPerCap2.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
pc30cp<- read.csv("data/Percapita/LevantPerCap2.csv")

#Subset to the period of the evolution of food production (Neolithic)
pc30cp2<-subset(pc30cp, calBP<11500 & calBP>7000)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(pc30cp2$MKDE-min(pc30cp2$MKDE))/(max((pc30cp2$MKDE)-min(pc30cp2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
pc30cp3<-cbind(StKDE,pc30cp2)

Cpc2 <- ggplot(pc30cp3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_continuous(breaks=c(0,.25, .5, .75, 1), limits=c(0,1.0))+
  scale_y_continuous(limits=c(-.06,0.05))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A. SW Asia KDE Per Capita Growth vs. Density")+
  geom_vline(xintercept = 0.20)+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
Cpc2

Cpt <- ggplot(pc30cp3,aes(x=(calBP), y=StKDE)) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(color=PeriodID), size=3.5) +
  geom_path(aes(y=StKDE),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  scale_y_continuous(limits=c(0,1))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. SW Asia Standardized Density vs. Time")+
  geom_hline(yintercept = 0.20)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
Cpt

#SW Asia paired plot===============================
FigSWAsia<-plot_grid(Cpc2, Cpt, ncol=2, align="hv", axis = "rl")
FigSWAsia

pdf("data/Figs/ExSWAsia1.pdf", width=20.55, height=14)
FigSWAsia
dev.off()

#####Mid continent of North America, Case ID 3==============================================
##===================================================
#Load radiocarbon ages from P3k data set.
SPD<-read.csv(file="data/RawP3Kc14.csv", header=T)
#subset the data to the specified region
boxz<- subset(SPD, Latitude>35 & Latitude<41 & Longitude>-93 & Longitude< -87.5)
#boxy<- subset(SPD, Latitude>35 & Latitude<41 & Longitude>-87.51 & Longitude< -80.5)

#remove NA's from the siteID column
boxz <- boxz[!is.na(boxz$SiteID), ]

#very quick and dirty map
counties<-map_data("state")

ArchGlobeMap<-ggplot() +
  geom_polygon(data = counties, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  #geom_polygon(data = canada, aes(x=long, y = lat, group = group),
  #    fill = "white", color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Longitude", y="Latitude", title = "US ArchaeoGlobe Regions and Radiocarbon")+
  #geom_point(data=boxy, aes(Longitude, Latitude, color=factor(Province)),
  #  inherit.aes = FALSE, alpha = 0.5, size = 2)+
  geom_point(data=boxz, aes(Longitude, Latitude, color=factor(Province)),
             inherit.aes = FALSE, alpha = 0.5, size = 2, shape=23)
ArchGlobeMap

#calibrate dates
CalMz <- calibrate(x = boxz$Age,  errors = boxz$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = boxz$SiteID, ages = SPD$Age, h = 100)

####Run SPD
spd.midus <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.midus, runm=200, xlim=c(8200,200), type="simple")

calBP<-spd.midus$grid$calBP
PrDens<-spd.midus$grid$PrDens

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
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

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/USDomestKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/USDomestKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
##Add in the 30 year bin dates
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)

write.table(sums, file = "data/Sumbin/USDomestSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("data/Sumbin/USDomestSumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

##Add ID variable for culture history periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 580, 1060, 1570, 2050, 3000, 6000, 8200),
                         labels=c('Oneota', 'Mississippian', 'Late Woodland', 'Middle Woodland','Early Woodland','Late Archic','Middle Archaic'))

write.table(pcgrowth, file = "data/Percapita/USDomestPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
us30pc<- read.csv("data/Percapita/USDomestPerCap.csv")

#Initial Domestication: cir. 4200
us30pc2<-subset(us30pc, calBP<4201 & calBP>499)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(us30pc2$MKDE-min(us30pc2$MKDE))/(max((us30pc2$MKDE)-min(us30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
us30pc3<-cbind(StKDE,us30pc2)

pcUSdom <- ggplot(us30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.25,0.2))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Mean KDE (density)", y="KDE per capita growth", title = "A. US MidContinent KDE Per Capita Growth vs. Density")+
  geom_vline(xintercept = .2)+
  geom_hline(yintercept = 0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcUSdom

uscpt <- ggplot(us30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Mean KDE", title = "B. US Midcontinent Density vs. Time")+
  geom_hline(yintercept = .2)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
uscpt

#paired plot===============================
Figmidcont<-plot_grid(pcUSdom, uscpt, ncol=2, align="hv", axis = "rl")
Figmidcont

pdf("data/Figs/Exmidcont.pdf", width=20.55, height=14)
Figmidcont
dev.off()

####Eastern China (Yangtze River Basin and points north) Case ID 20====================================================
##======================================================
box<- read.csv("data/RawP3Kc14.csv")
box2<- subset(box, Country %in% c("China", "Mongolia")& Longitude< 121 & Longitude> 111 
              & Latitude>28.5  & Latitude<37)
#box3<-subset(box, Country %in% c("China", "Mongolia")& Longitude< 121 & Longitude> 110 
# & Latitude>37  & Latitude<42)

world <- ne_countries(scale = "medium", returnclass = "sf")

# Define the list of Middle Eastern countries
eastasia <- c("China", "Mongolia", "North Korea", "South Korea")

# Filter the data
eastasia<- world[world$name %in% eastasia, ]

# Load rivers data (includes Yangtze)
rivers <- ne_download(scale = "medium", type = "rivers_lake_centerlines", category = "physical", returnclass = "sf")

# Filter for Yangtze River
yangtze <- rivers %>% filter(name == "Yangtze")

ggplot(data = eastasia) +
  geom_sf(fill = "NA", color = "black") +
  geom_sf(data = yangtze, color = "blue", size = 1) +  # Add Yangtze River
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Country)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)+
  #geom_map(map=middle_east)+
  # coord_map("moll")+
  #coord_sf(crs = "+proj=aeqd +lat_0=30 +lon_0=30") +
  theme_minimal() +
  labs(title = "Map of East Asia",
       subtitle = "",
       x = "Longitude",
       y = "Latitude")

###Calibrate the radiocarbon ages
cptcal <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$SiteName, ages = box2$Age, h = 100)


####Run analysis for component 3 logistic 3500 to 150
spd.echina <- spd(cptcal, bins=NA, runm=150, timeRange=c(15000,2000))
plot(spd.echina, runm=150, xlim=c(15000,2000), type="simple")

calBP<-spd.echina$grid$calBP
PrDens<-spd.echina$grid$PrDens

##KDE=============================================
##KDE
####make KDEs
US.randates = sampleDates(cptcal, bins=NA, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(15000,2000),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')
D.ckde$timeRange

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

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
#dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd, file = "data/KDEs/China1KDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/China1KDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')


###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)

calBP<-c(14970,14940,14910,14880,14850,14820,14790,14760,14730,14700,14670,14640,14610,14580,14550,14520,14490,14460,14430,14400,14370,14340,14310,14280,14250,14220,14190,14160,14130,14100,14070,14040,14010,13980,13950,13920,13890,13860,13830,
         13800,13770,13740,13710,13680,13650,13620,13590,13560,13530,13500,13470,13440,13410,13380,13350,13320,13290,13260,13230,13200,13170,13140,13110,13080,13050,13020,12990,12960,12930,12900,12870,12840,12810,12780,12750,12720,12690,12660,12630,12600,12570,12540,12510,12480,12450,12420,12390,12360,12330,12300,12270,
         12240,12210,12180, 12150,12120,12090,12060,12030,12000,11970,11940,11910,11880,11850,11820,11790,11760,11730,11700,11670,11640,11610,11580,11550,11520,11490,11460,11430,11400,11370,11340,11310,11280,11250,11220,11190,11160,11130,11100,11070,11040,11010,10980,10950,10920,10890,10860,10830,10800,10770,10740,10710,10680,10650,10620,10590,10560,10530,10500,10470,10440,10410,10380,10350,10320,10290,10260,10230,10200,10170,10140,
         10110,10080,10050,10020,9990,9960,9930,9900,9870,9840,9810,9780,9750,9720,9690,9660,9630,9600,9570,9540,9510,9480,9450,9420,9390,9360,9330,9300,9270,9240,9210,9180,9150,9120,9090,9060,9030,9000,8970,8940,8910,8880,8850,8820,
         8790,8760,8730,8700,8670,8640,8610,8580,8550,8520,8490,8460,8430,8400,8370,8340,8310,8280,8250,8220,8190,8160,8130,8100,8070,8040,8010,7980,7950,7920,7890,7860,7830,7800,7770,7740,7710,7680,7650,7620,7590,7560,7530,7500,7470,7440,7410,7380,7350,7320,7290,7260,7230,
         7200,7170,7140,7110,7080,7050,7020,6990,6960,6930,6900,6870,6840,6810,6780,6750,6720,6690,6660,6630,6600,6570,6540,6510,6480,6450,6420,6390,6360,
         6330,6300,6270,6240,6210,6180,6150,6120,6090,6060,6030,6000,5970,5940,5910,5880,5850,5820,5790,5760,5730,5700,5670,5640,5610,5580,5550,5520,5490,5460,5430,5400,5370,5340,5310,5280,5250,5220,5190,5160,5130,5100,5070,5040,5010,4980,4950,4920,4890,4860,4830,4800,4770,4740,4710,4680,4650,4620,4590,4560,4530,4500,
         4470,4440,4410,4380,4350,4320,4290,4260,4230,4200,4170,4140,4110,4080,4050,4020,3990,3960,3930,3900,3870,
         3840,3810,3780,3750,3720,3690,3660,3630,3600, 3570,3540,3510,3480,3450,3420,3390,3360,3330,3300,3270,3240,3210,3180,3150,3120,3090,3060,3030,3000,2970,2940,2910,2880,2850,2820,2790,2760,2730,2700,2670,2640,2610,2580,2550,2520,2490,2460,2430,
         2400,2370,2340,2310,2280,2250,2220,2190,2160,2130,2100,2070,2040,2010)

sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/China1Sumbin.csv", sep = ",", col.names=NA)


###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.
###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("data/Sumbin/China1Sumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
##Add ID variable for culture history periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(2000,3000, 4000, 6500, 8000, 10500, 11500, 15000),
                         labels=c('States', 'Bronze Age', 'Late Neolithic', 'Middle Neolithic','Early Neolithic','Early Holocence',' Pleistocene'))

write.table(pcgrowth, file = "data/Percapita/China1PerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
ch30pc<- read.csv("data/Percapita/China1PerCap.csv")

###Initial exploitation of rice at 10000
ch30pc2<-subset(ch30pc, calBP<11500 & calBP>7000)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(ch30pc2$MKDE-min(ch30pc2$MKDE))/(max((ch30pc2$MKDE)-min(ch30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
ch30pc3<-cbind(StKDE,ch30pc2)

pcch<- ggplot(ch30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.23,0.4))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE densit)", y="KDE per capita growth", title = "A. SE China KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcch

chcpt <- ggplot(ch30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years BP", y="Standardized KDE density", title = "B. SE China Density vs. Time")+
  geom_hline(yintercept=0.20)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
chcpt

Figch<-plot_grid(pcch, chcpt, ncol=2, align="hv", axis = "rl")
Figch

pdf("data/Figs/ExSEchina.pdf", width=20.55, height=14)
Figch
dev.off()

#####Exogenous Agriculture--------Start with Sonoran Desert
####Sonoran Desert, Case ID 5 ===========================================================
SPD<-read.csv(file="data/RawP3Kc14.csv", header=T)
boxsd<- subset(SPD, Latitude>23 & Latitude<35 & Longitude>-115 & Longitude< -108)
#write.table(boxsd, file = "NERDv4_0/SonoranDesert.csv", sep = ",", col.names=NA)
#boxsd<-read.csv(file="data/SonoranDesert.csv", header=T)
###MAP
counties<-map_data("state")

ArchGlobeMap<-ggplot() +
  geom_polygon(data = counties, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  #geom_polygon(data = canada, aes(x=long, y = lat, group = group),
  #    fill = "white", color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Longitude", y="Latitude", title = "US ArchaeoGlobe Regions and Radiocarbon")+
  geom_point(data=boxsd, aes(Longitude, Latitude, color=factor(Province)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)
ArchGlobeMap
#remove NA's from the siteID column
boxsd <- boxsd[!is.na(boxsd$SiteID), ]

CalMz <- calibrate(x = boxsd$Age,  errors = boxsd$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = boxsd$SiteID, ages = boxsd$Age, h = 100)

####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.mz, runm=200, xlim=c(8200,200), type="simple")

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')
D.ckde$timeRange

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

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/SonoranKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/SonoranKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)

calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)

write.table(sums, file = "data/Sumbin/SonoranSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("data/Sumbin/SonoranSumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

##Add ID variable for culture history periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200 ,700, 1780, 2800, 3130, 4050, 4500, 7500, 8200),
                         labels=c('Post-Hohokam','Hohokam', 'Cienega', 'San Pedro', 'Silver Bell','Late Archaic','Middle Archaic','Early Archaic'))


write.table(pcgrowth, file = "data/Percapita/SonoranPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
son30pc<- read.csv("data/Percapita/SonoranPerCap.csv")

son30pc2<-subset(son30pc, calBP<4900 & calBP>400)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(son30pc2$MKDE-min(son30pc2$MKDE))/(max((son30pc2$MKDE)-min(son30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
son30pc3<-cbind(StKDE, son30pc2)

pcSon <- ggplot(son30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.2,0.2))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A. Sonoran KDE Per Capita Growth vs. Density")+
  geom_vline(xintercept = 0.20)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcSon

Soncpt <- ggplot(son30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. Sonoran Density vs. Time")+
  geom_hline(yintercept = 0.20)
#geom_vline(xintercept = 2550)+
#geom_vline(xintercept = 2250)+
#geom_vline(xintercept = 830)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
Soncpt

Figsonora<-plot_grid(pcSon, Soncpt, ncol=1, align="hv", axis = "rl")
Figsonora

pdf("data/Figs/Sonora.pdf", width=17.55, height=15)
Figsonora
dev.off()

###Yucatan (Lowland Maya) Case ID:13=====================================================
#=====================================================================================
SPD.1<-read.csv(file="data/Yucatan.csv", header=T)

cptcal <- calibrate(x = SPD.1$Age,  errors = SPD.1$Error, calCurves = "intcal20",  normalised = FALSE)
bin<- binPrep(sites = SPD.1$SiteName, ages = SPD.1$Age, h = 100)

####Run SPD
spd.yuc <- spd(cptcal, bins=bin, runm=250, timeRange=c(8200,200))
plot(spd.yuc, runm=100, xlim=c(8200, 200), type="simple")

PrDens<-spd.yuc$grid$PrDens
calBP<-spd.yuc$grid$calBP

##KDE
####make KDEs
US.randates = sampleDates(cptcal, bins=bin, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')
D.ckde$timeRange

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

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-subset(dd, MKDE>0)
##Write the table
write.table(dd2, file = "data/KDEs/YucaKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/YucaKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)

calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)

write.table(sums, file = "data/Sumbin/YucaSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.
###Not working........
d <-read.csv("data/Sumbin/YucaSumbin.csv") 
#%>%
# dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
##Add ID variable for culture history periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 1200, 1710, 4800, 8200),
                         labels=c('Post-Classic','Classic','Formative','Archaic'))
write.table(pcgrowth, file = "data/Percapita/YucaPerCap.csv", sep = ",", col.names=NA)


###Plot mean KDE against the per capita growth rate in the North
yuc30pc<- read.csv("data/Percapita/YucaPerCap.csv")

yuc30pc2<-subset(yuc30pc, calBP<4700 & calBP>500)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(yuc30pc2$MKDE-min(yuc30pc2$MKDE))/(max((yuc30pc2$MKDE)-min(yuc30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
yuc30pc3<-cbind(StKDE, yuc30pc2)

pcYuc <- ggplot(yuc30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A. Yucatan KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcYuc

Yuccpt <- ggplot(yuc30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. Yucatan Density vs. Time")+
  geom_hline(yintercept=.2)
#geom_vline(xintercept = 3320)+
#geom_vline(xintercept = 2550)+
#geom_vline(xintercept = 2250)+
#geom_vline(xintercept = 830)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
Yuccpt

FigYuc<-plot_grid(pcYuc, Yuccpt, ncol=2, align="hv", axis = "rl")
FigYuc

pdf("data/Figs/ExYuc.pdf", width=20.55, height=14)
FigYuc
dev.off()

#####################Fremont CASE ID 7=======================================================

SPD<-read.csv(file="data/RawP3Kc14.csv", header=T)
boxsd<- subset(SPD, Latitude>38 & Latitude<42 & Longitude>-115 & Longitude< -109)
#write.table(boxsd, file = "data/Fremontdates.csv", sep = ",", col.names=NA)
#boxsd<-read.csv(file="data/Fremontdates.csv", header=T)
###MAP
counties<-map_data("state")

ArchGlobeMap<-ggplot() +
  geom_polygon(data = counties, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  #geom_polygon(data = canada, aes(x=long, y = lat, group = group),
  #    fill = "white", color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Longitude", y="Latitude", title = "US ArchaeoGlobe Regions and Radiocarbon")+
  geom_point(data=boxsd, aes(Longitude, Latitude, color=factor(Province)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)
ArchGlobeMap

boxsd <- boxsd[!is.na(boxsd$SiteID), ]

CalMz <- calibrate(x = boxsd$Age,  errors = boxsd$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = boxsd$SiteID, ages = boxsd$Age, h = 100)

####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.mz, runm=200, xlim=c(8200,200), type="simple")

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')
#D.ckde$timeRange

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

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<- subset(dd, MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/FremontKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/FremontKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/FremontSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read.csv("data/Sumbin/FremontSumbin.csv") 
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 700, 1850, 4800, 8200),
                         labels=c('Post-Fremont','Fremont','Late Archaic','Archaic'))

write.table(pcgrowth, file = "data/Percapita/FremontPerCap.csv", sep = ",", col.names=NA)


###Plot mean KDE against the per capita growth rate in the North
fre30pc<- read.csv("data/Percapita/FremontPerCap.csv")

fre30pc2<-subset(fre30pc, calBP<2700 & calBP>400)

StKDE<-(fre30pc2$MKDE-min(fre30pc2$MKDE))/(max((fre30pc2$MKDE)-min(fre30pc2$MKDE)))
fre30pc3<-cbind(StKDE,fre30pc2)

pcfre <- ggplot(fre30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.25,0.2))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="Per capita growth rate", title = "A. Fremont KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcfre

frecpt <- ggplot(fre30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Time (years cal BP)", y="Standardized KDE density", title = "B. Fremont Density vs. Time")+
  geom_hline(yintercept = 0.20)
#geom_vline(xintercept = 1400)+
#geom_vline(xintercept = 1100)+
#geom_vline(xintercept = 800)+
#annotate("text", x =2000, y = 1, label = "Period 1", size = 6)+
#annotate("text", x =1300, y = 1, label = "2", size = 6)+
#annotate("text", x =950, y = 1, label = "3", size = 6)+
#annotate("text", x =500, y = 1, label = "4", size = 6)
frecpt

Figfremont<-plot_grid(pcfre, frecpt, ncol=2, align="hv", axis = "rl")
Figfremont

pdf("data/figs/Exfremont.pdf", width=20.55, height=14)
Figfremont
dev.off()


#=Chihuahua Desert (Jornada US) Case ID 4===============================================================
#===========================================================
#Load radiocarbon
SPD<-read.csv(file="data/RawP3Kc14.csv", header=T)
#subset the data to the region of study
boxsd<- subset(SPD, Latitude>27 & Latitude<34 & Longitude>-108 & Longitude< -103)
#write.table(boxsd, file = "data/Jornadadates.csv", sep = ",", col.names=NA)
#boxsd<-read.csv(file="data/Jornadadates.csv", header=T)
###MAP
counties<-map_data("state")

ArchGlobeMap<-ggplot() +
  geom_polygon(data = counties, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  #geom_polygon(data = canada, aes(x=long, y = lat, group = group),
  #    fill = "white", color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Longitude", y="Latitude", title = "US ArchaeoGlobe Regions and Radiocarbon")+
  geom_point(data=boxsd, aes(Longitude, Latitude, color=factor(Province)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)
ArchGlobeMap

boxsd <- boxsd[!is.na(boxsd$SiteID), ]


CalMz <- calibrate(x = boxsd$Age,  errors = boxsd$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = boxsd$SiteID, ages = boxsd$Age, h = 100)

####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.mz, runm=200, xlim=c(8200,200), type="simple")

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')
D.ckde$timeRange

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

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<- subset(dd, MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/JornadaKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/JornadaKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)

write.table(sums, file = "data/Sumbin/JornadaSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read.csv("data/Sumbin/JornadaSumbin.csv") 
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

##Code cultural historical periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 650, 1000, 1800, 3800, 6000, 8200),
                         labels=c('Post-Pueblo','Pueblo','Formative','Early Agricultural Period','Middle Archaic','Early Archaic'))

write.table(pcgrowth, file = "data/Percapita/JornadaPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
jor30pc<- read.csv("data/Percapita/JornadaPerCap.csv")

jor30pc2<-subset(jor30pc, calBP<4000 & calBP>499)

StKDE<-(jor30pc2$MKDE-min(jor30pc2$MKDE))/(max((jor30pc2$MKDE)-min(jor30pc2$MKDE)))
jor30pc3<-cbind(StKDE, jor30pc2)

pcjor <- ggplot(jor30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.14,0.14))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A. Jornada KDE Per Capita Growth vs. Density")+
  geom_vline(xintercept = 0.20)+
  geom_hline(yintercept = 0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcjor

jorcpt <- ggplot(jor30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. Jornada Density vs. Time")+
  geom_hline(yintercept = 0.20)
#geom_vline(xintercept = 2250)+
#geom_vline(xintercept = 830)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
jorcpt

#paired plot===============================
Figjornada<-plot_grid(pcjor, jorcpt, ncol=2, align="hv", axis = "rl")
Figjornada

pdf("data/Figs/Exjornada.pdf", width=20.55, height=14)
Figjornada
dev.off()

###===Southern Colorado Plateau (Upland US Southwest) Case ID 6===========================

SPD<-read.csv(file="data/RawP3Kc14.csv", header=T)
boxsd<- subset(SPD, Latitude>34 & Latitude<38.01 & Longitude>-113 & Longitude< -105)
#write.table(boxsd, file = "data/ColoradoPlatdates.csv", sep = ",", col.names=NA)
#boxsd<-read.csv(file="data/ColoradoPlatdates.csv", header=T)
###MAP
counties<-map_data("state")

ArchGlobeMap<-ggplot() +
  geom_polygon(data = counties, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  #geom_polygon(data = canada, aes(x=long, y = lat, group = group),
  #    fill = "white", color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Longitude", y="Latitude", title = "US ArchaeoGlobe Regions and Radiocarbon")+
  geom_point(data=boxsd, aes(Longitude, Latitude, color=factor(Province)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)
ArchGlobeMap

boxsd <- boxsd[!is.na(boxsd$SiteID), ]

CalMz <- calibrate(x = boxsd$Age,  errors = boxsd$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = boxsd$SiteID, ages = boxsd$Age, h = 100)

####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.mz, runm=200, xlim=c(8200,200), type="simple")

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
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

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<- subset(dd, MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/ColoradoPlatKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/ColoradoPlatKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/ColoradoPlatSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read.csv("data/Sumbin/ColoradoPlatSumbin.csv") 
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

##Code cultural historical periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 320, 670, 820, 1050, 1210, 1455, 3430, 4000, 8200),
                         labels=c('PV','PIV','PIII','PII','PI','BMIII','BMII','Early Agricultural Period','Archaic'))

write.table(pcgrowth, file = "data/Percapita/ColoradoPlatPerCap.csv", sep = ",", col.names=NA)


###Plot mean KDE against the per capita growth rate in the North
col30pc<- read.csv("data/Percapita/ColoradoPlatPerCap.csv")

col30pc2<-subset(col30pc, calBP<4000 & calBP>500)

StKDE<-(col30pc2$MKDE-min(col30pc2$MKDE))/(max((col30pc2$MKDE)-min(col30pc2$MKDE)))
col30pc3<-cbind(StKDE, col30pc2)

pccol <- ggplot(col30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.3,0.15))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A. Upland SW US KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept = 0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pccol

colcpt <- ggplot(col30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. Upland SW US Density vs. Time")+
  geom_hline(yintercept = 0.20)
#geom_vline(xintercept = 3320)+
#geom_vline(xintercept = 2550)+
#geom_vline(xintercept = 2250)+
#geom_vline(xintercept = 830)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
colcpt

#paired plot===============================
FiguplandSW<-plot_grid(pccol, colcpt, ncol=2, align="hv", axis = "rl")
FiguplandSW

pdf("data/Figs/ExuplandSW.pdf", width=20.55, height=14)
FiguplandSW
dev.off()

###===Western Great Basin Case ID 11===============================
#=============================================================

SPD<-read.csv(file="data/RawP3Kc14.csv", header=T)
boxsd<- subset(SPD, Latitude>36 & Latitude<41 & Longitude>-119 & Longitude< -115)
#write.table(boxsd, file = "data/WesternBasindates.csv", sep = ",", col.names=NA)
#boxsd<-read.csv(file="data/WesternBasindates.csv", header=T)
###MAP
counties<-map_data("state")

ArchGlobeMap<-ggplot() +
  geom_polygon(data = counties, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  #geom_polygon(data = canada, aes(x=long, y = lat, group = group),
  #    fill = "white", color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Longitude", y="Latitude", title = "US ArchaeoGlobe Regions and Radiocarbon")+
  geom_point(data=boxsd, aes(Longitude, Latitude, color=factor(Province)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)
ArchGlobeMap

boxsd <- boxsd[!is.na(boxsd$SiteID), ]


CalMz <- calibrate(x = boxsd$Age,  errors = boxsd$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = boxsd$SiteID, ages = boxsd$Age, h = 100)

####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.mz, runm=200, xlim=c(8200,100), type="simple")

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')
#D.ckde$timeRange

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

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<- subset(dd, MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/WesternBasinKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/WesternBasinKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')
MKDE<-rowMeans(out50)
###Calculate the mean KDE of the 30 year sums of the 200 KDEs
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/WesternBasinSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read.csv("data/Sumbin/WesternBasinSumbin.csv") 
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

##Code cultural historical periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 800, 4000, 6970, 8200),
                         labels=c('Post Archaic','Late Archaic','Middle Archaic','Early Archaic'))
write.table(pcgrowth, file = "data/Percapita/WesternBasinPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
west30pc<- read.csv("data/Percapita/WesternBasinPerCap.csv")

west30pc2<-subset(west30pc, calBP<4700 & calBP>200)

StKDE<-(west30pc2$MKDE-min(west30pc2$MKDE))/(max((west30pc2$MKDE)-min(west30pc2$MKDE)))
west30pc3<-cbind(StKDE, west30pc2)

pcwest <- ggplot(west30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.19,.2))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A. Western Basin KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept = 0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcwest

westcpt <- ggplot(west30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. Western Great Basin Density vs. Time")+
  geom_hline(yintercept = 0.2)
#geom_vline(xintercept = 2550)+
#geom_vline(xintercept = 2250)+
#geom_vline(xintercept = 830)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
westcpt

#paired plot===============================
FigWGB<-plot_grid(pcwest, westcpt, ncol=2, align="hv", axis = "rl")
FigWGB

pdf("data/Figs/ExwestGB.pdf", width=20.55, height=14)
FigWGB
dev.off()

###Southern California Case ID 10=========================================

SPD<-read.csv(file="data/RawP3Kc14.csv", header=T)
boxsd<- subset(SPD, Latitude>27 & Latitude<36.01 & Longitude>-122 & Longitude< -116 &  Material != "SHELL")
#write.table(boxsd, file = "data/SoCaldates.csv", sep = ",", col.names=NA)
#boxsd<-read.csv(file="data/SoCaldates.csv", header=T)
###MAP
counties<-map_data("state")

ArchGlobeMap<-ggplot() +
  geom_polygon(data = counties, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  #geom_polygon(data = canada, aes(x=long, y = lat, group = group),
  #    fill = "white", color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Longitude", y="Latitude", title = "US ArchaeoGlobe Regions and Radiocarbon")+
  geom_point(data=boxsd, aes(Longitude, Latitude, color=factor(Province)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)
ArchGlobeMap

#remove NA's from the siteID column
boxsd <- boxsd[!is.na(boxsd$SiteID), ]


CalMz <- calibrate(x = boxsd$Age,  errors = boxsd$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = boxsd$SiteID, ages = boxsd$Age, h = 100)

####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(8200,100))
plot(spd.mz, runm=200, xlim=c(8200,100), type="simple")

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,100),bw=50, normalised = FALSE)
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

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<- subset(dd, MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/SoCalKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/SoCalKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
###Calculate the mean KDE of the 30 year sums of the 200 KDEs
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220,190,170,140,110)
sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/SoCalSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read.csv("data/Sumbin/SoCalSumbin.csv") 
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
##Code cultural historical periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 950, 2560, 6220, 8200),
                         labels=c('Late Precontact','Late Archaic','Middle Archaic','Early Archaic'))
write.table(pcgrowth, file = "data/Percapita/SoCalPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
SoCal30pc<- read.csv("data/Percapita/SoCalPerCap.csv")

SoCal30pc2<-subset(SoCal30pc, calBP<4700 & calBP>200)

StKDE<-(SoCal30pc2$MKDE-min(SoCal30pc2$MKDE))/(max((SoCal30pc2$MKDE)-min(SoCal30pc2$MKDE)))
SoCal30pc3<-cbind(StKDE, SoCal30pc2)


pcsocal <- ggplot(SoCal30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.15,0.15))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A. So. Cal. KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcsocal

socalcpt <- ggplot(SoCal30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. So. Cal. Density vs. Time")+
  geom_hline(yintercept=0.2)
#geom_vline(xintercept = 3320)+
#geom_vline(xintercept = 2550)+
#geom_vline(xintercept = 2250)+
#geom_vline(xintercept = 830)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
socalcpt

##Make paired plot
Figsocal<-plot_grid(pcsocal, socalcpt, ncol=2, align="hv", axis = "rl")
Figsocal

pdf("data/Figs/Exsocal.pdf", width=20.55, height=14)
Figsocal
dev.off()

##Med Basin---Iberia=================================================================

box<- read.csv("data/RawP3Kc14.csv")
box2<- subset(box, Country %in% c("Spain"))

#write.table(box2, file = "data/Spaindates.csv", sep = ",", col.names=NA)
#boxsd<-read.csv(file="data/SoCaldates.csv", header=T)

# Load world map data
world <- ne_countries(scale = "medium", returnclass = "sf")
europe <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")

# Define the list of countries
eu_countries <- c("Spain", "France", "Italy", "Greece")

# Filter the data for the Middle East
eu <- world[world$name %in% eu_countries, ]

ggplot(data = eu) +
  geom_sf(fill = "NA", color = "black") +
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Country)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)+
  #geom_map(map=middle_east)+
  # coord_map("moll")+
  #coord_sf(crs = "+proj=aeqd +lat_0=30 +lon_0=30") +
  theme_minimal() +
  labs(title = "Map of the Middle East",
       subtitle = "Countries in the Middle East",
       x = "Longitude",
       y = "Latitude")

CalMz <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$SiteName, ages = box2$Age, h = 100)

####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(15000,2000))
plot(spd.mz, runm=200, xlim=c(15000,2000), type="simple")

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(15000,2000),bw=50, normalised = FALSE)
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

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<- subset(dd, MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/SpainKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/SpainKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(14970,14940,14910,14880,14850,14820,14790,14760,14730,14700,14670,14640,14610,14580,14550,14520,14490,14460,14430,14400,14370,14340,14310,14280,14250,14220,14190,14160,14130,14100,14070,14040,14010,13980,13950,13920,13890,13860,13830,
         13800,13770,13740,13710,13680,13650,13620,13590,13560,13530,13500,13470,13440,13410,13380,13350,13320,13290,13260,13230,13200,13170,13140,13110,13080,13050,13020,12990,12960,12930,12900,12870,12840,12810,12780,12750,12720,12690,12660,12630,12600,12570,12540,12510,12480,12450,12420,12390,12360,12330,12300,12270,
         12240,12210,12180, 12150,12120,12090,12060,12030,12000,11970,11940,11910,11880,11850,11820,11790,11760,11730,11700,11670,11640,11610,11580,11550,11520,11490,11460,11430,11400,11370,11340,11310,11280,11250,11220,11190,11160,11130,11100,11070,11040,11010,10980,10950,10920,10890,10860,10830,10800,10770,10740,10710,10680,10650,10620,10590,10560,10530,10500,10470,10440,10410,10380,10350,10320,10290,10260,10230,10200,10170,10140,
         10110,10080,10050,10020,9990,9960,9930,9900,9870,9840,9810,9780,9750,9720,9690,9660,9630,9600,9570,9540,9510,9480,9450,9420,9390,9360,9330,9300,9270,9240,9210,9180,9150,9120,9090,9060,9030,9000,8970,8940,8910,8880,8850,8820,
         8790,8760,8730,8700,8670,8640,8610,8580,8550,8520,8490,8460,8430,8400,8370,8340,8310,8280,8250,8220,8190,8160,8130,8100,8070,8040,8010,7980,7950,7920,7890,7860,7830,7800,7770,7740,7710,7680,7650,7620,7590,7560,7530,7500,7470,7440,7410,7380,7350,7320,7290,7260,7230,
         7200,7170,7140,7110,7080,7050,7020,6990,6960,6930,6900,6870,6840,6810,6780,6750,6720,6690,6660,6630,6600,6570,6540,6510,6480,6450,6420,6390,6360,
         6330,6300,6270,6240,6210,6180,6150,6120,6090,6060,6030,6000,5970,5940,5910,5880,5850,5820,5790,5760,5730,5700,5670,5640,5610,5580,5550,5520,5490,5460,5430,5400,5370,5340,5310,5280,5250,5220,5190,5160,5130,5100,5070,5040,5010,4980,4950,4920,4890,4860,4830,4800,4770,4740,4710,4680,4650,4620,4590,4560,4530,4500,
         4470,4440,4410,4380,4350,4320,4290,4260,4230,4200,4170,4140,4110,4080,4050,4020,3990,3960,3930,3900,3870,
         3840,3810,3780,3750,3720,3690,3660,3630,3600, 3570,3540,3510,3480,3450,3420,3390,3360,3330,3300,3270,3240,3210,3180,3150,3120,3090,3060,3030,3000,2970,2940,2910,2880,2850,2820,2790,2760,2730,2700,2670,2640,2610,2580,2550,2520,2490,2460,2430,
         2400,2370,2340,2310,2280,2250,2220,2190,2160,2130,2100,2070,2040,2010)

sums<-cbind(calBP, MKDE, out50)

write.table(sums, file = "data/Sumbin/SpainSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read.csv("data/Sumbin/SpainSumbin.csv") 
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
##Code cultural historical periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(2000, 4000, 5500, 8050, 11500, 15000),
                         labels=c('Bronze Age','Calcolithic','Neolithic','Mesolithic','Paleolithic'))

write.table(pcgrowth, file = "data/Percapita/SpainPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
Spal30pc<- read.csv("data/Percapita/SpainPerCap.csv")

Spal30pc2<-subset(Spal30pc, calBP<8500 & calBP>5000)

StKDE<-(Spal30pc2$MKDE-min(Spal30pc2$MKDE))/(max((Spal30pc2$MKDE)-min(Spal30pc2$MKDE)))
Spal30pc3<-cbind(StKDE, Spal30pc2)

pcsp <- ggplot(Spal30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.1,0.2))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "C. Iberia KDE Per Capita Growth vs. Density")+
  geom_vline(xintercept = .20)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcsp

spcpt <- ggplot(Spal30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "D. Iberia Density vs. Time")+
  geom_hline(yintercept = 0.20)
#geom_smooth(se=FALSE)
#geom_vline(xintercept = 2550)+
#geom_vline(xintercept = 2250)+
#geom_vline(xintercept = 830)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
spcpt

#========== Med. Basin S. France Case ID======================================================
box<- read.csv("data/RawP3Kc14.csv")
box2<- subset(box, Latitude>40 & Latitude<46 & Longitude>1 & Longitude<9)

#write.table(box2, file = "data/Spaindates.csv", sep = ",", col.names=NA)
#boxsd<-read.csv(file="data/SoCaldates.csv", header=T)

# Load world map data
world <- ne_countries(scale = "medium", returnclass = "sf")
europe <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")

# Define the list of Middle Eastern countries
eu_countries <- c("Spain", "France", "Italy", "Greece")

# Filter the data for the Middle East
eu <- world[world$name %in% eu_countries, ]

ggplot(data = eu) +
  geom_sf(fill = "NA", color = "black") +
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Country)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)+
  #geom_map(map=middle_east)+
  # coord_map("moll")+
  #coord_sf(crs = "+proj=aeqd +lat_0=30 +lon_0=30") +
  theme_minimal() +
  labs(title = "Map of W Europe",
       subtitle = "France radiocarbon ages",
       x = "Longitude",
       y = "Latitude")

CalMz <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$SiteName, ages = box2$Age, h = 100)

####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(15000,2000))
plot(spd.mz, runm=200, xlim=c(15000,2000), type="simple")

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(15000,2000),bw=50, normalised = FALSE)
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

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<- subset(dd, MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/SFranceKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/SFranceKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(14970,14940,14910,14880,14850,14820,14790,14760,14730,14700,14670,14640,14610,14580,14550,14520,14490,14460,14430,14400,14370,14340,14310,14280,14250,14220,14190,14160,14130,14100,14070,14040,14010,13980,13950,13920,13890,13860,13830,
         13800,13770,13740,13710,13680,13650,13620,13590,13560,13530,13500,13470,13440,13410,13380,13350,13320,13290,13260,13230,13200,13170,13140,13110,13080,13050,13020,12990,12960,12930,12900,12870,12840,12810,12780,12750,12720,12690,12660,12630,12600,12570,12540,12510,12480,12450,12420,12390,12360,12330,12300,12270,
         12240,12210,12180, 12150,12120,12090,12060,12030,12000,11970,11940,11910,11880,11850,11820,11790,11760,11730,11700,11670,11640,11610,11580,11550,11520,11490,11460,11430,11400,11370,11340,11310,11280,11250,11220,11190,11160,11130,11100,11070,11040,11010,10980,10950,10920,10890,10860,10830,10800,10770,10740,10710,10680,10650,10620,10590,10560,10530,10500,10470,10440,10410,10380,10350,10320,10290,10260,10230,10200,10170,10140,
         10110,10080,10050,10020,9990,9960,9930,9900,9870,9840,9810,9780,9750,9720,9690,9660,9630,9600,9570,9540,9510,9480,9450,9420,9390,9360,9330,9300,9270,9240,9210,9180,9150,9120,9090,9060,9030,9000,8970,8940,8910,8880,8850,8820,
         8790,8760,8730,8700,8670,8640,8610,8580,8550,8520,8490,8460,8430,8400,8370,8340,8310,8280,8250,8220,8190,8160,8130,8100,8070,8040,8010,7980,7950,7920,7890,7860,7830,7800,7770,7740,7710,7680,7650,7620,7590,7560,7530,7500,7470,7440,7410,7380,7350,7320,7290,7260,7230,
         7200,7170,7140,7110,7080,7050,7020,6990,6960,6930,6900,6870,6840,6810,6780,6750,6720,6690,6660,6630,6600,6570,6540,6510,6480,6450,6420,6390,6360,
         6330,6300,6270,6240,6210,6180,6150,6120,6090,6060,6030,6000,5970,5940,5910,5880,5850,5820,5790,5760,5730,5700,5670,5640,5610,5580,5550,5520,5490,5460,5430,5400,5370,5340,5310,5280,5250,5220,5190,5160,5130,5100,5070,5040,5010,4980,4950,4920,4890,4860,4830,4800,4770,4740,4710,4680,4650,4620,4590,4560,4530,4500,
         4470,4440,4410,4380,4350,4320,4290,4260,4230,4200,4170,4140,4110,4080,4050,4020,3990,3960,3930,3900,3870,
         3840,3810,3780,3750,3720,3690,3660,3630,3600, 3570,3540,3510,3480,3450,3420,3390,3360,3330,3300,3270,3240,3210,3180,3150,3120,3090,3060,3030,3000,2970,2940,2910,2880,2850,2820,2790,2760,2730,2700,2670,2640,2610,2580,2550,2520,2490,2460,2430,
         2400,2370,2340,2310,2280,2250,2220,2190,2160,2130,2100,2070,2040,2010)

sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/SFranceSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read.csv("data/Sumbin/SFranceSumbin.csv") 
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

##Code cultural historical periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(2000, 4250, 8050, 11500, 15000),
                         labels=c('Bronze Age','Neolithic','Mesolithic','Paleolithic'))
write.table(pcgrowth, file = "data/Percapita/SFrancePerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
sfr30pc<- read.csv("data/Percapita/SFrancePerCap.csv")

sfr30pc2<-subset(sfr30pc, calBP<8500 & calBP>4500)

StKDE<-(sfr30pc2$MKDE-min(sfr30pc2$MKDE))/(max((sfr30pc2$MKDE)-min(sfr30pc2$MKDE)))
sfr30pc3<-cbind(StKDE, sfr30pc2)

pcsfr <- ggplot(sfr30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.1,0.15))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A. South France KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept = 0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcsfr

frcpt <- ggplot(sfr30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. South France Density vs. Time")+
  geom_hline(yintercept = 0.2)
#geom_vline(xintercept = 2550)+
#geom_vline(xintercept = 2250)+
#geom_vline(xintercept = 830)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
frcpt

Figsf<-plot_grid(pcsfr, frcpt, ncol=2, align="hv", axis = "rl")
Figsf

pdf("data/Figs/ExSfrance.pdf", width=20.55, height=14)
Figsf
dev.off()

#========== Med. Basin Italy======================================================
box<- read.csv("data/RawP3Kc14.csv")
box2<- subset(box, Latitude>36.5 & Latitude<44 & Longitude>10 & Longitude<17)

#write.table(box2, file = "data/Spaindates.csv", sep = ",", col.names=NA)
#boxsd<-read.csv(file="data/SoCaldates.csv", header=T)

# Load world map data
#world <- ne_countries(scale = "medium", returnclass = "sf")
europe <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")

# Define the list of Middle Eastern countries
eu_countries <- c("Spain", "France", "Italy", "Greece", "Albania", "Germany", "Turkey", "Croatia", "Slovinia")

# Filter the data for the Middle East
eu <- europe[europe$name %in% eu_countries, ]

ggplot(data = eu) +
  geom_sf(fill = "NA", color = "black") +
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Country)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)+
  #geom_map(map=middle_east)+
  # coord_map("moll")+
  #coord_sf(crs = "+proj=aeqd +lat_0=30 +lon_0=30") +
  theme_minimal() +
  labs(title = "Map of the Middle East",
       subtitle = "Countries in the Middle East",
       x = "Longitude",
       y = "Latitude")

CalMz <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$SiteName, ages = box2$Age, h = 100)

####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(15000,2000))
plot(spd.mz, runm=200, xlim=c(15000,2000), type="simple")

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(15000,2000),bw=50, normalised = FALSE)
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

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<- subset(dd, MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/ItalyKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/ItalyKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(14970,14940,14910,14880,14850,14820,14790,14760,14730,14700,14670,14640,14610,14580,14550,14520,14490,14460,14430,14400,14370,14340,14310,14280,14250,14220,14190,14160,14130,14100,14070,14040,14010,13980,13950,13920,13890,13860,13830,
         13800,13770,13740,13710,13680,13650,13620,13590,13560,13530,13500,13470,13440,13410,13380,13350,13320,13290,13260,13230,13200,13170,13140,13110,13080,13050,13020,12990,12960,12930,12900,12870,12840,12810,12780,12750,12720,12690,12660,12630,12600,12570,12540,12510,12480,12450,12420,12390,12360,12330,12300,12270,
         12240,12210,12180, 12150,12120,12090,12060,12030,12000,11970,11940,11910,11880,11850,11820,11790,11760,11730,11700,11670,11640,11610,11580,11550,11520,11490,11460,11430,11400,11370,11340,11310,11280,11250,11220,11190,11160,11130,11100,11070,11040,11010,10980,10950,10920,10890,10860,10830,10800,10770,10740,10710,10680,10650,10620,10590,10560,10530,10500,10470,10440,10410,10380,10350,10320,10290,10260,10230,10200,10170,10140,
         10110,10080,10050,10020,9990,9960,9930,9900,9870,9840,9810,9780,9750,9720,9690,9660,9630,9600,9570,9540,9510,9480,9450,9420,9390,9360,9330,9300,9270,9240,9210,9180,9150,9120,9090,9060,9030,9000,8970,8940,8910,8880,8850,8820,
         8790,8760,8730,8700,8670,8640,8610,8580,8550,8520,8490,8460,8430,8400,8370,8340,8310,8280,8250,8220,8190,8160,8130,8100,8070,8040,8010,7980,7950,7920,7890,7860,7830,7800,7770,7740,7710,7680,7650,7620,7590,7560,7530,7500,7470,7440,7410,7380,7350,7320,7290,7260,7230,
         7200,7170,7140,7110,7080,7050,7020,6990,6960,6930,6900,6870,6840,6810,6780,6750,6720,6690,6660,6630,6600,6570,6540,6510,6480,6450,6420,6390,6360,
         6330,6300,6270,6240,6210,6180,6150,6120,6090,6060,6030,6000,5970,5940,5910,5880,5850,5820,5790,5760,5730,5700,5670,5640,5610,5580,5550,5520,5490,5460,5430,5400,5370,5340,5310,5280,5250,5220,5190,5160,5130,5100,5070,5040,5010,4980,4950,4920,4890,4860,4830,4800,4770,4740,4710,4680,4650,4620,4590,4560,4530,4500,
         4470,4440,4410,4380,4350,4320,4290,4260,4230,4200,4170,4140,4110,4080,4050,4020,3990,3960,3930,3900,3870,
         3840,3810,3780,3750,3720,3690,3660,3630,3600, 3570,3540,3510,3480,3450,3420,3390,3360,3330,3300,3270,3240,3210,3180,3150,3120,3090,3060,3030,3000,2970,2940,2910,2880,2850,2820,2790,2760,2730,2700,2670,2640,2610,2580,2550,2520,2490,2460,2430,
         2400,2370,2340,2310,2280,2250,2220,2190,2160,2130,2100,2070,2040,2010)

sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/ItalySumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read.csv("data/Sumbin/ItalySumbin.csv") 
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

##Code cultural historical periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(2000, 4000, 5560, 8500, 11500, 15000),
                         labels=c('Bronze Age','Cooper Age','Neolithic','Mesolithic','Paleolithic'))

write.table(pcgrowth, file = "data/Percapita/ItalyPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
It30pc<- read.csv("data/Percapita/ItalyPerCap.csv")

It30pc2<-subset(It30pc, calBP<9000 & calBP>4000)

StKDE<-(It30pc2$MKDE-min(It30pc2$MKDE))/(max((It30pc2$MKDE)-min(It30pc2$MKDE)))
It30pc3<-cbind(StKDE, It30pc2)

pcIt <- ggplot(It30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.1,0.23))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A.Italy KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept = 0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcIt

Itcpt <- ggplot(It30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. Italy Density vs. Time")+
  geom_hline(yintercept = 0.2)
#geom_vline(xintercept = 3320)+
#geom_vline(xintercept = 2550)+
#geom_vline(xintercept = 2250)+
#geom_vline(xintercept = 830)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
Itcpt

##Plot Fig.
Figitaly<-plot_grid(pcIt, Itcpt, ncol=2, align="hv", axis = "rl")
Figitaly

pdf("data/Figs/ExItaly.pdf", width=20.55, height=14)
Figitaly
dev.off()

#========== Med. Basin Greece Case ID 23======================================================
box<- read.csv("data/RawP3Kc14.csv")
box2<- subset(box, Latitude>33 & Latitude<39 & Longitude>17 & Longitude<27)

#write.table(box2, file = "data/Spaindates.csv", sep = ",", col.names=NA)
#boxsd<-read.csv(file="data/SoCaldates.csv", header=T)

# Load world map data
#world <- ne_countries(scale = "medium", returnclass = "sf")
europe <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")

# Define the list of Middle Eastern countries
eu_countries <- c("Spain", "France", "Italy", "Greece", "Albania", "Germany", "Turkey", "Croatia", "Slovinia")

# Filter the data for the Middle East
eu <- europe[europe$name %in% eu_countries, ]

ggplot(data = eu) +
  geom_sf(fill = "NA", color = "black") +
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Country)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)+
  #geom_map(map=middle_east)+
  # coord_map("moll")+
  #coord_sf(crs = "+proj=aeqd +lat_0=30 +lon_0=30") +
  theme_minimal() +
  labs(title = "Map of the Middle East",
       subtitle = "Countries in the Middle East",
       x = "Longitude",
       y = "Latitude")

CalMz <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$SiteName, ages = box2$Age, h = 100)

####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(15000,2000))
plot(spd.mz, runm=200, xlim=c(15000,2000), type="simple")

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(15000,2000),bw=50, normalised = FALSE)
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

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<- subset(dd, MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/GreeceKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/GreeceKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........

# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(14970,14940,14910,14880,14850,14820,14790,14760,14730,14700,14670,14640,14610,14580,14550,14520,14490,14460,14430,14400,14370,14340,14310,14280,14250,14220,14190,14160,14130,14100,14070,14040,14010,13980,13950,13920,13890,13860,13830,
         13800,13770,13740,13710,13680,13650,13620,13590,13560,13530,13500,13470,13440,13410,13380,13350,13320,13290,13260,13230,13200,13170,13140,13110,13080,13050,13020,12990,12960,12930,12900,12870,12840,12810,12780,12750,12720,12690,12660,12630,12600,12570,12540,12510,12480,12450,12420,12390,12360,12330,12300,12270,
         12240,12210,12180, 12150,12120,12090,12060,12030,12000,11970,11940,11910,11880,11850,11820,11790,11760,11730,11700,11670,11640,11610,11580,11550,11520,11490,11460,11430,11400,11370,11340,11310,11280,11250,11220,11190,11160,11130,11100,11070,11040,11010,10980,10950,10920,10890,10860,10830,10800,10770,10740,10710,10680,10650,10620,10590,10560,10530,10500,10470,10440,10410,10380,10350,10320,10290,10260,10230,10200,10170,10140,
         10110,10080,10050,10020,9990,9960,9930,9900,9870,9840,9810,9780,9750,9720,9690,9660,9630,9600,9570,9540,9510,9480,9450,9420,9390,9360,9330,9300,9270,9240,9210,9180,9150,9120,9090,9060,9030,9000,8970,8940,8910,8880,8850,8820,
         8790,8760,8730,8700,8670,8640,8610,8580,8550,8520,8490,8460,8430,8400,8370,8340,8310,8280,8250,8220,8190,8160,8130,8100,8070,8040,8010,7980,7950,7920,7890,7860,7830,7800,7770,7740,7710,7680,7650,7620,7590,7560,7530,7500,7470,7440,7410,7380,7350,7320,7290,7260,7230,
         7200,7170,7140,7110,7080,7050,7020,6990,6960,6930,6900,6870,6840,6810,6780,6750,6720,6690,6660,6630,6600,6570,6540,6510,6480,6450,6420,6390,6360,
         6330,6300,6270,6240,6210,6180,6150,6120,6090,6060,6030,6000,5970,5940,5910,5880,5850,5820,5790,5760,5730,5700,5670,5640,5610,5580,5550,5520,5490,5460,5430,5400,5370,5340,5310,5280,5250,5220,5190,5160,5130,5100,5070,5040,5010,4980,4950,4920,4890,4860,4830,4800,4770,4740,4710,4680,4650,4620,4590,4560,4530,4500,
         4470,4440,4410,4380,4350,4320,4290,4260,4230,4200,4170,4140,4110,4080,4050,4020,3990,3960,3930,3900,3870,
         3840,3810,3780,3750,3720,3690,3660,3630,3600, 3570,3540,3510,3480,3450,3420,3390,3360,3330,3300,3270,3240,3210,3180,3150,3120,3090,3060,3030,3000,2970,2940,2910,2880,2850,2820,2790,2760,2730,2700,2670,2640,2610,2580,2550,2520,2490,2460,2430,
         2400,2370,2340,2310,2280,2250,2220,2190,2160,2130,2100,2070,2040,2010)

sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/GreeceSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read.csv("data/Sumbin/GreeceSumbin.csv") 
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
##Code cultural historical periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(2000, 4000, 5560, 8500, 11500, 15000),
                         labels=c('Bronze Age','Cooper Age','Neolithic','Mesolithic','Paleolithic'))
write.table(pcgrowth, file = "data/Percapita/GreecePerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
gree30pc<- read.csv("data/Percapita/GreecePerCap.csv")

gree30pc2<-subset(gree30pc, calBP<9500 & calBP>5000)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(gree30pc2$MKDE-min(gree30pc2$MKDE))/(max((gree30pc2$MKDE)-min(gree30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
gree30pc3<-cbind(StKDE,gree30pc2)

pcgree <- ggplot(gree30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.1,0.2))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A. Greece KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept = 0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcgree

greecpt <- ggplot(gree30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. Greece Density vs. Time")+
  geom_hline(yintercept = 0.2)+
  geom_vline(xintercept = 7250)
#geom_vline(xintercept = 2250)+
#geom_vline(xintercept = 830)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
greecpt

Figgree<-plot_grid(pcgree, greecpt, ncol=2, align="hv", axis = "rl")
Figgree

pdf("data/ExGreece.pdf", width=20.55, height=14)
Figgree
dev.off()

##N. France Case ID 29=====================================================
box<- read.csv("data/RawP3Kc14.csv")
box2<- subset(box, Latitude>46 & Latitude<49 & Longitude>-3 & Longitude<6)

#write.table(box2, file = "data/NorFrancedates.csv", sep = ",", col.names=NA)
box2<-read.csv(file="data/NorFrancedates.csv", header=T)

# Load map data
europe <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")

# Define the list of Middle Eastern countries
eu_countries <- c("Spain", "France", "Italy", "Greece", "Albania", "Germany", "United Kingdom", "Croatia", "Slovinia")

# Filter the data for the Middle East
eu <- europe[europe$name %in% eu_countries, ]

ggplot(data = eu) +
  geom_sf(fill = "NA", color = "black") +
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Country)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)+
  #geom_map(map=middle_east)+
  # coord_map("moll")+
  #coord_sf(crs = "+proj=aeqd +lat_0=30 +lon_0=30") +
  theme_minimal() +
  labs(title = "Map of the Middle East",
       subtitle = "Countries in the Middle East",
       x = "Longitude",
       y = "Latitude")


CalMz <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$SiteName, ages = box2$Age, h = 100)
#write.table(boxbins, file = "data/NorFranceboxbins.csv", sep = ",", col.names=NA)
####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(15000,2000))
plot(spd.mz, runm=200, xlim=c(15000,2000), type="simple")

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(15000,2000),bw=50, normalised = FALSE)
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

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<- subset(dd, MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/NorFranceKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/NorFranceKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(14970,14940,14910,14880,14850,14820,14790,14760,14730,14700,14670,14640,14610,14580,14550,14520,14490,14460,14430,14400,14370,14340,14310,14280,14250,14220,14190,14160,14130,14100,14070,14040,14010,13980,13950,13920,13890,13860,13830,
         13800,13770,13740,13710,13680,13650,13620,13590,13560,13530,13500,13470,13440,13410,13380,13350,13320,13290,13260,13230,13200,13170,13140,13110,13080,13050,13020,12990,12960,12930,12900,12870,12840,12810,12780,12750,12720,12690,12660,12630,12600,12570,12540,12510,12480,12450,12420,12390,12360,12330,12300,12270,
         12240,12210,12180, 12150,12120,12090,12060,12030,12000,11970,11940,11910,11880,11850,11820,11790,11760,11730,11700,11670,11640,11610,11580,11550,11520,11490,11460,11430,11400,11370,11340,11310,11280,11250,11220,11190,11160,11130,11100,11070,11040,11010,10980,10950,10920,10890,10860,10830,10800,10770,10740,10710,10680,10650,10620,10590,10560,10530,10500,10470,10440,10410,10380,10350,10320,10290,10260,10230,10200,10170,10140,
         10110,10080,10050,10020,9990,9960,9930,9900,9870,9840,9810,9780,9750,9720,9690,9660,9630,9600,9570,9540,9510,9480,9450,9420,9390,9360,9330,9300,9270,9240,9210,9180,9150,9120,9090,9060,9030,9000,8970,8940,8910,8880,8850,8820,
         8790,8760,8730,8700,8670,8640,8610,8580,8550,8520,8490,8460,8430,8400,8370,8340,8310,8280,8250,8220,8190,8160,8130,8100,8070,8040,8010,7980,7950,7920,7890,7860,7830,7800,7770,7740,7710,7680,7650,7620,7590,7560,7530,7500,7470,7440,7410,7380,7350,7320,7290,7260,7230,
         7200,7170,7140,7110,7080,7050,7020,6990,6960,6930,6900,6870,6840,6810,6780,6750,6720,6690,6660,6630,6600,6570,6540,6510,6480,6450,6420,6390,6360,
         6330,6300,6270,6240,6210,6180,6150,6120,6090,6060,6030,6000,5970,5940,5910,5880,5850,5820,5790,5760,5730,5700,5670,5640,5610,5580,5550,5520,5490,5460,5430,5400,5370,5340,5310,5280,5250,5220,5190,5160,5130,5100,5070,5040,5010,4980,4950,4920,4890,4860,4830,4800,4770,4740,4710,4680,4650,4620,4590,4560,4530,4500,
         4470,4440,4410,4380,4350,4320,4290,4260,4230,4200,4170,4140,4110,4080,4050,4020,3990,3960,3930,3900,3870,
         3840,3810,3780,3750,3720,3690,3660,3630,3600, 3570,3540,3510,3480,3450,3420,3390,3360,3330,3300,3270,3240,3210,3180,3150,3120,3090,3060,3030,3000,2970,2940,2910,2880,2850,2820,2790,2760,2730,2700,2670,2640,2610,2580,2550,2520,2490,2460,2430,
         2400,2370,2340,2310,2280,2250,2220,2190,2160,2130,2100,2070,2040,2010)

sums<-cbind(calBP, MKDE, out50)

write.table(sums, file = "data/Sumbin/NorFranceSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read.csv("data/Sumbin/NorFranceSumbin.csv") 
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(2000, 4320, 7800, 11500, 15000),
                         labels=c('Metal Age','Neolithic','Mesolithic','Paleolithic'))

write.table(pcgrowth, file = "data/Percapita/NorFrancePerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
nf30pc<- read.csv("data/Percapita/NorFrancePerCap.csv")

nf30pc2<-subset(nf30pc, calBP<8500 & calBP>4500)

StKDE<-(nf30pc2$MKDE-min(nf30pc2$MKDE))/(max((nf30pc2$MKDE)-min(nf30pc2$MKDE)))
nf30pc3<-cbind(StKDE, nf30pc2)

pcnfr <- ggplot(nf30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.1,0.2))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A. N. France KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept = 0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcnfr

nfrcpt <- ggplot(nf30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. N. France Density vs. Time")+
  geom_hline(yintercept = 0.2)
#geom_vline(xintercept = 3320)+
#geom_vline(xintercept = 2550)+
#geom_vline(xintercept = 2250)+
#geom_vline(xintercept = 830)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
nfrcpt

##Paired plots
Fignf<-plot_grid(pcnfr, nfrcpt, ncol=2, align="hv", axis = "rl")
Fignf

pdf("data/Figs/EXNFrance.pdf", width=20.55, height=14)
Fignf
dev.off()

##West Central Europe Case ID 22=====================================================
box<- read.csv("data/RawP3Kc14.csv")
box2<- subset(box, Latitude>47 & Latitude<52 & Longitude>7 & Longitude<14)

# Load map data
europe <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")

# Define the list of Middle Eastern countries
eu_countries <- c("Spain", "France", "Italy", "Greece", "Albania", "Germany", "United Kingdom", "Croatia", "Slovinia")

# Filter the data for the Middle East
eu <- europe[europe$name %in% eu_countries, ]

ggplot(data = eu) +
  geom_sf(fill = "NA", color = "black") +
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Country)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)+
  #geom_map(map=middle_east)+
  # coord_map("moll")+
  #coord_sf(crs = "+proj=aeqd +lat_0=30 +lon_0=30") +
  theme_minimal() +
  labs(title = "Map of Central Eurropean Ages",
       x = "Longitude",
       y = "Latitude")


CalMz <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$SiteName, ages = box2$Age, h = 100)
#write.table(boxbins, file = "data/NorFranceboxbins.csv", sep = ",", col.names=NA)
####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(15000,2000))
plot(spd.mz, runm=200, xlim=c(15000,2000), type="simple")

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(15000,2000),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')
#D.ckde$timeRange

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

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<- subset(dd, MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/WestCenEuropeKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/WestCenEuropeKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(14970,14940,14910,14880,14850,14820,14790,14760,14730,14700,14670,14640,14610,14580,14550,14520,14490,14460,14430,14400,14370,14340,14310,14280,14250,14220,14190,14160,14130,14100,14070,14040,14010,13980,13950,13920,13890,13860,13830,
         13800,13770,13740,13710,13680,13650,13620,13590,13560,13530,13500,13470,13440,13410,13380,13350,13320,13290,13260,13230,13200,13170,13140,13110,13080,13050,13020,12990,12960,12930,12900,12870,12840,12810,12780,12750,12720,12690,12660,12630,12600,12570,12540,12510,12480,12450,12420,12390,12360,12330,12300,12270,
         12240,12210,12180, 12150,12120,12090,12060,12030,12000,11970,11940,11910,11880,11850,11820,11790,11760,11730,11700,11670,11640,11610,11580,11550,11520,11490,11460,11430,11400,11370,11340,11310,11280,11250,11220,11190,11160,11130,11100,11070,11040,11010,10980,10950,10920,10890,10860,10830,10800,10770,10740,10710,10680,10650,10620,10590,10560,10530,10500,10470,10440,10410,10380,10350,10320,10290,10260,10230,10200,10170,10140,
         10110,10080,10050,10020,9990,9960,9930,9900,9870,9840,9810,9780,9750,9720,9690,9660,9630,9600,9570,9540,9510,9480,9450,9420,9390,9360,9330,9300,9270,9240,9210,9180,9150,9120,9090,9060,9030,9000,8970,8940,8910,8880,8850,8820,
         8790,8760,8730,8700,8670,8640,8610,8580,8550,8520,8490,8460,8430,8400,8370,8340,8310,8280,8250,8220,8190,8160,8130,8100,8070,8040,8010,7980,7950,7920,7890,7860,7830,7800,7770,7740,7710,7680,7650,7620,7590,7560,7530,7500,7470,7440,7410,7380,7350,7320,7290,7260,7230,
         7200,7170,7140,7110,7080,7050,7020,6990,6960,6930,6900,6870,6840,6810,6780,6750,6720,6690,6660,6630,6600,6570,6540,6510,6480,6450,6420,6390,6360,
         6330,6300,6270,6240,6210,6180,6150,6120,6090,6060,6030,6000,5970,5940,5910,5880,5850,5820,5790,5760,5730,5700,5670,5640,5610,5580,5550,5520,5490,5460,5430,5400,5370,5340,5310,5280,5250,5220,5190,5160,5130,5100,5070,5040,5010,4980,4950,4920,4890,4860,4830,4800,4770,4740,4710,4680,4650,4620,4590,4560,4530,4500,
         4470,4440,4410,4380,4350,4320,4290,4260,4230,4200,4170,4140,4110,4080,4050,4020,3990,3960,3930,3900,3870,
         3840,3810,3780,3750,3720,3690,3660,3630,3600, 3570,3540,3510,3480,3450,3420,3390,3360,3330,3300,3270,3240,3210,3180,3150,3120,3090,3060,3030,3000,2970,2940,2910,2880,2850,2820,2790,2760,2730,2700,2670,2640,2610,2580,2550,2520,2490,2460,2430,
         2400,2370,2340,2310,2280,2250,2220,2190,2160,2130,2100,2070,2040,2010)

sums<-cbind(calBP, MKDE, out50)

write.table(sums, file = "data/Sumbin/WestCenEuropeSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read.csv("data/Sumbin/WestCenEuropeSumbin.csv") 
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(2000, 3960, 8010, 11500, 15000),
                         labels=c('Metal Age','Neolithic','Mesolithic','Paleolithic'))
write.table(pcgrowth, file = "data/Percapita/WestCenEuropePerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
wce30pc<- read.csv("data/Percapita/WestCenEuropePerCap.csv")

wce30pc2<-subset(wce30pc, calBP<8500 & calBP>5950)

StKDE<-(wce30pc2$MKDE-min(wce30pc2$MKDE))/(max((wce30pc2$MKDE)-min(wce30pc2$MKDE)))
wce30pc3<-cbind(StKDE, wce30pc2)

pcwce <- ggplot(wce30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.15,0.15))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", 
       title = "A. West Central Europe KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept = 0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcwce

wcecpt <- ggplot(wce30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. West Central Europe Density vs. Time")+
  geom_hline(yintercept = 0.2)
#geom_vline(xintercept = 3320)+
#geom_vline(xintercept = 2550)+
#geom_vline(xintercept = 2250)+
#geom_vline(xintercept = 830)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
wcecpt

##Paired Plots
Figwce<-plot_grid(pcwce, wcecpt, ncol=2, align="hv", axis = "rl")
Figwce

pdf("data/ExWCEurope.pdf", width=20.55, height=14)
Figwce
dev.off()

####S. Britain Case ID 27=============================================
box<- read.csv("NERDv4_0/RawP3Kc14.csv")
box2<- subset(box, Latitude>51 & Latitude<55 & Longitude>-5 & Longitude<2)

#write.table(box2, file = "data/SouthBritdates.csv", sep = ",", col.names=NA)
#box2<-read.csv(file="data/SouthBritdates.csv", header=T)

# Load world map data
#world <- ne_countries(scale = "medium", returnclass = "sf")
europe <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")

# Define the list of Middle Eastern countries
eu_countries <- c("Spain", "France", "Italy", "Greece", "Albania", "Germany","Belgium", "Sweeden", "United Kingdom", "Croatia", "Slovinia", "Norway", "Netherlands")

# Filter the data for the Middle East
eu <- europe[europe$name %in% eu_countries, ]

ggplot(data = eu) +
  geom_sf(fill = "NA", color = "black") +
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Country)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)+
  #geom_map(map=middle_east)+
  # coord_map("moll")+
  #coord_sf(crs = "+proj=aeqd +lat_0=30 +lon_0=30") +
  theme_minimal() +
  labs(title = "Map of the Middle East",
       subtitle = "Countries in the Middle East",
       x = "Longitude",
       y = "Latitude")
#remove NA's from the siteID column
box2 <- box2[!is.na(box2$SiteName), ]

CalMz <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$SiteName, ages = box2$Age, h = 100)
#write.table(boxbins, file = "data/NorFranceboxbins.csv", sep = ",", col.names=NA)
####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(15000,2000))
plot(spd.mz, runm=200, xlim=c(15000,2000), type="simple")

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(15000,2000),bw=50, normalised = FALSE)
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

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<- subset(dd, MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/SouthBritKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/SouthBritKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)

calBP<-c(14970,14940,14910,14880,14850,14820,14790,14760,14730,14700,14670,14640,14610,14580,14550,14520,14490,14460,14430,14400,14370,14340,14310,14280,14250,14220,14190,14160,14130,14100,14070,14040,14010,13980,13950,13920,13890,13860,13830,
         13800,13770,13740,13710,13680,13650,13620,13590,13560,13530,13500,13470,13440,13410,13380,13350,13320,13290,13260,13230,13200,13170,13140,13110,13080,13050,13020,12990,12960,12930,12900,12870,12840,12810,12780,12750,12720,12690,12660,12630,12600,12570,12540,12510,12480,12450,12420,12390,12360,12330,12300,12270,
         12240,12210,12180, 12150,12120,12090,12060,12030,12000,11970,11940,11910,11880,11850,11820,11790,11760,11730,11700,11670,11640,11610,11580,11550,11520,11490,11460,11430,11400,11370,11340,11310,11280,11250,11220,11190,11160,11130,11100,11070,11040,11010,10980,10950,10920,10890,10860,10830,10800,10770,10740,10710,10680,10650,10620,10590,10560,10530,10500,10470,10440,10410,10380,10350,10320,10290,10260,10230,10200,10170,10140,
         10110,10080,10050,10020,9990,9960,9930,9900,9870,9840,9810,9780,9750,9720,9690,9660,9630,9600,9570,9540,9510,9480,9450,9420,9390,9360,9330,9300,9270,9240,9210,9180,9150,9120,9090,9060,9030,9000,8970,8940,8910,8880,8850,8820,
         8790,8760,8730,8700,8670,8640,8610,8580,8550,8520,8490,8460,8430,8400,8370,8340,8310,8280,8250,8220,8190,8160,8130,8100,8070,8040,8010,7980,7950,7920,7890,7860,7830,7800,7770,7740,7710,7680,7650,7620,7590,7560,7530,7500,7470,7440,7410,7380,7350,7320,7290,7260,7230,
         7200,7170,7140,7110,7080,7050,7020,6990,6960,6930,6900,6870,6840,6810,6780,6750,6720,6690,6660,6630,6600,6570,6540,6510,6480,6450,6420,6390,6360,
         6330,6300,6270,6240,6210,6180,6150,6120,6090,6060,6030,6000,5970,5940,5910,5880,5850,5820,5790,5760,5730,5700,5670,5640,5610,5580,5550,5520,5490,5460,5430,5400,5370,5340,5310,5280,5250,5220,5190,5160,5130,5100,5070,5040,5010,4980,4950,4920,4890,4860,4830,4800,4770,4740,4710,4680,4650,4620,4590,4560,4530,4500,
         4470,4440,4410,4380,4350,4320,4290,4260,4230,4200,4170,4140,4110,4080,4050,4020,3990,3960,3930,3900,3870,
         3840,3810,3780,3750,3720,3690,3660,3630,3600, 3570,3540,3510,3480,3450,3420,3390,3360,3330,3300,3270,3240,3210,3180,3150,3120,3090,3060,3030,3000,2970,2940,2910,2880,2850,2820,2790,2760,2730,2700,2670,2640,2610,2580,2550,2520,2490,2460,2430,
         2400,2370,2340,2310,2280,2250,2220,2190,2160,2130,2100,2070,2040,2010)

sums<-cbind(calBP, MKDE, out50)


write.table(sums, file = "data/Sumbin/SouthBritSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read.csv("data/Sumbin/SouthBritSumbin.csv") 
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(2000, 4200, 6060, 11500, 15000),
                         labels=c('Metal Age','Neolithic','Mesolithic','Paleolithic'))

write.table(pcgrowth, file = "data/Percapita/SouthBritPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
gr30pc<- read.csv("data/Percapita/SouthBritPerCap.csv")

gr30pc2<-subset(gr30pc, calBP<7001 & calBP>4499)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(gr30pc2$MKDE-min(gr30pc2$MKDE))/(max((gr30pc2$MKDE)-min(gr30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
gr30pc3<-cbind(StKDE,gr30pc2)

pcGr <- ggplot(gr30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.05,0.17))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A. S. Britain KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept = 0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcGr

Grcpt <- ggplot(gr30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  #scale_y_continuous(limits=c(-.05,0.17))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. S. Britain KDE density vs. Time")+
  geom_hline(yintercept = 0.20)
#geom_vline(xintercept = 2550)+
#geom_vline(xintercept = 2250)+
#geom_vline(xintercept = 830)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
Grcpt

#paired plots
Figsbrit<-plot_grid(pcGr, Grcpt, ncol=2, align="hv", axis = "rl")
Figsbrit

pdf("data/ExSBritain.pdf", width=20.55, height=14)
Figsbrit
dev.off()

####North Sea Basin Case ID=============================================
box<- read.csv("data/RawP3Kc14.csv")
#box2<- subset(box, Latitude>55 & Latitude<60 & Longitude>4 & Longitude<14)

#write.table(box2, file = "data/NorthSeadates.csv", sep = ",", col.names=NA)
###North Sea dates cleaned to remove underscores
box2<-read.csv(file="data/NorthSeadates.csv", header=T)

# Load world map data
#world <- ne_countries(scale = "medium", returnclass = "sf")
europe <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")

# Define the list of Middle Eastern countries
eu_countries <- c("Spain", "France", "Italy", "Greece", "Albania","Denmark", "Germany","Belgium", "Sweden", "United Kingdom", "Croatia", "Slovinia", "Norway", "Netherlands")

# Filter the data for the Middle East
eu <- europe[europe$name %in% eu_countries, ]

ggplot(data = eu) +
  geom_sf(fill = "NA", color = "black") +
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Country)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)+
  #geom_map(map=middle_east)+
  # coord_map("moll")+
  #coord_sf(crs = "+proj=aeqd +lat_0=30 +lon_0=30") +
  theme_minimal() +
  labs(title = "Map of the Middle East",
       subtitle = "Countries in the Middle East",
       x = "Longitude",
       y = "Latitude")

box2 <- box2[!is.na(box2$SiteName), ]

CalMz <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$SiteName, ages = box2$Age, h = 100)
#write.table(boxbins, file = "data/NorFranceboxbins.csv", sep = ",", col.names=NA)
####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(15000,2000))
plot(spd.mz, runm=200, xlim=c(15000,2000), type="simple")

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(15000,2000),bw=50, normalised = FALSE)
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

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<- subset(dd, MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/NorthSeaKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/NorthSeaKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(14970,14940,14910,14880,14850,14820,14790,14760,14730,14700,14670,14640,14610,14580,14550,14520,14490,14460,14430,14400,14370,14340,14310,14280,14250,14220,14190,14160,14130,14100,14070,14040,14010,13980,13950,13920,13890,13860,13830,
         13800,13770,13740,13710,13680,13650,13620,13590,13560,13530,13500,13470,13440,13410,13380,13350,13320,13290,13260,13230,13200,13170,13140,13110,13080,13050,13020,12990,12960,12930,12900,12870,12840,12810,12780,12750,12720,12690,12660,12630,12600,12570,12540,12510,12480,12450,12420,12390,12360,12330,12300,12270,
         12240,12210,12180, 12150,12120,12090,12060,12030,12000,11970,11940,11910,11880,11850,11820,11790,11760,11730,11700,11670,11640,11610,11580,11550,11520,11490,11460,11430,11400,11370,11340,11310,11280,11250,11220,11190,11160,11130,11100,11070,11040,11010,10980,10950,10920,10890,10860,10830,10800,10770,10740,10710,10680,10650,10620,10590,10560,10530,10500,10470,10440,10410,10380,10350,10320,10290,10260,10230,10200,10170,10140,
         10110,10080,10050,10020,9990,9960,9930,9900,9870,9840,9810,9780,9750,9720,9690,9660,9630,9600,9570,9540,9510,9480,9450,9420,9390,9360,9330,9300,9270,9240,9210,9180,9150,9120,9090,9060,9030,9000,8970,8940,8910,8880,8850,8820,
         8790,8760,8730,8700,8670,8640,8610,8580,8550,8520,8490,8460,8430,8400,8370,8340,8310,8280,8250,8220,8190,8160,8130,8100,8070,8040,8010,7980,7950,7920,7890,7860,7830,7800,7770,7740,7710,7680,7650,7620,7590,7560,7530,7500,7470,7440,7410,7380,7350,7320,7290,7260,7230,
         7200,7170,7140,7110,7080,7050,7020,6990,6960,6930,6900,6870,6840,6810,6780,6750,6720,6690,6660,6630,6600,6570,6540,6510,6480,6450,6420,6390,6360,
         6330,6300,6270,6240,6210,6180,6150,6120,6090,6060,6030,6000,5970,5940,5910,5880,5850,5820,5790,5760,5730,5700,5670,5640,5610,5580,5550,5520,5490,5460,5430,5400,5370,5340,5310,5280,5250,5220,5190,5160,5130,5100,5070,5040,5010,4980,4950,4920,4890,4860,4830,4800,4770,4740,4710,4680,4650,4620,4590,4560,4530,4500,
         4470,4440,4410,4380,4350,4320,4290,4260,4230,4200,4170,4140,4110,4080,4050,4020,3990,3960,3930,3900,3870,
         3840,3810,3780,3750,3720,3690,3660,3630,3600, 3570,3540,3510,3480,3450,3420,3390,3360,3330,3300,3270,3240,3210,3180,3150,3120,3090,3060,3030,3000,2970,2940,2910,2880,2850,2820,2790,2760,2730,2700,2670,2640,2610,2580,2550,2520,2490,2460,2430,
         2400,2370,2340,2310,2280,2250,2220,2190,2160,2130,2100,2070,2040,2010)

sums<-cbind(calBP, MKDE, out50)

write.table(sums, file = "data/Sumbin/NorthSeaSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read.csv("data/Sumbin/NorthSeaSumbin.csv") 
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(2000, 3630, 6060, 11500, 15000),
                         labels=c('Metal Age','Neolithic','Mesolithic','Paleolithic'))

write.table(pcgrowth, file = "data/Percapita/NorthSeaPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
ns30pc<- read.csv("data/Percapita/NorthSeaPerCap.csv")

ns30pc2<-subset(ns30pc, calBP<6800 & calBP>4300)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(ns30pc2$MKDE-min(ns30pc2$MKDE))/(max((ns30pc2$MKDE)-min(ns30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
ns30pc3<-cbind(StKDE,ns30pc2)

pcns <- ggplot(ns30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.1,0.1))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A.North Sea Coast KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcns

nscpt <- ggplot(ns30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. North Sea Coast Density vs. Time")+
  geom_hline(yintercept=.2)
#geom_vline(xintercept = 3320)+
#geom_vline(xintercept = 2550)+
#geom_vline(xintercept = 2250)+
#geom_vline(xintercept = 830)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
nscpt

##Plot example figs
Figns<-plot_grid(pcns, nscpt, ncol=2, align="hv", axis = "rl")
Figns

pdf("data/Figs/ExNSea.pdf", width=20.55, height=14)
Figns
dev.off()

####Central Texas Cycles===========================================================
box<- read.csv("data/FinalRCDTexas3.csv")
box2<- subset(box, Region=="CTx")

counties<-map_data("state")

ArchGlobeMap<-ggplot() +
  geom_polygon(data = counties, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  #geom_polygon(data = canada, aes(x=long, y = lat, group = group),
  #    fill = "white", color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Longitude", y="Latitude", title = "US ArchaeoGlobe Regions and Radiocarbon")+
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Region)),
             inherit.aes = FALSE, alpha = 0.5, size = 2, shape=23)
ArchGlobeMap

cptcal <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$Trinomial, ages = box2$Age, h = 100)

####Run analysis for component 3 logistic 3500 to 150
spd.CTx <- spd(cptcal, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.CTx, runm=200, xlim=c(8200,200), type="simple")

PrDens<-spd.CTx$grid$PrDens
calBP<-spd.CTx$grid$calBP

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(cptcal, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
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

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/CTexKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/CTexKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)

write.table(sums, file = "data/Sumbin/CTexSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("data/Sumbin/CTexSumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 500, 790, 1210, 4200, 8200),
                         labels=c('Historic','Toya','Austin','Late Archaic','Middle Archaic'))

write.table(pcgrowth, file = "data/Percapita/CTexPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
ct30pc<- read.csv("data/Percapita/CTexPerCap.csv")

ct30pc2<-subset(ct30pc, calBP<4000 & calBP>500)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(ct30pc2$MKDE-min(ct30pc2$MKDE))/(max((ct30pc2$MKDE)-min(ct30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
ct30pc3<-cbind(StKDE,ct30pc2)

pcctex <- ggplot(ct30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.1,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Mean KDE (density)", y="KDE per capita growth", title = "A. Central Texas Per Capita Growth vs. Density")+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcctex

ctexcpt <- ggplot(ct30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=factor(PeriodID)), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE Density", title = "B. Central Texas Standardized Density vs. Time")+
  geom_hline(yintercept=0.25)
ctexcpt

#paired plot===============================
Figctex<-plot_grid(pcctex, ctexcpt, ncol=2, align="hv", axis = "rl")
Figctex

pdf("data/Figs/Exctex.pdf", width=20.55, height=14)
Figctex
dev.off()

####Texas Coastal Plain Case ID 2===========================================================
box<- read.csv("data/FinalRCDTexas3.csv")
box2<- subset(box, Region=="TCP")

counties<-map_data("state")

ArchGlobeMap<-ggplot() +
  geom_polygon(data = counties, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  #geom_polygon(data = canada, aes(x=long, y = lat, group = group),
  #    fill = "white", color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Longitude", y="Latitude", title = "US ArchaeoGlobe Regions and Radiocarbon")+
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Region)),
             inherit.aes = FALSE, alpha = 0.5, size = 2, shape=23)
ArchGlobeMap

cptcal <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$Trinomial, ages = box2$Age, h = 100)

####Run analysis for component 3 logistic 3500 to 150
spd.CTx <- spd(cptcal, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.CTx, runm=200, xlim=c(8200,200), type="simple")

PrDens<-spd.CTx$grid$PrDens
calBP<-spd.CTx$grid$calBP

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(cptcal, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
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

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/TCPKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/TCPKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/TCPSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("data/Sumbin/TCPSumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 500, 790, 1210, 4000, 8200),
                         labels=c('Historic','Toya','Austin','Late Archaic','Middle Archaic'))
write.table(pcgrowth, file = "data/Percapita/TCPPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
txc30pc<- read.csv("data/Percapita/TCPPerCap.csv")

txc30pc2<-subset(txc30pc, calBP<5000 & calBP>500)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(txc30pc2$MKDE-min(txc30pc2$MKDE))/(max((txc30pc2$MKDE)-min(txc30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
txc30pc3<-cbind(StKDE, txc30pc2)

pctxcoast <- ggplot(txc30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=factor(PeriodID)), size=3.5) +
  geom_path(aes(),linewidth=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.15,0.17))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Mean KDE (density)", y="KDE per capita growth", title = "A. TCP KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pctxcoast

txcoastcpt <- ggplot(txc30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=factor(PeriodID)), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  #scale_x_reverse()+
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Mean KDE", title = "B. TCP Standardized Density vs. Time")+
  geom_hline(yintercept=0.25)
#geom_vline(xintercept = 770)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
txcoastcpt

#paired plot===============================
Figtcp<-plot_grid(pctxcoast, txcoastcpt, ncol=2, align="hv", axis = "rl")
Figtcp

pdf("data/ExTexCoast.pdf", width=20.55, height=14)
Figtcp
dev.off()

####SE US Case ID 14===========================================================
SPD<-read.csv(file="NERDv4_0/RawP3Kc14.csv", header=T)
box2<- subset(SPD, Latitude>30 & Latitude<35 & Longitude>-87.51 & Longitude< -82.5)
#write.table(box2, file = "data/Southeastdates.csv", sep = ",", col.names=NA)
#box2<-read.csv(file="data/Southeastdates.csv", header=T)
counties<-map_data("state")

ArchGlobeMap<-ggplot() +
  geom_polygon(data = counties, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  #geom_polygon(data = canada, aes(x=long, y = lat, group = group),
  #    fill = "white", color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Longitude", y="Latitude", title = "US ArchaeoGlobe Regions and Radiocarbon")+
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Province)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)
ArchGlobeMap

box2 <- box2[!is.na(box2$SiteID), ]

cptcal <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$SiteID, ages = box2$Age, h = 100)

####Run analysis for component 3 logistic 3500 to 150
spd.CTx <- spd(cptcal, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.CTx, runm=200, xlim=c(8200,200), type="simple")

PrDens<-spd.CTx$grid$PrDens
calBP<-spd.CTx$grid$calBP

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(cptcal, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
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

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/SoutheastKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/SoutheastKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/SoutheastSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("data/Sumbin/SoutheastSumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
##Add ID variable for culture history periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 580, 1060, 1570, 2050, 3000, 6000, 8200),
                         labels=c('Post-Mississippian', 'Mississippian', 'Late Woodland', 'Middle Woodland','Early Woodland','Late Archic','Middle Archaic'))
write.table(pcgrowth, file = "data/Percapita/SoutheastPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
se30pc<- read.csv("data/Percapita/SoutheastPerCap.csv")

se30pc2<-subset(se30pc, calBP<4900 & calBP>400)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(se30pc2$MKDE-min(se30pc2$MKDE))/(max((se30pc2$MKDE)-min(se30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
se30pc3<-cbind(StKDE, se30pc2)

pcseus <- ggplot(se30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),linewidth=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE Density", y="KDE per capita growth", title = "A. Southeast US KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept=0.0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcseus

seuscpt <- ggplot(se30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE Density", title = "B. Southeast US Density vs. Time")+
  geom_hline(yintercept=0.2)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
seuscpt

#Paired plots==================================
FigSEUS<-plot_grid(pcseus, seuscpt, ncol=2, align="hv", axis = "rl")
FigSEUS

pdf("data/ExSEus.pdf", width=20.55, height=14)
FigSEUS
dev.off()


####NE US Case ID 8===========================================================
SPD<-read.csv(file="data/RawP3Kc14.csv", header=T)
box2<- subset(SPD, Latitude>41 & Latitude<46 & Longitude>-81 & Longitude< -74)
#write.table(box2, file = "data/Northeastdates.csv", sep = ",", col.names=NA)
#box2<-read.csv(file="data/Northeastdates.csv", header=T)
counties<-map_data("state")

ArchGlobeMap<-ggplot() +
  geom_polygon(data = counties, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  #geom_polygon(data = canada, aes(x=long, y = lat, group = group),
  #    fill = "white", color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Longitude", y="Latitude", title = "US ArchaeoGlobe Regions and Radiocarbon")+
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Province)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)
ArchGlobeMap

box2 <- box2[!is.na(box2$SiteName), ]

cptcal <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$SiteName, ages = box2$Age, h = 100)

####Run analysis for component 3 logistic 3500 to 150
spd.CTx <- spd(cptcal, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.CTx, runm=200, xlim=c(8200,200), type="simple")

PrDens<-spd.CTx$grid$PrDens
calBP<-spd.CTx$grid$calBP

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(cptcal, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
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

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/NortheastKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/NortheastKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........

# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/NortheastSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("data/Sumbin/NortheastSumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

##Add ID variable for culture history periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 400, 1570, 2700, 6000, 8200),
                         labels=c('Historic' ,'Late', 'Woodland','Late Archaic','Middle Archaic'))
write.table(pcgrowth, file = "data/Percapita/NortheastPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
ne30pc<- read.csv("data/Percapita/NortheastPerCap.csv")

ne30pc2<-subset(ne30pc, calBP<2401 & calBP>299)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(ne30pc2$MKDE-min(ne30pc2$MKDE))/(max((ne30pc2$MKDE)-min(ne30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
ne30pc3<-cbind(StKDE, ne30pc2)


pcne <- ggplot(ne30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),linewidth=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A. Northeast US KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcne

necpt <- ggplot(ne30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. Northeast US Density vs. Time")+
  geom_vline(xintercept = 770)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
necpt

#paired plots
FigneUS<-plot_grid(pcne, necpt, ncol=2, align="hv", axis = "rl")
FigneUS

pdf("data/Figs/ExNortheastUS.pdf", width=20.55, height=14)
FigneUS
dev.off()


####NW Coast===========================================================
SPD<-read.csv(file="NERDv4_0/RawP3Kc14.csv", header=T)
box2<- subset(SPD, Latitude>42 & Latitude<52 & Longitude>-125 & Longitude< -121)
#write.table(box2, file = "data/NWCoastdates.csv", sep = ",", col.names=NA)
#box2<-read.csv(file="data/NWCoastdates.csv", header=T)
counties<-map_data("state")

ArchGlobeMap<-ggplot() +
  geom_polygon(data = counties, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  #geom_polygon(data = canada, aes(x=long, y = lat, group = group),
  #    fill = "white", color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Longitude", y="Latitude", title = "US ArchaeoGlobe Regions and Radiocarbon")+
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Province)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)
ArchGlobeMap

box2 <- box2[!is.na(box2$SiteName), ]

cptcal <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$SiteName, ages = box2$Age, h = 100)

####Run analysis for component 3 logistic 3500 to 150
spd.CTx <- spd(cptcal, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.CTx, runm=200, xlim=c(8200,200), type="simple")

PrDens<-spd.CTx$grid$PrDens
calBP<-spd.CTx$grid$calBP

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(cptcal, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
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

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/NWCoastKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/NWCoastKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/NWCoastSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("data/Sumbin/NWCoastSumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

##Add ID variable for culture history periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 325, 1750, 3750, 6350, 8200),
                         labels=c('Historic' ,'Late Pacific', 'Middle Pacific','Early Pacific','Archaic'))
write.table(pcgrowth, file = "data/Percapita/NWCoastPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
nwc30pc<- read.csv("data/Percapita/NWCoastPerCap.csv")

nwc30pc2<-subset(nwc30pc, calBP<5000 & calBP>500)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(nwc30pc2$MKDE-min(nwc30pc2$MKDE))/(max((nwc30pc2$MKDE)-min(nwc30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
nwc30pc3<-cbind(StKDE, nwc30pc2)

pcnwc <- ggplot(nwc30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),linewidth=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A. US NW Coast KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept = 0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcnwc

nwccpt <- ggplot(nwc30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. US NW Coast Density vs. Time")+
  geom_hline(yintercept = 0.2)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
nwccpt

#paired plots

Fignwc<-plot_grid(pcnwc, nwccpt, ncol=2, align="hv", axis = "rl")
Fignwc

pdf("data/Exnwc.pdf", width=20.55, height=14)
Fignwc
dev.off()

####SW Wyoming===========================================================
SPD<-read.csv(file="data/RawP3Kc14.csv", header=T)
box2<- subset(SPD, Latitude>41 & Latitude<44 & Longitude>-111 & Longitude< -106)
#write.table(box2, file = "data/SWWydates.csv", sep = ",", col.names=NA)
#box2<-read.csv(file="data/SWWydates.csv", header=T)
counties<-map_data("state")

ArchGlobeMap<-ggplot() +
  geom_polygon(data = counties, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  #geom_polygon(data = canada, aes(x=long, y = lat, group = group),
  #    fill = "white", color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Longitude", y="Latitude", title = "US ArchaeoGlobe Regions and Radiocarbon")+
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Province)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)
ArchGlobeMap

box2 <- box2[!is.na(box2$SiteID), ]

cptcal <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$SiteID, ages = box2$Age, h = 100)

####Run analysis for component 3 logistic 3500 to 150
spd.CTx <- spd(cptcal, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.CTx, runm=200, xlim=c(8200,200), type="simple")

PrDens<-spd.CTx$grid$PrDens
calBP<-spd.CTx$grid$calBP

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(cptcal, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
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

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/SWWyKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/SWWyKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/SWWySumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("data/Sumbin/SWWySumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

##Add ID variable for culture history periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 400, 1800, 3670, 8200),
                         labels=c('Historic', 'Late Prehistoric','Late Archaic','Early Archaic'))
write.table(pcgrowth, file = "data/Percapita/SWWyPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
wy30pc<- read.csv("data/Percapita/SWWyPerCap.csv")

wy30pc2<-subset(wy30pc, calBP<8000 & calBP>3500)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(wy30pc2$MKDE-min(wy30pc2$MKDE))/(max((wy30pc2$MKDE)-min(wy30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
wy30pc3<-cbind(StKDE,wy30pc2)

pcswwy <- ggplot(wy30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),linewidth=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.1,0.15))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A. SW Wyoming KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcswwy

swwycpt <- ggplot(wy30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. SW Wyoming Density vs. Time")+
  geom_hline(yintercept=.20)
#geom_vline(xintercept = 770)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
swwycpt

##Paired Plots=================

Figswwy<-plot_grid(pcswwy, swwycpt, ncol=2, align="hv", axis = "rl")
Figswwy

pdf("data/Exswwy.pdf", width=20.55, height=14)
Figswwy
dev.off()


#### East S. Africa==================================================================
SPD<-read.csv(file="data/RawP3Kc14.csv", header=T)
box2<- subset(SPD, Latitude>-35 & Latitude< -23 & Longitude> 25 & Longitude<35 )
#write.table(box2, file = "data/EastSAfricadates.csv", sep = ",", col.names=NA)
#box2<-read.csv(file="data/SWWydates.csv", header=T)

world <- ne_countries(scale = "medium", returnclass = "sf")

# Filter for South Africa
south_africa <- world[world$name == "South Africa", ]

ggplot(data = south_africa) +
  geom_sf(fill = "lightblue", color = "black") +
  labs(title = "Map of South Africa",
       subtitle = "Created with ggplot2 and rnaturalearth") +
  theme_minimal()+
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Province)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)

box2 <- box2[!is.na(box2$SiteName), ]

cptcal <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "shcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$SiteName, ages = box2$Age, h = 100)

####Run analysis for component 3 logistic 3500 to 150
spd.CTx <- spd(cptcal, bins=boxbins, runm=200, timeRange=c(8200,100))
plot(spd.CTx, runm=200, xlim=c(8200,100), type="simple")

PrDens<-spd.CTx$grid$PrDens
calBP<-spd.CTx$grid$calBP

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(cptcal, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,100),bw=50, normalised = FALSE)
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

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/ESouthAfKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/ESouthAfKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/ESouthAfSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("data/Sumbin/ESouthAfSumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
##Add ID variable for culture history periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 700, 2110, 8200),
                         labels=c('Polity','Agg','HG'))
write.table(pcgrowth, file = "data/Percapita/ESouthAfPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
esa30pc<- read.csv("data/Percapita/ESouthAfPerCap.csv")

esa30pc2<-subset(esa30pc, calBP<2101 & calBP>600)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(esa30pc2$MKDE-min(esa30pc2$MKDE))/(max((esa30pc2$MKDE)-min(esa30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
esa30pc3<-cbind(StKDE, esa30pc2)

pcesa <- ggplot(esa30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),linewidth=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A. East S. Africa KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcesa

esacpt <- ggplot(esa30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. East S. Africa Density vs. Time")+
  geom_hline(yintercept=.2)
#geom_vline(xintercept = 770)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
esacpt

#paired Plots
Figesa<-plot_grid(pcesa, esacpt, ncol=2, align="hv", axis = "rl")
Figesa

pdf("data/Exesa.pdf", width=20.55, height=14)
Figesa
dev.off()

#### West S. Africa==================================================================
SPD<-read.csv(file="data/RawP3Kc14.csv", header=T)
box2<- subset(SPD, Latitude>-35 & Latitude< -23 & Longitude> 16 & Longitude<25 )
#write.table(box2, file = "data/EastSAfricadates.csv", sep = ",", col.names=NA)
#box2<-read.csv(file="data/SWWydates.csv", header=T)

world <- ne_countries(scale = "medium", returnclass = "sf")

# Filter for South Africa
south_africa <- world[world$name == "South Africa", ]

ggplot(data = south_africa) +
  geom_sf(fill = "lightblue", color = "black") +
  labs(title = "Map of South Africa",
       subtitle = "Created with ggplot2 and rnaturalearth") +
  theme_minimal()+
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Province)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)

box2 <- box2[!is.na(box2$SiteName), ]

cptcal <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "shcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$SiteName, ages = box2$Age, h = 100)

####Run analysis for component 3 logistic 3500 to 150
spd.CTx <- spd(cptcal, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.CTx, runm=200, xlim=c(8200,200), type="simple")

PrDens<-spd.CTx$grid$PrDens
calBP<-spd.CTx$grid$calBP

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(cptcal, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
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

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/WSouthAfKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/WSouthAfKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/WSouthAfSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("data/Sumbin/WSouthAfSumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 700, 2110, 8200),
                         labels=c('Polity','Herding','HG'))
write.table(pcgrowth, file = "data/Percapita/WSouthAfPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
wsa30pc<- read.csv("data/Percapita/WSouthAfPerCap.csv")

wsa30pc2<-subset(wsa30pc, calBP<2700 & calBP>299)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(wsa30pc2$MKDE-min(wsa30pc2$MKDE))/(max((wsa30pc2$MKDE)-min(wsa30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
wsa30pc3<-cbind(StKDE, wsa30pc2)

pcwsa <- ggplot(wsa30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),linewidth=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A. West S. Africa KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcwsa

wsacpt <- ggplot(wsa30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. West S. Africa Density vs. Time")+
  geom_hline(yintercept=.2)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
wsacpt

##Paired plots

Figwsa<-plot_grid(pcwsa, wsacpt, ncol=2, align="hv", axis = "rl")
Figwsa

pdf("data/Exwsa.pdf", width=20.55, height=14)
Figwsa
dev.off()

#### NW Australia==================================================================
SPD<-read.csv(file="data/RawP3Kc14.csv", header=T)
box2<- subset(SPD, Latitude>-23 & Latitude< -10 & Longitude>125 & Longitude<140 )

australia_map <- map_data("world", region = "Australia")

# Create the map
ggplot(data = australia_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "blue", color = "black") +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Map of Australia")+
  geom_point(data=box2, aes(Longitude, Latitude),
             inherit.aes = FALSE, alpha = 0.5, size = 2)

box2 <- box2[!is.na(box2$SiteName), ]


cptcal <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "shcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$SiteName, ages = box2$Age, h = 100)

spd.CTx <- spd(cptcal, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.CTx, runm=200, xlim=c(8200,200), type="simple")

PrDens<-spd.CTx$grid$PrDens
calBP<-spd.CTx$grid$calBP

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(cptcal, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
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

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/NWAustKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/NWAustKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/NWAustSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("data/Sumbin/NWAustSumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 500, 2050, 8200),
                         labels=c('Late','Intensive HG','HG'))
write.table(pcgrowth, file = "data/Percapita/NWAustPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
nwau30pc<- read.csv("data/Percapita/NWAustPerCap.csv")

nwau30pc2<-subset(nwau30pc, calBP<5000 & calBP>500)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(nwau30pc2$MKDE-min(nwau30pc2$MKDE))/(max((nwau30pc2$MKDE)-min(nwau30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
nwau30pc3<-cbind(StKDE, nwau30pc2)

pcnwaus <- ggplot(nwau30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),linewidth=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.17,0.2))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A. NW Australia KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcnwaus

nwauscpt <- ggplot(nwau30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.15,0.17))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. NW Australia Density vs. Time")+
  geom_hline(yintercept=.2)+
  geom_vline(xintercept = 1000)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
nwauscpt

#paired Plots
Fignwaus<-plot_grid(pcnwaus, nwauscpt, ncol=2, align="hv", axis = "rl")
Fignwaus

pdf("data/Exnwestaus.pdf", width=20.55, height=14)
Fignwaus
dev.off()

#### SE Australia==================================================================
SPD<-read.csv(file="data/RawP3Kc14.csv", header=T)
#box2<- subset(SPD, Latitude>-40 & Latitude< -30 & Longitude>140 & Longitude<155 )
#write.table(box2, file = "data/SEAustdates.csv", sep = ",", col.names=NA)
box2<-read.csv(file="data/SEAustdates.csv", header=T)

australia_map <- map_data("world", region = "Australia")

# Create the map
ggplot(data = australia_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "blue", color = "black") +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Map of Australia")+
  geom_point(data=box2, aes(Longitude, Latitude),
             inherit.aes = FALSE, alpha = 0.5, size = 2)

box2 <- box2[!is.na(box2$SiteName), ]


cptcal <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "shcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$SiteName, ages = box2$Age, h = 100)

####Run analysis for component 3 logistic 3500 to 150
spd.CTx <- spd(cptcal, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.CTx, runm=200, xlim=c(8200,200), type="simple")

PrDens<-spd.CTx$grid$PrDens
calBP<-spd.CTx$grid$calBP

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(cptcal, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
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

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/SEAustKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/SEAustKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/SEAustSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("data/Sumbin/SEAustSumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 400, 2500, 8200),
                         labels=c('Late','Intensive HG','HG'))
write.table(pcgrowth, file = "data/Percapita/SEAustPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
seau30pc<- read.csv("data/Percapita/SEAustPerCap.csv")

seau30pc2<-subset(seau30pc, calBP<5000 & calBP>500)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(seau30pc2$MKDE-min(seau30pc2$MKDE))/(max((seau30pc2$MKDE)-min(seau30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
seau30pc3<-cbind(StKDE, seau30pc2)

pcseaus <- ggplot(seau30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),linewidth=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.07,0.15))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", 
       title = "A. SE Australia KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcseaus

seauscpt <- ggplot(seau30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. SE Australia Desnity vs. Time")+
  geom_vline(xintercept = 770)+
  geom_hline(yintercept=.2)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
seauscpt

Figseaus<-plot_grid(pcseaus, seauscpt, ncol=2, align="hv", axis = "rl")
Figseaus

pdf("data/Exseaus.pdf", width=20.55, height=14)
Figseaus
dev.off()

####Japan South ===================================================================
SPD2<-read.csv(file="data/JapanDates.csv", header=T)
box2<- subset(SPD2, Latitude>32 & Latitude< 36 & Longitude>129 & Longitude<134 )
#write.table(box2, file = "data/SEAustdates.csv", sep = ",", col.names=NA)
#box2<-read.csv(file="data/SEAustdates.csv", header=T)

japan_map <- map_data("world", region = "Japan")

# Create the map using ggplot2
ggplot(data = japan_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "lightblue", color = "black") +
  coord_fixed(1.3) +
  labs(title = "Map of Japan",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal()+
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Prefecture)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)

box2 <- box2[!is.na(box2$SiteNameEn), ]

cptcal <- calibrate(x = box2$CRA,  errors = box2$CRAError, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$SiteNameEn, ages = box2$CRA, h = 100)

####Run analysis for component 3 logistic 3500 to 150
spd.CTx <- spd(cptcal, bins=boxbins, runm=200, timeRange=c(15000,200))
plot(spd.CTx, runm=200, xlim=c(15000,200), type="simple")

PrDens<-spd.CTx$grid$PrDens
calBP<-spd.CTx$grid$calBP

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(cptcal, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(15000,200),bw=50, normalised = FALSE)
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

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/JapanSouthKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs//JapanSouthKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(14970,14940,14910,14880,14850,14820,14790,14760,14730,14700,14670,14640,14610,14580,14550,14520,14490,14460,14430,14400,14370,14340,14310,14280,14250,14220,14190,14160,14130,
         14100,14070,14040,14010,13980,13950,13920,13890,13860,13830,13800,13770,13740,13710,13680,13650,13620,13590,13560,13530,13500,13470,13440,13410,13380,13350,13320,13290,13260,
         13230,13200,13170,13140,13110,13080,13050,13020,12990,12960,12930,12900,12870,12840,12810,12780,12750,12720,12690,12660,12630,12600,12570,12540,12510,12480,12450,12420,12390,
         12360,12330,12300,12270,12240,12210,12180,12150,12120,12090,12060,12030,12000,11970,11940,11910,11880,11850,11820,11790,11760,11730,11700,11670,11640,11610,11580,11550,11520,
         11490,11460,11430,11400,11370,11340,11310,11280,11250,11220,11190,11160,11130,11100,11070,11040,11010,10980,10950,10920,10890,10860,10830,10800,10770,10740,10710,10680,10650,
         10620,10590,10560,10530,10500,10470,10440,10410,10380,10350,10320,10290,10260,10230,10200,10170,10140,10110,10080,10050,10020,9990,9960,9930,9900,9870,9840,9810,9780,9750,9720,
         9690,9660,9630,9600,9570,9540,9510,9480,9450,9420,9390,9360,9330,9300,9270,9240,9210,9180,9150,9120,9090,9060,9030,9000,8970,8940,8910,8880,8850,8820,8790,8760,8730,8700,8670,
         8640,8610,8580,8550,8520,8490,8460,8430,8400,8370,8340,8310,8280,8250,8220,8190,8160,8130,8100,8070,8040,8010,7980,7950,7920,7890,7860,7830,7800,7770,7740,7710,7680,7650,7620,
         7590,7560,7530,7500,7470,7440,7410,7380,7350,7320,7290,7260,7230,7200,7170,7140,7110,7080,7050,7020,6990,6960,6930,6900,6870,6840,6810,6780,6750,6720,6690,6660,6630,6600,6570,
         6540,6510,6480,6450,6420,6390,6360,6330,6300,6270,6240,6210,6180,6150,6120,6090,6060,6030,6000,5970,5940,5910,5880,5850,5820,5790,5760,5730,5700,5670,5640,5610,5580,5550,5520,
         5490,5460,5430,5400,5370,5340,5310,5280,5250,5220,5190,5160,5130,5100,5070,5040,5010,4980,4950,4920,4890,4860,4830,4800,4770,4740,4710,4680,4650,4620,4590,4560,4530,4500,4470,
         4440,4410,4380,4350,4320,4290,4260,4230,4200,4170,4140,4110,4080,4050,4020,3990,3960,3930,3900,3870,3840,3810,3780,3750,3720,3690,3660,3630,3600,3570,3540,3510,3480,3450,3420,
         3390,3360,3330,3300,3270,3240,3210,3180,3150,3120,3090,3060,3030,3000,2970,2940,2910,2880,2850,2820,2790,2760,2730,2700,2670,2640,2610,2580,2550,2520,2490,2460,2430,2400,2370,
         2340,2310,2280,2250,2220,2190,2160,2130,2100,2070,2040,2010,1980,1950,1920,1890,1860,1830,1800,1770,1740,1710,1680,1650,1620,1590,1560,1530,1500,1470,1440,1410,1380,1350,1320,
         1290,1260,1230,1200,1170,1140,1110,1080,1050,1020,990,960,930,900,870,840,810,780,750,720,690,660,630,600,570,540,510,480,450,420,390,360,330,300,270,240,210)
sums<-cbind(calBP, MKDE, out50)

write.table(sums, file = "data/Sumbin/JapanSouthSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("data/Sumbin/JapanSouthSumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 1350, 1590, 1770, 3100, 11500, 15000),
                         labels=c('Classic','Kofun','Shonai', 'Yayoi','Jamon','Paleolithic'))
write.table(pcgrowth, file = "data/Percapita/JapanSouthPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
sjp30pc<- read.csv("data/Percapita/JapanSouthPerCap.csv")

sjp30pc2<-subset(sjp30pc, calBP<3001 & calBP>999)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(sjp30pc2$MKDE-min(sjp30pc2$MKDE))/(max((sjp30pc2$MKDE)-min(sjp30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
sjp30pc3<-cbind(StKDE, sjp30pc2)

sjppc <- ggplot(sjp30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),linewidth=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.05,0.15))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth",
       title = "A. S. Japan KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
sjppc

sjpcpt <- ggplot(sjp30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. S. Japan Density vs. Time")+
  geom_hline(yintercept=.2)
#geom_vline(xintercept = 770)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
sjpcpt

#paired Plots

Figsjp<-plot_grid(sjppc, sjpcpt, ncol=2, align="hv", axis = "rl")
Figsjp

pdf("data/Exsjapan.pdf", width=20.55, height=14)
Figsjp
dev.off()

####Japan North ===================================================================
SPD2<-read.csv(file="data/JapanDates.csv", header=T)
box2<- subset(SPD2, Latitude>42 & Latitude< 46 & Longitude>138 & Longitude<146 )
#write.table(box2, file = "data/JapanNorthdates.csv", sep = ",", col.names=NA)
#box2<-read.csv(file="data/JapanNorthdates.csv", header=T)

japan_map <- map_data("world", region = "Japan")

# Create the map using ggplot2
ggplot(data = japan_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "lightblue", color = "black") +
  coord_fixed(1.3) +
  labs(title = "Map of Japan",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal()+
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Prefecture)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)

box2 <- box2[!is.na(box2$SiteNameEn), ]
box2 <- box2[!is.na(box2$CRA), ]

cptcal <- calibrate(x = box2$CRA,  errors = box2$CRAError, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$SiteNameEn, ages = box2$CRA, h = 100)

####Run analysis for component 3 logistic 3500 to 150
spd.CTx <- spd(cptcal, bins=boxbins, runm=200, timeRange=c(15000,200))
plot(spd.CTx, runm=200, xlim=c(15000,200), type="simple")

PrDens<-spd.CTx$grid$PrDens
calBP<-spd.CTx$grid$calBP

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(cptcal, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(15000,200),bw=50, normalised = FALSE)
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

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/JapanNorthKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/JapanNorthKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(14970,14940,14910,14880,14850,14820,14790,14760,14730,14700,14670,14640,14610,14580,14550,14520,14490,14460,14430,14400,14370,14340,14310,14280,14250,14220,14190,14160,14130,
         14100,14070,14040,14010,13980,13950,13920,13890,13860,13830,13800,13770,13740,13710,13680,13650,13620,13590,13560,13530,13500,13470,13440,13410,13380,13350,13320,13290,13260,
         13230,13200,13170,13140,13110,13080,13050,13020,12990,12960,12930,12900,12870,12840,12810,12780,12750,12720,12690,12660,12630,12600,12570,12540,12510,12480,12450,12420,12390,
         12360,12330,12300,12270,12240,12210,12180,12150,12120,12090,12060,12030,12000,11970,11940,11910,11880,11850,11820,11790,11760,11730,11700,11670,11640,11610,11580,11550,11520,
         11490,11460,11430,11400,11370,11340,11310,11280,11250,11220,11190,11160,11130,11100,11070,11040,11010,10980,10950,10920,10890,10860,10830,10800,10770,10740,10710,10680,10650,
         10620,10590,10560,10530,10500,10470,10440,10410,10380,10350,10320,10290,10260,10230,10200,10170,10140,10110,10080,10050,10020,9990,9960,9930,9900,9870,9840,9810,9780,9750,9720,
         9690,9660,9630,9600,9570,9540,9510,9480,9450,9420,9390,9360,9330,9300,9270,9240,9210,9180,9150,9120,9090,9060,9030,9000,8970,8940,8910,8880,8850,8820,8790,8760,8730,8700,8670,
         8640,8610,8580,8550,8520,8490,8460,8430,8400,8370,8340,8310,8280,8250,8220,8190,8160,8130,8100,8070,8040,8010,7980,7950,7920,7890,7860,7830,7800,7770,7740,7710,7680,7650,7620,
         7590,7560,7530,7500,7470,7440,7410,7380,7350,7320,7290,7260,7230,7200,7170,7140,7110,7080,7050,7020,6990,6960,6930,6900,6870,6840,6810,6780,6750,6720,6690,6660,6630,6600,6570,
         6540,6510,6480,6450,6420,6390,6360,6330,6300,6270,6240,6210,6180,6150,6120,6090,6060,6030,6000,5970,5940,5910,5880,5850,5820,5790,5760,5730,5700,5670,5640,5610,5580,5550,5520,
         5490,5460,5430,5400,5370,5340,5310,5280,5250,5220,5190,5160,5130,5100,5070,5040,5010,4980,4950,4920,4890,4860,4830,4800,4770,4740,4710,4680,4650,4620,4590,4560,4530,4500,4470,
         4440,4410,4380,4350,4320,4290,4260,4230,4200,4170,4140,4110,4080,4050,4020,3990,3960,3930,3900,3870,3840,3810,3780,3750,3720,3690,3660,3630,3600,3570,3540,3510,3480,3450,3420,
         3390,3360,3330,3300,3270,3240,3210,3180,3150,3120,3090,3060,3030,3000,2970,2940,2910,2880,2850,2820,2790,2760,2730,2700,2670,2640,2610,2580,2550,2520,2490,2460,2430,2400,2370,
         2340,2310,2280,2250,2220,2190,2160,2130,2100,2070,2040,2010,1980,1950,1920,1890,1860,1830,1800,1770,1740,1710,1680,1650,1620,1590,1560,1530,1500,1470,1440,1410,1380,1350,1320,
         1290,1260,1230,1200,1170,1140,1110,1080,1050,1020,990,960,930,900,870,840,810,780,750,720,690,660,630,600,570,540,510,480,450,420,390,360,330,300,270,240,210)
sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/JapanNorthSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("data/Sumbin/JapanNorthSumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 700, 1300, 2320, 3280,4000, 5000, 6000, 11500, 15000),
                         labels=c('Ainu','Satsumon','Epi Jamon','Final Jamon','Late Jamon','Middle Jamon', 'Early Jamon','Initial Jamon','Paleolithic'))
write.table(pcgrowth, file = "data/Percapita/JapanNorthPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
njp30pc<- read.csv("data/Percapita/JapanNorthPerCap.csv")

njp30pc2<-subset(njp30pc, calBP<6500 & calBP>3000)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(njp30pc2$MKDE-min(njp30pc2$MKDE))/(max((njp30pc2$MKDE)-min(njp30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
njp30pc3<-cbind(StKDE, njp30pc2)

pcnjp <- ggplot(njp30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),linewidth=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.1,0.2))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Mean KDE (density)", y="KDE per capita growth", title = "B. N. Japan KDE Per Capita Growth vs. Density")
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcnjp

njpcpt <- ggplot(njp30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Mean KDE", title = "B. N. Japan KDE vs. Time")
# geom_vline(xintercept =1020)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
njpcpt

###Paired Plots
Fignjp<-plot_grid(pcnjp, njpcpt, ncol=2, align="hv", axis = "rl")
Fignjp

pdf("data/Exnjapan.pdf", width=20.55, height=14)
Fignjp
dev.off()


##Northern Mendoza Eastern Andes=====================

###Load data
FullData<-read.csv(file="data/ArgentinaS1File.csv", header=T) ##(Supplemental Table 1)
SPD<-subset(FullData, Area=="North")

##calibrate dates from the North region of CW Argentina
CalMz <- calibrate(x = SPD$Age,  errors = SPD$Error, calCurves = "shcal20",  normalised = FALSE)
boxbins <- binPrep(sites = SPD$SiteName, ages = SPD$Age, h = 100)

####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.mz, runm=200, xlim=c(8200,200), type="simple")

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')
D.ckde$timeRange

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

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/NorthKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/NorthKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
##Add in the 30 year bin dates
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/North30Sumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("data/Sumbin/North30Sumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200,500, 1300, 2380, 3970, 8200),
                         labels=c('Phase 4','Phase 3', 'Phase 2','Phase 1','Archaic'))

write.table(pcgrowth, file = "data/Percapita/North30PerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
nmz30pc<- read.csv("data/Percapita/North30PerCap.csv")

nmz30pc2<-subset(nmz30pc, calBP<3000 & calBP>200)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(nmz30pc2$MKDE-min(nmz30pc2$MKDE))/(max((nmz30pc2$MKDE)-min(nmz30pc2$MKDE)))
##Add the standardized KDE to the Neolithic data frame
nmz30pc3<-cbind(StKDE, nmz30pc2)

pcnmdz <- ggplot(nmz30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),linewidth=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", 
       title = "A. N. Mendoza KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcnmdz

nmdzcpt <- ggplot(nmz30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. N. Mendoza Density vs. Time")+
  geom_hline(yintercept =.2)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
nmdzcpt

#paired Plots
Figmdz<-plot_grid(pcnmdz, nmdzcpt, ncol=2, align="hv", axis = "rl")
Figmdz

pdf("data/ExNmendoza.pdf", width=20.55, height=14)
Figmdz
dev.off()


##Northern Patagonia==========================================================================
FullData<-read.csv(file="data/ArgentinaS1File.csv", header=T) ##(Supplemental Table 1)
SPD<-subset(FullData, Area=="South")

##calibrate dates Southern Mendoza
CalMz <- calibrate(x = SPD$Age,  errors = SPD$Error, calCurves = "shcal20",  normalised = FALSE)
boxbins <- binPrep(sites = SPD$SiteName, ages = SPD$Age, h = 100)
####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.mz, runm=200, xlim=c(8200,200), type="simple")

##KDE
####Attempt KDE
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')
D.ckde$timeRange

##Write matrix of KDEs as a data frame
Checks<-as.data.frame(D.ckde$res.matrix)
#I then convert NAs at the end due to KDE bandwidth to zeros to easily cbind the data frames below
Checks2 <- replace(Checks, is.na(Checks), 0)
#Calculate the mean KDE of each time step of the 200 simulations
MKDE<-rowMeans(Checks2)
#calculate the 5th and 95th percentile of each time step of the 200 KDE simulations
lo<-apply( Checks2, #select columns
           1, # row-wise calcs
           quantile, probs=0.05) # give `quantile`
hi <-apply( Checks2, #select columns
            1, # row-wise calcs
            quantile, probs=0.95) # give `quantile`
calBPs<-spd.mz$grid$calBP
PrDenss<-spd.mz$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
ddc<-cbind(calBPs,PrDenss, Checks, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2c<-ddc %>%  filter(MKDE >0)
##Write the table
write.table(dd2c, file = "data/KDEs/SouthKDE50bin.csv", sep = ",", col.names=NA)

#load Center KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/SouthKDE50bin.csv") %>%
  dplyr::select(-X,-calBPs,-PrDenss, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps.
###Sum SPD by 30 year intervals
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
##Add in the 30 year bin dates
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)

write.table(sums, file = "Data/Sumbin/South30Sumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("data/Sumbin/South30Sumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200,500, 1300, 2380, 3970, 8200),
                         labels=c('Phase 4','Phase 3', 'Phase 2','Phase 1','Archaic'))
write.table(pcgrowth, file = "data/Percapita/South30PerCap.csv", sep = ",", col.names=NA)


###Plot mean KDE against the per capita growth rate in the North
npt30pc<- read.csv("data/Percapita/South30PerCap.csv")

npt30pc2<-subset(npt30pc, calBP<3000 & calBP>200)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(npt30pc2$MKDE-min(npt30pc2$MKDE))/(max((npt30pc2$MKDE)-min(npt30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
npt30pc3<-cbind(StKDE, npt30pc2)

pcnpat <- ggplot(npt30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=factor(PeriodID)), size=3.5) +
  geom_path(aes(),linewidth=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth",
       title = "A. N. Patagonia KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcnpat

npatcpt <- ggplot(npt30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=factor(PeriodID)), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. N. Patagonia Density vs. Time")+
  geom_hline(yintercept=.20)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
npatcpt

#Paired Plots

Fignpat<-plot_grid(pcnpat, npatcpt, ncol=2, align="hv", axis = "rl")
Fignpat

pdf("data/Figs/ExNpatagonia.pdf", width=20.55, height=14)
Fignpat
dev.off()


