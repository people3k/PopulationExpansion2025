
(1) Citation to the manuscript: This document describes the data files and R-code that accompany ``Innovations in food production and patterns of long-term population expansion" by Jacob Freeman and Erick Robinson. This README document describes the csv files and how to get started with the analysis. 


(2) Description of the files included in the repository and their relationship to the figures and tables in the manuscript
Files: 


(i) Radiocarbon files. These files are needed to run the code CaseStudyRadiocarbon.R For example: load data/FinalRCDTexas3.csv or data/nerd.csv.

These files contains all radiocarbon ages associated with archaeological remains collected from the archaeological regions. These files are essential to reproduce each file called RegionPerCap.csv and, ultimately, Figures 3 and S8-S38 in the main manuscript and the Supporting Material. These data are provided for researchers interested in reproducing the KDEs used in the study and engaging in their own analysis.

Typical Columns include:

Latitude-decimal degrees
Longitude-decimal degrees
SiteName--name of archaeological site
SiteID" Unique id number for an archaeological site.	
Trinomial--unique trinomial for each site	
Assay No.--unique id for lab and each radiocarbon sample	
Provenience--provenience of dated material if known
Feature Type--archaeological feature from which material was recovered
Material--material dated	
Taxa dated-species or genera of data material	
Human--yes or no. yes indicates human remains dated	
Age-- radiocarbon age	
Error--standard deviation of raw radiocarbon age	
Corrected/Raw	
Corrected/Normalized 14C age	 
Corrected/Normalized Age	
Delta 13C	
Delta 13C source	
Raw/Measured Age	 
Raw/Measured Age	
AMS or Radiometric	
Comments	
Reference
Region--Central Texas =CTx, Coastal Plain=TCP
Zone--Inland, riverine or coastal ecosystems


(ii) File path: data/RegionPerCap.csv 

These files contains the time-series of the per capita growth rates of all simulated KDEs, the mean KDE, and per capita growth rate of the mean KDE. The files are essential to reproduce population dynamics figures in the main manuscript and the Supporting Material. These files can be used by those people who do not want to reproduce the full analysis of archaeological radiocarbon based on the raw data.

Each regional file contains the same columns. For example, TCPPerCap.csv contains the following columns:

calBP--30 year time-steps in cal BP.
PeriodID--Cultural historical time period	
v1-v200: The per capita growth rate of each KDE in 30 year intervals
MKDE--The mean value of the 200 KDEs at each 30 year time-step
PerCap--The per capita growth rate of the mean KDE


(iii) File Path data/expansion.csv
This file contains all of the summary data necessary to replicate Figure 4 in the main text. This file contains the following columns

Case ID--unique ID number for each region.
Latitude-decimal degrees
Longitude-decimal degrees
AggID--States whether the sequence contains economies based on domesticated plants and or animals (Agg) or wild species (HG)
TimePeriod--the time period of study
Time-the length of the study sequence
cycles-number of counted short-term cycles
DomAnimal--0=no domesticated animals recorded. 1=Domesticated animals present.
Continent--Continent of each ase. Europe and Asia are treated as one continent (Eurasia).
Innovation--Innovation ranking from internal (1) to mixed (2), to external (3)
Innovation2--Same as Innovations just a string variable rather than numeric
CycleRate--Number of cycles divided by the Time variable
NumDates--Number of radiocarbon ages used to construct each sequence
ExpandK-The number of demographic transitions in a sequence
YearsperDT--The number of years per demographic transition
Kincrease--The increase in standardized population density per demographic transition in a sequence

(3)The code to analyze the data is contained in CaseStudyRadiocarbon.R, SumAnalysis.R, and PopExpandMap. The analysis was run in "R version 4.2.2 (2022-10-31 ucrt)" with the Rstudio integrated development environment version ``2023.06.0 Build 421." ChatGPT was used to help create the map.


(4) Getting started working with the project:

To begin working with the data, (i) open the R files. (ii) install and load all of the relevant packages at the beginning of the document, and (iii) follow the comments to begin reproducing the results reported in the manuscript. 
