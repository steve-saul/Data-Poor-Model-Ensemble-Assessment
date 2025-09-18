###################### Information Limited Fisheries Assessment Workflow ############################
######################################################################################
#Note: If you run the Catch MSY Portion of the code and are working out of a Dropbox or network syncing location....
#....you will need to "pause" syncing, or write to a non-network or cloud-based drive because the code is appending....
#....to the output file.

rm(list=ls(all=TRUE))
library(data.table)
library(bit64)
library(stats)
library(nloptr)
library(optimx)
library(fishmethods)
library(DBI)
library(RPostgreSQL)
library(RPostgres)
library(LBSPR)
library(zoo)
library(devtools)
#Note: need to first install Rtools from Rtools website. It is installed as an executable windows file, NOT installed as an R library.
#devtools::install_github("merrillrudd/LIME")
library(LIME)
#devtools::install_github("james-thorson/FishLife")   #This library is loaded later on in the code because it conflicts with rfishbase.
library(TropFishR)
library(ggplot2)
library(ggfan)
library(fanplot)
library(ggplot2)
library(ggspatial)
library(mapplots)
library(rgdal)
library(raster)
library(maps)
library(prettymapr)
library(rgeos)
library(MASS)

#################  Input Parameters (most will be dynamic from a parameter file once updated) ###############################################
#You need to manually create these two folders below on your computer first. The OtherInputs folder should be provided to you together with this script
PATH = "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/Assessment2024_2_23_2024"
#PATH = "C:/Users/sesaul/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/Assessment2023"
PATH_otherInput = paste(PATH,"OtherInputs",sep="/")
#Specify the location of your 7-zip software. You will need to download and install this Free software to extract the TNC data....
#...There is no work around for this because the file is password protected and R's build in upzip function doesn't work when the zip file is password protected.
PATH_7zip = paste(PATH,"7-Zip",sep="/")
#Here just specify the full paths and folder names for the remaining folders.
PATH_iFish = paste(PATH,"iFish",sep="/")
PATH_cpue = paste(PATH,"CPUE",sep="/")
PATH_output = paste(PATH,"Output",sep="/")
PATH_summaryTablesPaper = paste(PATH,"SummaryTablesPublication",sep="/")
PATH_histograms = paste(PATH,"Histograms",sep="/")
PATH_plots = paste(PATH,"Plots",sep="/")
PATH_historicalCatchPlots = paste(PATH,"HistoricalCatchPlots",sep="/")
PATH_catchMSY_diagnostics = paste(PATH,"CatchMSY_Diagnostics",sep="/")
PATH_catchMSY_plots = paste(PATH,"CatchMSY_Plots",sep="/")
PATH_saveProjectionLists = paste(PATH_output,"ProjectionResults",sep="/")
PATH_projectionPlots = paste(PATH,"ProjectionPlots",sep="/")
PATH_projectionPlots_byProjectType = paste(PATH,"ProjectionPlots_byProjectType",sep="/")
PATH_benchmarkPlots = paste(PATH,"BenchmarkPlots",sep="/")
PATH_benchmarkPlotsFamily = paste(PATH,"BenchmarkPlotsFamily",sep="/")
PATH_precentChangeF_histograms = paste(PATH,"PercentChangeF_histograms",sep="/")
PATH_summaryFanProjectionPlots = paste(PATH,"SummaryFanProjectionPlots",sep="/")
PATH_standardizedFanProjectionPlots = paste(PATH,"StandardizedFanProjectionPlots",sep="/")
PATH_highLevelSummaryPlots = paste(PATH,"HighLevelSummaryPlots",sep="/")
#Specify the species you want analyzed and the corresponding species codes.
GenusForAnalysis_names = c("Aphareus","Aprion","Argyrops","Atrobucca","Carangoides","Carangoides","Carangoides","Carangoides","Carangoides","Caranx","Caranx","Caranx","Caranx","Caranx","Cephalopholis","Cephalopholis","Cephalopholis","Cephalopholis","Cookeolus","Dentex","Diagramma","Diagramma","Elagatis","Epinephelus","Epinephelus","Epinephelus","Epinephelus","Epinephelus","Epinephelus","Epinephelus","Epinephelus",
  "Epinephelus","Epinephelus","Epinephelus","Epinephelus","Epinephelus","Epinephelus","Epinephelus","Epinephelus","Epinephelus","Epinephelus","Erythrocles","Etelis","Etelis","Etelis","Etelis","Glaucosoma","Gymnocranius","Gymnocranius","Hyporthodus","Lethrinus","Lethrinus","Lethrinus","Lethrinus","Lethrinus","Lethrinus","Lipocheilus","Lutjanus","Lutjanus","Lutjanus","Lutjanus","Lutjanus","Lutjanus","Lutjanus",
  "Lutjanus","Lutjanus","Lutjanus","Lutjanus","Lutjanus","Lutjanus","Lutjanus","Ostichthys","Paracaesio","Paracaesio","Paracaesio","Paracaesio","Parascolopsis","Pinjalo","Pinjalo","Plectropomus","Plectropomus","Pomadasys","Pristipomoides","Pristipomoides","Pristipomoides","Pristipomoides","Pristipomoides","Pristipomoides","Pristipomoides","Protonibea","Rachycentron","Saloptia","Seriola","Seriola","Sphyraena",
  "Sphyraena","Sphyraena","Symphorus","Variola","Wattsia","Other","Total")
#GenusForAnalysis_names = c("Lethrinus")
SpeciesForAnalysis_names = c("rutilans","virescens","spinifer","brevis","chrysophrys","coeruleopinnatus","fulvoguttatus","gymnostethus","malabaricus","bucculentus","ignobilis","lugubris","sexfasciatus",
  "tille","igarashiensis","miniata","sexmaculata","sonnerati","japonicus","carpenteri","labiosum","pictum","bipinnulata","amblycephalus","areolatus","bilobatus","bleekeri",
  "chlorostigma","coioides","epistictus","heniochus","latifasciatus","malabaricus","miliaris","morrhua","multinotatus","poecilonotus","radiatus","retouti","stictus",
  "undulosus","schlegelii","carbunculus","coruscans","radiosus","sp.","buergeri","grandoculis","griseus","octofasciatus","amboinensis","laticaudis","lentjan","nebulosus",
  "olivaceus","rubrioperculatus","carnolabrum","argentimaculatus","bitaeniatus","bohar","boutton","erythropterus","gibbus","johnii","lemniscatus","malabaricus","rivulatus",
  "russelli","sebae","timoriensis","vitta","japonicus","gonzalesi","kusakarii","stonei","xanthura","eriomma","lewisi","pinjalo","leopardus","maculatus","kaakan",
  "argyrogrammicus","filamentosus","flavipinnis","multidens","sieboldii","typus","zonatus","diacanthus","canadum","powelli","dumerili","rivoliana","barracuda",
  "forsteri","putnamae","nematophorus","albimarginata","mossambica","other","total")
#For Species Codes Below: use -999 as a placeholder for "Other" Species (i.e. all other biomass except the species we want separately).....
#...., and -9999 for "All Biomass", which provides an estimate of all reef fish together for biomass comparison later on.
SpeciesForAnalysis_codes = c("504404580790866828","504404580790866829","504404580790866917","1535446533334","1479214152054","1481081169755","1480666847385","1479214186452","1480666879099","1489993150718",
"504404580790866804","504404580790866805","1427691346501","1461212978732","504404580790866870","1483860901750","504404580790866871","504404580790866872","504404580790866869",
"504404580790866921","1457465274503","1457465391353","1461213079644","504404580790866873","504404580790866874","1427690470901","504404580790866875","504404580790866876",
"504404580790866877","504404580790866880","504404580790866883","504404580790866886","504404580790866889","504404580790866891","504404580790866892","504404580790866893",
"504404580790866894","504404580790866897","1427690629178","1427690796770","504404580790866903","504404580790866808","504404580790866832","504404580790866831","504404580790866833",
"504404580790866830","504404580790866810","504404580790866813","1427691198908","504404580790866905","504404580790866814","504404580790866817","504404580790866818","504404580790866821",
"504404580790866822","504404580790866823","504404580790866834","504404580790866835","1445981896175","504404580790866836","504404580790866837","504404580790866839",
"504404580790866840","504404580790866841","504404580790866842","504404580790866843","504404580790866845","504404580790866846","504404580790866847","504404580790866848",
"504404580790866849","504404580790866811","504404580790866851","504404580790866852","504404580790866853","504404580790866854","504404580790866868","504404580790866855",
"504404580790866856","504404580790866911","504404580790866912","1486438630050","1484033046946","504404580790866858","504404580790866859","504404580790866860","504404580790866861",
"504404580790866862","504404580790866863","1457465578947","1457465701638","1427690154905","504404580790866806","504404580790866807","1480666508070","1480666530132","1480666551916",
"504404580790866866","504404580790866914","504404580790866826","-999","-9999")
Species = data.frame(Species=paste(GenusForAnalysis_names,SpeciesForAnalysis_names))
Species$Codes = SpeciesForAnalysis_codes

#SpeciesForAnalysis_codes = c("504404580790866817")
FleetInventory_FileName = "FleetInventoryCapacity.csv"  #This file should be sent to you in the OtherInputs folder
#Specify the areas you want analyzed.
wholeCountry=TRUE       #If you also want the full EEZ analyzed with the pooled data, specify this as TRUE.
WPP = c(571,572,573,711,712,713,714,715,716,717,718)   #WPPs you want analyzed
#WPP = c(718,712)
#Flags to turn on and off different code sections
computeSPR_createInputFile=FALSE
setUpFilingSystem=FALSE   #WARNING: Setting this to "TRUE", will erase any output from a previous run of the code!!!!!
plotHistoricalCatches = FALSE
pullAndEstimateLifeHistoryParams = FALSE
plotFittedHistograms = FALSE
runCatchMSY = FALSE
plotCatchMSY_results = FALSE
plotCatchMSY_diagnostics = FALSE
ranCatchMSYonSeparateProcessors = FALSE
processCatchMSYoutput=FALSE
runLengthReconstruction = FALSE
ranLengthReconstructionSeparateProcessors=FALSE
devevlopLengthReconstructionSummmaryTable=FALSE
runProjections = FALSE
runRebuildingProjectionsRegardlessOfStatus = FALSE
createKobePlots = FALSE
plotProjections = FALSE
plotProjectionsByProjectionType = FALSE
createFanPlots = FALSE
createStandardizedFanPlots = FALSE
summarizeOutputForPublications = FALSE
fishingDays_fName = "fishing_days.csv"
CatchMSY_inputFileName = "StartingParametersCatchMSY.csv"         #Part of the OtherInputs folder. Note this file need to be filled out for each species and area you want analyzed.
histogramBinSizes_fName = "HistogramBinSizes.csv"
lifeHistoryValues_CODRS = "IFishEstimatedParameters_2023.csv"
#LandingsUncertaintyScenarios = c(0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8)  #These are used as multipliers for difference scenarios about what we might know about the landings
LandingsUncertaintyScenarios = c(1.0) #c(1.0,1.5,2.0)   #These are used as multipliers for difference scenarios about what we might know about the landings
CatchMSY_scenarioNames = c("Historical","Optimistic","Pessimistic","Constant")  # Don't use these other scenarios and deleted them from the code below "Historical_spr10","Optimistic_spr10","Pessimistic_spr10") #,"Historical_spr20","Optimistic_spr20","Pessimistic_spr20","Historical_spr30","Optimistic_spr30","Pessimistic_spr30")      #must match the length of the "CatchMSY_inputFileNames" vector above
CatchMSY_outfileName = "CatchMSY_Output.csv"    #Name of output file from CatchMSY
cMSY_optimisticScenario_yrsToCurrentDepletion = 75      #CatchMSY optimistic scenario timeframe in years 1950 start year
cMSY_pessimisticScenario_yrsToCurrentDepletion = 55     #CatchMSY pessimistic scenario timeframe in years 1970 start year
cMSY_constantScenario_yrsToCurrentDepletion = 15
populationReconstructionMethods = c("SolveBaranov","LBSPR","BevertonHoltInstantaneousMortality","LIME","CatchCurve")    #Do Not use TropFishR anymore - it is a VPA, so only useful if you have multiple years of length data over time
yearsToProject = 76
projectionYearsToPlot = 20
controlRuleSPR30_percent = 0.3
controlRuleSPR20_percent = 0.2
controlRuleSPR40_percent = 0.4
terminalYear = 2023
yearToStartProjection = 2024
WaterTempAssume_base = 20
Steepness_Scenarios = c(0.7,0.8,0.9)
ishingMortaltiyProjectionScenarios=c("CurrentF","F_equal_0","F_at_SPR20","F_at_SPR30","F_at_SPR40","MSY","OpenAccess","MEY") #,"F_half","F_double","ConstantCatch")
ColorsFor_ProjectionF_ScenarioPlots = c("red","royalblue","green","grey","yellow","purple","orange","brown","turquoise","pink","wheat2","gold3")   #must be same length as vector above plus one additonal color.



##################################### Read in CPUE Estimates and Read Into R ######################################
RawDataAverageCPUE_SpeciesComposition_kgPerLanding = read.csv(paste(PATH_cpue,"CPUE_Raw_kg_per_landing.csv",sep="/"),header=TRUE,sep=",")
AverageCPUE_bySpecies_allAreas = read.csv(paste(PATH_cpue,"AverageCPUE_bySpecies_allAreas.csv",sep="/"),header=TRUE,sep=",")
AverageCPUE_bySpeciesWPP = read.csv(paste(PATH_cpue,"AverageCPUE_bySpeciesWPP.csv",sep="/"),header=TRUE,sep=",")
AverageCPUE_WPP_AllSpeciesCombined = read.table(paste(PATH_cpue,"AverageCPUE_WPP_AllSpeciesCombined.csv",sep="/"),header=TRUE,sep=",")

####################################### Process IFish Data and Format for Analysis ###################

Lengths = read.table(paste(PATH_iFish,"ifish_sizing.txt",sep="/"),sep="$",header=TRUE,comment.char="",colClasses=c("character"),stringsAsFactors=FALSE,quote="",fill=TRUE)
Species = read.table(paste(PATH_iFish,"ifish_fish.txt",sep="/"),sep="$",header=TRUE,comment.char="",colClasses=c("character"),stringsAsFactors=FALSE,quote="",fill=TRUE)
Landings = read.table(paste(PATH_iFish,"ifish_deepslope.txt",sep="/"),sep="$",header=TRUE,comment.char="",colClasses=c("character"),stringsAsFactors=FALSE,quote="",fill=TRUE)
Boats = read.table(paste(PATH_iFish,"ifish_boat_pub.txt",sep="/"),sep="$",header=TRUE,comment.char="",colClasses=c("character"),stringsAsFactors=FALSE,quote="",fill=TRUE)
BoatsRAW = read.table(paste(PATH_iFish,"ifish_boat_pub.txt",sep="/"),sep="$",header=TRUE,comment.char="",colClasses=c("character"),stringsAsFactors=FALSE,quote="",fill=TRUE)
LocationPings = read.table(paste(PATH_iFish,"ifish_findmespot.txt",sep="/"),sep="$",header=TRUE,comment.char="",colClasses=c("character"),stringsAsFactors=FALSE,quote="",fill=TRUE)
BoatTracker = read.table(paste(PATH_iFish,"ifish_boat_tracker.txt",sep="/"),sep="$",header=TRUE,comment.char="",colClasses=c("character"),stringsAsFactors=FALSE,quote="",fill=TRUE)
TrackerList = read.table(paste(PATH_iFish,"ifish_tracker.txt",sep="/"),sep="$",header=TRUE,comment.char="",colClasses=c("character"),stringsAsFactors=FALSE,quote="",fill=TRUE)
iSnapperBoatSurvey2020 = read.csv(paste(PATH_otherInput,"IFishSnapperFleet_2023.csv",sep="/"),header=TRUE,sep=",",colClasses=c("character"),stringsAsFactors=FALSE)

#Remove eroneous duplicate record from species
Species$flag=0
Species$flag[Species$fish_species=="timorensis" & Species$common_name=="Yellowspotted grouper"]=1
Species = subset(Species,Species$flag!=1)
Species$flag=NULL

#Fix name of a species that is misspelled
Species$fish_species[Species$fish_species=="timorensis"]="timoriensis"

#Convert relevant variables to numeric
Lengths$cm = as.numeric(Lengths$cm)
Species$counter = as.numeric(Species$counter)
Species$lmat = as.numeric(Species$lmat)
Species$lopt = as.numeric(Species$lopt)
Species$linf = as.numeric(Species$linf)
Species$lmax = as.numeric(Species$lmax)
Species$lmatm = as.numeric(Species$lmatm)
Species$is_family_id = as.numeric(Species$is_family_id)
Species$species_id_number = as.numeric(Species$species_id_number)
Species$reported_trade_limit_weight = as.numeric(Species$reported_trade_limit_weight)
Species$var_a = as.numeric(Species$var_a)
Species$var_b = as.numeric(Species$var_b)
Species$converted_trade_limit_l = as.numeric(Species$converted_trade_limit_l)
Species$plotted_trade_limit_tl = as.numeric(Species$plotted_trade_limit_tl)
Species$conversion_factor_tl2fl = as.numeric(Species$conversion_factor_tl2fl)
Species$largest_specimen_cm = as.numeric(Species$largest_specimen_cm)
Landings$wpp1 = as.numeric(Landings$wpp1)
Landings$wpp2 = as.numeric(Landings$wpp2)
Landings$wpp3 = as.numeric(Landings$wpp3)
Landings$fishing_ledger = as.numeric(Landings$fishing_ledger)
Landings$total_catch = as.numeric(Landings$total_catch)
Boats$year_built = as.numeric(Boats$year_built)
Boats$length_of_boat = as.numeric(Boats$length_of_boat)
Boats$width_of_boat = as.numeric(Boats$width_of_boat)
Boats$height_of_boat = as.numeric(Boats$height_of_boat)
Boats$capacity_palka_kg = as.numeric(Boats$capacity_palka_kg)
Boats$gt_estimate = as.numeric(Boats$gt_estimate)
Boats$gt_declared = as.numeric(Boats$gt_declared)
Boats$engine_hp1 = as.numeric(Boats$engine_hp1)
Boats$engine_hp2 = as.numeric(Boats$engine_hp2)
Boats$engine_hp3 = as.numeric(Boats$engine_hp3)
Boats$number_of_crew = as.numeric(Boats$number_of_crew)
Boats$fishing_wpp1 = as.numeric(Boats$fishing_wpp1)
Boats$fishing_wpp2 = as.numeric(Boats$fishing_wpp2)
Boats$fishing_wpp3 = as.numeric(Boats$fishing_wpp3)
LocationPings$latitude = as.numeric(LocationPings$latitude)
LocationPings$longitude = as.numeric(LocationPings$longitude)
LocationPings$daily_avg_latitude = as.numeric(LocationPings$daily_avg_latitude)
LocationPings$daily_avg_longitude = as.numeric(LocationPings$daily_avg_longitude)

#################### Clean Data, Combine Components, and Check Sample Sizes #####################################
#Length data
Species = subset(Species,select=c(oid,fish_phylum,fish_class,fish_order,fish_family,fish_genus,
	fish_species,var_a,var_b,length_basis,conversion_factor_tl2fl))
Species = subset(Species,Species$fish_species!="")
names(Species)[names(Species)=="oid"]="fish_id"
Species$length_basis[Species$length_basis=="SL"]="FL"   #Make assumption that SL is FL because we don't have the conversion factors for SL to TL (SL is close to FL for many species) - best we can do
Species$var_a[Species$length_basis=="FL" & Species$conversion_factor_tl2fl > 0] = Species$var_a[Species$length_basis=="FL" & Species$conversion_factor_tl2fl > 0] * Species$conversion_factor_tl2fl[Species$length_basis=="FL" & Species$conversion_factor_tl2fl > 0]   #Convert fork length a and b length-weight parameters to TL
Species$length_basis[Species$length_basis=="FL" & Species$conversion_factor_tl2fl > 0]="TL"  #re-label those converted to TL as total length. NOTE: some species are in standard length 
Species$length_basis=NULL
Lengths = subset(Lengths,select=c(fish_id,cm,landing_id,length_type,codrs_picture_date))
Lengths = merge(Lengths,Species,all.x=TRUE,by=c("fish_id"))
Lengths = subset(Lengths,Lengths$fish_id!="0")   #Remove any records that don't have a species id associated
Lengths = subset(Lengths,!is.na(Lengths$fish_species))    #Remove records with no species name associated with them
Lengths$flag=0
Lengths$flag[Lengths$length_type=="FL" & Lengths$conversion_factor_tl2fl==0]=1
Lengths = subset(Lengths,Lengths$flag==0)
Lengths$flag=NULL
Lengths$length_type[Lengths$length_type==""]="TL"  #If length type is blank, assume it is in total length
Lengths$length_type[Lengths$length_type=="SL"]="FL"  #If length type is standard length, assume it is in fork length, since we don't have a conversion for SL to TL
Lengths$cm[Lengths$length_type=="FL" & Lengths$conversion_factor_tl2fl > 0] = Lengths$cm[Lengths$length_type=="FL" & Lengths$conversion_factor_tl2fl > 0] / Lengths$conversion_factor_tl2fl[Lengths$length_type=="FL" & Lengths$conversion_factor_tl2fl > 0]   #Convert fork length a and b length-weight parameters to TL
Lengths$length_type[Lengths$length_type=="FL" & Lengths$conversion_factor_tl2fl > 0]="TL"  #re-label those converted to TL as total length. NOTE: some species are in standard length 
#Boat Data
Boats = subset(Boats,select=c(oid,fishing_gear1,landing_port1,gt_estimate,category,fishing_wpp1,fishing_wpp2,fishing_wpp3))
names(Boats)[names(Boats)=="oid"]="boat_id"
names(Boats)[names(Boats)=="fishing_gear1"]="gear"
names(Boats)[names(Boats)=="landing_port1"]="landing_port"
Boats$WPP = Boats$fishing_wpp1
Boats$WPP[is.na(Boats$fishing_wpp1)]=Boats$fishing_wpp2[is.na(Boats$fishing_wpp1)]
Boats$WPP[is.na(Boats$fishing_wpp1) & is.na(Boats$fishing_wpp2)]=Boats$fishing_wpp3[is.na(Boats$fishing_wpp1) & is.na(Boats$fishing_wpp2)]
Boats$fishing_wpp1=NULL
Boats$fishing_wpp2=NULL
Boats$fishing_wpp3=NULL
Boats = unique(Boats)
#Vessel Survey Data - reconsile it with "Boats" data.frame and combine both into "Boats"
iSnapperBoatSurvey2020 = subset(iSnapperBoatSurvey2020,select=c(oid,fishing_gear1,landing_port1,gt_estimate,category,fishing_wpp1))
iSnapperBoatSurvey2020$gt_estimate = as.numeric(iSnapperBoatSurvey2020$gt_estimate)
iSnapperBoatSurvey2020$fishing_wpp1 = as.numeric(iSnapperBoatSurvey2020$fishing_wpp1)
names(iSnapperBoatSurvey2020)[names(iSnapperBoatSurvey2020)=="oid"]="boat_id"
names(iSnapperBoatSurvey2020)[names(iSnapperBoatSurvey2020)=="fishing_gear1"]="gear_boatSurvey"
names(iSnapperBoatSurvey2020)[names(iSnapperBoatSurvey2020)=="landing_port1"]="landing_port_boatSurvey"
names(iSnapperBoatSurvey2020)[names(iSnapperBoatSurvey2020)=="gt_estimate"]="gt_estimate_boatSurvey"
names(iSnapperBoatSurvey2020)[names(iSnapperBoatSurvey2020)=="category"]="category_boatSurvey"
names(iSnapperBoatSurvey2020)[names(iSnapperBoatSurvey2020)=="fishing_wpp1"]="WPP_boatSurvey"
Boats = merge(Boats,iSnapperBoatSurvey2020,all=TRUE,by=c("boat_id"))
Boats$gear[Boats$gear=="DriftOneLargeNatural"]="Other"
Boats$gear[Boats$gear=="DropMultipleSmallArtificial"]="Dropline"
Boats$gear[Boats$gear=="DropOneLargeArtificial"]="Dropline"
Boats$gear[Boats$gear=="DropOneLargeNatural"]="Dropline"
Boats$gear[Boats$gear=="DropOneSmallArtificial"]="Dropline"
Boats$gear[Boats$gear=="DropOneSmallNatural"]="Dropline"
Boats$gear[Boats$gear=="LongMultipleLargeNatural"]="Longline"
Boats$gear[Boats$gear=="PoleAndLine"]="Dropline"
Boats$gear[Boats$gear=="SurfaceOneLargeLive"]="Other"
Boats$gear[Boats$gear=="TrollMultipleLargeArtificial"]="Other"
Boats$gear[Boats$gear=="TrollMultipleSmallArtificial"]="Other"
Boats$gear[Boats$gear=="TrollOneLargeArtificial"]="Other"
Boats$gear[is.na(Boats$gear)]=Boats$gear_boatSurvey[is.na(Boats$gear)]
Boats$gear[Boats$gear==""]=Boats$gear_boatSurvey[Boats$gear==""]
Boats$gear_boatSurvey=NULL
Boats$landing_port[is.na(Boats$landing_port)]=Boats$landing_port_boatSurvey[is.na(Boats$landing_port)]
Boats$landing_port[Boats$landing_port==""]=Boats$landing_port_boatSurvey[Boats$landing_port==""]
Boats$landing_port_boatSurvey=NULL
Boats$gt_estimate[is.na(Boats$gt_estimate)]=Boats$gt_estimate_boatSurvey[is.na(Boats$gt_estimate)]
Boats$gt_estimate_boatSurvey=NULL
Boats$category[Boats$category==""]=NA
Boats$category[is.na(Boats$category)]=Boats$category_boatSurvey[is.na(Boats$category)]
Boats$category_boatSurvey=NULL
Boats$WPP[is.na(Boats$WPP)]=Boats$WPP_boatSurvey[is.na(Boats$WPP)]
Boats$WPP_boatSurvey=NULL
Boats$VesselCategory=""
Boats$VesselCategory[Boats$gt_estimate<5]="Nano"
Boats$VesselCategory[Boats$gt_estimate>=5 & Boats$gt_estimate<15]="Small"
Boats$VesselCategory[Boats$gt_estimate>=15 & Boats$gt_estimate<30]="Medium"
Boats$VesselCategory[Boats$gt_estimate>=30]="Large"
#Landings data
names(Landings)[names(Landings)=="oid"]="landing_id"
Landings$wpp2[is.na(Landings$wpp2)]=Landings$wpp3[is.na(Landings$wpp2)]
Landings$wpp1[is.na(Landings$wpp1)]=Landings$wpp2[is.na(Landings$wpp1)]
Landings$WPP_landings = Landings$wpp1
Landings = subset(Landings,select=c(landing_id,boat_id,WPP_landings,fishery_type))
names(Boats)[names(Boats)=="WPP"]="WPP_boats"   #reconsile WPP between Boats dataset and Landings
Landings = merge(Landings,Boats,all.x=TRUE,by=c("boat_id"))  #reconsile WPP between Boats dataset and Landings
Landings$WPP_landings[is.na(Landings$WPP_landings)] = Landings$WPP_boats[is.na(Landings$WPP_landings)]  #reconsile WPP between Boats dataset and Landings
Landings$WPP_boats=NULL
names(Landings)[names(Landings)=="WPP_landings"]="WPP"
Lengths = merge(Lengths,Landings,all.x=TRUE,by=c("landing_id"))
#separate date
Lengths = subset(Lengths,Lengths$codrs_picture_date!="\\N")
dateTimeVec_vec = lapply(X=Lengths$codrs_picture_date,FUN=function(x) {unlist(strsplit(as.character(x)," "))})
dateTimeVec_vec = do.call(rbind.data.frame,dateTimeVec_vec)
rownames(dateTimeVec_vec)=NULL
names(dateTimeVec_vec) = c("Date","Time")
dateOnlyVec_vec = lapply(X=dateTimeVec_vec$Date,FUN=function(x) {unlist(strsplit(as.character(x),"-"))})
dateOnlyVec_vec = do.call(rbind.data.frame,dateOnlyVec_vec)
rownames(dateOnlyVec_vec)=NULL
names(dateOnlyVec_vec) = c("Year_OutFishing","Month_OutFishing","Day_OutFishing")
dateOnlyVec_vec$Year_OutFishing = as.numeric(as.character(dateOnlyVec_vec$Year_OutFishing))
dateOnlyVec_vec$Month_OutFishing = as.numeric(as.character(dateOnlyVec_vec$Month_OutFishing))
dateOnlyVec_vec$Day_OutFishing = as.numeric(as.character(dateOnlyVec_vec$Day_OutFishing))
Lengths = cbind(Lengths,dateOnlyVec_vec)
Lengths$codrs_picture_date=NULL
Lengths = subset(Lengths,Lengths$Year_OutFishing<=terminalYear)  #Remove records from years after the terminal year as the year label is wrong! 
Lengths$species = paste(Lengths$fish_genus,Lengths$fish_species)
#Not necessary I think- keep here as placeholder: use only samples from the "snapper fishery"
#Lengths = subset(Lengths,Lengths$fishery_type=="BandaSnapper" |  Lengths$fishery_type=="Snapper" | Lengths$fishery_type=="SnapperFIP2B")

#Check sample sizes
SampleSizes = as.data.frame(as.data.frame.matrix(table(Lengths$species,Lengths$WPP)))
SampleSizes$EEZ = rowSums(SampleSizes)
CanAssess = as.data.frame(SampleSizes>150)
SampleSizes = subset(SampleSizes,CanAssess$EEZ==TRUE)
CanAssess = subset(CanAssess,CanAssess$EEZ==TRUE)
names(CanAssess) = paste(names(CanAssess),"_toAssess",sep="")
SampleSizesCanAssess = cbind(SampleSizes,CanAssess)
SampleSizesCanAssess$species = rownames(SampleSizesCanAssess)
rownames(SampleSizesCanAssess)=NULL
Spp_fishId = subset(Lengths,select=c(species,fish_id))
Spp_fishId = unique(Spp_fishId)
SampleSizesCanAssess = merge(SampleSizesCanAssess,Spp_fishId,all.x=TRUE,by=c("species"))

################################################################################################################
################# Set Up Filing System [Now that we know what species will be modeled] #########################
################################################################################################################
OptimalHistogramBreaks = read.table(paste(PATH_otherInput,histogramBinSizes_fName,sep="/"),header=TRUE,sep=",")
SpeciesFolderNames = gsub(" ","_",OptimalHistogramBreaks$species)
if(setUpFilingSystem==TRUE)
{
#  #Histogram plots of length data
#  if(dir.exists(PATH_histograms)==TRUE)
#  {
#    unlink(PATH_histograms,recursive=TRUE)
#  }
#  dir.create(PATH_histograms)
#  for(i in 1:length(SpeciesFolderNames))
#  {
#    if(dir.exists(paste(PATH_histograms,SpeciesFolderNames[i],sep="/"))==TRUE)
#    {
#      unlink(paste(PATH_histograms,SpeciesFolderNames[i],sep="/"),recursive=TRUE)
#    }
#    dir.create(paste(PATH_histograms,SpeciesFolderNames[i],sep="/"))
#  }
  #plots
  if(dir.exists(PATH_plots)==TRUE)
  {
    unlink(PATH_plots,recursive=TRUE)
  }
  dir.create(PATH_plots)
  for(i in 1:length(SpeciesFolderNames))
  {
    if(dir.exists(paste(PATH_plots,SpeciesFolderNames[i],sep="/"))==TRUE)
    {
      unlink(paste(PATH_plots,SpeciesFolderNames[i],sep="/"),recursive=TRUE)
    }
    dir.create(paste(PATH_plots,SpeciesFolderNames[i],sep="/"))
  }
  #Historical catch time series plots
  if(dir.exists(PATH_historicalCatchPlots)==TRUE)
  {
    unlink(PATH_historicalCatchPlots,recursive=TRUE)
  }
  dir.create(PATH_historicalCatchPlots)
  for(i in 1:length(SpeciesFolderNames))
  {
    if(dir.exists(paste(PATH_historicalCatchPlots,SpeciesFolderNames[i],sep="/"))==TRUE)
    {
      unlink(paste(PATH_historicalCatchPlots,SpeciesFolderNames[i],sep="/"),recursive=TRUE)
    }
    dir.create(paste(PATH_historicalCatchPlots,SpeciesFolderNames[i],sep="/"))
  }
  #CatchMSY Diagnostic Plots
  if(dir.exists(PATH_catchMSY_diagnostics)==TRUE)
  {
    unlink(PATH_catchMSY_diagnostics,recursive=TRUE)
  }
  dir.create(PATH_catchMSY_diagnostics)
 
  #CatchMSY Results Plots
  if(dir.exists(PATH_catchMSY_plots)==TRUE)
  {
    unlink(PATH_catchMSY_plots,recursive=TRUE)
  }
  dir.create(PATH_catchMSY_plots)
  for(i in 1:length(SpeciesFolderNames))
  {
    if(dir.exists(paste(PATH_catchMSY_plots,SpeciesFolderNames[i],sep="/"))==TRUE)
    {
      unlink(paste(PATH_catchMSY_plots,SpeciesFolderNames[i],sep="/"),recursive=TRUE)
    }
    dir.create(paste(PATH_catchMSY_plots,SpeciesFolderNames[i],sep="/"))
  }
#  #Projection spaghetti plots - traditional, by species and area with all projection scenarios on the same figure
#  if(dir.exists(PATH_projectionPlots)==TRUE)
#  {
#    unlink(PATH_projectionPlots,recursive=TRUE)
#  }
#  dir.create(PATH_projectionPlots)
#  for(i in 1:length(SpeciesFolderNames))
#  {
#    if(dir.exists(paste(PATH_projectionPlots,SpeciesFolderNames[i],sep="/"))==TRUE)
#    {
#      unlink(paste(PATH_projectionPlots,SpeciesFolderNames[i],sep="/"),recursive=TRUE)
#    }
#    dir.create(paste(PATH_projectionPlots,SpeciesFolderNames[i],sep="/"))
#  }
#  #Projection spaghetti plots - with one plot for each projection type, and the different scenarios on that figure
#  if(dir.exists(PATH_projectionPlots_byProjectType)==TRUE)
#  {
#    unlink(PATH_projectionPlots_byProjectType,recursive=TRUE)
#  }
#  dir.create(PATH_projectionPlots_byProjectType)
#  for(i in 1:length(SpeciesFolderNames))
#  {
#    if(dir.exists(paste(PATH_projectionPlots_byProjectType,SpeciesFolderNames[i],sep="/"))==TRUE)
#    {
#      unlink(paste(PATH_projectionPlots_byProjectType,SpeciesFolderNames[i],sep="/"),recursive=TRUE)
#    }
#    dir.create(paste(PATH_projectionPlots_byProjectType,SpeciesFolderNames[i],sep="/"))
#  }
  #Kobe Benchmark Plots
  if(dir.exists(PATH_benchmarkPlots)==TRUE)
  {
    unlink(PATH_benchmarkPlots,recursive=TRUE)
  }
  dir.create(PATH_benchmarkPlots)
  for(i in 1:length(SpeciesFolderNames))
  {
    if(dir.exists(paste(PATH_benchmarkPlots,SpeciesFolderNames[i],sep="/"))==TRUE)
    {
      unlink(paste(PATH_benchmarkPlots,SpeciesFolderNames[i],sep="/"),recursive=TRUE)
    }
    dir.create(paste(PATH_benchmarkPlots,SpeciesFolderNames[i],sep="/"))
  }
  #Percent Change in F Histograms
  if(dir.exists(PATH_precentChangeF_histograms)==TRUE)
  {
    unlink(PATH_precentChangeF_histograms,recursive=TRUE)
  }
  dir.create(PATH_precentChangeF_histograms)
  for(i in 1:length(SpeciesFolderNames))
  {
    if(dir.exists(paste(PATH_precentChangeF_histograms,SpeciesFolderNames[i],sep="/"))==TRUE)
    {
      unlink(paste(PATH_precentChangeF_histograms,SpeciesFolderNames[i],sep="/"),recursive=TRUE)
    }
    dir.create(paste(PATH_precentChangeF_histograms,SpeciesFolderNames[i],sep="/"))
  }
  #Summary Fan Plots
  if(dir.exists(PATH_summaryFanProjectionPlots)==TRUE)
  {
    unlink(PATH_summaryFanProjectionPlots,recursive=TRUE)
  }
  dir.create(PATH_summaryFanProjectionPlots)
  for(i in 1:length(SpeciesFolderNames))
  {
    if(dir.exists(paste(PATH_summaryFanProjectionPlots,SpeciesFolderNames[i],sep="/"))==TRUE)
    {
      unlink(paste(PATH_summaryFanProjectionPlots,SpeciesFolderNames[i],sep="/"),recursive=TRUE)
    }
    dir.create(paste(PATH_summaryFanProjectionPlots,SpeciesFolderNames[i],sep="/"))
  }
  #Standardized Fan Plots
  if(dir.exists(PATH_standardizedFanProjectionPlots)==TRUE)
  {
    unlink(PATH_standardizedFanProjectionPlots,recursive=TRUE)
  }
  dir.create(PATH_standardizedFanProjectionPlots)
  for(i in 1:length(SpeciesFolderNames))
  {
    if(dir.exists(paste(PATH_standardizedFanProjectionPlots,SpeciesFolderNames[i],sep="/"))==TRUE)
    {
      unlink(paste(PATH_standardizedFanProjectionPlots,SpeciesFolderNames[i],sep="/"),recursive=TRUE)
    }
    dir.create(paste(PATH_standardizedFanProjectionPlots,SpeciesFolderNames[i],sep="/"))
  }
  #High Level Summary Plots
  if(dir.exists(PATH_highLevelSummaryPlots)==TRUE)
  {
    unlink(PATH_highLevelSummaryPlots,recursive=TRUE)
  }
  dir.create(PATH_highLevelSummaryPlots)
}
#################################################################################################################
############################### Create Bins in Length Data ###################################################
#################################################################################################################
#NOTE: I Created histograms first for each species and WPP to see if there are any outliers....
#....and built a table that contains the histogram bin size for each species and WPP. 
#Spp="Wattsia mossambica"
#Wpp = "EEZ"
#Brk = 1
#if(Wpp!="EEZ")
#{
#	Sub = subset(Lengths,Lengths$species==Spp & Lengths$wpp1==Wpp)
#}
#if(Wpp=="EEZ")
#{
#	Sub = subset(Lengths,Lengths$species==Spp)
#}
#Max = (ceiling(max(Sub$cm)/Brk))*Brk
#windows(6,6)
#Hist = hist(Sub$cm,breaks=seq(0,Max,by=Brk),xlab="Length(cm)",main=Spp)
#mtext(paste(Wpp," | Breaks=",Brk,sep=""))

#Read in Pre-Determined Length Bins as Determine From Above Code
OptimalHistogramBreaks = read.table(paste(PATH_otherInput,histogramBinSizes_fName,sep="/"),header=TRUE,sep=",")

#Set up species names for life history parameter download and read in the CODRS life history parameters
GenusSpeciesForAnalysis_names = OptimalHistogramBreaks$species
GenusSpeciesSplit = transpose(as.data.frame(strsplit(OptimalHistogramBreaks$species," ")))
names(GenusSpeciesSplit) = c("Genus","Species")
GenusForAnalysis_names = GenusSpeciesSplit$Genus
SpeciesForAnalysis_names = GenusSpeciesSplit$Species
LifeHistoryParams_CODRS = read.table(paste(PATH_otherInput,lifeHistoryValues_CODRS,sep="/"),header=TRUE,sep=",")
LifeHistoryParams_CODRS$ID=NULL
LifeHistoryParams_CODRS$Family=NULL
LifeHistoryParams_CODRS$N.Last365Days=NULL
LifeHistoryParams_CODRS$Fmax=NULL   #NOTE: Don't use the Fmax value in the file because it is oddly too low!!!
LifeHistoryParams_CODRS$Mnat=NULL   #NOTE: Don't use the Mnat value in the file because it is oddly too low!!!
LifeHistoryParams_CODRS$VBG_K=NULL  #NOTE: Don't use the VBG_K value because it is the same value repeated over and over again
LifeHistoryParams_CODRS$SPR..=NULL  #NOTE: Don't use the SPR value because i will calculate it on my own in the code
LifeHistoryParams_CODRS$SPR.Risk.Level=NULL  #Drop the SPR risk level column
names(LifeHistoryParams_CODRS)[names(LifeHistoryParams_CODRS)=="Var_a"]="var_a"
names(LifeHistoryParams_CODRS)[names(LifeHistoryParams_CODRS)=="Var_b"]="var_b"
names(LifeHistoryParams_CODRS)[names(LifeHistoryParams_CODRS)=="Species.Name"]="species"
Species$species = paste(Species$fish_genus,Species$fish_species,sep=" ")
Species$var_a=NULL   #Use the var_a from the other life history FILE
Species$var_b=NULL   #Use the var_b from the other life history FILE
Species$fish_id=NULL  #no need for fish_id field if we use the genus and species
Species = merge(Species,LifeHistoryParams_CODRS,all=TRUE,by=c("species"))
Species = merge(data.frame(species=OptimalHistogramBreaks$species),Species,all.x=TRUE,by=c("species"))
Species$Lc[Species$Lc==0]=NA

#Create and Assign Group Common Names for Each Fish Family
Lengths$FamilyCommonName=""
Lengths$FamilyCommonName[Lengths$fish_family=="Carangidae"]="Jacks_Mackerels"
Lengths$FamilyCommonName[Lengths$fish_family=="Scombridae"]="Jacks_Mackerels"
Lengths$FamilyCommonName[Lengths$fish_family=="Emmelichthydae"]="Other"
Lengths$FamilyCommonName[Lengths$fish_family=="Emmelichthyidae"]="Other"
Lengths$FamilyCommonName[Lengths$fish_family=="Epinephelidae"]="Grouper"
Lengths$FamilyCommonName[Lengths$fish_family=="Glaucosomatidae"]="Other"
Lengths$FamilyCommonName[Lengths$fish_family=="Haemulidae"]="Grunts"
Lengths$FamilyCommonName[Lengths$fish_family=="Holocentridae"]="Other"
Lengths$FamilyCommonName[Lengths$fish_family=="Lethrinidae"]="Emperor"
Lengths$FamilyCommonName[Lengths$fish_family=="Lutjanidae"]="Snapper"
Lengths$FamilyCommonName[Lengths$fish_family=="Nemipteridae"]="SmallPelagic"
Lengths$FamilyCommonName[Lengths$fish_family=="Priacanthidae"]="Bigeyes"
Lengths$FamilyCommonName[Lengths$fish_family=="Rachycentridae"]="Cobia"
Lengths$FamilyCommonName[Lengths$fish_family=="Sciaenidae"]="Drums_Croaker"
Lengths$FamilyCommonName[Lengths$fish_family=="Sparidae"]="Porgies"
Lengths$FamilyCommonName[Lengths$fish_family=="Sphyraenidae"]="Barracuda"
Lengths$FamilyCommonName[Lengths$fish_family=="Acanthuridae"]="Surgeonfish"
Lengths$FamilyCommonName[Lengths$fish_family=="Albulidae"]="Bonefish"
Lengths$FamilyCommonName[Lengths$fish_family=="Malacanthidae"]="Other"
Lengths$FamilyCommonName[Lengths$fish_family=="Balistidae"]="Triggerfish"
Lengths$FamilyCommonName[Lengths$fish_family=="Priacanthidae"]="Bigeyes"
Lengths$FamilyCommonName[Lengths$fish_family=="Coryphaenidae"]="Dolphinfish"
Lengths$FamilyCommonName[Lengths$fish_family=="Serranidae"]="Grouper"
Lengths$FamilyCommonName[Lengths$fish_family=="Glaucosomatidae"]="Other"
Lengths$FamilyCommonName[Lengths$fish_family=="Lampridae"]="Other"
Lengths$FamilyCommonName[Lengths$fish_family=="Istiophoridae"]="LargePelagic"
Lengths$FamilyCommonName[Lengths$fish_family=="miscoded"]="Other"
Lengths$FamilyCommonName[Lengths$fish_family=="Ariidae"]="Other"
Lengths$FamilyCommonName[Lengths$fish_family=="Xiphiidae"]="LargePelagic"
SpeciesFamilyCommonName = subset(Lengths,select=c(species,FamilyCommonName))
SpeciesFamilyCommonName = unique(SpeciesFamilyCommonName)

#Merge with Species and SampleSizesCanAssess
Species = merge(Species,SpeciesFamilyCommonName,all.x=TRUE,by=c("species"))
SampleSizesCanAssess = merge(SampleSizesCanAssess,SpeciesFamilyCommonName,all.x=TRUE,by=c("species"))


MeanLengthEachYear = aggregate.data.frame(Lengths$cm,by=list(Lengths$species,Lengths$WPP,Lengths$Year_OutFishing),FUN=mean,na.rm=TRUE)
names(MeanLengthEachYear) = c("species","WPP","Year","meanLength")
MeanLengthEachYear = MeanLengthEachYear[order(MeanLengthEachYear$species,MeanLengthEachYear$WPP,MeanLengthEachYear$Year),]
write.table(MeanLengthEachYear,paste(PATH_output,"MeanLengthEachYearWPP.csv",sep="/"),sep="/",col.names=TRUE,row.names=FALSE)

############################################################################################
######### Plot Length Histograms With Ideal Bins and Lognormal Curve #######################
############################################################################################
if(plotFittedHistograms==TRUE)
{
	OptimalHistogramBreaks = read.table(paste(PATH_otherInput,histogramBinSizes_fName,sep="/"),header=TRUE,sep=",")
	names(OptimalHistogramBreaks) = gsub("X","",names(OptimalHistogramBreaks))
	WPP_EEZ = c(as.character(WPP),"EEZ")
	for(i in 1:length(GenusSpeciesForAnalysis_names))
	{
		for(j in 1:length(WPP_EEZ))
		{
			OptimalHistogramBreaks_i = subset(OptimalHistogramBreaks,OptimalHistogramBreaks$species==GenusSpeciesForAnalysis_names[i])
			Breaks_i = OptimalHistogramBreaks_i[,names(OptimalHistogramBreaks_i)==WPP_EEZ[j]]
			if(!is.na(Breaks_i))
			{
				if(WPP_EEZ[j]!="EEZ")
				{
					Sub = subset(Lengths,Lengths$species==GenusSpeciesForAnalysis_names[i] & as.character(Lengths$WPP)==WPP_EEZ[j])		
				}
				if(WPP_EEZ[j]=="EEZ")
				{			
					Sub = subset(Lengths,Lengths$species==GenusSpeciesForAnalysis_names[i])		
				}
				Breaks = seq(0,ceiling(max(Sub$cm+Breaks_i)),by=Breaks_i)	
				png(paste(PATH_histograms,paste("Histogram_",gsub(" ","_",GenusSpeciesForAnalysis_names[i]),"_",WPP_EEZ[j],".png",sep=""),sep="/"),units="px",width=3200,height=4200,res=600)
				Hist = hist(Sub$cm,xlab="Fish Length (cm)",main=GenusSpeciesForAnalysis_names[i],breaks=Breaks)
				mtext(paste("FMA=",WPP_EEZ[j],sep=""))
				DF = data.frame(Mids=Hist$mids,Count=Hist$counts)
				ExpandedDF = as.data.frame(lapply(DF,rep,DF$Count))
				ExpandedDF$Count=NULL			
				DistModel = fitdistr(ExpandedDF$Mids,densfun="lognormal")
				xfit=Hist$mids
				yfit=dlnorm(xfit,mean=DistModel[[1]][[1]],sd=DistModel[[1]][[2]])							
				multiplyer = (max(Hist$counts))/max(yfit)
				yfit = yfit*multiplyer
				lines(x=xfit,y=yfit, type="l",col="red", lwd=3)
				dev.off()
			}	
		}	
		if(i%%5==0)
		{
			print(paste("Working on item ",i," of ",length(GenusSpeciesForAnalysis_names),"...",sep=""))
			flush.console()			
		}
	}
}

#############################################################################################
############ Setup, Calculate and Pull Life History Parameters From Online Databases ########
#############################################################################################

if(pullAndEstimateLifeHistoryParams==TRUE)
{
	#First, setup generic, standardized data.frame to capture all the versions of life history parameters.
	#Then, put each data.frame into a list object. We create that empty list object here
	LifeHistoryScenariosList = list()
	#Set up generic data.frame with slots for each species to capture life history values 
	LifeHistoryParametersDataFrame_Generic = data.frame(species=GenusSpeciesForAnalysis_names,Lmax=-999,Linf=-999,Lmat=-999,Lopt=-999,Temperature_min=-999,Temperature_max=-999,Temperature_mean=-999,M=-999,
		VBK=-999,TrophicLevel=-999,TrophicLevelSE=-999,t0=-999,Age_max=-999,var_a=-999,var_b=-999,Winf=-999,
		h_steepness=-999,generation_time_years=-999,r_populationGrowthRate=-999,Fmsy_FishLife=-999,DepthLimit_deep=-999,DepthLimit_shallow=-999,DepthLimit_mean=-999,Vulnerability=-999,
		Lmax_imputed=FALSE,VBK_imputed=FALSE,M_imputed=FALSE,Linf_imputed=FALSE,Linf_assumedLmax=FALSE,Lmat_imputed=FALSE,
		Lopt_imputed=FALSE,Winf_imputed=FALSE,Age_max_imputed=FALSE,WaterTempAssumed=FALSE,WaterTempFromFishLife=FALSE,t0_assumed=FALSE,var_a_soureThisDataset=FALSE,var_b_soureThisDataset=FALSE,
		h_steepness_imputed=FALSE,generation_time_years_imputed=FALSE,r_populationGrowthRate_imputed=FALSE)
	FishBaseFecundityValues = data.frame(species=GenusSpeciesForAnalysis_names,FecundityMean=-999,FecundityMin=-999,
		FecundityMax=-999,LengthFecunMin=-999,LengthFecunMax=-999)
	LifeHistoryParametersDataFrame_Generic[LifeHistoryParametersDataFrame_Generic==-999]=NA
	#Scenario 1: Pull FishBase Values Directly Where They are Available
	library(rfishbase)
	LifeHistoryParametersDataFrame = LifeHistoryParametersDataFrame_Generic
	for(i in 1:dim(LifeHistoryParametersDataFrame)[1])
	{
		#popll has conversion from SL to TL and FL to TL, etc. Length1 is what you get if you have Length2. The formula is Length1 = a + b*Length2
		Length_Length = data.frame(popll(species_list=LifeHistoryParametersDataFrame$species[i]))
		Length_Length = subset(Length_Length,select=c(Length1,Length2,a,b))
		Length_Length$a = as.numeric(Length_Length$a)
		Length_Length$b = as.numeric(Length_Length$b)	
		nms = names(Length_Length)
		Length_Length = aggregate.data.frame(list(Length_Length$a,Length_Length$b),by=list(Length_Length$Length1,Length_Length$Length2),FUN=mean,na.rm=TRUE)
		names(Length_Length) = nms	
		#length_weight has values to calculate the length-weight relationship; Type is FL, TL or SL; 
		#IMPORTANT: a and b are Length-weight relationship values that correspond to the length type ("Type") listed!!!!
		LengthWeight = data.frame(length_weight(species_list=LifeHistoryParametersDataFrame$species[i]))
		LengthWeight = subset(LengthWeight,select=c(Type,a,aTL,b))
		#Contains depth, temperature, Lmat as Lm, Linf as Loo, Lopt, Lc, a and b weight-length conversion values
		PopLW = data.frame(poplw(species_list=LifeHistoryParametersDataFrame$species[i]))
		PopLW = subset(PopLW,select=c(Type,a,aTL,b))
		#Contains DepthMin, DepthMax,TempMin,TempMax,Lm (length at maturity), Loo, Lopt, Lc, Length-weight a and b parameters
		PopLf = data.frame(poplf(species_list=LifeHistoryParametersDataFrame$species[i]))
		PopLf = subset(PopLf,select=c(DepthMin,DepthMax,TempMin,TempMax,Lm,Loo,Lopt,Lc,a,b))
		PopLf_weightLength = data.frame(Type="TL",a=PopLf$a,aTL=NA,b=PopLf$b)
		LengthWeight = rbind(LengthWeight,PopLW,PopLf_weightLength)
		LengthWeight = unique(LengthWeight)
		LengthWeight = subset(LengthWeight,!is.na(LengthWeight$a) & !is.na(LengthWeight$b))	
		LengthWeight$a[!is.na(LengthWeight$aTL)]=LengthWeight$aTL[!is.na(LengthWeight$aTL)]
		LengthWeight$Type[!is.na(LengthWeight$aTL)]="TL"
		LengthWeight$aTL=NULL
		LengthWeight$Type[is.na(LengthWeight$Type)]="TL"  #If length type missing, assume it is Total Length!
		names(LengthWeight)[names(LengthWeight)=="Type"]="Type_LW"
		#Loo and TLinfinity are Linf, K is VBK, to is t0, Winfinity is Winf,Temperature which is mean temperature, M is natural mortality
		Popgrowth_fishBase = data.frame(popgrowth(species_list=LifeHistoryParametersDataFrame$species[i]))
		Popgrowth_fishBase = subset(Popgrowth_fishBase,select=c(Loo,Type,TLinfinity,K,to,Temperature,M))
		Popgrowth_fishBase$Loo[!is.na(Popgrowth_fishBase$TLinfinity)]=Popgrowth_fishBase$TLinfinity[!is.na(Popgrowth_fishBase$TLinfinity)]
		Popgrowth_fishBase$Type[!is.na(Popgrowth_fishBase$TLinfinity)]="TL"
		Popgrowth_fishBase_LooNA_TLinfinity_notNA = subset(Popgrowth_fishBase,is.na(Popgrowth_fishBase$Loo) & !is.na(Popgrowth_fishBase$TLinfinity))
		Popgrowth_fishBase$TLinfinity=NULL
		Popgrowth_fishBase_LooNA_TLinfinity_notNA$Loo=NULL
		names(Popgrowth_fishBase_LooNA_TLinfinity_notNA)[names(Popgrowth_fishBase_LooNA_TLinfinity_notNA)=="TLinfinity"]="Loo"
		Popgrowth_fishBase_LooNA_TLinfinity_notNA$TLinfinity=NULL	
		Popgrowth_fishBase_LooNA_TLinfinity_notNA = Popgrowth_fishBase_LooNA_TLinfinity_notNA[names(Popgrowth_fishBase)]
		Popgrowth_fishBase = rbind(Popgrowth_fishBase,Popgrowth_fishBase_LooNA_TLinfinity_notNA)
		#species pulls stuff from home page of FishBase. Length is Lmax and CommonLength is Lc. Weight is Wmax in grams
		SppHomePage = data.frame(species(species_list=LifeHistoryParametersDataFrame$species[i]))
		SppHomePage = subset(SppHomePage,select=c(DepthRangeShallow,DepthRangeDeep,Vulnerability,Length,LTypeMaxM,CommonLength,LTypeComM))
		#Get trophic level and combine with species home page information (both seem to only be one row so should combine seamlessly in theory...we'll see);
		Ecology = data.frame(ecology(species_list=LifeHistoryParametersDataFrame$species[i]))
		Ecology = subset(Ecology,select=c(FoodTroph,FoodSeTroph))
		SppHomePage = cbind(SppHomePage,Ecology)	
		#Wmax, TypeWeigt are units of Wmax, Lmax, Type are units of Lmax, tmax is longevity (Age_max) in years
		PopChar_fishBase = data.frame(popchar(species_list=LifeHistoryParametersDataFrame$species[i]))	
		PopChar_fishBase = subset(PopChar_fishBase,select=c(Lmax,Type,tmax))
		#Length at Maturity
		Mat_fishBase = data.frame(maturity(species_list=LifeHistoryParametersDataFrame$species[i]))
		Mat_fishBase = subset(Mat_fishBase,select=c(Lm,Type1,LengthMatMin,LengthMatMin2))
		Mat_fishBase$Lm[is.na(Mat_fishBase$Lm)]=((Mat_fishBase$LengthMatMin[is.na(Mat_fishBase$Lm)]+Mat_fishBase$LengthMatMin2[is.na(Mat_fishBase$Lm)])/2)
		Mat_fishBase$LengthMatMin=NULL
		Mat_fishBase$LengthMatMin2=NULL
		#Fecundity	
		Fecundity = data.frame(fecundity(species_list=LifeHistoryParametersDataFrame$species[i]))
		Fecundity = subset(Fecundity,select=c(FecundityMean,FecundityMin,FecundityMax))
		#Spawning Values (probably duplicated from Fecundity table, but just in case) 
		Spawning = data.frame(spawning(species_list=LifeHistoryParametersDataFrame$species[i]))
		Spawning = subset(Spawning,select=c(FecundityMin,FecundityMax,LengthFecunMin,LengthFecunMax))
		#Pull all parameters together into one data.frame	
		Parameters = data.frame(Loo=Popgrowth_fishBase$Loo,Type_Loo=Popgrowth_fishBase$Type,K=Popgrowth_fishBase$K,
			to=Popgrowth_fishBase$to,Temperature=Popgrowth_fishBase$Temperature,M=Popgrowth_fishBase$M)	
		#Add LengthWeight to Parameters
		rowDiff = dim(Parameters)[1]- dim(LengthWeight)[1]
		if(rowDiff<0)
		{
			temp = data.frame(matrix(nrow=abs(rowDiff),ncol=(dim(Parameters)[2]))) 
			names(temp) = names(Parameters)
			Parameters = rbind(Parameters,temp)
			Parameters = cbind(Parameters,LengthWeight)
		}
		if(rowDiff==0)
		{
			Parameters = cbind(Parameters,LengthWeight)	
		}
		if(rowDiff>0)
		{
			temp = data.frame(matrix(nrow=abs(rowDiff),ncol=(dim(LengthWeight)[2]))) 
			names(temp) = names(LengthWeight)
			LengthWeight = rbind(LengthWeight,temp)
			Parameters = cbind(Parameters,LengthWeight)
		}
		#Combine FishBase home page information with population characeristics, then add it to Parameters	
		names(SppHomePage)[names(SppHomePage)=="DepthRangeShallow"]="DepthMin"
		names(SppHomePage)[names(SppHomePage)=="DepthRangeDeep"]="DepthMax"
		names(SppHomePage)[names(SppHomePage)=="Length"]="Lmax"
		names(SppHomePage)[names(SppHomePage)=="LTypeMaxM"]="Type_Lmax"
		names(SppHomePage)[names(SppHomePage)=="CommonLength"]="Lc"
		names(SppHomePage)[names(SppHomePage)=="LTypeComM"]="Type_Lc"
		names(PopChar_fishBase)[names(PopChar_fishBase)=="Type"]="Type_Lmax"
		SppHomePage_Lmax = subset(SppHomePage,select=c(Lmax,Type_Lmax))
		SppHomePage_Lmax$tmax=NA
		PopChar_fishBase = rbind(PopChar_fishBase,SppHomePage_Lmax)
		SppHomePage$Lmax=NULL
		SppHomePage$Type_Lmax=NULL	
		rowDiff = dim(SppHomePage)[1]- dim(PopChar_fishBase)[1]
		if(rowDiff<0)
		{
			temp = data.frame(matrix(nrow=abs(rowDiff),ncol=(dim(SppHomePage)[2]))) 
			names(temp) = names(SppHomePage)
			SppHomePage = rbind(SppHomePage,temp)
			SppHomePage = cbind(SppHomePage,PopChar_fishBase)
		}
		if(rowDiff==0)
		{
			SppHomePage = cbind(SppHomePage,PopChar_fishBase)
		}
		if(rowDiff>0)
		{
			temp = data.frame(matrix(nrow=abs(rowDiff),ncol=(dim(PopChar_fishBase)[2]))) 
			names(temp) = names(PopChar_fishBase)
			PopChar_fishBase = rbind(PopChar_fishBase,temp)
			SppHomePage = cbind(SppHomePage,PopChar_fishBase)
		}		
		rowDiff = dim(SppHomePage)[1]- dim(Parameters)[1]
		if(rowDiff<0)
		{
			temp = data.frame(matrix(nrow=abs(rowDiff),ncol=(dim(SppHomePage)[2]))) 
			names(temp) = names(SppHomePage)
			SppHomePage = rbind(SppHomePage,temp)
			Parameters = cbind(Parameters,SppHomePage)
		}
		if(rowDiff==0)
		{
			Parameters = cbind(Parameters,SppHomePage)
		}
		if(rowDiff>0)
		{
			temp = data.frame(matrix(nrow=abs(rowDiff),ncol=(dim(Parameters)[2]))) 
			names(temp) = names(Parameters)
			Parameters = rbind(Parameters,temp)
			Parameters = cbind(Parameters,SppHomePage)
		}
		#Combine Length at maturity with Fecundity and Spawning, then add these to Parameters
		Fecundity$LengthFecunMin=NA
		Fecundity$LengthFecunMax=NA
		Spawning$FecundityMean=NA
		Spawning=Spawning[names(Fecundity)]
		SpawningFecundity = rbind(Spawning,Fecundity)
		rowDiff = dim(SpawningFecundity)[1]- dim(Mat_fishBase)[1]
		if(rowDiff<0)
		{
			temp = data.frame(matrix(nrow=abs(rowDiff),ncol=(dim(SpawningFecundity)[2]))) 
			names(temp) = names(SpawningFecundity)
			SpawningFecundity = rbind(SpawningFecundity,temp)
			Mat_fishBase = cbind(Mat_fishBase,SpawningFecundity)
		}
		if(rowDiff==0)
		{
			Mat_fishBase = cbind(Mat_fishBase,SpawningFecundity)
		}
		if(rowDiff>0)
		{
			temp = data.frame(matrix(nrow=abs(rowDiff),ncol=(dim(Mat_fishBase)[2]))) 
			names(temp) = names(Mat_fishBase)
			Mat_fishBase = rbind(Mat_fishBase,temp)
			Mat_fishBase = cbind(Mat_fishBase,SpawningFecundity)
		}
		SpawningFecundity$Lmat = Mat_fishBase$Lm
		SpawningFecundity$Type_Lmat = Mat_fishBase$Type1
		#Now add SpawningFecundity to Parameters 
		rowDiff = dim(Parameters)[1]- dim(SpawningFecundity)[1]
		if(rowDiff<0)
		{
			temp = data.frame(matrix(nrow=abs(rowDiff),ncol=(dim(Parameters)[2]))) 
			names(temp) = names(Parameters)
			Parameters = rbind(Parameters,temp)
			Parameters = cbind(Parameters,SpawningFecundity)
		}
		if(rowDiff==0)
		{
			Parameters = cbind(Parameters,SpawningFecundity)	
		}
		if(rowDiff>0)
		{
			temp = data.frame(matrix(nrow=abs(rowDiff),ncol=(dim(SpawningFecundity)[2]))) 
			names(temp) = names(SpawningFecundity)
			SpawningFecundity = rbind(SpawningFecundity,temp)
			Parameters = cbind(Parameters,SpawningFecundity)
		}
		Parameters$Loo = as.numeric(Parameters$Loo)
		Parameters$K = as.numeric(Parameters$K)
		Parameters$to = as.numeric(Parameters$to)
		Parameters$Temperature = as.numeric(Parameters$Temperature)
		Parameters$M = as.numeric(Parameters$M)
		Parameters$a = as.numeric(Parameters$a)
		Parameters$b = as.numeric(Parameters$b)
		Parameters$DepthMin = as.numeric(Parameters$DepthMin)
		Parameters$DepthMax = as.numeric(Parameters$DepthMax)
		Parameters$Vulnerability = as.numeric(Parameters$Vulnerability)
		Parameters$Lc = as.numeric(Parameters$Lc)
		Parameters$FoodTroph = as.numeric(Parameters$FoodTroph)
		Parameters$FoodSeTroph = as.numeric(Parameters$FoodSeTroph)
		Parameters$Lmax = as.numeric(Parameters$Lmax)
		Parameters$tmax = as.numeric(Parameters$tmax)
		Parameters$FecundityMean = as.numeric(Parameters$FecundityMean)
		Parameters$FecundityMin = as.numeric(Parameters$FecundityMin)
		Parameters$FecundityMax = as.numeric(Parameters$FecundityMax)
		Parameters$LengthFecunMin = as.numeric(Parameters$LengthFecunMin)
		Parameters$LengthFecunMax = as.numeric(Parameters$LengthFecunMax)
		Parameters$Lmat = as.numeric(Parameters$Lmat)
		#Set any length types that don't have the label TL, SL, FL or are NA to NA
		Parameters$Type_Loo[Parameters$Type_Loo!="TL" | Parameters$Type_Loo!="FL" | Parameters$Type_Loo!="SL" | !is.na(Parameters$Type_Loo)]=NA
		Parameters$Type_LW[Parameters$Type_LW!="TL" | Parameters$Type_LW!="FL" | Parameters$Type_LW!="SL" | !is.na(Parameters$Type_LW)]=NA
		Parameters$Type_Lc[Parameters$Type_Lc!="TL" | Parameters$Type_Lc!="FL" | Parameters$Type_Lc!="SL" | !is.na(Parameters$Type_Lc)]=NA
		Parameters$Type_Lmax[Parameters$Type_Lmax!="TL" | Parameters$Type_Lmax!="FL" | Parameters$Type_Lmax!="SL" | !is.na(Parameters$Type_Lmax)]=NA	
		#Pass Fecundity Values to Fecundity data.frame
		if(all(is.na(Parameters$FecundityMean))==FALSE)
		{
			FishBaseFecundityValues$FecundityMean[i] = mean(Parameters$FecundityMean,na.rm=TRUE)
		} else {
			FishBaseFecundityValues$FecundityMean[i] = NA
		}		
		if(all(is.na(Parameters$FecundityMin))==FALSE)
		{
			FishBaseFecundityValues$FecundityMin[i] = mean(Parameters$FecundityMin,na.rm=TRUE)
		} else {
			FishBaseFecundityValues$FecundityMin[i] = NA
		}		
		if(all(is.na(Parameters$FecundityMax))==FALSE)
		{
			FishBaseFecundityValues$FecundityMax[i] = mean(Parameters$FecundityMax,na.rm=TRUE)
		} else {
			FishBaseFecundityValues$FecundityMax[i] = NA
		}		
		if(all(is.na(Parameters$LengthFecunMin))==FALSE)
		{
			FishBaseFecundityValues$LengthFecunMin[i] = mean(Parameters$LengthFecunMin,na.rm=TRUE)
		} else {
			FishBaseFecundityValues$LengthFecunMin[i] = NA
		}		
			if(all(is.na(Parameters$LengthFecunMax))==FALSE)
		{
			FishBaseFecundityValues$LengthFecunMax[i] = mean(Parameters$LengthFecunMax,na.rm=TRUE)
		} else {
			FishBaseFecundityValues$LengthFecunMax[i] = NA
		}
		#Convert all length values to TL
		Length_Lenth_FL2TL = subset(Length_Length,Length_Length$Length2=="FL" & Length_Length$Length1=="TL")
		if(dim(Length_Lenth_FL2TL)[1]>0)
		{
			if(!is.na(Length_Lenth_FL2TL$a) & !is.na(Length_Lenth_FL2TL$b))
			{
				Parameters$Loo[Parameters$Type_Loo=="FL" & !is.na(Parameters$Loo) & !is.na(Parameters$Type_Loo)] = Length_Lenth_FL2TL$a + Length_Lenth_FL2TL$b * Parameters$Loo[Parameters$Type_Loo=="FL" & !is.na(Parameters$Loo) & !is.na(Parameters$Type_Loo)]
				Parameters$Type_Loo[Parameters$Type_Loo=="FL" & !is.na(Parameters$Loo) & !is.na(Parameters$Type_Loo)]="TL"
				Parameters$Lc[Parameters$Type_Lc=="FL" & !is.na(Parameters$Lc) & !is.na(Parameters$Type_Lc)] = Length_Lenth_FL2TL$a + Length_Lenth_FL2TL$b * Parameters$Lc[Parameters$Type_Lc=="FL" & !is.na(Parameters$Lc) & !is.na(Parameters$Type_Lc)]
				Parameters$Type_Lc[Parameters$Type_Lc=="FL" & !is.na(Parameters$Lc) & !is.na(Parameters$Type_Lc)]="TL"
				Parameters$Lmax[Parameters$Type_Lmax=="FL" & !is.na(Parameters$Lmax) & !is.na(Parameters$Type_Lmax)] = Length_Lenth_FL2TL$a + Length_Lenth_FL2TL$b * Parameters$Lmax[Parameters$Type_Lmax=="FL" & !is.na(Parameters$Lmax) & !is.na(Parameters$Type_Lmax)]
				Parameters$Type_Lmax[Parameters$Type_Lmax=="FL" & !is.na(Parameters$Lmax) & !is.na(Parameters$Type_Lmax)]="TL"
				Parameters$Lmat[Parameters$Type_Lmat=="FL" & !is.na(Parameters$Lmat) & !is.na(Parameters$Type_Lmat)] = Length_Lenth_FL2TL$a + Length_Lenth_FL2TL$b * Parameters$Lmat[Parameters$Type_Lmat=="FL" & !is.na(Parameters$Lmat) & !is.na(Parameters$Type_Lmat)]
				Parameters$Type_Lmat[Parameters$Type_Lmat=="FL" & !is.na(Parameters$Lmat) & !is.na(Parameters$Type_Lmat)]="TL"			
			}
		}
		Length_Lenth_TL2FL = subset(Length_Length,Length_Length$Length2=="TL" & Length_Length$Length1=="FL")
		if(dim(Length_Lenth_FL2TL)[1]==0 & dim(Length_Lenth_TL2FL)[1]>0)
		{
			if(!is.na(Length_Lenth_TL2FL$a) & !is.na(Length_Lenth_TL2FL$b))
			{
				Parameters$Loo[Parameters$Type_Loo=="FL" & !is.na(Parameters$Loo) & !is.na(Parameters$Type_Loo)] = (Parameters$Loo[Parameters$Type_Loo=="FL" & !is.na(Parameters$Loo) & !is.na(Parameters$Type_Loo)] - Length_Lenth_FL2TL$a)/Length_Lenth_FL2TL$b 
				Parameters$Type_Loo[Parameters$Type_Loo=="FL" & !is.na(Parameters$Loo) & !is.na(Parameters$Type_Loo)]="TL"
				Parameters$Lc[Parameters$Type_Lc=="FL" & !is.na(Parameters$Lc) & !is.na(Parameters$Type_Lc)] = (Parameters$Lc[Parameters$Type_Lc=="FL" & !is.na(Parameters$Lc) & !is.na(Parameters$Type_Lc)] - Length_Lenth_FL2TL$a)/Length_Lenth_FL2TL$b 
				Parameters$Type_Lc[Parameters$Type_Lc=="FL" & !is.na(Parameters$Lc) & !is.na(Parameters$Type_Lc)]="TL"
				Parameters$Lmax[Parameters$Type_Lmax=="FL" & !is.na(Parameters$Lmax) & !is.na(Parameters$Type_Lmax)] = (Parameters$Lmax[Parameters$Type_Lmax=="FL" & !is.na(Parameters$Lmax) & !is.na(Parameters$Type_Lmax)] - Length_Lenth_FL2TL$a)/Length_Lenth_FL2TL$b
				Parameters$Type_Lmax[Parameters$Type_Lmax=="FL" & !is.na(Parameters$Lmax) & !is.na(Parameters$Type_Lmax)]="TL"
				Parameters$Lmat[Parameters$Type_Lmat=="FL" & !is.na(Parameters$Lmat) & !is.na(Parameters$Type_Lmat)] = (Parameters$Lmat[Parameters$Type_Lmat=="FL" & !is.na(Parameters$Lmat) & !is.na(Parameters$Type_Lmat)] - Length_Lenth_FL2TL$a)/Length_Lenth_FL2TL$b
				Parameters$Type_Lmat[Parameters$Type_Lmat=="FL" & !is.na(Parameters$Lmat) & !is.na(Parameters$Type_Lmat)]="TL"
			}
		}
		Length_Lenth_SL2TL = subset(Length_Length,Length_Length$Length2=="SL" & Length_Length$Length1=="TL")
		if(dim(Length_Lenth_SL2TL)[1]>0)
		{
			if(!is.na(Length_Lenth_SL2TL$a) & !is.na(Length_Lenth_SL2TL$b))
			{
				Parameters$Loo[Parameters$Type_Loo=="SL" & !is.na(Parameters$Loo) & !is.na(Parameters$Type_Loo)] = Length_Lenth_SL2TL$a + Length_Lenth_SL2TL$b * Parameters$Loo[Parameters$Type_Loo=="SL" & !is.na(Parameters$Loo) & !is.na(Parameters$Type_Loo)]
				Parameters$Type_Loo[Parameters$Type_Loo=="SL" & !is.na(Parameters$Loo) & !is.na(Parameters$Type_Loo)]="TL"
				Parameters$Lc[Parameters$Type_Lc=="SL" & !is.na(Parameters$Lc) & !is.na(Parameters$Type_Lc)] = Length_Lenth_SL2TL$a + Length_Lenth_SL2TL$b * Parameters$Lc[Parameters$Type_Lc=="SL" & !is.na(Parameters$Lc) & !is.na(Parameters$Type_Lc)]
				Parameters$Type_Lc[Parameters$Type_Lc=="SL" & !is.na(Parameters$Lc) & !is.na(Parameters$Type_Lc)]="TL"
				Parameters$Lmax[Parameters$Type_Lmax=="SL" & !is.na(Parameters$Lmax) & !is.na(Parameters$Type_Lmax)] = Length_Lenth_SL2TL$a + Length_Lenth_SL2TL$b * Parameters$Lmax[Parameters$Type_Lmax=="SL" & !is.na(Parameters$Lmax) & !is.na(Parameters$Type_Lmax)]
				Parameters$Type_Lmax[Parameters$Type_Lmax=="SL" & !is.na(Parameters$Lmax) & !is.na(Parameters$Type_Lmax)]="TL"
				Parameters$Lmat[Parameters$Type_Lmat=="SL" & !is.na(Parameters$Lmat) & !is.na(Parameters$Type_Lmat)] = Length_Lenth_SL2TL$a+ Length_Lenth_SL2TL$b * Parameters$Lmat[Parameters$Type_Lmat=="SL" & !is.na(Parameters$Lmat) & !is.na(Parameters$Type_Lmat)]
				Parameters$Type_Lmat[Parameters$Type_Lmat=="SL" & !is.na(Parameters$Lmat) & !is.na(Parameters$Type_Lmat)]="TL"
			}
		}
		Length_Lenth_TL2SL = subset(Length_Length,Length_Length$Length2=="TL" & Length_Length$Length1=="SL")
		if(dim(Length_Lenth_SL2TL)[1]==0 & dim(Length_Lenth_TL2SL)[1]>0)
		{
			if(!is.na(Length_Lenth_TL2SL$a) & !is.na(Length_Lenth_TL2SL$b))
			{
				Parameters$Loo[Parameters$Type_Loo=="SL" & !is.na(Parameters$Loo) & !is.na(Parameters$Type_Loo)] = (Parameters$Loo[Parameters$Type_Loo=="SL" & !is.na(Parameters$Loo) & !is.na(Parameters$Type_Loo)] - Length_Lenth_TL2SL$a)/Length_Lenth_TL2SL$b
				Parameters$Type_Loo[Parameters$Type_Loo=="SL" & !is.na(Parameters$Loo) & !is.na(Parameters$Type_Loo)]="TL"
				Parameters$Lc[Parameters$Type_Lc=="SL" & !is.na(Parameters$Lc) & !is.na(Parameters$Type_Lc)] = (Parameters$Lc[Parameters$Type_Lc=="SL" & !is.na(Parameters$Lc) & !is.na(Parameters$Type_Lc)] - Length_Lenth_TL2SL$a)/Length_Lenth_TL2SL$b 
				Parameters$Type_Lc[Parameters$Type_Lc=="SL" & !is.na(Parameters$Lc) & !is.na(Parameters$Type_Lc)]="TL"
				Parameters$Lmax[Parameters$Type_Lmax=="SL" & !is.na(Parameters$Lmax) & !is.na(Parameters$Type_Lmax)] = (Parameters$Lmax[Parameters$Type_Lmax=="SL" & !is.na(Parameters$Lmax) & !is.na(Parameters$Type_Lmax)] - Length_Lenth_TL2SL$a)/Length_Lenth_TL2SL$b
				Parameters$Type_Lmax[Parameters$Type_Lmax=="SL" & !is.na(Parameters$Lmax) & !is.na(Parameters$Type_Lmax)]="TL"
				Parameters$Lmat[Parameters$Type_Lmat=="SL" & !is.na(Parameters$Lmat) & !is.na(Parameters$Type_Lmat)] = (Parameters$Lmat[Parameters$Type_Lmat=="SL" & !is.na(Parameters$Lmat) & !is.na(Parameters$Type_Lmat)] - Length_Lenth_TL2SL$a)/Length_Lenth_TL2SL$b
				Parameters$Type_Lmat[Parameters$Type_Lmat=="SL" & !is.na(Parameters$Lmat) & !is.na(Parameters$Type_Lmat)]="TL"
			}
		}	
		#Convert Length-Weight a and b parameters to all be expressed for Total Length
		Sub = subset(Lengths,Lengths$species==LifeHistoryParametersDataFrame$species[i])
		LengthTypeWeight = data.frame(Length=c(mean(Sub$cm,na.rm=TRUE),NA,NA),Type=c("TL","SL","FL"))
		Length_Length_SL2TL = subset(Length_Length,Length_Length$Length2=="TL" & Length_Length$Length1=="SL")
		if(dim(Length_Length_SL2TL)[1]>0)
		{
			LengthTypeWeight$Length[LengthTypeWeight$Type=="SL" & !is.na(LengthTypeWeight$Type)] = Length_Length_SL2TL$a + Length_Length_SL2TL$b * LengthTypeWeight$Length[LengthTypeWeight$Type=="TL" & !is.na(LengthTypeWeight$Type)]
		}
		Length_Length_TL2SL = subset(Length_Length,Length_Length$Length2=="SL" & Length_Length$Length1=="TL")
		if(dim(Length_Length_SL2TL)[1]==0 & dim(Length_Length_TL2SL)[1]>0)
		{	
			LengthTypeWeight$Length[LengthTypeWeight$Type=="SL" & !is.na(LengthTypeWeight$Type)] = (LengthTypeWeight$Length[LengthTypeWeight$Type=="TL" & !is.na(LengthTypeWeight$Type)] - Length_Length_TL2SL$a)/Length_Length_TL2SL$b
		}	
		Length_Length_FL2TL = subset(Length_Length,Length_Length$Length2=="TL" & Length_Length$Length1=="FL")
		if(dim(Length_Length_FL2TL)[1]>0)
		{
			LengthTypeWeight$Length[LengthTypeWeight$Type=="FL" & !is.na(LengthTypeWeight$Type)] = Length_Length_FL2TL$a + Length_Length_FL2TL$b * LengthTypeWeight$Length[LengthTypeWeight$Type=="TL" & !is.na(LengthTypeWeight$Type)]
		}
		Length_Length_TL2FL = subset(Length_Length,Length_Length$Length2=="FL" & Length_Length$Length1=="TL")
		if(dim(Length_Length_FL2TL)[1]==0 & dim(Length_Length_TL2FL)[1]>0)
		{	
			LengthTypeWeight$Length[LengthTypeWeight$Type=="FL" & !is.na(LengthTypeWeight$Type)] = (LengthTypeWeight$Length[LengthTypeWeight$Type=="TL" & !is.na(LengthTypeWeight$Type)] - Length_Length_TL2FL$a)/Length_Length_TL2FL$b
		}
		LengthWeight_ParamsAndTypes = subset(Parameters,select=c(Type_LW,a,b))	
		names(LengthWeight_ParamsAndTypes)[names(LengthWeight_ParamsAndTypes)=="Type_LW"]="Type"
		LengthWeight_ParamsAndTypes = merge(LengthWeight_ParamsAndTypes,LengthTypeWeight,all.x=TRUE,by=c("Type"))
		LengthWeight_ParamsAndTypes$W_estimated = LengthWeight_ParamsAndTypes$a + (LengthWeight_ParamsAndTypes$Length * LengthWeight_ParamsAndTypes$b)
		LengthWeight_ParamsAndTypes$a_new=NA
		#Remember: Length1 is what you get if you have Length2
		Length_Length_SL2TL = subset(Length_Length,Length_Length$Length2=="SL" & Length_Length$Length1=="TL")
		if(all(is.na(LengthWeight_ParamsAndTypes$Type))==FALSE)
		{
			if(dim(Length_Length_SL2TL)[1]>0)
			{
				LengthWeight_ParamsAndTypes$a_new[LengthWeight_ParamsAndTypes$Type=="SL" & !is.na(LengthWeight_ParamsAndTypes$Type)] = LengthWeight_ParamsAndTypes$W_estimated[LengthWeight_ParamsAndTypes$Type=="SL" & !is.na(LengthWeight_ParamsAndTypes$Type)]/
					(Length_Length_SL2TL$a + (Length_Length_SL2TL$b * mean(Sub$cm,na.rm=TRUE)))^LengthWeight_ParamsAndTypes$b[LengthWeight_ParamsAndTypes$Type=="SL" & !is.na(LengthWeight_ParamsAndTypes$Type)]
				LengthWeight_ParamsAndTypes$Type[LengthWeight_ParamsAndTypes$Type=="SL" & !is.na(LengthWeight_ParamsAndTypes$Type)]="TL"
			}
			Length_Length_TL2SL = subset(Length_Length,Length_Length$Length2=="TL" & Length_Length$Length1=="SL")
			if(dim(Length_Length_SL2TL)[1]==0 & dim(Length_Length_TL2SL)[1]>0)
			{	
				LengthWeight_ParamsAndTypes$a_new[LengthWeight_ParamsAndTypes$Type=="SL" & !is.na(LengthWeight_ParamsAndTypes$Type)] = LengthWeight_ParamsAndTypes$W_estimated[LengthWeight_ParamsAndTypes$Type=="SL" & !is.na(LengthWeight_ParamsAndTypes$Type)]/
					((mean(Sub$cm,na.rm=TRUE) - Length_Length_TL2SL$a)/Length_Length_TL2SL$b)
				LengthWeight_ParamsAndTypes$Type[LengthWeight_ParamsAndTypes$Type=="SL" & !is.na(LengthWeight_ParamsAndTypes$Type)]="TL"
			}		
			Length_Length_FL2TL = subset(Length_Length,Length_Length$Length2=="FL" & Length_Length$Length1=="TL")
			if(dim(Length_Length_FL2TL)[1]>0)
			{
				LengthWeight_ParamsAndTypes$a_new[LengthWeight_ParamsAndTypes$Type=="FL" & !is.na(LengthWeight_ParamsAndTypes$Type)] = LengthWeight_ParamsAndTypes$W_estimated[LengthWeight_ParamsAndTypes$Type=="FL" & !is.na(LengthWeight_ParamsAndTypes$Type)]/
					(Length_Length_FL2TL$a + (Length_Length_FL2TL$b * mean(Sub$cm,na.rm=TRUE)))^LengthWeight_ParamsAndTypes$b[LengthWeight_ParamsAndTypes$Type=="FL" & !is.na(LengthWeight_ParamsAndTypes$Type)]
				LengthWeight_ParamsAndTypes$Type[LengthWeight_ParamsAndTypes$Type=="FL" & !is.na(LengthWeight_ParamsAndTypes$Type)]="TL"
			}
			Length_Length_TL2FL = subset(Length_Length,Length_Length$Length2=="TL" & Length_Length$Length1=="FL")
			if(dim(Length_Length_FL2TL)[1]==0 & dim(Length_Length_TL2FL)[1]>0)
			{	
				LengthWeight_ParamsAndTypes$a_new[LengthTypeWeight$Type=="FL" & !is.na(LengthWeight_ParamsAndTypes$Type)] = LengthWeight_ParamsAndTypes$W_estimated[LengthWeight_ParamsAndTypes$Type=="FL" & !is.na(LengthWeight_ParamsAndTypes$Type)]/
					((mean(Sub$cm,na.rm=TRUE) - Length_Length_TL2FL$a)/Length_Length_TL2FL$b)
				LengthWeight_ParamsAndTypes$Type[LengthWeight_ParamsAndTypes$Type=="FL" & !is.na(LengthWeight_ParamsAndTypes$Type)]="TL"
			}
		}
		LengthWeight_ParamsAndTypes$a[!is.na(LengthWeight_ParamsAndTypes$a_new)] = LengthWeight_ParamsAndTypes$a_new[!is.na(LengthWeight_ParamsAndTypes$a_new)]
		LengthWeight_ParamsAndTypes$a_new=NULL
		Parameters$a = LengthWeight_ParamsAndTypes$a
		Parameters$b = LengthWeight_ParamsAndTypes$b
		Parameters$Type_LW = LengthWeight_ParamsAndTypes$Type
		#Now, Parameters data.frame is all in the same measurement type: Total Length! Pass to data.frame
		if(all(is.na(Parameters$Lmax))==FALSE)
		{
			LifeHistoryParametersDataFrame$Lmax[i] = max(Parameters$Lmax,na.rm=TRUE)
		} else {
			LifeHistoryParametersDataFrame$Lmax[i] = NA
		}
		LifeHistoryParametersDataFrame$Lmax_imputed[i]=FALSE
		LifeHistoryParametersDataFrame$Linf_assumedLmax[i]=FALSE
		if(all(is.na(Parameters$Loo))==FALSE)
		{
			LifeHistoryParametersDataFrame$Linf[i] = max(Parameters$Loo,na.rm=TRUE)
		} else {
			LifeHistoryParametersDataFrame$Linf[i] = NA
		}
		LifeHistoryParametersDataFrame$Linf_imputed[i]=FALSE
		if(is.na(LifeHistoryParametersDataFrame$Linf[i]) & !is.na(LifeHistoryParametersDataFrame$Lmax[i])) 
		{
			LifeHistoryParametersDataFrame$Linf[i] = LifeHistoryParametersDataFrame$Lmax[i]  #If Linf missing and Lmax present, set Linf=Lmax
		}
		if(all(is.na(Parameters$Lmat))==FALSE)
		{
			LifeHistoryParametersDataFrame$Lmat[i] = min(Parameters$Lmat,na.rm=TRUE)
		} else {
			LifeHistoryParametersDataFrame$Lmat[i] = NA
		}
		LifeHistoryParametersDataFrame$Lmat_imputed[i]=FALSE
		LifeHistoryParametersDataFrame$Lopt[i] = NA  #note: this value comes from the data directly and can only be calculated from the empirical data.class
		LifeHistoryParametersDataFrame$Lopt_imputed[i]=FALSE
		LifeHistoryParametersDataFrame$Temperature_min[i] = NA #fishbase only provides mean temperature
		LifeHistoryParametersDataFrame$Temperature_max[i] = NA #fishbase only provides mean temperature
		if(all(is.na(Parameters$Temperature))==FALSE)
		{	
			LifeHistoryParametersDataFrame$Temperature_mean[i] = mean(Parameters$Temperature,na.rm=TRUE)
		} else {
			LifeHistoryParametersDataFrame$Temperature_mean[i] = NA
		}
		LifeHistoryParametersDataFrame$WaterTempAssumed[i]=FALSE
		if(all(is.na(Parameters$M))==FALSE)
		{
			LifeHistoryParametersDataFrame$M[i] = min(Parameters$M,na.rm=TRUE)
		} else {
			LifeHistoryParametersDataFrame$M[i] = NA
		}
		LifeHistoryParametersDataFrame$M_imputed[i]=FALSE
		if(all(is.na(Parameters$K))==FALSE)
		{
			LifeHistoryParametersDataFrame$VBK[i] = min(Parameters$K,na.rm=TRUE)
		} else {
			LifeHistoryParametersDataFrame$VBK[i] = NA
		}
		if(all(is.na(Parameters$FoodTroph))==FALSE)
		{
			LifeHistoryParametersDataFrame$TrophicLevel[i] = mean(Parameters$FoodTroph,na.rm=TRUE)
		} else {
			LifeHistoryParametersDataFrame$TrophicLevel[i]=NA
		}
		if(all(is.na(Parameters$FoodSeTroph))==FALSE)
		{
			LifeHistoryParametersDataFrame$TrophicLevelSE[i] = mean(Parameters$FoodSeTroph,na.rm=TRUE)
		} else {
			LifeHistoryParametersDataFrame$TrophicLevelSE[i]=NA
		}	
		LifeHistoryParametersDataFrame$VBK_imputed[i]=FALSE
		if(all(is.na(Parameters$to))==FALSE)
		{
			LifeHistoryParametersDataFrame$t0[i] = mean(Parameters$to,na.rm=TRUE)
		} else {
			LifeHistoryParametersDataFrame$t0[i] = NA
		}
		LifeHistoryParametersDataFrame$t0_assumed[i]=FALSE
		if(all(is.na(Parameters$tmax))==FALSE)
		{
			LifeHistoryParametersDataFrame$Age_max[i] = max(Parameters$tmax,na.rm=TRUE)
		} else {
			LifeHistoryParametersDataFrame$Age_max[i] = NA
		}
		LifeHistoryParametersDataFrame$Age_max_imputed[i]=FALSE
		if(all(is.na(Parameters$a))==FALSE)
		{
			LifeHistoryParametersDataFrame$var_a[i] = mean(Parameters$a,na.rm=TRUE)
		} else {
			LifeHistoryParametersDataFrame$var_a[i] = NA
		}
		LifeHistoryParametersDataFrame$var_a_soureThisDataset[i]=TRUE
		if(all(is.na(Parameters$b))==FALSE)
		{
			LifeHistoryParametersDataFrame$var_b[i] = mean(Parameters$b,na.rm=TRUE)
		} else {
			LifeHistoryParametersDataFrame$var_b[i] = NA
		}
		LifeHistoryParametersDataFrame$var_b_soureThisDataset[i]=TRUE
		LifeHistoryParametersDataFrame$Winf[i] = NA  #This is sometimes provided by FishBase...,
		#...however here we are not going to use this parameter value from FishBase because we...
		#are taking an average of the length-weight parameter and this value will not properly...
		#correspond. If we need it, we will calculate from the a and b values we use. 
		LifeHistoryParametersDataFrame$Winf_imputed[i]=FALSE
		if(all(is.na(Parameters$DepthMax))==FALSE)
		{
			LifeHistoryParametersDataFrame$DepthLimit_deep[i] = max(Parameters$DepthMax,na.rm=TRUE)
		} else {
			LifeHistoryParametersDataFrame$DepthLimit_deep[i] = NA
		}
		if(all(is.na(Parameters$DepthMin))==FALSE)
		{
			LifeHistoryParametersDataFrame$DepthLimit_shallow[i] = min(Parameters$DepthMin,na.rm=TRUE)
		} else {
			LifeHistoryParametersDataFrame$DepthLimit_shallow[i] = NA
		}
		if(all(is.na(Parameters$DepthMax))==FALSE & all(is.na(Parameters$DepthMin))==FALSE)
		{
			LifeHistoryParametersDataFrame$DepthLimit_mean[i] = mean(c(Parameters$DepthMax,Parameters$DepthMin),na.rm=TRUE)
		} else {
			LifeHistoryParametersDataFrame$DepthLimit_mean[i] = NA
		}
		if(all(is.na(Parameters$Vulnerability))==FALSE)
		{
			LifeHistoryParametersDataFrame$Vulnerability[i] = mean(Parameters$Vulnerability,na.rm=TRUE)
		} else {
			LifeHistoryParametersDataFrame$Vulnerability[i] = NA
		}
		if(i%%10==0)
		{
			print(paste("FishBase: Working on species ",i," of ",dim(LifeHistoryParametersDataFrame)[1],"...",sep=""))
			flush.console()
		}
	}
	#Address missing length-weight conversion a and b parameter values: pull all the a and b parameters from FishBase by family and genus, and take the mean - use the mean as estimates for those species where these values are missing
	MissingLengthWeightParams = subset(LifeHistoryParametersDataFrame,is.na(LifeHistoryParametersDataFrame$var_a))
	if(dim(MissingLengthWeightParams)[1]>0)
	{
		MissingLengthWeightParams = subset(MissingLengthWeightParams,select=c(species))
		GenusSpeciesFamily = subset(Lengths,select=c(species,fish_genus,fish_family))
		GenusSpeciesFamily = unique(GenusSpeciesFamily)
		names(GenusSpeciesFamily)[names(GenusSpeciesFamily)=="fish_genus"]="Genus"
		MissingLengthWeightParams = merge(MissingLengthWeightParams,GenusSpeciesFamily,all.x=TRUE,by=c("species"))
		GenusFamily = subset(MissingLengthWeightParams,select=c(fish_family,Genus))
		GenusFamily = unique(GenusFamily)
		GenusFamily$a_lenWeight_estimate = NA
		GenusFamily$b_lenWeight_estimate = NA
		for(i in 1:dim(GenusFamily)[1])
		{
			GenusChar = GenusFamily$Genus[i]
			familyChar = GenusFamily$fish_family[i]
			LengthWeight_GenusFamily_i = data.frame(length_weight(species_list(Family=familyChar,Genus=GenusChar)))
			GenusFamily$a_lenWeight_estimate[i] = mean(LengthWeight_GenusFamily_i$aTL,na.rm=TRUE)
			GenusFamily$b_lenWeight_estimate[i] = mean(LengthWeight_GenusFamily_i$b,na.rm=TRUE)	
			if((i%%10)==0)
			{
				print(paste("Working on item ",i," of ",dim(GenusFamily)[1],"...",sep=""))
				flush.console()
			}
		}
		is.nan.data.frame = function(x)
		{
			do.call(cbind, lapply(x, is.nan))
		}
		GenusFamily[is.nan.data.frame(GenusFamily)]=NA
		for(i in 1:dim(GenusFamily)[1])
		{
			LifeHistoryParametersDataFrame$var_a[is.na(LifeHistoryParametersDataFrame$var_a) & LifeHistoryParametersDataFrame$Genus==GenusFamily$Genus[i]] = GenusFamily$a_lenWeight_estimate[i]
			LifeHistoryParametersDataFrame$var_b[is.na(LifeHistoryParametersDataFrame$var_b) & LifeHistoryParametersDataFrame$Genus==GenusFamily$Genus[i]] = GenusFamily$b_lenWeight_estimate[i]
		} 
	}
	LifeHistoryParametersDataFrame$var_a[is.na(LifeHistoryParametersDataFrame$var_a)]=mean(LifeHistoryParametersDataFrame$var_a,na.rm=TRUE)
	LifeHistoryParametersDataFrame$var_b[is.na(LifeHistoryParametersDataFrame$var_b)]=mean(LifeHistoryParametersDataFrame$var_b,na.rm=TRUE)
	#Convert weight-length a parameters to kg instead of grams
	LifeHistoryParametersDataFrame$var_a = LifeHistoryParametersDataFrame$var_a/1000
	LifeHistoryScenariosList$FromFishBase = LifeHistoryParametersDataFrame
	detach("package:rfishbase",unload=TRUE)
	gc()
	#Scenario 2: Pull Values from FishLife - start with this one because it has relatively accurate temperatures
	library(FishLife)
	LifeHistoryParametersDataFrame = LifeHistoryParametersDataFrame_Generic
	for(i in 1:dim(LifeHistoryParametersDataFrame)[1])
	{
		if(LifeHistoryParametersDataFrame$species[i]=="Monotaxis heterodon")
		{
			next
		}
		GenSppSplit = unlist(strsplit(LifeHistoryParametersDataFrame$species[i]," "))
		FishLifeResults = data.frame(Plot_taxa(Search_species(Genus=GenSppSplit[1],Species=GenSppSplit[2])$match_taxonomy, Database = FishLife::FishBase_and_RAM)[[1]]$Mean_pred)
		dev.off()
		names(FishLifeResults) = "ValueFromFishLife_log"
		FishLifeResults$Parameter=row.names(FishLifeResults)
		row.names(FishLifeResults)=NULL
		FishLifeResults$ValueFromFishLife=NA
		FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="Loo"]=exp(FishLifeResults$ValueFromFishLife_log[FishLifeResults$Parameter=="Loo"])  #cm
		FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="K"]=exp(FishLifeResults$ValueFromFishLife_log[FishLifeResults$Parameter=="K"])  #per year
		FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="Winfinity"]=exp(FishLifeResults$ValueFromFishLife_log[FishLifeResults$Parameter=="Winfinity"])  #grams
		FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="tmax"]=exp(FishLifeResults$ValueFromFishLife_log[FishLifeResults$Parameter=="tmax"])  #longevity in years
		FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="tm"]=exp(FishLifeResults$ValueFromFishLife_log[FishLifeResults$Parameter=="tm"])  #years to maturity
		FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="M"]=exp(FishLifeResults$ValueFromFishLife_log[FishLifeResults$Parameter=="M"])  #natural mortality
		FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="Lm"]=exp(FishLifeResults$ValueFromFishLife_log[FishLifeResults$Parameter=="Lm"])  #length at maturity in cm
		FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="Temperature"]=FishLifeResults$ValueFromFishLife_log[FishLifeResults$Parameter=="Temperature"]  #Temperature NOT in ln space, in celceus 
		FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="ln_var"]=NA  #conditional variance of recruitment devs - not needed for our study
		FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="rho"]=NA   #R1 parameter for recruitment devs  - not needed for our study
		FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="ln_MASPS"]=NA  #marginal annual spawners per recruit
		FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="ln_margsd"]=NA  #maringal standard deviation of recruitment deviations
		FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="h"]=FishLifeResults$ValueFromFishLife_log[FishLifeResults$Parameter=="h"]  #h steepness - not in log space
		FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="logitbound_h"]=NA  #internal value - don't need to use
		FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="ln_Fmsy_over_M"]=NA  #Fmsy divided by M - don't need
		FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="ln_Fmsy"]=exp(FishLifeResults$ValueFromFishLife_log[FishLifeResults$Parameter=="ln_Fmsy"])  #Fmsy - I think real MSY from biomass production type model - can be compared to what Catch-MSY estimtaes 
		FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="ln_r"]=NA  #we have the value not in log space - drop this one
		FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="r"]=FishLifeResults$ValueFromFishLife_log[FishLifeResults$Parameter=="r"]  #Intrinsic growth rate - not in log space
		FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="ln_G"]=NA  #log of generation time - drop becasue we have this not in log space
		FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="G"]=FishLifeResults$ValueFromFishLife_log[FishLifeResults$Parameter=="G"]  #generation time - not in log space - I don't know what the units are (days???)
		FishLifeResults = subset(FishLifeResults,!is.na(FishLifeResults$ValueFromFishLife))
		#Now pass parameters to data.frame
		LifeHistoryParametersDataFrame$Lmax[i] = NA   #No value from FishLife
		LifeHistoryParametersDataFrame$Lmax_imputed[i]=FALSE
		LifeHistoryParametersDataFrame$Linf_assumedLmax[i]=FALSE
		if(!is.na(FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="Loo"]))
		{
			LifeHistoryParametersDataFrame$Linf[i] = FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="Loo"]
		} else {
			LifeHistoryParametersDataFrame$Linf[i] = NA
		}
		LifeHistoryParametersDataFrame$Linf_imputed[i]=FALSE
		if(!is.na(FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="Lm"]))
		{
			LifeHistoryParametersDataFrame$Lmat[i] = FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="Lm"]
		} else {
			LifeHistoryParametersDataFrame$Lmat[i] = NA
		}
		LifeHistoryParametersDataFrame$Lmat_imputed[i]=FALSE	
		LifeHistoryParametersDataFrame$Lopt[i] = NA  #note: this value comes from the data directly and can only be calculated from the empirical data.class
		LifeHistoryParametersDataFrame$Lopt_imputed[i]=FALSE
		LifeHistoryParametersDataFrame$Temperature_min[i] = NA #fishbase only provides mean temperature
		LifeHistoryParametersDataFrame$Temperature_max[i] = NA #fishbase only provides mean temperature
		if(!is.na(FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="Temperature"]))
		{
			LifeHistoryParametersDataFrame$Temperature_mean[i] = FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="Temperature"]
		} else {
			LifeHistoryParametersDataFrame$Temperature_mean[i] = NA
		}
		LifeHistoryParametersDataFrame$WaterTempAssumed[i]=FALSE
		LifeHistoryParametersDataFrame$WaterTempFromFishLife[i]=TRUE
		if(!is.na(FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="M"]))
		{
			LifeHistoryParametersDataFrame$M[i] = FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="M"]
		} else {
			LifeHistoryParametersDataFrame$M[i] = NA
		}
		LifeHistoryParametersDataFrame$M_imputed[i]=FALSE
		if(!is.na(FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="K"]))
		{
			LifeHistoryParametersDataFrame$VBK[i] = FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="K"]
		} else {
			LifeHistoryParametersDataFrame$VBK[i] = NA
		}
		LifeHistoryParametersDataFrame$VBK_imputed[i]=FALSE
		LifeHistoryParametersDataFrame$t0[i] = 0  #Not provided by FishLife - assume it is zero
		LifeHistoryParametersDataFrame$t0_assumed[i]=TRUE
		if(!is.na(FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="tmax"]))
		{
			LifeHistoryParametersDataFrame$Age_max[i] = FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="tmax"]
		} else {
			LifeHistoryParametersDataFrame$Age_max[i] = NA
		}
		LifeHistoryParametersDataFrame$Age_max_imputed[i]=FALSE
		LifeHistoryParametersDataFrame$var_a[i] = LifeHistoryScenariosList$FromFishBase$var_a[i]  #From FishBase, Not provided in FishLife
		LifeHistoryParametersDataFrame$var_a_soureThisDataset[i]=FALSE
		LifeHistoryParametersDataFrame$var_b[i] = LifeHistoryScenariosList$FromFishBase$var_b[i]  #From FishBase, Not provided in FishLife
		LifeHistoryParametersDataFrame$var_b_soureThisDataset[i]=FALSE	
		#since the a and b weight-length parameters are not provided by FishLife, pull and use the Winfinity value
		if(!is.na(FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="Winfinity"]))
		{
			LifeHistoryParametersDataFrame$Winf[i] = FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="Winfinity"]
		} else {
			LifeHistoryParametersDataFrame$Winf[i] = NA
		}
		LifeHistoryParametersDataFrame$Winf_imputed[i]=FALSE
		if(!is.na(FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="h"]))
		{
			LifeHistoryParametersDataFrame$h_steepness[i] = FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="h"]
		} else {
			LifeHistoryParametersDataFrame$h_steepness[i] = NA
		}	
		LifeHistoryParametersDataFrame$h_steepness_imputed[i]=FALSE	
		if(!is.na(FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="G"]))
		{
			LifeHistoryParametersDataFrame$generation_time_years[i] = FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="G"]
		} else {
			LifeHistoryParametersDataFrame$generation_time_years[i] = NA
		}	
		LifeHistoryParametersDataFrame$generation_time_years_imputed[i]=FALSE
		if(!is.na(FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="r"]))
		{
			LifeHistoryParametersDataFrame$r_populationGrowthRate[i] = FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="r"]
		} else {
			LifeHistoryParametersDataFrame$r_populationGrowthRate[i] = NA
		}	
		LifeHistoryParametersDataFrame$r_populationGrowthRate_imputed[i]=FALSE
		if(!is.na(FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="ln_Fmsy"]))
		{
			LifeHistoryParametersDataFrame$Fmsy_FishLife[i] = FishLifeResults$ValueFromFishLife[FishLifeResults$Parameter=="ln_Fmsy"]
		} else {
			LifeHistoryParametersDataFrame$Fmsy_FishLife[i] = NA
		}	
		if(i%%20==0)
		{
			print(paste("FishLife: Working on species ",i," of ",dim(LifeHistoryParametersDataFrame)[1],"...",sep=""))
			flush.console()
		}
	}
	detach("package:FishLife")
	gc()
	LifeHistoryScenariosList$FromFishLife = LifeHistoryParametersDataFrame
	#Life History Scenario 3: Calculate Values from the Raw Data
	LifeHistoryParametersDataFrame = LifeHistoryParametersDataFrame_Generic
	for(i in 1:dim(LifeHistoryParametersDataFrame)[1])
	{
		Sub = subset(Lengths,Lengths$species==LifeHistoryParametersDataFrame$species[i])
		LifeHistoryParametersDataFrame$Lmax[i] = as.numeric(max(Sub$cm,na.rm=TRUE))
		LifeHistoryParametersDataFrame$Lmax_imputed[i]=TRUE
		LifeHistoryParametersDataFrame$Linf[i] = as.numeric(max(Sub$cm,na.rm=TRUE))  #Linf assumed to correspond to the largest fish - same value as Lmax
		LifeHistoryParametersDataFrame$Linf_imputed[i]=TRUE
		LifeHistoryParametersDataFrame$Linf_assumedLmax[i]=TRUE
		if(!is.na(LifeHistoryParametersDataFrame$Linf[i]))
		{
			if(Species$FamilyCommonName[Species$species==LifeHistoryParametersDataFrame$species[i]]=="Snapper")
			{
				LifeHistoryParametersDataFrame$Lmat[i] = LifeHistoryParametersDataFrame$Linf[i] * 0.59 
			}
			if(Species$FamilyCommonName[Species$species==LifeHistoryParametersDataFrame$species[i]]=="Grouper")
			{
				LifeHistoryParametersDataFrame$Lmat[i] = LifeHistoryParametersDataFrame$Linf[i] * 0.46 
			}
			if(Species$FamilyCommonName[Species$species==LifeHistoryParametersDataFrame$species[i]]!="Snapper" & Species$FamilyCommonName[Species$species==LifeHistoryParametersDataFrame$species[i]]!="Grouper")
			{
				LifeHistoryParametersDataFrame$Lmat[i] = LifeHistoryParametersDataFrame$Linf[i] * 0.50 
			}	
			LifeHistoryParametersDataFrame$Lmat_imputed[i]=TRUE
		}
		if(!is.na(LifeHistoryParametersDataFrame$Lmat[i]))
		{
			LifeHistoryParametersDataFrame$Lopt[i] = 1.33 * LifeHistoryParametersDataFrame$Lmat[i]
			LifeHistoryParametersDataFrame$Lopt_imputed[i]=TRUE
		}
		LifeHistoryParametersDataFrame$Temperature_min[i]=NA
		LifeHistoryParametersDataFrame$Temperature_max[i]=NA
		LifeHistoryParametersDataFrame$Temperature_mean[i] = LifeHistoryScenariosList$FromFishLife$Temperature_mean[i]   #Borrow water temperature from FishLife
		LifeHistoryParametersDataFrame$WaterTempFromFishLife[i]=TRUE
		LifeHistoryParametersDataFrame$WaterTempAssumed[i]=FALSE	
		if(is.na(LifeHistoryParametersDataFrame$Temperature_mean[i]))   #If temperature value missing from FishLife, use the average of all values from FishLife of the species that we are analyzing
		{
			LifeHistoryParametersDataFrame$Temperature_mean[i] = mean(LifeHistoryScenariosList$FromFishLife$Temperature_mean,na.rm=TRUE)
		}	
		LifeHistoryParametersDataFrame$Age_max[i] = LifeHistoryScenariosList$FromFishLife$Age_max[i]
		LifeHistoryParametersDataFrame$Age_max_imputed[i]=FALSE
		LifeHistoryParametersDataFrame$M[i] = -1*(log(0.015)/LifeHistoryParametersDataFrame$Age_max[i])  #Calculate M from the Dureuil and Froesie paper in Communications Biology
		LifeHistoryParametersDataFrame$M_imputed[i]=TRUE			
		#Calculate VBK from Pauly formula
		LifeHistoryParametersDataFrame$VBK[i] = (LifeHistoryParametersDataFrame$M[i]*LifeHistoryParametersDataFrame$Lopt[i])/(3*(LifeHistoryParametersDataFrame$Linf[i]-LifeHistoryParametersDataFrame$Lopt[i]))
		LifeHistoryParametersDataFrame$VBK_imputed[i]=TRUE
		LifeHistoryParametersDataFrame$t0[i] = 0  #Assume it is zero - no value provided
		LifeHistoryParametersDataFrame$var_a[i]=LifeHistoryScenariosList$FromFishBase$var_a[i] #borrow from FishLife - can not calculate these values from the data
		LifeHistoryParametersDataFrame$var_a_soureThisDataset[i]=FALSE
		LifeHistoryParametersDataFrame$var_b[i]=LifeHistoryScenariosList$FromFishBase$var_b[i] #borrow from FishLife - can not calculate these valeus from the data
		LifeHistoryParametersDataFrame$var_b_soureThisDataset[i]=FALSE
		LifeHistoryParametersDataFrame$Winf[i]=NA #can not calculate this without the var_a and var_b 
		LifeHistoryParametersDataFrame$Winf_imputed[i]=FALSE
	}
	LifeHistoryScenariosList$CalculatedFromData = LifeHistoryParametersDataFrame	

	#Life History Scenario 4: 5% Above Calculated Values from the Raw Data
	multFactor = 1.05
	LifeHistoryParametersDataFrame = LifeHistoryScenariosList$CalculatedFromData
	LifeHistoryParametersDataFrame$Lmax = LifeHistoryParametersDataFrame$Lmax * multFactor
	LifeHistoryParametersDataFrame$Linf = LifeHistoryParametersDataFrame$Linf * multFactor
	LifeHistoryParametersDataFrame$Lmat = LifeHistoryParametersDataFrame$Lmat * multFactor
	LifeHistoryParametersDataFrame$Lopt = LifeHistoryParametersDataFrame$Lopt * multFactor
	LifeHistoryParametersDataFrame$M = LifeHistoryParametersDataFrame$M * multFactor
	LifeHistoryParametersDataFrame$VBK = LifeHistoryParametersDataFrame$VBK * multFactor
	LifeHistoryParametersDataFrame$Age_max = LifeHistoryParametersDataFrame$Age_max * multFactor
	LifeHistoryScenariosList$CalculatedFromDataUp5 = LifeHistoryParametersDataFrame	
	#Life History Scenario 5: 5% Below Calculated Values from the Raw Data	
	multFactor = 0.95
	LifeHistoryParametersDataFrame = LifeHistoryScenariosList$CalculatedFromData
	LifeHistoryParametersDataFrame$Lmax = LifeHistoryParametersDataFrame$Lmax * multFactor
	LifeHistoryParametersDataFrame$Linf = LifeHistoryParametersDataFrame$Linf * multFactor
	LifeHistoryParametersDataFrame$Lmat = LifeHistoryParametersDataFrame$Lmat * multFactor
	LifeHistoryParametersDataFrame$Lopt = LifeHistoryParametersDataFrame$Lopt * multFactor
	LifeHistoryParametersDataFrame$M = LifeHistoryParametersDataFrame$M * multFactor
	LifeHistoryParametersDataFrame$VBK = LifeHistoryParametersDataFrame$VBK * multFactor
	LifeHistoryParametersDataFrame$Age_max = LifeHistoryParametersDataFrame$Age_max * multFactor
	LifeHistoryScenariosList$CalculatedFromDataDown5 = LifeHistoryParametersDataFrame	
	#FILL HOLES IN LIFE HISTORY INFORMATION
	#Fill holes in FishBase values using Calculated from Data Parameters
	#Lmax
	LifeHistoryScenariosList$FromFishBase$Lmax_imputed[is.na(LifeHistoryScenariosList$FromFishBase$Lmax)]=TRUE
	LifeHistoryScenariosList$FromFishBase$Lmax[is.na(LifeHistoryScenariosList$FromFishBase$Lmax)] = LifeHistoryScenariosList$CalculatedFromData$Lmax[is.na(LifeHistoryScenariosList$FromFishBase$Lmax)]
	LifeHistoryScenariosList$FromFishBase$Lmax_imputed[is.na(LifeHistoryScenariosList$FromFishBase$Lmax)]=FALSE
	#Linf
	LifeHistoryScenariosList$FromFishBase$Linf_imputed[is.na(LifeHistoryScenariosList$FromFishBase$Lmax)]=TRUE
	LifeHistoryScenariosList$FromFishBase$Linf[is.na(LifeHistoryScenariosList$FromFishBase$Linf)] = LifeHistoryScenariosList$CalculatedFromData$Linf[is.na(LifeHistoryScenariosList$FromFishBase$Linf)]
	LifeHistoryScenariosList$FromFishBase$Linf_imputed[is.na(LifeHistoryScenariosList$FromFishBase$Lmax)]=FALSE
	#Lmat
	LifeHistoryScenariosList$FromFishBase$Lmat_imputed[is.na(LifeHistoryScenariosList$FromFishBase$Lmat)]=TRUE
	LifeHistoryScenariosList$FromFishBase$Lmat[is.na(LifeHistoryScenariosList$FromFishBase$Lmat)] = LifeHistoryScenariosList$CalculatedFromData$Lmat[is.na(LifeHistoryScenariosList$FromFishBase$Lmat)]
	LifeHistoryScenariosList$FromFishBase$Lmat_imputed[is.na(LifeHistoryScenariosList$FromFishBase$Lmat)]=FALSE
	#Lopt
	LifeHistoryScenariosList$FromFishBase$Lopt_imputed[is.na(LifeHistoryScenariosList$FromFishBase$Lopt)] = TRUE
	LifeHistoryScenariosList$FromFishBase$Lopt[is.na(LifeHistoryScenariosList$FromFishBase$Lopt)] = LifeHistoryScenariosList$CalculatedFromData$Lopt[is.na(LifeHistoryScenariosList$FromFishBase$Lopt)]
	LifeHistoryScenariosList$FromFishBase$Lopt_imputed[is.na(LifeHistoryScenariosList$FromFishBase$Lopt)] = FALSE
	#Temperature
	LifeHistoryScenariosList$FromFishBase$WaterTempFromFishLife[is.na(LifeHistoryScenariosList$FromFishBase$Temperature_mean)]=TRUE
	LifeHistoryScenariosList$FromFishBase$Temperature_mean[is.na(LifeHistoryScenariosList$FromFishBase$Temperature_mean)] = LifeHistoryScenariosList$FromFishLife$Temperature_mean[is.na(LifeHistoryScenariosList$FromFishBase$Temperature_mean)]
	LifeHistoryScenariosList$FromFishBase$WaterTempFromFishLife[is.na(LifeHistoryScenariosList$FromFishBase$Temperature_mean)]=FALSE
	#M
	LifeHistoryScenariosList$FromFishBase$M_imputed[is.na(LifeHistoryScenariosList$FromFishBase$M)]=TRUE
	LifeHistoryScenariosList$FromFishBase$M[is.na(LifeHistoryScenariosList$FromFishBase$M)] = LifeHistoryScenariosList$CalculatedFromData$M[is.na(LifeHistoryScenariosList$FromFishBase$M)]
	LifeHistoryScenariosList$FromFishBase$M_imputed[is.na(LifeHistoryScenariosList$FromFishBase$M)]=FALSE
	#VBK
	LifeHistoryScenariosList$FromFishBase$VBK_imputed[is.na(LifeHistoryScenariosList$FromFishBase$VBK)]=TRUE
	LifeHistoryScenariosList$FromFishBase$VBK[is.na(LifeHistoryScenariosList$FromFishBase$VBK)] = LifeHistoryScenariosList$CalculatedFromData$VBK[is.na(LifeHistoryScenariosList$FromFishBase$VBK)]
	LifeHistoryScenariosList$FromFishBase$VBK_imputed[is.na(LifeHistoryScenariosList$FromFishBase$VBK)]=FALSE
	#t0
	LifeHistoryScenariosList$FromFishBase$t0_assumed[is.na(LifeHistoryScenariosList$FromFishBase$t0)] = TRUE
	LifeHistoryScenariosList$FromFishBase$t0[is.na(LifeHistoryScenariosList$FromFishBase$t0)] = 0 
	LifeHistoryScenariosList$FromFishBase$t0_assumed[is.na(LifeHistoryScenariosList$FromFishBase$t0)] = FALSE
	#Remove trophic levels that are outlyers
	LifeHistoryScenariosList$FromFishBase$TrophicLevel[LifeHistoryScenariosList$FromFishBase$TrophicLevel > 5]=NA
	#Fill holes in FishLife database using calculated from Data Parameters
	#Lmax
	LifeHistoryScenariosList$FromFishLife$Lmax_imputed[is.na(LifeHistoryScenariosList$FromFishLife$Lmax)]=TRUE
	LifeHistoryScenariosList$FromFishLife$Lmax[is.na(LifeHistoryScenariosList$FromFishLife$Lmax)] = LifeHistoryScenariosList$CalculatedFromData$Lmax[is.na(LifeHistoryScenariosList$FromFishLife$Lmax)]
	LifeHistoryScenariosList$FromFishLife$Lmax_imputed[is.na(LifeHistoryScenariosList$FromFishLife$Lmax)]=FALSE
	#Lopt
	LifeHistoryScenariosList$FromFishLife$Lopt_imputed[is.na(LifeHistoryScenariosList$FromFishLife$Lopt)]=TRUE
	LifeHistoryScenariosList$FromFishLife$Lopt[is.na(LifeHistoryScenariosList$FromFishLife$Lopt)] = LifeHistoryScenariosList$CalculatedFromData$Lopt[is.na(LifeHistoryScenariosList$FromFishLife$Lopt)]
	LifeHistoryScenariosList$FromFishLife$Lopt_imputed[is.na(LifeHistoryScenariosList$FromFishLife$Lopt)]=FALSE
	#VBK (there is only one observation missing)
	LifeHistoryScenariosList$FromFishLife$VBK_imputed[is.na(LifeHistoryScenariosList$FromFishLife$VBK)]=TRUE
	LifeHistoryScenariosList$FromFishLife$VBK[is.na(LifeHistoryScenariosList$FromFishLife$VBK)]=LifeHistoryScenariosList$CalculatedFromData$VBK[is.na(LifeHistoryScenariosList$FromFishLife$VBK)]
	LifeHistoryScenariosList$FromFishLife$VBK_imputed[is.na(LifeHistoryScenariosList$FromFishLife$VBK)]=FALSE
	#t0 (there is only one observation missing)
	LifeHistoryScenariosList$FromFishLife$t0_assumed[is.na(LifeHistoryScenariosList$FromFishLife$t0)]=TRUE
	LifeHistoryScenariosList$FromFishLife$t0[is.na(LifeHistoryScenariosList$FromFishLife$t0)]=0
	#weight-length variables a and b take from FishBase
	LifeHistoryScenariosList$FromFishLife$var_a = LifeHistoryScenariosList$FromFishBase$var_a
	LifeHistoryScenariosList$FromFishLife$var_b = LifeHistoryScenariosList$FromFishBase$var_b
	#For FishBase and FishLife Values: check some of the life history parameters with the size of the samples in the data: make sure that Linf is not too small considering the size of the samples
	Linf_CompareTable = data.frame(species=GenusSpeciesForAnalysis_names,FromFishLife=0,FromFishBase=0,CalculatedFromData=0,CalculatedFromDataUp5=0,CalculatedFromDataDown5=0,MaxLengthInData=0,nSamples=0)
	Lmax_CompareTable = data.frame(species=GenusSpeciesForAnalysis_names,FromFishLife=0,FromFishBase=0,CalculatedFromData=0,CalculatedFromDataUp5=0,CalculatedFromDataDown5=0,MaxLengthInData=0,nSamples=0)
	Lmat_CompareTable = data.frame(species=GenusSpeciesForAnalysis_names,FromFishLife=0,FromFishBase=0,CalculatedFromData=0,CalculatedFromDataUp5=0,CalculatedFromDataDown5=0,MaxLengthInData=0,nSamples=0)
	nSamples_aboveLinf_df = data.frame(nSamplesAboveLinf=rep(0,length(GenusSpeciesForAnalysis_names)))
	nSamples_aboveLmax_df = data.frame(nSamplesAboveLmax=rep(0,length(GenusSpeciesForAnalysis_names)))
	nSamples_aboveLmat_df = data.frame(nSamplesAboveLmat=rep(0,length(GenusSpeciesForAnalysis_names)))
	for(i in 1:length(LifeHistoryScenariosList))
	{
		for(j in 1:length(GenusSpeciesForAnalysis_names))
		{	
			Sub = subset(Lengths,Lengths$species==GenusSpeciesForAnalysis_names[j])	
			Linf_CompareTable$MaxLengthInData[Linf_CompareTable$species==GenusSpeciesForAnalysis_names[j]]=max(Sub$cm,na.rm=TRUE)
			Lmax_CompareTable$MaxLengthInData[Lmax_CompareTable$species==GenusSpeciesForAnalysis_names[j]]=max(Sub$cm,na.rm=TRUE)
			Lmat_CompareTable$MaxLengthInData[Lmat_CompareTable$species==GenusSpeciesForAnalysis_names[j]]=max(Sub$cm,na.rm=TRUE)
			Linf_CompareTable$nSamples[Linf_CompareTable$species==GenusSpeciesForAnalysis_names[j]]=dim(Sub)[1]
			Lmax_CompareTable$nSamples[Lmax_CompareTable$species==GenusSpeciesForAnalysis_names[j]]=dim(Sub)[1]
			Lmat_CompareTable$nSamples[Lmat_CompareTable$species==GenusSpeciesForAnalysis_names[j]]=dim(Sub)[1]
			LifeHistoryScenarios_i_j = subset(LifeHistoryScenariosList[[i]],LifeHistoryScenariosList[[i]]$species==GenusSpeciesForAnalysis_names[j])
			nSamples_aboveLinf_df$nSamplesAboveLinf[j] = dim(subset(Sub,Sub$cm > LifeHistoryScenarios_i_j$Linf))[1]
			nSamples_aboveLmax_df$nSamplesAboveLmax[j] = dim(subset(Sub,Sub$cm > LifeHistoryScenarios_i_j$Lmax))[1]
			nSamples_aboveLmat_df$nSamplesAboveLmat[j] = dim(subset(Sub,Sub$cm > LifeHistoryScenarios_i_j$Lmat))[1]
			Linf_CompareTable$nSamplesAbove[Linf_CompareTable$species==GenusSpeciesForAnalysis_names[j]]=nSamples_aboveLinf_df$nSamplesAboveLinf[j]
			Lmax_CompareTable$nSamplesAbove[Lmax_CompareTable$species==GenusSpeciesForAnalysis_names[j]]=nSamples_aboveLmax_df$nSamplesAboveLmax[j]
			Lmat_CompareTable$nSamplesAbove[Lmat_CompareTable$species==GenusSpeciesForAnalysis_names[j]]=nSamples_aboveLmat_df$nSamplesAboveLmat[j]
			Linf_CompareTable[Linf_CompareTable$species==GenusSpeciesForAnalysis_names[j],names(LifeHistoryScenariosList)[i]] = LifeHistoryScenarios_i_j$Linf
			Lmax_CompareTable[Lmax_CompareTable$species==GenusSpeciesForAnalysis_names[j],names(LifeHistoryScenariosList)[i]] = LifeHistoryScenarios_i_j$Lmax
			Lmat_CompareTable[Lmat_CompareTable$species==GenusSpeciesForAnalysis_names[j],names(LifeHistoryScenariosList)[i]] = LifeHistoryScenarios_i_j$Lmat
			if(j%%5==0)
			{
				print(paste("Working on ",i," , ",j,"...",sep=""))
				flush.console()
			}
		}
		names(nSamples_aboveLinf_df) = paste(names(nSamples_aboveLinf_df),names(LifeHistoryScenariosList)[i],sep="_")
		names(nSamples_aboveLmax_df) = paste(names(nSamples_aboveLmax_df),names(LifeHistoryScenariosList)[i],sep="_")
		names(nSamples_aboveLmat_df) = paste(names(nSamples_aboveLmat_df),names(LifeHistoryScenariosList)[i],sep="_")
		Linf_CompareTable = cbind(Linf_CompareTable,nSamples_aboveLinf_df)
		Lmax_CompareTable = cbind(Lmax_CompareTable,nSamples_aboveLmax_df)
		Lmat_CompareTable = cbind(Lmat_CompareTable,nSamples_aboveLmat_df)	
		nSamples_aboveLinf_df = data.frame(nSamplesAboveLinf=rep(0,length(GenusSpeciesForAnalysis_names)))
		nSamples_aboveLmax_df = data.frame(nSamplesAboveLmax=rep(0,length(GenusSpeciesForAnalysis_names)))
		nSamples_aboveLmat_df = data.frame(nSamplesAboveLmat=rep(0,length(GenusSpeciesForAnalysis_names)))		
	}	
	Linf_CompareTable$FractionSamplesAbove_FishLife = Linf_CompareTable$nSamplesAboveLinf_FromFishLife/Linf_CompareTable$nSamples
	Linf_CompareTable$FractionSamplesAbove_FishBase = Linf_CompareTable$nSamplesAboveLinf_FromFishBase/Linf_CompareTable$nSamples
	Linf_CompareTable$Linf_NEW_FishLife = 0
	Linf_CompareTable$Linf_NEW_FishBase = 0
	findLinfAtFivePercData = function(LenVecIn,start_Linf)
	{
		LenVec = data.frame(Len=LenVecIn,flag=0)
		LenVec$flag[LenVec$Len>start_Linf]=1
		Proportion = sum(LenVec$flag)/dim(LenVec)[1]
		while(Proportion>0.05)
		{
			start_Linf = start_Linf+1
			LenVec = data.frame(Len=LenVecIn,flag=0)
			LenVec$flag[LenVec$Len>start_Linf]=1
			Proportion = sum(LenVec$flag)/dim(LenVec)[1]		
		}
		return(start_Linf)
	}	
	for(i in 1:dim(Linf_CompareTable)[1])
	{
		if(Linf_CompareTable$FractionSamplesAbove_FishLife[i]>0.05)
		{
			Sub = subset(Lengths,Lengths$species==Linf_CompareTable$species[i])
			Linf_CompareTable$Linf_NEW_FishLife[i] = findLinfAtFivePercData(Sub$cm,0)
		}
		if(Linf_CompareTable$FractionSamplesAbove_FishLife[i]<=0.05)
		{
			Linf_CompareTable$Linf_NEW_FishLife[i] = Linf_CompareTable$FromFishLife[i]		
		}
		LifeHistoryScenariosList$FromFishLife$Linf[LifeHistoryScenariosList$FromFishLife$species==Linf_CompareTable$species[i]]=Linf_CompareTable$Linf_NEW_FishLife[i]
		if(Linf_CompareTable$FractionSamplesAbove_FishBase[i]>0.05)
		{
			Sub = subset(Lengths,Lengths$species==Linf_CompareTable$species[i])
			Linf_CompareTable$Linf_NEW_FishBase[i] = findLinfAtFivePercData(Sub$cm,0)			
		}
		if(Linf_CompareTable$FractionSamplesAbove_FishBase[i]<=0.05)
		{
			Linf_CompareTable$Linf_NEW_FishBase[i] = Linf_CompareTable$FromFishBase[i]		
		}
		LifeHistoryScenariosList$FromFishBase$Linf[LifeHistoryScenariosList$FromFishBase$species==Linf_CompareTable$species[i]]=Linf_CompareTable$Linf_NEW_FishBase[i]
		if(i%%5==0)
		{
			print(paste("Working on species ",i,"....",sep=""))
			flush.console()
		}	
	}
	#If Lmax is greater than Linf, set Lmax equal to Linf [NOTE: not an issue for CalculatedFromData
	LifeHistoryScenariosList$FromFishBase$Lmax[LifeHistoryScenariosList$FromFishBase$Lmax<LifeHistoryScenariosList$FromFishBase$Linf]=LifeHistoryScenariosList$FromFishBase$Linf[LifeHistoryScenariosList$FromFishBase$Lmax<LifeHistoryScenariosList$FromFishBase$Linf]
	LifeHistoryScenariosList$FromFishLife$Lmax[LifeHistoryScenariosList$FromFishLife$Lmax<LifeHistoryScenariosList$FromFishLife$Linf]=LifeHistoryScenariosList$FromFishLife$Linf[LifeHistoryScenariosList$FromFishLife$Lmax<LifeHistoryScenariosList$FromFishLife$Linf]
	#If Lmat in FishBase is very different from Lmat in FishLife and is far away from the ratio between Lmat and Linf in FishBase, then use the FishLife ratio and apply it to the Linf in FishBase to get a corrected Lmat that is more realistic
	LifeHistoryScenariosList$FromFishBase$Ratio_Lmat_Linf = LifeHistoryScenariosList$FromFishBase$Lmat/LifeHistoryScenariosList$FromFishBase$Linf
	LifeHistoryScenariosList$FromFishLife$Ratio_Lmat_Linf = LifeHistoryScenariosList$FromFishLife$Lmat/LifeHistoryScenariosList$FromFishLife$Linf
	LifeHistoryScenariosList$FromFishBase$Ratio_Lmat_Linf_FromFishLife = LifeHistoryScenariosList$FromFishLife$Ratio_Lmat_Linf
	LifeHistoryScenariosList$FromFishBase$Diff_Ratio_Lmat_Linf_FishBaseMinusFishLife = LifeHistoryScenariosList$FromFishBase$Ratio_Lmat_Linf - LifeHistoryScenariosList$FromFishBase$Ratio_Lmat_Linf_FromFishLife
	LifeHistoryScenariosList$FromFishBase$Lmat[abs(LifeHistoryScenariosList$FromFishBase$Diff_Ratio_Lmat_Linf_FishBaseMinusFishLife)>0.1] = LifeHistoryScenariosList$FromFishBase$Linf[abs(LifeHistoryScenariosList$FromFishBase$Diff_Ratio_Lmat_Linf_FishBaseMinusFishLife)>0.1] * LifeHistoryScenariosList$FromFishBase$Ratio_Lmat_Linf_FromFishLife[abs(LifeHistoryScenariosList$FromFishBase$Diff_Ratio_Lmat_Linf_FishBaseMinusFishLife)>0.1]
	#Save all four scenarios to an RDS file for checkpointing later on 
	saveRDS(LifeHistoryScenariosList,file=paste(PATH_output,"LifeHistoryScenarios.rds",sep="/"))
	for(i in 1:length(LifeHistoryScenariosList))
	{
		write.table(LifeHistoryScenariosList[[i]],paste(PATH_output,"/LifeHistoryTable_",names(LifeHistoryScenariosList[i]),".csv",sep=""),col.names=TRUE,row.names=FALSE,sep=",")
	}
	#Handle Imputing missing Fecundity values from FishBase and write out file
	write.table(FishBaseFecundityValues,paste(PATH_output,"FishBaseFecundityValues_RAW.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)


}


LifeHistoryScenariosList = readRDS(file=paste(PATH_output,"LifeHistoryScenarios.rds",sep="/"))

######################### Compute CPUE and Days Fished Using Peter and Wawan's Method ########################################################
#Filter Data As Per Wawan's SQL Code for Computing CPUE
#Species_cpue = Species[Species$species_id_number>0,]
#Species_cpue = Species_cpue[Species_cpue$fish_species!="",]
#Species_cpue = subset(Species_cpue,select=c(oid,fish_phylum,fish_class,fish_order,fish_family,fish_genus,fish_species,lmat,lopt,linf,lmax,lmatm,var_a,var_b,conversion_factor_tl2fl))
#names(Species_cpue)[names(Species_cpue)=="oid"]="fish_id"
#names(Species_cpue)[names(Species_cpue)=="fish_family"]="Family"
#names(Species_cpue)[names(Species_cpue)=="fish_genus"]="Genus"
#names(Species_cpue)[names(Species_cpue)=="fish_species"]="Species"

#Boats_cpue = Boats[Boats$gt_estimate>0,]
#names(Boats_cpue)[names(Boats_cpue)=="oid"]="boat_id"
#names(Boats_cpue)[names(Boats_cpue)=="fishing_gear1"]="fishing_gear_boatDataSet"
#names(Boats_cpue)[names(Boats_cpue)=="fishing_wpp1"]="fishing_area1"
#Boats_cpue = subset(Boats_cpue,select=c(boat_id,program_type,fishing_gear_boatDataSet,registration_port,year_built,length_of_boat,gt_estimate,category,fishing_area1))
#Boats_cpue = unique(Boats_cpue)
#Lengths = merge(Lengths,Boats_cpue,all.x=TRUE,by=c("boat_id"))

#Add boat characteristics (i.e. boat gt size, category, WPP fishing in, etc.) to data
#VesselCharacteristics = subset(Lengths,select=c(boat_id,fishing_gear_boatDataSet,length_of_boat,gt_estimate,category,fishing_area1))
#VesselCharacteristics = unique(VesselCharacteristics)
#names(VesselCharacteristics)[names(VesselCharacteristics)=="fishing_area1"]="WPP"
#Work with pre-calculated read in CPUE values from the iFish Community Server
AverageCPUE_WPP_AllSpeciesCombined$X=NULL
AverageCPUE_WPP_AllSpeciesCombined$n=NULL
names(AverageCPUE_WPP_AllSpeciesCombined)[names(AverageCPUE_WPP_AllSpeciesCombined)=="wpp"]="WPP"
names(AverageCPUE_WPP_AllSpeciesCombined)[names(AverageCPUE_WPP_AllSpeciesCombined)=="gear"]="fishing_gear"
names(AverageCPUE_WPP_AllSpeciesCombined)[names(AverageCPUE_WPP_AllSpeciesCombined)=="size"]="VesselCategory"
names(AverageCPUE_WPP_AllSpeciesCombined)[names(AverageCPUE_WPP_AllSpeciesCombined)=="cpue"]="kg_gt_day"
AverageCPUE_EEZ_AllSpeciesCombined = aggregate.data.frame(AverageCPUE_WPP_AllSpeciesCombined$kg_gt_day,by=list(AverageCPUE_WPP_AllSpeciesCombined$fishing_gear,AverageCPUE_WPP_AllSpeciesCombined$VesselCategory),FUN=mean,na.rm=TRUE)
names(AverageCPUE_EEZ_AllSpeciesCombined) = c("fishing_gear","VesselCategory","kg_gt_day")
AverageCPUE_EEZ_AllSpeciesCombined$WPP="EEZ"
AverageCPUE_WPP_AllSpeciesCombined$WPP = as.character(AverageCPUE_WPP_AllSpeciesCombined$WPP)
AverageCPUE_EEZ_AllSpeciesCombined=AverageCPUE_EEZ_AllSpeciesCombined[names(AverageCPUE_WPP_AllSpeciesCombined)]
AverageCPUE_WPP_AllSpeciesCombined = rbind(AverageCPUE_WPP_AllSpeciesCombined,AverageCPUE_EEZ_AllSpeciesCombined)
#Add static CPUE values for missing or too few observations as per Peter and Wawan's suggestion.
StaticValues_withCategory = data.frame(fishing_gear=c("Gillnet","Gillnet","Gillnet","Gillnet","Traps","Traps","Traps","Traps","Mixgears"),
    VesselCategory=c("Nano","Small","Medium","Large","Nano","Small","Medium","Large","Large"),
    kg_gt_day_Static=c(11,11,11,1,7,7,7,7,7),stringsAsFactors = FALSE)
StaticValues_withCategory = merge(StaticValues_withCategory,data.frame(WPP=c(as.character(WPP),"EEZ")))
CpueData = merge(AverageCPUE_WPP_AllSpeciesCombined,StaticValues_withCategory,all.x=TRUE,by=c("fishing_gear","VesselCategory","WPP"))
CpueData$kg_gt_day[is.na(CpueData$kg_gt_day)==TRUE & !is.na(CpueData$kg_gt_day_Static)]=CpueData$kg_gt_day_Static[is.na(CpueData$kg_gt_day)==TRUE & !is.na(CpueData$kg_gt_day_Static)]
CpueData$kg_gt_day_Static=NULL
#Fill in missing cpue values for trap boats
AverageCPUEboatSize = aggregate.data.frame(CpueData$kg_gt_day,by=list(CpueData$VesselCategory,CpueData$WPP),FUN=mean,na.rm=TRUE)
names(AverageCPUEboatSize) = c("VesselCategory","WPP","kg_gt_day_toFill")
CpueData = merge(CpueData,AverageCPUEboatSize,all.x=TRUE,by=c("VesselCategory","WPP"))
CpueData$kg_gt_day[is.na(CpueData$kg_gt_day)]=CpueData$kg_gt_day_toFill[is.na(CpueData$kg_gt_day)]
CpueData$kg_gt_day_toFill=NULL
AverageCPUEboatSize_allAreas = aggregate.data.frame(CpueData$kg_gt_day,by=list(CpueData$VesselCategory),FUN=mean,na.rm=TRUE)
names(AverageCPUEboatSize_allAreas) = c("VesselCategory","kg_gt_day_toFill")
CpueData = merge(CpueData,AverageCPUEboatSize_allAreas,all.x=TRUE,by=c("VesselCategory"))
CpueData$kg_gt_day[is.na(CpueData$kg_gt_day)]=CpueData$kg_gt_day_toFill[is.na(CpueData$kg_gt_day)]
CpueData$kg_gt_day_toFill=NULL
write.table(CpueData,paste(PATH_output,"CPUEdata.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
#Calculate Average Days Fishing
BoatTracker = subset(BoatTracker,select=c(boat_id,tracker_id))
BoatTracker = unique(BoatTracker)
LocationPings = merge(LocationPings,BoatTracker,all.x=TRUE,by=c("tracker_id"))
FishingDays = subset(LocationPings,select=c(boat_id,date_time))
dateTimeVec_vec = lapply(X=FishingDays$date_time,FUN=function(x) {unlist(strsplit(as.character(x)," "))})
dateTimeVec_vec = do.call(rbind.data.frame,dateTimeVec_vec)
rownames(dateTimeVec_vec)=NULL
names(dateTimeVec_vec) = c("Date","Time")
dateOnlyVec_vec = lapply(X=dateTimeVec_vec$Date,FUN=function(x) {unlist(strsplit(as.character(x),"-"))})
dateOnlyVec_vec = do.call(rbind.data.frame,dateOnlyVec_vec)
rownames(dateOnlyVec_vec)=NULL
names(dateOnlyVec_vec) = c("Year_OutFishing","Month_OutFishing","Day_OutFishing")
dateOnlyVec_vec$Year_OutFishing = as.numeric(as.character(dateOnlyVec_vec$Year_OutFishing))
dateOnlyVec_vec$Month_OutFishing = as.numeric(as.character(dateOnlyVec_vec$Month_OutFishing))
dateOnlyVec_vec$Day_OutFishing = as.numeric(as.character(dateOnlyVec_vec$Day_OutFishing))
FishingDays = cbind(FishingDays,dateOnlyVec_vec)
FishingDays$date_time=NULL
FishingDays = unique(FishingDays)
FishingDays = data.frame(table(FishingDays$boat_id,FishingDays$Year_OutFishing))
names(FishingDays) = c("boat_id","Year_OutFishing","DaysFishing")
FishingDays$boat_id = as.character(FishingDays$boat_id)
FishingDays$Year_OutFishing = as.character(FishingDays$Year_OutFishing)
FishingDays = subset(FishingDays,FishingDays$DaysFishing>0)
FishingDays = merge(FishingDays,Boats,all.x=TRUE,by=c("boat_id"))
#Not all of the records for each vessel in "Boats" has data attached to it such as gross tons, etc. Therefore, develop scaling factors for the unlabeled vessels, to account for them and properly scale up the catch to the amount per year it should be
FishingDays$CategoryNew = paste(FishingDays$category,FishingDays$VesselCategory)
FishingDays$CategoryNew[is.na(FishingDays$category) | is.na(FishingDays$VesselCategory)]="NO CATEGORY ASSIGNED"
ScalingFactor = as.data.frame.matrix(table(FishingDays$Year_OutFishing,FishingDays$CategoryNew))
ScalingFactor$TOTAL = rowSums(ScalingFactor)
ScalingFactor$Year = rownames(ScalingFactor)
rownames(ScalingFactor)=NULL
ScalingFactor$RatioVesselsNoCharacteristics = ScalingFactor$"NO CATEGORY ASSIGNED"/ScalingFactor$TOTAL
ScalingFactor = subset(ScalingFactor,select=c(Year,RatioVesselsNoCharacteristics))  #Use this a little below to adjust estimates to include all vessels
#Continue Calculating Average Days Fished per Fleet Segment  
FishingDays$YearNew = FishingDays$Year_OutFishing  #Create "YearNew" variable to combine some years becuase very few observations in some years to best be able to estimate average days fishing
FishingDays$YearNew[FishingDays$Year_OutFishing=="2014"]="2014_2015_2016_2017" 
FishingDays$YearNew[FishingDays$Year_OutFishing=="2015"]="2014_2015_2016_2017" 
FishingDays$YearNew[FishingDays$Year_OutFishing=="2016"]="2014_2015_2016_2017"
FishingDays$YearNew[FishingDays$Year_OutFishing=="2017"]="2014_2015_2016_2017"
FishingDays$YearNew[FishingDays$Year_OutFishing=="2021"]="2021_2022_2023"
FishingDays$YearNew[FishingDays$Year_OutFishing=="2022"]="2021_2022_2023"
FishingDays$YearNew[FishingDays$Year_OutFishing=="2023"]="2021_2022_2023"
Year_YearNew_Convert = subset(FishingDays,select=c(Year_OutFishing,YearNew))
Year_YearNew_Convert = unique(Year_YearNew_Convert)
names(Year_YearNew_Convert) = c("Year","YearNew")

Samples_FishingDays = as.data.frame(table(FishingDays$YearNew,FishingDays$category,FishingDays$gear,FishingDays$VesselCategory))
names(Samples_FishingDays) = c("YearNew","category","gear","VesselSize","numSamples")
Samples_FishingDays = subset(Samples_FishingDays,Samples_FishingDays$numSamples>0)
AvgDaysFishedPerYear_withCategory_mean = aggregate.data.frame(FishingDays$DaysFishing,by=list(FishingDays$YearNew,FishingDays$category,FishingDays$gear,FishingDays$VesselCategory),FUN=mean,na.rm=TRUE)
names(AvgDaysFishedPerYear_withCategory_mean) = c("YearNew","category","gear","VesselSize","AvgDaysFished")
AvgDaysFishedPerYear_withCategory_max = aggregate.data.frame(FishingDays$DaysFishing,by=list(FishingDays$YearNew,FishingDays$category,FishingDays$gear,FishingDays$VesselCategory),FUN=max,na.rm=TRUE)
names(AvgDaysFishedPerYear_withCategory_max) = c("YearNew","category","gear","VesselSize","MaxDaysFished")
AvgDaysFishedPerYear_withCategory_sd = aggregate.data.frame(FishingDays$DaysFishing,by=list(FishingDays$YearNew,FishingDays$category,FishingDays$gear,FishingDays$VesselCategory),FUN=sd,na.rm=TRUE)
names(AvgDaysFishedPerYear_withCategory_sd) = c("YearNew","category","gear","VesselSize","StdevDaysFished")
AvgDaysFishedPerYear_withCategory = merge(AvgDaysFishedPerYear_withCategory_mean,AvgDaysFishedPerYear_withCategory_sd,all.x=TRUE,by=c("YearNew","category","gear","VesselSize"))
AvgDaysFishedPerYear_withCategory = merge(AvgDaysFishedPerYear_withCategory,AvgDaysFishedPerYear_withCategory_max,all.x=TRUE,by=c("YearNew","category","gear","VesselSize"))
AvgDaysFishedPerYear_withCategory = merge(AvgDaysFishedPerYear_withCategory,Samples_FishingDays,all.x=TRUE,by=c("YearNew","category","gear","VesselSize"))
AvgDaysFishedPerYear_withCategory = merge(AvgDaysFishedPerYear_withCategory,Year_YearNew_Convert,all=TRUE,by=c("YearNew"))
AvgDaysFishedPerYear_withCategory = AvgDaysFishedPerYear_withCategory[order(AvgDaysFishedPerYear_withCategory$category,AvgDaysFishedPerYear_withCategory$gear,AvgDaysFishedPerYear_withCategory$VesselSize,AvgDaysFishedPerYear_withCategory$Year),]
YearsCategoryGearVesselSize = expand.grid(Year=unique(AvgDaysFishedPerYear_withCategory$Year),category=unique(AvgDaysFishedPerYear_withCategory$category),
	gear=unique(AvgDaysFishedPerYear_withCategory$gear),VesselSize=unique(AvgDaysFishedPerYear_withCategory$VesselSize))
AvgDaysFishedPerYear_withCategory = merge(AvgDaysFishedPerYear_withCategory,YearsCategoryGearVesselSize,all.y=TRUE,by=c("Year","category","gear","VesselSize"))
#***********************************************************************************************************************************************
#************************** Temporary Hard Coded Area With Wawan's Numbers: I don't know how he calculated these numbers ***********************
#***********************************************************************************************************************************************
gear2 <- c("Dropline", "Dropline", "Dropline", "Dropline", "Dropline", "Dropline", "Dropline", "Dropline",
           "Longline", "Longline", "Longline", "Longline", "Longline", "Longline", "Longline", "Longline",
           "Gillnet", "Gillnet", "Gillnet", "Gillnet", "Gillnet", "Gillnet", "Gillnet", "Gillnet",
           "Trap", "Trap", "Trap", "Trap", "Trap", "Trap", "Trap", "Trap",
           "Mixgears", "Mixgears", "Mixgears", "Mixgears", "Mixgears", "Mixgears", "Mixgears", "Mixgears")
size2 <- c("Nano Dedicated", "Nano Seasonal", "Small Dedicated", "Small Seasonal", "Medium Dedicated", "Medium Seasonal", "Large Dedicated", "Large Seasonal",
           "Nano Dedicated", "Nano Seasonal", "Small Dedicated", "Small Seasonal", "Medium Dedicated", "Medium Seasonal", "Large Dedicated", "Large Seasonal",
           "Nano Dedicated", "Nano Seasonal", "Small Dedicated", "Small Seasonal", "Medium Dedicated", "Medium Seasonal", "Large Dedicated", "Large Seasonal",
           "Nano Dedicated", "Nano Seasonal", "Small Dedicated", "Small Seasonal", "Medium Dedicated", "Medium Seasonal", "Large Dedicated", "Large Seasonal",
           "Nano Dedicated", "Nano Seasonal", "Small Dedicated", "Small Seasonal", "Medium Dedicated", "Medium Seasonal", "Large Dedicated", "Large Seasonal")
fdays <- c(188, 94,  196, 98,  178, 89,  178, 89,
           202, 101, 202, 101, 242, 121, 244, 122,
           196, 98,  199, 99,  282, 141, 211, 105,
           197, 98,  199, 99,  214, 107, 211, 105,
           204, 102, 201, 100, 155, 77,  211, 105)
HardCodedAvgDaysFishedPerYear_withCategory = data.frame(gear=gear2,VesselSize=size2,AvgDaysFished_hardCoded=fdays)
dateTimeVec_vec = lapply(X=HardCodedAvgDaysFishedPerYear_withCategory$VesselSize,FUN=function(x) {unlist(strsplit(as.character(x)," "))})
dateTimeVec_vec = do.call(rbind.data.frame,dateTimeVec_vec)
rownames(dateTimeVec_vec)=NULL
names(dateTimeVec_vec) = c("VesselSize","category")
HardCodedAvgDaysFishedPerYear_withCategory$VesselSize=NULL
HardCodedAvgDaysFishedPerYear_withCategory = cbind(HardCodedAvgDaysFishedPerYear_withCategory,dateTimeVec_vec)
AvgDaysFishedPerYear_withCategory = merge(AvgDaysFishedPerYear_withCategory,HardCodedAvgDaysFishedPerYear_withCategory,all.x=TRUE,by=c("gear","VesselSize","category"))
AvgDaysFishedPerYear_withCategory = AvgDaysFishedPerYear_withCategory[order(AvgDaysFishedPerYear_withCategory$category,AvgDaysFishedPerYear_withCategory$gear,AvgDaysFishedPerYear_withCategory$VesselSize,AvgDaysFishedPerYear_withCategory$Year),]
AvgDaysFishedPerYear_withCategory$AvgDaysFished[is.na(AvgDaysFishedPerYear_withCategory$AvgDaysFished)] = AvgDaysFishedPerYear_withCategory$AvgDaysFished_hardCoded[is.na(AvgDaysFishedPerYear_withCategory$AvgDaysFished)]
AvgDaysFishedPerYear_withCategory = subset(AvgDaysFishedPerYear_withCategory,!is.na(AvgDaysFishedPerYear_withCategory$AvgDaysFished))
AvgDaysFishedPerYear_withCategory$PlusOneStandardDeviation = AvgDaysFishedPerYear_withCategory$AvgDaysFished + AvgDaysFishedPerYear_withCategory$StdevDaysFished
#To clean the data: calculate one standard deviation from the "AvgDaysFished" value, flag any values outside of HALF OF ONE standard deviation, and replace them with the "AvgDaysFished" value plus one standard deviation ("PlusOneStandardDeviation")
AvgDaysFishedPerYear_withCategory$Difference_AvgDaysFished_PlusOneStdev = AvgDaysFishedPerYear_withCategory$AvgDaysFished - AvgDaysFishedPerYear_withCategory$AvgDaysFished_hardCoded
AbsDifferencesWithoutZero = abs(AvgDaysFishedPerYear_withCategory$Difference_AvgDaysFished_PlusOneStdev)
AbsDifferencesWithoutZero = subset(AbsDifferencesWithoutZero,AbsDifferencesWithoutZero>0)
sd_Difference_AvgDaysFished = sd(AbsDifferencesWithoutZero)   #This represents the standard deviation distance from zero both positive and negative that we should look to include or exclude estimated values
AvgDaysFishedPerYear_withCategory$keepCalculatedValue = FALSE
AvgDaysFishedPerYear_withCategory$keepCalculatedValue[(AvgDaysFishedPerYear_withCategory$Difference_AvgDaysFished_PlusOneStdev>(sd_Difference_AvgDaysFished*-1) & AvgDaysFishedPerYear_withCategory$Difference_AvgDaysFished_PlusOneStdev<sd_Difference_AvgDaysFished)]=TRUE
AvgDaysFishedPerYear_withCategory$keepCalculatedValue[is.na(AvgDaysFishedPerYear_withCategory$AvgDaysFished_hardCoded)]=TRUE
AvgDaysFishedPerYear_withCategory$AvgDaysFished[AvgDaysFishedPerYear_withCategory$keepCalculatedValue==FALSE] = AvgDaysFishedPerYear_withCategory$PlusOneStandardDeviation[AvgDaysFishedPerYear_withCategory$keepCalculatedValue==FALSE]  #where the value was one standard deviation away, then use the PlusOneStandardDeviation value and replace with it
AvgDaysFishedPerYear_withCategory$Difference_AvgDaysFished_PlusOneStdev = AvgDaysFishedPerYear_withCategory$AvgDaysFished - AvgDaysFishedPerYear_withCategory$AvgDaysFished_hardCoded  #recalculate the difference now that we have done some replacements 
AvgDaysFishedPerYear_withCategory$keepCalculatedValue = FALSE  #check again to see if any of the values replaced with the PlusOneStandardDeviation are still outside th standard deviation, and if so, replace those with the fixed value. 
AvgDaysFishedPerYear_withCategory$keepCalculatedValue[(AvgDaysFishedPerYear_withCategory$Difference_AvgDaysFished_PlusOneStdev>(sd_Difference_AvgDaysFished*-1) & AvgDaysFishedPerYear_withCategory$Difference_AvgDaysFished_PlusOneStdev<sd_Difference_AvgDaysFished)]=TRUE
AvgDaysFishedPerYear_withCategory$keepCalculatedValue[is.na(AvgDaysFishedPerYear_withCategory$AvgDaysFished_hardCoded)]=TRUE
AvgDaysFishedPerYear_withCategory$AvgDaysFished[AvgDaysFishedPerYear_withCategory$keepCalculatedValue==FALSE] = AvgDaysFishedPerYear_withCategory$AvgDaysFished_hardCoded[AvgDaysFishedPerYear_withCategory$keepCalculatedValue==FALSE]
AvgDaysFishedPerYear_withCategory = subset(AvgDaysFishedPerYear_withCategory,select=c(gear,VesselSize,category,Year,AvgDaysFished))

#***********************************************************************************************************************************************
#***********************************************************************************************************************************************
#CPUE * Days Fished * gross ton across entire fleet to estimate total cathch, then scale it for vessels without gross ton estimates 
names(Boats)[names(Boats)=="WPP_boats"]="WPP"
names(CpueData)[names(CpueData)=="fishing_gear"]="gear"
FullFleet = merge(Boats,CpueData,all.x=TRUE,by=c("VesselCategory","gear","WPP"))
FullFleet = subset(FullFleet,!is.na(FullFleet$gt_estimate))
names(AvgDaysFishedPerYear_withCategory)[names(AvgDaysFishedPerYear_withCategory)=="VesselSize"]="VesselCategory"
FullFleet = merge(FullFleet,AvgDaysFishedPerYear_withCategory,all.x=TRUE,by=c("gear","VesselCategory","category"))
FullFleet$TotalCatch_kg = FullFleet$kg_gt_day*FullFleet$AvgDaysFished*FullFleet$gt_estimate
FullFleet$TotalCatch_mT = FullFleet$TotalCatch_kg*0.001
CatchPerWPP = aggregate.data.frame(FullFleet$TotalCatch_mT,by=list(FullFleet$Year,FullFleet$WPP),FUN=sum,na.rm=TRUE)
names(CatchPerWPP) = c("Year","WPP","TotalCatch_mT")
CatchPerWPP$yearsCombined=FALSE
CatchPerWPP$yearsCombined[grep("_",CatchPerWPP$Year)]=TRUE
CatchPerWPP_noYearsCombined = CatchPerWPP[1,]
CatchPerWPP_noYearsCombined = CatchPerWPP_noYearsCombined[-1,]
for(i in 1:dim(CatchPerWPP)[1])
{
	if(CatchPerWPP$yearsCombined[i]==TRUE)
	{
		yrsSplit = as.numeric(unlist(strsplit(CatchPerWPP$Year[i],"_")))
		for(j in 1:length(yrsSplit))
		{
			temp = CatchPerWPP[i,]
			temp$Year=yrsSplit[j]
			CatchPerWPP_noYearsCombined = rbind(CatchPerWPP_noYearsCombined,temp)
		}	
	}
	if(CatchPerWPP$yearsCombined[i]==FALSE)
	{
		CatchPerWPP_noYearsCombined = rbind(CatchPerWPP_noYearsCombined,CatchPerWPP[i,])
	}
}
#Did scaling in the days out fishing (effort) estimations
#CatchPerWPP = merge(CatchPerWPP_noYearsCombined,ScalingFactor,all.x=TRUE,by=c("Year"))
#CatchPerWPP$TotalCatch_mT = CatchPerWPP$TotalCatch_mT/(1-CatchPerWPP$RatioVesselsNoCharacteristics)
#CatchPerWPP$RatioVesselsNoCharacteristics=NULL

################## Compute Amount of Catch By Fish Family and Species in Each WPP #############################

#Calculate the Proportion of Each Family Within Each WPP
Lengths$Weight_kg = Lengths$var_a*(Lengths$cm^Lengths$var_b)
Lengths_WPPnotMissing = Lengths[!is.na(Lengths$WPP) & Lengths$FamilyCommonName!="",]
FamilySampledPerWPP = aggregate.data.frame(Lengths_WPPnotMissing$Weight_kg,by=list(Lengths_WPPnotMissing$FamilyCommonName,Lengths_WPPnotMissing$WPP),FUN=sum,na.rm=TRUE)
names(FamilySampledPerWPP) = c("FamilyCommonName","WPP","Weight_kg")
FamilySampled = aggregate.data.frame(Lengths_WPPnotMissing$Weight_kg,by=list(Lengths_WPPnotMissing$WPP),FUN=sum,na.rm=TRUE)
names(FamilySampled) = c("WPP","Weight_kg_WPP")
FamilySampledPerWPP = merge(FamilySampledPerWPP,FamilySampled,all.x=TRUE,by=c("WPP"))
FamilySampledPerWPP$ProportionOfFamilyInWPP = FamilySampledPerWPP$Weight_kg/FamilySampledPerWPP$Weight_kg_WPP
FamilySampledPerWPP$Weight_kg_WPP=NULL
FamilySampledPerWPP$Weight_kg=NULL
FamilySampledPerWPP = reshape(FamilySampledPerWPP,v.names="ProportionOfFamilyInWPP",idvar="FamilyCommonName",timevar="WPP",direction="wide")
FamilySampledPerWPP[is.na(FamilySampledPerWPP)]=0
#Apply the Proportion of Each Family Within Each WPP to Sea Around Us Data
SeaAroundUsIndo = read.csv(paste(PATH_otherInput,"SeaAroundUs/SAU EEZ 937,361,938 v50-0/SAU EEZ 937,361,938 v50-0.csv",sep="/"),header=TRUE,sep=",",stringsAsFactors=FALSE)
SeaAroundUs_TNC_FamilyMapping = read.table(paste(PATH_otherInput,"SeaAroundUs/SeaAroundUs_Family_mapTo_FamilyGrouping.csv",sep="/"),header=TRUE,sep=",",stringsAsFactors=FALSE)
names(SeaAroundUsIndo)[names(SeaAroundUsIndo)=="common_name"]="SeaAroundUsFamily"
SeaAroundUsIndo = merge(SeaAroundUsIndo,SeaAroundUs_TNC_FamilyMapping,all.x=TRUE,by=c("SeaAroundUsFamily"))
SeaAroundUsIndo_byFamily = aggregate.data.frame(SeaAroundUsIndo$tonnes,by=list(SeaAroundUsIndo$FamilyCommonName,SeaAroundUsIndo$year),FUN=sum)
names(SeaAroundUsIndo_byFamily) = c("FamilyCommonName","Year","tonnes")
SeaAroundUsIndo_byFamily = merge(SeaAroundUsIndo_byFamily,FamilySampledPerWPP,all.x=TRUE,by=c("FamilyCommonName"))
SeaAroundUsIndo_byFamily = SeaAroundUsIndo_byFamily[!is.na(SeaAroundUsIndo_byFamily$ProportionOfFamilyInWPP.571),]
SeaAroundUsIndo_OnlyFractions = SeaAroundUsIndo_byFamily[,grep("ProportionOfFamilyInWPP",names(SeaAroundUsIndo_byFamily))]
SeaAroundUsIndo_RowInfo = SeaAroundUsIndo_byFamily[,grep("ProportionOfFamilyInWPP",names(SeaAroundUsIndo_byFamily),invert=TRUE)]
tonnes_rep = do.call("cbind", replicate((dim(SeaAroundUsIndo_OnlyFractions)[2]),data.frame(SeaAroundUsIndo_byFamily$tonnes),simplify = FALSE))
SeaAroundUsIndo_OnlyCatchPerWPP = SeaAroundUsIndo_OnlyFractions*tonnes_rep
SeaAroundUsIndo_byFamily = cbind(SeaAroundUsIndo_RowInfo,SeaAroundUsIndo_OnlyCatchPerWPP)
SeaAroundUsIndo_byFamily = SeaAroundUsIndo_byFamily[order(SeaAroundUsIndo_byFamily$FamilyCommonName,SeaAroundUsIndo_byFamily$Year),]
names(SeaAroundUsIndo_byFamily) = gsub("ProportionOfFamilyInWPP","Catch_mT",names(SeaAroundUsIndo_byFamily))
#Calculate the fraction of each species, in each family, in each WPP
FractionSpeciesInFamilyAndWPP = aggregate.data.frame(Lengths_WPPnotMissing$Weight_kg,by=list(Lengths_WPPnotMissing$species,Lengths_WPPnotMissing$FamilyCommonName,Lengths_WPPnotMissing$WPP),FUN=sum)
names(FractionSpeciesInFamilyAndWPP) = c("species","FamilyCommonName","WPP","Weight_kg")
FractionFamilyAndWPP = aggregate.data.frame(Lengths_WPPnotMissing$Weight_kg,by=list(Lengths_WPPnotMissing$FamilyCommonName,Lengths_WPPnotMissing$WPP),FUN=sum)
names(FractionFamilyAndWPP) = c("FamilyCommonName","WPP","Weight_kg_FamilyWPP")
FractionSpeciesInFamilyAndWPP = merge(FractionSpeciesInFamilyAndWPP,FractionFamilyAndWPP,all.x=TRUE,by=c("FamilyCommonName","WPP"))
FractionSpeciesInFamilyAndWPP$ProportionOfSpeciesInFamilyAndWPP = FractionSpeciesInFamilyAndWPP$Weight_kg/FractionSpeciesInFamilyAndWPP$Weight_kg_FamilyWPP
FractionSpeciesInFamilyAndWPP$Weight_kg=NULL
FractionSpeciesInFamilyAndWPP$Weight_kg_FamilyWPP=NULL
FractionSpeciesInFamilyAndWPP = reshape(FractionSpeciesInFamilyAndWPP,v.names="ProportionOfSpeciesInFamilyAndWPP",idvar="species",timevar="WPP",direction="wide")
FractionSpeciesInFamilyAndWPP[is.na(FractionSpeciesInFamilyAndWPP)]=0
#Calculate the historical catch by species and WPP by applying the fractions to the Sea Around Us data.
SeaAroundUsIndo_bySpecies = merge(SeaAroundUsIndo_byFamily,FractionSpeciesInFamilyAndWPP,by=c("FamilyCommonName"))
FamilyCatchMatrix = SeaAroundUsIndo_bySpecies[,grep("Catch_mT",names(SeaAroundUsIndo_bySpecies))]
SpeciesProportionMatrix = SeaAroundUsIndo_bySpecies[,grep("ProportionOfSpeciesInFamilyAndWPP",names(SeaAroundUsIndo_bySpecies))]
names(FamilyCatchMatrix) = gsub("Catch_mT.","",names(FamilyCatchMatrix))
names(SpeciesProportionMatrix) = gsub("ProportionOfSpeciesInFamilyAndWPP.","",names(SpeciesProportionMatrix))
SpeciesProportionMatrix = SpeciesProportionMatrix[,names(FamilyCatchMatrix)]
FamilyCatchMatrix = FamilyCatchMatrix[,as.character(WPP)]
SpeciesProportionMatrix = SpeciesProportionMatrix[,as.character(WPP)]
SpeciesCatchMatrix_mT = FamilyCatchMatrix*SpeciesProportionMatrix
RowInfo = subset(SeaAroundUsIndo_bySpecies,select=c(FamilyCommonName,Year,species))
SpeciesCatchMatrix_mT = cbind(RowInfo,SpeciesCatchMatrix_mT)
#Group those species that will not be inclued in the analysis as "Other"
SpeciesCatchMatrix_mT = merge(SpeciesCatchMatrix_mT,data.frame(species=GenusSpeciesForAnalysis_names,flag=1),all.x=TRUE,by=c("species"))
SpeciesCatchMatrix_mT$species[is.na(SpeciesCatchMatrix_mT$flag)]="Other"
SpeciesCatchMatrix_mT$flag=NULL
SpeciesCatchMatrix_mT$FamilyCommonName[SpeciesCatchMatrix_mT$species=="Other"]="Other"
SpeciesCatchMatrix_mT = aggregate.data.frame(SpeciesCatchMatrix_mT[,4:dim(SpeciesCatchMatrix_mT)[2]],by=list(SpeciesCatchMatrix_mT$species,SpeciesCatchMatrix_mT$FamilyCommonName,SpeciesCatchMatrix_mT$Year),FUN=sum,na.rm=TRUE)
names(SpeciesCatchMatrix_mT)[1:3]=c("species","FamilyCommonName","Year")
SpeciesCatchMatrix_mT = SpeciesCatchMatrix_mT[order(SpeciesCatchMatrix_mT$species,SpeciesCatchMatrix_mT$Year),]
#Calculate the modern catch by family and species since CODRS started: proportion of each family in each WPP applied to "FullFleet" catch as estimated using the CODRS Fleet Survey that Peter and TNC conduct
FamilySampled = aggregate.data.frame(Lengths_WPPnotMissing$Weight_kg,by=list(Lengths_WPPnotMissing$FamilyCommonName),FUN=sum,na.rm=TRUE)
names(FamilySampled) = c("FamilyCommonName","Weight_kg_Family")
FamilySampled$ProportionSampled_Family = FamilySampled$Weight_kg_Family/sum(FamilySampled$Weight_kg_Family)
FullFleet_byYear = aggregate.data.frame(FullFleet$TotalCatch_mT,by=list(FullFleet$Year),FUN=sum,na.rm=TRUE)
names(FullFleet_byYear) = c("Year","TotalCatch_mT")
FamilySampled = merge(FamilySampled,FullFleet_byYear)
FamilySampled$Catch_Family_mT = FamilySampled$ProportionSampled_Family*FamilySampled$TotalCatch_mT
FamilySampled$Weight_kg_Family=NULL
FamilySampled$ProportionSampled_Family=NULL
FamilySampled$TotalCatch_mT=NULL
#Calculate the Proportion of catch within each Family for each WPP
Lengths$Weight_kg = Lengths$var_a*(Lengths$cm^Lengths$var_b)
Lengths_WPPnotMissing = Lengths[!is.na(Lengths$WPP) & Lengths$FamilyCommonName!="",]
FamilySampledWithinWPP = aggregate.data.frame(Lengths_WPPnotMissing$Weight_kg,by=list(Lengths_WPPnotMissing$FamilyCommonName,Lengths_WPPnotMissing$WPP),FUN=sum,na.rm=TRUE)
names(FamilySampledWithinWPP) = c("FamilyCommonName","WPP","Weight_kg")
FamilySampledPerFamily = aggregate.data.frame(Lengths_WPPnotMissing$Weight_kg,by=list(Lengths_WPPnotMissing$FamilyCommonName),FUN=sum,na.rm=TRUE)
names(FamilySampledPerFamily) = c("FamilyCommonName","Weight_kg_Family")
FamilySampledWithinWPP = merge(FamilySampledWithinWPP,FamilySampledPerFamily,all.x=TRUE,by=c("FamilyCommonName"))
FamilySampledWithinWPP$ProportionOfFamilyInWPP = FamilySampledWithinWPP$Weight_kg/FamilySampledWithinWPP$Weight_kg_Family
FamilySampledWithinWPP$Weight_kg_Family=NULL
FamilySampledWithinWPP$Weight_kg=NULL
FamilySampledWithinWPP = reshape(FamilySampledWithinWPP,v.names="ProportionOfFamilyInWPP",idvar="FamilyCommonName",timevar="WPP",direction="wide")
FamilySampledWithinWPP[is.na(FamilySampledWithinWPP)]=0
CatchPerFamilyWithinWPP = merge(FamilySampledWithinWPP,FamilySampled,all=TRUE,by=c("FamilyCommonName"))
CatchPerFamilyWithinWPP[,grep("ProportionOfFamilyInWPP",names(CatchPerFamilyWithinWPP))] = CatchPerFamilyWithinWPP[,grep("ProportionOfFamilyInWPP",names(CatchPerFamilyWithinWPP))] * CatchPerFamilyWithinWPP$Catch_Family_mT
names(CatchPerFamilyWithinWPP) = gsub("ProportionOfFamilyInWPP.","Catch_mT_",names(CatchPerFamilyWithinWPP))
CatchPerFamilyWithinWPP$Catch_Family_mT=NULL
#For the modern time period of catches/landings fractioin out the portion of species within each fish family for each year and WPP
CatchPerYearWPP_Species = merge(CatchPerFamilyWithinWPP,FractionSpeciesInFamilyAndWPP,all=TRUE,by=c("FamilyCommonName"))
CatchFamily_matrix = CatchPerYearWPP_Species[,grep("Catch_mT",names(CatchPerYearWPP_Species))]
PropSppInFamilyWPP = CatchPerYearWPP_Species[,grep("ProportionOfSpeciesInFamilyAndWPP",names(CatchPerYearWPP_Species))]
names(CatchFamily_matrix) = gsub("Catch_mT_","",names(CatchFamily_matrix))
CatchFamily_matrix = CatchFamily_matrix[,as.character(WPP)]
names(PropSppInFamilyWPP) = gsub("ProportionOfSpeciesInFamilyAndWPP.","",names(PropSppInFamilyWPP))
PropSppInFamilyWPP = PropSppInFamilyWPP[,as.character(WPP)]
CatchSpecies_Matrix = CatchFamily_matrix * PropSppInFamilyWPP
CatchSpecies_Modern = cbind(data.frame(Year=CatchPerYearWPP_Species$Year,FamilyCommonName=CatchPerYearWPP_Species$FamilyCommonName,species=CatchPerYearWPP_Species$species),CatchSpecies_Matrix)
#Group those species in modern catch that will not be inclued in the analysis as "Other"
CatchSpecies_Modern = merge(CatchSpecies_Modern,data.frame(species=GenusSpeciesForAnalysis_names,flag=1),all.x=TRUE,by=c("species"))
CatchSpecies_Modern$species[is.na(CatchSpecies_Modern$flag)]="Other"
CatchSpecies_Modern$flag=NULL
CatchSpecies_Modern$FamilyCommonName[CatchSpecies_Modern$species=="Other"]="Other"
CatchSpecies_Modern = aggregate.data.frame(CatchSpecies_Modern[,4:dim(CatchSpecies_Modern)[2]],by=list(CatchSpecies_Modern$species,CatchSpecies_Modern$FamilyCommonName,CatchSpecies_Modern$Year),FUN=sum,na.rm=TRUE)
names(CatchSpecies_Modern)[1:3]=c("species","FamilyCommonName","Year")
CatchSpecies_Modern = CatchSpecies_Modern[order(CatchSpecies_Modern$species,CatchSpecies_Modern$Year),]


################ Scale Historical Data to Current Split and Re-Plot Reconstructed Landings #######################################
#The reason that this historical scaling is probably needed is because  Sea Around Us Data likely also has shallow water species included, and we are only modeling the deep slope seamount species;
Years_Modern = data.frame(Year=sort(unique(CatchSpecies_Modern$Year)))
Years_Historical = data.frame(Year=sort(unique(SpeciesCatchMatrix_mT$Year)))
Years_Scaling = merge(Years_Modern,Years_Historical,by=c("Year"))
Catches_Scaling_Modern = subset(CatchSpecies_Modern,CatchSpecies_Modern$Year >= min(Years_Scaling$Year) & CatchSpecies_Modern$Year <= max(Years_Scaling$Year))
Catches_Scaling_Historical = subset(SpeciesCatchMatrix_mT,SpeciesCatchMatrix_mT$Year >= min(Years_Scaling$Year) & SpeciesCatchMatrix_mT$Year <= max(Years_Scaling$Year))
Scaling = Catches_Scaling_Modern[,4:dim(Catches_Scaling_Modern)[2]]/Catches_Scaling_Historical[,4:dim(Catches_Scaling_Historical)[2]]
Scaling = cbind(Catches_Scaling_Modern[,1:3],Scaling)
for(i in 4:dim(Catches_Scaling_Modern)[2])   #Fill in missing scaling values by (1) year for any species where just some year are missing, and then by (2) fish family for those species where all years are missing scaling
{
	Sub = subset(Scaling,select=c("species","FamilyCommonName","Year",names(Scaling)[i]))
	Sub = subset(Sub,is.na(Sub[,names(Scaling)[i]]))
	if(dim(Sub)[1]>0)   #If there are not any missing scaling values, than skip that WPP
	{
		#Check to see which species within a family may be missing scaling and take average for that family and apply it to those missing
		Scaling_Families = merge(Scaling,data.frame(FamilyCommonName=sort(unique(Sub$FamilyCommonName)),flag=1),all.x=TRUE,by=c("FamilyCommonName"))
		Scaling_Families = subset(Scaling_Families,Scaling_Families$flag==1)
		Scaling_Families = subset(Scaling_Families,select=c("species","FamilyCommonName","Year",names(Scaling_Families)[i]))
		Scaling_Families = aggregate.data.frame(Scaling_Families[,names(Scaling_Families)[4]],by=list(Scaling_Families$FamilyCommonName),FUN=mean,na.rm=TRUE)
		names(Scaling_Families) = c("FamilyCommonName","AvgScaler")
		Scaling = merge(Scaling,Scaling_Families,all.x=TRUE,by=c("FamilyCommonName"))
		Scaling[,names(Scaling)[i]][is.na(Scaling[,names(Scaling)[i]])] = Scaling$AvgScaler[is.na(Scaling[,names(Scaling)[i]])]
		Scaling$AvgScaler=NULL
	}
}
for(i in 1:dim(Scaling)[1])
{
	if(anyNA(Scaling[i,5:dim(Scaling)[2]])==TRUE)
	{
		Mean = mean(as.numeric(Scaling[i,4:dim(Scaling)[2]]),na.rm=TRUE)
		Scaling[i,4:dim(Scaling)[2]][is.na(Scaling[i,4:dim(Scaling)[2]])] = Mean
	}
}	
Scaling_meanAcrossYears = aggregate.data.frame(Scaling[,4:dim(Scaling)[2]],by=list(Scaling$FamilyCommonName,Scaling$species),FUN=mean,na.rm=TRUE)
names(Scaling_meanAcrossYears)[names(Scaling_meanAcrossYears)=="Group.1"]="FamilyCommonName"
names(Scaling_meanAcrossYears)[names(Scaling_meanAcrossYears)=="Group.2"]="species"
#Add Zeros to historical data to cover years where some species were not caught, then scale historical catch
#NOTE: If this is not done, Catch MSY will mess up
fish_id_analyzed = unique(SpeciesCatchMatrix_mT$species)
year_start = min(SpeciesCatchMatrix_mT$Year)
year_end = max(SpeciesCatchMatrix_mT$Year)
for(i in 1:length(fish_id_analyzed))
{
  Sub = SpeciesCatchMatrix_mT[SpeciesCatchMatrix_mT$species==fish_id_analyzed[i],]
  fish_id_i = unique(Sub$species)
  FamilyCommonName_i = unique(Sub$FamilyCommonName)
  Genus_i = unique(Sub$Genus)
  Species_i = unique(Sub$Species)
  toMerge = data.frame(Year=seq(year_start,year_end,by=1))
  Sub = merge(toMerge,Sub,all.x=TRUE,by=c("Year"))
  Sub$species[is.na(Sub$species)]=fish_id_i
  Sub$FamilyCommonName[is.na(Sub$FamilyCommonName)]=FamilyCommonName_i
  Sub[is.na(Sub)]=0
  if(i==1)
  {
    Answer = Sub
  }
  if(i>1)
  {
    Answer = rbind(Answer,Sub)
  }
  if(i%%10==0)
  {
    print(paste("Working on species ",i," of ",dim(fish_id_analyzed)[1],"...",sep=""))
    flush.console()
  }
}
SpeciesCatchMatrix_mT = Answer
#Now apply scaling values to historical data
YrsToUseMeanScaling = Years_Historical$Year[Years_Historical$Year %in% Years_Scaling$Year==FALSE]
Scaling_meanAcrossYears = do.call("rbind", replicate(length(YrsToUseMeanScaling),Scaling_meanAcrossYears,simplify = FALSE))
YrsToUseMeanScaling = data.frame(Year=sort(rep(YrsToUseMeanScaling,length(unique(Scaling_meanAcrossYears$species)))))
Scaling_meanAcrossYears= cbind(YrsToUseMeanScaling,Scaling_meanAcrossYears)
Scaling_meanAcrossYears = Scaling_meanAcrossYears[,names(Scaling)]
Scaling = rbind(Scaling_meanAcrossYears,Scaling)
Scaling = Scaling[,names(SpeciesCatchMatrix_mT)]
Scaling = Scaling[order(Scaling$species,Scaling$Year),]
SpeciesCatchMatrix_mT = SpeciesCatchMatrix_mT[order(SpeciesCatchMatrix_mT$species,SpeciesCatchMatrix_mT$Year),]
HistoricalCatch_mT_scaled = SpeciesCatchMatrix_mT[,4:dim(SpeciesCatchMatrix_mT)[2]] * Scaling[,4:dim(Scaling)[2]]
HistoricalCatch_mT_scaled = cbind(SpeciesCatchMatrix_mT[,1:3],HistoricalCatch_mT_scaled)
HistoricalCatch_mT_scaled_YearsNotInModern = subset(HistoricalCatch_mT_scaled,HistoricalCatch_mT_scaled$Year<min(Years_Modern$Year))
HistoricalCatch_mT_scaled_YearsNotInModern = HistoricalCatch_mT_scaled_YearsNotInModern[order(HistoricalCatch_mT_scaled_YearsNotInModern$species,HistoricalCatch_mT_scaled_YearsNotInModern$Year),]
CatchSpecies_Modern = CatchSpecies_Modern[order(CatchSpecies_Modern$species,CatchSpecies_Modern$Year),]
SpeciesCatch_mT_Scaled = rbind(HistoricalCatch_mT_scaled_YearsNotInModern,CatchSpecies_Modern)
SpeciesCatch_mT_Scaled = SpeciesCatch_mT_Scaled[order(SpeciesCatch_mT_Scaled$species,SpeciesCatch_mT_Scaled$Year),]
row.names(SpeciesCatch_mT_Scaled) = NULL
SpeciesCatch_mT_Scaled = subset(SpeciesCatch_mT_Scaled,SpeciesCatch_mT_Scaled$Year<=terminalYear)
write.table(SpeciesCatch_mT_Scaled,paste(PATH_output,"ExtrapolatedLandings.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")

################################# Plot Scaled Historical Landings #########################################
if(plotHistoricalCatches==TRUE)
{  
  SppToPlot = sort(unique(SpeciesCatch_mT_Scaled$species))
  Colors = c("blue","red","green","orange","purple","black","brown","pink","gray","gold","yellowgreen")
  for(i in 1:length(SppToPlot))
  {
	Sub = SpeciesCatch_mT_Scaled[SpeciesCatch_mT_Scaled$species==SppToPlot[i],]
	Sub$FamilyCommonName=NULL
	Sub$species=NULL
	FolderNameSpp = gsub(" ","_",SppToPlot[i])
	png(paste(PATH_historicalCatchPlots,paste(FolderNameSpp,"_Historical_Scaled.png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
	for(j in 1:length(WPP))
	{
	  if(j==1)
	  {
		plot(Sub$Year,Sub[,as.character(WPP[j])],xlab="Year",ylab="Metric Tons",type="l",lwd=2,col=Colors[j],ylim=c(0,max(Sub[,2:length(Sub)])),cex.axis=1.1,cex.lab=1.1)
		title("Reconstructed Landings")
		mtext(SppToPlot[i])
	  }
	  if(j>1)
	  {
		points(Sub$Year,Sub[,as.character(WPP[j])],type="l",lwd=2,col=Colors[j])
	  }
	}
	legend("topleft",legend=WPP,col=Colors,lwd=2)
	dev.off()
  }
}


######################## Compute Updated SPR Values and Develop New Input File #############################################
if(computeSPR_createInputFile==TRUE)
{
	#Read fecundity values
	counter=1
	WPP_withEEZ = c(as.character(WPP),"EEZ")
	for(i in 1:length(GenusSpeciesForAnalysis_names))
	{
		for(j in 1:length(WPP_withEEZ))
		{
			Params_i = LifeHistoryScenariosList$CalculatedFromData[LifeHistoryScenariosList$CalculatedFromData$species==GenusSpeciesForAnalysis_names[i],]
			histBreaks=NA
			if(WPP_withEEZ[j]!="EEZ")
			{
				Sub = Lengths[Lengths$species==GenusSpeciesForAnalysis_names[i] & Lengths$WPP==as.numeric(WPP_withEEZ[j]),]
				histBreaks = subset(OptimalHistogramBreaks,OptimalHistogramBreaks$species==GenusSpeciesForAnalysis_names[i])[,paste("X",WPP_withEEZ[j],sep="")]
			}
			if(WPP_withEEZ[j]=="EEZ")
			{
				Sub = Lengths[Lengths$species==GenusSpeciesForAnalysis_names[i],]
				histBreaks = subset(OptimalHistogramBreaks,OptimalHistogramBreaks$species==GenusSpeciesForAnalysis_names[i])[,WPP_withEEZ[j]]
			}
			if(is.na(histBreaks))
			{
				next   #Skip that WPP for that species if there are not enough samples to do the analysis; we will know this by the fact that the OptimalHistogramBreaks value for that species and WPP will be NA
			}
			#Use LBSPR Instead to calculate initial F-Based
			MyPars = new("LB_pars")
			MyPars@Species = GenusSpeciesForAnalysis_names[i]
			MyPars@Linf = Params_i$Linf
			MyPars@L50 = Params_i$Lmat #length at 50% maturity
			MyPars@L95 = Params_i$Lmat * 1.27 #length at 95% maturity - we don't know this. From looking at Gulf species, it seems to be approximately 1.27 times the Length at 50% maturity. Apply this assumption here. 
			MyPars@MK = Params_i$M/Params_i$VBK
			MyPars@M = Params_i$M
			MyPars@L_units = "cm"
			MyPars@BinMin=(floor(min(Sub$cm,na.rm=TRUE)/histBreaks))*histBreaks
			MyPars@BinWidth=histBreaks
			MyPars@BinMax=ceiling(max(c(Params_i$Lmax,Params_i$Linf,max(Sub$cm,na.rm=TRUE)),na.rm=TRUE)/histBreaks)*histBreaks
			Len1 = new("LB_lengths",LB_pars=MyPars,file=Sub$cm,sataType="raw")
			myFit1 = LBSPRfit(MyPars, Len1)
			F_est = as.numeric(myFit1@Ests[,"FM"])*Params_i$M
			SPR = as.numeric(myFit1@Ests[,"SPR"])
			if(counter==1)
			{
				Answer = data.frame(species=GenusSpeciesForAnalysis_names[i],WPP=WPP_withEEZ[j],spr_current=SPR,F_est=F_est,M=Params_i$M)
			}
			if(counter>1)
			{
				Answer = rbind(Answer,data.frame(species=GenusSpeciesForAnalysis_names[i],WPP=WPP_withEEZ[j],spr_current=SPR,F_est=F_est,M=Params_i$M))
			}
			counter = counter + 1
			if(i%%5==0)
			{
				print(paste("Working on species ",i," of ",length(GenusSpeciesForAnalysis_names),"...",sep=""))
				flush.console()
			}
		}
	}
	Answer$spr_current[Answer$spr_current>1]=1
	InputFile = merge(data.frame(Scenario=CatchMSY_scenarioNames),Answer)
	InputFile = InputFile[c("Scenario","species","WPP","spr_current","F_est","M")]
	InputFile$l0_low=0.8
	InputFile$l0_up=1.0  
	InputFile$l0_low[InputFile$Scenario=="Constant"]=InputFile$spr_current[InputFile$Scenario=="Constant"]-0.1
	InputFile$l0_up[InputFile$Scenario=="Constant"]=InputFile$spr_current[InputFile$Scenario=="Constant"]+0.1
	InputFile$l0_low[InputFile$l0_low<0]=0
	InputFile$l0_up[InputFile$l0_up>1]=1  
	InputFile$l0_step=0
	InputFile$lt_low=InputFile$spr_current-0.1
	InputFile$lt_low[InputFile$lt_low<0]=0
	InputFile$lt_up=InputFile$spr_current+0.1
	InputFile$lt_up[InputFile$lt_up>1]=1
	InputFile$lt_refyr=max(SpeciesCatch_mT_Scaled$Year)
	InputFile$sigv=0
	InputFile$r_dist="unif"
	InputFile$r_low=0.05
	InputFile$r_up=1
	InputFile$r_mean=0
	InputFile$r_sd=0
	InputFile$nsims=30000
	InputFile$grout=2
	InputFile = InputFile[order(InputFile$species,InputFile$Scenario,InputFile$WPP),]
	write.table(InputFile,paste(PATH_otherInput,"StartingParametersCatchMSY.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
}



#######################################################################################################
########################### Compute Biomass and MSY Using Catch-MSY ####################################
########################################################################################################
#Read in The functions we need:
#Schaefer Logistic Growth Model
.schaefer	= function(theta)
{
  with(as.list(theta), {  ## for all combinations of ri & ki
    bt=vector()
    ell = 0  ## initialize ell
    for (j in startbt)
    {
      if(ell == 0)
      {
        bt[1]=j*k*exp(rnorm(1,0, sigR))  ## set biomass in first year
        for(i in 1:nyr) ## for all years in the time series
        {
          xt=rnorm(1,0, sigR)
          bt[i+1]=(bt[i]+r*bt[i]*(1-bt[i]/k)-ct[i])*exp(xt) ## calculate biomass as function of previous year's biomass plus net production minus catch
        }
        #Bernoulli likelihood, assign 0 or 1 to each combination of r and k
        ell = 0
        if(bt[nyr+1]/k>=lam1 && bt[nyr+1]/k <=lam2 && min(bt) > 0 && max(bt) <=k && bt[which(yr==interyr)]/k>=interbio[1] && bt[which(yr==interyr)]/k<=interbio[2])
          ell = 1
      }
    }
    return(list(ell=ell,bt=bt))
  })
}
#Stock Reduction Analysis for N trials, with argument "theta" - a list object containing: r (lower and upper bounds for r), K (lower and upper bounds for k), and lambda (limits for current depletion)
sraMSY = function(theta, N)
{
  with(as.list(theta),
  {
     ri = exp(runif(N, log(r[1]), log(r[2])))  ## get N values between r[1] and r[2], assign to ri
     ki = exp(runif(N, log(k[1]), log(k[2])))  ## get N values between k[1] and k[2], assing to ki
     itheta=cbind(r=ri,k=ki, lam1=lambda[1],lam2=lambda[2], sigR=sigR) ## assign ri, ki, and final biomass range to itheta
     M = apply(itheta,1,.schaefer) ## call Schaefer function with parameters in itheta
     i=1:N
     ## prototype objective function
     get.ell=function(i) M[[i]]$ell
     ell = sapply(i, get.ell)
     get.bt=function(i) M[[i]]$bt
     bt = sapply(i, get.bt)
     return(list(r=ri,k=ki,ell=ell,bt=bt))
  })
}
fileLineWriteCounter=1
StartValues = read.table(paste(PATH_otherInput,CatchMSY_inputFileName,sep="/"),sep=",",header=TRUE,stringsAsFactors=FALSE)







if(runCatchMSY==TRUE)
{

	CatchMSY_outfileName_original = CatchMSY_outfileName
	i=4     #1-4
	G=1   
	CatchMSY_outfileName = paste("CatchMSY_Output_",CatchMSY_scenarioNames[i],".csv",sep="")


  #for(i in 1:length(CatchMSY_scenarioNames))
  {
    if(CatchMSY_scenarioNames[i]=="Historical")
    {
      Catch = SpeciesCatch_mT_Scaled
    }
    if(CatchMSY_scenarioNames[i]=="Optimistic")
    {
	  Years = data.frame(Year=seq((max(as.numeric(SpeciesCatch_mT_Scaled$Year))-cMSY_optimisticScenario_yrsToCurrentDepletion),max(SpeciesCatch_mT_Scaled$Year),by=1))
	  Catch_modern = subset(SpeciesCatch_mT_Scaled,SpeciesCatch_mT_Scaled$Year>min(Years_Modern$Year))
	  Catch_max = aggregate.data.frame(Catch_modern[,4:dim(Catch_modern)[2]],by=list(Catch_modern$species,Catch_modern$FamilyCommonName),FUN=max,na.rm=TRUE)
	  names(Catch_max)[names(Catch_max)=="Group.1"]="species"
	  names(Catch_max)[names(Catch_max)=="Group.2"]="FamilyCommonName"
	  Catch = merge(Years,Catch_max)
    }
    if(CatchMSY_scenarioNames[i]=="Pessimistic")
    {
	  Years = data.frame(Year=seq((max(as.numeric(SpeciesCatch_mT_Scaled$Year))-cMSY_pessimisticScenario_yrsToCurrentDepletion),max(SpeciesCatch_mT_Scaled$Year),by=1))
	  Catch_modern = subset(SpeciesCatch_mT_Scaled,SpeciesCatch_mT_Scaled$Year>min(Years_Modern$Year))
	  Catch_max = aggregate.data.frame(Catch_modern[,4:dim(Catch_modern)[2]],by=list(Catch_modern$species,Catch_modern$FamilyCommonName),FUN=max,na.rm=TRUE)
	  names(Catch_max)[names(Catch_max)=="Group.1"]="species"
	  names(Catch_max)[names(Catch_max)=="Group.2"]="FamilyCommonName"
	  Catch = merge(Years,Catch_max)
    }
    if(CatchMSY_scenarioNames[i]=="Constant")
    {
	  Years = data.frame(Year=seq((max(as.numeric(SpeciesCatch_mT_Scaled$Year))-cMSY_constantScenario_yrsToCurrentDepletion),max(SpeciesCatch_mT_Scaled$Year),by=1))
	  Catch_modern = subset(SpeciesCatch_mT_Scaled,SpeciesCatch_mT_Scaled$Year>min(Years_Modern$Year))
	  Catch_max = aggregate.data.frame(Catch_modern[,4:dim(Catch_modern)[2]],by=list(Catch_modern$species,Catch_modern$FamilyCommonName),FUN=max,na.rm=TRUE)
	  names(Catch_max)[names(Catch_max)=="Group.1"]="species"
	  names(Catch_max)[names(Catch_max)=="Group.2"]="FamilyCommonName"
	  Catch = merge(Years,Catch_max)
    }
    for(G in 1:length(LandingsUncertaintyScenarios))
    {
      Catch[as.character(WPP)] = Catch[as.character(WPP)]*LandingsUncertaintyScenarios[G]   #Scale catch to various sensitivity levels
      yr = sort(unique(Catch$Year))  #years in time series
      nyr = length(yr)    #number of years in the time series
      for(j in 1:length(GenusSpeciesForAnalysis_names))
      {
         Catch_j = Catch[Catch$species==GenusSpeciesForAnalysis_names[j],]
         for(p in 1:(length(WPP)+1))
         {
            if(p!=(length(WPP)+1))
            {
              currentSPRvalue = (StartValues[StartValues$WPP==as.character(WPP[p]) & StartValues$species==GenusSpeciesForAnalysis_names[j] & StartValues$Scenario==CatchMSY_scenarioNames[i],])$spr_current
            }
            if(p==(length(WPP)+1))
            {
              currentSPRvalue = (StartValues[StartValues$WPP=="EEZ" & StartValues$species==GenusSpeciesForAnalysis_names[j] & StartValues$Scenario==CatchMSY_scenarioNames[i],])$spr_current
            }
            if(length(currentSPRvalue)==0)
            {
              next
            }
            if(currentSPRvalue<0.1 | currentSPRvalue>=0.1)
            {
              if(p!=(length(WPP)+1))
              {
                ct = Catch_j[,as.character(WPP[p])] #assumes that catch is given in tonnes
                Params_j = StartValues[StartValues$Scenario==CatchMSY_scenarioNames[i] &
                           StartValues$species==GenusSpeciesForAnalysis_names[j] & StartValues$WPP==as.character(WPP[p]),]
              }
              if(p==(length(WPP)+1))
              {
                 CatchAllWpp = subset(Catch_j,select=-c(Year,species,FamilyCommonName))
                 ct = rowSums(CatchAllWpp)
                Params_j = StartValues[StartValues$Scenario==CatchMSY_scenarioNames[i] &
                           StartValues$species==GenusSpeciesForAnalysis_names[j] & StartValues$WPP=="EEZ",]
              }
              if(sum(ct)==0)
              {
                sigR=0.05
                output = data.frame(CatchMSY_scenarioNames[i],LandingsUncertaintyScenarios[G],GenusSpeciesForAnalysis_names[j],
                    as.character(WPP[p]),sigR,min(yr),max(yr),max(ct),ct[1],
                    ct[nyr],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"No Catch")
                 if(fileLineWriteCounter==1)
                 {
                    colNames = c("Scenario","LandingsSensitivity","species","WPP","sigR","min_yr","max_yr","max_catch","catch_t0","catch_t","NumSims","BiomassMean_Year1",
                       "BiomassUpper_Year1","BiomassLower_Year1","BiomassMean_LastYear","BiomassUpper_LastYear","BiomassLower_LastYear",
                       "r_mean","r_stdev","r_upper","r_lower", "k_mean","k_stdev","k_upper","k_lower","msy_mean","msy_stdev",
                       "msy_lower","msy_upper","Notes")
                    write.table(transpose(data.frame(colNames)), file = paste(PATH_output,CatchMSY_outfileName,sep="/"),sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
                 }
                 write.table(output, file = paste(PATH_output,CatchMSY_outfileName,sep="/"), append = TRUE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
                 fileLineWriteCounter=fileLineWriteCounter+1
              }
              if(sum(ct)>0)
              {
                start_r  = c(Params_j$r_low,Params_j$r_up)  #starting R values
                start_k = c(max(ct),50*max(ct)) #using defaults for lower and upper k (NOT using the file values) e.g. 100 * max catch as per using Martell's paper guesses
                startbio = c(Params_j$l0_low,Params_j$l0_up)
                #settings that ignore intermediate year
                interyr = yr[2]   ## ignore!
                interbio = c(0,1) ## ignore!
                finalbio = c(Params_j$lt_low,Params_j$lt_up)
                #parameters about the simulation
                n = Params_j$nsims
                sigR = 0.05      ## process error; 0 if deterministic model; 0.05 reasonable value? 0.2 is too high
                startbt = seq(startbio[1],startbio[2],by=0.05) ## apply range of start biomass in steps of 0.05
                parbound = list(r=start_r,k=start_k,lambda=finalbio,sigR=sigR)
                #Stock reduction analysis
                R1 = sraMSY(parbound,n)
                ## Get statistics on r, k, MSY and determine new bounds for r and k
                r1 = R1$r[R1$ell==1]
                k1 = R1$k[R1$ell==1]
                msy1 = r1*k1/4
                mean_msy1 = exp(mean(log(msy1)))
                max_k1a = min(k1[r1<1.1*parbound$r[1]]) ## smallest k1 near initial lower bound of r
                max_k1b = max(k1[r1*k1/4<mean_msy1]) ## largest k1 that gives mean MSY
                max_k1 = if(max_k1a < max_k1b) {max_k1a} else {max_k1b}
                if(length(r1)<10)
                {
                  output = data.frame(CatchMSY_scenarioNames[i],LandingsUncertaintyScenarios[G],GenusSpeciesForAnalysis_names[j],
                    as.character(WPP[p]),sigR,min(yr),max(yr),max(ct),ct[1],
                    ct[nyr],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,paste("Too few (", length(r1), ") possible r-k combinations, check input parameters",sep=""))
                  if(fileLineWriteCounter==1)
                  {
                    colNames = c("Scenario","LandingsSensitivity","species","WPP","sigR","min_yr","max_yr","max_catch","catch_t0","catch_t","NumSims","BiomassMean_Year1",
                       "BiomassUpper_Year1","BiomassLower_Year1","BiomassMean_LastYear","BiomassUpper_LastYear","BiomassLower_LastYear",
                       "r_mean","r_stdev","r_upper","r_lower", "k_mean","k_stdev","k_upper","k_lower","msy_mean","msy_stdev",
                       "msy_lower","msy_upper","Notes")
                    write.table(transpose(data.frame(colNames)), file = paste(PATH_output,CatchMSY_outfileName,sep="/"),sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
                  }
                  write.table(output, file = paste(PATH_output,CatchMSY_outfileName,sep="/"), append = TRUE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
                  fileLineWriteCounter=fileLineWriteCounter+1
                  cat("Too few (", length(r1), ") possible r-k combinations, check input parameters","\n")
                  flush.console()
                }
                if(length(r1)>=10)
                {
                  ## set new upper bound of r to 1.2 max r1
                  parbound$r[2] = 1.2*max(r1)
                  ## set new lower bound for k to 0.9 min k1 and upper bound to max_k1
                  parbound$k = c(0.9 * min(k1), max_k1)
                  cat("First MSY =", format(mean_msy1, digits=3),"\n")
                  cat("First r =", format(exp(mean(log(r1))), digits=3),"\n")
                  cat("New upper bound for r =", format(parbound$r[2],digits=2),"\n")
                  cat("New range for k =", format(parbound$k[1], digits=3), "-", format(parbound$k[2],digits=3),"\n")
                  #Repeat analysis with new r-k bounds
                  R1 = sraMSY(parbound, n)
                  #Get statistics on r, k and msy
                  r = R1$r[R1$ell==1]
                  r_mean = exp(mean(log(r)))
                  r_sd = exp(sd(log(r)))
                  r_upper = r_mean + r_sd
                  r_lower = r_mean - r_sd
                  k = R1$k[R1$ell==1]
                  k_mean = exp(mean(log(k)))
                  k_sd = exp(sd(log(k)))
                  k_upper = k_mean + k_sd
                  k_lower = k_mean - k_sd
                  msy = r * k / 4
                  msy_mean = exp(mean(log(msy)))
                  msy_sd = exp(sd(log(msy)))
                  msy_upper = msy_mean + msy_sd
                  msy_lower = msy_mean - msy_sd
                  mean_ln_msy = mean(log(msy))
                  bt = data.frame(t(R1$bt),ell=R1$ell)
                  bt = bt[bt$ell==1,]
                  bt$ell=NULL
                  bt_mean = colMeans(bt)
                  bt_stdev = apply(bt,2,sd)
                  bt_upper = bt_mean + bt_stdev
                  bt_lower = bt_mean - bt_stdev
                  BiomassEstimates = data.frame(Year=c(0,yr),Mean=bt_mean,Upper=bt_upper,Lower=bt_lower)
                  row.names(BiomassEstimates)=NULL
                  write.table(bt,paste(PATH_catchMSY_diagnostics,paste(paste("BiomassEstimates_RAW",GenusSpeciesForAnalysis_names[j],CatchMSY_scenarioNames[i],as.character(WPP[p]),paste("LandingSens",LandingsUncertaintyScenarios[G],sep=""),sep="_"),".csv",sep=""),sep="/"),
                      row.names=FALSE,col.names=TRUE,sep=",")
                  write.table(BiomassEstimates,paste(PATH_catchMSY_diagnostics,paste(paste("BiomassEstimates",GenusSpeciesForAnalysis_names[j],CatchMSY_scenarioNames[i],as.character(WPP[p]),paste("LandingSens",LandingsUncertaintyScenarios[G],sep=""),sep="_"),".csv",sep=""),sep="/"),
                      row.names=FALSE,col.names=TRUE,sep=",")
                  ## plot MSY over catch data
                  if(plotCatchMSY_diagnostics==TRUE)
                  {
                    FolderNameSpp = gsub(" ","_",GenusSpeciesForAnalysis_names[j])				
                    png(paste(PATH_catchMSY_diagnostics,paste(paste("Cmsy",FolderNameSpp,CatchMSY_scenarioNames[i],as.character(WPP[p]),paste("LandingSens",LandingsUncertaintyScenarios[G],sep=""),sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
                    par(mfcol=c(2,3))
                    plot(yr, ct, type="l", ylim = c(0, max(ct,exp(mean_ln_msy + 2 * sd(log(msy))))), xlab = "Year", ylab = "Catch (mT)", main = paste(GenusForAnalysis_names[j],SpeciesForAnalysis_names[j],sep=" "),cex.axis=1.1,cex.lab=1.1)
                    if(p!=(length(WPP)+1))
                    {
                      mtext(paste(CatchMSY_scenarioNames[i],", WPP=",WPP[p],sep=""),cex=0.7)
                    }
                    if(p==(length(WPP)+1))
                    {
                      mtext(paste(CatchMSY_scenarioNames[i],", WPP=All",sep=""),cex=0.7)
                    }
                    abline(h=exp(mean(log(msy))),col="red",lwd=2)
                    abline(h=exp(mean_ln_msy - 2 * sd(log(msy))),col="red",lty=2)
                    abline(h=exp(mean_ln_msy + 2 * sd(log(msy))),col="red",lty=2)
                    hist(r, freq=F, xlim=c(0, max(c(1.2*max(r),exp(mean(log(r))+2*sd(log(r)))))), main = paste(GenusForAnalysis_names[j],SpeciesForAnalysis_names[j],sep=" "),cex.axis=1.1,cex.lab=1.1)
                    if(p!=(length(WPP)+1))
                    {
                      mtext(paste(CatchMSY_scenarioNames[i],", WPP=",WPP[p],sep=""),cex=0.7)
                    }
                    if(p==(length(WPP)+1))
                    {
                      mtext(paste(CatchMSY_scenarioNames[i],", WPP=All",sep=""),cex=0.7)
                    }
                    abline(v=exp(mean(log(r))),col="red",lwd=2)
                    abline(v=exp(mean(log(r))-2*sd(log(r))),col="red",lty=2)
                    abline(v=exp(mean(log(r))+2*sd(log(r))),col="red",lty=2)
                    plot(r1, k1, xlim = start_r, ylim = start_k, xlab="r", ylab="k (mT)",cex.axis=1.1,cex.lab=1.1)
                    title(GenusSpeciesForAnalysis_names[j])
                    if(p!=(length(WPP)+1))
                    {
                      mtext(paste(CatchMSY_scenarioNames[i],", WPP=",WPP[p],sep=""),cex=0.7)
                    }
                    if(p==(length(WPP)+1))
                    {
                      mtext(paste(CatchMSY_scenarioNames[i],", WPP=All",sep=""),cex=0.7)
                    }
                    hist(k, freq=F, xlim=c(0, max(c(1.2*max(k),exp(mean(log(k))+2*sd(log(k)))))), xlab="k (mT)", main = paste(GenusForAnalysis_names[j],SpeciesForAnalysis_names[j],sep=" "),cex.axis=1.1,cex.lab=1.1)
                    if(p!=(length(WPP)+1))
                    {
                      mtext(paste(CatchMSY_scenarioNames[i],", WPP=",WPP[p],sep=""),cex=0.7)
                    }
                    if(p==(length(WPP)+1))
                    {
                      mtext(paste(CatchMSY_scenarioNames[i],", WPP=All",sep=""),cex=0.7)
                    }
                    abline(v=exp(mean(log(k))),col="red", lwd=2)
                    abline(v=exp(mean(log(k))-2*sd(log(k))),col="red",lty=2)
                    abline(v=exp(mean(log(k))+2*sd(log(k))),col="red",lty=2)
                    plot(log(r), log(k),xlab="ln(r)",ylab="ln(k)",cex.axis=1.1,cex.lab=1.1)
                    title(paste(GenusForAnalysis_names[j],SpeciesForAnalysis_names[j],sep=" "))
                    if(p!=(length(WPP)+1))
                    {
                      mtext(paste(CatchMSY_scenarioNames[i],", WPP=",WPP[p],sep=""),cex=0.7)
                    }
                    if(p==(length(WPP)+1))
                    {
                      mtext(paste(CatchMSY_scenarioNames[i],", WPP=All",sep=""),cex=0.7)
                    }
                    abline(v=mean(log(r)))
                    abline(h=mean(log(k)))
                    abline(mean(log(msy))+log(4),-1, col="red",lwd=2)
                    abline(mean(log(msy))-2*sd(log(msy))+log(4),-1, col="red")
                    abline(mean(log(msy))+2*sd(log(msy))+log(4),-1, col="red")
                    hist(msy, freq=F, xlim=c(0, 1.2 * max(msy)), xlab="MSY (mT)",main = GenusSpeciesForAnalysis_names[j],cex.axis=1.1,cex.lab=1.1)
                    if(p!=(length(WPP)+1))
                    {
                      mtext(paste(CatchMSY_scenarioNames[i],", WPP=",WPP[p],sep=""),cex=0.7)
                    }
                    if(p==(length(WPP)+1))
                    {
                      mtext(paste(CatchMSY_scenarioNames[i],", WPP=All",sep=""),cex=0.7)
                    }
                    abline(v=exp(mean(log(msy))),col="red", lwd=2)
                    abline(v=exp(mean_ln_msy - 2 * sd(log(msy))),col="red")
                    abline(v=exp(mean_ln_msy + 2 * sd(log(msy))),col="red")
                    dev.off()
                    cat("Possible combinations = ", length(r),"\n")
                    cat("geom. mean r =", format(exp(mean(log(r))),digits=3), "\n")
                    cat("r +/- 2 SD =", format(exp(mean(log(r))-2*sd(log(r))),digits=3),"-",format(exp(mean(log(r))+2*sd(log(r))),digits=3), "\n")
                    cat("geom. mean k =", format(exp(mean(log(k))),digits=3), "\n")
                    cat("k +/- 2 SD =", format(exp(mean(log(k))-2*sd(log(k))),digits=3),"-",format(exp(mean(log(k))+2*sd(log(k))),digits=3), "\n")
                    cat("geom. mean MSY =", format(exp(mean(log(msy))),digits=3),"\n")
                    cat("MSY +/- 2 SD =", format(exp(mean_ln_msy - 2 * sd(log(msy))),digits=3), "-", format(exp(mean_ln_msy + 2 * sd(log(msy))),digits=3), "\n")
                  }
                  ## Write results into outfile, in append mode
                  if(p!=(length(WPP)+1))
                  {
                    output = data.frame(CatchMSY_scenarioNames[i],LandingsUncertaintyScenarios[G],GenusSpeciesForAnalysis_names[j],
                      as.character(WPP[p]),sigR,min(yr),max(yr),max(ct),ct[1],
                      ct[nyr],length(r),BiomassEstimates$Mean[1],BiomassEstimates$Upper[1],BiomassEstimates$Lower[1],
                      BiomassEstimates$Mean[dim(BiomassEstimates)[1]],BiomassEstimates$Upper[dim(BiomassEstimates)[1]],BiomassEstimates$Lower[dim(BiomassEstimates)[1]],
                      r_mean,r_sd,r_upper,r_lower,k_mean,k_sd,k_upper,k_lower,msy_mean,msy_sd,msy_upper,msy_lower,"Normal Run")
                  }
                  if(p==(length(WPP)+1))
                  {
                    output = data.frame(CatchMSY_scenarioNames[i],LandingsUncertaintyScenarios[G],GenusSpeciesForAnalysis_names[j],
                      "EEZ",sigR,min(yr),max(yr),max(ct),ct[1],
                      ct[nyr],length(r),BiomassEstimates$Mean[1],BiomassEstimates$Upper[1],BiomassEstimates$Lower[1],
                      BiomassEstimates$Mean[dim(BiomassEstimates)[1]],BiomassEstimates$Upper[dim(BiomassEstimates)[1]],BiomassEstimates$Lower[dim(BiomassEstimates)[1]],
                      r_mean,r_sd,r_upper,r_lower,k_mean,k_sd,k_upper,k_lower,msy_mean,msy_sd,msy_upper,msy_lower,"Normal Run")
                  }
                  if(fileLineWriteCounter==1)
                  {
                      colNames = c("Scenario","LandingsSensitivity","species","WPP","sigR","min_yr","max_yr","max_catch","catch_t0","catch_t","NumSims","BiomassMean_Year1",
                         "BiomassUpper_Year1","BiomassLower_Year1","BiomassMean_LastYear","BiomassUpper_LastYear","BiomassLower_LastYear",
                         "r_mean","r_stdev","r_upper","r_lower", "k_mean","k_stdev","k_upper","k_lower","msy_mean","msy_stdev",
                         "msy_lower","msy_upper","Notes")
                       write.table(transpose(data.frame(colNames)), file = paste(PATH_output,CatchMSY_outfileName,sep="/"),sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
                  }
                  write.table(output, file = paste(PATH_output,CatchMSY_outfileName,sep="/"), append = TRUE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
                  fileLineWriteCounter=fileLineWriteCounter+1
                }
              }
           }
        }
      }
    }
  }
  
  
  
  
}






if(ranCatchMSYonSeparateProcessors==TRUE)
{
	CatchMSY_Output_Historical = read.table(paste(PATH_output,"CatchMSY_Output_Historical.csv",sep="/"),header=TRUE,sep=",")
	CatchMSY_Output_Pessimistic = read.table(paste(PATH_output,"CatchMSY_Output_Pessimistic.csv",sep="/"),header=TRUE,sep=",")
	CatchMSY_Output_Optimistic = read.table(paste(PATH_output,"CatchMSY_Output_Optimistic.csv",sep="/"),header=TRUE,sep=",")
	CatchMSY_Output_Constant = read.table(paste(PATH_output,"CatchMSY_Output_Constant.csv",sep="/"),header=TRUE,sep=",")
	CatchMSY_Output = rbind(CatchMSY_Output_Historical,CatchMSY_Output_Pessimistic,CatchMSY_Output_Optimistic,CatchMSY_Output_Constant)
	write.table(CatchMSY_Output,paste(PATH_output,"CatchMSY_Output.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
}

if(processCatchMSYoutput==TRUE)
{
	#Need to modify the code below
	CmsyResults = read.table(paste(PATH_output,"CatchMSY_Output.csv",sep="/"),header=TRUE,sep=",")
	BiomassCarryingCapacityMSY_mean = subset(CmsyResults,select=c(Scenario,species,LandingsSensitivity,WPP,catch_t,BiomassMean_LastYear,r_mean,k_mean,msy_mean,Notes))
	BiomassCarryingCapacityMSY_mean$Bound="mean"
	BiomassCarryingCapacityMSY_lower = subset(CmsyResults,select=c(Scenario,species,LandingsSensitivity,WPP,catch_t,BiomassLower_LastYear,r_lower,k_lower,msy_lower,Notes))
	BiomassCarryingCapacityMSY_lower$Bound="lower"
	BiomassCarryingCapacityMSY_upper = subset(CmsyResults,select=c(Scenario,species,LandingsSensitivity,WPP,catch_t,BiomassUpper_LastYear,r_upper,k_upper,msy_upper,Notes))
	BiomassCarryingCapacityMSY_upper$Bound="upper"
	names(BiomassCarryingCapacityMSY_lower) = names(BiomassCarryingCapacityMSY_mean)
	names(BiomassCarryingCapacityMSY_upper) = names(BiomassCarryingCapacityMSY_mean)
	BiomassCarryingCapacityMSY = rbind(BiomassCarryingCapacityMSY_mean,BiomassCarryingCapacityMSY_lower,BiomassCarryingCapacityMSY_upper)
#	BiomassCarryingCapacityMSY$Genus = as.character(BiomassCarryingCapacityMSY$Genus)
#	BiomassCarryingCapacityMSY$Species = as.character(BiomassCarryingCapacityMSY$Species)
#	BiomassCarryingCapacityMSY$Scenario = as.character(BiomassCarryingCapacityMSY$Scenario)
#	BiomassCarryingCapacityMSY$WPP = as.character(BiomassCarryingCapacityMSY$WPP)
#	BiomassCarryingCapacityMSY$Notes = as.character(BiomassCarryingCapacityMSY$Notes)
	BiomassCarryingCapacityMSY = BiomassCarryingCapacityMSY[order(BiomassCarryingCapacityMSY$species,BiomassCarryingCapacityMSY$WPP,BiomassCarryingCapacityMSY$Scenario,BiomassCarryingCapacityMSY$Bound),]
	BiomassCarryingCapacityMSY = BiomassCarryingCapacityMSY[BiomassCarryingCapacityMSY$catch_t>0,]
	BiomassCarryingCapacityMSY = BiomassCarryingCapacityMSY[BiomassCarryingCapacityMSY$BiomassMean_LastYear>0,]
	BiomassCarryingCapacityMSY$Scenario_New = paste(BiomassCarryingCapacityMSY$Scenario,BiomassCarryingCapacityMSY$Bound,BiomassCarryingCapacityMSY$LandingsSensitivity,sep="-")
	BiomassCarryingCapacityMSY = subset(BiomassCarryingCapacityMSY,BiomassCarryingCapacityMSY$species!="")
	write.table(BiomassCarryingCapacityMSY,paste(PATH_output,"BiomassCarryingCapacityMSY.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
}

################################ Make Plots of Catch MSY Results #############################################
if(plotCatchMSY_results==TRUE)
{
  BiomassCarryingCapacityMSY = read.table(paste(PATH_output,"BiomassCarryingCapacityMSY.csv",sep="/"),header=TRUE,sep=",")
  SpeciesAreas_ForPlots = subset(BiomassCarryingCapacityMSY,select=c(species,WPP))
  SpeciesAreas_ForPlots = unique(SpeciesAreas_ForPlots)
  for(i in 1:dim(SpeciesAreas_ForPlots)[1])
  {
     BiomassCarryingCapacityMSY_i_j = BiomassCarryingCapacityMSY[BiomassCarryingCapacityMSY$species==SpeciesAreas_ForPlots$species[i] & BiomassCarryingCapacityMSY$WPP==SpeciesAreas_ForPlots$WPP[i],]
     BiomassCarryingCapacityMSY_i_j = BiomassCarryingCapacityMSY_i_j[order(BiomassCarryingCapacityMSY_i_j$Scenario_New),]
     FolderNameSpp = gsub(" ","_",SpeciesAreas_ForPlots$species[i]) 
     png(paste(paste(PATH_catchMSY_plots,FolderNameSpp,sep="/"),paste(paste("CatchMSY_Results",gsub(" ","_",SpeciesAreas_ForPlots$species[i]),SpeciesAreas_ForPlots$WPP[i],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=4200,res=600)
     par(mfrow=c(2,2),mar=c(9, 4, 4, 2) + 0.1) 
     plot(as.factor(c(1:dim(BiomassCarryingCapacityMSY_i_j)[1])),BiomassCarryingCapacityMSY_i_j$BiomassMean_LastYear,ylab="Current Biomass (mT)",main=SpeciesAreas_ForPlots$species[i],xaxt="n")
     axis(1,at=1:dim(BiomassCarryingCapacityMSY_i_j)[1],labels=BiomassCarryingCapacityMSY_i_j$Scenario_New,las=2,cex.axis=0.8)
     mtext(SpeciesAreas_ForPlots$WPP[i])
     plot(as.factor(c(1:dim(BiomassCarryingCapacityMSY_i_j)[1])),BiomassCarryingCapacityMSY_i_j$r_mean,ylab="r (population growth rate)",main=SpeciesAreas_ForPlots$species[i],xaxt="n")
     axis(1,at=1:dim(BiomassCarryingCapacityMSY_i_j)[1],labels=BiomassCarryingCapacityMSY_i_j$Scenario_New,las=2,cex.axis=0.8)
     mtext(SpeciesAreas_ForPlots$WPP[i])
     plot(as.factor(c(1:dim(BiomassCarryingCapacityMSY_i_j)[1])),BiomassCarryingCapacityMSY_i_j$k_mean,ylab="Carrying Capacity (mT)",main=SpeciesAreas_ForPlots$species[i],xaxt="n")
     axis(1,at=1:dim(BiomassCarryingCapacityMSY_i_j)[1],labels=BiomassCarryingCapacityMSY_i_j$Scenario_New,las=2,cex.axis=0.8)
     mtext(SpeciesAreas_ForPlots$WPP[i])
     plot(as.factor(c(1:dim(BiomassCarryingCapacityMSY_i_j)[1])),BiomassCarryingCapacityMSY_i_j$msy_mean,ylab="MSY (mT)",main=SpeciesAreas_ForPlots$species[i],xaxt="n")
     axis(1,at=1:dim(BiomassCarryingCapacityMSY_i_j)[1],labels=BiomassCarryingCapacityMSY_i_j$Scenario_New,las=2,cex.axis=0.8)
     abline(h=BiomassCarryingCapacityMSY_i_j$catch_t,col="red",lwd=2,lty=2)
     mtext(SpeciesAreas_ForPlots$WPP[i])
     dev.off()
     if(i%%5==0)
     {
        print(paste("Working on item ",i," of ",dim(SpeciesAreas_ForPlots)[1],"...",sep=""))
        flush.console()
     }
   }
}


############# Function: Beverton and Holt Recruitment as in Stock Synthesis (no phi) #############
BevHoltRecruitment = function(R0,h_steepness,SSB,SSB_virgin)
{
   Recruits =  (4*h_steepness*R0*SSB)/(SSB_virgin*(1-h_steepness)+((5*h_steepness)-1)*SSB)
   return(Recruits)
}

############## Function: Length-Based Logistic Selectivity #################################
SelectivityLogistic_Length = function(SizeBinMidPointsVector,Selex_Logistic_start_a,Selex_Logistic_start_b)
{
  Selex = 1/(1+exp(-1*log(19)*((SizeBinMidPointsVector-Selex_Logistic_start_a)/Selex_Logistic_start_b)))
  return(Selex)
}
SelectivityLogistic_Age = function(AgeVector,Selex_Logistic_start_a,Selex_Logistic_start_b)
{
  Selex = 1/(1+exp(-1*log(19)*((AgeVector-Selex_Logistic_start_a)/Selex_Logistic_start_b)))
  return(Selex)
}
BaranovCatch = function(AgeVector,LengthVector,PopulationVector,Mnat,F_mort,Selex)
{
  CatchVector = PopulationVector*(1-exp(-1*((F_mort*Selex)+Mnat)))*((F_mort*Selex)/((F_mort*Selex)+Mnat))
  return(CatchVector)
}

############## Function: Compute Fecundity ########################
#ComputeFecundity = function(Weight_at_size_vector,Weight_at_size_vector_min,Weight_at_size_vector_max,FecundityMin,FecundityMax)
#{
#  Fec = data.frame(Weight=c(Weight_at_size_vector_min,Weight_at_size_vector_max),Fecundity=c(FecundityMin,FecundityMax))
#  LinMod = lm(Fec$Fecundity ~ Fec$Weight)
#  Fecundity = (coef(LinMod)[[1]] + Weight_at_size_vector*coef(LinMod)[[2]])
#  Fecundity[Fecundity<0]=0
#  return(Fecundity)
#}

######################### Function: Find R0 ####################################
findR0 = function(M,Linf,VBK,a,b,K,maxAge,matureAge,init_min_GT,init_max_GT,incrementBy)
{
  age = c(0:maxAge)
  data = data.frame(age=age)
  data$maturity = 1
  data$maturity[1:matureAge] = 0 #maturity is knife-edge. only age mature=3 and above are matured
  natmort = exp(-M) # this is per year, we can make it per month
  data$length = (1-exp(-data$age*VBK))*Linf #this is the age-length relationship   - length is in cm
  data$weight = (a*(data$length^b)) #this is the length-weight relationship - weight is in kilograms
  runtime = 1000 #iterations to get equilibrium
  Nmat = maxAge+1 #matriz size includes the zero age
  PopINI = matrix(100, ncol = 1, nrow = Nmat) #initial population with 100 individuals per age class
  Leslie = matrix(0, ncol = Nmat, nrow = Nmat)
  x = rep(natmort,Nmat)
  y = diag(x)
  Leslie[c(2:Nmat),] = y[c(1:Nmat-1),]
  Leslie[1,] = 0 #this is the number of eggs produced per age class. This will not matter as we will adjust GK.
  # ----GET GK. Gk is the number of recruits needed to bring the population to K, assuming that there will be no fishing
  Biomass = vector()
  SSB = vector()
  Pop = PopINI
  #the range here needs a bit of guestimate
  tryValuesR0 = function(Leslie,Pop,init_min_GT,init_max_GT,incrementBy)
  {
    for(GK in seq(init_min_GT,init_max_GT,by=incrementBy))
    {
      count<-0
      for (t in 1:runtime)
      {
        count<-count+1
        Pop2<-floor(Leslie%*%Pop)
        Pop2[1,1]<-GK
        Pop<-Pop2
        Biomass[count]<-sum(Pop*data$weight)/(1000) #from kg to mt
        SSB[count]<-(sum(floor(data$maturity*Pop)*data$weight))/(1000) #from kg to mt;
      }
      indicator<-tail(Biomass,n=1)
      if (indicator >K){break}
    }
    return(GK)
  }
  GK =  tryValuesR0(Leslie,Pop,init_min_GT,init_max_GT,incrementBy)
  while(incrementBy>1)
  {
    init_min_GT = GK-incrementBy
    init_max_GT = GK+incrementBy
    incrementBy = round(incrementBy/10)
    if(incrementBy==0)  {incrementBy=1}
    GK =  tryValuesR0(Leslie,Pop,init_min_GT,init_max_GT,incrementBy)
  }
  #K #this is the carrying capacity from catch-MSY
  #GK #this is the number of recruits needed to reach K
  #indicator #if we have a constant recruitment assumption where recruitment == GK, this is the biomass of the system in the absense of fishing
  return(GK)
}
computeLengthSPR = function(Mnat,Linf,VBK,var_a,var_b,Lmat,F_estimated,LengthVec)
{
     age = c(1:100)
     popdf = data.frame(age)
     popdf$popF0[1] = 1000*exp(-Mnat)
     for(i in 2:100)
     {
       popdf$popF0[i] = popdf$popF0[i-1]*exp(-Mnat)
     }
     popdf$TL = Linf*(1-exp(-VBK*age))
     popdf$BW = var_a*popdf$TL^var_b
     popdf$BiomassF0 = popdf$BW*popdf$popF0
     popdf$Mature = ifelse(popdf$TL < Lmat,0,1)
     popdf$SpBiomassF0 = popdf$BiomassF0*popdf$Mature
     SpBiomassF0 = sum(popdf$SpBiomassF0)
     tmp = table(as.vector(LengthVec))
     lc = min(as.numeric(names(tmp)[tmp==max(tmp)]))
     lc_sub = quantile(LengthVec, c(.01), names=FALSE)
     popdf$Fact = ifelse(popdf$TL < lc_sub,0,
                    ifelse(popdf$TL < lc & popdf$TL>= lc_sub,(popdf$TL*F_estimated/(lc-lc_sub))-(lc_sub*F_estimated/(lc-lc_sub)),F_estimated))
     popdf$popFact[1] = 1000*exp(-(Mnat+popdf$Fact[1]))
     for(i in 2:100)
     {
       popdf$popFact[i] = popdf$popFact[i-1]*exp(-(Mnat+popdf$Fact[i])*1)
     }
     popdf$BiomassFact = popdf$BW*popdf$popFact
     popdf$SpBiomassFact = popdf$BiomassFact*popdf$Mature
     SpBiomassFact = sum(popdf$SpBiomassFact)
     BiomassFact = sum(popdf$BiomassFact)
     SPR = SpBiomassFact/SpBiomassF0
     return(SPR)
}
computeLengthSPR_selex = function(Mnat,Linf,VBK,var_a,var_b,Lmat,F_estimated,Selex_Logistic_a_age,Selex_Logistic_b_age)
{
     age = c(1:100)
     popdf = data.frame(age)
     popdf$popF0[1] = 1000*exp(-Mnat)
     for(i in 2:100)
     {
       popdf$popF0[i] = popdf$popF0[i-1]*exp(-Mnat)
     }
     popdf$TL = Linf*(1-exp(-VBK*age))
     popdf$BW = var_a*popdf$TL^var_b
     popdf$BiomassF0 = popdf$BW*popdf$popF0
     popdf$Mature = ifelse(popdf$TL < Lmat,0,1)
     popdf$SpBiomassF0 = popdf$BiomassF0*popdf$Mature
     SpBiomassF0 = sum(popdf$SpBiomassF0)
     popdf$Fact = SelectivityLogistic_Age(popdf$age,Selex_Logistic_a_age,Selex_Logistic_b_age)*F_estimated
     popdf$popFact[1] = 1000*exp(-(Mnat+popdf$Fact[1]))
     for(i in 2:100)
     {
       popdf$popFact[i] = popdf$popFact[i-1]*exp(-(Mnat+popdf$Fact[i])*1)
     }
     popdf$BiomassFact = popdf$BW*popdf$popFact
     popdf$SpBiomassFact = popdf$BiomassFact*popdf$Mature
     SpBiomassFact = sum(popdf$SpBiomassFact)
     BiomassFact = sum(popdf$BiomassFact)
     SPR = SpBiomassFact/SpBiomassF0
     return(SPR)
}

gc()





#################################################################################################################
############# Reconstruct Length Demographics & Calculate Current F #############################################
#################################################################################################################

#4743 scenarios, 591 on 8 processors: 1-591, 592-1183, 1184-1775, 1776-2367, 2368-2959, 2960-3551, 3552-4143, 4144-4734



if(runLengthReconstruction==TRUE)
{

	fName_ListOfPopulations = "ListOfPopulations8.RData"

	LifeHistoryScenariosList$FromFishBase$t0=0
  BiomassCarryingCapacityMSY = read.table(paste(PATH_output,"BiomassCarryingCapacityMSY.csv",sep="/"),header=TRUE,sep=",")
  ListOfPopulations=list()
  names(OptimalHistogramBreaks) = gsub("X","",names(OptimalHistogramBreaks))
  CanAssessSppWPP = CanAssess
  names(CanAssessSppWPP) = gsub("_toAssess","",names(CanAssessSppWPP))
  CanAssessSppWPP$species = row.names(CanAssessSppWPP)
  row.names(CanAssessSppWPP)=NULL
  
  #IMPORTANT NOTE: set this to "FALSE", because someone in the ministry of fisheries is doing their PhD work on this spp in this area. 
 CanAssessSppWPP[CanAssessSppWPP$species=="Lutjanus malabaricus",]["713"]=FALSE

  #for(i in 1:dim(BiomassCarryingCapacityMSY)[1])   
  #There are 4740 scenarios from CatchMSY
#  for(i in 1:592)      
#  for(i in 593:1185)   
#  for(i in 1186:1778)   
#  for(i in 1779:2371)   
#  for(i in 2372:2964)   
#  for(i in 2965:3557)   
#  for(i in 3558:4150)   
  for(i in 4151:4740)   
  {
    if(BiomassCarryingCapacityMSY$WPP[i]!="EEZ")
    {
      Sub = Lengths[Lengths$species==BiomassCarryingCapacityMSY$species[i] & Lengths$WPP==BiomassCarryingCapacityMSY$WPP[i],]
      StartValues_i_j_o_p = StartValues[StartValues$Scenario==BiomassCarryingCapacityMSY$Scenario[i] & StartValues$species==BiomassCarryingCapacityMSY$species[i] & StartValues$WPP==BiomassCarryingCapacityMSY$WPP[i],]
    }
    if(BiomassCarryingCapacityMSY$WPP[i]=="EEZ")
    {
      Sub = Lengths[Lengths$species==BiomassCarryingCapacityMSY$species[i],]
      StartValues_i_j_o_p = StartValues[StartValues$Scenario==BiomassCarryingCapacityMSY$Scenario[i] & StartValues$species==BiomassCarryingCapacityMSY$species[i] & StartValues$WPP=="All",]
    }
    if(subset(CanAssessSppWPP,CanAssessSppWPP$species==BiomassCarryingCapacityMSY$species[i])[,BiomassCarryingCapacityMSY$WPP[i]]==FALSE)     #Only run if sample size greater than 150 samples as determined much higher in script
    {
      next
    }
	histBreaks = subset(OptimalHistogramBreaks,OptimalHistogramBreaks$species==BiomassCarryingCapacityMSY$species[i])[,BiomassCarryingCapacityMSY$WPP[i]]
    for(p in 1:length(LifeHistoryScenariosList))
    {
	  LifeHistoryScenarios_p = subset(LifeHistoryScenariosList[[p]],LifeHistoryScenariosList[[p]]$species==BiomassCarryingCapacityMSY$species[i])
	  #Start things length-based
      Sub$Age = suppressWarnings(log(1-(Sub$cm/LifeHistoryScenarios_p$Linf))/(-1*LifeHistoryScenarios_p$VBK))
      Sub = Sub[Sub$Age!=Inf,]
      Sub = Sub[!is.na(Sub$Age),]
      Sub$Age_Int = floor(Sub$Age)
	  Sub$Length_Bin_lower = (floor(Sub$cm/histBreaks))*histBreaks
	  Sub$Length_Bin_upper = Sub$Length_Bin_lower + histBreaks
	  Sub$LengthBin_Midpoint = (Sub$Length_Bin_lower + Sub$Length_Bin_upper)/2
	  NumAtLength_i_j = data.frame(table(Sub$LengthBin_Midpoint))
	  names(NumAtLength_i_j) = c("Midpoint","NumSampled")
	  NumAtLength_i_j$Midpoint = as.numeric(as.character(NumAtLength_i_j$Midpoint))	  
	  Midpoints = data.frame(Midpoint=seq((histBreaks/2),((ceiling(max(c(LifeHistoryScenarios_p$Lmax,LifeHistoryScenarios_p$Linf),na.rm=TRUE)/histBreaks)*histBreaks)-(histBreaks/2)),by=histBreaks))
	  NumAtLength_i_j = merge(Midpoints,NumAtLength_i_j,all.x=TRUE,by=c("Midpoint"))
	  NumAtLength_i_j[is.na(NumAtLength_i_j)]=0
	  NumAtLength_i_j$Length_min = NumAtLength_i_j$Midpoint-(histBreaks/2)
	  NumAtLength_i_j$Length_max = NumAtLength_i_j$Midpoint+(histBreaks/2)
	  NumAtLength_i_j$Age_min = suppressWarnings((log(1-(NumAtLength_i_j$Length_min/(LifeHistoryScenarios_p$Linf)))/(-1*LifeHistoryScenarios_p$VBK))+LifeHistoryScenarios_p$t0)
	  NumAtLength_i_j$Age_max = suppressWarnings((log(1-(NumAtLength_i_j$Length_max/(LifeHistoryScenarios_p$Linf)))/(-1*LifeHistoryScenarios_p$VBK))+LifeHistoryScenarios_p$t0)
	  NumAtLength_i_j$Weight_kg_min = (LifeHistoryScenarios_p$var_a*(NumAtLength_i_j$Length_min^LifeHistoryScenarios_p$var_b))
	  NumAtLength_i_j$Weight_kg_max = (LifeHistoryScenarios_p$var_a*(NumAtLength_i_j$Length_max^LifeHistoryScenarios_p$var_b))
	  NumAtLength_i_j$Weight_kg = (LifeHistoryScenarios_p$var_a*(NumAtLength_i_j$Midpoint^LifeHistoryScenarios_p$var_b))
	  NumAtLength_i_j$AgeDifference_Years = NumAtLength_i_j$Age_max - NumAtLength_i_j$Age_min
	  NumAtLength_i_j$AgeDifference_Years[NumAtLength_i_j$AgeDifference_Years==Inf]=NA
	  library(zoo)
	  NumAtLength_i_j$AgeDifference_Years = na.locf(NumAtLength_i_j$AgeDifference_Years)
      for(u in 1:length(NumAtLength_i_j$AgeDifference_Years))
      {
        if(is.nan(NumAtLength_i_j$Age_max[u]) & !is.nan(NumAtLength_i_j$Age_min[u]))
        {
          NumAtLength_i_j$Age_max[u] = NumAtLength_i_j$Age_min[u] + NumAtLength_i_j$AgeDifference_Years[u]
        }
        if(is.nan(NumAtLength_i_j$Age_max[u]) & is.nan(NumAtLength_i_j$Age_min[u]))
        {
          NumAtLength_i_j$Age_min[u] = NumAtLength_i_j$Age_max[u-1]
          NumAtLength_i_j$Age_max[u] = NumAtLength_i_j$Age_min[u] + NumAtLength_i_j$AgeDifference_Years[u]
        }
        if(NumAtLength_i_j$Age_min[u]==Inf)
        {
          NumAtLength_i_j$Age_min[u] = NumAtLength_i_j$Age_max[u-1]
        }
        if(NumAtLength_i_j$Age_max[u]==Inf)
        {
          NumAtLength_i_j$Age_max[u] = NumAtLength_i_j$Age_min[u] + NumAtLength_i_j$AgeDifference_Years[u]
        }
      }
      NumAtLength_i_j$Age = (NumAtLength_i_j$Age_min+NumAtLength_i_j$Age_max)/2
      NumAtLength_i_j$AgeDifference_Days = round(NumAtLength_i_j$AgeDifference_Years*365)
      #Set up Num At Age
	   maxAge = max(floor(NumAtLength_i_j$Age))
      NumAtAge_i_j = data.frame(Age=seq(0,maxAge,by=1))
      NumAtAge_i_j$Length = LifeHistoryScenarios_p$Linf*(1-exp(-LifeHistoryScenarios_p$VBK*(NumAtAge_i_j$Age)))
      NumAtAge_i_j$Weight_kg = LifeHistoryScenarios_p$var_a*(NumAtAge_i_j$Length^LifeHistoryScenarios_p$var_b)
      NumAtAge_i_j$Maturity=0
      NumAtAge_i_j$Maturity[NumAtAge_i_j$Length>LifeHistoryScenarios_p$Lmat]=1
      #NumAtAge_i_j$Fecundity = ComputeFecundity(NumAtAge_i_j$Weight_kg,min(NumAtAge_i_j$Weight_kg[NumAtAge_i_j$Maturity==1]),max(NumAtAge_i_j$Weight_kg[NumAtAge_i_j$Maturity==1]),StartValues_i_j_o_p$Fecundity_min,StartValues_i_j_o_p$Fecundity_max)
      NumAtAge_i_j$Survival_M = exp(-LifeHistoryScenarios_p$M*(NumAtAge_i_j$Age))
      NumAtAge_i_j$Survival_M_rescaled = NumAtAge_i_j$Survival_M/sum(NumAtAge_i_j$Survival_M)
      NumAtAge_i_j$WeightFractionVirgin = NumAtAge_i_j$Weight_kg*NumAtAge_i_j$Survival_M_rescaled
      Weight_kg_oneFish = sum(NumAtAge_i_j$WeightFractionVirgin)
      NumFish_Virgin = (BiomassCarryingCapacityMSY$k_mean[i]*1000)/Weight_kg_oneFish
      NumAtAge_i_j$PopulationVirgin = NumAtAge_i_j$Survival_M_rescaled*NumFish_Virgin
      NumSampledAtAge = data.frame(table(Sub$Age_Int))
      names(NumSampledAtAge) = c("Age","NumSampled")
      NumAtAge_i_j = merge(NumAtAge_i_j,NumSampledAtAge,all.x=TRUE,by=c("Age"))
      NumAtAge_i_j[is.na(NumAtAge_i_j)]=0
      NumAtAge_i_j$SampledFraction = NumAtAge_i_j$NumSampled/sum(NumAtAge_i_j$NumSampled)
      CatchNumber = (BiomassCarryingCapacityMSY$catch_t[i]*1000)/Weight_kg_oneFish
      NumAtAge_i_j$ExtrapolatedCatchInNumber = NumAtAge_i_j$SampledFraction*CatchNumber
      NumAtAge_i_j$SSB_virginAtAge = NumAtAge_i_j$PopulationVirgin*(NumAtAge_i_j$Weight_kg*0.001)*NumAtAge_i_j$Maturity
      SSB_virgin = sum(NumAtAge_i_j$SSB_virginAtAge)
      AgeMat = round(log(1-(LifeHistoryScenarios_p$Lmat/LifeHistoryScenarios_p$Linf))/(-1*LifeHistoryScenarios_p$VBK))
      print("Finding R0...")
      flush.console()
      init_max_GT=((BiomassCarryingCapacityMSY$k_mean[i]*1000)/Weight_kg_oneFish)
      digits = nchar(round(init_max_GT,0))
      incrementBy = 10^(digits-3)
      if(digits<=3) {incrementBy=1}
      R0 = findR0(M=LifeHistoryScenarios_p$M,Linf=LifeHistoryScenarios_p$Linf,VBK=LifeHistoryScenarios_p$VBK,a=LifeHistoryScenarios_p$var_a,b=LifeHistoryScenarios_p$var_b,K=BiomassCarryingCapacityMSY$k_mean[i],maxAge=ceiling(max(NumAtAge_i_j$Age)),matureAge=AgeMat,init_min_GT=1,init_max_GT=init_max_GT,incrementBy=incrementBy)          #lower and upper starting values for seach are 10 and 50  percent of the total population number as extrapolated from catch MSY
      print("R0 found.")
      flush.console()
      #Compute Logistic Selectivity at Age
      NumAtAge_i_j$Selex_Empirical = NumAtAge_i_j$SampledFraction/max(NumAtAge_i_j$SampledFraction)
      booleanVec=FALSE
      for(g in 1:length(NumAtAge_i_j$Selex_Empirical))
      {
         if(NumAtAge_i_j$Selex_Empirical[g]==1) {booleanVec=TRUE}
         if(booleanVec==TRUE) {NumAtAge_i_j$Selex_Empirical[g]=1}
      }
      #Develop Catch Age-Length Key
      Catch_ALK = as.data.frame.matrix(table(Sub$LengthBin_Midpoint,Sub$Age_Int))
      Catch_ALK$Midpoint = row.names(Catch_ALK)
      row.names(Catch_ALK)=NULL
      Midpoints = subset(NumAtLength_i_j,select=c("Midpoint"))
      Catch_ALK = merge(Midpoints,Catch_ALK,all.x=TRUE,by=c("Midpoint"))
      row.names(Catch_ALK) = Catch_ALK$Midpoint
      Catch_ALK$Midpoint=NULL
      Ages=as.character(seq(0,maxAge,by=1))
      for(g in 1:length(Ages))
      {
        isInNames=FALSE
        for(gg in 1:length(names(Catch_ALK)))
        {
          if(as.character(Ages[g])==names(Catch_ALK)[gg])
          {
            isInNames=TRUE
          }
        }
        if(isInNames==FALSE)
        {
           toAdd = data.frame(rep(0,times=length(Catch_ALK[,1])))
           names(toAdd) = Ages[g]
           Catch_ALK = cbind(Catch_ALK,toAdd)
        }
      }
      Catch_ALK = subset(Catch_ALK,select=as.character(Ages))
      Catch_ALK[is.na(Catch_ALK)]=0
      Catch_ALK_Sum_rep = do.call("rbind", replicate(length(Catch_ALK[,1]),colSums(Catch_ALK),simplify = FALSE))
      Catch_ALK = Catch_ALK/Catch_ALK_Sum_rep
      Catch_ALK[is.na(Catch_ALK)]=0
      #Develop Population Age-Length Key: Length to Age
      ALK_Len2Age = data.frame(matrix(nrow=length(NumAtLength_i_j$Midpoint),ncol=(maxAge+1),by=1))
      names(ALK_Len2Age) = as.character(seq(0,maxAge,by=1))
      row.names(ALK_Len2Age) = NumAtLength_i_j$Midpoint
      ALK_Len2Age[is.na(ALK_Len2Age)]=0
      for(u in 1:length(NumAtLength_i_j$Midpoint))
      {
         IntAge_min = floor(NumAtLength_i_j$Age_min[u])
         IntAge_max = floor(NumAtLength_i_j$Age_max[u])
         if(IntAge_max>maxAge)
         {
           IntAge_max=maxAge
         }
         if(IntAge_min==IntAge_max)
         {
            ALK_Len2Age[as.character(NumAtLength_i_j$Midpoint[u]),as.character(IntAge_max)] = 1
         }
         if(IntAge_min!=IntAge_max)
         {
            IntAge_max = IntAge_max
            Diff = seq(IntAge_min,IntAge_max,by=1)
            Proportion = c()
            for(x in 1:length(Diff))
            {
              if(x==1)
              {
                 Prop = Diff[x+1]-NumAtLength_i_j$Age_min[u]
                 Proportion = c(Proportion,Prop)
              }
              if(x>1 & x<length(Diff))
              {
                Prop=1
                Proportion = c(Proportion,Prop)
              }
              if(x==length(Diff))
              {
                 Prop = NumAtLength_i_j$Age_max[u]-Diff[x]
                 Proportion = c(Proportion,Prop)
              }
            }
            Proportion = Proportion/sum(Proportion)
            ALK_Len2Age[as.character(NumAtLength_i_j$Midpoint[u]),as.character(Diff)] = Proportion
         }
      }
      ALK_Len2Age[is.na(ALK_Len2Age)]=0
      #Develop Age-Length Key: Age to Length
      ALK_Age2Len = data.frame(matrix(nrow=length(NumAtLength_i_j$Midpoint),ncol=(maxAge+1),by=1))
      names(ALK_Age2Len) = as.character(seq(0,maxAge,by=1))
      row.names(ALK_Age2Len) = NumAtLength_i_j$Midpoint
      ALK_Age2Len[is.na(ALK_Age2Len)]=0
      for(u in 1:length(NumAtLength_i_j$Midpoint))
      {
         IntAge_min = floor(NumAtLength_i_j$Age_min[u])
         IntAge_max = floor(NumAtLength_i_j$Age_max[u])
         if(IntAge_max>maxAge)
         {
           IntAge_max=maxAge
         }
         if(IntAge_min==IntAge_max)
         {
            ALK_Age2Len[as.character(NumAtLength_i_j$Midpoint[u]),as.character(IntAge_max)] = 1
         }
         if(IntAge_min!=IntAge_max)
         {
            Diff = seq(IntAge_min,IntAge_max,by=1)
            Proportion = c()
            for(x in 1:length(Diff))
            {
              if(x==1)
              {
                 Prop = ceiling(NumAtLength_i_j$Age_min[u]) - NumAtLength_i_j$Age_min[u]
					if(Prop==0)
					{
						Prop=1   #If both values are zero, then set this equal to 1 so that all age 0 fish are assigned to the smallest size group. Otherwise, no fish are assigned to age zero. 
					}
                 Proportion = c(Proportion,Prop)
              }
              if(x>1 & x<length(Diff))
              {
                Prop=1
                Proportion = c(Proportion,Prop)
              }
              if(x==length(Diff))
              {
                 Prop = NumAtLength_i_j$Age_max[u]-floor(NumAtLength_i_j$Age_max[u])
                 Proportion = c(Proportion,Prop)
              }
            }
            ALK_Age2Len[as.character(NumAtLength_i_j$Midpoint[u]),as.character(Diff)] = Proportion
         }
      }
	  #For the last age class and the last length bin, if it does not equal to one, set it so that it does	  
	  LastAgeClass = ALK_Age2Len[,as.character(max(as.numeric(names(ALK_Age2Len))))]
	  LastAgeClass_sum = sum(LastAgeClass)
	  if(LastAgeClass_sum!=1)
	  {
			LastAgeClass[length(LastAgeClass)] = 1-sum(LastAgeClass[-length(LastAgeClass)])
			ALK_Age2Len[,as.character(max(as.numeric(names(ALK_Age2Len))))] = LastAgeClass
	  }
      ALK_Age2Len[is.na(ALK_Age2Len)]=0
      ColSums_ALK_Age2Len = transpose(data.frame(colSums(ALK_Age2Len)))
      ColSums_ALK_Age2Len = do.call("rbind", replicate(length(ALK_Age2Len$"0"),ColSums_ALK_Age2Len,simplify = FALSE))
      ALK_Age2Len = ALK_Age2Len/ColSums_ALK_Age2Len
      #Compute virgin population size at length using age-length key
      PopVirginAtAge_rep = do.call("rbind", replicate(length(ALK_Age2Len[,1]),transpose(data.frame(NumAtAge_i_j$PopulationVirgin)),simplify = FALSE))
      NumAtLength_Virgin = PopVirginAtAge_rep*ALK_Age2Len
      NumAtLength_Virgin = rowSums(NumAtLength_Virgin)
      NumAtLength_Virgin = data.frame(Midpoint = NumAtLength_i_j$Midpoint,PopulationVirgin=NumAtLength_Virgin)
      NumAtLength_i_j = merge(NumAtLength_i_j,NumAtLength_Virgin,all.x=TRUE,by=c("Midpoint"))
      ExtrapCatch_Rep = transpose(do.call("cbind", replicate(dim(Catch_ALK)[1],data.frame(NumAtAge_i_j$ExtrapolatedCatchInNumber),simplify = FALSE)))
      ExtrapolatedCatchAtLength = Catch_ALK*ExtrapCatch_Rep
      ExtrapolatedCatchAtLength = data.frame(Midpoint=row.names(ExtrapolatedCatchAtLength),ExtrapolatedCatch=rowSums(ExtrapolatedCatchAtLength))
      NumAtLength_i_j = merge(NumAtLength_i_j,ExtrapolatedCatchAtLength,all.x=TRUE,by=c("Midpoint"))
      #Compute Logistic Selectivity at Length
      NumAtLength_i_j$Selex_Empirical = NumAtLength_i_j$NumSampled/max(NumAtLength_i_j$NumSampled)
      booleanVec=FALSE
      for(g in 1:length(NumAtLength_i_j$Selex_Empirical))
      {
         if(NumAtLength_i_j$Selex_Empirical[g]==1) {booleanVec=TRUE}
         if(booleanVec==TRUE) {NumAtLength_i_j$Selex_Empirical[g]=1}
      }
      #Fit Selectivity at age to empirical data
      Selex_Logistic_a_age = 5
      Selex_Logistic_b_age = 10
      Params = c(Selex_Logistic_a_age,Selex_Logistic_b_age)
      min.RSS = function(Params)
      {
        Selex = SelectivityLogistic_Age(NumAtAge_i_j$Age,Params[1],Params[2])
        return(sum((Selex - NumAtAge_i_j$Selex_Empirical)^2))
      }
      Model = optim(par=Params,fn=min.RSS,control=list(trace=0,pgtol=0.001,maxit=10000))
      Selex_Logistic_a_age = Model$par[1]
      Selex_Logistic_b_age = Model$par[2]
      NumAtAge_i_j$Selex_Fitted = SelectivityLogistic_Age(NumAtAge_i_j$Age,Selex_Logistic_a_age,Selex_Logistic_b_age)
      FolderNameSpp = gsub(" ","_",BiomassCarryingCapacityMSY$species[i])
      png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("SelectivityAtAge",gsub(" ","_",BiomassCarryingCapacityMSY$species[i]),BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
      plot(NumAtAge_i_j$Age,NumAtAge_i_j$Selex_Empirical,pch=16,xlab="Age",ylab="Selectivity",cex.axis=1.1,cex.lab=1.1)
      points(NumAtAge_i_j$Age,NumAtAge_i_j$Selex_Fitted,type="l",lwd=2,col="red")
      title(BiomassCarryingCapacityMSY$species[i])
      mtext(paste(as.character(BiomassCarryingCapacityMSY$WPP[i]),paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],sep=", "))
      dev.off()
      #Fit Selectivity at length to empirical data
      Selex_Logistic_a_length = 33
      Selex_Logistic_b_length = 8
      Params = c(Selex_Logistic_a_length,Selex_Logistic_b_length)
      min.RSS = function(Params)
      {
        Selex = SelectivityLogistic_Length(NumAtLength_i_j$Midpoint,Params[1],Params[2])
        return(sum((Selex - NumAtLength_i_j$Selex_Empirical)^2))
      }
      Model = optim(par=Params,fn=min.RSS,control=list(trace=0,pgtol=0.001,maxit=10000))
      Selex_Logistic_a_length = Model$par[1]
      Selex_Logistic_b_length = Model$par[2]
      NumAtLength_i_j$Selex_Fitted = SelectivityLogistic_Length(NumAtLength_i_j$Midpoint,Selex_Logistic_a_length,Selex_Logistic_b_length)
      FolderNameSpp = gsub(" ","_",BiomassCarryingCapacityMSY$species[i])
      png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("SelectivityAtLength",gsub(" ","_",BiomassCarryingCapacityMSY$species[i]),BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
      plot(NumAtLength_i_j$Midpoint,NumAtLength_i_j$Selex_Empirical,pch=16,xlab="Length (cm)",ylab="Selectivity",cex.axis=1.1,cex.lab=1.1)
      points(NumAtLength_i_j$Midpoint,NumAtLength_i_j$Selex_Fitted,type="l",lwd=2,col="red")
      title(BiomassCarryingCapacityMSY$species[i])
      mtext(paste(as.character(BiomassCarryingCapacityMSY$WPP[i]),paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],sep=", "))
      dev.off()
      #Create Refreshed numbers at age data.frames for each scenario
      NumAtLength_v = NumAtLength_i_j
      NumAtAge_v = NumAtAge_i_j
      for(v in 1:length(populationReconstructionMethods))
      {
        if(populationReconstructionMethods[v]=="LBSPR")
        {
          #Use LBSPR
          NumAtLength_i_j = NumAtLength_v
          NumAtAge_i_j = NumAtAge_v
          print("Reconstructing population size using LBSPR")
          flush.console()
          MyPars = new("LB_pars")
          MyPars@Species = BiomassCarryingCapacityMSY$species[i]
          MyPars@Linf = LifeHistoryScenarios_p$Linf
          MyPars@L50 = LifeHistoryScenarios_p$Lmat #length at 50% maturity
          MyPars@L95 = LifeHistoryScenarios_p$Lmat * 1.27  #length at 95% maturity - From looking at Gulf species, this is the approximate relationship
          MyPars@MK = LifeHistoryScenarios_p$M/LifeHistoryScenarios_p$VBK
          MyPars@M = LifeHistoryScenarios_p$M
          MyPars@L_units = "cm"
          MyPars@BinMin=min(NumAtLength_i_j$Midpoint)-(histBreaks/2)
          MyPars@BinWidth=histBreaks
          MyPars@BinMax=ceiling(max(c(LifeHistoryScenarios_p$Lmax,LifeHistoryScenarios_p$Linf,max(Sub$cm,na.rm=TRUE)),na.rm=TRUE)/histBreaks)*histBreaks
          Len1 = new("LB_lengths",LB_pars=MyPars,file=Sub$cm,sataType="raw")
          myFit1 = LBSPRfit(MyPars, Len1)
			F_estimated = as.numeric(myFit1@Ests[,"FM"])*LifeHistoryScenarios_p$M 
			NumAtAge_i_j$Survival_F_M = exp(-1*(LifeHistoryScenarios_p$M+(F_estimated*NumAtAge_i_j$Selex_Empirical))*NumAtAge_i_j$Age)     #"pseudosurvival"
			NumAtAge_i_j$Survival_F_M_Rescaled = NumAtAge_i_j$Survival_F_M/sum(NumAtAge_i_j$Survival_F_M)
			EstNumFish_Population = (BiomassCarryingCapacityMSY$BiomassMean_LastYear[i]*1000)/Weight_kg_oneFish
			NumAtAge_i_j$Population_Extrapolated = NumAtAge_i_j$Survival_F_M_Rescaled*EstNumFish_Population
			NumAtAge_i_j$EstimatedCatch_Baranov = BaranovCatch(NumAtAge_i_j$Age,NumAtAge_i_j$Length,NumAtAge_i_j$Population_Extrapolated,LifeHistoryScenarios_p$M,F_estimated,NumAtAge_i_j$Selex_Empirical)
			NumAtAge_i_j$Residuals = NumAtAge_i_j$ExtrapolatedCatchInNumber-NumAtAge_i_j$EstimatedCatch_Baranov
          FolderNameSpp = gsub(" ","_",BiomassCarryingCapacityMSY$species[i])
          png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("CatchAtLength",FolderNameSpp,BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
          plot(NumAtLength_i_j$Midpoint,NumAtLength_i_j$ExtrapolatedCatch,pch=16,xlab="Length (cm)",ylab="Catch in Number",cex.axis=1.1,cex.lab=1.1,ylim=c(0,max(NumAtLength_i_j$EstimatedCatch_Baranov,NumAtLength_i_j$ExtrapolatedCatch)))
          points(NumAtLength_i_j$Midpoint,NumAtLength_i_j$EstimatedCatch_Baranov,type="l",lwd=2,col="red")
          title(BiomassCarryingCapacityMSY$species[i])
          mtext(paste(BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],populationReconstructionMethods[v],sep=", "))
          dev.off()
			NumAtAge_i_j$SSB_CurrentAtAge = NumAtAge_i_j$Population_Extrapolated*(NumAtAge_i_j$Weight_kg)*NumAtAge_i_j$Maturity
			SSB_current = sum(NumAtAge_i_j$SSB_CurrentAtAge)
			SPR_current = SSB_current/SSB_virgin
			NumAgeCurrent = do.call("rbind",replicate(dim(ALK_Age2Len)[1],transpose(data.frame(NumAtAge_i_j$Population_Extrapolated)),simplify = FALSE))
			NumLengthCurrent = NumAgeCurrent*ALK_Age2Len
			NumAtLength_i_j$Population_Extrapolated = rowSums(NumLengthCurrent)
			CatchCurrent = do.call("rbind",replicate(dim(ALK_Age2Len)[1],transpose(data.frame(NumAtAge_i_j$EstimatedCatch_Baranov)),simplify = FALSE))
			CatchCurrent = CatchCurrent*ALK_Age2Len
			NumAtLength_i_j$EstimatedCatch_Baranov = rowSums(CatchCurrent)
          FolderNameSpp = gsub(" ","_",BiomassCarryingCapacityMSY$species[i])
          png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("CatchAtAge",FolderNameSpp,BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
          plot(NumAtAge_i_j$Age,NumAtAge_i_j$ExtrapolatedCatchInNumber,pch=16,xlab="Age",ylab="Catch in Number",cex.axis=1.1,cex.lab=1.1,ylim=c(0,max(NumAtAge_i_j$EstimatedCatch_Baranov,NumAtAge_i_j$ExtrapolatedCatchInNumber)))
          points(NumAtAge_i_j$Age,NumAtAge_i_j$EstimatedCatch_Baranov,type="l",lwd=2,col="red")
          title(BiomassCarryingCapacityMSY$species[i])
          mtext(paste(BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],populationReconstructionMethods[v],sep=", "))
          dev.off()
		   SPR_lengthBased = as.numeric(myFit1@Ests[,"SPR"])
        }
        #Find population that reconstructs the population at the current SPR level as compared with virgin population
        if(populationReconstructionMethods[v]=="SolveBaranov")
        {
          NumAtLength_i_j = NumAtLength_v
          NumAtAge_i_j = NumAtAge_v
          print("Reconstructing population size solving Baranov's Catch Equation")
          flush.console()
          F_start = 1.0
          min.RSS = function(F_start)
          {
            NumAtAge_i_j$Survival_F_M = exp(-1*(LifeHistoryScenarios_p$M+(F_start*NumAtAge_i_j$Selex_Empirical))*NumAtAge_i_j$Age)
            NumAtAge_i_j$Survival_F_M_Rescaled = NumAtAge_i_j$Survival_F_M/sum(NumAtAge_i_j$Survival_F_M)
            EstNumFish_Population = (BiomassCarryingCapacityMSY$BiomassMean_LastYear[i]*1000)/Weight_kg_oneFish
            NumAtAge_i_j$Population_Extrapolated = NumAtAge_i_j$Survival_F_M_Rescaled*EstNumFish_Population
            NumAtAge_i_j$SSB_CurrentAtAge = NumAtAge_i_j$Population_Extrapolated*(NumAtAge_i_j$Weight_kg*0.001)*NumAtAge_i_j$Maturity
            SSB_current = sum(NumAtAge_i_j$SSB_CurrentAtAge)
            SPR_current = SSB_current/SSB_virgin
            return(sum((SPR_current - StartValues_i_j_o_p$spr_current)^2))
          }
          Model = optimize(min.RSS,c(0,10),tol=0.001)
          NumAtAge_i_j$Survival_F_M = exp(-1*(LifeHistoryScenarios_p$M+(Model$minimum[1]*NumAtAge_i_j$Selex_Empirical))*NumAtAge_i_j$Age)     #"pseudosurvival"
          NumAtAge_i_j$Survival_F_M_Rescaled = NumAtAge_i_j$Survival_F_M/sum(NumAtAge_i_j$Survival_F_M)
          EstNumFish_Population = (BiomassCarryingCapacityMSY$BiomassMean_LastYear[i]*1000)/Weight_kg_oneFish
          NumAtAge_i_j$Population_Extrapolated = NumAtAge_i_j$Survival_F_M_Rescaled*EstNumFish_Population
          NumAtAge_i_j$Survival_F_M=NULL
          NumAtAge_i_j$Survival_F_M_Rescaled=NULL
          F_start = 1.0
          min.RSS = function(F_start)
          {
            EstimatedCatch_Baranov = BaranovCatch(NumAtAge_i_j$Age,NumAtAge_i_j$Length,NumAtAge_i_j$Population_Extrapolated,LifeHistoryScenarios_p$M,F_start,NumAtAge_i_j$Selex_Empirical)
            return(sum((EstimatedCatch_Baranov - NumAtAge_i_j$ExtrapolatedCatchInNumber)^2))
          }
          Model = optimize(min.RSS,c(0,10),tol=0.001)
          F_estimated = Model$minimum[1]
		  F_estimated = F_estimated
          NumAtAge_i_j$EstimatedCatch_Baranov = BaranovCatch(NumAtAge_i_j$Age,NumAtAge_i_j$Length,NumAtAge_i_j$Population_Extrapolated,LifeHistoryScenarios_p$M,F_estimated,NumAtAge_i_j$Selex_Empirical)
          NumAtAge_i_j$Residuals = NumAtAge_i_j$ExtrapolatedCatchInNumber-NumAtAge_i_j$EstimatedCatch_Baranov
          NumAtAge_i_j$SSB_CurrentAtAge = NumAtAge_i_j$Population_Extrapolated*(NumAtAge_i_j$Weight_kg*0.001)*NumAtAge_i_j$Maturity
          SSB_current = sum(NumAtAge_i_j$SSB_CurrentAtAge)
          SPR_current = SSB_current/SSB_virgin
          NumAgeCurrent = do.call("rbind",replicate(dim(ALK_Age2Len)[1],transpose(data.frame(NumAtAge_i_j$Population_Extrapolated)),simplify = FALSE))
          NumLengthCurrent = NumAgeCurrent*ALK_Age2Len
          NumAtLength_i_j$Population_Extrapolated = rowSums(NumLengthCurrent)
          CatchCurrent = do.call("rbind",replicate(dim(ALK_Age2Len)[1],transpose(data.frame(NumAtAge_i_j$EstimatedCatch_Baranov)),simplify = FALSE))
          CatchCurrent = CatchCurrent*ALK_Age2Len
          NumAtLength_i_j$EstimatedCatch_Baranov = rowSums(CatchCurrent)
          FolderNameSpp = gsub(" ","_",BiomassCarryingCapacityMSY$species[i])
          png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("CatchAtLength",FolderNameSpp,BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
          plot(NumAtLength_i_j$Midpoint,NumAtLength_i_j$ExtrapolatedCatch,pch=16,xlab="Length (cm)",ylab="Catch in Number",cex.axis=1.1,cex.lab=1.1,ylim=c(0,max(NumAtLength_i_j$EstimatedCatch_Baranov,NumAtLength_i_j$ExtrapolatedCatch)))
          points(NumAtLength_i_j$Midpoint,NumAtLength_i_j$EstimatedCatch_Baranov,type="l",lwd=2,col="red")
          title(BiomassCarryingCapacityMSY$species[i])
          mtext(paste(BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],populationReconstructionMethods[v],sep=", "))
          dev.off()
          png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("CatchAtAge",FolderNameSpp,BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
          plot(NumAtAge_i_j$Age,NumAtAge_i_j$ExtrapolatedCatchInNumber,pch=16,xlab="Age",ylab="Catch in Number",cex.axis=1.1,cex.lab=1.1,ylim=c(0,max(NumAtAge_i_j$EstimatedCatch_Baranov,NumAtAge_i_j$ExtrapolatedCatchInNumber)))
          points(NumAtAge_i_j$Age,NumAtAge_i_j$EstimatedCatch_Baranov,type="l",lwd=2,col="red")
          title(BiomassCarryingCapacityMSY$species[i])
          mtext(paste(BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],populationReconstructionMethods[v],sep=", "))
          dev.off()
		  SPR_lengthBased = computeLengthSPR_selex(LifeHistoryScenarios_p$M,LifeHistoryScenarios_p$Linf,LifeHistoryScenarios_p$VBK,LifeHistoryScenarios_p$var_a,LifeHistoryScenarios_p$var_b,LifeHistoryScenarios_p$Lmat,F_estimated,Selex_Logistic_a_age,Selex_Logistic_b_age)					
        }
        if(populationReconstructionMethods[v]=="BevertonHoltInstantaneousMortality")
        {
          NumAtLength_i_j = NumAtLength_v
          NumAtAge_i_j = NumAtAge_v
          print("Reconstructing population size using Beverton-Holt Instantaneous Mortality...")
          flush.console()
   		   HistForBevHolt = hist(Sub$cm,breaks=seq(0,(ceiling(max(Sub$cm,na.rm=TRUE)/histBreaks)*histBreaks),by=histBreaks),plot=FALSE)
		   HistDataBreaks = data.frame(Length=HistForBevHolt$mids,Counts=HistForBevHolt$counts)
		   HistDataBreaks = subset(HistDataBreaks,HistDataBreaks$Counts==max(HistDataBreaks$Counts))
		   BevHoltEqOut = bheq(len=Sub$cm,Linf=LifeHistoryScenarios_p$Linf,K=LifeHistoryScenarios_p$VBK,Lc=min(Sub$cm),La=LifeHistoryScenarios_p$Linf,nboot = 200)
		   Z = BevHoltEqOut$z[2]
		   F_estimated = Z - LifeHistoryScenarios_p$M
		   F_estimated = F_estimated
		   NumAtAge_i_j$Survival_F_M = exp(-1*(LifeHistoryScenarios_p$M+(F_estimated*NumAtAge_i_j$Selex_Empirical))*NumAtAge_i_j$Age)
          NumAtAge_i_j$Survival_F_M_Rescaled = NumAtAge_i_j$Survival_F_M/sum(NumAtAge_i_j$Survival_F_M)
          EstNumFish_Population = (BiomassCarryingCapacityMSY$BiomassMean_LastYear[i]*1000)/Weight_kg_oneFish
          NumAtAge_i_j$Population_Extrapolated = NumAtAge_i_j$Survival_F_M_Rescaled*EstNumFish_Population
          NumAtAge_i_j$SSB_CurrentAtAge = NumAtAge_i_j$Population_Extrapolated*(NumAtAge_i_j$Weight_kg*0.001)*NumAtAge_i_j$Maturity
          SSB_current = sum(NumAtAge_i_j$SSB_CurrentAtAge)
          SPR_current = SSB_current/SSB_virgin
          NumAtAge_i_j$EstimatedCatch_Baranov = BaranovCatch(NumAtAge_i_j$Age,NumAtAge_i_j$Length,NumAtAge_i_j$Population_Extrapolated,LifeHistoryScenarios_p$M,F_estimated,NumAtAge_i_j$Selex_Empirical)
          NumAtAge_i_j$Residuals = NumAtAge_i_j$ExtrapolatedCatchInNumber-NumAtAge_i_j$EstimatedCatch_Baranov
          NumAtAge_i_j$SSB_CurrentAtAge = NumAtAge_i_j$Population_Extrapolated*(NumAtAge_i_j$Weight_kg*0.001)*NumAtAge_i_j$Maturity
          SSB_current = sum(NumAtAge_i_j$SSB_CurrentAtAge)
          SPR_current = SSB_current/SSB_virgin
          NumAgeCurrent = do.call("rbind",replicate(dim(ALK_Age2Len)[1],transpose(data.frame(NumAtAge_i_j$Population_Extrapolated)),simplify = FALSE))
          NumLengthCurrent = NumAgeCurrent*ALK_Age2Len
          NumAtLength_i_j$Population_Extrapolated = rowSums(NumLengthCurrent)
          CatchCurrent = do.call("rbind",replicate(dim(ALK_Age2Len)[1],transpose(data.frame(NumAtAge_i_j$EstimatedCatch_Baranov)),simplify = FALSE))
          CatchCurrent = CatchCurrent*ALK_Age2Len
          NumAtLength_i_j$EstimatedCatch_Baranov = rowSums(CatchCurrent)
          FolderNameSpp = gsub(" ","_",BiomassCarryingCapacityMSY$species[i])
          png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("CatchAtLength",FolderNameSpp,BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
          plot(NumAtLength_i_j$Midpoint,NumAtLength_i_j$ExtrapolatedCatch,pch=16,xlab="Length (cm)",ylab="Catch in Number",cex.axis=1.1,cex.lab=1.1,ylim=c(0,max(NumAtLength_i_j$EstimatedCatch_Baranov,NumAtLength_i_j$ExtrapolatedCatch)))
          points(NumAtLength_i_j$Midpoint,NumAtLength_i_j$EstimatedCatch_Baranov,type="l",lwd=2,col="red")
          title(BiomassCarryingCapacityMSY$species[i])
          mtext(paste(BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],populationReconstructionMethods[v],sep=", "))
          dev.off()
          png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("CatchAtAge",FolderNameSpp,BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
          plot(NumAtAge_i_j$Age,NumAtAge_i_j$ExtrapolatedCatchInNumber,pch=16,xlab="Age",ylab="Catch in Number",cex.axis=1.1,cex.lab=1.1,ylim=c(0,max(NumAtAge_i_j$EstimatedCatch_Baranov,NumAtAge_i_j$ExtrapolatedCatchInNumber)))
          points(NumAtAge_i_j$Age,NumAtAge_i_j$EstimatedCatch_Baranov,type="l",lwd=2,col="red")
          title(BiomassCarryingCapacityMSY$species[i])
          mtext(paste(BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],populationReconstructionMethods[v],sep=", "))
          dev.off()
		  SPR_lengthBased = computeLengthSPR_selex(LifeHistoryScenarios_p$M,LifeHistoryScenarios_p$Linf,LifeHistoryScenarios_p$VBK,LifeHistoryScenarios_p$var_a,LifeHistoryScenarios_p$var_b,LifeHistoryScenarios_p$Lmat,F_estimated,Selex_Logistic_a_age,Selex_Logistic_b_age)					
        }
        if(populationReconstructionMethods[v]=="LIME")  #The condition that only the mean is used for LIME is becasue LIME estiamtes its own biomass levels (unlike the other 4 methods which use biomass from Catch-MSY). As a result, because we have upper and lower CatchMSY bounds as well, if we don't include this condition, we will accidentially run LIME three times and get the same results, whcih will go into the averaging and bias the averaged values.
        {
          NumAtLength_i_j = NumAtLength_v
          NumAtAge_i_j = NumAtAge_v
          print("Reconstructing population size using LIME...")
          flush.console()
          Len50_selex = max(Sub$cm,na.rm=TRUE)/2
          min.RSS = function(Len50_selex)
          {
            Sel = SelectivityLogistic_Length(Len50_selex,Selex_Logistic_a_length,Selex_Logistic_b_length)
            return(sum((Sel - 0.5)^2))
          }
          Model = optimize(min.RSS,c(0,max(Sub$cm)),tol=0.001)
          Len50_selex = Model$minimum[1]
          Len95_selex = max(Sub$cm,na.rm=TRUE)/2
          min.RSS = function(Len95_selex)
          {
            Sel = SelectivityLogistic_Length(Len95_selex,Selex_Logistic_a_length,Selex_Logistic_b_length)
            return(sum((Sel - 0.95)^2))
          }
          Model = optimize(min.RSS,c(0,max(Sub$cm)),tol=0.001)
          Len95_selex = Model$minimum[1]
		  #NOTE: the steepness h value is NOT used in the estimation of F, but needs to be included as an input otherwise the model will not work. 
          lh = create_lh_list(vbk=LifeHistoryScenarios_p$VBK,linf=LifeHistoryScenarios_p$Linf,lwa=(LifeHistoryScenarios_p$var_a*1000),lwb=LifeHistoryScenarios_p$var_b,
              S50=Len50_selex,S95=Len95_selex,selex_input="length",selex_type=c("logistic"), #Starting selectivity values and function specified
              M50=LifeHistoryScenarios_p$Lmat,maturity_input="length",M=LifeHistoryScenarios_p$M,binwidth=histBreaks,R0=R0, #R0 is starting value - will be found by LIME
              rho=0,h=0.8,  #steepness (h) and cruitment autocorrelation (rho)
              nseasons=1,nfleets=1,CVlen=0.2,SigmaR=0.1,AgeMax=max(NumAtAge_i_j$Age)) #if we could get MONTHLY effort, and split catches per fleet that would go here
          LF=matrix(NumAtLength_i_j$NumSampled,nrow=1)
          colnames(LF)=NumAtLength_i_j$Length_max
          rownames(LF)=max(SpeciesCatch_mT_Scaled$Year)
          Catch_matrix = as.matrix(BiomassCarryingCapacityMSY$catch_t[i]*1000*1000)   #Convert from metric tons to kilograms to grams - LIME wants in grams
          Catch_matrix = Catch_matrix[,length(Catch_matrix)]
          data_all=list(years=max(SpeciesCatch_mT_Scaled$Year),LF=LF,C_ft=as.matrix(Catch_matrix))   #list of years, length distribution (as a matrix of one row), catches and number of observed catches
          inputs_all=create_inputs(lh=lh,input_data=data_all)
          vals_selex_ft = matrix(NumAtLength_i_j$Selex_Fitted,nrow=1)
          Results_LIME = run_LIME(modpath=NULL,input=inputs_all,data_avail="Catch_LC1",C_type=2,est_selex_f=FALSE,vals_selex_ft=vals_selex_ft,derive_quants=TRUE,est_totalF=TRUE,Fpen=0)       #C_type is catch in weight, not numbers
          gradient = Results_LIME$opt$max_gradient <= 0.001
          #hessian = Results_LIME$Sdreport$pdHess    #Comment these out - otherwise, it will stop the code from working if model does not fit. I can filter out these models that don't fit later on. 
          #Converge = (hessian==TRUE & gradient == TRUE)
          #print(paste("gradient should be below 0.001, it is",Results_LIME$opt$max_gradient))
          #print(paste("does the hessian test pass? ",Results_LIME$Sdreport$pdHess))
          #print(paste("SPR estimated as ",Results_LIME$Report$SPR_t,sep=""))
          #flush.console()
		  if(length(Results_LIME$Report)==1)
		  {
			next
		  }
         F_estimated = Results_LIME$Report$F_t
		  NumAtAge_i_j$Survival_F_M = exp(-1*(LifeHistoryScenarios_p$M+(F_estimated*NumAtAge_i_j$Selex_Empirical))*NumAtAge_i_j$Age)
		  NumAtAge_i_j$Survival_F_M_Rescaled = NumAtAge_i_j$Survival_F_M/sum(NumAtAge_i_j$Survival_F_M)
		  EstNumFish_Population = (BiomassCarryingCapacityMSY$BiomassMean_LastYear[i]*1000)/Weight_kg_oneFish
		  NumAtAge_i_j$Population_Extrapolated = NumAtAge_i_j$Survival_F_M_Rescaled*EstNumFish_Population
		  NumAtAge_i_j$SSB_CurrentAtAge = NumAtAge_i_j$Population_Extrapolated*(NumAtAge_i_j$Weight_kg)*NumAtAge_i_j$Maturity
		  SSB_current = sum(NumAtAge_i_j$SSB_CurrentAtAge)
		  SPR_current = SSB_current/SSB_virgin
		  NumAtAge_i_j$EstimatedCatch_Baranov = BaranovCatch(NumAtAge_i_j$Age,NumAtAge_i_j$Length,NumAtAge_i_j$Population_Extrapolated,LifeHistoryScenarios_p$M,F_estimated,NumAtAge_i_j$Selex_Empirical)
		  NumAtAge_i_j$Residuals = NumAtAge_i_j$ExtrapolatedCatchInNumber-NumAtAge_i_j$EstimatedCatch_Baranov
		  NumAtAge_i_j$SSB_CurrentAtAge = NumAtAge_i_j$Population_Extrapolated*(NumAtAge_i_j$Weight_kg)*NumAtAge_i_j$Maturity
		  SSB_current = sum(NumAtAge_i_j$SSB_CurrentAtAge)
		  SPR_current = SSB_current/SSB_virgin
		  NumAgeCurrent = do.call("rbind",replicate(dim(ALK_Age2Len)[1],transpose(data.frame(NumAtAge_i_j$Population_Extrapolated)),simplify = FALSE))
		  NumLengthCurrent = NumAgeCurrent*ALK_Age2Len
		  NumAtLength_i_j$Population_Extrapolated = rowSums(NumLengthCurrent)
		  CatchCurrent = do.call("rbind",replicate(dim(ALK_Age2Len)[1],transpose(data.frame(NumAtAge_i_j$EstimatedCatch_Baranov)),simplify = FALSE))
		  CatchCurrent = CatchCurrent*ALK_Age2Len
		  NumAtLength_i_j$EstimatedCatch_Baranov = rowSums(CatchCurrent)		  
	      FolderNameSpp = gsub(" ","_",BiomassCarryingCapacityMSY$species[i])
          png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("CatchAtLength",FolderNameSpp,BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
          plot(NumAtLength_i_j$Midpoint,NumAtLength_i_j$ExtrapolatedCatch,pch=16,xlab="Length (cm)",ylab="Catch in Number",cex.axis=1.1,cex.lab=1.1,ylim=c(0,max(NumAtLength_i_j$EstimatedCatch_Baranov,NumAtLength_i_j$ExtrapolatedCatch)))
          points(NumAtLength_i_j$Midpoint,NumAtLength_i_j$EstimatedCatch_Baranov,type="l",lwd=2,col="red")
          title(BiomassCarryingCapacityMSY$species[i])
          mtext(paste(BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],populationReconstructionMethods[v],sep=", "))
          dev.off()
          png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("CatchAtAge",FolderNameSpp,BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
          plot(NumAtAge_i_j$Age,NumAtAge_i_j$ExtrapolatedCatchInNumber,pch=16,xlab="Age",ylab="Catch in Number",cex.axis=1.1,cex.lab=1.1,ylim=c(0,max(NumAtAge_i_j$EstimatedCatch_Baranov,NumAtAge_i_j$ExtrapolatedCatchInNumber)))
          points(NumAtAge_i_j$Age,NumAtAge_i_j$EstimatedCatch_Baranov,type="l",lwd=2,col="red")
          title(BiomassCarryingCapacityMSY$species[i])
          mtext(paste(BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],populationReconstructionMethods[v],sep=", "))
          dev.off()
		  SPR_lengthBased = computeLengthSPR_selex(LifeHistoryScenarios_p$M,LifeHistoryScenarios_p$Linf,LifeHistoryScenarios_p$VBK,LifeHistoryScenarios_p$var_a,LifeHistoryScenarios_p$var_b,LifeHistoryScenarios_p$Lmat,F_estimated,Selex_Logistic_a_age,Selex_Logistic_b_age)					
        }		
		if(populationReconstructionMethods[v]=="CatchCurve")
		{
			NumAtLength_i_j = NumAtLength_v
			NumAtAge_i_j = NumAtAge_v
			print("Reconstructing population size using Catch Curve analysis")
			flush.console()
			DataForCatchCurve = data.frame(LogCatch=log(NumAtAge_i_j$ExtrapolatedCatchInNumber[seq(which.max(NumAtAge_i_j$ExtrapolatedCatchInNumber),length(NumAtAge_i_j$ExtrapolatedCatchInNumber),by=1)]),Age=NumAtAge_i_j$Age[seq(which.max(NumAtAge_i_j$ExtrapolatedCatchInNumber),length(NumAtAge_i_j$ExtrapolatedCatchInNumber),by=1)])
			DataForCatchCurve$LogCatch[is.infinite(DataForCatchCurve$LogCatch)]=0
			DataForCatchCurve = subset(DataForCatchCurve,DataForCatchCurve$LogCatch>0)
			if(dim(DataForCatchCurve)[1]<2)
			{
				print(paste("Can not fit Catch Curve for ",BiomassCarryingCapacityMSY$species[i],"...Moving on to next model....",sep=""))		  
				flush.console()
				F_estimated=-999
				SPR_lengthBased=-999
				next
			}
			if(dim(DataForCatchCurve)[1]==2)
			{
				CatchCurveAnalysis = lm(DataForCatchCurve$LogCatch ~ DataForCatchCurve$Age)
				F_estimated = (coef(CatchCurveAnalysis)[2]*-1) - LifeHistoryScenarios_p$M			
			}
			if(dim(DataForCatchCurve)[1]>2)
			{
				DataForCatchCurve$AgeDiff=NA
				for(u in 1:dim(DataForCatchCurve)[1])
				{	
					if(u>1)
					{
						DataForCatchCurve$AgeDiff[u] = DataForCatchCurve$Age[u] - DataForCatchCurve$Age[u-1]
					}
				}
				DataForCatchCurve$AgeDiff[is.na(DataForCatchCurve$AgeDiff)]=1
				DataForCatchCurve = subset(DataForCatchCurve,DataForCatchCurve$AgeDiff==1)
				DataForCatchCurve$AgeDiff=NULL
			}
			if(dim(DataForCatchCurve)[1]<2)
			{
				print(paste("Can not fit Catch Curve for ",BiomassCarryingCapacityMSY$species[i],"...Moving on to next model....",sep=""))		  
				flush.console()
				F_estimated=-999
				SPR_lengthBased=-999
				next
			}
			if(dim(DataForCatchCurve)[1]==2)
			{
				CatchCurveAnalysis = lm(DataForCatchCurve$LogCatch ~ DataForCatchCurve$Age)
				F_estimated = (coef(CatchCurveAnalysis)[2]*-1) - LifeHistoryScenarios_p$M			
			}
			if(dim(DataForCatchCurve)[1]>2)
			{
				DataForCatchCurve$Difference=NA
				for(u in 1:dim(DataForCatchCurve)[1])
				{	
					if(u>1)
					{
						DataForCatchCurve$Difference[u] = DataForCatchCurve$LogCatch[u-1] - DataForCatchCurve$LogCatch[u]
					}
				}
				DataForCatchCurve = subset(DataForCatchCurve,DataForCatchCurve$Difference>0)
				DataForCatchCurve$Difference=NULL
			}
			if(dim(DataForCatchCurve)[1]<2)
			{
				print(paste("Can not fit Catch Curve for ",BiomassCarryingCapacityMSY$species[i],"...Moving on to next model....",sep=""))	  
				flush.console()
				F_estimated=-999
				SPR_lengthBased=-999
				next
			}
			if(dim(DataForCatchCurve)[1]>=2)
			{
				CatchCurveAnalysis = lm(DataForCatchCurve$LogCatch ~ DataForCatchCurve$Age)
				F_estimated = (coef(CatchCurveAnalysis)[2]*-1) - LifeHistoryScenarios_p$M			
			}
			names(F_estimated)=NULL
			if(F_estimated<0)
			{
				print(paste("Can not fit Catch Curve for ",BiomassCarryingCapacityMSY$species[i],"...Moving on to next model....",sep=""))			
				flush.console()
				F_estimated=-999
				SPR_lengthBased=-999
				next
			}
			F_estimated = F_estimated
			if(F_estimated>=0)
			{
				NumAtAge_i_j$Survival_F_M = exp(-1*(LifeHistoryScenarios_p$M+(F_estimated*NumAtAge_i_j$Selex_Empirical))*NumAtAge_i_j$Age)
				NumAtAge_i_j$Survival_F_M_Rescaled = NumAtAge_i_j$Survival_F_M/sum(NumAtAge_i_j$Survival_F_M)
				EstNumFish_Population = (BiomassCarryingCapacityMSY$BiomassMean_LastYear[i]*1000)/Weight_kg_oneFish
				NumAtAge_i_j$Population_Extrapolated = NumAtAge_i_j$Survival_F_M_Rescaled*EstNumFish_Population
				NumAtAge_i_j$SSB_CurrentAtAge = NumAtAge_i_j$Population_Extrapolated*(NumAtAge_i_j$Weight_kg)*NumAtAge_i_j$Maturity
				SSB_current = sum(NumAtAge_i_j$SSB_CurrentAtAge)
				SPR_current = SSB_current/SSB_virgin
				NumAtAge_i_j$EstimatedCatch_Baranov = BaranovCatch(NumAtAge_i_j$Age,NumAtAge_i_j$Length,NumAtAge_i_j$Population_Extrapolated,LifeHistoryScenarios_p$M,F_estimated,NumAtAge_i_j$Selex_Empirical)
				NumAtAge_i_j$Residuals = NumAtAge_i_j$ExtrapolatedCatchInNumber-NumAtAge_i_j$EstimatedCatch_Baranov
				NumAtAge_i_j$SSB_CurrentAtAge = NumAtAge_i_j$Population_Extrapolated*(NumAtAge_i_j$Weight_kg)*NumAtAge_i_j$Maturity
				SSB_current = sum(NumAtAge_i_j$SSB_CurrentAtAge)
				SPR_current = SSB_current/SSB_virgin
				NumAgeCurrent = do.call("rbind",replicate(dim(ALK_Age2Len)[1],transpose(data.frame(NumAtAge_i_j$Population_Extrapolated)),simplify = FALSE))
				NumLengthCurrent = NumAgeCurrent*ALK_Age2Len
				NumAtLength_i_j$Population_Extrapolated = rowSums(NumLengthCurrent)
				CatchCurrent = do.call("rbind",replicate(dim(ALK_Age2Len)[1],transpose(data.frame(NumAtAge_i_j$EstimatedCatch_Baranov)),simplify = FALSE))
				CatchCurrent = CatchCurrent*ALK_Age2Len
				NumAtLength_i_j$EstimatedCatch_Baranov = rowSums(CatchCurrent)		  
			    FolderNameSpp = gsub(" ","_",BiomassCarryingCapacityMSY$species[i])
			    png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("CatchAtLength",FolderNameSpp,BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
			    plot(NumAtLength_i_j$Midpoint,NumAtLength_i_j$ExtrapolatedCatch,pch=16,xlab="Length (cm)",ylab="Catch in Number",cex.axis=1.1,cex.lab=1.1,ylim=c(0,max(NumAtLength_i_j$EstimatedCatch_Baranov,NumAtLength_i_j$ExtrapolatedCatch)))
			    points(NumAtLength_i_j$Midpoint,NumAtLength_i_j$EstimatedCatch_Baranov,type="l",lwd=2,col="red")
			    title(BiomassCarryingCapacityMSY$species[i])
			    mtext(paste(BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],populationReconstructionMethods[v],sep=", "))
			    dev.off()
			    png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("CatchAtAge",FolderNameSpp,BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
			    plot(NumAtAge_i_j$Age,NumAtAge_i_j$ExtrapolatedCatchInNumber,pch=16,xlab="Age",ylab="Catch in Number",cex.axis=1.1,cex.lab=1.1,ylim=c(0,max(NumAtAge_i_j$EstimatedCatch_Baranov,NumAtAge_i_j$ExtrapolatedCatchInNumber)))
			    points(NumAtAge_i_j$Age,NumAtAge_i_j$EstimatedCatch_Baranov,type="l",lwd=2,col="red")
			    title(BiomassCarryingCapacityMSY$species[i])
			    mtext(paste(BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),names(LifeHistoryScenariosList)[p],populationReconstructionMethods[v],sep=", "))
			    dev.off()
				SPR_lengthBased = computeLengthSPR_selex(LifeHistoryScenarios_p$M,LifeHistoryScenarios_p$Linf,LifeHistoryScenarios_p$VBK,LifeHistoryScenarios_p$var_a,LifeHistoryScenarios_p$var_b,LifeHistoryScenarios_p$Lmat,F_estimated,Selex_Logistic_a_age,Selex_Logistic_b_age)					
			}
		}
		ResultsList_i_j_o_p = list(species=BiomassCarryingCapacityMSY$species[i],WPP=BiomassCarryingCapacityMSY$WPP[i],
			CatchMSY_Scenario=BiomassCarryingCapacityMSY$Scenario[i],LandingsSensitivity=BiomassCarryingCapacityMSY$LandingsSensitivity[i],
			Scenario_New=BiomassCarryingCapacityMSY$Scenario_New[i],LifeHistory_Scenario=names(LifeHistoryScenariosList)[p],
			populationReconstructionMethod=populationReconstructionMethods[v],Lmat=LifeHistoryScenarios_p$Lmat,Linf=LifeHistoryScenarios_p$Linf,
			VBK=LifeHistoryScenarios_p$VBK,a=LifeHistoryScenarios_p$var_a,b=LifeHistoryScenarios_p$var_b,M=LifeHistoryScenarios_p$M,
			F_estimated=F_estimated,SSB_virgin=SSB_virgin,R0=R0,Fecundity_min=StartValues_i_j_o_p$Fecundity_min,Fecundity_max=StartValues_i_j_o_p$Fecundity_max,
			Selex_Logistic_a_length=Selex_Logistic_a_length,Selex_Logistic_b_length=Selex_Logistic_b_length,Selex_Logistic_a_age=Selex_Logistic_a_age,
			Selex_Logistic_b_age=Selex_Logistic_b_age,ALK_Len2Age=ALK_Len2Age,ALK_Age2Len=ALK_Age2Len,NumAtLength=NumAtLength_i_j,NumAtAge=NumAtAge_i_j,
			BoundFromCatchMSY=BiomassCarryingCapacityMSY$Bound[i],SPR_lengthBased=SPR_lengthBased)
        ListOfPopulations[[length(ListOfPopulations)+1]] = ResultsList_i_j_o_p
        if(i%%50==0)
        {
          print(paste("Working on iteration ",i," of ",dim(BiomassCarryingCapacityMSY)[1],"...",sep=""))
          flush.console()
        }
      }
    }
  }
  #saveRDS(ListOfPopulations,file=paste(PATH_output,"ListOfPopulations.RData",sep="/"),ascii=FALSE,version=NULL,compress=TRUE,refhook=NULL)
  saveRDS(ListOfPopulations,file=paste(PATH_output,fName_ListOfPopulations,sep="/"),ascii=FALSE,version=NULL,compress=TRUE,refhook=NULL)



}


if(ranLengthReconstructionSeparateProcessors==TRUE
{
	#If multi-threaded the above code, will need to read in all of the individual ".RData" results and combine them into one list that is saved and read back in below
	ListOfPopulations1 = readRDS(file=paste(PATH_output,"ListOfPopulations1.RData",sep="/"),refhook=NULL)
	ListOfPopulations2 = readRDS(file=paste(PATH_output,"ListOfPopulations2.RData",sep="/"),refhook=NULL)
	ListOfPopulations3 = readRDS(file=paste(PATH_output,"ListOfPopulations3.RData",sep="/"),refhook=NULL)
	ListOfPopulations4 = readRDS(file=paste(PATH_output,"ListOfPopulations4.RData",sep="/"),refhook=NULL)
	ListOfPopulations5 = readRDS(file=paste(PATH_output,"ListOfPopulations5.RData",sep="/"),refhook=NULL)
	ListOfPopulations6 = readRDS(file=paste(PATH_output,"ListOfPopulations6.RData",sep="/"),refhook=NULL)
	ListOfPopulations7 = readRDS(file=paste(PATH_output,"ListOfPopulations7.RData",sep="/"),refhook=NULL)
	ListOfPopulations8 = readRDS(file=paste(PATH_output,"ListOfPopulations8.RData",sep="/"),refhook=NULL)
	#ListOfPopulations9 = readRDS(file=paste(PATH_output,"ListOfPopulations9.RData",sep="/"),refhook=NULL)
	#ListOfPopulations10 = readRDS(file=paste(PATH_output,"ListOfPopulations10.RData",sep="/"),refhook=NULL)
	#ListOfPopulations11 = readRDS(file=paste(PATH_output,"ListOfPopulations11.RData",sep="/"),refhook=NULL)
	#ListOfPopulations12 = readRDS(file=paste(PATH_output,"ListOfPopulations12.RData",sep="/"),refhook=NULL)

	#NOTE: Need to combine all sub-lists into one master list!
	TotalLength = length(ListOfPopulations1)+length(ListOfPopulations2)+length(ListOfPopulations3)+length(ListOfPopulations4)+
		length(ListOfPopulations5)+length(ListOfPopulations6)+length(ListOfPopulations7)+length(ListOfPopulations8)#+length(ListOfPopulations9)+
		#length(ListOfPopulations10)+length(ListOfPopulations11)+length(ListOfPopulations12)
	ListOfPopulations = list(length=TotalLength)
	cumLen = length(ListOfPopulations1)
	ListOfPopulations[1:cumLen] = ListOfPopulations1
	cumLen_prev = cumLen+1
	cumLen = cumLen + length(ListOfPopulations2)
	ListOfPopulations[cumLen_prev:cumLen] = ListOfPopulations2
	cumLen_prev = cumLen+1
	cumLen = cumLen + length(ListOfPopulations3)
	ListOfPopulations[cumLen_prev:cumLen] = ListOfPopulations3
	cumLen_prev = cumLen+1
	cumLen = cumLen + length(ListOfPopulations4)
	ListOfPopulations[cumLen_prev:cumLen] = ListOfPopulations4
	cumLen_prev = cumLen+1
	cumLen = cumLen + length(ListOfPopulations5)
	ListOfPopulations[cumLen_prev:cumLen] = ListOfPopulations5
	cumLen_prev = cumLen+1
	cumLen = cumLen + length(ListOfPopulations6)
	ListOfPopulations[cumLen_prev:cumLen] = ListOfPopulations6
	cumLen_prev = cumLen+1
	cumLen = cumLen + length(ListOfPopulations7)
	ListOfPopulations[cumLen_prev:cumLen] = ListOfPopulations7
	cumLen_prev = cumLen+1
	cumLen = cumLen + length(ListOfPopulations8)
	ListOfPopulations[cumLen_prev:cumLen] = ListOfPopulations8
	cumLen_prev = cumLen+1
	#cumLen = cumLen + length(ListOfPopulations9)
	#ListOfPopulations[cumLen_prev:cumLen] = ListOfPopulations9
	#cumLen_prev = cumLen+1
	#cumLen = cumLen + length(ListOfPopulations10)
	#ListOfPopulations[cumLen_prev:cumLen] = ListOfPopulations10
	#cumLen_prev = cumLen+1
	#cumLen = cumLen + length(ListOfPopulations11)
	#ListOfPopulations[cumLen_prev:cumLen] = ListOfPopulations11
	#cumLen_prev = cumLen+1
	#cumLen = cumLen + length(ListOfPopulations12)
	#ListOfPopulations[cumLen_prev:cumLen] = ListOfPopulations12
	saveRDS(ListOfPopulations,file=paste(PATH_output,"ListOfPopulations.RData",sep="/"),ascii=FALSE,version=NULL,compress=TRUE,refhook=NULL)
}



###################### Create Summary Table of Fitted Values #################################

if(devevlopLengthReconstructionSummmaryTable==TRUE)
{
	ListOfPopulations = readRDS(file=paste(PATH_output,"ListOfPopulations.RData",sep="/"),refhook=NULL)
	SummaryTableFittedValues = data.frame(species="",WPP="",CatchMSY_Scenario="",LifeHistory_Scenario="",populationReconstructionMethod="",
		Lmat=0,Linf=0,VBK=0,a=0,b=0,Mnat=0,SSB_virgin=0,F_estimated=0,R0=0,Selex_Logistic_a_length=0,Selex_Logistic_b_length=0,Selex_Logistic_a_age=0,Selex_Logistic_b_age=0,SPR_lengthBased=SPR_lengthBased,ListIndexNumber=0)
	SummaryTableFittedValues = SummaryTableFittedValues[-1,]
	for(i in 1:length(ListOfPopulations))
	{
	  temp = data.frame(species=ListOfPopulations[[i]]$species,WPP=ListOfPopulations[[i]]$WPP,CatchMSY_Scenario=ListOfPopulations[[i]]$Scenario_New,LifeHistory_Scenario=ListOfPopulations[[i]]$LifeHistory_Scenario,
		populationReconstructionMethod=ListOfPopulations[[i]]$populationReconstructionMethod,Lmat=ListOfPopulations[[i]]$Lmat,Linf=ListOfPopulations[[i]]$Linf,VBK=ListOfPopulations[[i]]$VBK,
		a=ListOfPopulations[[i]]$a,b=ListOfPopulations[[i]]$b,Mnat=ListOfPopulations[[i]]$M,SSB_virgin=ListOfPopulations[[i]]$SSB_virgin, 
		F_estimated=ListOfPopulations[[i]]$F_estimated,R0=ListOfPopulations[[i]]$R0,Selex_Logistic_a_length=ListOfPopulations[[i]]$Selex_Logistic_a_length,Selex_Logistic_b_length=ListOfPopulations[[i]]$Selex_Logistic_b_length,
		Selex_Logistic_a_age=ListOfPopulations[[i]]$Selex_Logistic_a_age,Selex_Logistic_b_age=ListOfPopulations[[i]]$Selex_Logistic_b_age,SPR_lengthBased=ListOfPopulations[[i]]$SPR_lengthBased,ListIndexNumber=i)
	  SummaryTableFittedValues = rbind(SummaryTableFittedValues,temp)
	  if(i%%50==0)
	  {
		print(paste("Working on iteration ",i," of ",length(ListOfPopulations),"...",sep=""))
		flush.console()
	  }
	}
	names(BiomassCarryingCapacityMSY)[names(BiomassCarryingCapacityMSY)=="Scenario_New"]="CatchMSY_Scenario"
	names(BiomassCarryingCapacityMSY)[names(BiomassCarryingCapacityMSY)=="Bound"]="CatchMSY_Bound"
	SummaryTableFittedValues = merge(SummaryTableFittedValues,BiomassCarryingCapacityMSY,all.x=TRUE,by=c("species","WPP","CatchMSY_Scenario"))
	write.table(SummaryTableFittedValues,paste(PATH_output,"SummaryTableFittedValues.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
}


##################################### Population Forecasting Simulator #######################
#Setup and Read F-Based Projection Function
projectionFunction = function(Scenario_i,yearsToProject,yearToStartProjection,F_project,h_steepness)
{
  Years = seq(yearToStartProjection,(yearToStartProjection+yearsToProject),by=1)
  NumAtAge = Scenario_i$NumAtAge
  #Project Forward in time
  for(m in 1:length(Years))
  {
    if(m==1)
    {
      Population_end = NumAtAge$Population_Extrapolated
      Population_start = rep(0,length(Population_end))
      for(o in 1:length(Population_end))
      {
        if(o==1)
        {
          Population_start[o] = 0
        }
        if(o>1)
        {
          Population_start[o] = Population_end[o-1]
        }
      }
      SSB_m = sum(Population_end*(NumAtAge$Weight_kg*0.001)*NumAtAge$Maturity)
      Recruitment_m = BevHoltRecruitment(Scenario_i$R0,h_steepness,SSB_m,Scenario_i$SSB_virgin)
      Population_start[1] = Recruitment_m 
      SelexToApply = c(0,NumAtAge$Selex_Empirical[-length(NumAtAge$Selex_Empirical)])
      Population_end = Population_start*exp(-1*(Scenario_i$M+(F_project*SelexToApply)))
      TotalDeaths_m = Population_start - Population_end
      Catch_m = ((F_project*SelexToApply)/(Scenario_i$M+(F_project*SelexToApply)))*Population_start*(1-exp(-(Scenario_i$M+(F_project*SelexToApply))))
      Catch_m_kg = Catch_m*Scenario_i$NumAtAge$Weight_kg
      NaturalDeaths_m = TotalDeaths_m - Catch_m
      Population_end[Population_end<0]=0   #can't have negative numbers of fish
      Population = data.frame(Population_end)
      names(Population) = as.character(Years[m])
      row.names(Population) = as.character(NumAtAge$Age)
      Catch = data.frame(Catch_m)
      names(Catch) = as.character(Years[m])
      row.names(Catch) = as.character(NumAtAge$Age)
      Catch_kg = data.frame(Catch_m_kg)
      names(Catch_kg) = as.character(Years[m])
      row.names(Catch_kg) = as.character(NumAtAge$Age)
      NaturalDeaths = data.frame(NaturalDeaths_m)
      names(NaturalDeaths) = as.character(Years[m])
      row.names(NaturalDeaths) = as.character(NumAtAge$Age)
      SSB = c(SSB_m)
      names(SSB) = as.character(Years[m])
      Recruitment = c(Recruitment_m)
      names(Recruitment) = as.character(Years[m])
      NumberOfFish = c(sum(Population_end))
      names(NumberOfFish) = as.character(Years[m])
      Biomass = c(sum(Population_end*NumAtAge$Weight_kg))
      names(Biomass) = as.character(Years[m])
      Population_sum = sum(Population_end)
      names(Population_sum) = as.character(Years[m])
      Population_Year = Population_sum
      Catch_Year = sum(Catch_m)
      names(Catch_Year) = as.character(Years[m])
      Catch_Year_kg = sum(Catch_m_kg)
      names(Catch_Year_kg) = as.character(Years[m])
      NaturalDeathsYear = sum(NaturalDeaths_m)
      names(NaturalDeathsYear) = as.character(Years[m])
    }
    if(m>1)
    {
      for(o in 1:length(Population_end))
      {
        if(o==1)
        {
          Population_start[o] = 0
        }
        if(o>1)
        {
          Population_start[o] = Population_end[o-1]
        }
      }
      SSB_m = sum(Population_end*(NumAtAge$Weight_kg*0.001)*NumAtAge$Maturity)
      Recruitment_m = BevHoltRecruitment(Scenario_i$R0,h_steepness,SSB_m,Scenario_i$SSB_virgin)
	  Population_start[1] = Recruitment_m 
      SelexToApply = c(0,NumAtAge$Selex_Empirical[-length(NumAtAge$Selex_Empirical)])
      Population_end = Population_start*exp(-1*(Scenario_i$M+(F_project*SelexToApply)))
      TotalDeaths_m = Population_start - Population_end
      Catch_m = ((F_project*SelexToApply)/(Scenario_i$M+(F_project*SelexToApply)))*Population_start*(1-exp(-(Scenario_i$M+(F_project*SelexToApply))))
      Catch_m_kg = Catch_m*Scenario_i$NumAtAge$Weight_kg
      NaturalDeaths_m = TotalDeaths_m - Catch_m
      Population_end[Population_end<0]=0   #can't have negative numbers of fish
      Population_temp = data.frame(Population_end)
      names(Population_temp) = as.character(Years[m])
      row.names(Population_temp) = as.character(NumAtAge$Age)
      Population = cbind(Population,Population_temp)
      Catch_temp = data.frame(Catch_m)
      names(Catch_temp) = as.character(Years[m])
      row.names(Catch_temp) = as.character(NumAtAge$Age)
      Catch = cbind(Catch,Catch_temp)
      Catch_temp_kg = data.frame(Catch_m_kg)
      names(Catch_temp_kg) = as.character(Years[m])
      row.names(Catch_temp_kg) = as.character(NumAtAge$Age)
      Catch_kg = cbind(Catch_kg,Catch_temp_kg)
      NaturalDeaths_temp = data.frame(NaturalDeaths_m)
      names(NaturalDeaths_temp) = as.character(Years[m])
      row.names(NaturalDeaths_temp) = as.character(NumAtAge$Age)
      NaturalDeaths = cbind(NaturalDeaths,NaturalDeaths_temp)
      names(SSB_m) = as.character(Years[m])
      SSB = c(SSB,SSB_m)
      names(Recruitment_m) = as.character(Years[m])
      Recruitment = c(Recruitment,Recruitment_m)
      NumberOfFish_temp = c(sum(Population_end))
      names(NumberOfFish_temp) = as.character(Years[m])
      NumberOfFish = c(NumberOfFish,NumberOfFish_temp)
      Biomass_temp = c(sum(Population_end*NumAtAge$Weight_kg))
      names(Biomass_temp) = as.character(Years[m])
      Biomass = c(Biomass,Biomass_temp)
      Population_sum = sum(Population_end)
      names(Population_sum) = as.character(Years[m])
      Population_Year = c(Population_Year,Population_sum)
      Catch_Year_sum = sum(Catch_temp)
      names(Catch_Year_sum) = as.character(Years[m])
      Catch_Year = c(Catch_Year,Catch_Year_sum)
      Catch_Year_sum_kg = sum(Catch_m_kg)
      names(Catch_Year_sum_kg) = as.character(Years[m])
      Catch_Year_kg = c(Catch_Year_kg,Catch_Year_sum_kg)  
      NaturalDeathsYear_temp = sum(NaturalDeaths_temp)
      names(NaturalDeathsYear_temp) = as.character(Years[m])
      NaturalDeathsYear = c(NaturalDeathsYear,NaturalDeathsYear_temp)
    }
  }
  Scenario_i$F_project = F_project
  Scenario_i$ProjectedPopulationAtAge = Population
  Scenario_i$ProjectedCatchAtAge = Catch
  Scenario_i$ProjectedCatchAtAge_kg = Catch_kg
  Scenario_i$ProjectedNaturalDeathsAtAge = NaturalDeaths
  Scenario_i$SSB = SSB
  Scenario_i$Recruitment = Recruitment
  Scenario_i$ProjectedPopulationYear = Population_Year
  Scenario_i$ProjectedBiomassYear = Biomass
  Scenario_i$ProjectedCatchYear = Catch_Year
  Scenario_i$ProjectedCatchYear_kg = Catch_Year_kg
  Scenario_i$ProjectedNaturalDeathsYear = NaturalDeathsYear
  return(Scenario_i)
}

FindSprCorrespondingF = function(correspondingSPR,Scenario_i,Steepness_Scenario,projectYrs)
{	
	F_project=0
	#Do the initial condition, with only natural mortality and F=0 and SPR 100%. Then incrementally increase F until you get the population to the desired SPR level. 
	Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenario)
	SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin	
	while(SPR_lastYear>correspondingSPR)
	{
		F_project = F_project + 0.1
		Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenario)
		SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
		if(F_project>10)
		{
			F_project = -999
			break
		}
	}
	if(F_project != -999)
	{
		F_project = F_project - 0.1
		if(F_project<0) {F_project=0}
		Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenario)
		SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
		while(SPR_lastYear>correspondingSPR)
		{
			F_project = F_project + 0.01
			Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenario)
			SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
		}
		F_project = F_project - 0.01
		if(F_project<0) {F_project=0}
		Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenario)
		SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
		while(SPR_lastYear>correspondingSPR)
		{
			F_project = F_project + 0.001
			Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenario)
			SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
		}
		F_project = F_project - 0.001
		if(F_project<0) {F_project=0}
		Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenario)
		SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
		while(SPR_lastYear>correspondingSPR)
		{
			F_project = F_project + 0.0001
			Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenario)
			SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
		}
	}
	return(F_project)
}

FindSprCorrespondingF_inFixedNumOfYears = function(correspondingSPR,Scenario_i,Steepness_Scenario,fixedYears)
{	
	F_project=0
	#Do the initial condition, with only natural mortality and F=0 and SPR 100%. Then incrementally increase F until you get the population to the desired SPR level. 
	Projection = projectionFunction(Scenario_i,50,yearToStartProjection,F_project,Steepness_Scenario)
	SPR_lastYear = Projection$SSB[fixedYears]/Projection$SSB_virgin	
	while(SPR_lastYear>correspondingSPR)
	{
		F_project = F_project + 0.1
		Projection = projectionFunction(Scenario_i,50,yearToStartProjection,F_project,Steepness_Scenario)
		SPR_lastYear = Projection$SSB[fixedYears]/Projection$SSB_virgin
		if(F_project>10)
		{
			F_project = -999
			break
		}
	}
	if(F_project != -999)
	{	
		F_project = F_project - 0.1
		if(F_project<0) {F_project=0}
		Projection = projectionFunction(Scenario_i,50,yearToStartProjection,F_project,Steepness_Scenario)
		SPR_lastYear = Projection$SSB[fixedYears]/Projection$SSB_virgin
		while(SPR_lastYear>correspondingSPR)
		{
			F_project = F_project + 0.01
			Projection = projectionFunction(Scenario_i,50,yearToStartProjection,F_project,Steepness_Scenario)
			SPR_lastYear = Projection$SSB[fixedYears]/Projection$SSB_virgin
		}
		F_project = F_project - 0.01
		if(F_project<0) {F_project=0}
		Projection = projectionFunction(Scenario_i,50,yearToStartProjection,F_project,Steepness_Scenario)
		SPR_lastYear = Projection$SSB[fixedYears]/Projection$SSB_virgin
		while(SPR_lastYear>correspondingSPR)
		{
			F_project = F_project + 0.001
			Projection = projectionFunction(Scenario_i,50,yearToStartProjection,F_project,Steepness_Scenario)
			SPR_lastYear = Projection$SSB[fixedYears]/Projection$SSB_virgin
		}
		F_project = F_project - 0.001
		if(F_project<0) {F_project=0}
		Projection = projectionFunction(Scenario_i,50,yearToStartProjection,F_project,Steepness_Scenario)
		SPR_lastYear = Projection$SSB[fixedYears]/Projection$SSB_virgin
		while(SPR_lastYear>correspondingSPR)
		{
			F_project = F_project + 0.0001
			Projection = projectionFunction(Scenario_i,50,yearToStartProjection,F_project,Steepness_Scenario)
			SPR_lastYear = Projection$SSB[fixedYears]/Projection$SSB_virgin
		}
	}
	return(F_project)
}

findFmsy = function(Scenario_i,Steepness_Scenario)
{
	F_project=0    #Initial condition 
	Projection = projectionFunction(Scenario_i,100,yearToStartProjection,F_project,Steepness_Scenarios[j])
	Catch_kg_i = Projection$ProjectedCatchYear_kg[length(Projection$ProjectedCatchYear_kg)]
	Catch_kg_last = -1
	while(Catch_kg_i>Catch_kg_last)
	{		
		Catch_kg_last = Catch_kg_i
		F_project = F_project + 0.1
		Projection = projectionFunction(Scenario_i,100,yearToStartProjection,F_project,Steepness_Scenario)
		Catch_kg_i = Projection$ProjectedCatchYear_kg[length(Projection$ProjectedCatchYear_kg)]
		if(F_project>10)
		{
			F_project = -999
			break
		}
	}
	if(F_project != -999)
	{	
		F_project = F_project - 0.1
		Catch_kg_last = Catch_kg_i
		Projection = projectionFunction(Scenario_i,100,yearToStartProjection,F_project,Steepness_Scenario)
		Catch_kg_i = Projection$ProjectedCatchYear_kg[length(Projection$ProjectedCatchYear_kg)]
		while(Catch_kg_i>Catch_kg_last)
		{		
			Catch_kg_last = Catch_kg_i
			F_project = F_project + 0.01
			Projection = projectionFunction(Scenario_i,100,yearToStartProjection,F_project,Steepness_Scenario)
			Catch_kg_i = Projection$ProjectedCatchYear_kg[length(Projection$ProjectedCatchYear_kg)]
		}
		F_project = F_project - 0.01
		Catch_kg_last = Catch_kg_i
		Projection = projectionFunction(Scenario_i,100,yearToStartProjection,F_project,Steepness_Scenario)
		Catch_kg_i = Projection$ProjectedCatchYear_kg[length(Projection$ProjectedCatchYear_kg)]
		while(Catch_kg_i>Catch_kg_last)
		{		
			Catch_kg_last = Catch_kg_i
			F_project = F_project + 0.001
			Projection = projectionFunction(Scenario_i,100,yearToStartProjection,F_project,Steepness_Scenario)
			Catch_kg_i = Projection$ProjectedCatchYear_kg[length(Projection$ProjectedCatchYear_kg)]
		}
	}
	return(F_project)
}

findFmey = function(Scenario_i,Steepness_Scenario,Cbar,Fmsy,MSY_kg,PricePer_spp)
{
	F_project=0    #Initial condition 
	Projection = projectionFunction(Scenario_i,100,yearToStartProjection,F_project,Steepness_Scenario)
	Fyear_atMSY = suppressWarnings(log(1-(Projection$ProjectedCatchYear_kg/Projection$ProjectedBiomassYear)))*-1			
	F_over_Fmsy_atMSY = Fyear_atMSY[length(Fyear_atMSY)]/Fmsy
	B_over_Bmsy_atMSY = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)]/Bmsy_kg
	Revenue_atMSY = PricePer_spp * F_over_Fmsy_atMSY * B_over_Bmsy_atMSY * MSY_kg
	Cost_atMSY = Cbar * (Fmsy * F_over_Fmsy_atMSY)^beta
	Profit_atMSY_i = Revenue_atMSY - Cost_atMSY
	Profit_atMSY_last = -1
	while(Profit_atMSY_i>Profit_atMSY_last)
	{		
		Profit_atMSY_last = Profit_atMSY_i
		F_project = F_project + 0.1
		Projection = projectionFunction(Scenario_i,100,yearToStartProjection,F_project,Steepness_Scenario)
		Fyear_atMSY = suppressWarnings(log(1-(Projection$ProjectedCatchYear_kg/Projection$ProjectedBiomassYear)))*-1			
		if(is.na(Fyear_atMSY[length(Fyear_atMSY)]))
		{
			break
		}
		F_over_Fmsy_atMSY = Fyear_atMSY[length(Fyear_atMSY)]/Fmsy
		B_over_Bmsy_atMSY = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)]/Bmsy_kg
		Revenue_atMSY = PricePer_spp * F_over_Fmsy_atMSY * B_over_Bmsy_atMSY * MSY_kg
		Cost_atMSY = Cbar * (Fmsy * F_over_Fmsy_atMSY)^beta
		Profit_atMSY_i = Revenue_atMSY - Cost_atMSY
		if(F_project>10)
		{
			F_project = -999
			break
		}
	}
	if(F_project != -999)
	{
		F_project = F_project - 0.1
		Profit_atMSY_last = Profit_atMSY_i
		Projection = projectionFunction(Scenario_i,100,yearToStartProjection,F_project,Steepness_Scenario)
		Fyear_atMSY = suppressWarnings(log(1-(Projection$ProjectedCatchYear_kg/Projection$ProjectedBiomassYear)))*-1			
		F_over_Fmsy_atMSY = Fyear_atMSY[length(Fyear_atMSY)]/Fmsy
		B_over_Bmsy_atMSY = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)]/Bmsy_kg
		Revenue_atMSY = PricePer_spp * F_over_Fmsy_atMSY * B_over_Bmsy_atMSY * MSY_kg
		Cost_atMSY = Cbar * (Fmsy * F_over_Fmsy_atMSY)^beta
		Profit_atMSY_i = Revenue_atMSY - Cost_atMSY
		while(Profit_atMSY_i>Profit_atMSY_last)
		{		
			Profit_atMSY_last = Profit_atMSY_i
			F_project = F_project + 0.01
			Projection = projectionFunction(Scenario_i,100,yearToStartProjection,F_project,Steepness_Scenario)
			Fyear_atMSY = suppressWarnings(log(1-(Projection$ProjectedCatchYear_kg/Projection$ProjectedBiomassYear)))*-1			
			if(is.na(Fyear_atMSY[length(Fyear_atMSY)]))
			{
				break
			}
			F_over_Fmsy_atMSY = Fyear_atMSY[length(Fyear_atMSY)]/Fmsy
			B_over_Bmsy_atMSY = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)]/Bmsy_kg
			Revenue_atMSY = PricePer_spp * F_over_Fmsy_atMSY * B_over_Bmsy_atMSY * MSY_kg
			Cost_atMSY = Cbar * (Fmsy * F_over_Fmsy_atMSY)^beta
			Profit_atMSY_i = Revenue_atMSY - Cost_atMSY
		}
		F_project = F_project - 0.01
		Profit_atMSY_last = Profit_atMSY_i
		Projection = projectionFunction(Scenario_i,100,yearToStartProjection,F_project,Steepness_Scenario)
		Fyear_atMSY = suppressWarnings(log(1-(Projection$ProjectedCatchYear_kg/Projection$ProjectedBiomassYear)))*-1			
		F_over_Fmsy_atMSY = Fyear_atMSY[length(Fyear_atMSY)]/Fmsy
		B_over_Bmsy_atMSY = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)]/Bmsy_kg
		Revenue_atMSY = PricePer_spp * F_over_Fmsy_atMSY * B_over_Bmsy_atMSY * MSY_kg
		Cost_atMSY = Cbar * (Fmsy * F_over_Fmsy_atMSY)^beta
		Profit_atMSY_i = Revenue_atMSY - Cost_atMSY
		while(Profit_atMSY_i>Profit_atMSY_last)
		{		
			Profit_atMSY_last = Profit_atMSY_i
			F_project = F_project + 0.001
			Projection = projectionFunction(Scenario_i,100,yearToStartProjection,F_project,Steepness_Scenario)
			Fyear_atMSY = suppressWarnings(log(1-(Projection$ProjectedCatchYear_kg/Projection$ProjectedBiomassYear)))*-1			
			if(is.na(Fyear_atMSY[length(Fyear_atMSY)]))
			{
				F_project = F_project-0.001
				break
			}
			F_over_Fmsy_atMSY = Fyear_atMSY[length(Fyear_atMSY)]/Fmsy
			B_over_Bmsy_atMSY = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)]/Bmsy_kg
			Revenue_atMSY = PricePer_spp * F_over_Fmsy_atMSY * B_over_Bmsy_atMSY * MSY_kg
			Cost_atMSY = Cbar * (Fmsy * F_over_Fmsy_atMSY)^beta
			Profit_atMSY_i = Revenue_atMSY - Cost_atMSY
		}
	}
	return(F_project)
}

findF_openAccess = function(Scenario_i,Steepness_Scenario,Bmsy)
{
	F_project=0    #Initial condition 
	Projection = projectionFunction(Scenario_i,100,yearToStartProjection,F_project,Steepness_Scenarios[j])
	B_kg_i = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)]
	B_OA_ratio_i = B_kg_i/Bmsy
	Catch_kg_i = 0 
	while(B_OA_ratio_i>0.3)
	{		
		Catch_kg_last = Catch_kg_i
		F_project = F_project + 0.1
		Projection = projectionFunction(Scenario_i,100,yearToStartProjection,F_project,Steepness_Scenario)
		B_kg_i = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)]
		B_OA_ratio_i = B_kg_i/Bmsy
		if(F_project>10)
		{
			F_project = -999
			break
		}
	}
	if(F_project != -999)
	{
		F_project = F_project - 0.1
		Projection = projectionFunction(Scenario_i,100,yearToStartProjection,F_project,Steepness_Scenarios[j])
		B_kg_i = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)]
		B_OA_ratio_i = B_kg_i/Bmsy
		while(B_OA_ratio_i>0.3)
		{		
			Catch_kg_last = Catch_kg_i
			F_project = F_project + 0.01
			Projection = projectionFunction(Scenario_i,100,yearToStartProjection,F_project,Steepness_Scenario)
			B_kg_i = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)]
			B_OA_ratio_i = B_kg_i/Bmsy
		}
		F_project = F_project - 0.01
		Projection = projectionFunction(Scenario_i,100,yearToStartProjection,F_project,Steepness_Scenarios[j])
		B_kg_i = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)]
		B_OA_ratio_i = B_kg_i/Bmsy
		while(B_OA_ratio_i>0.3)
		{		
			Catch_kg_last = Catch_kg_i
			F_project = F_project + 0.001
			Projection = projectionFunction(Scenario_i,100,yearToStartProjection,F_project,Steepness_Scenario)
			B_kg_i = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)]
			B_OA_ratio_i = B_kg_i/Bmsy
		}
	}
	return(F_project)
}



################################ TEMP CODE ############################# 
#Figure out which species were not projected and forecast those
#SummaryTableFittedValues = read.table(paste(PATH_output,"SummaryTableFittedValues.csv",sep="/"),header=TRUE,sep=",")
#SpeciesToProject = sort(unique(SummaryTableFittedValues$species))
#fileNamesProjectionResults_list = list.files(path=PATH_saveProjectionLists,pattern=".RData")
#speciesAssessed = as.data.frame.matrix(t(as.data.frame(strsplit(fileNamesProjectionResults_list,"ListOfPopulations_WithProjections_",fixed=TRUE))))
#speciesAssessed = as.data.frame.matrix(t(as.data.frame(strsplit(speciesAssessed$V2,".",fixed=TRUE))))
#speciesAssessed = speciesAssessed$V1
#speciesAssessed = gsub("_"," ",speciesAssessed)
#SpeciesToProject = data.frame(species=SpeciesToProject)
#speciesAssessed = data.frame(species=speciesAssessed,flag=1)
#SpeciesThatStillNeedForecasting = merge(SpeciesToProject,speciesAssessed,all.x=TRUE,by=c("species"))
#SpeciesThatStillNeedForecasting = subset(SpeciesThatStillNeedForecasting,is.na(SpeciesThatStillNeedForecasting$flag))
#SpeciesThatStillNeedForecasting$flag=NULL
#write.table(SpeciesThatStillNeedForecasting,paste(PATH_output,"SpeciesThatStillNeedForecasting.csv",sep=""),sep=",",col.names=TRUE,row.names=FALSE)

#####################################################################################################


#CHANGE PARAMS STARTING HERE

#Read in price data - NOTE: the price data is from the North Maluku work but should be the same here - NO price data is available in national level data
PricePerKg = read.table(paste(PATH_otherInput,"PricePerKg.csv",sep="/"),header=TRUE,sep=",")

#PROJECTION LOOP
#if(runProjections==TRUE)
{
	ListOfPopulations = readRDS(file=paste(PATH_output,"ListOfPopulations.RData",sep="/"),refhook=NULL)
	SummaryTableFittedValues = read.table(paste(PATH_output,"SummaryTableFittedValues.csv",sep="/"),header=TRUE,sep=",")
	Years = seq(yearToStartProjection,(yearToStartProjection+yearsToProject),by=1)
	SpeciesToProject = sort(unique(SummaryTableFittedValues$species))
	
	############################################## TEMP CODE #######################################################
	
	#SpeciesThatStillNeedForecasting = read.table(paste(PATH_output,"SpeciesThatStillNeedForecasting.csv",sep=""),sep=",",header=TRUE)
	#SpeciesToProject = SpeciesThatStillNeedForecasting$species
	#spp=9
	
	#################################################################################################################

	
	
	#for(spp in 1:length(SpeciesToProject))
	#There are 43 species in total - run about 2 species each on 21 processors/R sessions [last one will have 1 spp]
	#for(spp in 1:2)		#(1)
	#for(spp in 3:4)		#(2)
	#for(spp in 5:6)		#(3)
	#for(spp in 7:8)		#(4)
	#for(spp in 9:10)	#(5)
	#for(spp in 11:12)	#(6)
	#for(spp in 13:14)	#(7)
	#for(spp in 15:16)	#(8)
	#for(spp in 17:18)	#(9)
	#for(spp in 19:20)	#(10)
	#for(spp in 21:22)	#(11)
	#for(spp in 23:24)	#(12)
	
	
	#for(spp in 25:26)	#(13)
	#for(spp in 27:28)	#(14)
	#for(spp in 29:30)	#(15)
	#for(spp in 31:32)	#(16)
	#for(spp in 33:34)	#(17)
	#for(spp in 35:36)	#(18)
	#for(spp in 37:38)	#(19)
	#for(spp in 39:40)	#(20)
	#for(spp in 41:42)	#(21)
	for(spp in 43:43)	#(22)
	{
		ListOfPopulations_WithProjections = list()
		gc()   #garbage collect here to remove the old list so that the computer has enough RAM!!!!!
		SpeciesScenariosToRemove = data.frame(species="",CatchMSY_Scenario="",LifeHistory_Scenario="",populationReconstructionMethod="")
		SpeciesScenariosToRemove=SpeciesScenariosToRemove[-1,]
		ListOfPopulations_WithProjections_fName = paste("ListOfPopulations_WithProjections_",gsub(" ","_",SpeciesToProject[spp]),".RData",sep="")
		SpeciesScenariosToRemove_fName = paste("SpeciesScenariosToRemove_",gsub(" ","_",SpeciesToProject[spp]),".csv",sep="")
		SummaryTableFittedValues_spp = subset(SummaryTableFittedValues,SummaryTableFittedValues$species==SpeciesToProject[spp]) 
	
		#Calculate Price and set beta cost function economic variable
		beta = 1.3
		PricePerKg_sub = PricePerKg[PricePerKg$species==SpeciesToProject[spp],]
		if(dim(PricePerKg_sub)[1]==0)
		{
			PricePer_spp = mean(PricePerKg$price.per.kg,na.rm=TRUE)
		}
		if(dim(PricePerKg_sub)[1]>=1)
		{
			PricePer_spp = mean(PricePerKg_sub$price.per.kg,na.rm=TRUE)
		}
		#Loop through each scenario for each species 		
		for(i in 1:dim(SummaryTableFittedValues_spp)[1])
		{
			#Pull out the population from the list
			Scenario_i = ListOfPopulations[[SummaryTableFittedValues_spp$ListIndexNumber[i]]]
			#If we were not able to estimate current fishing mortality in the assessment/population demographic reconstruction portion of the code.....
			#.....then do that project that scenario and skip to the next one. 
			if(Scenario_i$F_estimated==-999 | sum(Scenario_i$NumAtAge$SSB_virginAtAge,na.rm=TRUE)==0)
			{
				SpeciesScenariosToRemove_i = data.frame(species=Scenario_i$species,CatchMSY_Scenario=Scenario_i$CatchMSY_Scenario,LifeHistory_Scenario=Scenario_i$LifeHistory_Scenario,populationReconstructionMethod=Scenario_i$populationReconstructionMethod)
				SpeciesScenariosToRemove = rbind(SpeciesScenariosToRemove,SpeciesScenariosToRemove_i)
				next
			}	
			NumAtAge = Scenario_i$NumAtAge
			LifeHistParam_i_p = LifeHistoryScenariosList[[Scenario_i$LifeHistory_Scenario]]   #NOTE: Replace Species_i with this value
			LifeHistParam_i_p = subset(LifeHistParam_i_p,LifeHistParam_i_p$species==Scenario_i$species)
			SummaryTableFittedValues_i = subset(SummaryTableFittedValues_spp,SummaryTableFittedValues_spp$WPP==Scenario_i$WPP & SummaryTableFittedValues_spp$CatchMSY_Scenario==Scenario_i$Scenario_New & SummaryTableFittedValues_spp$LifeHistory_Scenario==Scenario_i$LifeHistory_Scenario &
				SummaryTableFittedValues_spp$populationReconstructionMethod==Scenario_i$populationReconstructionMethod)
			#Try three different steepness scenarios: 0.7, 0.8, and 0.9 as specified at the top of the script
			for(j in 1:length(Steepness_Scenarios))
			{
				#Find Reference points and rebuilding times at MSY, open access, MEY, SPR20, SPR30, and spr40
				#MSY
				Fmsy = findFmsy(Scenario_i,Steepness_Scenarios[j])
				if(Fmsy==-999)  {next}
				Projection = projectionFunction(Scenario_i,100,yearToStartProjection,Fmsy,Steepness_Scenarios[j])
				MSY_kg = Projection$ProjectedCatchYear_kg[length(Projection$ProjectedCatchYear_kg)] 
				Bmsy_kg = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)] 
				SPR_at_MSY = Projection$SPR_lengthBased 
				SSB_at_MSY = Projection$SSB[length(Projection$SSB)]
				Recruitment_at_MSY = Projection$Recruitment[length(Projection$Recruitment)]
				MSY_N = Projection$ProjectedPopulationYear[length(Projection$ProjectedPopulationYear)]
				NaturalDeaths_at_MSY = Projection$ProjectedNaturalDeathsYear[length(Projection$ProjectedNaturalDeathsYear)]
				BiomassRelativeToBenchmark_MSY = Projection$ProjectedBiomassYear/max(Projection$ProjectedBiomassYear)
				RebuildTime_MSY = (min(as.numeric(names(which(BiomassRelativeToBenchmark_MSY>=0.95)))) - min(as.numeric(names(BiomassRelativeToBenchmark_MSY))))
				Projection_MSY = Projection
				#Open Access
				F_OA = findF_openAccess(Scenario_i,Steepness_Scenarios[j],Bmsy_kg)
				if(F_OA==-999)  {next}
				Projection = projectionFunction(Scenario_i,100,yearToStartProjection,F_OA,Steepness_Scenarios[j])
				Yield_OA_kg = Projection$ProjectedCatchYear_kg[length(Projection$ProjectedCatchYear_kg)] 
				B_OA_kg = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)] 
				SPR_OA = Projection$SPR_lengthBased 
				SSB_OA = Projection$SSB[length(Projection$SSB)]
				Recruitment_OA = Projection$Recruitment[length(Projection$Recruitment)]
				Yield_N_OA = Projection$ProjectedPopulationYear[length(Projection$ProjectedPopulationYear)]
				NaturalDeaths_OA = Projection$ProjectedNaturalDeathsYear[length(Projection$ProjectedNaturalDeathsYear)]
				Projection_OA = Projection
				#Now that we know where "Open Access" and "MSY" are, calculate the three Upside Economic variables we need
				Fbar = F_OA/Fmsy  
				Bbar = 0.3    
				Cbar = (PricePer_spp*Fbar*Bbar*MSY_kg)/((Fmsy*Fbar)^2)
				#SPR20
				Fspr20 = FindSprCorrespondingF(controlRuleSPR20_percent,Scenario_i,Steepness_Scenarios[j],100)
				if(Fspr20==-999)  {next}
				Projection = projectionFunction(Scenario_i,100,yearToStartProjection,Fspr20,Steepness_Scenarios[j])
				Yield_spr20_kg = Projection$ProjectedCatchYear_kg[length(Projection$ProjectedCatchYear_kg)] 
				B_spr20_kg = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)] 
				SPR_spr20 = Projection$SPR_lengthBased 
				SSB_spr20 = Projection$SSB[length(Projection$SSB)]
				Recruitment_spr20 = Projection$Recruitment[length(Projection$Recruitment)]
				Yield_N_spr20 = Projection$ProjectedPopulationYear[length(Projection$ProjectedPopulationYear)]
				NaturalDeaths_spr20 = Projection$ProjectedNaturalDeathsYear[length(Projection$ProjectedNaturalDeathsYear)]
				BiomassRelativeToBenchmark_SPR20 = Projection$ProjectedBiomassYear/max(Projection$ProjectedBiomassYear)
				RebuildTime_spr20 = (min(as.numeric(names(which(BiomassRelativeToBenchmark_SPR20>=0.95)))) - min(as.numeric(names(BiomassRelativeToBenchmark_SPR20))))
				Projection_SPR20 = Projection
				#SPR30
				Fspr30 = FindSprCorrespondingF(controlRuleSPR30_percent,Scenario_i,Steepness_Scenarios[j],100)
				if(Fspr30==-999)  {next}
				Projection = projectionFunction(Scenario_i,100,yearToStartProjection,Fspr30,Steepness_Scenarios[j])
				Yield_spr30_kg = Projection$ProjectedCatchYear_kg[length(Projection$ProjectedCatchYear_kg)] 
				B_spr30_kg = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)] 
				SPR_spr30 = Projection$SPR_lengthBased 
				SSB_spr30 = Projection$SSB[length(Projection$SSB)]
				Recruitment_spr30 = Projection$Recruitment[length(Projection$Recruitment)]
				Yield_N_spr30 = Projection$ProjectedPopulationYear[length(Projection$ProjectedPopulationYear)]
				NaturalDeaths_spr30 = Projection$ProjectedNaturalDeathsYear[length(Projection$ProjectedNaturalDeathsYear)]
				BiomassRelativeToBenchmark_SPR30 = Projection$ProjectedBiomassYear/max(Projection$ProjectedBiomassYear)
				RebuildTime_spr30 = (min(as.numeric(names(which(BiomassRelativeToBenchmark_SPR30>=0.95)))) - min(as.numeric(names(BiomassRelativeToBenchmark_SPR30))))
				Projection_SPR30 = Projection
				#SPR40
				Fspr40 = FindSprCorrespondingF(controlRuleSPR40_percent,Scenario_i,Steepness_Scenarios[j],100)
				if(Fspr40==-999)  {next}
				Projection = projectionFunction(Scenario_i,100,yearToStartProjection,Fspr40,Steepness_Scenarios[j])
				Yield_spr40_kg = Projection$ProjectedCatchYear_kg[length(Projection$ProjectedCatchYear_kg)] 
				B_spr40_kg = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)] 
				SPR_spr40 = Projection$SPR_lengthBased 
				SSB_spr40 = Projection$SSB[length(Projection$SSB)]
				Recruitment_spr40 = Projection$Recruitment[length(Projection$Recruitment)]
				Yield_N_spr40 = Projection$ProjectedPopulationYear[length(Projection$ProjectedPopulationYear)]
				NaturalDeaths_spr40 = Projection$ProjectedNaturalDeathsYear[length(Projection$ProjectedNaturalDeathsYear)]
				BiomassRelativeToBenchmark_SPR40 = Projection$ProjectedBiomassYear/max(Projection$ProjectedBiomassYear)
				RebuildTime_spr40 = (min(as.numeric(names(which(BiomassRelativeToBenchmark_SPR40>=0.95)))) - min(as.numeric(names(BiomassRelativeToBenchmark_SPR40))))
				Projection_SPR40 = Projection
				#Find other rebuild times relative to other biomass benchmarks:
				#how long it takes to rebuild the population to BatSPR20 when fishing at Fspr30
				BiomassRelativeToBenchmark = Projection_SPR30$ProjectedBiomassYear/B_spr20_kg
				RebuildTime_to_spr30_when_Fspr20 = (min(as.numeric(names(which(BiomassRelativeToBenchmark>=0.95)))) - min(as.numeric(names(BiomassRelativeToBenchmark))))
				#how long it takes to rebuild the population to BatSPR20 when fishing at Fspr40
				BiomassRelativeToBenchmark = Projection_SPR40$ProjectedBiomassYear/B_spr20_kg
				RebuildTime_to_spr40_when_Fspr20 = (min(as.numeric(names(which(BiomassRelativeToBenchmark>=0.95)))) - min(as.numeric(names(BiomassRelativeToBenchmark))))
				#how long it takes to rebuild the population to BatSPR30 when fishing at Fspr40
				BiomassRelativeToBenchmark = Projection_SPR40$ProjectedBiomassYear/B_spr30_kg
				RebuildTime_to_spr40_when_Fspr30 = (min(as.numeric(names(which(BiomassRelativeToBenchmark>=0.95)))) - min(as.numeric(names(BiomassRelativeToBenchmark))))
				#Find F at MEY (Maximum Economic Yield)
				Fmey = findFmey(Scenario_i,Steepness_Scenarios[j],Cbar,Fmsy,MSY_kg,PricePer_spp)
				Projection = projectionFunction(Scenario_i,100,yearToStartProjection,Fmey,Steepness_Scenarios[j])
				Yield_MEY_kg = Projection$ProjectedCatchYear_kg[length(Projection$ProjectedCatchYear_kg)] 
				B_MEY_kg = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)] 
				SPR_MEY = Projection$SPR_lengthBased 
				SSB_MEY = Projection$SSB[length(Projection$SSB)]
				Recruitment_MEY = Projection$Recruitment[length(Projection$Recruitment)]
				Yield_N_MEY = Projection$ProjectedPopulationYear[length(Projection$ProjectedPopulationYear)]
				NaturalDeaths_MEY = Projection$ProjectedNaturalDeathsYear[length(Projection$ProjectedNaturalDeathsYear)]
				BiomassRelativeToBenchmark_MEY = Projection$ProjectedBiomassYear/max(Projection$ProjectedBiomassYear)
				RebuildTime_MEY = (min(as.numeric(names(which(BiomassRelativeToBenchmark_MEY>=0.95)))) - min(as.numeric(names(BiomassRelativeToBenchmark_MEY))))
				Projection_MEY = Projection
				#Will one of the benchmarks not fit or does its estimation hit boundary conditions, if so, then skip this scenario and don't use it because it means the Fcurrent estimate also hit boundary conditions and we will not want to use this scenario anyway
				if(Fmsy==-999 | F_OA==-999 | Fspr20==-999 | Fspr30==-999 | Fspr40==-999 | Fmey==-999)
				{
					next
				}
				#Is Overfished or experiencing overfishing?
				Overfished = FALSE
				Overfishing = FALSE				
				if((SummaryTableFittedValues_i$BiomassMean_LastYear*1000) < B_spr20_kg)
				{
					Overfished = TRUE
				}
				if(Scenario_i$F_estimated > Fspr40)
				{
					Overfishing = TRUE
				}				
				#Find the Fishing mortality for each fishing mortality projection scenario
				for(v in 1:length(FishingMortaltiyProjectionScenarios))
				{
					F_project=NA
					RebuildingTime=NA
					willNotRebuild = NA
					if(FishingMortaltiyProjectionScenarios[v]=="CurrentF")				{F_project = Scenario_i$F_estimated}
					if(FishingMortaltiyProjectionScenarios[v]=="F_equal_0")				{F_project = 0}
					if(FishingMortaltiyProjectionScenarios[v]=="F_at_SPR20")			{F_project = Fspr20}
					if(FishingMortaltiyProjectionScenarios[v]=="F_at_SPR30")			{F_project = Fspr30}
					if(FishingMortaltiyProjectionScenarios[v]=="F_at_SPR40")			{F_project = Fspr40}
					if(FishingMortaltiyProjectionScenarios[v]=="MSY")					{F_project = Fmsy}
					if(FishingMortaltiyProjectionScenarios[v]=="OpenAccess")			{F_project = F_OA}
					if(FishingMortaltiyProjectionScenarios[v]=="MEY")					{F_project = Fmey}
					#Regardless of whether or not it is overfished, do the regular F-based projection with Beverton-Holt Recruitment for each fishing mortality scenario
					if(FishingMortaltiyProjectionScenarios[v]=="CurrentF")				
					{
						Projection = projectionFunction(Scenario_i,yearsToProject,yearToStartProjection,F_project,Steepness_Scenarios[j])
						Projection$FishingMortaltiyProjectionScenario = FishingMortaltiyProjectionScenarios[v]
						Projection$SteepnessAssumed = Steepness_Scenarios[j]
						Projection$ProjectionType = "RegularF"
						Projection$F_project=F_project								
						Yield_kg = Projection$ProjectedCatchYear_kg[length(Projection$ProjectedCatchYear_kg)] 
						B_kg = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)] 
						SPR = Projection$SPR_lengthBased 
						SSB = Projection$SSB[length(Projection$SSB)]
						Recruitment = Projection$Recruitment[length(Projection$Recruitment)]
						Yield_N = Projection$ProjectedPopulationYear[length(Projection$ProjectedPopulationYear)]
						NaturalDeaths = Projection$ProjectedNaturalDeathsYear[length(Projection$ProjectedNaturalDeathsYear)]
						Projection$Overfished = Overfished
						Projection$Overfishing = Overfishing
						Projection$RebuildingTime = NA
						Projection$willNotRebuild=NA
						Projection$WPP = Scenario_i$WPP
						#Calculate Profit, Cost, and Revenue (Upside Model)
						Fbar = F_OA/Fmsy  
						Bbar = 0.3
						Cbar = (PricePer_spp*Fbar*Bbar*MSY_kg)/((Fmsy*Fbar)^2)
						Fyear_atOA = (log(1-(Projection$ProjectedCatchYear_kg/Projection$ProjectedBiomassYear)))*-1			
						F_over_Fmsy = Fyear_atOA/Fmsy
						B_over_Bmsy = Projection$ProjectedBiomassYear/Bmsy_kg
						Revenue = PricePer_spp * F_over_Fmsy * B_over_Bmsy * MSY_kg
						Cost = Cbar * (Fmsy * F_over_Fmsy)^beta
						Profit = Revenue - Cost
						Projection$Revenue = Revenue
						Projection$Cost = Cost
						Projection$Profit = Profit	
						#Attach this other stuff ONLY to this scenario to easily build the benchmark table later on. ONLY attach it to the F_current projections. There will be a secondary loop in the projection summary part that will pull these values into a separate status file. 								
						Projection$Fspr40 = Fspr40
						Projection$Yield_spr40_kg = Yield_spr40_kg								
						Projection$B_spr40_kg = B_spr40_kg
						Projection$SSB_spr40 = SSB_spr40
						Projection$Recruitment_spr40 = Recruitment_spr40
						Projection$Yield_N_spr40 = Yield_N_spr40
						Projection$NaturalDeaths_spr40 = NaturalDeaths_spr40
						Projection$RebuildTime_spr40 = RebuildTime_spr40
						Projection$Fspr30 = Fspr30
						Projection$Yield_spr30_kg = Yield_spr30_kg
						Projection$B_spr30_kg = B_spr30_kg
						Projection$SSB_spr30 = SSB_spr30
						Projection$Recruitment_spr30 = Recruitment_spr30
						Projection$Yield_N_spr30 = Yield_N_spr30
						Projection$NaturalDeaths_spr30 = NaturalDeaths_spr30								
						Projection$RebuildTime_spr30 = RebuildTime_spr30
						Projection$Fspr20 = Fspr20
						Projection$Yield_spr20_kg = Yield_spr20_kg
						Projection$B_spr20_kg = B_spr20_kg
						Projection$SSB_spr20 = SSB_spr20
						Projection$Recruitment_spr20 = Recruitment_spr20
						Projection$Yield_N_spr20 = Yield_N_spr20
						Projection$NaturalDeaths_spr20 = NaturalDeaths_spr20
						Projection$RebuildTime_spr20 = RebuildTime_spr20
						Projection$Fmsy = Fmsy
						Projection$MSY_kg = MSY_kg
						Projection$Bmsy_kg = Bmsy_kg
						Projection$SPR_at_MSY = SPR_at_MSY
						Projection$SSB_at_MSY = SSB_at_MSY
						Projection$Recruitment_at_MSY = Recruitment_at_MSY
						Projection$MSY_N = MSY_N
						Projection$NaturalDeaths_at_MSY = NaturalDeaths_at_MSY
						Projection$RebuildTime_MSY = RebuildTime_MSY
						Projection$F_OA = F_OA
						Projection$Yield_OA_kg = Yield_OA_kg
						Projection$B_OA_kg = B_OA_kg
						Projection$SPR_OA = SPR_OA
						Projection$SSB_OA = SSB_OA
						Projection$Recruitment_OA = Recruitment_OA
						Projection$Yield_N_OA = Yield_N_OA
						Projection$NaturalDeaths_OA = NaturalDeaths_OA
						Projection$Fmey = Fmey
						Projection$Yield_MEY_kg = Yield_MEY_kg
						Projection$B_MEY_kg = B_MEY_kg
						Projection$SPR_MEY = SPR_MEY
						Projection$SSB_MEY = SSB_MEY
						Projection$Recruitment_MEY = Recruitment_MEY
						Projection$Yield_N_MEY = Yield_N_MEY
						Projection$NaturalDeaths_MEY = NaturalDeaths_MEY
						Projection$RebuildTime_MEY = RebuildTime_MEY
						Projection$RebuildTime_to_spr30_when_Fspr20 = RebuildTime_to_spr30_when_Fspr20
						Projection$RebuildTime_to_spr40_when_Fspr20 = RebuildTime_to_spr40_when_Fspr20
						Projection$RebuildTime_to_spr40_when_Fspr30 = RebuildTime_to_spr40_when_Fspr30
						#Now add this object to the list
						ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection
					}
					if(FishingMortaltiyProjectionScenarios[v]!="CurrentF")
					{
						Projection = projectionFunction(Scenario_i,yearsToProject,yearToStartProjection,F_project,Steepness_Scenarios[j])
						Projection$FishingMortaltiyProjectionScenario = FishingMortaltiyProjectionScenarios[v]
						Projection$SteepnessAssumed = Steepness_Scenarios[j]
						Projection$ProjectionType = "RegularF"
						Projection$F_project=F_project								
						Projection$Yield_Equilibrium = Projection$ProjectedCatchYear_kg[length(Projection$ProjectedCatchYear_kg)] 
						Projection$B_kg_Equilibrium = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)] 
						Projection$SPR_Equilibrium = Projection$SPR_lengthBased 
						Projection$SSB_Equilibrium = Projection$SSB[length(Projection$SSB)]
						Projection$Recruitment_Equilibrium = Projection$Recruitment[length(Projection$Recruitment)]
						Projection$Yield_N_Equilibrium = Projection$ProjectedPopulationYear[length(Projection$ProjectedPopulationYear)]
						Projection$NaturalDeaths_Equilibrium = Projection$ProjectedNaturalDeathsYear[length(Projection$ProjectedNaturalDeathsYear)]
						Projection$Overfished = Overfished
						Projection$Overfishing = Overfishing								
						BiomassRelativeToBenchmark = Projection$ProjectedBiomassYear/B_spr40_kg
						if(length(which(BiomassRelativeToBenchmark>=0.95))==0)
						{
							Projection$RebuildingTime = NA
							Projection$willNotRebuild=TRUE								
						}
						if(length(which(BiomassRelativeToBenchmark>=0.95))>0)
						{
							Projection$RebuildingTime = (min(as.numeric(names(which(BiomassRelativeToBenchmark>=0.95)))) - min(as.numeric(names(BiomassRelativeToBenchmark))))
							Projection$willNotRebuild=FALSE
						}
						Projection$WPP = Scenario_i$WPP
						#Calculate Profit, Cost, and Revenue (Upside Model)
						Fbar = F_OA/Fmsy  
						Bbar = 0.3    
						Cbar = (PricePer_spp*Fbar*Bbar*MSY_kg)/((Fmsy*Fbar)^2)
						Fyear_atOA = (log(1-(Projection$ProjectedCatchYear_kg/Projection$ProjectedBiomassYear)))*-1			
						F_over_Fmsy = Fyear_atOA/Fmsy
						B_over_Bmsy = Projection$ProjectedBiomassYear/Bmsy_kg
						Revenue = PricePer_spp * F_over_Fmsy * B_over_Bmsy * MSY_kg
						Cost = Cbar * (Fmsy * F_over_Fmsy)^beta
						Profit = Revenue - Cost
						Projection$Revenue = Revenue
						Projection$Cost = Cost
						Projection$Profit = Profit	
						ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection							
					}
				}	#this closes the bracket for the regular fishing mortality projections 
				#If it is overfished, then you need to run rebuilding scenarios........
				if(Overfished==TRUE)
				{
					#Will it / does it have the potential to rebuild even when F=0?
					ProjectionAtF0 = projectionFunction(Scenario_i,yearsToProject,yearToStartProjection,0,Steepness_Scenarios[j])
					BiomassRelativeToBenchmarkSPR40_atF0 = ProjectionAtF0$ProjectedBiomassYear/B_spr40_kg
					BiomassRelativeToBenchmarkSPR40_atF0_subset = subset(BiomassRelativeToBenchmarkSPR40_atF0,BiomassRelativeToBenchmarkSPR40_atF0>=0.95)
					if(length(BiomassRelativeToBenchmarkSPR40_atF0_subset)==0)
					{
						willNotRebuild=TRUE
					}
					if(length(BiomassRelativeToBenchmarkSPR40_atF0_subset)>0)
					{
						RebuildTimeTotargetSPR40_atF0 = (min(as.numeric(names(which(BiomassRelativeToBenchmarkSPR40_atF0>=0.95)))) - min(as.numeric(names(BiomassRelativeToBenchmarkSPR40_atF0))))+1
						if(is.na(RebuildTimeTotargetSPR40_atF0) | RebuildTimeTotargetSPR40_atF0>20)
						{
							willNotRebuild=TRUE
						} else {
							willNotRebuild=FALSE
						}
					}
					if(willNotRebuild==TRUE)
					{
						Projection = ProjectionAtF0
						Projection$FishingMortaltiyProjectionScenario = "WillNotRebuild"
						Projection$SteepnessAssumed = Steepness_Scenarios[j]
						Projection$ProjectionType = "RegularF"
						Projection$F_project=NA
						Projection$Yield_Equilibrium = NA
						Projection$B_kg_Equilibrium = NA
						Projection$SPR_Equilibrium = NA
						Projection$SSB_Equilibrium = NA
						Projection$Recruitment_Equilibrium = NA
						Projection$Yield_N_Equilibrium = NA
						Projection$NaturalDeaths_Equilibrium = NA
						Projection$Overfished = Overfished
						Projection$Overfishing = Overfishing								
						Projection$RebuildingTime = NA
						Projection$willNotRebuild=TRUE
						Projection$WPP = Scenario_i$WPP
						Projection$Revenue = NA
						Projection$Cost = NA
						Projection$Profit = NA	
						ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection
					}
					if(willNotRebuild==FALSE)
					{
						#Run Three Rebuilding Scenarios for each of the climate scenarios: 
						#(1) Find the F that rebuilds the population in 20 years back to spr 40%
						projectYrs20 = 20
						F_project = FindSprCorrespondingF_inFixedNumOfYears(controlRuleSPR40_percent,Scenario_i,Steepness_Scenarios[j],projectYrs20)
						Projection = projectionFunction(Scenario_i,yearsToProject,yearToStartProjection,F_project,Steepness_Scenarios[j])
						Projection$FishingMortaltiyProjectionScenario = "F_rebuild_20YrsMSC_max"
						Projection$SteepnessAssumed = Steepness_Scenarios[j]
						Projection$ProjectionType = "RegularF"
						Projection$F_project=F_project								
						Projection$Yield_Equilibrium = Projection$ProjectedCatchYear_kg[length(Projection$ProjectedCatchYear_kg)] 
						Projection$B_kg_Equilibrium = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)] 
						Projection$SPR_Equilibrium = Projection$SPR_lengthBased 
						Projection$SSB_Equilibrium = Projection$SSB[length(Projection$SSB)]
						Projection$Recruitment_Equilibrium = Projection$Recruitment[length(Projection$Recruitment)]
						Projection$Yield_N_Equilibrium = Projection$ProjectedPopulationYear[length(Projection$ProjectedPopulationYear)]
						Projection$NaturalDeaths_Equilibrium = Projection$ProjectedNaturalDeathsYear[length(Projection$ProjectedNaturalDeathsYear)]
						Projection$Overfished = Overfished
						Projection$Overfishing = Overfishing								
						BiomassRelativeToBenchmark = Projection$ProjectedBiomassYear/B_spr40_kg
						Projection$RebuildingTime = projectYrs20
						Projection$willNotRebuild=FALSE								
						Projection$WPP = Scenario_i$WPP
						#Calculate Profit, Cost, and Revenue (Upside Model)
						Fbar = F_OA/Fmsy  
						Bbar = 0.3    
						Cbar = (PricePer_spp*Fbar*Bbar*MSY_kg)/((Fmsy*Fbar)^2)
						Fyear_atOA = (log(1-(Projection$ProjectedCatchYear_kg/Projection$ProjectedBiomassYear)))*-1			
						F_over_Fmsy = Fyear_atOA/Fmsy
						B_over_Bmsy = Projection$ProjectedBiomassYear/Bmsy_kg
						Revenue = PricePer_spp * F_over_Fmsy * B_over_Bmsy * MSY_kg
						Cost = Cbar * (Fmsy * F_over_Fmsy)^beta
						Profit = Revenue - Cost
						Projection$Revenue = Revenue
						Projection$Cost = Cost
						Projection$Profit = Profit	
						ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection							
						#(2) Calculate the Tmin plus one generation time number of years and find the F that rebuilds the population in that time period back to spr 40%
						#one generation time calculation
						oneGenTime_numerator = Scenario_i$NumAtAge$Survival_M*Scenario_i$NumAtAge$Age*Scenario_i$NumAtAge$Maturity
						oneGenTime_denominator = Scenario_i$NumAtAge$Survival_M*Scenario_i$NumAtAge$Maturity
						oneGenerationTime = ceiling(sum(oneGenTime_numerator)/sum(oneGenTime_denominator))
						YearCounter=1
						SPR_lastYear=sum(Scenario_i$NumAtAge$Population_Extrapolated*Scenario_i$NumAtAge$Weight_kg*Scenario_i$NumAtAge$Maturity)/sum(Scenario_i$NumAtAge$PopulationVirgin*Scenario_i$NumAtAge$Weight_kg*Scenario_i$NumAtAge$Maturity)
						F_project=0
						while(SPR_lastYear<controlRuleSPR40_percent)
						{
							if(YearCounter==1)
							{
								Population_end = NumAtAge$Population_Extrapolated
								Population_start = rep(0,length(Population_end))
								for(o in 1:length(Population_end))
								{
									if(o==1)
									{
										Population_start[o] = 0
									}
									if(o>1)
									{
										Population_start[o] = Population_end[o-1]
									}
								}
								Catch_m = ((F_project*NumAtAge$Selex_Fitted)/(Scenario_i$M+(F_project*NumAtAge$Selex_Fitted)))*Population_start*(1-exp(-(Scenario_i$M+(F_project*NumAtAge$Selex_Fitted))))
								NaturalDeaths_m = (Scenario_i$M/(Scenario_i$M+(F_project*NumAtAge$Selex_Fitted)))*Population_start*(1-exp(-(Scenario_i$M+(F_project*NumAtAge$Selex_Fitted))))
								Population_end = Population_start - (Catch_m + NaturalDeaths_m)
								Population_end[Population_end<0]=0   #can't have negative numbers of fish
								SSB_m = sum(NumAtAge$SSB_CurrentAtAge)
								Recruitment_m =  BevHoltRecruitment(Scenario_i$R0,Steepness_Scenarios[j],SSB_m,Scenario_i$SSB_virgin)
								SPR_lastYear = SSB_m/Scenario_i$SSB_virgin
								YearCounter = YearCounter + 1
							}
							if(YearCounter>1)
							{
								for(o in 1:length(Population_end))
								{
									if(o==1)
									{
										Population_start[o] = 0
									}
									if(o>1)
									{
										Population_start[o] = Population_end[o-1]
									}
								}
								Catch_m = ((F_project*NumAtAge$Selex_Fitted)/(Scenario_i$M+(F_project*NumAtAge$Selex_Fitted)))*Population_start*(1-exp(-(Scenario_i$M+(F_project*NumAtAge$Selex_Fitted))))
								NaturalDeaths_m = (Scenario_i$M/(Scenario_i$M+(F_project*NumAtAge$Selex_Fitted)))*Population_start*(1-exp(-(Scenario_i$M+(F_project*NumAtAge$Selex_Fitted))))
								Population_end = Population_start - (Catch_m + NaturalDeaths_m)
								Population_end[Population_end<0]=0   #can't have negative numbers of fish
								SSB_m = sum(Population_end*(NumAtAge$Weight_kg*0.001)*NumAtAge$Maturity)
								Recruitment_m =  BevHoltRecruitment(Scenario_i$R0,Steepness_Scenarios[j],SSB_m,Scenario_i$SSB_virgin)
								Population_end[1] = Recruitment_m
								SPR_lastYear = SSB_m/Scenario_i$SSB_virgin
								YearCounter = YearCounter + 1
							}
						}
						timeToRebuld_noFishing = YearCounter
						#Add the time periods together, find the F that builds in that time, and project  
						projectYrs_oneGenTmim = oneGenerationTime + timeToRebuld_noFishing
						F_project = FindSprCorrespondingF_inFixedNumOfYears(controlRuleSPR40_percent,Scenario_i,Steepness_Scenarios[j],projectYrs_oneGenTmim)
						Projection = projectionFunction(Scenario_i,yearsToProject,yearToStartProjection,F_project,Steepness_Scenarios[j])
						Projection$FishingMortaltiyProjectionScenario = "F_rebuild_Tmin_oneGeneration"
						Projection$SteepnessAssumed = Steepness_Scenarios[j]
						Projection$ProjectionType = "RegularF"
						Projection$F_project=F_project								
						Projection$Yield_Equilibrium = Projection$ProjectedCatchYear_kg[length(Projection$ProjectedCatchYear_kg)] 
						Projection$B_kg_Equilibrium = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)] 
						Projection$SPR_Equilibrium = Projection$SPR_lengthBased 
						Projection$SSB_Equilibrium = Projection$SSB[length(Projection$SSB)]
						Projection$Recruitment_Equilibrium = Projection$Recruitment[length(Projection$Recruitment)]
						Projection$Yield_N_Equilibrium = Projection$ProjectedPopulationYear[length(Projection$ProjectedPopulationYear)]
						Projection$NaturalDeaths_Equilibrium = Projection$ProjectedNaturalDeathsYear[length(Projection$ProjectedNaturalDeathsYear)]
						Projection$Overfished = Overfished
						Projection$Overfishing = Overfishing								
						BiomassRelativeToBenchmark = Projection$ProjectedBiomassYear/B_spr40_kg
						Projection$RebuildingTime = projectYrs_oneGenTmim
						Projection$willNotRebuild=FALSE								
						Projection$WPP = Scenario_i$WPP
						#Calculate Profit, Cost, and Revenue (Upside Model)
						Fbar = F_OA/Fmsy  
						Bbar = 0.3    
						Cbar = (PricePer_spp*Fbar*Bbar*MSY_kg)/((Fmsy*Fbar)^2)
						Fyear_atOA = (log(1-(Projection$ProjectedCatchYear_kg/Projection$ProjectedBiomassYear)))*-1			
						F_over_Fmsy = Fyear_atOA/Fmsy
						B_over_Bmsy = Projection$ProjectedBiomassYear/Bmsy_kg
						Revenue = PricePer_spp * F_over_Fmsy * B_over_Bmsy * MSY_kg
						Cost = Cbar * (Fmsy * F_over_Fmsy)^beta
						Profit = Revenue - Cost
						Projection$Revenue = Revenue
						Projection$Cost = Cost
						Projection$Profit = Profit	
						ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection							
						#(3) Find the F that rebuilds the population in 10 years back to spr 40%
						projectYrs10 = 10
						F_project = FindSprCorrespondingF_inFixedNumOfYears(controlRuleSPR40_percent,Scenario_i,Steepness_Scenarios[j],projectYrs10)
						Projection = projectionFunction(Scenario_i,yearsToProject,yearToStartProjection,F_project,Steepness_Scenarios[j])
						Projection$FishingMortaltiyProjectionScenario = "F_rebuild_10YrsUS"
						Projection$SteepnessAssumed = Steepness_Scenarios[j]
						Projection$ProjectionType = "RegularF"
						Projection$F_project=F_project								
						Projection$Yield_Equilibrium = Projection$ProjectedCatchYear_kg[length(Projection$ProjectedCatchYear_kg)] 
						Projection$B_kg_Equilibrium = Projection$ProjectedBiomassYear[length(Projection$ProjectedBiomassYear)] 
						Projection$SPR_Equilibrium = Projection$SPR_lengthBased 
						Projection$SSB_Equilibrium = Projection$SSB[length(Projection$SSB)]
						Projection$Recruitment_Equilibrium = Projection$Recruitment[length(Projection$Recruitment)]
						Projection$Yield_N_Equilibrium = Projection$ProjectedPopulationYear[length(Projection$ProjectedPopulationYear)]
						Projection$NaturalDeaths_Equilibrium = Projection$ProjectedNaturalDeathsYear[length(Projection$ProjectedNaturalDeathsYear)]
						Projection$Overfished = Overfished
						Projection$Overfishing = Overfishing								
						BiomassRelativeToBenchmark = Projection$ProjectedBiomassYear/B_spr40_kg
						Projection$RebuildingTime = projectYrs10
						Projection$willNotRebuild=FALSE								
						Projection$WPP = Scenario_i$WPP
						#Calculate Profit, Cost, and Revenue (Upside Model)
						Fbar = F_OA/Fmsy  
						Bbar = 0.3    
						Cbar = (PricePer_spp*Fbar*Bbar*MSY_kg)/((Fmsy*Fbar)^2)
						Fyear_atOA = (log(1-(Projection$ProjectedCatchYear_kg/Projection$ProjectedBiomassYear)))*-1			
						F_over_Fmsy = Fyear_atOA/Fmsy
						B_over_Bmsy = Projection$ProjectedBiomassYear/Bmsy_kg
						Revenue = PricePer_spp * F_over_Fmsy * B_over_Bmsy * MSY_kg
						Cost = Cbar * (Fmsy * F_over_Fmsy)^beta
						Profit = Revenue - Cost
						Projection$Revenue = Revenue
						Projection$Cost = Cost
						Projection$Profit = Profit	
						ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection							
					}
				}  #is overfished bracket - rebulding scenarios
			}  #this bracket is for steepness scenarios
			#if(i%%50==0)
			{
				print(paste("Projecting population scenario ",i," for species ",SpeciesToProject[spp]," of ",dim(SummaryTableFittedValues_spp)[1],"...",sep=""))
				flush.console()
			}
		}   #this bracket is for the list of scenarios within one species
		saveRDS(ListOfPopulations_WithProjections,file=paste(PATH_saveProjectionLists,ListOfPopulations_WithProjections_fName,sep="/"),ascii=FALSE,version=NULL,compress=TRUE,refhook=NULL)
		write.table(SpeciesScenariosToRemove,paste(PATH_saveProjectionLists,SpeciesScenariosToRemove_fName,sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
	

	}     #this bracket is for the species that is being worked on


}        #this bracket is for the boolean as to whether or not to run projections 




if(processProjectionOutputCreateSummaryFile==TRUE)
{

	SummaryProjectionScenarios_fName = "SummaryProjectionScenario8.csv"
	Benchmarks_fName = "Benchmarks8.csv"

	fileNamesProjectionResults_list = list.files(path=PATH_saveProjectionLists,pattern=".RData")
	fileNamesProjectionResults_NoFit = list.files(path=PATH_saveProjectionLists,pattern=".csv")
	SummaryProjectionScenarios = data.frame(species="",WPP="",CatchMSY_Scenario="",LifeHistory_Scenario="",populationReconstructionMethod="",FishingMortaltiyProjectionScenario="",
		ProjectionType="",ListFileName="",ListIndexNumber=0,
		Lmat=0,Linf=0,VBK=0,a=0,
		b=0,Mnat=0,F_estimated=0,SSB_virgin=0,
		SSB_current=0,R0=0,Yield_estimated=0,N_estimated=0,
		B_estimated=0,Recruits_estimated=0,Selex_Logistic_a_length=0,
		Selex_Logistic_b_length=0,Selex_Logistic_a_age=0,Selex_Logistic_b_age=0,
		SteepnessAssumed=0,F_project=0,NumFishAtEnd=0,BiomassAtEnd=0,
		RecruitmentAtEnd=0,SSBatEnd=0,CatchAtEnd=0,SPR_lengthBased=0,
		Overfished=NA,Overfishing=NA,RebuildingTime=0,willNotRebuild=NA,
		RevenueAtEnd=0,CostAtEnd=0,
		ProfitAtEnd=0)
	SummaryProjectionScenarios = SummaryProjectionScenarios[-1,]
	Benchmarks = data.frame(species="",WPP="",CatchMSY_Scenario="",LifeHistory_Scenario="",populationReconstructionMethod="",SteepnessAssumed=0,
		ListFileName="",ListIndexNumber="",F_current=0,B_current=0,SSB_current=0,Yield_current=0,SPR_current=0,SSB_virgin=0,R0=0,
		Fspr40=0,Yield_spr40_kg=0,B_spr40_kg=0,SSB_spr40=0,Recruitment_spr40=0,Yield_N_spr40=0,NaturalDeaths_spr40=0,RebuildTime_spr40=0,
		Fspr30=0,Yield_spr30_kg=0,B_spr30_kg=0,SSB_spr30=0,Recruitment_spr30=0,Yield_N_spr30=0,NaturalDeaths_spr30=0,RebuildTime_spr30=0,
		Fspr20=0,Yield_spr20_kg=0,B_spr20_kg=0,SSB_spr20=0,Recruitment_spr20=0,Yield_N_spr20=0,NaturalDeaths_spr20=0,RebuildTime_spr20=0,                 
		Fmsy=0,MSY_kg=0,Bmsy_kg=0,SPR_at_MSY=0,SSB_at_MSY=0,Recruitment_at_MSY=0,MSY_N=0,NaturalDeaths_at_MSY=0,RebuildTime_MSY=0,Yield_OA_kg=0, 
		B_OA_kg=0,SPR_OA=0,SSB_OA=0,Recruitment_OA=0,Yield_N_OA=0,NaturalDeaths_OA=0,Fmey=0,Yield_MEY_kg=0,B_MEY_kg=0,SPR_MEY=0,Recruitment_MEY=0,                                           
		Yield_N_MEY=0,NaturalDeaths_MEY=0,RebuildTime_MEY=0) 
	Benchmarks = Benchmarks[-1,]
	#for(spp in 1:length(fileNamesProjectionResults_list))
	#for(spp in 1:5)
	#for(spp in 6:10)
	#for(spp in 11:15)
	#for(spp in 16:20)
	#for(spp in 21:25)
	#for(spp in 26:30)
	#for(spp in 31:35)
	for(spp in 36:43)
	{
		ProjectionList_thisSpecies = readRDS(paste(PATH_saveProjectionLists,fileNamesProjectionResults_list[spp],sep="/"),refhook=FALSE)
		Projections_NoFit = read.table(paste(PATH_saveProjectionLists,fileNamesProjectionResults_NoFit[spp],sep="/"),sep=",",header=TRUE)
		SummaryProjectionScenarios_spp = data.frame(species=rep("",length(ProjectionList_thisSpecies)),WPP=rep("",length(ProjectionList_thisSpecies)),CatchMSY_Scenario=rep("",length(ProjectionList_thisSpecies)),
			LifeHistory_Scenario=rep("",length(ProjectionList_thisSpecies)),populationReconstructionMethod=rep("",length(ProjectionList_thisSpecies)),FishingMortaltiyProjectionScenario=rep("",length(ProjectionList_thisSpecies)),
			ProjectionType=rep("",length(ProjectionList_thisSpecies)),ListFileName=rep("",length(ProjectionList_thisSpecies)),ListIndexNumber=rep(0,length(ProjectionList_thisSpecies)),
			Lmat=rep(0,length(ProjectionList_thisSpecies)),Linf=rep(0,length(ProjectionList_thisSpecies)),VBK=rep(0,length(ProjectionList_thisSpecies)),a=rep(0,length(ProjectionList_thisSpecies)),
			b=rep(0,length(ProjectionList_thisSpecies)),Mnat=rep(0,length(ProjectionList_thisSpecies)),F_estimated=rep(0,length(ProjectionList_thisSpecies)),SSB_virgin=rep(0,length(ProjectionList_thisSpecies)),
			SSB_current=rep(0,length(ProjectionList_thisSpecies)),R0=rep(0,length(ProjectionList_thisSpecies)),Yield_estimated=rep(0,length(ProjectionList_thisSpecies)),N_estimated=rep(0,length(ProjectionList_thisSpecies)),
			B_estimated=rep(0,length(ProjectionList_thisSpecies)),Recruits_estimated=rep(0,length(ProjectionList_thisSpecies)),Selex_Logistic_a_length=rep(0,length(ProjectionList_thisSpecies)),
			Selex_Logistic_b_length=rep(0,length(ProjectionList_thisSpecies)),Selex_Logistic_a_age=rep(0,length(ProjectionList_thisSpecies)),Selex_Logistic_b_age=rep(0,length(ProjectionList_thisSpecies)),
			SteepnessAssumed=rep(0,length(ProjectionList_thisSpecies)),F_project=rep(0,length(ProjectionList_thisSpecies)),NumFishAtEnd=rep(0,length(ProjectionList_thisSpecies)),BiomassAtEnd=rep(0,length(ProjectionList_thisSpecies)),
			RecruitmentAtEnd=rep(0,length(ProjectionList_thisSpecies)),SSBatEnd=rep(0,length(ProjectionList_thisSpecies)),CatchAtEnd=rep(0,length(ProjectionList_thisSpecies)),SPR_lengthBased=rep(0,length(ProjectionList_thisSpecies)),
			Overfished=rep(NA,length(ProjectionList_thisSpecies)),Overfishing=rep(NA,length(ProjectionList_thisSpecies)),RebuildingTime=rep(0,length(ProjectionList_thisSpecies)),willNotRebuild=rep(NA,length(ProjectionList_thisSpecies)),
			RevenueAtEnd=rep(0,length(ProjectionList_thisSpecies)),CostAtEnd=rep(0,length(ProjectionList_thisSpecies)),
			ProfitAtEnd=rep(0,length(ProjectionList_thisSpecies)))
		Benchmarks_spp = data.frame(species="",WPP="",CatchMSY_Scenario="",LifeHistory_Scenario="",populationReconstructionMethod="",SteepnessAssumed=0,
			ListFileName="",ListIndexNumber="",F_current=0,B_current=0,SSB_current=0,Yield_current=0,SPR_current=0,SSB_virgin=0,R0=0,
			Fspr40=0,Yield_spr40_kg=0,B_spr40_kg=0,SSB_spr40=0,Recruitment_spr40=0,Yield_N_spr40=0,NaturalDeaths_spr40=0,RebuildTime_spr40=0,
			Fspr30=0,Yield_spr30_kg=0,B_spr30_kg=0,SSB_spr30=0,Recruitment_spr30=0,Yield_N_spr30=0,NaturalDeaths_spr30=0,RebuildTime_spr30=0,
			Fspr20=0,Yield_spr20_kg=0,B_spr20_kg=0,SSB_spr20=0,Recruitment_spr20=0,Yield_N_spr20=0,NaturalDeaths_spr20=0,RebuildTime_spr20=0,                 
			Fmsy=0,MSY_kg=0,Bmsy_kg=0,SPR_at_MSY=0,SSB_at_MSY=0,Recruitment_at_MSY=0,MSY_N=0,NaturalDeaths_at_MSY=0,RebuildTime_MSY=0,Yield_OA_kg=0, 
			B_OA_kg=0,SPR_OA=0,SSB_OA=0,Recruitment_OA=0,Yield_N_OA=0,NaturalDeaths_OA=0,Fmey=0,Yield_MEY_kg=0,B_MEY_kg=0,SPR_MEY=0,Recruitment_MEY=0,                                           
			Yield_N_MEY=0,NaturalDeaths_MEY=0,RebuildTime_MEY=0)
			Benchmarks_spp = Benchmarks_spp[-1,]
		for(i in 1:length(ProjectionList_thisSpecies))
		{
			SummaryProjectionScenarios_spp$species[i]=ProjectionList_thisSpecies[[i]]$species
			SummaryProjectionScenarios_spp$WPP[i]=ProjectionList_thisSpecies[[i]]$WPP
			SummaryProjectionScenarios_spp$CatchMSY_Scenario[i]=ProjectionList_thisSpecies[[i]]$Scenario_New
			SummaryProjectionScenarios_spp$LifeHistory_Scenario[i]=ProjectionList_thisSpecies[[i]]$LifeHistory_Scenario
			SummaryProjectionScenarios_spp$populationReconstructionMethod[i]=ProjectionList_thisSpecies[[i]]$populationReconstructionMethod
			SummaryProjectionScenarios_spp$FishingMortaltiyProjectionScenario[i]=ProjectionList_thisSpecies[[i]]$FishingMortaltiyProjectionScenario
			SummaryProjectionScenarios_spp$ProjectionType[i]=ProjectionList_thisSpecies[[i]]$ProjectionType
			SummaryProjectionScenarios_spp$ListFileName[i]=fileNamesProjectionResults_list[spp]
			SummaryProjectionScenarios_spp$ListIndexNumber[i]=i
			SummaryProjectionScenarios_spp$Lmat[i]=ProjectionList_thisSpecies[[i]]$Lmat
			SummaryProjectionScenarios_spp$Linf[i]=ProjectionList_thisSpecies[[i]]$Linf
			SummaryProjectionScenarios_spp$VBK[i]=ProjectionList_thisSpecies[[i]]$VBK
			SummaryProjectionScenarios_spp$a[i]=ProjectionList_thisSpecies[[i]]$a
			SummaryProjectionScenarios_spp$b[i]=ProjectionList_thisSpecies[[i]]$b
			SummaryProjectionScenarios_spp$Mnat[i]=ProjectionList_thisSpecies[[i]]$M
			SummaryProjectionScenarios_spp$F_estimated[i]=ProjectionList_thisSpecies[[i]]$F_estimated
			SummaryProjectionScenarios_spp$SSB_virgin[i]=sum(ProjectionList_thisSpecies[[i]]$NumAtAge$SSB_virginAtAge)
			SummaryProjectionScenarios_spp$SSB_current[i]=sum(ProjectionList_thisSpecies[[i]]$NumAtAge$SSB_CurrentAtAge)
			SummaryProjectionScenarios_spp$R0[i]=ProjectionList_thisSpecies[[i]]$R0
			SummaryProjectionScenarios_spp$Yield_estimated[i]=sum(ProjectionList_thisSpecies[[i]]$NumAtAge$ExtrapolatedCatchInNumber*ProjectionList_thisSpecies[[i]]$NumAtAge$Weight_kg)
			SummaryProjectionScenarios_spp$N_estimated[i]=sum(ProjectionList_thisSpecies[[i]]$NumAtAge$Population_Extrapolated)
			SummaryProjectionScenarios_spp$B_estimated[i]=sum(ProjectionList_thisSpecies[[i]]$NumAtAge$Population_Extrapolated*ProjectionList_thisSpecies[[i]]$NumAtAge$Weight_kg)
			SummaryProjectionScenarios_spp$Recruits_estimated[i]=ProjectionList_thisSpecies[[i]]$Recruitment[1]
			SummaryProjectionScenarios_spp$Selex_Logistic_a_length[i]=ProjectionList_thisSpecies[[i]]$Selex_Logistic_a_length
			SummaryProjectionScenarios_spp$Selex_Logistic_b_length[i]=ProjectionList_thisSpecies[[i]]$Selex_Logistic_b_length
			SummaryProjectionScenarios_spp$Selex_Logistic_a_age[i]=ProjectionList_thisSpecies[[i]]$Selex_Logistic_a_age
			SummaryProjectionScenarios_spp$Selex_Logistic_b_age[i]=ProjectionList_thisSpecies[[i]]$Selex_Logistic_b_age
			SummaryProjectionScenarios_spp$SteepnessAssumed[i]=ProjectionList_thisSpecies[[i]]$SteepnessAssumed
			SummaryProjectionScenarios_spp$F_project[i]=ProjectionList_thisSpecies[[i]]$F_project
			SummaryProjectionScenarios_spp$NumFishAtEnd[i]=ProjectionList_thisSpecies[[i]]$ProjectedPopulationYear[length(ProjectionList_thisSpecies[[i]]$ProjectedPopulationYear)]
			SummaryProjectionScenarios_spp$BiomassAtEnd[i]=ProjectionList_thisSpecies[[i]]$ProjectedBiomassYear[length(ProjectionList_thisSpecies[[i]]$ProjectedBiomassYear)]
			SummaryProjectionScenarios_spp$RecruitmentAtEnd[i]=ProjectionList_thisSpecies[[i]]$Recruitment[length(ProjectionList_thisSpecies[[i]]$Recruitment)]
			SummaryProjectionScenarios_spp$SSBatEnd[i]=ProjectionList_thisSpecies[[i]]$SSB[length(ProjectionList_thisSpecies[[i]]$SSB)]
			SummaryProjectionScenarios_spp$CatchAtEnd[i]=sum(ProjectionList_thisSpecies[[i]]$ProjectedCatchAtAge[,length(ProjectionList_thisSpecies[[i]]$ProjectedCatchAtAge)]*ProjectionList_thisSpecies[[i]]$NumAtAge$Weight_kg)
			SummaryProjectionScenarios_spp$SPR_lengthBased[i]=ProjectionList_thisSpecies[[i]]$SPR_lengthBased
			SummaryProjectionScenarios_spp$Overfished[i] = ProjectionList_thisSpecies[[i]]$Overfished
			SummaryProjectionScenarios_spp$Overfishing[i] = ProjectionList_thisSpecies[[i]]$Overfishing
			SummaryProjectionScenarios_spp$RebuildingTime[i] = ProjectionList_thisSpecies[[i]]$RebuildingTime
			SummaryProjectionScenarios_spp$willNotRebuild[i] = ProjectionList_thisSpecies[[i]]$willNotRebuild
			SummaryProjectionScenarios_spp$RevenueAtEnd[i] = ProjectionList_thisSpecies[[i]]$Revenue[length(ProjectionList_thisSpecies[[i]]$Revenue)]
			SummaryProjectionScenarios_spp$CostAtEnd[i] = ProjectionList_thisSpecies[[i]]$Cost[length(ProjectionList_thisSpecies[[i]]$Cost)]
			SummaryProjectionScenarios_spp$ProfitAtEnd[i] = ProjectionList_thisSpecies[[i]]$Profit[length(ProjectionList_thisSpecies[[i]]$Profit)]
			if(ProjectionList_thisSpecies[[i]]$FishingMortaltiyProjectionScenario=="CurrentF")
			{
				Benchmarks_spp_i = data.frame(species=ProjectionList_thisSpecies[[i]]$species,WPP=ProjectionList_thisSpecies[[i]]$WPP,CatchMSY_Scenario=ProjectionList_thisSpecies[[i]]$Scenario_New,LifeHistory_Scenario=ProjectionList_thisSpecies[[i]]$LifeHistory_Scenario,
					populationReconstructionMethod=ProjectionList_thisSpecies[[i]]$populationReconstructionMethod,SteepnessAssumed=ProjectionList_thisSpecies[[i]]$SteepnessAssumed,ListFileName=fileNamesProjectionResults_list[spp],ListIndexNumber=i,F_current=ProjectionList_thisSpecies[[i]]$F_estimated,
					B_current=sum(ProjectionList_thisSpecies[[i]]$NumAtAge$Population_Extrapolated*ProjectionList_thisSpecies[[i]]$NumAtAge$Weight_kg),SSB_current=sum(ProjectionList_thisSpecies[[i]]$NumAtAge$SSB_CurrentAtAge),
					Yield_current=sum(ProjectionList_thisSpecies[[i]]$NumAtAge$ExtrapolatedCatchInNumber*ProjectionList_thisSpecies[[i]]$NumAtAge$Weight_kg),SPR_current=ProjectionList_thisSpecies[[i]]$SPR_lengthBased,SSB_virgin=ProjectionList_thisSpecies[[i]]$SSB_virgin,Recruitment_current_Equilibrum=ProjectionList_thisSpecies[[i]]$Recruitment[length(ProjectionList_thisSpecies[[i]]$Recruitment)],  
					R0=ProjectionList_thisSpecies[[i]]$R0,Fspr40=ProjectionList_thisSpecies[[i]]$Fspr40,Yield_spr40_kg=ProjectionList_thisSpecies[[i]]$Yield_spr40_kg,B_spr40_kg=ProjectionList_thisSpecies[[i]]$B_spr40_kg,SSB_spr40=ProjectionList_thisSpecies[[i]]$SSB_spr40,
					Recruitment_spr40=ProjectionList_thisSpecies[[i]]$Recruitment_spr40,Yield_N_spr40=ProjectionList_thisSpecies[[i]]$Yield_N_spr40,NaturalDeaths_spr40=ProjectionList_thisSpecies[[i]]$NaturalDeaths_spr40,
					RebuildTime_spr40=ProjectionList_thisSpecies[[i]]$RebuildTime_spr40,Fspr30=ProjectionList_thisSpecies[[i]]$Fspr30,Yield_spr30_kg=ProjectionList_thisSpecies[[i]]$Yield_spr30_kg,
					B_spr30_kg=ProjectionList_thisSpecies[[i]]$B_spr30_kg,SSB_spr30=ProjectionList_thisSpecies[[i]]$SSB_spr30,Recruitment_spr30=ProjectionList_thisSpecies[[i]]$Recruitment_spr30,Yield_N_spr30=ProjectionList_thisSpecies[[i]]$Yield_N_spr30,
					NaturalDeaths_spr30=ProjectionList_thisSpecies[[i]]$NaturalDeaths_spr30,RebuildTime_spr30=ProjectionList_thisSpecies[[i]]$RebuildTime_spr30,Fspr20=ProjectionList_thisSpecies[[i]]$Fspr20,
					Yield_spr20_kg=ProjectionList_thisSpecies[[i]]$Yield_spr20_kg,B_spr20_kg=ProjectionList_thisSpecies[[i]]$B_spr20_kg,SSB_spr20=ProjectionList_thisSpecies[[i]]$SSB_spr20,Recruitment_spr20=ProjectionList_thisSpecies[[i]]$Recruitment_spr20,Yield_N_spr20=ProjectionList_thisSpecies[[i]]$Yield_N_spr20,
					NaturalDeaths_spr20=ProjectionList_thisSpecies[[i]]$NaturalDeaths_spr20,RebuildTime_spr20=ProjectionList_thisSpecies[[i]]$RebuildTime_spr20,Fmsy=ProjectionList_thisSpecies[[i]]$Fmsy,
					MSY_kg=ProjectionList_thisSpecies[[i]]$MSY_kg,Bmsy_kg=ProjectionList_thisSpecies[[i]]$Bmsy_kg,SPR_at_MSY=ProjectionList_thisSpecies[[i]]$SPR_at_MSY,
					SSB_at_MSY=ProjectionList_thisSpecies[[i]]$SSB_at_MSY,Recruitment_at_MSY=ProjectionList_thisSpecies[[i]]$Recruitment_at_MSY,MSY_N=ProjectionList_thisSpecies[[i]]$MSY_N,
					NaturalDeaths_at_MSY=ProjectionList_thisSpecies[[i]]$NaturalDeaths_at_MSY,RebuildTime_MSY=ProjectionList_thisSpecies[[i]]$RebuildTime_MSY,Yield_OA_kg=ProjectionList_thisSpecies[[i]]$Yield_OA_kg, 
					B_OA_kg=ProjectionList_thisSpecies[[i]]$B_OA_kg,SPR_OA=ProjectionList_thisSpecies[[i]]$SPR_OA,SSB_OA=ProjectionList_thisSpecies[[i]]$SSB_OA,Recruitment_OA=ProjectionList_thisSpecies[[i]]$Recruitment_OA,
					Yield_N_OA=ProjectionList_thisSpecies[[i]]$Yield_N_OA,NaturalDeaths_OA=ProjectionList_thisSpecies[[i]]$NaturalDeaths_OA,Fmey=ProjectionList_thisSpecies[[i]]$Fmey,Yield_MEY_kg=ProjectionList_thisSpecies[[i]]$Yield_MEY_kg,
					B_MEY_kg=ProjectionList_thisSpecies[[i]]$B_MEY_kg,SPR_MEY=ProjectionList_thisSpecies[[i]]$SPR_MEY,Recruitment_MEY=ProjectionList_thisSpecies[[i]]$Recruitment_MEY,Yield_N_MEY=ProjectionList_thisSpecies[[i]]$Yield_N_MEY,
					NaturalDeaths_MEY=ProjectionList_thisSpecies[[i]]$NaturalDeaths_MEY,RebuildTime_MEY=ProjectionList_thisSpecies[[i]]$RebuildTime_MEY)
				Benchmarks_spp = rbind(Benchmarks_spp,Benchmarks_spp_i)			
			}
			if(i%%50==0)
			{
				print(paste("For species ",spp,": Summarizing iteration ",i," of ",length(ProjectionList_thisSpecies),"...",sep=""))
				flush.console()
			}
		}
		SummaryProjectionScenarios = rbind(SummaryProjectionScenarios,SummaryProjectionScenarios_spp)
		Benchmarks = rbind(Benchmarks,Benchmarks_spp)
		gc()
	}
	write.table(SummaryProjectionScenarios,paste(PATH_output,SummaryProjectionScenarios_fName,sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
	write.table(Benchmarks,paste(PATH_output,Benchmarks_fName,sep="/"),col.names=TRUE,row.names=FALSE,sep=",")


}


if(randOnMultipleProcessors==TRUE)
{
	SummaryProjectionScenarios1 = read.table(paste(PATH_output,"SummaryProjectionScenario1.csv",sep="/"),header=TRUE,sep=",")
	SummaryProjectionScenarios2 = read.table(paste(PATH_output,"SummaryProjectionScenario2.csv",sep="/"),header=TRUE,sep=",")
	SummaryProjectionScenarios3 = read.table(paste(PATH_output,"SummaryProjectionScenario3.csv",sep="/"),header=TRUE,sep=",")
	SummaryProjectionScenarios4 = read.table(paste(PATH_output,"SummaryProjectionScenario4.csv",sep="/"),header=TRUE,sep=",")
	SummaryProjectionScenarios5 = read.table(paste(PATH_output,"SummaryProjectionScenario5.csv",sep="/"),header=TRUE,sep=",")
	SummaryProjectionScenarios6 = read.table(paste(PATH_output,"SummaryProjectionScenario6.csv",sep="/"),header=TRUE,sep=",")
	SummaryProjectionScenarios7 = read.table(paste(PATH_output,"SummaryProjectionScenario7.csv",sep="/"),header=TRUE,sep=",")
	SummaryProjectionScenarios8 = read.table(paste(PATH_output,"SummaryProjectionScenario8.csv",sep="/"),header=TRUE,sep=",")
	SummaryProjectionScenarios = rbind(SummaryProjectionScenarios1,SummaryProjectionScenarios2,SummaryProjectionScenarios3,SummaryProjectionScenarios4,
		SummaryProjectionScenarios5,SummaryProjectionScenarios6,SummaryProjectionScenarios7,SummaryProjectionScenarios8)
	write.table(SummaryProjectionScenarios,paste(PATH_output,"SummaryProjectionScenarios.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
	Benchmarks1 = read.table(paste(PATH_output,"Benchmarks1.csv",sep="/"),header=TRUE,sep=",")
	Benchmarks2 = read.table(paste(PATH_output,"Benchmarks2.csv",sep="/"),header=TRUE,sep=",")
	Benchmarks3 = read.table(paste(PATH_output,"Benchmarks3.csv",sep="/"),header=TRUE,sep=",")
	Benchmarks4 = read.table(paste(PATH_output,"Benchmarks4.csv",sep="/"),header=TRUE,sep=",")
	Benchmarks5 = read.table(paste(PATH_output,"Benchmarks5.csv",sep="/"),header=TRUE,sep=",")
	Benchmarks6 = read.table(paste(PATH_output,"Benchmarks6.csv",sep="/"),header=TRUE,sep=",")
	Benchmarks7 = read.table(paste(PATH_output,"Benchmarks7.csv",sep="/"),header=TRUE,sep=",")
	Benchmarks8 = read.table(paste(PATH_output,"Benchmarks8.csv",sep="/"),header=TRUE,sep=",")
	Benchmarks = rbind(Benchmarks1,Benchmarks2,Benchmarks3,Benchmarks4,
		Benchmarks5,Benchmarks6,Benchmarks7,Benchmarks8)
	write.table(Benchmarks,paste(PATH_output,"Benchmarks.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
}



if(calculateRebuildingDifferentWay==TRUE)
{
	SummaryProjectionScenarios = read.table(paste(PATH_output,"SummaryProjectionScenarios.csv",sep="/"),header=TRUE,sep=",")
	Benchmarks = read.table(paste(PATH_output,"Benchmarks.csv",sep="/"),header=TRUE,sep=",")
	SummaryProjectionScenarios$RebuildTime_to_spr20_New=NA
	firstTimeLoop=TRUE
	#total is 3109584, thus across 8 processors, 388698 on each processor - do it via subsetting
	#SummaryProjectionScenarios = SummaryProjectionScenarios[1:388698,]
	#SummaryProjectionScenarios = SummaryProjectionScenarios[388699:777397,]
	#SummaryProjectionScenarios = SummaryProjectionScenarios[777398:1166096,]	
	#SummaryProjectionScenarios = SummaryProjectionScenarios[1166097:1554795,]	
	SummaryProjectionScenarios = SummaryProjectionScenarios[1554796:1943494,]	
	#SummaryProjectionScenarios = SummaryProjectionScenarios[1943495:2332193,]	
	#SummaryProjectionScenarios = SummaryProjectionScenarios[2332194:2720892,]	
	#SummaryProjectionScenarios = SummaryProjectionScenarios[2720893:3109584,]
	for(i in 1:dim(SummaryProjectionScenarios)[1])
	{
		if(firstTimeLoop==TRUE)
		{
			ProjectionList_thisSpecies = readRDS(paste(PATH_saveProjectionLists,SummaryProjectionScenarios$ListFileName[i],sep="/"),refhook=FALSE)		
		}
		if(firstTimeLoop==FALSE)
		{
			if(SummaryProjectionScenarios$ListFileName[i]!=SummaryProjectionScenarios$ListFileName[i-1])
			{
				ProjectionList_thisSpecies = readRDS(paste(PATH_saveProjectionLists,SummaryProjectionScenarios$ListFileName[i],sep="/"),refhook=FALSE)		
			}
		}
		firstTimeLoop=FALSE
		Projection_i = ProjectionList_thisSpecies[[SummaryProjectionScenarios$ListIndexNumber[i]]]
		Benchmark_i = subset(Benchmarks,Benchmarks$species==SummaryProjectionScenarios$species[i] & Benchmarks$WPP==SummaryProjectionScenarios$WPP[i] &
			Benchmarks$CatchMSY_Scenario==SummaryProjectionScenarios$CatchMSY_Scenario[i] & Benchmarks$LifeHistory_Scenario==SummaryProjectionScenarios$LifeHistory_Scenario[i] & 
			Benchmarks$populationReconstructionMethod==SummaryProjectionScenarios$populationReconstructionMethod[i] & Benchmarks$SteepnessAssumed==SummaryProjectionScenarios$SteepnessAssumed[i]) 
		#Rebuilding times for each "FishingMortaltiyProjectionScenario"
		BiomassRelativeToBenchmark = Projection_i$ProjectedBiomassYear/Benchmark_i$B_spr20_kg
		RebuildTime=NA
		RebuildTime = suppressWarnings(min(as.numeric(names(which(BiomassRelativeToBenchmark>=0.95)))) - min(as.numeric(names(BiomassRelativeToBenchmark))))
		if(is.infinite(RebuildTime)==TRUE | is.na(RebuildTime)) 
		{
			RebuildTime=NA
		}
		SummaryProjectionScenarios$RebuildTime_to_spr20_New[i]=RebuildTime
		if(i%%500==0)
		{
			print(paste("Working on row ",i," of ",dim(SummaryProjectionScenarios)[1],"...",sep=""))
			flush.console()	
		}
	}
	#write.table(SummaryProjectionScenarios,paste(PATH_output,"SummaryProjectionScenarios_rebuildingToSPR20_1.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
	#write.table(SummaryProjectionScenarios,paste(PATH_output,"SummaryProjectionScenarios_rebuildingToSPR20_2.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
	#write.table(SummaryProjectionScenarios,paste(PATH_output,"SummaryProjectionScenarios_rebuildingToSPR20_3.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
	#write.table(SummaryProjectionScenarios,paste(PATH_output,"SummaryProjectionScenarios_rebuildingToSPR20_4.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
	write.table(SummaryProjectionScenarios,paste(PATH_output,"SummaryProjectionScenarios_rebuildingToSPR20_5.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
	#write.table(SummaryProjectionScenarios,paste(PATH_output,"SummaryProjectionScenarios_rebuildingToSPR20_6.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
	#write.table(SummaryProjectionScenarios,paste(PATH_output,"SummaryProjectionScenarios_rebuildingToSPR20_7.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
	#write.table(SummaryProjectionScenarios,paste(PATH_output,"SummaryProjectionScenarios_rebuildingToSPR20_8.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")


############################# Read in tables and combine ############################

	SummaryProjectionScenarios1 = read.table(paste(PATH_output,"SummaryProjectionScenarios_rebuildingToSPR20_1.csv",sep="/"),header=TRUE,sep=",")
	SummaryProjectionScenarios2 = read.table(paste(PATH_output,"SummaryProjectionScenarios_rebuildingToSPR20_2.csv",sep="/"),header=TRUE,sep=",")
	SummaryProjectionScenarios3 = read.table(paste(PATH_output,"SummaryProjectionScenarios_rebuildingToSPR20_3.csv",sep="/"),header=TRUE,sep=",")
	SummaryProjectionScenarios4 = read.table(paste(PATH_output,"SummaryProjectionScenarios_rebuildingToSPR20_4.csv",sep="/"),header=TRUE,sep=",")
	SummaryProjectionScenarios5 = read.table(paste(PATH_output,"SummaryProjectionScenarios_rebuildingToSPR20_5.csv",sep="/"),header=TRUE,sep=",")
	SummaryProjectionScenarios6 = read.table(paste(PATH_output,"SummaryProjectionScenarios_rebuildingToSPR20_6.csv",sep="/"),header=TRUE,sep=",")
	SummaryProjectionScenarios7 = read.table(paste(PATH_output,"SummaryProjectionScenarios_rebuildingToSPR20_7.csv",sep="/"),header=TRUE,sep=",")
	SummaryProjectionScenarios8 = read.table(paste(PATH_output,"SummaryProjectionScenarios_rebuildingToSPR20_8.csv",sep="/"),header=TRUE,sep=",")
	SummaryProjectionScenarios = rbind(SummaryProjectionScenarios1,SummaryProjectionScenarios2,SummaryProjectionScenarios3,SummaryProjectionScenarios4,
		SummaryProjectionScenarios5,SummaryProjectionScenarios6,SummaryProjectionScenarios7,SummaryProjectionScenarios8)
	write.table(SummaryProjectionScenarios,paste(PATH_output,"SummaryProjectionScenarios_rebuildingToSPR20.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")

}



#############################################################################################################################################
########################## Calculate Ratios, Filter Output and Remote Outlyer Runs and Estimates ############################################
#############################################################################################################################################

#Remove Projection and Estimation Outliers
if(removeProjectionOutliers==TRUE)
{
	#Use geometric mean and geometric sd for biomass and some other items depending on how they are distributed. 
	geoMean = function(x)
	{
		return(exp(mean(log(x))))	
	}
	geoStDev = function(x)
	{
		gm1 = mean(log(x), na.rm = T)
		cil = exp(gm1-(1.96*(sd(log(x), na.rm = T)/sqrt(length(x)))))
		ciupp = exp(gm1+(1.96*(sd(log(x), na.rm = T)/sqrt(length(x)))))
		vec = c(round(cil,2), round(ciupp,2))
		return((vec[2]-vec[1])/2)
	}	
	geoMarginError = function(x)
	{
		return(exp((1.96*sd(log(x)))/length(x)))	
	}	
	geoCV = function(x)
	{
		return(sqrt(exp((sd(log(x)))^2)-1))
	}
	Benchmarks = read.table(paste(PATH_output,"Benchmarks.csv",sep="/"),header=TRUE,sep=",")
	SummaryProjectionScenarios = read.table(paste(PATH_output,"SummaryProjectionScenarios_rebuildingToSPR20.csv",sep="/"),header=TRUE,sep=",")
	#Pull the new rebuilding times from the SummaryProjectionScenarios object and merge them into the Benchmarks data.frame
	RebuildTimes_New = subset(SummaryProjectionScenarios,select=c(species,WPP,CatchMSY_Scenario,LifeHistory_Scenario,
		populationReconstructionMethod,SteepnessAssumed,FishingMortaltiyProjectionScenario,RebuildTime_to_spr20_New))
	RebuildTimes_New$idVariable = paste(RebuildTimes_New$species,RebuildTimes_New$WPP,RebuildTimes_New$CatchMSY_Scenario,
		RebuildTimes_New$LifeHistory_Scenario,RebuildTimes_New$populationReconstructionMethod,RebuildTimes_New$SteepnessAssumed,sep="_")
	RebuildTimes_New_reshaped = reshape(RebuildTimes_New,v.names="RebuildTime_to_spr20_New",idvar="idVariable",timevar="FishingMortaltiyProjectionScenario",direction="wide")
	RebuildTimes_New_reshaped$idVariable=NULL
	Benchmarks = merge(Benchmarks,RebuildTimes_New_reshaped,all.x=TRUE,by=c("species","WPP","CatchMSY_Scenario","LifeHistory_Scenario","populationReconstructionMethod","SteepnessAssumed"))
############################################### TEMP CODE: Use Biomass from Catch-MSY at Biomass_current ###############################################
	#Merge biomass estimates made using Catch-MSY into this data frame to compare with Biomass estimates at Fcurrent calculated at equilibrium
	BiomassCarryingCapacityMSY = read.table(paste(PATH_output,"BiomassCarryingCapacityMSY.csv",sep="/"),header=TRUE,sep=",")
	BiomassCarryingCapacityMSY = subset(BiomassCarryingCapacityMSY,select=c(species,WPP,Scenario_New,BiomassMean_LastYear,catch_t))
	names(BiomassCarryingCapacityMSY)[names(BiomassCarryingCapacityMSY)=="Scenario_New"]="CatchMSY_Scenario"
	names(BiomassCarryingCapacityMSY)[names(BiomassCarryingCapacityMSY)=="BiomassMean_LastYear"]="BiomassCurrent_CatchMSY_mT"
	names(BiomassCarryingCapacityMSY)[names(BiomassCarryingCapacityMSY)=="catch_t"]="Catch_CatchMSY_mT"
	Benchmarks = merge(Benchmarks,BiomassCarryingCapacityMSY,all.x=TRUE,by=c("species","WPP","CatchMSY_Scenario"))
	Benchmarks$B_current = Benchmarks$BiomassCurrent_CatchMSY_mT*1000
	Benchmarks$Yield_current = Benchmarks$Catch_CatchMSY_mT*1000	
############################################## END TEMP CODE ############################################################################################	
	#Calculate Ratios
	Benchmarks$F_over_Fspr20 = Benchmarks$F_current/Benchmarks$Fspr20
	Benchmarks$B_over_Bspr20 = Benchmarks$B_current/Benchmarks$B_spr20_kg
	Benchmarks$SSB_over_SSBspr20 = Benchmarks$SSB_current/Benchmarks$SSB_spr20
	Benchmarks$Yield_over_Yieldspr20 = Benchmarks$Yield_current/Benchmarks$Yield_spr20_kg
	Benchmarks$Recruit_over_RecruitSpr20 = Benchmarks$Recruitment_current_Equilibrum/Benchmarks$Recruitment_spr20
	Benchmarks$F_over_Fspr30 = Benchmarks$F_current/Benchmarks$Fspr30
	Benchmarks$B_over_Bspr30 = Benchmarks$B_current/Benchmarks$B_spr30_kg
	Benchmarks$SSB_over_SSBspr30 = Benchmarks$SSB_current/Benchmarks$SSB_spr30
	Benchmarks$Yield_over_Yieldspr30 = Benchmarks$Yield_current/Benchmarks$Yield_spr30_kg
	Benchmarks$Recruit_over_RecruitSpr30 = Benchmarks$Recruitment_current_Equilibrum/Benchmarks$Recruitment_spr30
	Benchmarks$F_over_Fspr40 = Benchmarks$F_current/Benchmarks$Fspr40
	Benchmarks$B_over_Bspr40 = Benchmarks$B_current/Benchmarks$B_spr40_kg
	Benchmarks$SSB_over_SSBspr40 = Benchmarks$SSB_current/Benchmarks$SSB_spr40
	Benchmarks$Yield_over_Yieldspr40 = Benchmarks$Yield_current/Benchmarks$Yield_spr40_kg
	Benchmarks$Recruit_over_RecruitSpr40 = Benchmarks$Recruitment_current_Equilibrum/Benchmarks$Recruitment_spr40
	Benchmarks$F_over_Fmsy = Benchmarks$F_current/Benchmarks$Fmsy
	Benchmarks$B_over_Bmsy = Benchmarks$B_current/Benchmarks$Bmsy_kg
	Benchmarks$SSB_over_SSBmsy = Benchmarks$SSB_current/Benchmarks$SSB_at_MSY
	Benchmarks$Yield_over_YieldMSY = Benchmarks$Yield_current/Benchmarks$MSY_kg
	Benchmarks$Recruit_over_RecruitMSY = Benchmarks$Recruitment_current_Equilibrum/Benchmarks$Recruitment_at_MSY
	Benchmarks$F_over_Fmey = Benchmarks$F_current/Benchmarks$Fmey
	Benchmarks$B_over_Bmey = Benchmarks$B_current/Benchmarks$B_MEY_kg
	Benchmarks$Yield_over_YieldMEY = Benchmarks$Yield_current/Benchmarks$Yield_MEY_kg
	Benchmarks$Recruit_over_RecruitMEY = Benchmarks$Recruitment_current_Equilibrum/Benchmarks$Recruitment_MEY
	#Use Median and Median Absolute Deviation (see function below) to estimate the central tendancy and error since the distributions are so skewed with a long tail of large values, even when truncated
	MedianAbsoluteDeviation = function(x)
	{
		x = subset(x,!is.na(x))
		Med = median(x)
		return(median(abs(Med-x)))
	}	
	#Remove outlyers: Fishing Mortality
	MinQuant_F = 0     #set quantiles: truncation takes care of entire lower limit
	MaxQuant_F = 0.95   #set quantiles: still need to bound and rein in the long tail of few large estimated values that pull central tendancy artificially higher 
	#Set up variables to flag records to remove 
	Benchmarks$F_current_FLAG=0
	Benchmarks$Fspr20_FLAG=0
	Benchmarks$Fspr30_FLAG=0
	Benchmarks$Fspr40_FLAG=0
	Benchmarks$Fmsy_FLAG=0
	Benchmarks$Fmey_FLAG=0	
	Benchmarks$Fratio20_FLAG=0
	Benchmarks$Fratio30_FLAG=0
	Benchmarks$Fratio40_FLAG=0	
	Benchmarks$Fmsy_ratio_FLAG=0
	Benchmarks$Fmey_ratio_FLAG=0
	#first pass, remove records that are not plausable, i.e. can't have F larger than 5 or less than 0
	Benchmarks$flag=1
	Benchmarks$F_current_FLAG[Benchmarks$F_current > 4 | Benchmarks$F_current < 0.01]=1
	Benchmarks$Fspr20_FLAG[Benchmarks$Fspr20 > 5 | Benchmarks$Fspr20 < 0]=1
	Benchmarks$Fspr30_FLAG[Benchmarks$Fspr30 > 5 | Benchmarks$Fspr30 < 0]=1
	Benchmarks$Fspr40_FLAG[Benchmarks$Fspr40 > 5 | Benchmarks$Fspr40 < 0]=1
	Benchmarks$Fmsy_FLAG[Benchmarks$Fmsy > 4 | Benchmarks$Fmsy < 0]=1
	Benchmarks$Fmey_FLAG[Benchmarks$Fmey > 4 | Benchmarks$Fmey < 0]=1
	Benchmarks$Fratio20_FLAG[Benchmarks$F_over_Fspr20 > 5 | Benchmarks$F_over_Fspr20 < 0]=1
	Benchmarks$Fratio20_FLAG[is.infinite(Benchmarks$F_over_Fspr20)]=1
	Benchmarks$Fratio30_FLAG[Benchmarks$F_over_Fspr30 > 5 | Benchmarks$F_over_Fspr30 < 0]=1
	Benchmarks$Fratio30_FLAG[is.infinite(Benchmarks$F_over_Fspr30)]=1
	Benchmarks$Fratio40_FLAG[Benchmarks$F_over_Fspr40 > 5 | Benchmarks$F_over_Fspr40 < 0]=1
	Benchmarks$Fratio40_FLAG[is.infinite(Benchmarks$F_over_Fspr40)]=1
	Benchmarks$Fmsy_ratio_FLAG[Benchmarks$F_over_Fmsy > 5 | Benchmarks$F_over_Fmsy < 0]=1
	Benchmarks$Fmsy_ratio_FLAG[is.infinite(Benchmarks$F_over_Fmsy)]=1
	Benchmarks$Fmey_ratio_FLAG[Benchmarks$F_over_Fmey > 5 | Benchmarks$F_over_Fmey < 0]=1
	Benchmarks$Fmey_ratio_FLAG[is.infinite(Benchmarks$F_over_Fmey)]=1
	Benchmarks$F_current_ln = suppressWarnings(log(Benchmarks$F_current))  
	Benchmarks$Fspr20_ln = suppressWarnings(log(Benchmarks$Fspr20))  
	Benchmarks$Fspr30_ln = suppressWarnings(log(Benchmarks$Fspr30))  
	Benchmarks$Fspr40_ln = suppressWarnings(log(Benchmarks$Fspr40))  
	Benchmarks$Fmsy_ln = suppressWarnings(log(Benchmarks$Fmsy))  
	Benchmarks$Fmey_ln = suppressWarnings(log(Benchmarks$Fmey))  
	Benchmarks$F_over_Fspr20_ln = suppressWarnings(log(Benchmarks$F_over_Fspr20))
	Benchmarks$F_over_Fspr30_ln = suppressWarnings(log(Benchmarks$F_over_Fspr30))
	Benchmarks$F_over_Fspr40_ln = suppressWarnings(log(Benchmarks$F_over_Fspr40))
	Benchmarks$F_over_Fmsy_ln = suppressWarnings(log(Benchmarks$F_over_Fmsy))
	Benchmarks$F_over_Fmey_ln = suppressWarnings(log(Benchmarks$F_over_Fmey))
	SppWpp = subset(Benchmarks,select=c(species,WPP))
	SppWpp = unique(SppWpp)
	SppWpp = SppWpp[order(SppWpp$species,SppWpp$WPP),]
	for(i in 1:dim(SppWpp)[1])
	{
		#F_current: use quantiles to refine the natural log of the distribution and remove last outliers
		thisSpp = Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i] & Benchmarks$F_current_FLAG==0,]		
		Quant = quantile(thisSpp$F_current_ln,probs=c(MinQuant_F,MaxQuant_F),na.rm=TRUE,type=8)
		Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$F_current_FLAG[Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$F_current_ln < Quant[1] | Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$F_current_ln > Quant[2]]=1	
		#Fspr20: use quantiles to refine the natural log of the distribution and remove last outliers
		thisSpp = Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i] & Benchmarks$Benchmarks$Fspr20_FLAG==0,]		
		Quant = quantile(thisSpp$Fspr20_ln,probs=c(MinQuant_F,MaxQuant_F),na.rm=TRUE,type=8)
		Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Fspr20_FLAG[Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Fspr20_ln < Quant[1] | Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Fspr20_ln > Quant[2]]=1	
		#Fspr30: use quantiles to refine the natural log of the distribution and remove last outliers
		thisSpp = Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i] & Benchmarks$Benchmarks$Fspr30_FLAG==0,]		
		Quant = quantile(thisSpp$Fspr30_ln,probs=c(MinQuant_F,MaxQuant_F),na.rm=TRUE,type=8)
		Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Fspr30_FLAG[Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Fspr30_ln < Quant[1] | Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Fspr30_ln > Quant[2]]=1	
		#Fspr40: use quantiles to refine the natural log of the distribution and remove last outliers
		thisSpp = Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i] & Benchmarks$Benchmarks$Fspr40_FLAG==0,]		
		Quant = quantile(thisSpp$Fspr40_ln,probs=c(MinQuant_F,MaxQuant_F),na.rm=TRUE,type=8)
		Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Fspr40_FLAG[Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Fspr40_ln < Quant[1] | Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Fspr40_ln > Quant[2]]=1	
		#Fmsy: use quantiles to refine the natural log of the distribution and remove last outliers
		thisSpp = Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i] & Benchmarks$Benchmarks$Fmsy_FLAG==0,]		
		Quant = quantile(thisSpp$Fmsy_ln,probs=c(MinQuant_F,MaxQuant_F),na.rm=TRUE,type=8)
		Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Fmsy_FLAG[Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Fmsy_ln < Quant[1] | Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Fmsy_ln > Quant[2]]=1	
		#Fmey: use quantiles to refine the natural log of the distribution and remove last outliers
		thisSpp = Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i] & Benchmarks$Benchmarks$Fmey_FLAG==0,]		
		Quant = quantile(thisSpp$Fmey_ln,probs=c(MinQuant_F,MaxQuant_F),na.rm=TRUE,type=8)
		Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Fmey_FLAG[Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Fmey_ln < Quant[1] | Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Fmey_ln > Quant[2]]=1	
		#F_over_Fspr20: Use Quantiles now to refine the natrual log of the distribution
		thisSpp = Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i] & Benchmarks$Benchmarks$Fratio20_FLAG==0,]	
		Quant = quantile(thisSpp$F_over_Fspr20_ln,probs=c(MinQuant_F,MaxQuant_F),na.rm=TRUE,type=8)
		Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Fratio20_FLAG[Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$F_over_Fspr20_ln < Quant[1] | Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$F_over_Fspr20_ln > Quant[2]]=1
		#F_over_Fspr30: Use Quantiles now to refine the natrual log of the distribution
		thisSpp = Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i] & Benchmarks$Benchmarks$Fratio30_FLAG==0,]	
		Quant = quantile(thisSpp$F_over_Fspr30_ln,probs=c(MinQuant_F,MaxQuant_F),na.rm=TRUE,type=8)
		Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Fratio30_FLAG[Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$F_over_Fspr30_ln < Quant[1] | Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$F_over_Fspr30_ln > Quant[2]]=1
		#F_over_Fspr40: Use Quantiles now to refine the natrual log of the distribution
		thisSpp = Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i] & Benchmarks$Benchmarks$Fratio40_FLAG==0,]	
		Quant = quantile(thisSpp$F_over_Fspr40_ln,probs=c(MinQuant_F,MaxQuant_F),na.rm=TRUE,type=8)
		Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Fratio40_FLAG[Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$F_over_Fspr40_ln < Quant[1] | Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$F_over_Fspr40_ln > Quant[2]]=1
		#F_over_Fmsy: Use Quantiles now to refine the natrual log of the distribution
		thisSpp = Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i] & Benchmarks$Benchmarks$Fmsy_ratio_FLAG==0,]	
		Quant = quantile(thisSpp$F_over_Fmsy_ln,probs=c(MinQuant_F,MaxQuant_F),na.rm=TRUE,type=8)
		Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Fmsy_ratio_FLAG[Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$F_over_Fmsy_ln < Quant[1] | Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$F_over_Fmsy_ln > Quant[2]]=1
		#F_over_Fmey: Use Quantiles now to refine the natrual log of the distribution
		thisSpp = Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i] & Benchmarks$Benchmarks$Fmey_ratio_FLAG==0,]	
		Quant = quantile(thisSpp$F_over_Fmey_ln,probs=c(MinQuant_F,MaxQuant_F),na.rm=TRUE,type=8)
		Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Fmey_ratio_FLAG[Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$F_over_Fmey_ln < Quant[1] | Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$F_over_Fmey_ln > Quant[2]]=1
		if(i%%10==0)
		{
			print(paste("Working on species ",SppWpp$species[i]," and WPP ",SppWpp$WPP[i],"...",sep=""))
			flush.console()	
		}
	}
	#Fcurrent
	Bench_F_current = subset(Benchmarks,Benchmarks$F_current_FLAG==0) # & Benchmarks$F_current_FLAG==0)
	Bench_F_current_mean = aggregate.data.frame(Bench_F_current$F_current,by=list(Bench_F_current$species,Bench_F_current$WPP),FUN=median,na.rm=TRUE)   #median   geoMean   geoStDev
	names(Bench_F_current_mean) = c("species","WPP","F_current_median")
	Samples = aggregate.data.frame(list(Benchmarks$F_current_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_F_current = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_F_current = Samples$ScenariosIncluded_F_current/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_F_current_mean = merge(Bench_F_current_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_F_current_mean$ObsNumForLower = ceiling((Bench_F_current_mean$ScenariosIncluded_F_current*0.5) - (qnorm(0.99)*sqrt(Bench_F_current_mean$ScenariosIncluded_F_current*0.5)))
	Bench_F_current_mean$ObsNumForUpper = ceiling((Bench_F_current_mean$ScenariosIncluded_F_current*0.5) + (qnorm(0.99)*sqrt(Bench_F_current_mean$ScenariosIncluded_F_current*0.5)))
	ItemsToLoop = data.frame(species=Bench_F_current$species,WPP=Bench_F_current$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$F_current_UpperCI=NA
	ItemsToLoop$F_current_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_F_current,Bench_F_current$species==ItemsToLoop$species[i] & Bench_F_current$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$F_current),]
		Bench_F_current_mean_Sub = subset(Bench_F_current_mean,Bench_F_current_mean$species==ItemsToLoop$species[i] & Bench_F_current_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$F_current_UpperCI[i]=Sub$F_current[Bench_F_current_mean_Sub$ObsNumForUpper]
		ItemsToLoop$F_current_LowerCI[i]=Sub$F_current[Bench_F_current_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_F_current = merge(Bench_F_current_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_F_current$ScenariosIncluded_F_current=NULL
	Bench_F_current$FractionScenariosIncluded_F_current=NULL
	Bench_F_current$ObsNumForLower=NULL
	Bench_F_current$ObsNumForUpper=NULL
	#Fspr20
	Bench_Fspr20 = subset(Benchmarks,Benchmarks$Fspr20_FLAG==0)  
	Bench_Fspr20_mean = aggregate.data.frame(Bench_Fspr20$Fspr20,by=list(Bench_Fspr20$species,Bench_Fspr20$WPP),FUN=median,na.rm=TRUE)
	names(Bench_Fspr20_mean) = c("species","WPP","Fspr20_median")
	Samples = aggregate.data.frame(list(Benchmarks$Fspr20_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_Fspr20 = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_Fspr20 = Samples$ScenariosIncluded_Fspr20/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_Fspr20_mean = merge(Bench_Fspr20_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_Fspr20_mean$ObsNumForLower = ceiling((Bench_Fspr20_mean$ScenariosIncluded_Fspr20*0.5) - (qnorm(0.99)*sqrt(Bench_Fspr20_mean$ScenariosIncluded_Fspr20*0.5)))
	Bench_Fspr20_mean$ObsNumForUpper = ceiling((Bench_Fspr20_mean$ScenariosIncluded_Fspr20*0.5) + (qnorm(0.99)*sqrt(Bench_Fspr20_mean$ScenariosIncluded_Fspr20*0.5)))
	ItemsToLoop = data.frame(species=Bench_Fspr20_mean$species,WPP=Bench_Fspr20_mean$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$Fspr20_UpperCI=NA
	ItemsToLoop$Fspr20_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_Fspr20,Bench_Fspr20$species==ItemsToLoop$species[i] & Bench_Fspr20$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$Fspr20),]
		Bench_Fspr20_mean_Sub = subset(Bench_Fspr20_mean,Bench_Fspr20_mean$species==ItemsToLoop$species[i] & Bench_Fspr20_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$Fspr20_UpperCI[i]=Sub$Fspr20[Bench_Fspr20_mean_Sub$ObsNumForUpper]
		ItemsToLoop$Fspr20_LowerCI[i]=Sub$Fspr20[Bench_Fspr20_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_Fspr20 = merge(Bench_Fspr20_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_Fspr20$ScenariosIncluded_Fspr20=NULL
	Bench_Fspr20$FractionScenariosIncluded_Fspr20=NULL
	Bench_Fspr20$ObsNumForLower=NULL
	Bench_Fspr20$ObsNumForUpper=NULL	
	#Fspr30
	Bench_Fspr30 = subset(Benchmarks,Benchmarks$Fspr30_FLAG==0)  
	Bench_Fspr30_mean = aggregate.data.frame(Bench_Fspr30$Fspr30,by=list(Bench_Fspr30$species,Bench_Fspr30$WPP),FUN=median,na.rm=TRUE)
	names(Bench_Fspr30_mean) = c("species","WPP","Fspr30_median")
	Samples = aggregate.data.frame(list(Benchmarks$Fspr30_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_Fspr30 = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_Fspr30 = Samples$ScenariosIncluded_Fspr30/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_Fspr30_mean = merge(Bench_Fspr30_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_Fspr30_mean$ObsNumForLower = ceiling((Bench_Fspr30_mean$ScenariosIncluded_Fspr30*0.5) - (qnorm(0.99)*sqrt(Bench_Fspr30_mean$ScenariosIncluded_Fspr30*0.5)))
	Bench_Fspr30_mean$ObsNumForUpper = ceiling((Bench_Fspr30_mean$ScenariosIncluded_Fspr30*0.5) + (qnorm(0.99)*sqrt(Bench_Fspr30_mean$ScenariosIncluded_Fspr30*0.5)))
	ItemsToLoop = data.frame(species=Bench_Fspr30_mean$species,WPP=Bench_Fspr30_mean$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$Fspr30_UpperCI=NA
	ItemsToLoop$Fspr30_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_Fspr30,Bench_Fspr30$species==ItemsToLoop$species[i] & Bench_Fspr30$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$Fspr30),]
		Bench_Fspr30_mean_Sub = subset(Bench_Fspr30_mean,Bench_Fspr30_mean$species==ItemsToLoop$species[i] & Bench_Fspr30_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$Fspr30_UpperCI[i]=Sub$Fspr30[Bench_Fspr30_mean_Sub$ObsNumForUpper]
		ItemsToLoop$Fspr30_LowerCI[i]=Sub$Fspr30[Bench_Fspr30_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_Fspr30 = merge(Bench_Fspr30_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_Fspr30$ScenariosIncluded_Fspr30=NULL
	Bench_Fspr30$FractionScenariosIncluded_Fspr30=NULL
	Bench_Fspr30$ObsNumForLower=NULL
	Bench_Fspr30$ObsNumForUpper=NULL	
	#Fspr40
	Bench_Fspr40 = subset(Benchmarks,Benchmarks$Fspr40_FLAG==0)  
	Bench_Fspr40_mean = aggregate.data.frame(Bench_Fspr40$Fspr40,by=list(Bench_Fspr40$species,Bench_Fspr40$WPP),FUN=median,na.rm=TRUE)
	names(Bench_Fspr40_mean) = c("species","WPP","Fspr40_median")
	Samples = aggregate.data.frame(list(Benchmarks$Fspr40_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_Fspr40 = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_Fspr40 = Samples$ScenariosIncluded_Fspr40/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_Fspr40_mean = merge(Bench_Fspr40_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_Fspr40_mean$ObsNumForLower = ceiling((Bench_Fspr40_mean$ScenariosIncluded_Fspr40*0.5) - (qnorm(0.99)*sqrt(Bench_Fspr40_mean$ScenariosIncluded_Fspr40*0.5)))
	Bench_Fspr40_mean$ObsNumForUpper = ceiling((Bench_Fspr40_mean$ScenariosIncluded_Fspr40*0.5) + (qnorm(0.99)*sqrt(Bench_Fspr40_mean$ScenariosIncluded_Fspr40*0.5)))
	ItemsToLoop = data.frame(species=Bench_Fspr40_mean$species,WPP=Bench_Fspr40_mean$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$Fspr40_UpperCI=NA
	ItemsToLoop$Fspr40_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_Fspr40,Bench_Fspr40$species==ItemsToLoop$species[i] & Bench_Fspr40$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$Fspr40),]
		Bench_Fspr40_mean_Sub = subset(Bench_Fspr40_mean,Bench_Fspr40_mean$species==ItemsToLoop$species[i] & Bench_Fspr40_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$Fspr40_UpperCI[i]=Sub$Fspr40[Bench_Fspr40_mean_Sub$ObsNumForUpper]
		ItemsToLoop$Fspr40_LowerCI[i]=Sub$Fspr40[Bench_Fspr40_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_Fspr40 = merge(Bench_Fspr40_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_Fspr40$ScenariosIncluded_Fspr40=NULL 
	Bench_Fspr40$FractionScenariosIncluded_Fspr40=NULL 
	Bench_Fspr40$ObsNumForLower=NULL 
	Bench_Fspr40$ObsNumForUpper=NULL	
	#Fmsy
	Bench_Fmsy = subset(Benchmarks,Benchmarks$Fmsy_FLAG==0)  
	Bench_Fmsy_mean = aggregate.data.frame(Bench_Fmsy$Fmsy,by=list(Bench_Fmsy$species,Bench_Fmsy$WPP),FUN=median,na.rm=TRUE)
	names(Bench_Fmsy_mean) = c("species","WPP","Fmsy_median")
	Samples = aggregate.data.frame(list(Benchmarks$Fmsy_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_Fmsy = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_Fmsy = Samples$ScenariosIncluded_Fmsy/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_Fmsy_mean = merge(Bench_Fmsy_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_Fmsy_mean$ObsNumForLower = ceiling((Bench_Fmsy_mean$ScenariosIncluded_Fmsy*0.5) - (qnorm(0.99)*sqrt(Bench_Fmsy_mean$ScenariosIncluded_Fmsy*0.5)))
	Bench_Fmsy_mean$ObsNumForUpper = ceiling((Bench_Fmsy_mean$ScenariosIncluded_Fmsy*0.5) + (qnorm(0.99)*sqrt(Bench_Fmsy_mean$ScenariosIncluded_Fmsy*0.5)))
	ItemsToLoop = data.frame(species=Bench_Fmsy_mean$species,WPP=Bench_Fmsy_mean$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$Fmsy_UpperCI=NA
	ItemsToLoop$Fmsy_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_Fmsy,Bench_Fmsy$species==ItemsToLoop$species[i] & Bench_Fmsy$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$Fmsy),]
		Bench_Fmsy_mean_Sub = subset(Bench_Fmsy_mean,Bench_Fmsy_mean$species==ItemsToLoop$species[i] & Bench_Fmsy_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$Fmsy_UpperCI[i]=Sub$Fmsy[Bench_Fmsy_mean_Sub$ObsNumForUpper]
		ItemsToLoop$Fmsy_LowerCI[i]=Sub$Fmsy[Bench_Fmsy_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_Fmsy = merge(Bench_Fmsy_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_Fmsy$ScenariosIncluded_Fmsy=NULL 
	Bench_Fmsy$FractionScenariosIncluded_Fmsy=NULL
	Bench_Fmsy$ObsNumForLower=NULL
	Bench_Fmsy$ObsNumForUpper=NULL	
	#Fmey
	Bench_Fmey = subset(Benchmarks,Benchmarks$Fmey_FLAG==0)  
	Bench_Fmey_mean = aggregate.data.frame(Bench_Fmey$Fmey,by=list(Bench_Fmey$species,Bench_Fmey$WPP),FUN=median,na.rm=TRUE)
	names(Bench_Fmey_mean) = c("species","WPP","Fmey_median")
	Samples = aggregate.data.frame(list(Benchmarks$Fmey_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_Fmey = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_Fmey = Samples$ScenariosIncluded_Fmey/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_Fmey_mean = merge(Bench_Fmey_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_Fmey_mean$ObsNumForLower = ceiling((Bench_Fmey_mean$ScenariosIncluded_Fmey*0.5) - (qnorm(0.99)*sqrt(Bench_Fmey_mean$ScenariosIncluded_Fmey*0.5)))
	Bench_Fmey_mean$ObsNumForUpper = ceiling((Bench_Fmey_mean$ScenariosIncluded_Fmey*0.5) + (qnorm(0.99)*sqrt(Bench_Fmey_mean$ScenariosIncluded_Fmey*0.5)))
	ItemsToLoop = data.frame(species=Bench_Fmey_mean$species,WPP=Bench_Fmey_mean$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$Fmey_UpperCI=NA
	ItemsToLoop$Fmey_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_Fmey,Bench_Fmey$species==ItemsToLoop$species[i] & Bench_Fmey$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$Fmey),]
		Bench_Fmey_mean_Sub = subset(Bench_Fmey_mean,Bench_Fmey_mean$species==ItemsToLoop$species[i] & Bench_Fmey_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$Fmey_UpperCI[i]=Sub$Fmey[Bench_Fmey_mean_Sub$ObsNumForUpper]
		ItemsToLoop$Fmey_LowerCI[i]=Sub$Fmey[Bench_Fmey_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_Fmey = merge(Bench_Fmey_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_Fmey$ScenariosIncluded_Fmey=NULL 
	Bench_Fmey$FractionScenariosIncluded_Fmey=NULL 
	Bench_Fmey$ObsNumForLower=NULL 
	Bench_Fmey$ObsNumForUpper=NULL
	#FratioSPR20
	Bench_Fratio20 = subset(Benchmarks,Benchmarks$Fratio20_FLAG==0)  
	Bench_Fratio20_mean = aggregate.data.frame(Bench_Fratio20$F_over_Fspr20,by=list(Bench_Fratio20$species,Bench_Fratio20$WPP),FUN=median,na.rm=TRUE)
	names(Bench_Fratio20_mean) = c("species","WPP","Fratio20_median")
	Samples = aggregate.data.frame(list(Benchmarks$Fratio20_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_Fratio20 = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_Fratio20 = Samples$ScenariosIncluded_Fratio20/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_Fratio20_mean = merge(Bench_Fratio20_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_Fratio20_mean$ObsNumForLower = ceiling((Bench_Fratio20_mean$ScenariosIncluded_Fratio20*0.5) - (qnorm(0.99)*sqrt(Bench_Fratio20_mean$ScenariosIncluded_Fratio20*0.5)))
	Bench_Fratio20_mean$ObsNumForUpper = ceiling((Bench_Fratio20_mean$ScenariosIncluded_Fratio20*0.5) + (qnorm(0.99)*sqrt(Bench_Fratio20_mean$ScenariosIncluded_Fratio20*0.5)))
	ItemsToLoop = data.frame(species=Bench_Fratio20_mean$species,WPP=Bench_Fratio20_mean$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$Fratio20_UpperCI=NA
	ItemsToLoop$Fratio20_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_Fratio20,Bench_Fratio20$species==ItemsToLoop$species[i] & Bench_Fratio20$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$F_over_Fspr20),]
		Bench_Fratio20_mean_Sub = subset(Bench_Fratio20_mean,Bench_Fratio20_mean$species==ItemsToLoop$species[i] & Bench_Fratio20_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$Fratio20_UpperCI[i]=Sub$F_over_Fspr20[Bench_Fratio20_mean_Sub$ObsNumForUpper]
		ItemsToLoop$Fratio20_LowerCI[i]=Sub$F_over_Fspr20[Bench_Fratio20_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_Fratio20 = merge(Bench_Fratio20_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_Fratio20$ScenariosIncluded_Fratio20=NULL 
	Bench_Fratio20$FractionScenariosIncluded_Fratio20=NULL 
	Bench_Fratio20$ObsNumForLower=NULL 
	Bench_Fratio20$ObsNumForUpper=NULL	
	#FratioSPR30
	Bench_Fratio30 = subset(Benchmarks,Benchmarks$Fratio30_FLAG==0)  
	Bench_Fratio30_mean = aggregate.data.frame(Bench_Fratio30$F_over_Fspr30,by=list(Bench_Fratio30$species,Bench_Fratio30$WPP),FUN=median,na.rm=TRUE)
	names(Bench_Fratio30_mean) = c("species","WPP","Fratio30_median")
	Samples = aggregate.data.frame(list(Benchmarks$Fratio30_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_Fratio30 = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_Fratio30 = Samples$ScenariosIncluded_Fratio30/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_Fratio30_mean = merge(Bench_Fratio30_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_Fratio30_mean$ObsNumForLower = ceiling((Bench_Fratio30_mean$ScenariosIncluded_Fratio30*0.5) - (qnorm(0.99)*sqrt(Bench_Fratio30_mean$ScenariosIncluded_Fratio30*0.5)))
	Bench_Fratio30_mean$ObsNumForUpper = ceiling((Bench_Fratio30_mean$ScenariosIncluded_Fratio30*0.5) + (qnorm(0.99)*sqrt(Bench_Fratio30_mean$ScenariosIncluded_Fratio30*0.5)))
	ItemsToLoop = data.frame(species=Bench_Fratio30_mean$species,WPP=Bench_Fratio30_mean$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$Fratio30_UpperCI=NA
	ItemsToLoop$Fratio30_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_Fratio30,Bench_Fratio30$species==ItemsToLoop$species[i] & Bench_Fratio30$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$F_over_Fspr30),]
		Bench_Fratio30_mean_Sub = subset(Bench_Fratio30_mean,Bench_Fratio30_mean$species==ItemsToLoop$species[i] & Bench_Fratio30_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$Fratio30_UpperCI[i]=Sub$F_over_Fspr30[Bench_Fratio30_mean_Sub$ObsNumForUpper]
		ItemsToLoop$Fratio30_LowerCI[i]=Sub$F_over_Fspr30[Bench_Fratio30_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_Fratio30 = merge(Bench_Fratio30_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_Fratio30$ScenariosIncluded_Fratio30=NULL 
	Bench_Fratio30$FractionScenariosIncluded_Fratio30=NULL 
	Bench_Fratio30$ObsNumForLower=NULL 
	Bench_Fratio30$ObsNumForUpper=NULL
	#FratioSPR40
	Bench_Fratio40 = subset(Benchmarks,Benchmarks$Fratio40_FLAG==0)  
	Bench_Fratio40_mean = aggregate.data.frame(Bench_Fratio40$F_over_Fspr40,by=list(Bench_Fratio40$species,Bench_Fratio40$WPP),FUN=median,na.rm=TRUE)
	names(Bench_Fratio40_mean) = c("species","WPP","Fratio40_median")
	Samples = aggregate.data.frame(list(Benchmarks$Fratio40_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_Fratio40 = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_Fratio40 = Samples$ScenariosIncluded_Fratio40/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_Fratio40_mean = merge(Bench_Fratio40_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_Fratio40_mean$ObsNumForLower = ceiling((Bench_Fratio40_mean$ScenariosIncluded_Fratio40*0.5) - (qnorm(0.99)*sqrt(Bench_Fratio40_mean$ScenariosIncluded_Fratio40*0.5)))
	Bench_Fratio40_mean$ObsNumForUpper = ceiling((Bench_Fratio40_mean$ScenariosIncluded_Fratio40*0.5) + (qnorm(0.99)*sqrt(Bench_Fratio40_mean$ScenariosIncluded_Fratio40*0.5)))
	ItemsToLoop = data.frame(species=Bench_Fratio40_mean$species,WPP=Bench_Fratio40_mean$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$Fratio40_UpperCI=NA
	ItemsToLoop$Fratio40_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_Fratio40,Bench_Fratio40$species==ItemsToLoop$species[i] & Bench_Fratio40$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$F_over_Fspr40),]
		Bench_Fratio40_mean_Sub = subset(Bench_Fratio40_mean,Bench_Fratio40_mean$species==ItemsToLoop$species[i] & Bench_Fratio40_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$Fratio40_UpperCI[i]=Sub$F_over_Fspr40[Bench_Fratio40_mean_Sub$ObsNumForUpper]
		ItemsToLoop$Fratio40_LowerCI[i]=Sub$F_over_Fspr40[Bench_Fratio40_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_Fratio40 = merge(Bench_Fratio40_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_Fratio40$ScenariosIncluded_Fratio40=NULL 
	Bench_Fratio40$FractionScenariosIncluded_Fratio40=NULL 
	Bench_Fratio40$ObsNumForLower=NULL 
	Bench_Fratio40$ObsNumForUpper=NULL	
	#Fmsy_ratio
	Bench_FratioMSY = subset(Benchmarks,Benchmarks$Fmsy_ratio_FLAG==0)  
	Bench_FratioMSY_mean = aggregate.data.frame(Bench_FratioMSY$F_over_Fmsy,by=list(Bench_FratioMSY$species,Bench_FratioMSY$WPP),FUN=median,na.rm=TRUE)
	names(Bench_FratioMSY_mean) = c("species","WPP","FratioMSY_median")
	Samples = aggregate.data.frame(list(Benchmarks$Fmsy_ratio_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_FratioMSY = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_FratioMSY = Samples$ScenariosIncluded_FratioMSY/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_FratioMSY_mean = merge(Bench_FratioMSY_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_FratioMSY_mean$ObsNumForLower = ceiling((Bench_FratioMSY_mean$ScenariosIncluded_FratioMSY*0.5) - (qnorm(0.99)*sqrt(Bench_FratioMSY_mean$ScenariosIncluded_FratioMSY*0.5)))
	Bench_FratioMSY_mean$ObsNumForUpper = ceiling((Bench_FratioMSY_mean$ScenariosIncluded_FratioMSY*0.5) + (qnorm(0.99)*sqrt(Bench_FratioMSY_mean$ScenariosIncluded_FratioMSY*0.5)))
	ItemsToLoop = data.frame(species=Bench_FratioMSY_mean$species,WPP=Bench_FratioMSY_mean$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$FratioMSY_UpperCI=NA
	ItemsToLoop$FratioMSY_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_FratioMSY,Bench_FratioMSY$species==ItemsToLoop$species[i] & Bench_FratioMSY$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$F_over_Fmsy),]
		Bench_FratioMSY_mean_Sub = subset(Bench_FratioMSY_mean,Bench_FratioMSY_mean$species==ItemsToLoop$species[i] & Bench_FratioMSY_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$FratioMSY_UpperCI[i]=Sub$F_over_Fmsy[Bench_FratioMSY_mean_Sub$ObsNumForUpper]
		ItemsToLoop$FratioMSY_LowerCI[i]=Sub$F_over_Fmsy[Bench_FratioMSY_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_FratioMSY = merge(Bench_FratioMSY_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_FratioMSY$ScenariosIncluded_FratioMSY=NULL 
	Bench_FratioMSY$FractionScenariosIncluded_FratioMSY=NULL 
	Bench_FratioMSY$ObsNumForLower=NULL 
	Bench_FratioMSY$ObsNumForUpper=NULL		
	#Fmey_Ratio
	Bench_FratioMEY = subset(Benchmarks,Benchmarks$Fmey_ratio_FLAG==0)  
	Bench_FratioMEY_mean = aggregate.data.frame(Bench_FratioMEY$F_over_Fmey,by=list(Bench_FratioMEY$species,Bench_FratioMEY$WPP),FUN=median,na.rm=TRUE)
	names(Bench_FratioMEY_mean) = c("species","WPP","FratioMEY_median")
	Samples = aggregate.data.frame(list(Benchmarks$Fmey_ratio_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_FratioMEY = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_FratioMEY = Samples$ScenariosIncluded_FratioMEY/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_FratioMEY_mean = merge(Bench_FratioMEY_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_FratioMEY_mean$ObsNumForLower = ceiling((Bench_FratioMEY_mean$ScenariosIncluded_FratioMEY*0.5) - (qnorm(0.99)*sqrt(Bench_FratioMEY_mean$ScenariosIncluded_FratioMEY*0.5)))
	Bench_FratioMEY_mean$ObsNumForUpper = ceiling((Bench_FratioMEY_mean$ScenariosIncluded_FratioMEY*0.5) + (qnorm(0.99)*sqrt(Bench_FratioMEY_mean$ScenariosIncluded_FratioMEY*0.5)))
	ItemsToLoop = data.frame(species=Bench_FratioMEY_mean$species,WPP=Bench_FratioMEY_mean$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$FratioMEY_UpperCI=NA
	ItemsToLoop$FratioMEY_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_FratioMEY,Bench_FratioMEY$species==ItemsToLoop$species[i] & Bench_FratioMEY$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$F_over_Fmey),]
		Bench_FratioMEY_mean_Sub = subset(Bench_FratioMEY_mean,Bench_FratioMEY_mean$species==ItemsToLoop$species[i] & Bench_FratioMEY_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$FratioMEY_UpperCI[i]=Sub$F_over_Fmey[Bench_FratioMEY_mean_Sub$ObsNumForUpper]
		ItemsToLoop$FratioMEY_LowerCI[i]=Sub$F_over_Fmey[Bench_FratioMEY_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_FratioMEY = merge(Bench_FratioMEY_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_FratioMEY$ScenariosIncluded_FratioMEY=NULL 
	Bench_FratioMEY$FractionScenariosIncluded_FratioMEY=NULL 
	Bench_FratioMEY$ObsNumForLower=NULL 
	Bench_FratioMEY$ObsNumForUpper=NULL	
	#Remove Outliers: Biomass
	Benchmarks$B_current_FLAG=0
	Benchmarks$Bspr20_FLAG=0
	Benchmarks$Bspr30_FLAG=0
	Benchmarks$Bspr40_FLAG=0
	Benchmarks$Bmsy_FLAG=0
	Benchmarks$Bmey_FLAG=0
	Benchmarks$Bratio20_FLAG=0
	Benchmarks$Bratio30_FLAG=0
	Benchmarks$Bratio40_FLAG=0
	Benchmarks$BratioMSY_FLAG=0
	Benchmarks$BratioMEY_FLAG=0
	Benchmarks$Bratio20_FLAG[Benchmarks$B_over_Bspr20 > 5 | Benchmarks$B_over_Bspr20 < 0.1]=1
	Benchmarks$Bratio20_FLAG[is.infinite(Benchmarks$B_over_Bspr20)]=1
	Benchmarks$Bratio30_FLAG[Benchmarks$B_over_Bspr30 > 5 | Benchmarks$B_over_Bspr30 < 0.1]=1
	Benchmarks$Bratio30_FLAG[is.infinite(Benchmarks$B_over_Bspr30)]=1
	Benchmarks$Bratio40_FLAG[Benchmarks$B_over_Bspr40 > 5 | Benchmarks$B_over_Bspr40 < 0.1]=1
	Benchmarks$Bratio40_FLAG[is.infinite(Benchmarks$B_over_Bspr40)]=1	
	Benchmarks$BratioMSY_FLAG[Benchmarks$B_over_Bmsy > 5 | Benchmarks$B_over_Bmsy < 0.1]=1
	Benchmarks$BratioMSY_FLAG[is.infinite(Benchmarks$B_over_Bmsy)]=1	
	Benchmarks$BratioMEY_FLAG[Benchmarks$B_over_Bmey > 5 | Benchmarks$B_over_Bmey < 0.1]=1
	Benchmarks$BratioMEY_FLAG[is.infinite(Benchmarks$B_over_Bmey)]=1
	#Set biomass tolerance limits to identify outliers in estimated biomass
	MinQuant_B = 0.025
	MaxQuant_B = 0.975
	Spp = sort(unique(Benchmarks$species))
	Benchmarks$B_current_ln = log(Benchmarks$B_current)
	Benchmarks$Bspr20_ln = log(Benchmarks$B_spr20_kg)
	Benchmarks$Bspr30_ln = log(Benchmarks$B_spr30_kg)
	Benchmarks$Bspr40_ln = log(Benchmarks$B_spr40_kg)
	Benchmarks$Bmsy_ln = log(Benchmarks$Bmsy_kg)
	Benchmarks$Bmey_ln = log(Benchmarks$B_MEY_kg)
	Benchmarks$B_over_Bspr20_ln = log(Benchmarks$B_over_Bspr20)
	Benchmarks$B_over_Bspr30_ln = log(Benchmarks$B_over_Bspr30)
	Benchmarks$B_over_Bspr40_ln = log(Benchmarks$B_over_Bspr40)
	Benchmarks$B_over_Bmsy_ln = log(Benchmarks$B_over_Bmsy)
	Benchmarks$B_over_Bmey_ln = log(Benchmarks$B_over_Bmey)
	for(i in 1:dim(SppWpp)[1])
	{
#		#Determine Cut-off values from distribution for B_current
#		thisSpp = Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]
#		Hist = hist(thisSpp$B_current,plot=FALSE)
#		HistDat = data.frame(Mids=Hist$mids,Counts=Hist$counts)
#		HistDat_sub = subset(HistDat,HistDat$Counts==0)
#		Median = median(thisSpp$B_current)
#		cutOff = HistDat_sub$Mids[1]
#		while(dim(HistDat_sub)[1]>0)
#		{
#			if(cutOff > Median)
#			{
#				Benchmarks$B_current_FLAG[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i] & Benchmarks$B_current > cutOff]=1
#			}
#			if(cutOff < Median)
#			{
#				thisSpp = subset(thisSpp,thisSpp$B_current>cutOff)
#				Benchmarks$B_current_FLAG[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i] & Benchmarks$B_current < cutOff]=1
#			}
#			Hist = hist(thisSpp$B_current,plot=FALSE)
#			HistDat = data.frame(Mids=Hist$mids,Counts=Hist$counts)
#			HistDat_sub = subset(HistDat,HistDat$Counts==0)
#			Median = median(thisSpp$B_current)
#			cutOff = HistDat_sub$Mids[1]
#		}
		
		
		#B_estiamted: Use Quantiles now to refine the natrual log of the distribution
#		thisSpp = Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i] & Benchmarks$B_current_FLAG==0,]	
#		Quant = quantile(thisSpp$B_current_ln,probs=c(MinQuant_B,MaxQuant_B),na.rm=TRUE,type=8)
#		Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$B_current_FLAG[Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$B_current_ln < Quant[1] | Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$B_current_ln > Quant[2]]=1
		#BratioSPR20
		thisSpp = Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i] & Benchmarks$Bratio20_FLAG==0,]	
		Quant = quantile(thisSpp$B_over_Bspr20_ln,probs=c(MinQuant_B,MaxQuant_B),na.rm=TRUE,type=8)
		Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Bspr20_FLAG[Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$B_over_Bspr20_ln < Quant[1] | Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$B_over_Bspr20_ln > Quant[2]]=1
		#BratioSPR30
		thisSpp = Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i] & Benchmarks$Bratio30_FLAG==0,]	
		Quant = quantile(thisSpp$B_over_Bspr30_ln,probs=c(MinQuant_B,MaxQuant_B),na.rm=TRUE,type=8)
		Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Bspr30_FLAG[Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$B_over_Bspr30_ln < Quant[1] | Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$B_over_Bspr30_ln > Quant[2]]=1
		#BratioSPR40
		thisSpp = Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i] & Benchmarks$Bratio40_FLAG==0,]	
		Quant = quantile(thisSpp$B_over_Bspr40_ln,probs=c(MinQuant_B,MaxQuant_B),na.rm=TRUE,type=8)
		Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$Bspr40_FLAG[Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$B_over_Bspr40_ln < Quant[1] | Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$B_over_Bspr40_ln > Quant[2]]=1
		#BratioMSY
		thisSpp = Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i] & Benchmarks$BratioMSY_FLAG==0,]	
		Quant = quantile(thisSpp$B_over_Bmsy_ln,probs=c(MinQuant_B,MaxQuant_B),na.rm=TRUE,type=8)
		Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$BratioMSY_FLAG[Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$B_over_Bmsy_ln < Quant[1] | Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$B_over_Bmsy_ln > Quant[2]]=1
		#BratioMEY
		thisSpp = Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i] & Benchmarks$BratioMEY_FLAG==0,]	
		Quant = quantile(thisSpp$B_over_Bmey_ln,probs=c(MinQuant_B,MaxQuant_B),na.rm=TRUE,type=8)
		Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$BratioMEY_FLAG[Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$B_over_Bmsy_ln < Quant[1] | Benchmarks[Benchmarks$species==SppWpp$species[i] & Benchmarks$WPP==SppWpp$WPP[i],]$B_over_Bmsy_ln > Quant[2]]=1
		if(i%%10==0)
		{
			print(paste("Working on item ",i," of ",dim(SppWpp)[1],"...",sep=""))
			flush.console()	
		}
	}
	write.table(Benchmarks,paste(PATH_output,"Benchmarks_withFlags.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)

	#Bcurrent
	#Bench_B_current = subset(Benchmarks,Benchmarks$B_current_FLAG==0) 
	Bench_B_current = Benchmarks   #Do NOT filter any observations if using the Catch-MSY calculated biomass!
	Bench_B_current_mean = aggregate.data.frame(Bench_B_current$B_current,by=list(Bench_B_current$species,Bench_B_current$WPP),FUN=geoMean)  #geoMean
	names(Bench_B_current_mean) = c("species","WPP","B_current_geoMean")
	Bench_B_current_stdev = aggregate.data.frame(Bench_B_current$B_current,by=list(Bench_B_current$species,Bench_B_current$WPP),FUN=geoStDev)  #geoStDev
	names(Bench_B_current_stdev) = c("species","WPP","B_current_geoStdev")
	Bench_B_current = merge(Bench_B_current_mean,Bench_B_current_stdev,all.x=TRUE,by=c("species","WPP"))
	Samples = aggregate.data.frame(list(Benchmarks$B_current_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_B_current = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_B_current = Samples$ScenariosIncluded_B_current/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_B_current = merge(Bench_B_current,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_B_current$B_current_UpperCI = Bench_B_current$B_current_geoMean + ((Bench_B_current$B_current_geoStdev/sqrt(Bench_B_current$ScenariosIncluded_B_current))*1.96)
	Bench_B_current$B_current_LowerCI = Bench_B_current$B_current_geoMean - ((Bench_B_current$B_current_geoStdev/sqrt(Bench_B_current$ScenariosIncluded_B_current))*1.96)
	Bench_B_current$ScenariosIncluded_B_current=NULL 
	Bench_B_current$FractionScenariosIncluded_B_current=NULL 
	#SPRcurrent
	Bench_SPR_current = subset(Benchmarks,Benchmarks$F_current_FLAG==0)
	Bench_mean_SPR_current = aggregate.data.frame(Bench_SPR_current$SPR_current,by=list(Bench_SPR_current$species,Bench_SPR_current$WPP),FUN=median,na.rm=TRUE)
	names(Bench_mean_SPR_current) = c("species","WPP","SPR_median")
	Samples = aggregate.data.frame(list(Benchmarks$F_current_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_SPR_current = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_SPR_current = Samples$ScenariosIncluded_SPR_current/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_mean_SPR_current = merge(Bench_mean_SPR_current,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_mean_SPR_current$ObsNumForLower = ceiling((Bench_mean_SPR_current$ScenariosIncluded_SPR_current*0.5) - (qnorm(0.99)*sqrt(Bench_mean_SPR_current$ScenariosIncluded_SPR_current*0.5)))
	Bench_mean_SPR_current$ObsNumForUpper = ceiling((Bench_mean_SPR_current$ScenariosIncluded_SPR_current*0.5) + (qnorm(0.99)*sqrt(Bench_mean_SPR_current$ScenariosIncluded_SPR_current*0.5)))
	ItemsToLoop = data.frame(species=Bench_SPR_current$species,WPP=Bench_SPR_current$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$SPR_current_UpperCI=NA
	ItemsToLoop$SPR_current_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_SPR_current,Bench_SPR_current$species==ItemsToLoop$species[i] & Bench_SPR_current$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$SPR_current),]
		Bench_SPR_current_mean_Sub = subset(Bench_mean_SPR_current,Bench_mean_SPR_current$species==ItemsToLoop$species[i] & Bench_mean_SPR_current$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$SPR_current_UpperCI[i]=Sub$SPR_current[Bench_SPR_current_mean_Sub$ObsNumForUpper]
		ItemsToLoop$SPR_current_LowerCI[i]=Sub$SPR_current[Bench_SPR_current_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_SPR_current = merge(Bench_mean_SPR_current,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_SPR_current$ScenariosIncluded_SPR_current=NULL 
	Bench_SPR_current$FractionScenariosIncluded_SPR_current=NULL 
	Bench_SPR_current$ObsNumForLower=NULL 
	Bench_SPR_current$ObsNumForUpper=NULL



	#Yield_current
	Bench_Yield_current = Benchmarks   #Do NOT filter any observations if using the Catch-MSY calculated biomass!
	Bench_Yield_current_mean = aggregate.data.frame(Bench_Yield_current$Yield_current,by=list(Bench_Yield_current$species,Bench_Yield_current$WPP),FUN=geoMean)  #geoMean
	names(Bench_Yield_current_mean) = c("species","WPP","Yield_current_geoMean")
	Bench_Yield_current_stdev = aggregate.data.frame(Bench_Yield_current$Yield_current,by=list(Bench_Yield_current$species,Bench_Yield_current$WPP),FUN=geoStDev)  #geoStDev
	names(Bench_Yield_current_stdev) = c("species","WPP","Yield_current_geoStdev")
	Bench_Yield_current = merge(Bench_Yield_current_mean,Bench_Yield_current_stdev,all.x=TRUE,by=c("species","WPP"))
	Samples = aggregate.data.frame(list(Benchmarks$B_current_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_B_current = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_B_current = Samples$ScenariosIncluded_B_current/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_Yield_current = merge(Bench_Yield_current,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_Yield_current$Yield_current_UpperCI = Bench_Yield_current$Yield_current_geoMean + ((Bench_Yield_current$Yield_current_geoStdev/sqrt(Bench_Yield_current$ScenariosIncluded_B_current))*1.96)
	Bench_Yield_current$Yield_current_LowerCI = Bench_Yield_current$Yield_current_geoMean - ((Bench_Yield_current$Yield_current_geoStdev/sqrt(Bench_Yield_current$ScenariosIncluded_B_current))*1.96)
	Bench_Yield_current$ScenariosIncluded_B_current=NULL 
	Bench_Yield_current$FractionScenariosIncluded_B_current=NULL 
	#BratioSPR20
	Bench_Bratio20 = subset(Benchmarks,Benchmarks$Bratio20_FLAG==0)
	Bench_Bratio20_mean = aggregate.data.frame(Bench_Bratio20$B_over_Bspr20,by=list(Bench_Bratio20$species,Bench_Bratio20$WPP),FUN=median,na.rm=TRUE)   #geoMean   MedianAbsoluteDeviation
	names(Bench_Bratio20_mean) = c("species","WPP","B_over_Bspr20_median")
	Samples = aggregate.data.frame(list(Benchmarks$Bratio20_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_Bratio20_current = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_Bratio20_current = Samples$ScenariosIncluded_Bratio20_current/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_Bratio20_mean = merge(Bench_Bratio20_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_Bratio20_mean$ObsNumForLower = ceiling((Bench_Bratio20_mean$ScenariosIncluded_Bratio20_current*0.5) - (qnorm(0.99)*sqrt(Bench_Bratio20_mean$ScenariosIncluded_Bratio20_current*0.5)))
	Bench_Bratio20_mean$ObsNumForUpper = ceiling((Bench_Bratio20_mean$ScenariosIncluded_Bratio20_current*0.5) + (qnorm(0.99)*sqrt(Bench_Bratio20_mean$ScenariosIncluded_Bratio20_current*0.5)))
	ItemsToLoop = data.frame(species=Bench_Bratio20$species,WPP=Bench_Bratio20$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$Bratio20_UpperCI=NA
	ItemsToLoop$Bratio20_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_Bratio20,Bench_Bratio20$species==ItemsToLoop$species[i] & Bench_Bratio20$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$B_over_Bspr20),]
		Bench_mean_Bratio20_current_Sub = subset(Bench_Bratio20_mean,Bench_Bratio20_mean$species==ItemsToLoop$species[i] & Bench_Bratio20_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$Bratio20_UpperCI[i]=Sub$B_over_Bspr20[Bench_mean_Bratio20_current_Sub$ObsNumForUpper]
		ItemsToLoop$Bratio20_LowerCI[i]=Sub$B_over_Bspr20[Bench_mean_Bratio20_current_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_Bratio20 = merge(Bench_Bratio20_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_Bratio20$ScenariosIncluded_Bratio20_current=NULL 
	Bench_Bratio20$FractionScenariosIncluded_Bratio20_current=NULL 
	Bench_Bratio20$ObsNumForLower=NULL 
	Bench_Bratio20$ObsNumForUpper=NULL




	#Biomass at SPR 20
	Bench_BiomassSpr20 = subset(Benchmarks,Benchmarks$Bspr20_FLAG==0)
	Bench_BiomassSpr20_mean = aggregate.data.frame(Bench_BiomassSpr20$B_spr20_kg,by=list(Bench_BiomassSpr20$species,Bench_BiomassSpr20$WPP),FUN=median,na.rm=TRUE)
	names(Bench_BiomassSpr20_mean) = c("species","WPP","Bspr20_median")
	Samples = aggregate.data.frame(list(Benchmarks$Bspr20_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_Bspr20_current = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_Bspr20_current = Samples$ScenariosIncluded_Bspr20_current/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_BiomassSpr20_mean = merge(Bench_BiomassSpr20_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_BiomassSpr20_mean$ObsNumForLower = ceiling((Bench_BiomassSpr20_mean$ScenariosIncluded_Bspr20_current*0.5) - (qnorm(0.99)*sqrt(Bench_BiomassSpr20_mean$ScenariosIncluded_Bspr20_current*0.5)))
	Bench_BiomassSpr20_mean$ObsNumForUpper = ceiling((Bench_BiomassSpr20_mean$ScenariosIncluded_Bspr20_current*0.5) + (qnorm(0.99)*sqrt(Bench_BiomassSpr20_mean$ScenariosIncluded_Bspr20_current*0.5)))
	ItemsToLoop = data.frame(species=Bench_BiomassSpr20$species,WPP=Bench_BiomassSpr20$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$BiomassSpr20_UpperCI=NA
	ItemsToLoop$BiomassSpr20_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_BiomassSpr20,Bench_BiomassSpr20$species==ItemsToLoop$species[i] & Bench_BiomassSpr20$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$B_spr20_kg),]
		Bench_BiomassSpr20_mean_Sub = subset(Bench_BiomassSpr20_mean,Bench_BiomassSpr20_mean$species==ItemsToLoop$species[i] & Bench_BiomassSpr20_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$BiomassSpr20_UpperCI[i]=Sub$B_spr20_kg[Bench_BiomassSpr20_mean_Sub$ObsNumForUpper]
		ItemsToLoop$BiomassSpr20_LowerCI[i]=Sub$B_spr20_kg[Bench_BiomassSpr20_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_BiomassSpr20 = merge(Bench_BiomassSpr20_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_BiomassSpr20$ScenariosIncluded_Bspr20_current=NULL 
	Bench_BiomassSpr20$FractionScenariosIncluded_Bspr20_current=NULL 
	Bench_BiomassSpr20$ObsNumForLower=NULL 
	Bench_BiomassSpr20$ObsNumForUpper=NULL	
	#Yield at SPR 20
	Bench_YieldSpr20 = subset(Benchmarks,Benchmarks$Bspr20_FLAG==0)
	Bench_YieldSpr20_mean = aggregate.data.frame(Bench_YieldSpr20$Yield_spr20_kg,by=list(Bench_YieldSpr20$species,Bench_YieldSpr20$WPP),FUN=median,na.rm=TRUE)
	names(Bench_YieldSpr20_mean) = c("species","WPP","YieldSpr20_median")
	Samples = aggregate.data.frame(list(Benchmarks$Bspr20_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_YieldSpr20 = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_YieldSpr20 = Samples$ScenariosIncluded_YieldSpr20/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_YieldSpr20_mean = merge(Bench_YieldSpr20_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_YieldSpr20_mean$ObsNumForLower = ceiling((Bench_YieldSpr20_mean$ScenariosIncluded_YieldSpr20*0.5) - (qnorm(0.99)*sqrt(Bench_YieldSpr20_mean$ScenariosIncluded_YieldSpr20*0.5)))
	Bench_YieldSpr20_mean$ObsNumForUpper = ceiling((Bench_YieldSpr20_mean$ScenariosIncluded_YieldSpr20*0.5) + (qnorm(0.99)*sqrt(Bench_YieldSpr20_mean$ScenariosIncluded_YieldSpr20*0.5)))
	ItemsToLoop = data.frame(species=Bench_YieldSpr20$species,WPP=Bench_YieldSpr20$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$YieldSpr20_UpperCI=NA
	ItemsToLoop$YieldSpr20_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_YieldSpr20,Bench_YieldSpr20$species==ItemsToLoop$species[i] & Bench_YieldSpr20$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$Yield_spr20_kg),]
		Bench_YieldSpr20_mean_Sub = subset(Bench_YieldSpr20_mean,Bench_YieldSpr20_mean$species==ItemsToLoop$species[i] & Bench_YieldSpr20_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$YieldSpr20_UpperCI[i]=Sub$Yield_spr20_kg[Bench_YieldSpr20_mean_Sub$ObsNumForUpper]
		ItemsToLoop$YieldSpr20_LowerCI[i]=Sub$Yield_spr20_kg[Bench_YieldSpr20_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_YieldSpr20 = merge(Bench_YieldSpr20_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_YieldSpr20$ScenariosIncluded_YieldSpr20=NULL 
	Bench_YieldSpr20$FractionScenariosIncluded_YieldSpr20=NULL 
	Bench_YieldSpr20$ObsNumForLower=NULL 
	Bench_YieldSpr20$ObsNumForUpper=NULL	
	#Rebuiling Time to spr20
	Bench_RebuildSPR20 = subset(Benchmarks,Benchmarks$F_current_FLAG==0 & Benchmarks$Fratio20_FLAG==0)
	Bench_mean_RebuildSPR20 = aggregate.data.frame(Bench_RebuildSPR20$RebuildTime_spr20,by=list(Bench_RebuildSPR20$species,Bench_RebuildSPR20$WPP),FUN=median,na.rm=TRUE)
	names(Bench_mean_RebuildSPR20) = c("species","WPP","RebuildToSPR20_median")
	Samples = aggregate.data.frame(list(Benchmarks$Fratio20_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_RebuildSPR20 = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_RebuildSPR20 = Samples$ScenariosIncluded_RebuildSPR20/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_mean_RebuildSPR20 = merge(Bench_mean_RebuildSPR20,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_mean_RebuildSPR20$ObsNumForLower = ceiling((Bench_mean_RebuildSPR20$ScenariosIncluded_RebuildSPR20*0.5) - (qnorm(0.99)*sqrt(Bench_mean_RebuildSPR20$ScenariosIncluded_RebuildSPR20*0.5)))
	Bench_mean_RebuildSPR20$ObsNumForUpper = ceiling((Bench_mean_RebuildSPR20$ScenariosIncluded_RebuildSPR20*0.5) + (qnorm(0.99)*sqrt(Bench_mean_RebuildSPR20$ScenariosIncluded_RebuildSPR20*0.5)))
	ItemsToLoop = data.frame(species=Bench_RebuildSPR20$species,WPP=Bench_RebuildSPR20$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$RebuildSPR20_UpperCI=NA
	ItemsToLoop$RebuildSPR20_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_RebuildSPR20,Bench_RebuildSPR20$species==ItemsToLoop$species[i] & Bench_RebuildSPR20$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$RebuildTime_spr20),]
		Bench_mean_RebuildSPR20_Sub = subset(Bench_mean_RebuildSPR20,Bench_mean_RebuildSPR20$species==ItemsToLoop$species[i] & Bench_mean_RebuildSPR20$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$RebuildSPR20_UpperCI[i]=Sub$RebuildTime_spr20[Bench_mean_RebuildSPR20_Sub$ObsNumForUpper]
		ItemsToLoop$RebuildSPR20_LowerCI[i]=Sub$RebuildTime_spr20[Bench_mean_RebuildSPR20_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_RebuildSPR20 = merge(Bench_mean_RebuildSPR20,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_RebuildSPR20$ScenariosIncluded_RebuildSPR20=NULL 
	Bench_RebuildSPR20$FractionScenariosIncluded_RebuildSPR20=NULL 
	Bench_RebuildSPR20$ObsNumForLower=NULL 
	Bench_RebuildSPR20$ObsNumForUpper=NULL
	#BratioSPR30
	Bench_Bratio30 = subset(Benchmarks,Benchmarks$Bratio30_FLAG==0)
	Bench_Bratio30_mean = aggregate.data.frame(Bench_Bratio30$B_over_Bspr30,by=list(Bench_Bratio30$species,Bench_Bratio30$WPP),FUN=median,na.rm=TRUE)   #geoMean   MedianAbsoluteDeviation
	names(Bench_Bratio30_mean) = c("species","WPP","B_over_Bspr30_median")
	Samples = aggregate.data.frame(list(Benchmarks$Bratio30_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_Bratio30_current = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_Bratio30_current = Samples$ScenariosIncluded_Bratio30_current/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_Bratio30_mean = merge(Bench_Bratio30_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_Bratio30_mean$ObsNumForLower = ceiling((Bench_Bratio30_mean$ScenariosIncluded_Bratio30_current*0.5) - (qnorm(0.99)*sqrt(Bench_Bratio30_mean$ScenariosIncluded_Bratio30_current*0.5)))
	Bench_Bratio30_mean$ObsNumForUpper = ceiling((Bench_Bratio30_mean$ScenariosIncluded_Bratio30_current*0.5) + (qnorm(0.99)*sqrt(Bench_Bratio30_mean$ScenariosIncluded_Bratio30_current*0.5)))
	ItemsToLoop = data.frame(species=Bench_Bratio30$species,WPP=Bench_Bratio30$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$Bratio30_UpperCI=NA
	ItemsToLoop$Bratio30_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_Bratio30,Bench_Bratio30$species==ItemsToLoop$species[i] & Bench_Bratio30$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$B_over_Bspr30),]
		Bench_mean_Bratio30_current_Sub = subset(Bench_Bratio30_mean,Bench_Bratio30_mean$species==ItemsToLoop$species[i] & Bench_Bratio30_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$Bratio30_UpperCI[i]=Sub$B_over_Bspr30[Bench_mean_Bratio30_current_Sub$ObsNumForUpper]
		ItemsToLoop$Bratio30_LowerCI[i]=Sub$B_over_Bspr30[Bench_mean_Bratio30_current_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_Bratio30 = merge(Bench_Bratio30_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_Bratio30$ScenariosIncluded_Bratio30_current=NULL 
	Bench_Bratio30$FractionScenariosIncluded_Bratio30_current=NULL 
	Bench_Bratio30$ObsNumForLower=NULL 
	Bench_Bratio30$ObsNumForUpper=NULL
	#Biomass at SPR 30
	Bench_BiomassSpr30 = subset(Benchmarks,Benchmarks$Bspr30_FLAG==0)
	Bench_BiomassSpr30_mean = aggregate.data.frame(Bench_BiomassSpr30$B_spr30_kg,by=list(Bench_BiomassSpr30$species,Bench_BiomassSpr30$WPP),FUN=median,na.rm=TRUE)
	names(Bench_BiomassSpr30_mean) = c("species","WPP","Bspr30_median")
	Samples = aggregate.data.frame(list(Benchmarks$Bspr30_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_Bspr30_current = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_Bspr30_current = Samples$ScenariosIncluded_Bspr30_current/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_BiomassSpr30_mean = merge(Bench_BiomassSpr30_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_BiomassSpr30_mean$ObsNumForLower = ceiling((Bench_BiomassSpr30_mean$ScenariosIncluded_Bspr30_current*0.5) - (qnorm(0.99)*sqrt(Bench_BiomassSpr30_mean$ScenariosIncluded_Bspr30_current*0.5)))
	Bench_BiomassSpr30_mean$ObsNumForUpper = ceiling((Bench_BiomassSpr30_mean$ScenariosIncluded_Bspr30_current*0.5) + (qnorm(0.99)*sqrt(Bench_BiomassSpr30_mean$ScenariosIncluded_Bspr30_current*0.5)))
	ItemsToLoop = data.frame(species=Bench_BiomassSpr30$species,WPP=Bench_BiomassSpr30$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$BiomassSpr30_UpperCI=NA
	ItemsToLoop$BiomassSpr30_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_BiomassSpr30,Bench_BiomassSpr30$species==ItemsToLoop$species[i] & Bench_BiomassSpr30$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$B_spr30_kg),]
		Bench_BiomassSpr30_mean_Sub = subset(Bench_BiomassSpr30_mean,Bench_BiomassSpr30_mean$species==ItemsToLoop$species[i] & Bench_BiomassSpr30_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$BiomassSpr30_UpperCI[i]=Sub$B_spr30_kg[Bench_BiomassSpr30_mean_Sub$ObsNumForUpper]
		ItemsToLoop$BiomassSpr30_LowerCI[i]=Sub$B_spr30_kg[Bench_BiomassSpr30_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_BiomassSpr30 = merge(Bench_BiomassSpr30_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_BiomassSpr30$ScenariosIncluded_Bspr30_current=NULL 
	Bench_BiomassSpr30$FractionScenariosIncluded_Bspr30_current=NULL 
	Bench_BiomassSpr30$ObsNumForLower=NULL 
	Bench_BiomassSpr30$ObsNumForUpper=NULL	
	#Yield at SPR 30
	Bench_YieldSpr30 = subset(Benchmarks,Benchmarks$Bspr30_FLAG==0)
	Bench_YieldSpr30_mean = aggregate.data.frame(Bench_YieldSpr30$Yield_spr30_kg,by=list(Bench_YieldSpr30$species,Bench_YieldSpr30$WPP),FUN=median,na.rm=TRUE)
	names(Bench_YieldSpr30_mean) = c("species","WPP","YieldSpr30_median")
	Samples = aggregate.data.frame(list(Benchmarks$Bspr30_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_YieldSpr30 = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_YieldSpr30 = Samples$ScenariosIncluded_YieldSpr30/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_YieldSpr30_mean = merge(Bench_YieldSpr30_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_YieldSpr30_mean$ObsNumForLower = ceiling((Bench_YieldSpr30_mean$ScenariosIncluded_YieldSpr30*0.5) - (qnorm(0.99)*sqrt(Bench_YieldSpr30_mean$ScenariosIncluded_YieldSpr30*0.5)))
	Bench_YieldSpr30_mean$ObsNumForUpper = ceiling((Bench_YieldSpr30_mean$ScenariosIncluded_YieldSpr30*0.5) + (qnorm(0.99)*sqrt(Bench_YieldSpr30_mean$ScenariosIncluded_YieldSpr30*0.5)))
	ItemsToLoop = data.frame(species=Bench_YieldSpr30$species,WPP=Bench_YieldSpr30$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$YieldSpr30_UpperCI=NA
	ItemsToLoop$YieldSpr30_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_YieldSpr30,Bench_YieldSpr30$species==ItemsToLoop$species[i] & Bench_YieldSpr30$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$Yield_spr30_kg),]
		Bench_YieldSpr30_mean_Sub = subset(Bench_YieldSpr30_mean,Bench_YieldSpr30_mean$species==ItemsToLoop$species[i] & Bench_YieldSpr30_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$YieldSpr30_UpperCI[i]=Sub$Yield_spr30_kg[Bench_YieldSpr30_mean_Sub$ObsNumForUpper]
		ItemsToLoop$YieldSpr30_LowerCI[i]=Sub$Yield_spr30_kg[Bench_YieldSpr30_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_YieldSpr30 = merge(Bench_YieldSpr30_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_YieldSpr30$ScenariosIncluded_YieldSpr30=NULL 
	Bench_YieldSpr30$FractionScenariosIncluded_YieldSpr30=NULL 
	Bench_YieldSpr30$ObsNumForLower=NULL 
	Bench_YieldSpr30$ObsNumForUpper=NULL	
	#Rebuiling Time to spr30
	Bench_RebuildSPR30 = subset(Benchmarks,Benchmarks$F_current_FLAG==0 & Benchmarks$Fratio30_FLAG==0)
	Bench_mean_RebuildSPR30 = aggregate.data.frame(Bench_RebuildSPR30$RebuildTime_spr30,by=list(Bench_RebuildSPR30$species,Bench_RebuildSPR30$WPP),FUN=median,na.rm=TRUE)
	names(Bench_mean_RebuildSPR30) = c("species","WPP","RebuildToSPR30_median")
	Samples = aggregate.data.frame(list(Benchmarks$Fratio30_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_RebuildSPR30 = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_RebuildSPR30 = Samples$ScenariosIncluded_RebuildSPR30/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_mean_RebuildSPR30 = merge(Bench_mean_RebuildSPR30,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_mean_RebuildSPR30$ObsNumForLower = ceiling((Bench_mean_RebuildSPR30$ScenariosIncluded_RebuildSPR30*0.5) - (qnorm(0.99)*sqrt(Bench_mean_RebuildSPR30$ScenariosIncluded_RebuildSPR30*0.5)))
	Bench_mean_RebuildSPR30$ObsNumForUpper = ceiling((Bench_mean_RebuildSPR30$ScenariosIncluded_RebuildSPR30*0.5) + (qnorm(0.99)*sqrt(Bench_mean_RebuildSPR30$ScenariosIncluded_RebuildSPR30*0.5)))
	ItemsToLoop = data.frame(species=Bench_RebuildSPR30$species,WPP=Bench_RebuildSPR30$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$RebuildSPR30_UpperCI=NA
	ItemsToLoop$RebuildSPR30_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_RebuildSPR30,Bench_RebuildSPR30$species==ItemsToLoop$species[i] & Bench_RebuildSPR30$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$RebuildTime_spr30),]
		Bench_mean_RebuildSPR30_Sub = subset(Bench_mean_RebuildSPR30,Bench_mean_RebuildSPR30$species==ItemsToLoop$species[i] & Bench_mean_RebuildSPR30$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$RebuildSPR30_UpperCI[i]=Sub$RebuildTime_spr30[Bench_mean_RebuildSPR30_Sub$ObsNumForUpper]
		ItemsToLoop$RebuildSPR30_LowerCI[i]=Sub$RebuildTime_spr30[Bench_mean_RebuildSPR30_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_RebuildSPR30 = merge(Bench_mean_RebuildSPR30,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_RebuildSPR30$ScenariosIncluded_RebuildSPR30=NULL 
	Bench_RebuildSPR30$FractionScenariosIncluded_RebuildSPR30=NULL 
	Bench_RebuildSPR30$ObsNumForLower=NULL 
	Bench_RebuildSPR30$ObsNumForUpper=NULL
	#BratioSPR40
	Bench_Bratio40 = subset(Benchmarks,Benchmarks$Bratio40_FLAG==0)
	Bench_Bratio40_mean = aggregate.data.frame(Bench_Bratio40$B_over_Bspr40,by=list(Bench_Bratio40$species,Bench_Bratio40$WPP),FUN=median,na.rm=TRUE)   #geoMean   MedianAbsoluteDeviation
	names(Bench_Bratio40_mean) = c("species","WPP","B_over_Bspr40_median")
	Samples = aggregate.data.frame(list(Benchmarks$Bratio40_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_Bratio40_current = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_Bratio40_current = Samples$ScenariosIncluded_Bratio40_current/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_Bratio40_mean = merge(Bench_Bratio40_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_Bratio40_mean$ObsNumForLower = ceiling((Bench_Bratio40_mean$ScenariosIncluded_Bratio40_current*0.5) - (qnorm(0.99)*sqrt(Bench_Bratio40_mean$ScenariosIncluded_Bratio40_current*0.5)))
	Bench_Bratio40_mean$ObsNumForUpper = ceiling((Bench_Bratio40_mean$ScenariosIncluded_Bratio40_current*0.5) + (qnorm(0.99)*sqrt(Bench_Bratio40_mean$ScenariosIncluded_Bratio40_current*0.5)))
	ItemsToLoop = data.frame(species=Bench_Bratio40$species,WPP=Bench_Bratio40$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$Bratio40_UpperCI=NA
	ItemsToLoop$Bratio40_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_Bratio40,Bench_Bratio40$species==ItemsToLoop$species[i] & Bench_Bratio40$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$B_over_Bspr40),]
		Bench_mean_Bratio40_current_Sub = subset(Bench_Bratio40_mean,Bench_Bratio40_mean$species==ItemsToLoop$species[i] & Bench_Bratio40_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$Bratio40_UpperCI[i]=Sub$B_over_Bspr40[Bench_mean_Bratio40_current_Sub$ObsNumForUpper]
		ItemsToLoop$Bratio40_LowerCI[i]=Sub$B_over_Bspr40[Bench_mean_Bratio40_current_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_Bratio40 = merge(Bench_Bratio40_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_Bratio40$ScenariosIncluded_Bratio40_current=NULL 
	Bench_Bratio40$FractionScenariosIncluded_Bratio40_current=NULL 
	Bench_Bratio40$ObsNumForLower=NULL 
	Bench_Bratio40$ObsNumForUpper=NULL
	#Biomass at SPR 40
	Bench_BiomassSpr40 = subset(Benchmarks,Benchmarks$Bspr40_FLAG==0)
	Bench_BiomassSpr40_mean = aggregate.data.frame(Bench_BiomassSpr40$B_spr40_kg,by=list(Bench_BiomassSpr40$species,Bench_BiomassSpr40$WPP),FUN=median,na.rm=TRUE)
	names(Bench_BiomassSpr40_mean) = c("species","WPP","Bspr40_median")
	Samples = aggregate.data.frame(list(Benchmarks$Bspr40_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_Bspr40_current = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_Bspr40_current = Samples$ScenariosIncluded_Bspr40_current/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_BiomassSpr40_mean = merge(Bench_BiomassSpr40_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_BiomassSpr40_mean$ObsNumForLower = ceiling((Bench_BiomassSpr40_mean$ScenariosIncluded_Bspr40_current*0.5) - (qnorm(0.99)*sqrt(Bench_BiomassSpr40_mean$ScenariosIncluded_Bspr40_current*0.5)))
	Bench_BiomassSpr40_mean$ObsNumForUpper = ceiling((Bench_BiomassSpr40_mean$ScenariosIncluded_Bspr40_current*0.5) + (qnorm(0.99)*sqrt(Bench_BiomassSpr40_mean$ScenariosIncluded_Bspr40_current*0.5)))
	ItemsToLoop = data.frame(species=Bench_BiomassSpr40$species,WPP=Bench_BiomassSpr40$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$BiomassSpr40_UpperCI=NA
	ItemsToLoop$BiomassSpr40_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_BiomassSpr40,Bench_BiomassSpr40$species==ItemsToLoop$species[i] & Bench_BiomassSpr40$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$B_spr40_kg),]
		Bench_BiomassSpr40_mean_Sub = subset(Bench_BiomassSpr40_mean,Bench_BiomassSpr40_mean$species==ItemsToLoop$species[i] & Bench_BiomassSpr40_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$BiomassSpr40_UpperCI[i]=Sub$B_spr40_kg[Bench_BiomassSpr40_mean_Sub$ObsNumForUpper]
		ItemsToLoop$BiomassSpr40_LowerCI[i]=Sub$B_spr40_kg[Bench_BiomassSpr40_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_BiomassSpr40 = merge(Bench_BiomassSpr40_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_BiomassSpr40$ScenariosIncluded_Bspr40_current=NULL 
	Bench_BiomassSpr40$FractionScenariosIncluded_Bspr40_current=NULL 
	Bench_BiomassSpr40$ObsNumForLower=NULL 
	Bench_BiomassSpr40$ObsNumForUpper=NULL	
	#Yield at SPR 40
	Bench_YieldSpr40 = subset(Benchmarks,Benchmarks$Bspr40_FLAG==0)
	Bench_YieldSpr40_mean = aggregate.data.frame(Bench_YieldSpr40$Yield_spr40_kg,by=list(Bench_YieldSpr40$species,Bench_YieldSpr40$WPP),FUN=median,na.rm=TRUE)
	names(Bench_YieldSpr40_mean) = c("species","WPP","YieldSpr40_median")
	Samples = aggregate.data.frame(list(Benchmarks$Bspr40_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_YieldSpr40 = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_YieldSpr40 = Samples$ScenariosIncluded_YieldSpr40/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_YieldSpr40_mean = merge(Bench_YieldSpr40_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_YieldSpr40_mean$ObsNumForLower = ceiling((Bench_YieldSpr40_mean$ScenariosIncluded_YieldSpr40*0.5) - (qnorm(0.99)*sqrt(Bench_YieldSpr40_mean$ScenariosIncluded_YieldSpr40*0.5)))
	Bench_YieldSpr40_mean$ObsNumForUpper = ceiling((Bench_YieldSpr40_mean$ScenariosIncluded_YieldSpr40*0.5) + (qnorm(0.99)*sqrt(Bench_YieldSpr40_mean$ScenariosIncluded_YieldSpr40*0.5)))
	ItemsToLoop = data.frame(species=Bench_YieldSpr40$species,WPP=Bench_YieldSpr40$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$YieldSpr40_UpperCI=NA
	ItemsToLoop$YieldSpr40_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_YieldSpr40,Bench_YieldSpr40$species==ItemsToLoop$species[i] & Bench_YieldSpr40$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$Yield_spr40_kg),]
		Bench_YieldSpr40_mean_Sub = subset(Bench_YieldSpr40_mean,Bench_YieldSpr40_mean$species==ItemsToLoop$species[i] & Bench_YieldSpr40_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$YieldSpr40_UpperCI[i]=Sub$Yield_spr40_kg[Bench_YieldSpr40_mean_Sub$ObsNumForUpper]
		ItemsToLoop$YieldSpr40_LowerCI[i]=Sub$Yield_spr40_kg[Bench_YieldSpr40_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_YieldSpr40 = merge(Bench_YieldSpr40_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_YieldSpr40$ScenariosIncluded_YieldSpr40=NULL 
	Bench_YieldSpr40$FractionScenariosIncluded_YieldSpr40=NULL 
	Bench_YieldSpr40$ObsNumForLower=NULL 
	Bench_YieldSpr40$ObsNumForUpper=NULL	
	#Rebuiling Time to spr40
	Bench_RebuildSPR40 = subset(Benchmarks,Benchmarks$F_current_FLAG==0 & Benchmarks$Fratio40_FLAG==0)
	Bench_mean_RebuildSPR40 = aggregate.data.frame(Bench_RebuildSPR40$RebuildTime_spr40,by=list(Bench_RebuildSPR40$species,Bench_RebuildSPR40$WPP),FUN=median,na.rm=TRUE)
	names(Bench_mean_RebuildSPR40) = c("species","WPP","RebuildToSPR40_median")
	Samples = aggregate.data.frame(list(Benchmarks$Fratio40_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_RebuildSPR40 = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_RebuildSPR40 = Samples$ScenariosIncluded_RebuildSPR40/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_mean_RebuildSPR40 = merge(Bench_mean_RebuildSPR40,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_mean_RebuildSPR40$ObsNumForLower = ceiling((Bench_mean_RebuildSPR40$ScenariosIncluded_RebuildSPR40*0.5) - (qnorm(0.99)*sqrt(Bench_mean_RebuildSPR40$ScenariosIncluded_RebuildSPR40*0.5)))
	Bench_mean_RebuildSPR40$ObsNumForUpper = ceiling((Bench_mean_RebuildSPR40$ScenariosIncluded_RebuildSPR40*0.5) + (qnorm(0.99)*sqrt(Bench_mean_RebuildSPR40$ScenariosIncluded_RebuildSPR40*0.5)))
	ItemsToLoop = data.frame(species=Bench_RebuildSPR40$species,WPP=Bench_RebuildSPR40$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$RebuildSPR40_UpperCI=NA
	ItemsToLoop$RebuildSPR40_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_RebuildSPR40,Bench_RebuildSPR40$species==ItemsToLoop$species[i] & Bench_RebuildSPR40$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$RebuildTime_spr40),]
		Bench_mean_RebuildSPR40_Sub = subset(Bench_mean_RebuildSPR40,Bench_mean_RebuildSPR40$species==ItemsToLoop$species[i] & Bench_mean_RebuildSPR40$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$RebuildSPR40_UpperCI[i]=Sub$RebuildTime_spr40[Bench_mean_RebuildSPR40_Sub$ObsNumForUpper]
		ItemsToLoop$RebuildSPR40_LowerCI[i]=Sub$RebuildTime_spr40[Bench_mean_RebuildSPR40_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_RebuildSPR40 = merge(Bench_mean_RebuildSPR40,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_RebuildSPR40$ScenariosIncluded_RebuildSPR40=NULL 
	Bench_RebuildSPR40$FractionScenariosIncluded_RebuildSPR40=NULL 
	Bench_RebuildSPR40$ObsNumForLower=NULL 
	Bench_RebuildSPR40$ObsNumForUpper=NULL
	#BratioMSY
	Bench_BratioMSY = subset(Benchmarks,Benchmarks$BratioMSY_FLAG==0)
	Bench_BratioMSY_mean = aggregate.data.frame(Bench_BratioMSY$B_over_Bmsy,by=list(Bench_BratioMSY$species,Bench_BratioMSY$WPP),FUN=median,na.rm=TRUE)   #geoMean   MedianAbsoluteDeviation
	names(Bench_BratioMSY_mean) = c("species","WPP","B_over_Bmsy_median")
	Samples = aggregate.data.frame(list(Benchmarks$BratioMSY_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_BratioMSY_current = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_BratioMSY_current = Samples$ScenariosIncluded_BratioMSY_current/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_BratioMSY_mean = merge(Bench_BratioMSY_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_BratioMSY_mean$ObsNumForLower = ceiling((Bench_BratioMSY_mean$ScenariosIncluded_BratioMSY_current*0.5) - (qnorm(0.99)*sqrt(Bench_BratioMSY_mean$ScenariosIncluded_BratioMSY_current*0.5)))
	Bench_BratioMSY_mean$ObsNumForUpper = ceiling((Bench_BratioMSY_mean$ScenariosIncluded_BratioMSY_current*0.5) + (qnorm(0.99)*sqrt(Bench_BratioMSY_mean$ScenariosIncluded_BratioMSY_current*0.5)))
	ItemsToLoop = data.frame(species=Bench_BratioMSY$species,WPP=Bench_BratioMSY$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$BratioMSY_UpperCI=NA
	ItemsToLoop$BratioMSY_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_BratioMSY,Bench_BratioMSY$species==ItemsToLoop$species[i] & Bench_BratioMSY$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$B_over_Bmsy),]
		Bench_mean_BratioMSY_current_Sub = subset(Bench_BratioMSY_mean,Bench_BratioMSY_mean$species==ItemsToLoop$species[i] & Bench_BratioMSY_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$BratioMSY_UpperCI[i]=Sub$B_over_Bmsy[Bench_mean_BratioMSY_current_Sub$ObsNumForUpper]
		ItemsToLoop$BratioMSY_LowerCI[i]=Sub$B_over_Bmsy[Bench_mean_BratioMSY_current_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_BratioMSY = merge(Bench_BratioMSY_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_BratioMSY$ScenariosIncluded_BratioMSY_current=NULL 
	Bench_BratioMSY$FractionScenariosIncluded_BratioMSY_current=NULL 
	Bench_BratioMSY$ObsNumForLower=NULL 
	Bench_BratioMSY$ObsNumForUpper=NULL
	#Biomass at MSY
	Bench_BiomassMSY = subset(Benchmarks,Benchmarks$Bmsy_FLAG==0)
	Bench_BiomassMSY_mean = aggregate.data.frame(Bench_BiomassMSY$Bmsy_kg,by=list(Bench_BiomassMSY$species,Bench_BiomassMSY$WPP),FUN=median,na.rm=TRUE)
	names(Bench_BiomassMSY_mean) = c("species","WPP","Bmsy_median")
	Samples = aggregate.data.frame(list(Benchmarks$Bmsy_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_BMSY_current = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_BMSY_current = Samples$ScenariosIncluded_BMSY_current/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_BiomassMSY_mean = merge(Bench_BiomassMSY_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_BiomassMSY_mean$ObsNumForLower = ceiling((Bench_BiomassMSY_mean$ScenariosIncluded_BMSY_current*0.5) - (qnorm(0.99)*sqrt(Bench_BiomassMSY_mean$ScenariosIncluded_BMSY_current*0.5)))
	Bench_BiomassMSY_mean$ObsNumForUpper = ceiling((Bench_BiomassMSY_mean$ScenariosIncluded_BMSY_current*0.5) + (qnorm(0.99)*sqrt(Bench_BiomassMSY_mean$ScenariosIncluded_BMSY_current*0.5)))
	ItemsToLoop = data.frame(species=Bench_BiomassMSY$species,WPP=Bench_BiomassMSY$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$BiomassMSY_UpperCI=NA
	ItemsToLoop$BiomassMSY_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_BiomassMSY,Bench_BiomassMSY$species==ItemsToLoop$species[i] & Bench_BiomassMSY$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$Bmsy_kg),]
		Bench_BiomassMSY_mean_Sub = subset(Bench_BiomassMSY_mean,Bench_BiomassMSY_mean$species==ItemsToLoop$species[i] & Bench_BiomassMSY_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$BiomassMSY_UpperCI[i]=Sub$Bmsy_kg[Bench_BiomassMSY_mean_Sub$ObsNumForUpper]
		ItemsToLoop$BiomassMSY_LowerCI[i]=Sub$Bmsy_kg[Bench_BiomassMSY_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_BiomassMSY = merge(Bench_BiomassMSY_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_BiomassMSY$ScenariosIncluded_BMSY_current=NULL 
	Bench_BiomassMSY$FractionScenariosIncluded_BMSY_current=NULL 
	Bench_BiomassMSY$ObsNumForLower=NULL 
	Bench_BiomassMSY$ObsNumForUpper=NULL	
	#Yield at MSY
	Bench_YieldMSY = subset(Benchmarks,Benchmarks$Bmsy_FLAG==0)
	Bench_YieldMSY_mean = aggregate.data.frame(Bench_YieldMSY$MSY_kg,by=list(Bench_YieldMSY$species,Bench_YieldMSY$WPP),FUN=median,na.rm=TRUE)
	names(Bench_YieldMSY_mean) = c("species","WPP","YieldMSY_median")
	Samples = aggregate.data.frame(list(Benchmarks$Bmsy_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_YieldMSY = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_YieldMSY = Samples$ScenariosIncluded_YieldMSY/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_YieldMSY_mean = merge(Bench_YieldMSY_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_YieldMSY_mean$ObsNumForLower = ceiling((Bench_YieldMSY_mean$ScenariosIncluded_YieldMSY*0.5) - (qnorm(0.99)*sqrt(Bench_YieldMSY_mean$ScenariosIncluded_YieldMSY*0.5)))
	Bench_YieldMSY_mean$ObsNumForUpper = ceiling((Bench_YieldMSY_mean$ScenariosIncluded_YieldMSY*0.5) + (qnorm(0.99)*sqrt(Bench_YieldMSY_mean$ScenariosIncluded_YieldMSY*0.5)))
	ItemsToLoop = data.frame(species=Bench_YieldMSY$species,WPP=Bench_YieldMSY$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$YieldMSY_UpperCI=NA
	ItemsToLoop$YieldMSY_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_YieldMSY,Bench_YieldMSY$species==ItemsToLoop$species[i] & Bench_YieldMSY$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$MSY_kg),]
		Bench_YieldMSY_mean_Sub = subset(Bench_YieldMSY_mean,Bench_YieldMSY_mean$species==ItemsToLoop$species[i] & Bench_YieldMSY_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$YieldMSY_UpperCI[i]=Sub$MSY_kg[Bench_YieldMSY_mean_Sub$ObsNumForUpper]
		ItemsToLoop$YieldMSY_LowerCI[i]=Sub$MSY_kg[Bench_YieldMSY_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_YieldMSY = merge(Bench_YieldMSY_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_YieldMSY$ScenariosIncluded_YieldMSY=NULL 
	Bench_YieldMSY$FractionScenariosIncluded_YieldMSY=NULL 
	Bench_YieldMSY$ObsNumForLower=NULL 
	Bench_YieldMSY$ObsNumForUpper=NULL	
	#Rebuiling Time to MSY
	Bench_RebuildMSY = subset(Benchmarks,Benchmarks$F_current_FLAG==0 & Benchmarks$Fmsy_ratio_FLAG==0)
	Bench_mean_RebuildMSY = aggregate.data.frame(Bench_RebuildMSY$RebuildTime_MSY,by=list(Bench_RebuildMSY$species,Bench_RebuildMSY$WPP),FUN=median,na.rm=TRUE)
	names(Bench_mean_RebuildMSY) = c("species","WPP","RebuildToMSY_median")
	Samples = aggregate.data.frame(list(Benchmarks$Fmsy_ratio_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_RebuildMSY = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_RebuildMSY = Samples$ScenariosIncluded_RebuildMSY/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_mean_RebuildMSY = merge(Bench_mean_RebuildMSY,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_mean_RebuildMSY$ObsNumForLower = ceiling((Bench_mean_RebuildMSY$ScenariosIncluded_RebuildMSY*0.5) - (qnorm(0.99)*sqrt(Bench_mean_RebuildMSY$ScenariosIncluded_RebuildMSY*0.5)))
	Bench_mean_RebuildMSY$ObsNumForUpper = ceiling((Bench_mean_RebuildMSY$ScenariosIncluded_RebuildMSY*0.5) + (qnorm(0.99)*sqrt(Bench_mean_RebuildMSY$ScenariosIncluded_RebuildMSY*0.5)))
	ItemsToLoop = data.frame(species=Bench_RebuildMSY$species,WPP=Bench_RebuildMSY$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$RebuildMSY_UpperCI=NA
	ItemsToLoop$RebuildMSY_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_RebuildMSY,Bench_RebuildMSY$species==ItemsToLoop$species[i] & Bench_RebuildMSY$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$RebuildTime_MSY),]
		Bench_mean_RebuildMSY_Sub = subset(Bench_mean_RebuildMSY,Bench_mean_RebuildMSY$species==ItemsToLoop$species[i] & Bench_mean_RebuildMSY$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$RebuildMSY_UpperCI[i]=Sub$RebuildTime_MSY[Bench_mean_RebuildMSY_Sub$ObsNumForUpper]
		ItemsToLoop$RebuildMSY_LowerCI[i]=Sub$RebuildTime_MSY[Bench_mean_RebuildMSY_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_RebuildMSY = merge(Bench_mean_RebuildMSY,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_RebuildMSY$ScenariosIncluded_RebuildMSY=NULL 
	Bench_RebuildMSY$FractionScenariosIncluded_RebuildMSY=NULL 
	Bench_RebuildMSY$ObsNumForLower=NULL 
	Bench_RebuildMSY$ObsNumForUpper=NULL
	#BratioMEY
	Bench_BratioMEY = subset(Benchmarks,Benchmarks$BratioMEY_FLAG==0)
	Bench_BratioMEY_mean = aggregate.data.frame(Bench_BratioMEY$B_over_Bmey,by=list(Bench_BratioMEY$species,Bench_BratioMEY$WPP),FUN=median,na.rm=TRUE)   #geoMean   MedianAbsoluteDeviation
	names(Bench_BratioMEY_mean) = c("species","WPP","B_over_Bmey_median")
	Samples = aggregate.data.frame(list(Benchmarks$BratioMEY_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_BratioMEY_current = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_BratioMEY_current = Samples$ScenariosIncluded_BratioMEY_current/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_BratioMEY_mean = merge(Bench_BratioMEY_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_BratioMEY_mean$ObsNumForLower = ceiling((Bench_BratioMEY_mean$ScenariosIncluded_BratioMEY_current*0.5) - (qnorm(0.99)*sqrt(Bench_BratioMEY_mean$ScenariosIncluded_BratioMEY_current*0.5)))
	Bench_BratioMEY_mean$ObsNumForUpper = ceiling((Bench_BratioMEY_mean$ScenariosIncluded_BratioMEY_current*0.5) + (qnorm(0.99)*sqrt(Bench_BratioMEY_mean$ScenariosIncluded_BratioMEY_current*0.5)))
	ItemsToLoop = data.frame(species=Bench_BratioMEY$species,WPP=Bench_BratioMEY$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$BratioMEY_UpperCI=NA
	ItemsToLoop$BratioMEY_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_BratioMEY,Bench_BratioMEY$species==ItemsToLoop$species[i] & Bench_BratioMEY$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$B_over_Bmey),]
		Bench_mean_BratioMEY_current_Sub = subset(Bench_BratioMEY_mean,Bench_BratioMEY_mean$species==ItemsToLoop$species[i] & Bench_BratioMEY_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$BratioMEY_UpperCI[i]=Sub$B_over_Bmey[Bench_mean_BratioMEY_current_Sub$ObsNumForUpper]
		ItemsToLoop$BratioMEY_LowerCI[i]=Sub$B_over_Bmey[Bench_mean_BratioMEY_current_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_BratioMEY = merge(Bench_BratioMEY_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_BratioMEY$ScenariosIncluded_BratioMEY_current=NULL 
	Bench_BratioMEY$FractionScenariosIncluded_BratioMEY_current=NULL 
	Bench_BratioMEY$ObsNumForLower=NULL 
	Bench_BratioMEY$ObsNumForUpper=NULL
	#Biomass at MEY
	Bench_BiomassMEY = subset(Benchmarks,Benchmarks$Bmey_FLAG==0)
	Bench_BiomassMEY_mean = aggregate.data.frame(Bench_BiomassMEY$B_MEY_kg,by=list(Bench_BiomassMEY$species,Bench_BiomassMEY$WPP),FUN=median,na.rm=TRUE)
	names(Bench_BiomassMEY_mean) = c("species","WPP","Bmey_median")
	Samples = aggregate.data.frame(list(Benchmarks$Bmey_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_BMEY_current = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_BMEY_current = Samples$ScenariosIncluded_BMEY_current/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_BiomassMEY_mean = merge(Bench_BiomassMEY_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_BiomassMEY_mean$ObsNumForLower = ceiling((Bench_BiomassMEY_mean$ScenariosIncluded_BMEY_current*0.5) - (qnorm(0.99)*sqrt(Bench_BiomassMEY_mean$ScenariosIncluded_BMEY_current*0.5)))
	Bench_BiomassMEY_mean$ObsNumForUpper = ceiling((Bench_BiomassMEY_mean$ScenariosIncluded_BMEY_current*0.5) + (qnorm(0.99)*sqrt(Bench_BiomassMEY_mean$ScenariosIncluded_BMEY_current*0.5)))
	ItemsToLoop = data.frame(species=Bench_BiomassMEY$species,WPP=Bench_BiomassMEY$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$BiomassMEY_UpperCI=NA
	ItemsToLoop$BiomassMEY_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_BiomassMEY,Bench_BiomassMEY$species==ItemsToLoop$species[i] & Bench_BiomassMEY$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$B_MEY_kg),]
		Bench_BiomassMEY_mean_Sub = subset(Bench_BiomassMEY_mean,Bench_BiomassMEY_mean$species==ItemsToLoop$species[i] & Bench_BiomassMEY_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$BiomassMEY_UpperCI[i]=Sub$B_MEY_kg[Bench_BiomassMEY_mean_Sub$ObsNumForUpper]
		ItemsToLoop$BiomassMEY_LowerCI[i]=Sub$B_MEY_kg[Bench_BiomassMEY_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_BiomassMEY = merge(Bench_BiomassMEY_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_BiomassMEY$ScenariosIncluded_BMEY_current=NULL 
	Bench_BiomassMEY$FractionScenariosIncluded_BMEY_current=NULL 
	Bench_BiomassMEY$ObsNumForLower=NULL 
	Bench_BiomassMEY$ObsNumForUpper=NULL	
	#Yield at MEY
	Bench_YieldMEY = subset(Benchmarks,Benchmarks$Bmey_FLAG==0)
	Bench_YieldMEY_mean = aggregate.data.frame(Bench_YieldMEY$Yield_MEY_kg,by=list(Bench_YieldMEY$species,Bench_YieldMEY$WPP),FUN=median,na.rm=TRUE)
	names(Bench_YieldMEY_mean) = c("species","WPP","YieldMEY_median")
	Samples = aggregate.data.frame(list(Benchmarks$Bmey_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_YieldMEY = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_YieldMEY = Samples$ScenariosIncluded_YieldMEY/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_YieldMEY_mean = merge(Bench_YieldMEY_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_YieldMEY_mean$ObsNumForLower = ceiling((Bench_YieldMEY_mean$ScenariosIncluded_YieldMEY*0.5) - (qnorm(0.99)*sqrt(Bench_YieldMEY_mean$ScenariosIncluded_YieldMEY*0.5)))
	Bench_YieldMEY_mean$ObsNumForUpper = ceiling((Bench_YieldMEY_mean$ScenariosIncluded_YieldMEY*0.5) + (qnorm(0.99)*sqrt(Bench_YieldMEY_mean$ScenariosIncluded_YieldMEY*0.5)))
	ItemsToLoop = data.frame(species=Bench_YieldMEY$species,WPP=Bench_YieldMEY$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$YieldMEY_UpperCI=NA
	ItemsToLoop$YieldMEY_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_YieldMEY,Bench_YieldMEY$species==ItemsToLoop$species[i] & Bench_YieldMEY$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$Yield_MEY_kg),]
		Bench_YieldMEY_mean_Sub = subset(Bench_YieldMEY_mean,Bench_YieldMEY_mean$species==ItemsToLoop$species[i] & Bench_YieldMEY_mean$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$YieldMEY_UpperCI[i]=Sub$Yield_MEY_kg[Bench_YieldMEY_mean_Sub$ObsNumForUpper]
		ItemsToLoop$YieldMEY_LowerCI[i]=Sub$Yield_MEY_kg[Bench_YieldMEY_mean_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_YieldMEY = merge(Bench_YieldMEY_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_YieldMEY$ScenariosIncluded_YieldMEY=NULL 
	Bench_YieldMEY$FractionScenariosIncluded_YieldMEY=NULL 
	Bench_YieldMEY$ObsNumForLower=NULL 
	Bench_YieldMEY$ObsNumForUpper=NULL	
	#Rebuiling Time to MEY
	Bench_RebuildMEY = subset(Benchmarks,Benchmarks$F_current_FLAG==0 & Benchmarks$Fmey_ratio_FLAG==0)
	Bench_mean_RebuildMEY = aggregate.data.frame(Bench_RebuildMEY$RebuildTime_MEY,by=list(Bench_RebuildMEY$species,Bench_RebuildMEY$WPP),FUN=median,na.rm=TRUE)
	names(Bench_mean_RebuildMEY) = c("species","WPP","RebuildToMEY_median")
	Samples = aggregate.data.frame(list(Benchmarks$Fmey_ratio_FLAG,Benchmarks$flag),by=list(Benchmarks$species,Benchmarks$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","ScenariosExcluded","TotalScenarios")
	Samples$ScenariosIncluded_RebuildMEY = Samples$TotalScenarios - Samples$ScenariosExcluded
	Samples$FractionScenariosIncluded_RebuildMEY = Samples$ScenariosIncluded_RebuildMEY/Samples$TotalScenarios
	Samples$ScenariosExcluded=NULL
	Samples$TotalScenarios=NULL
	Bench_mean_RebuildMEY = merge(Bench_mean_RebuildMEY,Samples,all.x=TRUE,by=c("species","WPP"))
	Bench_mean_RebuildMEY$ObsNumForLower = ceiling((Bench_mean_RebuildMEY$ScenariosIncluded_RebuildMEY*0.5) - (qnorm(0.99)*sqrt(Bench_mean_RebuildMEY$ScenariosIncluded_RebuildMEY*0.5)))
	Bench_mean_RebuildMEY$ObsNumForUpper = ceiling((Bench_mean_RebuildMEY$ScenariosIncluded_RebuildMEY*0.5) + (qnorm(0.99)*sqrt(Bench_mean_RebuildMEY$ScenariosIncluded_RebuildMEY*0.5)))
	ItemsToLoop = data.frame(species=Bench_RebuildMEY$species,WPP=Bench_RebuildMEY$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$RebuildMEY_UpperCI=NA
	ItemsToLoop$RebuildMEY_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		Sub = subset(Bench_RebuildMEY,Bench_RebuildMEY$species==ItemsToLoop$species[i] & Bench_RebuildMEY$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$RebuildTime_MEY),]
		Bench_mean_RebuildMEY_Sub = subset(Bench_mean_RebuildMEY,Bench_mean_RebuildMEY$species==ItemsToLoop$species[i] & Bench_mean_RebuildMEY$WPP==ItemsToLoop$WPP[i])
		ItemsToLoop$RebuildMEY_UpperCI[i]=Sub$RebuildTime_MEY[Bench_mean_RebuildMEY_Sub$ObsNumForUpper]
		ItemsToLoop$RebuildMEY_LowerCI[i]=Sub$RebuildTime_MEY[Bench_mean_RebuildMEY_Sub$ObsNumForLower]
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	Bench_RebuildMEY = merge(Bench_mean_RebuildMEY,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	Bench_RebuildMEY$ScenariosIncluded_RebuildMEY=NULL 
	Bench_RebuildMEY$FractionScenariosIncluded_RebuildMEY=NULL 
	Bench_RebuildMEY$ObsNumForLower=NULL 
	Bench_RebuildMEY$ObsNumForUpper=NULL
	#Put Benchmarks together 
	BenchmarkSummary = merge(Bench_F_current,Bench_Fspr20,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_Fspr30,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_Fspr40,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_Fmsy,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_Fmey,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_Fratio20,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_Fratio30,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_Fratio40,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_FratioMSY,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_FratioMEY,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_B_current ,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_SPR_current,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_Yield_current,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_Bratio20,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_BiomassSpr20,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_YieldSpr20,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_RebuildSPR20,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_Bratio30,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_BiomassSpr30,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_YieldSpr30,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_RebuildSPR30,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_Bratio40,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_BiomassSpr40,all=TRUE,by=c("species","WPP"))	
	BenchmarkSummary = merge(BenchmarkSummary,Bench_YieldSpr40,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_RebuildSPR40,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_BratioMSY,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_BiomassMSY,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_YieldMSY,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_RebuildMSY,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_BratioMEY,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_BiomassMEY,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_YieldMEY,all=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,Bench_RebuildMEY,all=TRUE,by=c("species","WPP"))	
	#Now, pull cost, revenue, and profit for different projection scenarios: F current, SPR 20, 30 and 40, MSY, and MEY
	SummaryProjCostRevProfit = subset(SummaryProjectionScenarios,SummaryProjectionScenarios$FishingMortaltiyProjectionScenario!="OpenAccess" & SummaryProjectionScenarios$FishingMortaltiyProjectionScenario!="willNotRebuild")
	BenchmarkFlags = subset(Benchmarks,select=c(species,WPP,CatchMSY_Scenario,LifeHistory_Scenario,populationReconstructionMethod,F_current_FLAG,Fspr20_FLAG,Fspr30_FLAG,Fspr40_FLAG,
		Fmsy_FLAG,Fmey_FLAG,Fratio20_FLAG,Fratio30_FLAG,Fratio40_FLAG,Fmsy_ratio_FLAG,Fmey_ratio_FLAG))
	BenchmarkFlags = unique(BenchmarkFlags)
	SummaryProjCostRevProfit = merge(SummaryProjCostRevProfit,BenchmarkFlags,all.x=TRUE,by=c("species","WPP","CatchMSY_Scenario","LifeHistory_Scenario","populationReconstructionMethod"))
	SummaryProjCostRevProfit$flag=1
	SummaryProjCostRevProfit = subset(SummaryProjCostRevProfit,!is.na(SummaryProjCostRevProfit$RevenueAtEnd))
	#CurrentF - Cost, Revenue and Profit 
	ProjCostRevProfit_Fcurrent = subset(SummaryProjCostRevProfit,SummaryProjCostRevProfit$F_current_FLAG==0 & SummaryProjCostRevProfit$FishingMortaltiyProjectionScenario=="CurrentF")
	ProjCostRevProfit_Fcurrent_mean = aggregate.data.frame(list(ProjCostRevProfit_Fcurrent$RevenueAtEnd,ProjCostRevProfit_Fcurrent$CostAtEnd),by=list(ProjCostRevProfit_Fcurrent$species,ProjCostRevProfit_Fcurrent$WPP),FUN=median,na.rm=TRUE)
	names(ProjCostRevProfit_Fcurrent_mean) = c("species","WPP","Revenue_Fcurrent","Cost_Fcurrent")
	Samples = aggregate.data.frame(list(ProjCostRevProfit_Fcurrent$flag),by=list(ProjCostRevProfit_Fcurrent$species,ProjCostRevProfit_Fcurrent$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","TotalScenarios")
	ProjCostRevProfit_Fcurrent_mean = merge(ProjCostRevProfit_Fcurrent_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	ProjCostRevProfit_Fcurrent_mean$ObsNumForLower = ceiling((ProjCostRevProfit_Fcurrent_mean$TotalScenarios*0.5) - (qnorm(0.99)*sqrt(ProjCostRevProfit_Fcurrent_mean$TotalScenarios*0.5)))
	ProjCostRevProfit_Fcurrent_mean$ObsNumForUpper = ceiling((ProjCostRevProfit_Fcurrent_mean$TotalScenarios*0.5) + (qnorm(0.99)*sqrt(ProjCostRevProfit_Fcurrent_mean$TotalScenarios*0.5)))
	ItemsToLoop = data.frame(species=ProjCostRevProfit_Fcurrent$species,WPP=ProjCostRevProfit_Fcurrent$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$Revenue_Fcurrent_UpperCI=NA
	ItemsToLoop$Revenue_Fcurrent_LowerCI=NA
	ItemsToLoop$Cost_Fcurrent_UpperCI=NA
	ItemsToLoop$Cost_Fcurrent_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		ProjCostRevProfit_Fcurrent_mean_Sub = subset(ProjCostRevProfit_Fcurrent_mean,ProjCostRevProfit_Fcurrent_mean$species==ItemsToLoop$species[i] & ProjCostRevProfit_Fcurrent_mean$WPP==ItemsToLoop$WPP[i])
		Sub = subset(ProjCostRevProfit_Fcurrent,ProjCostRevProfit_Fcurrent$species==ItemsToLoop$species[i] & ProjCostRevProfit_Fcurrent$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$RevenueAtEnd),]
		ItemsToLoop$Revenue_Fcurrent_UpperCI[i]=Sub$RevenueAtEnd[ProjCostRevProfit_Fcurrent_mean_Sub$ObsNumForUpper]
		ItemsToLoop$Revenue_Fcurrent_LowerCI[i]=Sub$RevenueAtEnd[ProjCostRevProfit_Fcurrent_mean_Sub$ObsNumForLower]
		Sub = Sub[order(Sub$CostAtEnd),]
		ItemsToLoop$Cost_Fcurrent_UpperCI[i]=Sub$CostAtEnd[ProjCostRevProfit_Fcurrent_mean_Sub$ObsNumForUpper]
		ItemsToLoop$Cost_Fcurrent_LowerCI[i]=Sub$CostAtEnd[ProjCostRevProfit_Fcurrent_mean_Sub$ObsNumForLower]		
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	ProjCostRevProfit_Fcurrent = merge(ProjCostRevProfit_Fcurrent_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	ProjCostRevProfit_Fcurrent$TotalScenarios=NULL 
	ProjCostRevProfit_Fcurrent$ObsNumForLower=NULL 
	ProjCostRevProfit_Fcurrent$ObsNumForUpper=NULL	
	#F_at_SPR20 - Cost, Revenue and Profit 
	ProjCostRevProfit_Fspr20 = subset(SummaryProjCostRevProfit,SummaryProjCostRevProfit$Fspr20_FLAG==0 & SummaryProjCostRevProfit$FishingMortaltiyProjectionScenario=="F_at_SPR20")
	ProjCostRevProfit_Fspr20_mean = aggregate.data.frame(list(ProjCostRevProfit_Fspr20$RevenueAtEnd,ProjCostRevProfit_Fspr20$CostAtEnd),by=list(ProjCostRevProfit_Fspr20$species,ProjCostRevProfit_Fspr20$WPP),FUN=median,na.rm=TRUE)
	names(ProjCostRevProfit_Fspr20_mean) = c("species","WPP","Revenue_Fspr20","Cost_Fspr20")
	Samples = aggregate.data.frame(list(ProjCostRevProfit_Fspr20$flag),by=list(ProjCostRevProfit_Fspr20$species,ProjCostRevProfit_Fspr20$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","TotalScenarios")
	ProjCostRevProfit_Fspr20_mean = merge(ProjCostRevProfit_Fspr20_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	ProjCostRevProfit_Fspr20_mean$ObsNumForLower = ceiling((ProjCostRevProfit_Fspr20_mean$TotalScenarios*0.5) - (qnorm(0.99)*sqrt(ProjCostRevProfit_Fspr20_mean$TotalScenarios*0.5)))
	ProjCostRevProfit_Fspr20_mean$ObsNumForUpper = ceiling((ProjCostRevProfit_Fspr20_mean$TotalScenarios*0.5) + (qnorm(0.99)*sqrt(ProjCostRevProfit_Fspr20_mean$TotalScenarios*0.5)))
	ItemsToLoop = data.frame(species=ProjCostRevProfit_Fspr20$species,WPP=ProjCostRevProfit_Fspr20$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$Revenue_Fspr20_UpperCI=NA
	ItemsToLoop$Revenue_Fspr20_LowerCI=NA
	ItemsToLoop$Cost_Fspr20_UpperCI=NA
	ItemsToLoop$Cost_Fspr20_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		ProjCostRevProfit_Fspr20_mean_Sub = subset(ProjCostRevProfit_Fspr20_mean,ProjCostRevProfit_Fspr20_mean$species==ItemsToLoop$species[i] & ProjCostRevProfit_Fspr20_mean$WPP==ItemsToLoop$WPP[i])
		Sub = subset(ProjCostRevProfit_Fspr20,ProjCostRevProfit_Fspr20$species==ItemsToLoop$species[i] & ProjCostRevProfit_Fspr20$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$RevenueAtEnd),]
		ItemsToLoop$Revenue_Fspr20_UpperCI[i]=Sub$RevenueAtEnd[ProjCostRevProfit_Fspr20_mean_Sub$ObsNumForUpper]
		ItemsToLoop$Revenue_Fspr20_LowerCI[i]=Sub$RevenueAtEnd[ProjCostRevProfit_Fspr20_mean_Sub$ObsNumForLower]
		Sub = Sub[order(Sub$CostAtEnd),]
		ItemsToLoop$Cost_Fspr20_UpperCI[i]=Sub$CostAtEnd[ProjCostRevProfit_Fspr20_mean_Sub$ObsNumForUpper]
		ItemsToLoop$Cost_Fspr20_LowerCI[i]=Sub$CostAtEnd[ProjCostRevProfit_Fspr20_mean_Sub$ObsNumForLower]		
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	ProjCostRevProfit_Fspr20 = merge(ProjCostRevProfit_Fspr20_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	ProjCostRevProfit_Fspr20$TotalScenarios=NULL 
	ProjCostRevProfit_Fspr20$ObsNumForLower=NULL 
	ProjCostRevProfit_Fspr20$ObsNumForUpper=NULL
	#F_at_SPR30 - Cost, Revenue and Profit 
	ProjCostRevProfit_Fspr30 = subset(SummaryProjCostRevProfit,SummaryProjCostRevProfit$Fspr30_FLAG==0 & SummaryProjCostRevProfit$FishingMortaltiyProjectionScenario=="F_at_SPR30")
	ProjCostRevProfit_Fspr30_mean = aggregate.data.frame(list(ProjCostRevProfit_Fspr30$RevenueAtEnd,ProjCostRevProfit_Fspr30$CostAtEnd),by=list(ProjCostRevProfit_Fspr30$species,ProjCostRevProfit_Fspr30$WPP),FUN=median,na.rm=TRUE)
	names(ProjCostRevProfit_Fspr30_mean) = c("species","WPP","Revenue_Fspr30","Cost_Fspr30")
	Samples = aggregate.data.frame(list(ProjCostRevProfit_Fspr30$flag),by=list(ProjCostRevProfit_Fspr30$species,ProjCostRevProfit_Fspr30$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","TotalScenarios")
	ProjCostRevProfit_Fspr30_mean = merge(ProjCostRevProfit_Fspr30_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	ProjCostRevProfit_Fspr30_mean$ObsNumForLower = ceiling((ProjCostRevProfit_Fspr30_mean$TotalScenarios*0.5) - (qnorm(0.99)*sqrt(ProjCostRevProfit_Fspr30_mean$TotalScenarios*0.5)))
	ProjCostRevProfit_Fspr30_mean$ObsNumForUpper = ceiling((ProjCostRevProfit_Fspr30_mean$TotalScenarios*0.5) + (qnorm(0.99)*sqrt(ProjCostRevProfit_Fspr30_mean$TotalScenarios*0.5)))
	ItemsToLoop = data.frame(species=ProjCostRevProfit_Fspr30$species,WPP=ProjCostRevProfit_Fspr30$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$Revenue_Fspr30_UpperCI=NA
	ItemsToLoop$Revenue_Fspr30_LowerCI=NA
	ItemsToLoop$Cost_Fspr30_UpperCI=NA
	ItemsToLoop$Cost_Fspr30_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		ProjCostRevProfit_Fspr30_mean_Sub = subset(ProjCostRevProfit_Fspr30_mean,ProjCostRevProfit_Fspr30_mean$species==ItemsToLoop$species[i] & ProjCostRevProfit_Fspr30_mean$WPP==ItemsToLoop$WPP[i])
		Sub = subset(ProjCostRevProfit_Fspr30,ProjCostRevProfit_Fspr30$species==ItemsToLoop$species[i] & ProjCostRevProfit_Fspr30$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$RevenueAtEnd),]
		ItemsToLoop$Revenue_Fspr30_UpperCI[i]=Sub$RevenueAtEnd[ProjCostRevProfit_Fspr30_mean_Sub$ObsNumForUpper]
		ItemsToLoop$Revenue_Fspr30_LowerCI[i]=Sub$RevenueAtEnd[ProjCostRevProfit_Fspr30_mean_Sub$ObsNumForLower]
		Sub = Sub[order(Sub$CostAtEnd),]
		ItemsToLoop$Cost_Fspr30_UpperCI[i]=Sub$CostAtEnd[ProjCostRevProfit_Fspr30_mean_Sub$ObsNumForUpper]
		ItemsToLoop$Cost_Fspr30_LowerCI[i]=Sub$CostAtEnd[ProjCostRevProfit_Fspr30_mean_Sub$ObsNumForLower]		
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	ProjCostRevProfit_Fspr30 = merge(ProjCostRevProfit_Fspr30_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	ProjCostRevProfit_Fspr30$TotalScenarios=NULL 
	ProjCostRevProfit_Fspr30$ObsNumForLower=NULL 
	ProjCostRevProfit_Fspr30$ObsNumForUpper=NULL
	#F_at_SPR40 - Cost, Revenue and Profit 
	ProjCostRevProfit_Fspr40 = subset(SummaryProjCostRevProfit,SummaryProjCostRevProfit$Fspr40_FLAG==0 & SummaryProjCostRevProfit$FishingMortaltiyProjectionScenario=="F_at_SPR40")
	ProjCostRevProfit_Fspr40_mean = aggregate.data.frame(list(ProjCostRevProfit_Fspr40$RevenueAtEnd,ProjCostRevProfit_Fspr40$CostAtEnd),by=list(ProjCostRevProfit_Fspr40$species,ProjCostRevProfit_Fspr40$WPP),FUN=median,na.rm=TRUE)
	names(ProjCostRevProfit_Fspr40_mean) = c("species","WPP","Revenue_Fspr40","Cost_Fspr40")
	Samples = aggregate.data.frame(list(ProjCostRevProfit_Fspr40$flag),by=list(ProjCostRevProfit_Fspr40$species,ProjCostRevProfit_Fspr40$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","TotalScenarios")
	ProjCostRevProfit_Fspr40_mean = merge(ProjCostRevProfit_Fspr40_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	ProjCostRevProfit_Fspr40_mean$ObsNumForLower = ceiling((ProjCostRevProfit_Fspr40_mean$TotalScenarios*0.5) - (qnorm(0.99)*sqrt(ProjCostRevProfit_Fspr40_mean$TotalScenarios*0.5)))
	ProjCostRevProfit_Fspr40_mean$ObsNumForUpper = ceiling((ProjCostRevProfit_Fspr40_mean$TotalScenarios*0.5) + (qnorm(0.99)*sqrt(ProjCostRevProfit_Fspr40_mean$TotalScenarios*0.5)))
	ItemsToLoop = data.frame(species=ProjCostRevProfit_Fspr40$species,WPP=ProjCostRevProfit_Fspr40$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$Revenue_Fspr40_UpperCI=NA
	ItemsToLoop$Revenue_Fspr40_LowerCI=NA
	ItemsToLoop$Cost_Fspr40_UpperCI=NA
	ItemsToLoop$Cost_Fspr40_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		ProjCostRevProfit_Fspr40_mean_Sub = subset(ProjCostRevProfit_Fspr40_mean,ProjCostRevProfit_Fspr40_mean$species==ItemsToLoop$species[i] & ProjCostRevProfit_Fspr40_mean$WPP==ItemsToLoop$WPP[i])
		Sub = subset(ProjCostRevProfit_Fspr40,ProjCostRevProfit_Fspr40$species==ItemsToLoop$species[i] & ProjCostRevProfit_Fspr40$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$RevenueAtEnd),]
		ItemsToLoop$Revenue_Fspr40_UpperCI[i]=Sub$RevenueAtEnd[ProjCostRevProfit_Fspr40_mean_Sub$ObsNumForUpper]
		ItemsToLoop$Revenue_Fspr40_LowerCI[i]=Sub$RevenueAtEnd[ProjCostRevProfit_Fspr40_mean_Sub$ObsNumForLower]
		Sub = Sub[order(Sub$CostAtEnd),]
		ItemsToLoop$Cost_Fspr40_UpperCI[i]=Sub$CostAtEnd[ProjCostRevProfit_Fspr40_mean_Sub$ObsNumForUpper]
		ItemsToLoop$Cost_Fspr40_LowerCI[i]=Sub$CostAtEnd[ProjCostRevProfit_Fspr40_mean_Sub$ObsNumForLower]		
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	ProjCostRevProfit_Fspr40 = merge(ProjCostRevProfit_Fspr40_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	ProjCostRevProfit_Fspr40$ScenariosIncluded_Fspr40_RevProfit=NULL 
	ProjCostRevProfit_Fspr40$FractionScenariosIncluded_Fspr40_RevProfit=NULL 
	ProjCostRevProfit_Fspr40$ObsNumForLower=NULL 
	ProjCostRevProfit_Fspr40$ObsNumForUpper=NULL
	#MSY - Cost, Revenue and Profit 
	ProjCostRevProfit_msy = subset(SummaryProjCostRevProfit,SummaryProjCostRevProfit$Fmsy_FLAG==0 & SummaryProjCostRevProfit$FishingMortaltiyProjectionScenario=="MSY")
	ProjCostRevProfit_msy_mean = aggregate.data.frame(list(ProjCostRevProfit_msy$RevenueAtEnd,ProjCostRevProfit_msy$CostAtEnd),by=list(ProjCostRevProfit_msy$species,ProjCostRevProfit_msy$WPP),FUN=median,na.rm=TRUE)
	names(ProjCostRevProfit_msy_mean) = c("species","WPP","Revenue_MSY","Cost_MSY")
	Samples = aggregate.data.frame(list(ProjCostRevProfit_msy$flag),by=list(ProjCostRevProfit_msy$species,ProjCostRevProfit_msy$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","TotalScenarios")
	ProjCostRevProfit_msy_mean = merge(ProjCostRevProfit_msy_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	ProjCostRevProfit_msy_mean$ObsNumForLower = ceiling((ProjCostRevProfit_msy_mean$TotalScenarios*0.5) - (qnorm(0.99)*sqrt(ProjCostRevProfit_msy_mean$TotalScenarios*0.5)))
	ProjCostRevProfit_msy_mean$ObsNumForUpper = ceiling((ProjCostRevProfit_msy_mean$TotalScenarios*0.5) + (qnorm(0.99)*sqrt(ProjCostRevProfit_msy_mean$TotalScenarios*0.5)))
	ItemsToLoop = data.frame(species=ProjCostRevProfit_msy$species,WPP=ProjCostRevProfit_msy$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$Revenue_msy_UpperCI=NA
	ItemsToLoop$Revenue_msy_LowerCI=NA
	ItemsToLoop$Cost_msy_UpperCI=NA
	ItemsToLoop$Cost_msy_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		ProjCostRevProfit_msy_mean_Sub = subset(ProjCostRevProfit_msy_mean,ProjCostRevProfit_msy_mean$species==ItemsToLoop$species[i] & ProjCostRevProfit_msy_mean$WPP==ItemsToLoop$WPP[i])
		Sub = subset(ProjCostRevProfit_msy,ProjCostRevProfit_msy$species==ItemsToLoop$species[i] & ProjCostRevProfit_msy$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$RevenueAtEnd),]
		ItemsToLoop$Revenue_msy_UpperCI[i]=Sub$RevenueAtEnd[ProjCostRevProfit_msy_mean_Sub$ObsNumForUpper]
		ItemsToLoop$Revenue_msy_LowerCI[i]=Sub$RevenueAtEnd[ProjCostRevProfit_msy_mean_Sub$ObsNumForLower]
		Sub = Sub[order(Sub$CostAtEnd),]
		ItemsToLoop$Cost_msy_UpperCI[i]=Sub$CostAtEnd[ProjCostRevProfit_msy_mean_Sub$ObsNumForUpper]
		ItemsToLoop$Cost_msy_LowerCI[i]=Sub$CostAtEnd[ProjCostRevProfit_msy_mean_Sub$ObsNumForLower]		
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	ProjCostRevProfit_msy = merge(ProjCostRevProfit_msy_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	ProjCostRevProfit_msy$TotalScenarios=NULL 
	ProjCostRevProfit_msy$ObsNumForLower=NULL 
	ProjCostRevProfit_msy$ObsNumForUpper=NULL
	#MEY - Cost, Revenue and Profit 	
	ProjCostRevProfit_mey = subset(SummaryProjCostRevProfit,SummaryProjCostRevProfit$Fmey_FLAG==0 & SummaryProjCostRevProfit$FishingMortaltiyProjectionScenario=="MEY")
	ProjCostRevProfit_mey_mean = aggregate.data.frame(list(ProjCostRevProfit_mey$RevenueAtEnd,ProjCostRevProfit_mey$CostAtEnd),by=list(ProjCostRevProfit_mey$species,ProjCostRevProfit_mey$WPP),FUN=median,na.rm=TRUE)
	names(ProjCostRevProfit_mey_mean) = c("species","WPP","Revenue_MEY","Cost_MEY")
	Samples = aggregate.data.frame(list(ProjCostRevProfit_mey$flag),by=list(ProjCostRevProfit_mey$species,ProjCostRevProfit_mey$WPP),FUN=sum,na.rm=TRUE)
	names(Samples) = c("species","WPP","TotalScenarios")
	ProjCostRevProfit_mey_mean = merge(ProjCostRevProfit_mey_mean,Samples,all.x=TRUE,by=c("species","WPP"))
	ProjCostRevProfit_mey_mean$ObsNumForLower = ceiling((ProjCostRevProfit_mey_mean$TotalScenarios*0.5) - (qnorm(0.99)*sqrt(ProjCostRevProfit_mey_mean$TotalScenarios*0.5)))
	ProjCostRevProfit_mey_mean$ObsNumForUpper = ceiling((ProjCostRevProfit_mey_mean$TotalScenarios*0.5) + (qnorm(0.99)*sqrt(ProjCostRevProfit_mey_mean$TotalScenarios*0.5)))
	ItemsToLoop = data.frame(species=ProjCostRevProfit_mey$species,WPP=ProjCostRevProfit_mey$WPP)
	ItemsToLoop = unique(ItemsToLoop)
	ItemsToLoop$Revenue_mey_UpperCI=NA
	ItemsToLoop$Revenue_mey_LowerCI=NA
	ItemsToLoop$Cost_mey_UpperCI=NA
	ItemsToLoop$Cost_mey_LowerCI=NA
	for(i in 1:dim(ItemsToLoop)[1])
	{
		ProjCostRevProfit_mey_mean_Sub = subset(ProjCostRevProfit_mey_mean,ProjCostRevProfit_mey_mean$species==ItemsToLoop$species[i] & ProjCostRevProfit_mey_mean$WPP==ItemsToLoop$WPP[i])
		Sub = subset(ProjCostRevProfit_mey,ProjCostRevProfit_mey$species==ItemsToLoop$species[i] & ProjCostRevProfit_mey$WPP==ItemsToLoop$WPP[i])
		Sub = Sub[order(Sub$RevenueAtEnd),]
		ItemsToLoop$Revenue_mey_UpperCI[i]=Sub$RevenueAtEnd[ProjCostRevProfit_mey_mean_Sub$ObsNumForUpper]
		ItemsToLoop$Revenue_mey_LowerCI[i]=Sub$RevenueAtEnd[ProjCostRevProfit_mey_mean_Sub$ObsNumForLower]
		Sub = Sub[order(Sub$CostAtEnd),]
		ItemsToLoop$Cost_mey_UpperCI[i]=Sub$CostAtEnd[ProjCostRevProfit_mey_mean_Sub$ObsNumForUpper]
		ItemsToLoop$Cost_mey_LowerCI[i]=Sub$CostAtEnd[ProjCostRevProfit_mey_mean_Sub$ObsNumForLower]		
		if(i%%50==0)
		{
			print(paste("Working on ",i," of ",dim(ItemsToLoop)[1],"...",sep=""))
			flush.console()
		}
	}
	ProjCostRevProfit_mey = merge(ProjCostRevProfit_mey_mean,ItemsToLoop,all.x=TRUE,by=c("species","WPP"))
	ProjCostRevProfit_mey$TotalScenarios=NULL 
	ProjCostRevProfit_mey$ObsNumForLower=NULL 
	ProjCostRevProfit_mey$ObsNumForUpper=NULL
	#Add the economic information to the output
	BenchmarkSummary = merge(BenchmarkSummary,ProjCostRevProfit_Fcurrent,all.x=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,ProjCostRevProfit_Fspr20,all.x=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,ProjCostRevProfit_Fspr30,all.x=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,ProjCostRevProfit_Fspr40,all.x=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,ProjCostRevProfit_msy,all.x=TRUE,by=c("species","WPP"))
	BenchmarkSummary = merge(BenchmarkSummary,ProjCostRevProfit_mey,all.x=TRUE,by=c("species","WPP"))
	write.table(BenchmarkSummary,paste(PATH_output,"BenchmarkSummary.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
}





#######################################################################################################################################################
#######################################################################################################################################################
################# MAKE PLOTS AND DEVELOP OUTPUT PRODUCTS ##############################################################################################
#######################################################################################################################################################
#######################################################################################################################################################


#################################### Create Kobe Plots By Species #####################################################
if(createKobePlots==TRUE)
{
  ExtrapolatedLandings = read.table(paste(PATH_output,"ExtrapolatedLandings.csv",sep="/"),header=TRUE,sep=",")
  SpeciesFamily = subset(ExtrapolatedLandings,select=c(species,FamilyCommonName))
  SpeciesFamily = unique(SpeciesFamily)
  CurrentStatusValues = read.table(paste(PATH_output,"Benchmarks_withFlags.csv",sep="/"),header=TRUE,sep=",",as.is=TRUE)
  CurrentStatusValues = subset(CurrentStatusValues,select=c(species,WPP,CatchMSY_Scenario,LifeHistory_Scenario,populationReconstructionMethod,F_over_Fspr40,B_over_Bspr20,Fratio40_FLAG,Bratio20_FLAG))
  CurrentStatusValues = subset(CurrentStatusValues,CurrentStatusValues$Fratio40_FLAG==0 & CurrentStatusValues$Bratio20_FLAG==0)
  CurrentStatusValues = merge(CurrentStatusValues,SpeciesFamily,all.x=TRUE,by=c("species"))
  CurrentStatusValues$F_ratioAboveOne = 0
  CurrentStatusValues$F_ratioAboveOne[CurrentStatusValues$F_over_Fspr40>=1]=1
  CurrentStatusValues$B_ratioAboveOne = 0
  CurrentStatusValues$B_ratioAboveOne[CurrentStatusValues$B_over_Bspr20>=1]=1
  CurrentStatusValues$Quadrant=0      #Quadrant 1=Overfished and Overfishing; 2=Overfishing, not overfished; 3=not overfishing, not overfished; 4=not overfishing, overfished
  CurrentStatusValues$Quadrant[CurrentStatusValues$F_over_Fspr40>=1 & CurrentStatusValues$B_over_Bspr20<1]=1
  CurrentStatusValues$Quadrant[CurrentStatusValues$F_over_Fspr40>=1 & CurrentStatusValues$B_over_Bspr20>=1]=2
  CurrentStatusValues$Quadrant[CurrentStatusValues$F_over_Fspr40<1 & CurrentStatusValues$B_over_Bspr20>=1]=3
  CurrentStatusValues$Quadrant[CurrentStatusValues$F_over_Fspr40<1 & CurrentStatusValues$B_over_Bspr20<1]=4
  StatusTable = data.frame(table(CurrentStatusValues$species,CurrentStatusValues$WPP,CurrentStatusValues$Quadrant))
  names(StatusTable) = c("species","WPP","Quadrant","Freq")
  StatusTable = reshape(StatusTable,v.names="Freq",idvar=c("species","WPP"),timevar="Quadrant",direction="wide")
  StatusTable$Total = StatusTable$Freq.1 + StatusTable$Freq.2 + StatusTable$Freq.3 + StatusTable$Freq.4
  StatusTable = subset(StatusTable,StatusTable$Total>0)
  StatusTable$Percent.1 = round((StatusTable$Freq.1/StatusTable$Total)*100)
  StatusTable$Percent.2 = round((StatusTable$Freq.2/StatusTable$Total)*100)
  StatusTable$Percent.3 = round((StatusTable$Freq.3/StatusTable$Total)*100)
  StatusTable$Percent.4 = round((StatusTable$Freq.4/StatusTable$Total)*100)
  StatusTable[is.na(StatusTable)]=0
  StatusTable$species = as.character(StatusTable$species)
  StatusTable$WPP = as.character(StatusTable$WPP)
  write.table(StatusTable,paste(PATH_output,"StatusTable.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
  #Now make the Kobe plots
  BenchmarksToPlot = subset(CurrentStatusValues,select=c(species,WPP))
  BenchmarksToPlot = unique(BenchmarksToPlot)
  CurrentStatusValues$pch_arg=0
  CurrentStatusValues$pch_arg[CurrentStatusValues$populationReconstructionMethod=="LBSPR"]=15
  CurrentStatusValues$pch_arg[CurrentStatusValues$populationReconstructionMethod=="BevertonHoltInstantaneousMortality"]=16
  CurrentStatusValues$pch_arg[CurrentStatusValues$populationReconstructionMethod=="LIME"]=17
  CurrentStatusValues$pch_arg[CurrentStatusValues$populationReconstructionMethod=="TropFishR"]=18
  CurrentStatusValues$pch_arg[CurrentStatusValues$populationReconstructionMethod=="SolveBaranov"]=19
  ColorsFor_benchmarkPlot_CmsyScenarios = data.frame(CatchMSY_Scenario=sort(unique(CurrentStatusValues$CatchMSY_Scenario)),Colors=c("red","royalblue","green","purple","orange","brown"),stringsAsFactors=FALSE)
  LengthReconstructionSymbols = unique(subset(CurrentStatusValues,select=c(populationReconstructionMethod,pch_arg)))
  CurrentStatusValues = merge(CurrentStatusValues,ColorsFor_benchmarkPlot_CmsyScenarios,all.x=TRUE,by=c("CatchMSY_Scenario"))
  for(i in 1:dim(BenchmarksToPlot)[1])
  {
    thisIteration = CurrentStatusValues[CurrentStatusValues$species==BenchmarksToPlot$species[i] & CurrentStatusValues$WPP==BenchmarksToPlot$WPP[i],]
    FolderNameSpp = gsub(" ","_",BenchmarksToPlot$species[i])
    png(paste(paste(PATH_benchmarkPlots,FolderNameSpp,sep="/"),paste(paste(gsub(" ","_",BenchmarksToPlot$species[i]),BenchmarksToPlot$WPP[i],sep="_"),"_Benchmarks.png",sep=""),sep="/"),units="px",width=5800,height=3200,res=600)
    layout(matrix(c(1,2), nrow = 1), widths = c(1,0.5))
    par(mar = c(5, 4, 4, 2) + 0.1)
    maxBratio = max(thisIteration$B_over_Bspr20)
    maxFratio = max(thisIteration$F_over_Fspr40)
    if(maxBratio<2) {maxBratio=2}
    if(maxBratio>10) {maxBratio=10}
    if(maxFratio<2) {maxFratio=2}
    if(maxFratio>10) {maxFratio=10}
    plot(NULL,xlab="Bcurrent/B@SPR20",ylab="Fcurrent/F@SPR40",xlim=c(0,ceiling(maxBratio)),ylim=c(0,ceiling(maxFratio)),xaxs="i",yaxs="i")
    rect(xleft=0,ybottom=0,xright=1,ytop=1,par("usr")[3] - 1, y0 + 3, par("usr")[4] + 1,border="khaki",col="khaki")
    rect(xleft=1,ybottom=0,xright=ceiling(maxBratio),ytop=1,par("usr")[3] - 1, y0 + 3, par("usr")[4] + 1,border="lightgreen",col="lightgreen")
    rect(xleft=0,ybottom=1,xright=1,ytop=ceiling(maxFratio),par("usr")[3] - 1, y0 + 3, par("usr")[4] + 1,border="plum1",col="plum1")
    rect(xleft=1,ybottom=1,xright=ceiling(maxBratio),ytop=ceiling(maxFratio),par("usr")[3] - 1, y0 + 3, par("usr")[4] + 1,border="bisque",col="bisque") #tan1
    for(j in 1:dim(thisIteration)[1])
    {
      if(j==1)
      {
        points(thisIteration$B_over_Bspr20[j],thisIteration$F_over_Fspr40[j],xlab="Bcurrent/B@SPR20",ylab="Fcurrent/F@SPR40",pch=thisIteration$pch_arg[j],col=thisIteration$Colors[j])
      }
      if(j>1)
      {
        points(thisIteration$B_over_Bspr20[j],thisIteration$F_over_Fspr40[j],xlab="Bcurrent/B@SPR20",ylab="Fcurrent/F@SPR40",pch=thisIteration$pch_arg[j],col=thisIteration$Colors[j])
      }
    }
    if(ceiling(maxFratio)>3)
    {
      yMultiplier = ceiling(maxFratio)/18
    } else {
      yMultiplier = 0.1
    }
    if(ceiling(maxBratio)>3)
    {
      xMultiplier = ceiling(maxBratio)/18
    } else {
      xMultiplier = 0.1
    }
    StatusPercents = StatusTable[StatusTable$species==BenchmarksToPlot$species[i] & StatusTable$WPP==BenchmarksToPlot$WPP[i],]
    rect(xleft=0.5-xMultiplier,ybottom=(((ceiling(maxFratio)-1)/2)-yMultiplier)+1,xright=0.5+xMultiplier,ytop=(((ceiling(maxFratio)-1)/2)+yMultiplier)+1,par("usr")[3] - 1, y0 + 3, par("usr")[4] + 1,border="black",col="white",lty=1,lwd=1)
    text(x=0.5,y=((ceiling(maxFratio)-1)/2)+1,labels=paste(StatusPercents$Percent.1,"%",sep=""))
    rect(xleft=(((ceiling(maxBratio)-1)/2)-xMultiplier)+1,ybottom=(((ceiling(maxFratio)-1)/2)-yMultiplier)+1,xright=(((ceiling(maxBratio)-1)/2)+xMultiplier)+1,ytop=(((ceiling(maxFratio)-1)/2)+yMultiplier)+1,par("usr")[3] - 1, y0 + 3, par("usr")[4] + 1,border="black",col="white",lty=1,lwd=1)
    text(x=((ceiling(maxBratio)-1)/2)+1,y=((ceiling(maxFratio)-1)/2)+1,labels=paste(StatusPercents$Percent.2,"%",sep=""))
    rect(xleft=(((ceiling(maxBratio)-1)/2)-xMultiplier)+1,ybottom=0.5-yMultiplier,xright=(((ceiling(maxBratio)-1)/2)+xMultiplier)+1,ytop=0.5+yMultiplier,par("usr")[3] - 1, y0 + 3, par("usr")[4] + 1,border="black",col="white",lty=1,lwd=1)
    text(x=((ceiling(maxBratio)-1)/2)+1,y=0.5,labels=paste(StatusPercents$Percent.3,"%",sep=""))
    rect(xleft=0.5-xMultiplier,ybottom=0.5-yMultiplier,xright=0.5+xMultiplier,ytop=0.5+yMultiplier,par("usr")[3] - 1, y0 + 3, par("usr")[4] + 1,border="black",col="white",lty=1,lwd=1)
    text(x=0.5,y=0.5,labels=paste(StatusPercents$Percent.4,"%",sep=""))
    abline(v=1,col="black",lty=3,lwd=1)
    abline(h=1,col="black",lty=3,lwd=1)
    title(paste("Current Status: ",BenchmarksToPlot$species[i],sep=" "))
    mtext(paste("WPP = ",BenchmarksToPlot$WPP[i],sep=""))
    par(mar = c(5, 0, 4, 2))
    plot(thisIteration$B_over_Bspr20[dim(thisIteration)[1]],thisIteration$F_over_Fspr40[dim(thisIteration)[1]],type="n",axes=FALSE,ann=FALSE,xlim=c(0,maxBratio),ylim=c(0,maxFratio))
    LengthReconstructionSymbols_forLegend = LengthReconstructionSymbols
    LengthReconstructionSymbols_forLegend$populationReconstructionMethod = as.character(LengthReconstructionSymbols_forLegend$populationReconstructionMethod)
    LengthReconstructionSymbols_forLegend$populationReconstructionMethod[LengthReconstructionSymbols_forLegend$populationReconstructionMethod=="BevertonHoltInstantaneousMortality"]="BevHoltF"
    legend("topleft",legend=c(as.character(LengthReconstructionSymbols_forLegend$populationReconstructionMethod),as.character(ColorsFor_benchmarkPlot_CmsyScenarios$CatchMSY_Scenario)),lty=c(rep(NA,times=dim(LengthReconstructionSymbols_forLegend)[1]),rep(1,times=dim(ColorsFor_benchmarkPlot_CmsyScenarios)[1])),lwd=c(rep(NA,times=dim(LengthReconstructionSymbols_forLegend)[1]),rep(2,times=dim(ColorsFor_benchmarkPlot_CmsyScenarios)[1])),col=c(rep("black",times=dim(LengthReconstructionSymbols_forLegend)[1]),ColorsFor_benchmarkPlot_CmsyScenarios$Colors),pch=c(LengthReconstructionSymbols_forLegend$pch_arg,rep(NA,times=dim(ColorsFor_benchmarkPlot_CmsyScenarios)[1])),cex=0.8)
    dev.off()
  }

}



#################################### Create Kobe Plots By Fish Family #####################################################
if(createKobePlots==TRUE)
{
  ExtrapolatedLandings = read.table(paste(PATH_output,"ExtrapolatedLandings.csv",sep="/"),header=TRUE,sep=",")
  SpeciesFamily = subset(ExtrapolatedLandings,select=c(species,FamilyCommonName))
  SpeciesFamily = unique(SpeciesFamily)
  CurrentStatusValues = read.table(paste(PATH_output,"Benchmarks_withFlags.csv",sep="/"),header=TRUE,sep=",",as.is=TRUE)
  CurrentStatusValues = subset(CurrentStatusValues,select=c(species,WPP,CatchMSY_Scenario,LifeHistory_Scenario,populationReconstructionMethod,F_over_Fspr40,B_over_Bspr20,Fratio40_FLAG,Bratio20_FLAG))
  CurrentStatusValues = subset(CurrentStatusValues,CurrentStatusValues$Fratio40_FLAG==0 & CurrentStatusValues$Bratio20_FLAG==0)
  CurrentStatusValues = merge(CurrentStatusValues,SpeciesFamily,all.x=TRUE,by=c("species"))
  CurrentStatusValues$F_ratioAboveOne = 0
  CurrentStatusValues$F_ratioAboveOne[CurrentStatusValues$F_over_Fspr40>=1]=1
  CurrentStatusValues$B_ratioAboveOne = 0
  CurrentStatusValues$B_ratioAboveOne[CurrentStatusValues$B_over_Bspr20>=1]=1
  CurrentStatusValues$Quadrant=0      #Quadrant 1=Overfished and Overfishing; 2=Overfishing, not overfished; 3=not overfishing, not overfished; 4=not overfishing, overfished
  CurrentStatusValues$Quadrant[CurrentStatusValues$F_over_Fspr40>=1 & CurrentStatusValues$B_over_Bspr20<1]=1
  CurrentStatusValues$Quadrant[CurrentStatusValues$F_over_Fspr40>=1 & CurrentStatusValues$B_over_Bspr20>=1]=2
  CurrentStatusValues$Quadrant[CurrentStatusValues$F_over_Fspr40<1 & CurrentStatusValues$B_over_Bspr20>=1]=3
  CurrentStatusValues$Quadrant[CurrentStatusValues$F_over_Fspr40<1 & CurrentStatusValues$B_over_Bspr20<1]=4
  StatusTable = data.frame(table(CurrentStatusValues$FamilyCommonName,CurrentStatusValues$WPP,CurrentStatusValues$Quadrant))
  names(StatusTable) = c("FamilyCommonName","WPP","Quadrant","Freq")
  StatusTable = reshape(StatusTable,v.names="Freq",idvar=c("FamilyCommonName","WPP"),timevar="Quadrant",direction="wide")
  StatusTable$Total = StatusTable$Freq.1 + StatusTable$Freq.2 + StatusTable$Freq.3 + StatusTable$Freq.4
  StatusTable = subset(StatusTable,StatusTable$Total>0)
  StatusTable$Percent.1 = round((StatusTable$Freq.1/StatusTable$Total)*100)
  StatusTable$Percent.2 = round((StatusTable$Freq.2/StatusTable$Total)*100)
  StatusTable$Percent.3 = round((StatusTable$Freq.3/StatusTable$Total)*100)
  StatusTable$Percent.4 = round((StatusTable$Freq.4/StatusTable$Total)*100)
  StatusTable[is.na(StatusTable)]=0
  StatusTable$FamilyCommonName = as.character(StatusTable$FamilyCommonName)
  StatusTable$WPP = as.character(StatusTable$WPP)
  write.table(StatusTable,paste(PATH_output,"StatusTableFamily.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
  #Now make the Kobe plots
  BenchmarksToPlot = subset(CurrentStatusValues,select=c(FamilyCommonName,WPP))
  BenchmarksToPlot = unique(BenchmarksToPlot)
  CurrentStatusValues$pch_arg=0
  CurrentStatusValues$pch_arg[CurrentStatusValues$populationReconstructionMethod=="LBSPR"]=15
  CurrentStatusValues$pch_arg[CurrentStatusValues$populationReconstructionMethod=="BevertonHoltInstantaneousMortality"]=16
  CurrentStatusValues$pch_arg[CurrentStatusValues$populationReconstructionMethod=="LIME"]=17
  CurrentStatusValues$pch_arg[CurrentStatusValues$populationReconstructionMethod=="TropFishR"]=18
  CurrentStatusValues$pch_arg[CurrentStatusValues$populationReconstructionMethod=="SolveBaranov"]=19
  ColorsFor_benchmarkPlot_CmsyScenarios = data.frame(CatchMSY_Scenario=sort(unique(CurrentStatusValues$CatchMSY_Scenario)),Colors=c("red","royalblue","green","purple","orange","brown"),stringsAsFactors=FALSE)
  LengthReconstructionSymbols = unique(subset(CurrentStatusValues,select=c(populationReconstructionMethod,pch_arg)))
  CurrentStatusValues = merge(CurrentStatusValues,ColorsFor_benchmarkPlot_CmsyScenarios,all.x=TRUE,by=c("CatchMSY_Scenario"))
  for(i in 1:dim(BenchmarksToPlot)[1])
  {
    thisIteration = CurrentStatusValues[CurrentStatusValues$FamilyCommonName==BenchmarksToPlot$FamilyCommonName[i] & CurrentStatusValues$WPP==BenchmarksToPlot$WPP[i],]
    png(paste(PATH_benchmarkPlotsFamily,paste(paste(BenchmarksToPlot$FamilyCommonName[i],BenchmarksToPlot$WPP[i],sep="_"),"_Benchmarks.png",sep=""),sep="/"),units="px",width=5800,height=3200,res=600)
    layout(matrix(c(1,2), nrow = 1), widths = c(1,0.5))
    par(mar = c(5, 4, 4, 2) + 0.1)
    maxBratio = max(thisIteration$B_over_Bspr20)
    maxFratio = max(thisIteration$F_over_Fspr40)
    if(maxBratio<2) {maxBratio=2}
    if(maxBratio>10) {maxBratio=10}
    if(maxFratio<2) {maxFratio=2}
    if(maxFratio>10) {maxFratio=10}
    plot(NULL,xlab="Bcurrent/B@SPR20",ylab="Fcurrent/F@SPR40",xlim=c(0,ceiling(maxBratio)),ylim=c(0,ceiling(maxFratio)),xaxs="i",yaxs="i")
    rect(xleft=0,ybottom=0,xright=1,ytop=1,par("usr")[3] - 1, y0 + 3, par("usr")[4] + 1,border="khaki",col="khaki")
    rect(xleft=1,ybottom=0,xright=ceiling(maxBratio),ytop=1,par("usr")[3] - 1, y0 + 3, par("usr")[4] + 1,border="lightgreen",col="lightgreen")
    rect(xleft=0,ybottom=1,xright=1,ytop=ceiling(maxFratio),par("usr")[3] - 1, y0 + 3, par("usr")[4] + 1,border="plum1",col="plum1")
    rect(xleft=1,ybottom=1,xright=ceiling(maxBratio),ytop=ceiling(maxFratio),par("usr")[3] - 1, y0 + 3, par("usr")[4] + 1,border="bisque",col="bisque") #tan1
    for(j in 1:dim(thisIteration)[1])
    {
      if(j==1)
      {
        points(thisIteration$B_over_Bspr20[j],thisIteration$F_over_Fspr40[j],xlab="Bcurrent/B@SPR20",ylab="Fcurrent/F@SPR40",pch=thisIteration$pch_arg[j],col=thisIteration$Colors[j])
      }
      if(j>1)
      {
        points(thisIteration$B_over_Bspr20[j],thisIteration$F_over_Fspr40[j],xlab="Bcurrent/B@SPR20",ylab="Fcurrent/F@SPR40",pch=thisIteration$pch_arg[j],col=thisIteration$Colors[j])
      }
    }
    if(ceiling(maxFratio)>3)
    {
      yMultiplier = ceiling(maxFratio)/18
    } else {
      yMultiplier = 0.1
    }
    if(ceiling(maxBratio)>3)
    {
      xMultiplier = ceiling(maxBratio)/18
    } else {
      xMultiplier = 0.1
    }
    StatusPercents = StatusTable[StatusTable$FamilyCommonName==BenchmarksToPlot$FamilyCommonName[i] & StatusTable$WPP==BenchmarksToPlot$WPP[i],]
    rect(xleft=0.5-xMultiplier,ybottom=(((ceiling(maxFratio)-1)/2)-yMultiplier)+1,xright=0.5+xMultiplier,ytop=(((ceiling(maxFratio)-1)/2)+yMultiplier)+1,par("usr")[3] - 1, y0 + 3, par("usr")[4] + 1,border="black",col="white",lty=1,lwd=1)
    text(x=0.5,y=((ceiling(maxFratio)-1)/2)+1,labels=paste(StatusPercents$Percent.1,"%",sep=""))
    rect(xleft=(((ceiling(maxBratio)-1)/2)-xMultiplier)+1,ybottom=(((ceiling(maxFratio)-1)/2)-yMultiplier)+1,xright=(((ceiling(maxBratio)-1)/2)+xMultiplier)+1,ytop=(((ceiling(maxFratio)-1)/2)+yMultiplier)+1,par("usr")[3] - 1, y0 + 3, par("usr")[4] + 1,border="black",col="white",lty=1,lwd=1)
    text(x=((ceiling(maxBratio)-1)/2)+1,y=((ceiling(maxFratio)-1)/2)+1,labels=paste(StatusPercents$Percent.2,"%",sep=""))
    rect(xleft=(((ceiling(maxBratio)-1)/2)-xMultiplier)+1,ybottom=0.5-yMultiplier,xright=(((ceiling(maxBratio)-1)/2)+xMultiplier)+1,ytop=0.5+yMultiplier,par("usr")[3] - 1, y0 + 3, par("usr")[4] + 1,border="black",col="white",lty=1,lwd=1)
    text(x=((ceiling(maxBratio)-1)/2)+1,y=0.5,labels=paste(StatusPercents$Percent.3,"%",sep=""))
    rect(xleft=0.5-xMultiplier,ybottom=0.5-yMultiplier,xright=0.5+xMultiplier,ytop=0.5+yMultiplier,par("usr")[3] - 1, y0 + 3, par("usr")[4] + 1,border="black",col="white",lty=1,lwd=1)
    text(x=0.5,y=0.5,labels=paste(StatusPercents$Percent.4,"%",sep=""))
    abline(v=1,col="black",lty=3,lwd=1)
    abline(h=1,col="black",lty=3,lwd=1)
    title(paste("Current Status: ",BenchmarksToPlot$FamilyCommonName[i],sep=" "))
    mtext(paste("WPP = ",BenchmarksToPlot$WPP[i],sep=""))
    par(mar = c(5, 0, 4, 2))
    plot(thisIteration$B_over_Bspr20[dim(thisIteration)[1]],thisIteration$F_over_Fspr40[dim(thisIteration)[1]],type="n",axes=FALSE,ann=FALSE,xlim=c(0,maxBratio),ylim=c(0,maxFratio))
    LengthReconstructionSymbols_forLegend = LengthReconstructionSymbols
    LengthReconstructionSymbols_forLegend$populationReconstructionMethod = as.character(LengthReconstructionSymbols_forLegend$populationReconstructionMethod)
    LengthReconstructionSymbols_forLegend$populationReconstructionMethod[LengthReconstructionSymbols_forLegend$populationReconstructionMethod=="BevertonHoltInstantaneousMortality"]="BevHoltF"
    legend("topleft",legend=c(as.character(LengthReconstructionSymbols_forLegend$populationReconstructionMethod),as.character(ColorsFor_benchmarkPlot_CmsyScenarios$CatchMSY_Scenario)),lty=c(rep(NA,times=dim(LengthReconstructionSymbols_forLegend)[1]),rep(1,times=dim(ColorsFor_benchmarkPlot_CmsyScenarios)[1])),lwd=c(rep(NA,times=dim(LengthReconstructionSymbols_forLegend)[1]),rep(2,times=dim(ColorsFor_benchmarkPlot_CmsyScenarios)[1])),col=c(rep("black",times=dim(LengthReconstructionSymbols_forLegend)[1]),ColorsFor_benchmarkPlot_CmsyScenarios$Colors),pch=c(LengthReconstructionSymbols_forLegend$pch_arg,rep(NA,times=dim(ColorsFor_benchmarkPlot_CmsyScenarios)[1])),cex=0.8)
    dev.off()
  }
}


################## Create Box Plots of F_over_Fspr40 and B_over_Bspr20 for each species and Family Nationally (EEZ) #######################
if(createBarPlots==TRUE)
{
	library(ggplot2)
	library(tidyverse)
	ExtrapolatedLandings = read.table(paste(PATH_output,"ExtrapolatedLandings.csv",sep="/"),header=TRUE,sep=",")
	SpeciesFamily = subset(ExtrapolatedLandings,select=c(species,FamilyCommonName))
	SpeciesFamily = unique(SpeciesFamily)
	CurrentStatusValues_notFiltered = read.table(paste(PATH_output,"Benchmarks_withFlags.csv",sep="/"),header=TRUE,sep=",",as.is=TRUE)
	CurrentStatusValues_notFiltered = merge(CurrentStatusValues_notFiltered,SpeciesFamily,all.x=TRUE,by=c("species"))
	WPPs = sort(unique(CurrentStatusValues_notFiltered$WPP))
	#Overfishing Box Plots by Species for Each WPP and the EEZ with the Target set at F@SPR40%. 
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio40_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (species), y = F_over_Fspr40)) +      #fct_infreq
		  geom_abline(intercept = 1, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("F/F@SPR40")+
		  xlab(NULL)+
		  guides(fill="none")+
		  ggtitle(label="Fcurrent/F@SPR40%",subtitle=SubTitles)+		  
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18),plot.title=element_text(family='', face='bold', colour='black', size=26,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("Overfishing_Fratio40_Species_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Overfished Box Plots by Species for Each WPP and the EEZ with the Limit set at B@SPR20% 
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio40_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (species), y = B_over_Bspr20)) +
		  geom_abline(intercept = 1, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("B/B@SPR20")+
		  xlab(NULL)+
		  guides(fill="none")+
		  ggtitle(label="Bcurrent/B@SPR20%",subtitle=SubTitles)+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18),plot.title=element_text(family='', face='bold', colour='black', size=26,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("Overfished_Bratio20_Species_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#SPR Box Plots by Species for Each WPP and the EEZ
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio40_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (species), y = SPR_current)) +
		  geom_abline(intercept = 0.4, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Current SPR")+
		  xlab(NULL)+
		  guides(fill="none")+
		  ggtitle(label="Current Length-Based SPR Estimates",subtitle=SubTitles)+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18),plot.title=element_text(family='', face='bold', colour='black', size=26,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("CurrentSPR_Species_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Percent Reduction in F by Species Needed to Rebuild to SPR 40%
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio40_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues$PercentReductionF_toSPR40 = (((CurrentStatusValues$F_current - CurrentStatusValues$Fspr40)/CurrentStatusValues$F_current)*100)*-1
	CurrentStatusValues$PercentReductionF_toSPR40[CurrentStatusValues$PercentReductionF_toSPR40>0]=NA
	CurrentStatusValues$PercentReductionF_toSPR40[CurrentStatusValues$B_over_Bspr20>=1]=NA
	CurrentStatusValues = subset(CurrentStatusValues,!is.na(CurrentStatusValues$PercentReductionF_toSPR40))
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (species), y = PercentReductionF_toSPR40)) +
		  geom_abline(intercept = 0, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = -50, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Percent Reduction in F")+
		  xlab(NULL)+
		  guides(fill="none")+
		  ggtitle(label="Percent Reduction When F=F@SPR40% to Rebuild to Limit (B@SPR20%)",subtitle=SubTitles)+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18), plot.title=element_text(family='', face='bold', colour='black', size=26,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("PercentReductionInFrebuildSPR40_Species_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Percent Reduction in F by Species Needed to Rebuild to SPR 30%
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio30_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues$PercentReductionF_toSPR30 = (((CurrentStatusValues$F_current - CurrentStatusValues$Fspr30)/CurrentStatusValues$F_current)*100)*-1
	CurrentStatusValues$PercentReductionF_toSPR30[CurrentStatusValues$PercentReductionF_toSPR30>0]=NA
	CurrentStatusValues$PercentReductionF_toSPR30[CurrentStatusValues$B_over_Bspr20>=1]=NA
	CurrentStatusValues = subset(CurrentStatusValues,!is.na(CurrentStatusValues$PercentReductionF_toSPR30))
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}		
		ggplot(CurrentStatusValues_Sub, aes(x = (species), y = PercentReductionF_toSPR30)) +
		  geom_abline(intercept = 0, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = -50, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Percent Reduction in F")+
		  xlab(NULL)+
		  ggtitle(label="Percent Reduction When F=F@SPR30% to Rebuild to Limit (B@SPR20%)",subtitle=SubTitles)+
		  guides(fill="none")+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18), plot.title=element_text(family='', face='bold', colour='black', size=22,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("PercentReductionInFrebuildSPR30_Species_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Percent Reduction in F by Species Needed to Rebuild to SPR 20%
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio20_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues$PercentReductionF_toSPR20 = (((CurrentStatusValues$F_current - CurrentStatusValues$Fspr20)/CurrentStatusValues$F_current)*100)*-1
	CurrentStatusValues$PercentReductionF_toSPR20[CurrentStatusValues$PercentReductionF_toSPR20>0]=NA
	CurrentStatusValues$PercentReductionF_toSPR20[CurrentStatusValues$B_over_Bspr20>=1]=NA
	CurrentStatusValues = subset(CurrentStatusValues,!is.na(CurrentStatusValues$PercentReductionF_toSPR20))
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (species), y = PercentReductionF_toSPR20)) +
		  geom_abline(intercept = 0, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = -50, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Percent Reduction in F")+
		  xlab(NULL)+
		  ggtitle(label="Percent Reduction When F=F@SPR20% to Rebuild to Limit (B@SPR20%)",subtitle=SubTitles)+
		  guides(fill="none")+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18), plot.title=element_text(family='', face='bold', colour='black', size=22,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))			
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("PercentReductionInFrebuildSPR20_Species_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Percent Reduction in F by Species Needed to Rebuild to MSY
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fmsy_ratio_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues$PercentReductionF_toMSY = (((CurrentStatusValues$F_current - CurrentStatusValues$Fmsy)/CurrentStatusValues$F_current)*100)*-1
	CurrentStatusValues$PercentReductionF_toMSY[CurrentStatusValues$PercentReductionF_toMSY>0]=NA
	CurrentStatusValues$PercentReductionF_toMSY[CurrentStatusValues$B_over_Bspr20>=1]=NA
	CurrentStatusValues = subset(CurrentStatusValues,!is.na(CurrentStatusValues$PercentReductionF_toMSY))
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (species), y = PercentReductionF_toMSY)) +
		  geom_abline(intercept = 0, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = -50, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Percent Reduction in F")+
		  xlab(NULL)+
		  ggtitle(label="Percent Reduction When F=F@MSY to Rebuild to Limit (B@SPR20%)",subtitle=SubTitles)+
		  guides(fill="none")+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18), plot.title=element_text(family='', face='bold', colour='black', size=22,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("PercentReductionInFrebuildMSY_Species_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Percent Reduction in F by Species Needed to Rebuild to MEY
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fmey_ratio_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues$PercentReductionF_toMEY = (((CurrentStatusValues$F_current - CurrentStatusValues$Fmey)/CurrentStatusValues$F_current)*100)*-1
	CurrentStatusValues$PercentReductionF_toMEY[CurrentStatusValues$PercentReductionF_toMEY>0]=NA
	CurrentStatusValues$PercentReductionF_toMEY[CurrentStatusValues$B_over_Bspr20>=1]=NA
	CurrentStatusValues = subset(CurrentStatusValues,!is.na(CurrentStatusValues$PercentReductionF_toMEY))
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (species), y = PercentReductionF_toMEY)) +
		  geom_abline(intercept = 0, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = -50, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Percent Reduction in F")+
		  xlab(NULL)+
		  ggtitle(label="Percent Reduction When F=F@MEY to Rebuild to Limit (B@SPR20%)",subtitle=SubTitles)+
		  guides(fill="none")+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18), plot.title=element_text(family='', face='bold', colour='black', size=22,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("PercentReductionInFrebuildMEY_Species_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Rebuilding Time to B@SPR20% by Species When Fishing At F@SPR40% - NOTE: the variable name "RebuildingTimeToSPR40" is not really accurate and is misleading for this.
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio40_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues$RebuildingTimeToSPR40 = CurrentStatusValues$RebuildTime_to_spr20_New.F_at_SPR40
	CurrentStatusValues$RebuildingTimeToSPR40[CurrentStatusValues$RebuildingTimeToSPR40<=0]=NA
	CurrentStatusValues$RebuildingTimeToSPR40[CurrentStatusValues$B_over_Bspr20>=1]=NA
	CurrentStatusValues = subset(CurrentStatusValues,!is.na(CurrentStatusValues$RebuildingTimeToSPR40))
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (species), y = RebuildingTimeToSPR40)) +
		  geom_abline(intercept = 40, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 30, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 20, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 10, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Rebuilding Time (Years)")+
		  xlab(NULL)+
		  ggtitle("Rebuilding Time to Limit (B@SPR20%) When F=F@SPR40%",subtitle=SubTitles)+
		  guides(fill="none")+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18), plot.title=element_text(family='', face='bold', colour='black', size=22,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("RebuildingTimeFishAtSPR40_Species_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Rebuilding Time to B@SPR20% by Species When Fishing At F@SPR30% - NOTE: the variable name "RebuildingTimeToSPR40" is not really accurate and is misleading for this.
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio30_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues$RebuildingTimeToSPR40 = CurrentStatusValues$RebuildTime_to_spr20_New.F_at_SPR30
	CurrentStatusValues$RebuildingTimeToSPR40[CurrentStatusValues$RebuildingTimeToSPR30<=0]=NA
	CurrentStatusValues$RebuildingTimeToSPR40[CurrentStatusValues$B_over_Bspr20>=1]=NA
	CurrentStatusValues = subset(CurrentStatusValues,!is.na(CurrentStatusValues$RebuildingTimeToSPR40))
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (species), y = RebuildingTimeToSPR40)) +
		  geom_abline(intercept = 40, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 30, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 20, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 10, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Rebuilding Time (Years)")+
		  xlab(NULL)+
		  ggtitle("Rebuilding Time to Limit (B@SPR20%) When F=F@SPR30%",subtitle=SubTitles)+
		  guides(fill="none")+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18), plot.title=element_text(family='', face='bold', colour='black', size=22,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("RebuildingTimeFishAtSPR30_Species_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Rebuilding Time to B@SPR20% by Species When Fishing At F@SPR20% - NOTE: the variable name "RebuildingTimeToSPR40" is not really accurate and is misleading for this.
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio20_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues$RebuildingTimeToSPR40 = CurrentStatusValues$RebuildTime_to_spr20_New.F_at_SPR20
	CurrentStatusValues$RebuildingTimeToSPR40[CurrentStatusValues$RebuildingTimeToSPR20<=0]=NA
	CurrentStatusValues$RebuildingTimeToSPR40[CurrentStatusValues$B_over_Bspr20>=1]=NA
	CurrentStatusValues = subset(CurrentStatusValues,!is.na(CurrentStatusValues$RebuildingTimeToSPR40))
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (species), y = RebuildingTimeToSPR40)) +
		  geom_abline(intercept = 40, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 30, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 20, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 10, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Rebuilding Time (Years)")+
		  xlab(NULL)+
		  ggtitle("Rebuilding Time to Limit (B@SPR20%) When F=F@SPR20%",subtitle=SubTitles)+
		  guides(fill="none")+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18), plot.title=element_text(family='', face='bold', colour='black', size=22,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("RebuildingTimeFishAtSPR20_Species_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Rebuilding Time to B@SPR20% by Species When Fishing At Fmey - NOTE: the variable name "RebuildingTimeToSPR40" is not really accurate and is misleading for this.
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fmey_ratio_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues$RebuildingTimeToSPR40 = CurrentStatusValues$RebuildTime_to_spr20_New.MEY
	CurrentStatusValues$RebuildingTimeToSPR40[CurrentStatusValues$RebuildingTimeToSPR20<=0]=NA
	CurrentStatusValues$RebuildingTimeToSPR40[CurrentStatusValues$B_over_Bspr20>=1]=NA
	CurrentStatusValues = subset(CurrentStatusValues,!is.na(CurrentStatusValues$RebuildingTimeToSPR40))
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (species), y = RebuildingTimeToSPR40)) +
		  geom_abline(intercept = 40, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 30, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 20, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 10, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Rebuilding Time (Years)")+
		  xlab(NULL)+
		  ggtitle("Rebuilding Time to Limit (B@SPR20%) When F=F@MEY",subtitle=SubTitles)+
		  guides(fill="none")+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18), plot.title=element_text(family='', face='bold', colour='black', size=22,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("RebuildingTimeFishAtMEY_Species_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Rebuilding Time to B@SPR20% by Species When Fishing At Fmsy - NOTE: the variable name "RebuildingTimeToSPR40" is not really accurate and is misleading for this.
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fmsy_ratio_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues$RebuildingTimeToSPR40 = CurrentStatusValues$RebuildTime_to_spr20_New.MSY
	CurrentStatusValues$RebuildingTimeToSPR40[CurrentStatusValues$RebuildingTimeToSPR20<=0]=NA
	CurrentStatusValues$RebuildingTimeToSPR40[CurrentStatusValues$B_over_Bspr20>=1]=NA
	CurrentStatusValues = subset(CurrentStatusValues,!is.na(CurrentStatusValues$RebuildingTimeToSPR40))
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (species), y = RebuildingTimeToSPR40)) +
		  geom_abline(intercept = 40, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 30, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 20, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 10, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Rebuilding Time (Years)")+
		  xlab(NULL)+
		  ggtitle("Rebuilding Time to Limit (B@SPR20%) When F=F@MSY",subtitle=SubTitles)+
		  guides(fill="none")+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18), plot.title=element_text(family='', face='bold', colour='black', size=22,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("RebuildingTimeFishAtMSY_Species_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	######## All of the above now but by FishFamily
	#Overfishing Box Plots by Fish Family for Each WPP and the EEZ with the Target set at F@SPR40%. 
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio40_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (FamilyCommonName), y = F_over_Fspr40)) +
		  geom_abline(intercept = 1, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("F/F@SPR40")+
		  xlab(NULL)+
		  guides(fill="none")+
		  ggtitle(label="Fcurrent/F@SPR40%",subtitle=SubTitles)+		  
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18),plot.title=element_text(family='', face='bold', colour='black', size=26,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("Overfishing_Fratio40_FamilyCommonName_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Overfished Box Plots by Fish Family for Each WPP and the EEZ with the Limit set at B@SPR20% 
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio40_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (FamilyCommonName), y = B_over_Bspr20)) +
		  geom_abline(intercept = 1, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("B/B@SPR20")+
		  xlab(NULL)+
		  guides(fill="none")+
		  ggtitle(label="Bcurrent/B@SPR20%",subtitle=SubTitles)+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18),plot.title=element_text(family='', face='bold', colour='black', size=26,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("Overfished_Bratio20_FamilyCommonName_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#SPR Box Plots by Fish Family for Each WPP and the EEZ
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio40_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (FamilyCommonName), y = SPR_current)) +
		  geom_abline(intercept = 0.4, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Current SPR")+
		  xlab(NULL)+
		  guides(fill="none")+
		  ggtitle(label="Current Length-Based SPR Estimates",subtitle=SubTitles)+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18),plot.title=element_text(family='', face='bold', colour='black', size=26,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("CurrentSPR_FamilyCommonName_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Percent Reduction in F by Fish Family Needed to Rebuild to SPR 40%
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio40_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues$PercentReductionF_toSPR40 = (((CurrentStatusValues$F_current - CurrentStatusValues$Fspr40)/CurrentStatusValues$F_current)*100)*-1
	CurrentStatusValues$PercentReductionF_toSPR40[CurrentStatusValues$PercentReductionF_toSPR40>0]=NA
	CurrentStatusValues$PercentReductionF_toSPR40[CurrentStatusValues$B_over_Bspr20>=1]=NA	
	CurrentStatusValues = subset(CurrentStatusValues,!is.na(CurrentStatusValues$PercentReductionF_toSPR40))
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (FamilyCommonName), y = PercentReductionF_toSPR40)) +
		  geom_abline(intercept = 0, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = -50, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Percent Reduction in F")+
		  xlab(NULL)+
		  guides(fill="none")+
		  ggtitle(label="Percent Reduction When F=F@SPR40% to Rebuild to Limit (B@SPR20%)",subtitle=SubTitles)+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18), plot.title=element_text(family='', face='bold', colour='black', size=26,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("PercentReductionInFrebuildSPR40_FamilyCommonName_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Percent Reduction in F by Fish Family Needed to Rebuild to SPR 30%
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio30_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues$PercentReductionF_toSPR30 = (((CurrentStatusValues$F_current - CurrentStatusValues$Fspr30)/CurrentStatusValues$F_current)*100)*-1
	CurrentStatusValues$PercentReductionF_toSPR30[CurrentStatusValues$PercentReductionF_toSPR30>0]=NA
	CurrentStatusValues$PercentReductionF_toSPR30[CurrentStatusValues$B_over_Bspr20>=1]=NA
	CurrentStatusValues = subset(CurrentStatusValues,!is.na(CurrentStatusValues$PercentReductionF_toSPR30))
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (FamilyCommonName), y = PercentReductionF_toSPR30)) +
		  geom_abline(intercept = 0, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = -50, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Percent Reduction in F")+
		  xlab(NULL)+
		  ggtitle(label="Percent Reduction When F=F@SPR30% to Rebuild to Limit (B@SPR20%)",subtitle=SubTitles)+
		  guides(fill="none")+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18), plot.title=element_text(family='', face='bold', colour='black', size=22,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("PercentReductionInFrebuildSPR30_FamilyCommonName_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Percent Reduction in F by Fish Family Needed to Rebuild to SPR 20%
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio20_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues$PercentReductionF_toSPR20 = (((CurrentStatusValues$F_current - CurrentStatusValues$Fspr20)/CurrentStatusValues$F_current)*100)*-1
	CurrentStatusValues$PercentReductionF_toSPR20[CurrentStatusValues$PercentReductionF_toSPR20>0]=NA
	CurrentStatusValues$PercentReductionF_toSPR20[CurrentStatusValues$B_over_Bspr20>=1]=NA
	CurrentStatusValues = subset(CurrentStatusValues,!is.na(CurrentStatusValues$PercentReductionF_toSPR20))
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (FamilyCommonName), y = PercentReductionF_toSPR20)) +
		  geom_abline(intercept = 0, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = -50, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Percent Reduction in F")+
		  xlab(NULL)+
		  ggtitle(label="Percent Reduction When F=F@SPR20% to Rebuild to Limit (B@SPR20%)",subtitle=SubTitles)+
		  guides(fill="none")+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18), plot.title=element_text(family='', face='bold', colour='black', size=22,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))			
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("PercentReductionInFrebuildSPR20_FamilyCommonName_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Percent Reduction in F by Fish Family Needed to Rebuild to MSY
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fmsy_ratio_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues$PercentReductionF_toMSY = (((CurrentStatusValues$F_current - CurrentStatusValues$Fmsy)/CurrentStatusValues$F_current)*100)*-1
	CurrentStatusValues$PercentReductionF_toMSY[CurrentStatusValues$PercentReductionF_toMSY>0]=NA
	CurrentStatusValues$PercentReductionF_toMSY[CurrentStatusValues$B_over_Bspr20>=1]=NA
	CurrentStatusValues = subset(CurrentStatusValues,!is.na(CurrentStatusValues$PercentReductionF_toMSY))
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (FamilyCommonName), y = PercentReductionF_toMSY)) +
		  geom_abline(intercept = 0, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = -50, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Percent Reduction in F")+
		  xlab(NULL)+
		  ggtitle(label="Percent Reduction When F=F@MSY to Rebuild to Limit (B@SPR20%)",subtitle=SubTitles)+
		  guides(fill="none")+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18), plot.title=element_text(family='', face='bold', colour='black', size=22,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("PercentReductionInFrebuildMSY_FamilyCommonName_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Percent Reduction in F by Fish Family Needed to Rebuild to MEY
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fmey_ratio_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues$PercentReductionF_toMEY = (((CurrentStatusValues$F_current - CurrentStatusValues$Fmey)/CurrentStatusValues$F_current)*100)*-1
	CurrentStatusValues$PercentReductionF_toMEY[CurrentStatusValues$PercentReductionF_toMEY>0]=NA
	CurrentStatusValues$PercentReductionF_toMEY[CurrentStatusValues$B_over_Bspr20>=1]=NA
	CurrentStatusValues = subset(CurrentStatusValues,!is.na(CurrentStatusValues$PercentReductionF_toMEY))
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (FamilyCommonName), y = PercentReductionF_toMEY)) +
		  geom_abline(intercept = 0, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = -50, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Percent Reduction in F")+
		  xlab(NULL)+
		  ggtitle(label="Percent Reduction When F=F@MEY to Rebuild to Limit (B@SPR20%)",subtitle=SubTitles)+
		  guides(fill="none")+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18), plot.title=element_text(family='', face='bold', colour='black', size=22,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("PercentReductionInFrebuildMEY_FamilyCommonName_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Rebuilding Time to B@SPR20% by Fish Family When Fishing At F@SPR40% - NOTE: the variable name "RebuildingTimeToSPR40" is not really accurate and is misleading for this.
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio40_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues$RebuildingTimeToSPR40 = CurrentStatusValues$RebuildTime_to_spr20_New.F_at_SPR40
	CurrentStatusValues$RebuildingTimeToSPR40[CurrentStatusValues$RebuildingTimeToSPR40<=0]=NA
	CurrentStatusValues$RebuildingTimeToSPR40[CurrentStatusValues$B_over_Bspr20>=1]=NA
	CurrentStatusValues = subset(CurrentStatusValues,!is.na(CurrentStatusValues$RebuildingTimeToSPR40))
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (FamilyCommonName), y = RebuildingTimeToSPR40)) +
		  geom_abline(intercept = 40, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 30, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 20, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 10, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Rebuilding Time (Years)")+
		  xlab(NULL)+
		  ggtitle("Rebuilding Time to Limit (B@SPR20%) When F=F@SPR40%",subtitle=SubTitles)+
		  guides(fill="none")+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18), plot.title=element_text(family='', face='bold', colour='black', size=22,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("RebuildingTimeFishAtSPR40_FamilyCommonName_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Rebuilding Time to B@SPR20% by Fish Family When Fishing At F@SPR30% - NOTE: the variable name "RebuildingTimeToSPR40" is not really accurate and is misleading for this.
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio30_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues$RebuildingTimeToSPR40 = CurrentStatusValues$RebuildTime_to_spr20_New.F_at_SPR30
	CurrentStatusValues$RebuildingTimeToSPR40[CurrentStatusValues$RebuildingTimeToSPR30<=0]=NA
	CurrentStatusValues$RebuildingTimeToSPR40[CurrentStatusValues$B_over_Bspr20>=1]=NA
	CurrentStatusValues = subset(CurrentStatusValues,!is.na(CurrentStatusValues$RebuildingTimeToSPR40))
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (FamilyCommonName), y = RebuildingTimeToSPR40)) +
		  geom_abline(intercept = 40, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 30, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 20, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 10, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Rebuilding Time (Years)")+
		  xlab(NULL)+
		  ggtitle("Rebuilding Time to Limit (B@SPR20%) When F=F@SPR30%",subtitle=SubTitles)+
		  guides(fill="none")+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18), plot.title=element_text(family='', face='bold', colour='black', size=22,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("RebuildingTimeFishAtSPR30_FamilyCommonName_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Rebuilding Time to B@SPR20% by Fish Family When Fishing At F@SPR20% - NOTE: the variable name "RebuildingTimeToSPR40" is not really accurate and is misleading for this.
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio20_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues$RebuildingTimeToSPR40 = CurrentStatusValues$RebuildTime_to_spr20_New.F_at_SPR20
	CurrentStatusValues$RebuildingTimeToSPR40[CurrentStatusValues$RebuildingTimeToSPR20<=0]=NA
	CurrentStatusValues$RebuildingTimeToSPR40[CurrentStatusValues$B_over_Bspr20>=1]=NA
	CurrentStatusValues = subset(CurrentStatusValues,!is.na(CurrentStatusValues$RebuildingTimeToSPR40))
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (FamilyCommonName), y = RebuildingTimeToSPR40)) +
		  geom_abline(intercept = 40, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 30, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 20, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 10, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Rebuilding Time (Years)")+
		  xlab(NULL)+
		  ggtitle("Rebuilding Time to Limit (B@SPR20%) When F=F@SPR20%",subtitle=SubTitles)+
		  guides(fill="none")+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18), plot.title=element_text(family='', face='bold', colour='black', size=22,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("RebuildingTimeFishAtSPR20_FamilyCommonName_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Rebuilding Time to B@SPR20% by Fish Family When Fishing At Fmey - NOTE: the variable name "RebuildingTimeToSPR40" is not really accurate and is misleading for this.
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fmey_ratio_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues$RebuildingTimeToSPR40 = CurrentStatusValues$RebuildTime_to_spr20_New.MEY
	CurrentStatusValues$RebuildingTimeToSPR40[CurrentStatusValues$RebuildingTimeToSPR20<=0]=NA
	CurrentStatusValues$RebuildingTimeToSPR40[CurrentStatusValues$B_over_Bspr20>=1]=NA
	CurrentStatusValues = subset(CurrentStatusValues,!is.na(CurrentStatusValues$RebuildingTimeToSPR40))
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (FamilyCommonName), y = RebuildingTimeToSPR40)) +
		  geom_abline(intercept = 40, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 30, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 20, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 10, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Rebuilding Time (Years)")+
		  xlab(NULL)+
		  ggtitle("Rebuilding Time to Limit (B@SPR20%) When F=F@MEY",subtitle=SubTitles)+
		  guides(fill="none")+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18), plot.title=element_text(family='', face='bold', colour='black', size=22,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("RebuildingTimeFishAtMEY_FamilyCommonName_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Rebuilding Time to B@SPR20% by Fish Family When Fishing At Fmsy - NOTE: the variable name "RebuildingTimeToSPR40" is not really accurate and is misleading for this.
	CurrentStatusValues = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fmsy_ratio_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues$RebuildingTimeToSPR40 = CurrentStatusValues$RebuildTime_to_spr20_New.MSY
	CurrentStatusValues$RebuildingTimeToSPR40[CurrentStatusValues$RebuildingTimeToSPR20<=0]=NA
	CurrentStatusValues$RebuildingTimeToSPR40[CurrentStatusValues$B_over_Bspr20>=1]=NA
	CurrentStatusValues = subset(CurrentStatusValues,!is.na(CurrentStatusValues$RebuildingTimeToSPR40))
	for(i in 1:length(WPPs))
	{
		CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}	
		ggplot(CurrentStatusValues_Sub, aes(x = (FamilyCommonName), y = RebuildingTimeToSPR40)) +
		  geom_abline(intercept = 40, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 30, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 20, slope = 0, colour = "dark red", linetype = 4)+
		  geom_abline(intercept = 10, slope = 0, colour = "dark red", linetype = 4)+
		  geom_boxplot(fill = "#487eb0", alpha = 0.3) +
		  theme_classic()+
		  ylab("Rebuilding Time (Years)")+
		  xlab(NULL)+
		  ggtitle("Rebuilding Time to Limit (B@SPR20%) When F=F@MSY",subtitle=SubTitles)+
		  guides(fill="none")+
		  theme(axis.title = element_text(size=16), axis.text = element_text(size = 18), 
				axis.text.x = element_text(angle = 90,  hjust=1), legend.title=element_text(size=16), 
				legend.text=element_text(size=18), plot.title=element_text(family='', face='bold', colour='black', size=22,hjust=0.5),
				plot.subtitle=element_text(family='',face="italic",colour="black",size=22,hjust=0.5))
		ggsave(file = paste(PATH_highLevelSummaryPlots,paste("RebuildingTimeFishAtMSY_FamilyCommonName_",WPPs[i],".jpg",sep=""),sep="/"), width = 18, height = 8)
	}
	#Rebuilding Time Trade-Off Plot by fish Family
	CurrentStatusValues_SPR20 = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio20_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues_SPR20 = subset(CurrentStatusValues_SPR20,select=c(species,FamilyCommonName,WPP,CatchMSY_Scenario,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,RebuildTime_to_spr20_New.F_at_SPR20))
	CurrentStatusValues_SPR20$BenchmarkType = "SPR20%"
	names(CurrentStatusValues_SPR20)[names(CurrentStatusValues_SPR20)=="RebuildTime_to_spr20_New.F_at_SPR20"]="RebuildingTime"
	CurrentStatusValues_SPR20$RebuildingTime[CurrentStatusValues_SPR20$RebuildingTime<=0]=NA
	CurrentStatusValues_SPR20$RebuildingTime[CurrentStatusValues_SPR20$B_over_Bspr20>=1]=NA
	CurrentStatusValues_SPR20 = subset(CurrentStatusValues_SPR20,!is.na(CurrentStatusValues_SPR20$RebuildingTime))
	CurrentStatusValues_SPR20$RebuildTime_to_spr20_New.F_at_SPR20=NULL
	CurrentStatusValues_SPR30 = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio30_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues_SPR30 = subset(CurrentStatusValues_SPR30,select=c(species,FamilyCommonName,WPP,CatchMSY_Scenario,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,RebuildTime_to_spr20_New.F_at_SPR30))
	CurrentStatusValues_SPR30$BenchmarkType = "SPR30%"
	names(CurrentStatusValues_SPR30)[names(CurrentStatusValues_SPR30)=="RebuildTime_to_spr20_New.F_at_SPR30"]="RebuildingTime"
	CurrentStatusValues_SPR30$RebuildingTime[CurrentStatusValues_SPR30$RebuildingTime<=0]=NA
	CurrentStatusValues_SPR30$RebuildingTime[CurrentStatusValues_SPR30$B_over_Bspr20>=1]=NA
	CurrentStatusValues_SPR30 = subset(CurrentStatusValues_SPR30,!is.na(CurrentStatusValues_SPR30$RebuildingTime))
	CurrentStatusValues_SPR30$RebuildTime_to_spr20_New.F_at_SPR30=NULL
	CurrentStatusValues_SPR40 = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio40_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues_SPR40 = subset(CurrentStatusValues_SPR40,select=c(species,FamilyCommonName,WPP,CatchMSY_Scenario,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,RebuildTime_to_spr20_New.F_at_SPR40))
	CurrentStatusValues_SPR40$BenchmarkType = "SPR40%"
	names(CurrentStatusValues_SPR40)[names(CurrentStatusValues_SPR40)=="RebuildTime_to_spr20_New.F_at_SPR40"]="RebuildingTime"
	CurrentStatusValues_SPR40$RebuildingTime[CurrentStatusValues_SPR40$RebuildingTime<=0]=NA
	CurrentStatusValues_SPR40$RebuildingTime[CurrentStatusValues_SPR40$B_over_Bspr20>=1]=NA
	CurrentStatusValues_SPR40 = subset(CurrentStatusValues_SPR40,!is.na(CurrentStatusValues_SPR40$RebuildingTime))
	CurrentStatusValues_SPR40$RebuildTime_to_spr20_New.F_at_SPR40=NULL
	CurrentStatusValues_MSY = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fmsy_ratio_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues_MSY = subset(CurrentStatusValues_MSY,select=c(species,FamilyCommonName,WPP,CatchMSY_Scenario,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,RebuildTime_to_spr20_New.MSY))
	CurrentStatusValues_MSY$BenchmarkType = "MSY"
	names(CurrentStatusValues_MSY)[names(CurrentStatusValues_MSY)=="RebuildTime_to_spr20_New.MSY"]="RebuildingTime"
	CurrentStatusValues_MSY$RebuildingTime[CurrentStatusValues_MSY$RebuildingTime<=0]=NA
	CurrentStatusValues_MSY$RebuildingTime[CurrentStatusValues_MSY$B_over_Bspr20>=1]=NA
	CurrentStatusValues_MSY = subset(CurrentStatusValues_MSY,!is.na(CurrentStatusValues_MSY$RebuildingTime))
	CurrentStatusValues_MSY$RebuildTime_to_spr20_New.MSY=NULL
	CurrentStatusValues_MEY = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fmey_ratio_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues_MEY = subset(CurrentStatusValues_MEY,select=c(species,FamilyCommonName,WPP,CatchMSY_Scenario,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,RebuildTime_to_spr20_New.MEY))
	CurrentStatusValues_MEY$BenchmarkType = "MEY"
	names(CurrentStatusValues_MEY)[names(CurrentStatusValues_MEY)=="RebuildTime_to_spr20_New.MEY"]="RebuildingTime"
	CurrentStatusValues_MEY$RebuildingTime[CurrentStatusValues_MEY$RebuildingTime<=0]=NA
	CurrentStatusValues_MEY$RebuildingTime[CurrentStatusValues_MEY$B_over_Bspr20>=1]=NA
	CurrentStatusValues_MEY = subset(CurrentStatusValues_MEY,!is.na(CurrentStatusValues_MEY$RebuildingTime))
	CurrentStatusValues_MEY$RebuildTime_to_spr20_New.MEY=NULL
	CurrentStatusValues = rbind(CurrentStatusValues_SPR20,CurrentStatusValues_SPR30,CurrentStatusValues_SPR40,
		CurrentStatusValues_MSY,CurrentStatusValues_MEY)
	CheckForZeros = data.frame(table(CurrentStatusValues$FamilyCommonName,CurrentStatusValues$BenchmarkType,CurrentStatusValues$WPP))
	names(CheckForZeros) = c("FamilyCommonName","BenchmarkType","WPP","Freq")
	CheckForZeros$Freq[CheckForZeros$Freq>0]=1
	for(i in 1:length(WPPs))
	{
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}
		CheckForZeros_i = subset(CheckForZeros,CheckForZeros$WPP==WPPs[i])
		Zeros = subset(CheckForZeros_i,CheckForZeros_i$Freq==0)
		if(dim(Zeros)[1]>0)
		{
			Zeros$species=""
			Zeros$CatchMSY_Scenario=""
			Zeros$LifeHistory_Scenario=""
			Zeros$populationReconstructionMethod=""
			Zeros$SteepnessAssumed=NA
			names(Zeros)[names(Zeros)=="Freq"]="RebuildingTime"
			CurrentStatusValues_Sub = subset(CurrentStatusValues,CurrentStatusValues$WPP==WPPs[i])
			Zeros = Zeros[,names(CurrentStatusValues_Sub)]
			CurrentStatusValues_Sub = rbind(CurrentStatusValues_Sub,Zeros)
		}
		png(paste(PATH_highLevelSummaryPlots,paste("RebuildingTime_TradeOff_",WPPs[i],".png",sep=""),sep="/"),units="px",width=5200,height=3200,res=600)
		par(mar=c(10, 5, 5, 5), xpd=TRUE)	
		cols = c("royalblue2","steelblue2","turquoise1","powderblue","aquamarine") # rainbow(5,s=0.5)
		boxplot(RebuildingTime ~ BenchmarkType + FamilyCommonName, data = CurrentStatusValues_Sub,
				at = c(1:5,7:11,13:17,19:23,25:29,31:35), col = cols,range=0.5,outline=FALSE,xaxt="n",
				xlab="Family Common Name",ylab="Rebuilding Time (years)",ylim=c(0,25))
		axis(1,at=c(3,9,15,21,27,33),labels=c("Drums/Croaker","Emperor","Grouper","Grunts","Other","Snapper"))
		legend(3.5,-12, fill = cols, legend = c("MEY","MSY","SPR20%","SPR30%","SPR40%"), horiz = TRUE)
		title("Rebuilding Time Trade-Offs")
		mtext(SubTitles)
		dev.off()
	}
	#Percent Reduction in F Trade-Off Plot by Fish Family
	CurrentStatusValues_reductionToFspr20 = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio20_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues_reductionToFspr20$PercentReductionF_toSPR20 = (((CurrentStatusValues_reductionToFspr20$F_current - CurrentStatusValues_reductionToFspr20$Fspr20)/CurrentStatusValues_reductionToFspr20$F_current)*100)*-1
	CurrentStatusValues_reductionToFspr20$PercentReductionF_toSPR20[CurrentStatusValues_reductionToFspr20$PercentReductionF_toSPR20>=0]=NA
	CurrentStatusValues_reductionToFspr20$PercentReductionF_toSPR20[CurrentStatusValues_reductionToFspr20$B_over_Bspr20>=1]=NA
	CurrentStatusValues_reductionToFspr20 = subset(CurrentStatusValues_reductionToFspr20,select=c(species,FamilyCommonName,WPP,CatchMSY_Scenario,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,PercentReductionF_toSPR20))
	CurrentStatusValues_reductionToFspr20$BenchmarkType="SPR20%"
	names(CurrentStatusValues_reductionToFspr20)[names(CurrentStatusValues_reductionToFspr20)=="PercentReductionF_toSPR20"]="PercentReductionF"
	CurrentStatusValues_reductionToFspr20 = subset(CurrentStatusValues_reductionToFspr20,!is.na(CurrentStatusValues_reductionToFspr20$PercentReductionF))
	CurrentStatusValues_reductionToFspr30 = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio30_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues_reductionToFspr30$PercentReductionF_toSPR30 = (((CurrentStatusValues_reductionToFspr30$F_current - CurrentStatusValues_reductionToFspr30$Fspr30)/CurrentStatusValues_reductionToFspr30$F_current)*100)*-1
	CurrentStatusValues_reductionToFspr30$PercentReductionF_toSPR30[CurrentStatusValues_reductionToFspr30$PercentReductionF_toSPR30>=0]=NA
	CurrentStatusValues_reductionToFspr30$PercentReductionF_toSPR30[CurrentStatusValues_reductionToFspr30$B_over_Bspr20>=1]=NA
	CurrentStatusValues_reductionToFspr30 = subset(CurrentStatusValues_reductionToFspr30,select=c(species,FamilyCommonName,WPP,CatchMSY_Scenario,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,PercentReductionF_toSPR30))
	CurrentStatusValues_reductionToFspr30$BenchmarkType="SPR30%"
	names(CurrentStatusValues_reductionToFspr30)[names(CurrentStatusValues_reductionToFspr30)=="PercentReductionF_toSPR30"]="PercentReductionF"
	CurrentStatusValues_reductionToFspr30 = subset(CurrentStatusValues_reductionToFspr30,!is.na(CurrentStatusValues_reductionToFspr30$PercentReductionF))
	CurrentStatusValues_reductionToFspr40 = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fratio40_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues_reductionToFspr40$PercentReductionF_toSPR40 = (((CurrentStatusValues_reductionToFspr40$F_current - CurrentStatusValues_reductionToFspr40$Fspr40)/CurrentStatusValues_reductionToFspr40$F_current)*100)*-1
	CurrentStatusValues_reductionToFspr40$PercentReductionF_toSPR40[CurrentStatusValues_reductionToFspr40$PercentReductionF_toSPR40>=0]=NA
	CurrentStatusValues_reductionToFspr40$PercentReductionF_toSPR40[CurrentStatusValues_reductionToFspr40$B_over_Bspr20>=1]=NA
	CurrentStatusValues_reductionToFspr40 = subset(CurrentStatusValues_reductionToFspr40,select=c(species,FamilyCommonName,WPP,CatchMSY_Scenario,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,PercentReductionF_toSPR40))
	CurrentStatusValues_reductionToFspr40$BenchmarkType="SPR40%"
	names(CurrentStatusValues_reductionToFspr40)[names(CurrentStatusValues_reductionToFspr40)=="PercentReductionF_toSPR40"]="PercentReductionF"
	CurrentStatusValues_reductionToFspr40 = subset(CurrentStatusValues_reductionToFspr40,!is.na(CurrentStatusValues_reductionToFspr40$PercentReductionF))
	CurrentStatusValues_reductionToMSY = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fmsy_ratio_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues_reductionToMSY$PercentReductionF_MSY = (((CurrentStatusValues_reductionToMSY$F_current - CurrentStatusValues_reductionToMSY$Fmsy)/CurrentStatusValues_reductionToMSY$F_current)*100)*-1
	CurrentStatusValues_reductionToMSY$PercentReductionF_MSY[CurrentStatusValues_reductionToMSY$PercentReductionF_MSY>=0]=NA
	CurrentStatusValues_reductionToMSY$PercentReductionF_MSY[CurrentStatusValues_reductionToMSY$B_over_Bspr20>=1]=NA
	CurrentStatusValues_reductionToMSY = subset(CurrentStatusValues_reductionToMSY,select=c(species,FamilyCommonName,WPP,CatchMSY_Scenario,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,PercentReductionF_MSY))
	CurrentStatusValues_reductionToMSY$BenchmarkType="MSY"
	names(CurrentStatusValues_reductionToMSY)[names(CurrentStatusValues_reductionToMSY)=="PercentReductionF_MSY"]="PercentReductionF"
	CurrentStatusValues_reductionToMSY = subset(CurrentStatusValues_reductionToMSY,!is.na(CurrentStatusValues_reductionToMSY$PercentReductionF))
	CurrentStatusValues_reductionToMEY = subset(CurrentStatusValues_notFiltered,CurrentStatusValues_notFiltered$Fmey_ratio_FLAG==0 & CurrentStatusValues_notFiltered$Bratio20_FLAG==0)
	CurrentStatusValues_reductionToMEY$PercentReductionF_MEY = (((CurrentStatusValues_reductionToMEY$F_current - CurrentStatusValues_reductionToMEY$Fmey)/CurrentStatusValues_reductionToMEY$F_current)*100)*-1
	CurrentStatusValues_reductionToMEY$PercentReductionF_MEY[CurrentStatusValues_reductionToMEY$PercentReductionF_MEY>=0]=NA
	CurrentStatusValues_reductionToMEY$PercentReductionF_MEY[CurrentStatusValues_reductionToMEY$B_over_Bspr20>=1]=NA
	CurrentStatusValues_reductionToMEY = subset(CurrentStatusValues_reductionToMEY,select=c(species,FamilyCommonName,WPP,CatchMSY_Scenario,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,PercentReductionF_MEY))
	CurrentStatusValues_reductionToMEY$BenchmarkType="MEY"
	names(CurrentStatusValues_reductionToMEY)[names(CurrentStatusValues_reductionToMEY)=="PercentReductionF_MEY"]="PercentReductionF"
	CurrentStatusValues_reductionToMEY = subset(CurrentStatusValues_reductionToMEY,!is.na(CurrentStatusValues_reductionToMSY$PercentReductionF))
	CurrentStatusValues_reductionInF = rbind(CurrentStatusValues_reductionToFspr20,CurrentStatusValues_reductionToFspr30,
		CurrentStatusValues_reductionToFspr40,CurrentStatusValues_reductionToMSY,CurrentStatusValues_reductionToMEY)
	CheckForZeros = data.frame(table(CurrentStatusValues_reductionInF$FamilyCommonName,CurrentStatusValues_reductionInF$BenchmarkType,CurrentStatusValues_reductionInF$WPP))
	names(CheckForZeros) = c("FamilyCommonName","BenchmarkType","WPP","Freq")
	CheckForZeros$Freq[CheckForZeros$Freq>0]=1
	for(i in 1:length(WPPs))
	{
		if(WPPs[i]!="EEZ")
		{
			SubTitles = paste("Fisheries Management Area:",WPPs[i],sep=" ")
		}
		if(WPPs[i]=="EEZ")
		{
			SubTitles = "All FMAs Combined"
		}
		CheckForZeros_i = subset(CheckForZeros,CheckForZeros$WPP==WPPs[i])
		Zeros = subset(CheckForZeros_i,CheckForZeros_i$Freq==0)
		if(dim(Zeros)[1]>0)
		{
			Zeros$species=""
			Zeros$CatchMSY_Scenario=""
			Zeros$LifeHistory_Scenario=""
			Zeros$populationReconstructionMethod=""
			Zeros$SteepnessAssumed=NA
			names(Zeros)[names(Zeros)=="Freq"]="PercentReductionF"
			CurrentStatusValues_reductionInF_Sub = subset(CurrentStatusValues_reductionInF,CurrentStatusValues_reductionInF$WPP==WPPs[i])
			Zeros = Zeros[,names(CurrentStatusValues_reductionInF_Sub)]
			CurrentStatusValues_reductionInF_Sub = rbind(CurrentStatusValues_reductionInF_Sub,Zeros)
		}
		png(paste(PATH_highLevelSummaryPlots,paste("PercentReductionF_TradeOff_",WPPs[i],".png",sep=""),sep="/"),units="px",width=5200,height=3200,res=600)
		par(mar=c(10, 5, 5, 5), xpd=TRUE)	
		cols = c("royalblue2","steelblue2","turquoise1","powderblue","aquamarine") # rainbow(5,s=0.5)
		BoxPlotFigData = boxplot(PercentReductionF ~ BenchmarkType + FamilyCommonName, data = CurrentStatusValues_reductionInF_Sub,
				at = c(1:5,7:11,13:17,19:23,25:29,31:35), col = cols,range=0.5,outline=FALSE,xaxt="n",
				xlab="Family Common Name",ylab="Percent Reduction in F",ylim=c(-80,0))
		axis(1,at=c(3,9,15,21,27,33),labels=c("Drums/Croaker","Emperor","Grouper","Grunts","Other","Snapper"))
		legend(3.5,-125, fill = cols, legend = c("MEY","MSY","SPR20%","SPR30%","SPR40%"), horiz = TRUE)
		title("Percent Reduction in F Trade-Offs")
		mtext(SubTitles)
		dev.off()
	}
	
}

################# Develop Output Tables For Paper of (1) Current Status, and (2) Tradeoffs #######################
if(createFinalSummaryTables==TRUE)
{
	BenchmarkSummary = read.table(paste(PATH_output,"BenchmarkSummary.csv",sep="/"),header=TRUE,sep=",")
	#(1a) Current Values 
	Current_F_B = subset(BenchmarkSummary,select=c(species,WPP,F_current_median,F_current_UpperCI,
		F_current_LowerCI,B_current_geoMean,B_current_UpperCI,B_current_LowerCI))
	Current_F_B$F_current_range = paste(round(Current_F_B$F_current_LowerCI,2),round(Current_F_B$F_current_UpperCI,2),sep="-")
	Current_F_B$B_current_range = paste(prettyNum(round(Current_F_B$B_current_LowerCI,0),big.mark=",",scientific=FALSE),prettyNum(round(Current_F_B$B_current_UpperCI,0),big.mark=",",scientific=FALSE),sep=" - ")
	Current_F_B = subset(Current_F_B,select=c(species,WPP,F_current_median,F_current_range,B_current_geoMean,B_current_range))
	Current_F_B$F_current_median = round(Current_F_B$F_current_median,2)
	Current_F_B$B_current_geoMean = prettyNum(round(Current_F_B$B_current_geoMean,0),big.mark=",",scientific=FALSE)
	names(Current_F_B)[names(Current_F_B)=="species"]="Species"
	names(Current_F_B)[names(Current_F_B)=="WPP"]="FMA"
	write.table(Current_F_B,paste(PATH_highLevelSummaryPlots,"Current_F_B.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
	#(1b) Status Tables
	StatusRatios = subset(BenchmarkSummary,select=c(species,WPP,Fratio40_median,Fratio40_UpperCI,
		Fratio40_LowerCI,B_over_Bspr20_median,Bratio20_UpperCI,Bratio20_LowerCI))
	StatusRatios$Fratio40_range = paste(round(StatusRatios$Fratio40_LowerCI,2),round(StatusRatios$Fratio40_UpperCI,2),sep=" - ")
	StatusRatios$Bratio20_range = paste(prettyNum(round(StatusRatios$Bratio20_LowerCI,2),big.mark=",",scientific=FALSE),prettyNum(round(StatusRatios$Bratio20_UpperCI,2),big.mark=",",scientific=FALSE),sep=" - ")
	StatusRatios = subset(StatusRatios,select=c(species,WPP,Fratio40_median,Fratio40_range,B_over_Bspr20_median,Bratio20_range))
	StatusRatios$Fratio40_median = round(StatusRatios$Fratio40_median,2)
	StatusRatios$B_over_Bspr20_median = prettyNum(round(StatusRatios$B_over_Bspr20_median,2),big.mark=",",scientific=FALSE)
	names(StatusRatios)[names(StatusRatios)=="species"]="Species"
	names(StatusRatios)[names(StatusRatios)=="WPP"]="FMA"
	write.table(StatusRatios,paste(PATH_highLevelSummaryPlots,"StatusRatios.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
	#(2) Create F-ratio Tradeoff Tables 
	TradeOffs_FishingMortality = subset(BenchmarkSummary,select=c(species,WPP,Fratio20_median,Fratio30_median,Fratio40_median,FratioMSY_median,FratioMEY_median))
	TradeOffs_FishingMortality$Fratio20_median = round(TradeOffs_FishingMortality$Fratio20_median,2)
	TradeOffs_FishingMortality$Fratio30_median = round(TradeOffs_FishingMortality$Fratio30_median,2)
	TradeOffs_FishingMortality$Fratio40_median = round(TradeOffs_FishingMortality$Fratio40_median,2)
	TradeOffs_FishingMortality$FratioMSY_median = round(TradeOffs_FishingMortality$FratioMSY_median,2)
	TradeOffs_FishingMortality$FratioMEY_median = round(TradeOffs_FishingMortality$FratioMEY_median,2)
	names(TradeOffs_FishingMortality)[names(TradeOffs_FishingMortality)=="Fratio20_median"]="Fcurrent/F@SPR20%"
	names(TradeOffs_FishingMortality)[names(TradeOffs_FishingMortality)=="Fratio30_median"]="Fcurrent/F@SPR30%"
	names(TradeOffs_FishingMortality)[names(TradeOffs_FishingMortality)=="Fratio40_median"]="Fcurrent/F@SPR40%"
	names(TradeOffs_FishingMortality)[names(TradeOffs_FishingMortality)=="FratioMSY_median"]="Fcurrent/F@MSY"
	names(TradeOffs_FishingMortality)[names(TradeOffs_FishingMortality)=="FratioMEY_median"]="Fcurrent/F@MEY"
	names(TradeOffs_FishingMortality)[names(TradeOffs_FishingMortality)=="species"]="Species"
	names(TradeOffs_FishingMortality)[names(TradeOffs_FishingMortality)=="WPP"]="FMA"
	write.table(TradeOffs_FishingMortality,paste(PATH_highLevelSummaryPlots,"TradeOffs_FishingMortality.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
	#(3) Create Biomass-ratio Tradeoff Tables
	TradeOffs_Biomass = subset(BenchmarkSummary,select=c(species,WPP,B_over_Bspr20_median,B_over_Bspr30_median,B_over_Bspr40_median,
		B_over_Bmsy_median,B_over_Bmey_median))	
	TradeOffs_Biomass$B_over_Bspr20_median = round(TradeOffs_Biomass$B_over_Bspr20_median,2)
	TradeOffs_Biomass$B_over_Bspr30_median = round(TradeOffs_Biomass$B_over_Bspr30_median,2)	
	TradeOffs_Biomass$B_over_Bspr40_median = round(TradeOffs_Biomass$B_over_Bspr40_median,2)	
	TradeOffs_Biomass$B_over_Bmsy_median = round(TradeOffs_Biomass$B_over_Bmsy_median,2)	
	TradeOffs_Biomass$B_over_Bmey_median = round(TradeOffs_Biomass$B_over_Bmey_median,2)	
	names(TradeOffs_Biomass)[names(TradeOffs_Biomass)=="B_over_Bspr20_median"]="Bcurrent/B@SPR20%"
	names(TradeOffs_Biomass)[names(TradeOffs_Biomass)=="B_over_Bspr30_median"]="Bcurrent/B@SPR30%"
	names(TradeOffs_Biomass)[names(TradeOffs_Biomass)=="B_over_Bspr40_median"]="Bcurrent/B@SPR40%"
	names(TradeOffs_Biomass)[names(TradeOffs_Biomass)=="B_over_Bmsy_median"]="Bcurrent/B@MSY"
	names(TradeOffs_Biomass)[names(TradeOffs_Biomass)=="B_over_Bmey_median"]="Bcurrent/B@MEY"
	names(TradeOffs_Biomass)[names(TradeOffs_Biomass)=="species"]="Species"
	names(TradeOffs_Biomass)[names(TradeOffs_Biomass)=="WPP"]="FMA"
	write.table(TradeOffs_Biomass,paste(PATH_highLevelSummaryPlots,"TradeOffs_Biomass.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
	#(4) Create yield tradeoff tables
	TradeOffs_Yield = subset(BenchmarkSummary,select=c(species,WPP,Yield_current_geoMean,YieldSpr20_median,YieldSpr30_median,YieldSpr40_median,YieldMSY_median,YieldMEY_median))	
	TradeOffs_Yield$Yield_current_geoMean = prettyNum(round(TradeOffs_Yield$Yield_current_geoMean),big.mark=",",scientific=FALSE)
	TradeOffs_Yield$YieldSpr20_median = prettyNum(round(TradeOffs_Yield$YieldSpr20_median),big.mark=",",scientific=FALSE)
	TradeOffs_Yield$YieldSpr30_median = prettyNum(round(TradeOffs_Yield$YieldSpr30_median),big.mark=",",scientific=FALSE)	
	TradeOffs_Yield$YieldSpr40_median = prettyNum(round(TradeOffs_Yield$YieldSpr40_median),big.mark=",",scientific=FALSE)	
	TradeOffs_Yield$YieldMSY_median = prettyNum(round(TradeOffs_Yield$YieldMSY_median),big.mark=",",scientific=FALSE)	
	TradeOffs_Yield$YieldMEY_median = prettyNum(round(TradeOffs_Yield$YieldMEY_median),big.mark=",",scientific=FALSE)	
	names(TradeOffs_Yield)[names(TradeOffs_Yield)=="Yield_current_geoMean"]="Yield@Fcurrent"
	names(TradeOffs_Yield)[names(TradeOffs_Yield)=="YieldSpr20_median"]="Yield@SPR20%"
	names(TradeOffs_Yield)[names(TradeOffs_Yield)=="YieldSpr30_median"]="Yield@SPR30%"
	names(TradeOffs_Yield)[names(TradeOffs_Yield)=="YieldSpr40_median"]="Yield@SPR40%"
	names(TradeOffs_Yield)[names(TradeOffs_Yield)=="YieldMSY_median"]="Yield@MSY"
	names(TradeOffs_Yield)[names(TradeOffs_Yield)=="YieldMEY_median"]="Yield@MEY"
	names(TradeOffs_Yield)[names(TradeOffs_Yield)=="species"]="Species"
	names(TradeOffs_Yield)[names(TradeOffs_Yield)=="WPP"]="FMA"
	write.table(TradeOffs_Yield,paste(PATH_highLevelSummaryPlots,"TradeOffs_Yield.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
	#(5) Create biomass tradeoff tables
	TradeOffs_Biomass = subset(BenchmarkSummary,select=c(species,WPP,B_current_geoMean,Bspr20_median,Bspr30_median,Bspr40_median,Bmsy_median,Bmey_median))	
	TradeOffs_Biomass$B_current_geoMean = prettyNum(round(TradeOffs_Biomass$B_current_geoMean),big.mark=",",scientific=FALSE)
	TradeOffs_Biomass$Bspr20_median = prettyNum(round(TradeOffs_Biomass$Bspr20_median),big.mark=",",scientific=FALSE)
	TradeOffs_Biomass$Bspr30_median = prettyNum(round(TradeOffs_Biomass$Bspr30_median),big.mark=",",scientific=FALSE)	
	TradeOffs_Biomass$Bspr40_median = prettyNum(round(TradeOffs_Biomass$Bspr40_median),big.mark=",",scientific=FALSE)	
	TradeOffs_Biomass$Bmsy_median = prettyNum(round(TradeOffs_Biomass$Bmsy_median),big.mark=",",scientific=FALSE)	
	TradeOffs_Biomass$Bmey_median = prettyNum(round(TradeOffs_Biomass$Bmey_median),big.mark=",",scientific=FALSE)	
	names(TradeOffs_Biomass)[names(TradeOffs_Biomass)=="B_current_geoMean"]="Biomass@Fcurrent"
	names(TradeOffs_Biomass)[names(TradeOffs_Biomass)=="Bspr20_median"]="Biomass@SPR20%"
	names(TradeOffs_Biomass)[names(TradeOffs_Biomass)=="Bspr30_median"]="Biomass@SPR30%"
	names(TradeOffs_Biomass)[names(TradeOffs_Biomass)=="Bspr40_median"]="Biomass@SPR40%"
	names(TradeOffs_Biomass)[names(TradeOffs_Biomass)=="Bmsy_median"]="Biomass@MSY"
	names(TradeOffs_Biomass)[names(TradeOffs_Biomass)=="Bmey_median"]="Biomass@MEY"
	names(TradeOffs_Biomass)[names(TradeOffs_Biomass)=="species"]="Species"
	names(TradeOffs_Biomass)[names(TradeOffs_Biomass)=="WPP"]="FMA"
	write.table(TradeOffs_Biomass,paste(PATH_highLevelSummaryPlots,"TradeOffs_Biomass.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
	#(6) Create fishing mortality rate tradeoff tables
	TradeOffs_F = subset(BenchmarkSummary,select=c(species,WPP,F_current_median,Fspr20_median,Fspr30_median,Fspr40_median,Fmsy_median,Fmey_median))	
	TradeOffs_F$F_current_median = round(TradeOffs_F$F_current_median,3)
	TradeOffs_F$Fspr20_median = round(TradeOffs_F$Fspr20_median,3)
	TradeOffs_F$Fspr30_median = round(TradeOffs_F$Fspr30_median,3)
	TradeOffs_F$Fspr40_median = round(TradeOffs_F$Fspr40_median,3)
	TradeOffs_F$Fmsy_median = round(TradeOffs_F$Fmsy_median,3)
	TradeOffs_F$Fmey_median = round(TradeOffs_F$Fmey_median,3)
	names(TradeOffs_F)[names(TradeOffs_F)=="F_current_median"]="Fcurrent"
	names(TradeOffs_F)[names(TradeOffs_F)=="Fspr20_median"]="F@SPR20%"
	names(TradeOffs_F)[names(TradeOffs_F)=="Fspr30_median"]="F@SPR30%"
	names(TradeOffs_F)[names(TradeOffs_F)=="Fspr40_median"]="F@SPR40%"
	names(TradeOffs_F)[names(TradeOffs_F)=="Fmsy_median"]="F@MSY"
	names(TradeOffs_F)[names(TradeOffs_F)=="Fmey_median"]="F@MEY"
	names(TradeOffs_F)[names(TradeOffs_F)=="species"]="Species"
	names(TradeOffs_F)[names(TradeOffs_F)=="WPP"]="FMA"
	write.table(TradeOffs_F,paste(PATH_highLevelSummaryPlots,"TradeOffs_F.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
	#(7) Create profit tradeoff tables
	BenchmarkSummary$Profit_current_median = BenchmarkSummary$Revenue_Fcurrent  - BenchmarkSummary$Cost_Fcurrent
	BenchmarkSummary$Profit_Fspr20_median = BenchmarkSummary$Revenue_Fspr20 - BenchmarkSummary$Cost_Fspr20
	BenchmarkSummary$Profit_Fspr30_median = BenchmarkSummary$Revenue_Fspr30 - BenchmarkSummary$Cost_Fspr30
	BenchmarkSummary$Profit_Fspr40_median = BenchmarkSummary$Revenue_Fspr40 - BenchmarkSummary$Cost_Fspr40
	BenchmarkSummary$Profit_MSY_median = BenchmarkSummary$Revenue_MSY - BenchmarkSummary$Cost_MSY
	BenchmarkSummary$Profit_MEY_median = BenchmarkSummary$Revenue_MEY - BenchmarkSummary$Cost_MEY
	TradeOffs_Profit = subset(BenchmarkSummary,select=c(species,WPP,Profit_current_median,Profit_Fspr20_median,
		Profit_Fspr30_median,Profit_Fspr40_median,Profit_MSY_median,Profit_MEY_median))
	TradeOffs_Profit$Profit_current_median = prettyNum(round(TradeOffs_Profit$Profit_current_median),big.mark=",",scientific=FALSE)
	TradeOffs_Profit$Profit_Fspr20_median = prettyNum(round(TradeOffs_Profit$Profit_Fspr20_median),big.mark=",",scientific=FALSE)
	TradeOffs_Profit$Profit_Fspr30_median = prettyNum(round(TradeOffs_Profit$Profit_Fspr30_median),big.mark=",",scientific=FALSE)	
	TradeOffs_Profit$Profit_Fspr40_median = prettyNum(round(TradeOffs_Profit$Profit_Fspr40_median),big.mark=",",scientific=FALSE)	
	TradeOffs_Profit$Profit_MSY_median = prettyNum(round(TradeOffs_Profit$Profit_MSY_median),big.mark=",",scientific=FALSE)	
	TradeOffs_Profit$Profit_MEY_median = prettyNum(round(TradeOffs_Profit$Profit_MEY_median),big.mark=",",scientific=FALSE)	
	names(TradeOffs_Profit)[names(TradeOffs_Profit)=="Profit_current_median"]="Profit@Fcurrent"
	names(TradeOffs_Profit)[names(TradeOffs_Profit)=="Profit_Fspr20_median"]="Profit@SPR20%"
	names(TradeOffs_Profit)[names(TradeOffs_Profit)=="Profit_Fspr30_median"]="Profit@SPR30%"
	names(TradeOffs_Profit)[names(TradeOffs_Profit)=="Profit_Fspr40_median"]="Profit@SPR40%"
	names(TradeOffs_Profit)[names(TradeOffs_Profit)=="Profit_MSY_median"]="Profit@MSY"
	names(TradeOffs_Profit)[names(TradeOffs_Profit)=="Profit_MEY_median"]="Profit@MEY"
	names(TradeOffs_Profit)[names(TradeOffs_Profit)=="species"]="Species"
	names(TradeOffs_Profit)[names(TradeOffs_Profit)=="WPP"]="FMA"
	write.table(TradeOffs_Profit,paste(PATH_highLevelSummaryPlots,"TradeOffs_Profit.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
	#(8) Calculate Tradeoff Rebuilding Times
	ExtrapolatedLandings = read.table(paste(PATH_output,"ExtrapolatedLandings.csv",sep="/"),header=TRUE,sep=",")
	SpeciesFamily = subset(ExtrapolatedLandings,select=c(species,FamilyCommonName))
	SpeciesFamily = unique(SpeciesFamily)
	CurrentStatusValues = read.table(paste(PATH_output,"Benchmarks_withFlags.csv",sep="/"),header=TRUE,sep=",",as.is=TRUE)
	CurrentStatusValues = subset(CurrentStatusValues,CurrentStatusValues$Fratio40_FLAG==0 & CurrentStatusValues$Bratio20_FLAG==0)
	CurrentStatusValues = merge(CurrentStatusValues,SpeciesFamily,all.x=TRUE,by=c("species"))
	WPPs = sort(unique(CurrentStatusValues$WPP))
	RebuildTime_Fcurrent = aggregate.data.frame(list(CurrentStatusValues$B_over_Bspr20,CurrentStatusValues$RebuildTime_to_spr20_New.CurrentF),by=list(CurrentStatusValues$species,CurrentStatusValues$WPP),FUN=median,na.rm=TRUE)
	names(RebuildTime_Fcurrent) = c("species","WPP","B_over_Bspr20","RebuildTimeToSPR20_currentF")
	RebuildTime_Fcurrent$RebuildTimeToSPR20_currentF[RebuildTime_Fcurrent$B_over_Bspr20>=1]=NA
	RebuildTime_Fcurrent$B_over_Bspr20=NULL
	RebuildTime_Fspr20 = aggregate.data.frame(list(CurrentStatusValues$B_over_Bspr20,CurrentStatusValues$RebuildTime_to_spr20_New.F_at_SPR20),by=list(CurrentStatusValues$species,CurrentStatusValues$WPP),FUN=median,na.rm=TRUE)
	names(RebuildTime_Fspr20) = c("species","WPP","B_over_Bspr20","RebuildTimeToSPR20_Fspr20")
	RebuildTime_Fspr20$RebuildTimeToSPR20_Fspr20[RebuildTime_Fspr20$B_over_Bspr20>=1]=NA
	RebuildTime_Fspr20$B_over_Bspr20=NULL
	RebuildTime_Fspr30 = aggregate.data.frame(list(CurrentStatusValues$B_over_Bspr20,CurrentStatusValues$RebuildTime_to_spr20_New.F_at_SPR30),by=list(CurrentStatusValues$species,CurrentStatusValues$WPP),FUN=median,na.rm=TRUE)
	names(RebuildTime_Fspr30) = c("species","WPP","B_over_Bspr20","RebuildTimeToSPR20_Fspr30")
	RebuildTime_Fspr30$RebuildTimeToSPR20_Fspr30[RebuildTime_Fspr30$B_over_Bspr20>=1]=NA
	RebuildTime_Fspr30$B_over_Bspr20=NULL	
	RebuildTime_Fspr40 = aggregate.data.frame(list(CurrentStatusValues$B_over_Bspr20,CurrentStatusValues$RebuildTime_to_spr20_New.F_at_SPR40),by=list(CurrentStatusValues$species,CurrentStatusValues$WPP),FUN=median,na.rm=TRUE)
	names(RebuildTime_Fspr40) = c("species","WPP","B_over_Bspr20","RebuildTimeToSPR20_Fspr40")
	RebuildTime_Fspr40$RebuildTimeToSPR20_Fspr40[RebuildTime_Fspr40$B_over_Bspr20>=1]=NA
	RebuildTime_Fspr40$B_over_Bspr20=NULL	
	RebuildTime_Fmsy = aggregate.data.frame(list(CurrentStatusValues$B_over_Bspr20,CurrentStatusValues$RebuildTime_to_spr20_New.MSY),by=list(CurrentStatusValues$species,CurrentStatusValues$WPP),FUN=median,na.rm=TRUE)
	names(RebuildTime_Fmsy) = c("species","WPP","B_over_Bspr20","RebuildTimeToSPR20_Fmsy")
	RebuildTime_Fmsy$RebuildTimeToSPR20_Fmsy[RebuildTime_Fmsy$B_over_Bspr20>=1]=NA
	RebuildTime_Fmsy$B_over_Bspr20=NULL	
	RebuildTime_Fmey = aggregate.data.frame(list(CurrentStatusValues$B_over_Bspr20,CurrentStatusValues$RebuildTime_to_spr20_New.MEY),by=list(CurrentStatusValues$species,CurrentStatusValues$WPP),FUN=median,na.rm=TRUE)
	names(RebuildTime_Fmey) = c("species","WPP","B_over_Bspr20","RebuildTimeToSPR20_Fmey")
	RebuildTime_Fmey$RebuildTimeToSPR20_Fmey[RebuildTime_Fmey$B_over_Bspr20>=1]=NA
	RebuildTime_Fmey$B_over_Bspr20=NULL	
	TradeOffs_RebuildTime = merge(RebuildTime_Fcurrent,RebuildTime_Fspr20,all.x=TRUE,by=c("species","WPP"))
	TradeOffs_RebuildTime = merge(TradeOffs_RebuildTime,RebuildTime_Fspr30,all.x=TRUE,by=c("species","WPP"))
	TradeOffs_RebuildTime = merge(TradeOffs_RebuildTime,RebuildTime_Fspr40,all.x=TRUE,by=c("species","WPP"))
	TradeOffs_RebuildTime = merge(TradeOffs_RebuildTime,RebuildTime_Fmsy,all.x=TRUE,by=c("species","WPP"))
	TradeOffs_RebuildTime = merge(TradeOffs_RebuildTime,RebuildTime_Fmey,all.x=TRUE,by=c("species","WPP"))
	names(TradeOffs_RebuildTime)[names(TradeOffs_RebuildTime)=="RebuildTimeToSPR20_currentF"]="RebuildTime@Fcurrent"
	names(TradeOffs_RebuildTime)[names(TradeOffs_RebuildTime)=="RebuildTimeToSPR20_Fspr20"]="RebuildTime@SPR20%"
	names(TradeOffs_RebuildTime)[names(TradeOffs_RebuildTime)=="RebuildTimeToSPR20_Fspr30"]="RebuildTime@SPR30%"
	names(TradeOffs_RebuildTime)[names(TradeOffs_RebuildTime)=="RebuildTimeToSPR20_Fspr40"]="RebuildTime@SPR40%"
	names(TradeOffs_RebuildTime)[names(TradeOffs_RebuildTime)=="RebuildTimeToSPR20_Fmsy"]="RebuildTime@MSY"
	names(TradeOffs_RebuildTime)[names(TradeOffs_RebuildTime)=="RebuildTimeToSPR20_Fmey"]="RebuildTime@MEY"
	names(TradeOffs_RebuildTime)[names(TradeOffs_RebuildTime)=="species"]="Species"
	names(TradeOffs_RebuildTime)[names(TradeOffs_RebuildTime)=="WPP"]="FMA"
	write.table(TradeOffs_RebuildTime,paste(PATH_highLevelSummaryPlots,"TradeOffs_RebuildTime.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
	#Percent F Reduction Tradeoff Table
	BenchmarkSummary$PercentReductionF_toSPR20 = round((((BenchmarkSummary$F_current_median - BenchmarkSummary$Fspr20_median)/BenchmarkSummary$F_current_median)*100)*-1)
	BenchmarkSummary$PercentReductionF_toSPR20[BenchmarkSummary$PercentReductionF_toSPR20>0]=NA
	BenchmarkSummary$PercentReductionF_toSPR20[BenchmarkSummary$B_over_Bspr20_median>=1]=NA	
	BenchmarkSummary$PercentReductionF_toSPR30 = round((((BenchmarkSummary$F_current_median - BenchmarkSummary$Fspr30_median)/BenchmarkSummary$F_current_median)*100)*-1)
	BenchmarkSummary$PercentReductionF_toSPR30[BenchmarkSummary$PercentReductionF_toSPR30>0]=NA
	BenchmarkSummary$PercentReductionF_toSPR30[BenchmarkSummary$B_over_Bspr30_median>=1]=NA	
	BenchmarkSummary$PercentReductionF_toSPR40 = round((((BenchmarkSummary$F_current_median - BenchmarkSummary$Fspr40_median)/BenchmarkSummary$F_current_median)*100)*-1)
	BenchmarkSummary$PercentReductionF_toSPR40[BenchmarkSummary$PercentReductionF_toSPR40>0]=NA
	BenchmarkSummary$PercentReductionF_toSPR40[BenchmarkSummary$B_over_Bspr40_median>=1]=NA	
	BenchmarkSummary$PercentReductionF_toMSY = round((((BenchmarkSummary$F_current_median - BenchmarkSummary$Fmsy_median)/BenchmarkSummary$F_current_median)*100)*-1)
	BenchmarkSummary$PercentReductionF_toMSY[BenchmarkSummary$PercentReductionF_toMSY>0]=NA
	BenchmarkSummary$PercentReductionF_toMSY[BenchmarkSummary$B_over_Bmsy_median>=1]=NA	
	BenchmarkSummary$PercentReductionF_toMEY = round((((BenchmarkSummary$F_current_median - BenchmarkSummary$Fmey_median)/BenchmarkSummary$F_current_median)*100)*-1)
	BenchmarkSummary$PercentReductionF_toMEY[BenchmarkSummary$PercentReductionF_toMEY>0]=NA
	BenchmarkSummary$PercentReductionF_toMEY[BenchmarkSummary$B_over_Bmey_median>=1]=NA	
	TradeOffs_PercentReductionF = subset(BenchmarkSummary,select=c(species,WPP,PercentReductionF_toSPR20,PercentReductionF_toSPR30,PercentReductionF_toSPR40,PercentReductionF_toMSY,PercentReductionF_toMEY))
	names(TradeOffs_PercentReductionF)[names(TradeOffs_PercentReductionF)=="PercentReductionF_toSPR20"]="PercReduceF@SPR20%"
	names(TradeOffs_PercentReductionF)[names(TradeOffs_PercentReductionF)=="RebuildTimeToSPR20_Fspr30"]="PercReduceF@SPR30%"
	names(TradeOffs_PercentReductionF)[names(TradeOffs_PercentReductionF)=="RebuildTimeToSPR20_Fspr40"]="PercReduceF@SPR40%"
	names(TradeOffs_PercentReductionF)[names(TradeOffs_PercentReductionF)=="RebuildTimeToSPR20_Fmsy"]="PercReduceF@MSY"
	names(TradeOffs_PercentReductionF)[names(TradeOffs_PercentReductionF)=="RebuildTimeToSPR20_Fmey"]="PercReduceF@MEY"
	names(TradeOffs_PercentReductionF)[names(TradeOffs_PercentReductionF)=="species"]="Species"
	names(TradeOffs_PercentReductionF)[names(TradeOffs_PercentReductionF)=="WPP"]="FMA"
	write.table(TradeOffs_PercentReductionF,paste(PATH_highLevelSummaryPlots,"TradeOffs_PercentReductionF.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)

}


############## Maps with Percent Overfished, Overfishing, and Uncertain: Three Pie Partitions ####################
#Read shape file and setup locations
if(createMapsWithPieCharts==TRUE)
{
	WPPshapeFile = readOGR("F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/SpatialDistributionBiomass/GIS_Data",layer="WPP_Indo_PKP18_2014_pl_geo")
	LocationsWPP = data.frame(ShapeFileName=as.character(WPPshapeFile$Descriptio),stringsAsFactors=FALSE)
	rowsToKeep=c()
	for(i in 1:length(WPP))
	{
	   rowsToKeep = c(rowsToKeep,grep(pattern=as.character(WPP[i]),x=LocationsWPP$ShapeFileName))
	}
	LocationsWPP = LocationsWPP[rowsToKeep,]
	LocationsWPP = data.frame(LocationsWPP_ShapeFile=LocationsWPP,WPP=WPP)
	Midpoints = data.frame(WPP=c(571,572,573,711,712,713,714,715,716,717,718),
	  Longitude=c(97.5,98.59572,114.72720,109.47087,111.40423,118.53350,126.44820,125.84402,122.70231,135.99417,135.45041),
	  Latitude=c(6,-1.4917099,-8.3130337,5.5395008,-3.3282201,-2.6460877,-3.8004656,0.1873852,2.6535561,1.4991783,-6.1616931))
	LocationsWPP = merge(LocationsWPP,Midpoints,all.x=TRUE,by=c("WPP"))
	LocationsWPP = subset(LocationsWPP,LocationsWPP$WPP!="717")
	#Read BenchmarkSummary file
	BenchmarkSummary = read.table(paste(PATH_output,"BenchmarkSummary.csv",sep="/"),header=TRUE,sep=",")
	BenchmarkSummary$TargetStatus = "Above or Below Target"
	BenchmarkSummary$TargetStatus[BenchmarkSummary$Fratio40_UpperCI<1]="Below Target"
	BenchmarkSummary$TargetStatus[BenchmarkSummary$Fratio40_LowerCI>1]="Above Target"
	BenchmarkSummary$LimitStatus = "Above or Below Limit"
	BenchmarkSummary$LimitStatus[BenchmarkSummary$Bratio20_UpperCI<1]="Below Limit"
	BenchmarkSummary$LimitStatus[BenchmarkSummary$Bratio20_LowerCI>1]="Above Limit"
	TargetStatus = as.data.frame(table(BenchmarkSummary$WPP,BenchmarkSummary$TargetStatus))
	names(TargetStatus) = c("WPP","TargetStatusLabel","NumTargetStatus")
	TargetStatusTotal = aggregate.data.frame(TargetStatus$NumTargetStatus,by=list(TargetStatus$WPP),FUN=sum)
	names(TargetStatusTotal) = c("WPP","TargetStatusTotal")
	TargetStatus = merge(TargetStatus,TargetStatusTotal,all.x=TRUE,by=c("WPP"))
	TargetStatus$TargetStatusPerc = TargetStatus$NumTargetStatus/TargetStatus$TargetStatusTotal
	TargetStatus = merge(TargetStatus,LocationsWPP,all.x=TRUE,by=c("WPP"))
	TargetStatus$LocationsWPP_ShapeFile=NULL
	TargetStatus = subset(TargetStatus,TargetStatus$WPP!="EEZ")
	LimitStatus = as.data.frame(table(BenchmarkSummary$WPP,BenchmarkSummary$LimitStatus))
	names(LimitStatus) = c("WPP","LimitStatusLabel","NumLimitStatus")
	LimitStatusTotal = aggregate.data.frame(LimitStatus$NumLimitStatus,by=list(LimitStatus$WPP),FUN=sum)
	names(LimitStatusTotal) = c("WPP","LimitStatusTotal")
	LimitStatus = merge(LimitStatus,LimitStatusTotal,all.x=TRUE,by=c("WPP"))
	LimitStatus$LimitStatusPerc = LimitStatus$NumLimitStatus/LimitStatus$LimitStatusTotal
	LimitStatus = merge(LimitStatus,LocationsWPP,all.x=TRUE,by=c("WPP"))
	LimitStatus$LocationsWPP_ShapeFile=NULL
	LimitStatus = subset(LimitStatus,LimitStatus$WPP!="EEZ")
	#Pie Chart Sample Size Text Locations for Each Pie
	TextSampleSizeLocations = data.frame(x=c(100,100.46866,117.3,107.47709,113.94177,118.77517,126.68987,127.05237,122.58147,138.7923,135.75250),
	  y=c(7,-3.0658615,-8.3130337,4.2277078,-3.3282201,-0.8620492,-5.4795607,1.6467065,4.6474816,1.012066,-7.9982033),
	  WPP=c(571,572,573,711,712,713,714,715,716,717,718))
	TextSampleSizeLocations = subset(TextSampleSizeLocations,TextSampleSizeLocations$WPP!="717")  
	#Plot Overfished Status: SPR20%
	xyz = make.xyz(LimitStatus$Longitude,LimitStatus$Latitude,LimitStatus$LimitStatusPerc,LimitStatus$LimitStatusLabel)
	LatLongXYZorder = data.frame(Longitude=xyz$x,Latitude=xyz$y)
	LimitStatus = merge(LatLongXYZorder,LimitStatus,all.x=TRUE,by=c("Longitude","Latitude"))
	SampleSize = subset(LimitStatus,select=c("Longitude","Latitude","LimitStatusTotal"))
	SampleSize = unique(SampleSize)
	png(paste(PATH_highLevelSummaryPlots,"OverfishedMapSPR20.png",sep="/"),units="px",width=5200,height=3200,res=600)
	par(mai = c(0.5, 0.5, 0.35, 0.2), omi = c(0.25, 0.5, 0, 0),mgp = c(2.5, 0.5, 0))
	prettymap(plot(WPPshapeFile,col=c("darkgreen",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),xlab="Longitude",ylab="Latitude"), oma=c(4,4,1,0),drawscale=TRUE,drawbox=FALSE, scale.plotunit="latlon", drawarrow=TRUE)
	axis(side=1,outer=TRUE,line=-2.5)
	axis(side=2,outer=TRUE,line=-2)
	mtext("Longitude",side=1,line=2,cex=1.2)
	mtext("Latitude",side=2,line=2.25,cex=1.2)
	title("Population Status: Limit Reference Point (Biomass at SPR 20%)")
	mtext("Indonesia Snapper-Grouper Fishery",font=3,line=-0.5)
	for(i in 1:(dim(SampleSize)[1]))
	{
	  #col=c("grey70","grey30","grey3")
	  add.pie(xyz$z[i,],xyz$x[i],xyz$y[i],radius=1.2,col=c("forestgreen","gold","tomato"),labels=NA) #,density=c(40,40,40),angle=c(45,90,135))  #,labels=c(NA,paste("n=",Status_WPP$N[i],sep="")),label.dist=1.3)
	}
	legend.pie(100,-10,labels=c("Above Limit","Above or Below Limit","Below Limit"),radius=1,bty="n",col=c("forestgreen","gold","tomato"),cex=0.9,label.dist=1.6)
	for(i in 1:(dim(SampleSize)[1]))
	{	
		text(TextSampleSizeLocations$x[i],TextSampleSizeLocations$y[i],paste("n=",SampleSize$LimitStatusTotal[i],sep=""))
	}
	dev.off()
	#Plot Overfishing Status at SPR 40%
	xyz = make.xyz(TargetStatus$Longitude,TargetStatus$Latitude,TargetStatus$TargetStatusPerc,TargetStatus$TargetStatusLabel)
	LatLongXYZorder = data.frame(Longitude=xyz$x,Latitude=xyz$y)
	TargetStatus = merge(LatLongXYZorder,TargetStatus,all.x=TRUE,by=c("Longitude","Latitude"))
	SampleSize = subset(TargetStatus,select=c("Longitude","Latitude","TargetStatusTotal"))
	SampleSize = unique(SampleSize) 
	png(paste(PATH_highLevelSummaryPlots,"OverfishingMapSPR40.png",sep="/"),units="px",width=5200,height=3200,res=600)
	par(mai = c(0.5, 0.5, 0.35, 0.2), omi = c(0.25, 0.5, 0, 0),mgp = c(2.5, 0.5, 0))
	prettymap(plot(WPPshapeFile,col=c("darkgreen",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),xlab="Longitude",ylab="Latitude"), oma=c(4,4,1,0),drawscale=TRUE,drawbox=FALSE, scale.plotunit="latlon", drawarrow=TRUE)
	axis(side=1,outer=TRUE,line=-2.5)
	axis(side=2,outer=TRUE,line=-2)
	mtext("Longitude",side=1,line=2,cex=1.2)
	mtext("Latitude",side=2,line=2.25,cex=1.2)
	title("Population Status: Target Reference Point (Fishing Mortality at SPR 40%)")
	mtext("Indonesia Snapper-Grouper Fishery",font=3,line=-0.5)
	for(i in 1:(dim(SampleSize)[1]))
	{
	  add.pie(xyz$z[i,],xyz$x[i],xyz$y[i],radius=1.2,col=c("gold","tomato","forestgreen"),labels=NA)  #,labels=c(NA,paste("n=",Status_WPP$N[i],sep="")),label.dist=1.3)
	}
	legend.pie(100,-10,labels=c("Above Target","Above or Below Target","Below Target"),radius=1,bty="n",col=c("tomato","gold","forestgreen"),cex=0.9,label.dist=1.6)
	for(i in 1:(dim(SampleSize)[1]))
	{
		text(TextSampleSizeLocations$x[i],TextSampleSizeLocations$y[i],paste("n=",SampleSize$TargetStatusTotal[i],sep=""))	 
	}
	dev.off()
}


########################## Nutritional Calculations for TradeOff Analysis ####################################
if(calculateNutritionalTradeOffs==TRUE)
{
	TradeOffs_Yield = read.table(paste(PATH_highLevelSummaryPlots,"TradeOffs_Yield.csv",sep="/"),sep=",",header=TRUE)
	for(i in 1:dim(TradeOffs_Yield)[1])
	{
		TradeOffs_Yield$Yield.Fcurrent[i] = gsub(",","",TradeOffs_Yield$Yield.Fcurrent[i])
		TradeOffs_Yield$Yield.SPR20.[i] = gsub(",","",TradeOffs_Yield$Yield.SPR20.[i])
		TradeOffs_Yield$Yield.SPR30.[i] = gsub(",","",TradeOffs_Yield$Yield.SPR30.[i])		
		TradeOffs_Yield$Yield.SPR40.[i] = gsub(",","",TradeOffs_Yield$Yield.SPR40.[i])		
		TradeOffs_Yield$Yield.MSY[i] = gsub(",","",TradeOffs_Yield$Yield.MSY[i])	
		TradeOffs_Yield$Yield.MEY[i] = gsub(",","",TradeOffs_Yield$Yield.MEY[i])	
	}
	TradeOffs_Yield$Yield.Fcurrent = as.numeric(TradeOffs_Yield$Yield.Fcurrent)
	TradeOffs_Yield$Yield.SPR20. = as.numeric(TradeOffs_Yield$Yield.SPR20.)
	TradeOffs_Yield$Yield.SPR30. = as.numeric(TradeOffs_Yield$Yield.SPR30.)
	TradeOffs_Yield$Yield.SPR40. = as.numeric(TradeOffs_Yield$Yield.SPR40.)
	TradeOffs_Yield$Yield.MSY = as.numeric(TradeOffs_Yield$Yield.MSY)	
	TradeOffs_Yield$Yield.MEY = as.numeric(TradeOffs_Yield$Yield.MEY)	
	names(TradeOffs_Yield)[names(TradeOffs_Yield)=="Species"]="species"
	names(TradeOffs_Yield)[names(TradeOffs_Yield)=="FMA"]="WPP"	
	#Units of nutrition are mg per 100 grams of fish 
	SpeciesNutritionPredictions = read.table(paste(PATH,"NutriCast","Species_Nutrient_Predictions.csv",sep="/"),sep=",",header=TRUE)
	#Fill holes for species not in the nutritional database with very similar species that have the same genus / phylogeny 
	#for "Atrobucca brevis" use "Atrobucca nibe"
	Substitute = subset(SpeciesNutritionPredictions,SpeciesNutritionPredictions$species=="Atrobucca nibe")
	Substitute$species = "Atrobucca brevis"
	SpeciesNutritionPredictions = rbind(SpeciesNutritionPredictions,Substitute)
	#for "Erythrocles schlegelii" use "Erythrocles monodi"
	Substitute = subset(SpeciesNutritionPredictions,SpeciesNutritionPredictions$species=="Erythrocles monodi")
	Substitute$species = "Erythrocles schlegelii"
	SpeciesNutritionPredictions = rbind(SpeciesNutritionPredictions,Substitute)
	#for "Glaucosoma buergeri" use average of "Glaucosoma hebraicum","Glaucosoma magnificum",and "Glaucosoma scapulare"
	Substitute = subset(SpeciesNutritionPredictions,SpeciesNutritionPredictions$species=="Glaucosoma hebraicum" | SpeciesNutritionPredictions$species=="Glaucosoma magnificum" | SpeciesNutritionPredictions$species=="Glaucosoma scapulare")
	Substitute$species = "Glaucosoma buergeri"
	Substitute = suppressWarnings(aggregate.data.frame(Substitute,by=list(Substitute$species),FUN=mean))
	Substitute$Group.1=NULL
	Substitute$species = "Glaucosoma buergeri"
	SpeciesNutritionPredictions = rbind(SpeciesNutritionPredictions,Substitute)
	#for "Paracaesio gonzalesi" use average of "Paracaesio kusakarii","Paracaesio sordida","Paracaesio stonei","Paracaesio xanthura"
	Substitute = subset(SpeciesNutritionPredictions,SpeciesNutritionPredictions$species=="Paracaesio kusakarii" | SpeciesNutritionPredictions$species=="Paracaesio sordida" | SpeciesNutritionPredictions$species=="Paracaesio stonei" | SpeciesNutritionPredictions$species=="Paracaesio xanthura")
	Substitute$species = "Paracaesio gonzalesi"
	Substitute = suppressWarnings(aggregate.data.frame(Substitute,by=list(Substitute$species),FUN=mean))
	Substitute$Group.1=NULL
	Substitute$species = "Paracaesio gonzalesi"
	SpeciesNutritionPredictions = rbind(SpeciesNutritionPredictions,Substitute)
	#for "Pinjalo lewisi" use "Pinjalo pinjalo"
	Substitute = subset(SpeciesNutritionPredictions,SpeciesNutritionPredictions$species=="Pinjalo pinjalo")
	Substitute$species = "Pinjalo lewisi"
	SpeciesNutritionPredictions = rbind(SpeciesNutritionPredictions,Substitute)
	SpeciesNutritionPredictions = subset(SpeciesNutritionPredictions,select=c(species,Selenium_mu,Zinc_mu,Protein_mu,Omega_3_mu,Calcium_mu,Iron_mu,Vitamin_A_mu))
	#Daily adult nutritional needs:
	#Selenium is 55 micrograms (1 microgram = 0.001 milligrams)
	#Zinc is 11 mg men and 8 mg women
	#Protein is about 60 grams per day (0.8 grams per kg per human body weight)
	#Omega 3 is about 450 to 500 mg per day
	#Calcium is about 1000 mg per day
	#Iron is 8 mg men 18 mg women (pregnant is 27 mg)
	#Vitamin A is 900 micrograms for men and 700 micrograms for women (1 microgram = 0.001 milligrams)
	DailyRequirements = data.frame(Nutrient=c("Selenium","Zinc","Protein","Omega3","Calcium","Iron","VitaminA"),DailyValue_mg=c(0.055,9.5,60000,475,1000,13,0.8))
	TradeOffs_Yield_Labels = subset(TradeOffs_Yield,select=c(species,WPP))
	TradeOffs_Yield_NoLabels = subset(TradeOffs_Yield,select=-c(species,WPP))
	NutritionSpeciesWPP = merge(TradeOffs_Yield_Labels,SpeciesNutritionPredictions,all.x=TRUE,by=c("species"))
	names(NutritionSpeciesWPP) = gsub("_mu","",names(NutritionSpeciesWPP))
	#Current Yield
	Nutrition_currentYield = subset(TradeOffs_Yield,select=c(species,WPP,Yield.Fcurrent))
	Nutrition_currentYield = merge(NutritionSpeciesWPP,Nutrition_currentYield,all.x=TRUE,by=c("species","WPP"))
	Nutrition_currentYield$Selenium = Nutrition_currentYield$Selenium*(Nutrition_currentYield$Yield.Fcurrent*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_currentYield$Zinc = Nutrition_currentYield$Zinc*(Nutrition_currentYield$Yield.Fcurrent*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_currentYield$Protein = Nutrition_currentYield$Protein*(Nutrition_currentYield$Yield.Fcurrent*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_currentYield$Omega_3 = Nutrition_currentYield$Omega_3*(Nutrition_currentYield$Yield.Fcurrent*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_currentYield$Calcium = Nutrition_currentYield$Calcium*(Nutrition_currentYield$Yield.Fcurrent*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_currentYield$Iron = Nutrition_currentYield$Iron*(Nutrition_currentYield$Yield.Fcurrent*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_currentYield$Vitamin_A = Nutrition_currentYield$Vitamin_A*(Nutrition_currentYield$Yield.Fcurrent*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_currentYield$Yield.Fcurrent=NULL
	Nutrition_currentYield = aggregate.data.frame(list(Nutrition_currentYield$Selenium,Nutrition_currentYield$Zinc,Nutrition_currentYield$Protein,Nutrition_currentYield$Omega_3,Nutrition_currentYield$Calcium,Nutrition_currentYield$Iron,Nutrition_currentYield$Vitamin_A),by=list(Nutrition_currentYield$WPP),FUN=sum)
	names(Nutrition_currentYield) = c("WPP","Selenium","Zinc","Protein","Omega_3","Calcium","Iron","Vitamin_A")
	#Yield When Managed at SPR20%
	Nutrition_YieldSpr20 = subset(TradeOffs_Yield,select=c(species,WPP,Yield.SPR20.))
	Nutrition_YieldSpr20 = merge(NutritionSpeciesWPP,Nutrition_YieldSpr20,all.x=TRUE,by=c("species","WPP"))
	Nutrition_YieldSpr20$Selenium = Nutrition_YieldSpr20$Selenium*(Nutrition_YieldSpr20$Yield.SPR20.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr20$Zinc = Nutrition_YieldSpr20$Zinc*(Nutrition_YieldSpr20$Yield.SPR20.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr20$Protein = Nutrition_YieldSpr20$Protein*(Nutrition_YieldSpr20$Yield.SPR20.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr20$Omega_3 = Nutrition_YieldSpr20$Omega_3*(Nutrition_YieldSpr20$Yield.SPR20.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr20$Calcium = Nutrition_YieldSpr20$Calcium*(Nutrition_YieldSpr20$Yield.SPR20.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr20$Iron = Nutrition_YieldSpr20$Iron*(Nutrition_YieldSpr20$Yield.SPR20.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr20$Vitamin_A = Nutrition_YieldSpr20$Vitamin_A*(Nutrition_YieldSpr20$Yield.SPR20.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr20$Yield.SPR20.=NULL
	Nutrition_YieldSpr20 = aggregate.data.frame(list(Nutrition_YieldSpr20$Selenium,Nutrition_YieldSpr20$Zinc,Nutrition_YieldSpr20$Protein,Nutrition_YieldSpr20$Omega_3,Nutrition_YieldSpr20$Calcium,Nutrition_YieldSpr20$Iron,Nutrition_YieldSpr20$Vitamin_A),by=list(Nutrition_YieldSpr20$WPP),FUN=sum)
	names(Nutrition_YieldSpr20) = c("WPP","Selenium","Zinc","Protein","Omega_3","Calcium","Iron","Vitamin_A")
	#Yield When Managed at SPR30%
	Nutrition_YieldSpr30 = subset(TradeOffs_Yield,select=c(species,WPP,Yield.SPR30.))
	Nutrition_YieldSpr30 = merge(NutritionSpeciesWPP,Nutrition_YieldSpr30,all.x=TRUE,by=c("species","WPP"))
	Nutrition_YieldSpr30$Selenium = Nutrition_YieldSpr30$Selenium*(Nutrition_YieldSpr30$Yield.SPR30.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr30$Zinc = Nutrition_YieldSpr30$Zinc*(Nutrition_YieldSpr30$Yield.SPR30.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr30$Protein = Nutrition_YieldSpr30$Protein*(Nutrition_YieldSpr30$Yield.SPR30.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr30$Omega_3 = Nutrition_YieldSpr30$Omega_3*(Nutrition_YieldSpr30$Yield.SPR30.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr30$Calcium = Nutrition_YieldSpr30$Calcium*(Nutrition_YieldSpr30$Yield.SPR30.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr30$Iron = Nutrition_YieldSpr30$Iron*(Nutrition_YieldSpr30$Yield.SPR30.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr30$Vitamin_A = Nutrition_YieldSpr30$Vitamin_A*(Nutrition_YieldSpr30$Yield.SPR30.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr30$Yield.SPR30.=NULL
	Nutrition_YieldSpr30 = aggregate.data.frame(list(Nutrition_YieldSpr30$Selenium,Nutrition_YieldSpr30$Zinc,Nutrition_YieldSpr30$Protein,Nutrition_YieldSpr30$Omega_3,Nutrition_YieldSpr30$Calcium,Nutrition_YieldSpr30$Iron,Nutrition_YieldSpr30$Vitamin_A),by=list(Nutrition_YieldSpr30$WPP),FUN=sum)
	names(Nutrition_YieldSpr30) = c("WPP","Selenium","Zinc","Protein","Omega_3","Calcium","Iron","Vitamin_A")
	#Yield When Managed at SPR40%
	Nutrition_YieldSpr40 = subset(TradeOffs_Yield,select=c(species,WPP,Yield.SPR40.))
	Nutrition_YieldSpr40 = merge(NutritionSpeciesWPP,Nutrition_YieldSpr40,all.x=TRUE,by=c("species","WPP"))
	Nutrition_YieldSpr40$Selenium = Nutrition_YieldSpr40$Selenium*(Nutrition_YieldSpr40$Yield.SPR40.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr40$Zinc = Nutrition_YieldSpr40$Zinc*(Nutrition_YieldSpr40$Yield.SPR40.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr40$Protein = Nutrition_YieldSpr40$Protein*(Nutrition_YieldSpr40$Yield.SPR40.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr40$Omega_3 = Nutrition_YieldSpr40$Omega_3*(Nutrition_YieldSpr40$Yield.SPR40.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr40$Calcium = Nutrition_YieldSpr40$Calcium*(Nutrition_YieldSpr40$Yield.SPR40.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr40$Iron = Nutrition_YieldSpr40$Iron*(Nutrition_YieldSpr40$Yield.SPR40.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr40$Vitamin_A = Nutrition_YieldSpr40$Vitamin_A*(Nutrition_YieldSpr40$Yield.SPR40.*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldSpr40$Yield.SPR40.=NULL
	Nutrition_YieldSpr40 = aggregate.data.frame(list(Nutrition_YieldSpr40$Selenium,Nutrition_YieldSpr40$Zinc,Nutrition_YieldSpr40$Protein,Nutrition_YieldSpr40$Omega_3,Nutrition_YieldSpr40$Calcium,Nutrition_YieldSpr40$Iron,Nutrition_YieldSpr40$Vitamin_A),by=list(Nutrition_YieldSpr40$WPP),FUN=sum)
	names(Nutrition_YieldSpr40) = c("WPP","Selenium","Zinc","Protein","Omega_3","Calcium","Iron","Vitamin_A")
	#Yield When Managed at MSY%
	Nutrition_YieldMSY = subset(TradeOffs_Yield,select=c(species,WPP,Yield.MSY))
	Nutrition_YieldMSY = merge(NutritionSpeciesWPP,Nutrition_YieldMSY,all.x=TRUE,by=c("species","WPP"))
	Nutrition_YieldMSY$Selenium = Nutrition_YieldMSY$Selenium*(Nutrition_YieldMSY$Yield.MSY*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldMSY$Zinc = Nutrition_YieldMSY$Zinc*(Nutrition_YieldMSY$Yield.MSY*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldMSY$Protein = Nutrition_YieldMSY$Protein*(Nutrition_YieldMSY$Yield.MSY*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldMSY$Omega_3 = Nutrition_YieldMSY$Omega_3*(Nutrition_YieldMSY$Yield.MSY*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldMSY$Calcium = Nutrition_YieldMSY$Calcium*(Nutrition_YieldMSY$Yield.MSY*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldMSY$Iron = Nutrition_YieldMSY$Iron*(Nutrition_YieldMSY$Yield.MSY*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldMSY$Vitamin_A = Nutrition_YieldMSY$Vitamin_A*(Nutrition_YieldMSY$Yield.MSY*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldMSY$Yield.MSY=NULL
	Nutrition_YieldMSY = aggregate.data.frame(list(Nutrition_YieldMSY$Selenium,Nutrition_YieldMSY$Zinc,Nutrition_YieldMSY$Protein,Nutrition_YieldMSY$Omega_3,Nutrition_YieldMSY$Calcium,Nutrition_YieldMSY$Iron,Nutrition_YieldMSY$Vitamin_A),by=list(Nutrition_YieldMSY$WPP),FUN=sum)
	names(Nutrition_YieldMSY) = c("WPP","Selenium","Zinc","Protein","Omega_3","Calcium","Iron","Vitamin_A")
	#Yield When Managed at MEY%
	Nutrition_YieldMEY = subset(TradeOffs_Yield,select=c(species,WPP,Yield.MEY))
	Nutrition_YieldMEY = merge(NutritionSpeciesWPP,Nutrition_YieldMEY,all.x=TRUE,by=c("species","WPP"))
	Nutrition_YieldMEY$Selenium = Nutrition_YieldMEY$Selenium*(Nutrition_YieldMEY$Yield.MEY*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldMEY$Zinc = Nutrition_YieldMEY$Zinc*(Nutrition_YieldMEY$Yield.MEY*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldMEY$Protein = Nutrition_YieldMEY$Protein*(Nutrition_YieldMEY$Yield.MEY*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldMEY$Omega_3 = Nutrition_YieldMEY$Omega_3*(Nutrition_YieldMEY$Yield.MEY*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldMEY$Calcium = Nutrition_YieldMEY$Calcium*(Nutrition_YieldMEY$Yield.MEY*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldMEY$Iron = Nutrition_YieldMEY$Iron*(Nutrition_YieldMEY$Yield.MEY*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldMEY$Vitamin_A = Nutrition_YieldMEY$Vitamin_A*(Nutrition_YieldMEY$Yield.MEY*1000/100)   #covert kg to grams by multiplying by 1000, then divide by 100 since nutrition is in per 100 grams of fish
	Nutrition_YieldMEY$Yield.MEY=NULL
	Nutrition_YieldMEY = aggregate.data.frame(list(Nutrition_YieldMEY$Selenium,Nutrition_YieldMEY$Zinc,Nutrition_YieldMEY$Protein,Nutrition_YieldMEY$Omega_3,Nutrition_YieldMEY$Calcium,Nutrition_YieldMEY$Iron,Nutrition_YieldMEY$Vitamin_A),by=list(Nutrition_YieldMEY$WPP),FUN=sum)
	names(Nutrition_YieldMEY) = c("WPP","Selenium","Zinc","Protein","Omega_3","Calcium","Iron","Vitamin_A")
	#Now Group by Mineral/Vitamin, and Calculate the Percent Change from Status Quo For Each
	Selenium = data.frame(WPP = Nutrition_currentYield$WPP,CurrentYield = Nutrition_currentYield$Selenium,Spr20=Nutrition_YieldSpr20$Selenium,
		Spr30=Nutrition_YieldSpr30$Selenium,Spr40=Nutrition_YieldSpr40$Selenium,MSY=Nutrition_YieldMSY$Selenium,MEY=Nutrition_YieldMEY$Selenium)
	Zinc = data.frame(WPP = Nutrition_currentYield$WPP,CurrentYield = Nutrition_currentYield$Zinc,Spr20=Nutrition_YieldSpr20$Zinc,
		Spr30=Nutrition_YieldSpr30$Zinc,Spr40=Nutrition_YieldSpr40$Zinc,MSY=Nutrition_YieldMSY$Zinc,MEY=Nutrition_YieldMEY$Zinc)
	Protein = data.frame(WPP = Nutrition_currentYield$WPP,CurrentYield = Nutrition_currentYield$Protein,Spr20=Nutrition_YieldSpr20$Protein,
		Spr30=Nutrition_YieldSpr30$Protein,Spr40=Nutrition_YieldSpr40$Protein,MSY=Nutrition_YieldMSY$Protein,MEY=Nutrition_YieldMEY$Protein)
	Omega_3 = data.frame(WPP = Nutrition_currentYield$WPP,CurrentYield = Nutrition_currentYield$Omega_3,Spr20=Nutrition_YieldSpr20$Omega_3,
		Spr30=Nutrition_YieldSpr30$Omega_3,Spr40=Nutrition_YieldSpr40$Omega_3,MSY=Nutrition_YieldMSY$Omega_3,MEY=Nutrition_YieldMEY$Omega_3)
	Calcium = data.frame(WPP = Nutrition_currentYield$WPP,CurrentYield = Nutrition_currentYield$Calcium,Spr20=Nutrition_YieldSpr20$Calcium,
		Spr30=Nutrition_YieldSpr30$Calcium,Spr40=Nutrition_YieldSpr40$Calcium,MSY=Nutrition_YieldMSY$Calcium,MEY=Nutrition_YieldMEY$Calcium)
	Iron = data.frame(WPP = Nutrition_currentYield$WPP,CurrentYield = Nutrition_currentYield$Iron,Spr20=Nutrition_YieldSpr20$Iron,
		Spr30=Nutrition_YieldSpr30$Iron,Spr40=Nutrition_YieldSpr40$Iron,MSY=Nutrition_YieldMSY$Iron,MEY=Nutrition_YieldMEY$Iron)
	Vitamin_A = data.frame(WPP = Nutrition_currentYield$WPP,CurrentYield = Nutrition_currentYield$Vitamin_A,Spr20=Nutrition_YieldSpr20$Vitamin_A,
		Spr30=Nutrition_YieldSpr30$Vitamin_A,Spr40=Nutrition_YieldSpr40$Vitamin_A,MSY=Nutrition_YieldMSY$Vitamin_A,MEY=Nutrition_YieldMEY$Vitamin_A)
	write.table(Selenium,paste(PATH_highLevelSummaryPlots,"Selenium.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
	write.table(Zinc,paste(PATH_highLevelSummaryPlots,"Zinc.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
	write.table(Protein,paste(PATH_highLevelSummaryPlots,"Protein.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
	write.table(Omega_3,paste(PATH_highLevelSummaryPlots,"Omega_3.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
	write.table(Calcium,paste(PATH_highLevelSummaryPlots,"Calcium.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
	write.table(Iron,paste(PATH_highLevelSummaryPlots,"Iron.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
	write.table(Vitamin_A,paste(PATH_highLevelSummaryPlots,"Vitamin_A.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
	#Composition change - more simple table for paper EEZ
	NationalNutrientTotal = data.frame(Nutrient=c("Selenium","Zinc","Protein","Omega_3","Calcium","Iron","Vitamin_A"),
		CurrentYield = round(as.numeric(t(subset(Nutrition_currentYield,Nutrition_currentYield$WPP=="EEZ"))[-1,])),
		Spr20 = round(as.numeric(t(subset(Nutrition_YieldSpr20,Nutrition_YieldSpr20$WPP=="EEZ"))[-1,])),
		Spr30 = round(as.numeric(t(subset(Nutrition_YieldSpr30,Nutrition_YieldSpr30$WPP=="EEZ"))[-1,])),
		Spr40 = round(as.numeric(t(subset(Nutrition_YieldSpr40,Nutrition_YieldSpr40$WPP=="EEZ"))[-1,])),
		MSY = round(as.numeric(t(subset(Nutrition_YieldMSY,Nutrition_YieldMSY$WPP=="EEZ"))[-1,])),
		MEY = round(as.numeric(t(subset(Nutrition_YieldMEY,Nutrition_YieldMEY$WPP=="EEZ"))[-1,])))
	write.table(NationalNutrientTotal,paste(PATH_highLevelSummaryPlots,"NationalNutrientTotal.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
}

#See Excel Sheet "TradeOffPlots.xlsx" for more summary plots for the paper!!!!

if(plotStatusByTrophicLevel==TRUE)
{
	SpeciesTrophicLevel = read.table(paste(PATH_output,"LifeHistoryTable_FromFishBase.csv",sep="/"),header=TRUE,sep=",")
	SpeciesTrophicLevel = subset(SpeciesTrophicLevel,select=c(species,TrophicLevel,Vulnerability))
	SpeciesTrophicLevel$TrophicLevel[SpeciesTrophicLevel$species=="Atrobucca brevis"]=3.6
	SpeciesTrophicLevel$TrophicLevel[SpeciesTrophicLevel$species=="Epinephelus amblycephalus"]=3.8
	SpeciesTrophicLevel$TrophicLevel[SpeciesTrophicLevel$species=="Epinephelus bleekeri"]=3.9
	SpeciesTrophicLevel$TrophicLevel[SpeciesTrophicLevel$species=="Glaucosoma buergeri"]=4.2
	SpeciesTrophicLevel$TrophicLevel[SpeciesTrophicLevel$species=="Paracaesio gonzalesi"]=3.4
	SpeciesTrophicLevel$TrophicLevel[SpeciesTrophicLevel$species=="Plectropomus leopardus"]=4.4
	BenchmarkSummary = read.table(paste(PATH_output,"BenchmarkSummary.csv",sep="/"),header=TRUE,sep=",")
	BenchmarkSummary = merge(BenchmarkSummary,SpeciesTrophicLevel,all.x=TRUE,by=c("species"))
	png(paste(PATH_highLevelSummaryPlots,"Overfished_TrophicLevel_ScatterPlot.png",sep="/"),units="px",width=3200,height=3200,res=600)
	plot(BenchmarkSummary$TrophicLevel,BenchmarkSummary$B_over_Bspr20_median,xlab="Trophic Level",ylab="B/B@SPR20%")
	title("Trophic Level and LRP Stock Status")
	dev.off()
	cor(BenchmarkSummary$TrophicLevel,BenchmarkSummary$B_over_Bspr20_median,use="complete.obs") 
	png(paste(PATH_highLevelSummaryPlots,"Overfishing_TrophicLevel_ScatterPlot.png",sep="/"),units="px",width=3200,height=3200,res=600)
	plot(BenchmarkSummary$TrophicLevel,BenchmarkSummary$Fratio40_median,xlab="Trophic Level",ylab="F/F@SPR40%")
	title("Trophic Level and TRP Stock Status")
	dev.off()
	cor(BenchmarkSummary$TrophicLevel,BenchmarkSummary$Fratio40_median,use="complete.obs")
}












