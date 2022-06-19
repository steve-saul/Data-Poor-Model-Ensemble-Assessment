###################### Develop Poseidon Biological Inputs ############################
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
library(LBSPR)
library(zoo)
library(devtools)
#Note: need to first install Rtools from Rtools website. It is installed as an executable windows file, NOT installed as an R library.
#devtools::install_github("merrillrudd/LIME")
library(LIME)
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


#################  Input Parameters (most will be dynamic from a parameter file once updated) ###############################################
#You need to manually create these two folders below on your computer first. The OtherInputs folder should be provided to you together with this script
PATH_otherInput = "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/BiologicalModel/OtherInputs"
#Specify the location of your 7-zip software. You will need to download and install this Free software to extract the TNC data....
#...There is no work around for this because the file is password protected and R's build in upzip function doesn't work when the zip file is password protected.
PATH_7zip = "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/BiologicalModel/7-Zip"
#Here just specify the full paths and folder names for the remaining folders.
PATH_iFish = "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/BiologicalModel/iFish"
PATH_cpue =  "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/BiologicalModel/CPUE"
PATH_output = "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/BiologicalModel/Output"
PATH_summaryTablesPaper = "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/BiologicalModel/SummaryTablesPublication"
PATH_plots = "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/BiologicalModel/Plots"
PATH_historicalCatchPlots = "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/BiologicalModel/HistoricalCatchPlots"
PATH_catchMSY_diagnostics = "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/BiologicalModel/CatchMSY_Diagnostics"
PATH_catchMSY_plots = "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/BiologicalModel/CatchMSY_Plots"
PATH_projectionPlots = "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/BiologicalModel/ProjectionPlots"
PATH_projectionPlots_byProjectType = "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/BiologicalModel/ProjectionPlots_byProjectType"
PATH_benchmarkPlots = "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/BiologicalModel/BenchmarkPlots"
PATH_precentChangeF_histograms = "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/BiologicalModel/PercentChangeF_histograms"
PATH_summaryFanProjectionPlots = "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/BiologicalModel/SummaryFanProjectionPlots"
PATH_standardizedFanProjectionPlots = "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/BiologicalModel/StandardizedFanProjectionPlots"
PATH_highLevelSummaryPlots = "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/BiologicalModel/HighLevelSummaryPlots"
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
  "russelli","sebae","timorensis","vitta","japonicus","gonzalesi","kusakarii","stonei","xanthura","eriomma","lewisi","pinjalo","leopardus","maculatus","kaakan",
  "argyrogrammicus","filamentosus","flavipinnis","multidens","sieboldii","typus","zonatus","diacanthus","canadum","powelli","dumerili","rivoliana","barracuda",
  "forsteri","putnamae","nematophorus","albimarginata","mossambica","other","total")
#SpeciesForAnalysis_names = c("laticaudis")
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
#SpeciesForAnalysis_codes = c("504404580790866817")
FleetInventory_FileName = "FleetInventoryCapacity.csv"  #This file should be sent to you in the OtherInputs folder
#Specify the areas you want analyzed.
wholeCountry=TRUE       #If you also want the full EEZ analyzed with the pooled data, specify this as TRUE.
WPP = c(571,572,573,711,712,713,714,715,716,718)   #WPPs you want analyzed
#WPP = c(718,712)
#Flags to turn on and off different code sections
computeSPR_createInputFile=TRUE
setUpFilingSystem=TRUE   #WARNING: Setting this to "TRUE", will erase any output from a previous run of the code!!!!!
pullIfishData = FALSE
pullCPUEestimates = FALSE
runCatchMSY = TRUE
plotHistoricalCatches = FALSE
plotCatchMSY_results = TRUE
plotCatchMSY_diagnostics = TRUE
runLengthReconstruction = TRUE
runProjections = TRUE
runRebuildingProjectionsRegardlessOfStatus = TRUE
createKobePlots = TRUE
plotProjections = TRUE
plotProjectionsByProjectionType = TRUE
createFanPlots = TRUE
createStandardizedFanPlots = TRUE
summarizeOutputForPublications = TRUE
#Address for Indonesia database - 
#THESE ARE PULLED FROM ONLINE GOVERNMENT DATABASE. LOCATION REMOVED TO PROTECT DATA ACCESS. Please contact Ministry of Fisheries, Republic of Indonesia for data access
Sql_URL_File_name = ""
password_sqlFile=""  
sql_user = ""       
sql_password = ""  
#Set other variables 
fishingDays_fName = "fishing_days.csv"
CatchMSY_inputFileName = "StartingParametersCatchMSY.csv"         #Part of the OtherInputs folder. Note this file need to be filled out for each species and area you want analyzed.
LandingsUncertaintyScenarios = c(0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8)  #These are used as multipliers for difference scenarios about what we might know about the landings
CatchMSY_scenarioNames = c("Historical","Optimistic","Pessimistic","Historical_spr10","Optimistic_spr10","Pessimistic_spr10") #,"Historical_spr20","Optimistic_spr20","Pessimistic_spr20","Historical_spr30","Optimistic_spr30","Pessimistic_spr30")      #must match the length of the "CatchMSY_inputFileNames" vector above
CatchMSY_outfileName = "CatchMSY_Output.csv"    #Name of output file from CatchMSY
cMSY_optimisticScenario_yrsToCurrentDepletion = 20      #CatchMSY optimistic scenario timeframe
cMSY_pessimisticScenario_yrsToCurrentDepletion = 10     #CatchMSY pessimistic scenario timeframe
populationReconstructionMethods = c("SolveBaranov","LBSPR","BevertonHoltInstantaneousMortality","LIME","TropFishR")
#LifeHistoryScenarios = c("Linf_90Perc_lmax","Linf_equal_lmax")
LifeHistoryScenarios = c("Linf_90Perc_lmax")
yearsToProject = 50
projectionYearsToPlot = 20
ControlRule_SPR = 0.3
Limit_SPR = 0.2
ControlRule40_SPR = 0.4
yearToStartProjection = 2019
Steepness_Scenarios = c(0.7,0.8,0.9)
FishingMortaltiyProjectionScenarios=c("CurrentF","F_equal_0","F_at_SPR20","F_at_SPR30","F_at_75percSPR30","F_at_SPR40","F_rebuild_20YrsMSC_max","F_rebuild_Tmin_oneGeneration","F_rebuild_10YrsUS") #,"F_half","F_double","ConstantCatch")
ColorsFor_ProjectionF_ScenarioPlots = c("red","royalblue","green","grey","yellow","purple","orange","brown","turquoise","pink","wheat2","gold3")   #must be same length as vector above plus one additonal color.

######################## Setup Filing System for Outputs ####################################################################################
#make first letter of species names capital
CapStr = function(y)
{
  c = strsplit(y, " ")[[1]]
  paste(toupper(substring(c,1,1)),substring(c,2),sep="",collapse=" ")
}
SpeciesForAnalysis_names_Cap = sapply(SpeciesForAnalysis_names,CapStr)
names(SpeciesForAnalysis_names_Cap)=NULL
SpeciesFolderNames = paste(GenusForAnalysis_names,SpeciesForAnalysis_names_Cap,sep="")
#general plot folder for  catch at age, catch at legnth, and selectivity plots
if(setUpFilingSystem==TRUE)
{
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
  for(i in 1:length(SpeciesFolderNames))
  {
    if(dir.exists(paste(PATH_catchMSY_diagnostics,SpeciesFolderNames[i],sep="/"))==TRUE)
    {
      unlink(paste(PATH_catchMSY_diagnostics,SpeciesFolderNames[i],sep="/"),recursive=TRUE)
    }
    dir.create(paste(PATH_catchMSY_diagnostics,SpeciesFolderNames[i],sep="/"))
  }
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
  #Projection spaghetti plots - traditional, by species and area with all projection scenarios on the same figure
  if(dir.exists(PATH_projectionPlots)==TRUE)
  {
    unlink(PATH_projectionPlots,recursive=TRUE)
  }
  dir.create(PATH_projectionPlots)
  for(i in 1:length(SpeciesFolderNames))
  {
    if(dir.exists(paste(PATH_projectionPlots,SpeciesFolderNames[i],sep="/"))==TRUE)
    {
      unlink(paste(PATH_projectionPlots,SpeciesFolderNames[i],sep="/"),recursive=TRUE)
    }
    dir.create(paste(PATH_projectionPlots,SpeciesFolderNames[i],sep="/"))
  }
  #Projection spaghetti plots - with one plot for each projection type, and the different scenarios on that figure
  if(dir.exists(PATH_projectionPlots_byProjectType)==TRUE)
  {
    unlink(PATH_projectionPlots_byProjectType,recursive=TRUE)
  }
  dir.create(PATH_projectionPlots_byProjectType)
  for(i in 1:length(SpeciesFolderNames))
  {
    if(dir.exists(paste(PATH_projectionPlots_byProjectType,SpeciesFolderNames[i],sep="/"))==TRUE)
    {
      unlink(paste(PATH_projectionPlots_byProjectType,SpeciesFolderNames[i],sep="/"),recursive=TRUE)
    }
    dir.create(paste(PATH_projectionPlots_byProjectType,SpeciesFolderNames[i],sep="/"))
  }
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


####################### iFish Data: Pull From Server or Read From Local Disk ###############################################################
if(pullIfishData==TRUE)
{
   if(dir.exists(PATH_iFish)==TRUE)
   {
     unlink(PATH_iFish,recursive=TRUE)
   }
   dir.create(PATH_iFish)
   url_split = strsplit(Sql_URL_File_name,"/")[[1]]
   Sql_fileName = url_split[length(url_split)]
   Sql_fileName_unzipped = strsplit(Sql_fileName,".zip")[[1]]
   if(file.exists(paste(PATH_iFish,Sql_fileName,sep="/")))
   {
      file.remove(paste(PATH_iFish,Sql_fileName,sep="/"))
   }
   if(file.exists(paste(PATH_iFish,Sql_fileName_unzipped,sep="/")))
   {
      file.remove(paste(PATH_iFish,Sql_fileName_unzipped,sep="/"))
   }
   download.file(Sql_URL_File_name,paste(PATH_iFish,Sql_fileName,sep="/"))
   system(paste0(paste(PATH_7zip,"7z",sep="/")," x ",paste(PATH_iFish,Sql_fileName,sep="/")," -o",PATH_iFish," -p",password_sqlFile))
   drv = dbDriver("PostgreSQL")
   SQL_codeToRun = readLines(paste(PATH_iFish,Sql_fileName_unzipped,sep="/"),n=-1)
   SQL_codeToRun = paste0(SQL_codeToRun,collapse="")
   con = dbConnect(drv, user=sql_user, password=sql_password, dbname="postgres",host="localhost")
   SqlExecute = dbSendQuery(con,SQL_codeToRun)
   TableNames = dbGetQuery(con, "SELECT table_name FROM information_schema.tables WHERE table_schema='public'")
   Headers = list()
   for(i in 1:length(TableNames$table_name))
   {
      Headers[[TableNames$table_name[i]]] = names(dbReadTable(con,TableNames$table_name[i]))
      SqlExecute_saveTable = dbSendQuery(con,paste("COPY ",TableNames$table_name[i]," TO '",paste(paste(PATH_iFish,TableNames$table_name[i],sep="/"),".txt",sep=""),"';",sep=""))
      CharVec = rep("",length(names(dbReadTable(con,TableNames$table_name[i]))))
      Table_i = scan(paste(paste(PATH_iFish,TableNames$table_name[i],sep="/"),".txt",sep=""),sep="\t",fill=TRUE,what=as.list(CharVec))
      Table_i_df = data.frame(Table_i,stringsAsFactors=FALSE)
      names(Table_i_df) = names(dbReadTable(con,TableNames$table_name[i]))
      write.table(Table_i_df,paste(paste(PATH_iFish,TableNames$table_name[i],sep="/"),".txt",sep=""),col.names=TRUE,row.names=FALSE,sep="\t")
   }
   RPostgreSQL::dbDisconnect(con)
}

##################################### Fetch CPUE Estimates ######################################
#THESE ARE PULLED FROM ONLINE GOVERNMENT DATABASE. LOCATION REMOVED TO PROTECT DATA ACCESS. Please contact Ministry of Fisheries, Republic of Indonesia for data access
if(pullCPUEestimates==TRUE)
{
  RawDataAverageCPUE_SpeciesComposition_kgPerLanding = read.table("",header=TRUE,sep=",")
  write.table(RawDataAverageCPUE_SpeciesComposition_kgPerLanding,paste(PATH_cpue,"CPUE_Raw_kg_per_landing.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
  AverageCPUE_bySpecies_allAreas = read.table("",header=TRUE,sep=",")
  write.table(AverageCPUE_bySpecies_allAreas,paste(PATH_cpue,"AverageCPUE_bySpecies_allAreas.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
  AverageCPUE_bySpeciesWPP = read.table("",header=TRUE,sep=",")
  write.table(AverageCPUE_bySpeciesWPP,paste(PATH_cpue,"AverageCPUE_bySpeciesWPP.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
  AverageCPUE_WPP_AllSpeciesCombined = read.table("",header=TRUE,sep=",")
  write.table(AverageCPUE_WPP_AllSpeciesCombined,paste(PATH_cpue,"AverageCPUE_WPP_AllSpeciesCombined.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
}


####################################### Process IFish Data and Format for Analysis ###################

Lengths = read.table(paste(PATH_iFish,"ifish_sizing.txt",sep="/"),sep="\t",header=TRUE,as.is=TRUE,colClasses=c("character"),stringsAsFactors=FALSE)
Species = read.table(paste(PATH_iFish,"ifish_fish.txt",sep="/"),sep="\t",header=TRUE,as.is=TRUE,colClasses=c("character"),stringsAsFactors=FALSE)
Landings = read.table(paste(PATH_iFish,"ifish_deepslope.txt",sep="/"),sep="\t",header=TRUE,as.is=TRUE,colClasses=c("character"),stringsAsFactors=FALSE)
Boats = read.table(paste(PATH_iFish,"ifish_boat_pub.txt",sep="/"),sep="\t",header=TRUE,as.is=TRUE,colClasses=c("character"),stringsAsFactors=FALSE)
LocationPings = read.table(paste(PATH_iFish,"ifish_findmespot.txt",sep="/"),sep="\t",header=TRUE,as.is=TRUE,colClasses=c("character"),stringsAsFactors=FALSE)
BoatTracker = read.table(paste(PATH_iFish,"ifish_boat_tracker.txt",sep="/"),sep="\t",header=TRUE,as.is=TRUE,colClasses=c("character"),stringsAsFactors=FALSE)
AverageCPUE_WPP_AllSpeciesCombined = read.table(paste(PATH_cpue,"AverageCPUE_WPP_AllSpeciesCombined.csv",sep="/"),header=TRUE,sep=",",as.is=TRUE)

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
Boats$capacity_palka_m3 = as.numeric(Boats$capacity_palka_m3)
Boats$gt_estimate = as.numeric(Boats$gt_estimate)
Boats$gt_declared = as.numeric(Boats$gt_declared)
Boats$engine_hp1 = as.numeric(Boats$engine_hp1)
Boats$engine_hp2 = as.numeric(Boats$engine_hp2)
Boats$engine_hp3 = as.numeric(Boats$engine_hp3)
Boats$number_of_crew = as.numeric(Boats$number_of_crew)
Boats$fishing_area1 = as.numeric(Boats$fishing_area1)
Boats$fishing_area2 = as.numeric(Boats$fishing_area2)
Boats$fishing_area3 = as.numeric(Boats$fishing_area3)
LocationPings$latitude = as.numeric(LocationPings$latitude)
LocationPings$longitude = as.numeric(LocationPings$longitude)
LocationPings$daily_avg_latitude = as.numeric(LocationPings$daily_avg_latitude)
LocationPings$daily_avg_longitude = as.numeric(LocationPings$daily_avg_longitude)


######################### Compute CPUE and Days Fished Using Peter and Wawan's Method ########################################################
#Filter Data As Per Wawan's SQL Code for Computing CPUE
Species_cpue = Species[Species$species_id_number>0,]
Species_cpue = Species_cpue[Species_cpue$fish_species!="",]
Species_cpue = subset(Species_cpue,select=c(oid,fish_phylum,fish_class,fish_order,fish_family,fish_genus,fish_species,lmat,lopt,linf,lmax,lmatm,var_a,var_b,conversion_factor_tl2fl))
names(Species_cpue)[names(Species_cpue)=="oid"]="fish_id"
names(Species_cpue)[names(Species_cpue)=="fish_family"]="Family"
names(Species_cpue)[names(Species_cpue)=="fish_genus"]="Genus"
names(Species_cpue)[names(Species_cpue)=="fish_species"]="Species"
Lengths_cpue = Lengths[Lengths$data_quality==1,]
Landings$landing_date = as.Date(Landings$landing_date)
Landings_cpue = Landings[Landings$fishery_type=="Snapper" & Landings$doc_status=="Posted" & Landings$data_status=="Complete" & Landings$landing_date >= as.Date("2015-01-01"),]
Boats_cpue = Boats[Boats$gt_estimate>0,]
#Combine data files and estimate CPUE
names(Landings_cpue)[names(Landings_cpue)=="oid"]="landing_id"
Landings_cpue = subset(Landings_cpue,select=c(landing_id,landing_date,wpp1,fishing_gear,boat_id))
Landings_cpue = unique(Landings_cpue)
dateOnlyVec_vec = lapply(X=Landings_cpue$landing_date,FUN=function(x) {unlist(strsplit(as.character(x),"-"))})
dateOnlyVec_vec = do.call(rbind.data.frame,dateOnlyVec_vec)
rownames(dateOnlyVec_vec)=NULL
names(dateOnlyVec_vec) = c("Year_Landing","Month_Landing","Day_Landing")
dateOnlyVec_vec$Year_Landing = as.numeric(as.character(dateOnlyVec_vec$Year_Landing))
dateOnlyVec_vec$Month_Landing = as.numeric(as.character(dateOnlyVec_vec$Month_Landing))
dateOnlyVec_vec$Day_Landing = as.numeric(as.character(dateOnlyVec_vec$Day_Landing))
Landings_cpue = cbind(Landings_cpue,dateOnlyVec_vec)
Landings_cpue$landing_date=NULL
Lengths_cpue = subset(Lengths_cpue,select=c(landing_id,fish_id,codrs_picture_date,cm,length_type))
dateTimeVec_vec = lapply(X=Lengths_cpue$codrs_picture_date,FUN=function(x) {unlist(strsplit(as.character(x)," "))})
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
Lengths_cpue = cbind(Lengths_cpue,dateOnlyVec_vec)
Lengths_cpue$codrs_picture_date=NULL
Lengths_cpue = merge(Lengths_cpue,Landings_cpue,all.x=TRUE,by=c("landing_id"))
Lengths_cpue = Lengths_cpue[!is.na(Lengths_cpue$boat_id),]    #This merge drops everythign that was screened out of the landings data file as per Wawan's sql code.
Lengths_cpue = Lengths_cpue[!is.na(Lengths_cpue$Year_OutFishing),]
Lengths_cpue = merge(Lengths_cpue,Species_cpue,all.x=TRUE,by=c("fish_id"))
Lengths_cpue = Lengths_cpue[!is.na(Lengths_cpue$var_a),]  #Since above merge is not totally clean, get rid of records that didn't merge
Lengths_cpue$Weight_kg = Lengths_cpue$var_a*((Lengths_cpue$cm*Lengths_cpue$conversion_factor_tl2fl)^Lengths_cpue$var_b)*0.001          #Note also that most are in total length - very few are not, and it is qustionable what type of measurement they are, so I'm treating everything as if it were total length;
names(Boats_cpue)[names(Boats_cpue)=="oid"]="boat_id"
names(Boats_cpue)[names(Boats_cpue)=="fishing_gear"]="fishing_gear_boatDataSet"
Boats_cpue = subset(Boats_cpue,select=c(boat_id,program_type,fishing_gear_boatDataSet,registration_port,year_built,length_of_boat,gt_estimate,category,fishing_area1))
Boats_cpue = unique(Boats_cpue)
Lengths_cpue = merge(Lengths_cpue,Boats_cpue,all.x=TRUE,by=c("boat_id"))

#Add boat characteristics (i.e. boat gt size, category, WPP fishing in, etc.) to data
VesselCharacteristics = subset(Lengths_cpue,select=c(boat_id,fishing_gear_boatDataSet,length_of_boat,gt_estimate,category,fishing_area1))
VesselCharacteristics = unique(VesselCharacteristics)
names(VesselCharacteristics)[names(VesselCharacteristics)=="fishing_area1"]="WPP"
AverageCPUE_WPP_AllSpeciesCombined1 = AverageCPUE_WPP_AllSpeciesCombined
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
FishingDays = FishingDays[FishingDays$Year_OutFishing==(as.numeric(strsplit(date()," ")[[1]][5])-1),]
FishingDays = FishingDays[order(FishingDays$boat_id,FishingDays$Year_OutFishing,FishingDays$Month_OutFishing,FishingDays$Day_OutFishing),]
FishingDays = data.frame(table(FishingDays$boat_id))
names(FishingDays) = c("boat_id","DaysFishingLastYear")
VesselCharacteristics = subset(Boats,select=c(oid,fishing_gear,gt_estimate,category))
VesselCharacteristics = unique(VesselCharacteristics)
names(VesselCharacteristics)[names(VesselCharacteristics)=="oid"]="boat_id"
VesselCharacteristics$VesselCategory=""
VesselCharacteristics$VesselCategory[VesselCharacteristics$gt_estimate<5]="Nano"
VesselCharacteristics$VesselCategory[VesselCharacteristics$gt_estimate>=5 & VesselCharacteristics$gt_estimate<15]="Small"
VesselCharacteristics$VesselCategory[VesselCharacteristics$gt_estimate>=15 & VesselCharacteristics$gt_estimate<30]="Medium"
VesselCharacteristics$VesselCategory[VesselCharacteristics$gt_estimate>=30]="Large"
FishingDays = merge(FishingDays,VesselCharacteristics,all.x=TRUE,by=c("boat_id"))
AvgDaysFishedPerYear_withCategory = aggregate.data.frame(FishingDays$DaysFishingLastYear,by=list(FishingDays$category,FishingDays$fishing_gear,FishingDays$VesselCategory),FUN=mean,na.rm=TRUE)
names(AvgDaysFishedPerYear_withCategory) = c("category","fishing_gear","VesselSize","AvgDaysFished")
AvgDaysFishedPerYear_noCategory = aggregate.data.frame(FishingDays$DaysFishingLastYear,by=list(FishingDays$fishing_gear,FishingDays$VesselCategory),FUN=mean,na.rm=TRUE)
names(AvgDaysFishedPerYear_noCategory) = c("fishing_gear","VesselSize","AvgDaysFished")
AvgDaysFishedPerYear_withCategory_Table = AvgDaysFishedPerYear_withCategory
AvgDaysFishedPerYear_withCategory_Table$Size_Category = paste(AvgDaysFishedPerYear_withCategory_Table$VesselSize,AvgDaysFishedPerYear_withCategory_Table$category,sep=" - ")
AvgDaysFishedPerYear_withCategory_Table = reshape(AvgDaysFishedPerYear_withCategory_Table,v.names="AvgDaysFished",idvar="Size_Category",timevar="fishing_gear",direction="wide")
#***********************************************************************************************************************************************
#************************** Temporary Hard Coded Area With Wawan's Numbers: I don't know how he calculated these numbers ***********************
#***********************************************************************************************************************************************
gear2 <- c("Dropline", "Dropline", "Dropline", "Dropline", "Dropline", "Dropline", "Dropline", "Dropline",
           "Longline", "Longline", "Longline", "Longline", "Longline", "Longline", "Longline", "Longline",
           "Gillnet", "Gillnet", "Gillnet", "Gillnet", "Gillnet", "Gillnet", "Gillnet", "Gillnet",
           "Traps", "Traps", "Traps", "Traps", "Traps", "Traps", "Traps", "Traps",
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
AvgDaysFishedPerYear_withCategory = data.frame(fishing_gear=gear2,VesselSize=size2,FishingDays=fdays)
dateTimeVec_vec = lapply(X=AvgDaysFishedPerYear_withCategory$VesselSize,FUN=function(x) {unlist(strsplit(as.character(x)," "))})
dateTimeVec_vec = do.call(rbind.data.frame,dateTimeVec_vec)
rownames(dateTimeVec_vec)=NULL
names(dateTimeVec_vec) = c("VesselSize","category")
AvgDaysFishedPerYear_withCategory$VesselSize=NULL
AvgDaysFishedPerYear_withCategory = cbind(AvgDaysFishedPerYear_withCategory,dateTimeVec_vec)

#***********************************************************************************************************************************************
#***********************************************************************************************************************************************
################ Use CPUE and Days Fished With Complete Fleet Census Data and Average Dats Fished To Determine Total Catch #########################################

Boats = Boats[Boats$program_type=="Snapper",]
FullFleet = subset(Boats,select=c(oid,fishing_gear,gt_estimate,category,fishing_area1,fishing_area2,fishing_area3))
names(FullFleet)[names(FullFleet)=="oid"]="boat_id"
FullFleet$WPP = FullFleet$fishing_area1
FullFleet$WPP[is.na(FullFleet$fishing_area1)]=FullFleet$fishing_area2[is.na(FullFleet$fishing_area1)]
FullFleet$WPP[is.na(FullFleet$fishing_area1) & is.na(FullFleet$fishing_area2)]=FullFleet$fishing_area3[is.na(FullFleet$fishing_area1) & is.na(FullFleet$fishing_area2)]
FullFleet$fishing_area1=NULL
FullFleet$fishing_area2=NULL
FullFleet$fishing_area3=NULL
FullFleet = unique(FullFleet)
FullFleet$VesselCategory=""
FullFleet$VesselCategory[FullFleet$gt_estimate<5]="Nano"
FullFleet$VesselCategory[FullFleet$gt_estimate>=5 & FullFleet$gt_estimate<15]="Small"
FullFleet$VesselCategory[FullFleet$gt_estimate>=15 & FullFleet$gt_estimate<30]="Medium"
FullFleet$VesselCategory[FullFleet$gt_estimate>=30]="Large"
FullFleet = merge(FullFleet,CpueData,all.x=TRUE,by=c("VesselCategory","fishing_gear","WPP"))
names(AvgDaysFishedPerYear_withCategory)[names(AvgDaysFishedPerYear_withCategory)=="VesselSize"]="VesselCategory"
FullFleet = merge(FullFleet,AvgDaysFishedPerYear_withCategory,all.x=TRUE,by=c("fishing_gear","VesselCategory","category"))
FullFleet$TotalCatch_kg = FullFleet$kg_gt_day*FullFleet$FishingDays*FullFleet$gt_estimate
FullFleet$TotalCatch_mT = FullFleet$TotalCatch_kg*0.001
CatchPerWPP = aggregate.data.frame(FullFleet$TotalCatch_mT,by=list(FullFleet$WPP),FUN=sum,na.rm=TRUE)
names(CatchPerWPP) = c("WPP","TotalCatch_mT")

################## Compute Amount of Catch By Fish Family and Species in Each WPP #############################
#Create group common names for each fish family
Lengths_cpue$FamilyCommonName=""
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Carangidae"]="Jacks_Mackerels"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Scombridae"]="Jacks_Mackerels"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Emmelichthydae"]="Other"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Emmelichthyidae"]="Other"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Epinephelidae"]="Grouper"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Glaucosomatidae"]="Other"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Haemulidae"]="Grunts"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Holocentridae"]="Other"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Lethrinidae"]="Emperor"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Lutjanidae"]="Snapper"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Nemipteridae"]="SmallPelagic"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Priacanthidae"]="Bigeyes"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Rachycentridae"]="Cobia"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Sciaenidae"]="Drums_Croaker"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Sparidae"]="Porgies"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Sphyraenidae"]="Barracuda"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Acanthuridae"]="Surgeonfish"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Albulidae"]="Bonefish"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Malacanthidae"]="Other"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Balistidae"]="Triggerfish"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Priacanthidae"]="Bigeyes"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Coryphaenidae"]="Dolphinfish"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Serranidae"]="Grouper"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Glaucosomatidae"]="Other"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Lampridae"]="Other"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Istiophoridae"]="LargePelagic"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="miscoded"]="Other"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Ariidae"]="Other"
Lengths_cpue$FamilyCommonName[Lengths_cpue$Family=="Xiphiidae"]="LargePelagic"

#Calculate the Proportion of Each Family Within Each WPP
Lengths_cpue_WPPnotMissing = Lengths_cpue[!is.na(Lengths_cpue$wpp1) & Lengths_cpue$FamilyCommonName!="",]
FamilySampledPerWPP = aggregate.data.frame(Lengths_cpue_WPPnotMissing$Weight_kg,by=list(Lengths_cpue_WPPnotMissing$FamilyCommonName,Lengths_cpue_WPPnotMissing$wpp1),FUN=sum,na.rm=TRUE)
names(FamilySampledPerWPP) = c("FamilyCommonName","WPP","Weight_kg")
FamilySampled = aggregate.data.frame(Lengths_cpue_WPPnotMissing$Weight_kg,by=list(Lengths_cpue_WPPnotMissing$wpp1),FUN=sum,na.rm=TRUE)
names(FamilySampled) = c("WPP","Weight_kg_WPP")
FamilySampledPerWPP = merge(FamilySampledPerWPP,FamilySampled,all.x=TRUE,by=c("WPP"))
FamilySampledPerWPP$ProportionOfFamilyInWPP = FamilySampledPerWPP$Weight_kg/FamilySampledPerWPP$Weight_kg_WPP
FamilySampledPerWPP$Weight_kg_WPP=NULL
FamilySampledPerWPP$Weight_kg=NULL
FamilySampledPerWPP = reshape(FamilySampledPerWPP,v.names="ProportionOfFamilyInWPP",idvar="FamilyCommonName",timevar="WPP",direction="wide")
FamilySampledPerWPP[is.na(FamilySampledPerWPP)]=0
#Apply the Proportion of Each Family Within Each WPP to Sea Around Us Data
SeaAroundUsIndo = read.table(paste(PATH_otherInput,"SeaAroundUs/SAU EEZ 361,937,938 v47-1/SAU EEZ 361,937,938 v47-1.csv",sep="/"),header=TRUE,sep=",",stringsAsFactors=FALSE)
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
FractionSpeciesInFamilyAndWPP = aggregate.data.frame(Lengths_cpue_WPPnotMissing$Weight_kg,by=list(Lengths_cpue_WPPnotMissing$fish_id,Lengths_cpue_WPPnotMissing$FamilyCommonName,Lengths_cpue_WPPnotMissing$wpp1),FUN=sum)
names(FractionSpeciesInFamilyAndWPP) = c("fish_id","FamilyCommonName","WPP","Weight_kg")
FractionFamilyAndWPP = aggregate.data.frame(Lengths_cpue_WPPnotMissing$Weight_kg,by=list(Lengths_cpue_WPPnotMissing$FamilyCommonName,Lengths_cpue_WPPnotMissing$wpp1),FUN=sum)
names(FractionFamilyAndWPP) = c("FamilyCommonName","WPP","Weight_kg_FamilyWPP")
FractionSpeciesInFamilyAndWPP = merge(FractionSpeciesInFamilyAndWPP,FractionFamilyAndWPP,all.x=TRUE,by=c("FamilyCommonName","WPP"))
FractionSpeciesInFamilyAndWPP$ProportionOfSpeciesInFamilyAndWPP = FractionSpeciesInFamilyAndWPP$Weight_kg/FractionSpeciesInFamilyAndWPP$Weight_kg_FamilyWPP
FractionSpeciesInFamilyAndWPP$Weight_kg=NULL
FractionSpeciesInFamilyAndWPP$Weight_kg_FamilyWPP=NULL
FractionSpeciesInFamilyAndWPP = reshape(FractionSpeciesInFamilyAndWPP,v.names="ProportionOfSpeciesInFamilyAndWPP",idvar="fish_id",timevar="WPP",direction="wide")
FractionSpeciesInFamilyAndWPP[is.na(FractionSpeciesInFamilyAndWPP)]=0
SeaAroundUsIndo_bySpecies = merge(SeaAroundUsIndo_byFamily,FractionSpeciesInFamilyAndWPP,by=c("FamilyCommonName"))
FamilyCatchMatrix = SeaAroundUsIndo_bySpecies[,grep("Catch_mT",names(SeaAroundUsIndo_bySpecies))]
SpeciesProportionMatrix = SeaAroundUsIndo_bySpecies[,grep("ProportionOfSpeciesInFamilyAndWPP",names(SeaAroundUsIndo_bySpecies))]
names(FamilyCatchMatrix) = gsub("Catch_mT.","",names(FamilyCatchMatrix))
names(SpeciesProportionMatrix) = gsub("ProportionOfSpeciesInFamilyAndWPP.","",names(SpeciesProportionMatrix))
SpeciesProportionMatrix = SpeciesProportionMatrix[,names(FamilyCatchMatrix)]
SpeciesCatchMatrix_mT = FamilyCatchMatrix*SpeciesProportionMatrix
RowInfo = subset(SeaAroundUsIndo_bySpecies,select=c(FamilyCommonName,Year,fish_id))
SpeciesCatchMatrix_mT = cbind(RowInfo,SpeciesCatchMatrix_mT)
#Estimate Current Catch by Family, Species, and WPP - use this constant value for 2015 throgh 2018 four years
FamilySampled = aggregate.data.frame(Lengths_cpue_WPPnotMissing$Weight_kg,by=list(Lengths_cpue_WPPnotMissing$FamilyCommonName),FUN=sum,na.rm=TRUE)
names(FamilySampled) = c("FamilyCommonName","Weight_kg_Family")
FamilySampled$ProportionSampled_Family = FamilySampled$Weight_kg_Family/sum(FamilySampled$Weight_kg_Family)
FamilySampled$Catch_Family_mT = FamilySampled$ProportionSampled_Family*sum(FullFleet$TotalCatch_mT,na.rm=TRUE)
FamilySampled$Weight_kg_Family=NULL
FamilySampled$ProportionSampled_Family=NULL
CatchPerWPP = merge(FamilySampledPerWPP,FamilySampled,all.x=TRUE,by=c("FamilyCommonName"))
CompleteCatches = aggregate.data.frame(FullFleet$TotalCatch_mT,by=list(FullFleet$WPP),FUN=sum,na.rm=TRUE)
names(CompleteCatches) = c("WPP","TotalCatch_mT")
CompleteCatches = merge(CompleteCatches,data.frame(WPP,flag=1),all.x=TRUE,by=c("WPP"))
CompleteCatches = CompleteCatches[!is.na(CompleteCatches$flag),]
CompleteCatches$flag=NULL
CatchPerWPP_ProportionFamilyMatrix = CatchPerWPP[,grep("ProportionOfFamilyInWPP",names(CatchPerWPP))]
WPPinCODRS_data = names(CatchPerWPP_ProportionFamilyMatrix)
tonnesFamily_rep = do.call("rbind", replicate((dim(CatchPerWPP_ProportionFamilyMatrix)[1]),data.frame(t(data.frame(CompleteCatches$TotalCatch_mT))),simplify = FALSE))
CatchPerWPP_FamilyMatrix = CatchPerWPP_ProportionFamilyMatrix*tonnesFamily_rep
names(CatchPerWPP_FamilyMatrix) = gsub("ProportionOfFamilyInWPP.","",names(CatchPerWPP_FamilyMatrix))
CatchPerWPP_Family = cbind(data.frame(FamilyCommonName=CatchPerWPP$FamilyCommonName),CatchPerWPP_FamilyMatrix)
CatchPerWPP_Species = merge(FractionSpeciesInFamilyAndWPP,CatchPerWPP_Family,all.x=TRUE,by=c("FamilyCommonName"))
CatchPerWPP_ProportionSpeciesInFamilyMatrix = CatchPerWPP_Species[,grep("ProportionOfSpeciesInFamilyAndWPP",names(CatchPerWPP_Species))]
names(CatchPerWPP_ProportionSpeciesInFamilyMatrix) = gsub("ProportionOfSpeciesInFamilyAndWPP.","",names(CatchPerWPP_ProportionSpeciesInFamilyMatrix))
CatchPerWPP_FamilyMatrix = CatchPerWPP_Species[,names(CatchPerWPP_ProportionSpeciesInFamilyMatrix)]
CatchPerWPP_Species = CatchPerWPP_ProportionSpeciesInFamilyMatrix*CatchPerWPP_FamilyMatrix
CatchPerWPP_Species = cbind(subset(FractionSpeciesInFamilyAndWPP,select=c(FamilyCommonName,fish_id)),CatchPerWPP_Species)
Years = data.frame(Year=c(2015,2016,2017,2018))
CatchPerWPP_Species = merge(CatchPerWPP_Species,Years)
CatchPerWPP_Species = CatchPerWPP_Species[,names(SpeciesCatchMatrix_mT)]
SpeciesCatchMatrix_mT = rbind(SpeciesCatchMatrix_mT,CatchPerWPP_Species)


################################# Plot Reconstructed Landings - RAW #######################################
if(plotHistoricalCatches==TRUE)
{
  if(wholeCountry==TRUE)
  {
    WPP = c(571,712,713,718,572,573,711,715,714,716)
  }
  Colors = c("blue","red","green","orange","purple","black","brown","pink","gray","gold")
  for(i in 1:length(SpeciesForAnalysis_codes))
  {
    if(SpeciesForAnalysis_codes[i]!=-999 & SpeciesForAnalysis_codes[i]!=-9999)
    {
      Sub = SpeciesCatchMatrix_mT[SpeciesCatchMatrix_mT$fish_id==SpeciesForAnalysis_codes[i],]
      Sub$FamilyCommonName=NULL
      Sub$fish_id=NULL
      FolderNameSpp = paste(as.character(GenusForAnalysis_names[i]),CapStr(as.character(SpeciesForAnalysis_names[i])),sep="")
      png(paste(paste(PATH_historicalCatchPlots,FolderNameSpp,sep="/"),paste(paste(GenusForAnalysis_names[i],SpeciesForAnalysis_names[i],sep="_"),"_Historical_RAW.png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
      for(j in 1:length(WPP))
      {
        if(j==1)
        {
          plot(Sub$Year,Sub[,as.character(WPP[j])],xlab="Year",ylab="Metric Tons",type="l",lwd=2,col=Colors[j],ylim=c(0,max(Sub[,2:length(Sub)])),cex.axis=1.1,cex.lab=1.1)
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
}


################ Scale Historical Data to Current Split and Re-Plot Reconstructed Landings #######################################
#The reason that this historical scaling is probably needed is because the Sea Around Us Data likely also has shallow water species included, and we are only modeling the deep slope seamount species;
LastYearHistorical = SpeciesCatchMatrix_mT[SpeciesCatchMatrix_mT$Year==2014,]
LastYearHistorical$FamilyCommonName=NULL
LastYearHistorical$Year=NULL
LastYearHistorical = LastYearHistorical[order(LastYearHistorical$fish_id),]
FirstYearModern = SpeciesCatchMatrix_mT[SpeciesCatchMatrix_mT$Year==2015,]
FirstYearModern$FamilyCommonName=NULL
FirstYearModern$Year=NULL
HistoricalDataSpecies = data.frame(fish_id=LastYearHistorical$fish_id)
HistoricalDataSpecies$flag=1
FirstYearModern = merge(FirstYearModern,HistoricalDataSpecies,all.x=TRUE,by=c("fish_id"))
FirstYearModern = FirstYearModern[FirstYearModern$flag==1,]
FirstYearModern$flag=NULL
FirstYearModern = FirstYearModern[!is.na(FirstYearModern$fish_id),]
FirstYearModern = FirstYearModern[order(FirstYearModern$fish_id),]
Scalers = FirstYearModern[,2:length(FirstYearModern)]/LastYearHistorical[,2:length(LastYearHistorical)]
Scalers[Scalers==NaN]=0
Scalers[is.na(Scalers)]=0
names(Scalers) = paste("Scaling",names(Scalers),sep="_")
Scalers = cbind(data.frame(fish_id=FirstYearModern$fish_id,Scalers))    #note to self: need to add in the fish_id variable back here, then merge with the historical data, sort both, pull apart, and apply the correction factors. Then re-plot
HistoricalCatch = SpeciesCatchMatrix_mT[SpeciesCatchMatrix_mT$Year<=2014,]
HistoricalCatch = merge(HistoricalCatch,Scalers,all.x=TRUE,by=c("fish_id"))
HistoricalCatch = HistoricalCatch[!is.na(HistoricalCatch$Scaling_712),]
HistoricalCatch_Matrix = HistoricalCatch[,as.character(WPP)]
Scalers_Matrix = HistoricalCatch[,grep("Scaling_",names(HistoricalCatch))]
names(Scalers_Matrix) = gsub("Scaling_","",names(Scalers_Matrix))
Scalers_Matrix = Scalers_Matrix[,names(HistoricalCatch_Matrix)]
HistoricalCatch_Matrix_Scaled = HistoricalCatch_Matrix*Scalers_Matrix
HistoricalCatch_Scaled = cbind(subset(HistoricalCatch,select=c(fish_id,FamilyCommonName,Year)),HistoricalCatch_Matrix_Scaled)
HistoricalCatch_Scaled = HistoricalCatch_Scaled[order(HistoricalCatch_Scaled$fish_id,HistoricalCatch_Scaled$Year),]
ModernCatch = SpeciesCatchMatrix_mT[SpeciesCatchMatrix_mT$Year>2014,]
ModernCatch = ModernCatch[,names(HistoricalCatch_Scaled)]
SpeciesCatch_mT_Scaled = rbind(HistoricalCatch_Scaled,ModernCatch)
SpeciesCatch_mT_Scaled = SpeciesCatch_mT_Scaled[order(SpeciesCatch_mT_Scaled$fish_id,SpeciesCatch_mT_Scaled$Year),]

#Add groups for all species per year and then all others species not broken out and add to the master dataset of catches
SpeciesBrokenOut = data.frame(fish_id=as.character(SpeciesForAnalysis_codes),stringsAsFactors=FALSE)
SpeciesBrokenOut = data.frame(fish_id=SpeciesBrokenOut[as.numeric(SpeciesBrokenOut$fish_id)>0,])
SpeciesBrokenOut$isBrokenOut = 1
SpeciesCatch_mT_Scaled = merge(SpeciesCatch_mT_Scaled,SpeciesBrokenOut,all.x=TRUE,by=c("fish_id"))
SpeciesCatch_mT_Scaled_Other = subset(SpeciesCatch_mT_Scaled,is.na(SpeciesCatch_mT_Scaled$isBrokenOut))
SpeciesCatch_mT_Scaled$isBrokenOut=NULL
SpeciesCatch_mT_Scaled_Other$isBrokenOut=NULL
SpeciesCatch_mT_Scaled_Other$FamilyCommonName=NULL
SpeciesCatch_mT_Scaled_Other$fish_id=NULL
if(dim(SpeciesCatch_mT_Scaled_Other)[1]!=0)
{
  SpeciesCatch_mT_Scaled_Other = aggregate.data.frame(SpeciesCatch_mT_Scaled_Other,by=list(SpeciesCatch_mT_Scaled_Other$Year),FUN=sum)
  SpeciesCatch_mT_Scaled_Other$Year=NULL
  names(SpeciesCatch_mT_Scaled_Other)[names(SpeciesCatch_mT_Scaled_Other)=="Group.1"]="Year"
  SpeciesCatch_mT_Scaled_Other$fish_id=as.character(-999)
  SpeciesCatch_mT_Scaled_Other$FamilyCommonName="Other"
  SpeciesCatch_mT_Scaled_Other = SpeciesCatch_mT_Scaled_Other[,names(SpeciesCatch_mT_Scaled)]
}
SpeciesCatch_mT_Scaled_All = SpeciesCatch_mT_Scaled
SpeciesCatch_mT_Scaled_All$FamilyCommonName=NULL
SpeciesCatch_mT_Scaled_All$fish_id=NULL
SpeciesCatch_mT_Scaled_All = aggregate.data.frame(SpeciesCatch_mT_Scaled_All,by=list(SpeciesCatch_mT_Scaled_All$Year),FUN=sum)
SpeciesCatch_mT_Scaled_All$Year=NULL
names(SpeciesCatch_mT_Scaled_All)[names(SpeciesCatch_mT_Scaled_All)=="Group.1"]="Year"
SpeciesCatch_mT_Scaled_All$fish_id=as.character(-9999)
SpeciesCatch_mT_Scaled_All$FamilyCommonName="All"
SpeciesCatch_mT_Scaled_All = SpeciesCatch_mT_Scaled_All[,names(SpeciesCatch_mT_Scaled)]
SpeciesCatch_mT_Scaled = rbind(SpeciesCatch_mT_Scaled,SpeciesCatch_mT_Scaled_Other,SpeciesCatch_mT_Scaled_All)
SpeciesCatch_mT_Scaled = SpeciesCatch_mT_Scaled[order(SpeciesCatch_mT_Scaled$fish_id,SpeciesCatch_mT_Scaled$Year),]
GenusSpeciesID = subset(Species_cpue,select=c(fish_id,Genus,Species))
GenusSpeciesID = unique(GenusSpeciesID)
SpeciesCatch_mT_Scaled = merge(SpeciesCatch_mT_Scaled,GenusSpeciesID,all.x=TRUE,by=c("fish_id"))
write.table(SpeciesCatch_mT_Scaled,paste(PATH_output,"ExtrapolatedLandings.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")

######################## Add Zeros to Years Where Some Species Were Not Caught ############################
#If this is not done, Catch MSY will mess up
fish_id_analyzed = unique(SpeciesCatch_mT_Scaled$fish_id)
year_start = min(SpeciesCatch_mT_Scaled$Year)
year_end = max(SpeciesCatch_mT_Scaled$Year)
for(i in 1:length(fish_id_analyzed))
{
  Sub = SpeciesCatch_mT_Scaled[SpeciesCatch_mT_Scaled$fish_id==fish_id_analyzed[i],]
  fish_id_i = unique(Sub$fish_id)
  FamilyCommonName_i = unique(Sub$FamilyCommonName)
  Genus_i = unique(Sub$Genus)
  Species_i = unique(Sub$Species)
  toMerge = data.frame(Year=seq(year_start,year_end,by=1))
  Sub = merge(toMerge,Sub,all.x=TRUE,by=c("Year"))
  Sub$fish_id[is.na(Sub$fish_id)]=fish_id_i
  Sub$FamilyCommonName[is.na(Sub$FamilyCommonName)]=FamilyCommonName_i
  Sub$Genus[is.na(Sub$Genus)]=Genus_i
  Sub$Species[is.na(Sub$Species)]=Species_i
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
SpeciesCatch_mT_Scaled = Answer
write.table(SpeciesCatch_mT_Scaled,paste(PATH_output,"ExtrapolatedLandings.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")

################################# Plot Scaled Historical Landings #########################################
if(plotHistoricalCatches==TRUE)
{
  if(wholeCountry==TRUE)
  {
    WPP = c(571,712,713,718,572,573,711,715,714,716)
  }
  Colors = c("blue","red","green","orange","purple","black","brown","pink","gray","gold")
  for(i in 1:length(SpeciesForAnalysis_codes))
  {
    Sub = SpeciesCatch_mT_Scaled[SpeciesCatch_mT_Scaled$fish_id==SpeciesForAnalysis_codes[i],]
    Sub$FamilyCommonName=NULL
    Sub$fish_id=NULL
    Sub$Genus=NULL
    Sub$Species=NULL
    FolderNameSpp = paste(as.character(GenusForAnalysis_names[i]),CapStr(as.character(SpeciesForAnalysis_names[i])),sep="")
    png(paste(paste(PATH_historicalCatchPlots,FolderNameSpp,sep="/"),paste(paste(GenusForAnalysis_names[i],SpeciesForAnalysis_names[i],sep="_"),"_Historical_Scaled.png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
    for(j in 1:length(WPP))
    {
      if(j==1)
      {
        plot(Sub$Year,Sub[,as.character(WPP[j])],xlab="Year",ylab="Metric Tons",type="l",lwd=2,col=Colors[j],ylim=c(0,max(Sub[,2:length(Sub)])),cex.axis=1.1,cex.lab=1.1)
        title("Reconstructed Landings")
        mtext(paste(GenusForAnalysis_names[i],SpeciesForAnalysis_names[i],sep=" "))
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
  #Function to compute length-based SPR
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
  #Read fecundity values
  FecundityValues = read.table(paste(PATH_otherInput,"FecundityValues_filledIn_top50spp.csv",sep="/"),header=TRUE,sep=",")
  #Read and prepare last year (or most recent 150 observations) of length data for input
  YearMonthDay = subset(Lengths_cpue,select=c(Year_Landing,Month_Landing,Day_Landing))
  YearMonthDay = unique(YearMonthDay)
  YearMonthDay = YearMonthDay[order(YearMonthDay$Year_Landing,YearMonthDay$Month_Landing,YearMonthDay$Day_Landing),]
  YearMonthDay$Date = as.Date(paste(YearMonthDay$Month_Landing,YearMonthDay$Day_Landing,YearMonthDay$Year_Landing,sep="/"),format="%m/%d/%Y")
  LastDate = YearMonthDay$Date[dim(YearMonthDay)[1]]
  YearMonthDay$DaysPassed = as.numeric(LastDate - YearMonthDay$Date)
  Lengths_cpue = merge(Lengths_cpue,YearMonthDay,all.x=TRUE,by=c("Year_Landing","Month_Landing","Day_Landing"))
  counter=1
  for(i in 1:length(SpeciesForAnalysis_codes))
  {
    Params_i = Species_cpue[Species_cpue$fish_id==SpeciesForAnalysis_codes[i],]
    for(j in 1:length(WPP))
    {
       Sub = Lengths_cpue[Lengths_cpue$fish_id==SpeciesForAnalysis_codes[i] & Lengths_cpue$wpp1==WPP[j],]
       if(dim(Sub)[1]<150)
       {
          next
       }
       if(dim(Sub)[1]>=150)
       {
          Sub_lastYear = Sub[Sub$DaysPassed<=365,]
          if(dim(Sub_lastYear)[1]<150)
          {
              Sub$Counter = seq(dim(Sub)[1],1)
              Sub_lastYear = Sub[Sub$Counter<=150,]
          }
          Mnat = 10^(0.566-(0.718*log10(Params_i$linf))+0.02*20)
          VBG_K = (Mnat*Params_i$lopt)/(3*(Params_i$linf-Params_i$lopt))
          Sub_lastYear$cm_bin5 = 5*round(Sub_lastYear$cm/5)
          tmp = table(as.vector(Sub_lastYear$cm_bin5))
          lc = as.numeric(names(tmp)[tmp==max(tmp)])
          BevHoltEqOut = bheq(len=Sub_lastYear$cm,Linf=Params_i$linf,K=VBG_K,Lc=lc[1],La=Params_i$linf,nboot = 200)
          Z = BevHoltEqOut$z[2]
          F_est = Z - Mnat
          LenVec = Sub_lastYear$cm
          LenVec = LenVec[!is.na(LenVec)]
          SPR = computeLengthSPR(Mnat,Params_i$linf,VBG_K,Params_i$var_a,Params_i$var_b,Params_i$lmat,F_est,LenVec)
          if(counter==1)
          {
            Answer = data.frame(fish_id=SpeciesForAnalysis_codes[i],WPP=as.character(WPP[j]),spr_current=SPR)
          }
          if(counter>1)
          {
            Answer = rbind(Answer,data.frame(fish_id=SpeciesForAnalysis_codes[i],WPP=as.character(WPP[j]),spr_current=SPR))
          }
          counter = counter + 1
       }
    }
  #Now for all WPP combined
   Sub = Lengths_cpue[Lengths_cpue$fish_id==SpeciesForAnalysis_codes[i],]
   if(dim(Sub)[1]<150)
   {
      next
   }
   if(dim(Sub)[1]>=150)
   {
      Sub_lastYear = Sub[Sub$DaysPassed<=365,]
      if(dim(Sub_lastYear)[1]<150)
      {
          Sub$Counter = seq(dim(Sub)[1],1)
          Sub_lastYear = Sub[Sub$Counter<=150,]
      }
      Mnat = 10^(0.566-(0.718*log10(Params_i$linf))+0.02*20)
      VBG_K = (Mnat*Params_i$lopt)/(3*(Params_i$linf-Params_i$lopt))
      Sub_lastYear$cm_bin5 = 5*round(Sub_lastYear$cm/5)
      tmp = table(as.vector(Sub_lastYear$cm_bin5))
      lc = as.numeric(names(tmp)[tmp==max(tmp)])
      BevHoltEqOut = bheq(len=Sub_lastYear$cm,Linf=Params_i$linf,K=VBG_K,Lc=lc[1],La=Params_i$linf,nboot = 200)
      Z = BevHoltEqOut$z[2]
      F_est = Z - Mnat
      LenVec = Sub_lastYear$cm
      LenVec = LenVec[!is.na(LenVec)]
      SPR = computeLengthSPR(Mnat,Params_i$linf,VBG_K,Params_i$var_a,Params_i$var_b,Params_i$lmat,F_est,LenVec)
      if(counter==1)
      {
        Answer = data.frame(fish_id=SpeciesForAnalysis_codes[i],WPP="All",spr_current=SPR)
      }
      if(counter>1)
      {
        Answer = rbind(Answer,data.frame(fish_id=SpeciesForAnalysis_codes[i],WPP="All",spr_current=SPR))
      }
      counter = counter + 1
      if(i%%5==0)
      {
        print(paste("Working on species ",i," of ",length(SpeciesForAnalysis_codes),"...",sep=""))
        flush.console()
      }
   }
  }
  Answer$spr_current[Answer$spr_current>1]=1
  IdGenusSpecies = subset(Species_cpue,select=c(fish_id,Genus,Species))
  Answer = merge(Answer,IdGenusSpecies,all.x=TRUE,by=c("fish_id"))
  Answer$fish_id=NULL
  InputFile = merge(data.frame(Scenario=CatchMSY_scenarioNames),Answer)
  InputFile = InputFile[c("Scenario","Genus","Species","WPP","spr_current")]
  InputFile$l0_low=0.8
  InputFile$l0_up=1.0
  InputFile$l0_step=0
  InputFile$lt_low[InputFile$Scenario=="Historical"]=InputFile$spr_current[InputFile$Scenario=="Historical"]-0.1
  InputFile$lt_up[InputFile$Scenario=="Historical"]=InputFile$spr_current[InputFile$Scenario=="Historical"]+0.1
  InputFile$lt_low[InputFile$Scenario=="Optimistic"]=InputFile$spr_current[InputFile$Scenario=="Optimistic"]-0.1
  InputFile$lt_up[InputFile$Scenario=="Optimistic"]=InputFile$spr_current[InputFile$Scenario=="Optimistic"]+0.1
  InputFile$lt_low[InputFile$Scenario=="Pessimistic"]=InputFile$spr_current[InputFile$Scenario=="Pessimistic"]-0.1
  InputFile$lt_up[InputFile$Scenario=="Pessimistic"]=InputFile$spr_current[InputFile$Scenario=="Pessimistic"]+0.1
  InputFile$lt_low[InputFile$Scenario=="Historical_spr10"]=0.05
  InputFile$lt_up[InputFile$Scenario=="Historical_spr10"]=0.15
  InputFile$lt_low[InputFile$Scenario=="Optimistic_spr10"]=0.05
  InputFile$lt_up[InputFile$Scenario=="Optimistic_spr10"]=0.15
  InputFile$lt_low[InputFile$Scenario=="Pessimistic_spr10"]=0.05
  InputFile$lt_up[InputFile$Scenario=="Pessimistic_spr10"]=0.15
  InputFile$l0_low[InputFile$Scenario=="Optimistic"]=0.8
  InputFile$l0_low[InputFile$Scenario=="Pessimistic"]=0.8
  InputFile$lt_low[InputFile$lt_low<0]=0
  InputFile$lt_refyr=2018
  InputFile$sigv=0
  InputFile$r_dist="unif"
  InputFile$r_low=0.05
  InputFile$r_up=1
  InputFile$r_mean=0
  InputFile$r_sd=0
  InputFile$nsims=30000
  InputFile$grout=2
  InputFile = merge(InputFile,FecundityValues,all.x=TRUE,by=c("Genus","Species"))
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
  for(i in 1:length(CatchMSY_scenarioNames))
  {
    if(CatchMSY_scenarioNames[i]=="Historical" | CatchMSY_scenarioNames[i]=="Historical_spr10" | CatchMSY_scenarioNames[i]=="Historical_spr20" | CatchMSY_scenarioNames[i]=="Historical_spr30")
    {
      Catch = SpeciesCatch_mT_Scaled
    }
    if(CatchMSY_scenarioNames[i]=="Optimistic" | CatchMSY_scenarioNames[i]=="Optimistic_spr10" | CatchMSY_scenarioNames[i]=="Optimistic_spr20" | CatchMSY_scenarioNames[i]=="Optimistic_spr30")
    {
      Catch = SpeciesCatch_mT_Scaled[SpeciesCatch_mT_Scaled$Year==max(SpeciesCatch_mT_Scaled$Year),]
      Catch$Year=NULL
      Years = data.frame(Year=seq((max(SpeciesCatch_mT_Scaled$Year)-cMSY_optimisticScenario_yrsToCurrentDepletion),max(SpeciesCatch_mT_Scaled$Year),by=1))
      Catch = merge(Catch,Years)
    }
    if(CatchMSY_scenarioNames[i]=="Pessimistic" | CatchMSY_scenarioNames[i]=="Pessimistic_spr10" | CatchMSY_scenarioNames[i]=="Pessimistic_spr20" | CatchMSY_scenarioNames[i]=="Pessimistic_spr30")
    {
      Catch = SpeciesCatch_mT_Scaled[SpeciesCatch_mT_Scaled$Year==max(SpeciesCatch_mT_Scaled$Year),]
      Catch$Year=NULL
      Years = data.frame(Year=seq((max(SpeciesCatch_mT_Scaled$Year)-cMSY_pessimisticScenario_yrsToCurrentDepletion),max(SpeciesCatch_mT_Scaled$Year),by=1))
      Catch = merge(Catch,Years)
    }
    if(CatchMSY_scenarioNames[i]=="Constant" | CatchMSY_scenarioNames[i]=="Constant_spr10" | CatchMSY_scenarioNames[i]=="Constant_spr20" | CatchMSY_scenarioNames[i]=="Constant_spr30")
    {
      Catch = SpeciesCatch_mT_Scaled[SpeciesCatch_mT_Scaled$Year==max(SpeciesCatch_mT_Scaled$Year),]
      Catch$Year=NULL
      Years = data.frame(Year=seq((max(SpeciesCatch_mT_Scaled$Year)-cMSY_pessimisticScenario_yrsToCurrentDepletion),max(SpeciesCatch_mT_Scaled$Year),by=1))
      Catch = merge(Catch,Years)
    }
    for(G in 1:length(LandingsUncertaintyScenarios))
    {
      Catch[as.character(WPP)] = Catch[as.character(WPP)]*LandingsUncertaintyScenarios[G]   #Scale catch to various sensitivity levels
      yr = sort(unique(Catch$Year))  #years in time series
      nyr = length(yr)    #number of years in the time series
      for(j in 1:length(SpeciesForAnalysis_codes))
      {
         Catch_j = Catch[Catch$fish_id==SpeciesForAnalysis_codes[j],]
         for(p in 1:(length(WPP)+1))
         {
            #Don't process those where the number of length samples are too few (less than 150 fish)
            if(p!=(length(WPP)+1))
            {
              Sub = Lengths_cpue[Lengths_cpue$fish_id==SpeciesForAnalysis_codes[j] & Lengths_cpue$wpp1==WPP[p],]
            }
            if(p==(length(WPP)+1))
            {
              Sub = Lengths_cpue[Lengths_cpue$fish_id==SpeciesForAnalysis_codes[j],]
            }
            if(dim(Sub)[1]<150)     #Only run if sample size greater than 150 samples
            {
              next
            }
            if(p!=(length(WPP)+1))
            {
              currentSPRvalue = (StartValues[StartValues$WPP==as.character(WPP[p]) & StartValues$Genus==GenusForAnalysis_names[j] & StartValues$Species==SpeciesForAnalysis_names[j] & StartValues$Scenario==CatchMSY_scenarioNames[i],])$spr_current
            }
            if(p==(length(WPP)+1))
            {
              currentSPRvalue = (StartValues[StartValues$WPP=="All" & StartValues$Genus==GenusForAnalysis_names[j] & StartValues$Species==SpeciesForAnalysis_names[j] & StartValues$Scenario==CatchMSY_scenarioNames[i],])$spr_current
            }
            if(length(currentSPRvalue)==0)
            {
              next
            }
            if(currentSPRvalue<0.1 | (currentSPRvalue>=0.1 & (CatchMSY_scenarioNames[i]!="Historical_spr10" & CatchMSY_scenarioNames[i]!="Optimistic_spr10" & CatchMSY_scenarioNames[i]!="Pessimistic_spr10")))
            {
              if(p!=(length(WPP)+1))
              {
                ct = Catch_j[,as.character(WPP[p])] #assumes that catch is given in tonnes

                Params_j = StartValues[StartValues$Scenario==CatchMSY_scenarioNames[i] &
                           StartValues$Genus==GenusForAnalysis_names[j] & StartValues$Species==SpeciesForAnalysis_names[j] & StartValues$WPP==as.character(WPP[p]),]
              }
              if(p==(length(WPP)+1))
              {
                 CatchAllWpp = subset(Catch_j,select=-c(fish_id,FamilyCommonName,Year,Genus,Species))
                 ct = rowSums(CatchAllWpp)
                Params_j = StartValues[StartValues$Scenario==CatchMSY_scenarioNames[i] &
                           StartValues$Genus==GenusForAnalysis_names[j] & StartValues$Species==SpeciesForAnalysis_names[j] & StartValues$WPP=="All",]
              }
              if(sum(ct)==0)
              {
                sigR=0.05
                output = data.frame(CatchMSY_scenarioNames[i],LandingsUncertaintyScenarios[G],GenusForAnalysis_names[j],SpeciesForAnalysis_names[j],
                    as.character(WPP[p]),sigR,min(yr),max(yr),max(ct),ct[1],
                    ct[nyr],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"No Catch")
                 if(fileLineWriteCounter==1)
                 {
                    colNames = c("Scenario","LandingsSensitivity","Genus","Species","WPP","sigR","min_yr","max_yr","max_catch","catch_t0","catch_t","NumSims","BiomassMean_Year1",
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
                  output = data.frame(CatchMSY_scenarioNames[i],LandingsUncertaintyScenarios[G],GenusForAnalysis_names[j],SpeciesForAnalysis_names[j],
                    as.character(WPP[p]),sigR,min(yr),max(yr),max(ct),ct[1],
                    ct[nyr],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,paste("Too few (", length(r1), ") possible r-k combinations, check input parameters",sep=""))
                  if(fileLineWriteCounter==1)
                  {
                    colNames = c("Scenario","LandingsSensitivity","Genus","Species","WPP","sigR","min_yr","max_yr","max_catch","catch_t0","catch_t","NumSims","BiomassMean_Year1",
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
                  write.table(bt,paste(PATH_catchMSY_diagnostics,paste(paste("BiomassEstimates_RAW",GenusForAnalysis_names[j],
                      SpeciesForAnalysis_names[j],CatchMSY_scenarioNames[i],as.character(WPP[p]),paste("LandingSens",LandingsUncertaintyScenarios[G],sep=""),sep="_"),".csv",sep=""),sep="/"),
                      row.names=FALSE,col.names=TRUE,sep=",")
                  write.table(BiomassEstimates,paste(PATH_catchMSY_diagnostics,paste(paste("BiomassEstimates",GenusForAnalysis_names[j],
                      SpeciesForAnalysis_names[j],CatchMSY_scenarioNames[i],as.character(WPP[p]),paste("LandingSens",LandingsUncertaintyScenarios[G],sep=""),sep="_"),".csv",sep=""),sep="/"),
                      row.names=FALSE,col.names=TRUE,sep=",")
                  ## plot MSY over catch data
                  if(plotCatchMSY_diagnostics==TRUE)
                  {
                    FolderNameSpp = paste(as.character(GenusForAnalysis_names[i]),CapStr(as.character(SpeciesForAnalysis_names[i])),sep="")
                    png(paste(paste(PATH_catchMSY_diagnostics,FolderNameSpp,sep="/"),paste(paste("Cmsy",GenusForAnalysis_names[j],SpeciesForAnalysis_names[j],CatchMSY_scenarioNames[i],as.character(WPP[p]),paste("LandingSens",LandingsUncertaintyScenarios[G],sep=""),sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
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
                    title(paste(GenusForAnalysis_names[j],SpeciesForAnalysis_names[j],sep=" "))
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
                    hist(msy, freq=F, xlim=c(0, 1.2 * max(msy)), xlab="MSY (mT)",main = paste(GenusForAnalysis_names[j],SpeciesForAnalysis_names[j],sep=" "),cex.axis=1.1,cex.lab=1.1)
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
                    output = data.frame(CatchMSY_scenarioNames[i],LandingsUncertaintyScenarios[G],GenusForAnalysis_names[j],SpeciesForAnalysis_names[j],
                      as.character(WPP[p]),sigR,min(yr),max(yr),max(ct),ct[1],
                      ct[nyr],length(r),BiomassEstimates$Mean[1],BiomassEstimates$Upper[1],BiomassEstimates$Lower[1],
                      BiomassEstimates$Mean[dim(BiomassEstimates)[1]],BiomassEstimates$Upper[dim(BiomassEstimates)[1]],BiomassEstimates$Lower[dim(BiomassEstimates)[1]],
                      r_mean,r_sd,r_upper,r_lower,k_mean,k_sd,k_upper,k_lower,msy_mean,msy_sd,msy_upper,msy_lower,"Normal Run")
                  }
                  if(p==(length(WPP)+1))
                  {
                    output = data.frame(CatchMSY_scenarioNames[i],LandingsUncertaintyScenarios[G],GenusForAnalysis_names[j],SpeciesForAnalysis_names[j],
                      "EEZ",sigR,min(yr),max(yr),max(ct),ct[1],
                      ct[nyr],length(r),BiomassEstimates$Mean[1],BiomassEstimates$Upper[1],BiomassEstimates$Lower[1],
                      BiomassEstimates$Mean[dim(BiomassEstimates)[1]],BiomassEstimates$Upper[dim(BiomassEstimates)[1]],BiomassEstimates$Lower[dim(BiomassEstimates)[1]],
                      r_mean,r_sd,r_upper,r_lower,k_mean,k_sd,k_upper,k_lower,msy_mean,msy_sd,msy_upper,msy_lower,"Normal Run")
                  }
                  if(fileLineWriteCounter==1)
                  {
                      colNames = c("Scenario","LandingsSensitivity","Genus","Species","WPP","sigR","min_yr","max_yr","max_catch","catch_t0","catch_t","NumSims","BiomassMean_Year1",
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


#Need to modify the code below
CmsyResults = read.table(paste(PATH_output,CatchMSY_outfileName,sep="/"),header=TRUE,sep=",")
BiomassCarryingCapacityMSY_mean = subset(CmsyResults,select=c(Scenario,Genus,Species,WPP,catch_t,BiomassMean_LastYear,r_mean,k_mean,msy_mean,Notes))
BiomassCarryingCapacityMSY_mean$Bound="mean"
BiomassCarryingCapacityMSY_lower = subset(CmsyResults,select=c(Scenario,Genus,Species,WPP,catch_t,BiomassLower_LastYear,r_lower,k_lower,msy_lower,Notes))
BiomassCarryingCapacityMSY_lower$Bound="lower"
BiomassCarryingCapacityMSY_upper = subset(CmsyResults,select=c(Scenario,Genus,Species,WPP,catch_t,BiomassUpper_LastYear,r_upper,k_upper,msy_upper,Notes))
BiomassCarryingCapacityMSY_upper$Bound="upper"
names(BiomassCarryingCapacityMSY_lower) = names(BiomassCarryingCapacityMSY_mean)
names(BiomassCarryingCapacityMSY_upper) = names(BiomassCarryingCapacityMSY_mean)
BiomassCarryingCapacityMSY = rbind(BiomassCarryingCapacityMSY_mean,BiomassCarryingCapacityMSY_lower,BiomassCarryingCapacityMSY_upper)
GenusSpeciesFishID = subset(Species_cpue,select=c(Genus,Species,fish_id))
BiomassCarryingCapacityMSY = merge(BiomassCarryingCapacityMSY,GenusSpeciesFishID,all.x=TRUE,by=c("Genus","Species"))
BiomassCarryingCapacityMSY$Genus = as.character(BiomassCarryingCapacityMSY$Genus)
BiomassCarryingCapacityMSY$Species = as.character(BiomassCarryingCapacityMSY$Species)
BiomassCarryingCapacityMSY$Scenario = as.character(BiomassCarryingCapacityMSY$Scenario)
BiomassCarryingCapacityMSY$WPP = as.character(BiomassCarryingCapacityMSY$WPP)
BiomassCarryingCapacityMSY$Notes = as.character(BiomassCarryingCapacityMSY$Notes)
BiomassCarryingCapacityMSY = BiomassCarryingCapacityMSY[order(BiomassCarryingCapacityMSY$Genus,BiomassCarryingCapacityMSY$Species,BiomassCarryingCapacityMSY$WPP,BiomassCarryingCapacityMSY$Scenario,BiomassCarryingCapacityMSY$Bound),]
BiomassCarryingCapacityMSY = BiomassCarryingCapacityMSY[BiomassCarryingCapacityMSY$catch_t>0,]
BiomassCarryingCapacityMSY = BiomassCarryingCapacityMSY[BiomassCarryingCapacityMSY$BiomassMean_LastYear>0,]
BiomassCarryingCapacityMSY$Scenario_New = paste(BiomassCarryingCapacityMSY$Scenario,BiomassCarryingCapacityMSY$Bound,sep="-")
write.table(BiomassCarryingCapacityMSY,paste(PATH_output,"BiomassCarryingCapacityMSY.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")

################################ Make Plots of Catch MSY Results #############################################
if(plotCatchMSY_results==TRUE)
{
  BiomassCarryingCapacityMSY = read.table(paste(PATH_output,"BiomassCarryingCapacityMSY.csv",sep="/"),header=TRUE,sep=",")
  SpeciesAreas_ForPlots = subset(BiomassCarryingCapacityMSY,select=c(Genus,Species,WPP))
  SpeciesAreas_ForPlots = unique(SpeciesAreas_ForPlots)
  for(i in 1:dim(SpeciesAreas_ForPlots)[1])
  {
     BiomassCarryingCapacityMSY_i_j = BiomassCarryingCapacityMSY[BiomassCarryingCapacityMSY$Genus==SpeciesAreas_ForPlots$Genus[i] & BiomassCarryingCapacityMSY$Species==SpeciesAreas_ForPlots$Species[i] & BiomassCarryingCapacityMSY$WPP==SpeciesAreas_ForPlots$WPP[i],]
     BiomassCarryingCapacityMSY_i_j = BiomassCarryingCapacityMSY_i_j[order(BiomassCarryingCapacityMSY_i_j$Scenario_New),]
     FolderNameSpp = paste(as.character(SpeciesAreas_ForPlots$Genus[i]),CapStr(as.character(SpeciesAreas_ForPlots$Species[i])),sep="")
     png(paste(paste(PATH_catchMSY_plots,FolderNameSpp,sep="/"),paste(paste("CatchMSY_Results",SpeciesAreas_ForPlots$Genus[i],SpeciesAreas_ForPlots$Species[i],SpeciesAreas_ForPlots$WPP[i],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=4200,res=600)
     par(mfrow=c(2,2),mar=c(9, 4, 4, 2) + 0.1) #HERE1 HERE1
     plot(as.factor(c(1:dim(BiomassCarryingCapacityMSY_i_j)[1])),BiomassCarryingCapacityMSY_i_j$BiomassMean_LastYear,ylab="Current Biomass (mT)",main=paste(GenusForAnalysis_names[i],SpeciesForAnalysis_names[i],sep=" "),xaxt="n")
     axis(1,at=1:dim(BiomassCarryingCapacityMSY_i_j)[1],labels=BiomassCarryingCapacityMSY_i_j$Scenario_New,las=2,cex.axis=0.8)
     mtext(SpeciesAreas_ForPlots$WPP[i])
     plot(as.factor(c(1:dim(BiomassCarryingCapacityMSY_i_j)[1])),BiomassCarryingCapacityMSY_i_j$r_mean,ylab="r (population growth rate)",main=paste(GenusForAnalysis_names[i],SpeciesForAnalysis_names[i],sep=" "),xaxt="n")
     axis(1,at=1:dim(BiomassCarryingCapacityMSY_i_j)[1],labels=BiomassCarryingCapacityMSY_i_j$Scenario_New,las=2,cex.axis=0.8)
     mtext(SpeciesAreas_ForPlots$WPP[i])
     plot(as.factor(c(1:dim(BiomassCarryingCapacityMSY_i_j)[1])),BiomassCarryingCapacityMSY_i_j$k_mean,ylab="Carrying Capacity (mT)",main=paste(GenusForAnalysis_names[i],SpeciesForAnalysis_names[i],sep=" "),xaxt="n")
     axis(1,at=1:dim(BiomassCarryingCapacityMSY_i_j)[1],labels=BiomassCarryingCapacityMSY_i_j$Scenario_New,las=2,cex.axis=0.8)
     mtext(SpeciesAreas_ForPlots$WPP[i])
     plot(as.factor(c(1:dim(BiomassCarryingCapacityMSY_i_j)[1])),BiomassCarryingCapacityMSY_i_j$msy_mean,ylab="MSY (mT)",main=paste(GenusForAnalysis_names[i],SpeciesForAnalysis_names[i],sep=" "),xaxt="n")
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

################################ Compute Natural Mortaltiy and Linf #########################################
# Calculate Von Bertalanffy growth paramenter K using L infinity and L opt (Froese and Binohlan 2000) Assuming an average water temp of 20C, see Froese & Pauly 2000 in Martinez-Andrade 2003
Species_cpue$Mnat_current = 10^(0.566-(0.718*log10(Species_cpue$linf))+0.02*20)   #Where Linf is 90% of Lmax
Species_cpue$Mnat_Lmax = 10^(0.566-(0.718*log10(Species_cpue$lmax))+0.02*20)   #Where Linf = Lmax
Species_cpue$VBG_K_current = (Species_cpue$Mnat_current*Species_cpue$lopt)/(3*(Species_cpue$linf-Species_cpue$lopt))
Species_cpue$VBG_K_Lmax = (Species_cpue$Mnat_Lmax*Species_cpue$lopt)/(3*(Species_cpue$lmax-Species_cpue$lopt))
Species_cpue = Species_cpue[order(Species_cpue$Genus,Species_cpue$Species),]
write.table(Species_cpue,paste(PATH_output,"SpeciesLifeHistoryParameters.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)

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
ComputeFecundity = function(Weight_at_size_vector,Weight_at_size_vector_min,Weight_at_size_vector_max,FecundityMin,FecundityMax)
{
  Fec = data.frame(Weight=c(Weight_at_size_vector_min,Weight_at_size_vector_max),Fecundity=c(FecundityMin,FecundityMax))
  LinMod = lm(Fec$Fecundity ~ Fec$Weight)
  Fecundity = (coef(LinMod)[[1]] + Weight_at_size_vector*coef(LinMod)[[2]])
  Fecundity[Fecundity<0]=0
  return(Fecundity)
}

######################### Function: Find R0 ####################################
findR0 = function(M,Linf,VBK,a,b,K,maxAge,matureAge,init_min_GT,init_max_GT,incrementBy)
{
  age = c(0:maxAge)
  data = data.frame(age=age)
  data$maturity = 1
  data$maturity[1:matureAge] = 0 #maturity is knife-edge. only age mature=3 and above are matured
  natmort = exp(-M) # this is per year, we can make it per month
  data$length = (1-exp(-data$age*VBK))*Linf #this is the age-length relationship   - length is in cm
  data$weight = (a*(data$length^b))*0.001 #this is the length-weight relationship - weight is in kilograms
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

############ Compute Sample Sizes of Lengths by Species and WPP - Do NOT Move Forward with WPPs that have too Few Samples #############

TableSampleSize = as.data.frame.matrix(table(Lengths_cpue$fish_id,Lengths_cpue$wpp1))
TableSampleSize$Total = rowSums(TableSampleSize)
TableSampleSize$fish_id = row.names(TableSampleSize)
row.names(TableSampleSize)=NULL
TableSampleSize_EEZ = subset(TableSampleSize,select=c(fish_id,Total))
TableSampleSize_EEZ$WPP="EEZ"
names(TableSampleSize_EEZ)[names(TableSampleSize_EEZ)=="Total"]="Freq"
TableSampleSize_long = data.frame(table(Lengths_cpue$fish_id,Lengths_cpue$wpp1))
names(TableSampleSize_long) = c("fish_id","WPP","Freq")
TableSampleSize_EEZ = TableSampleSize_EEZ[names(TableSampleSize_long)]
TableSampleSize_long = rbind(TableSampleSize_long,TableSampleSize_EEZ)
TableSampleSize_long$Analyze=TRUE
TableSampleSize_long$Analyze[TableSampleSize_long$Freq<100]=FALSE
BiomassCarryingCapacityMSY = merge(BiomassCarryingCapacityMSY,TableSampleSize_long,all.x=TRUE,by=c("fish_id","WPP"))
BiomassCarryingCapacityMSY = BiomassCarryingCapacityMSY[BiomassCarryingCapacityMSY$Analyze==TRUE,]

########### Set up a List Object of Data.frames at Each Length Bin for Each Species, Scenario, and WPP #######
StartValues = read.table(paste(PATH_otherInput,CatchMSY_inputFileName,sep="/"),sep=",",header=TRUE,as.is=TRUE)
Lengths_cpue$cm_converted = Lengths_cpue$cm*Lengths_cpue$conversion_factor_tl2fl
Lengths_cpue$LengthBin_Lower = floor(Lengths_cpue$cm_converted/5)*5
Lengths_cpue$LengthBin_Upper = ceiling(Lengths_cpue$cm_converted/5)*5
Lengths_cpue$LengthBin_Midpoint = Lengths_cpue$LengthBin_Lower + 2.5
SpeciesForAnalysis_codes_OnlySpecies = SpeciesForAnalysis_codes[SpeciesForAnalysis_codes!=as.character(-999) & SpeciesForAnalysis_codes!=as.character(-9999)]

if(runLengthReconstruction==TRUE)
{
  ListOfPopulations=list()
  for(i in 1:dim(BiomassCarryingCapacityMSY)[1])   #753 items, 150 on each of 5 processors: 1-150, 151-300, 301-450, 451-600, 601-753    for(i in 601:dim(BiomassCarryingCapacityMSY)[1])
  {
    Species_i = Species_cpue[Species_cpue$fish_id==BiomassCarryingCapacityMSY$fish_id[i],]
    if(BiomassCarryingCapacityMSY$WPP[i]!="EEZ")
    {
      Sub = Lengths_cpue[Lengths_cpue$fish_id==BiomassCarryingCapacityMSY$fish_id[i] & Lengths_cpue$wpp1==BiomassCarryingCapacityMSY$WPP[i],]
      StartValues_i_j_o_p = StartValues[StartValues$Scenario==BiomassCarryingCapacityMSY$Scenario[i] & StartValues$Genus==BiomassCarryingCapacityMSY$Genus[i] & StartValues$Species==BiomassCarryingCapacityMSY$Species[i] & StartValues$WPP==BiomassCarryingCapacityMSY$WPP[i],]
    }
    if(BiomassCarryingCapacityMSY$WPP[i]=="EEZ")
    {
      Sub = Lengths_cpue[Lengths_cpue$fish_id==BiomassCarryingCapacityMSY$fish_id[i],]
      StartValues_i_j_o_p = StartValues[StartValues$Scenario==BiomassCarryingCapacityMSY$Scenario[i] & StartValues$Genus==BiomassCarryingCapacityMSY$Genus[i] & StartValues$Species==BiomassCarryingCapacityMSY$Species[i] & StartValues$WPP=="All",]
    }
    if(dim(Sub)[1]<150)     #Only run if sample size greater than 150 samples
    {
      next
    }
    for(p in 1:length(LifeHistoryScenarios))
    {
      if(LifeHistoryScenarios[p]=="Linf_90Perc_lmax")
      {
          Linf = Species_i$linf
          Mnat = Species_i$Mnat_current
          VBK = Species_i$VBG_K_current
      }
      if(LifeHistoryScenarios[p]=="Linf_equal_lmax")
      {
          Linf = Species_i$lmax
          Mnat = Species_i$Mnat_Lmax
          VBK = Species_i$VBG_K_Lmax
      }
      #Start things length-based
      Sub$Age = log(1-(Sub$cm_converted/Linf))/(-1*VBK)
      Sub = Sub[Sub$Age!=Inf,]
      Sub = Sub[!is.na(Sub$Age),]
      Sub$Age_Int = floor(Sub$Age)
      Sub$Age_LowerLengthBin = log(1-(Sub$LengthBin_Lower/Linf))/(-1*VBK)
      Sub$Age_LowerLengthBin_Int = floor(Sub$Age_LowerLengthBin)
      Sub$Age_UpperLengthBin = log(1-(Sub$LengthBin_Upper/Linf))/(-1*VBK)
      Sub$Age_UpperLengthBin_Int = floor(Sub$Age_UpperLengthBin)
      NumAtLength_i_j = data.frame(table(Sub$LengthBin_Midpoint))
      names(NumAtLength_i_j) = c("Midpoint","NumSampled")
      NumAtLength_i_j$Midpoint = as.numeric(as.character(NumAtLength_i_j$Midpoint))
      Midpoints = data.frame(Midpoint=seq(2.5,((ceiling(Species_i$lmax/5)*5)-2.5),by=5))
      NumAtLength_i_j = merge(Midpoints,NumAtLength_i_j,all.x=TRUE,by=c("Midpoint"))
      NumAtLength_i_j[is.na(NumAtLength_i_j)]=0
      NumAtLength_i_j$Length_min = NumAtLength_i_j$Midpoint-2.5
      NumAtLength_i_j$Length_max = NumAtLength_i_j$Midpoint+2.5
      NumAtLength_i_j$Age_min = log(1-(NumAtLength_i_j$Length_min/Linf))/(-1*VBK)
      NumAtLength_i_j$Age_max = log(1-(NumAtLength_i_j$Length_max/Linf))/(-1*VBK)
      NumAtLength_i_j$Weight_kg_min = (Species_i$var_a*(NumAtLength_i_j$Length_min^Species_i$var_b))/1000
      NumAtLength_i_j$Weight_kg_max = (Species_i$var_a*(NumAtLength_i_j$Length_max^Species_i$var_b))/1000
      NumAtLength_i_j$Weight_kg = (Species_i$var_a*(NumAtLength_i_j$Midpoint^Species_i$var_b))/1000
      NumAtLength_i_j$AgeDifference_Years = NumAtLength_i_j$Age_max - NumAtLength_i_j$Age_min
      NumAtLength_i_j$AgeDifference_Years[NumAtLength_i_j$AgeDifference_Years==Inf]=NA
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
      NumAtAge_i_j$Length = Linf*(1-exp(-VBK*(NumAtAge_i_j$Age)))
      NumAtAge_i_j$Weight_kg = (Species_i$var_a*(NumAtAge_i_j$Length^Species_i$var_b))/1000
      NumAtAge_i_j$Maturity=0
      NumAtAge_i_j$Maturity[NumAtAge_i_j$Length>Species_i$lmat]=1
      NumAtAge_i_j$Fecundity = ComputeFecundity(NumAtAge_i_j$Weight_kg,min(NumAtAge_i_j$Weight_kg[NumAtAge_i_j$Maturity==1]),max(NumAtAge_i_j$Weight_kg[NumAtAge_i_j$Maturity==1]),StartValues_i_j_o_p$Fecundity_min,StartValues_i_j_o_p$Fecundity_max)
      NumAtAge_i_j$Survival_M = exp(-Mnat*(NumAtAge_i_j$Age))
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
      NumAtAge_i_j$SSB_virginAtAge = NumAtAge_i_j$PopulationVirgin*(NumAtAge_i_j$Weight_kg*0.001)*NumAtAge_i_j$Maturity*NumAtAge_i_j$Fecundity
      SSB_virgin = sum(NumAtAge_i_j$SSB_virginAtAge)
      AgeMat = round(log(1-(Species_i$lmat/Linf))/(-1*VBK))
      print("Finding R0...")
      flush.console()
      init_max_GT=((BiomassCarryingCapacityMSY$k_mean[i]*1000)/Weight_kg_oneFish)
      digits = nchar(round(init_max_GT,0))
      incrementBy = 10^(digits-3)
      if(digits<=3) {incrementBy=1}
      R0 = findR0(M=Mnat,Linf=Linf,VBK=VBK,a=Species_i$var_a,b=Species_i$var_b,K=BiomassCarryingCapacityMSY$k_mean[i],maxAge=ceiling(max(NumAtAge_i_j$Age)),matureAge=AgeMat,init_min_GT=1,init_max_GT=init_max_GT,incrementBy=incrementBy)          #lower and upper starting values for seach are 10 and 50  percent of the total population number as extrapolated from catch MSY
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
            ALK_Len2Age[as.character(NumAtLength_i_j$Midpoint[u]),as.character(IntAge_min)] = 1
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
            ALK_Age2Len[as.character(NumAtLength_i_j$Midpoint[u]),as.character(IntAge_min)] = NumAtLength_i_j$AgeDifference_Years[u]
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
      FolderNameSpp = paste(as.character(BiomassCarryingCapacityMSY$Genus[i]),CapStr(as.character(BiomassCarryingCapacityMSY$Species[i])),sep="")
      png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("SelectivityAtAge",BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
      plot(NumAtAge_i_j$Age,NumAtAge_i_j$Selex_Empirical,pch=16,xlab="Age",ylab="Selectivity",cex.axis=1.1,cex.lab=1.1)
      points(NumAtAge_i_j$Age,NumAtAge_i_j$Selex_Fitted,type="l",lwd=2,col="red")
      title(paste(BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],sep=" "))
      mtext(paste(as.character(BiomassCarryingCapacityMSY$WPP[i]),paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],sep=", "))
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
      FolderNameSpp = paste(as.character(BiomassCarryingCapacityMSY$Genus[i]),CapStr(as.character(BiomassCarryingCapacityMSY$Species[i])),sep="")
      png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("SelectivityAtLength",BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
      plot(NumAtLength_i_j$Midpoint,NumAtLength_i_j$Selex_Empirical,pch=16,xlab="Length (cm)",ylab="Selectivity",cex.axis=1.1,cex.lab=1.1)
      points(NumAtLength_i_j$Midpoint,NumAtLength_i_j$Selex_Fitted,type="l",lwd=2,col="red")
      title(paste(BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],sep=" "))
      mtext(paste(as.character(BiomassCarryingCapacityMSY$WPP[i]),paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],sep=", "))
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
          MyPars@Species = Species_i$Species
          MyPars@Linf = Linf
          MyPars@L50 = Species_i$lmat #length at 50% maturity
          MyPars@L95 = Species_i$lopt #length at 95% maturity - we don't know this so use lopt as a proxy. We don't know this because we are assuming knife edge maturity and don't have a study to know what the real maturity schedule looks like
          MyPars@MK = Mnat/VBK
          MyPars@L_units = "cm"
          MyPars@BinMin=min(NumAtLength_i_j$Midpoint)-2.5
          MyPars@BinWidth=NumAtLength_i_j$Midpoint[2]-NumAtLength_i_j$Midpoint[1]
          MyPars@BinMax=ceiling(Species_i$lmax/5)*5
          Len1 = new("LB_lengths",LB_pars=MyPars,file=Sub$cm_converted,sataType="raw")
          myFit1 = LBSPRfit(MyPars, Len1)
          myFit1@Ests
          plotSize(myFit1)
          plotMat(myFit1)
          #simulate population
          SimPars=MyPars
          SimPars@SL50 = Species_i$lmat #length at 50% maturity
          SimPars@SL95 = Species_i$lopt #length at 95% maturity - we don't know this so use lopt as a proxy. We don't know this because we are assuming knife edge maturity and don't have a study to know what the real maturity schedule looks like
          SimPars@SPR = StartValues_i_j_o_p$spr_current   #Use Peter's SPR estimate, NOT LBSPR computed value
          SimPars@R0=R0
          MySim = LBSPRsim(SimPars,Control = list("modtype"="absel"))
          FolderNameSpp = paste(as.character(BiomassCarryingCapacityMSY$Genus[i]),CapStr(as.character(BiomassCarryingCapacityMSY$Species[i])),sep="")
          png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("LBSPR_diagnostics",BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
          plotSim(MySim)
          plotSim(MySim, lf.type="pop",perRec = TRUE)
          dev.off()
          F_estimated = MySim@FM*Mnat
          NumAtLength_i_j$Survival_LengthBased_F_M_Rescaled = data.frame(MySim@pLPop)$PopF
          EstNumFish_Population = (BiomassCarryingCapacityMSY$BiomassMean_LastYear[i]*1000)/Weight_kg_oneFish
          NumAtLength_i_j$Population_Extrapolated = data.frame(MySim@pLPop)$PopF*EstNumFish_Population
          NumAtLength_i_j$EstimatedCatch_Baranov = MySim@pLCatch*CatchNumber
          NumAtLength_i_j$Residuals = NumAtLength_i_j$ExtrapolatedCatch-NumAtLength_i_j$EstimatedCatch_Baranov
          FolderNameSpp = paste(as.character(BiomassCarryingCapacityMSY$Genus[i]),CapStr(as.character(BiomassCarryingCapacityMSY$Species[i])),sep="")
          png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("CatchAtLength",BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
          plot(NumAtLength_i_j$Midpoint,NumAtLength_i_j$ExtrapolatedCatch,pch=16,xlab="Length (cm)",ylab="Catch in Number",cex.axis=1.1,cex.lab=1.1,ylim=c(0,max(NumAtLength_i_j$EstimatedCatch_Baranov,NumAtLength_i_j$ExtrapolatedCatch)))
          points(NumAtLength_i_j$Midpoint,NumAtLength_i_j$EstimatedCatch_Baranov,type="l",lwd=2,col="red")
          title(paste(BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],sep=" "))
          mtext(paste(BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep=", "))
          dev.off()
          PopAtLen = do.call("cbind",replicate(dim(ALK_Len2Age)[2],data.frame(NumAtLength_i_j$Population_Extrapolated),simplify = FALSE))
          PopAtLen = PopAtLen*ALK_Len2Age
          NumAtAge_i_j$Population_Extrapolated = colSums(PopAtLen)
          NumAtAge_i_j$SSB_CurrentAtAge = NumAtAge_i_j$Population_Extrapolated*(NumAtAge_i_j$Weight_kg*0.001)*NumAtAge_i_j$Maturity*NumAtAge_i_j$Fecundity
          CatchAtLen = do.call("cbind",replicate(dim(ALK_Len2Age)[2],data.frame(NumAtLength_i_j$EstimatedCatch_Baranov),simplify = FALSE))
          CatchAtLen = CatchAtLen*ALK_Len2Age
          NumAtAge_i_j$EstimatedCatch_Baranov = colSums(CatchAtLen)
          FolderNameSpp = paste(as.character(BiomassCarryingCapacityMSY$Genus[i]),CapStr(as.character(BiomassCarryingCapacityMSY$Species[i])),sep="")
          png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("CatchAtAge",BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
          plot(NumAtAge_i_j$Age,NumAtAge_i_j$ExtrapolatedCatchInNumber,pch=16,xlab="Age",ylab="Catch in Number",cex.axis=1.1,cex.lab=1.1,ylim=c(0,max(NumAtAge_i_j$EstimatedCatch_Baranov,NumAtAge_i_j$ExtrapolatedCatchInNumber)))
          points(NumAtAge_i_j$Age,NumAtAge_i_j$EstimatedCatch_Baranov,type="l",lwd=2,col="red")
          title(paste(BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],sep=" "))
          mtext(paste(BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep=", "))
          dev.off()
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
            NumAtAge_i_j$Survival_F_M = exp(-1*(Mnat+(F_start*NumAtAge_i_j$Selex_Empirical))*NumAtAge_i_j$Age)
            NumAtAge_i_j$Survival_F_M_Rescaled = NumAtAge_i_j$Survival_F_M/sum(NumAtAge_i_j$Survival_F_M)
            EstNumFish_Population = (BiomassCarryingCapacityMSY$BiomassMean_LastYear[i]*1000)/Weight_kg_oneFish
            NumAtAge_i_j$Population_Extrapolated = NumAtAge_i_j$Survival_F_M_Rescaled*EstNumFish_Population
            NumAtAge_i_j$SSB_CurrentAtAge = NumAtAge_i_j$Population_Extrapolated*(NumAtAge_i_j$Weight_kg*0.001)*NumAtAge_i_j$Maturity*NumAtAge_i_j$Fecundity
            SSB_current = sum(NumAtAge_i_j$SSB_CurrentAtAge)
            SPR_current = SSB_current/SSB_virgin
            return(sum((SPR_current - StartValues_i_j_o_p$spr_current)^2))
          }
          Model = optimize(min.RSS,c(0,10),tol=0.001)
          NumAtAge_i_j$Survival_F_M = exp(-1*(Mnat+(Model$minimum[1]*NumAtAge_i_j$Selex_Empirical))*NumAtAge_i_j$Age)     #"pseudosurvival"
          NumAtAge_i_j$Survival_F_M_Rescaled = NumAtAge_i_j$Survival_F_M/sum(NumAtAge_i_j$Survival_F_M)
          EstNumFish_Population = (BiomassCarryingCapacityMSY$BiomassMean_LastYear[i]*1000)/Weight_kg_oneFish
          NumAtAge_i_j$Population_Extrapolated = NumAtAge_i_j$Survival_F_M_Rescaled*EstNumFish_Population
          NumAtAge_i_j$Survival_F_M=NULL
          NumAtAge_i_j$Survival_F_M_Rescaled=NULL
          F_start = 1.0
          min.RSS = function(F_start)
          {
            EstimatedCatch_Baranov = BaranovCatch(NumAtAge_i_j$Age,NumAtAge_i_j$Length,NumAtAge_i_j$Population_Extrapolated,Mnat,F_start,NumAtAge_i_j$Selex_Empirical)
            return(sum((EstimatedCatch_Baranov - NumAtAge_i_j$ExtrapolatedCatchInNumber)^2))
          }
          Model = optimize(min.RSS,c(0,10),tol=0.001)
          F_estimated = Model$minimum[1]
          NumAtAge_i_j$EstimatedCatch_Baranov = BaranovCatch(NumAtAge_i_j$Age,NumAtAge_i_j$Length,NumAtAge_i_j$Population_Extrapolated,Mnat,F_estimated,NumAtAge_i_j$Selex_Empirical)
          NumAtAge_i_j$Residuals = NumAtAge_i_j$ExtrapolatedCatchInNumber-NumAtAge_i_j$EstimatedCatch_Baranov
          NumAtAge_i_j$SSB_CurrentAtAge = NumAtAge_i_j$Population_Extrapolated*(NumAtAge_i_j$Weight_kg*0.001)*NumAtAge_i_j$Maturity*NumAtAge_i_j$Fecundity
          SSB_current = sum(NumAtAge_i_j$SSB_CurrentAtAge)
          SPR_current = SSB_current/SSB_virgin
          NumAgeCurrent = do.call("rbind",replicate(dim(ALK_Age2Len)[1],transpose(data.frame(NumAtAge_i_j$Population_Extrapolated)),simplify = FALSE))
          NumLengthCurrent = NumAgeCurrent*ALK_Age2Len
          NumAtLength_i_j$Population_Extrapolated = rowSums(NumLengthCurrent)
          CatchCurrent = do.call("rbind",replicate(dim(ALK_Age2Len)[1],transpose(data.frame(NumAtAge_i_j$EstimatedCatch_Baranov)),simplify = FALSE))
          CatchCurrent = CatchCurrent*ALK_Age2Len
          NumAtLength_i_j$EstimatedCatch_Baranov = rowSums(CatchCurrent)
          FolderNameSpp = paste(as.character(BiomassCarryingCapacityMSY$Genus[i]),CapStr(as.character(BiomassCarryingCapacityMSY$Species[i])),sep="")
          png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("CatchAtLength",BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
          plot(NumAtLength_i_j$Midpoint,NumAtLength_i_j$ExtrapolatedCatch,pch=16,xlab="Length (cm)",ylab="Catch in Number",cex.axis=1.1,cex.lab=1.1,ylim=c(0,max(NumAtLength_i_j$EstimatedCatch_Baranov,NumAtLength_i_j$ExtrapolatedCatch)))
          points(NumAtLength_i_j$Midpoint,NumAtLength_i_j$EstimatedCatch_Baranov,type="l",lwd=2,col="red")
          title(paste(BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],sep=" "))
          mtext(paste(BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep=", "))
          dev.off()
          png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("CatchAtAge",BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
          plot(NumAtAge_i_j$Age,NumAtAge_i_j$ExtrapolatedCatchInNumber,pch=16,xlab="Age",ylab="Catch in Number",cex.axis=1.1,cex.lab=1.1,ylim=c(0,max(NumAtAge_i_j$EstimatedCatch_Baranov,NumAtAge_i_j$ExtrapolatedCatchInNumber)))
          points(NumAtAge_i_j$Age,NumAtAge_i_j$EstimatedCatch_Baranov,type="l",lwd=2,col="red")
          title(paste(BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],sep=" "))
          mtext(paste(BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep=", "))
          dev.off()
        }
        if(populationReconstructionMethods[v]=="BevertonHoltInstantaneousMortality")
        {
          NumAtLength_i_j = NumAtLength_v
          NumAtAge_i_j = NumAtAge_v
          print("Reconstructing population size using Baranov-Holt Instantaneous Mortality...")
          flush.console()
          F_start = 1.0
          min.RSS = function(F_start)
          {
            LenVec = Sub$cm_converted
            LenVec = LenVec[!is.na(LenVec)]
            SPR = computeLengthSPR(Mnat,Linf,VBK,Species_i$var_a,Species_i$var_b,Species_i$lmat,F_start,LenVec)
            return(sum((SPR - StartValues_i_j_o_p$spr_current)^2))
          }
          Model = optimize(min.RSS,c(0,10),tol=0.001)
          F_estimated = Model$minimum[1]
          NumAtAge_i_j$Survival_F_M = exp(-1*(Mnat+(F_estimated*NumAtAge_i_j$Selex_Empirical))*NumAtAge_i_j$Age)
          NumAtAge_i_j$Survival_F_M_Rescaled = NumAtAge_i_j$Survival_F_M/sum(NumAtAge_i_j$Survival_F_M)
          EstNumFish_Population = (BiomassCarryingCapacityMSY$BiomassMean_LastYear[i]*1000)/Weight_kg_oneFish
          NumAtAge_i_j$Population_Extrapolated = NumAtAge_i_j$Survival_F_M_Rescaled*EstNumFish_Population
          NumAtAge_i_j$SSB_CurrentAtAge = NumAtAge_i_j$Population_Extrapolated*(NumAtAge_i_j$Weight_kg*0.001)*NumAtAge_i_j$Maturity*NumAtAge_i_j$Fecundity
          SSB_current = sum(NumAtAge_i_j$SSB_CurrentAtAge)
          SPR_current = SSB_current/SSB_virgin
          NumAtAge_i_j$EstimatedCatch_Baranov = BaranovCatch(NumAtAge_i_j$Age,NumAtAge_i_j$Length,NumAtAge_i_j$Population_Extrapolated,Mnat,F_estimated,NumAtAge_i_j$Selex_Empirical)
          NumAtAge_i_j$Residuals = NumAtAge_i_j$ExtrapolatedCatchInNumber-NumAtAge_i_j$EstimatedCatch_Baranov
          NumAtAge_i_j$SSB_CurrentAtAge = NumAtAge_i_j$Population_Extrapolated*(NumAtAge_i_j$Weight_kg*0.001)*NumAtAge_i_j$Maturity*NumAtAge_i_j$Fecundity
          SSB_current = sum(NumAtAge_i_j$SSB_CurrentAtAge)
          SPR_current = SSB_current/SSB_virgin
          NumAgeCurrent = do.call("rbind",replicate(dim(ALK_Age2Len)[1],transpose(data.frame(NumAtAge_i_j$Population_Extrapolated)),simplify = FALSE))
          NumLengthCurrent = NumAgeCurrent*ALK_Age2Len
          NumAtLength_i_j$Population_Extrapolated = rowSums(NumLengthCurrent)
          CatchCurrent = do.call("rbind",replicate(dim(ALK_Age2Len)[1],transpose(data.frame(NumAtAge_i_j$EstimatedCatch_Baranov)),simplify = FALSE))
          CatchCurrent = CatchCurrent*ALK_Age2Len
          NumAtLength_i_j$EstimatedCatch_Baranov = rowSums(CatchCurrent)
          FolderNameSpp = paste(as.character(BiomassCarryingCapacityMSY$Genus[i]),CapStr(as.character(BiomassCarryingCapacityMSY$Species[i])),sep="")
          png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("CatchAtLength",BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
          plot(NumAtLength_i_j$Midpoint,NumAtLength_i_j$ExtrapolatedCatch,pch=16,xlab="Length (cm)",ylab="Catch in Number",cex.axis=1.1,cex.lab=1.1,ylim=c(0,max(NumAtLength_i_j$EstimatedCatch_Baranov,NumAtLength_i_j$ExtrapolatedCatch)))
          points(NumAtLength_i_j$Midpoint,NumAtLength_i_j$EstimatedCatch_Baranov,type="l",lwd=2,col="red")
          title(paste(BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],sep=" "))
          mtext(paste(BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep=", "))
          dev.off()
          png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("CatchAtAge",BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
          plot(NumAtAge_i_j$Age,NumAtAge_i_j$ExtrapolatedCatchInNumber,pch=16,xlab="Age",ylab="Catch in Number",cex.axis=1.1,cex.lab=1.1,ylim=c(0,max(NumAtAge_i_j$EstimatedCatch_Baranov,NumAtAge_i_j$ExtrapolatedCatchInNumber)))
          points(NumAtAge_i_j$Age,NumAtAge_i_j$EstimatedCatch_Baranov,type="l",lwd=2,col="red")
          title(paste(BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],sep=" "))
          mtext(paste(BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep=", "))
          dev.off()
        }
        if(populationReconstructionMethods[v]=="LIME" & BiomassCarryingCapacityMSY$Bound[i]=="mean")  #The condition that only the mean is used for LIME is becasue LIME estiamtes its own biomass levels (unlike the other 4 methods which use biomass from Catch-MSY). As a result, because we have upper and lower CatchMSY bounds as well, if we don't include this condition, we will accidentially run LIME three times and get the same results, whcih will go into the averaging and bias the averaged values.
        {
          NumAtLength_i_j = NumAtLength_v
          NumAtAge_i_j = NumAtAge_v
          print("Reconstructing population size using LIME...")
          flush.console()
          Len50_selex = 30
          min.RSS = function(Len50_selex)
          {
            Sel = SelectivityLogistic_Length(Len50_selex,Selex_Logistic_a_length,Selex_Logistic_b_length)
            return(sum((Sel - 0.5)^2))
          }
          Model = optimize(min.RSS,c(0,200),tol=0.001)
          Len50_selex = Model$minimum[1]
          Len95_selex = 50
          min.RSS = function(Len95_selex)
          {
            Sel = SelectivityLogistic_Length(Len95_selex,Selex_Logistic_a_length,Selex_Logistic_b_length)
            return(sum((Sel - 0.95)^2))
          }
          Model = optimize(min.RSS,c(0,200),tol=0.001)
          Len95_selex = Model$minimum[1]
          lh = create_lh_list(vbk=VBK,linf=Linf,lwa=Species_i$var_a,lwb=Species_i$var_b,
              S50=Len50_selex,S95=Len95_selex,selex_input="length",selex_type=c("logistic"), #Starting selectivity values and function specified
              M50=Species_i$lmat,maturity_input="length",M=Mnat,binwidth=5,R0=R0, #R0 is starting value - will be found by LIME
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
          hessian = Results_LIME$Sdreport$pdHess
          Converge = (hessian==TRUE & gradient == TRUE)
          print(paste("gradient should be below 0.001, it is",Results_LIME$opt$max_gradient))
          print(paste("does the hessian test pass? ",Results_LIME$Sdreport$pdHess))
          print(paste("SPR estimated as ",Results_LIME$Report$SPR_t,sep=""))
          flush.console()
          F_estimated = Results_LIME$Report$F_t
          NumAtAge_i_j$Population_Extrapolated = as.vector(Results_LIME$Report$N_ta)
          NumAtAge_i_j$EstimatedCatch_Baranov = as.vector(Results_LIME$Report$Cn_ta)
          NumAtAge_rep = do.call("cbind",replicate(dim(Results_LIME$Report$plba)[2],data.frame(NumAtAge_i_j$Population_Extrapolated),simplify = FALSE))
          NumAtLength = colSums(Results_LIME$Report$plba*NumAtAge_rep)
          NumAtLength_i_j$Survival_LengthBased_F_M_Rescaled = NumAtLength/sum(NumAtLength)
          NumAtLength_i_j$Population_Extrapolated = NumAtLength
          CatchAtAge_rep = do.call("cbind",replicate(dim(Results_LIME$Report$plba)[2],data.frame(NumAtAge_i_j$EstimatedCatch_Baranov),simplify = FALSE))
          CatchAtLength = colSums(Results_LIME$Report$plba*CatchAtAge_rep)
          NumAtLength_i_j$EstimatedCatch_Baranov = CatchAtLength
          NumAtLength_i_j$Residuals = NumAtLength_i_j$ExtrapolatedCatch-NumAtLength_i_j$EstimatedCatch_Baranov
          FolderNameSpp = paste(as.character(BiomassCarryingCapacityMSY$Genus[i]),CapStr(as.character(BiomassCarryingCapacityMSY$Species[i])),sep="")
          png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("CatchAtLength",BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
          plot(NumAtLength_i_j$Midpoint,NumAtLength_i_j$ExtrapolatedCatch,pch=16,xlab="Length (cm)",ylab="Catch in Number",cex.axis=1.1,cex.lab=1.1,ylim=c(0,max(NumAtLength_i_j$EstimatedCatch_Baranov,NumAtLength_i_j$ExtrapolatedCatch)))
          points(NumAtLength_i_j$Midpoint,NumAtLength_i_j$EstimatedCatch_Baranov,type="l",lwd=2,col="red")
          title(paste(BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],sep=" "))
          mtext(paste(BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep=", "))
          dev.off()
          png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("CatchAtAge",BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
          plot(NumAtAge_i_j$Age,NumAtAge_i_j$ExtrapolatedCatchInNumber,pch=16,xlab="Age",ylab="Catch in Number",cex.axis=1.1,cex.lab=1.1,ylim=c(0,max(NumAtAge_i_j$EstimatedCatch_Baranov,NumAtAge_i_j$ExtrapolatedCatchInNumber)))
          points(NumAtAge_i_j$Age,NumAtAge_i_j$EstimatedCatch_Baranov,type="l",lwd=2,col="red")
          title(paste(BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],sep=" "))
          mtext(paste(BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep=", "))
          dev.off()
        }
        if(populationReconstructionMethods[v]=="TropFishR")
        {
          NumAtLength_i_j = NumAtLength_v
          NumAtAge_i_j = NumAtAge_v
          print("Starting Length Reconstruction using VPA...")
          flush.console()
          catches_tropfish = lfqCreate(data.frame(length=Sub$cm_converted,date=as.Date("2018-01-02")),Lname="length",Dname="date",bin_size=5,species=paste(BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],sep=""),stock="wpp",Lmin=0)
          catches_tropfish$lmat = Species_i$lmat
          catches_tropfish$catch = catches_tropfish$catch[,1]
          IndexToRemove = c()
          for(u in length(catches_tropfish$catch):1)
          {
              if(catches_tropfish$catch[u]<=0)
              {
                IndexToRemove = c(IndexToRemove,u)
              }
              if(catches_tropfish$catch[u]>0)
              {
                break
              }
          }
          if(!is.null(IndexToRemove))
          {
            catches_tropfish$midLengths = catches_tropfish$midLengths[-IndexToRemove]
            catches_tropfish$catch = catches_tropfish$catch[-IndexToRemove]
          }
          catches_tropfish$Linf = Linf
          catches_tropfish$K = VBK
          scaling = 1
          DataForCatchCurve = data.frame(LogCatch=log(NumAtAge_i_j$ExtrapolatedCatch[seq(which.max(NumAtAge_i_j$ExtrapolatedCatch),length(NumAtAge_i_j$ExtrapolatedCatch),by=1)]),Age=NumAtAge_i_j$Age[seq(which.max(NumAtAge_i_j$ExtrapolatedCatch),length(NumAtAge_i_j$ExtrapolatedCatch),by=1)])
          DataForCatchCurve = DataForCatchCurve[!is.na(DataForCatchCurve$LogCatch) & !is.infinite(DataForCatchCurve$LogCatch),]
          CatchCurveAnalysis = lm(DataForCatchCurve$LogCatch ~ DataForCatchCurve$Age)
          F_estimated = (coef(CatchCurveAnalysis)[2]*-1) - Mnat
          names(F_estimated)=NULL
          catches_tropfish$a = Species_i$var_a/1000
          catches_tropfish$b = Species_i$var_b
          catches_tropfish$M = Mnat
          catches_tropfish$t0=0
          population_TropFish = VPA(param=catches_tropfish,terminalF=F_estimated,analysis_type="CA",plot=FALSE,catch_corFac=scaling,catch_unit=NA)
          Results_TropFish_byLength = data.frame(Midpoint=population_TropFish$classes.num,Population_Extrapolated=as.vector(population_TropFish$annualMeanNr),EstimatedCatch_Baranov=as.vector(population_TropFish$catch_numbers))
          NumAtLength_i_j = merge(NumAtLength_i_j,Results_TropFish_byLength,all.x=TRUE,by=c("Midpoint"))
          NumAtLength_i_j[is.na(NumAtLength_i_j)]=0
          NumAtLength_i_j$Survival_LengthBased_F_M_Rescaled = NumAtLength_i_j$Population_Extrapolated/sum(NumAtLength_i_j$Population_Extrapolated)
          NumAtLength_i_j$Residuals = NumAtLength_i_j$ExtrapolatedCatch-NumAtLength_i_j$EstimatedCatch_Baranov
          #convert results to at age
          PopAtLen = do.call("cbind",replicate(dim(ALK_Len2Age)[2],data.frame(NumAtLength_i_j$Population_Extrapolated),simplify = FALSE))
          PopAtLen = PopAtLen*ALK_Len2Age
          NumAtAge_i_j$Population_Extrapolated = colSums(PopAtLen)
          NumAtAge_i_j$SSB_current = NumAtAge_i_j$Population_Extrapolated*(NumAtAge_i_j$Weight_kg*0.001)*NumAtAge_i_j$Maturity*NumAtAge_i_j$Fecundity
          #Here, since VPA does NOT use any of the catchMSY information, I need to re-find the population at virgin, SSB at virgin, R0, etc. based on the current SPR level and current population as estimated by VPA. Then, replace all these values when the is written out to the large list
          SPR = sum(NumAtAge_i_j$SSB_current)/sum(NumAtAge_i_j$SSB_virginAtAge)
          N_virgin = sum(NumAtAge_i_j$PopulationVirgin)
          N_virgin_start = N_virgin
          min.RSS = function(N_virgin)
          {
            NumAtAge_i_j$PopulationVirgin = NumAtAge_i_j$Survival_M_rescaled*N_virgin
            NumAtAge_i_j$SSB_virginAtAge = NumAtAge_i_j$PopulationVirgin*(NumAtAge_i_j$Weight_kg*0.001)*NumAtAge_i_j$Maturity*NumAtAge_i_j$Fecundity
            SPR = sum(NumAtAge_i_j$SSB_current)/sum(NumAtAge_i_j$SSB_virginAtAge)
           return(sum((SPR - StartValues_i_j_o_p$spr_current)^2))
          }
          Model = optimize(min.RSS,c(0,(N_virgin_start*100000000)),tol=0.001)
          N_virgin = Model$minimum[1]
          NumAtAge_i_j$PopulationVirgin = NumAtAge_i_j$Survival_M_rescaled*N_virgin
          NumAtAge_i_j$SSB_virginAtAge = NumAtAge_i_j$PopulationVirgin*(NumAtAge_i_j$Weight_kg*0.001)*NumAtAge_i_j$Maturity*NumAtAge_i_j$Fecundity
          SSB_virgin=sum(NumAtAge_i_j$SSB_virginAtAge)
          B_virgin_mT = sum(NumAtAge_i_j$PopulationVirgin*(NumAtAge_i_j$Weight_kg*0.001))
          init_max_GT=N_virgin
          digits = nchar(round(init_max_GT,0))
          incrementBy = 10^(digits-3)
          if(digits<=3) {incrementBy=1}
          R0 = findR0(M=Mnat,Linf=Linf,VBK=VBK,a=Species_i$var_a,b=Species_i$var_b,K=B_virgin_mT,maxAge=ceiling(max(NumAtAge_i_j$Age)),matureAge=AgeMat,init_min_GT=1,init_max_GT=init_max_GT,incrementBy=incrementBy)          #lower and upper starting values for seach are 10 and 50  percent of the total population number as extrapolated from catch MSY
          PopVirginAgeCurrent = do.call("rbind",replicate(dim(ALK_Age2Len)[1],transpose(data.frame(NumAtAge_i_j$PopulationVirgin)),simplify = FALSE))
          PopVirginLengthCurrent = PopVirginAgeCurrent*ALK_Age2Len
          NumAtLength_i_j$PopulationVirgin = rowSums(PopVirginLengthCurrent)
          FolderNameSpp = paste(as.character(BiomassCarryingCapacityMSY$Genus[i]),CapStr(as.character(BiomassCarryingCapacityMSY$Species[i])),sep="")
          png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("CatchAtLength",BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
          plot(NumAtLength_i_j$Midpoint,NumAtLength_i_j$ExtrapolatedCatch,pch=16,xlab="Length (cm)",ylab="Catch in Number",cex.axis=1.1,cex.lab=1.1,ylim=c(0,max(NumAtLength_i_j$EstimatedCatch_Baranov,NumAtLength_i_j$ExtrapolatedCatch)))
          points(NumAtLength_i_j$Midpoint,NumAtLength_i_j$EstimatedCatch_Baranov,type="l",lwd=2,col="red")
          title(paste(BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],sep=" "))
          mtext(paste(BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep=", "))
          dev.off()
          png(paste(paste(PATH_plots,FolderNameSpp,sep="/"),paste(paste("CatchAtAge",BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
          plot(NumAtAge_i_j$Age,NumAtAge_i_j$ExtrapolatedCatchInNumber,pch=16,xlab="Age",ylab="Catch in Number",cex.axis=1.1,cex.lab=1.1,ylim=c(0,max(NumAtAge_i_j$EstimatedCatch_Baranov,NumAtAge_i_j$ExtrapolatedCatchInNumber)))
          points(NumAtAge_i_j$Age,NumAtAge_i_j$EstimatedCatch_Baranov,type="l",lwd=2,col="red")
          title(paste(BiomassCarryingCapacityMSY$Genus[i],BiomassCarryingCapacityMSY$Species[i],sep=" "))
          mtext(paste(BiomassCarryingCapacityMSY$WPP[i],paste(BiomassCarryingCapacityMSY$Scenario[i],BiomassCarryingCapacityMSY$Bound[i],sep="-"),LifeHistoryScenarios[p],populationReconstructionMethods[v],sep=", "))
          dev.off()
        }
        ResultsList_i_j_o_p = list(fish_id=BiomassCarryingCapacityMSY$fish_id[i],Genus=BiomassCarryingCapacityMSY$Genus[i],Species=BiomassCarryingCapacityMSY$Species[i],WPP=BiomassCarryingCapacityMSY$WPP[i],CatchMSY_Scenario=BiomassCarryingCapacityMSY$Scenario[i],LifeHistory_Scenario=LifeHistoryScenarios[p],populationReconstructionMethod=populationReconstructionMethods[v],Lmat=Species_i$lmat,Linf=Linf,VBK=VBK,a=Species_i$var_a,b=Species_i$var_b,M=Mnat,F_estimated=F_estimated,SSB_virgin=SSB_virgin,R0=R0,Fecundity_min=StartValues_i_j_o_p$Fecundity_min,Fecundity_max=StartValues_i_j_o_p$Fecundity_max,Selex_Logistic_a_length=Selex_Logistic_a_length,Selex_Logistic_b_length=Selex_Logistic_b_length,Selex_Logistic_a_age=Selex_Logistic_a_age,Selex_Logistic_b_age=Selex_Logistic_b_age,ALK_Len2Age=ALK_Len2Age,ALK_Age2Len=ALK_Age2Len,NumAtLength=NumAtLength_i_j,NumAtAge=NumAtAge_i_j,BoundFromCatchMSY=BiomassCarryingCapacityMSY$Bound[i])
        ListOfPopulations[[length(ListOfPopulations)+1]] = ResultsList_i_j_o_p
        if(i%%50==0)
        {
          print(paste("Working on iteration ",i," of ",dim(BiomassCarryingCapacityMSY)[1],"...",sep=""))
          flush.console()
        }
      }
    }
  }
  saveRDS(ListOfPopulations,file=paste(PATH_output,"ListOfPopulations.RData",sep="/"),ascii=FALSE,version=NULL,compress=TRUE,refhook=NULL)
}
ListOfPopulations = readRDS(file=paste(PATH_output,"ListOfPopulations.RData",sep="/"),refhook=NULL)


###################### Create Summary Table of Fitted Values #################################
SummaryTableFittedValues = data.frame(Genus="",Species="",WPP="",CatchMSY_Scenario="",LifeHistory_Scenario="",populationReconstructionMethod="",F_estimated=0,R0=0,Selex_Logistic_a_length=0,Selex_Logistic_b_length=0,Selex_Logistic_a_age=0,Selex_Logistic_b_age=0,ListIndexNumber=0)
SummaryTableFittedValues = SummaryTableFittedValues[-1,]
for(i in 1:length(ListOfPopulations))
{
  temp = data.frame(Genus=ListOfPopulations[[i]]$Genus,Species=ListOfPopulations[[i]]$Species,WPP=ListOfPopulations[[i]]$WPP,CatchMSY_Scenario=ListOfPopulations[[i]]$CatchMSY_Scenario,LifeHistory_Scenario=ListOfPopulations[[i]]$LifeHistory_Scenario,populationReconstructionMethod=ListOfPopulations[[i]]$populationReconstructionMethod,F_estimated=ListOfPopulations[[i]]$F_estimated,R0=ListOfPopulations[[i]]$R0,Selex_Logistic_a_length=ListOfPopulations[[i]]$Selex_Logistic_a_length,Selex_Logistic_b_length=ListOfPopulations[[i]]$Selex_Logistic_b_length,Selex_Logistic_a_age=ListOfPopulations[[i]]$Selex_Logistic_a_age,Selex_Logistic_b_age=ListOfPopulations[[i]]$Selex_Logistic_b_age,Mnat=ListOfPopulations[[i]]$M,ListIndexNumber=i,CatchMSY_Bound=ListOfPopulations[[i]]$BoundFromCatchMSY)
  SummaryTableFittedValues = rbind(SummaryTableFittedValues,temp)
  if(i%%50==0)
  {
    print(paste("Working on iteration ",i," of ",length(ListOfPopulations),"...",sep=""))
    flush.console()
  }
}
names(BiomassCarryingCapacityMSY)[names(BiomassCarryingCapacityMSY)=="Scenario"]="CatchMSY_Scenario"
names(BiomassCarryingCapacityMSY)[names(BiomassCarryingCapacityMSY)=="Bound"]="CatchMSY_Bound"
SummaryTableFittedValues = merge(SummaryTableFittedValues,BiomassCarryingCapacityMSY,all.x=TRUE,by=c("Genus","Species","WPP","CatchMSY_Scenario","CatchMSY_Bound"))
write.table(SummaryTableFittedValues,paste(PATH_output,"SummaryTableFittedValues.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")



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
      SSB_m = sum(Population_end*(NumAtAge$Weight_kg*0.001)*NumAtAge$Maturity*NumAtAge$Fecundity)
      Recruitment_m = BevHoltRecruitment(Scenario_i$R0,h_steepness,SSB_m,Scenario_i$SSB_virgin)
      SelexToApply = c(0,NumAtAge$Selex_Empirical[-length(NumAtAge$Selex_Empirical)])
      Population_end = Population_start*exp(-1*(Scenario_i$M+(F_project*SelexToApply)))
      TotalDeaths_m = Population_start - Population_end
      Catch_m = ((F_project*SelexToApply)/(Scenario_i$M+(F_project*SelexToApply)))*Population_start*(1-exp(-(Scenario_i$M+(F_project*SelexToApply))))
      NaturalDeaths_m = TotalDeaths_m - Catch_m
      Population_end[Population_end<0]=0   #can't have negative numbers of fish
      Population_end[1] = Recruitment_m
      Population = data.frame(Population_end)
      names(Population) = as.character(Years[m])
      row.names(Population) = as.character(NumAtAge$Age)
      Catch = data.frame(Catch_m)
      names(Catch) = as.character(Years[m])
      row.names(Catch) = as.character(NumAtAge$Age)
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
      SSB_m = sum(Population_end*(NumAtAge$Weight_kg*0.001)*NumAtAge$Maturity*NumAtAge$Fecundity)
      Recruitment_m = BevHoltRecruitment(Scenario_i$R0,h_steepness,SSB_m,Scenario_i$SSB_virgin)
      SelexToApply = c(0,NumAtAge$Selex_Empirical[-length(NumAtAge$Selex_Empirical)])
      Population_end = Population_start*exp(-1*(Scenario_i$M+(F_project*SelexToApply)))
      TotalDeaths_m = Population_start - Population_end
      Catch_m = ((F_project*SelexToApply)/(Scenario_i$M+(F_project*SelexToApply)))*Population_start*(1-exp(-(Scenario_i$M+(F_project*SelexToApply))))
      NaturalDeaths_m = TotalDeaths_m - Catch_m
      Population_end[Population_end<0]=0   #can't have negative numbers of fish
      Population_end[1] = Recruitment_m
      Population_temp = data.frame(Population_end)
      names(Population_temp) = as.character(Years[m])
      row.names(Population_temp) = as.character(NumAtAge$Age)
      Population = cbind(Population,Population_temp)
      Catch_temp = data.frame(Catch_m)
      names(Catch_temp) = as.character(Years[m])
      row.names(Catch_temp) = as.character(NumAtAge$Age)
      Catch = cbind(Catch,Catch_temp)
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
      NaturalDeathsYear_temp = sum(NaturalDeaths_temp)
      names(NaturalDeathsYear_temp) = as.character(Years[m])
      NaturalDeathsYear = c(NaturalDeathsYear,NaturalDeathsYear_temp)
    }
  }
  Scenario_i$F_project = F_project
  Scenario_i$ProjectedPopulationAtAge = Population
  Scenario_i$ProjectedCatchAtAge = Catch
  Scenario_i$ProjectedNaturalDeathsAtAge = NaturalDeaths
  Scenario_i$SSB = SSB
  Scenario_i$Recruitment = Recruitment
  Scenario_i$ProjectedPopulationYear = Population_Year
  Scenario_i$ProjectedBiomassYear = Biomass
  Scenario_i$ProjectedCatchYear = Catch_Year
  Scenario_i$ProjectedNaturalDeathsYear = NaturalDeathsYear
  return(Scenario_i)
}
#Setup and Read Constant Catch Based Projection Function
projectionFunction_constantCatch = function(Scenario_i,yearsToProject,yearToStartProjection,h_steepness)
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
      SSB_m = sum(Population_end*(NumAtAge$Weight_kg*0.001)*NumAtAge$Maturity*NumAtAge$Fecundity)
      Recruitment_m = BevHoltRecruitment(Scenario_i$R0,h_steepness,SSB_m,Scenario_i$SSB_virgin)
      SelexToApply = c(0,NumAtAge$Selex_Empirical[-length(NumAtAge$Selex_Empirical)])
      Population_end_afterM = Population_start*exp(-1*(Scenario_i$M))
      NaturalDeaths_m = Population_start - Population_end_afterM
      Catch_m = NumAtAge$ExtrapolatedCatchInNumber
      Population_end = Population_start - (NaturalDeaths_m + Catch_m)
      Population_end[Population_end<0]=0   #can't have negative numbers of fish
      Population_end[1] = Recruitment_m
      Population = data.frame(Population_end)
      names(Population) = as.character(Years[m])
      row.names(Population) = as.character(NumAtAge$Age)
      Catch = data.frame(Catch_m)
      names(Catch) = as.character(Years[m])
      row.names(Catch) = as.character(NumAtAge$Age)
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
      SSB_m = sum(Population_end*(NumAtAge$Weight_kg*0.001)*NumAtAge$Maturity*NumAtAge$Fecundity)
      Recruitment_m = BevHoltRecruitment(Scenario_i$R0,h_steepness,SSB_m,Scenario_i$SSB_virgin)
      SelexToApply = c(0,NumAtAge$Selex_Empirical[-length(NumAtAge$Selex_Empirical)])
      Population_end_afterM = Population_start*exp(-1*(Scenario_i$M))
      NaturalDeaths_m = Population_start - Population_end_afterM
      Catch_m = NumAtAge$ExtrapolatedCatchInNumber
      Population_end = Population_start - (NaturalDeaths_m + Catch_m)
      Population_end[Population_end<0]=0   #can't have negative numbers of fish
      Population_end[1] = Recruitment_m
      Population_temp = data.frame(Population_end)
      names(Population_temp) = as.character(Years[m])
      row.names(Population_temp) = as.character(NumAtAge$Age)
      Population = cbind(Population,Population_temp)
      Catch_temp = data.frame(Catch_m)
      names(Catch_temp) = as.character(Years[m])
      row.names(Catch_temp) = as.character(NumAtAge$Age)
      Catch = cbind(Catch,Catch_temp)
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
      NaturalDeathsYear_temp = sum(NaturalDeaths_temp)
      names(NaturalDeathsYear_temp) = as.character(Years[m])
      NaturalDeathsYear = c(NaturalDeathsYear,NaturalDeathsYear_temp)
    }
  }
  Scenario_i$F_project = NA
  Scenario_i$ProjectedPopulationAtAge = Population
  Scenario_i$ProjectedCatchAtAge = Catch
  Scenario_i$ProjectedNaturalDeathsAtAge = NaturalDeaths
  Scenario_i$SSB = SSB
  Scenario_i$Recruitment = Recruitment
  Scenario_i$ProjectedPopulationYear = Population_Year
  Scenario_i$ProjectedBiomassYear = Biomass
  Scenario_i$ProjectedCatchYear = Catch_Year
  Scenario_i$ProjectedNaturalDeathsYear = NaturalDeathsYear
  return(Scenario_i)
}
#Setup and Read Linear Recruitment Projection Function
projectionFunction_LinearRecruitment = function(Scenario_i,yearsToProject,yearToStartProjection,F_project)
{
  Years = seq(yearToStartProjection,(yearToStartProjection+yearsToProject),by=1)
  NumAtAge = Scenario_i$NumAtAge
  #Determine Linear Recruitment Function
  ConstantRecruitment = Scenario_i$R0
  LinearModel_recruitBetween0_20percent = lm(c(0,ConstantRecruitment) ~ c(0,0.2))
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
      SSB_m = sum(Population_end*(NumAtAge$Weight_kg*0.001)*NumAtAge$Maturity*NumAtAge$Fecundity)
      SPR = SSB_m/Scenario_i$SSB_virgin
      if(SPR>=0.2)
      {
        Recruitment_m =  ConstantRecruitment
      }
      if(SPR<0.2)
      {
        Recruitment_m = (SPR*LinearModel_recruitBetween0_20percent$coef[2]) + LinearModel_recruitBetween0_20percent$coef[1]
      }
      SelexToApply = c(0,NumAtAge$Selex_Empirical[-length(NumAtAge$Selex_Empirical)])
      Population_end = Population_start*exp(-1*(Scenario_i$M+(F_project*SelexToApply)))
      TotalDeaths_m = Population_start - Population_end
      Catch_m = ((F_project*SelexToApply)/(Scenario_i$M+(F_project*SelexToApply)))*Population_start*(1-exp(-(Scenario_i$M+(F_project*SelexToApply))))
      NaturalDeaths_m = TotalDeaths_m - Catch_m
      Population_end[Population_end<0]=0   #can't have negative numbers of fish
      Population_end[1] = Recruitment_m
      Population = data.frame(Population_end)
      names(Population) = as.character(Years[m])
      row.names(Population) = as.character(NumAtAge$Age)
      Catch = data.frame(Catch_m)
      names(Catch) = as.character(Years[m])
      row.names(Catch) = as.character(NumAtAge$Age)
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
      SSB_m = sum(Population_end*(NumAtAge$Weight_kg*0.001)*NumAtAge$Maturity*NumAtAge$Fecundity)
      SPR = SSB_m/Scenario_i$SSB_virgin
      if(SPR>=0.2)
      {
        Recruitment_m =  ConstantRecruitment
      }
      if(SPR<0.2)
      {
        Recruitment_m = (SPR*LinearModel_recruitBetween0_20percent$coef[2]) + LinearModel_recruitBetween0_20percent$coef[1]
      }
      SelexToApply = c(0,NumAtAge$Selex_Empirical[-length(NumAtAge$Selex_Empirical)])
      Population_end = Population_start*exp(-1*(Scenario_i$M+(F_project*SelexToApply)))
      TotalDeaths_m = Population_start - Population_end
      Catch_m = ((F_project*SelexToApply)/(Scenario_i$M+(F_project*SelexToApply)))*Population_start*(1-exp(-(Scenario_i$M+(F_project*SelexToApply))))
      NaturalDeaths_m = TotalDeaths_m - Catch_m
      Population_end[Population_end<0]=0   #can't have negative numbers of fish
      Population_end[1] = Recruitment_m      #MOVED ASSIGNMENT OF RECRUITMENT TO HERE, AND CHANGED IT TO be inserted in Population_end[1]
      Population_temp = data.frame(Population_end)
      names(Population_temp) = as.character(Years[m])
      row.names(Population_temp) = as.character(NumAtAge$Age)
      Population = cbind(Population,Population_temp)
      Catch_temp = data.frame(Catch_m)
      names(Catch_temp) = as.character(Years[m])
      row.names(Catch_temp) = as.character(NumAtAge$Age)
      Catch = cbind(Catch,Catch_temp)
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
      NaturalDeathsYear_temp = sum(NaturalDeaths_temp)
      names(NaturalDeathsYear_temp) = as.character(Years[m])
      NaturalDeathsYear = c(NaturalDeathsYear,NaturalDeathsYear_temp)
    }
  }
  Scenario_i$F_project = F_project
  Scenario_i$ProjectedPopulationAtAge = Population
  Scenario_i$ProjectedCatchAtAge = Catch
  Scenario_i$ProjectedNaturalDeathsAtAge = NaturalDeaths
  Scenario_i$Recruitment = Recruitment
  Scenario_i$SSB = SSB
  Scenario_i$ProjectedPopulationYear = Population_Year
  Scenario_i$ProjectedBiomassYear = Biomass
  Scenario_i$ProjectedCatchYear = Catch_Year
  Scenario_i$ProjectedNaturalDeathsYear = NaturalDeathsYear
  return(Scenario_i)
}
#Setup and Read Linear Recruitment and Constant Catch Projection Function
projectionFunction_LinearRecruitment_ConstantCatch = function(Scenario_i,yearsToProject,yearToStartProjection)
{
  Years = seq(yearToStartProjection,(yearToStartProjection+yearsToProject),by=1)
  NumAtAge = Scenario_i$NumAtAge
  #Determine Linear Recruitment Function
  ConstantRecruitment = Scenario_i$R0
  LinearModel_recruitBetween0_20percent = lm(c(0,ConstantRecruitment) ~ c(0,0.2))
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
      SSB_m = sum(Population_end*(NumAtAge$Weight_kg*0.001)*NumAtAge$Maturity*NumAtAge$Fecundity)
      SPR = SSB_m/Scenario_i$SSB_virgin
      if(SPR>=0.2)
      {
        Recruitment_m =  ConstantRecruitment
      }
      if(SPR<0.2)
      {
        Recruitment_m = (SPR*LinearModel_recruitBetween0_20percent$coef[2]) + LinearModel_recruitBetween0_20percent$coef[1]
      }
      SelexToApply = c(0,NumAtAge$Selex_Empirical[-length(NumAtAge$Selex_Empirical)])
      Population_end_afterM = Population_start*exp(-1*(Scenario_i$M))
      NaturalDeaths_m = Population_start - Population_end_afterM
      Catch_m = NumAtAge$ExtrapolatedCatchInNumber
      Population_end = Population_start - (NaturalDeaths_m + Catch_m)
      Population_end[Population_end<0]=0   #can't have negative numbers of fish
      Population_end[1] = Recruitment_m
      Population = data.frame(Population_end)
      names(Population) = as.character(Years[m])
      row.names(Population) = as.character(NumAtAge$Age)
      Catch = data.frame(Catch_m)
      names(Catch) = as.character(Years[m])
      row.names(Catch) = as.character(NumAtAge$Age)
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
      SSB_m = sum(Population_end*(NumAtAge$Weight_kg*0.001)*NumAtAge$Maturity*NumAtAge$Fecundity)
      SPR = SSB_m/Scenario_i$SSB_virgin
      if(SPR>=0.2)
      {
        Recruitment_m =  ConstantRecruitment
      }
      if(SPR<0.2)
      {
        Recruitment_m = (SPR*LinearModel_recruitBetween0_20percent$coef[2]) + LinearModel_recruitBetween0_20percent$coef[1]
      }
      SelexToApply = c(0,NumAtAge$Selex_Empirical[-length(NumAtAge$Selex_Empirical)])
      Population_end_afterM = Population_start*exp(-1*(Scenario_i$M))
      NaturalDeaths_m = Population_start - Population_end_afterM
      Catch_m = NumAtAge$ExtrapolatedCatchInNumber
      Population_end = Population_start - (NaturalDeaths_m + Catch_m)
      Population_end[Population_end<0]=0   #can't have negative numbers of fish
      Population_end[1] = Recruitment_m
      Population_temp = data.frame(Population_end)
      names(Population_temp) = as.character(Years[m])
      row.names(Population_temp) = as.character(NumAtAge$Age)
      Population = cbind(Population,Population_temp)
      Catch_temp = data.frame(Catch_m)
      names(Catch_temp) = as.character(Years[m])
      row.names(Catch_temp) = as.character(NumAtAge$Age)
      Catch = cbind(Catch,Catch_temp)
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
      NaturalDeathsYear_temp = sum(NaturalDeaths_temp)
      names(NaturalDeathsYear_temp) = as.character(Years[m])
      NaturalDeathsYear = c(NaturalDeathsYear,NaturalDeathsYear_temp)
    }
  }
  Scenario_i$F_project = NA
  Scenario_i$ProjectedPopulationAtAge = Population
  Scenario_i$ProjectedCatchAtAge = Catch
  Scenario_i$ProjectedNaturalDeathsAtAge = NaturalDeaths
  Scenario_i$SSB = SSB
  Scenario_i$Recruitment = Recruitment
  Scenario_i$ProjectedPopulationYear = Population_Year
  Scenario_i$ProjectedBiomassYear = Biomass
  Scenario_i$ProjectedCatchYear = Catch_Year
  Scenario_i$ProjectedNaturalDeathsYear = NaturalDeathsYear
  return(Scenario_i)
}



#PROJECTION LOOP
if(runProjections==TRUE)
{
  ListOfPopulations = readRDS(file=paste(PATH_output,"ListOfPopulations.RData",sep="/"),refhook=NULL)
  Years = seq(yearToStartProjection,(yearToStartProjection+yearsToProject),by=1)
  ListOfPopulations_WithProjections = list()
  for(i in 1:length(ListOfPopulations))
  {
#  for(i in 16444:18270)   #There are 18270 populations to project on 10 processors: 1-1827,1828-3654,3655-5481,5482-7308,7309-9135,9136-10962,10963-12789,12790-14616,14617-16443,16444-18270
#  {
    for(v in 1:length(FishingMortaltiyProjectionScenarios))
    {
      if(FishingMortaltiyProjectionScenarios[v]!="F_rebuild_20YrsMSC_max" & FishingMortaltiyProjectionScenarios[v]!="F_rebuild_Tmin_oneGeneration" & FishingMortaltiyProjectionScenarios[v]!="F_rebuild_10YrsUS")
      {
        #Pull out the population from the list
        Scenario_i = ListOfPopulations[[i]]
        NumAtAge = Scenario_i$NumAtAge
        Species_i = Species_cpue[Species_cpue$fish_id==Scenario_i$fish_id,]
        if(Scenario_i$WPP!="EEZ")
        {
          Sub = Lengths_cpue[Lengths_cpue$fish_id==Scenario_i$fish_id & Lengths_cpue$wpp1==as.numeric(Scenario_i$WPP),]
        }
        if(Scenario_i$WPP=="EEZ")
        {
          Sub = Lengths_cpue[Lengths_cpue$fish_id==Scenario_i$fish_id,]
        }
        #Find the Fishing mortality for each scenario
        if(FishingMortaltiyProjectionScenarios[v]=="CurrentF")
        {
          F_project = Scenario_i$F_estimated
        }
        if(FishingMortaltiyProjectionScenarios[v]=="F_equal_0")
        {
          F_project = 0
        }
        if(FishingMortaltiyProjectionScenarios[v]=="F_half")
        {
          F_project = Scenario_i$F_estimated/2
        }
        if(FishingMortaltiyProjectionScenarios[v]=="F_double")
        {
          F_project = Scenario_i$F_estimated*2
        }
        if(FishingMortaltiyProjectionScenarios[v]=="F_at_SPR30" | FishingMortaltiyProjectionScenarios[v]=="F_at_75percSPR30")
        {
          if(Scenario_i$populationReconstructionMethod=="LBSPR")
          {
            MyPars = new("LB_pars")
            MyPars@Species = Scenario_i$Species
            MyPars@Linf = Scenario_i$Linf
            MyPars@L50 = Scenario_i$Lmat #length at 50% maturity
            MyPars@L95 = Species_i$lopt #length at 95% maturity - we don't know this so use lopt as a proxy. We don't know this because we are assuming knife edge maturity and don't have a study to know what the real maturity schedule looks like
            MyPars@MK = Scenario_i$M/Scenario_i$VBK
            MyPars@L_units = "cm"
            MyPars@BinMin=min(Scenario_i$NumAtLength$Midpoint)-2.5
            MyPars@BinWidth=Scenario_i$NumAtLength$Midpoint[2]-Scenario_i$NumAtLength$Midpoint[1]
            MyPars@BinMax=ceiling(Species_i$lmax/5)*5
            Len1 = new("LB_lengths",LB_pars=MyPars,file=Sub$cm_converted,sataType="raw")
            myFit1 = LBSPRfit(MyPars, Len1)
            myFit1@Ests
            plotSize(myFit1)
            plotMat(myFit1)
            #simulate population
            SimPars=MyPars
            SimPars@SL50 = Scenario_i$Lmat #length at 50% maturity
            SimPars@SL95 = Species_i$lopt #length at 95% maturity - we don't know this so use lopt as a proxy. We don't know this because we are assuming knife edge maturity and don't have a study to know what the real maturity schedule looks like
            SimPars@SPR = ControlRule_SPR
            SimPars@R0=Scenario_i$R0
            MySim = LBSPRsim(SimPars,Control = list("modtype"="absel"))
            F_project = MySim@FM*Scenario_i$M
          }
          if(Scenario_i$populationReconstructionMethod=="SolveBaranov")
          {
            NumAtAge_thisScenario = subset(NumAtAge,select=-c(Population_Extrapolated,EstimatedCatch_Baranov,Residuals,SSB_CurrentAtAge))
            EstNumFish_Population = sum(NumAtAge_thisScenario$PopulationVirgin)*ControlRule_SPR
            F_start = 1.0
            min.RSS = function(F_start)
            {
              NumAtAge_thisScenario$Survival_F_M = exp(-1*(Scenario_i$M+(F_start*NumAtAge_thisScenario$Selex_Empirical))*NumAtAge_thisScenario$Age)
              NumAtAge_thisScenario$Survival_F_M_Rescaled = NumAtAge_thisScenario$Survival_F_M/sum(NumAtAge_thisScenario$Survival_F_M)
              NumAtAge_thisScenario$Population_Extrapolated = NumAtAge_thisScenario$Survival_F_M_Rescaled*EstNumFish_Population
              NumAtAge_thisScenario$SSB_CurrentAtAge = NumAtAge_thisScenario$Population_Extrapolated*(NumAtAge_thisScenario$Weight_kg*0.001)*NumAtAge_thisScenario$Maturity*NumAtAge_thisScenario$Fecundity
              SSB_current = sum(NumAtAge_thisScenario$SSB_CurrentAtAge)
              SPR_current = SSB_current/Scenario_i$SSB_virgin
              return(sum((SPR_current - ControlRule_SPR)^2))
            }
            Model = optimize(min.RSS,c(0,10),tol=0.001)
            NumAtAge_thisScenario$Survival_F_M = exp(-1*(Scenario_i$M+(Model$minimum[1]*NumAtAge_thisScenario$Selex_Empirical))*NumAtAge_thisScenario$Age)     #"pseudosurvival"
            NumAtAge_thisScenario$Survival_F_M_Rescaled = NumAtAge_thisScenario$Survival_F_M/sum(NumAtAge_thisScenario$Survival_F_M)
            NumAtAge_thisScenario$Population_Extrapolated = NumAtAge_thisScenario$Survival_F_M_Rescaled*EstNumFish_Population
            NumAtAge_thisScenario$Survival_F_M=NULL
            NumAtAge_thisScenario$Survival_F_M_Rescaled=NULL
            F_start = 1.0
            min.RSS = function(F_start)
            {
              EstimatedCatch_Baranov = BaranovCatch(NumAtAge_thisScenario$Age,NumAtAge_thisScenario$Length,NumAtAge_thisScenario$Population_Extrapolated,Scenario_i$M,F_start,NumAtAge_thisScenario$Selex_Empirical)
              return(sum((EstimatedCatch_Baranov - NumAtAge_thisScenario$ExtrapolatedCatchInNumber)^2))
            }
            Model = optimize(min.RSS,c(0,10),tol=0.001)
            F_project = Model$minimum[1]
          }
          if(Scenario_i$populationReconstructionMethod=="BevertonHoltInstantaneousMortality")
          {
            F_start = 1.0
            min.RSS = function(F_start)
            {
              LenVec = Sub$cm_converted
              LenVec = LenVec[!is.na(LenVec)]
              SPR = computeLengthSPR(Scenario_i$M,Scenario_i$Linf,Scenario_i$VBK,Species_i$var_a,Species_i$var_b,Species_i$lmat,F_start,LenVec)
              return(sum((SPR - ControlRule_SPR)^2))
            }
            Model = optimize(min.RSS,c(0,10),tol=0.001)
            F_project = Model$minimum[1]
          }
          if(Scenario_i$populationReconstructionMethod=="LIME")
          {
            Len50_selex = 30
            min.RSS = function(Len50_selex)
            {
              Sel = SelectivityLogistic_Length(Len50_selex,Scenario_i$Selex_Logistic_a_length,Scenario_i$Selex_Logistic_b_length)
              return(sum((Sel - 0.5)^2))
            }
            Model = optimize(min.RSS,c(0,200),tol=0.001)
            Len50_selex = Model$minimum[1]
            Len95_selex = 50
            min.RSS = function(Len95_selex)
            {
              Sel = SelectivityLogistic_Length(Len95_selex,Scenario_i$Selex_Logistic_a_length,Scenario_i$Selex_Logistic_b_length)
              return(sum((Sel - 0.95)^2))
            }
            Model = optimize(min.RSS,c(0,200),tol=0.001)
            Len95_selex = Model$minimum[1]
            lh = create_lh_list(vbk=Scenario_i$VBK,linf=Scenario_i$Linf,lwa=Species_i$var_a,lwb=Species_i$var_b,
                S50=Len50_selex,S95=Len95_selex,selex_input="length",selex_type=c("logistic"), #Starting selectivity values and function specified
                M50=Species_i$lmat,maturity_input="length",M=Scenario_i$M,binwidth=5,R0=Scenario_i$R0, #R0 is starting value - will be found by LIME
                rho=0,h=0.8,  #steepness (h) and cruitment autocorrelation (rho)
                nseasons=1,nfleets=1,CVlen=0.2,SigmaR=0.1,AgeMax=max(NumAtAge$Age)) #if we could get MONTHLY effort, and split catches per fleet that would go here
            LF=matrix(Scenario_i$NumAtLength$NumSampled,nrow=1)
            colnames(LF)=Scenario_i$NumAtLength$Length_max
            rownames(LF)=max(SpeciesCatch_mT_Scaled$Year)
            Catch_matrix = as.matrix(sum(Scenario_i$NumAtAge$ExtrapolatedCatchInNumber*Scenario_i$NumAtAge$Weight_kg*1000))
            Catch_matrix = Catch_matrix[,length(Catch_matrix)]
            data_all=list(years=max(SpeciesCatch_mT_Scaled$Year),LF=LF,C_ft=as.matrix(Catch_matrix))   #list of years, length distribution (as a matrix of one row), catches and number of observed catches
            inputs_all=create_inputs(lh=lh,input_data=data_all)
            vals_selex_ft = matrix(Scenario_i$NumAtLength$Selex_Fitted,nrow=1)
            Results_LIME = run_LIME(modpath=NULL,input=inputs_all,data_avail="Catch_LC1",C_type=2,est_selex_f=FALSE,vals_selex_ft=vals_selex_ft,derive_quants=TRUE,est_totalF=TRUE,Fpen=0)       #C_type is catch in weight, not numbers
            gradient = Results_LIME$opt$max_gradient <= 0.001
            hessian = Results_LIME$Sdreport$pdHess
            Converge = (hessian==TRUE & gradient == TRUE)
            print(paste("gradient should be below 0.001, it is",Results_LIME$opt$max_gradient))
            print(paste("does the hessian test pass? ",Results_LIME$Sdreport$pdHess))
            print(paste("SPR estimated as ",Results_LIME$Report$SPR_t,sep=""))
            flush.console()
            F_project = Results_LIME$Derived$F30
          }
          if(Scenario_i$populationReconstructionMethod=="TropFishR")
          {
            SPR_current = sum(Scenario_i$NumAtAge$Survival_M_F_rescaled*(Scenario_i$NumAtAge$Weight_kg*0.001)*Scenario_i$NumAtAge$Maturity*Scenario_i$NumAtAge$Fecundity)/sum(Scenario_i$NumAtAge$Survival_M_rescaled*(Scenario_i$NumAtAge$Weight_kg*0.001)*Scenario_i$NumAtAge$Maturity*Scenario_i$NumAtAge$Fecundity)
            F_start=0.5
            min.RSS = function(F_start)
            {
              Survival_F_M = exp(-1*(Scenario_i$M+(Scenario_i$NumAtAge$Selex_Empirical*F_start))*Scenario_i$NumAtAge$Age)
              Survival_F_M_rescaled = Survival_F_M/sum(Survival_F_M)
              SPR = sum(Survival_F_M_rescaled*(Scenario_i$NumAtAge$Weight_kg*0.001)*Scenario_i$NumAtAge$Maturity*Scenario_i$NumAtAge$Fecundity)/sum(Scenario_i$NumAtAge$Survival_M_rescaled*(Scenario_i$NumAtAge$Weight_kg*0.001)*Scenario_i$NumAtAge$Maturity*Scenario_i$NumAtAge$Fecundity)
              return(sum((SPR - ControlRule_SPR)^2))
            }
            Model = optimize(min.RSS,c(0,10),tol=0.001)
            F_project = Model$minimum[1]
          }
          if(FishingMortaltiyProjectionScenarios[v]=="F_at_75percSPR30")
          {
            F_project = F_project*0.75
          }
        }
        if(FishingMortaltiyProjectionScenarios[v]=="F_at_SPR20")
        {
          if(Scenario_i$populationReconstructionMethod=="LBSPR")
          {
            MyPars = new("LB_pars")
            MyPars@Species = Scenario_i$Species
            MyPars@Linf = Scenario_i$Linf
            MyPars@L50 = Scenario_i$Lmat #length at 50% maturity
            MyPars@L95 = Species_i$lopt #length at 95% maturity - we don't know this so use lopt as a proxy. We don't know this because we are assuming knife edge maturity and don't have a study to know what the real maturity schedule looks like
            MyPars@MK = Scenario_i$M/Scenario_i$VBK
            MyPars@L_units = "cm"
            MyPars@BinMin=min(Scenario_i$NumAtLength$Midpoint)-2.5
            MyPars@BinWidth=Scenario_i$NumAtLength$Midpoint[2]-Scenario_i$NumAtLength$Midpoint[1]
            MyPars@BinMax=ceiling(Species_i$lmax/5)*5
            Len1 = new("LB_lengths",LB_pars=MyPars,file=Sub$cm_converted,sataType="raw")
            myFit1 = LBSPRfit(MyPars, Len1)
            myFit1@Ests
            plotSize(myFit1)
            plotMat(myFit1)
            #simulate population
            SimPars=MyPars
            SimPars@SL50 = Scenario_i$Lmat #length at 50% maturity
            SimPars@SL95 = Species_i$lopt #length at 95% maturity - we don't know this so use lopt as a proxy. We don't know this because we are assuming knife edge maturity and don't have a study to know what the real maturity schedule looks like
            SimPars@SPR = Limit_SPR
            SimPars@R0=Scenario_i$R0
            MySim = LBSPRsim(SimPars,Control = list("modtype"="absel"))
            F_project = MySim@FM*Scenario_i$M
          }
          if(Scenario_i$populationReconstructionMethod=="SolveBaranov")
          {
            NumAtAge_thisScenario = subset(NumAtAge,select=-c(Population_Extrapolated,EstimatedCatch_Baranov,Residuals,SSB_CurrentAtAge))
            EstNumFish_Population = sum(NumAtAge_thisScenario$PopulationVirgin)*ControlRule_SPR
            F_start = 1.0
            min.RSS = function(F_start)
            {
              NumAtAge_thisScenario$Survival_F_M = exp(-1*(Scenario_i$M+(F_start*NumAtAge_thisScenario$Selex_Empirical))*NumAtAge_thisScenario$Age)
              NumAtAge_thisScenario$Survival_F_M_Rescaled = NumAtAge_thisScenario$Survival_F_M/sum(NumAtAge_thisScenario$Survival_F_M)
              NumAtAge_thisScenario$Population_Extrapolated = NumAtAge_thisScenario$Survival_F_M_Rescaled*EstNumFish_Population
              NumAtAge_thisScenario$SSB_CurrentAtAge = NumAtAge_thisScenario$Population_Extrapolated*(NumAtAge_thisScenario$Weight_kg*0.001)*NumAtAge_thisScenario$Maturity*NumAtAge_thisScenario$Fecundity
              SSB_current = sum(NumAtAge_thisScenario$SSB_CurrentAtAge)
              SPR_current = SSB_current/Scenario_i$SSB_virgin
              return(sum((SPR_current - Limit_SPR)^2))
            }
            Model = optimize(min.RSS,c(0,10),tol=0.001)
            NumAtAge_thisScenario$Survival_F_M = exp(-1*(Scenario_i$M+(Model$minimum[1]*NumAtAge_thisScenario$Selex_Empirical))*NumAtAge_thisScenario$Age)     #"pseudosurvival"
            NumAtAge_thisScenario$Survival_F_M_Rescaled = NumAtAge_thisScenario$Survival_F_M/sum(NumAtAge_thisScenario$Survival_F_M)
            NumAtAge_thisScenario$Population_Extrapolated = NumAtAge_thisScenario$Survival_F_M_Rescaled*EstNumFish_Population
            NumAtAge_thisScenario$Survival_F_M=NULL
            NumAtAge_thisScenario$Survival_F_M_Rescaled=NULL
            F_start = 1.0
            min.RSS = function(F_start)
            {
              EstimatedCatch_Baranov = BaranovCatch(NumAtAge_thisScenario$Age,NumAtAge_thisScenario$Length,NumAtAge_thisScenario$Population_Extrapolated,Scenario_i$M,F_start,NumAtAge_thisScenario$Selex_Empirical)
              return(sum((EstimatedCatch_Baranov - NumAtAge_thisScenario$ExtrapolatedCatchInNumber)^2))
            }
            Model = optimize(min.RSS,c(0,10),tol=0.001)
            F_project = Model$minimum[1]
          }
          if(Scenario_i$populationReconstructionMethod=="BevertonHoltInstantaneousMortality")
          {
            F_start = 1.0
            min.RSS = function(F_start)
            {
              LenVec = Sub$cm_converted
              LenVec = LenVec[!is.na(LenVec)]
              SPR = computeLengthSPR(Scenario_i$M,Scenario_i$Linf,Scenario_i$VBK,Species_i$var_a,Species_i$var_b,Species_i$lmat,F_start,LenVec)
              return(sum((SPR - Limit_SPR)^2))
            }
            Model = optimize(min.RSS,c(0,10),tol=0.001)
            F_project = Model$minimum[1]
          }
          if(Scenario_i$populationReconstructionMethod=="LIME")
          {
            Len50_selex = 30
            min.RSS = function(Len50_selex)
            {
              Sel = SelectivityLogistic_Length(Len50_selex,Scenario_i$Selex_Logistic_a_length,Scenario_i$Selex_Logistic_b_length)
              return(sum((Sel - 0.5)^2))
            }
            Model = optimize(min.RSS,c(0,200),tol=0.001)
            Len50_selex = Model$minimum[1]
            Len95_selex = 50
            min.RSS = function(Len95_selex)
            {
              Sel = SelectivityLogistic_Length(Len95_selex,Scenario_i$Selex_Logistic_a_length,Scenario_i$Selex_Logistic_b_length)
              return(sum((Sel - 0.95)^2))
            }
            Model = optimize(min.RSS,c(0,200),tol=0.001)
            Len95_selex = Model$minimum[1]
            lh = create_lh_list(vbk=Scenario_i$VBK,linf=Scenario_i$Linf,lwa=Species_i$var_a,lwb=Species_i$var_b,
                S50=Len50_selex,S95=Len95_selex,selex_input="length",selex_type=c("logistic"), #Starting selectivity values and function specified
                M50=Species_i$lmat,maturity_input="length",M=Scenario_i$M,binwidth=5,R0=Scenario_i$R0, #R0 is starting value - will be found by LIME
                rho=0,h=0.8,  #steepness (h) and cruitment autocorrelation (rho)
                nseasons=1,nfleets=1,CVlen=0.2,SigmaR=0.1,AgeMax=max(NumAtAge$Age)) #if we could get MONTHLY effort, and split catches per fleet that would go here
            LF=matrix(Scenario_i$NumAtLength$NumSampled,nrow=1)
            colnames(LF)=Scenario_i$NumAtLength$Length_max
            rownames(LF)=max(SpeciesCatch_mT_Scaled$Year)
            Catch_matrix = as.matrix(sum(Scenario_i$NumAtAge$ExtrapolatedCatchInNumber*Scenario_i$NumAtAge$Weight_kg*1000))
            Catch_matrix = Catch_matrix[,length(Catch_matrix)]
            data_all=list(years=max(SpeciesCatch_mT_Scaled$Year),LF=LF,C_ft=as.matrix(Catch_matrix))   #list of years, length distribution (as a matrix of one row), catches and number of observed catches
            inputs_all=create_inputs(lh=lh,input_data=data_all)
            vals_selex_ft = matrix(Scenario_i$NumAtLength$Selex_Fitted,nrow=1)
            Results_LIME = run_LIME(modpath=NULL,input=inputs_all,data_avail="Catch_LC1",C_type=2,est_selex_f=FALSE,vals_selex_ft=vals_selex_ft,derive_quants=TRUE,est_totalF=TRUE,Fpen=0)       #C_type is catch in weight, not numbers
            gradient = Results_LIME$opt$max_gradient <= 0.001
            hessian = Results_LIME$Sdreport$pdHess
            Converge = (hessian==TRUE & gradient == TRUE)
            print(paste("gradient should be below 0.001, it is",Results_LIME$opt$max_gradient))
            print(paste("does the hessian test pass? ",Results_LIME$Sdreport$pdHess))
            print(paste("SPR estimated as ",Results_LIME$Report$SPR_t,sep=""))
            flush.console()
            F_project = tryCatch(
            {
              uniroot(calc_ref,lower=0,upper=400,ages=Scenario_i$NumAtAge$Age,Mat_a=Scenario_i$NumAtAge$Maturity,W_a=Scenario_i$NumAtAge$Weight_kg*1000,M=Scenario_i$M,S_fa=t(as.matrix(Scenario_i$NumAtAge$Selex_Fitted)),ref = 0.2)$root
            },
            error=function(cond)
            {
              message("Can't find using 'uniroot' - will approximate another way...")
              DF = data.frame(SPRtarg=c(30,40),FatTarg=c(Results_LIME$Derived$F30,Results_LIME$Derived$F40))
              LinApprox = lm(DF$FatTarg ~ DF$SPRtarg)
              F_approx = (coef(LinApprox)[2]*20 + coef(LinApprox)[1])
              names(F_approx)=NULL
              return(F_approx)
            })
          }
          if(Scenario_i$populationReconstructionMethod=="TropFishR")
          {
            SPR_current = sum(Scenario_i$NumAtAge$Survival_M_F_rescaled*(Scenario_i$NumAtAge$Weight_kg*0.001)*Scenario_i$NumAtAge$Maturity*Scenario_i$NumAtAge$Fecundity)/sum(Scenario_i$NumAtAge$Survival_M_rescaled*(Scenario_i$NumAtAge$Weight_kg*0.001)*Scenario_i$NumAtAge$Maturity*Scenario_i$NumAtAge$Fecundity)
            F_start=0.5
            min.RSS = function(F_start)
            {
              Survival_F_M = exp(-1*(Scenario_i$M+(Scenario_i$NumAtAge$Selex_Empirical*F_start))*Scenario_i$NumAtAge$Age)
              Survival_F_M_rescaled = Survival_F_M/sum(Survival_F_M)
              SPR = sum(Survival_F_M_rescaled*(Scenario_i$NumAtAge$Weight_kg*0.001)*Scenario_i$NumAtAge$Maturity*Scenario_i$NumAtAge$Fecundity)/sum(Scenario_i$NumAtAge$Survival_M_rescaled*(Scenario_i$NumAtAge$Weight_kg*0.001)*Scenario_i$NumAtAge$Maturity*Scenario_i$NumAtAge$Fecundity)
              return(sum((SPR - Limit_SPR)^2))
            }
            Model = optimize(min.RSS,c(0,10),tol=0.001)
            F_project = Model$minimum[1]
          }
        }
        #Compute F to Project at SPR 40 percent
        if(FishingMortaltiyProjectionScenarios[v]=="F_at_SPR40")
        {
          if(Scenario_i$populationReconstructionMethod=="LBSPR")
          {
            MyPars = new("LB_pars")
            MyPars@Species = Scenario_i$Species
            MyPars@Linf = Scenario_i$Linf
            MyPars@L50 = Scenario_i$Lmat #length at 50% maturity
            MyPars@L95 = Species_i$lopt #length at 95% maturity - we don't know this so use lopt as a proxy. We don't know this because we are assuming knife edge maturity and don't have a study to know what the real maturity schedule looks like
            MyPars@MK = Scenario_i$M/Scenario_i$VBK
            MyPars@L_units = "cm"
            MyPars@BinMin=min(Scenario_i$NumAtLength$Midpoint)-2.5
            MyPars@BinWidth=Scenario_i$NumAtLength$Midpoint[2]-Scenario_i$NumAtLength$Midpoint[1]
            MyPars@BinMax=ceiling(Species_i$lmax/5)*5
            Len1 = new("LB_lengths",LB_pars=MyPars,file=Sub$cm_converted,sataType="raw")
            myFit1 = LBSPRfit(MyPars, Len1)
            myFit1@Ests
            plotSize(myFit1)
            plotMat(myFit1)
            #simulate population
            SimPars=MyPars
            SimPars@SL50 = Scenario_i$Lmat #length at 50% maturity
            SimPars@SL95 = Species_i$lopt #length at 95% maturity - we don't know this so use lopt as a proxy. We don't know this because we are assuming knife edge maturity and don't have a study to know what the real maturity schedule looks like
            SimPars@SPR = ControlRule40_SPR
            SimPars@R0=Scenario_i$R0
            MySim = LBSPRsim(SimPars,Control = list("modtype"="absel"))
            F_project = MySim@FM*Scenario_i$M
          }
          if(Scenario_i$populationReconstructionMethod=="SolveBaranov")
          {
            NumAtAge_thisScenario = subset(NumAtAge,select=-c(Population_Extrapolated,EstimatedCatch_Baranov,Residuals,SSB_CurrentAtAge))
            EstNumFish_Population = sum(NumAtAge_thisScenario$PopulationVirgin)*ControlRule_SPR
            F_start = 1.0
            min.RSS = function(F_start)
            {
              NumAtAge_thisScenario$Survival_F_M = exp(-1*(Scenario_i$M+(F_start*NumAtAge_thisScenario$Selex_Empirical))*NumAtAge_thisScenario$Age)
              NumAtAge_thisScenario$Survival_F_M_Rescaled = NumAtAge_thisScenario$Survival_F_M/sum(NumAtAge_thisScenario$Survival_F_M)
              NumAtAge_thisScenario$Population_Extrapolated = NumAtAge_thisScenario$Survival_F_M_Rescaled*EstNumFish_Population
              NumAtAge_thisScenario$SSB_CurrentAtAge = NumAtAge_thisScenario$Population_Extrapolated*(NumAtAge_thisScenario$Weight_kg*0.001)*NumAtAge_thisScenario$Maturity*NumAtAge_thisScenario$Fecundity
              SSB_current = sum(NumAtAge_thisScenario$SSB_CurrentAtAge)
              SPR_current = SSB_current/Scenario_i$SSB_virgin
              return(sum((SPR_current - ControlRule40_SPR)^2))
            }
            Model = optimize(min.RSS,c(0,10),tol=0.001)
            NumAtAge_thisScenario$Survival_F_M = exp(-1*(Scenario_i$M+(Model$minimum[1]*NumAtAge_thisScenario$Selex_Empirical))*NumAtAge_thisScenario$Age)     #"pseudosurvival"
            NumAtAge_thisScenario$Survival_F_M_Rescaled = NumAtAge_thisScenario$Survival_F_M/sum(NumAtAge_thisScenario$Survival_F_M)
            NumAtAge_thisScenario$Population_Extrapolated = NumAtAge_thisScenario$Survival_F_M_Rescaled*EstNumFish_Population
            NumAtAge_thisScenario$Survival_F_M=NULL
            NumAtAge_thisScenario$Survival_F_M_Rescaled=NULL
            F_start = 1.0
            min.RSS = function(F_start)
            {
              EstimatedCatch_Baranov = BaranovCatch(NumAtAge_thisScenario$Age,NumAtAge_thisScenario$Length,NumAtAge_thisScenario$Population_Extrapolated,Scenario_i$M,F_start,NumAtAge_thisScenario$Selex_Empirical)
              return(sum((EstimatedCatch_Baranov - NumAtAge_thisScenario$ExtrapolatedCatchInNumber)^2))
            }
            Model = optimize(min.RSS,c(0,10),tol=0.001)
            F_project = Model$minimum[1]
          }
          if(Scenario_i$populationReconstructionMethod=="BevertonHoltInstantaneousMortality")
          {
            F_start = 1.0
            min.RSS = function(F_start)
            {
              LenVec = Sub$cm_converted
              LenVec = LenVec[!is.na(LenVec)]
              SPR = computeLengthSPR(Scenario_i$M,Scenario_i$Linf,Scenario_i$VBK,Species_i$var_a,Species_i$var_b,Species_i$lmat,F_start,LenVec)
              return(sum((SPR - ControlRule40_SPR)^2))
            }
            Model = optimize(min.RSS,c(0,10),tol=0.001)
            F_project = Model$minimum[1]
          }
          if(Scenario_i$populationReconstructionMethod=="LIME")
          {
            Len50_selex = 30
            min.RSS = function(Len50_selex)
            {
              Sel = SelectivityLogistic_Length(Len50_selex,Scenario_i$Selex_Logistic_a_length,Scenario_i$Selex_Logistic_b_length)
              return(sum((Sel - 0.5)^2))
            }
            Model = optimize(min.RSS,c(0,200),tol=0.001)
            Len50_selex = Model$minimum[1]
            Len95_selex = 50
            min.RSS = function(Len95_selex)
            {
              Sel = SelectivityLogistic_Length(Len95_selex,Scenario_i$Selex_Logistic_a_length,Scenario_i$Selex_Logistic_b_length)
              return(sum((Sel - 0.95)^2))
            }
            Model = optimize(min.RSS,c(0,200),tol=0.001)
            Len95_selex = Model$minimum[1]
            lh = create_lh_list(vbk=Scenario_i$VBK,linf=Scenario_i$Linf,lwa=Species_i$var_a,lwb=Species_i$var_b,
                S50=Len50_selex,S95=Len95_selex,selex_input="length",selex_type=c("logistic"), #Starting selectivity values and function specified
                M50=Species_i$lmat,maturity_input="length",M=Scenario_i$M,binwidth=5,R0=Scenario_i$R0, #R0 is starting value - will be found by LIME
                rho=0,h=0.8,  #steepness (h) and cruitment autocorrelation (rho)
                nseasons=1,nfleets=1,CVlen=0.2,SigmaR=0.1,AgeMax=max(NumAtAge$Age)) #if we could get MONTHLY effort, and split catches per fleet that would go here
            LF=matrix(Scenario_i$NumAtLength$NumSampled,nrow=1)
            colnames(LF)=Scenario_i$NumAtLength$Length_max
            rownames(LF)=max(SpeciesCatch_mT_Scaled$Year)
            Catch_matrix = as.matrix(sum(Scenario_i$NumAtAge$ExtrapolatedCatchInNumber*Scenario_i$NumAtAge$Weight_kg*1000))
            Catch_matrix = Catch_matrix[,length(Catch_matrix)]
            data_all=list(years=max(SpeciesCatch_mT_Scaled$Year),LF=LF,C_ft=as.matrix(Catch_matrix))   #list of years, length distribution (as a matrix of one row), catches and number of observed catches
            inputs_all=create_inputs(lh=lh,input_data=data_all)
            vals_selex_ft = matrix(Scenario_i$NumAtLength$Selex_Fitted,nrow=1)
            Results_LIME = run_LIME(modpath=NULL,input=inputs_all,data_avail="Catch_LC1",C_type=2,est_selex_f=FALSE,vals_selex_ft=vals_selex_ft,derive_quants=TRUE,est_totalF=TRUE,Fpen=0)       #C_type is catch in weight, not numbers
            gradient = Results_LIME$opt$max_gradient <= 0.001
            hessian = Results_LIME$Sdreport$pdHess
            Converge = (hessian==TRUE & gradient == TRUE)
            print(paste("gradient should be below 0.001, it is",Results_LIME$opt$max_gradient))
            print(paste("does the hessian test pass? ",Results_LIME$Sdreport$pdHess))
            print(paste("SPR estimated as ",Results_LIME$Report$SPR_t,sep=""))
            flush.console()
            F_project = Results_LIME$Derived$F40
          }
          if(Scenario_i$populationReconstructionMethod=="TropFishR")
          {
            SPR_current = sum(Scenario_i$NumAtAge$Survival_M_F_rescaled*(Scenario_i$NumAtAge$Weight_kg*0.001)*Scenario_i$NumAtAge$Maturity*Scenario_i$NumAtAge$Fecundity)/sum(Scenario_i$NumAtAge$Survival_M_rescaled*(Scenario_i$NumAtAge$Weight_kg*0.001)*Scenario_i$NumAtAge$Maturity*Scenario_i$NumAtAge$Fecundity)
            F_start=0.5
            min.RSS = function(F_start)
            {
              Survival_F_M = exp(-1*(Scenario_i$M+(Scenario_i$NumAtAge$Selex_Empirical*F_start))*Scenario_i$NumAtAge$Age)
              Survival_F_M_rescaled = Survival_F_M/sum(Survival_F_M)
              SPR = sum(Survival_F_M_rescaled*(Scenario_i$NumAtAge$Weight_kg*0.001)*Scenario_i$NumAtAge$Maturity*Scenario_i$NumAtAge$Fecundity)/sum(Scenario_i$NumAtAge$Survival_M_rescaled*(Scenario_i$NumAtAge$Weight_kg*0.001)*Scenario_i$NumAtAge$Maturity*Scenario_i$NumAtAge$Fecundity)
              return(sum((SPR - ControlRule40_SPR)^2))
            }
            Model = optimize(min.RSS,c(0,10),tol=0.001)
            F_project = Model$minimum[1]
          }
        }
        #Projections Across Different Steepness Values for those scenarios with Bev. Holt. recruitment
        if(FishingMortaltiyProjectionScenarios[v]!="ConstantCatch")
        {
          for(j in 1:length(Steepness_Scenarios))
          {
            #F-based
            Projection = projectionFunction(Scenario_i,yearsToProject,yearToStartProjection,F_project,Steepness_Scenarios[j])
            Projection$FishingMortaltiyProjectionScenario = FishingMortaltiyProjectionScenarios[v]
            Projection$SteepnessAssumed = Steepness_Scenarios[j]
            Projection$ProjectionType = "RegularF"
            Projection$F_project = F_project
            ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection
          }
          #F Based Projections with Linear/Constant Recruitment
          Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,yearsToProject,yearToStartProjection,F_project)
          Projection_LinearRecruitment$FishingMortaltiyProjectionScenario = FishingMortaltiyProjectionScenarios[v]
          Projection_LinearRecruitment$SteepnessAssumed = NA
          Projection_LinearRecruitment$ProjectionType = "ConstantLinearRecruitment"
          Projection_LinearRecruitment$F_project = F_project
          ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection_LinearRecruitment
        }
        if(FishingMortaltiyProjectionScenarios[v]=="ConstantCatch")
        {
          for(j in 1:length(Steepness_Scenarios))
          {
            #Constant Catch
            Projection_ConstantCatch = projectionFunction_constantCatch(Scenario_i,yearsToProject,yearToStartProjection,Steepness_Scenarios[j])
            Projection_ConstantCatch$FishingMortaltiyProjectionScenario = FishingMortaltiyProjectionScenarios[v]
            Projection_ConstantCatch$SteepnessAssumed = Steepness_Scenarios[j]
            Projection_ConstantCatch$ProjectionType = "ConstantCatch"
            Projection_ConstantCatch$F_project = NA
            ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection_ConstantCatch
          }
          #Projection with Constant Catch and Linear/Constant Recruitment
          Projection_ConstantCatch_LinearRecruitment = projectionFunction_LinearRecruitment_ConstantCatch(Scenario_i,yearsToProject,yearToStartProjection)
          Projection_ConstantCatch_LinearRecruitment$FishingMortaltiyProjectionScenario = FishingMortaltiyProjectionScenarios[v]
          Projection_ConstantCatch_LinearRecruitment$SteepnessAssumed = NA
          Projection_ConstantCatch_LinearRecruitment$ProjectionType = "ConstantCatchAndConstantLinearRecruitment"
          Projection_ConstantCatch_LinearRecruitment$F_project = NA
          ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection_ConstantCatch_LinearRecruitment
        }
      }
      #If the population is overfished (i.e. SPR less than 30%) and we need to rebuild, then find the F value to rebuild in the timeframe, and do the rebuilding scenarios
      if(FishingMortaltiyProjectionScenarios[v]=="F_rebuild_20YrsMSC_max" | FishingMortaltiyProjectionScenarios[v]=="F_rebuild_Tmin_oneGeneration" | FishingMortaltiyProjectionScenarios[v]=="F_rebuild_10YrsUS")
      {
        SPR_current_thisScenario = sum(NumAtAge$SSB_CurrentAtAge)/Scenario_i$SSB_virgin
        if(runRebuildingProjectionsRegardlessOfStatus==TRUE | (runRebuildingProjectionsRegardlessOfStatus==FALSE & SPR_current_thisScenario>=ControlRule_SPR))
        {
          if(FishingMortaltiyProjectionScenarios[v]=="F_rebuild_20YrsMSC_max")
          {
            for(j in 1:length(Steepness_Scenarios))
            {
              #(1)Find F-Rebuild and Project Using the F Based Projections
              projectYrs = 20
              F_project=0
              Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
              SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
              if(SPR_lastYear<=ControlRule_SPR)
              {
                Projection$FishingMortaltiyProjectionScenario = FishingMortaltiyProjectionScenarios[v]
                Projection$SteepnessAssumed = Steepness_Scenarios[j]
                Projection$ProjectionType = "RegularF_willNotRebuildAtF=0"
                Projection$F_project=F_project
                ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection
              }
              if(SPR_lastYear>ControlRule_SPR)
              {
                 while(SPR_lastYear>ControlRule_SPR)
                 {
                    F_project = F_project + 0.1
                    Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                    SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 }
                 F_project = F_project - 0.1
                 if(F_project<0) {F_project=0}
                 Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                 SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 while(SPR_lastYear>ControlRule_SPR)
                 {
                    F_project = F_project + 0.01
                    Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                    SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 }
                 F_project = F_project - 0.01
                 if(F_project<0) {F_project=0}
                 Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                 SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 while(SPR_lastYear>ControlRule_SPR)
                 {
                    F_project = F_project + 0.001
                    Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                    SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 }
                 F_project = F_project - 0.001
                 if(F_project<0) {F_project=0}
                 Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                 SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 while(SPR_lastYear>ControlRule_SPR)
                 {
                    F_project = F_project + 0.0001
                    Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                    SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 }
                 Projection = projectionFunction(Scenario_i,yearsToProject,yearToStartProjection,F_project,Steepness_Scenarios[j])
                 Projection$FishingMortaltiyProjectionScenario = FishingMortaltiyProjectionScenarios[v]
                 Projection$SteepnessAssumed = Steepness_Scenarios[j]
                 Projection$ProjectionType = "RegularF"
                 Projection$F_project = F_project
                 ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection
              }
            }
            #(2)Find F-Rebuild and Project Using F Based Projections with Linear/Constant Recruitment
            projectYrs = 20
            F_project=0
            StartValues_species = subset(StartValues,select=c(Genus,Species,Fecundity_min,Fecundity_max))
            StartValues_species = unique(StartValues_species)
            Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
            SSB_lastYear = sum(Projection_LinearRecruitment$ProjectedPopulationAtAge[length(Projection_LinearRecruitment$ProjectedPopulationAtAge)]*(NumAtAge$Weight_kg*0.001)*NumAtAge$Maturity*NumAtAge$Fecundity)
            SPR_lastYear = SSB_lastYear/Projection_LinearRecruitment$SSB_virgin
            if(SPR_lastYear<=ControlRule_SPR)
            {
              Projection_LinearRecruitment$FishingMortaltiyProjectionScenario = FishingMortaltiyProjectionScenarios[v]
              Projection_LinearRecruitment$SteepnessAssumed = Steepness_Scenarios[j]
              Projection_LinearRecruitment$ProjectionType = "ConstantLinearRecruitment_willNotRebuildAtF=0"
              Projection_LinearRecruitment$F_project=F_project
              ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection_LinearRecruitment
            }
            if(SPR_lastYear>ControlRule_SPR)
            {
               while(SPR_lastYear>ControlRule_SPR)
               {
                  F_project = F_project + 0.1
                  Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
                  SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               }
               F_project = F_project - 0.1
               if(F_project<0) {F_project=0}
               Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
               SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               while(SPR_lastYear>ControlRule_SPR)
               {
                  F_project = F_project + 0.01
                  Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
                  SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               }
               F_project = F_project - 0.01
               if(F_project<0) {F_project=0}
               Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
               SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               while(SPR_lastYear>ControlRule_SPR)
               {
                  F_project = F_project + 0.001
                  Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
                  SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               }
               F_project = F_project - 0.001
               if(F_project<0) {F_project=0}
               Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
               SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               while(SPR_lastYear>ControlRule_SPR)
               {
                  F_project = F_project + 0.0001
                  Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
                  SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               }
               Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,yearsToProject,yearToStartProjection,F_project)
               Projection_LinearRecruitment$FishingMortaltiyProjectionScenario = FishingMortaltiyProjectionScenarios[v]
               Projection_LinearRecruitment$SteepnessAssumed = Steepness_Scenarios[j]
               Projection_LinearRecruitment$ProjectionType = "ConstantLinearRecruitment"
               Projection_LinearRecruitment$F_project=F_project
               ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection_LinearRecruitment
            }
          }
          if(FishingMortaltiyProjectionScenarios[v]=="F_rebuild_Tmin_oneGeneration")
          {
            oneGenTime_numerator = Scenario_i$NumAtAge$Survival_M*Scenario_i$NumAtAge$Age*Scenario_i$NumAtAge$Maturity
            oneGenTime_denominator = Scenario_i$NumAtAge$Survival_M*Scenario_i$NumAtAge$Maturity
            oneGenerationTime = ceiling(sum(oneGenTime_numerator)/sum(oneGenTime_denominator))
            for(j in 1:length(Steepness_Scenarios))
            {
                YearCounter=1
                SPR_lastYear=SPR_current_thisScenario
                F_project=0
                while(SPR_lastYear<ControlRule_SPR)
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
                    SSB_m = sum(Population_end*(NumAtAge$Weight_kg*0.001)*NumAtAge$Maturity*NumAtAge$Fecundity)
                    Recruitment_m =  BevHoltRecruitment(Scenario_i$R0,Steepness_Scenarios[j],SSB_m,Scenario_i$SSB_virgin)
                    Population_end[1] = Recruitment_m
                    SPR_lastYear = SSB_m/Scenario_i$SSB_virgin
                    YearCounter = YearCounter + 1
                 }
              }
              timeToRebuld_noFishing = YearCounter
              #Find the F value that rebuilds the stock in the time period we just calculated.
              projectYrs = oneGenerationTime + timeToRebuld_noFishing
              F_project=0
              Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
              SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
              if(SPR_lastYear<=ControlRule_SPR)
              {
                Projection$FishingMortaltiyProjectionScenario = FishingMortaltiyProjectionScenarios[v]
                Projection$SteepnessAssumed = Steepness_Scenarios[j]
                Projection$ProjectionType = "RegularF_willNotRebuildAtF=0"
                Projection$F_project = F_project
                ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection
              }
              if(SPR_lastYear>ControlRule_SPR)
              {
                 while(SPR_lastYear>ControlRule_SPR)
                 {
                    F_project = F_project + 0.1
                    Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                    SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 }
                 F_project = F_project - 0.1
                 if(F_project<0) {F_project=0}
                 Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                 SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 while(SPR_lastYear>ControlRule_SPR)
                 {
                    F_project = F_project + 0.01
                    Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                    SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 }
                 F_project = F_project - 0.01
                 if(F_project<0) {F_project=0}
                 Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                 SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 while(SPR_lastYear>ControlRule_SPR)
                 {
                    F_project = F_project + 0.001
                    Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                    SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 }
                 F_project = F_project - 0.001
                 if(F_project<0) {F_project=0}
                 Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                 SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 while(SPR_lastYear>ControlRule_SPR)
                 {
                    F_project = F_project + 0.0001
                    Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                    SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 }
                 Projection = projectionFunction(Scenario_i,yearsToProject,yearToStartProjection,F_project,Steepness_Scenarios[j])
                 Projection$FishingMortaltiyProjectionScenario = FishingMortaltiyProjectionScenarios[v]
                 Projection$SteepnessAssumed = Steepness_Scenarios[j]
                 Projection$ProjectionType = "RegularF"
                 Projection$F_project = F_project
                 Projection$timeToRebuild_noFishing = YearCounter
                 Projection$OneGenerationTime = oneGenerationTime
                 Projection$ProjectionTime = projectYrs
                 ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection
              }
            }
            #(2)Find F-Rebuild and Project Using F Based Projections with Linear/Constant Recruitment
            YearCounter=1
            SPR_lastYear=SPR_current_thisScenario
            F_project=0
            ConstantRecruitment = Scenario_i$R0
            LinearModel_recruitBetween0_20percent = lm(c(0,ConstantRecruitment) ~ c(0,0.2))
            while(SPR_lastYear<ControlRule_SPR)
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
                SSB_m = sum(NumAtAge$SSB_CurrentAtAge)
                SPR = SSB_m/Scenario_i$SSB_virgin
                if(SPR>=0.2)
                {
                  Recruitment_m =  ConstantRecruitment
                }
                if(SPR<0.2)
                {
                  Recruitment_m = (SPR*LinearModel_recruitBetween0_20percent$coef[2]) + LinearModel_recruitBetween0_20percent$coef[1]
                }
                Population_start[1] = Recruitment_m
                Catch_m = ((F_project*NumAtAge$Selex_Fitted)/(Scenario_i$M+(F_project*NumAtAge$Selex_Fitted)))*Population_start*(1-exp(-(Scenario_i$M+(F_project*NumAtAge$Selex_Fitted))))
                NaturalDeaths_m = (Scenario_i$M/(Scenario_i$M+(F_project*NumAtAge$Selex_Fitted)))*Population_start*(1-exp(-(Scenario_i$M+(F_project*NumAtAge$Selex_Fitted))))
                Population_end = Population_start - (Catch_m + NaturalDeaths_m)
                Population_end[Population_end<0]=0   #can't have negative numbers of fish
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
                SSB_m = sum(Population_end*(NumAtAge$Weight_kg*0.001)*NumAtAge$Maturity*NumAtAge$Fecundity)
                SPR = SSB_m/Scenario_i$SSB_virgin
                if(SPR>=0.2)
                {
                  Recruitment_m =  ConstantRecruitment
                }
                if(SPR<0.2)
                {
                  Recruitment_m = (SPR*LinearModel_recruitBetween0_20percent$coef[2]) + LinearModel_recruitBetween0_20percent$coef[1]
                }
                Population_start[1] = Recruitment_m
                Catch_m = ((F_project*NumAtAge$Selex_Fitted)/(Scenario_i$M+(F_project*NumAtAge$Selex_Fitted)))*Population_start*(1-exp(-(Scenario_i$M+(F_project*NumAtAge$Selex_Fitted))))
                NaturalDeaths_m = (Scenario_i$M/(Scenario_i$M+(F_project*NumAtAge$Selex_Fitted)))*Population_start*(1-exp(-(Scenario_i$M+(F_project*NumAtAge$Selex_Fitted))))
                Population_end = Population_start - (Catch_m + NaturalDeaths_m)
                Population_end[Population_end<0]=0   #can't have negative numbers of fish
                SPR_lastYear = SSB_m/Scenario_i$SSB_virgin
                YearCounter = YearCounter + 1
              }
            }
            timeToRebuld_noFishing = YearCounter
            #(1)Find F to rebuild in the time period we just calculated
            projectYrs = oneGenerationTime + timeToRebuld_noFishing
            F_project=0
            Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
            SSB_lastYear = sum(Projection_LinearRecruitment$ProjectedPopulationAtAge[length(Projection_LinearRecruitment$ProjectedPopulationAtAge)]*(NumAtAge$Weight_kg*0.001)*NumAtAge$Maturity*NumAtAge$Fecundity)
            SPR_lastYear = SSB_lastYear/Projection_LinearRecruitment$SSB_virgin
            if(SPR_lastYear<=ControlRule_SPR)
            {
              Projection_LinearRecruitment$FishingMortaltiyProjectionScenario = FishingMortaltiyProjectionScenarios[v]
              Projection_LinearRecruitment$SteepnessAssumed = Steepness_Scenarios[j]
              Projection_LinearRecruitment$ProjectionType = "ConstantLinearRecruitment_willNotRebuildAtF=0"
              Projection_LinearRecruitment$F_project = F_project
              ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection_LinearRecruitment
            }
            if(SPR_lastYear>ControlRule_SPR)
            {
               while(SPR_lastYear>ControlRule_SPR)
               {
                  F_project = F_project + 0.1
                  Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
                  SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               }
               F_project = F_project - 0.1
               if(F_project<0) {F_project=0}
               Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
               SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               while(SPR_lastYear>ControlRule_SPR)
               {
                  F_project = F_project + 0.01
                  Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
                  SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               }
               F_project = F_project - 0.01
               if(F_project<0) {F_project=0}
               Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
               SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               while(SPR_lastYear>ControlRule_SPR)
               {
                  F_project = F_project + 0.001
                  Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
                  SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               }
               F_project = F_project - 0.001
               if(F_project<0) {F_project=0}
               Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
               SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               while(SPR_lastYear>ControlRule_SPR)
               {
                  F_project = F_project + 0.0001
                  Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
                  SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               }
               Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,yearsToProject,yearToStartProjection,F_project)
               Projection_LinearRecruitment$FishingMortaltiyProjectionScenario = FishingMortaltiyProjectionScenarios[v]
               Projection_LinearRecruitment$SteepnessAssumed = NA
               Projection_LinearRecruitment$ProjectionType = "ConstantLinearRecruitment"
               Projection_LinearRecruitment$F_project = F_project
               Projection_LinearRecruitment$timeToRebuild_noFishing = YearCounter
               Projection_LinearRecruitment$OneGenerationTime = oneGenerationTime
               Projection_LinearRecruitment$ProjectionTime = projectYrs
               ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection_LinearRecruitment
            }
          }
          if(FishingMortaltiyProjectionScenarios[v]=="F_rebuild_10YrsUS")
          {
            for(j in 1:length(Steepness_Scenarios))
            {
              #(1)Find F-Rebuild and Project Using the F Based Projections
              projectYrs = 10
              F_project=0
              Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
              SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
              if(SPR_lastYear<=ControlRule_SPR)
              {
                Projection$FishingMortaltiyProjectionScenario = FishingMortaltiyProjectionScenarios[v]
                Projection$SteepnessAssumed = Steepness_Scenarios[j]
                Projection$ProjectionType = "RegularF_willNotRebuildAtF=0"
                Projection$F_project = F_project
                ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection
              }
              if(SPR_lastYear>ControlRule_SPR)
              {
                 while(SPR_lastYear>ControlRule_SPR)
                 {
                    F_project = F_project + 0.1
                    Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                    SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 }
                 F_project = F_project - 0.1
                 if(F_project<0) {F_project=0}
                 Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                 SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 while(SPR_lastYear>ControlRule_SPR)
                 {
                    F_project = F_project + 0.01
                    Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                    SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 }
                 F_project = F_project - 0.01
                 if(F_project<0) {F_project=0}
                 Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                 SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 while(SPR_lastYear>ControlRule_SPR)
                 {
                    F_project = F_project + 0.001
                    Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                    SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 }
                 F_project = F_project - 0.001
                 if(F_project<0) {F_project=0}
                 Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                 SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 while(SPR_lastYear>ControlRule_SPR)
                 {
                    F_project = F_project + 0.0001
                    Projection = projectionFunction(Scenario_i,projectYrs,yearToStartProjection,F_project,Steepness_Scenarios[j])
                    SPR_lastYear = Projection$SSB[length(Projection$SSB)]/Projection$SSB_virgin
                 }
                 Projection = projectionFunction(Scenario_i,yearsToProject,yearToStartProjection,F_project,Steepness_Scenarios[j])
                 Projection$FishingMortaltiyProjectionScenario = FishingMortaltiyProjectionScenarios[v]
                 Projection$SteepnessAssumed = Steepness_Scenarios[j]
                 Projection$ProjectionType = "RegularF"
                 Projection$F_project = F_project
                 ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection
              }
            }
            #(2)Find F-Rebuild and Project Using F Based Projections with Linear/Constant Recruitment
            projectYrs = 10
            F_project=0
            Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
            SSB_lastYear = sum(Projection_LinearRecruitment$ProjectedPopulationAtAge[length(Projection_LinearRecruitment$ProjectedPopulationAtAge)]*(NumAtAge$Weight_kg*0.001)*NumAtAge$Maturity*NumAtAge$Fecundity)
            SPR_lastYear = SSB_lastYear/Projection_LinearRecruitment$SSB_virgin
            if(SPR_lastYear<=ControlRule_SPR)
            {
              Projection_LinearRecruitment$FishingMortaltiyProjectionScenario = FishingMortaltiyProjectionScenarios[v]
              Projection_LinearRecruitment$SteepnessAssumed = NA
              Projection_LinearRecruitment$ProjectionType = "ConstantLinearRecruitment_willNotRebuildAtF=0"
              Projection_LinearRecruitment$F_project = F_project
              ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection_LinearRecruitment
            }
            if(SPR_lastYear>ControlRule_SPR)
            {
               while(SPR_lastYear>ControlRule_SPR)
               {
                  F_project = F_project + 0.1
                  Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
                  SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               }
               F_project = F_project - 0.1
               if(F_project<0) {F_project=0}
               Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
               SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               while(SPR_lastYear>ControlRule_SPR)
               {
                  F_project = F_project + 0.01
                  Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
                  SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               }
               F_project = F_project - 0.01
               if(F_project<0) {F_project=0}
               Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
               SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               while(SPR_lastYear>ControlRule_SPR)
               {
                  F_project = F_project + 0.001
                  Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
                  SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               }
               F_project = F_project - 0.001
               if(F_project<0) {F_project=0}
               Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
               SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               while(SPR_lastYear>ControlRule_SPR)
               {
                  F_project = F_project + 0.0001
                  Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,projectYrs,yearToStartProjection,F_project)
                  SPR_lastYear = Projection_LinearRecruitment$SSB[length(Projection_LinearRecruitment$SSB)]/Projection_LinearRecruitment$SSB_virgin
               }
               Projection_LinearRecruitment = projectionFunction_LinearRecruitment(Scenario_i,yearsToProject,yearToStartProjection,F_project)
               Projection_LinearRecruitment$FishingMortaltiyProjectionScenario = FishingMortaltiyProjectionScenarios[v]
               Projection_LinearRecruitment$SteepnessAssumed = NA
               Projection_LinearRecruitment$ProjectionType = "ConstantLinearRecruitment"
               Projection_LinearRecruitment$F_project = F_project
               ListOfPopulations_WithProjections[[length(ListOfPopulations_WithProjections)+1]] = Projection_LinearRecruitment
            }
          }
        }
      }
    }
    if(i%%50==0)
    {
      print(paste("Projecting population scenario ",i," of ",length(ListOfPopulations),"...",sep=""))
      flush.console()
    }
  }
  saveRDS(ListOfPopulations_WithProjections,file=paste(PATH_output,"ListOfPopulations_WithProjections.RData",sep="/"),ascii=FALSE,version=NULL,compress=TRUE,refhook=NULL)

  #Read the list of populations back in and interate through it to pull out key values and produce the look-up table file
  ListOfPopulations_WithProjections = readRDS(file=paste(PATH_output,"ListOfPopulations_WithProjections.RData",sep="/"),refhook=NULL)  ########################## Read in Projection RDS File and Create Summary #########################
  SummaryProjectionScenarios = data.frame(Genus=rep("",times=length(ListOfPopulations_WithProjections)),Species=rep("",times=length(ListOfPopulations_WithProjections)),
    WPP=rep("",times=length(ListOfPopulations_WithProjections)),CatchMSY_Scenario=rep("",times=length(ListOfPopulations_WithProjections)),CatchMSY_Scenario_Bound=rep("",times=length(ListOfPopulations_WithProjections)),
    LifeHistory_Scenario=rep("",times=length(ListOfPopulations_WithProjections)),populationReconstructionMethod=rep("",times=length(ListOfPopulations_WithProjections)),F_estimated=rep(0,times=length(ListOfPopulations_WithProjections)),
    SSB_estimated=rep(0,times=length(ListOfPopulations_WithProjections)),Yield_estimated=rep(0,times=length(ListOfPopulations_WithProjections)),N_estimated=rep(0,times=length(ListOfPopulations_WithProjections)),
    B_estimated=rep(0,times=length(ListOfPopulations_WithProjections)),Recruits_estimated=rep(0,times=length(ListOfPopulations_WithProjections)),R0=rep(0,times=length(ListOfPopulations_WithProjections)),
    Selex_Logistic_a_length=rep(0,times=length(ListOfPopulations_WithProjections)),Selex_Logistic_b_length=rep(0,times=length(ListOfPopulations_WithProjections)),Selex_Logistic_a_age=rep(0,times=length(ListOfPopulations_WithProjections)),
    Selex_Logistic_b_age=rep(0,times=length(ListOfPopulations_WithProjections)),Mnat=rep(0,times=length(ListOfPopulations_WithProjections)),ListIndexNumber=rep(0,times=length(ListOfPopulations_WithProjections)),
    FishingMortaltiyProjectionScenario=rep("",times=length(ListOfPopulations_WithProjections)),SteepnessAssumed=rep(0,times=length(ListOfPopulations_WithProjections)),ProjectionType=rep("",times=length(ListOfPopulations_WithProjections)),
    F_project=rep(0,times=length(ListOfPopulations_WithProjections)),NumFishAtEnd=rep(0,times=length(ListOfPopulations_WithProjections)),BiomassAtEnd=rep(0,times=length(ListOfPopulations_WithProjections)),
    RecruitmentAtEnd=rep(0,times=length(ListOfPopulations_WithProjections)),SSBatEnd=rep(0,times=length(ListOfPopulations_WithProjections)),CatchAtEnd=rep(0,times=length(ListOfPopulations_WithProjections)),Tmin=rep(NA,times=length(ListOfPopulations_WithProjections)),OneGenerationTime=rep(NA,times=length(ListOfPopulations_WithProjections)))
    i = sapply(SummaryProjectionScenarios, is.factor)
    SummaryProjectionScenarios[i] = lapply(SummaryProjectionScenarios[i], as.character)
  for(i in 1:length(ListOfPopulations_WithProjections))
  {
    SummaryProjectionScenarios$Genus[i]=ListOfPopulations_WithProjections[[i]]$Genus
    SummaryProjectionScenarios$Species[i]=ListOfPopulations_WithProjections[[i]]$Species
    SummaryProjectionScenarios$WPP[i]=ListOfPopulations_WithProjections[[i]]$WPP
    SummaryProjectionScenarios$CatchMSY_Scenario[i]=ListOfPopulations_WithProjections[[i]]$CatchMSY_Scenario
    SummaryProjectionScenarios$CatchMSY_Scenario_Bound[i]=ListOfPopulations_WithProjections[[i]]$BoundFromCatchMSY
    SummaryProjectionScenarios$LifeHistory_Scenario[i]=ListOfPopulations_WithProjections[[i]]$LifeHistory_Scenario
    SummaryProjectionScenarios$populationReconstructionMethod[i]=ListOfPopulations_WithProjections[[i]]$populationReconstructionMethod
    SummaryProjectionScenarios$F_estimated[i]=ListOfPopulations_WithProjections[[i]]$F_estimated
    SummaryProjectionScenarios$SSB_estimated[i]=sum(ListOfPopulations_WithProjections[[i]]$NumAtAge$SSB_CurrentAtAge)
    SummaryProjectionScenarios$Yield_estimated[i]=sum(ListOfPopulations_WithProjections[[i]]$NumAtAge$ExtrapolatedCatchInNumber*ListOfPopulations_WithProjections[[i]]$NumAtAge$Weight_kg)
    SummaryProjectionScenarios$N_estimated[i]=sum(ListOfPopulations_WithProjections[[i]]$NumAtAge$Population_Extrapolated)
    SummaryProjectionScenarios$B_estimated[i]=sum(ListOfPopulations_WithProjections[[i]]$NumAtAge$Population_Extrapolated*ListOfPopulations_WithProjections[[i]]$NumAtAge$Weight_kg)
    SummaryProjectionScenarios$Recruits_estimated[i]=ListOfPopulations_WithProjections[[i]]$Recruitment[1]
    SummaryProjectionScenarios$R0[i]=ListOfPopulations_WithProjections[[i]]$R0
    SummaryProjectionScenarios$Selex_Logistic_a_length[i]=ListOfPopulations_WithProjections[[i]]$Selex_Logistic_a_length
    SummaryProjectionScenarios$Selex_Logistic_b_length[i]=ListOfPopulations_WithProjections[[i]]$Selex_Logistic_b_length
    SummaryProjectionScenarios$Selex_Logistic_a_age[i]=ListOfPopulations_WithProjections[[i]]$Selex_Logistic_a_age
    SummaryProjectionScenarios$Selex_Logistic_b_age[i]=ListOfPopulations_WithProjections[[i]]$Selex_Logistic_b_age
    SummaryProjectionScenarios$Mnat[i]=ListOfPopulations_WithProjections[[i]]$M
    SummaryProjectionScenarios$ListIndexNumber[i]=i
    SummaryProjectionScenarios$FishingMortaltiyProjectionScenario[i]=ListOfPopulations_WithProjections[[i]]$FishingMortaltiyProjectionScenario
    if(ListOfPopulations_WithProjections[[i]]$ProjectionType=="ConstantLinearRecruitment" | ListOfPopulations_WithProjections[[i]]$ProjectionType=="ConstantCatchAndConstantLinearRecruitment")
    {
        SummaryProjectionScenarios$SteepnessAssumed[i]=NA
    } else  {
        SummaryProjectionScenarios$SteepnessAssumed[i]=ListOfPopulations_WithProjections[[i]]$SteepnessAssumed
    }
    SummaryProjectionScenarios$ProjectionType[i]=ListOfPopulations_WithProjections[[i]]$ProjectionType
    SummaryProjectionScenarios$F_project[i]=ListOfPopulations_WithProjections[[i]]$F_project
    SummaryProjectionScenarios$NumFishAtEnd[i]=ListOfPopulations_WithProjections[[i]]$ProjectedPopulationYear[length(ListOfPopulations_WithProjections[[i]]$ProjectedPopulationYear)]
    SummaryProjectionScenarios$BiomassAtEnd[i]=ListOfPopulations_WithProjections[[i]]$ProjectedBiomassYear[length(ListOfPopulations_WithProjections[[i]]$ProjectedBiomassYear)]
    SummaryProjectionScenarios$RecruitmentAtEnd[i]=ListOfPopulations_WithProjections[[i]]$Recruitment[length(ListOfPopulations_WithProjections[[i]]$Recruitment)]
    SummaryProjectionScenarios$SSBatEnd[i]=ListOfPopulations_WithProjections[[i]]$SSB[length(ListOfPopulations_WithProjections[[i]]$SSB)]
    SummaryProjectionScenarios$CatchAtEnd[i]=sum(ListOfPopulations_WithProjections[[i]]$ProjectedCatchAtAge[,length(ListOfPopulations_WithProjections[[i]]$ProjectedCatchAtAge)]*ListOfPopulations_WithProjections[[i]]$NumAtAge$Weight_kg)
    if(ListOfPopulations_WithProjections[[i]]$FishingMortaltiyProjectionScenario=="F_rebuild_Tmin_oneGeneration")
    {
      if(ListOfPopulations_WithProjections[[i]]$ProjectionType=="RegularF_willNotRebuildAtF=0")
      {
        SummaryProjectionScenarios$Tmin[i] = -999
        SummaryProjectionScenarios$OneGenerationTime[i] = -999
      } else {
        SummaryProjectionScenarios$Tmin[i] = ListOfPopulations_WithProjections[[i]]$timeToRebuild_noFishing
        SummaryProjectionScenarios$OneGenerationTime[i] = ListOfPopulations_WithProjections[[i]]$OneGenerationTime
      }
    }
    if(i%%1000==0)
    {
      print(paste("Summarizing iteration ",i," of ",length(ListOfPopulations_WithProjections),"...",sep=""))
      flush.console()
    }
  }
  SummaryProjectionScenarios1 = SummaryProjectionScenarios
  write.table(SummaryProjectionScenarios1,paste(PATH_output,"SummaryProjectionScenarios_RAW.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
  SummaryProjectionScenarios = read.table(paste(PATH_output,"SummaryProjectionScenarios_RAW.csv",sep="/"),sep=",",header=TRUE)
  BiomassCarryingCapacityMSY = read.table(paste(PATH_output,"BiomassCarryingCapacityMSY.csv",sep="/"),header=TRUE,sep=",")
  names(BiomassCarryingCapacityMSY)[names(BiomassCarryingCapacityMSY)=="Scenario"]="CatchMSY_Scenario"
  names(BiomassCarryingCapacityMSY)[names(BiomassCarryingCapacityMSY)=="Bound"]="CatchMSY_Scenario_Bound"
  names(BiomassCarryingCapacityMSY)[names(BiomassCarryingCapacityMSY)=="CatchMSY_Bound"]="CatchMSY_Scenario_Bound"
  SummaryProjectionScenarios = merge(SummaryProjectionScenarios,BiomassCarryingCapacityMSY,all.x=TRUE,by=c("Genus","Species","WPP","CatchMSY_Scenario","CatchMSY_Scenario_Bound"))
  CurrentStatusValues = subset(SummaryProjectionScenarios,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,
    populationReconstructionMethod,FishingMortaltiyProjectionScenario,ProjectionType,SteepnessAssumed,F_project,BiomassAtEnd,NumFishAtEnd,RecruitmentAtEnd,SSBatEnd,CatchAtEnd))
  CurrentStatusValues = CurrentStatusValues[CurrentStatusValues$ProjectionType=="RegularF" | CurrentStatusValues$ProjectionType=="ConstantLinearRecruitment",]
  CurrentStatusValues$ProjectionType=NULL
  BenchmarksSPR30 = CurrentStatusValues[CurrentStatusValues$FishingMortaltiyProjectionScenario=="F_at_SPR30",]
  BenchmarksSPR30_agg = aggregate.data.frame(list(BenchmarksSPR30$F_project,BenchmarksSPR30$NumFishAtEnd,BenchmarksSPR30$BiomassAtEnd,BenchmarksSPR30$RecruitmentAtEnd,BenchmarksSPR30$SSBatEnd,BenchmarksSPR30$CatchAtEnd),by=list(BenchmarksSPR30$Genus,BenchmarksSPR30$Species,BenchmarksSPR30$WPP,BenchmarksSPR30$CatchMSY_Scenario,BenchmarksSPR30$CatchMSY_Scenario_Bound,BenchmarksSPR30$LifeHistory_Scenario,BenchmarksSPR30$populationReconstructionMethod),FUN=mean)
  names(BenchmarksSPR30_agg) = c("Genus","Species","WPP","CatchMSY_Scenario","CatchMSY_Scenario_Bound","LifeHistory_Scenario","populationReconstructionMethod","Fspr30","Nspr30","Bspr30","Rspr30","SSBspr30","Yieldspr30")
  BenchmarksSPR20 = CurrentStatusValues[CurrentStatusValues$FishingMortaltiyProjectionScenario=="F_at_SPR20",]
  BenchmarksSPR20_agg = aggregate.data.frame(list(BenchmarksSPR20$F_project,BenchmarksSPR20$NumFishAtEnd,BenchmarksSPR20$BiomassAtEnd,BenchmarksSPR20$RecruitmentAtEnd,BenchmarksSPR20$SSBatEnd,BenchmarksSPR20$CatchAtEnd),by=list(BenchmarksSPR20$Genus,BenchmarksSPR20$Species,BenchmarksSPR20$WPP,BenchmarksSPR20$CatchMSY_Scenario,BenchmarksSPR20$CatchMSY_Scenario_Bound,BenchmarksSPR20$LifeHistory_Scenario,BenchmarksSPR20$populationReconstructionMethod),FUN=mean)
  names(BenchmarksSPR20_agg) = c("Genus","Species","WPP","CatchMSY_Scenario","CatchMSY_Scenario_Bound","LifeHistory_Scenario","populationReconstructionMethod","Fspr20","Nspr20","Bspr20","Rspr20","SSBspr20","Yieldspr20")
  BenchmarksSPR40 = CurrentStatusValues[CurrentStatusValues$FishingMortaltiyProjectionScenario=="F_at_SPR40",]
  BenchmarksSPR40_agg = aggregate.data.frame(list(BenchmarksSPR40$F_project,BenchmarksSPR40$NumFishAtEnd,BenchmarksSPR40$BiomassAtEnd,BenchmarksSPR40$RecruitmentAtEnd,BenchmarksSPR40$SSBatEnd,BenchmarksSPR40$CatchAtEnd),by=list(BenchmarksSPR40$Genus,BenchmarksSPR40$Species,BenchmarksSPR40$WPP,BenchmarksSPR40$CatchMSY_Scenario,BenchmarksSPR40$CatchMSY_Scenario_Bound,BenchmarksSPR40$LifeHistory_Scenario,BenchmarksSPR40$populationReconstructionMethod),FUN=mean)
  names(BenchmarksSPR40_agg) = c("Genus","Species","WPP","CatchMSY_Scenario","CatchMSY_Scenario_Bound","LifeHistory_Scenario","populationReconstructionMethod","Fspr40","Nspr40","Bspr40","Rspr40","SSBspr40","Yieldspr40")
  Benchmarks = merge(BenchmarksSPR20_agg,BenchmarksSPR30_agg,all.x=TRUE,by=c("Genus","Species","WPP","CatchMSY_Scenario","CatchMSY_Scenario_Bound","LifeHistory_Scenario","populationReconstructionMethod"))
  Benchmarks = merge(Benchmarks,BenchmarksSPR40_agg,all.x=TRUE,by=c("Genus","Species","WPP","CatchMSY_Scenario","CatchMSY_Scenario_Bound","LifeHistory_Scenario","populationReconstructionMethod"))
  write.table(Benchmarks,paste(PATH_output,"Benchmarks.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
  write.table(CurrentStatusValues,paste(PATH_output,"CurrentStatusValues_RAW.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
  write.table(BenchmarksSPR30,paste(PATH_output,"BenchmarksSPR30.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
  write.table(BenchmarksSPR20,paste(PATH_output,"BenchmarksSPR20.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
  write.table(BenchmarksSPR40,paste(PATH_output,"BenchmarksSPR40.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
  #Compute Ratio Values
  CurrentValues = subset(SummaryProjectionScenarios,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,
    populationReconstructionMethod,F_estimated,SSB_estimated,Yield_estimated,N_estimated,B_estimated))
  CurrentValues = unique(CurrentValues)
  CurrentValues = merge(CurrentValues,Benchmarks,by=c("Genus","Species","WPP","CatchMSY_Scenario","CatchMSY_Scenario_Bound","LifeHistory_Scenario","populationReconstructionMethod"),all.x=TRUE)
  CurrentValues$F_ratio30 = CurrentValues$F_estimated/CurrentValues$Fspr30
  CurrentValues$N_ratio30 = CurrentValues$N_estimated/CurrentValues$Nspr30
  CurrentValues$B_ratio30 = CurrentValues$B_estimated/CurrentValues$Bspr30
  CurrentValues$SSB_ratio30 = CurrentValues$SSB_estimated/CurrentValues$SSBspr30
  CurrentValues$Yield_ratio30 = CurrentValues$Yield_estimated/CurrentValues$Yieldspr30
  CurrentValues$F_ratio20 = CurrentValues$F_estimated/CurrentValues$Fspr20
  CurrentValues$N_ratio20 = CurrentValues$N_estimated/CurrentValues$Nspr20
  CurrentValues$B_ratio20 = CurrentValues$B_estimated/CurrentValues$Bspr20
  CurrentValues$SSB_ratio20 = CurrentValues$SSB_estimated/CurrentValues$SSBspr20
  CurrentValues$Yield_ratio20 = CurrentValues$Yield_estimated/CurrentValues$Yieldspr20
  CurrentValues$F_ratio40 = CurrentValues$F_estimated/CurrentValues$Fspr40
  CurrentValues$N_ratio40 = CurrentValues$N_estimated/CurrentValues$Nspr40
  CurrentValues$B_ratio40 = CurrentValues$B_estimated/CurrentValues$Bspr40
  CurrentValues$SSB_ratio40 = CurrentValues$SSB_estimated/CurrentValues$SSBspr40
  CurrentValues$Yield_ratio40 = CurrentValues$Yield_estimated/CurrentValues$Yieldspr40
  CurrentValues = unique(CurrentValues)
  write.table(CurrentValues,paste(PATH_output,"CurrentValuesBenchmarksAndRatios.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
  write.table(SummaryProjectionScenarios,paste(PATH_output,"SummaryProjectionScenarios.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
}






#######################################################################################################################################################
#######################################################################################################################################################
################# MAKE PLOTS AND DEVELOP OUTPUT PRODUCTS ##############################################################################################
#######################################################################################################################################################
#######################################################################################################################################################


##################### Create Kobe Plots and Percent Change F Table: Plot Ratio Current to Benchmarks ##########################
if(createKobePlots==TRUE)
{
  #Calculate Percent Change in F table
  CurrentStatusValues = read.table(paste(PATH_output,"CurrentValuesBenchmarksAndRatios.csv",sep="/"),header=TRUE,sep=",",as.is=TRUE)
  #TEMP CODE LINE: do NOT include SolveBaranov or TropFishR
  CurrentStatusValues = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="SolveBaranov",]   #TEMP CODE LINE: do NOT include SolveBaranov or TropFishR
  CurrentStatusValues = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="TropFishR",]      #TEMP CODE LINE: do NOT include SolveBaranov or TropFishR
  CurrentStatusValues = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="LIME",]
  CurrentStatusValues$GenusSpecies = paste(CurrentStatusValues$Genus,CurrentStatusValues$Species,sep="_")
  CurrentStatusValues$F_ratioAboveOne = 0
  CurrentStatusValues$F_ratioAboveOne[CurrentStatusValues$F_ratio30>=1]=1
  CurrentStatusValues$B_ratioAboveOne = 0
  CurrentStatusValues$B_ratioAboveOne[CurrentStatusValues$B_ratio30>=1]=1
  CurrentStatusValues$Quadrant=0      #Quadrant 1=Overfished and Overfishing; 2=Overfishing, not overfished; 3=not overfishing, not overfished; 4=not overfishing, overfished
  CurrentStatusValues$Quadrant[CurrentStatusValues$F_ratio30>=1 & CurrentStatusValues$B_ratio30<1]=1
  CurrentStatusValues$Quadrant[CurrentStatusValues$F_ratio30>=1 & CurrentStatusValues$B_ratio30>=1]=2
  CurrentStatusValues$Quadrant[CurrentStatusValues$F_ratio30<1 & CurrentStatusValues$B_ratio30>=1]=3
  CurrentStatusValues$Quadrant[CurrentStatusValues$F_ratio30<1 & CurrentStatusValues$B_ratio30<1]=4
  StatusTable = data.frame(table(CurrentStatusValues$GenusSpecies,CurrentStatusValues$WPP,CurrentStatusValues$Quadrant))
  names(StatusTable) = c("Species","WPP","Quadrant","Freq")
  StatusTable = reshape(StatusTable,v.names="Freq",idvar=c("Species","WPP"),timevar="Quadrant",direction="wide")
  StatusTable$Total = StatusTable$Freq.1 + StatusTable$Freq.2 + StatusTable$Freq.3 + StatusTable$Freq.4
  StatusTable$Percent.1 = round((StatusTable$Freq.1/StatusTable$Total)*100)
  StatusTable$Percent.2 = round((StatusTable$Freq.2/StatusTable$Total)*100)
  StatusTable$Percent.3 = round((StatusTable$Freq.3/StatusTable$Total)*100)
  StatusTable$Percent.4 = round((StatusTable$Freq.4/StatusTable$Total)*100)
  StatusTable[is.na(StatusTable)]=0
  StatusTable$Species = as.character(StatusTable$Species)
  StatusTable$WPP = as.character(StatusTable$WPP)
  write.table(StatusTable,paste(PATH_output,"StatusTable.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
  #Now make the Kobe plots
  CurrentStatusValues = read.table(paste(PATH_output,"CurrentValuesBenchmarksAndRatios.csv",sep="/"),header=TRUE,sep=",")
  #Remove runs that you don't want: here remove SolveBaranov and VPA and LIME
  CurrentStatusValues = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="SolveBaranov",]
  CurrentStatusValues = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="TropFishR",]
  CurrentStatusValues = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="LIME",]
  BenchmarksToPlot = subset(CurrentStatusValues,select=c(Genus,Species,WPP))
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
    thisIteration = CurrentStatusValues[CurrentStatusValues$Genus==BenchmarksToPlot$Genus[i] & CurrentStatusValues$Species==BenchmarksToPlot$Species[i] & CurrentStatusValues$WPP==BenchmarksToPlot$WPP[i],]
    FolderNameSpp = paste(as.character(BenchmarksToPlot$Genus[i]),CapStr(as.character(BenchmarksToPlot$Species[i])),sep="")
    png(paste(paste(PATH_benchmarkPlots,FolderNameSpp,sep="/"),paste(paste(BenchmarksToPlot$Genus[i],BenchmarksToPlot$Species[i],BenchmarksToPlot$WPP[i],sep="_"),"_Benchmarks.png",sep=""),sep="/"),units="px",width=5800,height=3200,res=600)
    layout(matrix(c(1,2), nrow = 1), widths = c(1,0.5))
    par(mar = c(5, 4, 4, 2) + 0.1)
    maxBratio = max(thisIteration$B_ratio30)
    maxFratio = max(thisIteration$F_ratio30)
    if(maxBratio<2) {maxBratio=2}
    if(maxBratio>10) {maxBratio=10}
    if(maxFratio<2) {maxFratio=2}
    if(maxFratio>10) {maxFratio=10}
    plot(NULL,xlab="Bcurrent/B@SPR30",ylab="Fcurrent/F@SPR30",xlim=c(0,ceiling(maxBratio)),ylim=c(0,ceiling(maxFratio)),xaxs="i",yaxs="i")
    rect(xleft=0,ybottom=0,xright=1,ytop=1,par("usr")[3] - 1, y0 + 3, par("usr")[4] + 1,border="khaki",col="khaki")
    rect(xleft=1,ybottom=0,xright=ceiling(maxBratio),ytop=1,par("usr")[3] - 1, y0 + 3, par("usr")[4] + 1,border="lightgreen",col="lightgreen")
    rect(xleft=0,ybottom=1,xright=1,ytop=ceiling(maxFratio),par("usr")[3] - 1, y0 + 3, par("usr")[4] + 1,border="plum1",col="plum1")
    rect(xleft=1,ybottom=1,xright=ceiling(maxBratio),ytop=ceiling(maxFratio),par("usr")[3] - 1, y0 + 3, par("usr")[4] + 1,border="bisque",col="bisque") #tan1
    for(j in 1:dim(thisIteration)[1])
    {
      if(j==1)
      {
        points(thisIteration$B_ratio30[j],thisIteration$F_ratio30[j],xlab="Bcurrent/B@SPR30",ylab="Fcurrent/F@SPR30",pch=thisIteration$pch_arg[j],col=thisIteration$Colors[j])
      }
      if(j>1)
      {
        points(thisIteration$B_ratio30[j],thisIteration$F_ratio30[j],xlab="Bcurrent/B@SPR30",ylab="Fcurrent/F@SPR30",pch=thisIteration$pch_arg[j],col=thisIteration$Colors[j])
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
    StatusPercents = StatusTable[StatusTable$Species==paste(as.character(BenchmarksToPlot$Genus[i]),as.character(BenchmarksToPlot$Species[i]),sep="_") & StatusTable$WPP==BenchmarksToPlot$WPP[i],]
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
    title(paste("Current Status: ",BenchmarksToPlot$Genus[i],BenchmarksToPlot$Species[i],sep=" "))
    mtext(paste("WPP = ",BenchmarksToPlot$WPP[i],sep=""))
    par(mar = c(5, 0, 4, 2))
    plot(thisIteration$B_ratio[dim(thisIteration)[1]],thisIteration$F_ratio[dim(thisIteration)[1]],type="n",axes=FALSE,ann=FALSE,xlim=c(0,maxBratio),ylim=c(0,maxFratio))
    LengthReconstructionSymbols_forLegend = LengthReconstructionSymbols
    LengthReconstructionSymbols_forLegend$populationReconstructionMethod = as.character(LengthReconstructionSymbols_forLegend$populationReconstructionMethod)
    LengthReconstructionSymbols_forLegend$populationReconstructionMethod[LengthReconstructionSymbols_forLegend$populationReconstructionMethod=="BevertonHoltInstantaneousMortality"]="BevHoltF"
    legend("topleft",legend=c(as.character(LengthReconstructionSymbols_forLegend$populationReconstructionMethod),as.character(ColorsFor_benchmarkPlot_CmsyScenarios$CatchMSY_Scenario)),lty=c(rep(NA,times=dim(LengthReconstructionSymbols_forLegend)[1]),rep(1,times=dim(ColorsFor_benchmarkPlot_CmsyScenarios)[1])),lwd=c(rep(NA,times=dim(LengthReconstructionSymbols_forLegend)[1]),rep(2,times=dim(ColorsFor_benchmarkPlot_CmsyScenarios)[1])),col=c(rep("black",times=dim(LengthReconstructionSymbols_forLegend)[1]),ColorsFor_benchmarkPlot_CmsyScenarios$Colors),pch=c(LengthReconstructionSymbols_forLegend$pch_arg,rep(NA,times=dim(ColorsFor_benchmarkPlot_CmsyScenarios)[1])),cex=0.8)
    dev.off()
  }
}

######################## Create Percentage Change in F Tables and Figures ###################################

Benchmarks = read.table(paste(PATH_output,"CurrentValuesBenchmarksAndRatios.csv",sep="/"),header=TRUE,sep=",")
#TEMP CODE LINE: do NOT include SolveBaranov or TropFishR or LIME
Benchmarks = Benchmarks[Benchmarks$populationReconstructionMethod!="TropFishR",]
Benchmarks = Benchmarks[Benchmarks$populationReconstructionMethod!="SolveBaranov",]
Benchmarks = Benchmarks[Benchmarks$populationReconstructionMethod!="LIME",]
Benchmarks$F_percentChangeToSPR30 = round(((Benchmarks$Fspr30 - Benchmarks$F_estimated)/Benchmarks$F_estimated)*100)
Benchmarks = subset(Benchmarks,select=c("CatchMSY_Scenario","Genus","Species","WPP","LifeHistory_Scenario","populationReconstructionMethod","F_percentChangeToSPR30"))
Benchmarks = unique(Benchmarks)
SpeciesWPP = subset(Benchmarks,select=c("Genus","Species","WPP"))
SpeciesWPP = unique(SpeciesWPP)
for(i in 1:dim(SpeciesWPP)[1])
{
  Sub = Benchmarks[Benchmarks$Genus==SpeciesWPP$Genus[i] & Benchmarks$Species==SpeciesWPP$Species[i] &  Benchmarks$WPP==SpeciesWPP$WPP[i],]
  FolderNameSpp = paste(as.character(SpeciesWPP$Genus[i]),CapStr(as.character(SpeciesWPP$Species[i])),sep="")
  png(paste(paste(PATH_precentChangeF_histograms,FolderNameSpp,sep="/"),paste(paste(SpeciesWPP$Genus[i],SpeciesWPP$Species[i],SpeciesWPP$WPP[i],sep="_"),"_F_percentChangeHistogram.png",sep=""),sep="/"),units="px",width=2200,height=2200,res=600)
  hist(Sub$F_percentChangeToSPR30,main=paste(SpeciesWPP$Genus[i],SpeciesWPP$Species[i],sep=" "),xlab="Percent Change")
  mtext(paste("WPP = ",SpeciesWPP$WPP[i],sep=""))
  dev.off()
}
                                         

########################### Plot Outputs by Projection Type #################################################
if(plotProjectionsByProjectionType==TRUE)
{
  ListOfPopulations_WithProjections = readRDS(file=paste(PATH_output,"ListOfPopulations_WithProjections.RData",sep="/"),refhook=NULL)
  SummaryProjectionScenarios = read.table(paste(PATH_output,"SummaryProjectionScenarios.csv",sep="/"),header=TRUE,sep=",",as.is=TRUE)
  #TEMP CODE LINE: do NOT include SolveBaranov or TropFishR or LIME
  SummaryProjectionScenarios = SummaryProjectionScenarios[SummaryProjectionScenarios$populationReconstructionMethod!="TropFishR",]
  SummaryProjectionScenarios = SummaryProjectionScenarios[SummaryProjectionScenarios$populationReconstructionMethod!="SolveBaranov",]
  SummaryProjectionScenarios = SummaryProjectionScenarios[SummaryProjectionScenarios$populationReconstructionMethod!="LIME",]
  #ONLY PLOT REGULAR F RUNS!!!!
  SummaryProjectionScenarios = SummaryProjectionScenarios[SummaryProjectionScenarios$ProjectionType=="RegularF" | SummaryProjectionScenarios$ProjectionType=="RegularF_willNotRebuildAtF=0",]
  PlotsToMake = subset(SummaryProjectionScenarios,select=c(Genus,Species,WPP,FishingMortaltiyProjectionScenario))
  PlotsToMake= unique(PlotsToMake)
#  PlotsToMake = subset(PlotsToMake,PlotsToMake$WPP=="712")    #TEMPORARY LINE OF CODE TO GET JUST WHAT WE WANT TO PLOT - NOT EVERYTHING YET
  WillNotRebuild = PlotsToMake
  WillNotRebuild$TotalRuns = dim(PlotsToMake)[1]
  WillNotRebuild$RunsThatDontRebuild = 0
  for(i in 1:dim(PlotsToMake)[1])
  {
    ScenariosToPlot = SummaryProjectionScenarios[SummaryProjectionScenarios$Genus==PlotsToMake$Genus[i] &
      SummaryProjectionScenarios$Species==PlotsToMake$Species[i] & SummaryProjectionScenarios$WPP==PlotsToMake$WPP[i] &
      SummaryProjectionScenarios$FishingMortaltiyProjectionScenario==PlotsToMake$FishingMortaltiyProjectionScenario[i],]
    PopReconstructionLineType = data.frame(populationReconstructionMethod=sort(unique(as.character(ScenariosToPlot$populationReconstructionMethod))),LineType=seq(1,length(unique(ScenariosToPlot$populationReconstructionMethod)),by=1))
    ScenariosToPlot = merge(ScenariosToPlot,PopReconstructionLineType,all.x=TRUE,by=c("populationReconstructionMethod"))
    ScenariosToPlot$Scenarios_Colors = paste(ScenariosToPlot$Scenario_New,ScenariosToPlot$SteepnessAssumed,sep="-")
    Scenario_Colors = data.frame(CatchMSY_Scenario=unique(ScenariosToPlot$CatchMSY_Scenario),Colors=ColorsFor_ProjectionF_ScenarioPlots[1:length(unique(ScenariosToPlot$CatchMSY_Scenario))])
    ScenariosToPlot =  merge(ScenariosToPlot,Scenario_Colors,all.x=TRUE,by=c("CatchMSY_Scenario"))
    ScenariosToPlot_noRebuild = ScenariosToPlot[ScenariosToPlot$ProjectionType=="RegularF_willNotRebuildAtF=0",]
    WillNotRebuild$RunsThatDontRebuild[i]=dim(ScenariosToPlot_noRebuild)[1]
    ScenariosToPlot = ScenariosToPlot[ScenariosToPlot$ProjectionType!="RegularF_willNotRebuildAtF=0",]
    if(dim(ScenariosToPlot)[1]==0)
    {
      next
    }
    for(j in 1:dim(ScenariosToPlot)[1])
    {
      if(j==1)
      {
        Yield = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$ProjectedCatchYear)
        names(Yield) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        SSB = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$SSB)
        names(SSB) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Abundance = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$ProjectedPopulationYear)
        names(Abundance) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Recruitment = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$Recruitment)
        names(Recruitment) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Biomass = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$ProjectedBiomassYear)
        names(Biomass) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
      }
      if(j>1)
      {
        Yield_j = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$ProjectedCatchYear)
        names(Yield_j) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Yield = cbind(Yield,Yield_j)
        SSB_j = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$SSB)
        names(SSB_j) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        SSB = cbind(SSB,SSB_j)
        Abundance_j = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$ProjectedPopulationYear)
        names(Abundance_j) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Abundance = cbind(Abundance,Abundance_j)
        Recruitment_j = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$Recruitment)
        names(Recruitment_j) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Recruitment = cbind(Recruitment,Recruitment_j)
        Biomass_j = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$ProjectedBiomassYear)
        names(Biomass_j) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Biomass = cbind(Biomass,Biomass_j)
      }
    }
    Yield$index = seq(1,dim(Yield)[1],by=1)
    Yield = Yield[Yield$index<projectionYearsToPlot,]
    Yield$index=NULL
    SSB$index = seq(1,dim(SSB)[1],by=1)
    SSB = SSB[SSB$index<projectionYearsToPlot,]
    SSB$index=NULL
    Abundance$index = seq(1,dim(Abundance)[1],by=1)
    Abundance = Abundance[Abundance$index<projectionYearsToPlot,]
    Abundance$index=NULL
    Recruitment$index = seq(1,dim(Recruitment)[1],by=1)
    Recruitment = Recruitment[Recruitment$index<projectionYearsToPlot,]
    Recruitment$index=NULL
    Biomass$index = seq(1,dim(Biomass)[1],by=1)
    Biomass = Biomass[Biomass$index<projectionYearsToPlot,]
    Biomass$index=NULL
    FolderNameSpp = paste(as.character(ScenariosToPlot$Genus[1]),CapStr(as.character(ScenariosToPlot$Species[1])),sep="")
    png(paste(paste(PATH_projectionPlots_byProjectType,FolderNameSpp,sep="/"),paste(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],
       ScenariosToPlot$WPP[1],ScenariosToPlot$CatchMSY_Scenario[1],ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep="_"),"_Projected_Yield.png",sep=""),sep="/"),units="px",width=5200,height=3200,res=600)
    layout(matrix(c(1,2), nrow = 1), widths = c(1,0.5))
    par(mar = c(5, 4, 4, 2) + 0.1)
    for(j in 1:dim(Yield)[2])
    {
      if(j==1)
      {
        plot(as.numeric(rownames(Yield)),Yield[,j],xlab="Year",ylab="Yield (numbers)",ylim=c(0,max(Yield)),type="l",lty=ScenariosToPlot$LineType[j],lwd=2,col=as.character(ScenariosToPlot$Colors[j]))
        title(paste(ScenariosToPlot$Genus[1]," ",ScenariosToPlot$Species[1]," (WPP ", ScenariosToPlot$WPP[1],"): ",ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep=""))
        mtext(paste(paste("CatchMSY Scenario = ",ScenariosToPlot$CatchMSY_Scenario[1],sep=""),"BevHoltRecruit",sep=" | "),cex=0.7)
      }
      if(j>1)
      {
        points(as.numeric(rownames(Yield)),Yield[,j],type="l",lty=ScenariosToPlot$LineType[j],lwd=2,col=as.character(ScenariosToPlot$Colors[j]))
      }
    }
    par(mar = c(5, 0, 4, 2))
    plot(as.numeric(rownames(Yield)),Yield[,j],type="n",axes=FALSE,ann=FALSE)
    MortalityNames=as.character(PopReconstructionLineType$populationReconstructionMethod)
    MortalityNames[MortalityNames=="BevertonHoltInstantaneousMortality"]="BevHoltF"
    legend("topleft",legend=c(MortalityNames,as.character(Scenario_Colors$CatchMSY_Scenario)),lty=c(PopReconstructionLineType$LineType,rep(1,times=dim(Scenario_Colors)[1])),lwd=2,col=c(rep("black",times=dim(PopReconstructionLineType)[1]),as.character(Scenario_Colors$Colors)),cex=1)
    dev.off()
    #Abundance Plots
    FolderNameSpp = paste(as.character(ScenariosToPlot$Genus[1]),CapStr(as.character(ScenariosToPlot$Species[1])),sep="")
    png(paste(paste(PATH_projectionPlots_byProjectType,FolderNameSpp,sep="/"),paste(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],
       ScenariosToPlot$WPP[1],ScenariosToPlot$CatchMSY_Scenario[1],ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep="_"),"_Projected_Abundance.png",sep=""),sep="/"),units="px",width=5200,height=3200,res=600)
    layout(matrix(c(1,2), nrow = 1), widths = c(1,0.5))
    par(mar = c(5, 4, 4, 2) + 0.1)
    for(j in 1:dim(Abundance)[2])
    {
      if(j==1)
      {
        plot(as.numeric(rownames(Abundance)),Abundance[,j],xlab="Year",ylab="Abundance (numbers)",ylim=c(0,max(Abundance)),type="l",lty=ScenariosToPlot$LineType[j],lwd=2,col=as.character(ScenariosToPlot$Colors[j]))
        title(paste(ScenariosToPlot$Genus[1]," ",ScenariosToPlot$Species[1]," (WPP ", ScenariosToPlot$WPP[1],"): ",ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep=""))
        mtext(paste(paste("CatchMSY Scenario = ",ScenariosToPlot$CatchMSY_Scenario[1],sep=""),"BevHoltRecruit",sep=" | "),cex=0.7)
      }
      if(j>1)
      {
        points(as.numeric(rownames(Abundance)),Abundance[,j],type="l",lty=ScenariosToPlot$LineType[j],lwd=2,col=as.character(ScenariosToPlot$Colors[j]))
      }
    }
    par(mar = c(5, 0, 4, 2))
    plot(as.numeric(rownames(Abundance)),Abundance[,j],type="n",axes=FALSE,ann=FALSE)
    MortalityNames=as.character(PopReconstructionLineType$populationReconstructionMethod)
    MortalityNames[MortalityNames=="BevertonHoltInstantaneousMortality"]="BevHoltF"
    legend("topleft",legend=c(MortalityNames,as.character(Scenario_Colors$CatchMSY_Scenario)),lty=c(PopReconstructionLineType$LineType,rep(1,times=dim(Scenario_Colors)[1])),lwd=2,col=c(rep("black",times=dim(PopReconstructionLineType)[1]),as.character(Scenario_Colors$Colors)),cex=1)
    dev.off()
    #SSB Plots
    FolderNameSpp = paste(as.character(ScenariosToPlot$Genus[1]),CapStr(as.character(ScenariosToPlot$Species[1])),sep="")
    png(paste(paste(PATH_projectionPlots_byProjectType,FolderNameSpp,sep="/"),paste(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],
       ScenariosToPlot$WPP[1],ScenariosToPlot$CatchMSY_Scenario[1],ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep="_"),"_Projected_SSB.png",sep=""),sep="/"),units="px",width=5200,height=3200,res=600)
    layout(matrix(c(1,2), nrow = 1), widths = c(1,0.5))
    par(mar = c(5, 4, 4, 2) + 0.1)
    for(j in 1:dim(SSB)[2])
    {
      if(j==1)
      {
        plot(as.numeric(rownames(SSB)),SSB[,j],xlab="Year",ylab="SSB",ylim=c(0,max(SSB)),type="l",lty=ScenariosToPlot$LineType[j],lwd=2,col=as.character(ScenariosToPlot$Colors[j]))
        title(paste(ScenariosToPlot$Genus[1]," ",ScenariosToPlot$Species[1]," (WPP ", ScenariosToPlot$WPP[1],"): ",ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep=""))
        mtext(paste(paste("CatchMSY Scenario = ",ScenariosToPlot$CatchMSY_Scenario[1],sep=""),"BevHoltRecruit",sep=" | "),cex=0.7)
      }
      if(j>1)
      {
        points(as.numeric(rownames(SSB)),SSB[,j],type="l",lty=ScenariosToPlot$LineType[j],lwd=2,col=as.character(ScenariosToPlot$Colors[j]))
      }
    }
    par(mar = c(5, 0, 4, 2))
    plot(as.numeric(rownames(SSB)),SSB[,j],type="n",axes=FALSE,ann=FALSE)
    MortalityNames=as.character(PopReconstructionLineType$populationReconstructionMethod)
    MortalityNames[MortalityNames=="BevertonHoltInstantaneousMortality"]="BevHoltF"
    legend("topleft",legend=c(MortalityNames,as.character(Scenario_Colors$CatchMSY_Scenario)),lty=c(PopReconstructionLineType$LineType,rep(1,times=dim(Scenario_Colors)[1])),lwd=2,col=c(rep("black",times=dim(PopReconstructionLineType)[1]),as.character(Scenario_Colors$Colors)),cex=1)
    dev.off()
    #Recruitment Plots
    FolderNameSpp = paste(as.character(ScenariosToPlot$Genus[1]),CapStr(as.character(ScenariosToPlot$Species[1])),sep="")
    png(paste(paste(PATH_projectionPlots_byProjectType,FolderNameSpp,sep="/"),paste(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],
       ScenariosToPlot$WPP[1],ScenariosToPlot$CatchMSY_Scenario[1],ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep="_"),"_Projected_Recruitment.png",sep=""),sep="/"),units="px",width=5200,height=3200,res=600)
    layout(matrix(c(1,2), nrow = 1), widths = c(1,0.5))
    par(mar = c(5, 4, 4, 2) + 0.1)
    for(j in 1:dim(Recruitment)[2])
    {
      if(j==1)
      {
        plot(as.numeric(rownames(Recruitment)),Recruitment[,j],xlab="Year",ylab="Recruitment",ylim=c(0,max(Recruitment)),type="l",lty=ScenariosToPlot$LineType[j],lwd=2,col=as.character(ScenariosToPlot$Colors[j]))
        title(paste(ScenariosToPlot$Genus[1]," ",ScenariosToPlot$Species[1]," (WPP ", ScenariosToPlot$WPP[1],"): ",ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep=""))
        mtext(paste(paste("CatchMSY Scenario = ",ScenariosToPlot$CatchMSY_Scenario[1],sep=""),"BevHoltRecruit",sep=" | "),cex=0.7)
      }
      if(j>1)
      {
        points(as.numeric(rownames(Recruitment)),Recruitment[,j],type="l",lty=ScenariosToPlot$LineType[j],lwd=2,col=as.character(ScenariosToPlot$Colors[j]))
      }
    }
    par(mar = c(5, 0, 4, 2))
    plot(as.numeric(rownames(Recruitment)),Recruitment[,j],type="n",axes=FALSE,ann=FALSE)
    MortalityNames=as.character(PopReconstructionLineType$populationReconstructionMethod)
    MortalityNames[MortalityNames=="BevertonHoltInstantaneousMortality"]="BevHoltF"
    legend("topleft",legend=c(MortalityNames,as.character(Scenario_Colors$CatchMSY_Scenario)),lty=c(PopReconstructionLineType$LineType,rep(1,times=dim(Scenario_Colors)[1])),lwd=2,col=c(rep("black",times=dim(PopReconstructionLineType)[1]),as.character(Scenario_Colors$Colors)),cex=1)
    dev.off()
    #Biomass Plots
    FolderNameSpp = paste(as.character(ScenariosToPlot$Genus[1]),CapStr(as.character(ScenariosToPlot$Species[1])),sep="")
    png(paste(paste(PATH_projectionPlots_byProjectType,FolderNameSpp,sep="/"),paste(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],
       ScenariosToPlot$WPP[1],ScenariosToPlot$CatchMSY_Scenario[1],ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep="_"),"_Projected_Biomass.png",sep=""),sep="/"),units="px",width=5200,height=3200,res=600)
    layout(matrix(c(1,2), nrow = 1), widths = c(1,0.5))
    par(mar = c(5, 4, 4, 2) + 0.1)
    for(j in 1:dim(Biomass)[2])
    {
      if(j==1)
      {
        plot(as.numeric(rownames(Biomass)),Biomass[,j],xlab="Year",ylab="Biomass",ylim=c(0,max(Biomass)),type="l",lty=ScenariosToPlot$LineType[j],lwd=2,col=as.character(ScenariosToPlot$Colors[j]))
        title(paste(ScenariosToPlot$Genus[1]," ",ScenariosToPlot$Species[1]," (WPP ", ScenariosToPlot$WPP[1],"): ",ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep=""))
        mtext(paste(paste("CatchMSY Scenario = ",ScenariosToPlot$CatchMSY_Scenario[1],sep=""),"BevHoltRecruit",sep=" | "),cex=0.7)
      }
      if(j>1)
      {
        points(as.numeric(rownames(Biomass)),Biomass[,j],type="l",lty=ScenariosToPlot$LineType[j],lwd=2,col=as.character(ScenariosToPlot$Colors[j]))
      }
    }
    par(mar = c(5, 0, 4, 2))
    plot(as.numeric(rownames(Biomass)),Biomass[,j],type="n",axes=FALSE,ann=FALSE)
    MortalityNames=as.character(PopReconstructionLineType$populationReconstructionMethod)
    MortalityNames[MortalityNames=="BevertonHoltInstantaneousMortality"]="BevHoltF"
    legend("topleft",legend=c(MortalityNames,as.character(Scenario_Colors$CatchMSY_Scenario)),lty=c(PopReconstructionLineType$LineType,rep(1,times=dim(Scenario_Colors)[1])),lwd=2,col=c(rep("black",times=dim(PopReconstructionLineType)[1]),as.character(Scenario_Colors$Colors)),cex=1)
    dev.off()
    if(i%%20==0)
    {
      print(paste("Ploting Scenario ",i," of ",dim(PlotsToMake)[1],".......",sep=""))
      flush.console()
    }
  }
  write.table(WillNotRebuild,paste(PATH_output,"RunsThatDoNotRebuild.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
}


########################### Create Fan Plots #################################################
if(createFanPlots==TRUE)
{
  ListOfPopulations_WithProjections = readRDS(file=paste(PATH_output,"ListOfPopulations_WithProjections.RData",sep="/"),refhook=NULL)
  SummaryProjectionScenarios = read.table(paste(PATH_output,"SummaryProjectionScenarios.csv",sep="/"),header=TRUE,sep=",",as.is=TRUE)
  #TEMP CODE LINE: do NOT include SolveBaranov or TropFishR
  SummaryProjectionScenarios = SummaryProjectionScenarios[SummaryProjectionScenarios$populationReconstructionMethod!="TropFishR",]
  SummaryProjectionScenarios = SummaryProjectionScenarios[SummaryProjectionScenarios$populationReconstructionMethod!="SolveBaranov",]
  SummaryProjectionScenarios = SummaryProjectionScenarios[SummaryProjectionScenarios$populationReconstructionMethod!="LIME",]
  #ONLY PLOT REGULAR F RUNS!!!!
  SummaryProjectionScenarios = SummaryProjectionScenarios[SummaryProjectionScenarios$ProjectionType=="RegularF" | SummaryProjectionScenarios$ProjectionType=="RegularF_willNotRebuildAtF=0",]
  PlotsToMake = subset(SummaryProjectionScenarios,select=c(Genus,Species,WPP,FishingMortaltiyProjectionScenario))
  PlotsToMake= unique(PlotsToMake)
#  PlotsToMake = subset(PlotsToMake,PlotsToMake$WPP=="712")    #TEMPORARY LINE OF CODE TO GET JUST WHAT WE WANT TO PLOT - NOT EVERYTHING YET
  for(i in 1:dim(PlotsToMake)[1])
  {
    ScenariosToPlot = SummaryProjectionScenarios[SummaryProjectionScenarios$Genus==PlotsToMake$Genus[i] &
      SummaryProjectionScenarios$Species==PlotsToMake$Species[i] & SummaryProjectionScenarios$WPP==PlotsToMake$WPP[i] &
      SummaryProjectionScenarios$FishingMortaltiyProjectionScenario==PlotsToMake$FishingMortaltiyProjectionScenario[i],]
    PopReconstructionLineType = data.frame(populationReconstructionMethod=sort(unique(as.character(ScenariosToPlot$populationReconstructionMethod))),LineType=seq(1,length(unique(ScenariosToPlot$populationReconstructionMethod)),by=1))
    ScenariosToPlot = merge(ScenariosToPlot,PopReconstructionLineType,all.x=TRUE,by=c("populationReconstructionMethod"))
    ScenariosToPlot$Scenarios_Colors = paste(ScenariosToPlot$Scenario_New,ScenariosToPlot$SteepnessAssumed,sep="-")
    Scenario_Colors = data.frame(CatchMSY_Scenario=unique(ScenariosToPlot$CatchMSY_Scenario),Colors=ColorsFor_ProjectionF_ScenarioPlots[1:length(unique(ScenariosToPlot$CatchMSY_Scenario))])
    ScenariosToPlot =  merge(ScenariosToPlot,Scenario_Colors,all.x=TRUE,by=c("CatchMSY_Scenario"))
    #Remove those scenarios that won't rebuild - they are recorded and written out to a file when the spaghetti plots are created in that code block
    ScenariosToPlot = ScenariosToPlot[ScenariosToPlot$ProjectionType!="RegularF_willNotRebuildAtF=0",]
    if(dim(ScenariosToPlot)[1]==0)
    {
      next
    }
    for(j in 1:dim(ScenariosToPlot)[1])
    {
      if(j==1)
      {
        Yield = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$ProjectedCatchYear)
        names(Yield) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        SSB = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$SSB)
        names(SSB) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Abundance = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$ProjectedPopulationYear)
        names(Abundance) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Recruitment = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$Recruitment)
        names(Recruitment) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Biomass = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$ProjectedBiomassYear)
        names(Biomass) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
      }
      if(j>1)
      {
        Yield_j = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$ProjectedCatchYear)
        names(Yield_j) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Yield = cbind(Yield,Yield_j)
        SSB_j = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$SSB)
        names(SSB_j) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        SSB = cbind(SSB,SSB_j)
        Abundance_j = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$ProjectedPopulationYear)
        names(Abundance_j) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Abundance = cbind(Abundance,Abundance_j)
        Recruitment_j = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$Recruitment)
        names(Recruitment_j) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Recruitment = cbind(Recruitment,Recruitment_j)
        Biomass_j = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$ProjectedBiomassYear)
        names(Biomass_j) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Biomass = cbind(Biomass,Biomass_j)
      }
    }
    Yield$index = seq(1,dim(Yield)[1],by=1)
    Yield = Yield[Yield$index<projectionYearsToPlot,]
    Yield$index=NULL
    SSB$index = seq(1,dim(SSB)[1],by=1)
    SSB = SSB[SSB$index<projectionYearsToPlot,]
    SSB$index=NULL
    Abundance$index = seq(1,dim(Abundance)[1],by=1)
    Abundance = Abundance[Abundance$index<projectionYearsToPlot,]
    Abundance$index=NULL
    Recruitment$index = seq(1,dim(Recruitment)[1],by=1)
    Recruitment = Recruitment[Recruitment$index<projectionYearsToPlot,]
    Recruitment$index=NULL
    Biomass$index = seq(1,dim(Biomass)[1],by=1)
    Biomass = Biomass[Biomass$index<projectionYearsToPlot,]
    Biomass$index=NULL
    #Yield Plots
    FolderNameSpp = paste(as.character(ScenariosToPlot$Genus[1]),CapStr(as.character(ScenariosToPlot$Species[1])),sep="")
    png(paste(paste(PATH_summaryFanProjectionPlots,FolderNameSpp,sep="/"),paste(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],
       ScenariosToPlot$WPP[1],ScenariosToPlot$CatchMSY_Scenario[1],ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep="_"),"_Fan_Yield.png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
    plot(NULL,xlim = c(min(as.numeric(row.names(Yield))),max(as.numeric(row.names(Yield)))),ylim = c(min(Yield),max(Yield)),ylab="Yield (N)",xlab="Year")
    Yield_thisReconstruction_m = t(Yield)
    fan(data = Yield_thisReconstruction_m,start=2019,type="percentile",fan.col=topo.colors,ln=NULL,rlab=NULL)
    title(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],sep=" "))
    mtext(paste("WPP=",ScenariosToPlot$WPP[1]," | Projection Type=",ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep=""))
    dev.off()
    #SSB Plots
    FolderNameSpp = paste(as.character(ScenariosToPlot$Genus[1]),CapStr(as.character(ScenariosToPlot$Species[1])),sep="")
    png(paste(paste(PATH_summaryFanProjectionPlots,FolderNameSpp,sep="/"),paste(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],
       ScenariosToPlot$WPP[1],ScenariosToPlot$CatchMSY_Scenario[1],ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep="_"),"_Fan_SSB.png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
    plot(NULL,xlim = c(min(as.numeric(row.names(SSB))),max(as.numeric(row.names(SSB)))),ylim = c(min(SSB),max(SSB)),ylab="SSB",xlab="Year")
    SSB_thisReconstruction_m = t(SSB)
    fan(data = SSB_thisReconstruction_m,start=2019,type="percentile",fan.col=topo.colors,ln=NULL,rlab=NULL)
    title(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],sep=" "))
    mtext(paste("WPP=",ScenariosToPlot$WPP[1]," | Projection Type=",ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep=""))
    dev.off()
    #Abundance Plots
    FolderNameSpp = paste(as.character(ScenariosToPlot$Genus[1]),CapStr(as.character(ScenariosToPlot$Species[1])),sep="")
    png(paste(paste(PATH_summaryFanProjectionPlots,FolderNameSpp,sep="/"),paste(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],
       ScenariosToPlot$WPP[1],ScenariosToPlot$CatchMSY_Scenario[1],ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep="_"),"_Fan_Abundance.png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
    plot(NULL,xlim = c(min(as.numeric(row.names(Abundance))),max(as.numeric(row.names(Abundance)))),ylim = c(min(Abundance),max(Abundance)),ylab="Abundance",xlab="Year")
    Abundance_thisReconstruction_m = t(Abundance)
    fan(data = Abundance_thisReconstruction_m,start=2019,type="percentile",fan.col=topo.colors,ln=NULL,rlab=NULL)
    title(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],sep=" "))
    mtext(paste("WPP=",ScenariosToPlot$WPP[1]," | Projection Type=",ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep=""))
    dev.off()
    #Recruitment Plots
    FolderNameSpp = paste(as.character(ScenariosToPlot$Genus[1]),CapStr(as.character(ScenariosToPlot$Species[1])),sep="")
    png(paste(paste(PATH_summaryFanProjectionPlots,FolderNameSpp,sep="/"),paste(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],
       ScenariosToPlot$WPP[1],ScenariosToPlot$CatchMSY_Scenario[1],ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep="_"),"_Fan_Recruitment.png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
    plot(NULL,xlim = c(min(as.numeric(row.names(Recruitment))),max(as.numeric(row.names(Recruitment)))),ylim = c(min(Recruitment),max(Recruitment)),ylab="Recruitment",xlab="Year")
    Recruitment_thisReconstruction_m = t(Recruitment)
    fan(data = Recruitment_thisReconstruction_m,start=2019,type="percentile",fan.col=topo.colors,ln=NULL,rlab=NULL)
    title(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],sep=" "))
    mtext(paste("WPP=",ScenariosToPlot$WPP[1]," | Projection Type=",ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep=""))
    dev.off()
    #Biomass Plots
    FolderNameSpp = paste(as.character(ScenariosToPlot$Genus[1]),CapStr(as.character(ScenariosToPlot$Species[1])),sep="")
    png(paste(paste(PATH_summaryFanProjectionPlots,FolderNameSpp,sep="/"),paste(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],
       ScenariosToPlot$WPP[1],ScenariosToPlot$CatchMSY_Scenario[1],ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep="_"),"_Fan_Biomass.png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
    plot(NULL,xlim = c(min(as.numeric(row.names(Biomass))),max(as.numeric(row.names(Biomass)))),ylim = c(min(Biomass),max(Biomass)),ylab="Biomass",xlab="Year")
    Biomass_thisReconstruction_m = t(Biomass)
    fan(data = Biomass_thisReconstruction_m,start=2019,type="percentile",fan.col=topo.colors,ln=NULL,rlab=NULL)
    title(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],sep=" "))
    mtext(paste("WPP=",ScenariosToPlot$WPP[1]," | Projection Type=",ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep=""))
    dev.off()
    if(i%%20==0)
    {
      print(paste("Plotting Fan Plot Scenario ",i," of ",dim(PlotsToMake)[1],".......",sep=""))
      flush.console()
    }
  }
}


########################### Create Standardized Fan Plots - FINAL #################################################
if(createStandardizedFanPlots==TRUE)
{
  ListOfPopulations_WithProjections = readRDS(file=paste(PATH_output,"ListOfPopulations_WithProjections.RData",sep="/"),refhook=NULL)
  SummaryProjectionScenarios = read.table(paste(PATH_output,"SummaryProjectionScenarios.csv",sep="/"),header=TRUE,sep=",",as.is=TRUE)
  #TEMP CODE LINE: do NOT include SolveBaranov or TropFishR
  SummaryProjectionScenarios = SummaryProjectionScenarios[SummaryProjectionScenarios$populationReconstructionMethod!="TropFishR",]
  SummaryProjectionScenarios = SummaryProjectionScenarios[SummaryProjectionScenarios$populationReconstructionMethod!="SolveBaranov",]
  SummaryProjectionScenarios = SummaryProjectionScenarios[SummaryProjectionScenarios$populationReconstructionMethod!="LIME",]
  #ONLY PLOT REGULAR F RUNS!!!!
  SummaryProjectionScenarios = SummaryProjectionScenarios[SummaryProjectionScenarios$ProjectionType=="RegularF" | SummaryProjectionScenarios$ProjectionType=="RegularF_willNotRebuildAtF=0",]
  #Obtain Benchmark values from summary file: NOTE: do NOT use the Benchmarks30 data.frame, instead grab the value again, because here I need the benchmark values by steepness value as well in order to properly standardize for the standardized fan plots.
  CurrentStatusValues = subset(SummaryProjectionScenarios,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,
    populationReconstructionMethod,FishingMortaltiyProjectionScenario,ProjectionType,SteepnessAssumed,F_project,BiomassAtEnd,NumFishAtEnd,RecruitmentAtEnd,SSBatEnd,CatchAtEnd))
  CurrentStatusValues = CurrentStatusValues[CurrentStatusValues$ProjectionType=="RegularF",]
  CurrentStatusValues$ProjectionType=NULL
  Benchmarks_atSpr30 = CurrentStatusValues[CurrentStatusValues$FishingMortaltiyProjectionScenario=="F_at_SPR30",]
  Benchmarks_atSpr30$FishingMortaltiyProjectionScenario=NULL
  names(Benchmarks_atSpr30)[names(Benchmarks_atSpr30)=="F_project"]="F_SPR30"
  names(Benchmarks_atSpr30)[names(Benchmarks_atSpr30)=="BiomassAtEnd"]="B_SPR30"
  names(Benchmarks_atSpr30)[names(Benchmarks_atSpr30)=="NumFishAtEnd"]="N_SPR30"
  names(Benchmarks_atSpr30)[names(Benchmarks_atSpr30)=="RecruitmentAtEnd"]="R_SPR30"
  names(Benchmarks_atSpr30)[names(Benchmarks_atSpr30)=="SSBatEnd"]="SSB_SPR30"
  names(Benchmarks_atSpr30)[names(Benchmarks_atSpr30)=="CatchAtEnd"]="Yield_SPR30"
  PlotsToMake = subset(SummaryProjectionScenarios,select=c(Genus,Species,WPP,FishingMortaltiyProjectionScenario))
  PlotsToMake = unique(PlotsToMake)
#  PlotsToMake = subset(PlotsToMake,PlotsToMake$WPP=="712")    #TEMPORARY LINE OF CODE TO GET JUST WHAT WE WANT TO PLOT - NOT EVERYTHING YET
  for(i in 1:dim(PlotsToMake)[1])
  {
    ScenariosToPlot = SummaryProjectionScenarios[SummaryProjectionScenarios$Genus==PlotsToMake$Genus[i] &
      SummaryProjectionScenarios$Species==PlotsToMake$Species[i] & SummaryProjectionScenarios$WPP==PlotsToMake$WPP[i] &
      SummaryProjectionScenarios$FishingMortaltiyProjectionScenario==PlotsToMake$FishingMortaltiyProjectionScenario[i],]
    PopReconstructionLineType = data.frame(populationReconstructionMethod=sort(unique(as.character(ScenariosToPlot$populationReconstructionMethod))),LineType=seq(1,length(unique(ScenariosToPlot$populationReconstructionMethod)),by=1))
    ScenariosToPlot = merge(ScenariosToPlot,PopReconstructionLineType,all.x=TRUE,by=c("populationReconstructionMethod"))
    ScenariosToPlot$Scenarios_Colors = paste(ScenariosToPlot$Scenario_New,ScenariosToPlot$SteepnessAssumed,sep="-")
    Scenario_Colors = data.frame(CatchMSY_Scenario=unique(ScenariosToPlot$CatchMSY_Scenario),Colors=ColorsFor_ProjectionF_ScenarioPlots[1:length(unique(ScenariosToPlot$CatchMSY_Scenario))])
    ScenariosToPlot =  merge(ScenariosToPlot,Scenario_Colors,all.x=TRUE,by=c("CatchMSY_Scenario"))
    #Don't include scenarios that won't rebuild. I record which of these don't rebuild in the non-standardized fan plot code area
    ScenariosToPlot = ScenariosToPlot[ScenariosToPlot$ProjectionType!="RegularF_willNotRebuildAtF=0",]
    ScenariosToPlot = merge(ScenariosToPlot,Benchmarks_atSpr30,all.x=TRUE,by=c("Genus","Species","WPP","CatchMSY_Scenario","CatchMSY_Scenario_Bound","LifeHistory_Scenario","populationReconstructionMethod","SteepnessAssumed"))
    if(dim(ScenariosToPlot)[1]==0)
    {
      next
    }
    for(j in 1:dim(ScenariosToPlot)[1])
    {
      if(j==1)
      {
        Yield = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$ProjectedCatchYear)
        Yield = Yield/ScenariosToPlot$Yield_SPR30[j]
        names(Yield) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        SSB = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$SSB)
        SSB = SSB/ScenariosToPlot$SSB_SPR30[j]
        names(SSB) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Abundance = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$ProjectedPopulationYear)
        Abundance = Abundance/ScenariosToPlot$N_SPR30[j]
        names(Abundance) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Recruitment = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$Recruitment)
        Recruitment = Recruitment/ScenariosToPlot$R_SPR30[j]
        names(Recruitment) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Biomass = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$ProjectedBiomassYear)
        Biomass = Biomass/ScenariosToPlot$B_SPR30[j]
        names(Biomass) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
      }
      if(j>1)
      {
        Yield_j = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$ProjectedCatchYear)
        Yield_j = Yield_j/ScenariosToPlot$Yield_SPR30[j]
        names(Yield_j) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Yield = cbind(Yield,Yield_j)
        SSB_j = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$SSB)
        SSB_j = SSB_j/ScenariosToPlot$SSB_SPR30[j]
        names(SSB_j) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        SSB = cbind(SSB,SSB_j)
        Abundance_j = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$ProjectedPopulationYear)
        Abundance_j = Abundance_j/ScenariosToPlot$N_SPR30[j]
        names(Abundance_j) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Abundance = cbind(Abundance,Abundance_j)
        Recruitment_j = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$Recruitment)
        Recruitment_j = Recruitment_j/ScenariosToPlot$R_SPR30[j]
        names(Recruitment_j) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Recruitment = cbind(Recruitment,Recruitment_j)
        Biomass_j = data.frame(ListOfPopulations_WithProjections[[ScenariosToPlot$ListIndexNumber[j]]]$ProjectedBiomassYear)
        Biomass_j = Biomass_j/ScenariosToPlot$B_SPR30[j]
        names(Biomass_j) = paste(ScenariosToPlot$LifeHistory_Scenario[j],ScenariosToPlot$populationReconstructionMethod[j],ScenariosToPlot$SteepnessAssumed[j],sep="-")
        Biomass = cbind(Biomass,Biomass_j)
      }
    }
    Yield$index = seq(1,dim(Yield)[1],by=1)
    Yield = Yield[Yield$index<projectionYearsToPlot,]
    Yield$index=NULL
    SSB$index = seq(1,dim(SSB)[1],by=1)
    SSB = SSB[SSB$index<projectionYearsToPlot,]
    SSB$index=NULL
    Abundance$index = seq(1,dim(Abundance)[1],by=1)
    Abundance = Abundance[Abundance$index<projectionYearsToPlot,]
    Abundance$index=NULL
    Recruitment$index = seq(1,dim(Recruitment)[1],by=1)
    Recruitment = Recruitment[Recruitment$index<projectionYearsToPlot,]
    Recruitment$index=NULL
    Biomass$index = seq(1,dim(Biomass)[1],by=1)
    Biomass = Biomass[Biomass$index<projectionYearsToPlot,]
    Biomass$index=NULL
    #Yield Plots
    FolderNameSpp = paste(as.character(ScenariosToPlot$Genus[1]),CapStr(as.character(ScenariosToPlot$Species[1])),sep="")
    png(paste(paste(PATH_standardizedFanProjectionPlots,FolderNameSpp,sep="/"),paste(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],
       ScenariosToPlot$WPP[1],ScenariosToPlot$CatchMSY_Scenario[1],ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep="_"),"_Fan_Standardized_Yield.png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
    plot(NULL,xlim = c(min(as.numeric(row.names(Yield))),max(as.numeric(row.names(Yield)))),ylim = c(min(Yield),max(Yield)),ylab="Yield/Yield@SPR30",xlab="Year")
    Yield_thisReconstruction_m = t(Yield)
    fan(data = Yield_thisReconstruction_m,start=2019,type="percentile",fan.col=topo.colors,ln=NULL,rlab=NULL)
    abline(h=1,col="red",lwd=2,lty=2)
    title(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],sep=" "))
    mtext(paste("WPP=",ScenariosToPlot$WPP[1]," | Projection Type=",ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep=""))
    dev.off()
    #SSB Plots
    FolderNameSpp = paste(as.character(ScenariosToPlot$Genus[1]),CapStr(as.character(ScenariosToPlot$Species[1])),sep="")
    png(paste(paste(PATH_standardizedFanProjectionPlots,FolderNameSpp,sep="/"),paste(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],
       ScenariosToPlot$WPP[1],ScenariosToPlot$CatchMSY_Scenario[1],ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep="_"),"_Fan_Standardized_SSB.png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
    plot(NULL,xlim = c(min(as.numeric(row.names(SSB))),max(as.numeric(row.names(SSB)))),ylim = c(min(SSB),max(SSB)),ylab="SSB/SSB@SPR30",xlab="Year")
    SSB_thisReconstruction_m = t(SSB)
    fan(data = SSB_thisReconstruction_m,start=2019,type="percentile",fan.col=topo.colors,ln=NULL,rlab=NULL)
    abline(h=1,col="red",lwd=2,lty=2)
    title(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],sep=" "))
    mtext(paste("WPP=",ScenariosToPlot$WPP[1]," | Projection Type=",ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep=""))
    dev.off()
    #Abundance Plots
    FolderNameSpp = paste(as.character(ScenariosToPlot$Genus[1]),CapStr(as.character(ScenariosToPlot$Species[1])),sep="")
    png(paste(paste(PATH_standardizedFanProjectionPlots,FolderNameSpp,sep="/"),paste(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],
       ScenariosToPlot$WPP[1],ScenariosToPlot$CatchMSY_Scenario[1],ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep="_"),"_Fan_Standardized_Abundance.png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
    plot(NULL,xlim = c(min(as.numeric(row.names(Abundance))),max(as.numeric(row.names(Abundance)))),ylim = c(min(Abundance),max(Abundance)),ylab="Abundance/Abundance@SPR30",xlab="Year")
    Abundance_thisReconstruction_m = t(Abundance)
    fan(data = Abundance_thisReconstruction_m,start=2019,type="percentile",fan.col=topo.colors,ln=NULL,rlab=NULL)
    abline(h=1,col="red",lwd=2,lty=2)
    title(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],sep=" "))
    mtext(paste("WPP=",ScenariosToPlot$WPP[1]," | Projection Type=",ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep=""))
    dev.off()
    #Recruitment Plots
    FolderNameSpp = paste(as.character(ScenariosToPlot$Genus[1]),CapStr(as.character(ScenariosToPlot$Species[1])),sep="")
    png(paste(paste(PATH_standardizedFanProjectionPlots,FolderNameSpp,sep="/"),paste(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],
       ScenariosToPlot$WPP[1],ScenariosToPlot$CatchMSY_Scenario[1],ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep="_"),"_Fan_Standardized_Recruitment.png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
    plot(NULL,xlim = c(min(as.numeric(row.names(Recruitment))),max(as.numeric(row.names(Recruitment)))),ylim = c(min(Recruitment),max(Recruitment)),ylab="Recruitment/Recruitment@SPR30",xlab="Year")
    Recruitment_thisReconstruction_m = t(Recruitment)
    fan(data = Recruitment_thisReconstruction_m,start=2019,type="percentile",fan.col=topo.colors,ln=NULL,rlab=NULL)
    abline(h=1,col="red",lwd=2,lty=2)
    title(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],sep=" "))
    mtext(paste("WPP=",ScenariosToPlot$WPP[1]," | Projection Type=",ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep=""))
    dev.off()
    #Biomass Plots
    FolderNameSpp = paste(as.character(ScenariosToPlot$Genus[1]),CapStr(as.character(ScenariosToPlot$Species[1])),sep="")
    png(paste(paste(PATH_standardizedFanProjectionPlots,FolderNameSpp,sep="/"),paste(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],
       ScenariosToPlot$WPP[1],ScenariosToPlot$CatchMSY_Scenario[1],ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep="_"),"_Fan_Standardized_Biomass.png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
    plot(NULL,xlim = c(min(as.numeric(row.names(Biomass))),max(as.numeric(row.names(Biomass)))),ylim = c(min(Biomass),max(Biomass)),ylab="Biomass/Biomass@SPR30",xlab="Year")
    Biomass_thisReconstruction_m = t(Biomass)
    fan(data = Biomass_thisReconstruction_m,start=2019,type="percentile",fan.col=topo.colors,ln=NULL,rlab=NULL)
    abline(h=1,col="red",lwd=2,lty=2)
    title(paste(ScenariosToPlot$Genus[1],ScenariosToPlot$Species[1],sep=" "))
    mtext(paste("WPP=",ScenariosToPlot$WPP[1]," | Projection Type=",ScenariosToPlot$FishingMortaltiyProjectionScenario[1],sep=""))
    dev.off()
    if(i%%20==0)
    {
      print(paste("Plotting Standardized Fan Plot Scenario ",i," of ",dim(PlotsToMake)[1],".......",sep=""))
      flush.console()
    }
  }
}





#######################################################################################################
################################## Get Specific Scenarios Out for Delivery ############################
#######################################################################################################

#SummaryProjectionScenarios = read.table(paste(PATH_output,"SummaryProjectionScenarios.csv",sep="/"),header=TRUE,sep=",",as.is=TRUE)
#SummaryProjectionScenarios_712 = subset(SummaryProjectionScenarios,SummaryProjectionScenarios$WPP=="712")
#SummaryProjectionScenarios_712 = subset(SummaryProjectionScenarios_712,SummaryProjectionScenarios_712$populationReconstructionMethod!="TropFishR")
#SummaryProjectionScenarios_712 = subset(SummaryProjectionScenarios_712,SummaryProjectionScenarios_712$populationReconstructionMethod!="SolveBaranov")
#SummaryProjectionScenarios_712 = subset(SummaryProjectionScenarios_712,SummaryProjectionScenarios_712$populationReconstructionMethod!="LIME")
#SummaryProjectionScenarios_712_currentF = subset(SummaryProjectionScenarios_712,SummaryProjectionScenarios_712$FishingMortaltiyProjectionScenario=="CurrentF")
#write.table(SummaryProjectionScenarios_712_currentF,paste(PATH_output,"SummaryProjectionScenarios_712_currentF_Poseidon.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
#
#
#Lengths_cpue_712_Erythropterus = subset(Lengths_cpue,Lengths_cpue$fish_id=="504404580790866839" & Lengths_cpue$wpp1==712)
#hist(Lengths_cpue_712_Erythropterus$cm)


#####################################################################################################
######################################### PLOTS FOR PAPER!!!!!!! ####################################
#####################################################################################################
#NOTE: Here, use the output from the run with 50 Species
PATH_output_old = PATH_output
PATH_output = paste(PATH_output_old,"Run50Species",sep="/")

PATH_highLevelSummaryPlots_old = PATH_highLevelSummaryPlots
PATH_highLevelSummaryPlots = paste(PATH_output,"HighLevelSummaryPlots",sep="/")

################# Maps with Percentage of overfished and overfishing Pie Charts ########################
StatusTable = read.table(paste(PATH_output,"StatusTable.csv",sep="/"),header=TRUE,sep=",")
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
CurrentStatusValues = read.table(paste(PATH_output,"CurrentValuesBenchmarksAndRatios.csv",sep="/"),header=TRUE,sep=",",as.is=TRUE)
#TEMP CODE LINE: do NOT include SolveBaranov or TropFishR
CurrentStatusValues = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="SolveBaranov",]   #TEMP CODE LINE: do NOT include SolveBaranov or TropFishR
CurrentStatusValues = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="TropFishR",]      #TEMP CODE LINE: do NOT include SolveBaranov or TropFishR
#CurrentStatusValues = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="LIME",]
#Identify and filter away the Runs that did not converge
CurrentStatusValues$F_ratio20_FLAG=0
CurrentStatusValues$F_ratio20_FLAG[CurrentStatusValues$F_ratio20 > 5 | CurrentStatusValues$F_ratio20 < 0.01]=1
CurrentStatusValues$F_ratio30_FLAG=0
CurrentStatusValues$F_ratio30_FLAG[CurrentStatusValues$F_ratio30 > 5 | CurrentStatusValues$F_ratio30 < 0.01]=1
CurrentStatusValues$F_ratio40_FLAG=0
CurrentStatusValues$F_ratio40_FLAG[CurrentStatusValues$F_ratio40 > 5 | CurrentStatusValues$F_ratio40 < 0.01]=1
CurrentStatusValues$B_ratio20_FLAG=0
CurrentStatusValues$B_ratio20_FLAG[CurrentStatusValues$B_ratio20 > 5 | CurrentStatusValues$B_ratio20 < 0.01]=1
CurrentStatusValues$B_ratio30_FLAG=0
CurrentStatusValues$B_ratio30_FLAG[CurrentStatusValues$B_ratio30 > 5 | CurrentStatusValues$B_ratio30 < 0.01]=1
CurrentStatusValues$B_ratio40_FLAG=0
CurrentStatusValues$B_ratio40_FLAG[CurrentStatusValues$B_ratio40 > 5 | CurrentStatusValues$B_ratio40 < 0.01]=1
#Combine genus_species name
CurrentStatusValues$GenusSpecies = paste(CurrentStatusValues$Genus,CurrentStatusValues$Species,sep="_")
#Setup count of which are overfishd and which not to create Kobe plot
CurrentStatusValues$F_ratioAboveOneSPR30 = 0
CurrentStatusValues$F_ratioAboveOneSPR30[CurrentStatusValues$F_ratio30>=1]=1
CurrentStatusValues$B_ratioBelowOneSPR30 = 0
CurrentStatusValues$B_ratioBelowOneSPR30[CurrentStatusValues$B_ratio30<1]=1
CurrentStatusValues$F_ratioAboveOneSPR20 = 0
CurrentStatusValues$F_ratioAboveOneSPR20[CurrentStatusValues$F_ratio20>=1]=1
CurrentStatusValues$B_ratioBelowOneSPR20 = 0
CurrentStatusValues$B_ratioBelowOneSPR20[CurrentStatusValues$B_ratio20<1]=1
CurrentStatusValues$F_ratioAboveOneSPR40 = 0
CurrentStatusValues$F_ratioAboveOneSPR40[CurrentStatusValues$F_ratio40>=1]=1
CurrentStatusValues$B_ratioBelowOneSPR40 = 0
CurrentStatusValues$B_ratioBelowOneSPR40[CurrentStatusValues$B_ratio40<1]=1
#F spr20
CurrentStatusValues_Fratio20 = subset(CurrentStatusValues,CurrentStatusValues$F_ratio20_FLAG==0)
Status_WPP_Fratio20 = aggregate.data.frame(list(CurrentStatusValues_Fratio20$F_ratioAboveOneSPR20),by=list(CurrentStatusValues_Fratio20$WPP),FUN=sum)
names(Status_WPP_Fratio20) = c("WPP","Num_Fratio20AboveOne")
Status_WPP_Fratio20_TOTAL = data.frame(table(CurrentStatusValues_Fratio20$WPP))
names(Status_WPP_Fratio20_TOTAL) = c("WPP","Total_Fratio20")
Status_WPP_Fratio20 = merge(Status_WPP_Fratio20,Status_WPP_Fratio20_TOTAL,all.x=TRUE,by=c("WPP"))
Status_WPP_Fratio20$Percent_OverfishingSPR20 = Status_WPP_Fratio20$Num_Fratio20AboveOne/Status_WPP_Fratio20$Total_Fratio20
Status_WPP_Fratio20$Num_Fratio20AboveOne=NULL
Status_WPP_Fratio20$Total_Fratio20=NULL
#F spr30
CurrentStatusValues_Fratio30 = subset(CurrentStatusValues,CurrentStatusValues$F_ratio30_FLAG==0)
Status_WPP_Fratio30 = aggregate.data.frame(list(CurrentStatusValues_Fratio30$F_ratioAboveOneSPR30),by=list(CurrentStatusValues_Fratio30$WPP),FUN=sum)
names(Status_WPP_Fratio30) = c("WPP","Num_Fratio30AboveOne")
Status_WPP_Fratio30_TOTAL = data.frame(table(CurrentStatusValues_Fratio30$WPP))
names(Status_WPP_Fratio30_TOTAL) = c("WPP","Total_Fratio30")
Status_WPP_Fratio30 = merge(Status_WPP_Fratio30,Status_WPP_Fratio30_TOTAL,all.x=TRUE,by=c("WPP"))
Status_WPP_Fratio30$Percent_OverfishingSPR30 = Status_WPP_Fratio30$Num_Fratio30AboveOne/Status_WPP_Fratio30$Total_Fratio30
Status_WPP_Fratio30$Num_Fratio30AboveOne=NULL
Status_WPP_Fratio30$Total_Fratio30=NULL
#F spr40
CurrentStatusValues_Fratio40 = subset(CurrentStatusValues,CurrentStatusValues$F_ratio40_FLAG==0)
Status_WPP_Fratio40 = aggregate.data.frame(list(CurrentStatusValues_Fratio40$F_ratioAboveOneSPR40),by=list(CurrentStatusValues_Fratio40$WPP),FUN=sum)
names(Status_WPP_Fratio40) = c("WPP","Num_Fratio40AboveOne")
Status_WPP_Fratio40_TOTAL = data.frame(table(CurrentStatusValues_Fratio40$WPP))
names(Status_WPP_Fratio40_TOTAL) = c("WPP","Total_Fratio40")
Status_WPP_Fratio40 = merge(Status_WPP_Fratio40,Status_WPP_Fratio40_TOTAL,all.x=TRUE,by=c("WPP"))
Status_WPP_Fratio40$Percent_OverfishingSPR40 = Status_WPP_Fratio40$Num_Fratio40AboveOne/Status_WPP_Fratio40$Total_Fratio40
Status_WPP_Fratio40$Num_Fratio40AboveOne=NULL
Status_WPP_Fratio40$Total_Fratio40=NULL
#B spr20
CurrentStatusValues_Bratio20 = subset(CurrentStatusValues,CurrentStatusValues$B_ratio20_FLAG==0)
Status_WPP_Bratio20 = aggregate.data.frame(list(CurrentStatusValues_Bratio20$B_ratioBelowOneSPR20),by=list(CurrentStatusValues_Bratio20$WPP),FUN=sum)
names(Status_WPP_Bratio20) = c("WPP","Num_Bratio20BelowOne")
Status_WPP_Bratio20_TOTAL = data.frame(table(CurrentStatusValues_Bratio20$WPP))
names(Status_WPP_Bratio20_TOTAL) = c("WPP","Total_Bratio20")
Status_WPP_Bratio20 = merge(Status_WPP_Bratio20,Status_WPP_Bratio20_TOTAL,all.x=TRUE,by=c("WPP"))
Status_WPP_Bratio20$Percent_OverfishedSPR20 = Status_WPP_Bratio20$Num_Bratio20BelowOne/Status_WPP_Bratio20$Total_Bratio20
Status_WPP_Bratio20$Num_Bratio20BelowOne=NULL
Status_WPP_Bratio20$Total_Bratio20=NULL
#B spr30
CurrentStatusValues_Bratio30 = subset(CurrentStatusValues,CurrentStatusValues$B_ratio30_FLAG==0)
Status_WPP_Bratio30 = aggregate.data.frame(list(CurrentStatusValues_Bratio30$B_ratioBelowOneSPR30),by=list(CurrentStatusValues_Bratio30$WPP),FUN=sum)
names(Status_WPP_Bratio30) = c("WPP","Num_Bratio30BelowOne")
Status_WPP_Bratio30_TOTAL = data.frame(table(CurrentStatusValues_Bratio30$WPP))
names(Status_WPP_Bratio30_TOTAL) = c("WPP","Total_Bratio30")
Status_WPP_Bratio30 = merge(Status_WPP_Bratio30,Status_WPP_Bratio30_TOTAL,all.x=TRUE,by=c("WPP"))
Status_WPP_Bratio30$Percent_OverfishedSPR30 = Status_WPP_Bratio30$Num_Bratio30BelowOne/Status_WPP_Bratio30$Total_Bratio30
Status_WPP_Bratio30$Num_Bratio30BelowOne=NULL
Status_WPP_Bratio30$Total_Bratio30=NULL
#B spr40
CurrentStatusValues_Bratio40 = subset(CurrentStatusValues,CurrentStatusValues$B_ratio40_FLAG==0)
Status_WPP_Bratio40 = aggregate.data.frame(list(CurrentStatusValues_Bratio40$B_ratioBelowOneSPR40),by=list(CurrentStatusValues_Bratio40$WPP),FUN=sum)
names(Status_WPP_Bratio40) = c("WPP","Num_Bratio40BelowOne")
Status_WPP_Bratio40_TOTAL = data.frame(table(CurrentStatusValues_Bratio40$WPP))
names(Status_WPP_Bratio40_TOTAL) = c("WPP","Total_Bratio40")
Status_WPP_Bratio40 = merge(Status_WPP_Bratio40,Status_WPP_Bratio40_TOTAL,all.x=TRUE,by=c("WPP"))
Status_WPP_Bratio40$Percent_OverfishedSPR40 = Status_WPP_Bratio40$Num_Bratio40BelowOne/Status_WPP_Bratio40$Total_Bratio40
Status_WPP_Bratio40$Num_Bratio40BelowOne=NULL
Status_WPP_Bratio40$Total_Bratio40=NULL
Status_WPP = merge(Status_WPP_Fratio20,Status_WPP_Fratio30,all.x=TRUE,by=c("WPP"))
Status_WPP = merge(Status_WPP,Status_WPP_Fratio40,all.x=TRUE,by=c("WPP"))
Status_WPP = merge(Status_WPP,Status_WPP_Bratio20,all.x=TRUE,by=c("WPP"))
Status_WPP = merge(Status_WPP,Status_WPP_Bratio30,all.x=TRUE,by=c("WPP"))
Status_WPP = merge(Status_WPP,Status_WPP_Bratio40,all.x=TRUE,by=c("WPP"))
Status_WPP = merge(Status_WPP,LocationsWPP,all.x=TRUE,by=c("WPP"))
Status_WPP = subset(Status_WPP,!is.na(Status_WPP$Longitude))
#Sample Sizes Determined
GenusSpeciesWPP = data.frame(GenusSpecies=CurrentStatusValues$GenusSpecies,WPP=CurrentStatusValues$WPP)
GenusSpeciesWPP = unique(GenusSpeciesWPP)
SampleSizes = data.frame(table(GenusSpeciesWPP$WPP))
names(SampleSizes) = c("WPP","N")
SampleSizes$WPP = as.numeric(as.character(SampleSizes$WPP))
SampleSizes = subset(SampleSizes,!is.na(SampleSizes$WPP))
Status_WPP = merge(Status_WPP,SampleSizes,all.x=TRUE,by=c("WPP"))

#Pie Chart Sample Size Text Locations for Each Pie
TextSampleSizeLocations = data.frame(x=c(100,100.46866,117.3,107.47709,113.94177,118.77517,126.68987,127.05237,122.58147,135.75250),
  y=c(7,-3.0658615,-8.3130337,4.2277078,-3.3282201,-0.8620492,-5.4795607,1.6467065,4.6474816,-7.9982033))

#Plot Overfished Status at SPR 30%
png(paste(PATH_highLevelSummaryPlots,"OverfishedMapSPR30.png",sep="/"),units="px",width=5200,height=3200,res=600)
par(mai = c(0.5, 0.5, 0.35, 0.2), omi = c(0.25, 0.5, 0, 0),mgp = c(2.5, 0.5, 0))
prettymap(plot(WPPshapeFile,col=c("green",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),xlab="Longitude",ylab="Latitude"), oma=c(4,4,1,0),drawscale=TRUE,drawbox=FALSE, scale.plotunit="latlon", drawarrow=TRUE)
axis(side=1,outer=TRUE,line=-2.5)
axis(side=2,outer=TRUE,line=-2)
mtext("Longitude",side=1,line=2,cex=1.2)
mtext("Latitude",side=2,line=2.25,cex=1.2)
title("Population Status: Overfished at SPR 30%")
mtext("Indonesia Deepwater Snapper-Grouper Fishery",font=3,line=-0.5)
formattedForPieCharts = subset(Status_WPP,select=c(Longitude,Latitude,Percent_OverfishedSPR30))
names(formattedForPieCharts)[names(formattedForPieCharts)=="Percent_OverfishedSPR30"]="Value"
formattedForPieCharts$Variable = "Overfished"
notOverfished = formattedForPieCharts
notOverfished$Variable = "NotOverfished"
notOverfished$Value = 1 - notOverfished$Value
formattedForPieCharts = rbind(formattedForPieCharts,notOverfished)
xyz = make.xyz(formattedForPieCharts$Longitude,formattedForPieCharts$Latitude,formattedForPieCharts$Value,formattedForPieCharts$Variable)
for(i in 1:(dim(Status_WPP)[1]))
{
  add.pie(xyz$z[i,],xyz$x[i],xyz$y[i],radius=1.2,col=c("blue","orange"),labels=NA)  #,labels=c(NA,paste("n=",Status_WPP$N[i],sep="")),label.dist=1.3)
}
legend.pie(99,-11,labels=c("Not Overfished","Overfished"),radius=1,bty="n",col=c("blue","orange"),cex=0.8,label.dist=1.3)
for(i in 1:(dim(Status_WPP)[1]))
{
  text(TextSampleSizeLocations$x[i],TextSampleSizeLocations$y[i],paste("n=",Status_WPP$N[i],sep=""))
}
dev.off()

#Plot Overfishing Status at SPR 30%
png(paste(PATH_highLevelSummaryPlots,"OverfishingMapSPR30.png",sep="/"),units="px",width=5200,height=3200,res=600)
par(mai = c(0.5, 0.5, 0.35, 0.2), omi = c(0.25, 0.5, 0, 0),mgp = c(2.5, 0.5, 0))
prettymap(plot(WPPshapeFile,col=c("green",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),xlab="Longitude",ylab="Latitude"), oma=c(4,4,1,0),drawscale=TRUE,drawbox=FALSE, scale.plotunit="latlon", drawarrow=TRUE)
axis(side=1,outer=TRUE,line=-2.5)
axis(side=2,outer=TRUE,line=-2)
mtext("Longitude",side=1,line=2,cex=1.2)
mtext("Latitude",side=2,line=2.25,cex=1.2)
title("Population Status: Overfishing at SPR 30%")
mtext("Indonesia Deepwater Snapper-Grouper Fishery",font=3,line=-0.5)
formattedForPieCharts = subset(Status_WPP,select=c(Longitude,Latitude,Percent_OverfishingSPR30))
names(formattedForPieCharts)[names(formattedForPieCharts)=="Percent_OverfishingSPR30"]="Value"
formattedForPieCharts$Variable = "Overfishing"
notOverfishing = formattedForPieCharts
notOverfishing$Variable = "NotOverfishing"
notOverfishing$Value = 1 - notOverfished$Value
formattedForPieCharts = rbind(formattedForPieCharts,notOverfishing)
xyz = make.xyz(formattedForPieCharts$Longitude,formattedForPieCharts$Latitude,formattedForPieCharts$Value,formattedForPieCharts$Variable)
for(i in 1:(dim(Status_WPP)[1]))
{
  add.pie(xyz$z[i,],xyz$x[i],xyz$y[i],radius=1.2,col=c("blue","orange"),labels=NA)  #,labels=c(NA,paste("n=",Status_WPP$N[i],sep="")),label.dist=1.3)
}
legend.pie(99,-11,labels=c("Not Overfishing","Overfishing"),radius=1,bty="n",col=c("blue","orange"),cex=0.8,label.dist=1.3)
for(i in 1:(dim(Status_WPP)[1]))
{
  text(TextSampleSizeLocations$x[i],TextSampleSizeLocations$y[i],paste("n=",Status_WPP$N[i],sep=""))
}
dev.off()

#Plot Overfished Status at SPR 40%
png(paste(PATH_highLevelSummaryPlots,"OverfishedMapSPR40.png",sep="/"),units="px",width=5200,height=3200,res=600)
par(mai = c(0.5, 0.5, 0.35, 0.2), omi = c(0.25, 0.5, 0, 0),mgp = c(2.5, 0.5, 0))
prettymap(plot(WPPshapeFile,col=c("green",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),xlab="Longitude",ylab="Latitude"), oma=c(4,4,1,0),drawscale=TRUE,drawbox=FALSE, scale.plotunit="latlon", drawarrow=TRUE)
axis(side=1,outer=TRUE,line=-2.5)
axis(side=2,outer=TRUE,line=-2)
mtext("Longitude",side=1,line=2,cex=1.2)
mtext("Latitude",side=2,line=2.25,cex=1.2)
title("Population Status: Overfished at SPR 40%")
mtext("Indonesia Deepwater Snapper-Grouper Fishery",font=3,line=-0.5)
formattedForPieCharts = subset(Status_WPP,select=c(Longitude,Latitude,Percent_OverfishedSPR40))
names(formattedForPieCharts)[names(formattedForPieCharts)=="Percent_OverfishedSPR40"]="Value"
formattedForPieCharts$Variable = "Overfished"
notOverfished = formattedForPieCharts
notOverfished$Variable = "NotOverfished"
notOverfished$Value = 1 - notOverfished$Value
formattedForPieCharts = rbind(formattedForPieCharts,notOverfished)
xyz = make.xyz(formattedForPieCharts$Longitude,formattedForPieCharts$Latitude,formattedForPieCharts$Value,formattedForPieCharts$Variable)
for(i in 1:(dim(Status_WPP)[1]))
{
  add.pie(xyz$z[i,],xyz$x[i],xyz$y[i],radius=1.2,col=c("blue","orange"),labels=NA)  #,labels=c(NA,paste("n=",Status_WPP$N[i],sep="")),label.dist=1.3)
}
legend.pie(99,-11,labels=c("Not Overfished","Overfished"),radius=1,bty="n",col=c("blue","orange"),cex=0.8,label.dist=1.3)
for(i in 1:(dim(Status_WPP)[1]))
{
  text(TextSampleSizeLocations$x[i],TextSampleSizeLocations$y[i],paste("n=",Status_WPP$N[i],sep=""))
}
dev.off()

#Plot Overfishing Status at SPR 40%
png(paste(PATH_highLevelSummaryPlots,"OverfishingMapSPR40.png",sep="/"),units="px",width=5200,height=3200,res=600)
par(mai = c(0.5, 0.5, 0.35, 0.2), omi = c(0.25, 0.5, 0, 0),mgp = c(2.5, 0.5, 0))
prettymap(plot(WPPshapeFile,col=c("green",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),xlab="Longitude",ylab="Latitude"), oma=c(4,4,1,0),drawscale=TRUE,drawbox=FALSE, scale.plotunit="latlon", drawarrow=TRUE)
axis(side=1,outer=TRUE,line=-2.5)
axis(side=2,outer=TRUE,line=-2)
mtext("Longitude",side=1,line=2,cex=1.2)
mtext("Latitude",side=2,line=2.25,cex=1.2)
title("Population Status: Overfishing at SPR 40%")
mtext("Indonesia Deepwater Snapper-Grouper Fishery",font=3,line=-0.5)
formattedForPieCharts = subset(Status_WPP,select=c(Longitude,Latitude,Percent_OverfishingSPR40))
names(formattedForPieCharts)[names(formattedForPieCharts)=="Percent_OverfishingSPR40"]="Value"
formattedForPieCharts$Variable = "Overfishing"
notOverfishing = formattedForPieCharts
notOverfishing$Variable = "NotOverfishing"
notOverfishing$Value = 1 - notOverfished$Value
formattedForPieCharts = rbind(formattedForPieCharts,notOverfishing)
xyz = make.xyz(formattedForPieCharts$Longitude,formattedForPieCharts$Latitude,formattedForPieCharts$Value,formattedForPieCharts$Variable)
for(i in 1:(dim(Status_WPP)[1]))
{
  add.pie(xyz$z[i,],xyz$x[i],xyz$y[i],radius=1.2,col=c("blue","orange"),labels=NA)  #,labels=c(NA,paste("n=",Status_WPP$N[i],sep="")),label.dist=1.3)
}
legend.pie(99,-11,labels=c("Not Overfishing","Overfishing"),radius=1,bty="n",col=c("blue","orange"),cex=0.8,label.dist=1.3)
for(i in 1:(dim(Status_WPP)[1]))
{
  text(TextSampleSizeLocations$x[i],TextSampleSizeLocations$y[i],paste("n=",Status_WPP$N[i],sep=""))
}
dev.off()

#Plot Overfished Status at SPR 20%
png(paste(PATH_highLevelSummaryPlots,"OverfishedMapSPR20.png",sep="/"),units="px",width=5200,height=3200,res=600)
par(mai = c(0.5, 0.5, 0.35, 0.2), omi = c(0.25, 0.5, 0, 0),mgp = c(2.5, 0.5, 0))
prettymap(plot(WPPshapeFile,col=c("green",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),xlab="Longitude",ylab="Latitude"), oma=c(4,4,1,0),drawscale=TRUE,drawbox=FALSE, scale.plotunit="latlon", drawarrow=TRUE)
axis(side=1,outer=TRUE,line=-2.5)
axis(side=2,outer=TRUE,line=-2)
mtext("Longitude",side=1,line=2,cex=1.2)
mtext("Latitude",side=2,line=2.25,cex=1.2)
title("Population Status: Overfished at SPR 20%")
mtext("Indonesia Deepwater Snapper-Grouper Fishery",font=3,line=-0.5)
formattedForPieCharts = subset(Status_WPP,select=c(Longitude,Latitude,Percent_OverfishedSPR20))
names(formattedForPieCharts)[names(formattedForPieCharts)=="Percent_OverfishedSPR20"]="Value"
formattedForPieCharts$Variable = "Overfished"
notOverfished = formattedForPieCharts
notOverfished$Variable = "NotOverfished"
notOverfished$Value = 1 - notOverfished$Value
formattedForPieCharts = rbind(formattedForPieCharts,notOverfished)
xyz = make.xyz(formattedForPieCharts$Longitude,formattedForPieCharts$Latitude,formattedForPieCharts$Value,formattedForPieCharts$Variable)
for(i in 1:(dim(Status_WPP)[1]))
{
  add.pie(xyz$z[i,],xyz$x[i],xyz$y[i],radius=1.2,col=c("blue","orange"),labels=NA)  #,labels=c(NA,paste("n=",Status_WPP$N[i],sep="")),label.dist=1.3)
}
legend.pie(99,-11,labels=c("Not Overfished","Overfished"),radius=1,bty="n",col=c("blue","orange"),cex=0.8,label.dist=1.3)
for(i in 1:(dim(Status_WPP)[1]))
{
  text(TextSampleSizeLocations$x[i],TextSampleSizeLocations$y[i],paste("n=",Status_WPP$N[i],sep=""))
}
dev.off()

#Plot Overfishing Status at SPR 20%
png(paste(PATH_highLevelSummaryPlots,"OverfishingMapSPR20.png",sep="/"),units="px",width=5200,height=3200,res=600)
par(mai = c(0.5, 0.5, 0.35, 0.2), omi = c(0.25, 0.5, 0, 0),mgp = c(2.5, 0.5, 0))
prettymap(plot(WPPshapeFile,col=c("green",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),xlab="Longitude",ylab="Latitude"), oma=c(4,4,1,0),drawscale=TRUE,drawbox=FALSE, scale.plotunit="latlon", drawarrow=TRUE)
axis(side=1,outer=TRUE,line=-2.5)
axis(side=2,outer=TRUE,line=-2)
mtext("Longitude",side=1,line=2,cex=1.2)
mtext("Latitude",side=2,line=2.25,cex=1.2)
title("Population Status: Overfishing at SPR 20%")
mtext("Indonesia Deepwater Snapper-Grouper Fishery",font=3,line=-0.5)
formattedForPieCharts = subset(Status_WPP,select=c(Longitude,Latitude,Percent_OverfishingSPR20))
names(formattedForPieCharts)[names(formattedForPieCharts)=="Percent_OverfishingSPR20"]="Value"
formattedForPieCharts$Variable = "Overfishing"
notOverfishing = formattedForPieCharts
notOverfishing$Variable = "NotOverfishing"
notOverfishing$Value = 1 - notOverfished$Value
formattedForPieCharts = rbind(formattedForPieCharts,notOverfishing)
xyz = make.xyz(formattedForPieCharts$Longitude,formattedForPieCharts$Latitude,formattedForPieCharts$Value,formattedForPieCharts$Variable)
for(i in 1:(dim(Status_WPP)[1]))
{
  add.pie(xyz$z[i,],xyz$x[i],xyz$y[i],radius=1.2,col=c("blue","orange"),labels=NA)  #,labels=c(NA,paste("n=",Status_WPP$N[i],sep="")),label.dist=1.3)
}
legend.pie(99,-11,labels=c("Not Overfishing","Overfishing"),radius=1,bty="n",col=c("blue","orange"),cex=0.8,label.dist=1.3)
for(i in 1:(dim(Status_WPP)[1]))
{
  text(TextSampleSizeLocations$x[i],TextSampleSizeLocations$y[i],paste("n=",Status_WPP$N[i],sep=""))
}
dev.off()


########### Three Pie Partitions: Second Version of Maps with Percentage of overfished and overfishing Pie Charts ####################
StatusTable = read.table(paste(PATH_output,"StatusTable.csv",sep="/"),header=TRUE,sep=",")
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
CurrentStatusValues = read.table(paste(PATH_output,"CurrentValuesBenchmarksAndRatios.csv",sep="/"),header=TRUE,sep=",",as.is=TRUE)
CurrentStatusValues$GenusSpecies = paste(CurrentStatusValues$Genus,CurrentStatusValues$Species,sep="_")

#TEMP CODE LINE: do NOT include SolveBaranov or TropFishR in both F and B ratio calculations; do NOT include LIME in B ratio calculations!!
CurrentStatusValues = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="SolveBaranov",]   #TEMP CODE LINE: do NOT include SolveBaranov or TropFishR
CurrentStatusValues = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="TropFishR",]      #TEMP CODE LINE: do NOT include SolveBaranov or TropFishR
CurrentStatusValues_noLIME = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="LIME",]
#Identify and filter away the Runs that did not converge
CurrentStatusValues$F_ratio20_FLAG=0
CurrentStatusValues$F_ratio20_FLAG[CurrentStatusValues$F_ratio20 > 5 | CurrentStatusValues$F_ratio20 < 0.01]=1
CurrentStatusValues$F_ratio30_FLAG=0
CurrentStatusValues$F_ratio30_FLAG[CurrentStatusValues$F_ratio30 > 5 | CurrentStatusValues$F_ratio30 < 0.01]=1
CurrentStatusValues$F_ratio40_FLAG=0
CurrentStatusValues$F_ratio40_FLAG[CurrentStatusValues$F_ratio40 > 5 | CurrentStatusValues$F_ratio40 < 0.01]=1
CurrentStatusValues$B_ratio20_FLAG=0
CurrentStatusValues$B_ratio20_FLAG[CurrentStatusValues$B_ratio20 > 5 | CurrentStatusValues$B_ratio20 < 0.01]=1
CurrentStatusValues$B_ratio30_FLAG=0
CurrentStatusValues$B_ratio30_FLAG[CurrentStatusValues$B_ratio30 > 5 | CurrentStatusValues$B_ratio30 < 0.01]=1
CurrentStatusValues$B_ratio40_FLAG=0
CurrentStatusValues$B_ratio40_FLAG[CurrentStatusValues$B_ratio40 > 5 | CurrentStatusValues$B_ratio40 < 0.01]=1
CurrentStatusValues_noLIME$F_ratio20_FLAG=0
CurrentStatusValues_noLIME$F_ratio20_FLAG[CurrentStatusValues_noLIME$F_ratio20 > 5 | CurrentStatusValues_noLIME$F_ratio20 < 0.01]=1
CurrentStatusValues_noLIME$F_ratio30_FLAG=0
CurrentStatusValues_noLIME$F_ratio30_FLAG[CurrentStatusValues_noLIME$F_ratio30 > 5 | CurrentStatusValues_noLIME$F_ratio30 < 0.01]=1
CurrentStatusValues_noLIME$F_ratio40_FLAG=0
CurrentStatusValues_noLIME$F_ratio40_FLAG[CurrentStatusValues_noLIME$F_ratio40 > 5 | CurrentStatusValues_noLIME$F_ratio40 < 0.01]=1
CurrentStatusValues_noLIME$B_ratio20_FLAG=0
CurrentStatusValues_noLIME$B_ratio20_FLAG[CurrentStatusValues_noLIME$B_ratio20 > 5 | CurrentStatusValues_noLIME$B_ratio20 < 0.01]=1
CurrentStatusValues_noLIME$B_ratio30_FLAG=0
CurrentStatusValues_noLIME$B_ratio30_FLAG[CurrentStatusValues_noLIME$B_ratio30 > 5 | CurrentStatusValues_noLIME$B_ratio30 < 0.01]=1
CurrentStatusValues_noLIME$B_ratio40_FLAG=0
CurrentStatusValues_noLIME$B_ratio40_FLAG[CurrentStatusValues_noLIME$B_ratio40 > 5 | CurrentStatusValues_noLIME$B_ratio40 < 0.01]=1
CurrentStatusValues_F_ratio20 = subset(CurrentStatusValues,CurrentStatusValues$F_ratio20_FLAG==0)
CurrentStatusValues_F_ratio30 = subset(CurrentStatusValues,CurrentStatusValues$F_ratio30_FLAG==0)
CurrentStatusValues_F_ratio40 = subset(CurrentStatusValues,CurrentStatusValues$F_ratio40_FLAG==0)
CurrentStatusValues_B_ratio20 = subset(CurrentStatusValues,CurrentStatusValues$B_ratio20_FLAG==0)
CurrentStatusValues_B_ratio30 = subset(CurrentStatusValues,CurrentStatusValues$B_ratio30_FLAG==0)
CurrentStatusValues_B_ratio40 = subset(CurrentStatusValues,CurrentStatusValues$B_ratio40_FLAG==0)
CurrentStatusValues_noLIME_F_ratio20 = subset(CurrentStatusValues_noLIME,CurrentStatusValues_noLIME$F_ratio20_FLAG==0)
CurrentStatusValues_noLIME_F_ratio30 = subset(CurrentStatusValues_noLIME,CurrentStatusValues_noLIME$F_ratio30_FLAG==0)
CurrentStatusValues_noLIME_F_ratio40 = subset(CurrentStatusValues_noLIME,CurrentStatusValues_noLIME$F_ratio40_FLAG==0)
CurrentStatusValues_noLIME_B_ratio20 = subset(CurrentStatusValues_noLIME,CurrentStatusValues_noLIME$B_ratio20_FLAG==0)
CurrentStatusValues_noLIME_B_ratio30 = subset(CurrentStatusValues_noLIME,CurrentStatusValues_noLIME$B_ratio30_FLAG==0)
CurrentStatusValues_noLIME_B_ratio40 = subset(CurrentStatusValues_noLIME,CurrentStatusValues_noLIME$B_ratio40_FLAG==0)
AvgRatios_F_ratio20 = aggregate.data.frame(list(CurrentStatusValues_F_ratio20$F_ratio20),by=list(CurrentStatusValues_F_ratio20$Genus,CurrentStatusValues_F_ratio20$Species,CurrentStatusValues_F_ratio20$WPP),FUN=mean,na.rm=TRUE)
names(AvgRatios_F_ratio20) = c("Genus","Species","WPP","F_ratio20_mean")
StDevRatios_F_ratio20 = aggregate.data.frame(list(CurrentStatusValues_F_ratio20$F_ratio20),by=list(CurrentStatusValues_F_ratio20$Genus,CurrentStatusValues_F_ratio20$Species,CurrentStatusValues_F_ratio20$WPP),FUN=sd,na.rm=TRUE)
names(StDevRatios_F_ratio20) = c("Genus","Species","WPP","F_ratio20_sd")
AvgRatios_F_ratio30 = aggregate.data.frame(list(CurrentStatusValues_F_ratio30$F_ratio30),by=list(CurrentStatusValues_F_ratio30$Genus,CurrentStatusValues_F_ratio30$Species,CurrentStatusValues_F_ratio30$WPP),FUN=mean,na.rm=TRUE)
names(AvgRatios_F_ratio30) = c("Genus","Species","WPP","F_ratio30_mean")
StDevRatios_F_ratio30 = aggregate.data.frame(list(CurrentStatusValues_F_ratio30$F_ratio30),by=list(CurrentStatusValues_F_ratio30$Genus,CurrentStatusValues_F_ratio30$Species,CurrentStatusValues_F_ratio30$WPP),FUN=sd,na.rm=TRUE)
names(StDevRatios_F_ratio30) = c("Genus","Species","WPP","F_ratio30_sd")
AvgRatios_F_ratio40 = aggregate.data.frame(list(CurrentStatusValues_F_ratio40$F_ratio40),by=list(CurrentStatusValues_F_ratio40$Genus,CurrentStatusValues_F_ratio40$Species,CurrentStatusValues_F_ratio40$WPP),FUN=mean,na.rm=TRUE)
names(AvgRatios_F_ratio40) = c("Genus","Species","WPP","F_ratio40_mean")
StDevRatios_F_ratio40 = aggregate.data.frame(list(CurrentStatusValues_F_ratio40$F_ratio40),by=list(CurrentStatusValues_F_ratio40$Genus,CurrentStatusValues_F_ratio40$Species,CurrentStatusValues_F_ratio40$WPP),FUN=sd,na.rm=TRUE)
names(StDevRatios_F_ratio40) = c("Genus","Species","WPP","F_ratio40_sd")
AvgRatios_noLIME_B_ratio20 = aggregate.data.frame(list(CurrentStatusValues_noLIME_B_ratio20$B_ratio20),by=list(CurrentStatusValues_noLIME_B_ratio20$Genus,CurrentStatusValues_noLIME_B_ratio20$Species,CurrentStatusValues_noLIME_B_ratio20$WPP),FUN=mean,na.rm=TRUE)
names(AvgRatios_noLIME_B_ratio20) = c("Genus","Species","WPP","B_ratio20_mean")
StDevRatios_noLIME_B_ratio20 = aggregate.data.frame(list(CurrentStatusValues_noLIME_B_ratio20$B_ratio20),by=list(CurrentStatusValues_noLIME_B_ratio20$Genus,CurrentStatusValues_noLIME_B_ratio20$Species,CurrentStatusValues_noLIME_B_ratio20$WPP),FUN=sd,na.rm=TRUE)
names(StDevRatios_noLIME_B_ratio20) = c("Genus","Species","WPP","B_ratio20_sd")
AvgRatios_noLIME_B_ratio30 = aggregate.data.frame(list(CurrentStatusValues_noLIME_B_ratio30$B_ratio30),by=list(CurrentStatusValues_noLIME_B_ratio30$Genus,CurrentStatusValues_noLIME_B_ratio30$Species,CurrentStatusValues_noLIME_B_ratio30$WPP),FUN=mean,na.rm=TRUE)
names(AvgRatios_noLIME_B_ratio30) = c("Genus","Species","WPP","B_ratio30_mean")
StDevRatios_noLIME_B_ratio30 = aggregate.data.frame(list(CurrentStatusValues_noLIME_B_ratio30$B_ratio30),by=list(CurrentStatusValues_noLIME_B_ratio30$Genus,CurrentStatusValues_noLIME_B_ratio30$Species,CurrentStatusValues_noLIME_B_ratio30$WPP),FUN=sd,na.rm=TRUE)
names(StDevRatios_noLIME_B_ratio30) = c("Genus","Species","WPP","B_ratio30_sd")
AvgRatios_noLIME_B_ratio40 = aggregate.data.frame(list(CurrentStatusValues_noLIME_B_ratio40$B_ratio40),by=list(CurrentStatusValues_noLIME_B_ratio40$Genus,CurrentStatusValues_noLIME_B_ratio40$Species,CurrentStatusValues_noLIME_B_ratio40$WPP),FUN=mean,na.rm=TRUE)
names(AvgRatios_noLIME_B_ratio40) = c("Genus","Species","WPP","B_ratio40_mean")
StDevRatios_noLIME_B_ratio40 = aggregate.data.frame(list(CurrentStatusValues_noLIME_B_ratio40$B_ratio40),by=list(CurrentStatusValues_noLIME_B_ratio40$Genus,CurrentStatusValues_noLIME_B_ratio40$Species,CurrentStatusValues_noLIME_B_ratio40$WPP),FUN=sd,na.rm=TRUE)
names(StDevRatios_noLIME_B_ratio40) = c("Genus","Species","WPP","B_ratio40_sd")
AvgRatios = merge(AvgRatios_F_ratio20,StDevRatios_F_ratio20,all.x=TRUE,by=c("Genus","Species","WPP"))
AvgRatios = merge(AvgRatios,AvgRatios_F_ratio30,all.x=TRUE,by=c("Genus","Species","WPP"))
AvgRatios = merge(AvgRatios,StDevRatios_F_ratio30,all.x=TRUE,by=c("Genus","Species","WPP"))
AvgRatios = merge(AvgRatios,AvgRatios_F_ratio40,all.x=TRUE,by=c("Genus","Species","WPP"))
AvgRatios = merge(AvgRatios,StDevRatios_F_ratio40,all.x=TRUE,by=c("Genus","Species","WPP"))
AvgRatios = merge(AvgRatios,AvgRatios_noLIME_B_ratio20,all.x=TRUE,by=c("Genus","Species","WPP"))
AvgRatios = merge(AvgRatios,StDevRatios_noLIME_B_ratio20,all.x=TRUE,by=c("Genus","Species","WPP"))
AvgRatios = merge(AvgRatios,AvgRatios_noLIME_B_ratio30,all.x=TRUE,by=c("Genus","Species","WPP"))
AvgRatios = merge(AvgRatios,StDevRatios_noLIME_B_ratio30,all.x=TRUE,by=c("Genus","Species","WPP"))
AvgRatios = merge(AvgRatios,AvgRatios_noLIME_B_ratio40,all.x=TRUE,by=c("Genus","Species","WPP"))
AvgRatios = merge(AvgRatios,StDevRatios_noLIME_B_ratio40,all.x=TRUE,by=c("Genus","Species","WPP"))
#Calcualte bounds
AvgRatios$F_ratio30_Lower = AvgRatios$F_ratio30_mean - AvgRatios$F_ratio30_sd
AvgRatios$F_ratio30_Upper = AvgRatios$F_ratio30_mean + AvgRatios$F_ratio30_sd
AvgRatios$B_ratio30_Lower = AvgRatios$B_ratio30_mean - AvgRatios$B_ratio30_sd
AvgRatios$B_ratio30_Upper = AvgRatios$B_ratio30_mean + AvgRatios$B_ratio30_sd
AvgRatios$F_ratio30_Lower[AvgRatios$F_ratio30_Lower<0]=0
AvgRatios$B_ratio30_Lower[AvgRatios$B_ratio30_Lower<0]=0
AvgRatios$OverfishingCategorySPR30=0
AvgRatios$OverfishingCategorySPR30[AvgRatios$OverfishingCategorySPR30==0 & AvgRatios$F_ratio30_Upper<1 & AvgRatios$F_ratio30_Lower<1]=1
AvgRatios$OverfishingCategorySPR30[AvgRatios$OverfishingCategorySPR30==0 & AvgRatios$F_ratio30_Upper>=1 & AvgRatios$F_ratio30_Lower<1]=2
AvgRatios$OverfishingCategorySPR30[AvgRatios$OverfishingCategorySPR30==0 & AvgRatios$F_ratio30_Upper>1 & AvgRatios$F_ratio30_Lower>=1]=3
AvgRatios$OverfishedCategorySPR30=0
AvgRatios$OverfishedCategorySPR30[AvgRatios$OverfishedCategorySPR30==0 & AvgRatios$B_ratio30_Lower>1 & AvgRatios$B_ratio30_Upper>1]=1
AvgRatios$OverfishedCategorySPR30[AvgRatios$OverfishedCategorySPR30==0 & AvgRatios$B_ratio30_Upper>=1 & AvgRatios$B_ratio30_Lower<1]=2
AvgRatios$OverfishedCategorySPR30[AvgRatios$OverfishedCategorySPR30==0 & AvgRatios$B_ratio30_Upper<=1 & AvgRatios$B_ratio30_Lower<1]=3
AvgRatios$F_ratio20_Lower = AvgRatios$F_ratio20_mean - AvgRatios$F_ratio20_sd
AvgRatios$F_ratio20_Upper = AvgRatios$F_ratio20_mean + AvgRatios$F_ratio20_sd
AvgRatios$B_ratio20_Lower = AvgRatios$B_ratio20_mean - AvgRatios$B_ratio20_sd
AvgRatios$B_ratio20_Upper = AvgRatios$B_ratio20_mean + AvgRatios$B_ratio20_sd
AvgRatios$F_ratio20_Lower[AvgRatios$F_ratio20_Lower<0]=0
AvgRatios$B_ratio20_Lower[AvgRatios$B_ratio20_Lower<0]=0
#If you can not compute the F or B multiplier at SPR 20 (convergence issue), borrow from those values at SPR 30, at least for the trend and "color" - i.e. overfished or experiencing overfishing 
#NOTE: this does not change the tabulated output or the actual values of these in the Kobe plots. 
AvgRatios$F_ratio20_mean[is.na(AvgRatios$F_ratio20_mean)]=AvgRatios$F_ratio30_mean[is.na(AvgRatios$F_ratio20_mean)]
AvgRatios$F_ratio20_sd[is.na(AvgRatios$F_ratio20_sd)]=AvgRatios$F_ratio30_sd[is.na(AvgRatios$F_ratio20_sd)]
AvgRatios$F_ratio20_Upper[is.na(AvgRatios$F_ratio20_Upper)]=AvgRatios$F_ratio30_Upper[is.na(AvgRatios$F_ratio20_Upper)]
AvgRatios$F_ratio20_Lower[is.na(AvgRatios$F_ratio20_Lower)]=AvgRatios$F_ratio30_Lower[is.na(AvgRatios$F_ratio20_Lower)]
AvgRatios$OverfishingCategorySPR20=0
AvgRatios$OverfishingCategorySPR20[AvgRatios$OverfishingCategorySPR20==0 & AvgRatios$F_ratio20_Upper<1 & AvgRatios$F_ratio20_Lower<1]=1
AvgRatios$OverfishingCategorySPR20[AvgRatios$OverfishingCategorySPR20==0 & AvgRatios$F_ratio20_Upper>=1 & AvgRatios$F_ratio20_Lower<1]=2
AvgRatios$OverfishingCategorySPR20[AvgRatios$OverfishingCategorySPR20==0 & AvgRatios$F_ratio20_Upper>1 & AvgRatios$F_ratio20_Lower>=1]=3
AvgRatios$OverfishedCategorySPR20=0
AvgRatios$OverfishedCategorySPR20[AvgRatios$OverfishedCategorySPR20==0 & AvgRatios$B_ratio20_Lower>1 & AvgRatios$B_ratio20_Upper>1]=1
AvgRatios$OverfishedCategorySPR20[AvgRatios$OverfishedCategorySPR20==0 & AvgRatios$B_ratio20_Upper>=1 & AvgRatios$B_ratio20_Lower<1]=2
AvgRatios$OverfishedCategorySPR20[AvgRatios$OverfishedCategorySPR20==0 & AvgRatios$B_ratio20_Upper<=1 & AvgRatios$B_ratio20_Lower<1]=3
AvgRatios$F_ratio40_Lower = AvgRatios$F_ratio40_mean - AvgRatios$F_ratio40_sd
AvgRatios$F_ratio40_Upper = AvgRatios$F_ratio40_mean + AvgRatios$F_ratio40_sd
AvgRatios$B_ratio40_Lower = AvgRatios$B_ratio40_mean - AvgRatios$B_ratio40_sd
AvgRatios$B_ratio40_Upper = AvgRatios$B_ratio40_mean + AvgRatios$B_ratio40_sd
AvgRatios$F_ratio40_Lower[AvgRatios$F_ratio40_Lower<0]=0
AvgRatios$B_ratio40_Lower[AvgRatios$B_ratio40_Lower<0]=0
#If you can not compute the F or B multiplier at SPR 40 (convergence issue), borrow from those values at SPR 30, at least for the trend and "color" - i.e. overfished or experiencing overfishing 
#NOTE: this does not change the tabulated output or the actual values of these in the Kobe plots. 
AvgRatios$F_ratio40_mean[is.na(AvgRatios$F_ratio40_mean)]=AvgRatios$F_ratio30_mean[is.na(AvgRatios$F_ratio40_mean)]
AvgRatios$F_ratio40_sd[is.na(AvgRatios$F_ratio40_sd)]=AvgRatios$F_ratio30_sd[is.na(AvgRatios$F_ratio40_sd)]
AvgRatios$F_ratio40_Upper[is.na(AvgRatios$F_ratio40_Upper)]=AvgRatios$F_ratio30_Upper[is.na(AvgRatios$F_ratio40_Upper)]
AvgRatios$F_ratio40_Lower[is.na(AvgRatios$F_ratio40_Lower)]=AvgRatios$F_ratio30_Lower[is.na(AvgRatios$F_ratio40_Lower)]
AvgRatios$OverfishingCategorySPR40=0
AvgRatios$OverfishingCategorySPR40[AvgRatios$OverfishingCategorySPR40==0 & AvgRatios$F_ratio40_Upper<1 & AvgRatios$F_ratio40_Lower<1]=1
AvgRatios$OverfishingCategorySPR40[AvgRatios$OverfishingCategorySPR40==0 & AvgRatios$F_ratio40_Upper>=1 & AvgRatios$F_ratio40_Lower<1]=2
AvgRatios$OverfishingCategorySPR40[AvgRatios$OverfishingCategorySPR40==0 & AvgRatios$F_ratio40_Upper>1 & AvgRatios$F_ratio40_Lower>=1]=3
AvgRatios$OverfishedCategorySPR40=0
AvgRatios$OverfishedCategorySPR40[AvgRatios$OverfishedCategorySPR40==0 & AvgRatios$B_ratio40_Lower>1 & AvgRatios$B_ratio40_Upper>1]=1
AvgRatios$OverfishedCategorySPR40[AvgRatios$OverfishedCategorySPR40==0 & AvgRatios$B_ratio40_Upper>=1 & AvgRatios$B_ratio40_Lower<1]=2
AvgRatios$OverfishedCategorySPR40[AvgRatios$OverfishedCategorySPR40==0 & AvgRatios$B_ratio40_Upper<=1 & AvgRatios$B_ratio40_Lower<1]=3
#If you can not compute the F or B multiplier at SPR 20 or SPR 40 (convergence issue), borrow from those values at SPR 30, at least for the trend and "color" - i.e. overfished or experiencing overfishing 
#NOTE: this does not change the tabulated output or the actual values of these in the Kobe plots. 
AvgRatios$F_ratio20_mean[is.na(AvgRatios$F_ratio20_mean)]=AvgRatios$F_ratio30_mean[is.na(AvgRatios$F_ratio20_mean)]
AvgRatios$F_ratio20_sd[is.na(AvgRatios$F_ratio20_sd)]=AvgRatios$F_ratio30_sd[is.na(AvgRatios$F_ratio20_sd)]
AvgRatios$F_ratio20_Upper[is.na(AvgRatios$F_ratio20_Upper)]=AvgRatios$F_ratio30_Upper[is.na(AvgRatios$F_ratio20_Upper)]
AvgRatios$F_ratio20_Lower[is.na(AvgRatios$F_ratio20_Lower)]=AvgRatios$F_ratio30_Lower[is.na(AvgRatios$F_ratio20_Lower)]
AvgRatios$F_ratio40_mean[is.na(AvgRatios$F_ratio40_mean)]=AvgRatios$F_ratio30_mean[is.na(AvgRatios$F_ratio40_mean)]
AvgRatios$F_ratio40_sd[is.na(AvgRatios$F_ratio40_sd)]=AvgRatios$F_ratio30_sd[is.na(AvgRatios$F_ratio40_sd)]
AvgRatios$F_ratio40_Upper[is.na(AvgRatios$F_ratio40_Upper)]=AvgRatios$F_ratio30_Upper[is.na(AvgRatios$F_ratio40_Upper)]
AvgRatios$F_ratio40_Lower[is.na(AvgRatios$F_ratio40_Lower)]=AvgRatios$F_ratio30_Lower[is.na(AvgRatios$F_ratio40_Lower)]
AvgRatios$B_ratio20_mean[is.na(AvgRatios$B_ratio20_mean)]=AvgRatios$B_ratio30_mean[is.na(AvgRatios$B_ratio20_mean)]
AvgRatios$B_ratio20_sd[is.na(AvgRatios$B_ratio20_sd)]=AvgRatios$B_ratio30_sd[is.na(AvgRatios$B_ratio20_sd)]
AvgRatios$B_ratio20_Upper[is.na(AvgRatios$B_ratio20_Upper)]=AvgRatios$B_ratio30_Upper[is.na(AvgRatios$B_ratio20_Upper)]
AvgRatios$B_ratio20_Lower[is.na(AvgRatios$B_ratio20_Lower)]=AvgRatios$B_ratio30_Lower[is.na(AvgRatios$B_ratio20_Lower)]
AvgRatios$B_ratio40_mean[is.na(AvgRatios$B_ratio40_mean)]=AvgRatios$B_ratio30_mean[is.na(AvgRatios$B_ratio40_mean)]
AvgRatios$B_ratio40_sd[is.na(AvgRatios$B_ratio40_sd)]=AvgRatios$B_ratio30_sd[is.na(AvgRatios$B_ratio40_sd)]
AvgRatios$B_ratio40_Upper[is.na(AvgRatios$B_ratio40_Upper)]=AvgRatios$B_ratio30_Upper[is.na(AvgRatios$B_ratio40_Upper)]
AvgRatios$B_ratio40_Lower[is.na(AvgRatios$B_ratio40_Lower)]=AvgRatios$B_ratio30_Lower[is.na(AvgRatios$B_ratio40_Lower)]
write.table(AvgRatios,paste(PATH_output,"PopulationRatios_pieChart.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
#spr30 case
Status_wpp_overfishingSPR30 = as.data.frame.matrix(table(AvgRatios$WPP,AvgRatios$OverfishingCategorySPR30))
names(Status_wpp_overfishingSPR30) = c("NoOverfishingSPR30","PossibleOverfishingSPR30","OverfishingSPR30")
Status_wpp_overfishingSPR30$Total = Status_wpp_overfishingSPR30$OverfishingSPR30 + Status_wpp_overfishingSPR30$PossibleOverfishingSPR30 + Status_wpp_overfishingSPR30$NoOverfishingSPR30
Status_wpp_overfishingSPR30$Overfishing_percSPR30 = Status_wpp_overfishingSPR30$OverfishingSPR30/Status_wpp_overfishingSPR30$Total
Status_wpp_overfishingSPR30$PossibleOverfishing_percSPR30 = Status_wpp_overfishingSPR30$PossibleOverfishingSPR30/Status_wpp_overfishingSPR30$Total
Status_wpp_overfishingSPR30$NoOverfishing_percSPR30 = Status_wpp_overfishingSPR30$NoOverfishingSPR30/Status_wpp_overfishingSPR30$Total
Status_wpp_overfishingSPR30$WPP = row.names(Status_wpp_overfishingSPR30)
row.names(Status_wpp_overfishingSPR30)=NULL
Status_wpp_overfishedSPR30 = as.data.frame.matrix(table(AvgRatios$WPP,AvgRatios$OverfishedCategorySPR30))
names(Status_wpp_overfishedSPR30) = c("NoOverfishedSPR30","PossibleOverfishedSPR30","OverfishedSPR30")
Status_wpp_overfishedSPR30$Total = Status_wpp_overfishedSPR30$OverfishedSPR30 + Status_wpp_overfishedSPR30$PossibleOverfishedSPR30 + Status_wpp_overfishedSPR30$NoOverfishedSPR30
Status_wpp_overfishedSPR30$Overfished_percSPR30 = Status_wpp_overfishedSPR30$OverfishedSPR30/Status_wpp_overfishedSPR30$Total
Status_wpp_overfishedSPR30$PossibleOverfished_percSPR30 = Status_wpp_overfishedSPR30$PossibleOverfishedSPR30/Status_wpp_overfishedSPR30$Total
Status_wpp_overfishedSPR30$NoOverfished_percSPR30 = Status_wpp_overfishedSPR30$NoOverfishedSPR30/Status_wpp_overfishedSPR30$Total
Status_wpp_overfishedSPR30$WPP = row.names(Status_wpp_overfishedSPR30)
row.names(Status_wpp_overfishedSPR30)=NULL
Status_wpp_overfishingSPR30$Total=NULL
Status_wpp_overfishedSPR30$Total=NULL
Status_wppSPR30 = merge(Status_wpp_overfishingSPR30,Status_wpp_overfishedSPR30,all.x=TRUE,by=c("WPP"))
#spr20 case
Status_wpp_overfishingSPR20 = as.data.frame.matrix(table(AvgRatios$WPP,AvgRatios$OverfishingCategorySPR20))
names(Status_wpp_overfishingSPR20) = c("NoOverfishingSPR20","PossibleOverfishingSPR20","OverfishingSPR20")
Status_wpp_overfishingSPR20$Total = Status_wpp_overfishingSPR20$OverfishingSPR20 + Status_wpp_overfishingSPR20$PossibleOverfishingSPR20 + Status_wpp_overfishingSPR20$NoOverfishingSPR20
Status_wpp_overfishingSPR20$Overfishing_percSPR20 = Status_wpp_overfishingSPR20$OverfishingSPR20/Status_wpp_overfishingSPR20$Total
Status_wpp_overfishingSPR20$PossibleOverfishing_percSPR20 = Status_wpp_overfishingSPR20$PossibleOverfishingSPR20/Status_wpp_overfishingSPR20$Total
Status_wpp_overfishingSPR20$NoOverfishing_percSPR20 = Status_wpp_overfishingSPR20$NoOverfishingSPR20/Status_wpp_overfishingSPR20$Total
Status_wpp_overfishingSPR20$WPP = row.names(Status_wpp_overfishingSPR20)
row.names(Status_wpp_overfishingSPR20)=NULL
Status_wpp_overfishedSPR20 = as.data.frame.matrix(table(AvgRatios$WPP,AvgRatios$OverfishedCategorySPR20))
names(Status_wpp_overfishedSPR20) = c("NoOverfishedSPR20","PossibleOverfishedSPR20","OverfishedSPR20")
Status_wpp_overfishedSPR20$Total = Status_wpp_overfishedSPR20$OverfishedSPR20 + Status_wpp_overfishedSPR20$PossibleOverfishedSPR20 + Status_wpp_overfishedSPR20$NoOverfishedSPR20
Status_wpp_overfishedSPR20$Overfished_percSPR20 = Status_wpp_overfishedSPR20$OverfishedSPR20/Status_wpp_overfishedSPR20$Total
Status_wpp_overfishedSPR20$PossibleOverfished_percSPR20 = Status_wpp_overfishedSPR20$PossibleOverfishedSPR20/Status_wpp_overfishedSPR20$Total
Status_wpp_overfishedSPR20$NoOverfished_percSPR20 = Status_wpp_overfishedSPR20$NoOverfishedSPR20/Status_wpp_overfishedSPR20$Total
Status_wpp_overfishedSPR20$WPP = row.names(Status_wpp_overfishedSPR20)
row.names(Status_wpp_overfishedSPR20)=NULL
Status_wpp_overfishingSPR20$Total=NULL
Status_wpp_overfishedSPR20$Total=NULL
Status_wppSPR20 = merge(Status_wpp_overfishingSPR20,Status_wpp_overfishedSPR20,all.x=TRUE,by=c("WPP"))
#spr40 case
Status_wpp_overfishingSPR40 = as.data.frame.matrix(table(AvgRatios$WPP,AvgRatios$OverfishingCategorySPR40))
names(Status_wpp_overfishingSPR40) = c("NoOverfishingSPR40","PossibleOverfishingSPR40","OverfishingSPR40")
Status_wpp_overfishingSPR40$Total = Status_wpp_overfishingSPR40$OverfishingSPR40 + Status_wpp_overfishingSPR40$PossibleOverfishingSPR40 + Status_wpp_overfishingSPR40$NoOverfishingSPR40
Status_wpp_overfishingSPR40$Overfishing_percSPR40 = Status_wpp_overfishingSPR40$OverfishingSPR40/Status_wpp_overfishingSPR40$Total
Status_wpp_overfishingSPR40$PossibleOverfishing_percSPR40 = Status_wpp_overfishingSPR40$PossibleOverfishingSPR40/Status_wpp_overfishingSPR40$Total
Status_wpp_overfishingSPR40$NoOverfishing_percSPR40 = Status_wpp_overfishingSPR40$NoOverfishingSPR40/Status_wpp_overfishingSPR40$Total
Status_wpp_overfishingSPR40$WPP = row.names(Status_wpp_overfishingSPR40)
row.names(Status_wpp_overfishingSPR40)=NULL
Status_wpp_overfishedSPR40 = as.data.frame.matrix(table(AvgRatios$WPP,AvgRatios$OverfishedCategorySPR40))
names(Status_wpp_overfishedSPR40) = c("NoOverfishedSPR40","PossibleOverfishedSPR40","OverfishedSPR40")
Status_wpp_overfishedSPR40$Total = Status_wpp_overfishedSPR40$OverfishedSPR40 + Status_wpp_overfishedSPR40$PossibleOverfishedSPR40 + Status_wpp_overfishedSPR40$NoOverfishedSPR40
Status_wpp_overfishedSPR40$Overfished_percSPR40 = Status_wpp_overfishedSPR40$OverfishedSPR40/Status_wpp_overfishedSPR40$Total
Status_wpp_overfishedSPR40$PossibleOverfished_percSPR40 = Status_wpp_overfishedSPR40$PossibleOverfishedSPR40/Status_wpp_overfishedSPR40$Total
Status_wpp_overfishedSPR40$NoOverfished_percSPR40 = Status_wpp_overfishedSPR40$NoOverfishedSPR40/Status_wpp_overfishedSPR40$Total
Status_wpp_overfishedSPR40$WPP = row.names(Status_wpp_overfishedSPR40)
row.names(Status_wpp_overfishedSPR40)=NULL
Status_wpp_overfishingSPR40$Total=NULL
Status_wpp_overfishedSPR40$Total=NULL
Status_wppSPR40 = merge(Status_wpp_overfishingSPR40,Status_wpp_overfishedSPR40,all.x=TRUE,by=c("WPP"))
Status_wpp = merge(Status_wppSPR30,Status_wppSPR20,by=c("WPP"))
Status_wpp = merge(Status_wpp,Status_wppSPR40,by=c("WPP"))
LocationsWPP$LocationsWPP_ShapeFile=NULL
Status_wpp = merge(Status_wpp,LocationsWPP,all.x=TRUE,by=c("WPP"))
Status_wpp = subset(Status_wpp,!is.na(Status_wpp$Longitude))
Status_wpp$Total = Status_wpp$NoOverfishingSPR30 + Status_wpp$PossibleOverfishingSPR30 + Status_wpp$OverfishingSPR30

#Pie Chart Sample Size Text Locations for Each Pie
TextSampleSizeLocations = data.frame(x=c(100,100.46866,117.3,107.47709,113.94177,118.77517,126.68987,127.05237,122.58147,135.75250),
  y=c(7,-3.0658615,-8.3130337,4.2277078,-3.3282201,-0.8620492,-5.4795607,1.6467065,4.6474816,-7.9982033))

#Plot Overfished Status at SPR 30%
png(paste(PATH_highLevelSummaryPlots,"OverfishedMapSPR30_uncertainty.png",sep="/"),units="px",width=5200,height=3200,res=600)
par(mai = c(0.5, 0.5, 0.35, 0.2), omi = c(0.25, 0.5, 0, 0),mgp = c(2.5, 0.5, 0))
prettymap(plot(WPPshapeFile,col=c("green",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),xlab="Longitude",ylab="Latitude"), oma=c(4,4,1,0),drawscale=TRUE,drawbox=FALSE, scale.plotunit="latlon", drawarrow=TRUE)
axis(side=1,outer=TRUE,line=-2.5)
axis(side=2,outer=TRUE,line=-2)
mtext("Longitude",side=1,line=2,cex=1.2)
mtext("Latitude",side=2,line=2.25,cex=1.2)
title("Population Status: Overfished at SPR 30%")
mtext("Indonesia Deepwater Snapper-Grouper Fishery",font=3,line=-0.5)
Overfished = subset(Status_wpp,select=c(Longitude,Latitude,Overfished_percSPR30))
names(Overfished)[names(Overfished)=="Overfished_percSPR30"]="Value"
Overfished$Variable = "Overfished"
PossibleOverfished = subset(Status_wpp,select=c(Longitude,Latitude,PossibleOverfished_percSPR30))
names(PossibleOverfished)[names(PossibleOverfished)=="PossibleOverfished_percSPR30"]="Value"
PossibleOverfished$Variable = "PossibleOverfished"
NoOverfished = subset(Status_wpp,select=c(Longitude,Latitude,NoOverfished_percSPR30))
names(NoOverfished)[names(NoOverfished)=="NoOverfished_percSPR30"]="Value"
NoOverfished$Variable = "NoOverfished"
formattedForPieCharts = rbind(Overfished,PossibleOverfished,NoOverfished)
xyz = make.xyz(formattedForPieCharts$Longitude,formattedForPieCharts$Latitude,formattedForPieCharts$Value,formattedForPieCharts$Variable)
for(i in 1:(dim(Status_wpp)[1]))
{
  add.pie(xyz$z[i,],xyz$x[i],xyz$y[i],radius=1.2,col=c("green3","yellow2","red3"),labels=NA)  #,labels=c(NA,paste("n=",Status_WPP$N[i],sep="")),label.dist=1.3)
}
legend.pie(100,-10,labels=c("Not Overfished","Possibly Overfished","Overfished"),radius=1,bty="n",col=c("green3","yellow2","red3"),cex=0.9,label.dist=1.6)
for(i in 1:(dim(Status_wpp)[1]))
{
  text(TextSampleSizeLocations$x[i],TextSampleSizeLocations$y[i],paste("n=",Status_wpp$Total[i],sep=""))
}
dev.off()

#Plot Overfishing Status at SPR 30%
png(paste(PATH_highLevelSummaryPlots,"OverfishingMapSPR30_uncertainty.png",sep="/"),units="px",width=5200,height=3200,res=600)
par(mai = c(0.5, 0.5, 0.35, 0.2), omi = c(0.25, 0.5, 0, 0),mgp = c(2.5, 0.5, 0))
prettymap(plot(WPPshapeFile,col=c("green",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),xlab="Longitude",ylab="Latitude"), oma=c(4,4,1,0),drawscale=TRUE,drawbox=FALSE, scale.plotunit="latlon", drawarrow=TRUE)
axis(side=1,outer=TRUE,line=-2.5)
axis(side=2,outer=TRUE,line=-2)
mtext("Longitude",side=1,line=2,cex=1.2)
mtext("Latitude",side=2,line=2.25,cex=1.2)
title("Population Status: Overfishing at SPR 30%")
mtext("Indonesia Deepwater Snapper-Grouper Fishery",font=3,line=-0.5)
Overfishing = subset(Status_wpp,select=c(Longitude,Latitude,Overfishing_percSPR30))
names(Overfishing)[names(Overfishing)=="Overfishing_percSPR30"]="Value"
Overfishing$Variable = "Overfishing"
PossibleOverfishing = subset(Status_wpp,select=c(Longitude,Latitude,PossibleOverfishing_percSPR30))
names(PossibleOverfishing)[names(PossibleOverfishing)=="PossibleOverfishing_percSPR30"]="Value"
PossibleOverfishing$Variable = "PossibleOverfishing"
NoOverfishing = subset(Status_wpp,select=c(Longitude,Latitude,NoOverfishing_percSPR30))
names(NoOverfishing)[names(NoOverfishing)=="NoOverfishing_percSPR30"]="Value"
NoOverfishing$Variable = "NoOverfishing"
formattedForPieCharts = rbind(Overfishing,PossibleOverfishing,NoOverfishing)
xyz = make.xyz(formattedForPieCharts$Longitude,formattedForPieCharts$Latitude,formattedForPieCharts$Value,formattedForPieCharts$Variable)
for(i in 1:(dim(Status_wpp)[1]))
{
  add.pie(xyz$z[i,],xyz$x[i],xyz$y[i],radius=1.2,col=c("green3","yellow2","red3"),labels=NA)  #,labels=c(NA,paste("n=",Status_WPP$N[i],sep="")),label.dist=1.3)
}
legend.pie(100,-10,labels=c("Not Overfishing","Possibly Overfishing","Overfishing"),radius=1,bty="n",col=c("green3","yellow2","red3"),cex=0.9,label.dist=1.6)
for(i in 1:(dim(Status_wpp)[1]))
{
  text(TextSampleSizeLocations$x[i],TextSampleSizeLocations$y[i],paste("n=",Status_wpp$Total[i],sep=""))
}
dev.off()

#Plot Overfished Status at SPR 20%
png(paste(PATH_highLevelSummaryPlots,"OverfishedMapSPR20_uncertainty.png",sep="/"),units="px",width=5200,height=3200,res=600)
par(mai = c(0.5, 0.5, 0.35, 0.2), omi = c(0.25, 0.5, 0, 0),mgp = c(2.5, 0.5, 0))
prettymap(plot(WPPshapeFile,col=c("green",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),xlab="Longitude",ylab="Latitude"), oma=c(4,4,1,0),drawscale=TRUE,drawbox=FALSE, scale.plotunit="latlon", drawarrow=TRUE)
axis(side=1,outer=TRUE,line=-2.5)
axis(side=2,outer=TRUE,line=-2)
mtext("Longitude",side=1,line=2,cex=1.2)
mtext("Latitude",side=2,line=2.25,cex=1.2)
title("Population Status: Overfished at SPR 20%")
mtext("Indonesia Deepwater Snapper-Grouper Fishery",font=3,line=-0.5)
Overfished = subset(Status_wpp,select=c(Longitude,Latitude,Overfished_percSPR20))
names(Overfished)[names(Overfished)=="Overfished_percSPR20"]="Value"
Overfished$Variable = "Overfished"
PossibleOverfished = subset(Status_wpp,select=c(Longitude,Latitude,PossibleOverfished_percSPR20))
names(PossibleOverfished)[names(PossibleOverfished)=="PossibleOverfished_percSPR20"]="Value"
PossibleOverfished$Variable = "PossibleOverfished"
NoOverfished = subset(Status_wpp,select=c(Longitude,Latitude,NoOverfished_percSPR20))
names(NoOverfished)[names(NoOverfished)=="NoOverfished_percSPR20"]="Value"
NoOverfished$Variable = "NoOverfished"
formattedForPieCharts = rbind(Overfished,PossibleOverfished,NoOverfished)
xyz = make.xyz(formattedForPieCharts$Longitude,formattedForPieCharts$Latitude,formattedForPieCharts$Value,formattedForPieCharts$Variable)
for(i in 1:(dim(Status_wpp)[1]))
{
  add.pie(xyz$z[i,],xyz$x[i],xyz$y[i],radius=1.2,col=c("green3","yellow2","red3"),labels=NA)  #,labels=c(NA,paste("n=",Status_WPP$N[i],sep="")),label.dist=1.3)
}
legend.pie(100,-10,labels=c("Not Overfished","Possibly Overfished","Overfished"),radius=1,bty="n",col=c("green3","yellow2","red3"),cex=0.9,label.dist=1.6)
for(i in 1:(dim(Status_wpp)[1]))
{
  text(TextSampleSizeLocations$x[i],TextSampleSizeLocations$y[i],paste("n=",Status_wpp$Total[i],sep=""))
}
dev.off()

#Plot Overfishing Status at SPR 20%
png(paste(PATH_highLevelSummaryPlots,"OverfishingMapSPR20_uncertainty.png",sep="/"),units="px",width=5200,height=3200,res=600)
par(mai = c(0.5, 0.5, 0.35, 0.2), omi = c(0.25, 0.5, 0, 0),mgp = c(2.5, 0.5, 0))
prettymap(plot(WPPshapeFile,col=c("green",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),xlab="Longitude",ylab="Latitude"), oma=c(4,4,1,0),drawscale=TRUE,drawbox=FALSE, scale.plotunit="latlon", drawarrow=TRUE)
axis(side=1,outer=TRUE,line=-2.5)
axis(side=2,outer=TRUE,line=-2)
mtext("Longitude",side=1,line=2,cex=1.2)
mtext("Latitude",side=2,line=2.25,cex=1.2)
title("Population Status: Overfishing at SPR 20%")
mtext("Indonesia Deepwater Snapper-Grouper Fishery",font=3,line=-0.5)
Overfishing = subset(Status_wpp,select=c(Longitude,Latitude,Overfishing_percSPR20))
names(Overfishing)[names(Overfishing)=="Overfishing_percSPR20"]="Value"
Overfishing$Variable = "Overfishing"
PossibleOverfishing = subset(Status_wpp,select=c(Longitude,Latitude,PossibleOverfishing_percSPR20))
names(PossibleOverfishing)[names(PossibleOverfishing)=="PossibleOverfishing_percSPR20"]="Value"
PossibleOverfishing$Variable = "PossibleOverfishing"
NoOverfishing = subset(Status_wpp,select=c(Longitude,Latitude,NoOverfishing_percSPR20))
names(NoOverfishing)[names(NoOverfishing)=="NoOverfishing_percSPR20"]="Value"
NoOverfishing$Variable = "NoOverfishing"
formattedForPieCharts = rbind(Overfishing,PossibleOverfishing,NoOverfishing)
xyz = make.xyz(formattedForPieCharts$Longitude,formattedForPieCharts$Latitude,formattedForPieCharts$Value,formattedForPieCharts$Variable)
for(i in 1:(dim(Status_wpp)[1]))
{
  add.pie(xyz$z[i,],xyz$x[i],xyz$y[i],radius=1.2,col=c("green3","yellow2","red3"),labels=NA)  #,labels=c(NA,paste("n=",Status_WPP$N[i],sep="")),label.dist=1.3)
}
legend.pie(100,-10,labels=c("Not Overfishing","Possibly Overfishing","Overfishing"),radius=1,bty="n",col=c("green3","yellow2","red3"),cex=0.9,label.dist=1.6)
for(i in 1:(dim(Status_wpp)[1]))
{
  text(TextSampleSizeLocations$x[i],TextSampleSizeLocations$y[i],paste("n=",Status_wpp$Total[i],sep=""))
}
dev.off()

#Plot Overfished Status at SPR 40%
png(paste(PATH_highLevelSummaryPlots,"OverfishedMapSPR40_uncertainty.png",sep="/"),units="px",width=5400,height=3400,res=600)
par(mai = c(0.5, 0.5, 0.35, 0.2), omi = c(0.25, 0.5, 0, 0),mgp = c(2.5, 0.5, 0))
prettymap(plot(WPPshapeFile,col=c("green",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),xlab="Longitude",ylab="Latitude"), oma=c(4,4,1,0),drawscale=TRUE,drawbox=FALSE, scale.plotunit="latlon", drawarrow=TRUE)
axis(side=1,outer=TRUE,line=-2.5)
axis(side=2,outer=TRUE,line=-2)
mtext("Longitude",side=1,line=2,cex=1.2)
mtext("Latitude",side=2,line=2.25,cex=1.2)
title("Population Status: Overfished at SPR 40%")
mtext("Indonesia Deepwater Snapper-Grouper Fishery",font=3,line=-0.5)
Overfished = subset(Status_wpp,select=c(Longitude,Latitude,Overfished_percSPR40))
names(Overfished)[names(Overfished)=="Overfished_percSPR40"]="Value"
Overfished$Variable = "Overfished"
PossibleOverfished = subset(Status_wpp,select=c(Longitude,Latitude,PossibleOverfished_percSPR40))
names(PossibleOverfished)[names(PossibleOverfished)=="PossibleOverfished_percSPR40"]="Value"
PossibleOverfished$Variable = "PossibleOverfished"
NoOverfished = subset(Status_wpp,select=c(Longitude,Latitude,NoOverfished_percSPR40))
names(NoOverfished)[names(NoOverfished)=="NoOverfished_percSPR40"]="Value"
NoOverfished$Variable = "NoOverfished"
formattedForPieCharts = rbind(Overfished,PossibleOverfished,NoOverfished)
xyz = make.xyz(formattedForPieCharts$Longitude,formattedForPieCharts$Latitude,formattedForPieCharts$Value,formattedForPieCharts$Variable)
for(i in 1:(dim(Status_wpp)[1]))
{
  add.pie(xyz$z[i,],xyz$x[i],xyz$y[i],radius=1.2,col=c("green3","yellow2","red3"),labels=NA)  #,labels=c(NA,paste("n=",Status_WPP$N[i],sep="")),label.dist=1.3)
}
legend.pie(100,-10,labels=c("Not Overfished","Possibly Overfished","Overfished"),radius=1,bty="n",col=c("green3","yellow2","red3"),cex=0.9,label.dist=1.6)
for(i in 1:(dim(Status_wpp)[1]))
{
  text(TextSampleSizeLocations$x[i],TextSampleSizeLocations$y[i],paste("n=",Status_wpp$Total[i],sep=""))
}
dev.off()

#Plot Overfishing Status at SPR 40%
png(paste(PATH_highLevelSummaryPlots,"OverfishingMapSPR40_uncertainty.png",sep="/"),units="px",width=5400,height=3400,res=600)
par(mai = c(0.5, 0.5, 0.35, 0.2), omi = c(0.25, 0.5, 0, 0),mgp = c(2.5, 0.5, 0))
prettymap(plot(WPPshapeFile,col=c("green",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),xlab="Longitude",ylab="Latitude"), oma=c(4,4,1,0),drawscale=TRUE,drawbox=FALSE, scale.plotunit="latlon", drawarrow=TRUE)
axis(side=1,outer=TRUE,line=-2.5)
axis(side=2,outer=TRUE,line=-2)
mtext("Longitude",side=1,line=2,cex=1.2)
mtext("Latitude",side=2,line=2.25,cex=1.2)
title("Population Status: Overfishing at SPR 40%")
mtext("Indonesia Deepwater Snapper-Grouper Fishery",font=3,line=-0.5)
Overfishing = subset(Status_wpp,select=c(Longitude,Latitude,Overfishing_percSPR40))
names(Overfishing)[names(Overfishing)=="Overfishing_percSPR40"]="Value"
Overfishing$Variable = "Overfishing"
PossibleOverfishing = subset(Status_wpp,select=c(Longitude,Latitude,PossibleOverfishing_percSPR40))
names(PossibleOverfishing)[names(PossibleOverfishing)=="PossibleOverfishing_percSPR40"]="Value"
PossibleOverfishing$Variable = "PossibleOverfishing"
NoOverfishing = subset(Status_wpp,select=c(Longitude,Latitude,NoOverfishing_percSPR40))
names(NoOverfishing)[names(NoOverfishing)=="NoOverfishing_percSPR40"]="Value"
NoOverfishing$Variable = "NoOverfishing"
formattedForPieCharts = rbind(Overfishing,PossibleOverfishing,NoOverfishing)
xyz = make.xyz(formattedForPieCharts$Longitude,formattedForPieCharts$Latitude,formattedForPieCharts$Value,formattedForPieCharts$Variable)
for(i in 1:(dim(Status_wpp)[1]))
{
  add.pie(xyz$z[i,],xyz$x[i],xyz$y[i],radius=1.2,col=c("green3","yellow2","red3"),labels=NA)  #,labels=c(NA,paste("n=",Status_WPP$N[i],sep="")),label.dist=1.3)
}
legend.pie(100,-10,labels=c("Not Overfishing","Possibly Overfishing","Overfishing"),radius=1,bty="n",col=c("green3","yellow2","red3"),cex=0.9,label.dist=1.6)
for(i in 1:(dim(Status_wpp)[1]))
{
  text(TextSampleSizeLocations$x[i],TextSampleSizeLocations$y[i],paste("n=",Status_wpp$Total[i],sep=""))
}
dev.off()


######################### Create Generic Indonesia Map Showing the WPPs ##############################
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
png(paste(PATH_highLevelSummaryPlots,"WPP_map.png",sep="/"),units="px",width=5200,height=3200,res=600)
par(mai = c(0.5, 0.5, 0.35, 0.2), omi = c(0.25, 0.5, 0, 0),mgp = c(2.5, 0.5, 0))
prettymap(plot(WPPshapeFile,col=c("green",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),xlab="Longitude",ylab="Latitude"), oma=c(4,4,1,0),drawscale=TRUE,drawbox=FALSE, scale.plotunit="latlon", drawarrow=TRUE)
axis(side=1,outer=TRUE,line=-2.5)
axis(side=2,outer=TRUE,line=-2)
mtext("Longitude",side=1,line=2,cex=1.2)
mtext("Latitude",side=2,line=2.25,cex=1.2)
title("Indonesian Fisheries Management Areas")
mtext("Wilayah Pengelolaan Perikanan (WPP)",font=3,line=-0.5)
WPPs = c("571","572","573","711","712","713","714","715","716","717","718")
MidpointsWPP = data.frame(x=c(97.5,98.59572,114.72720,109.47087,111.40423,118.53350,126.44820,125.84402,122.70231,135.99417,135.45041),
  y=c(6,-1.4917099,-8.3130337,5.5395008,-3.3282201,-2.6460877,-3.8004656,0.1873852,2.6535561,1.4991783,-6.1616931))
text(MidpointsWPP,WPPs,font=4)
dev.off()


################# Extract Information on Rebuilding Times #######################################
#and those that will NOT Rebuild

SummaryProjectionScenarios = read.table(paste(PATH_output,"SummaryProjectionScenarios.csv",sep="/"),header=TRUE,sep=",",as.is=TRUE)
ListOfPopulations_WithProjections = readRDS(file=paste(PATH_output,"ListOfPopulations_WithProjections.RData",sep="/"),refhook=NULL)

CurrentStatusValues = read.table(paste(PATH_output,"CurrentValuesBenchmarksAndRatios.csv",sep="/"),header=TRUE,sep=",",as.is=TRUE)
CurrentStatusValues$GenusSpecies = paste(CurrentStatusValues$Genus,CurrentStatusValues$Species,sep="_")

#TEMP CODE LINE: do NOT include SolveBaranov or TropFishR in both F and B ratio calculations; do NOT include LIME in B ratio calculations!!
CurrentStatusValues = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="SolveBaranov",]   #TEMP CODE LINE: do NOT include SolveBaranov or TropFishR
CurrentStatusValues = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="TropFishR",]      #TEMP CODE LINE: do NOT include SolveBaranov or TropFishR
CurrentStatusValues_noLIME = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="LIME",]

##Identify and filter away the Runs that did not converge
#CurrentStatusValues$F_ratio20_FLAG=0
#CurrentStatusValues$F_ratio20_FLAG[CurrentStatusValues$F_ratio20 > 5 | CurrentStatusValues$F_ratio20 < 0.01]=1
#CurrentStatusValues$F_ratio30_FLAG=0
#CurrentStatusValues$F_ratio30_FLAG[CurrentStatusValues$F_ratio30 > 5 | CurrentStatusValues$F_ratio30 < 0.01]=1
#CurrentStatusValues$F_ratio40_FLAG=0
#CurrentStatusValues$F_ratio40_FLAG[CurrentStatusValues$F_ratio40 > 5 | CurrentStatusValues$F_ratio40 < 0.01]=1
#CurrentStatusValues$B_ratio20_FLAG=0
#CurrentStatusValues$B_ratio20_FLAG[CurrentStatusValues$B_ratio20 > 5 | CurrentStatusValues$B_ratio20 < 0.01]=1
#CurrentStatusValues$B_ratio30_FLAG=0
#CurrentStatusValues$B_ratio30_FLAG[CurrentStatusValues$B_ratio30 > 5 | CurrentStatusValues$B_ratio30 < 0.01]=1
#CurrentStatusValues$B_ratio40_FLAG=0
#CurrentStatusValues$B_ratio40_FLAG[CurrentStatusValues$B_ratio40 > 5 | CurrentStatusValues$B_ratio40 < 0.01]=1
#CurrentStatusValues_toMerge = subset(CurrentStatusValues,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,populationReconstructionMethod,
#  F_ratio20_FLAG,F_ratio30_FLAG,F_ratio40_FLAG,B_ratio20_FLAG,B_ratio30_FLAG,B_ratio40_FLAG))
#CurrentStatusValues_toMerge = unique(CurrentStatusValues_toMerge)
#SummaryProjectionScenarios = merge(SummaryProjectionScenarios,CurrentStatusValues_toMerge,all.x=TRUE,
#  by=c("Genus","Species","WPP","CatchMSY_Scenario","CatchMSY_Scenario_Bound","LifeHistory_Scenario","populationReconstructionMethod"))
#
#create genus/species to fish family mapping and merge it in here.
GenusSpeciesCommonName = subset(Lengths_cpue,select=c(Genus,Species,FamilyCommonName))
GenusSpeciesCommonName = unique(GenusSpeciesCommonName)
SummaryProjectionScenarios = merge(SummaryProjectionScenarios,GenusSpeciesCommonName,all.x=TRUE,by=c("Genus","Species"))
SummaryProjectionScenarios_withRebuild = subset(SummaryProjectionScenarios,SummaryProjectionScenarios$FishingMortaltiyProjectionScenario=="F_rebuild_Tmin_oneGeneration")
SummaryProjectionScenarios_withRebuild$TimeToRebuild_NoFishing=0
SummaryProjectionScenarios_withRebuild$OneGenerationTime=0
SummaryProjectionScenarios_withRebuild$PopulationCanRebuild=TRUE

#Identify any runs that did not converge (as per above ratio FLAG parameters) and remove them here:
#SummaryProjectionScenarios_withRebuild$FLAG_SUM = SummaryProjectionScenarios_withRebuild$F_ratio20_FLAG +
#  SummaryProjectionScenarios_withRebuild$F_ratio30_FLAG + SummaryProjectionScenarios_withRebuild$F_ratio40_FLAG +
#  SummaryProjectionScenarios_withRebuild$B_ratio20_FLAG + SummaryProjectionScenarios_withRebuild$B_ratio30_FLAG + SummaryProjectionScenarios_withRebuild$B_ratio40_FLAG

for(i in 1:dim(SummaryProjectionScenarios_withRebuild)[1])
{
  if(ListOfPopulations_WithProjections[[SummaryProjectionScenarios_withRebuild$ListIndexNumber[i]]]$ProjectionType=="RegularF_willNotRebuildAtF=0")
  {
    SummaryProjectionScenarios_withRebuild$PopulationCanRebuild[i]=FALSE
  }
  if(ListOfPopulations_WithProjections[[SummaryProjectionScenarios_withRebuild$ListIndexNumber[i]]]$ProjectionType!="RegularF_willNotRebuildAtF=0")
  {
    SummaryProjectionScenarios_withRebuild$TimeToRebuild_NoFishing[i] = ListOfPopulations_WithProjections[[SummaryProjectionScenarios_withRebuild$ListIndexNumber[i]]]$timeToRebuild_noFishing
    SummaryProjectionScenarios_withRebuild$OneGenerationTime[i] = ListOfPopulations_WithProjections[[SummaryProjectionScenarios_withRebuild$ListIndexNumber[i]]]$OneGenerationTime
  }
  if(i%%5000==0)
  {
    print(paste("Working on item ",i," of ",dim(SummaryProjectionScenarios_withRebuild)[1],"...",sep=""))
    flush.console()
  }
}
SummaryProjectionScenarios_withRebuild$TotalRebuildTime = SummaryProjectionScenarios_withRebuild$TimeToRebuild_NoFishing + SummaryProjectionScenarios_withRebuild$OneGenerationTime
#Remove runs that took more than 30 years to rebuild, and obviously didn't converge properly.
SummaryProjectionScenarios_withRebuild = subset(SummaryProjectionScenarios_withRebuild,SummaryProjectionScenarios_withRebuild$TotalRebuildTime<30)
SummaryProjectionScenarios_noRebuild = subset(SummaryProjectionScenarios_withRebuild,SummaryProjectionScenarios_withRebuild$PopulationCanRebuild==FALSE)
NumRunsNoRebuild = data.frame(table(SummaryProjectionScenarios_noRebuild$Genus,SummaryProjectionScenarios_noRebuild$Species,SummaryProjectionScenarios_noRebuild$WPP))
names(NumRunsNoRebuild) = c("Genus","Species","WPP","Freq_noRebuild")
NumRunsNoRebuild = subset(NumRunsNoRebuild,NumRunsNoRebuild$Freq_noRebuild!=0)
TotalTriedToRebuild = data.frame(table(SummaryProjectionScenarios_withRebuild$Genus,SummaryProjectionScenarios_withRebuild$Species,SummaryProjectionScenarios_withRebuild$WPP))
names(TotalTriedToRebuild) = c("Genus","Species","WPP","TotalRuns")
TotalTriedToRebuild = subset(TotalTriedToRebuild,TotalTriedToRebuild$TotalRuns>0)
TotalTriedToRebuild = merge(TotalTriedToRebuild,NumRunsNoRebuild,all.x=TRUE,by=c("Genus","Species","WPP"))
TotalTriedToRebuild$FractionRunsNoRebuild = TotalTriedToRebuild$Freq_noRebuild/TotalTriedToRebuild$TotalRuns
TotalTriedToRebuild = subset(TotalTriedToRebuild,!is.na(TotalTriedToRebuild$FractionRunsNoRebuild))
TotalTriedToRebuild_WPP = data.frame(table(SummaryProjectionScenarios_noRebuild$WPP))
names(TotalTriedToRebuild_WPP) = c("WPP","TotalRuns_noRebuild")
TotalRuns_WPP = data.frame(table(SummaryProjectionScenarios_withRebuild$WPP))
names(TotalRuns_WPP) = c("WPP","TotalRuns")
RunsNotRebuilding_WPP = merge(TotalRuns_WPP,TotalTriedToRebuild_WPP,all.x=TRUE,by=c("WPP"))
RunsNotRebuilding_WPP[is.na(RunsNotRebuilding_WPP)]=0
RunsNotRebuilding_WPP$PercentNoRebuild = RunsNotRebuilding_WPP$TotalRuns_noRebuild/RunsNotRebuilding_WPP$TotalRuns
AvgRebuildingTimes = aggregate.data.frame(list(SummaryProjectionScenarios_withRebuild$TimeToRebuild_NoFishing,SummaryProjectionScenarios_withRebuild$OneGenerationTime,SummaryProjectionScenarios_withRebuild$TotalRebuildTime),by=list(SummaryProjectionScenarios_withRebuild$Genus,SummaryProjectionScenarios_withRebuild$Species,SummaryProjectionScenarios_withRebuild$WPP),FUN=mean)
names(AvgRebuildingTimes) = c("Genus","Species","WPP","TimeToRebuild_NoFishing","OneGenerationTime","TotalRebuilding")
AvgRebuildingTimes$GenusSpecies = paste(AvgRebuildingTimes$Genus,AvgRebuildingTimes$Species)
AvgRebuildingTimes_totalRebuild = subset(AvgRebuildingTimes,select=c(GenusSpecies,WPP,TotalRebuilding))
AvgRebuildingTimes_totalRebuild = reshape(AvgRebuildingTimes_totalRebuild,v.names="TotalRebuilding",idvar="GenusSpecies",timevar="WPP",direction="wide")
GenusSpecies = data.frame(t(data.frame((strsplit(AvgRebuildingTimes_totalRebuild$GenusSpecies," ")))))
names(GenusSpecies) = c("Genus","Species")
AvgRebuildingTimes_totalRebuild = cbind(AvgRebuildingTimes_totalRebuild,GenusSpecies)
AvgRebuildingTimes_totalRebuild$GenusSpecies=NULL
AvgRebuildingTimes_totalRebuild = AvgRebuildingTimes_totalRebuild[order(AvgRebuildingTimes_totalRebuild$Genus,AvgRebuildingTimes_totalRebuild$Species),]
write.table(AvgRebuildingTimes_totalRebuild,paste(PATH_output,"AvgRebuildingTimes.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
StandardDevRebuildingTimes = aggregate.data.frame(list(SummaryProjectionScenarios_withRebuild$TimeToRebuild_NoFishing,SummaryProjectionScenarios_withRebuild$OneGenerationTime,SummaryProjectionScenarios_withRebuild$TotalRebuildTime),by=list(SummaryProjectionScenarios_withRebuild$Genus,SummaryProjectionScenarios_withRebuild$Species,SummaryProjectionScenarios_withRebuild$WPP),FUN=sd)
names(StandardDevRebuildingTimes) = c("Genus","Species","WPP","TimeToRebuild_NoFishing","OneGenerationTime","TotalRebuilding")
StandardDevRebuildingTimes$GenusSpecies = paste(StandardDevRebuildingTimes$Genus,StandardDevRebuildingTimes$Species)
StandardDevRebuildingTimes_totalRebuild = subset(StandardDevRebuildingTimes,select=c(GenusSpecies,WPP,TotalRebuilding))
StandardDevRebuildingTimes_totalRebuild = reshape(StandardDevRebuildingTimes_totalRebuild,v.names="TotalRebuilding",idvar="GenusSpecies",timevar="WPP",direction="wide")
GenusSpecies = data.frame(t(data.frame((strsplit(StandardDevRebuildingTimes_totalRebuild$GenusSpecies," ")))))
names(GenusSpecies) = c("Genus","Species")
StandardDevRebuildingTimes_totalRebuild = cbind(StandardDevRebuildingTimes_totalRebuild,GenusSpecies)
StandardDevRebuildingTimes_totalRebuild$GenusSpecies=NULL
StandardDevRebuildingTimes_totalRebuild = StandardDevRebuildingTimes_totalRebuild[order(StandardDevRebuildingTimes_totalRebuild$Genus,StandardDevRebuildingTimes_totalRebuild$Species),]
write.table(StandardDevRebuildingTimes_totalRebuild,paste(PATH_output,"StandardDevRebuildingTimes_totalRebuild.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")

#Rebuilding Table and Plots by Fish Family and WPP
AvgRebuildingTimes_family = aggregate.data.frame(list(SummaryProjectionScenarios_withRebuild$TimeToRebuild_NoFishing,SummaryProjectionScenarios_withRebuild$OneGenerationTime,SummaryProjectionScenarios_withRebuild$TotalRebuildTime),by=list(SummaryProjectionScenarios_withRebuild$FamilyCommonName,SummaryProjectionScenarios_withRebuild$WPP),FUN=mean)
names(AvgRebuildingTimes_family) = c("Family","WPP","TimeToRebuild_NoFishing","OneGenerationTime","TotalRebuilding")
AvgRebuildingTimes_family_total = subset(AvgRebuildingTimes_family,select=c(Family,WPP,TotalRebuilding))
AvgRebuildingTimes_family_total = reshape(AvgRebuildingTimes_family_total,v.names="TotalRebuilding",idvar="Family",timevar="WPP",direction="wide")

#Histogram of Rebuilding Times
png(paste(PATH_highLevelSummaryPlots,"RebuildingTimes_Histogram.png",sep="/"),units="px",width=3200,height=3200,res=600)
hist(SummaryProjectionScenarios_withRebuild$TotalRebuildTime,freq=FALSE,xlab="Rebuilding Times (years)",main="Rebuilding Times",breaks=seq(0,ceiling(max(SummaryProjectionScenarios_withRebuild$TotalRebuildTime)+1),by=1))
dev.off()

FamilyNames = unique(SummaryProjectionScenarios_withRebuild$FamilyCommonName)
for(i in 1:length(FamilyNames))
{
  SummaryProjectionScenarios_withRebuild_family = subset(SummaryProjectionScenarios_withRebuild,SummaryProjectionScenarios_withRebuild$FamilyCommonName==FamilyNames[i])
  png(paste(PATH_highLevelSummaryPlots,paste("RebuildingTimes_Histogram_",FamilyNames[i],".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)
  hist(SummaryProjectionScenarios_withRebuild_family$TotalRebuildTime,freq=FALSE,xlab="Rebuilding Times (years)",main="Rebuilding Times",breaks=seq(0,ceiling(max(SummaryProjectionScenarios_withRebuild_family$TotalRebuildTime)+1),by=1))
  mtext(FamilyNames[i])
  dev.off()
}


################# Develop Large Scale Kobe Plot Summaries: ####################################
##******** Overfishing Defined as F below SPR 40% and overfished being defined as B below SPR 20% ******;
#Calculate Percent Change in F table
CurrentStatusValues = read.table(paste(PATH_output,"CurrentValuesBenchmarksAndRatios.csv",sep="/"),header=TRUE,sep=",",as.is=TRUE)
GenusSpeciesCommonName = subset(Lengths_cpue,select=c(Genus,Species,FamilyCommonName))
GenusSpeciesCommonName = unique(GenusSpeciesCommonName)
CurrentStatusValues = merge(CurrentStatusValues,GenusSpeciesCommonName,all.x=TRUE,by=c("Genus","Species"))

#Remove Jack and Mackerel Species and also L. malabaricus in 713 
toRemove = data.frame(Genus=c("Carangoides","Carangoides","Carangoides","Caranx","Caranx","Seriola"),Species=c("chrysophrys","coeruleopinnatus","gymnostethus","bucculentus","sexfasciatus","rivoliana"),flag=1)
CurrentStatusValues = merge(CurrentStatusValues,toRemove,all.x=TRUE,by=c("Genus","Species"))
CurrentStatusValues = subset(CurrentStatusValues,is.na(CurrentStatusValues$flag))
CurrentStatusValues$flag[CurrentStatusValues$Genus=="Lutjanus" & CurrentStatusValues$Species=="malabaricus" & CurrentStatusValues$WPP=="713"]=1
CurrentStatusValues = subset(CurrentStatusValues,is.na(CurrentStatusValues$flag))
CurrentStatusValues$flag=NULL                          

#Do NOT include SolveBaranov or TropFishR (or LIME for biomass benchmar estiamtes)
CurrentStatusValues = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="SolveBaranov",]   #TEMP CODE LINE: do NOT include SolveBaranov or TropFishR
CurrentStatusValues = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="TropFishR",]      #TEMP CODE LINE: do NOT include SolveBaranov or TropFishR
CurrentStatusValues_noLIME = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="LIME",]

#Identify and filter away the Runs that did not converge
CurrentStatusValues$F_ratio20_FLAG=0
CurrentStatusValues$F_ratio20_FLAG[CurrentStatusValues$F_ratio20 > 5 | CurrentStatusValues$F_ratio20 < 0.01]=1
CurrentStatusValues$F_ratio30_FLAG=0
CurrentStatusValues$F_ratio30_FLAG[CurrentStatusValues$F_ratio30 > 5 | CurrentStatusValues$F_ratio30 < 0.01]=1
CurrentStatusValues$F_ratio40_FLAG=0
CurrentStatusValues$F_ratio40_FLAG[CurrentStatusValues$F_ratio40 > 5 | CurrentStatusValues$F_ratio40 < 0.01]=1
CurrentStatusValues$B_ratio20_FLAG=0
CurrentStatusValues$B_ratio20_FLAG[CurrentStatusValues$B_ratio20 > 5 | CurrentStatusValues$B_ratio20 < 0.01]=1
CurrentStatusValues$B_ratio30_FLAG=0
CurrentStatusValues$B_ratio30_FLAG[CurrentStatusValues$B_ratio30 > 5 | CurrentStatusValues$B_ratio30 < 0.01]=1
CurrentStatusValues$B_ratio40_FLAG=0
CurrentStatusValues$B_ratio40_FLAG[CurrentStatusValues$B_ratio40 > 5 | CurrentStatusValues$B_ratio40 < 0.01]=1
CurrentStatusValues_noLIME$F_ratio20_FLAG=0
CurrentStatusValues_noLIME$F_ratio20_FLAG[CurrentStatusValues_noLIME$F_ratio20 > 5 | CurrentStatusValues_noLIME$F_ratio20 < 0.01]=1
CurrentStatusValues_noLIME$F_ratio30_FLAG=0
CurrentStatusValues_noLIME$F_ratio30_FLAG[CurrentStatusValues_noLIME$F_ratio30 > 5 | CurrentStatusValues_noLIME$F_ratio30 < 0.01]=1
CurrentStatusValues_noLIME$F_ratio40_FLAG=0
CurrentStatusValues_noLIME$F_ratio40_FLAG[CurrentStatusValues_noLIME$F_ratio40 > 5 | CurrentStatusValues_noLIME$F_ratio40 < 0.01]=1
CurrentStatusValues_noLIME$B_ratio20_FLAG=0
CurrentStatusValues_noLIME$B_ratio20_FLAG[CurrentStatusValues_noLIME$B_ratio20 > 5 | CurrentStatusValues_noLIME$B_ratio20 < 0.01]=1
CurrentStatusValues_noLIME$B_ratio30_FLAG=0
CurrentStatusValues_noLIME$B_ratio30_FLAG[CurrentStatusValues_noLIME$B_ratio30 > 5 | CurrentStatusValues_noLIME$B_ratio30 < 0.01]=1
CurrentStatusValues_noLIME$B_ratio40_FLAG=0
CurrentStatusValues_noLIME$B_ratio40_FLAG[CurrentStatusValues_noLIME$B_ratio40 > 5 | CurrentStatusValues_noLIME$B_ratio40 < 0.01]=1
CurrentStatusValues_F_ratio20 = subset(CurrentStatusValues,CurrentStatusValues$F_ratio20_FLAG==0)
CurrentStatusValues_F_ratio30 = subset(CurrentStatusValues,CurrentStatusValues$F_ratio30_FLAG==0)
CurrentStatusValues_F_ratio40 = subset(CurrentStatusValues,CurrentStatusValues$F_ratio40_FLAG==0)
CurrentStatusValues_B_ratio20 = subset(CurrentStatusValues,CurrentStatusValues$B_ratio20_FLAG==0)
CurrentStatusValues_B_ratio30 = subset(CurrentStatusValues,CurrentStatusValues$B_ratio30_FLAG==0)
CurrentStatusValues_B_ratio40 = subset(CurrentStatusValues,CurrentStatusValues$B_ratio40_FLAG==0)
CurrentStatusValues_noLIME_F_ratio20 = subset(CurrentStatusValues_noLIME,CurrentStatusValues_noLIME$F_ratio20_FLAG==0)
CurrentStatusValues_noLIME_F_ratio30 = subset(CurrentStatusValues_noLIME,CurrentStatusValues_noLIME$F_ratio30_FLAG==0)
CurrentStatusValues_noLIME_F_ratio40 = subset(CurrentStatusValues_noLIME,CurrentStatusValues_noLIME$F_ratio40_FLAG==0)
CurrentStatusValues_noLIME_B_ratio20 = subset(CurrentStatusValues_noLIME,CurrentStatusValues_noLIME$B_ratio20_FLAG==0)
CurrentStatusValues_noLIME_B_ratio30 = subset(CurrentStatusValues_noLIME,CurrentStatusValues_noLIME$B_ratio30_FLAG==0)
CurrentStatusValues_noLIME_B_ratio40 = subset(CurrentStatusValues_noLIME,CurrentStatusValues_noLIME$B_ratio40_FLAG==0)

CurrentStatusValues$GenusSpecies = paste(CurrentStatusValues$Genus,CurrentStatusValues$Species,sep="_")
CurrentStatusValues$F_ratioAboveOne = 0
CurrentStatusValues$F_ratioAboveOne[CurrentStatusValues$F_ratio40>=1]=1
CurrentStatusValues$B_ratioAboveOne = 0
CurrentStatusValues$B_ratioAboveOne[CurrentStatusValues$B_ratio20>=1]=1
CurrentStatusValues$Quadrant=0      #Quadrant 1=Overfished and Overfishing; 2=Overfishing, not overfished; 3=not overfishing, not overfished; 4=not overfishing, overfished
CurrentStatusValues$Quadrant[CurrentStatusValues$F_ratio40>=1 & CurrentStatusValues$B_ratio20<1]=1
CurrentStatusValues$Quadrant[CurrentStatusValues$F_ratio40>=1 & CurrentStatusValues$B_ratio20>=1]=2
CurrentStatusValues$Quadrant[CurrentStatusValues$F_ratio40<1 & CurrentStatusValues$B_ratio20>=1]=3
CurrentStatusValues$Quadrant[CurrentStatusValues$F_ratio40<1 & CurrentStatusValues$B_ratio20<1]=4
StatusTable = data.frame(table(CurrentStatusValues$FamilyCommonName,CurrentStatusValues$WPP,CurrentStatusValues$Quadrant))
names(StatusTable) = c("FamilyCommonName","WPP","Quadrant","Freq")
StatusTable$FamilyCommonName = as.character(StatusTable$FamilyCommonName)
StatusTable$WPP = as.character(StatusTable$WPP)
StatusTable$Quadrant = as.numeric(as.character(StatusTable$Quadrant))
StatusTable = reshape(StatusTable,v.names="Freq",idvar=c("FamilyCommonName","WPP"),timevar="Quadrant",direction="wide")
StatusTable$Total = StatusTable$Freq.1 + StatusTable$Freq.2 + StatusTable$Freq.3 + StatusTable$Freq.4
StatusTable$Percent.1 = round((StatusTable$Freq.1/StatusTable$Total)*100)
StatusTable$Percent.2 = round((StatusTable$Freq.2/StatusTable$Total)*100)
StatusTable$Percent.3 = round((StatusTable$Freq.3/StatusTable$Total)*100)
StatusTable$Percent.4 = round((StatusTable$Freq.4/StatusTable$Total)*100)
StatusTable[is.na(StatusTable)]=0
StatusTable$FamilyCommonName = as.character(StatusTable$FamilyCommonName)
StatusTable$WPP = as.character(StatusTable$WPP)
write.table(StatusTable,paste(PATH_output,"StatusTable_Family.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
#Now make the Kobe plots
#CurrentStatusValues = read.table(paste(PATH_output,"CurrentValuesBenchmarksAndRatios.csv",sep="/"),header=TRUE,sep=",")
#CurrentStatusValues = merge(CurrentStatusValues,GenusSpeciesCommonName,all.x=TRUE,by=c("Genus","Species"))
#Remove runs that you don't want: here remove SolveBaranov and VPA and LIME
#CurrentStatusValues = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="SolveBaranov",]
#CurrentStatusValues = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="TropFishR",]
#CurrentStatusValues = CurrentStatusValues[CurrentStatusValues$populationReconstructionMethod!="LIME",]
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
  thisIteration = CurrentStatusValues[CurrentStatusValues$FamilyCommonName==BenchmarksToPlot$FamilyCommonName[i] & as.character(CurrentStatusValues$WPP)==as.character(BenchmarksToPlot$WPP[i]),]
  png(paste(PATH_highLevelSummaryPlots,paste(paste(BenchmarksToPlot$FamilyCommonName[i],BenchmarksToPlot$WPP[i],sep="_"),"_Benchmarks.png",sep=""),sep="/"),units="px",width=5800,height=3200,res=600)
  layout(matrix(c(1,2), nrow = 1), widths = c(1,0.5))
  par(mar = c(5, 4, 4, 2) + 0.1)
  maxBratio = max(thisIteration$B_ratio20)
  maxFratio = max(thisIteration$F_ratio40)
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
      points(thisIteration$B_ratio20[j],thisIteration$F_ratio40[j],pch=thisIteration$pch_arg[j],col=thisIteration$Colors[j])
    }
    if(j>1)
    {
      points(thisIteration$B_ratio20[j],thisIteration$F_ratio40[j],pch=thisIteration$pch_arg[j],col=thisIteration$Colors[j])
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
  StatusPercents = StatusTable[StatusTable$FamilyCommonName==as.character(BenchmarksToPlot$FamilyCommonName[i]) & StatusTable$WPP==BenchmarksToPlot$WPP[i],]
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
  if(BenchmarksToPlot$WPP[i]=="EEZ")
  {
    mtext("Data Pooled Across EEZ")
  }
  if(BenchmarksToPlot$WPP[i]!="EEZ")
  {
    mtext(paste("FMA = ",BenchmarksToPlot$WPP[i],sep=""))
  }                                                        
  box()
  par(mar = c(5, 0, 4, 2))
  plot(thisIteration$B_ratio[dim(thisIteration)[1]],thisIteration$F_ratio[dim(thisIteration)[1]],type="n",axes=FALSE,ann=FALSE,xlim=c(0,maxBratio),ylim=c(0,maxFratio))
  LengthReconstructionSymbols_forLegend = LengthReconstructionSymbols
  LengthReconstructionSymbols_forLegend$populationReconstructionMethod = as.character(LengthReconstructionSymbols_forLegend$populationReconstructionMethod)
  LengthReconstructionSymbols_forLegend$populationReconstructionMethod[LengthReconstructionSymbols_forLegend$populationReconstructionMethod=="BevertonHoltInstantaneousMortality"]="BevHoltF"
  legend("topleft",legend=c(as.character(LengthReconstructionSymbols_forLegend$populationReconstructionMethod),as.character(ColorsFor_benchmarkPlot_CmsyScenarios$CatchMSY_Scenario)),lty=c(rep(NA,times=dim(LengthReconstructionSymbols_forLegend)[1]),rep(1,times=dim(ColorsFor_benchmarkPlot_CmsyScenarios)[1])),lwd=c(rep(NA,times=dim(LengthReconstructionSymbols_forLegend)[1]),rep(2,times=dim(ColorsFor_benchmarkPlot_CmsyScenarios)[1])),col=c(rep("black",times=dim(LengthReconstructionSymbols_forLegend)[1]),ColorsFor_benchmarkPlot_CmsyScenarios$Colors),pch=c(LengthReconstructionSymbols_forLegend$pch_arg,rep(NA,times=dim(ColorsFor_benchmarkPlot_CmsyScenarios)[1])),cex=0.8)
  dev.off()
}

########################## 5/25/2022: Calculate Forgone Catch - i.e. those that are under and over MSY #####################

SummaryProjectionScenarios = read.table(paste(PATH_output,"SummaryProjectionScenarios.csv",sep="/"),header=TRUE,sep=",",as.is=TRUE)
SummaryProjectionScenarios = SummaryProjectionScenarios[SummaryProjectionScenarios$populationReconstructionMethod!="SolveBaranov",]   #TEMP CODE LINE: do NOT include SolveBaranov or TropFishR
SummaryProjectionScenarios = SummaryProjectionScenarios[SummaryProjectionScenarios$populationReconstructionMethod!="TropFishR",]      #TEMP CODE LINE: do NOT include SolveBaranov or TropFishR
SummaryProjectionScenarios = SummaryProjectionScenarios[SummaryProjectionScenarios$ProjectionType=="RegularF",]  #Use only the regular F based projections to calculate benchmarks; the differnet recruitment assumptions made by Peter should not be included here - benchmarks can not be calculated from those;
#Remove Jack and Mackerel Species and also L. malabaricus in 713
toRemove = data.frame(Genus=c("Carangoides","Carangoides","Carangoides","Caranx","Caranx","Seriola"),Species=c("chrysophrys","coeruleopinnatus","gymnostethus","bucculentus","sexfasciatus","rivoliana"),flag=1)
SummaryProjectionScenarios = merge(SummaryProjectionScenarios,toRemove,all.x=TRUE,by=c("Genus","Species"))
SummaryProjectionScenarios = subset(SummaryProjectionScenarios,is.na(SummaryProjectionScenarios$flag))
SummaryProjectionScenarios$flag[SummaryProjectionScenarios$Genus=="Lutjanus" & SummaryProjectionScenarios$Species=="malabaricus" & SummaryProjectionScenarios$WPP=="713"]=1
SummaryProjectionScenarios = subset(SummaryProjectionScenarios,is.na(SummaryProjectionScenarios$flag))
SummaryProjectionScenarios$flag=NULL
#current biomass and catch
CatchBiomass_current = SummaryProjectionScenarios[SummaryProjectionScenarios$FishingMortaltiyProjectionScenario=="CurrentF",]
CatchBiomass_current_noLime = CatchBiomass_current[CatchBiomass_current$populationReconstructionMethod!="LIME",]
CatchBiomass_current_noLime = subset(CatchBiomass_current_noLime,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,ProjectionType,Yield_estimated,B_estimated))
names(CatchBiomass_current_noLime)[names(CatchBiomass_current_noLime)=="Yield_estimated"]="CurrentCatchEquil"
names(CatchBiomass_current_noLime)[names(CatchBiomass_current_noLime)=="B_estimated"]="CurrentB"
#biomass and catch at SPR20
CatchBiomass_atSPR20 = SummaryProjectionScenarios[SummaryProjectionScenarios$FishingMortaltiyProjectionScenario=="F_at_SPR20",]
CatchBiomass_atSPR20 = CatchBiomass_atSPR20[CatchBiomass_atSPR20$populationReconstructionMethod!="LIME",]
CatchBiomass_atSPR20 = subset(CatchBiomass_atSPR20,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,ProjectionType,CatchAtEnd,BiomassAtEnd))
names(CatchBiomass_atSPR20)[names(CatchBiomass_atSPR20)=="CatchAtEnd"]="CatchSPR20"
names(CatchBiomass_atSPR20)[names(CatchBiomass_atSPR20)=="BiomassAtEnd"]="Bspr20"
#biomass and catch at SPR30
CatchBiomass_atSPR30 = SummaryProjectionScenarios[SummaryProjectionScenarios$FishingMortaltiyProjectionScenario=="F_at_SPR30",]
CatchBiomass_atSPR30 = CatchBiomass_atSPR30[CatchBiomass_atSPR30$populationReconstructionMethod!="LIME",]
CatchBiomass_atSPR30 = subset(CatchBiomass_atSPR30,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,ProjectionType,CatchAtEnd,BiomassAtEnd))
names(CatchBiomass_atSPR30)[names(CatchBiomass_atSPR30)=="CatchAtEnd"]="CatchSPR30"
names(CatchBiomass_atSPR30)[names(CatchBiomass_atSPR30)=="BiomassAtEnd"]="Bspr30"
#biomass and catch at SPR40
CatchBiomass_atSPR40 = SummaryProjectionScenarios[SummaryProjectionScenarios$FishingMortaltiyProjectionScenario=="F_at_SPR40",]
CatchBiomass_atSPR40 = CatchBiomass_atSPR40[CatchBiomass_atSPR40$populationReconstructionMethod!="LIME",]
CatchBiomass_atSPR40 = subset(CatchBiomass_atSPR40,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,ProjectionType,CatchAtEnd,BiomassAtEnd))
names(CatchBiomass_atSPR40)[names(CatchBiomass_atSPR40)=="CatchAtEnd"]="CatchSPR40"
names(CatchBiomass_atSPR40)[names(CatchBiomass_atSPR40)=="BiomassAtEnd"]="Bspr40"
#Merge, calculate catch and biomass ratios, and filter the results to remove model runs that clearly did not converge; use the same criterian as used in the Excel spreadsheet
#NOTE: have this table written out and use it in the paper to replace the current table; have this as a standard output of the code. 
CatchBiomass = merge(CatchBiomass_current_noLime,CatchBiomass_atSPR20,all.x=TRUE,by=c("Genus","Species","WPP",
  "CatchMSY_Scenario","CatchMSY_Scenario_Bound","LifeHistory_Scenario","populationReconstructionMethod","SteepnessAssumed","ProjectionType"))
CatchBiomass = merge(CatchBiomass,CatchBiomass_atSPR30,all.x=TRUE,by=c("Genus","Species","WPP",
  "CatchMSY_Scenario","CatchMSY_Scenario_Bound","LifeHistory_Scenario","populationReconstructionMethod","SteepnessAssumed","ProjectionType"))
CatchBiomass = merge(CatchBiomass,CatchBiomass_atSPR40,all.x=TRUE,by=c("Genus","Species","WPP",
  "CatchMSY_Scenario","CatchMSY_Scenario_Bound","LifeHistory_Scenario","populationReconstructionMethod","SteepnessAssumed","ProjectionType"))
CatchBiomass$CatchRatioSPR20 = CatchBiomass$CurrentCatchEquil/CatchBiomass$CatchSPR20
CatchBiomass$CatchRatioSPR30 = CatchBiomass$CurrentCatchEquil/CatchBiomass$CatchSPR30
CatchBiomass$CatchRatioSPR40 = CatchBiomass$CurrentCatchEquil/CatchBiomass$CatchSPR40
CatchBiomass$CatchRatioSPR20_flag=0
CatchBiomass$CatchRatioSPR20_flag[CatchBiomass$CatchRatioSPR20>5 | CatchBiomass$CatchRatioSPR20<0.01]=1
CatchBiomass$CatchRatioSPR30_flag=0
CatchBiomass$CatchRatioSPR30_flag[CatchBiomass$CatchRatioSPR30>5 | CatchBiomass$CatchRatioSPR30<0.01]=1
CatchBiomass$CatchRatioSPR40_flag=0
CatchBiomass$CatchRatioSPR40_flag[CatchBiomass$CatchRatioSPR40>5 | CatchBiomass$CatchRatioSPR40<0.01]=1
CatchBiomass$BratioSPR20 = CatchBiomass$CurrentB/CatchBiomass$Bspr20
CatchBiomass$BratioSPR30 = CatchBiomass$CurrentB/CatchBiomass$Bspr30
CatchBiomass$BratioSPR40 = CatchBiomass$CurrentB/CatchBiomass$Bspr40
CatchBiomass$BratioSPR20_flag=0
CatchBiomass$BratioSPR20_flag[CatchBiomass$BratioSPR20>5 | CatchBiomass$BratioSPR20<0.01]=1
CatchBiomass$BratioSPR30_flag=0
CatchBiomass$BratioSPR30_flag[CatchBiomass$BratioSPR30>5 | CatchBiomass$BratioSPR30<0.01]=1
CatchBiomass$BratioSPR40_flag=0
CatchBiomass$BratioSPR40_flag[CatchBiomass$BratioSPR40>5 | CatchBiomass$BratioSPR40<0.01]=1 
CatchBiomass_catchRatioSPR20 = CatchBiomass[CatchBiomass$CatchRatioSPR20_flag==0,]
CatchBiomass_catchRatioSPR20_mean = aggregate.data.frame(CatchBiomass_catchRatioSPR20$CatchRatioSPR20,by=list(CatchBiomass_catchRatioSPR20$Genus,CatchBiomass_catchRatioSPR20$Species,CatchBiomass_catchRatioSPR20$WPP),FUN=mean,na.rm=TRUE)
names(CatchBiomass_catchRatioSPR20_mean) = c("Genus","Species","WPP","CatchRatioSPR20")
i = sapply(CatchBiomass_catchRatioSPR20_mean,is.factor)
CatchBiomass_catchRatioSPR20_mean[i] = lapply(CatchBiomass_catchRatioSPR20_mean[i],as.character)
CatchBiomass_catchRatioSPR20_sd = aggregate.data.frame(CatchBiomass_catchRatioSPR20$CatchRatioSPR20,by=list(CatchBiomass_catchRatioSPR20$Genus,CatchBiomass_catchRatioSPR20$Species,CatchBiomass_catchRatioSPR20$WPP),FUN=sd,na.rm=TRUE)
  names(CatchBiomass_catchRatioSPR20_sd) = c("Genus","Species","WPP","CatchRatioSPR20_sd")
i = sapply(CatchBiomass_catchRatioSPR20_sd,is.factor)
CatchBiomass_catchRatioSPR20_sd[i] = lapply(CatchBiomass_catchRatioSPR20_sd[i],as.character)
CatchBiomass_catchRatioSPR20 = merge(CatchBiomass_catchRatioSPR20_mean,CatchBiomass_catchRatioSPR20_sd,all.x=TRUE,by=c("Genus","Species","WPP"))
CatchBiomass_catchRatioSPR30 = CatchBiomass[CatchBiomass$CatchRatioSPR30_flag==0,]
CatchBiomass_catchRatioSPR30_mean = aggregate.data.frame(CatchBiomass_catchRatioSPR30$CatchRatioSPR30,by=list(CatchBiomass_catchRatioSPR30$Genus,CatchBiomass_catchRatioSPR30$Species,CatchBiomass_catchRatioSPR30$WPP),FUN=mean,na.rm=TRUE)
names(CatchBiomass_catchRatioSPR30_mean) = c("Genus","Species","WPP","CatchRatioSPR30")
i = sapply(CatchBiomass_catchRatioSPR30_mean,is.factor)
CatchBiomass_catchRatioSPR30_mean[i] = lapply(CatchBiomass_catchRatioSPR30_mean[i],as.character)
CatchBiomass_catchRatioSPR30_sd = aggregate.data.frame(CatchBiomass_catchRatioSPR30$CatchRatioSPR30,by=list(CatchBiomass_catchRatioSPR30$Genus,CatchBiomass_catchRatioSPR30$Species,CatchBiomass_catchRatioSPR30$WPP),FUN=sd,na.rm=TRUE)
names(CatchBiomass_catchRatioSPR30_sd) = c("Genus","Species","WPP","CatchRatioSPR30_sd")
i = sapply(CatchBiomass_catchRatioSPR30_sd,is.factor)
CatchBiomass_catchRatioSPR30_sd[i] = lapply(CatchBiomass_catchRatioSPR30_sd[i],as.character)
CatchBiomass_catchRatioSPR30 = merge(CatchBiomass_catchRatioSPR30_mean,CatchBiomass_catchRatioSPR30_sd,all.x=TRUE,by=c("Genus","Species","WPP"))
CatchBiomass_catchRatioSPR40 = CatchBiomass[CatchBiomass$CatchRatioSPR40_flag==0,]
CatchBiomass_catchRatioSPR40_mean = aggregate.data.frame(CatchBiomass_catchRatioSPR40$CatchRatioSPR40,by=list(CatchBiomass_catchRatioSPR40$Genus,CatchBiomass_catchRatioSPR40$Species,CatchBiomass_catchRatioSPR40$WPP),FUN=mean,na.rm=TRUE)
names(CatchBiomass_catchRatioSPR40_mean) = c("Genus","Species","WPP","CatchRatioSPR40")
i = sapply(CatchBiomass_catchRatioSPR40_mean,is.factor)
CatchBiomass_catchRatioSPR40_mean[i] = lapply(CatchBiomass_catchRatioSPR40_mean[i],as.character)
CatchBiomass_catchRatioSPR40_sd = aggregate.data.frame(CatchBiomass_catchRatioSPR40$CatchRatioSPR40,by=list(CatchBiomass_catchRatioSPR40$Genus,CatchBiomass_catchRatioSPR40$Species,CatchBiomass_catchRatioSPR40$WPP),FUN=sd,na.rm=TRUE)
names(CatchBiomass_catchRatioSPR40_sd) = c("Genus","Species","WPP","CatchRatioSPR40_sd")
i = sapply(CatchBiomass_catchRatioSPR40_sd,is.factor)
CatchBiomass_catchRatioSPR40_sd[i] = lapply(CatchBiomass_catchRatioSPR40_sd[i],as.character)
CatchBiomass_catchRatioSPR40 = merge(CatchBiomass_catchRatioSPR40_mean,CatchBiomass_catchRatioSPR40_sd,all.x=TRUE,by=c("Genus","Species","WPP"))
CatchBiomass_BratioSPR20 = CatchBiomass[CatchBiomass$BratioSPR20_flag==0,]
CatchBiomass_BratioSPR20_mean = aggregate.data.frame(CatchBiomass_BratioSPR20$BratioSPR20,by=list(CatchBiomass_BratioSPR20$Genus,CatchBiomass_BratioSPR20$Species,CatchBiomass_BratioSPR20$WPP),FUN=mean,na.rm=TRUE)
names(CatchBiomass_BratioSPR20_mean) = c("Genus","Species","WPP","BratioSPR20")
i = sapply(CatchBiomass_BratioSPR20_mean,is.factor)
CatchBiomass_BratioSPR20_mean[i] = lapply(CatchBiomass_BratioSPR20_mean[i],as.character)
CatchBiomass_BratioSPR20_sd = aggregate.data.frame(CatchBiomass_BratioSPR20$BratioSPR20,by=list(CatchBiomass_BratioSPR20$Genus,CatchBiomass_BratioSPR20$Species,CatchBiomass_BratioSPR20$WPP),FUN=sd,na.rm=TRUE)
names(CatchBiomass_BratioSPR20_sd) = c("Genus","Species","WPP","BratioSPR20_sd")
i = sapply(CatchBiomass_BratioSPR20_sd,is.factor)
CatchBiomass_BratioSPR20_sd[i] = lapply(CatchBiomass_BratioSPR20_sd[i],as.character)
CatchBiomass_BratioSPR20 = merge(CatchBiomass_BratioSPR20_mean,CatchBiomass_BratioSPR20_sd,all.x=TRUE,by=c("Genus","Species","WPP"))
CatchBiomass_BratioSPR30 = CatchBiomass[CatchBiomass$BratioSPR30_flag==0,]
CatchBiomass_BratioSPR30_mean = aggregate.data.frame(CatchBiomass_BratioSPR30$BratioSPR30,by=list(CatchBiomass_BratioSPR30$Genus,CatchBiomass_BratioSPR30$Species,CatchBiomass_BratioSPR30$WPP),FUN=mean,na.rm=TRUE)
names(CatchBiomass_BratioSPR30_mean) = c("Genus","Species","WPP","BratioSPR30")
i = sapply(CatchBiomass_BratioSPR30_mean,is.factor)
CatchBiomass_BratioSPR30_mean[i] = lapply(CatchBiomass_BratioSPR30_mean[i],as.character)
CatchBiomass_BratioSPR30_sd = aggregate.data.frame(CatchBiomass_BratioSPR30$BratioSPR30,by=list(CatchBiomass_BratioSPR30$Genus,CatchBiomass_BratioSPR30$Species,CatchBiomass_BratioSPR30$WPP),FUN=sd,na.rm=TRUE)
names(CatchBiomass_BratioSPR30_sd) = c("Genus","Species","WPP","BratioSPR30_sd")
i = sapply(CatchBiomass_BratioSPR30_sd,is.factor)
CatchBiomass_BratioSPR30_sd[i] = lapply(CatchBiomass_BratioSPR30_sd[i],as.character)
CatchBiomass_BratioSPR30 = merge(CatchBiomass_BratioSPR30_mean,CatchBiomass_BratioSPR30_sd,all.x=TRUE,by=c("Genus","Species","WPP"))
CatchBiomass_BratioSPR40 = CatchBiomass[CatchBiomass$BratioSPR40_flag==0,]
CatchBiomass_BratioSPR40_mean = aggregate.data.frame(CatchBiomass_BratioSPR40$BratioSPR40,by=list(CatchBiomass_BratioSPR40$Genus,CatchBiomass_BratioSPR40$Species,CatchBiomass_BratioSPR40$WPP),FUN=mean,na.rm=TRUE)
names(CatchBiomass_BratioSPR40_mean) = c("Genus","Species","WPP","BratioSPR40")
i = sapply(CatchBiomass_BratioSPR40_mean,is.factor)
CatchBiomass_BratioSPR40_mean[i] = lapply(CatchBiomass_BratioSPR40_mean[i],as.character)
CatchBiomass_BratioSPR40_sd = aggregate.data.frame(CatchBiomass_BratioSPR40$BratioSPR40,by=list(CatchBiomass_BratioSPR40$Genus,CatchBiomass_BratioSPR40$Species,CatchBiomass_BratioSPR40$WPP),FUN=sd,na.rm=TRUE)
names(CatchBiomass_BratioSPR40_sd) = c("Genus","Species","WPP","BratioSPR40_sd")
i = sapply(CatchBiomass_BratioSPR40_sd,is.factor)
CatchBiomass_BratioSPR40_sd[i] = lapply(CatchBiomass_BratioSPR40_sd[i],as.character)
CatchBiomass_BratioSPR40 = merge(CatchBiomass_BratioSPR40_mean,CatchBiomass_BratioSPR40_sd,all.x=TRUE,by=c("Genus","Species","WPP"))
CatchBiomass = merge(CatchBiomass_catchRatioSPR20,CatchBiomass_catchRatioSPR30,all.x=TRUE,by=c("Genus","Species","WPP"))
CatchBiomass = merge(CatchBiomass,CatchBiomass_catchRatioSPR40,all.x=TRUE,by=c("Genus","Species","WPP"))
CatchBiomass = merge(CatchBiomass,CatchBiomass_BratioSPR20,all.x=TRUE,by=c("Genus","Species","WPP"))
CatchBiomass = merge(CatchBiomass,CatchBiomass_BratioSPR30,all.x=TRUE,by=c("Genus","Species","WPP"))
CatchBiomass = merge(CatchBiomass,CatchBiomass_BratioSPR40,all.x=TRUE,by=c("Genus","Species","WPP"))
CatchBiomass$CatchRatioSPR20_Lower = CatchBiomass$CatchRatioSPR20 - CatchBiomass$CatchRatioSPR20_sd
CatchBiomass$CatchRatioSPR20_Lower[CatchBiomass$CatchRatioSPR20_Lower<0]=0
CatchBiomass$CatchRatioSPR20_Upper = CatchBiomass$CatchRatioSPR20 + CatchBiomass$CatchRatioSPR20_sd
CatchBiomass$CatchRatioSPR30_Lower = CatchBiomass$CatchRatioSPR30 - CatchBiomass$CatchRatioSPR30_sd
CatchBiomass$CatchRatioSPR30_Lower[CatchBiomass$CatchRatioSPR30_Lower<0]=0
CatchBiomass$CatchRatioSPR30_Upper = CatchBiomass$CatchRatioSPR30 + CatchBiomass$CatchRatioSPR30_sd
CatchBiomass$CatchRatioSPR40_Lower = CatchBiomass$CatchRatioSPR40 - CatchBiomass$CatchRatioSPR40_sd
CatchBiomass$CatchRatioSPR40_Lower[CatchBiomass$CatchRatioSPR40_Lower<0]=0
CatchBiomass$CatchRatioSPR40_Upper = CatchBiomass$CatchRatioSPR40 + CatchBiomass$CatchRatioSPR40_sd
CatchBiomass$BratioSPR20_Lower = CatchBiomass$BratioSPR20 - CatchBiomass$BratioSPR20_sd
CatchBiomass$BratioSPR20_Lower[CatchBiomass$BratioSPR20_Lower<0]=0
CatchBiomass$BratioSPR20_Upper = CatchBiomass$BratioSPR20 + CatchBiomass$BratioSPR20_sd
CatchBiomass$BratioSPR30_Lower = CatchBiomass$BratioSPR30 - CatchBiomass$BratioSPR30_sd
CatchBiomass$BratioSPR30_Lower[CatchBiomass$BratioSPR30_Lower<0]=0
CatchBiomass$BratioSPR30_Upper = CatchBiomass$BratioSPR30 + CatchBiomass$BratioSPR30_sd
CatchBiomass$BratioSPR40_Lower = CatchBiomass$BratioSPR40 - CatchBiomass$BratioSPR40_sd
CatchBiomass$BratioSPR40_Lower[CatchBiomass$BratioSPR40_Lower<0]=0
CatchBiomass$BratioSPR40_Upper = CatchBiomass$BratioSPR40 + CatchBiomass$BratioSPR40_sd
#current F
FishingMort_current = SummaryProjectionScenarios[SummaryProjectionScenarios$FishingMortaltiyProjectionScenario=="CurrentF",]
FishingMort_current = subset(FishingMort_current,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,ProjectionType,F_estimated))
names(FishingMort_current)[names(FishingMort_current)=="F_estimated"]="CurrentF"
#F at SPR20
FishingMort_atSPR20 = SummaryProjectionScenarios[SummaryProjectionScenarios$FishingMortaltiyProjectionScenario=="F_at_SPR20",]
FishingMort_atSPR20 = subset(FishingMort_atSPR20,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,ProjectionType,F_project))
names(FishingMort_atSPR20)[names(FishingMort_atSPR20)=="F_project"]="Fspr20"
#F at SPR30
FishingMort_atSPR30 = SummaryProjectionScenarios[SummaryProjectionScenarios$FishingMortaltiyProjectionScenario=="F_at_SPR30",]
FishingMort_atSPR30 = subset(FishingMort_atSPR30,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,ProjectionType,F_project))
names(FishingMort_atSPR30)[names(FishingMort_atSPR30)=="F_project"]="Fspr30"
#F at SPR40
FishingMort_atSPR40 = SummaryProjectionScenarios[SummaryProjectionScenarios$FishingMortaltiyProjectionScenario=="F_at_SPR40",]
FishingMort_atSPR40 = subset(FishingMort_atSPR40,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,ProjectionType,F_project))
names(FishingMort_atSPR40)[names(FishingMort_atSPR40)=="F_project"]="Fspr40"
#Merge, calculate F ratios, and filter the results to remove model runs that clearly did not converge; use the same criterian as used in the Excel spreadsheet
#NOTE: have this table written out and use it in the paper to replace the current table; have this as a standard output of the code.
FishingMort = merge(FishingMort_current,FishingMort_atSPR20,all.x=TRUE,by=c("Genus","Species","WPP",
  "CatchMSY_Scenario","CatchMSY_Scenario_Bound","LifeHistory_Scenario","populationReconstructionMethod","SteepnessAssumed","ProjectionType"))
FishingMort = merge(FishingMort,FishingMort_atSPR30,all.x=TRUE,by=c("Genus","Species","WPP",
  "CatchMSY_Scenario","CatchMSY_Scenario_Bound","LifeHistory_Scenario","populationReconstructionMethod","SteepnessAssumed","ProjectionType"))
FishingMort = merge(FishingMort,FishingMort_atSPR40,all.x=TRUE,by=c("Genus","Species","WPP",
  "CatchMSY_Scenario","CatchMSY_Scenario_Bound","LifeHistory_Scenario","populationReconstructionMethod","SteepnessAssumed","ProjectionType"))
FishingMort$SteepnessAssumed=NULL
FishingMort$ProjectionType=NULL
FishingMort = unique(FishingMort)
FishingMort$FratioSPR20 = FishingMort$CurrentF/FishingMort$Fspr20
FishingMort$FratioSPR30 = FishingMort$CurrentF/FishingMort$Fspr30
FishingMort$FratioSPR40 = FishingMort$CurrentF/FishingMort$Fspr40
FishingMort$FratioSPR20_flag=0
FishingMort$FratioSPR20_flag[FishingMort$FratioSPR20>5 | FishingMort$FratioSPR20<0.01]=1
FishingMort$FratioSPR30_flag=0
FishingMort$FratioSPR30_flag[FishingMort$FratioSPR30>5 | FishingMort$FratioSPR30<0.01]=1
FishingMort$FratioSPR40_flag=0
FishingMort$FratioSPR40_flag[FishingMort$FratioSPR40>5 | FishingMort$FratioSPR40<0.01]=1
#mean and sd for F@spr20
FishingMort_FratioSPR20 = FishingMort[FishingMort$FratioSPR20_flag==0,]
FishingMort_FratioSPR20_mean = aggregate.data.frame(FishingMort_FratioSPR20$FratioSPR20,by=list(FishingMort_FratioSPR20$Genus,FishingMort_FratioSPR20$Species,FishingMort_FratioSPR20$WPP),FUN=mean,na.rm=TRUE)
names(FishingMort_FratioSPR20_mean) = c("Genus","Species","WPP","FratioSPR20")
i = sapply(FishingMort_FratioSPR20_mean,is.factor)
FishingMort_FratioSPR20_mean[i] = lapply(FishingMort_FratioSPR20_mean[i],as.character)
FishingMort_FratioSPR20_sd = aggregate.data.frame(FishingMort_FratioSPR20$FratioSPR20,by=list(FishingMort_FratioSPR20$Genus,FishingMort_FratioSPR20$Species,FishingMort_FratioSPR20$WPP),FUN=sd,na.rm=TRUE)
names(FishingMort_FratioSPR20_sd) = c("Genus","Species","WPP","FratioSPR20_sd")
i = sapply(FishingMort_FratioSPR20_sd,is.factor)
FishingMort_FratioSPR20_sd[i] = lapply(FishingMort_FratioSPR20_sd[i],as.character)
FishingMort_FratioSPR20 = merge(FishingMort_FratioSPR20_mean,FishingMort_FratioSPR20_sd,all.x=TRUE,by=c("Genus","Species","WPP"))
#mean and sd for F@spr30
FishingMort_FratioSPR30 = FishingMort[FishingMort$FratioSPR30_flag==0,]
FishingMort_FratioSPR30_mean = aggregate.data.frame(FishingMort_FratioSPR30$FratioSPR30,by=list(FishingMort_FratioSPR30$Genus,FishingMort_FratioSPR30$Species,FishingMort_FratioSPR30$WPP),FUN=mean,na.rm=TRUE)
names(FishingMort_FratioSPR30_mean) = c("Genus","Species","WPP","FratioSPR30")
i = sapply(FishingMort_FratioSPR30_mean,is.factor)
FishingMort_FratioSPR30_mean[i] = lapply(FishingMort_FratioSPR30_mean[i],as.character)
FishingMort_FratioSPR30_sd = aggregate.data.frame(FishingMort_FratioSPR30$FratioSPR30,by=list(FishingMort_FratioSPR30$Genus,FishingMort_FratioSPR30$Species,FishingMort_FratioSPR30$WPP),FUN=sd,na.rm=TRUE)
names(FishingMort_FratioSPR30_sd) = c("Genus","Species","WPP","FratioSPR30_sd")
i = sapply(FishingMort_FratioSPR30_sd,is.factor)
FishingMort_FratioSPR30_sd[i] = lapply(FishingMort_FratioSPR30_sd[i],as.character)
FishingMort_FratioSPR30 = merge(FishingMort_FratioSPR30_mean,FishingMort_FratioSPR30_sd,all.x=TRUE,by=c("Genus","Species","WPP"))
#mean and sd for F@spr40
FishingMort_FratioSPR40 = FishingMort[FishingMort$FratioSPR40_flag==0,]
FishingMort_FratioSPR40_mean = aggregate.data.frame(FishingMort_FratioSPR40$FratioSPR40,by=list(FishingMort_FratioSPR40$Genus,FishingMort_FratioSPR40$Species,FishingMort_FratioSPR40$WPP),FUN=mean,na.rm=TRUE)
names(FishingMort_FratioSPR40_mean) = c("Genus","Species","WPP","FratioSPR40")
i = sapply(FishingMort_FratioSPR40_mean,is.factor)
FishingMort_FratioSPR40_mean[i] = lapply(FishingMort_FratioSPR40_mean[i],as.character)
FishingMort_FratioSPR40_sd = aggregate.data.frame(FishingMort_FratioSPR40$FratioSPR40,by=list(FishingMort_FratioSPR40$Genus,FishingMort_FratioSPR40$Species,FishingMort_FratioSPR40$WPP),FUN=sd,na.rm=TRUE)
names(FishingMort_FratioSPR40_sd) = c("Genus","Species","WPP","FratioSPR40_sd")
i = sapply(FishingMort_FratioSPR40_sd,is.factor)
FishingMort_FratioSPR40_sd[i] = lapply(FishingMort_FratioSPR40_sd[i],as.character)
FishingMort_FratioSPR40 = merge(FishingMort_FratioSPR40_mean,FishingMort_FratioSPR40_sd,all.x=TRUE,by=c("Genus","Species","WPP"))
FishingMort_Fratio = merge(FishingMort_FratioSPR20,FishingMort_FratioSPR30,all.x=TRUE,by=c("Genus","Species","WPP"))
FishingMort_Fratio = merge(FishingMort_Fratio,FishingMort_FratioSPR40,all.x=TRUE,by=c("Genus","Species","WPP"))
FishingMort_Fratio$FratioSPR20_Lower = FishingMort_Fratio$FratioSPR20 - FishingMort_Fratio$FratioSPR20_sd
FishingMort_Fratio$FratioSPR20_Lower[FishingMort_Fratio$FratioSPR20_Lower<0]=0
FishingMort_Fratio$FratioSPR20_Upper = FishingMort_Fratio$FratioSPR20 + FishingMort_Fratio$FratioSPR20_sd
FishingMort_Fratio$FratioSPR30_Lower = FishingMort_Fratio$FratioSPR30 - FishingMort_Fratio$FratioSPR30_sd
FishingMort_Fratio$FratioSPR30_Lower[FishingMort_Fratio$FratioSPR30_Lower<0]=0
FishingMort_Fratio$FratioSPR30_Upper = FishingMort_Fratio$FratioSPR30 + FishingMort_Fratio$FratioSPR30_sd
FishingMort_Fratio$FratioSPR40_Lower = FishingMort_Fratio$FratioSPR40 - FishingMort_Fratio$FratioSPR40_sd
FishingMort_Fratio$FratioSPR40_Lower[FishingMort_Fratio$FratioSPR40_Lower<0]=0
FishingMort_Fratio$FratioSPR40_Upper = FishingMort_Fratio$FratioSPR40 + FishingMort_Fratio$FratioSPR40_sd

#Merge biomass and fishing mortality ratios and determine colors for the figure
Ratios = merge(CatchBiomass,FishingMort_Fratio,all.x=TRUE,by=c("Genus","Species","WPP"))
Ratios$Color_BratioSPR20=""
Ratios$Color_BratioSPR20[Ratios$BratioSPR20_Lower >= 1 & Ratios$BratioSPR20_Upper >= 1]="green3"
Ratios$Color_BratioSPR20[Ratios$BratioSPR20_Lower < 1 & Ratios$BratioSPR20_Upper > 1]="yellow2"
Ratios$Color_BratioSPR20[Ratios$BratioSPR20_Lower < 1 & Ratios$BratioSPR20_Upper < 1]="red3"
Ratios$Color_BratioSPR30=""
Ratios$Color_BratioSPR30[Ratios$BratioSPR30_Lower >= 1 & Ratios$BratioSPR30_Upper >= 1]="green3"
Ratios$Color_BratioSPR30[Ratios$BratioSPR30_Lower < 1 & Ratios$BratioSPR30_Upper > 1]="yellow2"
Ratios$Color_BratioSPR30[Ratios$BratioSPR30_Lower < 1 & Ratios$BratioSPR30_Upper < 1]="red3"
Ratios$Color_BratioSPR40=""
Ratios$Color_BratioSPR40[Ratios$BratioSPR40_Lower >= 1 & Ratios$BratioSPR40_Upper >= 1]="green3"
Ratios$Color_BratioSPR40[Ratios$BratioSPR40_Lower < 1 & Ratios$BratioSPR40_Upper > 1]="yellow2"
Ratios$Color_BratioSPR40[Ratios$BratioSPR40_Lower < 1 & Ratios$BratioSPR40_Upper < 1]="red3"
Ratios$Color_FratioSPR20=""
Ratios$Color_FratioSPR20[Ratios$FratioSPR20_Lower >= 1 & Ratios$FratioSPR20_Upper >= 1]="green3"
Ratios$Color_FratioSPR20[Ratios$FratioSPR20_Lower < 1 & Ratios$FratioSPR20_Upper > 1]="yellow2"
Ratios$Color_FratioSPR20[Ratios$FratioSPR20_Lower < 1 & Ratios$FratioSPR20_Upper < 1]="red3"
Ratios$Color_FratioSPR30=""
Ratios$Color_FratioSPR30[Ratios$FratioSPR30_Lower >= 1 & Ratios$FratioSPR30_Upper >= 1]="green3"
Ratios$Color_FratioSPR30[Ratios$FratioSPR30_Lower < 1 & Ratios$FratioSPR30_Upper > 1]="yellow2"
Ratios$Color_FratioSPR30[Ratios$FratioSPR30_Lower < 1 & Ratios$FratioSPR30_Upper < 1]="red3"
Ratios$Color_FratioSPR40=""
Ratios$Color_FratioSPR40[Ratios$FratioSPR40_Lower >= 1 & Ratios$FratioSPR40_Upper >= 1]="green3"
Ratios$Color_FratioSPR40[Ratios$FratioSPR40_Lower < 1 & Ratios$FratioSPR40_Upper > 1]="yellow2"
Ratios$Color_FratioSPR40[Ratios$FratioSPR40_Lower < 1 & Ratios$FratioSPR40_Upper < 1]="red3"
#Range Tables: SPR 20 Overfished
SPR20_range_overfished = subset(Ratios,select=c(Genus,Species,WPP,BratioSPR20_Upper,BratioSPR20_Lower))
SPR20_range_overfished$Range = paste(round(SPR20_range_overfished$BratioSPR20_Lower,2)," - ",round(SPR20_range_overfished$BratioSPR20_Upper,2),sep="")
SPR20_range_overfished$GenusSpecies = paste(SPR20_range_overfished$Genus,"_",SPR20_range_overfished$Species,sep="")
SPR20_range_overfished_reshape = subset(SPR20_range_overfished,select=c(GenusSpecies,WPP,Range))
SPR20_range_overfished_reshape = reshape(SPR20_range_overfished_reshape,v.names="Range",idvar="GenusSpecies",timevar="WPP",direction="wide")
GenusSpecies = transpose(data.frame(strsplit(SPR20_range_overfished_reshape$GenusSpecies,split="_")))
names(GenusSpecies) = c("Genus","Species")
SPR20_range_overfished_reshape = cbind(GenusSpecies,SPR20_range_overfished_reshape)
SPR20_range_overfished_reshape$GenusSpecies=NULL
WPP_EEZ = c("Genus","Species",as.character(sort(WPP)),"EEZ")
names(SPR20_range_overfished_reshape) = gsub("Range.","",names(SPR20_range_overfished_reshape))
SPR20_range_overfished_reshape = SPR20_range_overfished_reshape[,WPP_EEZ]
write.table(SPR20_range_overfished_reshape,paste(PATH_output,"SPR20_range_overfishedTable_FINAL.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
SPR20_lower_overfished_reshape = subset(SPR20_range_overfished,select=c(GenusSpecies,WPP,BratioSPR20_Lower))
SPR20_lower_overfished_reshape = reshape(SPR20_lower_overfished_reshape,v.names="BratioSPR20_Lower",idvar="GenusSpecies",timevar="WPP",direction="wide")
GenusSpecies = transpose(data.frame(strsplit(SPR20_lower_overfished_reshape$GenusSpecies,split="_")))
names(GenusSpecies) = c("Genus","Species")
SPR20_lower_overfished_reshape = cbind(GenusSpecies,SPR20_lower_overfished_reshape)
SPR20_lower_overfished_reshape$GenusSpecies=NULL
WPP_EEZ = c("Genus","Species",as.character(sort(WPP)),"EEZ")
names(SPR20_lower_overfished_reshape) = gsub("BratioSPR20_Lower.","",names(SPR20_lower_overfished_reshape))
SPR20_lower_overfished_reshape = SPR20_lower_overfished_reshape[,WPP_EEZ]
write.table(SPR20_lower_overfished_reshape,paste(PATH_output,"SPR20_lower_overfishedTable_FINAL.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
SPR20_upper_overfished_reshape = subset(SPR20_range_overfished,select=c(GenusSpecies,WPP,BratioSPR20_Upper))
SPR20_upper_overfished_reshape = reshape(SPR20_upper_overfished_reshape,v.names="BratioSPR20_Upper",idvar="GenusSpecies",timevar="WPP",direction="wide")
GenusSpecies = transpose(data.frame(strsplit(SPR20_upper_overfished_reshape$GenusSpecies,split="_")))
names(GenusSpecies) = c("Genus","Species")
SPR20_upper_overfished_reshape = cbind(GenusSpecies,SPR20_upper_overfished_reshape)
SPR20_upper_overfished_reshape$GenusSpecies=NULL
WPP_EEZ = c("Genus","Species",as.character(sort(WPP)),"EEZ")
names(SPR20_upper_overfished_reshape) = gsub("BratioSPR20_Upper.","",names(SPR20_upper_overfished_reshape))
SPR20_upper_overfished_reshape = SPR20_upper_overfished_reshape[,WPP_EEZ]
write.table(SPR20_upper_overfished_reshape,paste(PATH_output,"SPR20_upper_overfishedTable_FINAL.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
#Range Tables: SPR 30 Overfished
SPR30_range_overfished = subset(Ratios,select=c(Genus,Species,WPP,BratioSPR30_Upper,BratioSPR30_Lower))
SPR30_range_overfished$Range = paste(round(SPR30_range_overfished$BratioSPR30_Lower,2)," - ",round(SPR30_range_overfished$BratioSPR30_Upper,2),sep="")
SPR30_range_overfished$GenusSpecies = paste(SPR30_range_overfished$Genus,"_",SPR30_range_overfished$Species,sep="")
SPR30_range_overfished_reshape = subset(SPR30_range_overfished,select=c(GenusSpecies,WPP,Range))
SPR30_range_overfished_reshape = reshape(SPR30_range_overfished_reshape,v.names="Range",idvar="GenusSpecies",timevar="WPP",direction="wide")
GenusSpecies = transpose(data.frame(strsplit(SPR30_range_overfished_reshape$GenusSpecies,split="_")))
names(GenusSpecies) = c("Genus","Species")
SPR30_range_overfished_reshape = cbind(GenusSpecies,SPR30_range_overfished_reshape)
SPR30_range_overfished_reshape$GenusSpecies=NULL
WPP_EEZ = c("Genus","Species",as.character(sort(WPP)),"EEZ")
names(SPR30_range_overfished_reshape) = gsub("Range.","",names(SPR30_range_overfished_reshape))
SPR30_range_overfished_reshape = SPR30_range_overfished_reshape[,WPP_EEZ]
write.table(SPR30_range_overfished_reshape,paste(PATH_output,"SPR30_range_overfishedTable_FINAL.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
SPR30_lower_overfished_reshape = subset(SPR30_range_overfished,select=c(GenusSpecies,WPP,BratioSPR30_Lower))
SPR30_lower_overfished_reshape = reshape(SPR30_lower_overfished_reshape,v.names="BratioSPR30_Lower",idvar="GenusSpecies",timevar="WPP",direction="wide")
GenusSpecies = transpose(data.frame(strsplit(SPR30_lower_overfished_reshape$GenusSpecies,split="_")))
names(GenusSpecies) = c("Genus","Species")
SPR30_lower_overfished_reshape = cbind(GenusSpecies,SPR30_lower_overfished_reshape)
SPR30_lower_overfished_reshape$GenusSpecies=NULL
WPP_EEZ = c("Genus","Species",as.character(sort(WPP)),"EEZ")
names(SPR30_lower_overfished_reshape) = gsub("BratioSPR30_Lower.","",names(SPR30_lower_overfished_reshape))
SPR30_lower_overfished_reshape = SPR30_lower_overfished_reshape[,WPP_EEZ]
write.table(SPR30_lower_overfished_reshape,paste(PATH_output,"SPR30_lower_overfishedTable_FINAL.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
SPR30_upper_overfished_reshape = subset(SPR30_range_overfished,select=c(GenusSpecies,WPP,BratioSPR30_Upper))
SPR30_upper_overfished_reshape = reshape(SPR30_upper_overfished_reshape,v.names="BratioSPR30_Upper",idvar="GenusSpecies",timevar="WPP",direction="wide")
GenusSpecies = transpose(data.frame(strsplit(SPR30_upper_overfished_reshape$GenusSpecies,split="_")))
names(GenusSpecies) = c("Genus","Species")
SPR30_upper_overfished_reshape = cbind(GenusSpecies,SPR30_upper_overfished_reshape)
SPR30_upper_overfished_reshape$GenusSpecies=NULL
WPP_EEZ = c("Genus","Species",as.character(sort(WPP)),"EEZ")
names(SPR30_upper_overfished_reshape) = gsub("BratioSPR30_Upper.","",names(SPR30_upper_overfished_reshape))
SPR30_upper_overfished_reshape = SPR30_upper_overfished_reshape[,WPP_EEZ]
write.table(SPR30_upper_overfished_reshape,paste(PATH_output,"SPR30_upper_overfishedTable_FINAL.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
#Range Tables: SPR 40 Overfished
SPR40_range_overfished = subset(Ratios,select=c(Genus,Species,WPP,BratioSPR40_Upper,BratioSPR40_Lower))
SPR40_range_overfished$Range = paste(round(SPR40_range_overfished$BratioSPR40_Lower,2)," - ",round(SPR40_range_overfished$BratioSPR40_Upper,2),sep="")
SPR40_range_overfished$GenusSpecies = paste(SPR40_range_overfished$Genus,"_",SPR40_range_overfished$Species,sep="")
SPR40_range_overfished_reshape = subset(SPR40_range_overfished,select=c(GenusSpecies,WPP,Range))
SPR40_range_overfished_reshape = reshape(SPR40_range_overfished_reshape,v.names="Range",idvar="GenusSpecies",timevar="WPP",direction="wide")
GenusSpecies = transpose(data.frame(strsplit(SPR40_range_overfished_reshape$GenusSpecies,split="_")))
names(GenusSpecies) = c("Genus","Species")
SPR40_range_overfished_reshape = cbind(GenusSpecies,SPR40_range_overfished_reshape)
SPR40_range_overfished_reshape$GenusSpecies=NULL
WPP_EEZ = c("Genus","Species",as.character(sort(WPP)),"EEZ")
names(SPR40_range_overfished_reshape) = gsub("Range.","",names(SPR40_range_overfished_reshape))
SPR40_range_overfished_reshape = SPR40_range_overfished_reshape[,WPP_EEZ]
write.table(SPR40_range_overfished_reshape,paste(PATH_output,"SPR40_range_overfishedTable_FINAL.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
SPR40_lower_overfished_reshape = subset(SPR40_range_overfished,select=c(GenusSpecies,WPP,BratioSPR40_Lower))
SPR40_lower_overfished_reshape = reshape(SPR40_lower_overfished_reshape,v.names="BratioSPR40_Lower",idvar="GenusSpecies",timevar="WPP",direction="wide")
GenusSpecies = transpose(data.frame(strsplit(SPR40_lower_overfished_reshape$GenusSpecies,split="_")))
names(GenusSpecies) = c("Genus","Species")
SPR40_lower_overfished_reshape = cbind(GenusSpecies,SPR40_lower_overfished_reshape)
SPR40_lower_overfished_reshape$GenusSpecies=NULL
WPP_EEZ = c("Genus","Species",as.character(sort(WPP)),"EEZ")
names(SPR40_lower_overfished_reshape) = gsub("BratioSPR40_Lower.","",names(SPR40_lower_overfished_reshape))
SPR40_lower_overfished_reshape = SPR40_lower_overfished_reshape[,WPP_EEZ]
write.table(SPR40_lower_overfished_reshape,paste(PATH_output,"SPR40_lower_overfishedTable_FINAL.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
SPR40_upper_overfished_reshape = subset(SPR40_range_overfished,select=c(GenusSpecies,WPP,BratioSPR40_Upper))
SPR40_upper_overfished_reshape = reshape(SPR40_upper_overfished_reshape,v.names="BratioSPR40_Upper",idvar="GenusSpecies",timevar="WPP",direction="wide")
GenusSpecies = transpose(data.frame(strsplit(SPR40_upper_overfished_reshape$GenusSpecies,split="_")))
names(GenusSpecies) = c("Genus","Species")
SPR40_upper_overfished_reshape = cbind(GenusSpecies,SPR40_upper_overfished_reshape)
SPR40_upper_overfished_reshape$GenusSpecies=NULL
WPP_EEZ = c("Genus","Species",as.character(sort(WPP)),"EEZ")
names(SPR40_upper_overfished_reshape) = gsub("BratioSPR40_Upper.","",names(SPR40_upper_overfished_reshape))
SPR40_upper_overfished_reshape = SPR40_upper_overfished_reshape[,WPP_EEZ]
write.table(SPR40_upper_overfished_reshape,paste(PATH_output,"SPR40_upper_overfishedTable_FINAL.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
#Range Tables: SPR 20 Overfishing
SPR20_range_overfishing = subset(Ratios,select=c(Genus,Species,WPP,FratioSPR20_Upper,FratioSPR20_Lower))
SPR20_range_overfishing$Range = paste(round(SPR20_range_overfishing$FratioSPR20_Lower,2)," - ",round(SPR20_range_overfishing$FratioSPR20_Upper,2),sep="")
SPR20_range_overfishing$GenusSpecies = paste(SPR20_range_overfishing$Genus,"_",SPR20_range_overfishing$Species,sep="")
SPR20_range_overfishing_reshape = subset(SPR20_range_overfishing,select=c(GenusSpecies,WPP,Range))
SPR20_range_overfishing_reshape = reshape(SPR20_range_overfishing_reshape,v.names="Range",idvar="GenusSpecies",timevar="WPP",direction="wide")
GenusSpecies = transpose(data.frame(strsplit(SPR20_range_overfishing_reshape$GenusSpecies,split="_")))
names(GenusSpecies) = c("Genus","Species")
SPR20_range_overfishing_reshape = cbind(GenusSpecies,SPR20_range_overfishing_reshape)
SPR20_range_overfishing_reshape$GenusSpecies=NULL
WPP_EEZ = c("Genus","Species",as.character(sort(WPP)),"EEZ")
names(SPR20_range_overfishing_reshape) = gsub("Range.","",names(SPR20_range_overfishing_reshape))
SPR20_range_overfishing_reshape = SPR20_range_overfishing_reshape[,WPP_EEZ]
write.table(SPR20_range_overfishing_reshape,paste(PATH_output,"SPR20_range_overfishingTable_FINAL.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
SPR20_lower_overfishing_reshape = subset(SPR20_range_overfishing,select=c(GenusSpecies,WPP,FratioSPR20_Lower))
SPR20_lower_overfishing_reshape = reshape(SPR20_lower_overfishing_reshape,v.names="FratioSPR20_Lower",idvar="GenusSpecies",timevar="WPP",direction="wide")
GenusSpecies = transpose(data.frame(strsplit(SPR20_lower_overfishing_reshape$GenusSpecies,split="_")))
names(GenusSpecies) = c("Genus","Species")
SPR20_lower_overfishing_reshape = cbind(GenusSpecies,SPR20_lower_overfishing_reshape)
SPR20_lower_overfishing_reshape$GenusSpecies=NULL
WPP_EEZ = c("Genus","Species",as.character(sort(WPP)),"EEZ")
names(SPR20_lower_overfishing_reshape) = gsub("FratioSPR20_Lower.","",names(SPR20_lower_overfishing_reshape))
SPR20_lower_overfishing_reshape = SPR20_lower_overfishing_reshape[,WPP_EEZ]
write.table(SPR20_lower_overfishing_reshape,paste(PATH_output,"SPR20_lower_overfishingTable_FINAL.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
SPR20_upper_overfishing_reshape = subset(SPR20_range_overfishing,select=c(GenusSpecies,WPP,FratioSPR20_Upper))
SPR20_upper_overfishing_reshape = reshape(SPR20_upper_overfishing_reshape,v.names="FratioSPR20_Upper",idvar="GenusSpecies",timevar="WPP",direction="wide")
GenusSpecies = transpose(data.frame(strsplit(SPR20_upper_overfishing_reshape$GenusSpecies,split="_")))
names(GenusSpecies) = c("Genus","Species")
SPR20_upper_overfishing_reshape = cbind(GenusSpecies,SPR20_upper_overfishing_reshape)
SPR20_upper_overfishing_reshape$GenusSpecies=NULL
WPP_EEZ = c("Genus","Species",as.character(sort(WPP)),"EEZ")
names(SPR20_upper_overfishing_reshape) = gsub("FratioSPR20_Upper.","",names(SPR20_upper_overfishing_reshape))
SPR20_upper_overfishing_reshape = SPR20_upper_overfishing_reshape[,WPP_EEZ]
write.table(SPR20_upper_overfishing_reshape,paste(PATH_output,"SPR20_upper_overfishingTable_FINAL.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
#Range Tables: SPR 30 Overfishing
SPR30_range_overfishing = subset(Ratios,select=c(Genus,Species,WPP,FratioSPR30_Upper,FratioSPR30_Lower))
SPR30_range_overfishing$Range = paste(round(SPR30_range_overfishing$FratioSPR30_Lower,2)," - ",round(SPR30_range_overfishing$FratioSPR30_Upper,2),sep="")
SPR30_range_overfishing$GenusSpecies = paste(SPR30_range_overfishing$Genus,"_",SPR30_range_overfishing$Species,sep="")
SPR30_range_overfishing_reshape = subset(SPR30_range_overfishing,select=c(GenusSpecies,WPP,Range))
SPR30_range_overfishing_reshape = reshape(SPR30_range_overfishing_reshape,v.names="Range",idvar="GenusSpecies",timevar="WPP",direction="wide")
GenusSpecies = transpose(data.frame(strsplit(SPR30_range_overfishing_reshape$GenusSpecies,split="_")))
names(GenusSpecies) = c("Genus","Species")
SPR30_range_overfishing_reshape = cbind(GenusSpecies,SPR30_range_overfishing_reshape)
SPR30_range_overfishing_reshape$GenusSpecies=NULL
WPP_EEZ = c("Genus","Species",as.character(sort(WPP)),"EEZ")
names(SPR30_range_overfishing_reshape) = gsub("Range.","",names(SPR30_range_overfishing_reshape))
SPR30_range_overfishing_reshape = SPR30_range_overfishing_reshape[,WPP_EEZ]
write.table(SPR30_range_overfishing_reshape,paste(PATH_output,"SPR30_range_overfishingTable_FINAL.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
SPR30_lower_overfishing_reshape = subset(SPR30_range_overfishing,select=c(GenusSpecies,WPP,FratioSPR30_Lower))
SPR30_lower_overfishing_reshape = reshape(SPR30_lower_overfishing_reshape,v.names="FratioSPR30_Lower",idvar="GenusSpecies",timevar="WPP",direction="wide")
GenusSpecies = transpose(data.frame(strsplit(SPR30_lower_overfishing_reshape$GenusSpecies,split="_")))
names(GenusSpecies) = c("Genus","Species")
SPR30_lower_overfishing_reshape = cbind(GenusSpecies,SPR30_lower_overfishing_reshape)
SPR30_lower_overfishing_reshape$GenusSpecies=NULL
WPP_EEZ = c("Genus","Species",as.character(sort(WPP)),"EEZ")
names(SPR30_lower_overfishing_reshape) = gsub("FratioSPR30_Lower.","",names(SPR30_lower_overfishing_reshape))
SPR30_lower_overfishing_reshape = SPR30_lower_overfishing_reshape[,WPP_EEZ]
write.table(SPR30_lower_overfishing_reshape,paste(PATH_output,"SPR30_lower_overfishingTable_FINAL.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
SPR30_upper_overfishing_reshape = subset(SPR30_range_overfishing,select=c(GenusSpecies,WPP,FratioSPR30_Upper))
SPR30_upper_overfishing_reshape = reshape(SPR30_upper_overfishing_reshape,v.names="FratioSPR30_Upper",idvar="GenusSpecies",timevar="WPP",direction="wide")
GenusSpecies = transpose(data.frame(strsplit(SPR30_upper_overfishing_reshape$GenusSpecies,split="_")))
names(GenusSpecies) = c("Genus","Species")
SPR30_upper_overfishing_reshape = cbind(GenusSpecies,SPR30_upper_overfishing_reshape)
SPR30_upper_overfishing_reshape$GenusSpecies=NULL
WPP_EEZ = c("Genus","Species",as.character(sort(WPP)),"EEZ")
names(SPR30_upper_overfishing_reshape) = gsub("FratioSPR30_Upper.","",names(SPR30_upper_overfishing_reshape))
SPR30_upper_overfishing_reshape = SPR30_upper_overfishing_reshape[,WPP_EEZ]
write.table(SPR30_upper_overfishing_reshape,paste(PATH_output,"SPR30_upper_overfishingTable_FINAL.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
#Range Tables: SPR 40 Overfishing
SPR40_range_overfishing = subset(Ratios,select=c(Genus,Species,WPP,FratioSPR40_Upper,FratioSPR40_Lower))
SPR40_range_overfishing$Range = paste(round(SPR40_range_overfishing$FratioSPR40_Lower,2)," - ",round(SPR40_range_overfishing$FratioSPR40_Upper,2),sep="")
SPR40_range_overfishing$GenusSpecies = paste(SPR40_range_overfishing$Genus,"_",SPR40_range_overfishing$Species,sep="")
SPR40_range_overfishing_reshape = subset(SPR40_range_overfishing,select=c(GenusSpecies,WPP,Range))
SPR40_range_overfishing_reshape = reshape(SPR40_range_overfishing_reshape,v.names="Range",idvar="GenusSpecies",timevar="WPP",direction="wide")
GenusSpecies = transpose(data.frame(strsplit(SPR40_range_overfishing_reshape$GenusSpecies,split="_")))
names(GenusSpecies) = c("Genus","Species")
SPR40_range_overfishing_reshape = cbind(GenusSpecies,SPR40_range_overfishing_reshape)
SPR40_range_overfishing_reshape$GenusSpecies=NULL
WPP_EEZ = c("Genus","Species",as.character(sort(WPP)),"EEZ")
names(SPR40_range_overfishing_reshape) = gsub("Range.","",names(SPR40_range_overfishing_reshape))
SPR40_range_overfishing_reshape = SPR40_range_overfishing_reshape[,WPP_EEZ]
write.table(SPR40_range_overfishing_reshape,paste(PATH_output,"SPR40_range_overfishingTable_FINAL.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
SPR40_lower_overfishing_reshape = subset(SPR40_range_overfishing,select=c(GenusSpecies,WPP,FratioSPR40_Lower))
SPR40_lower_overfishing_reshape = reshape(SPR40_lower_overfishing_reshape,v.names="FratioSPR40_Lower",idvar="GenusSpecies",timevar="WPP",direction="wide")
GenusSpecies = transpose(data.frame(strsplit(SPR40_lower_overfishing_reshape$GenusSpecies,split="_")))
names(GenusSpecies) = c("Genus","Species")
SPR40_lower_overfishing_reshape = cbind(GenusSpecies,SPR40_lower_overfishing_reshape)
SPR40_lower_overfishing_reshape$GenusSpecies=NULL
WPP_EEZ = c("Genus","Species",as.character(sort(WPP)),"EEZ")
names(SPR40_lower_overfishing_reshape) = gsub("FratioSPR40_Lower.","",names(SPR40_lower_overfishing_reshape))
SPR40_lower_overfishing_reshape = SPR40_lower_overfishing_reshape[,WPP_EEZ]
write.table(SPR40_lower_overfishing_reshape,paste(PATH_output,"SPR40_lower_overfishingTable_FINAL.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
SPR40_upper_overfishing_reshape = subset(SPR40_range_overfishing,select=c(GenusSpecies,WPP,FratioSPR40_Upper))
SPR40_upper_overfishing_reshape = reshape(SPR40_upper_overfishing_reshape,v.names="FratioSPR40_Upper",idvar="GenusSpecies",timevar="WPP",direction="wide")
GenusSpecies = transpose(data.frame(strsplit(SPR40_upper_overfishing_reshape$GenusSpecies,split="_")))
names(GenusSpecies) = c("Genus","Species")
SPR40_upper_overfishing_reshape = cbind(GenusSpecies,SPR40_upper_overfishing_reshape)
SPR40_upper_overfishing_reshape$GenusSpecies=NULL
WPP_EEZ = c("Genus","Species",as.character(sort(WPP)),"EEZ")
names(SPR40_upper_overfishing_reshape) = gsub("FratioSPR40_Upper.","",names(SPR40_upper_overfishing_reshape))
SPR40_upper_overfishing_reshape = SPR40_upper_overfishing_reshape[,WPP_EEZ]
write.table(SPR40_upper_overfishing_reshape,paste(PATH_output,"SPR40_upper_overfishingTable_FINAL.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)
#Establish Functions to draw custom circles half filled on the plot
upper.half.circle = function(x,y,r,nsteps=1000,...)
{
    rs <- seq(0,pi,len=nsteps)
    xc <- x+r*cos(rs)
    yc <- y+r*sin(rs)
    polygon(xc,yc,...)
}
lower.half.circle = function(x,y,r,nsteps=1000,...)
{
    rs <- seq(0,pi,len=nsteps)
    xc <- x-r*cos(rs)
    yc <- y-r*sin(rs)
    polygon(xc,yc,...)
}
#Foregone Catch and Status Plot for WPP 571
thisWPP = "571"
Sub_571 = Ratios[Ratios$WPP==thisWPP,]
Sub_571 = Sub_571[order(Sub_571$CatchRatioSPR40),]
GenSpp_thisWPP = unique(subset(Sub_571,select=c(Genus,Species)))
GenSpp_thisWPP$GenusSpecies = paste(GenSpp_thisWPP$Genus,GenSpp_thisWPP$Species,sep=" ")
png(paste(PATH_highLevelSummaryPlots,paste(paste("CatchRatios_OverfishingColors",thisWPP,sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)    #was 3200 and 3200
par(mar=c(14, 5, 4, 2) + 0.1,xpd=NA)
plot(1:dim(Sub_571)[1],Sub_571$CatchRatioSPR40,xaxt="n",yaxt="n",main="",xlab ="",ylab="Yield Relative to Y@SPR40",type="n",ylim=c(0,ceiling(max(Sub_571$CatchRatioSPR40_Upper)))) #pch=21,cex=1.5,col="black",) #,bg=Sub_571$Color_FratioSPR40)
axis(1,at=1:dim(Sub_571)[1],labels=GenSpp_thisWPP$GenusSpecies,las=2,cex.axis=0.8,tick=TRUE)
axis(2,at=0:ceiling(max(Sub_571$CatchRatioSPR40_Upper)),labels=seq(0,ceiling(max(Sub_571$CatchRatioSPR40_Upper)),by=1),las=2,cex.axis=0.8,tick=TRUE)
for(j in 1:dim(Sub_571)[1])
{
  upper.half.circle(j,Sub_571$CatchRatioSPR40[j],r=0.1,nsteps=10000,col=Sub_571$Color_FratioSPR40[j])   #for the larger areas like 712, the r was equal to 0.3
  lower.half.circle(j,Sub_571$CatchRatioSPR40[j],r=0.1,nsteps=10000,col=Sub_571$Color_BratioSPR20[j])
  points(rep(j,2),c(Sub_571$CatchRatioSPR40_Lower[j],Sub_571$CatchRatioSPR40_Upper[j]),col="gray",type="l",lwd=2)
}
xLineVal = seq((1-0.05),(dim(Sub_571)[1]+0.05),by=0.05)
yLineVal = rep(1,length(xLineVal))
points(xLineVal,yLineVal,type="l",lty=2,col="gray",lwd=2)
title("Yield Relative to Y@SPR40")
mtext(paste("FMA ",thisWPP,sep=""))
upper.half.circle(1.2,-(ceiling(max(Sub_571$CatchRatioSPR40_Upper))+1.0),r=0.1,nsteps=10000,col="green3")
text(1.3,-(ceiling(max(Sub_571$CatchRatioSPR40_Upper))+1.0),"No Overfishing",cex=0.8,pos=4)
upper.half.circle(1.2,-(ceiling(max(Sub_571$CatchRatioSPR40_Upper))+1.5),r=0.1,nsteps=10000,col="yellow2")
text(1.3,-(ceiling(max(Sub_571$CatchRatioSPR40_Upper))+1.5),"Possible Overfishing",cex=0.8,pos=4)
upper.half.circle(1.2,-(ceiling(max(Sub_571$CatchRatioSPR40_Upper))+2.0),r=0.1,nsteps=10000,col="red3")
text(1.3,-(ceiling(max(Sub_571$CatchRatioSPR40_Upper))+2.0),"Overfishing",cex=0.8,pos=4)
lower.half.circle(2.7,-(ceiling(max(Sub_571$CatchRatioSPR40_Upper))+1.0),r=0.1,nsteps=10000,col="green3")
text(2.8,-(ceiling(max(Sub_571$CatchRatioSPR40_Upper))+1.0),"Not Overfished",cex=0.8,pos=4)
lower.half.circle(2.7,-(ceiling(max(Sub_571$CatchRatioSPR40_Upper))+1.5),r=0.1,nsteps=10000,col="yellow2")
text(2.8,-(ceiling(max(Sub_571$CatchRatioSPR40_Upper))+1.5),"Possibly Overfished",cex=0.8,pos=4)
lower.half.circle(2.7,-(ceiling(max(Sub_571$CatchRatioSPR40_Upper))+2.0),r=0.1,nsteps=10000,col="red3")
text(2.8,-(ceiling(max(Sub_571$CatchRatioSPR40_Upper))+2.0),"Overfished",cex=0.8,pos=4)
dev.off()
#Foregone Catch and Status Plot for WPP 572
thisWPP = "572"
Sub_572 = Ratios[Ratios$WPP==thisWPP,]
Sub_572 = Sub_572[order(Sub_572$CatchRatioSPR40),]
GenSpp_thisWPP = unique(subset(Sub_572,select=c(Genus,Species)))
GenSpp_thisWPP$GenusSpecies = paste(GenSpp_thisWPP$Genus,GenSpp_thisWPP$Species,sep=" ")
png(paste(PATH_highLevelSummaryPlots,paste(paste("CatchRatios_OverfishingColors",thisWPP,sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)    #was 3200 and 3200
par(mar=c(14, 5, 4, 2) + 0.1,xpd=NA)
plot(1:dim(Sub_572)[1],Sub_572$CatchRatioSPR40,xaxt="n",yaxt="n",main="",xlab ="",ylab="Yield Relative to Y@SPR40",type="n",ylim=c(0,ceiling(max(Sub_572$CatchRatioSPR40_Upper)))) #pch=21,cex=1.5,col="black",) #,bg=Sub_572$Color_FratioSPR40)
axis(1,at=1:dim(Sub_572)[1],labels=GenSpp_thisWPP$GenusSpecies,las=2,cex.axis=0.8,tick=TRUE)
axis(2,at=0:ceiling(max(Sub_572$CatchRatioSPR40_Upper)),labels=seq(0,ceiling(max(Sub_572$CatchRatioSPR40_Upper)),by=1),las=2,cex.axis=0.8,tick=TRUE)
for(j in 1:dim(Sub_572)[1])
{
  upper.half.circle(j,Sub_572$CatchRatioSPR40[j],r=0.1,nsteps=10000,col=Sub_572$Color_FratioSPR40[j])   #for the larger areas like 712, the r was equal to 0.3
  lower.half.circle(j,Sub_572$CatchRatioSPR40[j],r=0.1,nsteps=10000,col=Sub_572$Color_BratioSPR20[j])
  points(rep(j,2),c(Sub_572$CatchRatioSPR40_Lower[j],Sub_572$CatchRatioSPR40_Upper[j]),col="gray",type="l",lwd=2)
}
xLineVal = seq((1-0.05),(dim(Sub_572)[1]+0.05),by=0.05)
yLineVal = rep(1,length(xLineVal))
points(xLineVal,yLineVal,type="l",lty=2,col="gray",lwd=2)
title("Yield Relative to Y@SPR40")
mtext(paste("FMA ",thisWPP,sep=""))
upper.half.circle(1.2,-(ceiling(max(Sub_572$CatchRatioSPR40_Upper))+ 0.4),r=0.1,nsteps=10000,col="green3")
text(1.3,-(ceiling(max(Sub_572$CatchRatioSPR40_Upper))+0.4),"No Overfishing",cex=0.8,pos=4)
upper.half.circle(1.2,-(ceiling(max(Sub_572$CatchRatioSPR40_Upper))+0.7),r=0.1,nsteps=10000,col="yellow2")
text(1.3,-(ceiling(max(Sub_572$CatchRatioSPR40_Upper))+0.7),"Possible Overfishing",cex=0.8,pos=4)
upper.half.circle(1.2,-(ceiling(max(Sub_572$CatchRatioSPR40_Upper))+1.0),r=0.1,nsteps=10000,col="red3")
text(1.3,-(ceiling(max(Sub_572$CatchRatioSPR40_Upper))+1.0),"Overfishing",cex=0.8,pos=4)
lower.half.circle(3.3,-(ceiling(max(Sub_572$CatchRatioSPR40_Upper))+0.4),r=0.1,nsteps=10000,col="green3")
text(3.4,-(ceiling(max(Sub_572$CatchRatioSPR40_Upper))+0.4),"Not Overfished",cex=0.8,pos=4)
lower.half.circle(3.3,-(ceiling(max(Sub_572$CatchRatioSPR40_Upper))+0.7),r=0.1,nsteps=10000,col="yellow2")
text(3.4,-(ceiling(max(Sub_572$CatchRatioSPR40_Upper))+0.7),"Possibly Overfished",cex=0.8,pos=4)
lower.half.circle(3.3,-(ceiling(max(Sub_572$CatchRatioSPR40_Upper))+1.0),r=0.1,nsteps=10000,col="red3")
text(3.4,-(ceiling(max(Sub_572$CatchRatioSPR40_Upper))+1.0),"Overfished",cex=0.8,pos=4)
dev.off()
#Foregone Catch and Status Plot for WPP 573
thisWPP = "573"
Sub_573 = Ratios[Ratios$WPP==thisWPP,]
Sub_573 = subset(Sub_573,!is.na(Sub_573$FratioSPR40))
Sub_573 = Sub_573[order(Sub_573$CatchRatioSPR40),]
GenSpp_thisWPP = unique(subset(Sub_573,select=c(Genus,Species)))
GenSpp_thisWPP$GenusSpecies = paste(GenSpp_thisWPP$Genus,GenSpp_thisWPP$Species,sep=" ")
png(paste(PATH_highLevelSummaryPlots,paste(paste("CatchRatios_OverfishingColors",thisWPP,sep="_"),".png",sep=""),sep="/"),units="px",width=5200,height=3200,res=600)    #was 3200 and 3200
par(mar=c(14, 5, 4, 2) + 0.1,xpd=NA)
plot(1:dim(Sub_573)[1],Sub_573$CatchRatioSPR40,xaxt="n",yaxt="n",main="",xlab ="",ylab="Yield Relative to Y@SPR40",type="n",ylim=c(0,ceiling(max(Sub_573$CatchRatioSPR40_Upper)))) #pch=21,cex=1.5,col="black",) #,bg=Sub_573$Color_FratioSPR40)
axis(1,at=1:dim(Sub_573)[1],labels=GenSpp_thisWPP$GenusSpecies,las=2,cex.axis=0.8,tick=TRUE)
axis(2,at=0:ceiling(max(Sub_573$CatchRatioSPR40_Upper)),labels=seq(0,ceiling(max(Sub_573$CatchRatioSPR40_Upper)),by=1),las=2,cex.axis=0.8,tick=TRUE)
for(j in 1:dim(Sub_573)[1])
{
  upper.half.circle(j,Sub_573$CatchRatioSPR40[j],r=0.25,nsteps=10000,col=Sub_573$Color_FratioSPR40[j])   #for the larger areas like 712, the r was equal to 0.3
  lower.half.circle(j,Sub_573$CatchRatioSPR40[j],r=0.25,nsteps=10000,col=Sub_573$Color_BratioSPR20[j])
  points(rep(j,2),c(Sub_573$CatchRatioSPR40_Lower[j],Sub_573$CatchRatioSPR40_Upper[j]),col="gray",type="l",lwd=2)
}
xLineVal = seq((1-0.05),(dim(Sub_573)[1]+0.05),by=0.05)
yLineVal = rep(1,length(xLineVal))
points(xLineVal,yLineVal,type="l",lty=2,col="gray",lwd=2)
title("Yield Relative to Y@SPR40")
mtext(paste("FMA ",thisWPP,sep=""))
upper.half.circle(1.0,-(ceiling(max(Sub_573$CatchRatioSPR40_Upper))+ 1.3),r=0.25,nsteps=10000,col="green3")
text(1.1,-(ceiling(max(Sub_573$CatchRatioSPR40_Upper))+1.3),"No Overfishing",cex=0.8,pos=4)
upper.half.circle(1.0,-(ceiling(max(Sub_573$CatchRatioSPR40_Upper))+1.8),r=0.25,nsteps=10000,col="yellow2")
text(1.1,-(ceiling(max(Sub_573$CatchRatioSPR40_Upper))+1.8),"Possible Overfishing",cex=0.8,pos=4)
upper.half.circle(1.0,-(ceiling(max(Sub_573$CatchRatioSPR40_Upper))+2.3),r=0.25,nsteps=10000,col="red3")
text(1.1,-(ceiling(max(Sub_573$CatchRatioSPR40_Upper))+2.3),"Overfishing",cex=0.8,pos=4)
lower.half.circle(9.1,-(ceiling(max(Sub_573$CatchRatioSPR40_Upper))+1.3),r=0.25,nsteps=10000,col="green3")
text(9.2,-(ceiling(max(Sub_573$CatchRatioSPR40_Upper))+1.3),"Not Overfished",cex=0.8,pos=4)
lower.half.circle(9.1,-(ceiling(max(Sub_573$CatchRatioSPR40_Upper))+1.8),r=0.25,nsteps=10000,col="yellow2")
text(9.2,-(ceiling(max(Sub_573$CatchRatioSPR40_Upper))+1.8),"Possibly Overfished",cex=0.8,pos=4)
lower.half.circle(9.1,-(ceiling(max(Sub_573$CatchRatioSPR40_Upper))+2.3),r=0.25,nsteps=10000,col="red3")
text(9.2,-(ceiling(max(Sub_573$CatchRatioSPR40_Upper))+2.3),"Overfished",cex=0.8,pos=4)
dev.off()
#Foregone Catch and Status Plot for WPP 711
thisWPP = "711"
Sub_711 = Ratios[Ratios$WPP==thisWPP,]
Sub_711 = subset(Sub_711,!is.na(Sub_711$FratioSPR40))
Sub_711 = subset(Sub_711,!is.na(Sub_711$BratioSPR20))
Sub_711 = Sub_711[order(Sub_711$CatchRatioSPR40),]
GenSpp_thisWPP = unique(subset(Sub_711,select=c(Genus,Species)))
GenSpp_thisWPP$GenusSpecies = paste(GenSpp_thisWPP$Genus,GenSpp_thisWPP$Species,sep=" ")
png(paste(PATH_highLevelSummaryPlots,paste(paste("CatchRatios_OverfishingColors",thisWPP,sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)    #was 3200 and 3200
par(mar=c(14, 5, 4, 2) + 0.1,xpd=NA)
plot(1:dim(Sub_711)[1],Sub_711$CatchRatioSPR40,xaxt="n",yaxt="n",main="",xlab ="",ylab="Yield Relative to Y@SPR40",type="n",ylim=c(0,ceiling(max(Sub_711$CatchRatioSPR40_Upper)))) #pch=21,cex=1.5,col="black",) #,bg=Sub_711$Color_FratioSPR40)
axis(1,at=1:dim(Sub_711)[1],labels=GenSpp_thisWPP$GenusSpecies,las=2,cex.axis=0.8,tick=TRUE)
axis(2,at=0:ceiling(max(Sub_711$CatchRatioSPR40_Upper)),labels=seq(0,ceiling(max(Sub_711$CatchRatioSPR40_Upper)),by=1),las=2,cex.axis=0.8,tick=TRUE)
for(j in 1:dim(Sub_711)[1])
{
  upper.half.circle(j,Sub_711$CatchRatioSPR40[j],r=0.15,nsteps=10000,col=Sub_711$Color_FratioSPR40[j])   #for the larger areas like 712, the r was equal to 0.3
  lower.half.circle(j,Sub_711$CatchRatioSPR40[j],r=0.15,nsteps=10000,col=Sub_711$Color_BratioSPR20[j])
  points(rep(j,2),c(Sub_711$CatchRatioSPR40_Lower[j],Sub_711$CatchRatioSPR40_Upper[j]),col="gray",type="l",lwd=2)
}
xLineVal = seq((1-0.05),(dim(Sub_711)[1]+0.05),by=0.05)
yLineVal = rep(1,length(xLineVal))
points(xLineVal,yLineVal,type="l",lty=2,col="gray",lwd=2)
title("Yield Relative to Y@SPR40")
mtext(paste("FMA ",thisWPP,sep=""))
upper.half.circle(0.7,-(ceiling(max(Sub_711$CatchRatioSPR40_Upper))+ 0.5),r=0.15,nsteps=10000,col="green3")
text(0.8,-(ceiling(max(Sub_711$CatchRatioSPR40_Upper))+0.5),"No Overfishing",cex=0.8,pos=4)
upper.half.circle(0.7,-(ceiling(max(Sub_711$CatchRatioSPR40_Upper))+0.8),r=0.15,nsteps=10000,col="yellow2")
text(0.8,-(ceiling(max(Sub_711$CatchRatioSPR40_Upper))+0.8),"Possible Overfishing",cex=0.8,pos=4)
upper.half.circle(0.7,-(ceiling(max(Sub_711$CatchRatioSPR40_Upper))+1.1),r=0.15,nsteps=10000,col="red3")
text(0.8,-(ceiling(max(Sub_711$CatchRatioSPR40_Upper))+1.1),"Overfishing",cex=0.8,pos=4)
lower.half.circle(2.5,-(ceiling(max(Sub_711$CatchRatioSPR40_Upper))+0.5),r=0.15,nsteps=10000,col="green3")
text(2.6,-(ceiling(max(Sub_711$CatchRatioSPR40_Upper))+0.5),"Not Overfished",cex=0.8,pos=4)
lower.half.circle(2.5,-(ceiling(max(Sub_711$CatchRatioSPR40_Upper))+0.8),r=0.15,nsteps=10000,col="yellow2")
text(2.6,-(ceiling(max(Sub_711$CatchRatioSPR40_Upper))+0.8),"Possibly Overfished",cex=0.8,pos=4)
lower.half.circle(2.5,-(ceiling(max(Sub_711$CatchRatioSPR40_Upper))+1.1),r=0.15,nsteps=10000,col="red3")
text(2.6,-(ceiling(max(Sub_711$CatchRatioSPR40_Upper))+1.1),"Overfished",cex=0.8,pos=4)
dev.off()
#Foregone Catch and Status Plot for WPP 712
thisWPP = "712"
Sub_712 = Ratios[Ratios$WPP==thisWPP,]
Sub_712 = subset(Sub_712,!is.na(Sub_712$FratioSPR40))
Sub_712 = Sub_712[order(Sub_712$CatchRatioSPR40),]
GenSpp_thisWPP = unique(subset(Sub_712,select=c(Genus,Species)))
GenSpp_thisWPP$GenusSpecies = paste(GenSpp_thisWPP$Genus,GenSpp_thisWPP$Species,sep=" ")
png(paste(PATH_highLevelSummaryPlots,paste(paste("CatchRatios_OverfishingColors",thisWPP,sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)    #was 3200 and 3200
par(mar=c(14, 5, 4, 2) + 0.1,xpd=NA)
plot(1:dim(Sub_712)[1],Sub_712$CatchRatioSPR40,xaxt="n",yaxt="n",main="",xlab ="",ylab="Yield Relative to Y@SPR40",type="n",ylim=c(0,ceiling(max(Sub_712$CatchRatioSPR40_Upper)))) #pch=21,cex=1.5,col="black",) #,bg=Sub_712$Color_FratioSPR40)
axis(1,at=1:dim(Sub_712)[1],labels=GenSpp_thisWPP$GenusSpecies,las=2,cex.axis=0.8,tick=TRUE)
axis(2,at=0:ceiling(max(Sub_712$CatchRatioSPR40_Upper)),labels=seq(0,ceiling(max(Sub_712$CatchRatioSPR40_Upper)),by=1),las=2,cex.axis=0.8,tick=TRUE)
for(j in 1:dim(Sub_712)[1])
{
  upper.half.circle(j,Sub_712$CatchRatioSPR40[j],r=0.25,nsteps=10000,col=Sub_712$Color_FratioSPR40[j])   #for the larger areas like 712, the r was equal to 0.3
  lower.half.circle(j,Sub_712$CatchRatioSPR40[j],r=0.25,nsteps=10000,col=Sub_712$Color_BratioSPR20[j])
  points(rep(j,2),c(Sub_712$CatchRatioSPR40_Lower[j],Sub_712$CatchRatioSPR40_Upper[j]),col="gray",type="l",lwd=2)
}
xLineVal = seq((1-0.05),(dim(Sub_712)[1]+0.05),by=0.05)
yLineVal = rep(1,length(xLineVal))
points(xLineVal,yLineVal,type="l",lty=2,col="gray",lwd=2)
title("Yield Relative to Y@SPR40")
mtext(paste("FMA ",thisWPP,sep=""))
upper.half.circle(1.0,-(ceiling(max(Sub_712$CatchRatioSPR40_Upper))+ 1.3),r=0.25,nsteps=10000,col="green3")
text(1.1,-(ceiling(max(Sub_712$CatchRatioSPR40_Upper))+1.3),"No Overfishing",cex=0.8,pos=4)
upper.half.circle(1.0,-(ceiling(max(Sub_712$CatchRatioSPR40_Upper))+1.8),r=0.25,nsteps=10000,col="yellow2")
text(1.1,-(ceiling(max(Sub_712$CatchRatioSPR40_Upper))+1.8),"Possible Overfishing",cex=0.8,pos=4)
upper.half.circle(1.0,-(ceiling(max(Sub_712$CatchRatioSPR40_Upper))+2.3),r=0.25,nsteps=10000,col="red3")
text(1.1,-(ceiling(max(Sub_712$CatchRatioSPR40_Upper))+2.3),"Overfishing",cex=0.8,pos=4)
lower.half.circle(8.1,-(ceiling(max(Sub_712$CatchRatioSPR40_Upper))+1.3),r=0.25,nsteps=10000,col="green3")
text(8.2,-(ceiling(max(Sub_712$CatchRatioSPR40_Upper))+1.3),"Not Overfished",cex=0.8,pos=4)
lower.half.circle(8.1,-(ceiling(max(Sub_712$CatchRatioSPR40_Upper))+1.8),r=0.25,nsteps=10000,col="yellow2")
text(8.2,-(ceiling(max(Sub_712$CatchRatioSPR40_Upper))+1.8),"Possibly Overfished",cex=0.8,pos=4)
lower.half.circle(8.1,-(ceiling(max(Sub_712$CatchRatioSPR40_Upper))+2.3),r=0.25,nsteps=10000,col="red3")
text(8.2,-(ceiling(max(Sub_712$CatchRatioSPR40_Upper))+2.3),"Overfished",cex=0.8,pos=4)
dev.off()
#Foregone Catch and Status Plot for WPP 713
thisWPP = "713"
Sub_713 = Ratios[Ratios$WPP==thisWPP,]
Sub_713 = subset(Sub_713,!is.na(Sub_713$FratioSPR40))
Sub_713 = Sub_713[order(Sub_713$CatchRatioSPR40),]
GenSpp_thisWPP = unique(subset(Sub_713,select=c(Genus,Species)))
GenSpp_thisWPP$GenusSpecies = paste(GenSpp_thisWPP$Genus,GenSpp_thisWPP$Species,sep=" ")
png(paste(PATH_highLevelSummaryPlots,paste(paste("CatchRatios_OverfishingColors",thisWPP,sep="_"),".png",sep=""),sep="/"),units="px",width=5200,height=3200,res=600)    #was 3200 and 3200
par(mar=c(14, 5, 4, 2) + 0.1,xpd=NA)
plot(1:dim(Sub_713)[1],Sub_713$CatchRatioSPR40,xaxt="n",yaxt="n",main="",xlab ="",ylab="Yield Relative to Y@SPR40",type="n",ylim=c(0,ceiling(max(Sub_713$CatchRatioSPR40_Upper)))) #pch=21,cex=1.5,col="black",) #,bg=Sub_713$Color_FratioSPR40)
axis(1,at=1:dim(Sub_713)[1],labels=GenSpp_thisWPP$GenusSpecies,las=2,cex.axis=0.8,tick=TRUE)
axis(2,at=0:ceiling(max(Sub_713$CatchRatioSPR40_Upper)),labels=seq(0,ceiling(max(Sub_713$CatchRatioSPR40_Upper)),by=1),las=2,cex.axis=0.8,tick=TRUE)
for(j in 1:dim(Sub_713)[1])
{
  upper.half.circle(j,Sub_713$CatchRatioSPR40[j],r=0.25,nsteps=10000,col=Sub_713$Color_FratioSPR40[j])   #for the larger areas like 713, the r was equal to 0.3
  lower.half.circle(j,Sub_713$CatchRatioSPR40[j],r=0.25,nsteps=10000,col=Sub_713$Color_BratioSPR20[j])
  points(rep(j,2),c(Sub_713$CatchRatioSPR40_Lower[j],Sub_713$CatchRatioSPR40_Upper[j]),col="gray",type="l",lwd=2)
}
xLineVal = seq((1-0.05),(dim(Sub_713)[1]+0.05),by=0.05)
yLineVal = rep(1,length(xLineVal))
points(xLineVal,yLineVal,type="l",lty=2,col="gray",lwd=2)
title("Yield Relative to Y@SPR40")
mtext(paste("FMA ",thisWPP,sep=""))
upper.half.circle(1.0,-(ceiling(max(Sub_713$CatchRatioSPR40_Upper))+ 1.3),r=0.25,nsteps=10000,col="green3")
text(1.1,-(ceiling(max(Sub_713$CatchRatioSPR40_Upper))+1.3),"No Overfishing",cex=0.8,pos=4)
upper.half.circle(1.0,-(ceiling(max(Sub_713$CatchRatioSPR40_Upper))+1.8),r=0.25,nsteps=10000,col="yellow2")
text(1.1,-(ceiling(max(Sub_713$CatchRatioSPR40_Upper))+1.8),"Possible Overfishing",cex=0.8,pos=4)
upper.half.circle(1.0,-(ceiling(max(Sub_713$CatchRatioSPR40_Upper))+2.3),r=0.25,nsteps=10000,col="red3")
text(1.1,-(ceiling(max(Sub_713$CatchRatioSPR40_Upper))+2.3),"Overfishing",cex=0.8,pos=4)
lower.half.circle(9.1,-(ceiling(max(Sub_713$CatchRatioSPR40_Upper))+1.3),r=0.25,nsteps=10000,col="green3")
text(9.2,-(ceiling(max(Sub_713$CatchRatioSPR40_Upper))+1.3),"Not Overfished",cex=0.8,pos=4)
lower.half.circle(9.1,-(ceiling(max(Sub_713$CatchRatioSPR40_Upper))+1.8),r=0.25,nsteps=10000,col="yellow2")
text(9.2,-(ceiling(max(Sub_713$CatchRatioSPR40_Upper))+1.8),"Possibly Overfished",cex=0.8,pos=4)
lower.half.circle(9.1,-(ceiling(max(Sub_713$CatchRatioSPR40_Upper))+2.3),r=0.25,nsteps=10000,col="red3")
text(9.2,-(ceiling(max(Sub_713$CatchRatioSPR40_Upper))+2.3),"Overfished",cex=0.8,pos=4)
dev.off()
#Foregone Catch and Status Plot for WPP 714
thisWPP = "714"
Sub_714 = Ratios[Ratios$WPP==thisWPP,]
Sub_714 = subset(Sub_714,!is.na(Sub_714$FratioSPR40))
Sub_714 = Sub_714[order(Sub_714$CatchRatioSPR40),]
GenSpp_thisWPP = unique(subset(Sub_714,select=c(Genus,Species)))
GenSpp_thisWPP$GenusSpecies = paste(GenSpp_thisWPP$Genus,GenSpp_thisWPP$Species,sep=" ")
png(paste(PATH_highLevelSummaryPlots,paste(paste("CatchRatios_OverfishingColors",thisWPP,sep="_"),".png",sep=""),sep="/"),units="px",width=4200,height=3200,res=600)    #was 3200 and 3200
par(mar=c(14, 5, 4, 2) + 0.1,xpd=NA)
plot(1:dim(Sub_714)[1],Sub_714$CatchRatioSPR40,xaxt="n",yaxt="n",main="",xlab ="",ylab="Yield Relative to Y@SPR40",type="n",ylim=c(0,ceiling(max(Sub_714$CatchRatioSPR40_Upper)))) #pch=21,cex=1.5,col="black",) #,bg=Sub_714$Color_FratioSPR40)
axis(1,at=1:dim(Sub_714)[1],labels=GenSpp_thisWPP$GenusSpecies,las=2,cex.axis=0.8,tick=TRUE)
axis(2,at=0:ceiling(max(Sub_714$CatchRatioSPR40_Upper)),labels=seq(0,ceiling(max(Sub_714$CatchRatioSPR40_Upper)),by=1),las=2,cex.axis=0.8,tick=TRUE)
for(j in 1:dim(Sub_714)[1])
{
  upper.half.circle(j,Sub_714$CatchRatioSPR40[j],r=0.25,nsteps=10000,col=Sub_714$Color_FratioSPR40[j])   #for the larger areas like 714, the r was equal to 0.3
  lower.half.circle(j,Sub_714$CatchRatioSPR40[j],r=0.25,nsteps=10000,col=Sub_714$Color_BratioSPR20[j])
  points(rep(j,2),c(Sub_714$CatchRatioSPR40_Lower[j],Sub_714$CatchRatioSPR40_Upper[j]),col="gray",type="l",lwd=2)
}
xLineVal = seq((1-0.05),(dim(Sub_714)[1]+0.05),by=0.05)
yLineVal = rep(1,length(xLineVal))
points(xLineVal,yLineVal,type="l",lty=2,col="gray",lwd=2)
title("Yield Relative to Y@SPR40")
mtext(paste("FMA ",thisWPP,sep=""))
upper.half.circle(1.0,-(ceiling(max(Sub_714$CatchRatioSPR40_Upper))+ 1.5),r=0.25,nsteps=10000,col="green3")
text(1.1,-(ceiling(max(Sub_714$CatchRatioSPR40_Upper))+1.5),"No Overfishing",cex=0.8,pos=4)
upper.half.circle(1.0,-(ceiling(max(Sub_714$CatchRatioSPR40_Upper))+2.0),r=0.25,nsteps=10000,col="yellow2")
text(1.1,-(ceiling(max(Sub_714$CatchRatioSPR40_Upper))+2.0),"Possible Overfishing",cex=0.8,pos=4)
upper.half.circle(1.0,-(ceiling(max(Sub_714$CatchRatioSPR40_Upper))+2.5),r=0.25,nsteps=10000,col="red3")
text(1.1,-(ceiling(max(Sub_714$CatchRatioSPR40_Upper))+2.5),"Overfishing",cex=0.8,pos=4)
lower.half.circle(12.1,-(ceiling(max(Sub_714$CatchRatioSPR40_Upper))+1.5),r=0.25,nsteps=10000,col="green3")
text(12.2,-(ceiling(max(Sub_714$CatchRatioSPR40_Upper))+1.5),"Not Overfished",cex=0.8,pos=4)
lower.half.circle(12.1,-(ceiling(max(Sub_714$CatchRatioSPR40_Upper))+2.0),r=0.25,nsteps=10000,col="yellow2")
text(12.2,-(ceiling(max(Sub_714$CatchRatioSPR40_Upper))+2.0),"Possibly Overfished",cex=0.8,pos=4)
lower.half.circle(12.1,-(ceiling(max(Sub_714$CatchRatioSPR40_Upper))+2.5),r=0.25,nsteps=10000,col="red3")
text(12.2,-(ceiling(max(Sub_714$CatchRatioSPR40_Upper))+2.5),"Overfished",cex=0.8,pos=4)
dev.off()
#Foregone Catch and Status Plot for WPP 715
thisWPP = "715"
Sub_715 = Ratios[Ratios$WPP==thisWPP,]
Sub_715 = subset(Sub_715,!is.na(Sub_715$FratioSPR40))
Sub_715 = Sub_715[order(Sub_715$CatchRatioSPR40),]
GenSpp_thisWPP = unique(subset(Sub_715,select=c(Genus,Species)))
GenSpp_thisWPP$GenusSpecies = paste(GenSpp_thisWPP$Genus,GenSpp_thisWPP$Species,sep=" ")
png(paste(PATH_highLevelSummaryPlots,paste(paste("CatchRatios_OverfishingColors",thisWPP,sep="_"),".png",sep=""),sep="/"),units="px",width=5200,height=3200,res=600)    #was 3200 and 3200
par(mar=c(14, 5, 4, 2) + 0.1,xpd=NA)
plot(1:dim(Sub_715)[1],Sub_715$CatchRatioSPR40,xaxt="n",yaxt="n",main="",xlab ="",ylab="Yield Relative to Y@SPR40",type="n",ylim=c(0,ceiling(max(Sub_715$CatchRatioSPR40_Upper)))) #pch=21,cex=1.5,col="black",) #,bg=Sub_715$Color_FratioSPR40)
axis(1,at=1:dim(Sub_715)[1],labels=GenSpp_thisWPP$GenusSpecies,las=2,cex.axis=0.8,tick=TRUE)
axis(2,at=0:ceiling(max(Sub_715$CatchRatioSPR40_Upper)),labels=seq(0,ceiling(max(Sub_715$CatchRatioSPR40_Upper)),by=1),las=2,cex.axis=0.8,tick=TRUE)
for(j in 1:dim(Sub_715)[1])
{
  upper.half.circle(j,Sub_715$CatchRatioSPR40[j],r=0.25,nsteps=10000,col=Sub_715$Color_FratioSPR40[j])   #for the larger areas like 715, the r was equal to 0.3
  lower.half.circle(j,Sub_715$CatchRatioSPR40[j],r=0.25,nsteps=10000,col=Sub_715$Color_BratioSPR20[j])
  points(rep(j,2),c(Sub_715$CatchRatioSPR40_Lower[j],Sub_715$CatchRatioSPR40_Upper[j]),col="gray",type="l",lwd=2)
}
xLineVal = seq((1-0.05),(dim(Sub_715)[1]+0.05),by=0.05)
yLineVal = rep(1,length(xLineVal))
points(xLineVal,yLineVal,type="l",lty=2,col="gray",lwd=2)
title("Yield Relative to Y@SPR40")
mtext(paste("FMA ",thisWPP,sep=""))
upper.half.circle(1.0,-(ceiling(max(Sub_715$CatchRatioSPR40_Upper))+ 1.5),r=0.25,nsteps=10000,col="green3")
text(1.1,-(ceiling(max(Sub_715$CatchRatioSPR40_Upper))+1.5),"No Overfishing",cex=0.8,pos=4)
upper.half.circle(1.0,-(ceiling(max(Sub_715$CatchRatioSPR40_Upper))+2.0),r=0.25,nsteps=10000,col="yellow2")
text(1.1,-(ceiling(max(Sub_715$CatchRatioSPR40_Upper))+2.0),"Possible Overfishing",cex=0.8,pos=4)
upper.half.circle(1.0,-(ceiling(max(Sub_715$CatchRatioSPR40_Upper))+2.5),r=0.25,nsteps=10000,col="red3")
text(1.1,-(ceiling(max(Sub_715$CatchRatioSPR40_Upper))+2.5),"Overfishing",cex=0.8,pos=4)
lower.half.circle(9.1,-(ceiling(max(Sub_715$CatchRatioSPR40_Upper))+1.5),r=0.25,nsteps=10000,col="green3")
text(9.2,-(ceiling(max(Sub_715$CatchRatioSPR40_Upper))+1.5),"Not Overfished",cex=0.8,pos=4)
lower.half.circle(9.1,-(ceiling(max(Sub_715$CatchRatioSPR40_Upper))+2.0),r=0.25,nsteps=10000,col="yellow2")
text(9.2,-(ceiling(max(Sub_715$CatchRatioSPR40_Upper))+2.0),"Possibly Overfished",cex=0.8,pos=4)
lower.half.circle(9.1,-(ceiling(max(Sub_715$CatchRatioSPR40_Upper))+2.5),r=0.25,nsteps=10000,col="red3")
text(9.2,-(ceiling(max(Sub_715$CatchRatioSPR40_Upper))+2.5),"Overfished",cex=0.8,pos=4)
dev.off()
#Foregone Catch and Status Plot for WPP 716
thisWPP = "716"
Sub_716 = Ratios[Ratios$WPP==thisWPP,]
Sub_716 = subset(Sub_716,!is.na(Sub_716$FratioSPR40))
Sub_716 = Sub_716[order(Sub_716$CatchRatioSPR40),]
GenSpp_thisWPP = unique(subset(Sub_716,select=c(Genus,Species)))
GenSpp_thisWPP$GenusSpecies = paste(GenSpp_thisWPP$Genus,GenSpp_thisWPP$Species,sep=" ")
png(paste(PATH_highLevelSummaryPlots,paste(paste("CatchRatios_OverfishingColors",thisWPP,sep="_"),".png",sep=""),sep="/"),units="px",width=3200,height=3200,res=600)    #was 3200 and 3200
par(mar=c(14, 5, 4, 2) + 0.1,xpd=NA)
plot(1:dim(Sub_716)[1],Sub_716$CatchRatioSPR40,xaxt="n",yaxt="n",main="",xlab ="",ylab="Yield Relative to Y@SPR40",type="n",ylim=c(0,ceiling(max(Sub_716$CatchRatioSPR40_Upper)))) #pch=21,cex=1.5,col="black",) #,bg=Sub_716$Color_FratioSPR40)
axis(1,at=1:dim(Sub_716)[1],labels=GenSpp_thisWPP$GenusSpecies,las=2,cex.axis=0.8,tick=TRUE)
axis(2,at=0:ceiling(max(Sub_716$CatchRatioSPR40_Upper)),labels=seq(0,ceiling(max(Sub_716$CatchRatioSPR40_Upper)),by=1),las=2,cex.axis=0.8,tick=TRUE)
for(j in 1:dim(Sub_716)[1])
{
  upper.half.circle(j,Sub_716$CatchRatioSPR40[j],r=0.1,nsteps=10000,col=Sub_716$Color_FratioSPR40[j])   #for the larger areas like 716, the r was equal to 0.3
  lower.half.circle(j,Sub_716$CatchRatioSPR40[j],r=0.1,nsteps=10000,col=Sub_716$Color_BratioSPR20[j])
  points(rep(j,2),c(Sub_716$CatchRatioSPR40_Lower[j],Sub_716$CatchRatioSPR40_Upper[j]),col="gray",type="l",lwd=2)
}
xLineVal = seq((1-0.05),(dim(Sub_716)[1]+0.05),by=0.05)
yLineVal = rep(1,length(xLineVal))
points(xLineVal,yLineVal,type="l",lty=2,col="gray",lwd=2)
title("Yield Relative to Y@SPR40")
mtext(paste("FMA ",thisWPP,sep=""))
upper.half.circle(1.0,-(ceiling(max(Sub_716$CatchRatioSPR40_Upper))+ 0.7),r=0.1,nsteps=10000,col="green3")
text(1.1,-(ceiling(max(Sub_716$CatchRatioSPR40_Upper))+0.7),"No Overfishing",cex=0.8,pos=4)
upper.half.circle(1.0,-(ceiling(max(Sub_716$CatchRatioSPR40_Upper))+1.0),r=0.1,nsteps=10000,col="yellow2")
text(1.1,-(ceiling(max(Sub_716$CatchRatioSPR40_Upper))+1.0),"Possible Overfishing",cex=0.8,pos=4)
upper.half.circle(1.0,-(ceiling(max(Sub_716$CatchRatioSPR40_Upper))+1.3),r=0.1,nsteps=10000,col="red3")
text(1.1,-(ceiling(max(Sub_716$CatchRatioSPR40_Upper))+1.3),"Overfishing",cex=0.8,pos=4)
lower.half.circle(3.1,-(ceiling(max(Sub_716$CatchRatioSPR40_Upper))+0.7),r=0.1,nsteps=10000,col="green3")
text(3.2,-(ceiling(max(Sub_716$CatchRatioSPR40_Upper))+0.7),"Not Overfished",cex=0.8,pos=4)
lower.half.circle(3.1,-(ceiling(max(Sub_716$CatchRatioSPR40_Upper))+1.0),r=0.1,nsteps=10000,col="yellow2")
text(3.2,-(ceiling(max(Sub_716$CatchRatioSPR40_Upper))+1.0),"Possibly Overfished",cex=0.8,pos=4)
lower.half.circle(3.1,-(ceiling(max(Sub_716$CatchRatioSPR40_Upper))+1.3),r=0.1,nsteps=10000,col="red3")
text(3.2,-(ceiling(max(Sub_716$CatchRatioSPR40_Upper))+1.3),"Overfished",cex=0.8,pos=4)
dev.off()
#Foregone Catch and Status Plot for WPP 718
thisWPP = "718"
Sub_718 = Ratios[Ratios$WPP==thisWPP,]
Sub_718 = subset(Sub_718,!is.na(Sub_718$FratioSPR40))
Sub_718 = Sub_718[order(Sub_718$CatchRatioSPR40),]
GenSpp_thisWPP = unique(subset(Sub_718,select=c(Genus,Species)))
GenSpp_thisWPP$GenusSpecies = paste(GenSpp_thisWPP$Genus,GenSpp_thisWPP$Species,sep=" ")
png(paste(PATH_highLevelSummaryPlots,paste(paste("CatchRatios_OverfishingColors",thisWPP,sep="_"),".png",sep=""),sep="/"),units="px",width=5200,height=3200,res=600)    #was 3200 and 3200
par(mar=c(14, 5, 4, 2) + 0.1,xpd=NA)
plot(1:dim(Sub_718)[1],Sub_718$CatchRatioSPR40,xaxt="n",yaxt="n",main="",xlab ="",ylab="Yield Relative to Y@SPR40",type="n",ylim=c(0,ceiling(max(Sub_718$CatchRatioSPR40_Upper)))) #pch=21,cex=1.5,col="black",) #,bg=Sub_718$Color_FratioSPR40)
axis(1,at=1:dim(Sub_718)[1],labels=GenSpp_thisWPP$GenusSpecies,las=2,cex.axis=0.8,tick=TRUE)
axis(2,at=0:ceiling(max(Sub_718$CatchRatioSPR40_Upper)),labels=seq(0,ceiling(max(Sub_718$CatchRatioSPR40_Upper)),by=1),las=2,cex.axis=0.8,tick=TRUE)
for(j in 1:dim(Sub_718)[1])
{
  upper.half.circle(j,Sub_718$CatchRatioSPR40[j],r=0.25,nsteps=10000,col=Sub_718$Color_FratioSPR40[j])   #for the larger areas like 718, the r was equal to 0.3
  lower.half.circle(j,Sub_718$CatchRatioSPR40[j],r=0.25,nsteps=10000,col=Sub_718$Color_BratioSPR20[j])
  points(rep(j,2),c(Sub_718$CatchRatioSPR40_Lower[j],Sub_718$CatchRatioSPR40_Upper[j]),col="gray",type="l",lwd=2)
}
xLineVal = seq((1-0.05),(dim(Sub_718)[1]+0.05),by=0.05)
yLineVal = rep(1,length(xLineVal))
points(xLineVal,yLineVal,type="l",lty=2,col="gray",lwd=2)
title("Yield Relative to Y@SPR40")
mtext(paste("FMA ",thisWPP,sep=""))
upper.half.circle(1.0,-(ceiling(max(Sub_718$CatchRatioSPR40_Upper))+ 1.5),r=0.25,nsteps=10000,col="green3")
text(1.1,-(ceiling(max(Sub_718$CatchRatioSPR40_Upper))+1.5),"No Overfishing",cex=0.8,pos=4)
upper.half.circle(1.0,-(ceiling(max(Sub_718$CatchRatioSPR40_Upper))+2.0),r=0.25,nsteps=10000,col="yellow2")
text(1.1,-(ceiling(max(Sub_718$CatchRatioSPR40_Upper))+2.0),"Possible Overfishing",cex=0.8,pos=4)
upper.half.circle(1.0,-(ceiling(max(Sub_718$CatchRatioSPR40_Upper))+2.5),r=0.25,nsteps=10000,col="red3")
text(1.1,-(ceiling(max(Sub_718$CatchRatioSPR40_Upper))+2.5),"Overfishing",cex=0.8,pos=4)
lower.half.circle(8.1,-(ceiling(max(Sub_718$CatchRatioSPR40_Upper))+1.5),r=0.25,nsteps=10000,col="green3")
text(8.2,-(ceiling(max(Sub_718$CatchRatioSPR40_Upper))+1.5),"Not Overfished",cex=0.8,pos=4)
lower.half.circle(8.1,-(ceiling(max(Sub_718$CatchRatioSPR40_Upper))+2.0),r=0.25,nsteps=10000,col="yellow2")
text(8.2,-(ceiling(max(Sub_718$CatchRatioSPR40_Upper))+2.0),"Possibly Overfished",cex=0.8,pos=4)
lower.half.circle(8.1,-(ceiling(max(Sub_718$CatchRatioSPR40_Upper))+2.5),r=0.25,nsteps=10000,col="red3")
text(8.2,-(ceiling(max(Sub_718$CatchRatioSPR40_Upper))+2.5),"Overfished",cex=0.8,pos=4)
dev.off()
#Foregone Catch and Status Plot for WPP EEZ
thisWPP = "EEZ"
Sub_EEZ = Ratios[Ratios$WPP==thisWPP,]
Sub_EEZ = subset(Sub_EEZ,!is.na(Sub_EEZ$FratioSPR40))
Sub_EEZ = Sub_EEZ[order(Sub_EEZ$CatchRatioSPR40),]
GenSpp_thisWPP = unique(subset(Sub_EEZ,select=c(Genus,Species)))
GenSpp_thisWPP$GenusSpecies = paste(GenSpp_thisWPP$Genus,GenSpp_thisWPP$Species,sep=" ")
png(paste(PATH_highLevelSummaryPlots,paste(paste("CatchRatios_OverfishingColors",thisWPP,sep="_"),".png",sep=""),sep="/"),units="px",width=5200,height=3200,res=600)    #was 3200 and 3200
par(mar=c(14, 5, 4, 2) + 0.1,xpd=NA)
plot(1:dim(Sub_EEZ)[1],Sub_EEZ$CatchRatioSPR40,xaxt="n",yaxt="n",main="",xlab ="",ylab="Yield Relative to Y@SPR40",type="n",ylim=c(0,ceiling(max(Sub_EEZ$CatchRatioSPR40_Upper)))) #pch=21,cex=1.5,col="black",) #,bg=Sub_EEZ$Color_FratioSPR40)
axis(1,at=1:dim(Sub_EEZ)[1],labels=GenSpp_thisWPP$GenusSpecies,las=2,cex.axis=0.8,tick=TRUE)
axis(2,at=0:ceiling(max(Sub_EEZ$CatchRatioSPR40_Upper)),labels=seq(0,ceiling(max(Sub_EEZ$CatchRatioSPR40_Upper)),by=1),las=2,cex.axis=0.8,tick=TRUE)
for(j in 1:dim(Sub_EEZ)[1])
{
  upper.half.circle(j,Sub_EEZ$CatchRatioSPR40[j],r=0.25,nsteps=10000,col=Sub_EEZ$Color_FratioSPR40[j])   #for the larger areas like EEZ, the r was equal to 0.3
  lower.half.circle(j,Sub_EEZ$CatchRatioSPR40[j],r=0.25,nsteps=10000,col=Sub_EEZ$Color_BratioSPR20[j])
  points(rep(j,2),c(Sub_EEZ$CatchRatioSPR40_Lower[j],Sub_EEZ$CatchRatioSPR40_Upper[j]),col="gray",type="l",lwd=2)
}
xLineVal = seq((1-0.05),(dim(Sub_EEZ)[1]+0.05),by=0.05)
yLineVal = rep(1,length(xLineVal))
points(xLineVal,yLineVal,type="l",lty=2,col="gray",lwd=2)
title("Yield Relative to Y@SPR40")
mtext("Data Aggregated Across All FMAs")
upper.half.circle(1.0,-(ceiling(max(Sub_EEZ$CatchRatioSPR40_Upper))+ 1.5),r=0.25,nsteps=10000,col="green3")
text(1.1,-(ceiling(max(Sub_EEZ$CatchRatioSPR40_Upper))+1.5),"No Overfishing",cex=0.8,pos=4)
upper.half.circle(1.0,-(ceiling(max(Sub_EEZ$CatchRatioSPR40_Upper))+2.0),r=0.25,nsteps=10000,col="yellow2")
text(1.1,-(ceiling(max(Sub_EEZ$CatchRatioSPR40_Upper))+2.0),"Possible Overfishing",cex=0.8,pos=4)
upper.half.circle(1.0,-(ceiling(max(Sub_EEZ$CatchRatioSPR40_Upper))+2.5),r=0.25,nsteps=10000,col="red3")
text(1.1,-(ceiling(max(Sub_EEZ$CatchRatioSPR40_Upper))+2.5),"Overfishing",cex=0.8,pos=4)
lower.half.circle(10.1,-(ceiling(max(Sub_EEZ$CatchRatioSPR40_Upper))+1.5),r=0.25,nsteps=10000,col="green3")
text(10.2,-(ceiling(max(Sub_EEZ$CatchRatioSPR40_Upper))+1.5),"Not Overfished",cex=0.8,pos=4)
lower.half.circle(10.1,-(ceiling(max(Sub_EEZ$CatchRatioSPR40_Upper))+2.0),r=0.25,nsteps=10000,col="yellow2")
text(10.2,-(ceiling(max(Sub_EEZ$CatchRatioSPR40_Upper))+2.0),"Possibly Overfished",cex=0.8,pos=4)
lower.half.circle(10.1,-(ceiling(max(Sub_EEZ$CatchRatioSPR40_Upper))+2.5),r=0.25,nsteps=10000,col="red3")
text(10.2,-(ceiling(max(Sub_EEZ$CatchRatioSPR40_Upper))+2.5),"Overfished",cex=0.8,pos=4)
dev.off()

########### Get Information on Rebuilding Times at SPR 20%, SPR 30% and SPR 40% ########################
debugRebuildTimes=FALSE
SummaryProjectionScenarios = read.table(paste(PATH_output,"SummaryProjectionScenarios.csv",sep="/"),header=TRUE,sep=",",as.is=TRUE)
SummaryProjectionScenarios = SummaryProjectionScenarios[SummaryProjectionScenarios$populationReconstructionMethod!="SolveBaranov",]   #TEMP CODE LINE: do NOT include SolveBaranov or TropFishR
SummaryProjectionScenarios = SummaryProjectionScenarios[SummaryProjectionScenarios$populationReconstructionMethod!="TropFishR",]      #TEMP CODE LINE: do NOT include SolveBaranov or TropFishR
ListOfPopulations_WithProjections = readRDS(file=paste(PATH_output,"ListOfPopulations_WithProjections.RData",sep="/"),refhook=NULL)
#Calculate Time to Rebuild the Stocks to SPR 20%
RunsSPR20 = SummaryProjectionScenarios[SummaryProjectionScenarios$FishingMortaltiyProjectionScenario=="F_at_SPR20",]
RunsSPR20 = RunsSPR20[RunsSPR20$ProjectionType=="RegularF",]
names(RunsSPR20)[names(RunsSPR20)=="BiomassAtEnd"]="BatSPR20"
BiomassEquilibriumAtFcurrent = SummaryProjectionScenarios[SummaryProjectionScenarios$FishingMortaltiyProjectionScenario=="CurrentF",]
BiomassEquilibriumAtFcurrent = subset(BiomassEquilibriumAtFcurrent,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,ProjectionType,BiomassAtEnd))
names(BiomassEquilibriumAtFcurrent)[names(BiomassEquilibriumAtFcurrent)=="BiomassAtEnd"]="BiomassCurrentEquilib"
RunsSPR20 = merge(RunsSPR20,BiomassEquilibriumAtFcurrent,all.x=TRUE,by=c("Genus","Species","WPP","CatchMSY_Scenario","CatchMSY_Scenario_Bound","LifeHistory_Scenario","populationReconstructionMethod","SteepnessAssumed","ProjectionType"))
RunsSPR20$B_ratioSPR20 = RunsSPR20$B_estimated/RunsSPR20$BatSPR20
RunsSPR20$timeToRebuildToSPR20=0
for(i in 1:dim(RunsSPR20)[1])
{
  timeToRebuildToSPR20=0
 if(RunsSPR20$B_ratioSPR20[i]>=1)
 {
   timeToRebuildToSPR20 = NA     #indicates that this particular run is not overfished and does not need to be rebuilt
 }
 if(RunsSPR20$B_ratioSPR20[i]<1)
 {
    RunsSPR20_this = RunsSPR20[i,]
    if(debugRebuildTimes==TRUE)
    {
      windows()
      plot(as.numeric(names(ListOfPopulations_WithProjections[[RunsSPR20_this$ListIndexNumber]]$ProjectedBiomassYear)),ListOfPopulations_WithProjections[[RunsSPR20_this$ListIndexNumber]]$ProjectedBiomassYear,pch=16,xlab="Year",ylab="Biomass (kg)")
      abline(h=RunsSPR20_this$BatSPR20,col="gray",lty=3,lwd=3)
    }
    #first, see if it is even going to recover. We will know if it is recovering if the slope of the log of biomass fit to a linear model is positive
    CheckSlope = lm(log(ListOfPopulations_WithProjections[[RunsSPR20_this$ListIndexNumber]]$ProjectedBiomassYear) ~ as.numeric(names(ListOfPopulations_WithProjections[[RunsSPR20_this$ListIndexNumber]]$ProjectedBiomassYear)))
    timeToRebuildToSPR20 = 0
    if(as.numeric(CheckSlope$coef[2])<0)
    {
      timeToRebuildToSPR20 = NA   #indicating that rebuilding is not needed (we are higher than SPR target)
    }
    if(as.numeric(CheckSlope$coef[2])>=0)
    {
      ProximityToTarget = ListOfPopulations_WithProjections[[RunsSPR20_this$ListIndexNumber]]$ProjectedBiomassYear/RunsSPR20_this$BatSPR20
      if(ProximityToTarget[length(ProximityToTarget)] < 0.995)
      {
        timeToRebuildToSPR20 = -999   #indicating that the stock will not rebuild
      }
      if(ProximityToTarget[length(ProximityToTarget)] >= 0.995)
      {
        ProximityToTarget_onceRebuilt = subset(ProximityToTarget,ProximityToTarget>=0.995)
        yearRebuilt = as.numeric(names(ProximityToTarget_onceRebuilt[1]))
        if(debugRebuildTimes==TRUE)
        {
          abline(v=yearRebuilt,col="gray",lty=3,lwd=3)
        }
        timeToRebuildToSPR20 = yearRebuilt - (min(as.numeric(names(ListOfPopulations_WithProjections[[RunsSPR20_this$ListIndexNumber]]$ProjectedBiomassYear)))-1)
      }
    }
  }
  RunsSPR20$timeToRebuildToSPR20[i] = timeToRebuildToSPR20
  if(i%%500==0)
  {
    print(paste("Working on state of nature for SPR20 ",i," of ",dim(RunsSPR20)[1],"...",sep=""))
    flush.console()
  }
}
write.table(RunsSPR20,paste(PATH_output,"RunsSPR20_withRebuildTimes.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")

#Calculate Time to Rebuild the Stocks to SPR 30%
RunsSPR30 = SummaryProjectionScenarios[SummaryProjectionScenarios$FishingMortaltiyProjectionScenario=="F_at_SPR30",]
RunsSPR30 = RunsSPR30[RunsSPR30$ProjectionType=="RegularF",]
names(RunsSPR30)[names(RunsSPR30)=="BiomassAtEnd"]="BatSPR30"
BiomassEquilibriumAtFcurrent = SummaryProjectionScenarios[SummaryProjectionScenarios$FishingMortaltiyProjectionScenario=="CurrentF",]
BiomassEquilibriumAtFcurrent = subset(BiomassEquilibriumAtFcurrent,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,ProjectionType,BiomassAtEnd))
names(BiomassEquilibriumAtFcurrent)[names(BiomassEquilibriumAtFcurrent)=="BiomassAtEnd"]="BiomassCurrentEquilib"
RunsSPR30 = merge(RunsSPR30,BiomassEquilibriumAtFcurrent,all.x=TRUE,by=c("Genus","Species","WPP","CatchMSY_Scenario","CatchMSY_Scenario_Bound","LifeHistory_Scenario","populationReconstructionMethod","SteepnessAssumed","ProjectionType"))
RunsSPR30$B_ratioSPR30 = RunsSPR30$B_estimated/RunsSPR30$BatSPR30
RunsSPR30$timeToRebuildToSPR30=0
for(i in 1:dim(RunsSPR30)[1])
{
  timeToRebuildToSPR30=0
 if(RunsSPR30$B_ratioSPR30[i]>=1)
 {
   timeToRebuildToSPR30 = NA     #indicates that this particular run is not overfished and does not need to be rebuilt
 }
 if(RunsSPR30$B_ratioSPR30[i]<1)
 {
    RunsSPR30_this = RunsSPR30[i,]
    if(debugRebuildTimes==TRUE)
    {
      windows()
      plot(as.numeric(names(ListOfPopulations_WithProjections[[RunsSPR30_this$ListIndexNumber]]$ProjectedBiomassYear)),ListOfPopulations_WithProjections[[RunsSPR30_this$ListIndexNumber]]$ProjectedBiomassYear,pch=16,xlab="Year",ylab="Biomass (kg)")
      abline(h=RunsSPR30_this$BatSPR30,col="gray",lty=3,lwd=3)
    }
    #first, see if it is even going to recover. We will know if it is recovering if the slope of the log of biomass fit to a linear model is positive
    CheckSlope = lm(log(ListOfPopulations_WithProjections[[RunsSPR30_this$ListIndexNumber]]$ProjectedBiomassYear) ~ as.numeric(names(ListOfPopulations_WithProjections[[RunsSPR30_this$ListIndexNumber]]$ProjectedBiomassYear)))
    timeToRebuildToSPR30 = 0
    if(as.numeric(CheckSlope$coef[2])<0)
    {
      timeToRebuildToSPR30 = NA   #indicating that rebuilding is not needed (we are higher than SPR target)
    }
    if(as.numeric(CheckSlope$coef[2])>=0)
    {
      ProximityToTarget = ListOfPopulations_WithProjections[[RunsSPR30_this$ListIndexNumber]]$ProjectedBiomassYear/RunsSPR30_this$BatSPR30
      if(ProximityToTarget[length(ProximityToTarget)] < 0.995)
      {
        timeToRebuildToSPR30 = -999   #indicating that the stock will not rebuild
      }
      if(ProximityToTarget[length(ProximityToTarget)] >= 0.995)
      {
        ProximityToTarget_onceRebuilt = subset(ProximityToTarget,ProximityToTarget>=0.995)
        yearRebuilt = as.numeric(names(ProximityToTarget_onceRebuilt[1]))
        if(debugRebuildTimes==TRUE)
        {
          abline(v=yearRebuilt,col="gray",lty=3,lwd=3)
        }
        timeToRebuildToSPR30 = yearRebuilt - (min(as.numeric(names(ListOfPopulations_WithProjections[[RunsSPR30_this$ListIndexNumber]]$ProjectedBiomassYear)))-1)
      }
    }
  }
  RunsSPR30$timeToRebuildToSPR30[i] = timeToRebuildToSPR30
  if(i%%500==0)
  {
    print(paste("Working on state of nature for SPR30 ",i," of ",dim(RunsSPR30)[1],"...",sep=""))
    flush.console()
  }
}
write.table(RunsSPR30,paste(PATH_output,"RunsSPR30_withRebuildTimes.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")

#Calculate Time to Rebuild the Stocks to SPR 40%
RunsSPR40 = SummaryProjectionScenarios[SummaryProjectionScenarios$FishingMortaltiyProjectionScenario=="F_at_SPR40",]
RunsSPR40 = RunsSPR40[RunsSPR40$ProjectionType=="RegularF",]
names(RunsSPR40)[names(RunsSPR40)=="BiomassAtEnd"]="BatSPR40"
BiomassEquilibriumAtFcurrent = SummaryProjectionScenarios[SummaryProjectionScenarios$FishingMortaltiyProjectionScenario=="CurrentF",]
BiomassEquilibriumAtFcurrent = subset(BiomassEquilibriumAtFcurrent,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,ProjectionType,BiomassAtEnd))
names(BiomassEquilibriumAtFcurrent)[names(BiomassEquilibriumAtFcurrent)=="BiomassAtEnd"]="BiomassCurrentEquilib"
RunsSPR40 = merge(RunsSPR40,BiomassEquilibriumAtFcurrent,all.x=TRUE,by=c("Genus","Species","WPP","CatchMSY_Scenario","CatchMSY_Scenario_Bound","LifeHistory_Scenario","populationReconstructionMethod","SteepnessAssumed","ProjectionType"))
RunsSPR40$B_ratioSPR40 = RunsSPR40$B_estimated/RunsSPR40$BatSPR40
RunsSPR40$timeToRebuildToSPR40=0
for(i in 1:dim(RunsSPR40)[1])
{
  timeToRebuildToSPR40=0
 if(RunsSPR40$B_ratioSPR40[i]>=1)
 {
   timeToRebuildToSPR40 = NA     #indicates that this particular run is not overfished and does not need to be rebuilt
 }
 if(RunsSPR40$B_ratioSPR40[i]<1)
 {
    RunsSPR40_this = RunsSPR40[i,]
    if(debugRebuildTimes==TRUE)
    {
      windows()
      plot(as.numeric(names(ListOfPopulations_WithProjections[[RunsSPR40_this$ListIndexNumber]]$ProjectedBiomassYear)),ListOfPopulations_WithProjections[[RunsSPR40_this$ListIndexNumber]]$ProjectedBiomassYear,pch=16,xlab="Year",ylab="Biomass (kg)")
      abline(h=RunsSPR40_this$BatSPR40,col="gray",lty=3,lwd=3)
    }
    #first, see if it is even going to recover. We will know if it is recovering if the slope of the log of biomass fit to a linear model is positive
    CheckSlope = lm(log(ListOfPopulations_WithProjections[[RunsSPR40_this$ListIndexNumber]]$ProjectedBiomassYear) ~ as.numeric(names(ListOfPopulations_WithProjections[[RunsSPR40_this$ListIndexNumber]]$ProjectedBiomassYear)))
    timeToRebuildToSPR40 = 0
    if(as.numeric(CheckSlope$coef[2])<0)
    {
      timeToRebuildToSPR40 = NA   #indicating that rebuilding is not needed (we are higher than SPR target)
    }
    if(as.numeric(CheckSlope$coef[2])>=0)
    {
      ProximityToTarget = ListOfPopulations_WithProjections[[RunsSPR40_this$ListIndexNumber]]$ProjectedBiomassYear/RunsSPR40_this$BatSPR40
      if(ProximityToTarget[length(ProximityToTarget)] < 0.995)
      {
        timeToRebuildToSPR40 = -999   #indicating that the stock will not rebuild
      }
      if(ProximityToTarget[length(ProximityToTarget)] >= 0.995)
      {
        ProximityToTarget_onceRebuilt = subset(ProximityToTarget,ProximityToTarget>=0.995)
        yearRebuilt = as.numeric(names(ProximityToTarget_onceRebuilt[1]))
        if(debugRebuildTimes==TRUE)
        {
          abline(v=yearRebuilt,col="gray",lty=3,lwd=3)
        }
        timeToRebuildToSPR40 = yearRebuilt - (min(as.numeric(names(ListOfPopulations_WithProjections[[RunsSPR40_this$ListIndexNumber]]$ProjectedBiomassYear)))-1)
      }
    }
  }
  RunsSPR40$timeToRebuildToSPR40[i] = timeToRebuildToSPR40
  if(i%%500==0)
  {
    print(paste("Working on state of nature for SPR40 ",i," of ",dim(RunsSPR40)[1],"...",sep=""))
    flush.console()
  }
}
write.table(RunsSPR40,paste(PATH_output,"RunsSPR40_withRebuildTimes.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")

#Percent Change in F at SPR30%
RunsSPR30 = read.table(paste(PATH_output,"RunsSPR30_withRebuildTimes.csv",sep="/"),header=TRUE,sep=",")
RunsSPR30$PercentChangeF = ((RunsSPR30$F_project - RunsSPR30$F_estimated)/RunsSPR30$F_estimated)*100
RunsSPR30$PercentChangeF[RunsSPR30$OverfishingCategorySPR30!=3]=NA
PercChangeF_SPR30 = subset(RunsSPR30,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,ProjectionType,PercentChangeF))
PercChangeF_SPR30 = PercChangeF_SPR30[!is.na(PercChangeF_SPR30$PercentChangeF) & PercChangeF_SPR30$PercentChangeF<0,]
PercChangeF_SPR30_mean = aggregate.data.frame(PercChangeF_SPR30$PercentChangeF,by=list(PercChangeF_SPR30$Genus,PercChangeF_SPR30$Species,PercChangeF_SPR30$WPP),FUN=mean)
names(PercChangeF_SPR30_mean) = c("Genus","Species","WPP","PercentChangeF_spr30_mean")
PercChangeF_SPR30_sd = aggregate.data.frame(PercChangeF_SPR30$PercentChangeF,by=list(PercChangeF_SPR30$Genus,PercChangeF_SPR30$Species,PercChangeF_SPR30$WPP),FUN=sd)
names(PercChangeF_SPR30_sd) = c("Genus","Species","WPP","PercentChangeF_spr30_sd")
PercChangeF_SPR30 = merge(PercChangeF_SPR30_mean,PercChangeF_SPR30_sd,all.x=TRUE,by=c("Genus","Species","WPP"))
PercChangeF_SPR30$Upper = PercChangeF_SPR30$PercentChangeF_spr30_mean + PercChangeF_SPR30$PercentChangeF_spr30_sd
PercChangeF_SPR30$Lower = PercChangeF_SPR30$PercentChangeF_spr30_mean - PercChangeF_SPR30$PercentChangeF_spr30_sd
PercChangeF_SPR30$UncertaintyRange = paste(round(PercChangeF_SPR30$Lower,0)," to ",round(PercChangeF_SPR30$Upper,0),sep="")
GenusSpeciesAreaStatusSPR30 = subset(RunsSPR30,select=c(Genus,Species,WPP))
GenusSpeciesAreaStatusSPR30 = unique(GenusSpeciesAreaStatusSPR30)
PercChangeF_SPR30 = merge(GenusSpeciesAreaStatusSPR30,PercChangeF_SPR30,all.x=TRUE,by=c("Genus","Species","WPP"))
PercChangeF_SPR30_mean = subset(PercChangeF_SPR30,select=c(Genus,Species,WPP,PercentChangeF_spr30_mean))
PercChangeF_SPR30_mean$GenusSpecies = paste(PercChangeF_SPR30_mean$Genus,PercChangeF_SPR30_mean$Species,sep="_")
GenusSpeciesCombineKey = subset(PercChangeF_SPR30_mean,select=c(Genus,Species,GenusSpecies))
GenusSpeciesCombineKey = unique(GenusSpeciesCombineKey)
PercChangeF_SPR30_mean$Genus=NULL
PercChangeF_SPR30_mean$Species=NULL
PercChangeF_SPR30_mean$PercentChangeF_spr30_mean = round(PercChangeF_SPR30_mean$PercentChangeF_spr30_mean,0)
PercChangeF_Table_mean_reshape_spr30 = reshape(PercChangeF_SPR30_mean,v.names="PercentChangeF_spr30_mean",idvar="GenusSpecies",timevar="WPP",direction="wide")
names(PercChangeF_Table_mean_reshape_spr30) = gsub("PercentChangeF_spr30_mean.","",names(PercChangeF_Table_mean_reshape_spr30))
PercChangeF_Table_mean_reshape_spr30 = merge(PercChangeF_Table_mean_reshape_spr30,GenusSpeciesCombineKey,all.x=TRUE,by=c("GenusSpecies"))
PercChangeF_Table_mean_reshape_spr30$GenusSpecies=NULL
write.table(PercChangeF_Table_mean_reshape_spr30,paste(PATH_highLevelSummaryPlots,"PercChangeF_Table_mean_reshape_spr30.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
PercChangeF_SPR30_range = subset(PercChangeF_SPR30,select=c(Genus,Species,WPP,UncertaintyRange))
PercChangeF_SPR30_range$GenusSpecies = paste(PercChangeF_SPR30_range$Genus,PercChangeF_SPR30_range$Species,sep="_")
PercChangeF_SPR30_range$Genus=NULL
PercChangeF_SPR30_range$Species=NULL
PercChangeF_Table_range_reshape_spr30 = reshape(PercChangeF_SPR30_range,v.names="UncertaintyRange",idvar="GenusSpecies",timevar="WPP",direction="wide")
names(PercChangeF_Table_range_reshape_spr30) = gsub("UncertaintyRange.","",names(PercChangeF_Table_range_reshape_spr30))
PercChangeF_Table_range_reshape_spr30 = merge(PercChangeF_Table_range_reshape_spr30,GenusSpeciesCombineKey,all.x=TRUE,by=c("GenusSpecies"))
PercChangeF_Table_range_reshape_spr30$GenusSpecies=NULL
write.table(PercChangeF_Table_range_reshape_spr30,paste(PATH_highLevelSummaryPlots,"PercChangeF_Table_range_reshape_spr30.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")

#Percent Change in F at SPR20%
RunsSPR20 = read.table(paste(PATH_output,"RunsSPR20_withRebuildTimes.csv",sep="/"),header=TRUE,sep=",")
RunsSPR20$PercentChangeF = ((RunsSPR20$F_project - RunsSPR20$F_estimated)/RunsSPR20$F_estimated)*100
RunsSPR20$PercentChangeF[RunsSPR20$OverfishingCategorySPR20!=3]=NA
PercChangeF_SPR20 = subset(RunsSPR20,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,ProjectionType,PercentChangeF))
PercChangeF_SPR20 = PercChangeF_SPR20[!is.na(PercChangeF_SPR20$PercentChangeF) & PercChangeF_SPR20$PercentChangeF<0,]
PercChangeF_SPR20_mean = aggregate.data.frame(PercChangeF_SPR20$PercentChangeF,by=list(PercChangeF_SPR20$Genus,PercChangeF_SPR20$Species,PercChangeF_SPR20$WPP),FUN=mean)
names(PercChangeF_SPR20_mean) = c("Genus","Species","WPP","PercentChangeF_SPR20_mean")
PercChangeF_SPR20_sd = aggregate.data.frame(PercChangeF_SPR20$PercentChangeF,by=list(PercChangeF_SPR20$Genus,PercChangeF_SPR20$Species,PercChangeF_SPR20$WPP),FUN=sd)
names(PercChangeF_SPR20_sd) = c("Genus","Species","WPP","PercentChangeF_SPR20_sd")
PercChangeF_SPR20 = merge(PercChangeF_SPR20_mean,PercChangeF_SPR20_sd,all.x=TRUE,by=c("Genus","Species","WPP"))
PercChangeF_SPR20$Upper = PercChangeF_SPR20$PercentChangeF_SPR20_mean + PercChangeF_SPR20$PercentChangeF_SPR20_sd
PercChangeF_SPR20$Lower = PercChangeF_SPR20$PercentChangeF_SPR20_mean - PercChangeF_SPR20$PercentChangeF_SPR20_sd
PercChangeF_SPR20$UncertaintyRange = paste(round(PercChangeF_SPR20$Lower,0)," to ",round(PercChangeF_SPR20$Upper,0),sep="")
GenusSpeciesAreaStatusSPR20 = subset(RunsSPR20,select=c(Genus,Species,WPP))
GenusSpeciesAreaStatusSPR20 = unique(GenusSpeciesAreaStatusSPR20)
PercChangeF_SPR20 = merge(GenusSpeciesAreaStatusSPR20,PercChangeF_SPR20,all.x=TRUE,by=c("Genus","Species","WPP"))
PercChangeF_SPR20_mean = subset(PercChangeF_SPR20,select=c(Genus,Species,WPP,PercentChangeF_SPR20_mean))
PercChangeF_SPR20_mean$GenusSpecies = paste(PercChangeF_SPR20_mean$Genus,PercChangeF_SPR20_mean$Species,sep="_")
GenusSpeciesCombineKey = subset(PercChangeF_SPR20_mean,select=c(Genus,Species,GenusSpecies))
GenusSpeciesCombineKey = unique(GenusSpeciesCombineKey)
PercChangeF_SPR20_mean$Genus=NULL
PercChangeF_SPR20_mean$Species=NULL
PercChangeF_SPR20_mean$PercentChangeF_SPR20_mean = round(PercChangeF_SPR20_mean$PercentChangeF_SPR20_mean,0)
PercChangeF_Table_mean_reshape_SPR20 = reshape(PercChangeF_SPR20_mean,v.names="PercentChangeF_SPR20_mean",idvar="GenusSpecies",timevar="WPP",direction="wide")
names(PercChangeF_Table_mean_reshape_SPR20) = gsub("PercentChangeF_SPR20_mean.","",names(PercChangeF_Table_mean_reshape_SPR20))
PercChangeF_Table_mean_reshape_SPR20 = merge(PercChangeF_Table_mean_reshape_SPR20,GenusSpeciesCombineKey,all.x=TRUE,by=c("GenusSpecies"))
PercChangeF_Table_mean_reshape_SPR20$GenusSpecies=NULL
write.table(PercChangeF_Table_mean_reshape_SPR20,paste(PATH_highLevelSummaryPlots,"PercChangeF_Table_mean_reshape_SPR20.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
PercChangeF_SPR20_range = subset(PercChangeF_SPR20,select=c(Genus,Species,WPP,UncertaintyRange))
PercChangeF_SPR20_range$GenusSpecies = paste(PercChangeF_SPR20_range$Genus,PercChangeF_SPR20_range$Species,sep="_")
PercChangeF_SPR20_range$Genus=NULL
PercChangeF_SPR20_range$Species=NULL
PercChangeF_Table_range_reshape_SPR20 = reshape(PercChangeF_SPR20_range,v.names="UncertaintyRange",idvar="GenusSpecies",timevar="WPP",direction="wide")
names(PercChangeF_Table_range_reshape_SPR20) = gsub("UncertaintyRange.","",names(PercChangeF_Table_range_reshape_SPR20))
PercChangeF_Table_range_reshape_SPR20 = merge(PercChangeF_Table_range_reshape_SPR20,GenusSpeciesCombineKey,all.x=TRUE,by=c("GenusSpecies"))
PercChangeF_Table_range_reshape_SPR20$GenusSpecies=NULL
write.table(PercChangeF_Table_range_reshape_SPR20,paste(PATH_highLevelSummaryPlots,"PercChangeF_Table_range_reshape_SPR20.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")

#Percent Change in F at SPR40%
RunsSPR40 = read.table(paste(PATH_output,"RunsSPR40_withRebuildTimes.csv",sep="/"),header=TRUE,sep=",")
RunsSPR40$PercentChangeF = ((RunsSPR40$F_project - RunsSPR40$F_estimated)/RunsSPR40$F_estimated)*100
RunsSPR40$PercentChangeF[RunsSPR40$OverfishingCategorySPR40!=3]=NA
PercChangeF_SPR40 = subset(RunsSPR40,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,ProjectionType,PercentChangeF))
PercChangeF_SPR40 = PercChangeF_SPR40[!is.na(PercChangeF_SPR40$PercentChangeF) & PercChangeF_SPR40$PercentChangeF<0,]
PercChangeF_SPR40_mean = aggregate.data.frame(PercChangeF_SPR40$PercentChangeF,by=list(PercChangeF_SPR40$Genus,PercChangeF_SPR40$Species,PercChangeF_SPR40$WPP),FUN=mean)
names(PercChangeF_SPR40_mean) = c("Genus","Species","WPP","PercentChangeF_SPR40_mean")
PercChangeF_SPR40_sd = aggregate.data.frame(PercChangeF_SPR40$PercentChangeF,by=list(PercChangeF_SPR40$Genus,PercChangeF_SPR40$Species,PercChangeF_SPR40$WPP),FUN=sd)
names(PercChangeF_SPR40_sd) = c("Genus","Species","WPP","PercentChangeF_SPR40_sd")
PercChangeF_SPR40 = merge(PercChangeF_SPR40_mean,PercChangeF_SPR40_sd,all.x=TRUE,by=c("Genus","Species","WPP"))
PercChangeF_SPR40$Upper = PercChangeF_SPR40$PercentChangeF_SPR40_mean + PercChangeF_SPR40$PercentChangeF_SPR40_sd
PercChangeF_SPR40$Lower = PercChangeF_SPR40$PercentChangeF_SPR40_mean - PercChangeF_SPR40$PercentChangeF_SPR40_sd
PercChangeF_SPR40$UncertaintyRange = paste(round(PercChangeF_SPR40$Lower,0)," to ",round(PercChangeF_SPR40$Upper,0),sep="")
GenusSpeciesAreaStatusSPR40 = subset(RunsSPR40,select=c(Genus,Species,WPP))
GenusSpeciesAreaStatusSPR40 = unique(GenusSpeciesAreaStatusSPR40)
PercChangeF_SPR40 = merge(GenusSpeciesAreaStatusSPR40,PercChangeF_SPR40,all.x=TRUE,by=c("Genus","Species","WPP"))
PercChangeF_SPR40_mean = subset(PercChangeF_SPR40,select=c(Genus,Species,WPP,PercentChangeF_SPR40_mean))
PercChangeF_SPR40_mean$GenusSpecies = paste(PercChangeF_SPR40_mean$Genus,PercChangeF_SPR40_mean$Species,sep="_")
GenusSpeciesCombineKey = subset(PercChangeF_SPR40_mean,select=c(Genus,Species,GenusSpecies))
GenusSpeciesCombineKey = unique(GenusSpeciesCombineKey)
PercChangeF_SPR40_mean$Genus=NULL
PercChangeF_SPR40_mean$Species=NULL
PercChangeF_SPR40_mean$PercentChangeF_SPR40_mean = round(PercChangeF_SPR40_mean$PercentChangeF_SPR40_mean,0)
PercChangeF_Table_mean_reshape_SPR40 = reshape(PercChangeF_SPR40_mean,v.names="PercentChangeF_SPR40_mean",idvar="GenusSpecies",timevar="WPP",direction="wide")
names(PercChangeF_Table_mean_reshape_SPR40) = gsub("PercentChangeF_SPR40_mean.","",names(PercChangeF_Table_mean_reshape_SPR40))
PercChangeF_Table_mean_reshape_SPR40 = merge(PercChangeF_Table_mean_reshape_SPR40,GenusSpeciesCombineKey,all.x=TRUE,by=c("GenusSpecies"))
PercChangeF_Table_mean_reshape_SPR40$GenusSpecies=NULL
write.table(PercChangeF_Table_mean_reshape_SPR40,paste(PATH_highLevelSummaryPlots,"PercChangeF_Table_mean_reshape_SPR40.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
PercChangeF_SPR40_range = subset(PercChangeF_SPR40,select=c(Genus,Species,WPP,UncertaintyRange))
PercChangeF_SPR40_range$GenusSpecies = paste(PercChangeF_SPR40_range$Genus,PercChangeF_SPR40_range$Species,sep="_")
PercChangeF_SPR40_range$Genus=NULL
PercChangeF_SPR40_range$Species=NULL
PercChangeF_Table_range_reshape_SPR40 = reshape(PercChangeF_SPR40_range,v.names="UncertaintyRange",idvar="GenusSpecies",timevar="WPP",direction="wide")
names(PercChangeF_Table_range_reshape_SPR40) = gsub("UncertaintyRange.","",names(PercChangeF_Table_range_reshape_SPR40))
PercChangeF_Table_range_reshape_SPR40 = merge(PercChangeF_Table_range_reshape_SPR40,GenusSpeciesCombineKey,all.x=TRUE,by=c("GenusSpecies"))
PercChangeF_Table_range_reshape_SPR40$GenusSpecies=NULL
write.table(PercChangeF_Table_range_reshape_SPR40,paste(PATH_highLevelSummaryPlots,"PercChangeF_Table_range_reshape_SPR40.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")

#Rebulding Times for overfished stocks at SPR30%
RebuildTimeSPR30 = subset(RunsSPR30,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,ProjectionType,timeToRebuildToSPR30))
RebuildTimeSPR30 = RebuildTimeSPR30[RebuildTimeSPR30$populationReconstructionMethod!="LIME",]    #Remove LIME from recovery time average calculation becuase it is based on biomass and I am not using the LIME runs for biomass related benchmark determination since LIME is often estimating unreasonable biomass
RebuildTimeSPR30_willNotRebuild = RebuildTimeSPR30[RebuildTimeSPR30$TimeToRebuildToSPR30==-999,]
RebuildTimeSPR30 = RebuildTimeSPR30[!is.na(RebuildTimeSPR30$timeToRebuildToSPR30) & RebuildTimeSPR30$timeToRebuildToSPR30!=-999,]
RebuildTimeSPR30_mean = aggregate.data.frame(RebuildTimeSPR30$timeToRebuildToSPR30,by=list(RebuildTimeSPR30$Genus,RebuildTimeSPR30$Species,RebuildTimeSPR30$WPP),FUN=mean)
names(RebuildTimeSPR30_mean) = c("Genus","Species","WPP","TimeToRebuid_spr30_mean")
RebuildTimeSPR30_sd = aggregate.data.frame(RebuildTimeSPR30$timeToRebuildToSPR30,by=list(RebuildTimeSPR30$Genus,RebuildTimeSPR30$Species,RebuildTimeSPR30$WPP),FUN=sd)
names(RebuildTimeSPR30_sd) = c("Genus","Species","WPP","TimeToRebuid_spr30_sd")
RebuildTimeSPR30 = merge(RebuildTimeSPR30_mean,RebuildTimeSPR30_sd,all.x=TRUE,by=c("Genus","Species","WPP"))
RebuildTimeSPR30$Upper = RebuildTimeSPR30$TimeToRebuid_spr30_mean + RebuildTimeSPR30$TimeToRebuid_spr30_sd
RebuildTimeSPR30$Lower = RebuildTimeSPR30$TimeToRebuid_spr30_mean - RebuildTimeSPR30$TimeToRebuid_spr30_sd
RebuildTimeSPR30$UncertaintyRange = paste(round(RebuildTimeSPR30$Lower,0)," to ",round(RebuildTimeSPR30$Upper,0),sep="")
GenusSpeciesAreaStatusSPR30 = subset(RunsSPR30,select=c(Genus,Species,WPP))
GenusSpeciesAreaStatusSPR30 = unique(GenusSpeciesAreaStatusSPR30)
RebuildTimeSPR30_mean = subset(RebuildTimeSPR30,select=c(Genus,Species,WPP,TimeToRebuid_spr30_mean))
RebuildTimeSPR30_mean$TimeToRebuid_spr30_mean = round(RebuildTimeSPR30_mean$TimeToRebuid_spr30_mean,0)
RebuildTimeSPR30_mean$GenusSpecies = paste(RebuildTimeSPR30_mean$Genus,RebuildTimeSPR30_mean$Species,sep="_")
RebuildTimeSPR30_mean$Genus=NULL
RebuildTimeSPR30_mean$Species=NULL
RebuildTimeSPR30_mean_reshape_SPR30 = reshape(RebuildTimeSPR30_mean,v.names="TimeToRebuid_spr30_mean",idvar="GenusSpecies",timevar="WPP",direction="wide")
names(RebuildTimeSPR30_mean_reshape_SPR30) = gsub("TimeToRebuid_spr30_mean.","",names(RebuildTimeSPR30_mean_reshape_SPR30))
RebuildTimeSPR30_mean_reshape_SPR30 = merge(RebuildTimeSPR30_mean_reshape_SPR30,GenusSpeciesCombineKey,all.x=TRUE,by=c("GenusSpecies"))
RebuildTimeSPR30_mean_reshape_SPR30$GenusSpecies=NULL
write.table(RebuildTimeSPR30_mean_reshape_SPR30,paste(PATH_highLevelSummaryPlots,"RebuildTimeSPR30_mean_reshape_SPR30.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
RebuildTimeSPR30_range = subset(RebuildTimeSPR30,select=c(Genus,Species,WPP,UncertaintyRange))
RebuildTimeSPR30_range$GenusSpecies = paste(RebuildTimeSPR30_range$Genus,RebuildTimeSPR30_range$Species,sep="_")
RebuildTimeSPR30_range$Genus=NULL
RebuildTimeSPR30_range$Species=NULL
RebuildTimeSPR30_range_reshape_SPR30 = reshape(RebuildTimeSPR30_range,v.names="UncertaintyRange",idvar="GenusSpecies",timevar="WPP",direction="wide")
names(RebuildTimeSPR30_range_reshape_SPR30) = gsub("UncertaintyRange.","",names(RebuildTimeSPR30_range_reshape_SPR30))
RebuildTimeSPR30_range_reshape_SPR30 = merge(RebuildTimeSPR30_range_reshape_SPR30,GenusSpeciesCombineKey,all.x=TRUE,by=c("GenusSpecies"))
RebuildTimeSPR30_range_reshape_SPR30$GenusSpecies=NULL
write.table(RebuildTimeSPR30_range_reshape_SPR30,paste(PATH_highLevelSummaryPlots,"RebuildTimeSPR30_range_reshape_SPR30.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")

#Rebulding Times for overfished stocks at SPR20%
RebuildTimeSPR20 = subset(RunsSPR20,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,ProjectionType,timeToRebuildToSPR20))
RebuildTimeSPR20 = RebuildTimeSPR20[RebuildTimeSPR20$populationReconstructionMethod!="LIME",]    #Remove LIME from recovery time average calculation becuase it is based on biomass and I am not using the LIME runs for biomass related benchmark determination since LIME is often estimating unreasonable biomass
RebuildTimeSPR20_willNotRebuild = RebuildTimeSPR20[RebuildTimeSPR20$timeToRebuildToSPR20==-999,]
RebuildTimeSPR20 = RebuildTimeSPR20[!is.na(RebuildTimeSPR20$timeToRebuildToSPR20) & RebuildTimeSPR20$timeToRebuildToSPR20!=-999,]
RebuildTimeSPR20_mean = aggregate.data.frame(RebuildTimeSPR20$timeToRebuildToSPR20,by=list(RebuildTimeSPR20$Genus,RebuildTimeSPR20$Species,RebuildTimeSPR20$WPP),FUN=mean)
names(RebuildTimeSPR20_mean) = c("Genus","Species","WPP","TimeToRebuid_SPR20_mean")
RebuildTimeSPR20_sd = aggregate.data.frame(RebuildTimeSPR20$timeToRebuildToSPR20,by=list(RebuildTimeSPR20$Genus,RebuildTimeSPR20$Species,RebuildTimeSPR20$WPP),FUN=sd)
names(RebuildTimeSPR20_sd) = c("Genus","Species","WPP","TimeToRebuid_SPR20_sd")
RebuildTimeSPR20 = merge(RebuildTimeSPR20_mean,RebuildTimeSPR20_sd,all.x=TRUE,by=c("Genus","Species","WPP"))
RebuildTimeSPR20$Upper = RebuildTimeSPR20$TimeToRebuid_SPR20_mean + RebuildTimeSPR20$TimeToRebuid_SPR20_sd
RebuildTimeSPR20$Lower = RebuildTimeSPR20$TimeToRebuid_SPR20_mean - RebuildTimeSPR20$TimeToRebuid_SPR20_sd
RebuildTimeSPR20$UncertaintyRange = paste(round(RebuildTimeSPR20$Lower,0)," to ",round(RebuildTimeSPR20$Upper,0),sep="")
GenusSpeciesAreaStatusSPR20 = subset(RunsSPR20,select=c(Genus,Species,WPP))
GenusSpeciesAreaStatusSPR20 = unique(GenusSpeciesAreaStatusSPR20)
RebuildTimeSPR20_mean = subset(RebuildTimeSPR20,select=c(Genus,Species,WPP,TimeToRebuid_SPR20_mean))
RebuildTimeSPR20_mean$TimeToRebuid_SPR20_mean = round(RebuildTimeSPR20_mean$TimeToRebuid_SPR20_mean,0)
RebuildTimeSPR20_mean$GenusSpecies = paste(RebuildTimeSPR20_mean$Genus,RebuildTimeSPR20_mean$Species,sep="_")
RebuildTimeSPR20_mean$Genus=NULL
RebuildTimeSPR20_mean$Species=NULL
RebuildTimeSPR20_mean_reshape_SPR20 = reshape(RebuildTimeSPR20_mean,v.names="TimeToRebuid_SPR20_mean",idvar="GenusSpecies",timevar="WPP",direction="wide")
names(RebuildTimeSPR20_mean_reshape_SPR20) = gsub("TimeToRebuid_SPR20_mean.","",names(RebuildTimeSPR20_mean_reshape_SPR20))
RebuildTimeSPR20_mean_reshape_SPR20 = merge(RebuildTimeSPR20_mean_reshape_SPR20,GenusSpeciesCombineKey,all.x=TRUE,by=c("GenusSpecies"))
RebuildTimeSPR20_mean_reshape_SPR20$GenusSpecies=NULL
write.table(RebuildTimeSPR20_mean_reshape_SPR20,paste(PATH_highLevelSummaryPlots,"RebuildTimeSPR20_mean_reshape_SPR20.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
RebuildTimeSPR20_range = subset(RebuildTimeSPR20,select=c(Genus,Species,WPP,UncertaintyRange))
RebuildTimeSPR20_range$GenusSpecies = paste(RebuildTimeSPR20_range$Genus,RebuildTimeSPR20_range$Species,sep="_")
RebuildTimeSPR20_range$Genus=NULL
RebuildTimeSPR20_range$Species=NULL
RebuildTimeSPR20_range_reshape_SPR20 = reshape(RebuildTimeSPR20_range,v.names="UncertaintyRange",idvar="GenusSpecies",timevar="WPP",direction="wide")
names(RebuildTimeSPR20_range_reshape_SPR20) = gsub("UncertaintyRange.","",names(RebuildTimeSPR20_range_reshape_SPR20))
RebuildTimeSPR20_range_reshape_SPR20 = merge(RebuildTimeSPR20_range_reshape_SPR20,GenusSpeciesCombineKey,all.x=TRUE,by=c("GenusSpecies"))
RebuildTimeSPR20_range_reshape_SPR20$GenusSpecies=NULL
write.table(RebuildTimeSPR20_range_reshape_SPR20,paste(PATH_highLevelSummaryPlots,"RebuildTimeSPR20_range_reshape_SPR20.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")

#Rebulding Times for overfished stocks at SPR40%
RebuildTimeSPR40 = subset(RunsSPR40,select=c(Genus,Species,WPP,CatchMSY_Scenario,CatchMSY_Scenario_Bound,LifeHistory_Scenario,populationReconstructionMethod,SteepnessAssumed,ProjectionType,timeToRebuildToSPR40))
RebuildTimeSPR40 = RebuildTimeSPR40[RebuildTimeSPR40$populationReconstructionMethod!="LIME",]    #Remove LIME from recovery time average calculation becuase it is based on biomass and I am not using the LIME runs for biomass related benchmark determination since LIME is often estimating unreasonable biomass
RebuildTimeSPR40_willNotRebuild = RebuildTimeSPR40[RebuildTimeSPR40$timeToRebuildToSPR40==-999,]
RebuildTimeSPR40 = RebuildTimeSPR40[!is.na(RebuildTimeSPR40$timeToRebuildToSPR40) & RebuildTimeSPR40$timeToRebuildToSPR40!=-999,]
RebuildTimeSPR40_mean = aggregate.data.frame(RebuildTimeSPR40$timeToRebuildToSPR40,by=list(RebuildTimeSPR40$Genus,RebuildTimeSPR40$Species,RebuildTimeSPR40$WPP),FUN=mean)
names(RebuildTimeSPR40_mean) = c("Genus","Species","WPP","TimeToRebuid_SPR40_mean")
RebuildTimeSPR40_sd = aggregate.data.frame(RebuildTimeSPR40$timeToRebuildToSPR40,by=list(RebuildTimeSPR40$Genus,RebuildTimeSPR40$Species,RebuildTimeSPR40$WPP),FUN=sd)
names(RebuildTimeSPR40_sd) = c("Genus","Species","WPP","TimeToRebuid_SPR40_sd")
RebuildTimeSPR40 = merge(RebuildTimeSPR40_mean,RebuildTimeSPR40_sd,all.x=TRUE,by=c("Genus","Species","WPP"))
RebuildTimeSPR40$Upper = RebuildTimeSPR40$TimeToRebuid_SPR40_mean + RebuildTimeSPR40$TimeToRebuid_SPR40_sd
RebuildTimeSPR40$Lower = RebuildTimeSPR40$TimeToRebuid_SPR40_mean - RebuildTimeSPR40$TimeToRebuid_SPR40_sd
RebuildTimeSPR40$UncertaintyRange = paste(round(RebuildTimeSPR40$Lower,0)," to ",round(RebuildTimeSPR40$Upper,0),sep="")
GenusSpeciesAreaStatusSPR40 = subset(RunsSPR40,select=c(Genus,Species,WPP))
GenusSpeciesAreaStatusSPR40 = unique(GenusSpeciesAreaStatusSPR40)
RebuildTimeSPR40_mean = subset(RebuildTimeSPR40,select=c(Genus,Species,WPP,TimeToRebuid_SPR40_mean))
RebuildTimeSPR40_mean$TimeToRebuid_SPR40_mean = round(RebuildTimeSPR40_mean$TimeToRebuid_SPR40_mean,0)
RebuildTimeSPR40_mean$GenusSpecies = paste(RebuildTimeSPR40_mean$Genus,RebuildTimeSPR40_mean$Species,sep="_")
RebuildTimeSPR40_mean$Genus=NULL
RebuildTimeSPR40_mean$Species=NULL
RebuildTimeSPR40_mean_reshape_SPR40 = reshape(RebuildTimeSPR40_mean,v.names="TimeToRebuid_SPR40_mean",idvar="GenusSpecies",timevar="WPP",direction="wide")
names(RebuildTimeSPR40_mean_reshape_SPR40) = gsub("TimeToRebuid_SPR40_mean.","",names(RebuildTimeSPR40_mean_reshape_SPR40))
RebuildTimeSPR40_mean_reshape_SPR40 = merge(RebuildTimeSPR40_mean_reshape_SPR40,GenusSpeciesCombineKey,all.x=TRUE,by=c("GenusSpecies"))
RebuildTimeSPR40_mean_reshape_SPR40$GenusSpecies=NULL
write.table(RebuildTimeSPR40_mean_reshape_SPR40,paste(PATH_highLevelSummaryPlots,"RebuildTimeSPR40_mean_reshape_SPR40.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")
RebuildTimeSPR40_range = subset(RebuildTimeSPR40,select=c(Genus,Species,WPP,UncertaintyRange))
RebuildTimeSPR40_range$GenusSpecies = paste(RebuildTimeSPR40_range$Genus,RebuildTimeSPR40_range$Species,sep="_")
RebuildTimeSPR40_range$Genus=NULL
RebuildTimeSPR40_range$Species=NULL
RebuildTimeSPR40_range_reshape_SPR40 = reshape(RebuildTimeSPR40_range,v.names="UncertaintyRange",idvar="GenusSpecies",timevar="WPP",direction="wide")
names(RebuildTimeSPR40_range_reshape_SPR40) = gsub("UncertaintyRange.","",names(RebuildTimeSPR40_range_reshape_SPR40))
RebuildTimeSPR40_range_reshape_SPR40 = merge(RebuildTimeSPR40_range_reshape_SPR40,GenusSpeciesCombineKey,all.x=TRUE,by=c("GenusSpecies"))
RebuildTimeSPR40_range_reshape_SPR40$GenusSpecies=NULL
write.table(RebuildTimeSPR40_range_reshape_SPR40,paste(PATH_highLevelSummaryPlots,"RebuildTimeSPR40_range_reshape_SPR40.csv",sep="/"),col.names=TRUE,row.names=FALSE,sep=",")

################# Histograms of Rebuilding Times for Paper #############################
png(paste(PATH_highLevelSummaryPlots,"RebuildingTimes_Histograms_allSpecies.png",sep="/"),units="px",width=3200,height=3200,res=600)
par(mfrow=c(2,2))
AvgRebuildingTimes_totalRebuild = read.table(paste(PATH_output,"AvgRebuildingTimes.csv",sep="/"),header=TRUE,sep=",")
AvgRebuildingTimes_totalRebuild$GenusSpecies = paste(AvgRebuildingTimes_totalRebuild$Genus,AvgRebuildingTimes_totalRebuild$Species,sep="_")
GenusSpeciesNameLink = subset(AvgRebuildingTimes_totalRebuild,select=c(Genus,Species,GenusSpecies))
GenusSpeciesNameLink = unique(GenusSpeciesNameLink)
AvgRebuildingTimes_totalRebuild$Genus=NULL
AvgRebuildingTimes_totalRebuild$Species=NULL
AvgRebuildingTimes_long = reshape(AvgRebuildingTimes_totalRebuild,varying=list(1:(dim(AvgRebuildingTimes_totalRebuild)[2]-1)),v.names="TotalRebuilding",idvar="GenusSpecies",direction="long")
AvgRebuildingTimes_long_noNA = AvgRebuildingTimes_long[!is.na(AvgRebuildingTimes_long$TotalRebuilding),]
hist(AvgRebuildingTimes_long_noNA$TotalRebuilding,xlab="Years",main="Time To Rebuild:\nTmin + One Generation",breaks=seq(2,52,by=2),freq=FALSE)
RunsSPR20 = read.table(paste(PATH_output,"RunsSPR20_withRebuildTimes.csv",sep="/"),header=TRUE,sep=",")
RunsSPR20 = subset(RunsSPR20,!is.na(RunsSPR20$timeToRebuildToSPR20) & RunsSPR20$timeToRebuildToSPR20!=-999)
hist(RunsSPR20$timeToRebuildToSPR20,xlab="Years",main="Time To Rebuild to SPR 20%",breaks=seq(0,52,by=2),freq=FALSE)
RunsSPR30 = read.table(paste(PATH_output,"RunsSPR30_withRebuildTimes.csv",sep="/"),header=TRUE,sep=",")
RunsSPR30 = subset(RunsSPR30,!is.na(RunsSPR30$timeToRebuildToSPR30) & RunsSPR30$timeToRebuildToSPR30!=-999)
hist(RunsSPR30$timeToRebuildToSPR30,xlab="Years",main="Time To Rebuild to SPR 30%",breaks=seq(0,52,by=2),freq=FALSE)
RunsSPR40 = read.table(paste(PATH_output,"RunsSPR40_withRebuildTimes.csv",sep="/"),header=TRUE,sep=",")
RunsSPR40 = subset(RunsSPR40,!is.na(RunsSPR40$timeToRebuildToSPR40) & RunsSPR40$timeToRebuildToSPR40!=-999)
hist(RunsSPR40$timeToRebuildToSPR40,xlab="Years",main="Time To Rebuild to SPR 40%",breaks=seq(0,52,by=2),freq=FALSE)
dev.off()

################# Snapper Family: Histograms of Rebuilding Times for Paper #############################
ExtrapolatedLandings = read.table(paste(PATH_output,"ExtrapolatedLandings.csv",sep="/"),header=TRUE,sep=",")
GenusSpeciesFamily = subset(ExtrapolatedLandings,select=c(Genus,Species,FamilyCommonName))
GenusSpeciesFamily = unique(GenusSpeciesFamily)
GenusSpeciesFamily$Genus = as.character(GenusSpeciesFamily$Genus)
GenusSpeciesFamily$Species = as.character(GenusSpeciesFamily$Species)
GenusSpeciesFamily$FamilyCommonName = as.character(GenusSpeciesFamily$FamilyCommonName)
GenusSpeciesFamily$GenusSpecies = paste(GenusSpeciesFamily$Genus,GenusSpeciesFamily$Species,sep="_")
row.names(GenusSpeciesFamily) = NULL
png(paste(PATH_highLevelSummaryPlots,"RebuildingTimes_Histograms_Snappers.png",sep="/"),units="px",width=3200,height=3200,res=600)
par(mfrow=c(2,2))
AvgRebuildingTimes_totalRebuild = read.table(paste(PATH_output,"AvgRebuildingTimes.csv",sep="/"),header=TRUE,sep=",")
AvgRebuildingTimes_totalRebuild$GenusSpecies = paste(AvgRebuildingTimes_totalRebuild$Genus,AvgRebuildingTimes_totalRebuild$Species,sep="_")
GenusSpeciesNameLink = subset(AvgRebuildingTimes_totalRebuild,select=c(Genus,Species,GenusSpecies))
GenusSpeciesNameLink = unique(GenusSpeciesNameLink)
AvgRebuildingTimes_totalRebuild$Genus=NULL
AvgRebuildingTimes_totalRebuild$Species=NULL
AvgRebuildingTimes_long = reshape(AvgRebuildingTimes_totalRebuild,varying=list(1:(dim(AvgRebuildingTimes_totalRebuild)[2]-1)),v.names="TotalRebuilding",idvar="GenusSpecies",direction="long")
AvgRebuildingTimes_long_noNA = AvgRebuildingTimes_long[!is.na(AvgRebuildingTimes_long$TotalRebuilding),]
row.names(AvgRebuildingTimes_long_noNA)=NULL
AvgRebuildingTimes_long_noNA = merge(AvgRebuildingTimes_long_noNA,GenusSpeciesFamily,all.x=TRUE,by=c("GenusSpecies"))
AvgRebuildingTimes_long_noNA_snapper = AvgRebuildingTimes_long_noNA[AvgRebuildingTimes_long_noNA$FamilyCommonName=="Snapper",]
hist(AvgRebuildingTimes_long_noNA_snapper$TotalRebuilding,xlab="Years",main="Time To Rebuild:\nTmin + One Generation",breaks=seq(2,52,by=2),freq=FALSE)
mtext("Snappers")
RunsSPR20 = read.table(paste(PATH_output,"RunsSPR20_withRebuildTimes.csv",sep="/"),header=TRUE,sep=",")
RunsSPR20 = merge(RunsSPR20,GenusSpeciesFamily,all.x=TRUE,by=c("Genus","Species"))
RunsSPR20 = subset(RunsSPR20,RunsSPR20$FamilyCommonName=="Snapper")
RunsSPR20 = subset(RunsSPR20,!is.na(RunsSPR20$timeToRebuildToSPR20) & RunsSPR20$timeToRebuildToSPR20!=-999)
hist(RunsSPR20$timeToRebuildToSPR20,xlab="Years",main="Time To Rebuild to SPR 20%",breaks=seq(0,52,by=2),freq=FALSE)
mtext("Snappers")
RunsSPR30 = read.table(paste(PATH_output,"RunsSPR30_withRebuildTimes.csv",sep="/"),header=TRUE,sep=",")
RunsSPR30 = merge(RunsSPR30,GenusSpeciesFamily,all.x=TRUE,by=c("Genus","Species"))
RunsSPR30 = subset(RunsSPR30,RunsSPR30$FamilyCommonName=="Snapper")
RunsSPR30 = subset(RunsSPR30,!is.na(RunsSPR30$timeToRebuildToSPR30) & RunsSPR30$timeToRebuildToSPR30!=-999)
hist(RunsSPR30$timeToRebuildToSPR30,xlab="Years",main="Time To Rebuild to SPR 30%",breaks=seq(0,52,by=2),freq=FALSE)
mtext("Snappers")
RunsSPR40 = read.table(paste(PATH_output,"RunsSPR40_withRebuildTimes.csv",sep="/"),header=TRUE,sep=",")
RunsSPR40 = merge(RunsSPR40,GenusSpeciesFamily,all.x=TRUE,by=c("Genus","Species"))
RunsSPR40 = subset(RunsSPR40,RunsSPR40$FamilyCommonName=="Snapper")
RunsSPR40 = subset(RunsSPR40,!is.na(RunsSPR40$timeToRebuildToSPR40) & RunsSPR40$timeToRebuildToSPR40!=-999)
hist(RunsSPR40$timeToRebuildToSPR40,xlab="Years",main="Time To Rebuild to SPR 40%",breaks=seq(0,52,by=2),freq=FALSE)
mtext("Snappers")
dev.off()

################# Grouper Family: Histograms of Rebuilding Times for Paper #############################
ExtrapolatedLandings = read.table(paste(PATH_output,"ExtrapolatedLandings.csv",sep="/"),header=TRUE,sep=",")
GenusSpeciesFamily = subset(ExtrapolatedLandings,select=c(Genus,Species,FamilyCommonName))
GenusSpeciesFamily = unique(GenusSpeciesFamily)
GenusSpeciesFamily$Genus = as.character(GenusSpeciesFamily$Genus)
GenusSpeciesFamily$Species = as.character(GenusSpeciesFamily$Species)
GenusSpeciesFamily$FamilyCommonName = as.character(GenusSpeciesFamily$FamilyCommonName)
GenusSpeciesFamily$GenusSpecies = paste(GenusSpeciesFamily$Genus,GenusSpeciesFamily$Species,sep="_")
row.names(GenusSpeciesFamily) = NULL
png(paste(PATH_highLevelSummaryPlots,"RebuildingTimes_Histograms_Groupers.png",sep="/"),units="px",width=3200,height=3200,res=600)
par(mfrow=c(2,2))
AvgRebuildingTimes_totalRebuild = read.table(paste(PATH_output,"AvgRebuildingTimes.csv",sep="/"),header=TRUE,sep=",")
AvgRebuildingTimes_totalRebuild$GenusSpecies = paste(AvgRebuildingTimes_totalRebuild$Genus,AvgRebuildingTimes_totalRebuild$Species,sep="_")
GenusSpeciesNameLink = subset(AvgRebuildingTimes_totalRebuild,select=c(Genus,Species,GenusSpecies))
GenusSpeciesNameLink = unique(GenusSpeciesNameLink)
AvgRebuildingTimes_totalRebuild$Genus=NULL
AvgRebuildingTimes_totalRebuild$Species=NULL
AvgRebuildingTimes_long = reshape(AvgRebuildingTimes_totalRebuild,varying=list(1:(dim(AvgRebuildingTimes_totalRebuild)[2]-1)),v.names="TotalRebuilding",idvar="GenusSpecies",direction="long")
AvgRebuildingTimes_long_noNA = AvgRebuildingTimes_long[!is.na(AvgRebuildingTimes_long$TotalRebuilding),]
row.names(AvgRebuildingTimes_long_noNA)=NULL
AvgRebuildingTimes_long_noNA = merge(AvgRebuildingTimes_long_noNA,GenusSpeciesFamily,all.x=TRUE,by=c("GenusSpecies"))
AvgRebuildingTimes_long_noNA_grouper = AvgRebuildingTimes_long_noNA[AvgRebuildingTimes_long_noNA$FamilyCommonName=="Grouper",]
hist(AvgRebuildingTimes_long_noNA_grouper$TotalRebuilding,xlab="Years",main="Time To Rebuild:\nTmin + One Generation",breaks=seq(2,52,by=2),freq=FALSE)
mtext("Groupers")
RunsSPR20 = read.table(paste(PATH_output,"RunsSPR20_withRebuildTimes.csv",sep="/"),header=TRUE,sep=",")
RunsSPR20 = merge(RunsSPR20,GenusSpeciesFamily,all.x=TRUE,by=c("Genus","Species"))
RunsSPR20 = subset(RunsSPR20,RunsSPR20$FamilyCommonName=="Grouper")
RunsSPR20 = subset(RunsSPR20,!is.na(RunsSPR20$timeToRebuildToSPR20) & RunsSPR20$timeToRebuildToSPR20!=-999)
hist(RunsSPR20$timeToRebuildToSPR20,xlab="Years",main="Time To Rebuild to SPR 20%",breaks=seq(0,52,by=2),freq=FALSE)
mtext("Groupers")
RunsSPR30 = read.table(paste(PATH_output,"RunsSPR30_withRebuildTimes.csv",sep="/"),header=TRUE,sep=",")
RunsSPR30 = merge(RunsSPR30,GenusSpeciesFamily,all.x=TRUE,by=c("Genus","Species"))
RunsSPR30 = subset(RunsSPR30,RunsSPR30$FamilyCommonName=="Grouper")
RunsSPR30 = subset(RunsSPR30,!is.na(RunsSPR30$timeToRebuildToSPR30) & RunsSPR30$timeToRebuildToSPR30!=-999)
hist(RunsSPR30$timeToRebuildToSPR30,xlab="Years",main="Time To Rebuild to SPR 30%",breaks=seq(0,52,by=2),freq=FALSE)
mtext("Groupers")
RunsSPR40 = read.table(paste(PATH_output,"RunsSPR40_withRebuildTimes.csv",sep="/"),header=TRUE,sep=",")
RunsSPR40 = merge(RunsSPR40,GenusSpeciesFamily,all.x=TRUE,by=c("Genus","Species"))
RunsSPR40 = subset(RunsSPR40,RunsSPR40$FamilyCommonName=="Grouper")
RunsSPR40 = subset(RunsSPR40,!is.na(RunsSPR40$timeToRebuildToSPR40) & RunsSPR40$timeToRebuildToSPR40!=-999)
hist(RunsSPR40$timeToRebuildToSPR40,xlab="Years",main="Time To Rebuild to SPR 40%",breaks=seq(0,52,by=2),freq=FALSE)
mtext("Groupers")
dev.off()

################# Jacks_Mackerels Family: Histograms of Rebuilding Times for Paper #############################
ExtrapolatedLandings = read.table(paste(PATH_output,"ExtrapolatedLandings.csv",sep="/"),header=TRUE,sep=",")
GenusSpeciesFamily = subset(ExtrapolatedLandings,select=c(Genus,Species,FamilyCommonName))
GenusSpeciesFamily = unique(GenusSpeciesFamily)
GenusSpeciesFamily$Genus = as.character(GenusSpeciesFamily$Genus)
GenusSpeciesFamily$Species = as.character(GenusSpeciesFamily$Species)
GenusSpeciesFamily$FamilyCommonName = as.character(GenusSpeciesFamily$FamilyCommonName)
GenusSpeciesFamily$GenusSpecies = paste(GenusSpeciesFamily$Genus,GenusSpeciesFamily$Species,sep="_")
row.names(GenusSpeciesFamily) = NULL
png(paste(PATH_highLevelSummaryPlots,"RebuildingTimes_Histograms_JacksMackerels.png",sep="/"),units="px",width=3200,height=3200,res=600)
par(mfrow=c(2,2))
AvgRebuildingTimes_totalRebuild = read.table(paste(PATH_output,"AvgRebuildingTimes.csv",sep="/"),header=TRUE,sep=",")
AvgRebuildingTimes_totalRebuild$GenusSpecies = paste(AvgRebuildingTimes_totalRebuild$Genus,AvgRebuildingTimes_totalRebuild$Species,sep="_")
GenusSpeciesNameLink = subset(AvgRebuildingTimes_totalRebuild,select=c(Genus,Species,GenusSpecies))
GenusSpeciesNameLink = unique(GenusSpeciesNameLink)
AvgRebuildingTimes_totalRebuild$Genus=NULL
AvgRebuildingTimes_totalRebuild$Species=NULL
AvgRebuildingTimes_long = reshape(AvgRebuildingTimes_totalRebuild,varying=list(1:(dim(AvgRebuildingTimes_totalRebuild)[2]-1)),v.names="TotalRebuilding",idvar="GenusSpecies",direction="long")
AvgRebuildingTimes_long_noNA = AvgRebuildingTimes_long[!is.na(AvgRebuildingTimes_long$TotalRebuilding),]
row.names(AvgRebuildingTimes_long_noNA)=NULL
AvgRebuildingTimes_long_noNA = merge(AvgRebuildingTimes_long_noNA,GenusSpeciesFamily,all.x=TRUE,by=c("GenusSpecies"))
AvgRebuildingTimes_long_noNA_jackMackerel = AvgRebuildingTimes_long_noNA[AvgRebuildingTimes_long_noNA$FamilyCommonName=="Jacks_Mackerels",]
hist(AvgRebuildingTimes_long_noNA_jackMackerel$TotalRebuilding,xlab="Years",main="Time To Rebuild:\nTmin + One Generation",breaks=seq(2,52,by=2),freq=FALSE)
mtext("Jacks/Mackerels")
RunsSPR20 = read.table(paste(PATH_output,"RunsSPR20_withRebuildTimes.csv",sep="/"),header=TRUE,sep=",")
RunsSPR20 = merge(RunsSPR20,GenusSpeciesFamily,all.x=TRUE,by=c("Genus","Species"))
RunsSPR20 = subset(RunsSPR20,RunsSPR20$FamilyCommonName=="Jacks_Mackerels")
RunsSPR20 = subset(RunsSPR20,!is.na(RunsSPR20$timeToRebuildToSPR20) & RunsSPR20$timeToRebuildToSPR20!=-999)
hist(RunsSPR20$timeToRebuildToSPR20,xlab="Years",main="Time To Rebuild to SPR 20%",breaks=seq(0,52,by=2),freq=FALSE)
mtext("Jacks/Mackerels")
RunsSPR30 = read.table(paste(PATH_output,"RunsSPR30_withRebuildTimes.csv",sep="/"),header=TRUE,sep=",")
RunsSPR30 = merge(RunsSPR30,GenusSpeciesFamily,all.x=TRUE,by=c("Genus","Species"))
RunsSPR30 = subset(RunsSPR30,RunsSPR30$FamilyCommonName=="Jacks_Mackerels")
RunsSPR30 = subset(RunsSPR30,!is.na(RunsSPR30$timeToRebuildToSPR30) & RunsSPR30$timeToRebuildToSPR30!=-999)
hist(RunsSPR30$timeToRebuildToSPR30,xlab="Years",main="Time To Rebuild to SPR 30%",breaks=seq(0,52,by=2),freq=FALSE)
mtext("Jacks/Mackerels")
RunsSPR40 = read.table(paste(PATH_output,"RunsSPR40_withRebuildTimes.csv",sep="/"),header=TRUE,sep=",")
RunsSPR40 = merge(RunsSPR40,GenusSpeciesFamily,all.x=TRUE,by=c("Genus","Species"))
RunsSPR40 = subset(RunsSPR40,RunsSPR40$FamilyCommonName=="Jacks_Mackerels")
RunsSPR40 = subset(RunsSPR40,!is.na(RunsSPR40$timeToRebuildToSPR40) & RunsSPR40$timeToRebuildToSPR40!=-999)
hist(RunsSPR40$timeToRebuildToSPR40,xlab="Years",main="Time To Rebuild to SPR 40%",breaks=seq(0,52,by=2),freq=FALSE)
mtext("Jacks/Mackerels")
dev.off()

################# Emperor Family: Histograms of Rebuilding Times for Paper #############################
ExtrapolatedLandings = read.table(paste(PATH_output,"ExtrapolatedLandings.csv",sep="/"),header=TRUE,sep=",")
GenusSpeciesFamily = subset(ExtrapolatedLandings,select=c(Genus,Species,FamilyCommonName))
GenusSpeciesFamily = unique(GenusSpeciesFamily)
GenusSpeciesFamily$Genus = as.character(GenusSpeciesFamily$Genus)
GenusSpeciesFamily$Species = as.character(GenusSpeciesFamily$Species)
GenusSpeciesFamily$FamilyCommonName = as.character(GenusSpeciesFamily$FamilyCommonName)
GenusSpeciesFamily$GenusSpecies = paste(GenusSpeciesFamily$Genus,GenusSpeciesFamily$Species,sep="_")
row.names(GenusSpeciesFamily) = NULL
png(paste(PATH_highLevelSummaryPlots,"RebuildingTimes_Histograms_Emperor.png",sep="/"),units="px",width=3200,height=3200,res=600)
par(mfrow=c(2,2))
AvgRebuildingTimes_totalRebuild = read.table(paste(PATH_output,"AvgRebuildingTimes.csv",sep="/"),header=TRUE,sep=",")
AvgRebuildingTimes_totalRebuild$GenusSpecies = paste(AvgRebuildingTimes_totalRebuild$Genus,AvgRebuildingTimes_totalRebuild$Species,sep="_")
GenusSpeciesNameLink = subset(AvgRebuildingTimes_totalRebuild,select=c(Genus,Species,GenusSpecies))
GenusSpeciesNameLink = unique(GenusSpeciesNameLink)
AvgRebuildingTimes_totalRebuild$Genus=NULL
AvgRebuildingTimes_totalRebuild$Species=NULL
AvgRebuildingTimes_long = reshape(AvgRebuildingTimes_totalRebuild,varying=list(1:(dim(AvgRebuildingTimes_totalRebuild)[2]-1)),v.names="TotalRebuilding",idvar="GenusSpecies",direction="long")
AvgRebuildingTimes_long_noNA = AvgRebuildingTimes_long[!is.na(AvgRebuildingTimes_long$TotalRebuilding),]
row.names(AvgRebuildingTimes_long_noNA)=NULL
AvgRebuildingTimes_long_noNA = merge(AvgRebuildingTimes_long_noNA,GenusSpeciesFamily,all.x=TRUE,by=c("GenusSpecies"))
AvgRebuildingTimes_long_noNA_emperor = AvgRebuildingTimes_long_noNA[AvgRebuildingTimes_long_noNA$FamilyCommonName=="Emperor",]
hist(AvgRebuildingTimes_long_noNA_emperor$TotalRebuilding,xlab="Years",main="Time To Rebuild:\nTmin + One Generation",breaks=seq(2,52,by=2),freq=FALSE)
mtext("Emperor")
RunsSPR20 = read.table(paste(PATH_output,"RunsSPR20_withRebuildTimes.csv",sep="/"),header=TRUE,sep=",")
RunsSPR20 = merge(RunsSPR20,GenusSpeciesFamily,all.x=TRUE,by=c("Genus","Species"))
RunsSPR20 = subset(RunsSPR20,RunsSPR20$FamilyCommonName=="Jacks_Mackerels")
RunsSPR20 = subset(RunsSPR20,!is.na(RunsSPR20$timeToRebuildToSPR20) & RunsSPR20$timeToRebuildToSPR20!=-999)
hist(RunsSPR20$timeToRebuildToSPR20,xlab="Years",main="Time To Rebuild to SPR 20%",breaks=seq(0,52,by=2),freq=FALSE)
mtext("Emperor")
RunsSPR30 = read.table(paste(PATH_output,"RunsSPR30_withRebuildTimes.csv",sep="/"),header=TRUE,sep=",")
RunsSPR30 = merge(RunsSPR30,GenusSpeciesFamily,all.x=TRUE,by=c("Genus","Species"))
RunsSPR30 = subset(RunsSPR30,RunsSPR30$FamilyCommonName=="Jacks_Mackerels")
RunsSPR30 = subset(RunsSPR30,!is.na(RunsSPR30$timeToRebuildToSPR30) & RunsSPR30$timeToRebuildToSPR30!=-999)
hist(RunsSPR30$timeToRebuildToSPR30,xlab="Years",main="Time To Rebuild to SPR 30%",breaks=seq(0,52,by=2),freq=FALSE)
mtext("Emperor")
RunsSPR40 = read.table(paste(PATH_output,"RunsSPR40_withRebuildTimes.csv",sep="/"),header=TRUE,sep=",")
RunsSPR40 = merge(RunsSPR40,GenusSpeciesFamily,all.x=TRUE,by=c("Genus","Species"))
RunsSPR40 = subset(RunsSPR40,RunsSPR40$FamilyCommonName=="Jacks_Mackerels")
RunsSPR40 = subset(RunsSPR40,!is.na(RunsSPR40$timeToRebuildToSPR40) & RunsSPR40$timeToRebuildToSPR40!=-999)
hist(RunsSPR40$timeToRebuildToSPR40,xlab="Years",main="Time To Rebuild to SPR 40%",breaks=seq(0,52,by=2),freq=FALSE)
mtext("Emperor")
dev.off()





####################################
PATH_mergeLifeHistVals = "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/BiologicalModel/Output/Run50Species"
iFishParamsUsed = read.csv(paste(PATH_mergeLifeHistVals,"LifeHistoryParamsUsed.csv",sep="/"),sep=",",header=TRUE)
SummaryTableFittedValues = read.csv(paste(PATH_mergeLifeHistVals,"SummaryTableFittedValues.csv",sep="/"),sep=",",header=TRUE)
iFishParamsUsed = subset(iFishParamsUsed,iFishParamsUsed$fish_genus!="")
iFishParamsUsed = subset(iFishParamsUsed,iFishParamsUsed$fish_species!="")
iFishParamsUsed = subset(iFishParamsUsed,select=c(fish_genus,fish_species,lmat,lopt,linf,lmax,lmatm,var_a,var_b))
iFishParamsUsed = unique(iFishParamsUsed)  
SummaryTableFittedValues = subset(SummaryTableFittedValues,SummaryTableFittedValues$CatchMSY_Bound=="mean")
SummaryTableFittedValues_selected = subset(SummaryTableFittedValues,select=c(Genus,Species,WPP,Selex_Logistic_a_length,Selex_Logistic_b_length,Selex_Logistic_a_age,Selex_Logistic_b_age,Mnat))
SummaryTableFittedValues_selected = unique(SummaryTableFittedValues_selected)
SummaryTableFvals = subset(SummaryTableFittedValues,select=c(Genus,Species,WPP,populationReconstructionMethod,F_estimated))
SummaryTableFvals = unique(SummaryTableFvals)
i = sapply(SummaryTableFvals,is.factor)
SummaryTableFvals[i] = lapply(SummaryTableFvals[i],as.character)
SummaryTableFvals_Identifier = subset(SummaryTableFvals,select=c(Genus,Species,WPP))
SummaryTableFvals_Identifier = unique(SummaryTableFvals_Identifier)
SummaryTableFvals_Identifier$Identifier = paste(SummaryTableFvals_Identifier$Genus,SummaryTableFvals_Identifier$Species,SummaryTableFvals_Identifier$WPP,sep="_")
SummaryTableFvals_Identifier = unique(SummaryTableFvals_Identifier)
SummaryTableFvals = merge(SummaryTableFvals,SummaryTableFvals_Identifier,all.x=TRUE,by=c("Genus","Species","WPP"))
SummaryTableFvals$Genus=NULL
SummaryTableFvals$Species=NULL
SummaryTableFvals$WPP=NULL
SummaryTableFvals$populationReconstructionMethod = as.character(SummaryTableFvals$populationReconstructionMethod)
SummaryTableFvals_reshape = reshape(SummaryTableFvals,v.names="F_estimated",idvar="Identifier",timevar="populationReconstructionMethod",direction="wide")
Identifiers = lapply(X=SummaryTableFvals_reshape$Identifier,FUN=function(x) {unlist(strsplit(as.character(x),"_"))})
Identifiers = do.call(rbind.data.frame,Identifiers)
names(Identifiers) = c("Genus","Species","WPP")
i = sapply(Identifiers,is.factor)
Identifiers[i] = lapply(Identifiers[i],as.character)
SummaryTableFvals_reshape = cbind(SummaryTableFvals_reshape,Identifiers)
SummaryTableFvals_reshape$Identifier=NULL
SummaryTableFittedValues_selected$Identifier=NULL
SummaryTableFittedValues_selected = merge(SummaryTableFittedValues_selected,SummaryTableFvals_reshape,all.x=TRUE,by=c("Genus","Species","WPP"))
names(iFishParamsUsed)[names(iFishParamsUsed)=="fish_genus"]="Genus"
names(iFishParamsUsed)[names(iFishParamsUsed)=="fish_species"]="Species"
ParametersUsedTable = merge(SummaryTableFittedValues_selected,iFishParamsUsed,all.x=TRUE,by=c("Genus","Species"))
ParametersUsedTable$F_estimated.SolveBaranov=NULL
ParametersUsedTable$F_estimated.TropFishR=NULL
write.table(ParametersUsedTable,paste(PATH_mergeLifeHistVals,"ParametersUsedTable.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)



################### WRITE OUT A TABLE WITH THE LIFE HISTORY VALUES AND OTHER PARAMETERS USED ##########################################
PATH_Top50 = "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/BiologicalModel/Output/Run50Species"
PATH_iFish = "F:/Dropbox/KeyMarineConsulting/OceansConservancy/Indonesia/BiologicalModel/iFish"
ListOfPopulations = readRDS(file=paste(PATH_Top50,"ListOfPopulations.RData",sep="/"),refhook=NULL)
SummaryTableFittedValues = read.table(paste(PATH_Top50,"SummaryTableFittedValues.csv",sep="/"),header=TRUE,sep=",")
SummaryTableFittedValues = subset(SummaryTableFittedValues,SummaryTableFittedValues$populationReconstructionMethod!="TropFishR")
R0_avg = aggregate.data.frame(SummaryTableFittedValues$R0,by=list(SummaryTableFittedValues$Genus,SummaryTableFittedValues$Species,SummaryTableFittedValues$WPP),FUN=mean,na.rm=TRUE)
names(R0_avg) = c("Genus","Species","WPP","R0")
SpeciesParams = subset(SummaryTableFittedValues,select=c(Genus,Species,WPP,Selex_Logistic_a_length,Selex_Logistic_b_length,Selex_Logistic_a_age,Selex_Logistic_b_age,fish_id))
SpeciesParams = unique(SpeciesParams)
SpeciesParams = merge(SpeciesParams,R0_avg,by=c("Genus","Species","WPP"),all.x=TRUE)
dim(SpeciesParams)
Species = read.table(paste(PATH_iFish,"ifish_fish.txt",sep="/"),sep="\t",header=TRUE,as.is=TRUE,colClasses=c("character"),stringsAsFactors=FALSE)
Species = subset(Species,Species$fish_species!="")
Species = subset(Species,select=c(fish_genus,fish_species,lmat,lopt,linf,lmax,var_a,var_b))
Species = unique(Species)
GenusSpecies = subset(SpeciesParams,select=c(Genus,Species))
GenusSpecies = unique(GenusSpecies)
names(Species)[names(Species)=="fish_genus"]="Genus"
names(Species)[names(Species)=="fish_species"]="Species"
GenusSpecies = merge(GenusSpecies,Species,all.x=TRUE,by=c("Genus","Species"))
write.table(GenusSpecies,paste(PATH_Top50,"PARAM_IN_PAPER.csv",sep="/"),sep=",",col.names=TRUE,row.names=FALSE)








