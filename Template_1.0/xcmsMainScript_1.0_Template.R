#XCMS Untargeted Analysis
#2023-05-11
#Analysis for RZcheckerboard in NEG mode
#For 2023-03-03 experiment set

BiocManager::install("xcms")
install.packages("RColorBrewer")

####PREP####

#Load your libraries
library(xcms)
library(manhattanly)
library(htmlwidgets)
library(MatrixGenerics)
library(RColorBrewer)
register(SerialParam())  #Prevents xcms errors on PCs. Other options will run faster on macs or clusters


#Set your working directory
setwd("C:/Users/RheeLab/Documents/R Pipelines/2023-05 XCMS Demo")

#Read in my custom functions
source("ES_XCMSfunctions_1.0.R")

#Name your experiment
exptkey = "RIFpH2_POS"  #Choose a short identifying tag for output files and data names. Don't start it with a number

#To load a previous session into R, run the command below.
#load(paste0(exptkey,"_xcms01.RData"))



####SETTINGS####
#Choose the options that match your instrument and method

#Load your experiment description tables
pdtable = read.csv("pdtable.csv", stringsAsFactors = FALSE)
condkey = read.csv("conditionkey.csv", header=TRUE, stringsAsFactors=FALSE)
tTestKey = read.csv("comparekey.csv", header=TRUE, stringsAsFactors=FALSE)

#Set ppm error and mode
seterr = 50  #for matching peak groups with molecules
setmode = "POS"  #Put NEG or POS here depending on mode.

#Retention time range:

#For ANP methods
rtrange = c(30, 1080)  #Minimum and maximum retention times to analyze, in seconds.

#For HILIC
#rtrange = c(30, 1440)  #Minimum and maximum retention times to analyze, in seconds.

#Peak finding parameters:

#For TOF3 and TOF2
# cwp <- CentWaveParam(
#   ppm = 25,                 #error for m/z slice within each file within each peak. 25 is generous
#   peakwidth = c(10, 90),    #min and max expected peak width, in seconds
#   noise = 1500,             #remove all reads with intensity below this threshold
#   snthresh = 10,            #signal-to-noise ratio cutoff. Not really sure what it does.
#   prefilter = c(3, 2500))   #Require at least x reads with intensity >= y to look for peaks in a m/z slice

# For QTOF2
cwp <- CentWaveParam(
  ppm = 25,                 #ppm error within each file within each peak. 25 is generous
  peakwidth = c(8, 30),    #min and max expected peak width, in seconds
  noise = 2000,             #remove all reads with intensity below this threshold
  snthresh = 10,            #signal-to-noise ratio cutoff. Not really sure what it does.
  prefilter = c(3, 3300))   #Require at least x reads with intensity >= y to look in a mass slice for peaks

#Peak merging parameters:
mpp <- MergeNeighboringPeaksParam(
  expandRt = 4,       #how close, in seconds, do two peaks need to be to count as overlapping?
  ppm = 10,           #how close, in m/z ppm, do the peaks need to be?
  minProp = 0.50)     #The lowest point between two peaks must be over this fraction of the lower peak's height to merge

#Retention time alignment parameters:

#If using QC samples
# owp <- ObiwarpParam(
#   binSize = 0.6,
#   subset = grep("qc", pdtable[,"sample_name"]), #where the QC samples are
#   subsetAdjust = "average") 

#If single run
owp <- ObiwarpParam(binSize = 0.6)

#Feature grouping parameters:
pdp <- PeakDensityParam(
  sampleGroups = pdtable[,"sample_group"],  
  minFraction = 0.8,   #Minimum fraction of samples in at least one sample group that have a peak 
  bw = 8,              #Adjusts the allowed retention time difference between peaks to group them in the same feature. Don't raise too high!
  binSize = 0.01)      #Adjusts the m/z difference allowed between peaks to get them grouped in the same feature. I'm not totally sure how it works.


#Peak filling parameters:
cpap = ChromPeakAreaParam()  #In short, using the defaults



####DATA READ-IN####

pdextra <- pdtable
lochold = lapply(condkey[,"label"], grep, x=pdtable[,"sample_name"])
for (i in c(1:length(lochold)) ){
  pdextra[lochold[[i]],"group_abbr"] = condkey[i,"label"]
  pdextra[lochold[[i]],"sample_color"] = condkey[i,"color"]
}
rownames(pdextra) = pdextra[,"sample_name"]

cdfs = paste0("./rawdata/", pdtable[,"file"])

raw_data <- readMSData(files = cdfs, 
                       pdata = new("NAnnotatedDataFrame", pdtable[,c("sample_name","sample_group")]),
                       mode = "onDisk") #2 minutes for 40 files

raw_data <- filterRt(raw_data, rtrange)  #Filter by retention time


####QUALITY CONTROL PART 1####

qcmols = read.csv("QCmols_ANP-POS.csv", header=TRUE, stringsAsFactors=FALSE)  #Choose your molecule set as appropriate for your method and mode
checkedconds = c("A", "E", "F", "J")  #Choose the subset of conditions to check
mzerr = 50  #Choose the chromatogram mz error in ppm 
RTplotrange = 180  #choose the time around the reference RT to plot, or enter "ALL"


qcmols = addmzrange(qcmols,setmode)
conddex = which(pdextra[,"group_abbr"] %in% checkedconds)
miniraw = filterFile(raw_data, conddex)

chr_raws <- chromatogram(miniraw, mz = as.matrix(qcmols[, c("mzmin","mzmax")]))  #5 minutes or so
chrs_wpeaks <- findChromPeaks(chr_raws, param = cwp)
chrs_merged <- refineChromPeaks(chrs_wpeaks, mpp)  

ChromsByConds(
  chrtoplot = chrs_merged,
  graphconds = checkedconds,
  RTrange = RTplotrange,
  molinfo = qcmols
)


####XCMS PROCESSES####
print(paste0("XCMS processes started at ", Sys.time()))

xdata <- findChromPeaks(raw_data, param = cwp)
print(paste0("Peak finding finished at ", Sys.time()))
#1hr8minutes for 40 samples

xdata_mrg <- refineChromPeaks(xdata, mpp)  
print(paste0("Peak merging finished at ", Sys.time()))
#1hr22 minutes for 40 samples

xdata_mrg <- adjustRtime(xdata_mrg, param = owp)
print(paste0("RT adjustment finished at ", Sys.time()))
#6 minutes for 40 samples.

xdata_ftg <- groupChromPeaks(xdata_mrg, param = pdp) 
print(paste0("Feature group finished at ", Sys.time()))
#3 minutes for 40 samples.

xdata_ftg <- fillChromPeaks(xdata_ftg, param = cpap)
print(paste0("Peak fill finished at ", Sys.time()))
#29 minutes for 30 samples

#save your progress
save.image(file= paste0(exptkey,"_xcms01.RData"))


####QUALITY CONTROL PART 2####

#View retention time adjustment
plotAdjustedRtime(xdata_mrg, col = pdextra[,"sample_color"])

#View retention time adjustment colored by batch
#Sorry this isn't automated - this is for my 3-batch run with 6 qc samples.
# qclocs = grep("qc", pdtable[,"sample_name"])
# batchcols = c(rep("red", qclocs[2]), rep("green3", qclocs[4]-qclocs[2]), rep("blue", qclocs[6]-qclocs[4]))
# plotAdjustedRtime(xdata_mrg, col = batchcols)




####UNTARGETED RESULTS####

#Save the table of all identified features
featdefs = featureDefinitions(xdata_ftg)
featvals = featureValues(xdata_ftg)
colnames(featvals) = paste0(exptkey, "_", pdextra[,"sample_name"])
allexport = cbind(rownames(featvals),subset(featdefs, select = c(mzmed, mzmin, mzmax, rtmed, rtmin, rtmax, npeaks)), featvals)
colnames(allexport)[1] = paste0(exptkey, "_FeatureID")
rownames(allexport)=rownames(featvals)
write.table(allexport, file = paste0(exptkey,"_AllFeatures.txt"), sep="\t", row.names = FALSE)

#To volcano plot all feature groups, with no formula info.

dir.create("./AllFeatVolcanos")
setwd("./AllFeatVolcanos")

volctest_all = MultiVolcano(dataset = xdata_ftg, 
                         comparekey = tTestKey, 
                         pthresh = 0.05,   #p-value threshold with BH multiple hypothesis correction
                         l2fcthresh = 1)   #log2 fold change threshold



#Save the feature information for all features that are significant in any one of your test categories.
colnames(volctest_all$ttests) = paste0("pval_", colnames(volctest_all$ttests))
colnames(volctest_all$l2fcs) = paste0("l2fc_", colnames(volctest_all$l2fcs))
broadfeat = sort(unique(unlist(volctest_all$sigfts)))
unfilt_sigfeatinfo = cbind(allexport[broadfeat,], 
                           volctest_all$ttests[broadfeat,],
                           volctest_all$l2fcs[broadfeat,])
write.table(unfilt_sigfeatinfo, file = paste0(exptkey,"_SigUnfiltered.txt"), sep="\t", row.names = FALSE)


#Print chromatograms of the top 50 features with highest foldchange in each test category

chromset_unfilt = vector(mode = "list", length = length(volctest_all$sigfts))
names(chromset_unfilt) = names(volctest_all$sigfts)

for (i in c(1:length(volctest_all$sigfts))) {
  dir.create(paste0("./", names(volctest_all$sigfts)[i]))
  setwd(paste0("./", names(volctest_all$sigfts)[i]))
  
  if (length(volctest_all$sigfts[[i]])>50){
    chromset_unfilt[[i]] <- featureChromatograms(xdata_ftg, features = volctest_all$sigfts[[i]][1:50], expandRt=60)
  } else {
  chromset_unfilt[[i]] <- featureChromatograms(xdata_ftg, features = volctest_all$sigfts[[i]], expandRt=60)
  }
  #Note that featureChromatograms() is very slow, especially if you have a lot of samples.
  
  ChromPrint(featchrom = chromset_unfilt[[i]], 
             condkey = condkey, 
             pdextra = pdextra)
  setwd("..")
  print(paste0("Chromatograms for ", names(volctest_all$sigfts)[i], " finished at ", Sys.time()))
}

setwd("..")


####FORMULA-TARGETED RESULTS####
#Attempt to match all features to a long list of formulae from KEGG, likely without retention times
#Then take successfully matched features and perform volcano thresholding
#The filtering helps with multiple hypothesis testing. 

mymols = read.table("PwaySet_ESpady01.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")

featfilter = FeatMatching(featdefs = featdefs,   #The feature definitions from xdata_ftg
                           molset = mymols,      #Formulae with masses to match
                           mode = setmode,    #Mode, either "POS" or "NEG"
                           mzerr = seterr,    #Error for m/z match, Default 50 ppm
                           rterr = 120)       #Error for retention time. Not used yet
#Takes 6 minutes or so

#Save matched molecules to a .txt file
matchexport = cbind(allexport[rownames(featfilter),], featfilter[,c("FormulaMatch", "ExactMass", "RT_known", "Label")])
write.table(matchexport, file = paste0(exptkey,"_FormMatchedFeatures.txt"), sep="\t", row.names = FALSE)

#Now we redo the volcano thresholding with only features that had a m/z we could 
#possibly identify

dir.create("./PwaySetVolcanos")
setwd("./PwaySetVolcanos")

volctest_PwaySet = MultiVolcano(dataset = xdata_ftg,            #Results as XCMSnExp class
                             ftsubset = rownames(featfilter),   #The features you matched
                             molmatchinfo = featfilter,         #Table with feature names and match candidates
                             comparekey = tTestKey,             #The comparison set from your csv
                             pthresh = 0.05,                    #p-value significance threshold
                             l2fcthresh = 1)                    #log2 fold change significance threshold

#Save the feature information for all features that are significant in any one of your test categories.
colnames(volctest_PwaySet$ttests) = paste0("pval_", colnames(volctest_PwaySet$ttests))
colnames(volctest_PwaySet$l2fcs) = paste0("l2fc_", colnames(volctest_PwaySet$l2fcs))
broadfeat = sort(unique(unlist(volctest_PwaySet$sigfts)))
unfilt_sigfeatinfo = cbind(allexport[broadfeat,], 
                           volctest_PwaySet$ttests[broadfeat,],
                           volctest_PwaySet$l2fcs[broadfeat,])
write.table(unfilt_sigfeatinfo, file = paste0(exptkey,"_SigPwaySet.txt"), sep="\t", row.names = FALSE)


#Print chromatograms of the top 50 features with highest foldchange in each test category

chromset_PwaySet = vector(mode = "list", length = length(volctest_PwaySet$sigfts))
names(chromset_PwaySet) = names(volctest_PwaySet$sigfts)

for (i in c(1:length(volctest_PwaySet$sigfts))) {
  dir.create(paste0("./", names(volctest_PwaySet$sigfts)[i]))
  setwd(paste0("./", names(volctest_PwaySet$sigfts)[i]))
  
  if (length(volctest_PwaySet$sigfts[[i]])>50){
    chromset_PwaySet[[i]] <- featureChromatograms(xdata_ftg, features = volctest_PwaySet$sigfts[[i]][1:50], expandRt=60)
  } else {
    chromset_PwaySet[[i]] <- featureChromatograms(xdata_ftg, features = volctest_PwaySet$sigfts[[i]], expandRt=60)
  }
  #Note that featureChromatograms() is very slow, especially if you have a lot of samples.
  
  ChromPrint(featchrom = chromset_PwaySet[[i]], 
             condkey = condkey, 
             pdextra = pdextra,
             molinfo = featfilter)
  setwd("..")
  print(paste0("Chromatograms for ", names(volctest_PwaySet$sigfts)[i], " finished at ", Sys.time()))
}

setwd("..")

####CUSTOM CHROMATOGRAMS####

#Make a chromatogram object for a set of features of interest

customfeats = c("FT1355", "FT0809", "FT1024")  #Put the features you want here

customftchroms = featureChromatograms(xdata_ftg,                #Dataset with the features 
                                      features = customfeats,   #your feature set 
                                      expandRt=60)              #time range to plot around the peaks

#Print chromatograms split into panels by chartgroup (from conditionkey.csv)

ChromPrint(featchrom = customftchroms, 
           condkey = condkey, 
           pdextra = pdextra,
           molinfo = featfilter)


#Make a chromatogram object for a set of (m/z, RT) ranges of interest.

EICinfo = read.csv("customEICinfo.csv", header=TRUE, stringsAsFactors=FALSE)  #Edit this csv for your ranges of interest

customEICs <- chromatogram(xdata_ftg,                                        #Your dataset
                           mz = as.matrix(EICinfo[, c("mzmin","mzmax")]),    #mzs pulled from csv
                           rt = as.matrix(EICinfo[, c("rtmin","rtmax")]))    #rts pulled from csv

#Print chromatograms split into panels by chartgroup (from conditionkey.csv)
#Sorry, under construction! Will need to write a new plotting loop function



###SAVE SYSTEM IMAGE####

#I recommend saving after slow processing steps
save.image(file= paste0(exptkey,"_xcms01.RData"))
