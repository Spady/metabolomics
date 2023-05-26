#2023-05-04 Functions for XCMS Targeted Analysis

#For streamlining the analysis pipeline as it stands

MplusHrange = function(extmass, ppmerr){
  acterr = (extmass+1.0073)*(ppmerr*1e-6)
  lowmass = extmass+1.0073-acterr
  himass = extmass+1.0073+acterr
  return(c(lowmass, himass))
}

MminusHrange = function(extmass, ppmerr){
  acterr = (extmass-1.0073)*(ppmerr*1e-6)
  lowmass = extmass-1.0073-acterr
  himass = extmass-1.0073+acterr
  return(c(lowmass, himass))
}


MolMatching = function(molset, featdefs, mode, mzerr, rterr, savefile = "matchedmols.csv"){
  # For troubleshooting:
  # molset = mymols
  # featdefs = gooddef
  # mode = "POS"
  # mzerr = 50
  # rterr = 120
  # rm(list = c("molset", "featdefs", "mode", "mzerr", "rterr"))
  
  if (mode == "POS"){
    molset = molset[which(molset[,"Mode"]=="POS"),]
    molset[,c("mzmin", "mzmax")] = t(apply(molset["ExactMass"], 1, MplusHrange, ppmerr=mzerr))
  } else if (mode == "NEG") {
    molset = molset[which(molset[,"Mode"]=="NEG"),]
    molset[,c("mzmin", "mzmax")] = t(apply(molset["ExactMass"], 1, MminusHrange, ppmerr=mzerr))
  } else if (mode == "HILIC") {
    molset = molset[which(molset[,"Mode"]=="HILIC"),]
    molset[,c("mzmin", "mzmax")] = t(apply(molset["ExactMass"], 1, MminusHrange, ppmerr=mzerr))
  } else {
    print("Mode not recognized. Enter POS, NEG, or HILIC")
  }
  
  for (i in c(1:nrow(molset)) ) {
    mzmatch = which ((featdefs[,"mzmed"] > molset[i,"mzmin"]) & (featdefs[,"mzmed"] < molset[i,"mzmax"]))
    
    if (is.na(molset[i,"RT_known"]) == TRUE){
      print(paste0("No RT given for ", molset[i,"Label"], " in row ", i))
      if (length(mzmatch>0)) {
        rtmatch = c(1:length(mzmatch))
      } else {
        rtmatch = mzmatch
      }
    } else {
      rtmatch = which ((featdefs[mzmatch,"rtmed"] > (molset[i,"RT_known"]-rterr) ) & (featdefs[mzmatch,"rtmed"] <  (molset[i,"RT_known"]+rterr) ) ) 
    }
    
    
    if (length(rtmatch) ==1) {
      molset[i,"FtMatch"] = rownames(featdefs)[mzmatch[rtmatch]]
    } else if (length(rtmatch) >1) {
      print(paste0("Multiple matches for ", molset[i,"Label"], " in row ", i))
      molset[i,"FtMatch"] = paste( rownames(featdefs)[mzmatch[rtmatch]], collapse = "," )
    } else {
      print(paste0("No match for ", molset[i,"Label"], " in row ", i))
    }
  }
  
  matchedfts = unlist(strsplit(paste( na.omit(molset[,"FtMatch"]), collapse = "," ), split = ","))
  
  write.csv(molset, file = savefile)
  matchinfo = list(molset, matchedfts)
  names(matchinfo) = c("molinfo", "featset")
  return(matchinfo)
}



#molinfo is an optional argument in this function. Put it in if you want names on your chromatograms

ChromPrint = function(featchrom, condkey, pdextra, molinfo){
  fdef = featureDefinitions(featchrom)
  ftset = rownames(fdef)
  ftheight = featureValues(featchrom, method="medret", value="maxo")

  sample_colors = pdextra[,"sample_color"]
  names(sample_colors) = rownames(pdextra)
  
  plotgroups = unique(condkey[,"plotgroup"])
  chartlocs = list()
  
  for (i in (1:length(plotgroups))) {
    conds = condkey[which(condkey[,"plotgroup"] == plotgroups[i]), "label"]
    chartlocs[[plotgroups[i]]] = which(pdextra[,"group_abbr"] %in% conds)
  }
  
  for (i in c(1:length(ftset) )){
    
    if (missing(molinfo)==FALSE) {
      if ("FtMatch"%in%colnames(molinfo) ){
        molloc = grep(ftset[i], molinfo[,"FtMatch"])
        imagename = substr(paste0(molinfo[molloc,"Label"], " ", ftset[i]), 1, 100)
        titletext = paste0(molinfo[molloc,"Label"], "  ", 
                           round(molinfo[molloc,"mzmin"], digits=4) , " - ", 
                           round(molinfo[molloc,"mzmax"], digits=4) , "  ", 
                           ftset[i] )
      } else if ("FormulaMatch"%in%colnames(molinfo)){
        imagename = substr(paste0(ftset[i], " ", molinfo[ftset[i],"Label"]), 1, 100)
        titletext = paste0(ftset[i], " mzmed:", 
                           round(molinfo[ftset[i],"mzmed"], digits=4) , "  ", 
                           molinfo[ftset[i],"FormulaMatch"], " ",
                           molinfo[ftset[i],"Label"] )
      }

    } else {
      imagename = ftset[i]
      titletext = paste0(ftset[i], "  ", 
                         round(fdef[ftset[i],"mzmin"], digits=4) , " - ", 
                         round(fdef[ftset[i],"mzmax"], digits=4) ) 
    }
    
    chtheight = max(ftheight[ftset[i],] , na.rm = TRUE)
    png(file = paste0(imagename, ".png"), width = 600*length(chartlocs), height = 600)
    par(mfcol = c(1, length(chartlocs)) )
    
    for (j in c(1:length(chartlocs))) {
      plot(featchrom[i,chartlocs[[j]]], 
           col = sample_colors[chartlocs[[j]]], peakType = "rectangle",
           peakCol = sample_colors[chartlocs[[j]]][chromPeaks(featchrom[i,chartlocs[[j]]])[, "column"]],
           peakBg = NA,
           ylim = c(0, chtheight) ,
           main = names(chartlocs)[j],
           cex.lab=2, cex.axis=2, cex.main=2)
    }
    
    mtext(titletext ,         
          side = 3,
          line = - 2,
          outer = TRUE,
          cex = 1)
    
    dev.off()
    
  }
  
  
}



MultiVolcano = function(dataset, ftsubset, molmatchinfo, comparekey, pthresh, l2fcthresh) {
  
  # For Troubleshooting
   # dataset=xdata_ftg
   # comparekey = tTestKey
   # pthresh = 0.05
   # l2fcthresh = 1
  # 
  # 
  # dir.create("./EditingMultiVolc")
  # setwd("./EditingMultiVolc")
  
  ftsort = featureValues(dataset, 
                           method="medret",   #Consider changing this
                           missing = "rowmin_half")
  colnames(ftsort) = dataset$sample_name
  if (missing(ftsubset) == FALSE) {
    ftsort = ftsort[ftsubset,]
  } 
  
  ttests = matrix(nrow = dim(ftsort)[1], ncol = dim(comparekey)[1] )
  rownames(ttests) = rownames(ftsort)
  colnames(ttests)= comparekey[,"compareID"]
  l2fcs = ttests
  
  for (i in c(1:dim(comparekey)[1])) {
    explocs = grep(comparekey[i,"expcond"], colnames(ftsort))
    nulllocs = grep(comparekey[i,"nullcond"], colnames(ftsort))
    for (j in c(1:dim(ftsort)[1]) ) {
      tryhold = try(t.test(ftsort[j,explocs], ftsort[j,nulllocs]), silent=TRUE)
      if (is(tryhold, "try-error")) {
        ttests[j,i] = 1
      } else {
        ttests[j,i] = tryhold$p.value
      }
      
      l2fcs[j,i] = log2(mean(ftsort[j,explocs]) /  mean(ftsort[j,nulllocs]))
    }
    
    #p-value multiple hypothesis correction
    ttests[,i] = p.adjust(ttests[,i], method="BH")

    #Plot the volcano here using the column i from ttests and l2fcs
    if (missing(molmatchinfo)) {
      volcframe = data.frame(FTName = rownames(ttests), 
                             P = ttests[,i], 
                             EFFECTSIZE = l2fcs[,i]
      )
      volcobj = volcanor(volcframe, snp = "FTName")
      
    } else {
      volcframe = data.frame(FTName = rownames(ttests), 
                             P = ttests[,i], 
                             EFFECTSIZE = l2fcs[,i], 
                             Formula = molmatchinfo[,"FormulaMatch"], 
                             CandidateID = molmatchinfo[,"Label"]
      )
      volcobj = volcanor(volcframe, snp = "FTName", gene = "Formula", annotation1 = "CandidateID")
    }
    

    volcwidg = volcanoly(volcobj, 
                         effect_size_line = c(l2fcthresh*-1, l2fcthresh),
                         genomewideline = -log10(pthresh),
                         xlab = "log2foldchange",
                         title = colnames(ttests)[i]  )
    
    htmlwidgets::saveWidget(
      widget = volcwidg, 
      file = paste0(colnames(ttests)[i], "_VolcanoPlot.html"), 
      selfcontained = TRUE #creates a single html file  
      )
  }
  
  cats = unique(comparekey[,"category"])
  tmins = matrix(nrow=nrow(ttests), ncol=length(cats))
  rownames(tmins) = rownames(ttests)
  fcabsmaxs <- tmins
  colnames(tmins) = paste0(cats,"_tmin")
  colnames(fcabsmaxs) = paste0(cats,"_fcabsmax")
  sigfts = vector(mode="list", length = length(cats))
  names(sigfts) = cats
  
  for (catname in cats) {
    print(catname)
    catlocs = grep(catname, comparekey[,"category"])
    if (length(catlocs) > 1) {
    tmins[,paste0(catname,"_tmin")] = rowMins(ttests[,catlocs])
    fcabsmaxs[,paste0(catname,"_fcabsmax")] = rowMaxs(abs(l2fcs[,catlocs]))
    } else {
      tmins[,paste0(catname,"_tmin")] = ttests[,catlocs]
      fcabsmaxs[,paste0(catname,"_fcabsmax")] = abs(l2fcs[,catlocs])
    }
    
    smallt = which(tmins[,paste0(catname,"_tmin")] < pthresh)
    bigfc = which(fcabsmaxs[,paste0(catname,"_fcabsmax")] > l2fcthresh)
    if ( (length(smallt) !=0) & (length(bigfc) !=0) ) {
      broadchems = intersect(names(smallt), names(bigfc))
      sigfts[[catname]] = broadchems[order(fcabsmaxs[broadchems, paste0(catname,"_fcabsmax")], decreasing=TRUE)]
      
    } else {
      sigfts[[catname]] = NULL
    }
  }
  
  ttests = cbind(ttests, tmins)
  l2fcs = cbind(l2fcs, fcabsmaxs)
  
  multivolcout = list(ttests, l2fcs, sigfts)
  names(multivolcout) = c("ttests", "l2fcs", "sigfts")
  return(multivolcout)
  
}


FeatMatching = function(featdefs, molset, mode, mzerr, rterr){
  # For troubleshooting:
   # molset = mymols
   # featdefs = gooddef
   # mode = "POS"
   # mzerr = 50
   # rterr = 120
  # rm(list = c("molset", "featdefs", "mode", "mzerr", "rterr"))
  
  if (mode == "POS"){
    molset[,c("mzmin", "mzmax")] = t(apply(molset["ExactMass"], 1, MplusHrange, ppmerr=mzerr))
  } else if (mode == "NEG") {
    molset[,c("mzmin", "mzmax")] = t(apply(molset["ExactMass"], 1, MminusHrange, ppmerr=mzerr))
  } else if (mode == "HILIC") {
    molset[,c("mzmin", "mzmax")] = t(apply(molset["ExactMass"], 1, MminusHrange, ppmerr=mzerr))
  } else {
    print("Mode not recognized. Enter POS, NEG, or HILIC")
  }
   
  for (i in c(1:nrow(featdefs)) ) {
    
    mzmatch = which ((featdefs[i,"mzmed"] > molset[,"mzmin"]) & (featdefs[i,"mzmed"] < molset[,"mzmax"]))
    
    if (is.na(molset[i,"RT_known"]) == TRUE){
      if (length(mzmatch>0)) {
        rtmatch = c(1:length(mzmatch))
      } else {
        rtmatch = mzmatch
      }
    } else {
      #dunno if this part works. I don't have RTs in my KEGG output stuff.
      rtmatch = which ((featdefs[i,"rtmed"] > (molset[mzmatch,"RT_known"]-rterr) ) & (featdefs[i,"rtmed"] <  (molset[mzmatch,"RT_known"]+rterr) ) ) 
    }
    
    if (length(rtmatch) ==1) {
      featdefs[i,"FormulaMatch"] = molset[,"Formula"][mzmatch[rtmatch]]
      featdefs[i,"ExactMass"] = molset[,"ExactMass"][mzmatch[rtmatch]]
      featdefs[i,"RT_known"] = molset[,"RT_known"][mzmatch[rtmatch]]
      featdefs[i,"Label"] = molset[,"Label"][mzmatch[rtmatch]]
    } else if (length(rtmatch) >1) {
      featdefs[i,"FormulaMatch"] = paste( molset[,"Formula"][mzmatch[rtmatch]], collapse = ", " )
      featdefs[i,"ExactMass"] = paste( molset[,"ExactMass"][mzmatch[rtmatch]], collapse = ", " )
      featdefs[i,"RT_known"] = paste( molset[,"RT_known"][mzmatch[rtmatch]], collapse = ", " )
      featdefs[i,"Label"] = paste( molset[,"Label"][mzmatch[rtmatch]], collapse = ", " )
      
    } 
  }
  
  featswmatch = featdefs[which(is.na(featdefs[,"FormulaMatch"]) == FALSE), ]
  return(featswmatch)
}


addmzrange = function(molset, mode){
  if (mode == "POS"){
    molset[,c("mzmin", "mzmax")] = t(apply(molset["ExactMass"], 1, MplusHrange, ppmerr=mzerr))
  } else if (mode == "NEG") {
    molset[,c("mzmin", "mzmax")] = t(apply(molset["ExactMass"], 1, MminusHrange, ppmerr=mzerr))
  } else {
    print("Mode not recognized. Enter POS or NEG")
  }
  return(molset)
}


#This plot function is going to have to assume you tried to find peaks on your chromatogram.
#Otherwise peakCol will fart out. 
#would have to take out peakType, peakCol and peakBg to make it work.
#Currently may throw warnings when there aren't any peaks in your plot range. 

ChromsByConds <- function(chrtoplot, graphconds, RTrange, molinfo){
  #Find which samples in the chrtoplot correspond to each graphconds entry
  chartlocs = list()
  
  for (k in c(1:length(graphconds))){
    chartlocs[[k]] = grep(graphconds[k], pData(chrtoplot)[,"sample_name"] )
  }
  
  #Make each image file
  for (i in c(1:nrow(molinfo))){
    
    imagename = paste0(substr(molinfo[i,"Label"], 1, 80), "_bycond")
    titletext = paste0(
      substr(molinfo[i,"Label"], 1, 80), "  ",
      round(molinfo[i,"mzmin"], digits=4) , " - ", 
      round(molinfo[i,"mzmax"], digits=4) 
    )
    
    xlims = NULL
    if (class(RTrange) == "numeric"){
      if (is.na(molinfo[i, "RT_known"]) == FALSE) {
        xlims = c(molinfo[i,"RT_known"]-RTrange, molinfo[i,"RT_known"]+RTrange) 
        if (xlims[1] < 0){
          xlims[1] <- 0
        }
      }
    } 
    
    png(file = paste0(imagename, ".png"), width = 1200, height = 400*length(chartlocs))
    par(mfrow = c(length(chartlocs), 1),
        oma = c(0, 0, 4, 0) )
    
    #Make each subplot
    for (j in c(1:length(chartlocs))) {
      arbcols = brewer.pal(n=length(chartlocs[[j]]), name="Set1")
      
      plot(chrtoplot[i,chartlocs[[j]]], 
           col = arbcols, 
           peakType = "rectangle",
           peakCol = arbcols[chromPeaks(chrtoplot[i,chartlocs[[j]]])[, "column"]],
           peakBg = NA,
           xlim = xlims, 
           main = paste0("Condition ",graphconds[j]),
           cex.lab=2, cex.axis=2, cex.main=2)
    }
    
    mtext(titletext ,         
          side = 3,
          line = 1,
          outer = TRUE,
          cex = 2)
    
    dev.off()
  }
}


#A vestigial part from trying to get the height to scale things
#chtheight = max(unlist(lapply(chrtoplot, intensity)), na.rm=TRUE)


