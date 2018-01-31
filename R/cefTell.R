# Process Pipeline to gain ion statistics

# read files into list of ceftell objects
cef.list = lapply(cef.fh[2], function(x) cefTell(x))


#Process isobars
cef.list = lapply(cef.list,function(x) getIsobarsCef(x))


# Get the compound list
cef.list = getCompoundList(cef.list)

# Define the ion types
cef.list$ionTypes = names(cef.list$CefFileData[[1]]$sumIonDataBySample)


# Get the ion frequencies
ionFrequencies = getCmpdFrequencyByIonType(cef.list)

# Get sum intensities per cmpd
ionIntensities = getCmpdSumIntensitiesByIonType(cef.list)


# Get mean intensities per cmpd
mean.ionIntensities = ionIntensities/ionFrequencies

# Assign values to object for export
cef.list$ionFrequencies = ionFrequencies
cef.list$sumIonIntensities = ionIntensities
cef.list$meanIonIntensities = ionIntensities/ionFrequencies
save(file="ceflist_final.rda",cef.list)

#Sum of square errors
SSQIntensities = getCmpdSSQByIonType(cef.list)

#A Test
ionFrequencies_1 = as.matrix(ionFrequencies-1)
ionFrequencies_1[which(ionFrequencies_1 < 2)] = NA

# Deviations and CV
SDIonIntensities = sqrt(SSQIntensities/(ionFrequencies_1))
CVIonIntensities = SDIonIntensities/mean.ionIntensities
#########################################################################

#save(file="ceflist.rda",cef.list)


## EXPORT DATA

cef.list$CMPDS
cef.list$ionFrequencies
cef.list$sumIonIntensities

cef.list$Mass = as.numeric(sapply(cef.list$CMPDS,function(x) strsplit(x,"@")[[1]][1]))
cef.list$RT = as.numeric(sapply(cef.list$CMPDS,function(x) strsplit(x,"@")[[1]][2]))

cef.df = data.frame("ID" = cef.list$CMPDS,
                   "Mass" = cef.list$Mass,
                   "RT" = cef.list$RT,
                   cef.list$meanIonIntensities)
