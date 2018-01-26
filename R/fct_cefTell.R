
#Cef Tell: functions to read multiple Profinder cef files and to obtain ion statistics. 
cefTell = function(fh,annot=F,frag=F){
  
  require(xml2)
  require(magrittr)
  
  # Read the cef xml file
  xml  = read_xml(fh)
  cat(paste("Processing file: ",fh,"\n"))

  # Gets the compound nodes
  Compounds = xml %>% xml_find_all("//Compound")
  cat(paste("    # Compounds: ",length(Compounds),"\n"))

  # Gets the location nodes mass and retention time
  Locations = Compounds %>% xml_find_all("Location")
  M = Locations  %>% xml_attr("m") %>% as.numeric()
  rt = Locations  %>% xml_attr("rt") %>% as.numeric()
  
  #ID = paste(round(M,4),"@",round(rt,8),sep="")
  ID = paste(M,"@",rt,sep="")
  
  Scores= Compounds %>% xml_find_all("CompoundScores")
  Spectra = Compounds %>% xml_find_all("Spectrum")
  
  MSPeaks  = Spectra %>% xml_find_all("MSPeaks")
  cat(paste("     Reading isotope peaks, this takes a while","\n"))
  Results = Compounds %>% xml_find_all("Results")
 
 
  
  
  
  peaks = lapply(MSPeaks, function(x) x %>% xml_find_all("p"))
  
  

  peakData = list(
    File = fh,
    "ID"= ID,
    mz = lapply(peaks,function(x) x %>% xml_attr("x") %>% as.numeric()),
    mzx =lapply(peaks,function(x) x %>% xml_attr("rx") %>% as.numeric()),
    area = lapply(peaks,function(x) x %>% xml_attr("y") %>% as.numeric()),
    z = lapply(peaks, function(x) x %>% xml_attr("z") %>% as.numeric()),
    ions  = sapply(peaks, function(x) x %>% xml_attr("s"))
  )
 if(frag==T){
  Fragments = Compounds %>% xml_find_all("FragmentIonConfirmation") %>% xml_find_all("MSDetails") %>% xml_find_all("FragmentIonPeaks")
  
  
  
  if(length(Fragments)> 0 ){
    fragmentpeaks = lapply(Fragments, function(x) x %>% xml_find_all("fp"))
    
    fragmentpeakData = list(
      File = fh,
      "ID"= ID,
      
      mz = lapply(fragmentpeaks,function(x) x %>% xml_attr("x") %>% as.numeric()),
      area = lapply(fragmentpeaks,function(x) x %>% xml_attr("y") %>% as.numeric()),
      sn = lapply(fragmentpeaks, function(x) x %>% xml_attr("sn") %>% as.numeric()),
      ce  = lapply(fragmentpeaks, function(x) x %>% xml_attr("ce") %>% as.numeric()),
      source = lapply(fragmentpeaks, function(x) x %>% xml_attr("source") %>% as.character())
    )
    fragment.rt  = parentmz = NULL
    for(i in 1:length(fragmentpeaks)){
      n.fragments = length(fragmentpeakData$mz[[i]]) 
      fragment.rt = c(fragment.rt,rep(rt[i],n.fragments))
      parentmz = c(parentmz,rep(M[i]+1.00728,n.fragments))
    }
    

    
  }

   peakData$Fragments = fragmentpeakData
  }
  if(annot == T){
    print("Retrieving annotations")
    Molecule = Results %>% xml_find_all("Molecule")
    Formula = gsub(" ","",Molecule %>% xml_attr("formula")   )
    Name = Molecule %>% xml_attr("name")
    
    MatchScores = Molecule %>% xml_find_all("MatchScores")
 
    Database = OverallScore = DBScore = HMDB=LMDB=KEGG=NULL
    
    for(i in 1:length(Compounds)){
      
      
      child  = xml_child(Molecule[[i]]) %>% xml_find_all("Match")
      algo = child %>% xml_attr("algo")
      m = match("overall",algo)
    
      if(!is.na(m)){
        
        OverallScore[i] = as.numeric(xml_attr(child[m],"score") )
      }else{  
        OverallScore[i] = NA  
      }
      
      
      m = match("db",algo)
      if(!is.na(m)){
         DBScore[i] = as.numeric(xml_attr(child[m],"score"))
      }else{
        DBScore[i] = NA
      }
      
      child  = xml_child(Results[[i]]) %>% xml_find_all("Database") %>% xml_find_all("Accession")
      if(length(child)>0){
        db = child %>% xml_attr("db")

        #HMDB
        m= match("HMP ID", db)
        if(!is.na(m)){
          HMDB[i]= xml_attr(child[m],"id") 
        }else{
          HMDB[i] = NA
        }


        #LipidMaps
        m= match("Lipid ID", db)
        if(!is.na(m)){
          LMDB[i]= xml_attr(child[m],"id") 
        }else{
          LMDB[i] = NA
        }

        #Kegg
        m= match("KEGG ID", db)
        if(!is.na(m)){
          KEGG[i]= xml_attr(child[m],"id") 
        }else{
          KEGG[i] = NA
        }

        
        

      }else{
        HMDB[i] =  LMDB[i]= KEGG[i]= NA
      }
      
    }
    
    peakData$Formula = Formula 
    peakData$Name = Name
    peakData$HMDB = HMDB
    peakData$LMDB = LMDB
    peakData$KEGG = KEGG 
    peakData$Score = OverallScore
  }
  print("Getting ion types")
  peakData = getIonTypes(peakData)
  print("Getting Ion Peak Data")
  if(frag==F){
  peakData = getIonPeakData(peakData)
  print("Getting Ratios by Adduct")
  peakData = getAdductRatios(peakData)
  }  
  peakData
}

getIonTypes = function(peakData){   
  
  ionTypes = unique(unlist(peakData$ions))
  peakData$isotopes = ionTypes
  g = grep("\\+[1-9]{1}$",ionTypes)
  
  isotopes = ionTypes[g]
  ionTypes = ionTypes[-g]
  
  cat(paste("    # Adduct types: ",length(ionTypes),"\n"))
  out = NULL
  for(ion in ionTypes){
    out[[ion]] = c(ion,isotopes[grep(   paste("^",gsub("\\+","\\\\+",ion),"\\+[1-9]{1}",sep="") ,isotopes) ])
  }
  
  peakData$ionTypes = out

  peakData
}

####################################################################
####################################################################

getAdductRatios = function(peakData){
    # Get adduct types
    adducts = colnames(peakData$sumIonDataByCompound$MZ)
    out.mono = out.ID = out.rats= matrix(nrow = length(peakData$ID), ncol = dim(peakData$sumIonDataByCompound$MZ)[2])
    
    #adduct = "^M\\+H$"
    k=1
    for(type in adducts){
       print(type)
       type = gsub("\\+","\\\\+",type)
       type = paste("^",type,sep="")
      
       out.ID = peakData$ID
       #out.ID = out.area = out.mz= NULL
    
       i=j=1
       for(cmpd in peakData$ions){
          
          g = grep(type,cmpd)

          if(length(g)>0){
            
             g = g[1]
             iso.topes ="\\+1"
             iso.adduct = paste(type,iso.topes,sep="")
             #iso.adduct = gsub("\\+","\\\\+",iso.adduct)
             
             
             g = c(g[1],grep(iso.adduct,cmpd))
        
             if(length(g)>1){
               
               out.mono[i,k] = peakData$mz[[i]][g][1]
               x  = peakData$area[[i]][g]
               out.rats[i,k] = (x[1]/x[1])/(x[2]/x[1])
               #print(out.mono[i,k])
               #print(out.rats[i,k])
             } 

          }
      
       i=i+1
       }
       
       #out.mono[,k] = sapply(out.mz, function(x) x[1])
       #out.rats[,k] = sapply(out.area,function(x) (x[1]/x[1])/(x[2]/x[1])) 
       k=k+1
    
    
    }
    
    
    
    
    out.M = sapply(out.ID,function(x) as.numeric(strsplit(split="@",x)[[1]][1]))
    out.RT= sapply(out.ID, function(x) as.numeric(strsplit(split="@",x)[[1]][2]))

    out = data.frame("M" = out.M,"RT" = out.RT,"mz" = out.mono, "ratio" =out.rats)
    peakData$ratios = out
    
    peakData
}







# Function gets the ion peak data from the 
getIonPeakData = function(peakData){
  cat(paste("     Processing adducts\n"))
  ionTypes = names(peakData$ionTypes)
  sumIonDataBySample = NULL
  sumIonDataByCompound = list("MZ"=NULL,"Intensities"=NULL);
  
  l = length(peakData$ions)
    
  for(i in 1:l){
    
    options(warn=-1)
    compoundIonData = ionIntensities = ionMZ = NULL
    
    for(ion in ionTypes){
      isotopes = peakData$ionTypes[ion]
      
      
      m = match(unlist(isotopes),peakData$ions[[i]])
      m = m[which(!is.na(m))]
      
      if(length(m)==0){
        ionIntensities[ion] = NA
        ionMZ[ion] = NA
      }else{ 
        ionIntensities[ion] = (sum(peakData$area[[i]][m],na.rm=T))
        ionMZ[ion] = peakData$mz[[i]][m][1]
      } 
    }
    
    sumIonDataByCompound$Intensities = rbind(sumIonDataByCompound$Intensities,ionIntensities)
    sumIonDataByCompound$MZ = rbind(sumIonDataByCompound$MZ,ionMZ)
    options(warn=0)
    
  }
  
  rownames(sumIonDataByCompound$MZ) = peakData$ID
  rownames(sumIonDataByCompound$Intensities) = peakData$ID
  peakData$sumIonDataByCompound = sumIonDataByCompound
  peakData$sumIonDataBySample = apply(sumIonDataByCompound$Intensities,2,sum,na.rm=T)
  peakData
}



getIsoBars = function(cef){
  
  print(cef$File)
  
  isobars = w = NULL;
  duplicates = which(duplicated(cef$ID)==T)
  
  for(isobar in duplicates){
    w = c(w,which(cef$ID == cef$ID[isobar]))
    
  }
  
  w 
}


getIsobarsCef = function(cef){
  w = getIsoBars(cef)
  cef.new=cef
  if(!is.null(w)){
    cef.isobars = list(
     "ID"  = cef$ID[w],
     "mz"  = cef$mz[w],
     "mzx" = cef$mzx[w],
     "area"  = cef$area[w],
     "z" = cef$z[w],
     "ions" = cef$ions[w],
     "sumIonDataByCompound" = list("MZ"  = cef$sumIonDataByCompound$MZ[w,],
                                  "Intensities" = cef$sumIonDataByCompound$Intensities[w,])
  );
  
    cef.new = list(
      "File" = cef$File,
      "ID"  = cef$ID[-w],
      "mz"  = cef$mz[-w],
      "mzx" = cef$mzx[-w],
      "area"  = cef$area[-w],
      "z" = cef$z[-w],
      "ions" = cef$ions[-w],
    
      "sumIonDataByCompound" = list("MZ" = cef$sumIonDataByCompound$MZ[-w,],
                                  "Intensities" = cef$sumIonDataByCompound$Intensities[-w,]),
      "sumIonDataBySample" = cef$sumIonDataBySample,
      "isobars" = cef.isobars
    )
  
  
  }

  cef.new  
  
}


getCompoundList = function(cef.list){
  ID = NULL;
  
  for(df in cef.list){
    ID = c(ID,df$ID)
  }
  cef.list = list("CMPDS"  =  unique(ID), "CefFileData" = cef.list)
}




getCmpdFrequency = function(cef.list){
  CMPDS = cef.list$CMPDS
  freq = array(data = 0,dim=length(CMPDS))
  i=1
  
  for(cmpd in CMPDS){
    
    for(cef in cef.list$CefFileData){
      w = which(cef$ID == cmpd)
      if(length(w) > 0){
        freq[i] = freq[i]+1
      } 
    }
    i=i+1  
  }  
}





getCmpdFrequencyByIonType = function(cef.list){
  
  CMPDS = cef.list$CMPDS
  ionFrequencies = matrix(nrow = length(CMPDS),ncol = length(cef.list$ionTypes))  
  
  j = 1

  for(ion in cef.list$ionTypes){
   print(ion)

    freq = array(data = 0,dim=length(CMPDS))
    i=1
    
    for(cmpd in CMPDS){
      
      
      for(cef in cef.list$CefFileData){

        w = which(cef$ID == cmpd)
        if(length(w) > 0){
          for(id in w)
          g = match(ion,cef$ions[[w]])
          if(!is.na(g) ){
             freq[i] = freq[i]+1
          }
        } 
      }
      i=i+1           

      
    }

    ionFrequencies[,j] = freq
    j=j+1
  }
   
  #ionFrequencies = as.data.frame(ionFrequencies)
  rownames(ionFrequencies)  = CMPDS
  colnames(ionFrequencies)  = cef.list$ionTypes
  ionFrequencies[which(ionFrequencies==0)] = NA
  ionFrequencies
  
}


getCmpdSumIntensitiesByIonType = function(cef.list){
  CMPDS = cef.list$CMPDS
  ionIntensities = matrix(nrow = length(CMPDS),ncol = length(cef.list$ionTypes))  

  j = 1
  CMPDS = cef.list$CMPDS 

  for(ion in cef.list$ionTypes){
 
    
    sumIntensity = array(data = 0,dim=length(CMPDS))
    i=1
    
    for(cmpd in CMPDS){
      for(cef in cef.list$CefFileData){
        
        w = which(cef$ID == cmpd)
        if(length(w) > 0){
          if(!is.na(cef$sumIonDataByCompound$Intensities[w,j])){
             sumIntensity[i] = sumIntensity[i]+cef$sumIonDataByCompound$Intensities[w,j]
          }
        } 
      }
      i=i+1           
      
      
    }
    ionIntensities[,j] = sumIntensity
    j=j+1
  }
  
  #ionIntensities = as.data.frame(ionIntensities)
  rownames(ionIntensities)  = CMPDS
  colnames(ionIntensities)  = cef.list$ionTypes
  ionIntensities[which(ionIntensities==0)] = NA
  ionIntensities
  
}


getCmpdSSQByIonType = function(cef.list){
  CMPDS = cef.list$CMPDS
  meanIonIntensities = cef.list$meanIonIntensities
  SSQIntensities = SSQRT = matrix(nrow = length(CMPDS),ncol = length(cef.list$ionTypes))  
  
  j = 1
  
  CMPDS = cef.list$CMPDS 
  
  for(ion in cef.list$ionTypes){
    
    
    SSQIonIntensity = array(data = 0,dim=length(CMPDS))
    
    i=1
    
    for(cmpd in CMPDS){
      
      for(cef in cef.list$CefFileData){
        
        w = which(cef$ID == cmpd)
        if(length(w) > 0){
          if(!is.na(cef$sumIonDataByCompound$Intensities[w,j])){
            #sumIntensity[i] = sumIntensity[i]+cef$sumIonDataByCompound$Intensities[w,j]
            SSQIonIntensity[i] = SSQIonIntensity[i] + ((cef$sumIonDataByCompound$Intensities[w,j]-meanIonIntensities[w,j])^2)
          }
        } 
      }
      i=i+1           
      
      
    }
    SSQIntensities[,j] = SSQIonIntensity
    j=j+1
  }
  
  #SSQIntensities = (SSQIntensities)
  rownames(SSQIntensities)  = CMPDS
  colnames(SSQIntensities)  = cef.list$ionTypes
  SSQIntensities[which(SSQIntensities==0)]=NA
  SSQIntensities  
}





