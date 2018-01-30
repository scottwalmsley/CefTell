
#Functions to parse annotations from hmdb and lmdb databases in flat file format
getAnnotations = function(data,outfile,ppmtol=10, mfg=70,db=70,maxRTdiff = 0.2){
  require(stringr)
  require(rcdk)
  #data = infile
  data  = read.delim(infile,check.names=F,sep="\t",comment="#")
  count.na = function(x) sum(!is.na(x))
  
  
  Mass  = as.numeric(as.character(data$Mass))
  rt    = as.numeric(as.character(data$"Retention Time"))
  
  re <- "\\[ ([^()]+) \\]"
  annot = as.character(data$Annotations)
  annot = (sapply(annot,function(x)str_extract_all(x, re)[[1]]))
  names(annot) = NULL
  l = sapply(annot,function(x) length(x))
  
  annot[which(l==0)]  =  ""
  annot = sub("\\[","",annot)
  annot = sub("\\]","",annot)
  
  cmpd =  as.character(data$Compound)
  formulae = sapply(annot,function(x) strsplit(x,",")[[1]][1])
  names(formulae) = NULL
  formulae = gsub(" ","",formulae)
  formulae[grep("<",formulae)] = NA
  formulae[grep(">",formulae)] = NA
  expected.mass = array(dim=dim(data)[1])
  i=1
  for(f in formulae){
    #print(i)
    if(!is.na(f) & !grepl("D",formulae[i])){ 
      m = get.formula(f,charge=0)@mass
    }else{
      m=NA
    }
    expected.mass[i] = m
    i=i+1
  }
  
  score = sapply(annot,function(x) strsplit(x,",")[[1]][3])
  names(score) = NULL
  score = sub("overall=","",score)
  score = sub(" ","",score)
  score = as.numeric(score)
  
  sub.score = sapply(annot,function(x) strsplit(x,",")[[1]][4])
  names(sub.score) = NULL
  
  mfg = array(dim=length(score))
  mfg[grep("mfg",sub.score)] = sub.score[grep("mfg",sub.score)]
  mfg = as.numeric(sub("mfg=","",mfg))
  
  db = array(dim=length(score))     
  db[grep("db",sub.score)] = sub.score[grep("db",sub.score)]
  db = as.numeric(sub("db=","",db))
  
  dbID = array(dim=length(score))
  dbID = sapply(annot,function(x) strsplit(x,",")[[1]][5])
  names(dbID)=NULL
  
  dbID.2 = array(dim=length(score))
  dbID.2 = sapply(annot,function(x) strsplit(x,",")[[1]][6])
  names(dbID.2)=NULL
  
  dbID.3 = array(dim=length(score))
  dbID.3 = sapply(annot,function(x) strsplit(x,",")[[1]][7])
  names(dbID.3)=NULL
  
  dbID.4 = array(dim=length(score))
  dbID.4 = sapply(annot,function(x) strsplit(x,",")[[1]][8])
  names(dbID.4)=NULL
  
  dbID.5 = array(dim=length(score))
  dbID.5 = sapply(annot,function(x) strsplit(x,",")[[1]][9])
  names(dbID.5)=NULL
  
  dbIDList = data.frame(dbID,dbID.2,dbID.3,dbID.4,dbID.5,stringsAsFactors=F)
  
  l = dim(data)[1]
  metlin = array(dim=l)
  kegg = data$`KEGG ID`      
  
  LipidMap = data$`LMP ID`
  hmp = data$`HMP ID`
  
  cas = data$`CAS Number`
  
  ###loop to get the id;s
  for(i in 1:5){
    g = grep("METLIN",dbIDList[,i])
    metlin[g] = dbIDList[g,i]
    #  g = grep("Lipid ID",dbIDList[,i])
    #  LipidMap[g] = dbIDList[g,i]
    # g = grep("CAS",dbIDList[,i])
    #  cas[g] = dbIDList[g,i]
    # g = grep("KEGG",dbIDList[,i]) 
    #  kegg[g] = dbIDList[g,i]
    # g = grep("HMP",dbIDList[,i])
    #hmp[g] = dbIDList[g,i]
    
  }
  
  #hmp = sub("HMDB","HMDB0",hmp)
  metlin = sub("METLIN ID=","",metlin)
  metlin = gsub(" ","",metlin)  
  metlin = as.numeric(metlin)
  
  #kegg = sub("KEGG ID=","",kegg)
  #kegg= gsub(" ","",kegg) 
  
  #LipidMap = sub("Lipid ID=","",LipidMap)
  #LipidMap = gsub(" ","",LipidMap) 
  
  #hmp = sub("HMP ID=","",hmp)
  #hmp = gsub(" ","",hmp) 
  
  #cas = sub("CAS ID=","",cas)
  #cas = gsub(" ","",cas) 
  
  
  #### annotation workflow
  PROTON = 1.007262
  
  massTol = function(M,ppm) {
    as.numeric((ppm * M) / 1e6)
  }
  
  getRange = function(M,ppm){
    tol = massTol(M,ppm)
    c(M-tol,M+tol)
  }
  
  ppm.err = function(M,Mhat){
    (M-Mhat)/Mhat * 1e6
  }
  
  options(warn=-1)
  hm.mass = as.numeric(as.character(hmdb$MW))
  lm.mass = as.numeric(as.character(lmdb$exactMass))
  hm.form = as.character(hmdb$FORMULA)
  lm.form = as.character(lmdb$formula)
  options(warn=0)
  
  #### search for near masses in 2 ppm tol window
  r.hmid = array(dim=dim(data)[1])  ## by mass
  r.lmid = array(dim=dim(data)[1])  
  r.hmf  = array(dim=dim(data)[1])  ## by formula
  r.lmf  = array(dim=dim(data)[1])
  cnt.hm = array(dim=dim(data)[1])
  cnt.lm = array(dim=dim(data)[1])
  
  i=1
  for(M in expected.mass){
    #for(M in data$Mass){
    
    ### get the score
    cpd = cmpd[i]
    cpd.scr = score[i]
    #cpd.mfg = mfg[i]
    #cpd.db = db[i]
    cpd.form = formulae[i]
    
    pr = getRange(M,ppmtol)
    
    w = which(hm.mass > pr[1] & hm.mass < pr[2])
    r.hmid[i]= (paste(hmdb$ACCESSION[w],collapse=","))
    cnt.hm[i] = length(w)
    w = which(lm.mass > pr[1] & lm.mass < pr[2])
    r.lmid[i]= (paste(lmdb$LMID[w],collapse=","))
    cnt.lm[i] = length(w)
    r.hmf[i] = paste(hmdb$ACCESSION[which(hm.form == cpd.form)],collapse=",")
    r.lmf[i] = paste(lmdb$LMID[which(lm.form == cpd.form)],collapse=",")
    i=i+1
  }
  
  
  #ion = array(dim=l)
  ### mass corrections
  Mass  = as.numeric(as.character(data$Mass))
  delta.mz = round(Mass-expected.mass,2)
  
  
  ppm = (Mass-expected.mass)/expected.mass *1e6
  ppm = round(ppm,2)  
  
  
  ##  now we are set to deine levels of annotation
  quality = array(dim=l)
  quality[which(mfg<70)] = 0                 # level 0 
  quality[which( is.na(mfg)&is.na(db))] = 0  # level 0   
  quality[which(mfg>70)] = 1                 # level 1
  quality[which(db <70)] = 1                 # level 1
  quality[which(db>70)]  = 2                 # level 2
  quality[which(db > 70 & cnt.hm ==1)]=3     # level 3
  ## add mslevel 2 at later time 
  quality[grep("D",formulae)]  = 5           # level 5
  
  
  m= match(as.character(hmp),as.character(hmdb$ACCESSION)  )
  w = which(!is.na(m)) 
  pathway = array(dim=l)
  origin = array(dim=l)
  class = superclass = subclass= array(dim=l)
  lipid = array(dim=l)
  
  #  pathway[w] = as.character(hmdb$pathway[m[w]])
  subclass[w] = as.character(hmdb$SUBCLASS[m[w]])
  superclass[w] = as.character(hmdb$SUPERCLASS[m[w]])
  class[w] = as.character(hmdb$CLASS[m[w]])
  origin[w] = as.character(hmdb$ENDOGENOUS[m[w]])
  
  
  m = match(LipidMap,lmdb$LMID)
  w = which(!is.na(m))
  lipid[w] = as.character(lmdb$category[m[w]])
  
  
  out = data.frame("COMPOUND" = cmpd,
                   #"N" = data$Frequency,
                   "FORMULA" = formulae,
                   "QUALITY" = quality,
                   #"OBSERVED ION" = ion,
                   "MASS" = Mass,
                   #"EXPECTED MASS" = expected.mass,
                   #"DELTA MASS" = delta.mz,
                   #"DELTAPPM" = ppm,  
                   "RETENTION TIME"=rt,
                   "SCOREDB"=db,
                   "SCOREMFG" = mfg,
                   "HMDB" = hmp,
                   "HMDB CLASS" = class,
                   "HMDB SUBCLASS" = subclass,
                   "HMDB SUPERCLASS" = superclass,
                   "HMDB ORIGIN" = origin,
                   "HMDB ADDITIONAL" = r.hmf,
                   "LMDB" = LipidMap,
                   "LMDB CATEGORY" = lipid,
                   "LMDB ADDITIONAL" = r.lmf,
                   
                   #"HMDB MatchedFormula" = r.hmf,
                   #"LMDB MatchedFormula" = r.lmf,
                   "KEGG" = kegg,
                   "Metlin" = metlin, 
                   "CAS" = cas,
                   stringsAsFactors=F,check.names=F)
  
  print(dim(out))
  #w = which(quality!=0)
  #out = out[w,]
  #delta.mz = delta.mz[w]   
  #print(dim(out))
  
  cmpd = as.character(out$COMPOUND)
  
  #### filter pesky  names
  
  #grep(" \\+ [0-9]*\\.[0-9]+",as.character(cmpd))
  nm = gsub("nullnullnull","+",as.character(cmpd))
  nm  = gsub(" \\+ [0-9]*\\.[0-9]+","", as.character(nm))
  g.nm = gsub("\\[","\\\\[",nm)
  g.nm = gsub("\\]","\\\\]",g.nm)
  g.nm = gsub(":[0-9]+$","",g.nm)
  tmp.nm = g.nm
  g.nm = unique(g.nm)
  
  i = 1;
  IDX = list(dim=length(g.nm))
  removal.IDX = NULL
  temp.dset = NULL
  
  
  for(n in g.nm){
    
    
    # l = grep(paste("^",as.character(n),"$",sep=""),as.character(tmp.nm))
    l = which(as.character(tmp.nm)==as.character(n))
    
    dset = out[l,]
    
    if(length(l) > 1 ){
      
      ### test for rt clusters
      d = dist(dset$"RETENTION TIME",method="euclidean")
      h = hclust(d)
      cut = cutree(h,h=maxRTdiff)
      tm = cbind(dset,"class" = cut)
      nC = length(unique(cut)) 
      
      if(length(l) == nC){
        temp.dset = rbind(temp.dset,dset)
        removal.IDX = c(removal.IDX,l)
      }
      ## now take care of redundant names
      
      if(length(l) > nC){
        ## filter for maximum quality score
        max.Q.i =  which(dset$QUALITY == max(dset$QUALITY))
        dset = dset[max.Q.i,]
        
        ## get the index to remove the bad  entrys
        removal.IDX = c(removal.IDX,l[-max.Q.i])
        l = l[max.Q.i]         
        
        ## now take the first by naming convention
        temp.dset = rbind(temp.dset,dset[1,])
        removal.IDX = c(removal.IDX,l)
        
      }
      
    }
    i=i+1
  }
  
  #out = out[-removal.IDX,]
  # out = rbind(out,temp.dset)
  print(dim(out))
  
  write.csv(file=outfile,out,row.names=F,na="")
  return(out)
}
