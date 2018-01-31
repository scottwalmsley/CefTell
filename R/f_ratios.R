

##########################################################
##### Filter on the predictedsmooth curve
##########################################################
ratio.list = ratios[,]
mz.list = ratios[,3:20]

#ratio.list = ratios[,21:38]
#mz.list = ratios[,3:20]
iontypes = colnames(cef.list[[1]]$sumIonDataByCompound$MZ)
out.ratio.list = list()
M = log(ratios$M)
out.ratio.mat = matrix(nrow = length(M),ncol=18)

t.M = mz.list


#### Load models:
## Dimer: 
load("lwupper2M.rda")
load("lwlower2M.rda")
dimer.upper = lw.upper
dimer.lower = lw.lower
#load("lwDIMER.rda")
#dimer.mod = lw

load("lwupperNH4.rda")
load("lwlowerNH4.rda")
NH4.upper = lw.upper
NH4.lower = lw.lower
#load("lwNH4.rda")
#nh4.mod = lw

load("lwupperK.rda")
load("lwlowerK.rda")
K.upper = lw.upper
K.lower = lw.lower
#load("lwK.rda")
#k.mod = lw


load("lwupperMH.rda")
load("lwlowerMH.rda")
M.upper = lw.upper
M.lower = lw.lower
#load("lwMH.rda")
#mh.mod = lw


for(k in 1:18){
  ion = iontypes[k]
  if(ion == "M+Na"){
      #t.M = M+22
      lw.upper = M.upper
      lw.lower = M.lower
      #lw = mh.mod
 
  }
  if(ion == "M+2H"){
      #t.M = M
      lw.upper = M.upper
      lw.lower = M.lower
      #lw = mh.mod      
  }
  if(ion == "2M+Na"){
      #t.M = (M*2)+22
      lw.upper = dimer.upper
      lw.lower = dimer.lower
      #lw = dimer.mod
  }
  if(ion == "M+3H"){
      #t.M = M
      lw.upper = M.upper
      lw.lower = M.lower
      #lw = mh.mod
  }
  if(ion == "M+H"){
      #t.M = M
      lw.upper = M.upper
      lw.lower = M.lower
      #lw = mh.mod
  }
  if(ion == "M+HNa"){
      #t.M = M
      lw.upper = M.upper
      lw.lower = M.lower
      #lw = mh.mod
  }
  if(ion == "M+K"){
      #t.M = M-39
      lw.upper = K.upper
      lw.lower = K.lower
      #lw = k.mod
  }
  if(ion == "2M+K"){
      #t.M = (M*2)+39
      lw.upper = dimer.upper
      lw.lower = dimer.lower
      #lw = dimer.mod
  }
  if(ion == "2M+H"){
      #t.M = M*2
      lw.upper = dimer.upper
      lw.lower = dimer.lower
      #lw = dimer.mod
  }
  if(ion == "M+NH4"){
      #t.M = M+18
      lw.upper = NH4.upper
      lw.lower = NH4.lower
      #lw = nh4.mod
  }
  if(ion == "2M+NH4"){
      #t.M = (M*2)+18
      lw.upper = dimer.upper
      lw.lower = dimer.lower
      #lw = dimer.mod
  }
  if(ion == "M+KNa"){
      #t.M = M+22
      lw.upper = K.upper
      lw.lower = K.lower
      #lw = k.mod
  }   
  if(ion == "M+HK"){
      #t.M = M+39
      lw.upper = K.upper
      lw.lower = K.lower
      #lw = k.mod
  }   
  if(ion == "M+H5N"){
      #t.M = M+10
      lw.upper = M.upper
      lw.lower = M.lower
      #lw = mh.mod
  } 
  if(ion == "M+H4NNa"){
      #t.M = M+22
      lw.upper = M.upper
      lw.lower = M.lower
      #lw = mh.mod
  } 
  if(ion == "M+NH4+NH4"){
      #t.M = M+18
      lw.upper = NH4.upper
      lw.lower = NH4.lower
      #lw = nh4.mod
  } 
  if(ion == "M+2Na"){
      #t.M = M+22
      lw.upper = M.upper
      lw.lower = M.lower
      #lw = mh.mod
  } 
  if(ion == "M+H4KN"){
      #t.M = M-56
      lw.upper = K.upper
      lw.lower = K.lower
      #lw = k.mod
  } 
  if(ion == "M+2K"){
      #t.M = M+39
      lw.upper = K.upper
      lw.lower = K.lower
      #lw = k.mod
  } 



  print(ion)
  test.out = NULL
  rat = ratio.list[,k]
  #t.M = log(mz.list[,k])


  for(i in 1:length(M)){
   # filter.upper = filter.lower = array(dim=4,data=0);
    filter =  array(dim=2,data=0)
    #m = t.M[i]
    #mz = ratios$mz.1[i]
    mz = log(mz.list[i,k])
    rt = ratios$RT[i]
    ratio = rat[i]
    m = log(mz.list[i,k])
    ### Upper bounds 
    
    if(!is.na(ratio)){
       pred.upper = predict(lw.upper,m)
       pred.lower = predict(lw.lower,m)
       #pred = predict(lw,m)    

       uy = pred.upper$y + (pred.upper$y*0.3)
       ly = pred.lower$y - (pred.lower$y*0.3)

       #uy = (pred$y+(pred$y*0.1))
       #ly = (pred$y-(pred$y*0.1))

      if( ratio > (uy) ){
        filter[2] = 1
      }
      if( ratio < (ly)){
        filter[1] = 1 
      }
      
    }
    if(is.na(ratio)){
      filter = array(dim=2)
    }
    line = c(exp(M[i]),mz,rt,ratio, ly,uy,filter)
    test.out = rbind(test.out,line)

  }
  
   colnames(test.out) = c("M","mz", "RT","Ratio","lowerLimit","upperLimit",
                          "l","u")

   
  rownames(test.out) = NULL
  
  
  result= as.data.frame(test.out)
  colnames(result)
  
  out.ratio.list[[k]] = result
  
}



uppers = lowers = matrix(nrow = length(M),ncol=18)
for(i in 1:18){
  uppers[,i] = out.ratio.list[[i]][,7]
  lowers[,i] = out.ratio.list[[i]][,8]
  
}
filt = cbind(lowers,uppers)



filt[which(filt==1)]=NA
filt[which(filt==0)]=1

iontypes

r = c(1:2,4,6,8,10)
#r = c(1:2,4,8)
r = 1:16
r = 1:16
remo = apply(filt[,r],1,sum,na.rm=T)

remo2 = apply(filt[,r+16],1,sum,na.rm=T)

#remo = apply(filt[,1:18],1,sum,na.rm=T)

wu = which(remo == 0)
length(remo[wu])
wu2 = which(remo2 == 0 ) 
length(remo2[wu2])

w = c(wu,wu2)
plot(result$RT, M,pch=19)
points(result$RT[w], M[w],col=2,pch=19,cex=0.75)


plot(x=result$RT[-w], y=result$M[-w],cex=0.75,pch=19,col=rgb(0.1,0.1,0.1,0.1),xlab="Retention Time",ylab="M/Z")
points(x=result$RT[w], y=result$M[w],cex=0.75,pch=19,col=rgb(0.8,0.1,0.1,0.82))





dim(result)

trim = cef.list[[1]]


dset = data.frame("ID" = trim$ID[-w],
                  "Mass"=trim$ratios$M[-w],
                  "RT"=trim$ratios$RT[-w],
                  "Formula"=trim$Formula[-w],
                  "Name"=trim$Name[-w],
                  "HMDB"=trim$HMDB[-w],
                  "LMDB"=trim$LMDB[-w])

colnames(trim$ratios) = c("Mass","RT",c(iontypes,iontypes))
rownames(trim$ratios) = paste(seq(1,12213,by=1),"_",trim$Name,sep="")
options(java.parameters = "-Xmx4096m")
write.xlsx2(trim$ratios[-w,1:20], file="final_BALDB.xlsx", sheetName="Observed Adduct MZ",
              col.names=TRUE, row.names=TRUE, append=FALSE)
write.xlsx2(trim$ratios[-w,c(1:2,21:38)], file="final_BALDB.xlsx", sheetName="Observed Adduct Ratio",
              col.names=TRUE, row.names=TRUE, append=TRUE)

save(file="trimmed_dset.rda",dset)
