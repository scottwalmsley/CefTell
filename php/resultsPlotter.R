# Routine to connect to SQL database and plot MSMS 'hit's for validataion

########################
## connect to mysql

EXPTAG = "baldb"  # experiment tag
library(RMySQL)



nistlib = read.delim("nistlib.tsv",sep="\t",header=F)
colnames(nistlib) = c("library","LIBID","names","formula","MW", "precursorMZ","casNO","ion","nPeaks","mzPeaks","intPeaks")


nistLIBID = as.numeric(nistlib$LIBID)
con <- dbConnect(MySQL(),
 user = 'root',    password = '######', host = '######',dbname='baldb');

sql = paste("select * from raw_data_",EXPTAG," where mslevel = 2",sep="")
raw = dbGetQuery(conn=con,sql)
#raw = read.delim(paste("raw_data_",EXPTAG,".tsv",sep=""),sep=" ")
write.table(raw,file=paste("raw_data_",EXPTAG,".tsv",sep=""))

##############################
## get filtered results table
sql = paste("select *, max(prob) as prob from ", EXPTAG ," where IDRANK=1 and dotP > 800 and abs(ppm) <50 and extraction like 'LIPID%' group by name,cev order by polarity,cev,name",
		sep="");

shortList = dbGetQuery(conn=con,sql)
names = shortList$name
fh = shortList$fh
fh = sub(".MGF.TSV",".txt",fh);
sample = sub(".MGF.TSV","",fh);
dim(shortList)

scanID = shortList$scan
libID = shortList$libID
libMz = shortList$libMZ
lib = shortList$lib
formula = shortList$formula
polarity = shortList$polarity
CeV = shortList$cev
ion = shortList$ion
rt = shortList$rt
mz = shortList$mz
extraction = shortList$extraction
score = shortList$score
dot = shortList$dotP
P = shortList$prob
index = seq(1,length(names), by=1)


shortList = cbind(index,shortList)
grep("Lipid",shortList$lib)
unique(shortList$lib)
grep("nist_msms2",shortList$lib)


#################################################################
#### loop to get data and print it
##################################################################
#pdf(file=paste(EXPTAG,"_SHORT_111117_HILIC_SPECLIB.pdf",sep=""))
pdf(file=paste(EXPTAG,"_SHORT_111117_LIPID_SPECLIB.pdf",sep=""))

filterMzbyRef = function(MZ,I,MZr,ppm){
   ret.MZ = array()
   ret.I = array()

   for(i in 1:length(MZr)){
      m = MZr[i] 
      mzTol = ppm*m/1e6
      w = which(MZ > m-mzTol & MZ < m+mzTol)
      
      if(length(w) ==0){
          ret.MZ = c(ret.MZ,0)
          ret.I = c(ret.I,0)
      }
      if( length(w)== 1){ 
	   ret.MZ = c(ret.MZ,MZ[w])
         ret.I = c(ret.I,I[w])
      }
      if(length(w) > 1){
          ret.MZ = c(ret.MZ, MZ[which(I[w] == max(I[w]))])
          ret.I = c(ret.I, I[which(I[w] == max(I[w]))])

      }
     
   }  
   ret.MZ = ret.MZ[-1]
   ret.I = ret.I[-1]
   ret = data.frame(mz = ret.MZ,I = ret.I)

}


par(mfrow=c(3,1))
par(mar=c(4,1,3,1))

IDX = array()
j=1;



for(i in 1:length(names)){


   ##############################
   ## get a single result SCAN


  # sql = paste("SELECT * from raw_data_",EXPTAG," where fh like '",fh[i],"' and scan = ",scanID[i],sep=""); 
   #res = dbGetQuery(conn=con,sql)

   res = raw[which(raw$fh ==fh[i] & raw$scan == scanID[i]),]

I = as.numeric(unlist(strsplit(as.character(res$intData), " ")))
MZ = as.numeric(unlist(strsplit(as.character(res$mzData), " ")))
I = I / max(I) * 999




if(grepl("Lipid",lib[i])){	
	#sql = paste("SELECT * from nistlib where library like '",lib[i],"' and LIBID = '",names[i],"'",sep="");
     
}else{

	#sql = paste("SELECT * from nistlib where library like '",lib[i],"' and LIBID = '",libID[i],"'",sep="");
	if(grepl("msms2",lib[i])){
          w = which(nistlib$library == lib[i] & nistlib$names == names[i] & nistlib$formula == shortList$formula[i])
          w = w[1]
	}else{
         w = which(nistlib$library == lib[i] &  nistlib$LIBID == libID[i])
       #w = which(nistlib$LIBID == libID[i])
     }

}


#LIBID = as.character(nistlib$LIBID)
#w = which(LIBID == libID[i])



#entry = dbGetQuery(conn=con,sql)
entry = nistlib[w,]



if(grepl("Lipid",lib[i])){
   #sql = paste("SELECT * from nistlib where library like '",lib[i],"' and LIBID = '",names[i],"'",sep="");
   #get  = dbGetQuery(conn=con,sql)
   libID[i] = shortList$libID[i]= as.character(entry$LIBID)
   formula[i] = shortList$formula[i] = as.character(entry$formula)
   lib[i] = shortList$lib[i] = as.character( entry$library)
   libMz[i] = shortList$libMZ[i] = as.numeric(entry$MW)

}



Ir = as.numeric(unlist(strsplit(as.character(entry$intPeaks), " ")))
MZr = as.numeric(unlist(strsplit(as.character(entry$mzPeaks), " ")))

err =  round((mz[i] - libMz[i]) / libMz[i] * 1e6,1)
errDa = round((mz[i] - libMz[i]),3)


getFilt = filterMzbyRef(MZ,I,MZr,100)

dMZ = getFilt$mz
dI = getFilt$I

distance = round(sqrt(sum((dI-Ir)^2)))

#MZ = dMZ
#I = dI


#######################
## now plot 3 / page


 
   IDX[j] = i; 
   j=j+1

   main =  paste(index[i],".  ",names[i],"\n",paste("Score=",score[i]," Dot=",dot[i]," prob=",P[i],sep=""))
   
   plot(main =main,ylim=c(-1100,1100),xlim=c((min(MZr)-100),mz[i]+10),xlab="",MZ,I,type="h",lwd=2,col=2,cex.main=1)

   abline(h=0,lwd=1)
   lines(MZr,-Ir,type="h",lwd=2,col=4)
 
   txt  = paste(sample[i],extraction[i],"\nmz = ",mz[i],";  rt=",rt[i], "\ncev=",CeV[i],"| polarity=", polarity[i],ion[i])
   txt2 = paste(entry$library,":",entry$LIBID,"\nCAS=",entry$casNO,"\nFormula:",formula[i],"\nlibMz= ",libMz[i])
   txt3 = paste("err = ",paste(err,"ppm , ",errDa," Da",sep=""))
   txt4 = paste("dist = ",distance,sep="")
   cx = 0.75
   text(min(MZr-100),700,txt,cex=cx,adj=c(0,0))
   text(min(MZr-100),-700,txt2,cex=cx,adj=c(0,0))
   text(min(MZr-100),-999,txt3,cex=cx,adj=c(0,0))
   text(min(MZr-100),100, txt4,cex=0.8,adj=c(0,0))

}


write.csv(shortList[IDX,],file=paste(EXPTAG,"_SPECLIB_LIPID_SHORTLIST.csv",sep=""))

dev.off()
dbDisconnect(con)



