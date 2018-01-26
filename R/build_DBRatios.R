###############################################
#### Build Theoratical Isootpe Ratio DB 
###############################################
#   This script builds the isotope ratio 
#   databases from the HMDB and LMDB databases.


library(xml2)
library(magrittr)
library(parallel)


#  Read the hmdb xml file:   download from www.hmdb.ca and unzip
xml  = read_xml("hmdb_metabolites.xml")
metabs = xml_children(xml)

cat(paste("    # Metabolites: ",length(metabs),"\n"))


formula = unique(unlist(formula))
formula = formula[-grep("R",formula)] #
formula = formula[-grep("\\)n",formula)]
formula = formula[-grep("Fe",formula)]
formula = formula[-grep("Hg",formula)]
formula

# Get the HMDB acessions
ID = lapply(metabs, function(x) xml_child(x,search=4))
ID = sapply(ID, function(x) xml_text(xml_contents(x)))

# Get the HMDB names
name = lapply(metabs, function(x) xml_child(x,search=6))
name = sapply(name, function(x) xml_text(xml_contents(x)))

# Get the HMDB formulae
formulae = lapply(metabs, function(x) xml_child(x,search=9))
formulae = sapply(formulae, function(x) xml_text(xml_contents(x)))
formulae[lapply(formulae, is.null)] = NA
formulae = unlist(formulae)


# Get the HMDB main class
mass = lapply(metabs, function(x) xml_child(x,search=11))
mass[sapply(mass,is.null)] = NA
mass = sapply(mass, function(x) xml_text(xml_contents(x)))


# Get the HMDB main class
class = lapply(metabs, function(x) xml_child(x,search=))
class = sapply(class, function(x) xml_text(xml_contents(x)))

# Clean up (large memory requirement)
rm(xml)
rm(metabs)
gc()



#Read the Lipid Maps data
lmdb = read.delim("LMDB.txt",header=T,check.names = F)
lmf = unique(as.character(lmdb$formula))

formula = c(lmf,formula)
formula = unique(formula)
formula = formula[-grep("D",formula)]


###### Run Emass to get masses and ratios: build target feature library
# you need to download emass.exe
start_time <- Sys.time()
K.DB = lapply(formula, function(x) shell(shQuote("emass "),input = paste(x,"K",sep=""),intern = T))
end_time <- Sys.time()
print(end_time - start_time)
save(K.DB,file="final_KDB.rda")

start_time <- Sys.time()
DIMER.DB = lapply(formula, function(x) shell(shQuote("emass "),input = paste("(",x,")2",sep=""),intern = T))
save(DIMER.DB,file="final_DIMERDB.rda")
end_time <- Sys.time()
print(end_time - start_time)

start_time <- Sys.time()
NH4.DB = lapply(formula, function(x) shell(shQuote("emass "),input = paste(x,"NH4",sep=""),intern = T))
save(NH4.DB,file="final_NH4DB.rda")
end_time <- Sys.time()
print(end_time - start_time)



#Fct to build libraries as R objects
prepDB = function(emass,maxPks){
  formula =  unlist(strsplit(emass[1]," "))[2]
  
  masses = array();
  intensities = array();
  for (i in 2:length(emass)) {
    line = as.numeric(unlist(strsplit(emass[i]," ")))
    masses[i - 1] = line[1]
    intensities[i - 1] = line[2]
    
  }
  intensities = intensities / sum(intensities)
  if(length(grep("Fe",formula) )>0){
    w = which(intensities == max(intensities))
    intensities = intensities[w:min(length(intensities),w+maxPks-1)]
    intensities = round(intensities,6)
  }else{
    intensities = intensities[1:min(length(intensities),maxPks)]
    intensities = round(intensities,6)
  }
  if(length(grep("Fe",formula))>0){
    masses = masses[w:min(length(masses),w+maxPks-1)]
    masses = round(masses,5)
  }else{
    masses = masses[1:min(length(masses),maxPks)]
    masses = round(masses,5)
  }
  delta = diff(masses)
  n = length(masses)
  #return(
  out = list(
    formula = as.character(formula),exactMass = masses[1], nPeaks = n,delta = delta,mass = masses,intensity = intensities
  )
  
  
  out
  #)
}

load("final_DB.rda")
load("final_KDB.rda")
load("final_DIMERDB.rda")
load("final_NH4DB.rda")

#Prepare the databases
DB_TFL = lapply(DB, function(x) prepDB(x,maxPks=2))
K_TFL = lapply(K.DB, function(x) prepDB(x,maxPks=2))
DIMER_TFL = lapply(DIMER.DB, function(x) prepDB(x,maxPks=2))
NH4_TFL = lapply(NH4.DB, function(x) prepDB(x,maxPks=2))

# Function to filter unwanted non organics or irrelevant atoms
filterAtoms = function(DB){
  #### Filter unwanted formulae
  formulae = unlist(lapply(DB,function(x) x$formula))
  g = grep("C[1-9]+H[1-9]+",formulae)
  formulae = formulae[g]
  DB = DB[g]
  
  #g = grep("D",formulae)
  #formulae = formulae[-g]
  #DB = DB[-g]

  g = grep("P[4-9]*",formulae)
  formulae = formulae[-g]
  DB = DB[-g]
  
  g = grep("F[4-9]{1}",formulae)
  formulae = formulae[-g]
  DB = DB[-g]
  
  g = grep("Br",formulae)
  formulae = formulae[-g]
  DB = DB[-g]
  
  g = grep("B",formulae)
  formulae = formulae[-g]
  DB = DB[-g]
  
  g = grep("Al",formulae)
  formulae = formulae[-g]
  DB = DB[-g]
  
  g = grep("As",formulae)
  formulae = formulae[-g]
  DB = DB[-g]
  
  g = grep("I",formulae)
  formulae = formulae[-g]
  DB = DB[-g]
  
  g = grep("Sb",formulae)
  formulae = formulae[-g]
  DB = DB[-g]
  
  g = grep("Se",formulae)
  formulae = formulae[-g]
  DB = DB[-g]
  DB
  
  g = grep("F17",formulae)
  formulae = formulae[-g]
  DB = DB[-g]
  DB
  
  g = grep("Cl",formulae)
  formulae = formulae[-g]
  DB = DB[-g]
  
}



DB_TFL = filterAtoms(DB_TFL)
KDB_TFL = filterAtoms(K_TFL)
DIMER_TFL = filterAtoms(DIMER_TFL)
NH4_TFL = filterAtoms(NH4_TFL)

set = DB_TFL
rats.DB = sapply(set, function(x) x$intensity[1]/x$intensity[2])
M.DB = sapply(set, function(x) x$exactMass)


# Visualize ratios
plot(log(M.DB),(rats.DB),xlab="Log Exact Mass",ylab= "M0/M1",main="M+H, Computed Ratios",pch=19,cex=0.75,
     col=rgb(0.0,0.1,0.1,0.55),xlim=c(3,10))


set = NH4_TFL
rats.DB = sapply(set, function(x) x$intensity[1]/x$intensity[2])
M.DB = sapply(set, function(x) x$exactMass)

points(log(M.DB),(rats.DB),xlab="Log Exact Mass",ylab= "M0/M1",main="M+NH4, Computed Ratios",pch=19,cex=0.25,
     col=rgb(0.1,0.8,0.1,0.25))



set = DIMER_TFL
rats.DB = sapply(set, function(x) x$intensity[1]/x$intensity[2])
M.DB = sapply(set, function(x) x$exactMass)
points(log(M.DB),(rats.DB),xlab="Log Exact Mass",ylab= "M0/M1",main="2M+H, Computed Ratios",pch=19,cex=0.25,
     col=rgb(0.8,0.1,0.1,0.25))


set = KDB_TFL
rats.DB = sapply(set, function(x) x$intensity[1]/x$intensity[2])
M.DB = sapply(set, function(x) x$exactMass)

points(log(M.DB),(rats.DB),pch=1,cex=0.25, col=rgb(0.1,0.8,0.8,0.25))



###########################################################
# Manually select points on boundary to build ratio limits.
set = DB_TFL
rats.DB = sapply(set, function(x) x$intensity[1]/x$intensity[2])
M.DB = sapply(set, function(x) x$exactMass)

plot(log(M.DB),(rats.DB),pch=19,cex=0.5, col=rgb(0.1,0.1,0.1,0.5))

identifyPch <- function(x, y = NULL, n = length(x), pch = 19,col, ...)
{
  xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
  sel <- rep(FALSE, length(x)); res <- integer(0)
  while(sum(sel) < n) {
    ans <- identify(x[!sel], y[!sel], n = 1, plot = FALSE, ...)
    if(!length(ans)) break
    ans <- which(!sel)[ans]
    points(x[ans], y[ans], pch = pch,col=col)
    sel[ans] <- TRUE
    res <- c(res, ans)
  }
  res
}
upperidx = identifyPch(log(M.DB),rats.DB,col=2)
loweridx = identifyPch(log(M.DB),(rats.DB),col=4)


#save(upperidx,file="upperidxNH4.rda")
#save(loweridx,file="loweridxNH4.rda")
#save(upperidx,file="upperidxK.rda")
#save(loweridx,file="loweridxK.rda")
#save(upperidx,file="upperidx2M.rda")
#save(loweridx,file="loweridx2M.rda")
#save(upperidx,file="upperidxMH.rda")
#save(loweridx,file="loweridxMH.rda")
#save(loweridx,file="loweridx2.rda")
