
#' clean.formulae
#'
#' Cleans Formulae of unwanted entries, and trims the list down to unique entries.
#'
#' @param formulae
#'
#' @return the character vector of clean formulae
#' @export
#'
#'
clean.formulae = function(formulae){
   formulae = as.character(formulae)
   formulae = unique(formulae)
   formulae = formulae[-grep("\\[",formulae)] ## remove one odd formula
   formulae = formulae[-grep("D",formulae)]  ## remove deuterated entries
   formulae = formulae[-grep("R",formulae)] ## remove entries with R group annotations
   formulae
}




#' run.emass
#'
#' Runs Emass for prediciton of isotopic abundances
#'
#' @param formula  The input Lewis formula
#' @param maxPks  Maximum number of desired peaks
#' @export
#' @return Returns isotopic (mass and intensity) values
#' @keywords internal
run.emass = function(formula, maxPks = 6) {
   print(formula)
   masses = array();
   intensities = array();

   emass = shell(shQuote("emass "),input = paste(formula,",0"),intern = T)


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


#' build.DB.emass
#'
#' Builds the mass and intensity database using emass.
#'
#' @param k number of isotopes
#' @param formula The input list of formulae
#' @export
#' @return The target feature library of emass predicted isotopic mass and intensities in DB object form.
build.DB.emass = function(k,formula) {
  #mclapply = getSys();
  tmp = lapply(formula, function(x)
    run.emass(x,k))

  ####cleanup:   get rid of small ions, wierd ions
  DB = list()
  i = 1;
  for (iso in tmp) {

    if ((iso$nPeaks > 2) & (iso$delta[1] < 1.2)) {
      DB[[i]] = iso
      i = i + 1
    }
  }
  DB
}




#' get.DB.dataframe
#'
#' A helper function to convert list of DB objects into a dataframe.
#'
#' @param DB The database (list of DB objects)
#' @param n Limit to n number of peaks in the output for each isotope
#'
#' @return  The database in dataframe format.
#' @export
#'
get.DB.dataframe = function(DB,n) {
  db.formula = unlist(sapply(DB,function(x)
    unlist(x$formula)))
  db.exactMass = unlist(sapply(DB,function(x)
    unlist(x$mass[1])))
  db.IsotopeIntensities = lapply(DB,function(x)
    unlist(x$intensity))
  db.IsotopeMasses = lapply(DB,function(x)
    unlist(x$mass))

  out = data.frame(formula = db.formula,exactMass = db.exactMass)

  int = data.frame(Io = unlist(sapply(db.IsotopeIntensities, function(x)
    unlist(x[1]))))
  mass = data.frame(Mo = unlist(sapply(db.IsotopeMasses, function(x)
    unlist(x[1]))))

  for (i in 2:n) {
    M = unlist(sapply(db.IsotopeMasses,function(x)
      unlist(x[i])))
    I = unlist(sapply(db.IsotopeIntensities,function(x)
      unlist(x[i])))
    int = cbind(int,I)
    mass = cbind(mass,M)
    colnames(mass)[i] = paste("M",i - 1,sep = "")
    colnames(int)[i] = paste("I",i - 1,sep = "")
  }

  cbind(out,mass,int)
}




#' write.DB
#'
#' Writes the Target Feature Library in tab deliminted form.
#'
#' @param file File name to save as
#' @param DB The list of DB objects.
#' @export
write.DB = function(file,DB) {
  out = "";
  for (obj in DB) {
    header = paste("#",obj$formula,obj$nPeaks);
    line = paste(
      obj$formula,obj$exactMass,obj$nPeaks,paste(obj$mass),paste(obj$intensity),sep =
        "\t"
    )
    out = c(out,header,line)
  }
  write(file = file,out,sep = "")
}



#' read.DB
#'
#' Reads the Target Feature Library from a tab delimited file.
#'
#' @param file
#'
#' @return Returns list of DB objects.
#' @export
read.DB  = function(file) {
  d = list()

  data = read.table(file);
  uf = unique(as.character(data[,1]))

  for (i in 1:length(uf)) {
    sub = data[as.character(data[,1]) == uf[i],]
    d[[i]] = list(
      formula = as.character(sub[1,1]),exactMass = sub[1,2],nPeaks = sub[1,3],mass = sub[,4],intensity = sub[,5]
    )

  }
  d
}


#' generate.decoy.DB
#'
#' Produce a decoy feature library by computing around mass (ppm) and %intensity.
#'
#' @param DB  The input database
#' @param jM  Jitter masses (ppm)
#' @param jI  Jitter intensities (0-1)
#' @param k   k number of decoys produced per input spectra
#' @export
#' @return  Returns the decoy feature library in DB object form.
#' @keywords internal
generate.decoy.DB = function(DB,jM,jI,k) {
  DFL = list()
  i = 1

  for (iso in DB) {
    for (n in 1:k) {
      formula = paste("decoy_",k,"_",iso$formula,sep = "")
      dec_mass = sapply(iso$mass, function(x)
        jitter(x, amount = ppmErr(x, jM)))
      dec_int = sapply(iso$intensity, function(x)
        jitter(x, amount = x * jI))
      dec_delta = diff(dec_mass)

      DFL[[i]] = list(
        formula = formula, nPeaks = iso$nPeaks, delta = dec_delta, mass = dec_mass,intensity = dec_int
      )
      i = i + 1
    }
  }
  DFL
}


#' generate.decoy.spectra
#'
#' Produce a decoy database by computing around mass (ppm) and %intensity.
#'
#' @param DB  The input target feature library
#' @param jM  Jitter the mass in ppm
#' @param dI  Second best delta intensities
#' @param k   Number of decoys produced per input spectra
#' @export
#' @return  Returns the decoy feature library in DB object form.
generate.decoy.spectra = function(DB,jM,dI,k) {

   DFL = list()
   i = 1

   # Generate for each isotope cluster encountered
   for (iso in DB) {

      for (n in 1:k) {  # k copies

         formula = paste("decoy_",k,"_",iso$formula,sep = "")


         ### Random uniform noise of mass deltas:  generate mass
         dec_mass = sapply(iso$mass, function(x)
            jitter(x, amount = ppmErr(x, jM)))


         ### Normal distributed intensities of isotope intensity deltas: generate intensities
         dec_int = NULL

         for(peak in iso$int){
            czech = -1;  #### check number for the while loop below
            while(czech == -1){
               jit.dI = peak+ rnorm(1,mean(dI),sd(dI));

               if((peak+jit.dI) > 0.00001){   ### a reasonable lower bounds for relative isotope intensities
                  dec_int = c(dec_int,peak+jit.dI);
                  czech = 0;
               }else{
                  jit.dI = rnorm(1,mean(dI),sd(dI));
                  czech=-1
               }
            }


         }

         dec_delta = diff(dec_mass)
         dec_int = dec_int/sum(dec_int)

         DFL[[i]] = list(
            formula = formula, nPeaks = iso$nPeaks, delta = dec_delta, mass = dec_mass,intensity = dec_int)

         i = i + 1
      }
   }

   DFL
}



#' ppmErr
#'
#' Compute the ppmErr for a given mass and expected ppm tolerance
#'
#' @param M the input mass
#' @param ppmTol the ppm tolerance window
#'
#' @return  the absolute ppm tolerance
#' @export
#' @keywords internal
ppmErr = function(M,ppmTol){
   return( ppmTol * M/1e6)
}

