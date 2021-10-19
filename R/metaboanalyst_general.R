# internal variables and functions not to be modified by users
# This is only for web version
.on.public.web <- FALSE; # only TRUE when on metaboanalyst web server

# note, this is usually used at the end of a function
# for local, return itself; for web, push to global environment
.set.mSet <- function(mSetObj=NA){
  if(.on.public.web){
    mSet <<- mSetObj;
    return (1);
  }
  return(mSetObj);
}

.get.mSet <- function(mSetObj=NA){
  if(.on.public.web){
    return(mSet)
  }else{
    return(mSetObj);
  }
}

#'Constructs a dataSet object for storing data
#'@description This functions handles the construction of a mSetObj object for storing data for further processing and analysis.
#'It is necessary to utilize this function to specify to MetaboAnalystR the type of data and the type of analysis you will perform.
#'@usage InitDataObjects(data.type, anal.type, paired=FALSE)
#'@param data.type The type of data, either list (Compound lists), conc (Compound concentration data),
#'specbin (Binned spectra data), pktable (Peak intensity table), nmrpeak (NMR peak lists), mspeak (MS peak lists),
#'or msspec (MS spectra data)
#'@param anal.type Indicate the analysis module to be performed: stat, pathora, pathqea, msetora, msetssp, msetqea, ts,
#'cmpdmap, smpmap, or pathinteg
#'@param paired Indicate if the data is paired or not. Logical, default set to FALSE
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import methods

InitDataObjects <- function(data.type, anal.type, paired=FALSE){

  dataSet <- list();
  dataSet$type <- data.type;
  dataSet$design.type <- "regular"; # one factor to two factor
  dataSet$cls.type <- "disc"; # default until specified otherwise
  dataSet$format <- "rowu";
  dataSet$paired <- paired;
  analSet <- list();
  analSet$type <- anal.type;
  mSetObj <- list();
  mSetObj$dataSet <- dataSet;
  mSetObj$analSet <- analSet;
  mSetObj$imgSet <- list();
  mSetObj$msgSet <- list(); # store various message during data processing
  mSetObj$msgSet$msg.vec <- vector(mode="character");     # store error messages
  mSetObj$cmdSet <- vector(mode="character"); # store R command

  .init.global.vars(anal.type);
  print("MetaboAnalyst R objects initialized ...");
  return(.set.mSet(mSetObj));
}

.init.global.vars <- function(anal.type){

  # other global variables
  msg.vec <<- "";
  err.vec <<- "";

  # for mummichog
  peakFormat <<- "mpt"
  mumRT.type <<- "NA";

  anal.type <<- anal.type;

  if(.on.public.web){
    # disable parallel prcessing for public server
    library(BiocParallel);
    register(SerialParam());
  }else{
    if("stat" %in% anal.type | "msetqea" %in% anal.type | "pathqea" %in% anal.type | "roc" %in% anal.type)
      # start Rserve engine for Rpackage
      load_Rserve();
  }

  # plotting required by all
  Cairo::CairoFonts(regular="Arial:style=Regular",bold="Arial:style=Bold",
                    italic="Arial:style=Italic",bolditalic = "Arial:style=Bold Italic",symbol = "Symbol")
}

#' Adds an error message
#'@description The error message will be printed in all cases.
#'Used in higher functions.
#'@param msg Error message to print
#'@export
AddErrMsg <- function(msg){
  err.vec <<- c(err.vec, msg);
  print(msg);
}

# read binary qs files
# sub.dir is sub folder, leave NULL is under main lib folder
.get.my.lib <- function(filenm, sub.dir=NULL){

  if(!is.null(sub.dir)){
    sub.dir <- paste0(sub.dir, "/");
  }
  if(.on.public.web){
    lib.path <- paste0("../../libs/", sub.dir, filenm);
    print(paste("loading library:", lib.path));
    return(qs::qread(lib.path));
  }

  lib.download <- FALSE;

  print(filenm)

  if(!file.exists(filenm)){
    lib.download <- TRUE;
  }else{
    # time <- file.info(filenm)
    # diff_time <- difftime(Sys.time(), time[,"mtime"], unit="days")
    # if(diff_time>90){
    #   lib.download <- TRUE;
    # }
  }

  lib.url <- paste0("https://www.metaboanalyst.ca/resources/libs/", sub.dir, filenm);
  # Deal with curl issues
  if(lib.download){
    tryCatch(
      {
        download.file(lib.url, destfile=filenm, method="curl")
      }, warning = function(w){ print() },
      error = function(e) {
        print("Download unsucceful. Ensure that curl is downloaded on your computer.")
        print("Attempting to re-try download using libcurl...")
        download.file(lib.url, destfile=filenm, method="libcurl")
      }
    )
  }

  lib.path <- filenm;

  # Deal w. corrupt downloaded files
  tryCatch({
    my.lib <- qs::qread(lib.path); # this is a returned value, my.lib never called outside this function, should not be in global env.
    print("Loaded files from MetaboAnalyst web-server.")
  },
  warning = function(w) { print() },
  error = function(err) {
    print("Reading data unsuccessful, attempting to re-download file...")
    tryCatch({
      download.file(lib.url, destfile=filenm, method="curl")
      my.lib <- qs::qread(lib.path);
      print("Loaded necessary files.")
    },
    warning = function(w) { print() },
    error = function(err) {
      print("Loading files from server unsuccessful. Ensure curl is downloaded on your computer.")
    }
    )
  })
  return(my.lib);
}

fast.write.csv <- function(dat, file, row.names=TRUE){
  tryCatch(
    {
      if(is.data.frame(dat)){
        # there is a rare bug in data.table (R 3.6) which kill the R process in some cases
        data.table::fwrite(dat, file, row.names=row.names);
      }else{
        write.csv(dat, file, row.names=row.names);
      }
    }, error=function(e){
      print(e);
      write.csv(dat, file, row.names=row.names);
    }, warning=function(w){
      print(w);
      write.csv(dat, file, row.names=row.names);
    });
}

GetShortNames<-function(nm.vec, max.len= 45){
  new.nms <- vector(mode="character", length=length(nm.vec));
  for(i in 1:length(nm.vec)){
    nm <- nm.vec[i];
    if(nchar(nm) <= max.len){
      new.nms[i] <- nm;
    }else{
      wrds <- strsplit(nm, "[[:space:]]+")[[1]];
      new.nm <- "";
      if(length(wrds)>1){
        for(m in 1:length(wrds)){
          wrd <- wrds[m];
          if(nchar(new.nm)+4+nchar(wrd) <= max.len){
            new.nm <- paste(new.nm, wrd);
          }else{
            new.nms[i] <- paste (new.nm, "...", sep="");
            break;
          }
        }
      }else{
        new.nms[i] <- paste (substr(nm, 0, 21), "...", sep="");
      }
    }
  }
  return(new.nms);
}

# return heat colors of given length
GetMyHeatCols <- function(len){
  if(len > 50){
    ht.col <- c(substr(heat.colors(50), 0, 7), rep("#FFFFFF", len-50));
  }else{
    # reduce to hex by remove the last character so HTML understand
    ht.col <- substr(heat.colors(len), 0, 7);
  }
}

########################
########################
#### Load libraries ####
########################
########################

# Load ggplot2
load_ggplot <- function(){
  suppressMessages(library(ggplot2))
}
