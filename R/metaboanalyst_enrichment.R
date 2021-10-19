#'Save compound name for mapping
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param qvec Input the vector to query
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
Setup.MapData <- function(mSetObj=NA, qvec){
  mSetObj <- .get.mSet(mSetObj);
  mSetObj$dataSet$cmpd <- qvec;
  return(.set.mSet(mSetObj));
}

#'Given a list of compound names or ids, find matched name or ids from selected databases
#'@description Given a list of compound names or ids
#'find matched name or IDs from selected databases
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects).
#'@param q.type Input the query type, "name" for compound names, "hmdb" for HMDB IDs, "kegg" for KEGG IDs, "pubchem"
#'for PubChem CIDs, "chebi" for ChEBI IDs, "metlin" for METLIN IDs, and "hmdb_kegg" for a both KEGG and HMDB IDs.
#'@param hmdb Logical, T to cross reference to HMDB, F to not.
#'@param pubchem Logical, T to cross reference to PubChem, F to not.
#'@param chebi Logical, T to cross reference to CheBI, F to not.
#'@param kegg Logical, T to cross reference to KEGG, F to not.
#'@param metlin Logical, T to cross reference to MetLin, F to not.
#'@param lipid Logical, if features are lipids (T), a different database will be used for
#'compound matching.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

CrossReferencing <- function(mSetObj=NA, q.type, hmdb=T, pubchem=T,
                             chebi=F, kegg=T, metlin=F, lipid=F, model = T){

  mSetObj <- .get.mSet(mSetObj);

  # record the filter for 8 major databases
  mSetObj$return.cols <- c(hmdb, pubchem, chebi, kegg, metlin);
  mSetObj$lipid.feats <- lipid

  # record all the data
  if(!exists("name.map", where = mSetObj)){
    mSetObj$name.map <- list();
  }

  # distribute job
  mSetObj$dataSet$q.type <- q.type;

  if(.on.public.web){
    .set.mSet(mSetObj);
    MetaboliteMappingExact(mSetObj, q.type, lipid, model);
    mSetObj <- .get.mSet(mSetObj);
  }else{
    mSetObj <- MetaboliteMappingExact(mSetObj, q.type, lipid, model);
  }

  # do some sanity check
  todo.inx <- which(is.na(mSetObj$name.map$hit.inx));
  if(length(mSetObj$name.map$hit.inx) == 0){
    mSetObj$msgSet$nmcheck.msg <- c(0, "No hits found for the given compound ID. Please make
                                    sure that correct compound IDs or common compound names are used.");
  }else if(length(todo.inx)/length(mSetObj$name.map$hit.inx) > 0.5){
    mSetObj$msgSet$nmcheck.msg <- c(0, "Over half of the compound IDs could not be matched to our database. Please make
                                    sure that correct compound IDs or common compound names are used.");
  }else if (length(todo.inx) > 15){
    mSetObj$msgSet$nmcheck.msg <- c(2, "There are >15 compounds without matches. You can either proceed or if necessary, update these compound IDs and upload again.");
  }else{
    mSetObj$msgSet$nmcheck.msg <- c(1, "Name matching OK, please inspect (and manual correct) the results then proceed.");
  }

  # if(!.on.public.web){
  #   print(mSetObj$msgSet$nmcheck.msg)
  #
  #   if(length(todo.inx) == length(mSetObj$name.map$hit.inx)){
  #     AddErrMsg("Name matching failed! Please make sure that correct standardized feature names are used!")
  #     return(0)
  #   }
  # }

  return(.set.mSet(mSetObj));
}

#'Mapping from different metabolite IDs
#'@description For compound names to other ids, can do exact or approximate matches
#'For other IDs, except HMDB ID, all others may return multiple/non-unique hits
#'Multiple hits or non-unique hits will allow users to manually select
#'@param mSetObj Input the name of the created mSetObj.
#'@param q.type Inpute the query-type, "name" for compound names, "hmdb" for HMDB IDs, "kegg" for KEGG IDs, "pubchem"
#'for PubChem CIDs, "chebi" for ChEBI IDs, "metlin" for METLIN IDs, and "hmdb_kegg" for a both KEGG and HMDB IDs.
#'@param lipid Boolean, if features are lipids, a different database will be used for
#'compound matching.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
MetaboliteMappingExact <- function(mSetObj=NA, q.type, lipid = F, model = T){

  mSetObj <- .get.mSet(mSetObj);

  qvec <- mSetObj$dataSet$cmpd;

  # variables to record results
  hit.inx <- vector(mode='numeric', length=length(qvec)); # record hit index, initial 0
  names(hit.inx) <- qvec;
  match.values <- vector(mode='character', length=length(qvec)); # the best matched values (hit names), initial ""
  match.state <- vector(mode='numeric', length=length(qvec));  # match status - 0, no match; 1, exact match; initial 0

  if(model){
    cmpd.db <- .get.my.lib("model_metabolites_master.qs");
  }else{
    # cmpd.db <- .get.my.lib("master_cpds_microbiomenet.qs");
    cmpd.db <- .get.my.lib("master_compound_db.qs");
  }

  if(q.type == "hmdb"){
    n <- 5 # Number of digits for V3 of HMDB
    hmdb.digits <- as.vector(sapply(cmpd.db$hmdb, function(x) strsplit(x, "HMDB", fixed=TRUE)[[1]][2]))
    hmdb.v3.ids <- paste0("HMDB", substr(hmdb.digits, nchar(hmdb.digits)-n+1, nchar(hmdb.digits)))
    hit.inx.v3 <- match(tolower(qvec), tolower(hmdb.v3.ids));
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$hmdb));
    hit.inx[is.na(hit.inx)] <- hit.inx.v3[is.na(hit.inx)]
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
  }else if(q.type == "pubchem"){
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$pubchem));
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
  }else if(q.type == "chebi"){
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$chebi));
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
  }else if(q.type == "metlin"){
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$metlin));
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
  }else if(q.type == "kegg"){
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$kegg));
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
  }else if(q.type == "name"){

    # first find exact match to the common compound names
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$name));
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;

    if(!model){

      syn.db <- .get.my.lib("master_syn_nms.qs")
      # syn.db <- .get.my.lib("master_syns_microbiomenet.qs")

      syns.list <-  syn.db$syns.list;
      todo.inx <- which(is.na(hit.inx));

      if(length(todo.inx) > 0){
        for(i in 1:length(syns.list)){
          syns <-  syns.list[[i]];
          hitInx <- match(tolower(qvec[todo.inx]), tolower(syns));

          hitPos <- which(!is.na(hitInx));
          if(length(hitPos)>0){
            # record matched ones
            orig.inx<-todo.inx[hitPos];
            hit.inx[orig.inx] <- i;
            # match.values[orig.inx] <- syns[hitInx[hitPos]];  # show matched synnames
            match.values[orig.inx] <- cmpd.db$name[i];    # show common name
            match.state[orig.inx] <- 1;

            # update unmatched list
            todo.inx<-todo.inx[is.na(hitInx)];
          }
          if(length(todo.inx) == 0) break;
        }
      }
    }
  }

  mSetObj$name.map$query.vec <- qvec;
  mSetObj$name.map$hit.inx <- hit.inx;
  mSetObj$name.map$hit.values <- match.values;
  mSetObj$name.map$match.state <- match.state;

  return(.set.mSet(mSetObj));
}

#'Creates the mapping result table
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@export
CreateMappingResultTable <- function(mSetObj=NA, model = T){

  mSetObj <- .get.mSet(mSetObj);

  lipid <- mSetObj$lipid.feats

  qvec <- mSetObj$dataSet$cmpd;

  if(is.null(qvec)){
    return();
  }
  # style for highlighted background for unmatched names
  pre.style<-NULL;
  post.style<-NULL;

  # style for no matches
  if(mSetObj$dataSet$q.type == "name"){
    no.prestyle<-"<strong style=\"background-color:yellow; font-size=125%; color=\"black\">";
    no.poststyle<-"</strong>";
  }else{
    no.prestyle<-"<strong style=\"background-color:red; font-size=125%; color=\"black\">";
    no.poststyle<-"</strong>";
  }

  hit.inx<-mSetObj$name.map$hit.inx;
  hit.values<-mSetObj$name.map$hit.values;
  match.state<-mSetObj$name.map$match.state;

  # construct the result table with cells wrapped in html tags
  # the unmatched will be highlighted in different background
  html.res <- matrix("", nrow=length(qvec), ncol=8);
  csv.res <- matrix("", nrow=length(qvec), ncol=9);
  colnames(csv.res) <- c("Query", "Match", "HMDB", "PubChem", "ChEBI", "KEGG", "METLIN", "SMILES", "Comment");

  if(model){

    cmpd.db <- .get.my.lib("model_metabolites_master.qs");

    for (i in 1:length(qvec)){
      if(match.state[i]==1){
        pre.style<-"";
        post.style="";
      }else{ # no matches
        pre.style<-no.prestyle;
        post.style<-no.poststyle;
      }
      hit <-cmpd.db[hit.inx[i], ,drop=F];
      html.res[i, ]<-c(paste(pre.style, qvec[i], post.style, sep=""),
                       paste(ifelse(match.state[i]==0, "", hit.values[i]), sep=""),
                       paste(ifelse(match.state[i]==0 || is.na(hit$HMDB) || hit$HMDB=="" || hit$HMDB=="NA","-", paste("<a href=http://www.hmdb.ca/metabolites/", hit$HMDB, " target='_blank'>",hit$HMDB,"</a>", sep="")),  sep=""),
                       paste(ifelse(match.state[i]==0 || is.na(hit$PubChem) || hit$PubChem=="" || hit$PubChem=="NA", "-", paste("<a href=http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=", hit$PubChem," target='_blank'>", hit$PubChem,"</a>", sep="")), sep=""),
                       paste(ifelse(match.state[i]==0 || is.na(hit$ChEBI) || hit$ChEBI==""|| hit$ChEBI=="NA","-", paste("<a href=http://www.ebi.ac.uk/chebi/searchId.do?chebiId=", hit$ChEBI, " target='_blank'>",hit$ChEBI,"</a>", sep="")), sep=""),
                       paste(ifelse(match.state[i]==0 || is.na(hit$KEGG) || hit$KEGG==""|| hit$KEGG=="NA","-",paste("<a href=http://www.genome.jp/dbget-bin/www_bget?", hit$KEGG, " target='_blank'>", hit$KEGG,"</a>", sep="")), sep=""),
                       paste(ifelse(match.state[i]==0 || is.na(hit$Metlin) || hit$Metlin==""|| hit$Metlin=="NA","-",paste("<a href=http://metlin.scripps.edu/metabo_info.php?molid=", hit$Metlin," target='_blank'>",hit$Metlin,"</a>", sep="")), sep=""),
                       ifelse(match.state[i]!=1,"View",""));
      csv.res[i, ]<-c(qvec[i],
                      ifelse(match.state[i]==0, "NA", hit.values[i]),
                      ifelse(match.state[i]==0, "NA", hit$HMDB),
                      ifelse(match.state[i]==0, "NA", hit$PubChem),
                      ifelse(match.state[i]==0, "NA", hit$ChEBI),
                      ifelse(match.state[i]==0, "NA", hit$KEGG),
                      ifelse(match.state[i]==0, "NA", hit$Metlin),
                      ifelse(match.state[i]==0, "NA", hit$Smile),
                      match.state[i]);
    }
  }else{

    # cmpd.db <- .get.my.lib("master_cpds_microbiomenet.qs");
    cmpd.db <- .get.my.lib("master_compound_db.qs");

    for (i in 1:length(qvec)){
      if(match.state[i]==1){
        pre.style<-"";
        post.style="";
      }else{ # no matches
        pre.style<-no.prestyle;
        post.style<-no.poststyle;
      }
      hit <-cmpd.db[hit.inx[i], ,drop=F];
      html.res[i, ]<-c(paste(pre.style, qvec[i], post.style, sep=""),
                       paste(ifelse(match.state[i]==0, "", hit.values[i]), sep=""),
                       paste(ifelse(match.state[i]==0 || is.na(hit$hmdb_id) || hit$hmdb_id=="" || hit$hmdb_id=="NA","-", paste("<a href=http://www.hmdb.ca/metabolites/", hit$hmdb_id, " target='_blank'>",hit$hmdb_id,"</a>", sep="")),  sep=""),
                       paste(ifelse(match.state[i]==0 || is.na(hit$pubchem_id) || hit$pubchem_id=="" || hit$pubchem_id=="NA", "-", paste("<a href=http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=", hit$pubchem_id," target='_blank'>", hit$pubchem_id,"</a>", sep="")), sep=""),
                       paste(ifelse(match.state[i]==0 || is.na(hit$chebi_id) || hit$chebi_id==""|| hit$chebi_id=="NA","-", paste("<a href=http://www.ebi.ac.uk/chebi/searchId.do?chebiId=", hit$chebi_id, " target='_blank'>",hit$chebi_id,"</a>", sep="")), sep=""),
                       paste(ifelse(match.state[i]==0 || is.na(hit$kegg_id) || hit$kegg_id==""|| hit$kegg_id=="NA","-",paste("<a href=http://www.genome.jp/dbget-bin/www_bget?", hit$kegg_id, " target='_blank'>", hit$kegg_id,"</a>", sep="")), sep=""),
                       paste(ifelse(match.state[i]==0 || is.na(hit$metlin_id) || hit$metlin_id==""|| hit$metlin_id=="NA","-",paste("<a href=http://metlin.scripps.edu/metabo_info.php?molid=", hit$metlin_id," target='_blank'>",hit$metlin_id,"</a>", sep="")), sep=""),
                       ifelse(match.state[i]!=1,"View",""));
      csv.res[i, ]<-c(qvec[i],
                      ifelse(match.state[i]==0, "NA", hit.values[i]),
                      ifelse(match.state[i]==0, "NA", hit$hmdb_id),
                      ifelse(match.state[i]==0, "NA", hit$pubchem_id),
                      ifelse(match.state[i]==0, "NA", hit$chebi_id),
                      ifelse(match.state[i]==0, "NA", hit$kegg_id),
                      ifelse(match.state[i]==0, "NA", hit$metlin_id),
                      ifelse(match.state[i]==0, "NA", hit$smiles),
                      match.state[i]);
    }
  }

  # return only columns user selected
  # add query and match columns at the the beginning, and 'Detail' at the end
  return.cols <- c(TRUE, TRUE, mSetObj$return.cols, TRUE);
  html.res <- html.res[,return.cols, drop=F];
  csv.res <- csv.res[,return.cols, drop=F];

  # store the value for report
  mSetObj$dataSet$map.table <- csv.res;
  fast.write.csv(csv.res, file="name_map.csv", row.names=F);

  if(.on.public.web){
    .set.mSet(mSetObj);
    return(as.vector(html.res));
  }else{
    return(.set.mSet(mSetObj));
  }
}

#'Set metabolome filter
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param TorF Input metabolome filter
#'@export
SetMetabolomeFilter<-function(mSetObj=NA, TorF){

  mSetObj <- .get.mSet(mSetObj);

  if(!.on.public.web){
    mSetObj$api$filter <- TorF
  }

  mSetObj$dataSet$use.metabo.filter <- TorF;
  return(.set.mSet(mSetObj));
}

#'Set current user selected metset library for search
#'@description if enrichment analysis, also prepare lib by
#'creating a list of metabolite sets
#'@usage SetCurrentMsetLib(mSetObj=NA, lib.type, excludeNum)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param lib.type Input user selected name of library, "self", "kegg_pathway",
#'"smpdb_pathway", "blood", "urine", "csf", "snp", "predicted", "location", and "drug".
#'@param excludeNum Users input the mimimum number compounds within selected metabolite sets (metabolitesets < excludeNum)
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@import qs
#'@export

SetCurrentMsetLib <- function(libname, excludeNum=0){

  destfile <- paste(libname, ".qs", sep = "");

  if(!exists("current.msetlib")) {

    if(.on.public.web){

      if(grepl("kegg", libname)){
        my.qs  <- paste("../../lib/kegg_libs/", destfile, sep="");
      }else{
        my.qs  <- paste("../../lib/metacyc_libs/", destfile, sep="");
      }

      current.msetlib <- qs::qread(my.qs);

    }else{

      if(grepl("kegg", libname)){
        my.qs <- paste("https://www.microbiomenet.com/resources/lib/kegg_libs/", destfile, sep="");
      }else{
        my.qs <- paste("https://www.microbiomenet.com/resources/lib/metacyc_libs/", destfile, sep="");
      }

      if(!file.exists(destfile)){
        download.file(my.qs, destfile);
      }
      current.msetlib <- qs::qread(destfile);
    }
  }else{
    current.msetlib <- qs::qread(destfile);
  }

  # create a named list, use the ids for list names
  ms.list <- strsplit(as.character(current.msetlib[,3]), "; ", fixed=TRUE);
  names(ms.list) <- current.msetlib[,2];

  if(excludeNum > 0){
    cmpd.count <- lapply(ms.list, length);
    sel.inx <- cmpd.count >= excludeNum;
    ms.list <- ms.list[sel.inx];
    current.msetlib <- current.msetlib[sel.inx, ]
  }

  # update current.mset and push to global env
  current.msetlib$member <- ms.list;

  current.msetlib <<- current.msetlib;
}

#'Return the final (after user selection) map as dataframe
#'@description Returns three columns: original name, HMDB name and KEGG ID,
#'for enrichment and pathway analysis, respectively
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
GetFinalNameMap <- function(mSetObj=NA, lipid = FALSE, model = T){

  mSetObj <- .get.mSet(mSetObj);

  lipid = mSetObj$lipid.feats

  if (is.null(lipid)) {
    lipid = FALSE
  }

  hit.inx <- mSetObj$name.map$hit.inx;
  hit.values <- mSetObj$name.map$hit.values;
  match.state <- mSetObj$name.map$match.state;

  qvec <- mSetObj$dataSet$cmpd;
  nm.mat <- matrix(nrow=length(qvec), ncol=5);
  colnames(nm.mat) <- c("query", "hmdb",  "kegg", "pubchem", "name");

  if(model){

    cmpd.db <- .get.my.lib("model_metabolites_master.qs");

    for (i in 1:length(qvec)){
      hit <-cmpd.db[hit.inx[i], ,drop=F];
      if(match.state[i]==0){
        hmdb.hit <- NA;
        pubchem.hit.id <- NA;
        kegg.hit <- NA;
        nm.hit <- NA
      }else{
        hmdb.hit <- ifelse(nchar(hit.values[i])==0, NA, hit.values[i]);
        pubchem.hit.id <- ifelse(nchar(hit$pubchem)==0, NA, hit$pubchem);
        kegg.hit <- ifelse(nchar(hit$kegg)==0, NA, hit$kegg);
        nm.hit <- ifelse(nchar(hit$name)==0, NA, hit$name);
      }
      nm.mat[i, ]<-c(qvec[i], hmdb.hit, kegg.hit, pubchem.hit.id, nm.hit);
    }
  }else{
    cmpd.db <- .get.my.lib("master_compound_db.qs");
    # cmpd.db <- .get.my.lib("master_cpds_microbiomenet.qs");

    for (i in 1:length(qvec)){
      hit <- cmpd.db[hit.inx[i], ,drop=F];
      if(match.state[i]==0){
        hmdb.hit <- NA;
        pubchem.hit.id <- NA;
        kegg.hit <- NA;
        nm.hit <- NA
      }else{
        hmdb.hit <- ifelse(nchar(hit.values[i])==0, NA, hit.values[i]);
        pubchem.hit.id <- ifelse(nchar(hit$pubchem_id)==0, NA, hit$pubchem_id);
        kegg.hit <- ifelse(nchar(hit$kegg_id)==0, NA, hit$kegg_id);
        nm.hit <- ifelse(nchar(hit$name)==0, NA, hit$name);
      }
      nm.mat[i, ]<-c(qvec[i], hmdb.hit, kegg.hit, pubchem.hit.id, nm.hit);
    }
  }
  return(as.data.frame(nm.mat));
}

#'Over-representation analysis using hypergeometric tests
#'@description Over-representation analysis using hypergeometric tests
#'The probability is calculated from obtaining equal or higher number
#'of hits using 1-phyper. Since phyper is a cumulative probability,
#'to get P(X>=hit.num) => P(X>(hit.num-1))
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
CalculatePermHyperScore <- function(mSetObj=NA, model = T, perm = F, library.name){

  mSetObj <- .get.mSet(mSetObj);

  # make a clean dataSet$cmpd data based on name mapping
  # only valid hmdb name will be used
  nm.map <- GetFinalNameMap(mSetObj, model = model);

  # KEGG libraries formatted with KEGG compound IDS
  if(grepl("kegg", library.name)){
    valid.inx <- !(is.na(nm.map$kegg)| duplicated(nm.map$kegg));
    ora.vec <- nm.map$kegg[valid.inx];
  }else if(grepl("mset", library.name)){
    valid.inx <- !(is.na(nm.map$name)| duplicated(nm.map$name));
    ora.vec <- nm.map$name[valid.inx];
  }else{
    # MetaCyc + RefMet libraries formatted with PubChem IDs
    valid.inx <- !(is.na(nm.map$pubchem)| duplicated(nm.map$pubchem));
    ora.vec <- as.character(nm.map$pubchem[valid.inx]);
  }

  q.size<-length(ora.vec);

  if(is.na(ora.vec) || q.size==0) {
    AddErrMsg("No valid compound names found!");
    return(0);
  }

  current.mset <- current.msetlib$member

  # make a clean metabilite set based on reference metabolome filtering
  # also need to update ora.vec to the updated mset
  if(mSetObj$dataSet$use.metabo.filter && !is.null(mSetObj$dataSet$metabo.filter.hmdb)){
    current.mset <- lapply(current.mset, function(x){x[x %in% mSetObj$dataSet$metabo.filter.hmdb]})
    mSetObj$dataSet$filtered.mset <- current.mset;
    ora.vec <- ora.vec[ora.vec %in% unique(unlist(current.mset, use.names = FALSE))]
    q.size <- length(ora.vec);
  }

  # total uniq cmpds in the current mset lib
  uniq.count <- length(unique(unlist(current.mset, use.names = FALSE)));

  set.size<-length(current.mset);

  if(set.size ==1){
    AddErrMsg("Cannot perform enrichment analysis on a single metabolite set!");
    return(0);
  }

  hits<-lapply(current.mset, function(x){x[x %in% ora.vec]});
  hit.num<-unlist(lapply(hits, function(x) length(x)), use.names = FALSE);

  if(sum(hit.num>0)==0){
    AddErrMsg("No match was found to the selected metabolite set library!");
    return(0);
  }

  set.num <- unlist(lapply(current.mset, length), use.names = FALSE);

  # prepare for the result table
  res.mat <- matrix(NA, nrow=set.size, ncol=1);
  rownames(res.mat) <- names(current.mset);

  for(i in 1:set.size){
    res.mat[i,1]<-phyper(hit.num[i]-1, set.num[i], uniq.count-set.num[i], q.size, lower.tail=F);
  }

  mSetObj$analSet$ora.mat <- res.mat

  return(.set.mSet(mSetObj));
}

#'Over-representation analysis using hypergeometric tests
#'@description Over-representation analysis using hypergeometric tests
#'The probability is calculated from obtaining equal or higher number
#'of hits using 1-phyper. Since phyper is a cumulative probability,
#'to get P(X>=hit.num) => P(X>(hit.num-1))
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
CalculateTrueHyperScore <- function(mSetObj=NA, model = T, perm = F,
                                    library.name, record){

  mSetObj <- .get.mSet(mSetObj);

  # make a clean dataSet$cmpd data based on name mapping
  # only valid hmdb name will be used
  nm.map <- GetFinalNameMap(mSetObj, model = model);

  # KEGG libraries formatted with KEGG compound IDS
  if(grepl("kegg", library.name)){
    valid.inx <- !(is.na(nm.map$kegg)| duplicated(nm.map$kegg));
    ora.vec <- nm.map$kegg[valid.inx];
  }else if(grepl("mset", library.name)){
    valid.inx <- !(is.na(nm.map$name)| duplicated(nm.map$name));
    ora.vec <- nm.map$name[valid.inx];
  }else{
    # MetaCyc + RefMet libraries formatted with PubChem IDs
    valid.inx <- !(is.na(nm.map$pubchem)| duplicated(nm.map$pubchem));
    ora.vec <- as.character(nm.map$pubchem[valid.inx]);
  }

  q.size<-length(ora.vec);

  if(is.na(ora.vec) || q.size==0) {
    AddErrMsg("No valid compound names found!");
    return(0);
  }

  current.mset <- current.msetlib$member

  # make a clean metabilite set based on reference metabolome filtering
  # also need to update ora.vec to the updated mset
  if(mSetObj$dataSet$use.metabo.filter && !is.null(mSetObj$dataSet$metabo.filter.hmdb)){
    current.mset <- lapply(current.mset, function(x){x[x %in% mSetObj$dataSet$metabo.filter.hmdb]})
    mSetObj$dataSet$filtered.mset <- current.mset;
    ora.vec <- ora.vec[ora.vec %in% unique(unlist(current.mset, use.names = FALSE))]
    q.size <- length(ora.vec);
  }

  # total uniq cmpds in the current mset lib
  uniq.count <- length(unique(unlist(current.mset, use.names = FALSE)));

  set.size<-length(current.mset);

  if(set.size ==1){
    AddErrMsg("Cannot perform enrichment analysis on a single metabolite set!");
    return(0);
  }

  hits <- lapply(current.mset, function(x){x[x %in% ora.vec]});
  hit.num <- unlist(lapply(hits, function(x) length(x)), use.names = FALSE);
  feat_vec <- sapply(hits, function(x) paste(x, collapse=";"))

  if(sum(hit.num>0)==0){
    AddErrMsg("No match was found to the selected metabolite set library!");
    return(0);
  }

  set.num<-unlist(lapply(current.mset, length), use.names = FALSE);

  # prepare for the result table
  res.mat<-matrix(NA, nrow=set.size, ncol=6);
  rownames(res.mat)<-names(current.mset);
  colnames(res.mat)<-c("total", "expected", "hits", "Raw p", "Holm p", "FDR");

  for(i in 1:set.size){
    res.mat[i,1] <- set.num[i];
    res.mat[i,2] <- q.size*(set.num[i]/uniq.count);
    res.mat[i,3] <- hit.num[i];
    res.mat[i,4] <- phyper(hit.num[i]-1, set.num[i], uniq.count-set.num[i], q.size, lower.tail=F);
    # negneg[[i]] <- uniq.count + hit.num[[i]] - set.num[i] - q.size;
  }
#
#   first <- unlist(lapply(hit.num, function(x) max(0, x-1)));
#   easematrix <- cbind(first, (set.num - hit.num + 1), (q.size - hit.num), unlist(negneg));

  # res.mat[,7] <- apply(easematrix, 1, function(x) fisher.test(matrix(x, nrow=2), alternative = "greater")$p.value);

  ##################################### from ms peaks

  # Gamma-adjusted p-values
#
#   perm_record <- qs::qread("permutation_hits.qs")
#   perm_record <- unlist(perm_record)
#
#   perm_minus <- abs(0.9999999999 - perm_record);
#   sigpvalue <-res.mat[,7]
#
#   sig_hits <- res.mat[,3]
#
#   if(length(sig_hits[sig_hits!=0]) < round(length(sig_hits)*0.05)){ # too few hits that can't calculate gamma dist!
#     if(!exists("adjustedp")){
#       adjustedp <- rep(NA, length = length(res.mat[,1]))
#     }
#     res.mat <- cbind(res.mat, Gamma=adjustedp);
#   }else{

    # tryCatch({
    #   fit.gamma <- fitdistrplus::fitdist(perm_minus, distr = "gamma", method = "mle", lower = c(0, 0), start = list(scale = 1, shape = 1));
    #   rawpval <- as.numeric(sigpvalue);
    #   adjustedp <- 1 - (pgamma(1-rawpval, shape = fit.gamma$estimate["shape"], rate = fit.gamma$estimate["scale"]));
    # }, error = function(e){
    #   if(!exists("adjustedp")){
    #     adjustedp <- rep(NA, length = length(res.mat[,1]))
    #   }
    #   res.mat <- cbind(res.mat, Gamma=adjustedp);
    #   print(e)
    # }, finally = {
    #   if(!exists("adjustedp")){
    #     adjustedp <- rep(NA, length = length(res.mat[,1]))
    #   }
    #   res.mat <- cbind(res.mat, Gamma=adjustedp);
    # })

#   }

  # calculate empirical p-values
  record <- record
  fisher.p <- as.numeric(res.mat[,4])

  #pathway in rows, perms in columns
  record_matrix <- do.call(cbind, do.call(cbind, record))
  num_perm <- ncol(record_matrix)

  #number of better hits for web
  better.hits <- sapply(seq_along(record_matrix[,1]), function(i) sum(record_matrix[i,] <= fisher.p[i])  )

  #account for a bias due to finite sampling - Davison and Hinkley (1997)
  emp.p <- sapply(seq_along(record_matrix[,1]), function(i) (sum(record_matrix[i,] <= fisher.p[i])/num_perm) )

  res.mat <- cbind(res.mat, Emp.Hits=better.hits, Empirical=emp.p, Cpd.Hits = feat_vec)

  ############### end of from ms peaks

  res.mat <- res.mat[hit.num>1, , drop=FALSE];

  # adjust for multiple testing problems
  res.mat[,5] <- p.adjust(res.mat[,4], "holm");
  res.mat[,6] <- p.adjust(res.mat[,4], "fdr");

  ord.inx <- order(as.numeric(res.mat[,4]));
  res.mat <- res.mat[ord.inx,]
  res.mat <- res.mat[!duplicated(rownames(res.mat)), ]

  mSetObj$analSet$ora.mat <- res.mat
  return(.set.mSet(mSetObj));
}

#'Over-representation analysis using hypergeometric tests
#'@description Over-representation analysis using hypergeometric tests
#'The probability is calculated from obtaining equal or higher number
#'of hits using 1-phyper. Since phyper is a cumulative probability,
#'to get P(X>=hit.num) => P(X>(hit.num-1))
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
CalculateHyperScore <- function(mSetObj=NA, model = T){

  mSetObj <- .get.mSet(mSetObj);

  # make a clean dataSet$cmpd data based on name mapping
  # only valid hmdb name will be used
  nm.map <- GetFinalNameMap(mSetObj, model = model);

  # KEGG libraries formatted with KEGG compound IDS
  if(grepl("kegg", mSetObj$analSet$msetlibname)){
    valid.inx <- !(is.na(nm.map$kegg)| duplicated(nm.map$kegg));
    ora.vec <- nm.map$kegg[valid.inx];
  }else if(grepl("mset", mSetObj$analSet$msetlibname)){
    valid.inx <- !(is.na(nm.map$name)| duplicated(nm.map$name));
    ora.vec <- nm.map$name[valid.inx];
  }else{
    # MetaCyc + RefMet libraries formatted with PubChem IDs
    valid.inx <- !(is.na(nm.map$pubchem)| duplicated(nm.map$pubchem));
    ora.vec <- as.character(nm.map$pubchem[valid.inx]);
  }

  q.size<-length(ora.vec);

  if(is.na(ora.vec) || q.size==0) {
    AddErrMsg("No valid compound names found!");
    return(0);
  }

  current.mset <- current.msetlib$member

  # make a clean metabilite set based on reference metabolome filtering
  # also need to update ora.vec to the updated mset
  if(mSetObj$dataSet$use.metabo.filter && !is.null(mSetObj$dataSet$metabo.filter.hmdb)){
    current.mset <- lapply(current.mset, function(x){x[x %in% mSetObj$dataSet$metabo.filter.hmdb]})
    mSetObj$dataSet$filtered.mset <- current.mset;
    ora.vec <- ora.vec[ora.vec %in% unique(unlist(current.mset, use.names = FALSE))]
    q.size <- length(ora.vec);
  }

  # total uniq cmpds in the current mset lib
  uniq.count <- length(unique(unlist(current.mset, use.names = FALSE)));

  set.size<-length(current.mset);

  if(set.size ==1){
    AddErrMsg("Cannot perform enrichment analysis on a single metabolite set!");
    return(0);
  }

  hits<-lapply(current.mset, function(x){x[x %in% ora.vec]});
  hit.num<-unlist(lapply(hits, function(x) length(x)), use.names = FALSE);

  if(sum(hit.num>0)==0){
    AddErrMsg("No match was found to the selected metabolite set library!");
    return(0);
  }

  set.num<-unlist(lapply(current.mset, length), use.names = FALSE);

  # prepare for the result table
  res.mat<-matrix(NA, nrow=set.size, ncol=6);
  rownames(res.mat)<-names(current.mset);
  colnames(res.mat)<-c("total", "expected", "hits", "Raw p", "Holm p", "FDR");

  for(i in 1:set.size){
    res.mat[i,1]<-set.num[i];
    res.mat[i,2]<-q.size*(set.num[i]/uniq.count);
    res.mat[i,3]<-hit.num[i];
    res.mat[i,4]<-phyper(hit.num[i]-1, set.num[i], uniq.count-set.num[i], q.size, lower.tail=F);
  }

  if(nrow(res.mat > 100)){

    if(sum(res.mat[,3] == 0)/nrow(res.mat) > 0.5 ){
      res.mat <- res.mat[hit.num>0,,drop=FALSE];
    }else{
      res.mat <- res.mat[hit.num>1,,drop=FALSE];
    }

  }else{
    res.mat <- res.mat[hit.num>0,,drop=FALSE];
  }

  # adjust for multiple testing problems
  res.mat[,5] <- p.adjust(res.mat[,4], "holm");
  res.mat[,6] <- p.adjust(res.mat[,4], "fdr");

  ord.inx<-order(res.mat[,4]);
  mSetObj$analSet$ora.mat <- signif(res.mat[ord.inx,],3);
  mSetObj$analSet$ora.hits <- hits;

  fast.write.csv(mSetObj$analSet$ora.mat, file="msea_ora_result.csv");
  return(.set.mSet(mSetObj));
}

##########################################################################################################
################# Plotting Functions
##########################################################################################################

#'Plot over-representation analysis (ORA)
#'@description Plot over-representation analysis (ORA)
#'@usage PlotORA(mSetObj=NA, imgName, imgOpt, format="png", dpi=72, width=NA)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param imgName Input a name for the plot
#'@param imgOpt "net"
#'@param format Select the image format, "png", or "pdf".
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images,
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.
#'@param width Input the width, there are 2 default widths, the first, width = NULL, is 10.5.
#'The second default is width = 0, where the width is 7.2. Otherwise users can input their own width.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PlotORA <- function(mbSetObj=NA, imgName, imgOpt,
                  format="png", dpi=72, width=NA){

  mbSetObj <- .get.mSet(mbSetObj);
  gemLib <- mbSetObj$mNet$gemLib

  if(gemLib == "both"){
    res <- c("msea_ora_result_a.csv", "msea_ora_result_c.csv")
    imgName <- c(paste(imgName, "dpi", dpi, "agora", ".", format, sep=""), paste(imgName, "dpi", dpi, "carveme", ".", format, sep=""))
  }else if(gemLib == "agora"){
    res <- "msea_ora_result_a.csv"
    imgName = paste(imgName, "dpi", dpi, "agora", ".", format, sep="");
  }else{
    res <- "msea_ora_result_c.csv"
    imgName = paste(imgName, "dpi", dpi, "carveme", ".", format, sep="");
  }

  mbSetObj$analSet$enrich.net <- list()

  for(i in seq_along(res)){

    ora.mat <- read.csv(res[[i]], row.names = 1)

    #calculate the enrichment fold change
    folds <- ora.mat[,3]/ora.mat[,2];
    names(folds) <- GetShortNames(rownames(ora.mat));
    pvals <- ora.mat[,4];

    if(is.na(width)){
      w <- 9;
    }else if(width == 0){
      w <- 7;
    }else{
      w <-width;
    }
    h <- w;

    Cairo::Cairo(file = imgName[[i]], unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");

    PlotMSEA.Overview(folds, pvals);
    dev.off();

    if(.on.public.web){
      mbSetObj$analSet$enrich.net[[i]] <- PlotEnrichNet.Overview(folds, pvals);
    }
  }
  return(.set.mbSetObj(mbSetObj))
}

#'Plot MSEA overview
#'@description Barplot height is enrichment fold change
#'color is based on p values, used in higher functions
#'@usage PlotMSEA.Overview(folds, pvals)
#'@param folds Input the fold-change values
#'@param pvals Input the p-values
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
PlotMSEA.Overview <- function(folds, pvals){

  # due to space limitation, plot top 50 if more than 50 were given
  title <- "Metabolite Sets Enrichment Overview";
  ht.col <- GetMyHeatCols(length(folds));
  if(length(folds) > 25){
    folds <- folds[1:25];
    pvals <- pvals[1:25];
    ht.col <- ht.col[1:25];
    title <- "Enrichment Overview (top 25)";
  }

  op <- par(mar=c(5,20,4,6), oma=c(0,0,0,4));

  barplot(rev(folds), horiz=T, col=rev(ht.col),
          xlab="Enrichment Ratio", las=1, cex.name=0.75, space=c(0.5, 0.5),
          main= title);

  minP <- min(pvals);
  maxP <- max(pvals);
  medP <- (minP+maxP)/2;

  axs.args <- list(at=c(minP, medP, maxP), labels=format(c(maxP, medP, minP), scientific=T, digit=1), tick = F);
  image.plot(legend.only=TRUE, zlim=c(minP, maxP), col=rev(ht.col),
             axis.args=axs.args, legend.shrink = 0.4, legend.lab="P value");
  par(op);
}

#'Plot MSEA Dot Plot
#'@description Dot plot of enrichment analysis results.
#'@usage PlotEnrichDotPlot(mSet)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param enrichType Input whether the enrichment analysis was over-respresentation
#'analysis (ora) or quantitative enrichment analysis (qea).
#'@param imgName Input a name for the plot
#'@param format Select the image format, "png", or "pdf".
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images,
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.
#'@param width Input the width, there are 2 default widths, the first, width = NULL, is 10.5.
#'The second default is width = 0, where the width is 7.2. Otherwise users can input their own width.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PlotEnrichDotPlot <- function(mbSetObj=NA, imgName,
                              format="png", dpi=72, width=NA){

  mbSetObj <- .get.mSet(mbSetObj);

  gemLib <- mbSetObj$mNet$gemLib

  if(.on.public.web){
    load_ggplot()
  }

  if(gemLib == "both"){
    res <- c("msea_ora_result_a.csv", "msea_ora_result_c.csv")
    imgName = c(paste(imgName, "dpi", dpi, "agora", ".", format, sep=""), paste(imgName, "dpi", dpi, "carveme", ".", format, sep=""))
  }else if(gemLib == "agora"){
    res <- "msea_ora_result_a.csv"
    imgName = paste(imgName, "dpi", dpi, "agora", ".", format, sep="");
  }else{
    res <- "msea_ora_result_c.csv"
    imgName = paste(imgName, "dpi", dpi, "carveme", ".", format, sep="");
  }

  for(i in seq_along(res)){
    results <- read.csv(res[[i]], row.names = 1)
    my.cols <- GetMyHeatCols(nrow(results));

    if(nrow(results) > 25){
      results <- results[1:25,]
      my.cols <- my.cols[1:25];
    }

    df <- data.frame(Name = factor(row.names(results), levels = rev(row.names(results))),
                     rawp = results[,4],
                     logp = -log10(results[,4]),
                     folds = results[,3]/results[,2])

    maxp <- max(df$rawp);

    if(is.na(width)){
      w <- 12;
      h <- 9;
    }else if(width == 0){
      h <- w <- 7;
    }else{
      h <- w <- width;
    }

    p <- ggplot(df,
                aes(x = logp, y = Name)) +
      geom_point(aes(size = folds, color = rawp)) + scale_size_continuous(range = c(2, 8)) +
      theme_bw(base_size = 14.5) +
      scale_colour_gradient(limits=c(0, maxp), low=my.cols[1], high = my.cols[length(my.cols)]) +
      ylab(NULL) + xlab("-log10 (p-value)") +
      ggtitle("Overview of Enriched Pathways") +
      theme(legend.text=element_text(size=14),
            legend.title=element_text(size=15))

    p$labels$colour <- "P-value"
    p$labels$size <- "Enrichment Ratio"

    ggsave(p, filename = imgName[[i]], dpi=dpi, width=w, height=h)

  }
  return(.set.mbSetObj(mbSetObj))
}

#'Barplot height is enrichment fold change
#'@description Used in higher functions, the color is based on p values
#'@param folds Input fold-change for bar plot
#'@param pvals Input p-values for bar plot
#'@param layoutOpt Input the layout option, default is set to layout.fruchterman.reingold
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import igraph
#'@import reshape

PlotEnrichNet.Overview <- function(folds, pvals, layoutOpt=layout.fruchterman.reingold){

  # due to space limitation, plot top 50 if more than 50 were given
  title <- "Enrichment Network Overview";

  if(length(folds) > 50){
    folds <- folds[1:50];
    pvals <- pvals[1:50];
    title <- "Enrichment Overview (top 50)";
  }

  if(.on.public.web){
    load_igraph()
    load_reshape()
  }

  pvalue <- pvals;
  id <- names(pvalue);
  geneSets <- current.msetlib$member;
  n <- length(pvalue);
  w <- matrix(NA, nrow=n, ncol=n);
  colnames(w) <- rownames(w) <- id;

  for (i in 1:n) {
    for (j in i:n) {
      w[i,j] = overlap_ratio(geneSets[id[i]], geneSets[id[j]])
    }
  }

  wd <- melt(w);
  wd <- wd[wd[,1] != wd[,2],];
  wd <- wd[!is.na(wd[,3]),];
  g <- graph.data.frame(wd[,-3], directed=F);
  E(g)$width <- sqrt(wd[,3]*20);
  g <- delete.edges(g, E(g)[wd[,3] < 0.25]);
  V(g)$color <- heat.colors(length(pvalue));

  cnt <- folds;
  names(cnt) <- id;
  V(g)$size <- cnt + 3;

  pos.xy <- layout.fruchterman.reingold(g,niter=500);

  # now create the json object
  nodes <- vector(mode="list");
  node.nms <- V(g)$name;
  node.sizes <- V(g)$size;
  node.cols <- V(g)$color;
  for(i in 1:length(node.sizes)){
    nodes[[i]] <- list(
      id = node.nms[i],
      label=node.nms[i],
      size=node.sizes[i],
      color=node.cols[i],
      x = pos.xy[i,1],
      y = pos.xy[i,2]
    );
  }

  edge.mat <- get.edgelist(g);
  edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2]);

  # covert to json
  netData <- list(nodes=nodes, edges=edge.mat);
  sink("msea_network.json");
  cat(RJSONIO::toJSON(netData));
  sink();

  return(g);
}

##################################################
############## image utilities ###################

# Plot a strip of color key beside a figure
# Adapted from the image.plot in fields package to correct label
# so that the small p value is bigger, located on top of the color key
# Jeff Xia, jeff.xia@mcgill.ca
# McGill University, Canada
# License: GNU GPL (>= 2)

image.plot <- function(..., add = FALSE, nlevel = 64,
                       horizontal = FALSE, legend.shrink = 0.9, legend.width = 1.2,
                       legend.mar = ifelse(horizontal, 3.1, 5.1), legend.lab = NULL,
                       graphics.reset = FALSE, bigplot = NULL, smallplot = NULL,
                       legend.only = FALSE, col = tim.colors(nlevel), lab.breaks = NULL,
                       axis.args = NULL, legend.args = NULL, midpoint = FALSE) {

  old.par <- par(no.readonly = TRUE)
  #  figure out zlim from passed arguments
  info <- image.plot.info(...)
  if (add) {
    big.plot <- old.par$plt
  }
  if (legend.only) {
    graphics.reset <- TRUE
  }
  if (is.null(legend.mar)) {
    legend.mar <- ifelse(horizontal, 3.1, 5.1)
  }
  #
  # figure out how to divide up the plotting real estate.
  #
  temp <- image.plot.plt(add = add, legend.shrink = legend.shrink,
                         legend.width = legend.width, legend.mar = legend.mar,
                         horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
  #
  # bigplot are plotting region coordinates for image
  # smallplot are plotting coordinates for legend
  smallplot <- temp$smallplot
  bigplot <- temp$bigplot
  #
  # draw the image in bigplot, just call the R base function
  # or poly.image for polygonal cells note logical switch
  # for poly.grid parsed out of call from image.plot.info
  if (!legend.only) {
    if (!add) {
      par(plt = bigplot)
    }
    if (!info$poly.grid) {
      image(..., add = add, col = col)
    }
    else {
      poly.image(..., add = add, col = col, midpoint = midpoint)
    }
    big.par <- par(no.readonly = TRUE)
  }
  ##
  ## check dimensions of smallplot
  if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
    par(old.par)
    stop("plot region too small to add legend\n")
  }
  # Following code draws the legend using the image function
  # and a one column image.
  # calculate locations for colors on legend strip
  ix <- 1
  minz <- info$zlim[1]
  maxz <- info$zlim[2]
  binwidth <- (maxz - minz)/nlevel
  midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
  iy <- midpoints
  iz <- matrix(iy, nrow = 1, ncol = length(iy))
  # extract the breaks from the ... arguments
  # note the breaks delineate intervals of common color
  breaks <- list(...)$breaks
  # draw either horizontal or vertical legends.
  # using either suggested breaks or not -- a total of four cases.
  #
  # next par call sets up a new plotting region just for the legend strip
  # at the smallplot coordinates
  par(new = TRUE, pty = "m", plt = smallplot, err = -1)
  # create the argument list to draw the axis
  #  this avoids 4 separate calls to axis and allows passing extra
  # arguments.
  # then add axis with specified lab.breaks at specified breaks
  if (!is.null(breaks) & !is.null(lab.breaks)) {
    # axis with labels at break points
    axis.args <- c(list(side = ifelse(horizontal, 1, 4),
                        mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2),
                        at = breaks, labels = lab.breaks), axis.args)
  }
  else {
    # If lab.breaks is not specified, with or without breaks, pretty
    # tick mark locations and labels are computed internally,
    # or as specified in axis.args at the function call
    axis.args <- c(list(side = ifelse(horizontal, 1, 4),
                        mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)),
                   axis.args)
  }
  #
  # draw color scales the four cases are horizontal/vertical breaks/no breaks
  # add a label if this is passed.
  if (!horizontal) {
    if (is.null(breaks)) {
      image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "",
            ylab = "", col = col)
    }
    else {
      image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "",
            ylab = "", col = col, breaks = breaks)
    }
  }
  else {
    if (is.null(breaks)) {
      image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "",
            ylab = "", col = col)
    }
    else {
      image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "",
            ylab = "", col = col, breaks = breaks)
    }
  }

  #
  # now add the axis to the legend strip.
  # notice how all the information is in the list axis.args
  #
  do.call("axis", axis.args)

  # add a box around legend strip
  box()

  #
  # add a label to the axis if information has been  supplied
  # using the mtext function. The arguments to mtext are
  # passed as a list like the drill for axis (see above)
  #
  if (!is.null(legend.lab)) {
    legend.args <- list(text = legend.lab, side = ifelse(horizontal,
                                                         1, 3), line = 1)
  }
  #
  # add the label using mtext function
  if (!is.null(legend.args)) {
    do.call(mtext, legend.args)
  }
  #
  #
  # clean up graphics device settings
  # reset to larger plot region with right user coordinates.
  mfg.save <- par()$mfg
  if (graphics.reset | add) {
    par(old.par)
    par(mfg = mfg.save, new = FALSE)
    invisible()
  }
  else {
    par(big.par)
    par(plt = big.par$plt, xpd = FALSE)
    par(mfg = mfg.save, new = FALSE)
    invisible()
  }
}


"image.plot.info" <- function(...) {
  temp <- list(...)
  #
  xlim <- NA
  ylim <- NA
  zlim <- NA
  poly.grid <- FALSE
  #
  # go through various cases of what these can be
  #
  ##### x,y,z list is first argument
  if (is.list(temp[[1]])) {
    xlim <- range(temp[[1]]$x, na.rm = TRUE)
    ylim <- range(temp[[1]]$y, na.rm = TRUE)
    zlim <- range(temp[[1]]$z, na.rm = TRUE)
    if (is.matrix(temp[[1]]$x) & is.matrix(temp[[1]]$y) &
        is.matrix(temp[[1]]$z)) {
      poly.grid <- TRUE
    }
  }
  ##### check for polygrid first three arguments should be matrices
  #####
  if (length(temp) >= 3) {
    if (is.matrix(temp[[1]]) & is.matrix(temp[[2]]) & is.matrix(temp[[3]])) {
      poly.grid <- TRUE
    }
  }
  #####  z is passed without an  x and y  (and not a poly.grid!)
  #####
  if (is.matrix(temp[[1]]) & !poly.grid) {
    xlim <- c(0, 1)
    ylim <- c(0, 1)
    zlim <- range(temp[[1]], na.rm = TRUE)
  }
  ##### if x,y,z have all been passed find their ranges.
  ##### holds if poly.grid or not
  #####
  if (length(temp) >= 3) {
    if (is.matrix(temp[[3]])) {
      xlim <- range(temp[[1]], na.rm = TRUE)
      ylim <- range(temp[[2]], na.rm = TRUE)
      zlim <- range(temp[[3]], na.rm = TRUE)
    }
  }
  #### parse x,y,z if they are  named arguments
  # determine if  this is polygon grid (x and y are matrices)
  if (is.matrix(temp$x) & is.matrix(temp$y) & is.matrix(temp$z)) {
    poly.grid <- TRUE
  }
  xthere <- match("x", names(temp))
  ythere <- match("y", names(temp))
  zthere <- match("z", names(temp))
  if (!is.na(zthere))
    zlim <- range(temp$z, na.rm = TRUE)
  if (!is.na(xthere))
    xlim <- range(temp$x, na.rm = TRUE)
  if (!is.na(ythere))
    ylim <- range(temp$y, na.rm = TRUE)
  # overwrite zlims with passed values
  if (!is.null(temp$zlim))
    zlim <- temp$zlim
  if (!is.null(temp$xlim))
    xlim <- temp$xlim
  if (!is.null(temp$ylim))
    ylim <- temp$ylim
  list(xlim = xlim, ylim = ylim, zlim = zlim, poly.grid = poly.grid)
}

# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

image.plot.plt <- function(x, add = FALSE, legend.shrink = 0.9,
                           legend.width = 1, horizontal = FALSE, legend.mar = NULL,
                           bigplot = NULL, smallplot = NULL, ...) {
  old.par <- par(no.readonly = TRUE)
  if (is.null(smallplot))
    stick <- TRUE
  else stick <- FALSE
  if (is.null(legend.mar)) {
    legend.mar <- ifelse(horizontal, 3.1, 5.1)
  }
  # compute how big a text character is
  char.size <- ifelse(horizontal, par()$cin[2]/par()$din[2],
                      par()$cin[1]/par()$din[1])
  # This is how much space to work with based on setting the margins in the
  # high level par command to leave between strip and big plot
  offset <- char.size * ifelse(horizontal, par()$mar[1], par()$mar[4])
  # this is the width of the legned strip itself.
  legend.width <- char.size * legend.width
  # this is room for legend axis labels
  legend.mar <- legend.mar * char.size
  # smallplot is the plotting region for the legend.
  if (is.null(smallplot)) {
    smallplot <- old.par$plt
    if (horizontal) {
      smallplot[3] <- legend.mar
      smallplot[4] <- legend.width + smallplot[3]
      pr <- (smallplot[2] - smallplot[1]) * ((1 - legend.shrink)/2)
      smallplot[1] <- smallplot[1] + pr
      smallplot[2] <- smallplot[2] - pr
    }
    else {
      smallplot[2] <- 1 - legend.mar
      smallplot[1] <- smallplot[2] - legend.width
      pr <- (smallplot[4] - smallplot[3]) * ((1 - legend.shrink)/2)
      smallplot[4] <- smallplot[4] - pr
      smallplot[3] <- smallplot[3] + pr
    }
  }
  if (is.null(bigplot)) {
    bigplot <- old.par$plt
    if (!horizontal) {
      bigplot[2] <- min(bigplot[2], smallplot[1] - offset)
    }
    else {
      bottom.space <- old.par$mar[1] * char.size
      bigplot[3] <- smallplot[4] + offset
    }
  }
  if (stick & (!horizontal)) {
    dp <- smallplot[2] - smallplot[1]
    smallplot[1] <- min(bigplot[2] + offset, smallplot[1])
    smallplot[2] <- smallplot[1] + dp
  }
  return(list(smallplot = smallplot, bigplot = bigplot))
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
