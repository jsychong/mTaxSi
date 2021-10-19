#' This function reads in the user's data.
#'
#' @param mbSetObj
#' @param fileName Character, input the name of the file containing the taxonomic signature.
#' @param idType Either 'asv', 'species', 'gg_otus', 'silva'.
readMNetTable <- function(mbSetObj, fileName, idType){

  mbSetObj <- .get.mbSetObj(mbSetObj);

  msg <- NULL;
  mydata <- .readDataTable(fileName);

  print(head(mydata))

  if(any(is.na(mydata)) || class(mydata) == "try-error"){
    AddErrMsg("Failed to read in the table! Please make sure that the table
              is in the right format and does not have empty cells or NA.");
    return(0);
  }

  # look for #NAME, store in a list
  col.nms <- tolower(colnames(mydata));

  pattern <- paste(c("#name", "#regulated"), collapse = "|")

  col.nms.inx <- grep(pattern, col.nms);

  if(length(col.nms.inx) < 2){
    AddErrMsg("Please make sure you have the labels #NAME
              and #REGULATED in your table!");
    return(0);
  }

  if(nrow(mydata)==1){
    AddErrMsg("Only one feature uploaded!")
    return(0);
  }

  # empty cells or NAs cannot be tolerated
  na.inx  <- is.na(mydata) | mydata == ""

  if(sum(na.inx) > 0){
    AddErrMsg(paste("A total of", sum(na.inx), "empty or NA values were
                    found in the uploaded table!"))
    return(0)
  }

  mbSetObj$mNet <- list()
  mbSetObj$mNet$idType <- idType

  qs::qsave(mydata, "user_data_orig.qs")

  if(.on.public.web){
    .set.mbSetObj(mbSetObj)
    return(1)
  }else{
    return(.set.mbSetObj(mbSetObj))
  }
}

#' This function does a sanity check on the user's
#' uploaded table
#'
#' @param mbSetObj
#' @param userPheno Character, input the user phenotype.
#'
#' @export

PerformMicrobiomeNetSanity <- function(mbSetObj, userPheno){

  mbSetObj <- .get.mbSetObj(mbSetObj);

  if(file.exists("user_data_orig.qs")){
    taxa_tbl <- qs::qread("user_data_orig.qs")
  }else{
    AddErrMsg("User data has not been successfully uploaded!")
    return(0)
  }

  # check that the phenotype label is valid
  pheno <- userPheno
  pheno <- trimws(unlist(strsplit(pheno, "vs.")))

  if(length(pheno) <2){
    AddErrMsg("Invalid phenotype!")
    return(0)
  }

  mbSetObj$mNet$pheno <- pheno

  # check no weird symbols in pheno labels
  # check for special characters
  if(sum(is.na(iconv(pheno)))>0){
    na.inx <- is.na(iconv(pheno));
    nms <- paste(pheno[na.inx], collapse="; ");
    AddErrMsg(paste("No special letters (i.e. Latin, Greek) are allowed in sample names!", nms, collapse=" "));
    return(0);
  }

  # double check no duplicate features
  feats <- taxa_tbl[1,]

  if(length(feats) > length(unique(feats))){
    AddErrMsg("Duplicate taxonomic features!")
    return(0)
  }

  # save as qs table
  colnames(taxa_tbl) <- c("Microbe", "Regulated")

  taxa_tbl[,1] <- gsub("\t", "", taxa_tbl[,1])

  mbSetObj$mNet$feat.tbl <- taxa_tbl;
  qs::qsave(taxa_tbl, "user_data_checked.qs")

  if(.on.public.web){
    .set.mbSetObj(mbSetObj)
    return(1)
  }else{
    return(.set.mbSetObj(mbSetObj))
  }
}

#' Getters for mNet sub-object
getPheno <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj)
  return(mbSetObj$mNet$pheno)
}

getNumFeats <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj)
  return(nrow(mbSetObj$mNet$feat.tbl))
}

setGemLibrary <- function(mbSetObj, gem.libs){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$mNet$gemLib <- tolower(gem.libs)
  return(.set.mbSetObj(mbSetObj))
}

#' This function is a wrapper function
#' to match user's input to our local
#' database of GEM models.
#'
#' @param mbSetObj
#' @param identity Numeric, only useful when matching ASVs. Default is 0.97 identity.
#'
#' @description Maps user's taxonomic
#' signatures to our database of GEM models.
#' @export

match2gems <- function(mbSetObj, identity = 0.97){

  mbSetObj <- .get.mbSetObj(mbSetObj)

  taxaType <- mbSetObj$mNet$idType
  userTable <- mbSetObj$mNet$feat.tbl
  gemLib <- mbSetObj$mNet$gemLib

  if(taxaType == "asv"){

    library(Biostrings)

    seqs = userTable$Microbe
    seqList = DNAStringSet(seqs)
    names(seqList) = paste0("seq", 1:length(seqList))

    for(s in names(seqList)) writeXStringSet(seqList[s], paste0(s, ".fasta"))

    len.gem.lib <- length(gemLib)

    fasta_files <- list.files(pattern = "\\.fasta$")

    if("agora" %in% gemLib | gemLib == "both"){

      if(.on.public.web){
        repSeqFilePath <- "../../lib/agora_ncbi_16s.udb"
      }else{
        repSeqFilePath = "~/Desktop/MicrobiomeNet/MicrobiomeNet_libs/agora_ncbi_16s.udb"
      }

      repSeqDir <- paste0(getwd(), "/")

      if(!file.exists("agora_ncbi_16s.udb")){
        file.copy(repSeqFilePath, repSeqDir)
      }

      resultFile = "vsearch_results_agora.txt"

      command_to_run = paste0("vsearch", " --usearch_global ", repSeqDir, fasta_files, " --db ",
                              repSeqDir, "agora_ncbi_16s.udb", " --id ", identity," --strand both --blast6out ", repSeqDir,
                              fasta_files, resultFile, " --maxaccepts 20 --maxrejects 500")

      run_vsearch(command_to_run)

      asv_df_a <- data.frame(asv = seqs, names = names(seqList))

      asv_matches <- add_matches(asv_df_a, "agora")
      asv_matches[sapply(asv_matches, is.null)] <- NA
      asv_df_a$matches_agora <- unlist(asv_matches)
      asv_df_a <- asv_df_a[, c("asv", "matches_agora")]
      colnames(asv_df_a) <- c("query", "matches")
      asv_df_a$db <- rep("agora")

      if(all(is.na(asv_df_a$matches))){
        AddErrMsg("No matches in GEMs!")
        return(0)
      }

      reg.inx <- match(asv_df_a$query, userTable$Microbe)
      asv_df_a$regulated <- userTable$Regulated[reg.inx]

      if(gemLib == "agora"){
        qs::qsave(asv_df_a, "gem_matches.qs")
      }
    }

    if("carveme" %in% gemLib | gemLib == "both"){

      if(.on.public.web){
        repSeqFilePath <- "../../lib/carveme_asv_16s.udb"
      }else{
        repSeqFilePath = "~/Desktop/MicrobiomeNet/MicrobiomeNet_libs/carveme_asv_16s.udb"
      }

      repSeqDir <- paste0(getwd(), "/")

      if(!file.exists("carveme_asv_16s.udb")){
        file.copy(repSeqFilePath, repSeqDir)
      }

      resultFile = "vsearch_results_carveme.txt"


      command_to_run = paste0("vsearch", " --usearch_global ", repSeqDir, fasta_files, " --db ",
                              repSeqDir, "carveme_asv_16s.udb", " --id ", identity," --strand both --blast6out ", repSeqDir,
                              fasta_files, resultFile, " --maxaccepts 20 --maxrejects 500")

      run_vsearch(command_to_run)

      asv_df_c <- data.frame(asv = seqs, names = names(seqList))

      asv_matches <- add_matches(asv_df_c, "carveme")
      asv_matches[sapply(asv_matches, is.null)] <- NA
      asv_df_c$matches_carveme <- unlist(asv_matches)

      if(.on.public.web){
        carveme_data = data.table::fread("../../lib/taxa_mapping/model_list_processed.txt")
      }else{
        carveme_data = data.table::fread("~/Desktop/Multi-Omics/Pediatric_IBD_Data/ASV_EMBL/data/embl_gems/model_list_processed.txt")
      }

      asv_df_c$matches_carveme <- gsub("_lcl.*", "", asv_df_c$matches_carveme)
      asv_df_c = merge(asv_df_c, carveme_data, by.x = "matches_carveme", by.y = "assembly_accession", all.x = T)
      asv_df_c <- asv_df_c[, c("asv", "ModelID")]
      colnames(asv_df_c) <- c("query", "matches")
      asv_df_c$db <- rep("carveme")

      if(all(is.na(asv_df_c$matches))){
        AddErrMsg("No matches in GEMs!")
        return(0)
      }

      reg.inx <- match(asv_df_c$query, userTable$Microbe)
      asv_df_c$regulated <- userTable$Regulated[reg.inx]

      if(gemLib == "carveme"){
        qs::qsave(asv_df_c, "gem_matches.qs")
      }
    }

    if(gemLib == "both"){
      # combine asv_df tables
      asv_df <- rbind(asv_df_a, asv_df_c)
      qs::qsave(asv_df, "gem_matches.qs")
    }

    # now fast delete the fasta files
    unlink('*fasta')
    unlink('*txt')

  }else if(taxaType == "species"){

    if(.on.public.web){
      gem_otus <- data.table::fread("../../lib/taxa_mapping/GEM_taxa.csv")
    }else{
      gem_otus <- data.table::fread("~/Desktop/MicrobiomeNet/MicrobiomeNet_libs/taxa_mapping/GEM_taxa.csv")
    }

    if("agora" %in% gemLib | gemLib == "both"){

      gem_otus.a <- gem_otus[gem_otus$data == "Agora", ]
      match.inx <- match(userTable$Microbe, gem_otus.a$species)

      # round 2 matching

      for(i in seq_along(match.inx)){
        if(is.na(match.inx[i])){
          match.inx[i] <- grep(userTable$Microbe[i], gem_otus.a$taxname)[1]
        }
      }

      userTable.a <- userTable[!is.na(match.inx), ]

      models <- gem_otus.a[na.omit(match.inx), ]

      if(nrow(models) == 0){
        AddErrMsg("No AGORA matches in GEMs!")
        return(0)
      }

      models <- models[, c(1,3)]
      models.df.a <- cbind(userTable.a[,1], models, userTable.a[,2])
      colnames(models.df.a) <- c("query", "matches", "db", "regulated")
      models.df.a$db <- tolower(models.df.a$db)

    }

    if("carveme" %in% gemLib| gemLib == "both"){
      gem_otus.c <- gem_otus[gem_otus$data == "CarveMe", ]
      match.inx <- match(userTable$Microbe, gem_otus.c$species)

      # round 2 matching

      for(i in seq_along(match.inx)){
        if(is.na(match.inx[i])){
          match.inx[i] <- grep(userTable$Microbe[i], gem_otus.c$taxname)[1]
        }
      }

      userTable.c <- userTable[!is.na(match.inx), ]

      models <- gem_otus.c[na.omit(match.inx), ]

      if(nrow(models) == 0){
        AddErrMsg("No CarveMe matches in GEMs!")
        return(0)
      }

      models <- models[, c(1,3)]
      models.df.c <- cbind(userTable.c[,1], models, userTable.c[,2])
      colnames(models.df.c) <- c("query", "matches", "db", "regulated")
      models.df.c$db <- tolower(models.df.c$db)
    }

    if(gemLib == "both"){
      models.df <- rbind(models.df.a, models.df.c)
    }else if(gemLib == "agora"){
      models.df <- models.df.a
    }else{
      models.df <- models.df.c
    }

    qs::qsave(models.df, "gem_matches.qs")

  }else if(taxaType == "gg_otus"){

    if(.on.public.web){
      gg_otus <- data.table::fread("../../lib/taxa_mapping/taxa_map_ggotu.csv")
    }else{
      gg_otus <- data.table::fread("~/Desktop/MicrobiomeNet/MicrobiomeNet_libs/taxa_mapping/taxa_map_ggotu.csv")
    }

    if(length(gemLib) == 1 & gemLib == "agora"){
      gg_otus <- gg_otus[gg_otus$data == "Agora", ]
    }else if(length(gemLib) == 1 & gemLib == "carveme"){
      gg_otus <- gg_otus[gg_otus$data == "CarveMe", ]
    }

    match.inx <- match(userTable$Microbe, gg_otus$gg_OTU)

    userTable <- userTable[!is.na(match.inx), ]

    models <- gg_otus[na.omit(match.inx), ]

    if(nrow(models) == 0){
      AddErrMsg("No matches in GEMs!")
      return(0)
    }

    models <- models[, c(1,3)]
    models.df <- cbind(userTable[,1], models, userTable[,2])
    colnames(models.df) <- c("query", "matches", "db", "regulated")
    models.df$db <- tolower(models.df$db)

    qs::qsave(models.df, "gem_matches.qs")

  }else if(taxaType == "silva"){

    if(.on.public.web){
      silva_otus_filt <- readRDS("../../lib/taxa_mapping/silva_otus_filt.rds")
    }else{
      silva_otus_filt <- readRDS("~/Desktop/MicrobiomeNet/MicrobiomeNet_libs/taxa_mapping/silva_otus_filt.rds")
    }

    silva_input <- userTable$Microbe

    if(gemLib == "agora" | gemLib == "both"){

      silva_otus_filt.a <- silva_otus_filt[silva_otus_filt$data == "Agora", ]

      matches <- match_silva(silva_input, silva_otus_filt.a)

      userTable.a <- userTable[!is.na(matches), ]

      # microbe_matches <- na.omit(silva_otus_filt$Model_Name[matches])
      microbe_matches_df <- silva_otus_filt.a[na.omit(matches),]

      matches.df.a <- cbind(userTable.a, microbe_matches_df)

      if(nrow(matches.df.a) == 0){
        AddErrMsg("No AGORA matches in GEMs!")
        return(0)
      }

      matches.df.a <- matches.df.a[,c("Microbe", "Model_Name", "data", "Regulated")]
      colnames(matches.df.a) <- c("query", "matches", "db", "regulated")
      matches.df.a$db <- tolower(matches.df.a$db)

    }

    if(gemLib == "carveme" | gemLib == "both"){

      silva_otus_filt.c <- silva_otus_filt[silva_otus_filt$data == "CarveMe", ]

      matches <- match_silva(silva_input, silva_otus_filt.c)

      userTable.c <- userTable[!is.na(matches), ]

      # microbe_matches <- na.omit(silva_otus_filt$Model_Name[matches])
      microbe_matches_df <- silva_otus_filt.c[na.omit(matches),]

      matches.df.c <- cbind(userTable.c, microbe_matches_df)

      if(nrow(matches.df.c) == 0){
        AddErrMsg("No CarveMe matches in GEMs!")
        return(0)
      }

      matches.df.c <- matches.df.c[,c("Microbe", "Model_Name", "data", "Regulated")]
      colnames(matches.df.c) <- c("query", "matches", "db", "regulated")
      matches.df.c$db <- tolower(matches.df.c$db)

    }

    if(gemLib == "both"){
      matches.df <- rbind(matches.df.a, matches.df.c)
    }else if(gemLib == "agora"){
      matches.df <- matches.df.a
    }else{
      matches.df <- matches.df.c
    }

    qs::qsave(matches.df, "gem_matches.qs")

  }else if(taxaType == "ncbi") {
    # fill in

  }else if(taxaType == "gg_tax"){
    #fill in

  }else{
    AddErrMsg("Invalid taxonomy label chosen!")
    return(0)
  }

  return(.set.mbSetObj(mbSetObj))
}

# utility function
# to map to SILVA taxonomy
match_silva <- function(silva_input, silva_mapping){

  silva_input <- gsub("p__Bacteroidetes", "p__Bacteroidota", silva_input)

  if(grepl("d__Bacteria", silva_input[1])){
    new_input <- silva_input
    match.inx <- match(silva_input, silva_mapping$silva_taxon)
  }else{
    new_input <- paste0("d__Bacteria", ";", silva_input)
    match.inx <- match(new_input, silva_mapping$silva_taxon)
  }

  round2 <- is.na(match.inx)

  if(sum(round2)>0){

    new_input <- gsub(".*;s__", "s__", new_input)

    # new_matches.inx <- lapply(new_input, function(x) grep(x, silva_mapping$silva_taxon)[1] )

    new_matches.inx <- sapply(new_input, function(x){

      no_match <- paste("unclassified", "uncultured", "human_gut_metagenome", collapse="|")

      if(grepl(no_match, x)){
        NA
      }else{

        new_inx <- grep(x, silva_mapping$silva_taxon)[1]

        if(is.na(new_inx)){
          x <- gsub("s__", "", x)
          new_inx <- match(x, silva_mapping$Model_Name)
        }
        new_inx
      }

    })
    match.inx <- new_matches.inx
  }

  names(match.inx) <- silva_input
  return(match.inx)
}

# Utility function for ASV matching
run_vsearch <- function(command_to_run){

  for(i in 1:length(command_to_run)){
    system(command_to_run[i])
  }

}

# Utility function for ASV matching
add_matches <- function(asv_df, lib){

  library(data.table)

  matches <- vector("list", length = nrow(asv_df))

  for(i in 1:nrow(asv_df)){

    if(lib == "agora"){
      results = data.table::fread(paste0(asv_df$names[i], ".fastavsearch_results_agora.txt"), header = F)
    }else{
      results = data.table::fread(paste0(asv_df$names[i], ".fastavsearch_results_carveme.txt"), header = F)
    }

    if(nrow(results) > 0){

      colnames(results)[1:6] <- c("seqID", "dbID", "matchPerc", "alnlen", "mism", "gapopens")
      results <- results[,1:6]

      if(nrow(results) > 1){
        max_id <- results$matchPerc == max(results$matchPerc)
        results_keep <- results[max_id,]
        longest_aln <- abs(results_keep$alnlen-max(results_keep$alnlen)) < 5
        results_keep <- results_keep[longest_aln,]
        matches[[i]] = results_keep$dbID[1]
      } else {
        matches[[i]] = results$dbID
      }
    }
  }
  return(matches)
}

GetGemMapRowNames <- function(mbSet){
  gem_matches <- qs::qread("gem_matches.qs")
  return(rownames(gem_matches))
}

GetGemMapCol <-function(mbSet, colInx){
  gem_matches <- qs::qread("gem_matches.qs")
  return(gem_matches[,colInx]);
}



