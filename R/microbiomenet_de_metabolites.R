#' Function to grab metabolits from matched GEMs
#'
#' @param mbSetObj
#' @param nperm Number of permutations, default is set to 100.
#'
#' @description This function creates a metabolite table
#' from matched GEMs as well as creates metabolite tables
#' for permutations (to be used by enrichment analysis).
#' Saves the permuted metabolite table in a list that is
#' saved as a .qs file
GetMetabolitesFromGems <- function(mbSetObj, nperm = 100){

  library(RSQLite)
  library(dplyr)

  mbSetObj <- .get.mbSetObj(mbSetObj)

  matches.df <- qs::qread("gem_matches.qs")
  db <- mbSetObj$mNet$gemLib

  if(db == "agora" | db == "both"){

    agora_matches <- matches.df[matches.df$db == "agora", ]
    a <- na.omit(agora_matches$matches)
    n.gems <- length(a)
    a <- unique(a)

    sqlite.path <- paste0(url.pre, "Agora_metabolites_mucin_trim.sqlite");
    con <- DBI::dbConnect(RSQLite::SQLite(), sqlite.path)
    tbl <- tbl(con, "Agora_metabolites_mucin_trim.sqlite")
    tbl2 <- filter(tbl, Model_Name %in% a) %>% collect()
    agora_gut <- splist2presabs(tbl2, "Model_Name", "Name")
    dbDisconnect(con)

    mbSetObj$mNet$gem.metabolites.a <- agora_gut

    ### code for permutations

    if(.on.public.web){
      gems <- data.table::fread("../../lib/taxa_mapping/GEM_taxa.csv")
    }else{
      gems <- data.table::fread("~/Desktop/MicrobiomeNet/MicrobiomeNet_libs/taxa_mapping/GEM_taxa.csv")
    }

    perm.a <- vector("list", length = nperm)

    gems <- gems[gems$data == "Agora", 1, drop = FALSE]

    perm.a <- lapply(perm.a, function(x) {
      unlist(dplyr::slice_sample(gems, n = n.gems))
    })

    sqlite.path <- paste0(url.pre, "Agora_metabolites_mucin_trim.sqlite");
    con <- DBI::dbConnect(RSQLite::SQLite(), sqlite.path)
    tbl <- tbl(con, "Agora_metabolites_mucin_trim.sqlite")

    perm.mets <- lapply(perm.a, function(x){

      tbl2 <- filter(tbl, Model_Name %in% x) %>% collect()
      agora_gut <- splist2presabs(tbl2, "Model_Name", "Name")

    })

    dbDisconnect(con)

    qs::qsave(perm.mets, "perm_mets_gems_a.qs")
  }

  if(db == "carveme" | db == "both"){

    carveme_matches <- matches.df[matches.df$db == "carveme", ]
    c <- na.omit(carveme_matches$matches)
    n.gems <- length(c)
    c <- unique(c)

    sqlite.path <- paste0(url.pre, "CarveMe_metabolites_gut_trim.sqlite");
    con <- DBI::dbConnect(RSQLite::SQLite(), sqlite.path)
    tbl <- tbl(con, "CarveMe_metabolites_gut_trim.sqlite")
    tbl2 <- filter(tbl, Model_Name %in% c) %>% collect()
    carveme_gut <- splist2presabs(tbl2, "Model_Name", "Name")
    dbDisconnect(con)

    if(nrow(carveme_gut) < length(c)) {

      # get non gut
      sqlite.path <- paste0(url.pre, "CarveMe_metabolites_nongut_trim.sqlite");
      con <- DBI::dbConnect(RSQLite::SQLite(), sqlite.path)
      tbl <- tbl(con, "CarveMe_metabolites_nongut_trim.sqlite")
      tbl3 <- filter(tbl, Model_Name %in% c) %>% collect()
      carveme_nongut <- splist2presabs(tbl3, "Model_Name", "Name")
      dbDisconnect(con)

      # carveme data
      carveme_df <- plyr::rbind.fill(carveme_gut, carveme_nongut)
    }else{
      carveme_df <- carveme_gut
    }

    mbSetObj$mNet$gem.metabolites.c <- carveme_df

    ### code for permutations

    if(.on.public.web){
      gems <- data.table::fread("../../lib/taxa_mapping/GEM_taxa.csv")
    }else{
      gems <- data.table::fread("~/Desktop/MicrobiomeNet/MicrobiomeNet_libs/taxa_mapping/GEM_taxa.csv")
    }

    perm.c <- vector("list", length = nperm)

    gems <- gems[gems$data == "CarveMe", ]
    gems <- gems[gems$type == "Human gut", 1, drop = FALSE]

    perm.c <- lapply(perm.c, function(x) {
      unlist(dplyr::slice_sample(gems, n = n.gems))
    })

    sqlite.path <- paste0(url.pre, "CarveMe_metabolites_gut_trim.sqlite");
    con <- DBI::dbConnect(RSQLite::SQLite(), sqlite.path)
    tbl <- tbl(con, "CarveMe_metabolites_gut_trim.sqlite")

    # Time difference of 16.6392 secs for 44 mets
    perm.mets <- lapply(perm.c, function(x){
      tbl2 <- filter(tbl, Model_Name %in% x) %>% collect()
      agora_gut <- splist2presabs(tbl2, "Model_Name", "Name")
    })

    dbDisconnect(con)

    qs::qsave(perm.mets, "perm_mets_gems_c.qs")

  }
  return(.set.mbSetObj(mbSetObj))
}

# converts long format to wide presence/absence table
# https://rdrr.io/rforge/fuzzySim/man/splist2presabs.html
splist2presabs <-
  function(data, sites.col, sp.col, keep.n = FALSE) {
    stopifnot(
      length(sites.col) == 1,
      length(sp.col) == 1,
      sites.col != sp.col,
      sites.col %in% 1:ncol(data) | sites.col %in% colnames(data),
      sp.col %in% 1:ncol(data) | sp.col %in% colnames(data),
      is.logical(keep.n)
    )
    presabs <- table(data[ , c(sites.col, sp.col)])
    presabs <- as.data.frame(unclass(presabs))
    if (!keep.n)  presabs[presabs > 1] <- 1
    presabs <- data.frame(row.names(presabs), presabs, check.names = F)
    names(presabs)[1] <- names(subset(data, select = sites.col))
    rownames(presabs) <- NULL
    return(presabs)
  }

#' Performs either tests of proportions or binomial tests
#' on the predicted metabolite table
#'
#' @param mbSetObj
#' @param test.type Either 'proportions' to do a test of proportions or
#' 'binomial' to do a binomial test.
#'
#' @export
GetDEMetabolites <- function(mbSetObj, test.type = "proportions"){

  library(tidyverse)

  mbSetObj <- .get.mbSetObj(mbSetObj);

  match.table <- qs::qread("gem_matches.qs")
  match.table <- match.table[!is.na(match.table$matches),]

  model.type <- tolower(mbSetObj$mNet$gemLib)

  # do it for agora then carveme

  if(model.type == "agora" | model.type == "both"){

    metabolites.table <- mbSetObj$mNet$gem.metabolites.a
    metabolites.table[is.na(metabolites.table)] <- 0

    match.table.a <- match.table[match.table$db == "agora", ]

    metabolites.table <- dplyr::left_join(match.table.a, metabolites.table, by = c("matches" = "Model_Name"))
    metabolites.table <- metabolites.table[,-1]
    metabolites.table_sum <- metabolites.table %>% distinct(.) %>%
      select(-c(matches, db)) %>% group_by(regulated) %>% summarise_all(list(sum))

    if(test.type == "tt"){

      tryCatch(
        {      test_res <- apply(metabolites.table[,-c(1:3)], 2, function(x){

          if(var(x) < 0.1){
            return(1)
          }else{
            tt <- t.test(x ~ metabolites.table$regulated)
            return(tt$p.value)
          }
        })}, error = function(e){
          if(!exists("test_res")){
            test_res <- apply(metabolites.table_sum[,-1], 2, function(x){

              if(var(x) < 0.1){
                return(1)
              }else{
                pt <- prop.test(x = x[1], n = nrow(metabolites.table))
                return(pt$p.value)
              }

            })
            test.type <- "proportions"
          }
        }, finally = {

          pval.cutoff <- GetDefaultPvalCutoff(test_res)
          pval_less <- test_res[test_res <= pval.cutoff]
          pval_eq <- test_res[test_res < pval.cutoff]

          pvals <- list(pval_less, pval_eq)
          lens <- c(length(pval_less), length(pval_eq))
          pval.inx <- which(abs(lens-150)==min(abs(lens-150)))
          pvalues <- pvals[[pval.inx[[1]]]]
        }
      )

    }else{

      if(test.type == "proportions"){
        test_prop <- apply(metabolites.table_sum[,-1], 2, function(x){

          if(var(x) < 0.1){
            return(1)
          }else{
            pt <- prop.test(x = x[1], n = nrow(metabolites.table))
            return(pt$p.value)
          }

        })

        pval.cutoff <- GetDefaultPvalCutoff(test_prop)
        pval_less <- test_prop[test_prop <= pval.cutoff]
        pval_eq <- test_prop[test_prop < pval.cutoff]

        pvals <- list(pval_less, pval_eq)
        lens <- c(length(pval_less), length(pval_eq))
        pval.inx <- which(abs(lens-150)==min(abs(lens-150)))
        pvalues <- pvals[[pval.inx[[1]]]]

      }else{

        test_binomial <- apply(metabolites.table_sum[,-1], 2, function(x){

          if(var(x) < 0.1){
            return(1)
          }else{
            bt <- binom.test(x = x[1], n = nrow(metabolites.table), alternative = "two.sided")
            return(bt$p.value)
          }

        })

        pval.cutoff <- GetDefaultPvalCutoff(test_binomial)
        pval_less <- test_binomial[test_binomial <= pval.cutoff]
        pval_eq <- test_binomial[test_binomial < pval.cutoff]

        pvals <- list(pval_less, pval_eq)
        lens <- c(length(pval_less), length(pval_eq))
        pval.inx <- which(abs(lens-150)==min(abs(lens-150)))
        pvalues <- pvals[[pval.inx[[1]]]]

      }
    }

    if(length(pvalues) == 0){
      AddErrMsg("No signficantly different metabolites found!")
      return(0)
    }

    qs::qsave(pvalues, "de_metabolites_a.qs")
    write.table(names(pvalues), file="de_metabolites_a.txt", quote = F, row.names = F, col.names = F)

    mbSetObj$mNet$metabolites.table_sum.a <- metabolites.table_sum
    mbSetObj$mNet$match.reg.a <- metabolites.table$regulated
    mbSetObj$mNet$pval.cutoff.a <- pval.cutoff
    mbSetObj$mNet$test.type <- test.type
  }

  if(model.type == "carveme" | model.type == "both"){

    metabolites.table <- mbSetObj$mNet$gem.metabolites.c
    metabolites.table[is.na(metabolites.table)] <- 0

    match.table.c <- match.table[match.table$db == "carveme", ]

    metabolites.table <- dplyr::left_join(match.table.c, metabolites.table, by = c("matches" = "Model_Name"))
    metabolites.table <- metabolites.table[,-1]
    metabolites.table_sum <- metabolites.table %>% distinct(.) %>%
      select(-c(matches, db)) %>% group_by(regulated) %>% summarise_all(list(sum))

    if(test.type == "tt"){

      tryCatch(
        {      test_res <- apply(metabolites.table[,-c(1:3)], 2, function(x){

          if(var(x) < 0.1){
            return(1)
          }else{
            tt <- t.test(x ~ metabolites.table$regulated)
            return(tt$p.value)
          }
        })}, error = function(e){
          if(!exists("test_res")){
            test_res <- apply(metabolites.table_sum[,-1], 2, function(x){

              if(length(x) == 0){
                return(1)
              }

              if(var(x) < 0.1){
                return(1)
              }else{
                pt <- prop.test(x = x[1], n = nrow(metabolites.table))
                return(pt$p.value)
              }

            })
            test.type <- "proportions"
          }
        }, finally = {

          pval.cutoff <- GetDefaultPvalCutoff(test_res)
          pval_less <- test_res[test_res <= pval.cutoff]
          pval_eq <- test_res[test_res < pval.cutoff]

          pvals <- list(pval_less, pval_eq)
          lens <- c(length(pval_less), length(pval_eq))
          pval.inx <- which(abs(lens-150)==min(abs(lens-150)))
          pvalues <- pvals[[pval.inx[[1]]]]
        }
      )

    }else{
      if(test.type == "proportions"){
        test_prop <- apply(metabolites.table_sum[,-1], 2, function(x){

          if(var(x) < 0.1){
            return(1)
          }else{
            pt <- prop.test(x = x[1], n = nrow(metabolites.table))
            return(pt$p.value)
          }

        })

        pval.cutoff <- GetDefaultPvalCutoff(test_prop)
        pval_less <- test_prop[test_prop <= pval.cutoff]
        pval_eq <- test_prop[test_prop < pval.cutoff]

        pvals <- list(pval_less, pval_eq)
        lens <- c(length(pval_less), length(pval_eq))
        pval.inx <- which(abs(lens-150)==min(abs(lens-150)))
        pvalues <- pvals[[pval.inx[[1]]]]

      }else{

        test_binomial <- apply(metabolites.table_sum[,-1], 2, function(x){

          if(var(x) < 0.1){
            return(1)
          }else{
            bt <- binom.test(x = x[1], n = nrow(metabolites.table), alternative = "two.sided")
            return(bt$p.value)
          }

        })

        pval.cutoff <- GetDefaultPvalCutoff(test_binomial)
        pval_less <- test_binomial[test_binomial <= pval.cutoff]
        pval_eq <- test_binomial[test_binomial < pval.cutoff]

        pvals <- list(pval_less, pval_eq)
        lens <- c(length(pval_less), length(pval_eq))
        pval.inx <- which(abs(lens-150)==min(abs(lens-150)))
        pvalues <- pvals[[pval.inx[[1]]]]

      }
    }

    if(length(pvalues) == 0){
      AddErrMsg("No signficantly different metabolites found!")
      return(0)
    }

    qs::qsave(pvalues, "de_metabolites_c.qs")
    write.table(names(pvalues), file="de_metabolites_c.txt", quote = F, row.names = F, col.names = F)

    mbSetObj$mNet$metabolites.table_sum.c <- metabolites.table_sum
    mbSetObj$mNet$match.reg.c <- metabolites.table$regulated
    mbSetObj$mNet$pval.cutoff.c <- pval.cutoff
  }

  return(.set.mbSetObj(mbSetObj))
}

#' Function that gets metabolites for permutations
#' Speed can be improved...
GetPermDEMetabolites <- function(mbSetObj){

  mbSetObj <- .get.mbSetObj(mbSetObj);

  if(mbSetObj$mNet$gemLib == "both"){

    perm_libs <- c("perm_mets_gems_a.qs", "perm_mets_gems_c.qs")
    perm_res <- c("permutation_metabolites_a.qs", "permutation_metabolites_c.qs")

    mets.reg <- list(mbSetObj$mNet$match.reg.a, mbSetObj$mNet$match.reg.c)
    pval.cutoff <- c(mbSetObj$mNet$pval.cutoff.a, mbSetObj$mNet$pval.cutoff.c)
    test.type <- c(mbSetObj$mNet$test.type.a, mbSetObj$mNet$test.type.c)

  }else if(mbSetObj$mNet$gemLib == "agora"){

    perm_libs <- "perm_mets_gems_a.qs"
    perm_res <- "permutation_metabolites_a.qs"

    mets.reg <- mbSetObj$mNet$match.reg.a
    pval.cutoff <- mbSetObj$mNet$pval.cutoff
    test.type <- mbSetObj$mNet$test.type

  }else{

    perm_libs <- "perm_mets_gems_c.qs"
    perm_res <- "permutation_metabolites_c.qs"

    mets.reg <- mbSetObj$mNet$match.reg.c
    pval.cutoff <- mbSetObj$mNet$pval.cutoff
    test.type <- mbSetObj$mNet$test.type

  }

  for(i in seq_along(perm_libs)){

    perm_mets <- qs::qread(perm_libs[i])
    num_perm <- length(perm_mets)

    set.seed(123)

    print(paste('Resampling, ', num_perm, 'permutations to estimate background ...', "round", i));
    permutation_mets <- vector("list", num_perm);

    permutation_mets <- lapply(seq_along(permutation_mets), function(x){

      metabolites.table <- perm_mets[[x]]

      if(class(mets.reg) == "list"){
        mets.reg <- mets.reg[[i]][1:nrow(metabolites.table)]
      }else{
        mets.reg <- mets.reg[1:nrow(metabolites.table)]
      }

      metabolites.table <- add_column(metabolites.table, Regulated = mets.reg, .after = "Model_Name")
      metabolites.table_sum <- metabolites.table[,-1] %>% group_by(Regulated) %>% summarise_all(list(sum))
      metabolites.table_sum[is.na(metabolites.table_sum)] <- 0

      tryCatch(
        {  test_res <- apply(metabolites.table[,-c(1:3)], 2, function(x){

          if(var(x) < 0.1){
            return(1)
          }else{
            tt <- t.test(x ~ metabolites.table$regulated)
            return(tt$p.value)
          }
        })}, error = function(e){
          if(!exists("test_res")){
            test_res <- apply(metabolites.table_sum[,-1], 2, function(x){

              tryCatch({
                if(var(x) < 0.1){
                  return(1)
                }else{
                  pt <- prop.test(x = x[1], n = nrow(metabolites.table))
                  return(pt$p.value)
                }
              }, error = function(e){
                return(1)
              })

            })
          }
        }, finally = {

          if(!exists("test_res")){
            test_res <- rep(1, nrow(metabolites.table))
          }

          pvalues <- test_res[test_res < pval.cutoff[[i]]]
        }
      )

      # if(test.type[[i]] == "proportions"){
      #   test_prop <- apply(metabolites.table_sum[,-1], 2, function(x){
      #
      #     if(var(x) < 0.1){
      #       return(1)
      #     }else{
      #       pt <- prop.test(x = x[1], n = nrow(metabolites.table))
      #       return(pt$p.value)
      #     }
      #
      #   })
      #   pvalues <- test_prop[test_prop <= pval.cutoff[[i]]]
      # }else{
      #
      #   test_binomial <- apply(metabolites.table_sum[,-1], 2, function(x){
      #
      #     if(var(x) < 0.1){
      #       return(1)
      #     }else{
      #       bt <- binom.test(x = x[1], n = nrow(metabolites.table), alternative = "two.sided")
      #       return(bt$p.value)
      #     }
      #
      #   })
      #   pvalues <- test_binomial[test_binomial <= pval.cutoff[[i]]]
      # }
    })

    qs::qsave(permutation_mets, perm_res[i])
  }

  return(.set.mbSetObj(mbSetObj))
}

#' Wrapper function that performs enrichment analysis
#'
#' @param mbSetObj
#' @param library.name Input the name of the pathway library, without extension.
#' @param model Either 'both', 'agora', or 'carveme'.
#' @param input.type
PerformTrueEnrichment <- function(mbSetObj, library.name,
                                  model = T, input.type = "name"){

  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$mNet$lib <- library.name

  gemLib <- mbSetObj$mNet$gemLib

  if(gemLib == "both"){
    text.file <- c("de_metabolites_a.txt", "de_metabolites_c.txt")
    record <- list(qs::qread("permutation_hits_a.qs"), qs::qread("permutation_hits_c.qs"))
  }else if(gemLib == "agora"){
    text.file <- "de_metabolites_a.txt"
    record <- list(qs::qread("permutation_hits_a.qs"))
  }else{
    text.file <- "de_metabolites_c.txt"
    record <- list(qs::qread("permutation_hits_c.qs"))
  }

  for(i in seq_along(text.file)){

    mSet <-InitDataObjects("conc", "msetora", FALSE)
    cmpd.vec <- unique(scan(text.file[[i]], character(), sep = "\n"))
    mSet <- Setup.MapData(mSet, cmpd.vec);
    mSet <- CrossReferencing(mSet, input.type, model = model);
    mSet <- SetMetabolomeFilter(mSet, F);
    mSet <- CalculateTrueHyperScore(mSet, model = model, library.name = library.name, record = record[[i]])

    if(text.file[[i]] == "de_metabolites_a.txt"){
      mbSetObj$mNet$resMat <- mSet$analSet$ora.mat
      fast.write.csv(mSet$analSet$ora.mat, file="msea_ora_result_a.csv");
    }else{
      mbSetObj$mNet$resMat <- mSet$analSet$ora.mat
      fast.write.csv(mSet$analSet$ora.mat, file="msea_ora_result_c.csv");
    }
  }

  return(.set.mbSetObj(mbSetObj))
}

GetMnetResRowNames <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(rownames(mbSetObj$mNet$resMat));
}

GetMnetResColNames <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(colnames(mbSetObj$mNet$resMat));
}

GetMnetResMat <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(as.matrix(mbSetObj$mNet$resMat));
}

#' Function that performs enrichment analysis
#' on permuted metabolites
#'
#' @param mbSetObj
#' @param library.name Input the name of the pathway library, without extension.
PerformPermEnrichment <- function(mbSetObj, library.name){

  mbSetObj <- .get.mbSetObj(mbSetObj);

  gemLib <- mbSetObj$mNet$gemLib

  if(gemLib == "both"){
    perm_mets <- list(qs::qread("permutation_metabolites_a.qs"), qs::qread("permutation_metabolites_c.qs"))
    perm_res <- c("permutation_hits_a.qs", "permutation_hits_c.qs")
  }else if(gemLib == "agora"){
    perm_mets <- list(qs::qread("permutation_metabolites_a.qs"))
    perm_res <- "permutation_hits_a.qs"
  }else{
    perm_mets <- list(qs::qread("permutation_metabolites_c.qs"))
    perm_res <- "permutation_hits_c.qs"
  }

  SetCurrentMsetLib(library.name, 2);

  for(i in seq_along(perm_mets)){

    num_perm <- length(perm_mets[[i]])
    permutation_hits <- vector("list", num_perm);

    for(j in seq_along(perm_mets[[i]])){

      cmpds <- names(perm_mets[[i]][[j]])

      if(length(cmpds) == 0){
        res <- list(matrix(1, nrow = nrow(current.msetlib), ncol = 1))
      }else{
        res <- .perform_perm_enrichment(cmpds, library.name)
      }

      permutation_hits[[j]] <- res
    }

    # end.for-start.for
    # Time difference of 6.294494 secs

    qs::qsave(permutation_hits, perm_res[[i]])
  }

  return(.set.mbSetObj(mbSetObj))
}

.perform_perm_enrichment <- function(cmpd.vec, library.name, model = T, input.type = "name"){

  mSet <- InitDataObjects("conc", "msetora", FALSE)
  mSet <- Setup.MapData(mSet, cmpd.vec);

  mSet <- tryCatch({
    CrossReferencing(mSet, input.type, model = model);
  }, error = function(e) {
    return(mSet)
  })

  if('name.map' %in% names(mbSet)){
    if(sum(is.na(mSet$name.map$hit.inx))/length(mSet$name.map$hit.inx) > 0.75){
      AddErrMsg("Too many missing IDs!")
      return(list(matrix(1, nrow = nrow(current.msetlib), ncol = 1)))
    }else{
      mSet <- SetMetabolomeFilter(mSet, F);
      mSet <- CalculatePermHyperScore(mSet, model = model, perm = T, library.name)
      res <- list(mSet$analSet$ora.mat)
      return(res)
    }
  }else{
    return(list(matrix(1, nrow = nrow(current.msetlib), ncol = 1)))
  }
}

###################################################################

#' Gets default p-val with reasonable cutoff
GetDefaultPvalCutoff <- function(test.pvals){

  pvals <- c(0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001,
             1.0e-6, 1.0e-7, 1.0e-8, 1.0e-9, 1.0e-10, 1.0e-11, 1.0e-12,
             1.0e-13, 1.0e-14, 1.0e-15, 1.0e-16, 1.0e-17, 1.0e-18, 1.0e-19, 1.0e-20)

  ndat <- test.pvals
  ndat <- ndat[order(ndat)]

  percP <- ifelse(0.1*length(ndat) > 200, 200/length(ndat), 0.1)
  n <- floor(percP*length(ndat))
  cutoff <- ndat[n+1]

  if(sum(test.pvals < cutoff) > 100){
    return(cutoff)
  }

  level <- (attr(regexpr("(?<=\\.)0+", cutoff, perl = TRUE), "match.length")) + 5

  maxp <- ceiling_dec(cutoff, level)

  if(sum(maxp) < 275){
    return(maxp)
  }else{
    return(0.05)
  }
}

ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)
