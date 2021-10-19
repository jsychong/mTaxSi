#' PCA and t-SNE plots of predicted metabolomics data
#'
#' @param mbSetObj
#' @param type Either 'pca' or 'tsne'.
#' @param imgName Character, the name of the plot without an extension.
#'
#' @export
global_visualization <- function(mbSetObj, type = "pca", imgName = imgName){

  mbSetObj <- .get.mbSetObj(mbSetObj);

  group.label <- mbSetObj$mNet$pheno
  group.label <- paste0(group.label, collapse = " vs. ")
  model.type <- mbSetObj$mNet$gemLib

  match.table <- qs::qread("gem_matches.qs")

  if(model.type == "agora" | model.type == "both"){

    metabolites.table <- mbSetObj$mNet$gem.metabolites.a
    metabolites.table[is.na(metabolites.table)] <- 0

    fast.write.csv(metabolites.table, "metabolite_table.csv", row.names = F)

    match.table.a <- match.table[match.table$db == "agora", ]
    match.table.inx <- match(metabolites.table$Model_Name, match.table.a$matches)
    match.table.a <- match.table.a[match.table.inx, ]

    if(type == "pca"){

      # first PCA score

      pca_matrix <- metabolites.table[, -1]
      rownames(pca_matrix) <- metabolites.table[,1]

      metabolites_pca <- prcomp(pca_matrix)
      df_out <- as.data.frame(metabolites_pca$x)
      df_out$Group <- factor(toupper(match.table.a$regulated), levels = c("DOWN", "UP"))
      df_out$Microbe <- rownames(df_out)

      var_importance <- summary(metabolites_pca)$importance[2,]
      var_importance <- paste0(names(var_importance), " (", var_importance*100, "%)")

      library(ggplot2)
      library(viridis)

      g.title <- c("PCA of AGORA Metabolic Content")
      g.title.save <- paste0(imgName, "agora", ".png")

      p <- ggplot(df_out, aes(x = PC1, y = PC2, color = Group, label = Microbe)) +
        geom_point(size=5) + theme_bw() + scale_color_manual(values=c("DOWN" = "#387590", "UP" = "#E64242")) +
        xlab(var_importance[1]) + ylab(var_importance[2]) + ggtitle(g.title) +
        labs(color=group.label) + theme(legend.text=element_text(size=9))

      ggsave(g.title.save, p, height = 7, width = 8)

      library(plotly)
      pp <- ggplotly(p)

      p.title.save <- paste0(imgName, "agora", ".html")

      htmlwidgets::saveWidget(pp, p.title.save)

      # loading plot
      loadings <- metabolites_pca$rotation
      loadings <- as.data.frame(loadings[,1:2])
      loadings$cpds <- rownames(loadings)

      library(ggplot2)
      library(viridis)

      g.title <- c("PCA Loadings of AGORA Metabolic Content")
      g.title.save <- paste0("pca_loading_", "agora", ".png")

      p <- ggplot(loadings, aes(x = PC1, y = PC2, label = cpds)) +
        geom_point(size=2.5) + theme_bw() +
        xlab("PC1") + ylab("PC2") + ggtitle(g.title) +
        theme(legend.text=element_text(size=9))

      ggsave(g.title.save, p, height = 7, width = 8)

      library(plotly)
      pp <- ggplotly(p)

      p.title.save <- paste0("pca_loading_", "agora", ".html")

      htmlwidgets::saveWidget(pp, p.title.save)


    }else{
      # t-Sne

      library(Rtsne)
      set.seed(123)

      rtsne_out <- tryCatch(Rtsne(as.matrix(pca_matrix), pca = FALSE, verbose = T,
                                  perplexity = floor(nrow(pca_matrix)/3)), error=function(e) NULL)

      if(is.null(rtsne_out)){
        current.msg <<- c("t-SNE incompatible with dataset!")
        return(0)
      }

      df_out <- as.data.frame(rtsne_out$Y)
      df_out$Group <- factor(toupper(match.table.a$regulated), levels = c("DOWN", "UP"))
      df_out$Microbe <- rownames(df_out)

      g.title <- c("tSNE of AGORA Metabolic Content")
      g.title.save <- paste0(imgName, "agora", ".png")

      pt <- ggplot(df_out, aes(x = V1, y = V2, color = Group, label = Microbe)) +
        geom_point(size=3) + theme_bw() + scale_color_manual(values=c("UP" = "#E64242", "DOWN" = "#387590")) +
        xlab("Dimension 1") + ylab("Dimension 2") + ggtitle(g.title) +
        labs(color=group.label) + theme(legend.text=element_text(size=9))

      ggsave(g.title.save, pt, height = 7, width = 8)

      ppt <- ggplotly(pt)

      p.title.save <- paste0(imgName, "agora", ".html")

      htmlwidgets::saveWidget(ppt, p.title.save)
    }

  }

  if(model.type == "carveme" | model.type == "both"){

    metabolites.table <- mbSetObj$mNet$gem.metabolites.c
    metabolites.table[is.na(metabolites.table)] <- 0

    fast.write.csv(metabolites.table, "metabolite_table.csv", row.names = F)

    match.table.c <- match.table[match.table$db == "carveme", ]
    match.table.inx <- match(metabolites.table$Model_Name, match.table.c$matches)
    match.table.c <- match.table.c[match.table.inx, ]

    if(type == "pca"){
      # first PCA

      pca_matrix <- metabolites.table[, -1]
      rownames(pca_matrix) <- metabolites.table[,1]

      metabolites_pca <- prcomp(pca_matrix)
      df_out <- as.data.frame(metabolites_pca$x)
      df_out$Group <- factor(toupper(match.table.c$regulated), levels = c("DOWN", "UP"))
      df_out$Microbe <- rownames(df_out)

      var_importance <- summary(metabolites_pca)$importance[2,]
      var_importance <- paste0(names(var_importance), " (", var_importance*100, "%)")

      library(ggplot2)
      library(viridis)

      g.title <- c("PCA of CarveMe Metabolic Content")
      g.title.save <- paste0(imgName, "carveme", ".png")

      p <- ggplot(df_out, aes(x = PC1, y = PC2, color = Group, label = Microbe)) +
        geom_point(size=5) + theme_bw() + scale_color_manual(values=c("UP" = "#E64242", "DOWN" = "#387590")) +
        xlab(var_importance[1]) + ylab(var_importance[2]) + ggtitle(g.title) +
        labs(color=group.label) + theme(legend.text=element_text(size=9))

      ggsave(g.title.save, p, height = 7, width = 8)

      library(plotly)
      pp <- ggplotly(p)

      p.title.save <- paste0(imgName, "carveme", ".html")

      htmlwidgets::saveWidget(pp, p.title.save)

      # loading plot
      loadings <- metabolites_pca$rotation
      loadings <- as.data.frame(loadings[,1:2])
      loadings$cpds <- rownames(loadings)

      library(ggplot2)
      library(viridis)

      g.title <- c("PCA Loadings of CarveMe Metabolic Content")
      g.title.save <- paste0("pca_loading_", "carveme", ".png")

      p <- ggplot(loadings, aes(x = PC1, y = PC2, label = cpds)) +
        geom_point(size=2.5) + theme_bw() +
        xlab("PC1") + ylab("PC2") + ggtitle(g.title) +
        theme(legend.text=element_text(size=9))

      ggsave(g.title.save, p, height = 7, width = 8)

      library(plotly)
      pp <- ggplotly(p)

      p.title.save <- paste0("pca_loading_", "carveme", ".html")

      htmlwidgets::saveWidget(pp, p.title.save)

    }else{
      # t-Sne

      library(Rtsne)

      set.seed(123)

      rtsne_out <- tryCatch(Rtsne(as.matrix(pca_matrix), pca = FALSE, verbose = T,
                                  perplexity = floor(nrow(pca_matrix)/3)), error=function(e) NULL)

      if(is.null(rtsne_out)){
        current.msg <<- c("t-SNE incompatible with dataset!")
        return(0)
      }

      df_out <- as.data.frame(rtsne_out$Y)
      df_out$Group <- factor(toupper(match.table.c$regulated), levels = c("DOWN", "UP"))
      df_out$Microbe <- rownames(df_out)

      g.title <- c("tSNE of CarveMe Metabolic Content")
      g.title.save <- paste0(imgName, "carveme", ".png")

      pt <- ggplot(df_out, aes(x = V1, y = V2, color = Group, label = Microbe)) +
        geom_point(size=3) + theme_bw() + scale_color_manual(values=c("UP" = "#E64242","DOWN" = "#387590")) +
        xlab("Dimension 1") + ylab("Dimension 2") + ggtitle(g.title) +
        labs(color=group.label) + theme(legend.text=element_text(size=9))

      ggsave(g.title.save, pt, height = 7, width = 8)

      ppt <- ggplotly(pt)

      p.title.save <- paste0(imgName, "carveme", ".html")

      htmlwidgets::saveWidget(ppt, p.title.save)
    }
  }

  return(.set.mbSetObj(mbSetObj))
}

#' Creates heatmap of bioactive compounds
#'
#' @param mbSetObj
#' @param metType Either 'bioactive' for bioactive metabolites or 'amino' for amino acids.
#' @param plotName Default is 'heatmap'.
#' @param order.by 'Hits' to order by the number of hits or else it will order the metabolites alphabetically.
#' @param width Default is 7.
#' @param height Default is 7.
#'
#' @description Heatmap can be of either known microbial-origin
#' bioactive compounds or amino acids.
#' @export
heatmap_bioactive_metabolites <- function(mbSetObj, metType = "bioactive",
                                          plotName = "heatmap", order.by = "hits",
                                          width = 7, height = 7){

  mbSetObj <- .get.mbSetObj(mbSetObj)

  if(metType == "bioactive"){
    amino.acids = FALSE
  }else{
    amino.acids = TRUE
  }

  gemLib <- mbSetObj$mNet$gemLib

  group.label <- mbSetObj$mNet$pheno
  group.label <- paste0(group.label, collapse = " vs. ")
  model.type <- mbSetObj$mNet$gemLib

  if(gemLib == "both"){
    metabolites.table_sum <- list(mbSetObj$mNet$metabolites.table_sum.a, mbSetObj$mNet$metabolites.table_sum.c)

    if(amino.acids){
      plot_name <- c("heatmap_amino_mets_agora.png", "heatmap_amino_mets_carveme.png")
    }else{
      plot_name <- c("heatmap_bioactive_mets_agora.png", "heatmap_bioactive_mets_carveme.png")
    }


  }else if(gemLib == "agora"){
    metabolites.table_sum <- list(mbSetObj$mNet$metabolites.table_sum.a)
    plot_name <- paste0(plotName, ".png")
  }else{
    metabolites.table_sum <- list(mbSetObj$mNet$metabolites.table_sum.c)
    plot_name <- paste0(plotName, ".png")
  }

  if(.on.public.web & amino.acids){
    bioactive.metabolites <- qs::qread("../../lib/plotting/bioactive_metabolites_amino.qs")
  }else if(.on.public.web & !amino.acids){
    bioactive.metabolites <- qs::qread("../../lib/plotting/bioactive_metabolites_no_amino.qs")
  }else if(!.on.public.web & amino.acids){
    bioactive.metabolites <- qs::qread("~/Desktop/MicrobiomeNet/MicrobiomeNet_libs/bioactive_metabolites_amino.qs")
  }else{
    bioactive.metabolites <- qs::qread("~/Desktop/MicrobiomeNet/MicrobiomeNet_libs/bioactive_metabolites_no_amino.qs")
  }

  bioactive.mets <- bioactive.metabolites[,1]

  for(i in seq_along(metabolites.table_sum)){

    match.table.inx <- match(bioactive.mets, colnames(metabolites.table_sum[[i]]))

    if(length(na.omit(match.table.inx)) < 5){
      print("Not enough metabolites matched to bioactive metabolites!")
      return(0)
    }

    metabolites.table_sum_matched <- metabolites.table_sum[[i]][, c(1, na.omit(match.table.inx))]
    metabolites.table_sum_matched <- t(metabolites.table_sum_matched)
    colnames(metabolites.table_sum_matched) <- metabolites.table_sum_matched[1,]
    metabolites.table_sum_matched <- metabolites.table_sum_matched[-1, ]
    metabolites.table.df <- reshape2::melt(metabolites.table_sum_matched)
    metabolites.table.df <- metabolites.table.df[order(as.numeric(trimws(metabolites.table.df$value)), decreasing = T), ]

    if(order.by == "hits"){
      metabolites.table.df$Var1 <- factor(metabolites.table.df$Var1, levels = rev(as.character(unique(metabolites.table.df$Var1))))
    }else{
      metabolites.table.df$Var1 <- factor(metabolites.table.df$Var1, levels = unique(bioactive.metabolites[,1]))
    }

    metabolites.table.df$value <- as.numeric(trimws(metabolites.table.df$value))

    group.inx <- match(as.character(metabolites.table.df$Var1), bioactive.metabolites[,1])

    anot.df <- data.frame(x = metabolites.table.df$Var1, y = "Class", group = bioactive.metabolites[group.inx, 2])

    library(ggplot2)
    library(RColorBrewer)
    library(ggnewscale)

    cols <- c("#003f5c", "#7F3C8D", "#11A579", "#E68310", "#FFFD7D", "#A5AA99")

    p <- ggplot(metabolites.table.df, aes(Var2, Var1)) +
      #add border white colour of line thickness 0.25
      geom_tile(aes(fill = value), colour="white", size = 0.25) +
      coord_equal() +
      xlab(paste0("\n", group.label)) + ylab("Metabolite\n") + ggtitle("Heatmap of\nBioactive Metabolites") +
      labs(fill = "Total\nCount") +
      scale_fill_distiller(palette = "RdBu") +

      # Set new scale fill after you've specified the scale for the heatmap
      new_scale_fill() +
      geom_tile(data = anot.df, aes(y, x, fill = group),
                width = 0.5, colour = "white") +
      scale_fill_manual(values = cols, name = "Metabolite\nClass") +    #remove extra space
      scale_y_discrete(expand=c(0,0)) +
      scale_x_discrete(expand=c(0,0), labels = c("", "Down", "Up")) +
      theme(
        legend.title = element_text(size=10, face="bold"),
        legend.text = element_text(size=10),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0.75, unit = 'cm'),
        axis.title = element_text(size = 10.5, face="bold"),
        axis.text = element_text(size = 10),
        #set thickness of axis ticks
        axis.ticks = element_blank(),
        #remove plot background
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        #remove plot border
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.6))

    ggsave(plot_name[[i]], p, width = width, height = height)

  }

  return(.set.mbSetObj(mbSetObj))
}

#' Creates a bubble plot of bioactive
#' metabolite presence or absence
#' from a group of metabolites
#'
#' @param mbSetObj
#' @param metType Either 'bioactive' for bioactive metabolites, 'amino' for amino acids,
#' or the pathway name from the KEGG pathways.
#' @param width Default is 12.
#' @param height Default is 7.
#' @param bubble_size Default is 5.
#'
#' @export
microbe_metabolite_bubble_chart <- function(mbSetObj, metType = "bioactive",
                                            width = 12, height = 7, bubble_size = 5){

  mbSetObj <- .get.mbSetObj(mbSetObj);

  gemLib <- mbSetObj$mNet$gemLib

  match.table <- qs::qread("gem_matches.qs")

  if(gemLib == "both"){
    metabolites.table <- list(mbSetObj$mNet$gem.metabolites.a, mbSetObj$mNet$gem.metabolites.c)
    match.table <- list(match.table[match.table$db == "agora", ], match.table[match.table$db == "carveme", ])
    plot_name <- c(paste0("bubble_", metType, "_metabolites_agora.png"), paste0("bubble_", metType, "_metabolites_carveme.png"))
  }else if(gemLib == "agora"){
    metabolites.table <- list(mbSetObj$mNet$gem.metabolites.a)
    match.table <- list(match.table)
    plot_name <- paste0("bubble_", metType, "_metabolites_agora.png")
  }else{
    metabolites.table <- list(mbSetObj$mNet$gem.metabolites.c)
    match.table <- list(match.table)
    plot_name <- paste0("bubble_", metType, "_metabolites_carveme.png")
  }

  if(metType %in% c("bioactive", "amino")){

    if(.on.public.web & metType == "amino"){
      bioactive.metabolites <- qs::qread("../../lib/plotting/bioactive_metabolites_amino.qs")
    }else if(.on.public.web & metType == "bioactive"){
      bioactive.metabolites <- qs::qread("../../lib/plotting/bioactive_metabolites_no_amino.qs")
    }else if(!.on.public.web & metType == "amino"){
      bioactive.metabolites <- qs::qread("~/Desktop/MicrobiomeNet/MicrobiomeNet_libs/bioactive_metabolites_amino.qs")
    }else{
      bioactive.metabolites <- qs::qread("~/Desktop/MicrobiomeNet/MicrobiomeNet_libs/bioactive_metabolites_no_amino.qs")
    }

    bioactive.mets <- bioactive.metabolites[,1]

  }else{

    model_mets <- qs::qread("model_metabolites_master.qs")
    kegg <- qs::qread("kegg_pathways_generic.qs")

    cpds <- unlist(strsplit(kegg$member[match(metType, kegg$name)] , "; "))
    bioactive.mets <- model_mets$name[na.omit(match(cpds, model_mets$kegg))]

    bioactive.metabolites <- data.frame(value = bioactive.mets,
                                        group = metType)

  }

  for(i in seq_along(metabolites.table)){

    metabolites.table.a <- metabolites.table[[i]]
    match.table.a <- match.table[[i]]
    metabolites.table.a[is.na(metabolites.table.a)] <- 0
    match.table.inx <- match(metabolites.table.a$Model_Name, match.table.a$matches)
    match.table.a <- match.table.a[match.table.inx, ]

    metabolites.table.a <- add_column(metabolites.table.a, Regulated = match.table.a$regulated, .after = "Model_Name")

    match.table.inx <- match(bioactive.mets, colnames(metabolites.table.a))

    if(length(na.omit(match.table.inx)) < 5){
      print("Not enough metabolites matched to bioactive metabolites!")
      return(0)
    }

    if(length(match.table.inx) < 12){
      height <- 7
    }else if(length(match.table.inx) < 18){
      height <- 8
    }else if(length(match.table.inx) < 24){
      height <- 9
    }else if(length(match.table.inx) < 32){
      height <- 10
    }else{
      height <- height
    }

    metabolites.table.a <- metabolites.table.a[, c(1, na.omit(match.table.inx))]
    metabolites.table.a <- t(metabolites.table.a)
    colnames(metabolites.table.a) <- metabolites.table.a[1,]
    metabolites.table.a <- metabolites.table.a[-1, ]
    metabolites.table.df <- reshape2::melt(metabolites.table.a)

    metabolites.table.df <- metabolites.table.df[(metabolites.table.df$value == "1"), ]
    group.inx <- match(as.character(metabolites.table.df$Var1), bioactive.metabolites[,1])
    metabolites.table.df$group = factor(bioactive.metabolites[group.inx, 2])
    metabolites.table.df <- metabolites.table.df[order(metabolites.table.df$group),]
    metabolites.table.df$regulated <- match.table.a$regulated[match(as.character(metabolites.table.df$Var2), match.table.a$matches)]
    metabolites.table.df <- metabolites.table.df[order(metabolites.table.df$regulated), ]
    metabolites.table.df$Var2 <- factor(metabolites.table.df$Var2, levels = unique(as.character(metabolites.table.df$Var2)))
    col_regulated <- match.table.a$regulated[match(levels(metabolites.table.df$Var2), match.table.a$matches)]
    col_regulated <- ifelse(tolower(col_regulated) == "down", "#000058", "#8b0000")

    groups <-  metabolites.table.df[,c(1,4)] %>%
      group_by(group) %>%
      slice(c(1, n())) %>%
      ungroup()

    groups_df <- data.frame(xmin = groups$Var1[c(TRUE, FALSE)],
                            xmax = groups$Var1[c(FALSE, TRUE)],
                            group = unique(groups$group))

    library(ggplot2)
    library(ggtext)
    library(RColorBrewer)

    cols <- c("#003f5c", "#7F3C8D", "#11A579", "#E68310", "#b8b633", "#575c4d")

    p <- ggplot(data = metabolites.table.df, mapping = aes(x = Var1, y = Var2,
                                                           fill = value, color = group)) +
      geom_point(alpha = 0.75, size = bubble_size) +
      scale_color_manual(values = cols, name = "Metabolite\nClass") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      theme(axis.text.y = element_markdown(colour = col_regulated)) +
      guides(fill = FALSE) +
      xlab("Bioactive Metabolite") + ylab("Microbes\n")

    ggsave(plot_name[[i]], p, width = width, height = height)
  }

  return(.set.mbSetObj(mbSetObj))
}

#' Function to create lollipop comparison chart
#' of pathway enrichment results from AGORA
#' versus CarveMe
#'
#' @param mbSetObj
#' @param plotName Character, the name of the plot.
#' @param anot Add annotations for the pathways.
#' @param width Default is 7.
#' @param height Default is 4.
#'
#' @export
lollipop_plot_enrichment <- function(mbSetObj, plotName,
                                     anot = T, width = 7, height = 4){

  library(data.table)
  library(tidyverse)
  library(ggplot2)

  mbSetObj <- .get.mbSetObj(mbSetObj);

  lib.type <- ifelse(grep("kegg", mbSet$mNet$lib), "kegg", "metacyc")

  gemLib <- mbSetObj$mNet$gemLib

  if(gemLib == "both"){
    AddErrMsg("This function is only intended for a single GEM library")
    return(0)
  }else if(gemLib == "agora"){
    super_generic <- read.csv("msea_ora_result_a.csv")
  }else{
    super_generic <- read.csv("msea_ora_result_c.csv")
  }

  # hits/expected = enrichment ratio
  super_generic$folds = super_generic[,3]/super_generic[,2]
  super_generic$logp = -log10(super_generic[,6])
  super_generic$logp <- scales::rescale(super_generic$logp, to = c(0, 5))

  if(nrow(super_generic) > 50){
    super_generic <- super_generic[1:50,]
  }

  super_generic <- super_generic[order(super_generic$logp), ]

  super_generic$Pathways <- factor(rownames(super_generic), levels = unique(rownames(super_generic)))
  super_generic <- super_generic %>% filter(logp > 0)

  if(anot){

    if(lib.type == "kegg"){

      if(.on.public.web){
        kegg_hier <- read.csv("../../lib/plotting/kegg_brite_hierarchy_pathways.csv")
      }else{
        kegg_hier <- read.csv("~/Desktop/MicrobiomeNet/pathways/kegg_brite_hierarchy_pathways.csv")
      }

      kegg_hier$Pathway_Name <- trimws(substring(kegg_hier$Pathway_Name, 6))
      hier.inx <- match(as.character(super_generic$Pathways), kegg_hier$Pathway_Name)
      anot.df <- data.frame(x = super_generic$Pathways, y = 0, hier = kegg_hier[hier.inx, 2] )

      cols <- c("#5F4690","#1D6996","#38A6A5","#0F8554",
                "#73AF48","#EDAD08","#E17C05","#CC503E","#94346E","#6F4070","#994E95","#666666")
    }else{

      if(.on.public.web){
        meta_hier <- readRDS("../../lib/plotting/meta_cyc_paths_hier.rds")
      }else{
        meta_hier <- readRDS("~/Desktop/MicrobiomeNet/pathways/MetaCyc/2021/meta_cyc_paths_hier.rds")
      }

      meta_hier$Super.Pathways <- gsub("\\s*\\([^\\)]+\\)","",as.character(meta_hier$Super.Pathways))
      meta_hier$Super.Pathways <- gsub("//.*","",as.character(meta_hier$Super.Pathways))
      meta_hier$Super.Pathways <- gsub("superpathway of","",as.character(meta_hier$Super.Pathways))
      meta_hier$Super.Pathways[meta_hier$Super.Pathways == ""] <- "Unclassified"
      meta_hier$Super.Pathways <- tools::toTitleCase(trimws(meta_hier$Super.Pathways))

      hier.inx <- match(as.character(super_generic$Pathways), meta_hier$Pathways)
      anot.df <- data.frame(x = super_generic$Pathways, y = -6, hier = meta_hier[hier.inx, 2] )

      if(length(unique(anot.df$hier)) <= 12){
        cols <- c("#5F4690","#1D6996","#38A6A5","#0F8554",
                  "#73AF48","#EDAD08","#E17C05","#CC503E","#94346E","#6F4070","#994E95","#666666")
      }else{
        library(RColorBrewer)
        nb.cols <- length(unique(anot.df$hier))
        cols <- colorRampPalette(c("#5F4690","#1D6996","#38A6A5","#0F8554",
                                   "#73AF48","#EDAD08","#E17C05","#CC503E","#94346E","#6F4070","#994E95","#666666"))(nb.cols)
      }
    }
  }

  plot_title <- "Lollipop Chart of Pathway Enrichment"
  plot_save <- paste0(plotName, ".png")

  # plot version 3
  # line plot, length corresponding to log p, circle at tip corresponding to enrichment ratio
  p <- ggplot(super_generic, aes(x = Pathways, y = logp, width = 0.1)) +
    geom_bar(stat = "identity", fill = "#24221e") +
    geom_point(aes(x = Pathways, y = logp, size = folds), colour = "#24221e") +
    new_scale_color() +
    geom_point(data = anot.df, aes(x, y, color = hier), shape = 15, size = 5, inherit.aes = FALSE) +
    scale_color_manual("Hierarchy", values = cols) +
    scale_size_continuous(range = c(0, 4)) +
    coord_flip() +
    labs(title = plot_title) +
    theme(plot.title = element_text(size=9)) +
    xlab("Pathway\n") + ylab("\n-log10 P-Value") + labs(size = "Enrichment Ratio") +
    theme_minimal()

  ggsave(plot_save, p, width = width, height = height)

  return(.set.mbSetObj(mbSetObj))

}
##################################################################################################

#' Function to create lollipop comparison chart
#' of pathway enrichment results from AGORA
#' versus CarveMe
#'
#' @param mbSetObj
#' @param plot.libname 'KEGG' or 'MetaCyc'.
#' @param save.libname Name of the pathway library to be used in the image name.
#' @param npath Default is set to 25.
#'
#' @export
lollipop_plot_comp <- function(mbSetObj, plot.libname, save.libname, npath = 25){

  library(data.table)
  library(tidyverse)
  library(ggplot2)

  mbSetObj <- .get.mbSetObj(mbSetObj);

  gemLib <- mbSetObj$mNet$gemLib

  pathResults <- vector("list")
  pathResults$agora_super <- read.csv("msea_ora_result_a.csv")
  pathResults$carveme_super <- read.csv("msea_ora_result_c.csv")

  path.names.all <- lapply(pathResults, rownames)
  path.intersect <- Reduce(intersect, path.names.all)

  # remove paths not found by all
  path.intersected <- lapply(pathResults, function(x) x[row.names(x) %in% path.intersect,])
  path.intersected <- lapply(path.intersected, function(x){
    if(nrow(x) > npath){
      x <- x[1:npath,]
    }
  })

  super_generic <- data.table::rbindlist(path.intersected, idcol = TRUE)

  # hits/expected = enrichment ratio
  super_generic$folds = super_generic[,5]/super_generic[,4]
  super_generic$logp = -log10(super_generic[,6])
  super_generic$logp <- scales::rescale(super_generic$logp, to = c(0, 5))

  super_generic <- super_generic %>% mutate(PlotLogP = ifelse(.id == "agora_super", -logp, logp))
  super_generic <- super_generic[order(super_generic$logp), ]
  super_generic$folds_x <- super_generic$PlotLogP/2

  super_generic$Pathways <- factor(super_generic$X, levels = unique(super_generic$X))
  super_generic <- super_generic %>% mutate(DB_ID = ifelse(.id == "agora_super", "AGORA", "CarveMe"))
  super_generic <- super_generic %>% filter(logp > 0)

  plot_title <- paste0("Enrichment of ", plot.libname, "\n(AGORA vs CarveMe)")
  plot_save <- paste0("enrichment_agora_carveme_", save.libname, ".png")

  # plot version 3
  # line plot, length corresponding to log p, circle at tip corresponding to enrichment ratio
  p <- ggplot(super_generic, aes(x = Pathways, y = PlotLogP, fill = DB_ID, width = 0.1)) +
    geom_bar(stat = "identity") +
    scale_fill_manual("GEM", values = c("AGORA" = "#ffa600",
                                        "CarveMe" = "#003f5c")) +
    geom_point(aes(x = Pathways, y = PlotLogP, size = folds, colour = DB_ID), alpha = 0.75) +
    scale_color_manual("GEM", values = c("AGORA" = "#ffa600",
                                         "CarveMe" = "#003f5c")) +
    scale_size_continuous(range = c(0, 4)) +
    scale_y_continuous(limits = c(-5, 5),
                       labels = paste0(as.character(c(seq(5, 0, -2.5), seq(2.5, 5, 2.5))))) +
    coord_flip() +
    labs(title = plot_title) +
    theme(plot.title = element_text(size=12)) +
    xlab("Pathway") + ylab("-log10 P-Value") + labs(size = "Enrichment Ratio") +
    theme_minimal()

  p

  ggsave(plot_save, p, width = 9, height = 7)

  return(.set.mbSetObj(mbSetObj))
}

#' Function to create lollipop comparison chart with asterisks
#' for significant pathways from real results.
#' @param agora_super
#' @param carveme_super
#' @param plot.libname
#' @param save.libname
#' @param asterisk_file Dataframe from csv containing enrichment results
#' from real omics data.
#' @export
# agora_super <- read.csv("~/Desktop/MicrobiomeNet/Case_Studies/Ped_IBD_CaseStudy/enrichment_results_gem_agora/metacyc_bacteria_pathways_enrichment_results/msea_ora_result.csv")
# carveme_super <- read.csv("~/Desktop/MicrobiomeNet/Case_Studies/Ped_IBD_CaseStudy/enrichment_results_gem_carveme/metacyc_bacteria_pathways_enrichment_results/msea_ora_result.csv")
# real_super <- read.csv("~/Desktop/MicrobiomeNet/Case_Studies/Ped_IBD_CaseStudy/mummichog_pathway_enrichment_metacyc_microbial.csv")
#
# lollipop_chart_real_res(agora_super, carveme_super, "MetaCyc Microbial", "metacyc_microbial_real_res", real_super, mum = T)

lollipop_chart_real_res <- function(mbSetObj, plot.libname, save.libname,
                                    asterisk_file, mum = F, width = 14, height = 7,
                                    lib.type = "kegg", anot = F, npath = 25){

  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(ggnewscale)

  mbSetObj <- .get.mbSetObj(mbSetObj);

  gemLib <- mbSetObj$mNet$gemLib

  pathResults <- vector("list")
  pathResults$agora_super <- read.csv("msea_ora_result_a.csv")
  pathResults$carveme_super <- read.csv("msea_ora_result_c.csv")
  asterisk_file <- read.csv(asterisk_file)

  path.names.all <- lapply(pathResults, rownames)
  path.intersect <- Reduce(intersect, path.names.all)

  # remove paths not found by all
  path.intersected <- lapply(pathResults, function(x) x[row.names(x) %in% path.intersect,])
  path.intersected <- lapply(path.intersected, function(x){
    if(nrow(x) > npath){
      x <- x[1:npath,]
    }
  })

  super_generic <- data.table::rbindlist(path.intersected, idcol = TRUE)
  # hits/expected = enrichment ratio
  super_generic$folds = super_generic[,5]/super_generic[,4]
  super_generic$logp = -log10(super_generic[,6])
  super_generic$logp <- scales::rescale(super_generic$logp, to = c(0, 4.75))

  super_generic <- super_generic %>% mutate(PlotLogP = ifelse(.id == "agora_super", -logp, logp))
  super_generic <- super_generic[order(super_generic$logp), ]
  super_generic$folds_x <- super_generic$PlotLogP/2

  super_generic$Pathways <- factor(super_generic$X, levels = unique(super_generic$X))
  super_generic <- super_generic %>% mutate(DB_ID = ifelse(.id == "agora_super", "AGORA", "CarveMe"))
  super_generic <- super_generic[(super_generic$logp > 0), ]

  if(mum){
    asterisk_file <- asterisk_file[asterisk_file$Gamma < 0.1, ]
  }else{
    asterisk_file <- asterisk_file[asterisk_file$FDR < 0.1, ]
  }

  super_generic$signif <- super_generic$X %in% asterisk_file$X
  super_generic$signif.coord <- ifelse(super_generic$PlotLogP < 0, super_generic$PlotLogP - 0.25, super_generic$PlotLogP + 0.25)

  plot_title <- paste0("Enrichment of ", plot.libname, "\n(AGORA vs CarveMe)")
  plot_save <- paste0("enrichment_agora_carveme_", save.libname, ".png")

  if(anot){

    if(lib.type == "kegg"){

      if(.on.public.web){
        kegg_hier <- read.csv("../../lib/plotting/kegg_brite_hierarchy_pathways.csv")
      }else{
        kegg_hier <- read.csv("~/Desktop/MicrobiomeNet/pathways/kegg_brite_hierarchy_pathways.csv")
      }

      kegg_hier$Pathway_Name <- trimws(substring(kegg_hier$Pathway_Name, 6))
      hier.inx <- match(as.character(super_generic$Pathways), kegg_hier$Pathway_Name)
      anot.df <- data.frame(x = super_generic$Pathways, y = -6, hier = kegg_hier[hier.inx, 2] )

      cols <- c("#5F4690","#1D6996","#38A6A5","#0F8554",
                "#73AF48","#EDAD08","#E17C05","#CC503E","#94346E","#6F4070","#994E95","#666666")
    }else{

      if(.on.public.web){
        meta_hier <- readRDS("../../lib/plotting/meta_cyc_paths_hier.rds")
      }else{
        meta_hier <- readRDS("~/Desktop/MicrobiomeNet/pathways/MetaCyc/2021/meta_cyc_paths_hier.rds")
      }

      meta_hier$Super.Pathways <- gsub("\\s*\\([^\\)]+\\)","",as.character(meta_hier$Super.Pathways))
      meta_hier$Super.Pathways <- gsub("//.*","",as.character(meta_hier$Super.Pathways))
      meta_hier$Super.Pathways <- gsub("superpathway of","",as.character(meta_hier$Super.Pathways))
      meta_hier$Super.Pathways[meta_hier$Super.Pathways == ""] <- "Unclassified"
      meta_hier$Super.Pathways <- tools::toTitleCase(trimws(meta_hier$Super.Pathways))

      hier.inx <- match(as.character(super_generic$Pathways), meta_hier$Pathways)
      anot.df <- data.frame(x = super_generic$Pathways, y = -6, hier = meta_hier[hier.inx, 2] )

      if(length(unique(anot.df$hier)) <= 12){
        cols <- c("#5F4690","#1D6996","#38A6A5","#0F8554",
                  "#73AF48","#EDAD08","#E17C05","#CC503E","#94346E","#6F4070","#994E95","#666666")
      }else{
        library(RColorBrewer)
        nb.cols <- length(unique(anot.df$hier))
        cols <- colorRampPalette(c("#5F4690","#1D6996","#38A6A5","#0F8554",
                                   "#73AF48","#EDAD08","#E17C05","#CC503E","#94346E","#6F4070","#994E95","#666666"))(nb.cols)
      }
    }
  }

  # plot version 3
  # line plot, length corresponding to log p, circle at tip corresponding to enrichment ratio
  p <- ggplot(super_generic, aes(x = Pathways, y = PlotLogP, fill = DB_ID, width = 0.1)) +
    geom_bar(stat = "identity") +
    scale_fill_manual("GEM", values = c("AGORA" = "#ffa600",
                                        "CarveMe" = "#003f5c")) +
    geom_point(aes(x = Pathways, y = PlotLogP, size = folds, colour = DB_ID), alpha = 0.75) +
    scale_color_manual("GEM", values = c("AGORA" = "#ffa600",
                                         "CarveMe" = "#003f5c")) +
    # try to add tile
    new_scale_color() +
    geom_point(data = anot.df, aes(x, y, color = hier), shape = 15, size = 5, inherit.aes = FALSE) +
    scale_color_manual("Hierarchy", values = cols) +
    scale_size_continuous(range = c(0, 4)) +
    scale_y_continuous(limits = c(-6, 5),
                       breaks = seq(-6, 5, by = 1),
                       labels = paste0(as.character(c("", seq(5, 0, -1), seq(1, 5, 1))))) +
    labs(title = plot_title) +
    theme(plot.title = element_text(size=12)) +
    xlab("Pathway") + ylab("-log10 P-Value") + labs(size = "Enrichment Ratio") +
    theme_minimal() +
    theme(legend.title = element_text(size=9.5),
          legend.text = element_text(size=8.5),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 9)) +
    coord_flip()

  if(any(super_generic$signif == T)){
    p <- p + geom_point(data = super_generic[super_generic$signif == T, ], aes(x = Pathways, y = signif.coord),
                        shape = "*", size = 8, show.legend = F)
  }

  ggsave(plot_save, p, width = width, height = height)

  return(.set.mbSetObj(mbSetObj))
}
