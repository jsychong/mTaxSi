#' Function to subset phylogenetic trees
#' @export
subset_tree <- function(tree.type, microbes){

  library(treeio)
  library(ape)

  if(tree.type == "agora"){
    tree_phylo <- ape::read.tree("~/Desktop/MicrobiomeNet/MicrobiomeNet_libs/trees/agora_phylo_tree.nwk")
    tree_rxn <- ape::read.tree("~/Desktop/MicrobiomeNet/MicrobiomeNet_libs/trees/agora_rxn_tree.nwk")
  }else{
    tree_phylo <- ape::read.tree("~/Desktop/MicrobiomeNet/MicrobiomeNet_libs/trees/carveme_phylo_tree.nwk")
    tree_rxn <- ape::read.tree("~/Desktop/MicrobiomeNet/MicrobiomeNet_libs/trees/carveme_rxn_tree.nwk")
  }

  phylo_tree <- keep.tip(tree_phylo, microbes)
  write.tree(phylo_tree, "subset_phylo_tree.nwk")

  rxn_tree <- keep.tip(tree_rxn, microbes)
  write.tree(rxn_tree, "subset_rxn_tree.nwk")

}
