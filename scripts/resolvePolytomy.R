library(shiny)
library(plotly)
library(optparse)
library(dplyr)
library(stringr)
library(ape)
library(testit)
library(phangorn)
library(phytools)
library(castor)
library(optparse)

source('./LingPhyloR/inferenceUtils.R')

option_list = list(
  make_option(c("-i", "--input-data"), type="character", default=NULL, 
              help="dataset file path", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

resolved_trees<-function(nexus_path) {
  trees = read.nexus(nexus_path, tree.names = NULL, force.multi = TRUE)
  # Since some of these trees are not binary, resolve them.
  all_resolved_trees = list()
  for (j in 1:length(trees)) {
    # print(c(j, length(trees), length(all_resolved_trees)), file=stderr())
    tree_j = trees[[j]]
    resolved = resolveAllNodesForUnrooted(tree_j)
    for (i in 1:length(resolved)) {
      has_found=FALSE
      tree_i = resolved[[i]]
      if (length(all_resolved_trees)>0) for (k in 1:length(all_resolved_trees)) {
        tree_k = all_resolved_trees[[k]]
        if (all.equal(tree_i, tree_k, use.edge.length =FALSE)) {
          has_found=TRUE
          break
        }
      }
      if (!has_found) all_resolved_trees[[length(all_resolved_trees)+1]] = tree_i
    }
  }
  
  class(all_resolved_trees) = 'multiPhylo'
  trees = all_resolved_trees
  for(t in trees) {
    # print(class(t))
    # print(t$tip.label[[1]])
    t_rooted = ape::root(t, outgroup=c(1), resolve.root=TRUE) # root at first leaf (arbitrary)
    write.tree(t_rooted, file=stdout())
  }
  # trees
}

resolved_trees(opt[['input-data']])

