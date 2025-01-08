
###
# The following functions are all part of the typical Spectronaut/FragPipe/MaxQuant => MSStats workflow
###

# These break down the Protein Identifier used by most peptide/protein summarization tools
SimplifyUniprot <- function(lst){
  return( unlist(tstrsplit(lst, "\\|", keep=2)) )
}

SimplifyCommon <- function(lst){
  return( unlist(tstrsplit(lst, "\\||_", keep=3)) )
}

# This takes group comparison data and protein level data to mark well-evidenced proteins in the group comparison data
# Specifically, whether all runs/replicates of a condition have a protein present
LabelComplete <- function(results, intens){
  
  infComplete <- intens[!is.na(LogIntensities), .N, by = .(GROUP, Protein)][, .(Protein, N, Complete = N == max(N)), by = GROUP]
  results[infComplete, PositiveComplete := i.Complete, on = .(PositiveCondition = GROUP, Protein = Protein)]
  results[infComplete, NegativeComplete := i.Complete, on = .(NegativeCondition = GROUP, Protein = Protein)]
  
  return(results)
}

# This marks Significant results passing both an adjusted pvalue and fold change threshold
# If the table has been modified by the above, infinites will be filtered appropriately
LabelSigEffect <- function( results, apThreshold = 0.05, l2fcThreshold = 1 ) {
  
  results[, Significant := F][adj.pvalue < apThreshold & abs(log2FC) > l2fcThreshold, Significant := T]
  infBool <- "PositiveComplete" %in% colnames(results)
  if (infBool){
    results[log2FC == Inf & PositiveComplete == F, Significant := F] [log2FC == -Inf & NegativeComplete == F, Significant := F]
  }
  
  results[, Effect := "insig"]
  results[is.na(log2FC), Effect := "missing"]
  results[sign(log2FC) == 1 & Significant == T, Effect := "overrepresented"]
  results[sign(log2FC) == -1 & Significant == T, Effect := "underrepresented"]
  
  if (infBool){
    results[log2FC == Inf & PositiveComplete == T, Effect := "exclusively"]
    results[log2FC == -Inf & NegativeComplete == T, Effect := "absent"]
  }
  
  return(results)
}

# hierarchical clustering function for use with Complex Heatmap 
rowClusterWithNA <- function(mat, corr = FALSE, na.value = 0, ...){
  # corr bool asserts whether to use euclidean distance or correlation 
  # euclidean distance is more consistent for similar magnitudes, but can overlook patterns in proteins with very different absolute magnitudes but similar kinetics
  if (corr){
    mat[is.na(mat)] <- na.value
    dst <- as.dist(1-cor(t(mat), use = "pairwise"))
    dst[is.na(dst)|!is.finite(dst)] <- na.value
    hclust(dst, method = "ward.D2", ...)
  } else {
    mat[is.na(mat)] <- na.value
    dst<-dist(mat)
    dst[is.na(dst)|!is.finite(dst)] <- na.value
    hclust(dst, method = "ward.D2", ...)
  }
}

# Following functions are especially useful for network analyses, ie network propagation or shortest paths

# Does what it says on the box--hardcoded filepath
LoadUniProtMap <- function( uniprotLst, idName="STRING", species = "HUMAN"){
  switch(toupper(species),
         HUMAN = {mapFile <- "~/Documents/Databases/HUMAN_9606_idmapping.dat.gz"}
  )
  
  map <- fread(mapFile, header = F)[V2 == idName, .(V1,V3)]
  setnames(map, c("V1","V3"), c("Uniprot", idName))
  return(map)
  #return( unlist(lapply( uniprotLst, function(uniProt){ map[V1 == uniProt, V3] } )) )
}

#
MakeIgraphFromEdges <- function(edges){
  net.graph <- igraph::graph_from_edgelist(as.matrix(edges ), directed= FALSE)
  components <- igraph::decompose(net.graph)
  componentSizes <- sapply (components, FUN = function(x)length(igraph::V(x)))
  
  dev.new()
  barplot(componentSizes)
  
  net.graph <- components[[which.max(componentSizes)]]
  print("Number of Nodes: " )
  print( as.character(length(igraph::V(net.graph))) )
  print( "Number of Edges: " )
  print( as.character(length(igraph::E(net.graph))) )
  return(net.graph)
}


GmtToList <- function(gmtDT, group = "ont", identifier = "gene"){
  
  groupNames <- unlist(unique(gmtDT[, group, with = F]))
  #
  groupLists <- pbapply::pblapply(groupNames, FUN = function(groupName, dt, groupCol, geneCol){
                        list(dt[ (dt[[groupCol]] %in% groupName), ][[geneCol]]) 
                       }, dt = gmtDT, groupCol = group, geneCol = identifier, cl = makeCluster(6))
  
  return( groupLists )
}

