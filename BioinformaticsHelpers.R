
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

# Forces a .gmt format to lists that are easier to handle in R
GmtToList <- function(gmtDT, group = "ont", identifier = "gene"){
  
  groupNames <- unlist(unique(gmtDT[, group, with = F]))
  #
  groupLists <- pbapply::pblapply(groupNames, FUN = function(groupName, dt, groupCol, geneCol){
                        list(dt[ (dt[[groupCol]] %in% groupName), ][[geneCol]]) 
                       }, dt = gmtDT, groupCol = group, geneCol = identifier, cl = makeCluster(6))
  
  return( groupLists )
}


### evidence.txt file handling functions using code from Yuan Zhou's "Fragpipe_to_SAINTexpress.R" file 
# Spectronaut/FragPipe peptide-wise "msstats.csv" => MaxQuant/artMS "evidence.txt"
formatEvidenceFile <- function(msstatsfile, bputilsPath = "~/Downloads/Code/bp_utils", outFilePath = "./"){
  
  if (is.character(msstatsFile)){
    peptides <- fread(msstatsFile)
  } else if(is.data.table(msstatsfile)){
    peptides <- msstatsfile
  } else if(is.data.frame(msstatsfile)){
    peptides <- data.table(msstatsfile)
  } else{
    stop("<msstatsfile> argument expects a valid file path to peptide-wise msstats.csv input file OR a data.table containing the same data.")
  }
  
  source (file.path(bputilsPath, "spectronautFile2ArtMS.R"))
  
  cf<- list()
  # normalization method FALSE = no normalization; default is global medians which you can get my removing/commenting out all normalization lines
  # cf$msstats$normalization_method = FALSE
  #cf$msstats$normalization_method = "globalStandards"
  #cf$msstats$normalization_reference <-  "P38398"
  # should artms attempt to annotate proteins 1 = yes; 0 = no
  cf$output_extras$annotate$enabled <- as.integer(1)
  # should artms do QC 1 = yes; 0= no
  cf$qc$extended <- as.integer(0)
  cf$qc$basic <- as.integer(0)
  # cf$output_extras$annotate$species <- "MOUSE"
  # make files in artMS format
  spectronautFile2ArtMS(spectronautPeptideFile, 
                        outFilePrefix = outFilePath, 
                        artmsConfig = cf, contrastPatterns  = contrastPatterns)
}

# Remove contaminants and subset evidence.txt file 
cleanupEvidenceFile <- function(evidence = NULL,
                                newEvidenceFileName = "evidence_sub.txt", save = T,
                                filePath = "./", 
                                contaminants = c("O77727", "P00698", "P00761", "P00883", "P02662", "P02663", "P02666", "P02668", "P02769")){
  if (!is.data.table(evidence)){
    print(paste("No <evidence> argument provided, cleanupEvidenceFile is reading evidence.txt from the directory:",filePath))
    evidence <- read.table(file = file.path(filePath,"evidence.txt"), header = T, sep = "\t", stringsAsFactors = F, check.names = F)
  }
  
  evidence_sub <- evidence[-which(is.na(evidence$Intensity)), ]
  # check Leading proteins format
  if(any(grepl("sp\\|", evidence_sub$`Leading proteins`)))
  {
    #evidence_sub$`Leading proteins` <- gsub("sp\\|", "", evidence_sub$`Leading proteins`)
    # evidence_sub$`Leading proteins`[grep("_HUMAN", evidence_sub$`Leading proteins`, invert = T)] <- paste("CON__", evidence_sub$`Leading proteins`[grep("_HUMAN", evidence_sub$`Leading proteins`, invert = T)], sep = "")
    #evidence_sub$`Leading proteins` <- gsub("\\|.*", "", evidence_sub$`Leading proteins`)
    evidence_sub$`Leading proteins` <- gsub('([a-z,0-9,A-Z,\\_]+\\|{1})([A-Z,0-9,\\_]+)(\\|[A-Z,a-z,0-9,_]+)',
                                            '\\2', evidence_sub$`Leading proteins`)
  }
  if(any(contaminates %in% evidence_sub$`Leading proteins`))
  {
    evidence_sub$`Leading proteins`[which(evidence_sub$`Leading proteins` %in% contaminates)] <- 
      paste("CON__", evidence_sub$`Leading proteins`[which(evidence_sub$`Leading proteins` %in% contaminates)], sep = "")
  }
  
  if (save & is.character(newEvidenceFileName))  write.table(evidence_sub, file.path(filePath,newEvidenceFileName), sep = "\t", row.names = F, quote = F)
  return(evidence_sub)
}

# spectralCounts should be data.table from reprint.spc.tsv
formatSpcTable <- function(spectralCounts){
  spcTable <- melt(spectralCounts[2:nrow(spectralCounts),], id.vars = c("PROTID","GENEID","PROTLEN"))
  setnames(spcTable, c("PROTID","GENEID","PROTLEN","variable","value"), c("Prey","PreyGene","PreyLength","RunName","SpCount"))
  spcTable[, c("Bait", "Batch", "Treatment", "Replicate") := tstrsplit(RunName, split="_")[1:4] ]
  return(spcTable)
}

# Run SAINTexpressR
runSaintExpressR_onSpcTable <- function(spcTable = NULL, pseudoControls = NULL){
    if (is.null(spcTable) ){
      spcTable <- formatSpcTable( fread("reprint.spc.tsv") ) }
  print(str(spcTable))
  preysTable <- unique( spcTable[, .(Prey, PreyLength, PreyGene)] )
  setnames(preysTable, c("Prey", "PreyLength", "PreyGene"), c("prey", "preyLength", "preyGene" ))
  baitsTable <- unique( spcTable[, .(RunName, Bait, Treatment)][, Treatment := "T"] )
  setnames(baitsTable, c("RunName", "Bait", "Treatment"), c("run", "bait", "treatment"))
  interactionsTable <- unique(spcTable[, .(RunName, Bait, Prey, as.numeric(SpCount))])
  setnames(interactionsTable, c("RunName", "Bait", "Prey", "V4"), c("run", "bait", "prey", "spc"))
  if (is.null(pseudoControls)){
    baitsTable[bait == "Control", treatment := "C"] #[, bait := tstrsplit(run, split = "_", keep = 2)[[1]]]
    saintOut <- SaintExpressR.SPC(interactionsTable, baitsTable, preysTable, fixedBeta1 = -4.6)
    #experimentalControlsTable <- spcTable[, .(RunName, Prey, Batch)] [spcTable[Bait == "Control", .(mean(as.numeric(SpCount)), var(as.numeric(SpCount))), by = .(Prey, Batch)], c("mean","omega") := .( i.V1, 1-sqrt(i.V1/i.V2)), on = .(Batch = Batch, Prey = Prey)] [!is.finite(omega), omega := 1.1] [, Batch := NULL]
    #setnames(experimentalControlsTable, c("run", "prey", "mean", "omega"))
    # Use Saint regularly with experimental controls per batch
    #saintReg <- SaintExpressR.SPC(interactionsTable, baitsTable, preysTable, experimentalControlsTable, fixedBeta1 = -4.6)
  } else {
    stop("No code for pseudocontrols yet, see B2AI_NMMFanalysis.Rmd")
  }
  return(saintOut)
}


# Create two file folders with int and spc files ready to run with saint

