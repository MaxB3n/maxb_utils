library(VennDiagram)
library(ggplot2)
library(data.table)
library(this.path)
library(ggrepel)

source(file.path(this.dir(), "GenericHelpers.R"))

### SIMPLE PLOTS
PlotVennFromDT <- function(dt, valueCol, groupCol, filename){
  
  groupNames <- unique(dt[, groupCol], with = F)
  
  valueVectors <- lapply(groupNames, FUN = function(group, dt, value){dt[grepl(group, dt[,groupCol, with = F ]), c(value)]}, dt = dt, value = valueCol )
  
  venn.diagram( x = valueVectors,
                category.names = groupNames,
                output = T,
                filename = filename)
}


### QC PLOTS (Bar, pca, correlation, heatmap)
PlotQCBar <- function( msstatsTable, plotType = "count", colorsManual = NA, groupLabel = "Treatment", dataName = "", showReplicateLegend = F, face = "mono"){
  
  # Expects GROUP (conditions), RUN (unique runs), Protein or PEPTIDE from MSStats output
  # Detects protein or feature table,  plotType should be provided as "count" or "intensity" or a column name
  msstTable <- copy(msstatsTable)
  
  if ("Protein" %in% colnames(msstTable)){
    dataCol <- "Protein"
  } else { 
    setnames(msstTable, c("PEPTIDE"), c("Peptide"))
    msstTable[, LogIntensities := log2(INTENSITY)]
    dataCol <- "Peptide" 
  }
  
  msstTable <- msstTable[!is.na(LogIntensities),]
  msstTable[, Rep := as.character(as.numeric(as.character(RUN)) %% (nun(RUN)/nun(GROUP)) +1)]
  
  barTable <- msstTable[, .(.N, "SumLogIntensities" = sum(LogIntensities)), by = .(Rep, GROUP)]
  #str(barTable)
  
  switch(plotType,
         count={qcPlot <- ggplot( barTable, aes(x=Rep, y=N, fill = Rep) ) +
           ylab(paste0("Counts of Unique ", dataCol, "s")) +
           ggtitle(paste(dataName, "Unique", dataCol, "Counts per Run"))},
         intensity={qcPlot <- ggplot( barTable, aes(x=Rep, y=SumLogIntensities, fill = Rep) ) +
                              ylab("Summed Log2Intensities") +
                              ggtitle(paste(dataName, dataCol, "Intensities per Run"))},
         {if (plotType %in% colnames(msstTable)){
           qcPlot <- ggplot( msstTable, aes(x=Rep, y=plotType, fill = Rep) ) +
             ylab(plotType) +
             ggtitle(paste( dataName, dataCol, plotType, "per Run"))
         } else{ stop("argument `plotType` must be \"count\", \"intensity\", or a valid column name of the provided table!")}
           
         }
  )
  
  
  qcPlot <- qcPlot +
    geom_bar(position = "dodge", stat="identity") +
    xlab(paste0("Runs Grouped by ", groupLabel))+
    theme_classic() +
    theme(text = element_text(family = face ), plot.title = element_text(size = 18, hjust = 0.5),
          axis.title.y = element_text(margin = margin(r = 15)), axis.title.x = element_text(margin = margin(t = 15))) +
    facet_wrap(~GROUP, nrow=1, strip.position = "bottom") +
    theme(strip.placement = "outside")
  
  if (!showReplicateLegend){
    qcPlot <- qcPlot + theme(legend.position = "none", axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }
  
  if (!is.na(colorsManual[[1]])){
    qcPlot <- qcPlot + scale_fill_manual(values = colorsManual)
  }
  
  return( qcPlot )
  
}

PlotQcHeatmap <- function(mat, scale = "Log2 Intensities", showrn = F, csplit = F, title = "", 
                          legendRange = colorRamp2(quantile(mat, probs = c(0.015, 0.5, 0.985),na.rm = T), c( "blue", "white", "red")), row_names_gp = gpar(fontfamily = face), face = "sans",
                          cluster_rows = rowClusterWithNA(mat)){
  
  if (is.character(csplit) & length(csplit) == ncol(mat) ){
    hmap <- Heatmap(mat, 
                    cluster_rows = cluster_rows,
                    name = scale,
                    col =  legendRange,
                    column_title = title,
                    row_title = NULL,
                    show_row_names = showrn,
                    show_column_names = T,
                    column_gap = unit(1, "mm"),
                    cluster_columns = F,
                    column_split = csplit,                    
                    column_names_rot = 75,
                    column_names_gp = gpar(fontfamily = face), row_names_gp = row_names_gp,
                    raster_by_magick = T) 
  } else {
    hmap <- Heatmap(mat, 
                    cluster_rows = cluster_rows,
                    name = scale,
                    col =  legendRange,
                    column_title = title,
                    row_title = NULL,
                    show_row_names = showrn,
                    show_column_names = T,
                    cluster_columns = T,
                    column_names_rot = 75,
                    column_names_gp = gpar(fontfamily = face), row_names_gp = row_names_gp, column_title_gp = gpar(fontfamily = face),
                    raster_by_magick = T) 
  }
    return(hmap)
}
  

PlotVolcano <- function(table, condition, l2fcT = 1, apValT = 0.05, repel = F, face = "mono"){
    
    v <- ggplot(table[Label == condition, ], aes(x = log2FC, y = -log10(adj.pvalue), color = Effect)) + 
    geom_point( alpha = 0.35, shape = 19, size = 0.8) + 
    scale_color_manual( values = c(overrepresented = "firebrick",
                                   exclusively = "firebrick",
                                   underrepresented = "navy",
                                   absent = "navy",
                                   insignificant = "grey") ) +
    geom_vline(xintercept = -l2fcT, linetype="dashed", 
               color = "grey24", size=0.5) +
    geom_vline(xintercept = l2fcT, linetype="dashed", 
               color = "grey24", size=0.5) +
    geom_hline(yintercept = -log10(apValT), linetype="dashed", 
               color = "grey24", size=0.5) +
    ggtitle(paste("Significant", sprintf("(adj.pval < %s, abs(l2fc) > %s)", apValT, l2fcT), "Proteins in", condition)) +
    theme(legend.position = "none",
          text = element_text(family = face ), plot.title = element_text(size = 18, hjust = 0.5),
          axis.title.y = element_text(margin = margin(r = 15)), axis.title.x = element_text(margin = margin(t = 15))) #+
    
    if (repel){
    v <- v + geom_label_repel(aes( x = log2FC, y = -log10(adj.pvalue), color = Effect, label = Gene),
                              box.padding = 0.3, max.overlaps = 15)
    }
    
    return(v)
}
  
  
PlotCorrelationMatrix <- function(){
  
  Heatmap(corMat, 
          name = "Spearman \nCorrelation",
          col =  colorRamp2(c(0.89, 0.99), c( "seashell", "red")),
          cluster_columns = F, cluster_rows = F,
          column_gap = unit(1, "mm"), row_gap = unit(1, "mm"),
          column_split = rep(as.character(1:8), each = 3),  row_split = rep(as.character(1:8), each = 3),
          row_title = rep("",8), column_title = c("Protein Intensity Correlation Matrix"),
          show_row_names = T, show_column_names = T, row_names_side = "left",
          column_names_rot = 75,
          column_names_gp = gpar(fontfamily = face), row_names_gp = gpar(fontfamily = face),
          column_title_gp = gpar(fontfamily = face),
          raster_by_magick = T,
          cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
            if(i >= j) {
              grid.rect(x, y, w, h, gp = gpar(fill = col, col = col))
              grid.text(round(corMat[i, j], 2), x, y, gp = gpar(fontfamily = face, fontsize = 8)) 
            } else{
              grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "white"))
            }
          })
  
}  
  
  
  