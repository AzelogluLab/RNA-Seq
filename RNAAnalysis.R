install.packages("rgl")


Benj_Hoch_Corr <- function(rna_params, counts){
  ##FindDiffGenes
  ##Description
  ##Takes the data produced from the exact test done on the counts data, then filters out the genes through the Benjamini Hochberg Correction
  ##Arguments: 
  ##rna_params:
  ##counts: 
  ##Returns: 
  ##counts: A counts table containing the genes that are filtered out by the Benjamini Hochberg Correction
  
  #Order the table by p - value
  counts <- counts[order(counts$PValue),]
  #Add a column with rank in correspondence with the p - values 
  counts$Rank <- 1:length(counts$PValue)
  #Add a column for the Benjamini Hochberg Correction Value
  counts$Benj_Hoch_Corr <- NA
  #Calculate Benjamini Hochberg Correction
  for (row in counts){
    counts$Benj_Hoch_Corr <- ((counts$Rank)/length(counts$Rank)) * as.numeric(rna_params[11])
  }
  
  stats <- counts[, c("logFC", "logCPM", "PValue", "Rank", "Benj_Hoch_Corr")]
  for(i in 1:length(rownames(stats))){
    if(stats[i,3] == stats[i,5]){
      cutoff <- i
      break
    }
    else if (stats[i,3] > stats[i,5]){
      cutoff <- i
      break
    }
  }
  cutoff <- cutoff - 1
  i <- length(colnames(counts))
  i2 <- i - 2
  CutoffCounts <- counts[, c(1, 2, 3, i - 1, i, 4:i2)]
  CutoffCounts <- CutoffCounts[1:cutoff, ]
  
  #Make counts table with the differentially expressed genes and remove the non - significant genes from the exactdge object 
  return(CutoffCounts)
}


DGEHeatmap <- function(counts){
  ##DGEHeatmap
  ##Description
  ##Takes the counts table that has been through the Benjamini Hochberg correction and the filtering of lowly expressed genes
  ##then it makes a heat map where the genes are put through a standard deviation accross the samples. 
  ##Arguments:
  ##counts: A counts table that contains the filtered genes 
  ##Returns: 
  ##None
  #Take the standard deviation of the rows
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  #read in the data
  counts <- counts[order(counts$logFC),]
  data <- counts[, 7:length(colnames(counts))]
  data <- t(apply(data, 1, function(x)(x-min(x))/(max(x)-min(x))))
  #genes <- rep(rownames(data), length(colnames(data)))
  col_fun <- colorRamp2(c(0, .5, 1), c("red", "white", "blue"))
  row_ha <- rowAnnotation(logFC = anno_barplot(counts$logFC))
  Heatmap(data, name = "Counts", col = col_fun, show_row_dend = FALSE, 
          clustering_distance_columns = "pearson", show_row_names = FALSE,
          right_annotation = row_ha, cluster_rows = FALSE)
  
}

CustomDGEHeatmap <- function(counts){
  ##CustomDGEHeatmap
  ##Description
  ##Takes a counts table and makes a heatmap
  ##Arguments:
  ##counts: A counts table containing genes and their corresponding count values
  ##Return:
  ##None
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  #read in the data
  counts <- counts[order(counts$logFC),]
  rows = counts$Gene
  data <- counts[, 5:length(colnames(counts))]
  data <- t(apply(data, 1, function(x)(x-min(x))/(max(x)-min(x))))
  rownames(data) = rows
  #genes <- rep(rownames(data), length(colnames(data)))
  col_fun <- colorRamp2(c(0, .5, 1), c("red", "white", "blue"))
  row_ha <- rowAnnotation(logFC = anno_barplot(counts$logFC))
  Heatmap(data, name = "Counts", col = col_fun, show_row_dend = FALSE, 
          clustering_distance_columns = "pearson", show_row_names = TRUE,
          right_annotation = row_ha, cluster_rows = FALSE, row_names_side = 'left')
  
}

GeneSearch <- function(counts, file){
  ##GeneSearch
  ##Description
  ##Reads a file and filters out a counts table based on the genes in the file
  ##Arguments
  ##counts: A counts table containing genes and their corresponding count values
  ##file: A string, a filename that hs gene ids inside
  ##Returns:
  ##counts: A new counts table with the genes specified in the file
  data <- read.delim(file = file)
  Index <- c()
  for (i in 1:length(data$Gene.Symbol)){
    Index <- c(Index, which(counts$Gene == data$Gene.Symbol[i]))
  }
  counts <- counts[Index, ]
  return(counts)
}

GeneBarGraph <- function(counts, gene){
  ##GeneBarGraph
  ##Description
  ##Searches through the counts table for the gene of interest 
  ##Arguments
  ##counts: A counts table containing genes and their corresponding counts values
  ##gene: A gene of interest
  ##Returns:
  ##None
  library(ggplot2)
  data <- counts[which(rownames(counts) == 'DDN'), ]
  group <- colnames(data)
  counts <- c()
  for (i in 1:length(colnames(data))){
    counts <- c(counts, data[1, i])
  }
  
  data <- data.frame("Group" = group, "Counts" = counts)
  print(data)
  ggplot(data, aes(x = Group, y = Counts)) + 
    geom_bar(fill = "#0073C2FF", stat = "identity") + 
    geom_text(aes(label = Counts), vjust = -0.3)
  
}
GeneBarGraph(counts = raw_counts, gene = 'DDN')
data <- GeneSearch(counts = exactCounts_results, file = "geneid_symbol_cytokines.txt")
CustomDGEHeatmap(counts = GeneSearch(counts = selected_filt_counts, file =  "geneid_symbol_cytokines.txt"))

FilterLowGenes <- function(rna_params, counts){
  ##FilterDiffGene
  ##Description
  ##Takes the counts table that only contain the differentially expressed genes and filter them further by removing lowly expressed genes
  ##"lowly expressed genes" are defined by the user in the params file
  ##Arguments: 
  ##params: A list of parameter values
  ##DGEcounts: A counts table that only contains the significantly differentiated genes 
  ### FILTER OUT LOW EXPRESSED GENES
  Filter <- c()
  for(row in 1:length(rownames(counts))){
    AddVector <- c()
    for(col in 1:length(colnames(counts))){
      AddVector <- c(AddVector, counts[row,col])
      if(col == length(colnames(counts))){
        median <- median(AddVector)
        if(median < as.numeric(rna_params[12])){
          Filter <- c(Filter, as.numeric(row))
        }
      }
    }
  }
  counts <- counts[-Filter, ]
  return(counts)
}

makeVolcano <- function(corrected_exactCounts, exactCounts_results, split_index){
  library(data.table)
  corrected_exactCounts <- corrected_exactCounts[order(corrected_exactCounts$logFC),]
  corrected_exactCounts <- setDT(corrected_exactCounts, keep.rownames = TRUE)[]
  exactCounts_results <- setDT(exactCounts_results, keep.rownames = TRUE)[]
  exactCounts_results <- exactCounts_results[order(exactCounts_results$logFC),]
  with(exactCounts_results, plot(logFC, -log10(PValue), pch = 20, main= "Volcano Plot", xlim = c(-4.5, 4.5), xlab = expression('log'[2]*Foldchange), ylab = expression('log'[10]*"(p-value)")))
  with(subset(exactCounts_results, exactCounts_results$Gene %in% exactCounts_results$Gene[1:split_index-1]), points(logFC, -log10(PValue), pch=20, col="red"))
  with(subset(exactCounts_results, exactCounts_results$Gene %in% exactCounts_results$Gene[split_index:length(exactCounts_results$Gene)]), points(logFC, -log10(PValue), pch=20, col="blue"))
  with(subset(exactCounts_results, -log10(PValue) < 5), points(logFC, -log10(PValue), pch = 20, col = "black"))
  legend("topright", c("Non - DE", "Up", "Down"), fill = c("black", "blue", "red"))
  #with(AllCountsList[c(insert_index), ], points(logFC, -log10(PValue), pch=20, col="orange"))
  #library(calibrate)
  #with(AllCountsList[c(insert_index), ], textxy(logFC, -log10(PValue), labs = Gene, cex = .8))
}

makePlotMA <- function(corrected_exactCounts, exactCounts_results){
  library(data.table)
  corrected_exactCounts <- corrected_exactCounts[order(corrected_exactCounts$logFC),]
  corrected_exactCounts <- setDT(corrected_exactCounts, keep.rownames = TRUE)[]
  exactCounts_results <- exactCounts_results[order(exactCounts_results$PValue), ]
  exactCounts_results <- setDT(exactCounts_results, keep.rownames = TRUE)[]
  with(exactCounts_results, plot(logCPM, logFC, pch = 20, main= paste("Plot MA (", length(corrected_exactCounts$PValue), " DE Genes)", sep = ""), xlab = "log CPM", ylab = "log fold-change"))
  with(subset(exactCounts_results, PValue < exactCounts_results$PValue[length(corrected_exactCounts$PValue) + 1] & logFC > 0), points(logCPM, logFC, pch=20, col="blue"))
  with(subset(exactCounts_results, PValue < exactCounts_results$PValue[length(corrected_exactCounts$PValue) + 1] & logFC < 0), points(logCPM, logFC, pch=20, col="red"))
  legend("topright", c("Non - DE", "Up", "Down"), fill = c("black", "blue", "red"))
}

SplitDEG <- function(csv, down_file, up_file){
  csv <- csv[order(csv$logFC),]
  count <- 0
  for(i in csv$logFC){
    count = count + 1
    if(i > 0){
      break
    }
  }
  write.csv(csv[1:count-1, ], file = down_file)
  write.csv(csv[count:length(rownames(csv)), ], file = up_file)
  return(count)
}


makePCA <- function(index1, index2, counts){
  library(ggfortify)
  counts <- counts[,c(index1$Index, index2$Index)]
  pcaList <- as.data.frame(counts)
  pcaList <- as.data.frame(t(pcaList))
  pcaList$Group <- c(rep(index1$Group, length(index1$Index)), rep(index2$Group, length(index2$Index)))
  autoplot(prcomp(pcaList[1:length(pcaList) - 1]), data = pcaList, colour = 'Group', label = TRUE)
}

makePCA_all <- function(counts, group_params){
  pcaList <- as.data.frame(counts)
  pcaList <- as.data.frame(t(pcaList))
  pcaList$Group <- group_params$GroupParameters$Groups
  autoplot(prcomp(pcaList[1:length(pcaList) - 1]), data = pcaList, colour = 'Group', label = TRUE)
}

make3D_PCA <- function(index1, index2, counts){
  library(rgl)
  library(pca3d)
  counts <- counts[,c(index1$Index, index2$Index)]
  pcaList <- as.data.frame(counts)
  pcaList <- as.data.frame(t(pcaList))
  pcaList$Group <- c(rep(index1$Group, length(index1$Index)), rep(index2$Group, length(index2$Index)))
  pca <- prcomp(pcaList[, 1:length(pcaList) - 1], scale. = TRUE)
  pca_2 <-princomp(pcaList[,1:length(pcaList)-1], cor=TRUE, scores=TRUE)
  plot3d(pca_2$scores[,1:length(pcaList)-2],col=pcaList$Group)
  gr <- factor(pcaList$Group)
  pca3d(pca, group = gr, show.labels = TRUE)
}
