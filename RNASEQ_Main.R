RNASEQ_Main_Sequencing <- function(){
  #Load in all of the functions
  source('ParameterSetup.R')
  source('STAR.R')
  source('QualityCheckSTAR.R')
  source('Subread.R')
  source('QualityCheckRsubread.R')
  source('EdgeR.R')
  source('GroupOrganize.R')
  source('RNAAnalysis.R')
  
  ###Setup Parameters###
  #Load in the csv files containing the RNA Sequencing and Group Parameters
  rna_params <- ParamSetup()
  group_params <- GroupParamSetup(rna_params = rna_params)
  
  ###STAR RUN###
  #Generate the Genome
  STARGenerateGenome(rna_params = rna_params)
  #unzip the files to be aligned
  gunzip_check(rna_params = rna_params)
  #Align the files to the genome
  STARAlign(rna_params = rna_params)
  
  ###Star Quality Check###
  #Quality Check functions to look at the percentage and number of reads per file
  QCSTAR_StackPerc(rna_params = rna_params, group_params = group_params)
  QCSTAR_StackNumreads(rna_params = rna_params, group_params = group_params)
  
  ###RUN Rsubread###
  #Run featurecounts function then save the counts table to a file
  featureCounts <- ftcounts(rna_params = rna_params)
  #Opens the file made by feature counts and reads it into a dataframe format
  raw_counts <- GetCountsTable(rna_params = rna_params, group_params = group_params, fc = featureCounts)
  seq_output <- list("rna_params" =  rna_params, "group_params" = group_params, "raw_counts" = raw_counts)
  return(seq_output)
}

RNASEQ_Main_Analysis <- function(seq_output){
  #Unpackage the list
  rna_params <- seq_output$rna_params
  group_params <- seq_output$group_params
  raw_counts <- seq_output$raw_counts
  #Filter out lowly expressed genes
  filter_counts <- FilterLowGenes(rna_params = rna_params, counts = raw_counts)
  #Normalize the raw counts 
  normalize_filter_counts <- NormalizeCounts(counts = filter_counts)

  
 ###Rsubread Quality Check###
  #Select groups to be compared in analysis
  group1_index <- SelectGroup1(group_params = group_params)
  group2_index <- SelectGroup2(group_params = group_params)
  
  
  #Take the selected groups and seperate them into their own dataframes
  group1_counts <- normalize_filter_counts[ ,group1_index$Index]
  group2_counts <- normalize_filter_counts[ ,group2_index$Index]
  
  #Create correlation heatmaps
  HeatmapGroup1 <- QCftcountsSample(counts = group1_counts)
  HeatmapGroup2 <- QCftcountsSample(counts = group2_counts)
  Heatmap_All_Groups <- QCftcountsSample(counts = normalize_filter_counts[,c(group1_index$Index, group2_index$Index)])
  
  #Create PCA
  makePCA(index1 = group1_index, index2 = group2_index, counts = normalize_filter_counts)
  make3D_PCA(index1 = group1_index, index2 = group2_index, counts = normalize_filter_counts)
  
  ###Run EdgeR###
  selected_counts <- raw_counts[ ,c(group1_index$Index, group2_index$Index)]
  selected_filt_counts <- FilterLowGenes(rna_params = rna_params, counts = selected_counts)
  selected_filt_norm_counts <- NormalizeCounts(counts = selected_filt_counts)
  
  dge <- Makedge(counts = selected_filt_counts, group1 = group1_index$Group, group2 = group2_index$Group, index1 = group1_index$Index, index2 = group2_index$Index) 
  exactdge <- ExactTestdge(dge = dge, group_params = group_params, group1_index = group1_index, group2_index = group2_index)
  exactResults <- exactdge$table
  exactCounts <- selected_filt_norm_counts[rownames(exactResults), ]
  exactCounts_results <- merge(exactResults, exactCounts, by = 0, all = TRUE)
  names(exactCounts_results)[1] <- "Gene"
  
  
  #DGE Heatmap
  corrected_exactCounts <- Benj_Hoch_Corr(rna_params = rna_params, counts = exactCounts_results)
  DGEHeatmap(counts = corrected_exactCounts) 
  
    #Volcano Plot
  SplitCount <- SplitDEG(csv = exactCounts_results, down_file = "main_downregulated_genes.csv", up_file = "main_upregulated_genes.csv")
  makeVolcano(corrected_exactCounts = corrected_exactCounts, exactCounts_results = exactCounts_results, split_index = SplitCount)
  makePlotMA(corrected_exactCounts = corrected_exactCounts, exactCounts_results = exactCounts_results)
  
  analysis_output <- list("group1_index" = group1_index, "group2_index" = group2_index)
}

RNASEQ_Second_Analysis <- function(seq_output, analysis_output){
  #Unpackage the lists
  rna_params <- seq_output$rna_params
  group_params <- seq_output$group_params
  raw_counts <- seq_output$raw_counts
  group1_index <- analysis_output$group1_index
  group2_index <- analysis_output$group2_index
  
  #Ask user which groups to remove
  removal_params <- SampleRemove(group_params = group_params, raw_counts = raw_counts)
  group_params <- removal_params$group_params
  raw_counts <- removal_params$raw_counts
  
  #Adjust group parameters
  adj_group1_index <- SelectGroup_secondary(group_params = group_params, group_index = group1_index)
  adj_group2_index <- SelectGroup_secondary(group_params = group_params, group_index = group2_index)
  
  #Filter out lowly expressed genes
  filter_counts <- FilterLowGenes(rna_params = rna_params, counts = raw_counts)
  #Normalize the raw counts 
  normalize_filter_counts <- NormalizeCounts(counts = filter_counts)
  
  #Take the selected groups and seperate them into their own dataframes
  group1_counts <- normalize_filter_counts[ ,adj_group1_index$Index]
  group2_counts <- normalize_filter_counts[ ,adj_group2_index$Index]
  
  #Create correlation heatmaps
  HeatmapGroup1 <- QCftcountsSample(counts = group1_counts)
  HeatmapGroup2 <- QCftcountsSample(counts = group2_counts)
  Heatmap_All_Groups <- QCftcountsSample(counts = normalize_filter_counts[,c(adj_group1_index$Index, adj_group2_index$Index)])
  
  #Create PCA
  makePCA(index1 = adj_group1_index, index2 = adj_group2_index, counts = normalize_filter_counts)
  make3D_PCA(index1 = adj_group1_index, index2 = adj_group2_index, counts = normalize_filter_counts)
  
  ###Run EdgeR###
  selected_counts <- raw_counts[ ,c(adj_group1_index$Index, adj_group2_index$Index)]
  selected_filt_counts <- FilterLowGenes(rna_params = rna_params, counts = selected_counts)
  selected_filt_norm_counts <- NormalizeCounts(counts = selected_filt_counts)
  
  dge <- Makedge(counts = selected_filt_counts, group1 = adj_group1_index$Group, group2 = adj_group2_index$Group, index1 = adj_group1_index$Index, index2 = adj_group2_index$Index) 
  exactdge <- ExactTestdge(dge = dge, group_params = group_params)
  exactResults <- exactdge$table
  exactCounts <- selected_filt_norm_counts[rownames(exactResults), ]
  exactCounts_results <- merge(exactResults, exactCounts, by = 0, all = TRUE)
  names(exactCounts_results)[1] <- "Gene"
  
  
  #DGE Heatmap
  corrected_exactCounts <- Benj_Hoch_Corr(rna_params = rna_params, counts = exactCounts_results)
  DGEHeatmap(counts = corrected_exactCounts) 
  
  #Volcano Plot
  SplitCount <- SplitDEG(csv = exactCounts_results, down_file = "second_downregulated_genes.csv", up_file = "second_upregulated_genes.csv")
  makeVolcano(corrected_exactCounts = corrected_exactCounts, exactCounts_results = exactCounts_results, split_index = SplitCount)
  makePlotMA(corrected_exactCounts = corrected_exactCounts, exactCounts_results = exactCounts_results)
}


###Run Code###
seq_output <- RNASEQ_Main_Sequencing()
analysis_output <- RNASEQ_Main_Analysis(seq_output = seq_output)
RNASEQ_Second_Analysis(seq_output = seq_output, analysis_output = analysis_output)  
