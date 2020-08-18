Makedge <- function(counts, group1, group2, index1, index2){
  ##Make Dge
  ##Description
  ##Makes A dge object 
  ##Arguments:
  ##counts: A dataframe, with the rows as genes and the columns as samples
  ##group1: The first selected group to be put through analysis, this group will serve as the baseline. 
  ##group2: The second selected group to be put through analysis. 
  ##index1: The index of the group1 group names in the group parameters
  ##index2: The index of the group2 group names in the group parameters
  ##Returns
  ##dge: A dge object created from a counts table to be used with the edgeR library
  library(edgeR)
  dge <- DGEList(counts = counts, group = c(rep(group1, length(index1)), rep(group2, length(index2))), genes = rownames(counts))
  return(dge)
}


ExactTestdge <- function(dge, group_params, group1_index, group2_index){
  ##ExactTestdge
  ##Description
  ##Takes the dge object and does an exact test on the data
  ##Arguments
  ##dge: A dge object created from a counts table to be used with the edgeR library
  ##group_params: A list, containing group parameter information
  ##group1_index: The index number of group1 being analyzed in the group_params
  ##group2_index: The index number of group2 being analyzed in the group_params
  ##Returns
  ##exactdge: A dge object that contains the results from the exact test
  library(edgeR)
  #create the design format for the data to be used when running the exact test function
  design <- model.matrix(~ 0 + factor(droplevels(group_params$GroupParameters$Groups[c(group1_index$Index, group2_index$Index)])))
  #Assign the group information to the design
  colnames(design) <- c(substr(colnames(design)[1], start = 99, stop = nchar(colnames(design)[1])), substr(colnames(design)[2], start = 99, stop = nchar(colnames(design)[2])))
  rownames(design) <- group_params$Names[c(group1_index$Index, group2_index$Index)] 
  dge$samples$group <- droplevels(group_params$GroupParameters$Groups[c(group1_index$Index, group2_index$Index)])
  #Make the dge and run the exact test on the data 
  dge <- estimateDisp(dge, design)
  exactdge <- exactTest(dge, pair = c(group1_index$Group, group2_index$Group) ,  dispersion = "common" )
  return(exactdge)
} 


NormalizeCounts <- function(counts){
  ##NormalizeCounts
  ##Description
  ##Takes a counts table and normalizes the data
  ##Arguments:
  ##counts: a counts table with the rows as genes and the columns as samples
  ##NormFactors: Normalization factors calculated from edgeR
  ##Returns: 
  ##counts: a counts table with the rows as genes and the columns as samples (normalized)
  ##Assigns a normalized counts table to the environmentl
  library(edgeR)
  dge <- DGEList(counts = counts, group = colnames(counts), genes = rownames(counts))
  dge <- calcNormFactors(dge, method = 'TMM')
  NormFactors <- dge$samples$norm.factors
  assign('NormFactors', NormFactors, envir = .GlobalEnv)
  for (row in 1:length(rownames(counts))){
    for (col in 1:length(colnames(counts))){
      counts[row, col] <- counts[row, col]/NormFactors[col]
    }
  }
  return(counts)
}