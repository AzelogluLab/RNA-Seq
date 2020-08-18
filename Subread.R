ftcounts <- function(rna_params){
  ##ftcounts
  ##Description
  ##Runs featureCounts function from Rsubread, and places the count table in a text file
  ##Arguments:
  ##rna_params: A list of parameter values  
  ##Returns:
  ##fc: A datatable containing data about the alignment results 
  library(Rsubread)
  #Change working directory
  setwd(rna_params[2])
  data_path <- "Results/AlignedFiles"
  
  #Add the filenames together
  subread_files <- list.files(path = data_path, all.file = FALSE, full.names = TRUE, recursive = FALSE)
  #Change featurecounts parameters based on the selections in the read file
  if (rna_params[9] == "N"){
    pairedEnd <- FALSE
  } 
  else if (rna_params[9] == "Y"){
    pairedEnd <-  TRUE
  }
  if(rna_params[10] == "N"){
    multMap <- FALSE
  } 
  else if (rna_params[10] == "Y"){
    multMap <- TRUE
  }
  dir.create("Results/ReadsReport")
  #Run featureCounts
  fc <- featureCounts(files = subread_files , annot.ext = rna_params[3], isGTFAnnotationFile = TRUE, GTF.featureType = "exon", GTF.attrType = "gene_id", countMultiMappingReads = multMap, allowMultiOverlap = multMap, isPairedEnd = pairedEnd, reportReads = "CORE", reportReadsPath = "Results/ReadsReport")
  return(fc)
}

GetCountsTable <- function(rna_params, group_params, fc){
  ##GetCountsTable
  ##Description
  ##Reads counts file and places it into a variable called counts, renames the rows and reorders them in alphabetical order.
  ##Arguments: 
  ##rna_params: A list containing parameters for the alignment and analysis
  ##group_params: A list containing parameters related to the sample grouping
  ##fc: A datatable containing data about the alignment results 
  ##Returns:
  ##counts : a data table of the counts, with the column names as the experiment name and the rows as the geneIDs, and the data being the counts
  #Writes the counts table to a file
  setwd(rna_params[2])
  dir.create("Results/DEG")
  write.table(x=data.frame(fc$annotation[,c("GeneID")], fc$counts, stringsAsFactors=FALSE), file = "Results/DEG/RAWCountsTable.txt", quote=FALSE, sep="\t", row.names=FALSE)
  #Reads the counts table from the file and assigns it to a variable
  counts <- read.delim("Results/DEG/RAWCountsTable.txt", check.names = FALSE, row.names=1)
  colnames(counts) <- group_params$Names
  return(counts)
}
