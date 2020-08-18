QCSTAR_StackPerc <- function(rna_params, group_params){
  ##QCSTAR_StackPerc
  ##Description
  ##Takes the data from the log files produced by STAR and creates a stacked bar graph containing the percent reads for each of the samples 
  ##Arguments:
  ##rna_params: A list containing parameters for the alignment and analysis
  ##group_params: A list containing parameters related to the sample grouping
  ##Returns:
  ##None
  setwd(rna_params[2])
  #Get filename from params list
  data_path <- "Results/STAR_metadata_finallog"
  filenames <- list.files(path = data_path, pattern = ".txt", full.names = TRUE, recursive = FALSE)
  #Make blank lists
  unique_mapped_num <- c()
  mult_loci_num <- c()
  too_many_loci_num <- c()
  unmapped_mismatch_num <- c()
  unmapped_tooshort_num <- c()
  unmapped_other_num <- c()
  chimeric_num <- c()
  #Read through files and extract the values to be use in the bar plot 
  for(file in filenames){
    logfile <- readLines(file)
    split_logfile <- strsplit(logfile, "[[:space:]]")
    unique_mapped <- split_logfile[[10]]
    unique_mapped_num <- c(unique_mapped_num, as.numeric(sub("%", "", unique_mapped[30])))
    mult_loci <- split_logfile[[25]]
    mult_loci_num <- c(mult_loci_num, as.numeric(sub("%", "", mult_loci[22])))
    too_many_loci <- split_logfile[[27]]
    too_many_loci_num <- c(too_many_loci_num, as.numeric(sub("%", "", too_many_loci[23])))
    unmapped_mismatch <- split_logfile[[29]]
    unmapped_mismatch_num <- c(unmapped_mismatch_num, as.numeric(sub("%", "", unmapped_mismatch[16])))
    unmapped_tooshort <- split_logfile[[30]]
    unmapped_tooshort_num <- c(unmapped_tooshort_num, as.numeric(sub("%", "", unmapped_tooshort[25])))
    unmapped_other <- split_logfile[[31]]
    unmapped_other_num <- c(unmapped_other_num, as.numeric(sub("%", "", unmapped_other[28])))
    chimeric <- split_logfile[[34]]
    chimeric_num <- c(chimeric_num, as.numeric(sub("%", "", chimeric[34])))
  }
  
  
  #Organize values to graph format
  SampleNames <- group_params$Names
  Legend <- c(rep("unique-mapped", length(SampleNames)), rep("multiple-loci", length(SampleNames)), rep("too-many-loci", length(SampleNames)), rep("unmapped-mismatch", length(SampleNames)), rep("unmapped-tooshort", length(SampleNames)), rep("unmapped-other", length(SampleNames)), rep("chimeric", length(SampleNames)))
  condition <- rep(c(SampleNames), 7)
  percent_reads <- c(unique_mapped_num, mult_loci_num, too_many_loci_num, unmapped_mismatch_num, unmapped_tooshort_num, unmapped_other_num, chimeric_num)
  bar_label <- c(unique_mapped_num, rep(NA, 6 * length(SampleNames))) 
  data <- data.frame(Legend, condition, percent_reads)
  
  library(ggplot2)
  library(reshape2)
  #Make graph
  ggplot(data, aes(fill=Legend, y=percent_reads, x=condition)) + 
    geom_bar(stat="identity") + 
    geom_text(aes(label = bar_label), vjust = 1.6, size = 3.5, position = "stack") +
    xlab("Samples") + ylab("Percent Reads (%)") +
    ggtitle("QC STAR (% of Reads)") 
}

QCSTAR_StackNumreads <- function(rna_params, group_params){
  ##QCSTAR_StackNumreads
  ##Description
  ##Takes the data from the log files produced by STAR and creates a stacked bar graph containing the number reads for each of the samples
  ##Arguments:
  ##rna_params: A list containing parameters for the alignment and analysis
  ##group_params: A list containing parameters related to the sample grouping
  ##Returns:
  ##None
  setwd(rna_params[2])
  data_path <- "Results/STAR_metadata_finallog"
  filenames <- list.files(path = data_path, pattern = ".txt", full.names = TRUE, recursive = FALSE)
  #Make blank lists
  #total_num <- c()
  unique_num <- c()
  mult_loci_num <- c()
  too_many_num <- c()
  chim_num <- c()
  #Read through files and extract the values for the bar pl ot
  for(file in filenames){
    logfile <- readLines(file)
    split_logfile <- strsplit(logfile, "[[:space:]]")
    unique <- split_logfile[[9]]
    unique_num <- c(unique_num, unique[25])
    mult_loci <- split_logfile[[24]]
    mult_loci_num <- c(mult_loci_num, mult_loci[17]) 
    too_many <- split_logfile[[26]]
    too_many_num <- c(too_many_num, too_many[18])
    chim <- split_logfile[[33]]
    chim_num <- c(chim_num, chim[[29]])
  }
  #Organize values into graph format
  SampleNames <- group_params$Names
  Legend <- c(rep("unique-mapped", length(SampleNames)), rep("multiple-loci", length(SampleNames)), rep("too many loci", length(SampleNames)), rep("chimeric", length(SampleNames)))
  condition <- rep(c(SampleNames), 4)
  num_reads <- as.numeric(c(unique_num, mult_loci_num, too_many_num, chim_num))
  bar_label <- c(unique_num, rep(NA, 3 * length(SampleNames)))
  data <- data.frame(Legend, condition, num_reads)
  
  
  library(ggplot2)
  library(reshape2)
  #Make Graph
  ggplot(data, aes(fill=Legend, y=num_reads, x=condition)) + 
    geom_bar(stat="identity") + 
    geom_text(aes(label = bar_label), vjust = 1.6, size = 3.5, position = "stack") +
    xlab("Samples") + ylab("Number of Reads") + 
    ggtitle("QC STAR (# of Reads)") 
}
    
