ParamSetup <- function(){
  ##Parameter Setup
  ##Description:
  ##Reads in a csv file containing paramter values for the analysis
  ##Arguments:
  ##None
  ##Returns:
  ##params
  ##A list of parameter values 
  
  #Asks for the input file
  inputFile <- readline(prompt = "Enter path to parameter file: ")
  #Reads in the file
  data <- read.csv(inputFile, header = FALSE, stringsAsFactors = FALSE) 
  #Extracts the data from the file into variables
  data <- data$V2
  threads <- trimws(data[2], which = "both")
  STARsource <- trimws(data[3], which = "both")
  annotation <- trimws(data[4], which = "both")
  FastaFile <- trimws(data[5], which = "both")
  readLength <- trimws(data[6], which = "both")
  fastqFile <- trimws(data[8], which = "both")
  isGzipped <- trimws(data[9], which = "both")
  GroupFile <- trimws(data[11], which = "both")
  isPairedEnd <- trimws(data[12], which = "both")
  isMulti <- trimws(data[13], which = "both")
  FDR <- trimws(data[15], which = "both")
  LowExpCounts <- trimws(data[16], which = "both")
  #Put the variables together in a list and return the list
  params <- c(threads, STARsource, annotation, FastaFile, readLength, fastqFile, isGzipped, GroupFile, isPairedEnd, isMulti, FDR, LowExpCounts)
  return(params)
}

GroupParamSetup <- function(rna_params){ 
  ##GroupParamSetup
  ##Description
  ##Reads a csv file containing parameter values related to the sample groups
  ##Arguments:
  ##rna_params: A list containing parameters
  ##Returns:
  ##csvParams: A dictionary containing perameters associated with the groups 

  #Reads the columns into variables
  group_csv <- read.csv(file = rna_params[8])
  Groups <- group_csv$Group
  GroupNum <- c()
  Group <- group_csv$Group.Number
  for (i in Group){
    if (length(i) == 1){
      Group <- paste("0", i, sep = "", collapse = NULL)  
    }
    GroupNum <- c(GroupNum, Group)
  }
  SecondID <- group_csv$Second.Identifier..Optional.
  SecondID <- ifelse(is.na(SecondID), " ", SecondID)
  
  #Put the variables in a list
  GroupParams <- data.frame(Groups, GroupNum, SecondID)
  
  SampleNames <- c()
  for (i in 1:length(GroupParams$Groups)){
    SampleNames <- c(SampleNames, trimws(paste(GroupParams$Groups[i],GroupParams$GroupNum[i],GroupParams$SecondID[i], sep = "", collapse = NULL), which = "both"))
  }
  
  Groups <- c() 
  for (i in GroupParams$Groups){
    if (! i %in% Groups){
      Groups <- c(Groups, i)
    }
  }
  return_list <- list("GroupParameters" = GroupParams, "Names" = SampleNames, "GroupNames" = Groups)
  return(return_list)
}

