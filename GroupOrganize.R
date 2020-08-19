SelectGroup1 <- function(group_params){
  ##SelectGroup1
  ##Description:
  ##The function asks the user their group of interest and returns Index values of the group within the parameters  
  ##Arguments:
  ##group_params: A list of parameter values
  ##Returns:
  ##return_list: Index values
  print("Which groups would you like to put through analysis?")
  #Ask the user for the first group
  Group1 <- readline(prompt = "Group 1: ")
  IndexCounter <- 0
  Index <- c()
  #Obtain the index numbers for the group selected
  for (group in group_params$GroupParameters$Groups){
    IndexCounter <- IndexCounter + 1
    if (Group1 == group){
      Index <- c(Index, IndexCounter)
    }
  }
  if (length(Index) == 0){
    readline(prompt = "ERROR:The Group you have selected does not match any group in your counts table")
  }
  return_list <- list("Group" = Group1, "Index" = Index)
  return(return_list)
}


SelectGroup2 <- function(group_params){
  ##SelectGroup2
  ##Description:
  ##The function asks the user their group of interest and returns Index values of the group within the parameters  
  ##Arguments:
  ##group_params: A list of parameter values
  ##Returns:
  ##return_list: Index values
  Index <- c()
  IndexCounter <- 0
  #Asks the user for the second group
  Group2 <- readline(prompt = "Group 2: ")
  #Obtain the index numbers for the second group selected
  for (group in group_params$GroupParameters$Groups){
    IndexCounter <- IndexCounter + 1
    if (Group2 == group){
      Index <- c(Index, IndexCounter)
    }
  }
  if (length(Index) == 0){
    readline(prompt = "ERROR:The Group you have selected does not match any group in your counts table")
  }
  return_list <- list("Group" = Group2, "Index" = Index)
  return(return_list)
}

SelectGroup_secondary <- function(group_params, group_index){
  ##SelectGroup_secondary
  ##Description:
  ##The function takes an index number 
  ##Arguments:
  ##group_params: A list of parameter values
  ##Returns:
  ##return_list: Index values
  Index <- c()
  IndexCounter <- 0
  #Obtain the index numbers for the second group selected
  for (group in group_params$GroupParameters$Groups){
    IndexCounter <- IndexCounter + 1
    if (toupper(group_index$Group) == toupper(group)){
      Index <- c(Index, IndexCounter)
    }
  }
  return_list <- list("Group" = group_index$Group, "Index" = Index)
  return(return_list)
}

SampleRemove <- function(group_params, raw_counts){
  ##Sample Check 
  ##Description
  ##A function which asks the user if they want to remove any samples from their data before continuing on to analysis
  ##Arguments:
  ##None
  ##Returns:
  ##return_list: a list of the group and counts associated with it 
  
  #Set the while loop condition
  library(dplyr)
  condition = TRUE
  index <- c()
  
  while(condition ==  TRUE){
    #Ask the user if they want to remove a group
    DumpGroup <- toupper(readline(prompt = "Do you want to dump a group? (If not, type NO): "))
    
    #Remove group specified from group_params
    index_counter <- 0
    for (name in group_params$Names){
      index_counter <- index_counter + 1
      
      if (toupper(name) == DumpGroup){
        index <- c(index, index_counter)
      }
      
      else if(DumpGroup == "NO"){     
        condition <- FALSE
      }
    }
    print(paste(DumpGroup, "group has been removed", sep = " ", collapse = NULL))
  }
  #Using the index, remove the groups from the dataset
  group_params$GroupParameters <- group_params$GroupParameters[-index, ]
  group_params$Names <- group_params$Names[-index]
  raw_counts <- raw_counts[ ,-index]
  colnames(raw_counts) <- group_params$Names
  return_list <- list("raw_counts" = raw_counts, "group_params" = group_params)
  return(return_list)
}