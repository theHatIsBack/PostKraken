#!/usr/bin/env Rscript

################################################################################
#
# filtering kraken output
# PostKraken
# Cameron Ferguson, 06-12-22
#
################################################################################

########################## preparing the workspace #############################

#importing the different library's 
library(data.table)
library(stringi)
library(argparse)


############################ creating the flags ################################

parser <- ArgumentParser(prog = 'filtering_kraken_output.R',
                         description = 'Extracting reads identified as a specific taxa(s)',
                         epilog = 'If you have any issues or sugestions for improvements please headover too:')

parser$add_argument('-k',
                    '--kraken_output',
                    dest = 'kraken_File',
                    type = 'character',
                    action = 'store',
                    required = TRUE,
                    help = 'The file name of the kraken output file you want to filter')

parser$add_argument('-o',
                    '--output_filename',
                    dest = 'output_FileName',
                    type = 'character',
                    required = TRUE,
                    help = 'The file name you want to give the filtered reads')

parser$add_argument('-m',
                    '--method',
                    dest = 'method',
                    choices = c('filter', 'exclude'),
                    nargs = '?',
                    default = 'filter',
                    type = 'character',
                    required = TRUE,
                    help = 'filter finds reads that match specified taxa IDs. exclude finds reads that do not match the specified taxa ID')

parser$add_argument('-t',
                    '--taxa',
                    dest = 'taxa_ID',
                    type = 'character',
                    required = TRUE,
                    help = 'The taxa ID you want to filter by')

parser$add_argument('--include_lower_taxa',
                    choices = c(T, F),
                    nargs = '?',
                    default = 'F',
                    type = 'logical',
                    help = 'whether to include lower taxa or not. Defualt is F if the T value is passed you will need to include the -r flag for it to work')

parser$add_argument('-r',
                    '--kraken_report',
                    dest = 'kraken_Report_Filename',
                    type = 'character',
                    help = 'The taxa ID you want to filter by')

args <- parser$parse_args()


########################## creating the functions ##############################

#creating a function to pull out all the reads that match the IDs
filtReads <- function(ID, dataframe){
  #creating the regular exspression search term
  regID <- paste('\\b', ID, '\\b', sep = '')
  
  #converting the list of kraken output to a data.table
  table <- as.data.table(matrix(krakenOut, 
                                ncol = 2, 
                                nrow = length(krakenOut)/2, 
                                byrow = T))
  
  #using lapply to loop in c and find the header lines that match our ID and 
  #return a vector of positions 
  index <- sort(unlist(lapply(regID, function(x){ grep(pattern = x, table$V1) })))
  
  #filtering the table for the sequences of interest
  fTable <- table[index, ,]
  
  #converting the data.table back to a list
  fReads <- unlist(as.list(t(fTable)))
  
  return(fReads)
}


#creating a function to exclude all the reads that match the IDs
exReads <- function(ID, dataframe){
  #ensuring the value passed is the correct datatype
  ID <- as.list(ID)
  
  #looping through the IDs and removing them from the seq file
  for (i in 1:length(ID)) {
    #creating the regular expression search term
    regID <- paste('\\b', ID[i], '\\b', sep = '')
    
    #finding lines that match our ID
    listOfMatches <- grep(pattern = regID, dataframe)
    
    #creating the temp variables
    index <- 1
    toExclude <- c()
    
    #creating a list of positions to exclude 
    for (seqLine in listOfMatches) {
      toExclude[index] <- seqLine #header
      toExclude[(index + 1)] <- (seqLine + 1) #sequence Data
      
      #adding to the iterator
      index <- index + 2
      
    }
    
    #excluding the matching reads 
    dataframe <- dataframe[-toExclude] 
  }
  
  return(dataframe)
}


#creating a function to pull out all taxonomic ID's bellow 
findIDsBellow <- function(ID, dataframe){
  #creating the regular exspression search term
  regID <- paste('\\b', ID, '\\b', sep = '')
  
  #removing any levels above the supplied taxa ID
  lower <- grep(pattern = regID, dataframe$V5)
  upper <- length(dataframe$V5)
  trimmedReport <- dataframe[lower:upper, , ]
  
  #creating a list of kraken taxa codes 
  taxaSymbols <- c('R', 'D', 'K', 'P', 'C', 'O', 'F', 'G', 'S')
  
  #pulling out the taxa symbol and splitting it up into characters 
  pat <- stri_split_boundaries(trimmedReport$V4[1], type = 'character')
  
  #identifying level of taxa code supplied in list 
  level <- grep(pattern = pat[[1]][1], taxaSymbols)
  
  #checking to make sure ID has levels bellow it and if it doesn't return the original ID
  if(level == 9){
    
    return(ID)
    
  } else {
    
    #removing the supplied taxa ID from the df so it doesn't return a null
    trimmedReport2 <- trimmedReport[2:length(V4), ,]
    
    #creating the list to populate with taxa codes
    lowerTaxa <- c()
    
    #the index serves two perpouses to act as a way to acess elements of the list and to keep 
    #track of the row of the data.table
    index <- 1
    
    #using a for loop to output a list of lower taxa ID's
    for (x in trimmedReport2$V4) {
      #pulling out the taxa symbol and splitting it up into characters 
      pat2 <- stri_split_boundaries(x, type = 'character')
      
      #checking if the position of the next taxa symbol is higher or lower then the supplied taxa
      #in the taxaSymbols list 
      
      if (level < grep(pattern = pat2[[1]][1], taxaSymbols)){
        results <- trimmedReport2[index, .(V5),]
        lowerTaxa[index] <- results$V5
        index <- index + 1
        
      } else if (level == grep(pattern = pat2[[1]][1], taxaSymbols) & is.na(pat2[[1]][2]) != is.na(pat[[1]][2])) {
        results <- trimmedReport2[index, .(V5),]
        lowerTaxa[index] <- results$V5
        index <- index + 1
        
      } else {
        break
      }
    }
    
    #combining the new list of taxa with the old supplied taxa
    listOfTaxa <- append(lowerTaxa, ID)
    return(listOfTaxa)
    
  }
  
}


############################## main workflow ###################################

main <- function(){
  #passing data from the flags to variables
  include_taxa_levels_bellow <- args$include_lower_taxa
  taxID <- args$taxa_ID
  krakenOut <- stri_read_lines(args$kraken_File)
  outputFilename <- args$output_FileName
  Method <- args$method
  
  #creating a list to collect the output
  filteredReads <- c()
  
  if (Method == 'filter') {
    #running the findIDsBellow function and filtering all ID's returned 
    if (include_taxa_levels_bellow == T) {
      #importing the data
      krakenReport <- fread(args$kraken_Report_Filename, header = F)
      
      #pulling out all the need ID's from the report
      taxaList <- findIDsBellow(taxID, krakenReport)
      
      #looping through all of the taxa ID's
      for (taxa in taxaList) {
        intReads <- filtReads(taxa, krakenOut)
        filteredReads <- append(filteredReads, intReads)
      }
      
      #just running the filtering for the ID provided 
    } else {
      taxaList <- c(taxID)
      filteredReads <- filtReads(taxaList[1], krakenOut)
      
    }
    
  } else if (Method == 'exclude'){
    #running the findIDsBellow function and removing all ID's returned
    if (include_taxa_levels_bellow == T) {
      #importing the data
      krakenReport <- fread(args$kraken_Report_Filename, header = F)
      
      #pulling out all the need ID's from the report
      taxaList <- findIDsBellow(taxID, krakenReport)
      
      #removing all of the ID's that where found 
      filteredReads <- exReads(taxaList, krakenOut)
       
      #just running the filtering for the ID provided 
    } else {
      taxaList <- c(taxID)
      filteredReads <- filtReads(taxaList, krakenOut)
      
    }
    
  }
  
  #writing the output of the filtering to a file 
  stri_write_lines(filteredReads, outputFilename, sep = '\n')
  
  #outputting info to the user
  print(paste('number of taxa IDs found:', length(taxaList)))
  print(paste('number of reads found:', length(filteredReads)/2))
  print(paste('percentage of reads retained:', (length(filteredReads)/length(krakenOut)*100), '%'))
  print(paste('filtered reads have been written to:', outputFilename))
  
}

main()

#write a findIDsAbove function for the package 
