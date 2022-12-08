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

parser <- ArgumentParser(description = 'Extracting reads identified as a specific taxa(s)',
                         epilog = 'If you have any issues or sugestions for improvements please headover too:')

parser$add_argument('-k',
                    '--kraken_output',
                    dest = 'krakenFile',
                    type = 'character',
                    action = 'store',
                    required = TRUE,
                    help = 'The file name of the kraken output file you want to filter')

parser$add_argument('-o',
                    '--filtered_reads_output',
                    dest = 'outputFileName',
                    type = 'character',
                    required = TRUE,
                    help = 'The file name you want to give the filtered reads')

parser$add_argument('-t',
                    '--taxa',
                    dest = 'taxaID',
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
                    dest = 'krakenReportFilename',
                    type = 'character',
                    help = 'The taxa ID you want to filter by')

args <- parser$parse_args()


########################## creating the functions ##############################

#creating a function to pull out all the reads that match the IDs
filtReads <- function(ID, dataframe){
  #finding lines that match our ID
  listOfMatches <- grep(pattern = ID, dataframe)
  
  #creating the temp variables
  i <- 1
  fReads <- c()
  
  for (seqLine in listOfMatches) {
    #pulling out the lines with matches and adding the strings to a list
    fReads[i] <- dataframe[seqLine] #header
    fReads[(i+1)] <- dataframe[(seqLine+1)] #sequence Data
    
    #adding to the iterator
    i <- i+2
    
  }
  
  return(fReads)
}


#creating a function to pull out all taxonomic ID's bellow 
findIDsBellow <- function(ID, dataframe){
  #removing any levels above the supplied taxa ID
  lower <- grep(dataframe$V5, pattern = ID)
  upper <- length(dataframe$V5)
  trimmedReport <- dataframe[lower:upper, , ]
  
  #creating a list of kraken taxa codes 
  taxaSymbols <- c('R', 'D', 'K', 'P', 'C', 'O', 'F', 'G', 'S')
  
  #pulling out the taxa symbol and splitting it up into characters 
  pat <- stri_split_boundaries(trimmedReport$V4[1], type='character')
  
  #identifying level of taxa code supplied in list 
  level <- grep(taxaSymbols, pattern = pat[1])
  
  #checking to make sure ID has levels bellow it
  if(level == 9){
    
    #if it doesn't return the original ID
    return(ID)
    
  } else {
    
    #removing the supplied taxa ID from the df so it doesn't return a null
    trimmedReport2 <- trimmedReport[2:length(V4), ,]
    
    #using lapply to loop in C and output a list of taxa ID's
    lowerTaxa <- lapply(trimmedReport2$V4, function(x){
      
      #pulling out the taxa symbol and splitting it up into characters 
      pat2 <- stri_split_boundaries(x, type='character')
      
      #checking if the position of the next taxa symbol is higher or lower then the supplied taxa
      #in the taxaSymbols list 
      
      if (level < grep(taxaSymbols, pattern = pat2[1])){
        results <- trimmedReport[V4 == x, .(V5),]
        results$V5
      }
      
    })
    
    #combining the new list of taxa with the old supplied taxa
    listOfTaxa <- append(lowerTaxa, ID)
    return(listOfTaxa)
    
  }
  
}


############################## main workflow ###################################

main <- function(){
  #passing data from the flags to variables
  include_taxa_levels_bellow <- args$include_lower_taxa
  taxID <- args$taxaID
  krakenOut <- stri_read_lines(args$krakenFile)
  outputFilename <- args$outputFileName
  
  #creating a list to collect the output
  filteredReads <- c()
  
  #running the findIDsBellow function and filtering all ID's returned 
  if (include_taxa_levels_bellow == T) {
    #importing the data
    krakenReport <- fread(args$krakenReportFilename, header = F)
    
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
  
  #writing the output of the filtering to a file 
  stri_write_lines(filteredReads, outputFilename, sep = '\n')
  
  #outputting info to the user
  print(paste('number of taxa IDs found:', length(taxaList)))
  print(paste('number of reads found:', length(filteredReads)/2))
  print(paste('percentage of reads retained:', (length(filteredReads)/length(krakenOut)*100), '%'))
  print(paste('filtered reads have been written to:', outputFilename))
  
}

main()

#write an exclude function and a findIDsAbove function for the package 
