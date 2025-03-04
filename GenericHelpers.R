nun <- function(lst){
  return( length(unique(lst)) )
}

substrRight <- function(string,start,stop = 0){
  #expects negative numbers
  l <- nchar(string)
  return(substr(string, l+start+1, l+stop))
}

substrLeft <- function(string,start,stop = 0){
  #expects negative numbers
  l <- nchar(string)
  return(substr(string, start, l+stop))
}


collapseVecPairwise <- function(vec, sep = "_"){
  n <- length(vec)
  if (n==2){
    return( paste(vec[1], vec[2], sep = sep) )
  } else {
    return( c( paste(vec[1], vec[2:n], sep = sep), collapseVecPairwise( vec[2:n], sep = sep )) )
  }
}


LoadPackagesFromCSV <- function(fpath, inclBioc = T){
  s
  packageTable <- read.csv(fpath)
  packageList <- packageTable$Package
  packagesToInstall <- packageList[!(list.of.packages %in% installed.packages()[,"Package"])]
  
  if(length(packagesToInstall)){
    
    if ( sum(grepl("Bioc", packageTable[packagesToInstall,]$Depends)) > 0){
      if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install( packageList[grepl("Bioc", packageTable[packagesToInstall,]$Depends)] )
      install.packages(packageList[!grepl("Bioc", packageTable[packagesToInstall,]$Depends)] )
    } else {
      install.packages( packagesToInstall )
    }
  }
}

#' FUNCTION: freadDirectory() 
#' This function reads all parse-able files into the parent environment as a data.table(s)
#' Some important functionalities include options to load files from all subdirectories, to supply names for each created table,
#' to limit the filetypes parsed, and to rbind all data into one data.table with the original file name (or supplied names) as a column. 
#' @param directory character vector that says where to look
#' @param Names list of character vectors to name created data.tables, applied in alphabetical order to files in directory. if this is   
#' shorter than the number of files to read, only the first files will be named and remaining files will use their filenames
#' @param recursive boolean to go into all subdirectories and look for files there too. this is all or nothing, so make sure only files you want to read are in there
#' @param patterns regex character vector to match the file extensions to parse. gzipped files matching this regex are also read. 
#' @param separate boolean to assign individual data.tables in parent envir or if FALSE, rbind into one data.table THAT IS RETURNED
#' @param verbose boolean to print names of files read and objects created
freadDirectory <- function(directory = ".", Names = NULL, recursive = T, patterns = "\\.csv|\\.tsv|\\.txt|\\.table", separate = TRUE, verbose = T){
  files <- list.files(path = directory, pattern = patterns, recursive = recursive)
  if (is.null(Names)) {
    Names <- gsub(paste0(patterns, "|\\.gz"),"", files)
  } else {
    if(length(Names) < length(files)){
      append( Names, files[length(Names):length(files)] )
    }
  }
  for (i in 1:length(files)){
    if (separate){
      assign(Names[i], data.table::fread(file.path(directory,files[i])), envir = parent.env(environment()))
      if (verbose) print(sprintf("Created data.table of %s named %s", files[i], Names[i]))
    } else {
      rbind()
    }
    
  }
}


