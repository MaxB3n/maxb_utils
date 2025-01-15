nun <- function(lst){
  return( length(unique(lst)) )
}

substrRight <- function(string,start,stop = 0){
  #expects negative numbers
  l <- nchar(string)
  return(substr(string, l+start+1, l+stop))
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

freadDirectory <- function(directory = ".", Names = NULL, separate = T, patterns = "\\.csv|\\.tsv|\\.txt|\\.table"){
  files <- list.files(path = directory, pattern = patterns)
  if (is.null(Names)) {
    Names <- gsub(paste0(patterns, "|\\.gz"),"", files)
  } else {
    if(length(Names) < length(files)){
      append( Names, files[length(Names):length(files)] )
    }
  }
  for (i in 1:length(files)){
    assign(Names[i], data.table::fread(file.path(directory,files[i])), envir = parent.env(environment()))
    print(sprintf("Created data.table of %s named %s", files[i], Names[i]))
  }
}
