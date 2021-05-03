#-------------------------------------------------------------------------------------------
## Build an obj that will save the previous objs
# Function
obj <- function(obj_list, x) {
  n <- length(obj_list)
  obj_list[[n+1]] <- x  
  params_list <<- obj_list
}
#-------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------
## Pipeline
# import packages 
library("dplyr")
source("~/Desktop/workstation/metaNanoPype_helper_funs.R")

# Import data
metananopype[["params"]] <- metananopype[["data"]] <- metananopype <- list()
metananopype$params[["tax_path"]] <- "~/Desktop/workstation" 
files2import <- list.files(path = metananopype$params$tax_path, 
                           pattern = ".report", full.names = TRUE)
metananopype$params[["samples"]] <- gsub(pattern = ".report", replacement = "", 
                                         x = basename(files2import))
metananopype$data[["kraken2"]] <- lapply(X = setNames(files2import, metananopype$params$samples), 
                                         function(x, y = "\t", z = FALSE) read.table(file = x, sep = y, stringsAsFactors = z), 
                                         y = "\t", z = FALSE)

# Parse table
metananopype$data[["kraken2_parsed"]] <- lapply(metananopype$data$kraken2, function(x) parse_kraken2_report(x))
metananopype$data[["kraken2_parsed"]] <- lapply(metananopype$data$kraken2_parsed, 
                                                function(x, cols2sel = c(1,3)) x[,cols2sel], cols2sel = c(1,3))
for ( x in names(metananopype$data$kraken2_parsed)) colnames(metananopype$data$kraken2_parsed[[x]])[2] <- x
metananopype$data[["tax_tbl"]] <- merge_tax_tbl(list_tbls = metananopype$data$kraken2_parsed, merge_by = "Taxonomy")

# Import into phyloseq 

# Plot taxonomy

# Plot alpha diversity

# Plot beta diversity

#-------------------------------------------------------------------------------------------