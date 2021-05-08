#---------------------------------------------------------------------------------------------------
## Function to parse kraken2 report tables 
# 
parse_kraken2_report <- function(kraken2_report) {
  
  # 'parse_kraken2_report()': function to parse kraken2
  #report tables. kraken2 report need to have 6 columns: 
  #(1st) percentage of reads rooted in the clade of each
  #taxon; (2nd) number of reads rooted in the clade of 
  #each taxon; (3rd) number of reads assigned directly 
  #to each taxon; (4th) taxonomic rank; (5th) NCBI taxonomic 
  #ID; (6th) taxon. The function checks the full taxonomic 
  #path/name by appending the following taxonomic ranks 
  #(separated by "; "): "D" (domain), "P" (phylum), "C" 
  #(class), "O" (order), "F" (family), "G" (genus), "S" 
  #(species). Other ranks are not considered. For each 
  #full taxonomic path/name. It returns a parsed table 
  #with the following columns: "Taxonomy" (full taxonomic 
  #path/name semi-colon plus space separated); "Percentage" 
  #(percentage of reads that belong to each taxon - recovered
  #from the 1st kraken2 report table); "Abundance" (number of 
  #reads that were assigned/classified to each taxon - 
  #recovered from the 2nd column).
  # 'kraken2_report' (mandatory): a data frame with the 
  #kraken2 report output file. 
  
  # Check input 
  stopifnot( is.data.frame(kraken2_report) & ncol(kraken2_report) == 6 )

  # Define taxonomic ranks to include in the output
  ranks <- 1:7 # rank number domain --> species | 1 --> 7
  names(ranks) <-  c("D", "P", "C", "O", "F", "G", "S")
  tax_parsed <- data.frame(stringsAsFactors = FALSE) # output table parsed
  rank_tax <- c() # the full tax path
  for ( t in 1:nrow(kraken2_report) ) { # loop over each kraken2 report table row
    rank <- kraken2_report[t,4] # current tax rank 
    if ( ! ( rank %in% c("U", "R", names(ranks)) ) ) next # discarding any rank that is not among the canonical ranks
    if ( rank %in% c("U", "R") ) { # if it is 'unclassified' or 'root' do:
      if ( rank == "U" ) { # if it is 'unclassified' save to output:
        tax_parsed[nrow(tax_parsed)+1,"Taxonomy"] <- kraken2_report[t,6]
        tax_parsed[nrow(tax_parsed),"Percentage"] <- kraken2_report[t,1]
        tax_parsed[nrow(tax_parsed),"Abundance"] <- kraken2_report[t,2]
      }
      next
    }
    rank_no <- ranks[rank] # curr rank no 
    if ( rank_no > length(rank_tax) ) { # if rank_no > than curr rank tax list, keep adding taxon:
      rank_tax <- append(rank_tax, trimws(kraken2_report[t,6])) # the full tax path
      next
    }
    if ( rank_no <= length(rank_tax) ) { # if rank_no <= than curr rank tax list, 
      #save output & update rank list to the curr rank no.:
      # save the previous tax
      tax_parsed[nrow(tax_parsed)+1,"Taxonomy"] <- paste(rank_tax, collapse = "; ")
      tax_parsed[nrow(tax_parsed),"Percentage"] <- kraken2_report[t-1,1]
      tax_parsed[nrow(tax_parsed),"Abundance"] <- kraken2_report[t-1,2]
      # parse the next tax
      rank_tax <- c(rank_tax[0:(rank_no-1)], trimws(kraken2_report[t,6])) # update the full tax path
    }
  }
  # save the last full tax path in the table
  tax_parsed[nrow(tax_parsed)+1,"Taxonomy"] <- paste(rank_tax, collapse = "; ")
  tax_parsed[nrow(tax_parsed),"Percentage"] <- kraken2_report[t-1,1]
  tax_parsed[nrow(tax_parsed),"Abundance"] <- kraken2_report[t-1,2]
  return(tax_parsed)
}

## Function to gather multiple tables based on full tax path/name
#
merge_tax_tbl <- function(list_tbls, merge_by) {
  
  # 'merge_tax_tbl()': merge tables (data.frame) 
  #including all values found in the left and 
  #right tables. 
  # 'list_tbls' (mandatory): list of data frames
  #to merge by 'merge_by' (all the tables need to 
  #include this column name). 
  # 'merge_by' (mandatory): column name to merge by
  #(character of length one).s
  
  # Check input 
  stopifnot(is.list(list_tbls) & all( unlist(lapply(list_tbls, function(x) is.data.frame(x))) ) )
  stopifnot( all( unlist(lapply(list_tbls, function(x) merge_by %in% colnames(x))) ) )
  stopifnot(is.character(merge_by) & length(merge_by)==1)
  
  # Merge tbls
  merged_tbl <- NULL
  for ( tbl in list_tbls ) { # loop over tables
    if ( is.null(merged_tbl) ) { # assign the 1st tbl to output 
      merged_tbl <- tbl
    } else { # merge the previous table to the current tbl
      merged_tbl <- merge(x = merged_tbl, y = tbl, by = merge_by, all = TRUE)
    }
  }
  merged_tbl[is.na(merged_tbl)] <- 0 # substitute NA by 0
  rownames(merged_tbl) <- paste0("seq_", 1:nrow(merged_tbl))
  return(merged_tbl)
}

## create function to parse tax into tax_table() for phyloseq
# 
tax_to_physeq <- function(full_tax, split_by = "; ", 
                          tax_ranks = c("Kingdom", "Phylum", 
                                        "Class", "Order", 
                                        "Family", "Genus", 
                                        "Species")) {
  
  # 'tax_to_physeq()': convert a vector with full taxonomic 
  #names into a data frame with the different taxonomic names 
  #separated by taxonomic rank. Named vector with seq id identification. 
  # 'full_tax': character vector with full taxonomic names, e.g., 
  #"Bacteria; Actinobacteriota; Acidimicrobiia; Microtrichales; Ilumatobacteraceae; uncultured".
  # 'split_by': character symbol to split the full taxonomic 
  #names. By default is '; '.
  # 'tax_ranks': character vector of taxonomic ranks to use 
  #to build the taxonomic data frame. If a rank does not exist
  #it will be filled in with NAs.
  
  # build data frame 
  tax_tbl <- data.frame(
    stringsAsFactors = FALSE 
  )
  tx_no <- length(tax_ranks)
  
  # split full taxonomy
  tax = lapply(full_tax, function(x) strsplit(x = x, split = split_by)[[1]])
  
  # loop over each full taxonomic name and fill the table
  ct <- 0
  for ( tx in tax ) {
    ct <- ct + 1 
    # if full taxonomic names do not have all the ranks fill with NAs
    if ( length(tx) < tx_no ) tx[(length(tx)+1):tx_no] <- NA
    tax_tbl[ct, tax_ranks] <- tx[1:tx_no]
  }
  row.names(tax_tbl) <- names(full_tax)
  tax_tbl <- phyloseq::tax_table(as.matrix(tax_tbl))
  return(tax_tbl)
}
#---------------------------------------------------------------------------------------------------
