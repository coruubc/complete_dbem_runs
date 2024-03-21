#' Double check DBEM runs
#'
#'This function double checks all runs from the DBEM are completed correctly
#'within Compute Canada. It cross-references the results files with
#'the settings and sbem inpout files. Specifically, the number of species ran,
#'that each species has the correct number of years, and will output that info
#'
#'
#' @author Juliano Palacios Abrantes | j.palacios@oceans.ubc
#' @param results_path path where results are if not on personal scracth folder
#' @param dbem_path path where the DBEM inout files exist
#' @param settings settings file for DBEM run. Preallocation Settings1F
#' @return
#'
#' @export
dbem_run_check <- function(results_path = NA, dbem_path = NA, settings = "Settings1F"){
  
  # ----------------#
  # Packages needed
  # ----------------#
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(stringr))
  
  #Extract global variables
  user <- Sys.info()[8] # Extracts user name for path
  
  # Get DBEM input info.
  if(is.na(dbem_path)){
    dbem_path <- paste0("/home/",user,"/projects/rrg-wailung/",user,"/Fortran/")
  }
  
  # Load settings file
  settings_file <- read.table(paste0(dbem_path,"Data/",settings,".txt"),
                              quote="\"", comment.char="")
  colnames(settings_file) <- c("variable","value")
  
  # Get result folder
  rpath <- settings_file$value[5]
  
  # Determines main data pathways if not personally included
  if(is.na(results_path)){
    results_path <- paste0("/home/",user,"/scratch/Results/",rpath)
  }
  
  # ---------- #
  ## Data checking protocol
  # ---------- #
  
  # Check how many species were ran
  spp_list <- list.files(results_path)
  taxon_ran <- length(spp_list)
  
  # Check number of years per species
  for(s in 1:taxon_ran){
    # Get info needed
    taxa <- spp_list[s]
    
    read_taxa <- list.files(paste0(results_path,"/",taxa), pattern = taxa,full.names = T)
    
    # Missing data
    if(length(read_taxa) ==0){
      missing_msg <- paste("No data for",taxa)
      next()
    }else{
      missing_msg <- paste("Data for",taxa,"good")
    }
    
    
    catch_runs <- str_subset(read_taxa, pattern = "Catch")
    abd_runs <- str_subset(read_taxa, pattern = "Abd")
    
    # Printing results
    n_files <- length(read_taxa)
    
    # Abundance information
    n_abd <- length(abd_runs)
    
    # Year one ran
    a_i_year <- str_sub(abd_runs[1],
                        str_count(abd_runs[1]) - 7,
                        str_count(abd_runs[1]) - 4
    )
    
    # Year N ran
    a_n_year <- str_sub(abd_runs[length(abd_runs)],
                        str_count(abd_runs[n_abd]) - 7,
                        str_count(abd_runs[n_abd]) - 4
    )
    
    # Catch information
    n_catch <- length(catch_runs)
    
    # Year one ran
    c_i_year <- str_sub(catch_runs[1],
                        str_count(catch_runs[1]) - 7,
                        str_count(catch_runs[1]) - 4
    )
    
    # Year N ran
    c_n_year <- str_sub(catch_runs[n_catch],
                        str_count(catch_runs[n_catch]) - 7,
                        str_count(catch_runs[n_catch]) - 4
    )
    
    # Return DF
    taxon_df <- data.frame(
      variable = c("nfiles" ,"n_abd_y","abd_i_year","abd_n_year","n_catch_y","catch_i_year","catch_n_year"),
      value = c(n_files,n_abd,a_i_year,a_n_year,n_catch,c_i_year,c_n_year)
    ) %>%
      mutate(taxon_key = taxa)
    
    if(s == 1){
      
      final_df <- taxon_df
      
    }else{
      
      final_df <- bind_rows(final_df,taxon_df)
    }
  }
  
  # Find mismatching taxa
  complete_df <- final_df %>%
    group_by(variable,value) %>%
    summarise(n = n(),
              taxa = paste(unique(taxon_key),collapse =  ";"),
              .groups = "keep"
    )
  
  
  
  # ---------- #
  ## Save check Message ##
  # ---------- #
  sink(paste0(results_path,"/0_dbem_check.txt"))
  
  # Summary
  summary_msg <- paste("Runs completed for", settings_file$value[7], "species in groups of", settings_file$value[1], "species each, using the",
                       settings_file$value[2],"ESM under SSP",settings_file$value[3],
                       "\n Fishing level allwed in HS:",settings_file$value[8],"\n Fishing level allowed in EEZ:",settings_file$value[9],
                       "\nMPA scenario:",settings_file$value[10],".\n Results saved at",rpath,". See bottom for settings file"
  )
  cat("# Summary\n", summary_msg, "\n\n")
  
  # Number of species ran
  n_species_msg <- paste("Ran",taxon_ran,"species with",settings_file$value[1],"species per run")
  cat("# Number of species ran\n", n_species_msg, "\n\n")
  
  
  # Missing data
  cat("# Missing data\n", missing_msg, "\n\n")
  
  # Potential miss runs
  cat("# Potential miss runs\n")
  
  if(nrow(complete_df) > 7){
    
    print(complete_df, n = 1000)
    
  }else{
    cat("All runs are OK\n\n")
  }
  
  # Settings file
  cat("# Settings file\n\n")
  
  print(settings_file)
  
  # Close the connection
  sink()
  
  # Return message
  ret_msg <- "Check up complete. Please see the dbem_check.txt file in the result's repository"
  
  return(ret_msg)
  
}

dbem_run_check()