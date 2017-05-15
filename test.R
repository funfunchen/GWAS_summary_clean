yu.files <- list.files("/net/twins/home/fangchen/TEST",pattern="sig_May_12$",full.names=T, recursive= F);
for (file in yu.files){ 
      col.name <- strsplit(basename(file), "[.]")[[1]][1]
      col.name <- paste0(col.name, ".")
      cols <- c("BETA", "PVALUE")
      cols <- setNames(cols, paste0(col.name, cols))
      # if the merged dataset doesn't exist, create it
      if (!exists("whole.file")){
              whole.file <- read_tsv(file) %>% rename_(.dots = cols)
       }
      # if the merged dataset does exist, append to it
      if (exists("whole.file")){
              temp_file <- read_tsv(file) %>% rename_(.dots = cols)
          whole.file <- left_join(whole.file, temp_file, by=c("CHROM"="CHROM", "POS"="POS"), copy=F) %>%
                   select()
        }
}
