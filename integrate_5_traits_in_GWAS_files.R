
yu.files <- list.files("/net/twins/home/fangchen/TEST",pattern="cut$",full.names=T, recursive= F);
## integrate all the smoking and drinking files to get the 5 traits
for (file in yu.files){ 
  col.name <- strsplit(basename(file), "[.]")[[1]][1]
  col.name <- paste0(col.name, ".")
  cols <- c("BETA", "PVALUE")
  cols <- setNames(cols, paste0(col.name, cols))
  # if the merged dataset doesn't exist, create it
  if (!exists("whole.file")){
    whole.file <- read_tsv(file) %>% rename_(.dots = cols) %>% arrange(CHROM, POS) ;
    next
  } 
  # if the merged dataset does exist, append to it
  if (exists("whole.file")){
    temp_file <- read_tsv(file) %>% rename_(.dots = cols) %>% arrange(CHROM, POS)
    whole.file <- left_join(whole.file, temp_file, by=c("CHROM"="CHROM", "POS"="POS"), copy=F) %>% na.omit
  }
}

na.omit(whole.file) %>% select(CHROM, POS, matches("BETA|PVALUE"))
test.file <- filter(test.file, cpd.PVALUE<=5*10^(-8)|dnd.PVALUE<=5*10^(-8)|dpw.PVALUE<=5*10^(-8)|sc.PVALUE<=5*10^(-8)|si.PVALUE<=5*10^(-8))

# Multivariate lm ---------------------------------------------------------

sig_file <- read_tsv("total_5_traits_beta.May15")

#library(lazyeval)
windowshift <- function(file, window.size, pvalue){
  # Get the most significant SNP signla in the window.size region
  #
  # Args:
  #    file: the GWAS/meta-analysis file
  #    window: the window.size, could be 10^6 or 2*10^6
  #    pvalue: the pvaule you want to use, especially for multi-traits, use the quotation mark("") 
  # 
  # Return: the sorted file 
  new.data <- file[0,]
  for (i in unique(file$CHROM)){
    file.chrom <- filter(file, CHROM == i) %>% arrange_(pvalue); #  use the standard evaluation versions of the dplyr functions (just append '_' to the function names, ie. group_by_ & summarise_) 
    new.data <- bind_rows(new.data, file.chrom[1,], .id = NULL) %>% filter(!is.na(CHROM));  ##find the smallest p value, pass to new.data
    sign.exist=TRUE;  ## get rid of the snps within the "window" of the snp we choose
    while(sign.exist){
      min.pos <- file.chrom[1,]$POS;
      file.chrom <- filter(file.chrom, !(POS>(min.pos-window.size/2) & POS<(min.pos+window.size/2))) %>% arrange_(pvalue);
      new.data <- bind_rows(new.data, file.chrom[1,], .id = NULL);
      if(length(file.chrom$POS)==0) {
        sign.exist=FALSE;
      }
    }
  }
  new.data <- filter(new.data, !is.na(CHROM)) %>% arrange(CHROM,POS);
}

mylist <- list(); 
for (file in  summary.files) {
  sum_data <- read_tsv(file) %>% select(CHROM, POS, Beta, pvalue) %>% rename(beta.GB=Beta, p.GB= pvalue)
  name <- basename(file);
  new.form <- left_join(sig_file, sum_data, by=c("CHROM"="CHROM", "POS"="POS"), copy=F) %>% na.omit; 
  cpd <- filter(new.form, cpd.PVALUE <= 5*10^(-8)) %>% windowshift(1000000, "cpd.PVALUE");
  dnd <- filter(new.form, dnd.PVALUE <= 5*10^(-8)) %>% windowshift(1000000, "dnd.PVALUE");
  dpw <- filter(new.form, dpw.PVALUE <= 5*10^(-8)) %>% windowshift(1000000, "dpw.PVALUE");
  sc <- filter(new.form, sc.PVALUE <= 5*10^(-8)) %>% windowshift(1000000, "sc.PVALUE");
  si <- filter(new.form, si.PVALUE <= 5*10^(-8)) %>% windowshift(1000000, "si.PVALUE");
  temp <- dplyr::union(cpd, dnd, dpw, si, sc) %>% na.omit;
  test.beta <- lm(beta.GB ~ dnd.BETA + dpw.BETA + cpd.BETA + sc.BETA + si.BETA, data = temp );
  result <- list(file=name, dnd.beta=summary(test.beta)$coefficients[2,1], dnd.pvalue=summary(test.beta)$coefficients[2,4], 
                 dpw.beta=summary(test.beta)$coefficients[3,1], dpw.pvalue=summary(test.beta)$coefficients[3,4], 
                 cpd.beta=summary(test.beta)$coefficients[4,1], cpd.pvalue=summary(test.beta)$coefficients[4,4], 
                 sc.beta=summary(test.beta)$coefficients[5,1], sc.pvalue=summary(test.beta)$coefficients[5,4],
                 si.beta=summary(test.beta)$coefficients[6,1], si.pvalue=summary(test.beta)$coefficients[6,4]);
  mylist <- rbind(mylist, result);  
}
my_res <- as.data.frame(mylist, row.names = F, stringsAsFactors=F) 
my_res <- as.data.frame(lapply(my_res,unlist))
my_res$file <- gsub(".clean$","",my_res$file)
write.table(my_res, "Multi_Overall_5traits_may_16.full", quote=F, row.names=F, sep="\t")
# my_res <- as_tibble(my_res)
