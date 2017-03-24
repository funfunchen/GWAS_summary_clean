yu.files <- list.files("/net/twins/home/fangchen/TEST/Gscan_Mar8",pattern="result$",full.names=T, recursive= F); ## read the GSCAN files in the folder
my.files <- list.files("/net/twins/home/fangchen/LD-hub/cleaned_data/MAGNETIC",pattern="clean$",full.names=T, recursive= F); ## read the GWAS summary files in the folder
for (yufile in yu.files){  ## loop in GSCAN files, first level loop
  yu <- read_delim(yufile, "\t")
  sig_yu <- yu %>% filter(PVALUE <= 5*10^(-8)) %>% select(CHROM, POS, REF, ALT, AF, PVALUE, BETA, SE, N) %>% arrange(CHROM, POS); ##find significant SNPs
  mylist <- list();  ## create empty list for the results
  outfilename <- basename(yufile); ## get the basemame of output files
  for (file in my.files) {  ## loop in GWAS files, second level loop
    my_data <- read_tsv(file);
    name <- basename(file);
    new.form <- left_join(sig_yu, my_data, by=c("CHROM"="CHROM", "POS"="POS"), copy=F) %>% rename(beta.GA=BETA, beta.GB=Beta) %>% filter(!is.na(SNP)) %>% 
      select(CHROM, POS, SNP, REF=REF.x, ALT=ALT.x, beta.GA.p=PVALUE, beta.GA, beta.GB, beta.GB.p=pvalue);  ## combine the GSCAN and summary files for the sig SNPs
    new.data <- new.form[0,];
    window.size <- 1000000; ## window shift
    for (i in unique(new.form$CHROM)){
      file.chrom <- filter(new.form,CHROM == i) %>% arrange(beta.GA.p);
      new.data <- bind_rows(new.data,file.chrom[1,], .id = NULL) %>% filter(!is.na(CHROM));  ##find the smallest p value, pass to new.data
      sign.exist=TRUE;  ## get rid of the snps within the "window" of the snp we choose
      while(sign.exist){
        min.pos <- file.chrom[1,]$POS;
        file.chrom <- filter(file.chrom,!(POS>(min.pos-window.size/2) & POS<(min.pos+window.size/2))) %>% arrange(beta.GA.p);
        new.data <- bind_rows(new.data,file.chrom[1,], .id = NULL);
        if(length(file.chrom$POS)==0) {
          sign.exist=FALSE;
        }
      }
    }
    new.data <- filter(new.data,!is.na(CHROM)) %>% arrange(CHROM,POS);
    if (nrow(new.data) == 0) next; ## make sure there exists beta in the data frame, otherwise move to the next file 
    test.beta <- lm(new.data$beta.GB ~ new.data$beta.GA);
    result <- list(file=name,beta=summary(test.beta)$coefficients[2,1], pvalue=summary(test.beta)$coefficients[2,4]);
    mylist <- rbind(mylist, result);  ## append the result
  }
  my_res <- as.data.frame(mylist, row.names = F, stringsAsFactors=F) 
  my_res <- as.data.frame(lapply(my_res,unlist))
  my_res$file <- gsub(".clean$","",my_res$file)
  my_res_sig <- my_res[which(as.numeric(my_res$pvalue)<0.001),]
  write.table(my_res, paste0(outfilename,".MAGNETIC.full"), quote=F, row.names=F, sep="\t")
  write.table(my_res_sig, paste0(outfilename,".MAGNETIC.sig"), quote=F, row.names=F, sep="\t")
}