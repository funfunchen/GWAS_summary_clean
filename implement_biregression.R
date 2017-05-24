#use with find /net/twins/home/fangchen/TEST/ -name "*cut"|xargs -n1 -P5 -I{} Rscript --vanilla 2way.R -t {}

library(tidyverse);
library(getopt);

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
    file.chrom <- filter(file, CHROM == i) %>% arrange_(pvalue); #  use the standard evaluation versions of the dplyr functions (just append '_' to the
    new.data <- bind_rows(new.data, file.chrom[1,], .id = NULL) %>% filter(!is.na(CHROM));  ##find the smallest p value, pass to new.data
    sign.exist=TRUE;  ## get rid of the snps within the "window" of the snp we choose
    while(sign.exist){
      min.pos <- file.chrom[1,]$POS;
      file.chrom <- filter(file.chrom, !(POS>(min.pos-window.size/2) & POS<(min.pos+window.size/2))) %>% arrange_(pvalue);
      new.data <- bind_rows(new.data, file.chrom[1,], .id = NULL);
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
  
  
  
  opt.spec <- matrix(c('trait.file', 't', 1, "character"),
                     ncol=4, byrow=TRUE);
  A <- getopt(opt.spec); ## use exposure file as an argument with Rscript command in shell; combined with xargs to Parallel
  
  trait.name <- strsplit(basename(A$trait.file), "[.]")[[1]][1]; # get the traits name
  yu.file <- read_tsv(A$trait.file);
  yu.sig.file <- paste0("/net/twins/home/fangchen/TEST/", trait.name, ".result.sig_May_12"); # already have the files which all SNP are significant
  yu.sig <- read_tsv(yu.sig.file);
  
  bi_regression <- function(file, expo.file, expo.sig){
    # Calculate the ratio in each file
    #
    # Args:
    #  file: the summaries file
    #
    # Result:
    #  a list consist of trait name(exposure), file name(outcome) and the r ratio
    name <- basename(file);
    my.data <- read_tsv(file);
    X.set <- left_join(expo.sig, my.data, by=c("CHROM"="CHROM", "POS"="POS"), copy=F) %>%
      select(CHROM, POS, SNP, X_x_beta=BETA, X_x_p=PVALUE, X_y_beta=Beta, X_y_p=pvalue) %>% na.omit; # the Trait X set(without info of Y)
    Y.set <- filter(my.data, pvalue<=5*10^(-8)) %>% left_join(expo.file,  by=c("CHROM"="CHROM", "POS"="POS"), copy=F) %>%
      select(CHROM, POS, SNP, Y_y_beta=Beta, Y_y_p=pvalue, Y_x_beta=BETA, Y_x_p=PVALUE) %>% na.omit; # the Trait Y set (without info of X)
    X.set <- windowshift(X.set, 1000000, "X_x_p");  # the effect sizes on Trait X (without info of Y)
    Y.set <- windowshift(Y.set, 1000000, "Y_y_p");  # the effect sizes on Trait Y (without info of X)
    beta.xx <- as.numeric(X.set$X_x_beta); beta.xy <- as.numeric(X.set$X_y_beta);
    beta.yy <- as.numeric(Y.set$Y_y_beta); beta.yx <- as.numeric(Y.set$Y_x_beta);
    n_x <- length(beta.xx);  n_y <- length(beta.yy)
    if(length(beta.xx)<10||length(beta.yy)<10){
      result <- list(trait=trait.name, file=name, expo_beta=NA, expo_p=0, out_beta=NA, out_p=NA, t=NA, t_p=NA, df=NA);
    }else {
      x <- lm(beta.xy ~ beta.xx); y <- lm(beta.yx ~ beta.yy);
      beta_x <- summary(x)$coefficients[2,1]; se_x <- summary(x)$coefficients[2,2]; p_x <- summary(x)$coefficients[2,4]
      beta_y <- summary(y)$coefficients[2,1]; se_y <- summary(y)$coefficients[2,2]; p_y <- summary(y)$coefficients[2,4]
      df <- (se_x^2/n_x + se_y^2/n_y)^2/((se_x^2/n_x)^2/(n_x - 1) + (se_y^2/n_y)^2/(n_y - 1));
      t <- abs(beta_x - beta_y)/sqrt(se_x^2 + se_y^2);
      p <- 2*pt(t, df, lower.tail=F);
      result <- list(trait=trait.name, file=name, expo_beta=beta_x, expo_p=p_x, out_beta=beta_y, out_p=p_y, t=t, p=p, df=df);
    }
  }
  
  my.files <- list.files("/net/twins/home/fangchen/LD-hub/cleaned_data/",pattern="clean$",full.names=T, recursive= F);
  my_res <- do.call(rbind, lapply(my.files, bi_regression, expo.file=yu.file, expo.sig=yu.sig));
  my_res <- as.data.frame(my_res);
  my_res$file <- gsub(".clean$", "", my_res$file);
  my_res_sig <- my_res[which(as.numeric(unlist(my_res$expo_p))<0.01),]
  write.table(format(my_res, digits = 3), paste0(trait.name,".biway_may24.full"), quote=F, row.names=F, sep="\t")
  write.table(format(my_res_sig, digits = 3), paste0(trait.name,".biway_may24.sig"), quote=F, row.names=F, sep="\t")