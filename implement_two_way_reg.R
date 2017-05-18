# 
# z_trans <- function(x, y){
#   rho<- cor(x, y, method = "spearman");
#   z <- 0.5 * (log(1+rho) - log(1-rho));
# }
# 
# z_x_h <- z_trans(beta_xx, beta_xy)
# z_y_h <- z_trans(beta_yy, beta_yx)
# 
# lh_M <- function(x_hat, x, n_x, y_hat, y, n_y){
#    dnorm(x_hat, x, sqrt(1/(length(n_x)-3)))*dnorm(y_hat, y, sqrt(1/(length(n_y)-3)))
# };
# lh_M1 <- lh_M(z_x_h, z_x_h, x, z_y_h, 0, y)
# lh_M2 <- lh_M(z_x_h, 0, x, z_y_h, z_y_h, y)
# lh_M3 <- lh_M(z_x_h, 0, x, z_y_h, 0, y)
# 
# lh_fun_M4 <- function(z){
#   -1*dnorm(z_x_h, z, sqrt(1/(length(x)-3)))*dnorm(z_y_h, z, sqrt(1/(length(y)-3)))
# };
# test <- nlminb(0.1, lh_fun_M4, gradient = NULL);
# lh_M4 <- -1*lh_fun_M4(test$par)
# 
# AIC_score_M1 <- 2-2*log(lh_M1)
# AIC_score_M2 <- 2-2*log(lh_M2)
# AIC_score_M3 <- 2-2*log(lh_M3)
# AIC_score_M4 <- 2-2*log(lh_M4)
# r <- exp((min(AIC_score_M1,AIC_score_M2)-min(AIC_score_M3,AIC_score_M4))/2)

# test --------------------------------------------------------------------

# /net/twins/home/fangchen/TEST/si.cut
# /net/twins/home/fangchen/LD-hub/cleaned_data/EduYears_Main.clean
# /net/twins/home/fangchen/TEST/si.result.sig_May_12
# /net/twins/home/fangchen/LD-hub/cleaned_data/body_fat_percentage_GWAS.clean
# 
# edu.file <- read_tsv("/net/twins/home/fangchen/LD-hub/cleaned_data/EduYears_Main.clean")
# si.cut <- read_tsv("/net/twins/home/fangchen/TEST/si.cut") %>% arrange(CHROM, POS)
# si_sig <- read_tsv("/net/twins/home/fangchen/TEST/si.result.sig_May_12")
# edu.set <- filter(edu.file, pvalue<=5*10^(-8)) %>% left_join(si.cut,  by=c("CHROM"="CHROM", "POS"="POS"), copy=F) %>% 
#   select(CHROM, POS, SNP, Y_edu_beta=Beta, Y_edu_p=pvalue, Y_si_beta=BETA, Y_si_p=PVALUE) %>% na.omit
# si.set <- left_join(si_sig, edu.file, by=c("CHROM"="CHROM", "POS"="POS"), copy=F) %>% 
#   select(CHROM, POS, SNP, X_si_beta=BETA, X_si_p=PVALUE, X_edu_beta=Beta, X_edu_p=pvalue) %>% na.omit
# edu.set <- windowshift(edu.set, 1000000, "Y_edu_p")
# si.set <- windowshift(si.set, 1000000, "X_si_p")
# z_x_h <- z_trans(si.set$X_si_beta, si.set$X_edu_beta)
# z_y_h <- z_trans(edu.set$Y_edu_beta, edu.set$Y_si_beta)
# x <- si.set$X_si_beta
# y <- edu.set$Y_edu_beta


#######################################################
# screeen test in the folder  ---------------------------------------------

library(tidyverse)

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


# Fisher’s Z-transformation -----------------------------------------------
z_trans <- function(x, y){
  # Compute the Fisher’s Z-transformation of each data set
  #
  #Args:
  #  x: vectors of effect sizes on Traint X 
  #  y: vectors of effect sizes on Traint Y
  #
  #Result:
  #  z: z hat
  rho<- cor(x, y, method = "spearman",, use = "pairwise.complete.obs");
  z <- 0.5 * (log(1+rho) - log(1-rho));
}


# likelyhood function -----------------------------------------------------
lh_M <- function(x_hat, x, n_x, y_hat, y, n_y){
  # Compute the Fisher’s Z-transformation of each data set
  #
  #Args:
  #  x_hat: z hat of trait X data set
  #  x:  trait X data set
  #  n_x: trait X data set (to get the length)
  #  y_hat: z hat of trait Y data set
  #  y: trait X data set
  #  n_y: trait Y data set (to get the length)
  #
  # Result: likelyhood
  dnorm(x_hat, x, sqrt(1/(length(n_x)-3)))*dnorm(y_hat, y, sqrt(1/(length(n_y)-3)))
};


yu.files <- list.files("/net/twins/home/fangchen/TEST",pattern="cut$",full.names=T, recursive= F); # the GWAS results of s/d 
my.files <- list.files("/net/twins/home/fangchen/LD-hub/cleaned_data/",pattern="clean$",full.names=T, recursive= F); # summaries data
for (yu in yu.files){
  trait.name <- strsplit(basename(yu), "[.]")[[1]][1]; # get the traits name
  yu.file <- read_tsv(yu);
  yu.sig.file <- paste0("/net/twins/home/fangchen/TEST/", trait.name, ".result.sig_May_12"); # already have the files which all SNP are significant 
  yu.sig <- read_tsv(yu.sig.file);
  mylist <- list();
  for (file in my.files){
    name <- basename(file);
    my.data <- read_tsv(file);
    X.set <- left_join(yu.sig, my.data, by=c("CHROM"="CHROM", "POS"="POS"), copy=F) %>% 
      select(CHROM, POS, SNP, X_x_beta=BETA, X_x_p=PVALUE, X_y_beta=Beta, X_y_p=pvalue) %>% na.omit; # the Trait X set(without info of Y)  
    Y.set <- filter(my.data, pvalue<=5*10^(-8)) %>% left_join(yu.file,  by=c("CHROM"="CHROM", "POS"="POS"), copy=F) %>% 
      select(CHROM, POS, SNP, Y_y_beta=Beta, Y_y_p=pvalue, Y_x_beta=BETA, Y_x_p=PVALUE) %>% na.omit; # the Trait Y set (without info of X) 
    X.set <- windowshift(X.set, 1000000, "X_x_p");  # the effect sizes on Trait X (without info of Y) 
    Y.set <- windowshift(Y.set, 1000000, "Y_y_p");  # the effect sizes on Trait Y (without info of X)
    beta.xx <- X.set$X_x_beta; beta.xy <- X.set$X_y_beta; 
    beta.yy <- Y.set$Y_y_beta; beta.yx <- Y.set$Y_x_beta;
    beta.xx <- na.omit(beta.xx);  beta.yy <- na.omit(beta.yy)
    if(length(beta.xx)<3||length(beta.yy)<3){
      result <- list(trait=trait.name, file=name, ratio=NA);
      mylist <- rbind(mylist, result); 
    }else {
    z_x_h <- z_trans(beta.xx, beta.xy);
    z_y_h <- z_trans(beta.yy, beta.yx);
    lh_M1 <- lh_M(z_x_h, z_x_h, beta.xx, z_y_h, 0, beta.yy);
    lh_M2 <- lh_M(z_x_h, 0, beta.xx, z_y_h, z_y_h, beta.yy);
    lh_M3 <- lh_M(z_x_h, 0, beta.xx, z_y_h, 0, beta.yy);
    AIC_score_M1 <- 2-2*log(lh_M1)
    AIC_score_M2 <- 2-2*log(lh_M2)
    AIC_score_M3 <- 2-2*log(lh_M3)
    ## get the mle of modle 4:
    lh_fun_M4 <- function(z){
      -1*dnorm(z_x_h, z, sqrt(1/(length(beta.xx)-3)))*dnorm(z_y_h, z, sqrt(1/(length(beta.yy)-3)))
    };
    test <- nlminb(0.1, lh_fun_M4, gradient = NULL);
    lh_M4 <- -1*lh_fun_M4(test$par);
    AIC_score_M4 <- 2-2*log(lh_M4);
    r <- exp((min(AIC_score_M1,AIC_score_M2)-min(AIC_score_M3,AIC_score_M4))/2);
    result <- list(trait=trait.name, file=name, ratio=r);
    mylist <- rbind(mylist, result);  ## append the result
    } ## end of else
  }
  my_res <- as.data.frame(mylist, row.names = F, stringsAsFactors=F);
  my_res <- as.data.frame(lapply(my_res,unlist));
  my_res$file <- gsub(".clean$", "", my_res$file);
  my_res_sig <- my_res[which(as.numeric(my_res$r)<0.01),]
  write.table(my_res, paste0(trait.name,".2way.full"), quote=F, row.names=F, sep="\t")
  write.table(my_res_sig, paste0(trait.name,".2way.sig"), quote=F, row.names=F, sep="\t")
}