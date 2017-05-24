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
  rho<- cor(x, y, method = "spearman", use = "pairwise.complete.obs");
  z <- 0.5 * (log(1+rho) - log(1-rho));
  output<-list(rho=rho, z=z)
  return(output)
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


yu.files <- list.files("/net/twins/home/fangchen/TEST",pattern="cut$",full.names=T, recursive= F);
my.files <- list.files("/net/twins/home/fangchen/LD-hub/cleaned_data/",pattern="clean$",full.names=T, recursive= F); 
for (yu in yu.files){
  trait.name <- strsplit(basename(yu), "[.]")[[1]][1];
  yu.file <- read_tsv(yu);
  yu.sig.file <- paste0(trait.name, ".result.sig_May_12");
  yu.sig <- read_tsv(yu.sig.file);
  mylist <- list();
  for (file in my.files){
    name <- basename(file);
    my.data <- read_tsv(file);
    X.set <- left_join(yu.sig, my.data, by=c("CHROM"="CHROM", "POS"="POS"), copy=F) %>% 
      select(CHROM, POS, SNP, X_x_beta=BETA, X_x_p=PVALUE, X_y_beta=Beta, X_y_p=pvalue) %>% na.omit;
    Y.set <- filter(my.data, pvalue<=5*10^(-8)) %>% left_join(yu.file,  by=c("CHROM"="CHROM", "POS"="POS"), copy=F) %>% 
      select(CHROM, POS, SNP, Y_y_beta=Beta, Y_y_p=pvalue, Y_x_beta=BETA, Y_x_p=PVALUE) %>% na.omit;
    X.set <- windowshift(X.set, 1000000, "X_x_p");
    Y.set <- windowshift(Y.set, 1000000, "Y_y_p");
    beta.xx <- X.set$X_x_beta; beta.xy <- X.set$X_y_beta;
    beta.yy <- Y.set$Y_y_beta; beta.yx <- Y.set$Y_x_beta;
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
  }
  my_res <- as.data.frame(mylist, row.names = F, stringsAsFactors=F);
  my_res <- as.data.frame(lapply(my_res,unlist));
  my_res$file <- gsub(".clean$", "", my_res$file);
  my_res_sig <- my_res[which(as.numeric(my_res$r)<0.01),]
  write.table(my_res, paste0(trait.name,".2way.full"), quote=F, row.names=F, sep="\t")
  write.table(my_res_sig, paste0(trait.name,".2way.sig"), quote=F, row.names=F, sep="\t")
}




cad <- read_tsv("/net/twins/home/fangchen/LD-hub/cleaned_data/cad.additive.Oct2015.clean")
dpw.cut <- read_tsv("/net/twins/home/fangchen/TEST/dpw.cut") %>% arrange(CHROM, POS)
dpw_sig <- read_tsv("/net/twins/home/fangchen/TEST/dpw.result.sig_May_12")
cad.set <- filter(cad, pvalue<=5*10^(-8)) %>% left_join(dpw.cut,  by=c("CHROM"="CHROM", "POS"="POS"), copy=F) %>%
  select(CHROM, POS, SNP, Y_cad_beta=Beta, Y_cad_p=pvalue, Y_dpw_beta=BETA, Y_dpw_p=PVALUE) %>% na.omit
dpw.set <- left_join(dpw_sig, cad, by=c("CHROM"="CHROM", "POS"="POS"), copy=F) %>%
  select(CHROM, POS, SNP, X_dpw_beta=BETA, X_dpw_p=PVALUE, X_cad_beta=Beta, X_cad_p=pvalue) %>% na.omit
cad.set <- windowshift(cad.set, 1000000, "Y_cad_p")
dpw.set <- windowshift(dpw.set, 1000000, "X_dpw_p")
z_x_h <- z_trans(dpw.set$X_dpw_beta, dpw.set$X_cad_beta)
z_y_h <- z_trans(cad.set$Y_cad_beta, cad.set$Y_dpw_beta)
x <- dpw.set$X_dpw_beta
y <- cad.set$Y_cad_beta

lh_M1 <- lh_M(z_x_h, z_x_h, x, z_y_h, 0, y);
lh_M2 <- lh_M(z_x_h, 0, x, z_y_h, z_y_h, y);
lh_M3 <- lh_M(z_x_h, 0, x, z_y_h, 0, y);
AIC_score_M1 <- 2-2*log(lh_M1)
AIC_score_M2 <- 2-2*log(lh_M2)
AIC_score_M3 <- 0-2*log(lh_M3)
lh_fun_M4 <- function(z, x_h, y_h, n_x, n_y){ 
  -1*dnorm(x_h, z, sqrt(1/(length(n_x)-3)))*dnorm(y_h, z, sqrt(1/(length(n_y)-3)))  
};
z_m4 <- nlminb(0.1, lh_fun_M4, x_h=z_x_h, y_h=z_y_h, n_x=x, n_y=y, gradient = NULL)$par;
lh_M4 <- lh_M(z_x_h, z_m4, x, z_y_h, z_m4, y);
AIC_score_M4 <- 2-2*log(lh_M4);
r <- exp((min(AIC_score_M1,AIC_score_M2)-min(AIC_score_M3,AIC_score_M4))/2);

test.beta <- lm(dpw.set$X_cad_beta ~ dpw.set$X_dpw_beta);
summary(test.beta)$coefficients


####
create sample

cad.clean <- cad[cad$SNP %in% cad.set$SNP, ]
dpw <- dpw.cut[dpw.cut$CHROM %in% dpw.set$CHROM & dpw.cut$POS %in% dpw.set$POS, ]




# trait.file <- "/Users/fangchen/My_Scripts/GWAS_summary_data_clean/GWAS_Data_Clean/GWAS_summary_clean/Sample/dpw.sampler.cut"
# trait.name <- strsplit(basename(trait.file), "[.]")[[1]][1]; # get the traits name
# yu.file <- read_tsv(trait.file);
# yu.sig <- read_tsv("/Users/fangchen/My_Scripts/GWAS_summary_data_clean/GWAS_Data_Clean/GWAS_summary_clean/Sample/dpw.result.sig_May_12");
# my.files <- "/Users/fangchen/My_Scripts/GWAS_summary_data_clean/GWAS_Data_Clean/GWAS_summary_clean/Sample/cad.sampler.clean";
# my_data <- read_tsv(my.files)

X.set <- read_tsv("~/My_Scripts/GWAS_summary_data_clean/GWAS_Data_Clean/GWAS_summary_clean/Sample/dpw.set")
Y.set <- read_tsv("~/My_Scripts/GWAS_summary_data_clean/GWAS_Data_Clean/GWAS_summary_clean/Sample/cad.set")


beta.xx <- as.numeric(X.set$X_dpw_beta); beta.xy <- as.numeric(X.set$X_cad_beta);
beta.yy <- as.numeric(Y.set$Y_cad_beta); beta.yx <- as.numeric(Y.set$Y_dpw_beta);
n_x <- length(beta.xx);  n_y <- length(beta.yy)
if(length(beta.xx)<10||length(beta.yy)<10){
  result <- list(trait="dpw", file="local", expo_beta=NA, expo_p=0, out_beta=NA, out_p=NA, t=NA, t_p=NA, df=NA);
}else {
  x <- lm(beta.xy ~ beta.xx); y <- lm(beta.yy ~ beta.yx);
  beta_x <- summary(x)$coefficients[2,1]; se_x <- summary(x)$coefficients[2,2]; p_x <- summary(x)$coefficients[2,4]
  beta_y <- summary(y)$coefficients[2,1]; se_y <- summary(y)$coefficients[2,2]; p_y <- summary(y)$coefficients[2,4]
  df <- (se_x^2/n_x + se_y^2/n_y)^2/((se_x^2/n_x)^2/(n_x - 1) + (se_y^2/n_y)^2/(n_y - 1));
  t <- abs(beta_x - beta_y)/sqrt(se_x^2 + se_y^2);
  p <- 2*pt(t, df, lower.tail=F);
  result <- list(trait="dpw", file="local", expo_beta=beta_x, expo_p=p_x, out_beta=beta_y, out_p=p_y, t=t, p=p, df=df);
}


my_res <- bi_regression(my.files, yu.file, yu.sig)
my_res <- na.omit(as.data.frame(my_res));
my_res$file <- gsub(".clean$", "", my_res$file);
my_res_sig <- my_res[which(as.numeric(unlist(my_res$expo_p))<0.01),]
write.table(format(my_res), paste0(trait.name,".biway_may24.full"), quote=F, row.names=F, sep="\t")
write.table(format(my_res_sig), paste0(trait.name,".biway_may24.sig"), quote=F, row.names=F, sep="\t")




setwd("/Users/fangchen/My_Scripts/GWAS_summary_data_clean/GWAS_Data_Clean/GWAS_summary_clean")
# Using the extracting data ----------------------------------------------------------

bi_regression <- function(X.file){ 
  # Calculate the ratio in each file
  #
  # Args:
  #  file: the summaries file
  #
  # Result:
  #  a list consist of trait name(exposure), file name(outcome) and the r ratio
  filename <- basename(X.file)
  # m <- regexpr("^[^_]+(?=_)", filename, perl = T)
  # trait.name <- regmatches(filename, m)
  trait.name <- strsplit(filename, split = "_")[[1]][1]
  Y.name <- gsub("^[^_]+(?=_)_","",filename,perl=TRUE)
  Y.file <- paste(Y.name, trait.name, sep = "_")
  X.set <- read_tsv(X.file);
  Y.set <- read_tsv(Y.file);
  beta.xx <- as.numeric(X.set$X_x_beta); beta.xy <- as.numeric(X.set$X_y_beta);
  rank.xx <- rank(beta.xx); rank.xy <- rank(beta.xy);
  beta.yy <- as.numeric(Y.set$Y_y_beta); beta.yx <- as.numeric(Y.set$Y_x_beta);
  rank.yy <- rank(beta.yy); rank.yx <- rank(beta.yx)
  n_x <- length(beta.xx);  n_y <- length(beta.yy)
  if (n_x == 0){
    result <- list(trait=trait.name, file=Y.name, expo_rho=NA, expo_p=NA, expo_se=NA, N_x=n_x, out_rho=NA, out_p=NA, out_se=NA, N_y=n_y, t=NA, t_p=NA, df=NA, AIC_r=NA, lm_beta=NA, lm_p=NA);
  } else if (length(beta.xx)<10 || length(beta.yy)<10){
    # calc the rank rho
    x <- cor.test(rank.xy, rank.xx);
    rho.x <- unname(x$estimate);
    se_x <- unname(sqrt((1 - x$estimate^2)/x$parameter)); p_x <- x$p.value;
    # calc the lm beta
    test.beta <- lm(beta.xy ~ beta.xx);
    beta=summary(test.beta)$coefficients[2,1]; pvalue=summary(test.beta)$coefficients[2,4];
    result <- list(trait=trait.name, file=Y.name, expo_rho=rho.x, expo_p=p_x, expo_se=se_x, N_x=n_x, out_rho=NA, out_p=NA, out_se=NA, N_y=n_y, t=NA, t_p=NA, df=NA, AIC_r=NA, lm_beta=beta, lm_p=pvalue);
  } else {
    #   x <- lm(beta.xy ~ beta.xx); y <- lm(beta.yy ~ beta.yx);
    #   beta_x <- summary(x)$coefficients[2,1]; se_x <- summary(x)$coefficients[2,2]; p_x <- summary(x)$coefficients[2,4]
    #   beta_y <- summary(y)$coefficients[2,1]; se_y <- summary(y)$coefficients[2,2]; p_y <- summary(y)$coefficients[2,4]
    x <- cor.test(rank.xy, rank.xx); y <- cor.test(rank.yy, rank.yx);
    rho.x <- unname(x$estimate); rho.y <- unname(y$estimate);
    se_x <- unname(sqrt((1 - x$estimate^2)/x$parameter)); p_x <- x$p.value;
    se_y <- unname(sqrt((1 - y$estimate^2)/y$parameter)); p_y <- y$p.value;
    df <- (se_x^2/n_x + se_y^2/n_y)^2/((se_x^2/n_x)^2/(n_x - 1) + (se_y^2/n_y)^2/(n_y - 1));
    t <- abs(rho.x - rho.y)/sqrt(se_x^2 + se_y^2);
    p <- 2*pt(t, df, lower.tail=F);
    test.beta <- lm(beta.xy ~ beta.xx);
    beta=summary(test.beta)$coefficients[2,1]; pvalue=summary(test.beta)$coefficients[2,4];
    ### calculate AIC-r
    z_x_h <- z_trans(beta.xx, beta.xy)$z; x_rho <- z_trans(beta.xx, beta.xy)$rho;
    z_y_h <- z_trans(beta.yy, beta.yx)$z; y_rho <- z_trans(beta.yy, beta.yx)$rho;
    lh_M1 <- lh_M(z_x_h, z_x_h, beta.xx, z_y_h, 0, beta.yy);
    lh_M2 <- lh_M(z_x_h, 0, beta.xx, z_y_h, z_y_h, beta.yy);
    lh_M3 <- lh_M(z_x_h, 0, beta.xx, z_y_h, 0, beta.yy);
    AIC_score_M1 <- 2-2*log(lh_M1);
    AIC_score_M2 <- 2-2*log(lh_M2);
    AIC_score_M3 <- 0-2*log(lh_M3);
    lh_fun_M4 <- function(z, x_h, y_h, n.x, n.y){
      -1*dnorm(x_h, z, sqrt(1/(length(n.x)-3)))*dnorm(y_h, z, sqrt(1/(length(n.y)-3)))
    };
    z_m4 <- nlminb(0.1, lh_fun_M4, x_h=z_x_h, y_h=z_y_h, n.x=beta.xx, n.y=beta.yy, gradient = NULL)$par;
    lh_M4 <- lh_M(z_x_h, z_m4, beta.xx, z_y_h, z_m4, beta.yy);
    AIC_score_M4 <- 2-2*log(lh_M4);
    r <- exp((min(AIC_score_M1,AIC_score_M2)-min(AIC_score_M3,AIC_score_M4))/2);
    result <- list(trait=trait.name, file=Y.name, expo_rho=rho.x, expo_p=p_x, expo_se=se_x, N_x=n_x, out_rho=rho.y, out_p=p_y, out_se=se_y, N_y=n_y, t=t, p=p, df=df, AIC_r=r, lm_beta=beta, lm_p=pvalue);
  }
}

# f1 <- "~/TEST/MR/Data_sets/dnd_set/dnd_body_fat_percentage_GWAS"
# bi_regression(f1)
set.files <- list.files(".",pattern="^sc",full.names=T, recursive= F);
my_res <- do.call(rbind, lapply(set.files, bi_regression));
my_res <- as.data.frame(apply(my_res, 2, unlist), stringsAsFactors = F); 
my_res[, c(3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15)] <- sapply(my_res[, c(3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15)], as.numeric)
write_excel_csv(format(my_res, digits=3), paste0("/Users/fangchen/TEST/MR/excel/", "sc", ".xls"), na = "NA")


#my_res_sig <-  my_res[which(as.numeric(my_res$expo_p)<0.05),]
# my_res <- my_res[!is.na(my_res$expo_rho),]
# my_res <- as_tibble(my_res,validate = TRUE);
# my_res <- filter(my_res, !is.na(expo_rho))








# plot test--------------------------------------------------------------------
Y.set <- read_tsv("~/TEST/MR/Data_sets/dpw_set/fn2stu_dpw")
X.set <- read_tsv("~/TEST/MR/Data_sets/dpw_set/dpw_fn2stu")
beta.xx <- as.numeric(X.set$X_x_beta); beta.xy <- as.numeric(X.set$X_y_beta);
beta.yy <- as.numeric(Y.set$Y_y_beta); beta.yx <- as.numeric(Y.set$Y_x_beta);
plot(beta.yx, beta.yy)
plot(beta.xx, beta.xy)
abline(lm(beta.xy ~ beta.xx))
beta.yx


# plot --------------------------------------------------------------------

library(grid)
library(gridExtra)

get_plot_pdf <- function(X.file){
  filename <- basename(X.file)
  trait.name <- strsplit(filename, split = "_")[[1]][1]
  Y.name <- gsub("^[^_]+(?=_)_","",filename,perl=TRUE)
  Y.file <- paste(Y.name, trait.name, sep = "_")
  X.set <- read_tsv(X.file);
  Y.set <- read_tsv(Y.file);
  # Y.set <- read_tsv("~/TEST/MR/Data_sets/cpd_set/Menarche_Nature2014_GWASMetaResults_17122014_cpd") 
  # X.set <- read_tsv("~/TEST/MR/Data_sets/cpd_set/cpd_Menarche_Nature2014_GWASMetaResults_17122014") 
  # beta.xx <- as.numeric(X.set$X_x_beta); beta.xy <- as.numeric(X.set$X_y_beta);
  # beta.yy <- as.numeric(Y.set$Y_y_beta); beta.yx <- as.numeric(Y.set$Y_x_beta);
  # plot(beta.yx, beta.yy)
  # plot(beta.xx, beta.xy)
  # abline(lm(beta.xy ~ beta.xx))
  p1 <- ggplot(X.set, aes(x=as.numeric(X.set$X_x_beta), y=as.numeric(X.set$X_y_beta))) + geom_point() + geom_smooth(method = 'lm')
  p2 <- ggplot(Y.set, aes(x=as.numeric(Y.set$Y_x_beta), y=as.numeric(Y.set$Y_y_beta))) + geom_point() + geom_smooth(method = 'lm')
  p3 <- try(grid.arrange(p1, p2, ncol=2, nrow =1, top = paste0(trait.name, "_", Y.name)))
  try(ggsave(paste0(trait.name, "_", Y.name,".pdf"),  plot=p3, device = "pdf"));
}

set.files <- list.files(".",pattern="^si",full.names=T, recursive= F);
sapply(set.files, get_plot_pdf)
