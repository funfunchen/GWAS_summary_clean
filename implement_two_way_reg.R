
z_trans <- function(x, y){
  rho<- cor(x, y, method = "spearman");
  z <- 0.5 * (log(1+rho) - log(1-rho));
}
         
lh_M1 <- function(z){
 -1*dnorm(z_x_h, z, sqrt(length(x)-3))*dnorm(z_y_h, 0, sqrt(length(y)-3))
}

lh_M2 <- function(z){
  -1*dnorm(z_x_h, 0, sqrt(length(x)-3))*dnorm(z_y_h, z, sqrt(length(y)-3))
}

lh_M4 <- function(z){
  -1*dnorm(z_x_h, z, sqrt(length(x)-3))*dnorm(z_y_h, z, sqrt(length(y)-3))
}

z_x_h <- z_trans(x, y)
z_y_h <- z_trans(y, x)

test <- nlminb(0.1, lh_M1, gradient = NULL)
test$par
AIC_score <- 2-2*log(test$par)
AIC_score


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