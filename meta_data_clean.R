###########################
### use tidyverse package 
### install.packages("tidyverse")
###########################
library(tidyverse)

#okg:1000 genomes data
ori.form <- read_delim("MAGIC_ln_HOMA-IR.txt", "\t", col_types = "cccdddd") %>%
  mutate(E_allele=toupper(effect_allele), O_allele=toupper(other_allele)) %>%
  select(-c(2,3))

# add ALT AF info to dbsnp data
okg.data <- read_delim("all.af", "\t",col_type="ciccc")
dbsnp.dic <- read_delim("dic_re.vcf","\t",col_type="ciccc")%>%
  left_join(okg.data, by=c("CHROM" = "CHROM","POS" = "POS"),suffix = c("2","1")) %>%
  select(-c(4,5))

## N shoul be caculated from maf or effect AF in original form
whole.form <- left_join(ori.form,dbsnp.dic,by=c("snp" = "ID"),copy=F) %>%
  mutate(ALT1=strsplit(ALT1,","),AF=strsplit(AF,",")) %>% filter(!is.na(AF)) %>%  
   unnest(ALT1,AF,.drop=F) %>% filter(nchar(REF1)==1 & nchar(ALT1)==1) %>% mutate(AF=as.numeric(AF)) %>%
        mutate(Otg=paste0(REF1,ALT1), Samp=paste0(O_allele,E_allele)) %>% 
          mutate(N=ceiling(1/(stderr*2*maf*(1 - maf))))

#######################
## define function for senarios for flip, stand and flip_strand
#######################
flip <- function(x, y){
  f_x=y;
  f_y=x;
  return(paste0(f_x,f_y));
}

strand <- function(x, y) {
  f_x=ifelse(x=="A","T",ifelse(x=="T","A",ifelse(x=="G","C","G")));
  f_y=ifelse(y=="A","T",ifelse(y=="T","A",ifelse(y=="G","C","G")));
  return(paste0(f_x,f_y));
}

flip_strand <- function(x,y){
  xy=unlist(strsplit(flip(x,y),""));
  strand(xy[1],xy[2]);
}
#####################

Set_one <- whole.form %>% filter(Otg %in% c("AG","GA","CT","TC")) %>%
  mutate(Flip=flip(REF1,ALT1), Strand=strand(REF1,ALT1), Stran_Flip=flip_strand(REF1,ALT1)) %>%
    mutate(Beta=ifelse((Samp==Otg|Samp==Strand),effect,ifelse((Samp==Flip|Samp==Stran_Flip),-1*effect,NA))) %>%
           select(CHROM,POS,SNP=snp,REF=REF1,ALT=ALT1,Alt_Freq=AF,Beta,N,SE=stderr,pvalue)


#############
#### if there are effect_allele frequency in original data
#############

# Set_two <- whole.form %>% filter(Otg %in% c("AT","TA","GC","CG")) %>% filter(abs(Freq-0.5)>=0.1) %>% mutate(Flip=flip(REF1,ALT1)) %>% 
#   mutate(Diff=ifelse(abs(Freq-AF)>abs(Freq+AF-1),"B","S")) %>% filter(!is.na(Diff)) %>% 
#   mutate(BETA=ifelse(((Samp==Otg & Diff=="S")|(Samp==Flip & Diff=="S")),Beta,ifelse((Samp==Flip & Diff=="B")|(Samp==Otg & Diff=="B"), -1*Beta, NA))) %>% 
#   select(CHROM,POS,SNP=rsID,REF=REF1,ALT=ALT1,Alt_Freq=AF,Beta=BETA,N,SE,pvalue=P)

# new.form <- bind_rows(Set_one,Set_two, .id=NULL) %>% arrange(CHROM,POS)

#############

new.form <- Set_one %>% arrange(CHROM,POS)
             
write.table(new.form, "new_summary", sep = "\t", quote = F, row.names =F, col.names =T)


