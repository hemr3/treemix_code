library(tidyverse)



add.col<-function(df, new.col) {n.row<-dim(df)[1]
length(new.col)<-n.row
cbind(df, new.col)


}

#merging 2 columns into one in a data frame
#below are the steps to make .snp data 'comprehensible' for treemix

d4_10_col = d4_10 %>% 
  subset(select=c(snp_id, tref, talt)) #subsetting out the columns i want

d4_10_un = d4_10_col %>%  
  transform(d4_10_col, deni4=paste(tref, talt, sep = ",")) #adding these columns together
  
d4_10_un = d4_10_un %>% 
  subset(select=c(snp_id, deni4)) #removing the tref and talt columns


d8_10_col = d8_10 %>% 
  subset(select=c(snp_id, tref, talt))

d8_10_un = d8_10_col %>% 
  transform(d8_10_col, deni8=paste(tref, talt, sep = ","))

d8_10_un = d8_10_un %>% 
  subset(select=c(snp_id, deni8))


chag10_col = chag10 %>% 
  subset(select=c(snp_id, tref, talt))

chag10_un = chag10_col %>% 
  transform(chag10_col, chag=paste(tref, talt, sep = ","))

chag10_un = chag10_un %>% 
  subset(select=c(snp_id, chag))


altai10_col = altai10 %>% 
  subset(select=c(snp_id, tref, talt))

altai10_un = altai10_col %>% 
  transform(altai10_col, altai=paste(tref, talt, sep = ","))

altai10_un = altai10_un %>% 
  subset(select=c(snp_id, altai))


vind10_col = vind10 %>% 
  subset(select=c(snp_id, tref, talt))

vind10_un = vind10_col %>% 
  transform(vind10_col, vind=paste(tref, talt, sep = ","))

vind10_un = vind10_un %>% 
  subset(select=c(snp_id, vind))


mes10_col = mes10 %>% 
  subset(select=c(snp_id, tref, talt))

mes10_un = mes10_col %>% 
  transform(mes10_col, mez1=paste(tref, talt, sep = ","))

mes10_un = mes10_un %>% 
  subset(select=c(snp_id, mez1))

#binding all the frames together; its super clunky but it works. 
#binding by snp_id here, but it might be more useful with 'pos' in the entire CHR length

comb_10 = full_join(d4_10_un, d8_10_un, by = "snp_id")

comb_10 = full_join(comb_10, altai10_un, by = "snp_id")

comb_10 = full_join(comb_10, chag10_un, by = "snp_id")

comb_10 = full_join(comb_10, mes10_un, by = "snp_id")

comb_10 = full_join(comb_10, vind10_un, by = "snp_id")

comb_10 = comb_10 %>% 
  subset(select=c(deni4, deni8, altai, chag, mez1, vind))

#replacing the NAs with 0,0 

comb_10[is.na(comb_10)] <- "0,0"

#writing out this to a snp file; needs to be compressed as a gzip
write_delim(comb_10, "total_chr10.snp")
