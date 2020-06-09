#pkg and wd pre
library(readxl)
library(tidyverse)
library(HGNChelper)
library(org.Hs.eg.db);library(org.Bt.eg.db)
setwd('/Users/liulihe95/Desktop/BrayanEnrich/snp_files')

# function pre - check validity of suspicious symbols in HGNC
check_symbol = function(notsure){
  tmp = checkGeneSymbols(notsure) %>% 
    mutate(Suggested.Symbol.Merge = ifelse(is.na(Suggested.Symbol),x,Suggested.Symbol))
  out = unlist(tmp$Suggested.Symbol.Merge,use.names = F)
  return(out)
}

# retrive bta gene information
Universe_id_bta = AnnotationDbi::select(org.Bt.eg.db, 
                                        as.character(AnnotationDbi::keys(org.Bt.eg.db,keytype = c("ENSEMBL"))),
                                        columns = c("ENTREZID","SYMBOL"),
                                        keytype = "ENSEMBL") %>% 
  dplyr::distinct(ENSEMBL,.keep_all= TRUE)

# Extract Gene and match identifiers - (ENS - Entrez - HGNC corrected symbol)
setwd('/Users/liulihe95/Desktop/BrayanEnrich/snp_files/')
GeneID_sig = c('sig_snps_L1.xlsx','sig_snps_L2.xlsx','sig_snps_L3.xlsx')
GeneID_all = c('total_snps_L1.xlsx','total_snps_L2.xlsx','total_snps_L3.xlsx')
sigGene_All = list()
for (i in seq_along(GeneFile_ID)){
  siglist = read_xlsx(GeneID_sig[i],col_names = 'Gene') %>% 
    unlist(use.names = F)
  sigGene_Each_Lac = read_xlsx(GeneID_all[i],col_names = 'Gene') %>% 
    dplyr::left_join(Universe_id_bta,by=c('Gene'='ENSEMBL')) %>% 
    mutate(SYMBOL_Suggested = check_symbol(SYMBOL)) %>% 
    mutate(Sig = ifelse(Gene %in% siglist,'1','0'))
  sigGene_All[[i]] = sigGene_Each_Lac
  names(sigGene_All)[i] = paste0('Lac',i)
}

# output examained identifiers
setwd('/Users/liulihe95/Desktop/BrayanEnrich')
save(sigGene_All,file = 'sigGene_All.rda')


