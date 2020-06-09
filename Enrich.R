# read in 
setwd('/Users/liulihe95/Desktop/BrayanEnrich')
sigGene_All = load('sigGene_All.rda')


#
HyperGEnrich = function(GeneSet,
                        Database = '', #'go','kegg,'interpro','mesh','msig','reactome'
                        minOverlap = 4,
                        pvalue_thres = 0.05,
                        adj_pvalue_thres = 1,
                        padj_method = "BH"#"bonferroni" or "hochberg" or "BH"
                        ){
  library(ggplot2);library(biomaRt);library(gage);library(magrittr)# load pkg
  # get db
  TestingSubsetNames = names(GeneSet)
  message("Total Number of subsets/module to check: ",length(TestingSubsetNames))
  if (Database == 'go'){
    DB_List = get(load('../testDB/GO_DB.rda'));IDtype = 1
  } else if (Database == 'kegg'){
    DB_List = get(load('../testDB/KEGG_DB.rda'));IDtype = 2
  } else if (Database == 'interpro'){
      DB_List = get(load('../testDB/Interpro_DB.rda'));IDtype = 1
  } else if (Database == 'mesh'){
        DB_List = get(load('../testDB/MeSH_DB.rda'));IDtype = 2
  } else if (Database == 'msig'){
          DB_List = get(load('../testDB/Msig_DB.rda'));IDtype = 2
  } else if (Database == 'reactome'){
            DB_List = get(load('../testDB/Reactome_DB.rda'));IDtype = 2
  } else {
              print('Please choose database: go,kegg,interpro,mesh,msig,reactome')
  }
  #
  GeneInDB = unique(unlist(DB_List,use.names = F))
  message('Database selected: ', Database,'.',
          '\nTotal number of terms to check ', length(DB_List),'.',
          '\nTotal number of genes showed up in any record ',length(GeneInDB),'.')
  #
  total_enrich = 0
  raw_pvalue_all = numeric()
  results = list()
  results_raw = list()
  
  #
  for (i in c(1:(length(TestingSubsetNames)))){
    Genelist_tmp = GeneSet[[i]]
    message("working on dataset #",i," - ",TestingSubsetNames[i])
    # get sig gene
    sig.genes_tmp = Genelist_tmp %>% data.frame() %>% 
      dplyr::filter(Sig == '1') %>% 
      dplyr::select(names(GeneSet[[i]])[IDtype]) %>% 
      unlist(use.names = F)
    sig.genes = sig.genes_tmp %>% na.omit()
    # get total gene
    total.genes_tmp = Genelist_tmp %>% data.frame() %>% 
      dplyr::filter(Sig %in% c('1','0')) %>% 
      dplyr::select(names(GeneSet[[i]])[IDtype]) %>% 
      unlist(use.names = F)
    total.genes = total.genes_tmp %>% na.omit()
    
    message('Total genes: ',length(total.genes),'. with ',length(total.genes_tmp)-length(total.genes),' loss due to identifier match',
            '\nSig Genes: ', length(sig.genes),'. with ',length(sig.genes_tmp)-length(sig.genes) ,' loss due to identifier match')
    
    # overlap with the database
    N = length(total.genes[total.genes %in% GeneInDB])
    S = length(sig.genes[sig.genes %in% GeneInDB]) #
    
    ExternalLoss_total = paste((length(total.genes) - N),round((length(total.genes) - N)/N,3),sep = "/")
    ExternalLoss_sig = paste((length(sig.genes) - S),round((length(sig.genes) - S)/S,3),sep = "/")
    # formatting out put
    out = data.frame(Term=character(),
                     totalG=numeric(),
                     sigG=numeric(),
                     Pvalue=numeric(),
                     ExternalLoss_total = character(),
                     ExternalLoss_sig = character(),
                     findG =  character())
    #
    for(j in 1:length(DB_List)){
      if (j%%1000 == 0) {message("tryingd on term ",j," - ",names(DB_List)[j])}
      gENEs = DB_List[[j]] # all gene in target GO #### note
      m = length(total.genes[total.genes %in% gENEs]) # genes from target GO and in our dataset
      findG = sig.genes[sig.genes %in% gENEs]
      s = length(findG)
      PastefindG = paste(findG,collapse="/")
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value,100)
      tmp = data.frame(Term= names(DB_List)[j],
                       totalG = m, 
                       sigG = s, 
                       pvalue = Pval, 
                       ExternalLoss_total = ExternalLoss_total,
                       ExternalLoss_sig = ExternalLoss_sig,
                       findG = PastefindG)
      out = rbind(out,tmp)
    }
     names(out) = c('Term','totalG','sigG','pvalue','ExternalLoss_total','ExternalLoss_sig','findG')
    # put all palues in a box
    #raw_pvalue_all = append(raw_pvalue_all,out$Pvalue,length(raw_pvalue_all))
    
    
    # raw
     "/-" <- function(x,y) ifelse(y==0,0,base:::"/"(x,y))
    final_raw = out %>% 
      arrange(pvalue) %>% 
      dplyr::mutate(hitsPerc = sigG*100 /- totalG) %>% 
      mutate(adj.pvalue = p.adjust(pvalue, method = padj_method))
    
    results_raw[[i]] = final_raw
    names(results_raw)[i] = paste(TestingSubsetNames[i],"with",dim(final_raw)[1],"terms tested - raw")
  
    # selection 
    final = final_raw %>% 
      dplyr::filter(totalG >= minOverlap) %>% 
      dplyr::filter(pvalue <= pvalue_thres) %>% 
      dplyr::filter(`adj.pvalue` <= adj_pvalue_thres)
    
    results[[i]] = final
    names(results)[i] = paste(TestingSubsetNames[i],"with",dim(final)[1],"enriched terms - selected")
    
    # selection ends
    message("Significant Enrichment Hits: ",nrow(final))
    total_enrich = total_enrich + nrow(final)
  }
  # raw_pvalue_index = seq(0.05,1,by=0.05)
  # raw_pvalue_sum = numeric()
  # for( z in seq_along(raw_pvalue_index)){raw_pvalue_sum[z] = length(which(raw_pvalue_all <= raw_pvalue_index[z]))}
  # raw_pvalue_distribution = data.frame(index = raw_pvalue_index,counts_GO = raw_pvalue_sum)
  #raw_pvalue_distribution
  
  save(results_raw,results,
       file = paste0(Database,'-Enrichment-',minOverlap,'-',pvalue.thres,'-',adj.pvalue.thre,'.rda'))
  
  message(total_enrich," significant terms found within ",
          length(TestingSubsetNames)," modules/subsets ", 
          "with pvalue and adj.pvalue set at ",pvalue.thres,' and ',adj.pvalue.thre)
}


HyperGEnrich(GeneSet = sigGene_All,
             Database = 'go',
             #Database = '', #'go','kegg,'interpro','mesh','msig','reactome'
             minOverlap = 4,
             pvalue_thres = 0.05,
             adj_pvalue_thres = 1,
             padj_method = "BH")

HyperGEnrich(GeneSet = sigGene_All,
             Database = 'kegg',
             #Database = '', #'go','kegg,'interpro','mesh','msig','reactome'
             minOverlap = 4,
             pvalue_thres = 0.05,
             adj_pvalue_thres = 1,
             padj_method = "BH")

HyperGEnrich(GeneSet = sigGene_All,
             Database = 'interpro',
             #Database = '', #'go','kegg,'interpro','mesh','msig','reactome'
             minOverlap = 4,
             pvalue_thres = 0.05,
             adj_pvalue_thres = 1,
             padj_method = "BH")

HyperGEnrich(GeneSet = sigGene_All,
             Database = 'mesh',
             #Database = '', #'go','kegg,'interpro','mesh','msig','reactome'
             minOverlap = 4,
             pvalue_thres = 0.05,
             adj_pvalue_thres = 1,
             padj_method = "BH")

HyperGEnrich(GeneSet = sigGene_All,
             Database = 'reactome',
             #Database = '', #'go','kegg,'interpro','mesh','msig','reactome'
             minOverlap = 4,
             pvalue_thres = 0.05,
             adj_pvalue_thres = 1,
             padj_method = "BH")

HyperGEnrich(GeneSet = sigGene_All,
             Database = 'msig',
             #Database = '', #'go','kegg,'interpro','mesh','msig','reactome'
             minOverlap = 4,
             pvalue_thres = 0.05,
             adj_pvalue_thres = 1,
             padj_method = "BH")
