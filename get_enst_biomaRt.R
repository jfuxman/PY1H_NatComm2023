library(dplyr) # Filtrar tablas data frames
library(ggplot2) # Realizar graficas
library(ggpubr) # Calculos estadisticos y graficas
library(stringr)
library(ShortRead)
library(ggridges)
library(Rsamtools)
library(seqinr)
library(phylotools) 
library(Biostrings)
library(RColorBrewer)
library(biomaRt)


df = readxl::read_xlsx('Het array with Ensembl Gene IDs.xlsx')
ensgs1 = df$`TF1 Ensembl Gene ID` %>% 
  na.omit() %>% 
  as.character() 

ensgs2 = df$`TF2 Ensembl Gene ID` %>% 
  na.omit() %>% 
  as.character()

ensgs = c(ensgs1, ensgs2) %>% unique()
'ENSG00000137504' %in% ensgs 

ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl")
biomartCacheClear()
biomartCacheInfo()
listAttributes(mart = ensembl)
listFilters(mart = ensembl)

df_ensts = getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id'),
                 filters = 'ensembl_gene_id',
                 values = ensgs, 
                 mart = ensembl,
                 useCache = FALSE)


df_dna = biomaRt::getSequence(id = df_ensts$ensembl_transcript_id,
                              type='ensembl_transcript_id',
                              seqType = 'coding',
                              mart= ensembl)

df_ensts$seq = df_dna$coding[match(df_ensts$ensembl_transcript_id, df_dna$ensembl_transcript_id)]
df_ensts[6,'seq'] == 'Sequence unavailable'

write.table(df_ensts, 'Gene_Transcript_Seq.txt', quote = F, row.names = F ,sep = '\t')
