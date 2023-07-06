library(dplyr) # Filtrar tablas data frames
library(ggplot2) # Realizar graficas
#library(ggpubr) # Calculos estadisticos y graficas
library(stringr)
library(ShortRead)
library(ggridges)
library(Rsamtools)
library(seqinr)


## Data for January
#setwd('D:/BostonLab/Anna_heterodimer_sequencing/new_data/221003_JFB10247_fastq/bwa_outputs')
setwd('bwa_outputs')
files = list.files(pattern = '.bam')
files = files[!str_detect(files, 'sorted')]

MetaDf = readxl::read_excel('../vtr-ht_changed.xlsx')
colnames(MetaDf) = c('Info',' TF1', 'Ensg1', 'Ensg2', 'Plate', 'Well')


PlateIndexs = c()
WellIndexs = c()
PlasmidNames = c()
AlignedReads = c()
CorrectReads = c()
TrimmedReads = c()
max_n_reads_names = c()
max_n_reads_s = c()
expectedSeqs = c()

for(bam_file in files){
  print(bam_file)
  
  param <- ScanBamParam(flag = scanBamFlag(),
                        tag = c('NM', 'MD', 'AS','XS', 'XA'),
                        what = scanBamWhat())

  
  TrueAln <- scanBam(bam_file, param = param)[[1]]
  
  
  TrueDf = data.frame(read_name = TrueAln$qname,
                      read_flag = TrueAln$flag,
                      read_seq = TrueAln$seq,
                      ref_name = TrueAln$rname,
                      read_length = TrueAln$qwidth,
                      #read_distance = TrueAln$tag$NM,
                      read_position = TrueAln$pos,
                      read_cigar = TrueAln$cigar,
                      read_primary_score = TrueAln$tag$AS,
                      read_sec_score = TrueAln$tag$XS)
  if (nrow(TrueDf)>1) {

    metavalues = TrueDf %>% filter(!is.na(ref_name) & read_primary_score/read_length >= 0.9) %>%
    group_by(read_name, ref_name) %>% 
    summarize(n = n()) %>% 
    arrange(desc(n)) %>%
    mutate(ref_gene = str_split_fixed(ref_name, '_', 2)[, 1]) %>% 
    ungroup() %>% 
    group_by(read_name, ref_gene) %>% 
    summarize(n = n()) %>% 
    arrange(ref_gene) %>% 
    ungroup() %>% 
    group_by(ref_gene) %>% 
    summarize(n = n()) %>% 
    arrange(desc(n))
  
   } else {
    metavalues = TrueDf %>% filter(!is.na(ref_name)) %>%
    group_by(read_name, ref_name) %>% 
    summarize(n = n()) %>% 
    arrange(desc(n)) %>%
    mutate(ref_gene = str_split_fixed(ref_name, '_', 2)[, 1]) %>% 
    ungroup() %>% 
    group_by(read_name, ref_gene) %>% 
    summarize(n = n()) %>% 
    arrange(ref_gene) %>% 
    ungroup() %>% 
    group_by(ref_gene) %>% 
    summarize(n = n()) %>% 
    arrange(desc(n))
  
   }
  name = str_split_fixed(bam_file, '\\.', 2)[,1]
  PlateIndex = str_split_fixed(name, '_', 3)[,1]
  WellIndex = str_split_fixed(name, '_', 3)[,2]
  PlasmidName = str_split_fixed(name, '_', 3)[,3]
  
  n_TrimmedRead = readFastq(paste('../',PlateIndex, '_clean_files_relaxed_prev/', PlateIndex, '_ROI_', WellIndex, '_', PlasmidName, '_F_trimmed.fastq', sep = ''))
  n_TrimmedRead = n_TrimmedRead@sread %>%  as.data.frame()
  n_TrimmedRead = nrow(n_TrimmedRead)
  
  if(PlasmidName == 'AD2U'){
    ensg_name = MetaDf %>% filter(Plate == PlateIndex, Well == WellIndex) %>% dplyr::select(Ensg2) %>% unlist()
    longitud = MetaDf %>% filter(Plate == PlateIndex, Well == WellIndex) %>% dplyr::select(Ensg2) %>% nrow()
  }else{
    ensg_name = MetaDf %>% filter(Plate == PlateIndex, Well == WellIndex) %>% dplyr::select(Ensg1) %>% unlist()
    longitud = MetaDf %>% filter(Plate == PlateIndex, Well == WellIndex) %>% dplyr::select(Ensg1) %>% nrow()
  }
   print(ensg_name)
   print(longitud)
   if(longitud == 0){
    ensg_name = 'NA'
    n_aligned_reads = sum(metavalues$n)
    n_correct_reads = 0
  } else if (is.na(ensg_name)) {
    ensg_name = 'NA'
    n_aligned_reads = sum(metavalues$n)
    n_correct_reads = 0
  } else if (ensg_name == 'empty'){
    ensg_name = 'NA'
    n_aligned_reads = sum(metavalues$n)
    n_correct_reads = 0    
  }
    else{
    n_aligned_reads = sum(metavalues$n)
    n_correct_reads = metavalues %>% filter(ref_gene == ensg_name)
    if(nrow(n_correct_reads) == 0) {
      n_correct_reads = 0
    } else{
      n_correct_reads = n_correct_reads %>% select(n) %>% unlist() %>% unname() 
    }
  }
  max_n_reads = metavalues[1,2] %>% unlist() %>% unname()
  max_n_reads_name = metavalues[1,1] %>% unlist() %>% unname()
  
  
  
  PlateIndexs = c(PlateIndexs, PlateIndex)
  WellIndexs = c(WellIndexs, WellIndex)
  PlasmidNames = c(PlasmidNames, PlasmidName)
  TrimmedReads = c(TrimmedReads, n_TrimmedRead)
  AlignedReads = c(AlignedReads, n_aligned_reads)
  CorrectReads = c(CorrectReads, n_correct_reads)
  max_n_reads_names = c(max_n_reads_names, max_n_reads_name)
  max_n_reads_s = c(max_n_reads_s, max_n_reads)
  expectedSeqs = c(expectedSeqs, ensg_name)
  }

results = data.frame(PlasmidNames, PlateIndexs, WellIndexs,TrimmedReads,AlignedReads,
                     CorrectReads, expectedSeqs, max_n_reads_s, max_n_reads_names)

results$PlateIndexs = as.numeric(results$PlateIndexs)
results$expected = ifelse(results$expectedSeqs == results$max_n_reads_names, 'Yes', ifelse(results$expectedSeqs =='NA', NA, 'No'))
results$letter = str_sub(results$WellIndexs, 1 ,1)
results$number = str_sub(results$WellIndexs, 2 ,-1)

for (Plate in unique(results$PlateIndexs)){
  tmp_results = results %>% filter(PlateIndexs == Plate)
  P1 = ggplot(tmp_results, aes(x = PlasmidNames, y = number, fill = log10(CorrectReads+1))) + 
    geom_tile(color = 'black')+
    theme_bw()+
    scale_fill_gradientn(colours = c('white', 'blue'))+
    #geom_text(aes(label = paste(expected , CorrectReads, sep = '-'))) + 
    geom_text(aes(label = paste(expected , round(CorrectReads/AlignedReads * 100,1), sep = '/'))) +  
    theme(legend.position = 'top') + facet_wrap(~letter) +
    ggtitle(Plate)
  
  ggsave(filename = paste(Plate, '.jpg', sep =''), plot = P1, height = 300, width = 250, units = 'mm')
}
getwd()
write.table(results, file = 'PGADT7_AD2U_results.txt', sep = '\t', row.names = F, quote = F)

