library(dplyr)
library(stringr)
library(TFMPvalue, lib = 'my_R_packages')
library(readxl)
library(matrixStats)
library(stringi)
library(Biostrings)


### Run as follows
## Rscript MotifAnalysis.R 12 1

args = commandArgs(trailingOnly = TRUE)
index.promoter = as.numeric(args[1]) # Index of the promoter (1-15)
n_elements = as.numeric(args[2]) 

dna_subsequences_fast <- function(dna_string, X) {
  n <- nchar(dna_string)
  dna_list <- stri_sub(dna_string, 1:(n-X+1), length=X)
  return(dna_list)
}

get_score <- function(dna_string, input.matrix1, input.matrix2) {
  # Split the DNA sequence into a list of individual nucleotides
  dna_list <- strsplit(dna_string, "")[[1]]
  
  # Transpose the PWM matrix and convert it to log odds form
  log_matrix <- log(input.matrix1 / 0.25)
  
  # Look up the corresponding log odds values for each nucleotide
  log_odds_values <- log_matrix[cbind(seq_along(dna_list), match(dna_list, c("A", "C", "G", "T")))]
  
  # Calculate the sum of log odds values using the colSums function from matrixStats
  log_sum <- sum(log_odds_values, na.rm =  TRUE)
  
  # pvalue
  pvalue = TFMsc2pv(mat = input.matrix2, score = log_sum, bg = bg, type = 'PWM')
  
  return(pvalue)
}



bg <- c(A=0.25, C=0.25, G=0.25, T=0.25)

# Events
events_df = read_excel('Final events list 6Feb23.xlsx')
tfs_events = c(events_df$TF1,events_df$TF2) %>% unique()

# Motif - TF relation table
df1 = read.table('TF_Information.txt',header = T, sep = '\t', fill = T) %>% as_tibble()

# Information content of TFs
df3 = read.table('Information_content.txt', header = T, sep = '\t', fill = T) %>% as_tibble()

# Useless motifs
useless_motifs = df3 %>% filter(is.na(start) & is.na(end)) %>% dplyr::select(ID) %>% unlist()

# Useful motifs
df1 = df1 %>% filter(!Motif_ID %in% useless_motifs)

# Remove empty motifs
df1 = df1 %>% filter(Motif_ID != '.')

# Remove motifs not beloning to usefuls TFS
df1 = df1 %>% filter(TF_Name %in% tfs_events)


# Promoter sequence
promoters <- readDNAStringSet("promoters.fasta")
IDs = promoters@ranges@NAMES
sequences = as.data.frame(promoters)
promoters.df = data.frame(ID = IDs, sequence = sequences)

all_positions = c()
for(i in index.promoter){
  id = promoters.df[i,1] %>% unlist()
  promoter = promoters.df[i,2] %>% unlist()
  
  for (j in c(n_elements:(n_elements+100))){
    
    print(paste(id, j))
    
    # Get file.name
    file.name = df1[j,'Motif_ID'] %>% unlist()

    # Read matrix
    if(file.name == 'NA' | is.na(file.name)){next}
    tmp.df = read.table(paste0('pwms_all_motifs/', file.name, '.txt'), header = T)

    # Remove first position column
    tmp.df = tmp.df %>% dplyr::select(-1) 

    # Define the PWM matrix
    pwm_matrix = as.matrix(tmp.df)

    # Save another matrix to get p-value
    row_pwm_matrix = t(pwm_matrix)

    # Remove names
    dimnames(pwm_matrix) = NULL

    # Split length
    length_seq = nrow(tmp.df)

    if(length_seq > 15){
      all_positions = c(all_positions, NA)
    } else{
      # Split in small sequences
      dna_list = dna_subsequences_fast(promoter, length_seq) 

      # Get p-values
      pvalues = lapply(dna_list, get_score, input.matrix1 = pwm_matrix ,input.matrix2 = row_pwm_matrix)

      # Get significant positions
      sig.positions = which(pvalues <= 0.0001)
      sig.positions = paste0(sig.positions, collapse =',')

      all_positions = c(all_positions, sig.positions)
    }

    }
}

all_positions = data.frame(all_positions)

write.table(all_positions, paste0('TFPM_motifs/all_positions_', index.promoter, '_', n_elements , '.txt'))

