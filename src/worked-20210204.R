#create a matrix with differences obtained from a set of sequences and output as a heatmap


# Step 1: FASTA to df
require(seqinr)
t <- read.fasta(file = "~/OneDrive - Kemri Wellcome Trust/variants-highlighter/data/dummy-data-02-20210204.fa")
tabl <- data.frame(seq_name = names(t), seqs = unlist(getSequence(t, as.string = T)))
tabl

# Step 2: convert to matrix
seqs <- matrix(unlist(strsplit(as.character(tabl$seqs), "")), ncol = 12, 
               byrow = T, dimnames = list(tabl$seq_name))
seqs

# Step 3: row-wise comparison of the sequences
compare = t(combn(nrow(seqs), 2, FUN = function(x)seqs[x[1],]==seqs[x[2],]))
samples = rownames(compare) = combn(nrow(seqs), 2, FUN = function(x)paste0("seqs", x[1], 
                                                                 "_seqs", x[2]))
compare
storage.mode(compare) = "integer"     # storing the logical as integers
compare
colnames(compare) <- paste0("S", seq(12))
compare

# Step 4: plotting with ggplot2
require(tidyverse)
require(ggplot2)


# 0 = mismatch, 1 = match
compare %>%
  as.data.frame() %>%
  rownames_to_column("seqs_name") %>%
  pivot_longer(-c(seqs_name), names_to = "Sites", values_to = "values") %>%
  ggplot(aes(x = Sites, y = seqs_name, fill = values)) +
  geom_tile() + 
  #geom_text(aes(label = values), color = "white") +
  scale_fill_gradient(low = '#f4a582', high = '#252525') +
  scale_color_gradient(low = '#f4a582', high = '#252525') +
  theme_bw() +
  xlab('sites') +
  ylab('sequences') 
 
  
#'#ca0020','#f4a582','#f7f7f7','#92c5de', '#4393c3', '#2166ac','#252525'

#==============================================================================

