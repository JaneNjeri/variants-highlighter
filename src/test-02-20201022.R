# Step 1: FASTA to df

require(seqinr)
t <- read.fasta(file = "/Users/janenjeri/Desktop/Project/Assignment/trial2.fa")
tabl <- data.frame(seq_name = names(t), seqs = unlist(getSequence(t, as.string = T)))
tabl

# Step 2: convert to matrix      #the no. of columns/sites, only works with fixed sequence length:-still working on this
seqs <- matrix(unlist(strsplit(as.character(tabl$seqs), "")), ncol = 10, 
               byrow = T, dimnames = list(tabl$seq_name))
seqs

# Step 3: row-wise comparison of the sequences
compare = t(combn(nrow(seqs), 2, FUN = function(x)seqs[x[1],]==seqs[x[2],]))
samples = rownames(compare) = combn(nrow(seqs), 2, FUN = function(x)paste0("seqs", x[1], 
                                                                 "_seqs", x[2]))
compare
storage.mode(compare) = "integer"     # storing the logical as integers
compare
colnames(compare) <- paste0("S", seq(10))
compare

# Step 4: plotting with ggplot2
require(tidyverse)
require(ggplot2)


# 0 = mismatch, 1 = match
compare %>%
  as.data.frame() %>%
  rownames_to_column("Seqs_name") %>%
  pivot_longer(-c(Seqs_name), names_to = "Sites", values_to = "values") %>%
  ggplot(aes(x = Sites, y = Seqs_name, fill = values)) +
  geom_tile() + 
  geom_text(aes(label = values), color = "white") +
  scale_fill_gradient(low = "red", high = "black") +
  scale_color_gradient(low = "red", high = "black")


#==============================================================================

