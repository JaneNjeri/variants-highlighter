#create a matrix with differences obtained from a set of sequences and output as a heatmap


# Step 1: FASTA to df
require(seqinr)
t <- read.fasta(file = "~/OneDrive - Kemri Wellcome Trust/variants-highlighter/data/test.fasta")
tabl <- data.frame(seq_name = names(t), seqs = unlist(getSequence(t, as.string = T)))
tabl

ref = read.fasta("~/OneDrive - Kemri Wellcome Trust/variants-highlighter/data/ref.fa")
ref.tab <- data.frame(seq_name = names(ref), seqs = unlist(getSequence(ref, as.string = T)))
ref.tab

# Step 2: convert to matrix
seqs <- matrix(unlist(strsplit(as.character(tabl$seqs), "")), ncol = 42,
               byrow = T, dimnames = list(tabl$seq_name))
seqs

seqref <- matrix(unlist(strsplit(as.character(ref.tab$seqs), "")), ncol = 42, byrow = T)
seqref

rownames(seqref) <- "ref"


# Step 3: row-wise comparison of the sequences
compare = t(combn(nrow(seqs), 2, FUN = function(x)seqs[x[1],]==seqs[x[2],]))
samples = rownames(compare) = combn(nrow(seqs), 2, FUN = function(x)paste0("sample", x[1],
                                                                 "_sample", x[2]))



apply(cbind(seqref, seqs), 1, function(x)isTRUE(all.equal(x[1:2],x[3:4])))

seqref==seqs

#-----------------------------------------------------------------------------------------

# Fancy way.
# similarity.matrix <- apply(seqs, 2, function(x)colSums(x==seqs))
# diag(similarity.matrix)<-0


# More understandable. But verbose.
#similarity.matrix <- matrix(nrow=ncol(seqs),ncol=ncol(seqs))

#for (col in 1:ncol(seqs)) {
  #matches <- seqs[,col]==seqs
  # match.counts <- colSums(matches)
  # match.counts[col] <- 0 # Set the same column comparison to zero.
  # similarity.matrix[,col] <- match.counts


#------------------------------------------------------------------------------


n = seq_len(ncol(seqs))
id <- expand.grid(n, n)
out <- matrix(colSums(seqs[, id[,1]] == seqs[, id[,2]]), ncol = length(n))
diag(out) <- 0

#------------------------------------------------------------------------------------------

compare
storage.mode(compare) = "integer"     # storing the logical as integers
compare
colnames(compare) <- paste0("S", seq(42))
compare

# Step 4: plotting with ggplot2
require(tidyverse)
require(ggplot2)


# 0 = mismatch, 1 = match
matches %>%
  as.data.frame() %>%
  rownames_to_column("seqs_name") %>%
  pivot_longer(-c(seqs_name), names_to = "Sites", values_to = "values") %>%
  ggplot(aes(x = Sites, y = seqs_name, fill = values)) +
  geom_tile() + 
  geom_text(aes(label = values), color = "white") +
  scale_fill_gradient(low = '#f4a582', high = '#252525') +
  scale_color_gradient(low = '#f4a582', high = '#252525') +
  theme_bw() +
  xlab('Alignment Position') +
  ylab('sequences') 


#'#ca0020','#f4a582','#f7f7f7','#92c5de', '#4393c3', '#2166ac','#252525'

#==============================================================================



