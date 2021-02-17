#------------------------------------------------------------------------------

require(seqinr)
require(ggplot2)
require(Biostrings)

t = read.fasta(file = "/Users/janenjeri/Desktop/Project/Assignment/trial.fa")
tbl = data.frame(seq_name=names(t), seqs=unlist(getSequence(t, as.string = T)))
tbl

compare <- function(infile = tbl, ref_name = "Test3_ref",
                    exclude = c("-"))
  {
  require(Biostrings)
  require(seqinr)
  seq0 = read.alignment(infile, format = format)
  s <- nucleotideSubstitutionMatrix(match = 2, mismatch = 0, baseOnly = TRUE)
  seq1 <- as.character(seq0[1:2,2])
  seq2 <- as.character(seq0[3,2])
  seq_diff <- list.string.diff(a = seq1, b = seq2, exclude = exclude)
  return(seq_diff)
  }


ggplot(tbl, aes(x = seq_name, y = seqs, fill = "blue")) + 
  geom_tile()

#------------------------------------------------------------------------------

library(tidyr)
library(dplyr)
f = read.fasta(file = "/Users/janenjeri/Desktop/Project/Assignment/trial.fa")
df = data.frame(seq_name=names(f), seqs=unlist(getSequence(f, as.string = T)))

df %>%
  separate_rows(seqs) %>%
  distinct %>%
  inner_join(., ., by = "seqs") %>%
  count(seq_name.x, seq_name.y) %>%
  complete(seq_name.x, seq_name.y)

Scan <- function(x) scan(text = x, what = "", sep = ",", quiet = TRUE)
countSame <- function(x, y) length(intersect(Scan(x), Scan(y)))
x <- setNames(df$seqs, df$seq_name)
outer(x, x, Vectorize(countSame))

#-----------------------------------------------------------------------------
# test2
install.packages("stringdist")
library(stringdist)

d = read.fasta(file = "/Users/janenjeri/Desktop/Project/Assignment/trial.fa")
dat = data.frame(seq_name=names(d), seqs=unlist(getSequence(t, as.string = T)))

# Distance methods available in stringdist
dist.methods = c("jw", "soundex")

# Try all the methods with the sample data
sapply(dist.methods, function(m) stringdist(dat[3,2],dat[1:2,2], method=m))

#------------------------------------------------------------------------------

#test
test<-data.frame(a=c("inf",1,"inf"),b=c("nan",3,"nan"))

############################################################################
# OTHER WORK

# Dummy data
x <- LETTERS[1:20]
y <- paste0("var", seq(1,20))
data <- expand.grid(X=x, Y=y)
data$Z <- runif(400, 0, 5)

# Heatmap 
ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile()

#----------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
require(Biostrings)


df <- data.frame(names=c("A ADAM", "S BEAN", "A APPLE", "J BOND", "J BOND"), 
                 v1=c("Test_a", "Test_b", "Test_a", "Test_b", "Test_b"), 
                 v2=c("Test_c", "Test_c", "Test_d", "Test_d", "Test_d"))
library(tidyverse)
map2(df, c('soundex', 'jw', 'jw'), ~stringdist::stringdistmatrix(.x, method = .y)) %>% 
  map_df(broom::tidy, .id = 'var') %>% 
  spread(var, distance)

#-----------------------------------------------------------------------------
s <- nucleotideSubstitutionMatrix(match = 2, mismatch = 0, baseOnly = TRUE)
s # Print out the matrix

compare <- function(fil){
  
}




#-----------------------------------------------------------------------------
# generate random data
set.seed(314)
df2 <- data.frame(id1 = sample(letters, size=5, replace=T))
df2$id2 <- sample(letters, size=5, replace=T)
df2$r.values <- runif(n=5)

# function to return pairwise matrix
xtabs(r.values ~ id1 + id2, data=df)

#-----------------------------------------------------------------------------
library(ggplot2)
BiocManager::install("ggtree")
library(ggtree)

set.seed(102)
tree <- rtree(60)
p <- ggtree(tree)
p1 <- p + geom_hilight(node=62) + geom_hilight(node=88, fill="red")
p1
dat <- data.frame(id=c(62, 88), type=c("A", "B"))
p2 <- p + geom_hilight(data=dat, mapping=aes(node=id, fill=type))
p2

#-----------------------------------------------------------------------------

