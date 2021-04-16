# from stackoverflow

dat <- data.frame(sequence=c(1,2,2,2), start=c(1,1,2001,8000), stop=c(9000,2000,7999,9000), type=c("mapped","mapped","deletion","mapped"))


library(ggplot2)

g <- ggplot(data=dat, mapping=aes(ymin=0, ymax=1, xmin=start, xmax=stop, fill=type)) +
  geom_rect() + facet_grid(sequence~., switch="y") +
  labs(x="Position (BP)", y="Sequence / Strain", title="Mapped regions for all sequences") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))

g



#----------------------------------------------------------------------------------------------------------
# not sure whats happening here lol

library(bio3d)

# Read alignment
aln <- read.fasta(system.file("~/OneDrive - Kemri Wellcome Trust/variants-highlighter/data/test.fasta",package="bio3d"))

## alignment plot
plot(aln, labels=basename.pdb(aln$id))

## Works also for a 'pdbs' object
attach(transducin)
plot(pdbs)

detach(transducin)

#-----------------------------------------------------------------------------------------------------------
# installing packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("msa")
install.packages("ggmsa")

# loading packages
library(msa)
library(ggmsa)
library(Biostrings)

# loading data
mySequences <- readDNAStringSet("data/test.fasta", format = "fasta")
mySequences

# performing msa
my_msa <- msa(mySequences)
my_msa

# write to file
rownames(my_msa)
detail(my_msa)

alignment2Fasta <- function(alignment, filename) {
  sink(filename)
  
  n <- length(rownames(alignment))
  for(i in seq(1, n)) {
    cat(paste0('>', rownames(alignment)[i]))
    cat('\n')
    the.sequence <- toString(unmasked(alignment)[[i]])
    cat(the.sequence)
    cat('\n')  
  }
  
  sink(NULL)
}

  
alignment2Fasta(my_msa, 'results/test_aln.fasta')


# plotting msa
file_align = readDNAMultipleAlignment("results/test_aln.fasta", format = "fasta")
ggmsa(file_align, start = 1, end = 66)

file_align = unmasked(file_align)
names(file_align)[1]
ref = file_align[1]


bm = sapply(1:length(file_align),function(i){
  as.numeric(as.matrix(file_align[i])==as.matrix(ref))
})

bm = t(bm)
rownames(bm) = names(file_align)


install.packages("pheatmap")
library(pheatmap)
pheatmap(bm[nrow(bm):1,1:66],cluster_rows=FALSE,cluster_cols=FALSE)


