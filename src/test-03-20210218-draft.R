# step 1: FASTA to df
require(seqinr)
fastafile <- read.fasta(file = "~/OneDrive - Kemri Wellcome Trust/variants-highlighter/data/dummy-data-02-20210204.fa")
tab <- data.frame(samples = names(fastafile), seqs = unlist(getSequence(fastafile, as.string = T)))
tab

reffile <- read.fasta(file = "~/OneDrive - Kemri Wellcome Trust/variants-highlighter/data/ref.fa")
refvector <- as.array(reffile)
refvector


typeof(tab)



# Step 2: convert to matrix
sequences <- matrix(unlist(strsplit(as.character(tab$seqs), "")), ncol = 12, 
               byrow = T, dimnames = list(tab$samples))
sequences


# step 2: compare or checking for similarity btwn sequences in a column
# if similar - print T
# if different - print F
# if missing - print N


for ( i in 1:length(tab$seqs)) {
  print(tab$seqs[[i]])
}




