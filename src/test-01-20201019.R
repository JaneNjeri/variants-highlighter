# FASTA to dataframe
# 1:
library(Biostrings)

fas_file <- readDNAStringSet("/Users/janenjeri/Desktop/Project/Assignment/trial.fa")
seq_name = names(fas_file)
sequence = paste(fas_file)
df <- data.frame(seq_name, sequence)
names(df) <- NULL


m <- as.matrix(df)
m

m2 <- matrix(m, ncol = ncol(df), dimnames = NULL)
m2



#ifelse(df[1,]==df[2,], 1, 0)

#------------------------------------------------------------------------------
# 2: 

fasta2dataframe=function(fastaFile){
  s = readDNAStringSet(fastaFile)
  RefSeqID = names(s)
  RefSeqID = sub(" .*", "", RefSeqID) 
  #erase all characters after the first space: regular expression matches a space followed by any sequence of characters and sub replaces that with a string having zero  characters 
  
  for (i in 1:length(s)){
    seq[i]=toString(s[i])
  }
  
  RefSeqID_seq=data.frame(RefSeqID,seq)
  return(RefSeqID_seq)
}

mydf = fasta2dataframe("/Users/janenjeri/Desktop/Project/Assignment/trial.fa")

#------------------------------------------------------------------------------
# 3: this function  does not put the NAs into consideration
ReadFasta<-function(file) {
  fasta<-readLines(file)          # Read the file line by line
  ind<-grep(">", fasta)           # Identify header lines
  s<-data.frame(ind=ind, from=ind+1, to=c((ind-1)[-1], length(fasta)))       # Identify the sequence lines
  seqs<-rep(NA, length(ind))      # Process sequence lines
  for(i in 1:length(ind)) {
    seqs[i]<-paste(fasta[s$from[i]:s$to[i]], collapse="")
  }
 
  DF<-data.frame(name=gsub(">", "", fasta[ind]), sequence=seqs)       # Create a data frame
  return(DF)                 # Return the data frame as a result object from the function
}

# Usage example
seqs<-ReadFasta("/Users/janenjeri/Desktop/Project/Assignment/trial.fa")
seqs

seqs1 <- matrix(unlist(strsplit(as.character(seqs$sequence), "")), ncol = 5,
                byrow = T, dimnames = list(seqs$name))
seqs1

compare = t(combn(nrow(seqs1), 2, FUN=function(x)seqs1[x[1],]==seqs1[x[2],]))
rownames(compare) = combn(nrow(seqs1), 2, FUN=function(x)paste0("sequence", x[1], 
                                                                "_sequence",x[2]))
compare

storage.mode(compare) = "integer"




# install.packages("ape")
# require(ape)
# reads <- read.dna("/Users/janenjeri/Desktop/Project/Assignment/trial.fa", 
#                   format = "fasta", as.character = F, as.matrix = T)
# reads
# 
# View(reads)
#-----------------------------------------------------------------------------
