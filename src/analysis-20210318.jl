# implementation of BioJulia for the analysis
pwd()

cd("$(homedir())/OneDrive - Kemri Wellcome Trust/variants-highlighter")

# adding necessary packages
using Pkg 


Pkg.add("BioSequences"); Pkg.add("BioAlignments"); Pkg.add("FastaIO")       # implementing biojulia: biosequences


using BioSequences; using BioAlignments; using FastaIO; 

# reading the fasta file

FastaReader("data/dummy-data-02-20210204.fa") do fr
    for (name, seq) in fr
        println("$name : $seq")
    end
end



file = readfasta("data/dummy-data-02-20210204.fa")
ref1 = readfasta("data/ref.fa")

for (name, seq) in FastaReader("data/dummy-data-02-20210204.fa")
    println(">$name", "\n", "$seq")
end


#-----------------------------------------------------------------------------------------------------

# pairwise comparison

function pairwise(x::String)
    return String.(split(x,""))
end


function pairwise(x::Array{Tuple{String, String},1})    #does not work with the readfasta 
    return String[2].(split(x,""))
end


function pairwise1(x::String, y::String, op = ==)
    broadcast(op, split(x,""), split(y,""))
end


function pairwise1(x::Array{String, 1}, y::Array{String, 1} = x, op = ==)
    arr = []
    for each in x
        tmp = Array{Int64,1}()
        for other in y
            push!(tmp, count(broadcast(op, split(each,""), split(other,""))))
        end
        push!(arr, tmp)
    end
    return arr
end

ref = "ACGTAAGCGACG"
seq1 = "ACGTTAGTCGTC"
sequences = ["ACGTTAGTCGTC", "A-GTT-GTCACG", "ACGTTGGTCATC", "ACGTAAGCGACG"]
sequences2 = ["ACGTTAGTCGTC", "A-GTT-GTCACG", "ACGTTGGTCATC", "ACGTAAGCGACG"]

pairwise(sequences)
pairwise(ref1)

pairwise1(ref, seq1)
pairwise1(ref, sequences)
pairwise1(sequences, sequences2)

#-----------------------------------------------------------------------------------------------------

