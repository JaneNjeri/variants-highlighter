# implementation of BioJulia for the analysis

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



readfasta("data/dummy-data-02-20210204.fa")


for (name, seq) in FastaReader("data/dummy-data-02-20210204.fa")
    println(">$name", "\n", "$seq")
end


