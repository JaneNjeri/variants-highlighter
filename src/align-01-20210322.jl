ref = "ACGTAAGCGACG"
seq1 = "ACGTTAGTCGTC"
sequences = ["ACGTTAGTCGTC", "A-GTT-GTCACG", "ACGTTGGTCATC", "ACGTAAGCGACG"]


let alphabetDNA = [1, 2, 3, 4, 5]
#                  A  C  G  T  ? 
    global letter2num
    function letter2num(c::Char)
        if c == 'A' || c == 'a'
            return 1
        elseif c == 'C' || c == 'c'
            return 2
        elseif c == 'G' || c == 'g'
            return 3
        elseif c == 'T' || c == 't' || c == 'U' || c == 'u'
            return 4
        else
            error("unknown base $c")
        end
    end
end   


function transformseq(seq::String)
    N = length(seq)
    vec = [letter2num(seq[i]) for i in 1:N]
    return vec
end



