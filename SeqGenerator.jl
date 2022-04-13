
using BioSequences, FASTX

"""
    PCRgenerator(n,p)

Generates a FASTA file with 3 records:

Random DNA sequence 'n' bp long  

Fwd primer 'p' bp long  

Rev primer 'p' bp long  

The Amplicon is located randomly on the sequence
"""
function PCRgenerator(n,p)

    # Sequence
    Seq = randseq(DNAAlphabet{4}(), n)
    CompSeq = complement(Seq)


    # Primers
    #FwdLoc = rand(1:Int(round(n/2)-p))
    FwdLoc = rand(1:(n-2p))

    #RevLoc = rand(1:Int(round(n/2)-p))
    RevLoc = rand(1:n-FwdLoc)


    Fwd = Seq[FwdLoc:FwdLoc+p]
    Rev = reverse(CompSeq[n-(RevLoc+p):n-RevLoc])

    SeqRecord = FASTA.Record("Sequnce", Seq)
    FwdRecord = FASTA.Record("Fwd primer", Fwd)
    RevRecord = FASTA.Record("Rev primer", Rev)

    writer = open(FASTA.Writer, "PCR.fasta")
    write(writer, SeqRecord)
    write(writer, FwdRecord)
    write(writer, RevRecord)
    close(writer)

end
