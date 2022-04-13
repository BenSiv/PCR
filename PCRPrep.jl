
cd(raw"C:\Users\admin\Documents\BenSiv_Equinom\PCR")

using Pkg
Pkg.activate(".")


include("SeqGenerator.jl")

n = parse(Int64, ARGS[1])
p = parse(Int64, ARGS[2])

PCRgenerator(n, p)

reader = open(FASTA.Reader, "PCR.fasta")
Reads = [FASTA.sequence(record) for record in reader]
Seq = Reads[1]
Fwd = Reads[2]
Rev = Reads[3]

FwdLoc = findall(String(Fwd), String(Seq))
RevLoc = findall(String(reverse(Rev)), String(complement(Seq)))

if all(length.([FwdLoc, RevLoc]) .== 1)
    FwdLoc = collect(FwdLoc...)
    RevLoc = collect(RevLoc...)
else
    print(length.([FwdLoc, RevLoc]))
end

AmpLoc = first(FwdLoc):last(RevLoc)

Amplicon = Seq[AmpLoc]

using DataStructures
FwdCount = DataStructures.counter(Fwd)
RevCount = DataStructures.counter(Rev)

# Calculate Annealing temprature (Tₘ ᵒC) => http://insilico.ehu.es/tm.php?formula=basic
# Tm = 64.9+41*(yG+zC-16.4)/(wA+xT+yG+zC)

FwdTₘ = 64.9+41*(FwdCount[DNA_G]+FwdCount[DNA_C]-16.4)/(FwdCount[DNA_A]+FwdCount[DNA_T]+FwdCount[DNA_G]+FwdCount[DNA_C]) |> round |> Int
RevTₘ = 64.9+41*(RevCount[DNA_G]+RevCount[DNA_C]-16.4)/(RevCount[DNA_A]+RevCount[DNA_T]+RevCount[DNA_G]+RevCount[DNA_C]) |> round |> Int

# If the Taq speed is 1kb per minute
ElongationTime = (length(Amplicon)/1000)*60 |> round |> Int

using DataFrames
Reaction = DataFrame(Temp = [95, 95, minimum([FwdTₘ,RevTₘ]), 72, 72, 8], 
                     Time = [1, 15, 15, ElongationTime, 2, "∞"], 
                     TimeScale = ["minute", "second", "second", "second", "minute", "-"], 
                     Cycles = [1, 30, 30, 30, 1, 1])
                     

using Statistics
using DataFrames
function ReactionPrep(Volume, DNAConc, PrimerConc = 10; Plasmid = false, Repeats = 1)
    """
    Requiers => Statistics, DataFrames

    Volumn in μl
    Concentration in ng/μl

    DNAConc calculation:
    5-50 ng of genomic DNA or 0.1-1 ng of plasmid DNA
    ~30                       ~0.5
    """

    if Plasmid
        DNAmass = 0.5
    else
        DNAmass = 30
    end

    taqmix = Volume/2
    primers = (0.25/PrimerConc)*Volume
    template = DNAmass/DNAConc
    water = Volume - sum([taqmix, primers, template])

    prep = DataFrame(TaqMix = [taqmix.*rep for rep in Repeats],
                     Primers = [primers.*rep for rep in Repeats],
                     Template = [template.*rep for rep in Repeats],
                     Water = [water.*rep for rep in Repeats],
                     Units = fill("μl", length(Repeats)))

    return prep
end

ReactionPrep(50, 15, 10, Plasmid = false, Repeats = [1,5,20])