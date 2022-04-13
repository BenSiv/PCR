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

using Plots

GeneRuller = [100, 300, 500, 750, 1000, 1500, 2000, 3000, 5000]

GelPlot = plot(xlimit = (0,5), ylimit = (0,6000), grid = false, legend = false, yaxis=(:log10, [10, :auto]))
for band in GeneRuller
    if band in [500, 1500]
        plot!([1,2], [band,band], color = :black, width = 3)
    else
        plot!([1,2], [band,band], color = :black)
    end
end

plot!([3,4], [length(Amplicon),length(Amplicon)], color = :black, width = 3)

savefig("PCR.png")

SeqPlot = plot([1,n], [44,44], xlimit = (0,n+1), ylimit = (0,100), grid = false, legend = false, color = :black, width = 3)
plot!([minimum(FwdLoc),maximum(FwdLoc)],[46,46], color = :black, width = 3)
plot!([minimum(RevLoc),maximum(RevLoc)],[46,46], color = :black, width = 3)

savefig("SeqPlot.png")