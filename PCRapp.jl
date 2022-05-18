cd(raw"C:\Users\Ben\Documents\GitHub\PCR")

using Pkg
Pkg.activate(".")

using Gtk

Window = GtkWindow("PCR preparation guide", 400, 600)

vBox = GtkBox(:v)
push!(Window, vBox)

Inputs = GtkLabel("Parameter inputs")
push!(vBox, Inputs)


# Sequence Input
SeqBox = GtkBox(:h)
push!(vBox, SeqBox)

SeqLabel = GtkLabel("Template sequence")
push!(SeqBox, SeqLabel)

Seq = GtkEntry()
push!(SeqBox, Seq)
# Expand the widget to the window size
set_gtk_property!(SeqBox, :expand, Seq, true)


# Forward Input
FwdBox = GtkBox(:h)
push!(vBox, FwdBox)

FwdLabel = GtkLabel("Forward primer")
push!(FwdBox, FwdLabel)

Fwd = GtkEntry()
push!(FwdBox, Fwd)
# Expand the widget to the window size
set_gtk_property!(FwdBox, :expand, Fwd, true)


# Reverse Input
RevBox = GtkBox(:h)
push!(vBox, RevBox)

RevLabel = GtkLabel("Reverse primer")
push!(RevBox, RevLabel)

Rev = GtkEntry()
push!(RevBox, Rev)
# Expand the widget to the window size
set_gtk_property!(RevBox, :expand, Rev, true)


TypeBox = GtkBox(:h)
push!(vBox, TypeBox)

TypeLabel = GtkLabel("Template type")
push!(TypeBox, TypeLabel)

# Initiate a DropDown widget
TempType = GtkComboBoxText()
choices = ["Genomic", "Plasmid"]
for choice in choices
  push!(TempType,choice)
end
# Lets set the active element to be "two"
set_gtk_property!(TempType,:active,1)

signal_connect(TempType, "changed") do widget, others...
  # get the active index
  idx = get_gtk_property(TempType, "active", Int)
  # get the active string 
  # We need to wrap the GAccessor call into a Gtk bytestring
  str = Gtk.bytestring( GAccessor.active_text(TempType) ) 
  println("Active element is \"$str\" at index $idx")
end

push!(TypeBox, TempType)


# Volume Input
VolumeBox = GtkBox(:h)
push!(vBox, VolumeBox)

VolumeLabel = GtkLabel("Reaction volume (μl)")
push!(VolumeBox, VolumeLabel)

Volume = GtkEntry()
push!(VolumeBox, Volume)


LadderBox = GtkBox(:h)
push!(vBox, LadderBox)

LadderLabel = GtkLabel("DNA ladder")
push!(LadderBox, LadderLabel)

# Initiate a DropDown widget
Ladder = GtkComboBoxText()
choices = ["GeneRuler Express", "GeneRuler Mix"]
for choice in choices
  push!(Ladder,choice)
end
# Lets set the active element to be "two"
set_gtk_property!(Ladder,:active,1)

signal_connect(Ladder, "changed") do widget, others...
  # get the active index
  idx = get_gtk_property(Ladder, "active", Int)
  # get the active string 
  # We need to wrap the GAccessor call into a Gtk bytestring
  str = Gtk.bytestring( GAccessor.active_text(Ladder) ) 
  println("Active element is \"$str\" at index $idx")
end

push!(LadderBox, Ladder)


# Repeats Input
RepBox = GtkBox(:h)
push!(vBox, RepBox)

RepLabel = GtkLabel("Repeats")
push!(RepBox, RepLabel)

Repeat1 = GtkEntry()
push!(RepBox, Repeat1)

Repeat2 = GtkEntry()
push!(RepBox, Repeat2)

Repeat3 = GtkEntry()
push!(RepBox, Repeat3)

# Primer Concentration Input
pConcBox = GtkBox(:h)
push!(vBox, pConcBox)

pConcLabel = GtkLabel("Primer Concentration (ng/μl)")
push!(pConcBox, pConcLabel)

pConcentration = GtkEntry()
push!(pConcBox, pConcentration)


# Cycles Input
CyclesBox = GtkBox(:h)
push!(vBox, CyclesBox)

CyclesLabel = GtkLabel("Cycles")
push!(CyclesBox, CyclesLabel)

Cycles = GtkEntry()
push!(CyclesBox, Cycles)


# DNA Concentration Input
dConcBox = GtkBox(:h)
push!(vBox, dConcBox)

dConcLabel = GtkLabel("DNA Concentration (ng/μl)")
push!(dConcBox, dConcLabel)

DNAConc = GtkEntry()
push!(dConcBox, DNAConc)


# Submit button
Button = GtkButton("Submit")
push!(vBox, Button)

using YAML, FASTX, Weave
function on_button_clicked(w)
    # println("The button has been clicked")
    Input = Dict("Seq" => get_gtk_property(Seq, :text, String),
                 "Fwd" => get_gtk_property(Fwd, :text, String),
                 "Rev" => get_gtk_property(Rev, :text, String),
                 "TempType" => Gtk.bytestring( GAccessor.active_text(TempType) ) ,
                 "Volume" =>  parse(Int64, get_gtk_property(Volume, :text, String)),
                 "DNAladder" => Gtk.bytestring( GAccessor.active_text(Ladder) ),
                 "Repeats" => parse.(Int64, [get_gtk_property(Repeat1, :text, String), get_gtk_property(Repeat2, :text, String), get_gtk_property(Repeat3, :text, String)]),
                 "PrimerConc" => parse(Int64, get_gtk_property(pConcentration, :text, String)),
                 "Cycles" => parse(Int64, get_gtk_property(Cycles, :text, String)),
                 "DNAConc" => parse(Int64, get_gtk_property(DNAConc, :text, String))
    )
    YAML.write_file("Input.yml", Input)
    Input = YAML.load_file("Input.yml")

    # Save the Sequences as a fasta file
    SeqRecord = FASTA.Record("Sequnce", Input["Seq"])
    FwdRecord = FASTA.Record("Fwd primer", Input["Fwd"])
    RevRecord = FASTA.Record("Rev primer", Input["Rev"])

    writer = open(FASTA.Writer, "PCR.fasta")
    write(writer, SeqRecord)
    write(writer, FwdRecord)
    write(writer, RevRecord)
    close(writer)

    PCRparams = Dict("FileName" => "PCR.fasta",
                    "Plasmid" => Input["TempType"] == "Plasmid",
                    "Repeats" => Input["Repeats"],
                    "Volume" => Input["Volume"],
                    "PrimerConc" => Input["PrimerConc"],
                    "Cycles" => Input["Cycles"],
                    "DNAladder" => Input["DNAladder"],
                    "DNAConc" => Input["DNAConc"]
                    )
    # Save the PCR parameters as a YAML file
    YAML.write_file("PCRparams.yml", PCRparams)

    weave("PCRprep.jmd")
end
signal_connect(on_button_clicked, Button, "clicked")  

showall(Window)
