using BioAlignments
using MIToS.MSA
using MIToS.SIFTS
using MIToS.PDB
using StatsPlots
using Statistics
using Printf
using DataFrames
using Plots
using Combinatorics
using CSV



function fasta2pir(aln)
    seq_pos = 0
    ref_pos = 0
    last_seq_pos = 0
    seq2ref = Dict{Int, Int}()
    for (seq_res, ref_res) in alignment(aln)
        if seq_res != '-'
            seq_pos += 1
        end
        if ref_res != '-'
            ref_pos += 1
        end
        if seq_pos != last_seq_pos
            seq2ref[seq_pos] = ref_pos
            last_seq_pos = seq_pos
        end
    end
    seq2ref
end



function get_exons(query_name)
    df = DataFrame(symbol = Char[], length = Int[])
    open("$(query_name)_$(iteration).exons_lengths.txt","r") do f
        global lines = readlines(f)
    end
    exons = [parse(Int, ss) for ss in split(lines[2])]
    exons_symbol = split(lines[1])[3]
    for i in 1:length(exons)-1
        push!(df, (exons_symbol[i], (exons[i+1]-exons[i])))
    end
    return df
end

function get_infos(exons, msa_total ,index_debut, alignment)
    new_f2p = Dict( v => k for (k, v) in f2p )
    fragment_templates= Int64[]
    fragment_identity_mean = Float64[]
    fragment_identity_max = Float64[]
    fragment_templates_names = []
    alltemplates = []
    first_pos_pir = Int64[]
    last_pos_pir = Int64[]
    index= 1
    loop_end = 0
    global i = 1
    loop_start = i
    goon = true
    global debut = 1
    while goon && i <= size(exons,1)
#         println(i)
        identityliste = Float64[]
        templateWellAligned = 0
        templatesIndex = []
#         println(debut)
#         println(debut+exons.length[i]-1)

        fin = debut+exons.length[i]-1
#         println(fin)
        
        # if returns a 0 element array, push 0 in everything
            fragment_index = [ new_f2p[i] for i in debut:fin  if i in keys(new_f2p) ]
#                     println("fragment index = $(fragment_index)")
        if isempty(fragment_index)
            push!(first_pos_pir, -1)
            push!(last_pos_pir, -1)
            push!(templatesIndex, "None")
            push!(alltemplates, "None")
            push!(fragment_templates_names, 0)
            push!(fragment_templates, 0)
            push!(fragment_identity_mean, 0)
            push!(fragment_identity_max, 0)
        else
#             println(fragment_index)
            index_debut = fragment_index[1]
            index_end = fragment_index[end]
#             println("index_debut = $(index_debut)")
#             println(" end = $(index_end)")
#         if index_end == size(msa_total, 2)
#             goon= false
#             break
#         end
#         if f2p[debut] == f2p[debut+exons.length[i]-1]
#             if f2p[debut] != size(msa_total, 2)
#                 println("go to the next iteration")
#                 loop_start = i+1
#             end
#         else
        
            fragment = msa_total[2:nsequences(msa_total),index_debut:index_end]
#             println(fragment)
            pir_positions = getcolumnmapping(fragment)
#             println("pir pos")
#             println(pir_positions)
            pir_first = pir_positions[1]
            pir_last = pir_positions[end]
            push!(first_pos_pir, pir_first)
            push!(last_pos_pir, pir_last)
            cov = 100 .* coverage(fragment) * size(fragment, 2) / exons.length[i]
#             println(cov)
            for k in 1:length(cov)
                if cov[k] >= 50
                    templateWellAligned += 1
                    push!(templatesIndex, names(cov)[1][k])
                    push!(alltemplates, names(cov)[1][k])
                end
            end
            push!(fragment_templates_names, templatesIndex)
            push!(fragment_templates, length(templatesIndex))
            for l in templatesIndex
                push!(identityliste, percentidentity(fragment[1,:], fragment[l,:]))
            end
            if(isempty(identityliste)) # to avoid LOADERROR when trying to do mean() or max() of an empty vector
                push!(fragment_identity_mean, 0)
            push!(fragment_identity_max, 0)
            else
                push!(fragment_identity_mean, mean(identityliste))
                push!(fragment_identity_max, maximum(identityliste))
            end
        end
        loop_end = i
        global debut = debut + exons.length[i]
        global i = i+1
    end
    return fragment_templates, fragment_identity_mean, fragment_identity_max, fragment_templates_names, alltemplates, loop_end, loop_start, first_pos_pir, last_pos_pir
end

# query_name = "./test_258682/ENSG00000058404_ENST00000258682/ENSG00000058404_ENST00000258682"#ARGS[1]
# iteration = 1#ARGS[2]
# query_name = "./FMR1/ENSG00000102081_ENST00000370475/ENSG00000102081_ENST00000370475"#ARGS[1]
# iteration = 4#ARGS[2]
query_name = ARGS[1]
iteration = ARGS[2]


gap=Residue('-')
open("$(query_name)_$(iteration).fa","r") do f
    global seq_lines = readlines(f)
end
sequence = seq_lines[2]
exons = DataFrame(symbol = Char[], length = Int[])
exons_filtered = DataFrame(symbol = Char[], length = Int[])
exons = get_exons(query_name)
# println(exons)
exons_filtered = filter(row -> row[:length] > 7, exons)
msa_total = read("$(query_name)_$(iteration).pir", PIR, deletefullgaps=false)
setreference!(msa_total, 1)
adjustreference!(msa_total)
debut_msa = msa_total[1,1:7]
debut_msa = join(string.(debut_msa))
costmodel = AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-1)
#aln = pairalign(SemiGlobalAlignment(), join(msa_total[1,:]),seq_lines[2], costmodel)
aln = pairalign(SemiGlobalAlignment(), join(msa_total[1,:]),seq_lines[2], costmodel)
f2p = fasta2pir(aln)
#println(debut_msa)
#println(sequence)
index_debut = findfirst(debut_msa, sequence)
# println(index_debut.start)
df = DataFrame(exonNumber = Int[],symbol= Char[], nbTemp = Int[], meanIdentity = Float64[], maxIdentity = Float64[],color = String[],  legend = String[], first_pos_pir = Int64[], last_pos_pir = Int64[], templatesnames = [])
fragments_identity_mean = Float64[]
fragments_identity_max = Float64[]
mean_rmsd = Float64[]
fragments_templates= Int64[]
fragments_templates_names = []
alltemplates = []
first_pos_pir = Int64[]
last_pos_pir = Int64[]
fragments_templates, fragments_identity_mean, fragments_identity_max, fragments_templates_names, alltemplates, loop_end , loop_start, first_pos_pir, last_pos_pir= get_infos(exons, msa_total, index_debut.start, f2p)
n = exons_filtered.symbol
global it = 0
for i in 1:loop_end-loop_start+1
#     println(i)
#      if exons.length[i] > 7
        if(fragments_identity_max[i] < 25)
            push!(df, (loop_start+it, Char(exons.symbol[loop_start+it]), fragments_templates[i], fragments_identity_mean[i], fragments_identity_max[i], "red","0-24.99%", first_pos_pir[i], last_pos_pir[i], fragments_templates_names[i]))#, mean_rmsd[i]))
        elseif ((fragments_identity_max[i] >=25) && (fragments_identity_max[i] <50))
            push!(df, (loop_start+it, Char(exons.symbol[loop_start+it]), fragments_templates[i], fragments_identity_mean[i], fragments_identity_max[i], "orange","25-49.99%", first_pos_pir[i], last_pos_pir[i], fragments_templates_names[i]))#, mean_rmsd[i]))
        elseif ((fragments_identity_max[i] >=50) && (fragments_identity_max[i] <75))
            push!(df, (loop_start+it, Char(exons.symbol[loop_start+it]), fragments_templates[i], fragments_identity_mean[i], fragments_identity_max[i], "blue","50-74.99%", first_pos_pir[i], last_pos_pir[i], fragments_templates_names[i]))#, mean_rmsd[i]))
        else
            push!(df, (loop_start+it, Char(exons.symbol[loop_start+it]), fragments_templates[i], fragments_identity_mean[i], fragments_identity_max[i], "green","75-100%", first_pos_pir[i], last_pos_pir[i], fragments_templates_names[i]))#, mean_rmsd[i]))
        end
        global it +=1
#      end
end
total_length = 0
iterative_number = DataFrame( symbol = Char[], number = Char[])
for i in 1:size(df, 1)
    if (exons.length[i] < 5)
        if (df.maxIdentity[i]) >25
            push!(iterative_number, (exons.symbol[loop_start+i-1], '1'))
        elseif (df.first_pos_pir[i] == -1)
            push!(iterative_number, (exons.symbol[loop_start+i-1], '0'))
        else
            push!(iterative_number, (exons.symbol[loop_start+i-1], '3'))
        end
    elseif (df.maxIdentity[i] > 25)
        push!(iterative_number, (exons.symbol[loop_start+i-1], '1'))
    elseif (df.first_pos_pir[i] == -1)
        push!(iterative_number, (exons.symbol[loop_start+i-1], '0'))
    else
        push!(iterative_number, (exons.symbol[loop_start+i-1], '2'))
    end
end
df2 = DataFrame(name = [], nb = [])
for i in sequencenames(msa_total)
    counted = count(c -> c == i, collect(alltemplates))
    push!(df2, (i, counted))
end
a = groupby(df2, 2, sort=true)
df3 =  DataFrame(numberOfExons = [], nbOfTemplates = [])
for i in 2:length(a)
    push!(df3, (a[i].nb[1], length(a[i].nb)))
end



open("$(query_name)_annotated.pir","r") do f
    global lines = readlines(f)
end
first = findfirst(isequal(iterative_number.symbol[1]), lines[2])
last = findlast(isequal(iterative_number.symbol[length(iterative_number.symbol)]), lines[2])

global exons_numbers = ""
for i in 1:length(iterative_number.symbol)
    global exons_numbers = exons_numbers * (iterative_number.number[i]^((findlast(isequal(iterative_number.symbol[i]), lines[2])+1)-findfirst(isequal(iterative_number.symbol[i]), lines[2])))
end
l = collect(lines[4])
global number = 1
for i in first:last
    # if already runned once but still not found, do not rerun
    if(l[i] == '2' && exons_numbers[number]=='2')
        l[i] = '3'
    # if runned once but not in the window search of hhblits, do not change the nuber
    elseif (l[i] == '2' && exons_numbers[number]=='0')
        l[i] = '2'
    else
        l[i] = exons_numbers[number]
    end
    global number = number + 1
end
                println(join(l))
open("$(query_name)_annotated.pir", "w") do f
    write(f, lines[1]*"\n")
    write(f, lines[2]*"\n")
    write(f, lines[3]*"\n")
    write(f, join(l))
end


# creation of every plots and saving it in <name query>.pdf
if !isempty(df)
    p1 = @df df bar(:nbTemp,color = :color, group = :legend, bar_width=0.9,
        title = "Number of templates by exons",
        xlabel = "Exons",
        xlim = (0,size(df,1)+1),
        xticks = (1:1:size(df,1)+2, :symbol),
        ylabel="Number of templates",
        ylim = (0,maximum(:nbTemp)+1),
        yticks = 0:5:maximum(:nbTemp),
        legend = :topright)

    p4 = @df df3 bar(:numberOfExons, :nbOfTemplates,
        title = "number of exons covered by the number of templates",
        xlabel = "Number of exons covered",
        ylabel="Number of templates",
        xlim = (0, loop_end+1),
        legend=false)
    p = plot(p1,p4,size = (1280, 720))
    savefig(p, "$(query_name)_$(iteration)_plots.pdf")

end
CSV.write("$(query_name)_$(iteration)_plot1.csv", df)
CSV.write("$(query_name)_$(iteration)_plot_2.csv", df2)

# p = plot(p1,p4,size = (1280, 720))