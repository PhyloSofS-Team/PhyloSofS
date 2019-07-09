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

function get_infos_from_msa_new(df, msa_total ,index_debut)
    fragment_templates= Int64[]
    fragment_identity_mean = Float64[]
    fragment_identity_max = Float64[]
    fragment_templates_names = []
    alltemplates = []
    index= 1
    loop_end = 0
    global pstart=1
    npos = length(msa_total[1,:])
    nposFrag = 1000
    s = 1
    global i = 1
    while s < index_debut 
        s = s + df.length[i]
        global i += 1
    end
    if s > index_debut
        global i -= 1
        nposFrag = s -index_debut
    end
    loop_start = i
    goon = true
    while goon && i <= size(df,1)
#        println(i)
        identityliste = Float64[]
        templateWellAligned = 0
        templatesIndex = []
        pend = pstart+df.length[i]-1
        if nposFrag < df.length[i]
            pend = pstart + nposFrag -1
        end
        if pend > npos 
            pend = npos
            goon = false
        end
        fragment = msa_total[2:nsequences(msa_total),pstart:pend]
        nposFrag = length(fragment[1,:])
#        println("gaps = ")
#        println()
        cov = 100 .* coverage(fragment) * nposFrag / df.length[i]
#        println(cov[1,:])
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
        loop_end = i
        global pstart = pstart + df.length[i]
        global i+= 1
        nposFrag = 1000
    end
#    println(loop_start)
    return fragment_templates, fragment_identity_mean, fragment_identity_max, fragment_templates_names, alltemplates, loop_end, loop_start
end

gap=Residue('-')
query_name = ARGS[1]
iteration = ARGS[2]
exons = DataFrame(symbol = Char[], length = Int[])
exons_filtered = DataFrame(symbol = Char[], length = Int[])
exons = get_exons(query_name)
# println(exons)
exons_filtered = filter(row -> row[:length] > 7, exons)
msa_total = read("$(query_name)_$(iteration).pir", PIR)
setreference!(msa_total, 1)
adjustreference!(msa_total)
debut_msa = msa_total[1,1:7]
debut_msa = join(string.(debut_msa))
#println(debut_msa)
open("$(query_name)_$(iteration).fa","r") do f
    global seq_lines = readlines(f)
end
sequence = seq_lines[2]
#println(sequence)
index_debut = findfirst(debut_msa, sequence)
#println(index_debut.start)
df = DataFrame(exonNumber = Int[],symbol= Char[], nbTemp = Int[], meanIdentity = Float64[], maxIdentity = Float64[],color = String[],  legend = String[], templatesnames = [])#, meanRmsd = Float64[])
fragments_identity_mean = Float64[]
fragments_identity_max = Float64[]
mean_rmsd = Float64[]
fragments_templates= Int64[]
fragments_templates_names = []
alltemplates = []
fragments_templates, fragments_identity_mean, fragments_identity_max, fragments_templates_names, alltemplates, loop_end , loop_start= get_infos_from_msa_new(exons, msa_total, index_debut.start)
n = exons_filtered.symbol
global it = 0
for i in 1:loop_end-loop_start+1
#      if exons.length[i] > 7
        if(fragments_identity_max[i] < 25)
            push!(df, (loop_start+it, Char(exons.symbol[loop_start+it]), fragments_templates[i], fragments_identity_mean[i], fragments_identity_max[i], "red","0-24.99%", fragments_templates_names[i]))#, mean_rmsd[i]))
        elseif ((fragments_identity_max[i] >=25) && (fragments_identity_max[i] <50))
            push!(df, (loop_start+it, Char(exons.symbol[loop_start+it]), fragments_templates[i], fragments_identity_mean[i], fragments_identity_max[i], "orange","25-49.99%", fragments_templates_names[i]))#, mean_rmsd[i]))
        elseif ((fragments_identity_max[i] >=50) && (fragments_identity_max[i] <75))
            push!(df, (loop_start+it, Char(exons.symbol[loop_start+it]), fragments_templates[i], fragments_identity_mean[i], fragments_identity_max[i], "blue","50-74.99%", fragments_templates_names[i]))#, mean_rmsd[i]))
        else
            push!(df, (loop_start+it, Char(exons.symbol[loop_start+it]), fragments_templates[i], fragments_identity_mean[i], fragments_identity_max[i], "green","75-100%", fragments_templates_names[i]))#, mean_rmsd[i]))
        end
        global it +=1
#      end
end
total_length = 0
iterative_number = DataFrame( symbol = Char[], number = Char[])
for i in 1:it
    if (fragments_identity_max[i] > 25)
        push!(iterative_number, (exons.symbol[loop_start+i-1], '1'))
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
    if(l[i] == '2' && exons_numbers[number]=='2')
        l[i] = '3'
    else
        l[i] = exons_numbers[number]
    end
    global number = number + 1
end
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