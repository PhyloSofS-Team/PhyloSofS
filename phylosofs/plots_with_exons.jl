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
    open("$(query_name).exons_lengths.txt","r") do f
        global lines = readlines(f)
    end
    exons = [parse(Int, ss) for ss in split(lines[2])]
    exons_symbol = split(lines[1])[3]
    for i in 1:length(exons)-1
        push!(df, (exons_symbol[i], (exons[i+1]-exons[i])))
    end
    return df
end

function get_infos_from_msa_new(df, msa_total)
    fragment_templates= Int64[]
    fragment_identity_mean = Float64[]
    fragment_identity_max = Float64[]
    fragment_templates_names = []
    alltemplates = []
    index= 1
    loop_end = 0
    global len=0
    stop = false
    for i in 1:size(exons,1)
        identityliste = Float64[]
        templateWellAligned = 0
        templatesIndex = []
        if stop == true
            break
        end
        global truc = 0
        for j in len+1:length(msa_total[1,:])
            if j == size(msa_total, 2)
                stop = true
                tostore = j
                break
            end
            if msa_total[1,j] != gap
                global truc =  truc + 1
            end
            if truc == exons.length[i]
                global tostore = j
                break
            end
        end
            fragment = msa_total[2:nsequences(msa_total),len+1:tostore]
            cov = 100 .* coverage(fragment)
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
        global len = tostore
    end
    return fragment_templates, fragment_identity_mean, fragment_identity_max, fragment_templates_names, alltemplates, loop_end
end

gap=Residue('-')
query_name = ARGS[1]
iteration = ARGS[2]
exons = DataFrame(symbol = Char[], length = Int[])
exons_filtered = DataFrame(symbol = Char[], length = Int[])
exons = get_exons(query_name)
exons_filtered = filter(row -> row[:length] > 7, exons)
msa_total = read("$(query_name)_$(iteration).pir", PIR,)
df = DataFrame(exonNumber = Int[], nbTemp = Int[], meanIdentity = Float64[], maxIdentity = Float64[],color = String[],  legend = String[], templatesnames = [])#, meanRmsd = Float64[])
fragments_identity_mean = Float64[]
fragments_identity_max = Float64[]
mean_rmsd = Float64[]
fragments_templates= Int64[]
fragments_templates_names = []
alltemplates = []
fragments_templates, fragments_identity_mean, fragments_identity_max, fragments_templates_names, alltemplates, loop_end = get_infos_from_msa_new(exons, msa_total)
n = exons_filtered.symbol
global it = 0
for i in 1:loop_end
    if exons.length[i] > 7
        if(fragments_identity_max[i] < 25)
            push!(df, (it+1, fragments_templates[i], fragments_identity_mean[i], fragments_identity_max[i], "red","0-24.99%", fragments_templates_names[i]))#, mean_rmsd[i]))
        elseif ((fragments_identity_max[i] >=25) && (fragments_identity_max[i] <50))
            push!(df, (it+1, fragments_templates[i], fragments_identity_mean[i], fragments_identity_max[i], "orange","25-49.99%", fragments_templates_names[i]))#, mean_rmsd[i]))
        elseif ((fragments_identity_max[i] >=50) && (fragments_identity_max[i] <75))
            push!(df, (it+1, fragments_templates[i], fragments_identity_mean[i], fragments_identity_max[i], "blue","50-74.99%", fragments_templates_names[i]))#, mean_rmsd[i]))
        else
            push!(df, (it+1, fragments_templates[i], fragments_identity_mean[i], fragments_identity_max[i], "green","75-100%", fragments_templates_names[i]))#, mean_rmsd[i]))
        end
        global it +=1
    end
end
total_length = 0
iterative_number = DataFrame( symbol = Char[], number = Char[])
for i in 1:it
    if (fragments_identity_max[i] > 25)
        push!(iterative_number, (exons.symbol[i], '1'))
    else
        push!(iterative_number, (exons.symbol[i], '2'))
    end
end
df2 = DataFrame(name = [], nb = [])
for i in sequencenames(msa_total)
    counted = count(c -> c == i, collect(alltemplates))
    # if the name is splitted, 2 alignments from the same template result in the same name
    # then messing up the plot
    # if(length(split(i, "_")) >=2)
    #     name = split(i, "_")[1] * "_" * split(i, "_")[2]
    # else
    #     name = split(i, "_")[1]
    # end
    # push!(df2, (name, counted))
    push!(df2, (i, counted))
end
a = groupby(df2, 2, sort=true)
df3 =  DataFrame(numberOfExons = [], nbOfTemplates = [])
for i in 2:length(a)
    push!(df3, (a[i].nb[1], length(a[i].nb)))
end
# creation of every plots and saving it in <name query>.pdf
p1 = @df df bar(:exonNumber, :nbTemp,color = :color, group = :legend, bar_width=0.9,
    title = "Number of templates by exons",
    xlabel = "Exons",
    xlim = (0,size(df,1)+1),
    xticks = (1:1:size(df,1)+3, n[1:size(df,1)]),
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
CSV.write("$(query_name)_$(iteration)_plot1.csv", df)
CSV.write("$(query_name)_$(iteration)_plot_2.csv", df2)
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