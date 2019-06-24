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

# functions for the RMSD
# functions needed to map the MSA to the pdb/cif files
is_pdb(name) = occursin(r"^\w{4}", name) #occursin(r"^\w{4}_\w$", name) gets only AB12_A and not AB12.
# since the query is UKNP, it is not taken

get_cif_file_name(pdb_name) = "$(split(pdb_name, '_')[1]).cif"

get_query_file_name(pdb_name) = "$(split(pdb_name, '_')[1]).B99990001.pdb"

get_pdb_file(pdb_name) = "$(split(pdb_name, '_')[1]).pdb.gz"

get_chain(pdb_name) = string(split(pdb_name, '_')[2])

function msa2pdbres(pdb_seq, msa_seq)
    costmodel = CostModel(match=0, mismatch=10, insertion=1, deletion=1)
    aln = pairalign(EditDistance(), String(pdb_seq), String(msa_seq), costmodel)
    msa2pdbres = Array{Int}(undef, length(msa_seq))
    pdb_pos=0
    msa_pos=0
    for (pdbres, msares) in alignment(aln)
        if msares != '-'
            msa_pos += 1
            if pdbres != '-'
                pdb_pos +=1
                msa2pdbres[msa_pos]=pdb_pos
            else
                msa2pdbres[msa_pos]=0
            end
        else
            if msares != '-'
                msa_pos += 1
            end
        end
    end
    msa2pdbres
end

# cut a sequence in 10 fragments stored in a vector
function get_10_num_from_int(a)
    vect = Int64[]
    push!(vect, 0)
    b = floor(a/10)
    for i in 1:10
        push!(vect, b*i)
    end
    vect
end

"""
Takes as an input a vector (length of the fragments of the msa) and returns
- the number of templates well aligned
- mean identity for each fragment for those aligned templates (as a vector)
- max identity for each fragment for those aligned templates (as a vector)
- the name of each templates selected by fragments
"""
function get_infos_from_msa(vect, msa)
    fragment_templates= Int64[]
    fragment_identity_mean = Float64[]
    fragment_identity_max = Float64[]
    fragment_templates_names = []
    alltemplates = []
    for i in 1:length(vect)-1
        identityliste = Float64[]
        templateWellAligned = 0
        templatesIndex = []
        fragment = msa_total[2:nsequences(msa_total),vect[i]+1:vect[i+1]]
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
    end
    return fragment_templates, fragment_identity_mean, fragment_identity_max, fragment_templates_names, alltemplates
end




"""
Takes as an input the vector containing the fragments length, the name of the pir file
And returns a vector containing the mean RMSD for all fragments, for all pairwise .PIR files
This functions works with 5 pairwise .pir files named each pir_name_1, pir_name_2...
"""
function pairwise_RMSD_calculation_by_fragments(vect, query_name)
    for l in 1:length(vect)-1
        rms = Float64[]
        for duo in comb
            go = true
            msa_1 = read("$(query_name)_$(duo[1]).pir", PIR, deletefullgaps=false)
            msa_2 = read("$(query_name)_$(duo[2]).pir", PIR, deletefullgaps=false)
            # when using HHmakemodel.py to create the .Pir file, if the cif file of the template is not found,
            # the script removes it from the alignment.
            # This can result in a pir file with only the query, messing up the next steps.
            try
                x = msa_1[2,:]
            catch y
                println("pir number $(duo[1]) does not have a template in it")
                go = false
            end
            try
                x = msa_2[2,:]
            catch y
                println("pir number $(duo[2]) does not have a template in it")
                go = false
            end
            if go
                seq1=Residue[]
                seq2=Residue[]
                t = []
                push!(t, length(msa_1[1,:]))
                push!(t, length(msa_2[1,:]))
                for i in 1:minimum(t)
                    if(msa_2[1,i]!=gap) #if query sequence[i] is not a gap
                        if(msa_1[2,i]!=gap && msa_2[2,i]!=gap)
                            push!(seq1,msa_1[2,i])
                            push!(seq2,msa_2[2,i])
                        end
                    end
                end

                open("new_pir_pairwise.pir", "w") do file
                    write(file, ">P1;$(sequencenames(msa_1)[2])\n")
                    write(file, "$(getannotsequence(msa_1, sequencenames(msa_1)[2], "Title"))\n")
                    write(file, "$(join(seq1))*\n")
                    write(file,">P1;$(sequencenames(msa_2)[2])\n")
                    write(file, "$(getannotsequence(msa_2, sequencenames(msa_2)[2], "Title"))\n")
                    write(file, "$(join(seq2))*")
                end
                msa_ref = read("new_pir_pairwise.pir", PIR,generatemapping=true, useidcoordinates=true)
                # test without some loops
                name_i = sequencenames(msa_ref)[1] #to erase
                name_j = sequencenames(msa_ref)[2] #to erase
                #msa_ref_i = gapstrip!(msa_1)
                #msa_ref_j = gapstrip!(msa_2)
                #name_i = sequencenames(msa_ref_i)[2]
                #name_j = sequencenames(msa_ref_j)[2]
                pdbfile_i = downloadpdb("$(split(name_i, '_')[1])", format=PDBFile)
                pdbfile_j = downloadpdb("$(split(name_j, '_')[1])", format=PDBFile)
                PDB_i = read(pdbfile_i, PDBFile, model = "1",group = "ATOM",chain = get_chain(name_i))
                PDB_j = read(pdbfile_j, PDBFile, model = "1",group = "ATOM",chain = get_chain(name_j))
                msa2mapping = Dict{String, Vector{Int}}()
                msa2pdb_res = Dict{String, Vector{PDBResidue}}()
                msa2pdb_res[name_i] = PDB_i
                msa2pdb_res[name_j] = PDB_j
                pdb_seq_i = [ three2residue(res.id.name) for res in PDB_i ]
                #msa_seq_i = [ res for res in getsequence(msa_ref_i, name_i) if res != GAP ]
                msa_seq_i = [ res for res in getsequence(msa_ref, name_i) if res != GAP ] #to erase
                msa2mapping[name_i] = msa2pdbres(pdb_seq_i, msa_seq_i)
                pdb_seq_j = [ three2residue(res.id.name) for res in PDB_j ]
                #msa_seq_j = [ res for res in getsequence(msa_ref_j, name_j) if res != GAP ]
                msa_seq_j = [ res for res in getsequence(msa_ref, name_i) if res != GAP ] #to erase
                msa2mapping[name_j] = msa2pdbres(pdb_seq_j, msa_seq_j)
                #msacol2msares_i = getsequencemapping(msa_ref_i, name_i)
                #msacol2msares_j = getsequencemapping(msa_ref_j, name_j)
                msacol2msares_i = getsequencemapping(msa_ref, name_i)  #to erase
                msacol2msares_j = getsequencemapping(msa_ref, name_j)  #to erase
                new_vect = get_10_num_from_int(length(msacol2msares_i))
                pdbres_i = PDBResidue[]
                pdbres_j = PDBResidue[]
                try #ugly, need to change this
                    for (msares_i, msares_j) in zip(msacol2msares_i[new_vect[l]+1:new_vect[l+1]], msacol2msares_j[new_vect[l]+1:new_vect[l+1]])
                        if msares_i != 0 &&  msares_j != 0
                            pdbres_i_index = msa2mapping[name_i][msares_i]
                            pdbres_j_index = msa2mapping[name_j][msares_j]
                            if pdbres_i_index != 0 && pdbres_j_index != 0
                                push!(pdbres_i, msa2pdb_res[name_i][pdbres_i_index])
                                push!(pdbres_j, msa2pdb_res[name_j][pdbres_j_index])
                            end
                        end
                    end
                    @assert length(pdbres_i) == length(pdbres_j)
                    superimposed_i, superimposed_j, RMSD = superimpose(pdbres_i, pdbres_j)
                    push!(rms, RMSD)
                catch x
                    #println(x)
                    #println(-1)
                    push!(rms, -1)
                end
            end
        end
        push!(mean_rmsd, mean(rms))
    end
    mean_rmsd
end

gap=Residue('-')
query_name = ARGS[1]
comb = collect(combinations(1:5,2))
df = DataFrame(exonNumber = Int[], nbTemp = Int[], meanIdentity = Float64[], maxIdentity = Float64[],color = String[],  legend = String[], templatesnames = [], meanRmsd = Float64[])
fragments_identity_mean = Float64[]
fragments_identity_max = Float64[]
mean_rmsd = Float64[]
fragments_templates= Int64[]
fragments_templates_names = []
alltemplates = []
# Part for coverage and sequence identity, using the .pir with every templates
msa = read("$(query_name).pir", PIR)
msa_total = adjustreference!(copy(msa))
fragments_msa = get_10_num_from_int(length(msa_total[1,:]))
println(fragments_msa)

fragments_templates, fragments_identity_mean, fragments_identity_max, fragments_templates_names, alltemplates = get_infos_from_msa(fragments_msa, msa_total)
mean_rmsd = pairwise_RMSD_calculation_by_fragments(fragments_msa, query_name)
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
df3 =  DataFrame(numberOfDeciles = [], nbOfTemplates = [])
for i in 2:length(a)
    push!(df3, (i, length(a[i].nb)))
end
# Put everything in the dataframe
for i in 1:10
    if(fragments_identity_max[i] < 25)
        push!(df, (i, fragments_templates[i], fragments_identity_mean[i], fragments_identity_max[i], "red","0-24.99%", fragments_templates_names[i], mean_rmsd[i]))
    elseif ((fragments_identity_max[i] >=25) && (fragments_identity_max[i] <50))
        push!(df, (i, fragments_templates[i], fragments_identity_mean[i], fragments_identity_max[i], "orange","25-49.99%", fragments_templates_names[i], mean_rmsd[i]))
    elseif ((fragments_identity_max[i] >=50) && (fragments_identity_max[i] <75))
        push!(df, (i, fragments_templates[i], fragments_identity_mean[i], fragments_identity_max[i], "blue","50-74.99%", fragments_templates_names[i], mean_rmsd[i]))
    else
        push!(df, (i, fragments_templates[i], fragments_identity_mean[i], fragments_identity_max[i], "green","75-100%", fragments_templates_names[i], mean_rmsd[i]))
    end
end
# creation of every plots and saving it in <name query>.pdf
p1 = @df df bar(:exonNumber, :nbTemp,color = :color, group = :legend, bar_width=0.9,
    title = "Number of templates by decile",
    xlabel = "decile",
    xlim = (0,15),
    xticks = 0:1:10,
    ylabel="number of templates",
    ylim = (0,maximum(:nbTemp)+1),
    yticks = 0:5:maximum(:nbTemp),
    legend = :topright)
p2 = @df df scatter(:meanIdentity, :meanRmsd, group = :exonNumber,
    title = "mean RMSD by mean % identity, by deciles",
    xlabel = "mean % identity",
    xlim = (-1,maximum(:meanIdentity)+5),
    xticks = -1:5:maximum(:meanIdentity),
    ylabel="mean RMSD",
    ylim = (-1,maximum(:meanRmsd)+1),
    yticks = -1:0.5:maximum(:meanRmsd),
    legend = :bottomleft)

p3 = @df df scatter(:exonNumber, :meanRmsd, group = :nbTemp,
    title = "mean RMSD by deciles, grouped by number of templates",
    xlabel = "deciles",
    xlim = (-1,13),
    xticks = -1:1:10,
    ylabel="mean RMSD",
    ylim = (-1,maximum(:meanRmsd)+1),
    yticks = -1:0.5:maximum(:meanRmsd),
    legend = :topright)
p4 = @df df3 bar(:numberOfDeciles, :nbOfTemplates,
    title = "",
    xlabel = "number of deciles covered",
    ylabel="number of templates",
    xlim = (0, maximum(:numberOfDeciles)+1),
    ylim = (0,maximum(:nbOfTemplates)+1),
    xticks = 0:1:maximum(:numberOfDeciles),
    yticks = 0:1:maximum(:nbOfTemplates),
    legend=false)
p = plot(p1,p4,p2,p3,layout=(2,2), size = (1500,1500))
savefig(p, "$(query_name)_plots.pdf")
CSV.write("$(query_name).csv", df)
# remove every pdb.gz files downloaded so that modeller doesn't try to use them
# instead of the cif files
files = readdir()
for i in files
    if endswith(i, "pdb.gz")
        rm(i)
    end
end
