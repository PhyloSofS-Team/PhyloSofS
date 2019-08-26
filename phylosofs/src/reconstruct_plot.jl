using DataFrames
using StatsPlots
using CSV
# using Plots


function get_number_of_df()
    open("number_of_loops.txt","r") do f
        global number = readlines(f)
    end
    global number = parse(Int,join(number))
    return number
end

function get_exons_symbols(query_name)
    df = DataFrame(symbol = Char[], length = Int[])
    open("$(query_name)_annotated.pir","r") do f
        global lines = readlines(f)
    end
    exons_symbol = split(lines[1])[3]
    return exons_symbol
end

function reconstruct_dataframe(query_name, exons, nb_of_dataframes)
    all_df = DataFrame(nb_df = Int[], exonNumber = Int[],symbol= String[], nbTemp = Int[], maxIdentity = Float64[],color = String[],  legend = String[], first_pos_pir = Int64[], last_pos_pir = Int64[])
    reconstructed_dataframe = DataFrame(exonNumber = Int[],symbol= String[], nbTemp = Int[], maxIdentity = Float64[],color = String[],  legend = String[], nb_df = Int64[], first_pos_pir = Int64[], last_pos_pir = Int64[])
    for i in 1:Int(nb_of_dataframes)
        global df = CSV.read("$(query_name)_$(i)_plot1.csv")
        for j in 1:size(df,1)
            push!(all_df, (i, df.exonNumber[j], string(df.symbol[j]), df.nbTemp[j], df.maxIdentity[j], df.color[j], df.legend[j], df.first_pos_pir[j], df.last_pos_pir[j]))
        end
    end
    subdf = groupby(all_df, :symbol, sort=false)
    for i in 1:length(subdf)
        row_to_take = findmax(subdf[i].maxIdentity)[2]
        # if every maxIdentity == 0, then take the one that has the most number of templates
        if findmax(subdf[i].maxIdentity)[1] ==0
            row_to_take = findmax(subdf[i].nbTemp)[2]
        end
        push!(reconstructed_dataframe, (i, subdf[i].symbol[row_to_take], subdf[i].nbTemp[row_to_take], subdf[i].maxIdentity[row_to_take], subdf[i].color[row_to_take], subdf[i].legend[row_to_take], subdf[i].nb_df[row_to_take], subdf[i].first_pos_pir[row_to_take], subdf[i].last_pos_pir[row_to_take]))
    end
    return reconstructed_dataframe
end


query_name = ARGS[1]
nb_of_dataframes = get_number_of_df()
exons = get_exons_symbols(query_name)
reconstructed_dataframe = reconstruct_dataframe(query_name, exons, nb_of_dataframes)

p1 = @df reconstructed_dataframe bar(:nbTemp,color = :color, group = :legend, bar_width=0.9,
        title = "Number of templates by exons",
        xlabel = "Exons",
        xlim = (0,size(reconstructed_dataframe,1)+1),
        xticks = (1:1:size(reconstructed_dataframe,1)+1, :symbol),
        ylabel="Number of templates",
        ylim = (0,maximum(:nbTemp)+1),
        yticks = 0:5:maximum(:nbTemp),
        legend = :topright)

savefig(p1, "$(query_name)_reconstructed_plots.pdf")
CSV.write("$(query_name)_reconstructed.csv", reconstructed_dataframe)