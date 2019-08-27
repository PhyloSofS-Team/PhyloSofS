# New version ?
# The fact that we create an aligignment for each exon mess up modeller and does weird things (the hub of CAMK2B is awful)
# but if we add everything on the same line everything is nice.
# the could do this:
# if df.nb_df is the same as the last iteration, append the line, else, start a new one

using BioAlignments
using CSV
using MIToS.MSA
using BioStructures


transcript_name = ARGS[1]

function get_seq_id(d, chain)
    seq_id = Int64[]
    for (i,v) in enumerate(d["_atom_site.label_atom_id"])
        if v =="CA"
            if d["_atom_site.label_asym_id"][i] == chain
                push!(seq_id, parse(Int64, d["_atom_site.label_seq_id"][i]))
            end
        end
    end
    # unique is there to prevent any issues caused by the field ["_atom_site.label_alt_id"]
    # sometimes, the CA is present 2 times for the same number, but the field ["_atom_site.label_alt_id"] differs
    # causing the number to be pushed 2 times
    unique(seq_id)
end

function get_index(seq_id, value)
    for (i,v) in enumerate(seq_id)
        if v == value
            return i
        end
    end
end


open("$(transcript_name)_reconstructed.pir", "w") do f
    df_number = -1


    df = CSV.read("$(transcript_name)_reconstructed.csv", copycols=true)
# println(df)


    open("$(transcript_name)_annotated.pir","r") do e
        global seq_lines = readlines(e)
    end


# we may need 2 loops :
# 1 that appends the query sequence
# 1 that writes each exons from the pir files on each lines

    global queryseq = ""
    global length_gaps = 0
    for i in 1:size(df, 1)
        # if there is no templates found (first_pos_pir == -1; then just add the fasta sequence)
        # add the size of the exon to last_pos_pir
        if df.first_pos_pir[i] ==-1
#         println("not found")
            fasta_start = findfirst(df.symbol[i],seq_lines[2]).start
            fasta_end = findlast(df.symbol[i], seq_lines[2]).stop
#         println(seq_lines[3][fasta_start:fasta_end])
            df.last_pos_pir[i] = fasta_end-fasta_start
            global queryseq = queryseq * seq_lines[3][fasta_start:fasta_end]
        else
            global pir_msa = read("$(transcript_name)_$(df.nb_df[i]).pir", PIR, deletefullgaps=false)
        # the position from the pir are the good ones
        # !!! keep the gaps, we need them, do no do adjustreference!
            exon_qseq = pir_msa[1:1,df.first_pos_pir[i]:df.last_pos_pir[i]]
            global queryseq = queryseq * (join(exon_qseq))
        end
    end
#here, get the query_sequence name + change the header to get the good positions
    query_name = sequencenames(pir_msa)[1]
    query_length = length([ i for i in queryseq  if i != '-'])
    query_header = split(getannotsequence(pir_msa, query_name, "Title"), ":")
    query_header[3] = "1"
    query_header[5] = string(query_length)



    write(f, ">P1;" * split(query_name, "_")[1] * "_" * split(query_name, "_")[2] * "\n")
    write(f, join(query_header, ":") * "\n")
    write(f, queryseq * '*')
    write(f, "\n")

    global templateseq = ""
    global name = ""
    global length_gaps = 0
# global template_header = ""
    first_pos_array = Int64[]
    last_pos_array = Int64[]
    for i in 1:size(df,1)
        if df.first_pos_pir[i] == -1
            global length_gaps = length_gaps + df.last_pos_pir[i]
        elseif df.nb_df[i] != df_number
            if !isempty(templateseq)
#             println("should print the sequence and start a new one")
                template_name = split(name[1], "_")[1]
                template_chain = split(name[1], "_")[2]
                d = MMCIFDict("$(template_name).cif")
                seq_id = get_seq_id(d, template_chain)
                if length(replace(templateseq, "-" => "")) == 0
                    println("template $(template_name) should not be the one taken, it is only gaps")
                else
#             println(seq_id)
#             println(get_index(seq_id, minimum(first_pos_array)))
#             println(length(replace(templateseq, "-" => "")))
                    template_header[5] = string(seq_id[get_index(seq_id, minimum(first_pos_array)) + length(replace(templateseq, "-"=> ""))-1])
                    end_gaps = length(queryseq) - (length_gaps + length(templateseq))
                    template_header[3] = string(minimum(first_pos_array))
                    write(f,">P1;" * name[1] * "\n")
                    write(f,join(template_header, ":")* "\n")
                    write(f,"-" ^ length_gaps * templateseq * "-" ^ end_gaps * "*")
                    write(f, "\n")
                end
                global length_gaps = length_gaps + length(templateseq)
                global templateseq = ""
                global pir_msa = read("$(transcript_name)_$(df.nb_df[i]).pir", PIR, deletefullgaps=false)
                fragment = pir_msa[2:2,df.first_pos_pir[i]:df.last_pos_pir[i]]
                full_seq = replace(join(getsequence(pir_msa, 2)), "-" => "")
                fragment_seq = replace(join(getsequence(fragment, 1)), "-" => "")
                name = sequencenames(fragment)
                first_pos_array = []
                last_pos_array = []
                pos_start = findfirst(fragment_seq, full_seq).start
                pos_end = findfirst(fragment_seq, full_seq).stop
                global template_header = split(getannotsequence(pir_msa, name[1], "Title"), ":")
                starting_position = parse(Int64, template_header[3])
                push!(first_pos_array, starting_position + pos_start-1)
                push!(last_pos_array, starting_position + pos_end-1)
                global templateseq = templateseq * join(fragment)
            else
#             println("get the right exons etc.")
                global pir_msa = read("$(transcript_name)_$(df.nb_df[i]).pir", PIR, deletefullgaps=false)
                fragment = pir_msa[2:2,df.first_pos_pir[i]:df.last_pos_pir[i]]
                full_seq = replace(join(getsequence(pir_msa, 2)), "-" => "")
                fragment_seq = replace(join(getsequence(fragment, 1)), "-" => "")
                name = sequencenames(fragment)
                pos_start = findfirst(fragment_seq, full_seq).start
                pos_end = findfirst(fragment_seq, full_seq).stop
                template_header = split(getannotsequence(pir_msa, name[1], "Title"), ":")
                starting_position = parse(Int64, template_header[3])
                push!(first_pos_array, starting_position + pos_start-1)
                push!(last_pos_array, starting_position + pos_end-1)
                global templateseq = templateseq * join(fragment)
            end
            df_number = df.nb_df[i]
        else
#         println("get the right exons etc.")
            global pir_msa = read("$(transcript_name)_$(df.nb_df[i]).pir", PIR, deletefullgaps=false)
            fragment = pir_msa[2:2,df.first_pos_pir[i]:df.last_pos_pir[i]]
            full_seq = replace(join(getsequence(pir_msa, 2)), "-" => "")
            fragment_seq = replace(join(getsequence(fragment, 1)), "-" => "")
            name = sequencenames(fragment)
            pos_start = findfirst(fragment_seq, full_seq).start
            pos_end = findfirst(fragment_seq, full_seq).stop
            template_header = split(getannotsequence(pir_msa, name[1], "Title"), ":")
            starting_position = parse(Int64, template_header[3])
            push!(first_pos_array, starting_position + pos_start-1)
            push!(last_pos_array, starting_position + pos_end-1)
            global templateseq = templateseq * join(fragment)
        end
    end
    template_name = split(name[1], "_")[1]
    template_chain = split(name[1], "_")[2]
    d = MMCIFDict("$(template_name).cif")
    seq_id = get_seq_id(d, template_chain)
# println(seq_id)
# println(get_index(seq_id, minimum(first_pos_array)))
# println(length(replace(templateseq, "-" => "")))
    if length(replace(templateseq, "-" => "")) == 0
            println("template $(template_name) should not be the one taken, it is only gaps")
    else
        write(f,">P1;" * name[1]* "\n")
        template_header[3] = string(minimum(first_pos_array))
        template_header[5] = string(seq_id[get_index(seq_id, minimum(first_pos_array)) + length(replace(templateseq, "-"=> ""))-1])
        end_gaps = length(queryseq) - (length_gaps + length(templateseq))
        write(f,join(template_header, ":")* "\n")
        write(f,"-" ^ length_gaps * templateseq * "-" ^ end_gaps * "*")
    end
end
