#!/usr/bin/env julia

using ArgParse
using BioStructures


function parse_commandline()
    settings = ArgParseSettings(description="""
    This script reads the `hhr` file returned by `hhsearch` using the pdb70
    database and download the needed PDB files into the indicated folder.
    """,)

    @add_arg_table settings begin
        "--input", "-i"
            help = "path to the `hhr` file returned by `hhsearch`."
            required = true
        "--pdb", "-p"
            help = """
            path where the pdb files in mmCIF format are going to be downloaded.
            """
            default = "."
    end

    return parse_args(settings)
end


"""
Return the list of pdb codes in the `hhr` file returned by `hhsearch`.

For example, a `hhr` file containing the following list in the header:

```

 No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
  1 2FE0_A small myristoylated pro  19.8      73  0.0012   23.4   0.0   16    1-16      6-21  (136)
  2 4MCB_B tRNA (guanine-N(1)-)-me  10.2   2E+02  0.0031   23.8   0.0   10   60-69     52-61  (246)

```

Is going to return:

```
10-element Array{String,1}:
 "2FE0"
 "4MCB"
```
"""
function get_pdb_list(hhr_file)
    open(hhr_file, "r") do file
        parse = false
        pdb_list = String[]
        for line in eachline(file)
            if occursin("No Hit", line)
                parse = true
                continue
            end
            if parse
                beginning = match(r"^\s*[0-9]+\s(\S+)_\S+\s", line)
                if beginning === nothing
                    return pdb_list
                end
                push!(pdb_list, beginning.captures[1])
            end
        end
        pdb_list
    end
end


"""
Download each pdb in pdb_list into the output_path in mmCIF format.
It will retry each pdb until 5 times if errors occur.
"""
function download_pdbs(pdb_list, output_path)
    for pdb in pdb_list
        retry(downloadpdb, delays=ExponentialBackOff(n=5))(
            pdb, file_format=MMCIF, pdb_dir=output_path
        )
    end
end


function main()
    parsed_args = parse_commandline()

    pdb_list = get_pdb_list(parsed_args["input"])

    output_path = abspath(parsed_args["pdb"])
    mkpath(output_path)

    download_pdbs(pdb_list, output_path)
end

main()
