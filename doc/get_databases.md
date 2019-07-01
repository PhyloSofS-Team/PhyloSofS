# Download the databases for HH-suite and modeller
-----
### These databases are needed for the structural part of PhyloSofS and are put after the arguments --hhdb, --structdb and --allpdb

Please make sure to have at least 300Go of free disk space before downloading the databases.

The first 2 databases can be found on the [hhsuite documentation](https://github.com/soedinglab/hh-suite) and can be found [here](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/)

The HHDB are archive called Uniprot/Uniclust  
The STRUCTDB files are all named pdb70...  
You need to decompress them.

The ALLPDB database consists of downloading every PDB contained in PDB. The steps needed to create this database are explained below:

First, download every .cif files from PDB. This can be done with [rsync](https://doc.ubuntu-fr.org/rsync). The command for rsync was taken from [the pdb website](https://www.rcsb.org/pages/download/ftp)

```
rsync -rlpt -v -z --delete rsync.ebi.ac.uk::pub/databases/pdb/data/structures/divided/mmCIF/ ./mmCIF
```
There is other way to download them if you don't want to use rsync, but they takes a lot more time. You can for example create a liste of every PDB ID from the PDB report section [list of all PDB entries as of 16/04/2019](http://www.rcsb.org/pdb/resultsV2/sids.jsp?qrid=2042F96E), then loop through it and wget every file.
```
#!/bin/bash
for i in `cat $1`;
do
    wget https://files.rcsb.org/download/$i.cif
	sleep .1
done
```

-----
## If you used rsync
Once you have downloaded the cif files with rsync, you need to decompress them and put them in the same folder
```
for i in `ls`;
do
    for j in `ls $i`;
    do 
        gunzip ./$i/$j && mv  $i/${j%???} ./../allpdb/;
    done;
done
```

-----
## If you used wget
If you used wget, no need to decompress and move the files, as they should be in the same folder already. You need to change the names of the files to lowercase. This can be done with [mmv](https://ss64.com/bash/mmv.html)
```
 mmv \*.cif\* \#l1.cif
```
(to check the output of the command before launching it, please add -n to the mmv command)


----
## use with PhyloSofS
-----
To use these databases with phylosofs, put them after the arguments:
    * --db  
    * --structdb  
    * --allpdb  
Be careful, for HHDB and STRUCTDB, you need to provide the path for the database, but also the basename of the files in the database.

The database pdb70_from_mmcif_latest contains files starting with pdb70. To use it with PhyloSofS, the argument should be:
> --structdb path_to_the_folder/pdb70

It is the same with HHDB. For ALLPDB, you need to give the path to the folder only, without the ending "/".
