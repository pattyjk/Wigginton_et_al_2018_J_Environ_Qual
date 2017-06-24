#N-Cycling_Wastewater_Wigginton

## Moving things around
All the work flow of this project is done on a Linux-based machine using QIIME 1.9.1 (thus far). 

First, need to move all the files into the same folder so we can join the paired end reads. We'll use the find command paired with the copy command. The basic structure of this command is:

find -type f -name "*.gz" -exec cp {} /directory_where_files_are_going\;

The '*' is a wildcard in Linux so basically anything with a *.gz extension are found here. So we'll navigate to the directory (folder) where the files are located and execute a handful of commands. First, we'll create a directory for each set of sequences with 'mkdir'.  

```
mkdir SK97

mkdir SKW1
```
Next we'll change into the directory where the files are and then copy them to our created folders.

```
cd SKW1-96-36782746-003

find -type f -name "*.gz" -exec cp {} /Desktop/N-Cycling_wastewater_wigginton/SKW1\;
```
`Hop back a level in our file structure.

```
cd ..
```
And change into one of our folders with the reads in it.

```
cd SK97-186-37112077
```
Find and copy the files of interest. 

```
find -type f -name "*.gz" -exec cp {} /Desktop/N-Cycling_wastewater_wigginton/SK97\;
```
Hop back a level in our file structure.

```
cd ..
```
## Joining Paired Ends and Demultiplexing

Now that our files are moved we can join our paired ends with 'multiple_join_paired_ends.py' in QIIME. We'll run 2 commands to accomlpish this where '-i' is our folder with our reads in it and '-o' is a output directory for our joined reads. 

```
mmultiple_join_paired_ends.py -i SKW1 -o SKW1_joined

multiple_join_paired_ends.py -i SK97 -o SK97_joined
```

Next we want to remove any reads that could not be joined (from both forward and reverse reads). We'll do this using a shell command similar to the one used above to do this for each of our folders. We're targeting files ending in 'fastqjoin.un1.fastq'  and 'fastqjoin.un2.fastq' like so.

```
find SKW1_joined -type f -name 'fastqjoin.un2.fastq' -delete
find SKW1_joined -type f -name 'fastqjoin.un1.fastq' -delete
ofind SK97_joined -type f -name 'fastqjoin.un2.fastq' -delete
find SK97_joined -type f -name 'fastqjoin.un1.fastq' -delete
```

With our reads joined and unjoined reads purged we want to demultiplex the files which trims and quality filters our sequuences leaving us with a fasta file for each gene. We'll use "multiple_split_libraries_fastq.py". 

```
multiple_split_libraries_fastq.py -i SKW1_joined -o SKW1_split_nos -m sampleid_by_file --include_input_dir_path -p demult_param.txt

multiple_split_libraries_fastq.py -i SK97_joined -o SK97_split_amo -m sampleid_by_file --include_input_dir_path -p demult_param.txt
i```

Where '-i' is our directory with the reads in it, '-o' is a output directory, '-m' is telling QIIME to use the sampleid (in lieu of a barcode) to figure out sample names, '-p' is a parameter file that tells QIIME to use a quality score threshold of 30 (1 in 1000 chance the base isn't the correct base), and '--include_input_dir_path' which uses the directory name for the sample name.

Now that our samples are demultiplexed we can proceed with downstream analyses. 

## Chimera detection

As with any amplicon-based sequencing project we need to check our sequences for chimeric sequences that derive from PCR errors. We're going to be using USEARCH in de novo mode to identify chimeric sequences and then discard them. De novo mode uses ab abundance-based approach to ID chimeras and is particularly suited for functional genes that typicall have meh databases (in relation to environmental sequences). 

b```
cd SKW1_split_nos
identify_chimeric_seqs.py -m usearch61 --suppress_usearch61_ref -i seqs.fna -o chimeric

#we can see how many chimeras we have by using the shell command 'grep'
grep -c "SKW" chimeric/chimeras.txt

#we'll take the output and filter out chimeric sequences. 
filter_fasta.py -n -o chimera_free_split.fna -f seqs.fna -s chimeric/chimeras.txt



cd ..

cd SK97_split_amo
identify_chimeric_seqs.py -m usearch61 --suppress_usearch61_ref -i skw.txt -o chimeric
grep -c "SKW" chimeric/chimeras.txt

filter_fasta.py -n -o chimera_free_split.fna -f seqs.fna -s chimeric/chimeras.txt

cd ..

```This can be sped up (parallelized) by developing a refernce database from NCBI or FunGene.

## Picking operational taxonomic units (OTUs), picking a rep set, & making an OTU table

To collapse our amplicons to species level units we'll cluster the reads into OTUs. We will use 85% sequence identity for amoA (Pester et al. 2011, Envrion. Microbiol.) and 90% identity for nosZ (Bowen et al. 2013, Front. Micro.; this reference is for nirS but sample principal applies).

To pick OTUs we'll use the 'pick_otus.py' command paired with the OTU picking method swarm (Mahe et al. 2013, PeerJ). This is the most robust and accurate de novo OTU picking method.

```
cd SK97_split_amo
pick_otus.py -m swarm -i chimera_free_split.fna -o swarm_otus -s 0.85 
cd ..

cd SKW1_split_nos
pick_otus.py -m swarm -i chimera_free_split.fna -o swarm_otus -s 0.90
cd ..
```
Now that we have our OTUs we want to create a OTU table (species by sample matrix) with the command 'make_otu_table.py' and pick a set of representative sequences with 'pick_rep_set.py'. 

```
cd SK97_split_amo

pick_rep_set.py -f chimera_free_split.fna -o amo_rep_set.fna -i swarm_otus/
make_otu_table.py -o amo_otu_table.biom -i swarm_otus/

#you can convert this table to a .txt file with 'biom convert' command, the nuts and bolts of this command varies depending on your version of the biom package (mine is older). This OTU table can be used in any statistical package (i.e. R). 

biom convert -i amo_otu_table.biom -o amo_otu_table.txt -b --table-type="OTU table"
biom summarize-table -i amo_otu_table.biom -o amo_sum.txt 

#see our sequencing depth
less nos_sum.txt

cd ..

cd SKW1_split_nos

pick_rep_set.py -f chimera_free_split.fna -o nos_rep_set.fna -i swarm_otus/
make_otu_table.py -o nos_otu_table.biom -i swarm_otus/chimera_free_split_otus.txt 
biom convert -i nos_otu_table.biom -o nos_otu_table.txt -b --table-type="OTU table"
biom summarize-table -i nos_otu_table.biom -o nos_sum.txt

less nos_sum.txt

#go back home
cd ..
```

With our rep set and OTU table we can perform downstream analyses (beta and alpha diversity and assigning taxonomy).

## Diversity

This can be done in any statistical package but we'll calculate alpha and beta diversity in QIIME because it is easy. 

```
cd SKW1_split_nos

single_rarefaction.py -i nos_otu_table.biom -o nos_rare.biom -d
#normalizing the OTU table to teh lowest sequencing depth
#-d should be the lowest sequencing depth from the 'biom summarize-table command'.

beta_diversity.py -m bray_curtis -i nos_rare.biom -o nos_beta
#calculate pairwise dissimilarity with bray curtis. To see all options run 'beta_diversity.py -s'.

principal_coordinates.py -i nos_beta -o nos_beta/coords.txt
#creates 2D coordinates for a principal coordinates analysis (PCoA).

cd ..

cd SK97_split_amo

single_rarefaction.py -i amo_otu_table.biom -o amo_rare.biom -d
#-d should be the lowest sequencing depth from the 'biom summarize-table command'.

beta_diversity.py -m bray_curtis -i amo_rare.biom -o amo_beta
#calculate pairwise dissimilarity with bray curtis. To see all options run 'beta_diversity.py -s'.

principal_coordinates.py -i amo_beta -o amo_beta/coords.txt
#creates 2D coordinates for a principal coordinates analysis (PCoA).

cd ..
``` 

