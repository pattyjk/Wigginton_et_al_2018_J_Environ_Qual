# N-Cycling_Wastewater_Wigginton

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
```

Where '-i' is our directory with the reads in it, '-o' is a output directory, '-m' is telling QIIME to use the sampleid (in lieu of a barcode) to figure out sample names, '-p' is a parameter file that tells QIIME to use a quality score threshold of 30 (1 in 1000 chance the base isn't the correct base), and '--include_input_dir_path' which uses the directory name for the sample name.

Now that our samples are demultiplexed we can proceed with downstream analyses. 

## Chimera detection

As with any amplicon-based sequencing project we need to check our sequences for chimeric sequences that derive from PCR errors. We're going to be using USEARCH in de novo mode to identify chimeric sequences and then discard them. De novo mode uses ab abundance-based approach to ID chimeras and is particularly suited for functional genes that typicall have meh databases (in relation to environmental sequences). 

```
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
```

This can be sped up (parallelized) by developing a refernce database from NCBI or FunGene.The final consideration we should take into account if framshift mutations. We can take care of these with Framebot from the RDPstaff (http://fungene.cme.msu.edu/FunGenePipeline/framebot/form.spr) and can use their online utility (the command line program is a pain to get working). You then take the resulting file (in the framebot folder) and use that in downstream analyses. 

## Picking operational taxonomic units (OTUs), picking a rep set, & making an OTU table

To collapse our amplicons to species level units we'll cluster the reads into OTUs. We will use 85% sequence identity for amoA (Pester et al. 2011, Envrion. Microbiol.) and 90% identity for nosZ (Bowen et al. 2013, Front. Micro.; this reference is for nirS but sample principal applies).

To pick OTUs we'll use the 'pick_otus.py' command paired with the de novo OTU picking method swarm (Mahe et al. 2013, PeerJ). This is the most robust and accurate de novo OTU picking method. Reference-based OTU picking is common with 16S datasets, however given the large proportion of unknown taxa de novo OTU picking will likely be better for functional genes. 

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

#biom convert -i amo_otu_table.biom -o amo_otu_table.txt -b --table-type="OTU table"
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



## BLASTn

To get identify our rep set and make sure we have all amoA and nosZ sequneces we'll use BLAST+ as, unlike 16S sequences, the databases aren't the best (FunGene is alright). You have to downlaod and install the software from NCBI. The first thing you have to do is get the nt database. We can do that with a simple script.

```
#ignore these, these are MSU HPCC specific. 
#ssh dev-intel14
#module load RMBlast/2.2.28

mkdir nt
cd nt

update_blastdb.pl nt
#downloads a 48 file database.

tar -xf *.gz
#extracts all the files to the current directory (nt)
```

We also need to get a list of NCBI identifiers to ignore as, for this round, we only want known cultured denitrifiers. We can use this link:

https://www.ncbi.nlm.nih.gov/nuccore/?term=%22environmental+samples%22%5Borganism%5D+OR+metagenomes%5Borgn%5D+OR+%22Unidentified%22+OR+%22clone%22

Then: Send to > File > Format > GI list

Now we can BLAST our rep set. This will return a tab delimited file that has the OTU the closested cultured rep and the percent identity to that cultured organism. 

```
blastn -query nos_rep_set.fna -max_target_seqs 1 -outfmt "6 qseqid sacc stitle pident evalue" -out nos_results_out -negative_gilist sequence.gi -db nt

blastn -query amo_rep_set.fna -max_target_seqs 1 -outfmt "6 qseqid sacc stitle pident evalue" -out amo_results_out -negative_gilist sequence.gi -db nt
```

The most important part of BLAST is ensuring we have the gene of interest and only the gene of interest. For amoA we're looking for Nitrosomonas, Nitrosococcus, Nitrospira, and ppossibly Methylobacterium (pmoA-methane oxidizing gene- is very similar to amoA). For nosZ we're looking for things like Pseudomonas, Bradyrhizobium (any rhizobia), Cheatococcus, and Dechloromonas among many others (nosZ is very widespread). The nosZ sequnences look good, however the amoA sequences give us some pause as many are hitting weird things like Streptococcus. 

We'll blast this one more time and get the closest clone or envrionmental sequence (we don't technically need it but it's helpful for inferring drivers). This will require another (bigger) list of NBCI ascession numbers:

https://www.ncbi.nlm.nih.gov/nuccore/?term=plasmid%5Borgn%5D+OR+%22genome%22+OR+%22complete%22+OR+%22strain%22+OR+%22ATCC%22

```
blastn -query nos_rep_set.fna -max_target_seqs 1 -outfmt "6 qseqid sacc stitle pident evalue" -out nos_results_out_envir -negative_gilist cultured_gi.gi -db nt
#cultured_gi.gi is our file from the last link.

#rerun later. 
#blastn -query amo_rep_set.fna -max_target_seqs 1 -outfmt "6 qseqid sacc stitle pident evalue" -out amo_results_out_envir -negative_gilist cultured_gi.gi -db nt
```
## What's wrong with amoA?

Something isn't quite right with the amoA data. There's several ways to go about fixing this but we'll do this the easiest way in QIIME. 

First, we're going to download the AOB seed sequences from FunGene: http://fungene.cme.msu.edu/hmm_detail.spr?hmm_id=325

Select Seed Sequences > Begin Analysis > Nucleotide Download > Download

Next we want to do a closed reference OTU pick with our AOB seed sequences as a reference database.

```
pick_closed_reference_otus.py -i chimera_free_split2.fna -o closed_pick -r fungene_9.1_amoA_AOB_38_aligned_nucleotide_seqs.fa -p closed_param.txt 
```

This produces a single folder (closed_pick) and it will make a list of OTUs that hit the database within 75% similarity and those that failed to hit the database within 80% similarity. We can create an OTU table, summarize it, and create a rep set to see what we have for depth and to BLAST a sequence or two at NCBI (https://blast.ncbi.nlm.nih.gov/Blast.cgi).

```
make_otu_table.py -i closed_pick/chimera_free_split2_otus.txt -o closed_pick/otu_table.biom
biom summarize-table -i closed_pick/otu_table.biom -o closed_pick/otu_sum.txt
pick_rep_set.py -i closed_pick/chimera_free_split2_otus.txt -o closed_pick/rep_set.fna -f chimera_free_split2.fna
count_seqs.py -i closed_pick/rep_set.fna -o closed_pick/count.txt
```

We see we've lost a some sequences and we can see that our standard deviation on our sequence count (count.txt) is large (over 100bp). In the rep set that we have several sequences that are tiny <~50bp in length). However when we hand BLAST a few long sequences they all come back as amoA or amoA-containing organisms. When we blast the short sequences they come back as the odd ball sequences (i.e. Steptomyces genome). 

From here we'll take our OTU file and filter our fasta file like we did for chimera filtering except we'll be using the failures files from our closed reference OTU pick. Further, we'll purge all sequences <200bp in length with a shell script using the 'awk' command. We could have scrapped our amoA data and gone back to demultiplexing and thrown out sequences <200bp in length as well (though, its supposed to do that by default). 

```
filter_fasta.py -n -o chimera_free_split2_filtered.fna -f chimera_free_split2.fna -s closed_pick/uclust_ref_picked_otus/chimera_free_split2_failures.txt

awk '!/^>/ { next } { getline seq } length(seq) >= 200 { print $0 "\n" seq }' chimera_free_split2_filtered.fna > chimera_free_split2_filtered2.fna
#removes any sequence <200bp in length. Can be modified by changing the number or the widget '>' to exclude sequences of any length.

#chimera_free_split2_filtered2.fna is the file we'll throw into OTU picking
```

Now we can repeat the processing we did above (de novo OTU picking, rep set generation, BLASTn, and making an OTU table)

```
pick_otus.py -i chimera_free_split2_filtered2.fna -m swarm -o swarm_otus2 -s 0.85
pick_rep_set.py -f chimera_free_split2_filtered2.fna -o amo_rep_set2.fna -i swarm_otus2/chimera_free_split2_filtered2_otus.txt
make_otu_table.py -o amo_otu_table2.biom -i swarm_otus2/chimera_free_split2_filtered2_otus.txt
biom convert -b -i amo_otu_table2.biom -o amo_otu_table2.txt --table-type='OTU table'
biom summarize-table -i amo_otu_table2.biom -o amo_otu_table2_sum.txt

blastn -query amo_rep_set2.fna -max_target_seqs 1 -outfmt "6 qseqid sacc stitle pident evalue" -out amo_results_out2 -negative_gilist sequence.gi -db ntcp amo_

blastn -query amo_rep_set2.fna -max_target_seqs 1 -outfmt "6 qseqid sacc stitle pident evalue" -out amo_results_out_envir2 -negative_gilist cultured_gi.gi -db nt
```

## Mapping file
For our diversity to meaningful we need to put together a mapping file for each gene that has all the information (metadata) for each sample. QIIME requires a pretty specific format for their mapping files: http://qiime.org/documentation/file_formats.html#mapping-file-overview

We need the header to look like this in tab delimited format (CSV or any other kind won't work):

#SampleID BarcodeSequence LinkerPrimerSequence Description

All other information has to be squeezed as columns in between "LinkerPrimerSequence" and "Description". Additionall "BarcodeSequence" and "LinkerPrimerSequence" can be left blank for use because of the way we demultiplexed our samples (via the MiSeq versus by barcode in QIIME). Additionally, description can be anything you want (such as a sentence describing the sample or study).

Both the mapping files I have look like this:

```
head -1 amo_map.txt
#SampleID       BarcodeSequence LinkerPrimerSequence    Gene    Fragment_length _bp     Field_ID        Type_of _Wastewater Month_Collected  Sampling_Location       Description
```

## Diversity

This can be done in any statistical package but we'll calculate alpha and beta diversity in QIIME because it's easy. 

```
cd SKW1_split_nos
#change to our nosZ folder

single_rarefaction.py -i nos_otu_table.biom -o nos_rare.biom -d 2028
#normalizing the OTU table to the lowest sequencing depth, I chose 2028 (the second most shallow sample to maximize diversity, lowest was <50l0)
#-d should be the lowest sequencing depth from the 'biom summarize-table command'.

beta_diversity.py -m bray_curtis -i nos_rare.biom -o nos_beta
#calculate pairwise dissimilarity with bray curtis. To see all options run 'beta_diversity.py -s'.

principal_coordinates.py -i nos_beta -o nos_beta/coords.txt
#creates 2D coordinates for a principal coordinates analysis (PCoA).

make_2d_plots.py -m nos_map.txt -i nos_beta/coords.txt -o nos_beta/2D_plot
#makes a 2D PCoA plot. There's options for a 3D plot using emperor if you so desired. 

cd ..
#jump back a level.

cd SK97_split_amo
#change into our amoA directory.

single_rarefaction.py -i amo_otu_table2.biom -o amo_rare.biom -d 1431q
#-d should be the lowest sequencing depth from the 'biom summarize-table command'. I threw out the poorly sequenced sample here (depth=200) as its not sufficient for describing diversity here. 

beta_diversity.py -m bray_curtis -i amo_rare.biom -o amo_beta
#calculate pairwise dissimilarity with bray curtis. To see all options run 'beta_diversity.py -s'.

principal_coordinates.py -i amo_beta -o amo_beta/coords.txt
#creates 2D coordinates for a principal coordinates analysis (PCoA).

make_2d_plots.py -m amo_map.txt -o amo_beta/2d_plot -i amo_beta/coords.txt
#makes a 2D PCoA plot. There's options for a 3D plot using emperor if you so desired. 

cd ..
``` 
We have our first appreciable piece of data! We ca open the folder 2D_plot and click the HTML file so see the ordination of our Bray-Curtis similairty values. We have a few overarching patterns here that we can come back to in a bit.

Next we want to calculate our alpha diversity. We first need to re-rarefy our data (thousands of times, instead of just a single type) becasue alpha diversity is more sensitive to errors associated with rarefaction. So by doing it a bunch of times we can reduce that error. 

I am going to paralleize this script to take advantage of the HPCC at MSU but the guts of the script  will be the same. The only funky thing here is the '-m' option which is the start point of rarefaction (I always set it to 10). Also, "-O 16" calls 16 cores. 

```
parallel_multiple_rarefactions.py -i nos_otu_table.biom -o nos_mult_rare -m 10 -x 2028 -s 100 -n 1000 -O 16

parallel_multiple_rarefactions.py -i amo_otu_table.biom -o amo_mult_rare -m 10 -x 7627 -s 100 -n 1000 -O 16
```

This produces a folder for each gene with lots of OTU tables for each iteration of the rarefaction. I did 1k iterations with ~200-700 steps, which might be overkill (both computationally and diversity estimate wise) but I'd rather be sure I've converged on the true mean. However, if you have enough replication you can drop the number of iterations. 

Next, we want to calculate the alhpa diversity. Again, I'm going to parallize this but the guts are exactly the same as the resular script (alpha_diversity.py).

```
parallel_alpha_diversity.py -i nos_mult_rare/ -o nos_alpha -m observed_species,shannon,singles,goods_coverage,simpson_e -O 16

#going to calculate shannon diversity, observed species, the number of singeltons (i.e. OTUs present only once), coverage (people sometimes ask for this in papers), and eveness.

parallel_alpha_diversity.py -i amo_mult_rare/ -o amo_alpha -m observed_otus,shannon,singles,goods_coverage,pielou_e -O 16
```

Next we collate our alpha diversity results into a file and we'll append it to our mapping files (for ease of use in plotting in R or other packages). 

```

```



## Significance testing


