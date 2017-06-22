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
Hop back a level in our file structure.

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
multiple_join_paired_ends.py -i SKW1 -o SKW1_joined

multiple_join_paired_ends.py -i SK97 -o SK97_joined
```

Next we want to remove any reads that could not be joined (from both forward and reverse reads). We'll do this using a shell command similar to the one used above to do this for each of our folders. We're targeting files ending in 'fastqjoin.un1.fastq'  and 'fastqjoin.un2.fastq' like so.

```
find SKW1_joined -type f -name 'fastqjoin.un2.fastq' -delete
find SKW1_joined -type f -name 'fastqjoin.un1.fastq' -delete
find SK97_joined -type f -name 'fastqjoin.un2.fastq' -delete
find SK97_joined -type f -name 'fastqjoin.un1.fastq' -delete
```

With our reads joined and unjoined reads purged we want to demultiplex the files which trims and quality filters our sequuences leaving us with a fasta file for each gene. We'll use "multiple_split_libraries_fastq.py". 

```
multiple_split_libraries_fastq.py -i SKW1_joined -o SKW1_split -m sampleid_by_file --include_input_dir_path -p demult_param.txt

multiple_split_libraries_fastq.py -i SK97_joined -o SK97_split -m sampleid_by_file --include_input_dir_path -p demult_param.txt
```

Where '-i' is our directory with the reads in it, '-o' is a output directory, '-m' is telling QIIME to use the sampleid (in lieu of a barcode) to figure out sample names, '-p' is a parameter file that tells QIIME to use a quality score threshold of 30 (1 in 1000 chance the base isn't the correct base), and '--include_input_dir_path' which uses the directory name for the sample name.

Now that our samples are demultiplexed we can proceed with downstream analyses. 
