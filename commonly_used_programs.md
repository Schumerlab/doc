# Commonly used programs

## Contents
1 [Example usage of commonly used programs](#example-usage-of-commonly-used-programs)
- 1.1 [Programs to overlap files & more](#programs-to-overlap-files--more)
- 1.2 [PLINK: the genomics Swiss army knife](#plink-the-genomics-swiss-army-knife)
- 1.3 [Running blast from the command line](#running-blast-from-the-command-line)
- 1.4 [Aligning fasta sequences](#aligning-fasta-sequences)
    - 1.4.1 [Aligning long sequences](#aligning-long-sequences)
    - 1.4.2 [Aligning shorter sequences, genes or proteins](#aligning-shorter-sequences-genes-or-proteins)
- 1.5 [Programs to adapter and quality trim fastq files](#programs-to-adapter-and-quality-trim-fastq-files)
- 1.6 [Programs to map reads and call variants](#programs-to-map-reads-and-call-variants)
- 1.7 [Programs to build phylogenies](#programs-to-build-phylogenies)
- 1.8 [Quick insight into species contributing to unknown hybrids](#quick-insight-into-species-contributing-to-unknown-hybrids)
- 1.9 [General use programs](#general-use-programs)

## Example usage of commonly used programs
### Programs to overlap files & more

bedtools is an incredibly useful program for overlap files of many common format types. Check out the bedtools
documentation here [[1] (https://web.archive.org/web/20230327074412/https://bedtools.readthedocs.io/en/latest/)]

Here are some examples we commonly use:

1) overlap a bed file with another bed file, or a large number of other possible formats (gff, gtf, vcf, etc):

`/home/groups/schumer/shared_bin/bedtools2/bin/intersectBed -a infile1.bed -b
infile2.vcf -wo > outfile.bed`

see bedtools documentation for the many options in overlapping files [[2] (https://web.archive.org/web/202303270
74412/https://bedtools.readthedocs.io/en/latest/)]

2) calculate coverage in a particular region of a bam file:

`/home/groups/schumer/shared_bin/bedtools2/bin/multiBamCov -bams file.bam -bed
region_of_interest.bed -q 30 > results`

3) convert a bam file to a fastq file:

`/home/groups/schumer/shared_bin/bedtools2/bin/bamToFastq -i file.bam -fq
file.fq`

### PLINK: the genomics Swiss army knife
If you want to calculate a basic pop-gen statistic or analysis, chances are plink can do it. Here is a small subset of
the things plink can do, see the plink homepage for details [[3] (https://web.archive.org/web/20230327074412/http
s://www.cog-genomics.org/plink/1.9/)]

Determine pairwise LD between SNPs or AIMs:

`/home/groups/schumer/shared_bin/plink infile-name-stem --ld-window-kb 5000 --
ld-window 20000 --r2 --hardy --hwe 0.001 --out outfil-name-stem-ld_decay --
ld-window-r2 0 --allow-extra-chr`

*Make sure to check out the command line parameters. If you want to calculate r or D instead, replace the r2 flag*

plink can convert between a large number of file formats, perform HWE filtering, LD pruning, identify mendelian
errors, and more

### Running blast from the command line

BLAST is a powerful algorithm for searching for alignments between sequences of varying levels of divergence. It
is most efficient for search for one of a few sequences in a large database (i.e. a genome).

On sherlock, load the blast modules:

`module load biology module load ncbi-blast+`

Format your database for blast search:

`makeblastdb -in genome.fa -dbtype nucl`

Run your blast job:

`blastn -db genome.fa -query my_gene.fa -task blastn -out
my_gene_blast_results -evalue 1e-20 -outfmt 6`

*Larger e-values will result in more alignments*

The unnamed columns in the output file are as follows:

```
query_sequence
reference_sequence
percent_identity
alignment_length
mismatches
gaps
query_start
query_end
reference_start
reference_end
e-value
bit_score
```

### Aligning fasta sequences

#### Aligning long sequences
For aligning long (i.e. ~chromosome length) sequences, mummer [[4] (https://web.archive.org/web/202303270744
12/https://mummer4.github.io/)] works great! Remember to read the documentation and tune your parameters
based on divergence between the sequences your aligning.

1) align sequences

`/home/groups/schumer/shared_bin/mummer-4.0.0rc1/nucmer fasta_chr1_species1.fa
fasta_chr1_species2.fa -l 100 -c 200 --prefix=chr-name-prefix`

2) mark and filter alignments

`/home/groups/schumer/shared_bin/mummer-4.0.0rc1/delta-filter -m chr-name-
prefix.delta > chr-name-prefix.delta.m`

3) plot results

`/home/groups/schumer/shared_bin/mummer-4.0.0rc1/mummerplot -large -layout
chr-name-prefix.delta.m --png -p chr-name-prefix`

if mummerplot in step 3 doesn't work (https://github.com/mummer4/mummer/issues/36), you can always plot the
*.delta.m file in R using this workflow (https://jmonlong.github.io/Hippocamplus/2017/09/19/mummerplots-with-
ggplot2/)

#### Aligning shorter sequences, genes or proteins
For aligning short sequences, particularly genes and proteins, clustal omega is very useful. **Note: clustal omega
will not work if the two sequences you're aligning aren't the same orientation.** The input file is a single fasta
containing all of the sequences you would like to align.

`/home/groups/schumer/shared_bin/clustalo-1.2.4-Ubuntu-x86_64 -i
infile_gene1_allspecies.fa -o clustal_align_gene1_allspecies.aln --
outfmt=clustal -v --force`

### Programs to adapter and quality trim fastq files

Trimgalore is a great program for trimming adapter sequences from your reads.

To trim single-end reads:

`/home/groups/schumer/shared_bin/TrimGalore/trim_galore --path_to_cutadapt
/home/groups/schumer/shared_bin/cutadapt --phred33 --quality 30 -e 0.001 --
stringency 1 --length 32 SE_read.fq.gz`

To trim paired-end reads:

`/home/groups/schumer/shared_bin/TrimGalore/trim_galore --path_to_cutadapt
/home/groups/schumer/shared_bin/cutadapt --phred33 --quality 30 -e 0.001 --
stringency 1 --length 32 --paired -retain_unpaired PE_read1.fq.gz
PE_read2.fq.gz`

### Programs to map reads and call variants

For genomic data, see Mapping and variant calling

For RNAseq data, see Mapping and RNA quantification

### Programs to build phylogenies

One of the most commonly use programs for phylogenetic reconstruction in the lab is RAxML
The RAxML program can be found here:

`/home/groups/schumer/shared_bin/raxmlHPC-PTHREADS`

1) Preparing fasta files:

RAxML requires as input a multi-fasta alignment file. RAxML does not accommodate polymorphic basepairs so
either mask them:

`perl -pi -e 's/[RYSWKM]/N/g' myfasta.fa`

or use a coin flip approach to convert them to non-IUPAC codes:

`perl coinflip_fasta_poly.pl infile.fa > outfile_coinflip.fa`

2) Optionally, trim aligned fasta files to variable sites to run RAxML quickly:

```
module load biology py-biopython
python Variable_sites_extractor.py -v myfasta.fa -o myvariablefa.fa
```

3) Run RAxML

See the RAxML documentation for specific applications (i.e. for coding sequences). For an alignment that is mainly
non-coding, here is a commonly used command to generate a maximum likelihood tree with 100 bootstraps:

`/home/groups/schumer/shared_bin/raxmlHPC-PTHREADS -f a -m GTRGAMMA -p 12345 -
x 12345 -s myvariablefa.fa -# 100 -n myresults_name`

You can add multi threading with the -T option but make sure your number of requested cpus matches or it will run
much more slowly than without multi threading

You can visualize your results with desktop applications like FigTree (https://web.archive.org/web/20230327074412/http://tree.bio.ed.ac.uk/software/figtree)

### Quick insight into species contributing to unknown hybrids

This can be done with sppIDer (a comparative mapping and and read depth analysis pipeline.)
(https://github.com/GLBRC/sppIDer)

This workflow has been tweaked to work on sherlock.

The first step is to create a combination reference genome of all the species of interest:

```
ml py-biopython/1.70_py27
ml bwa
ml samtools
python sppIDer_combineRefGenomes_sherlock.py --out combinedRef_XXX.fasta --key refGenomeKey.txt --trim ###
```

--out will be the final combined reference genome

--key is a key to the reference genomes to be combined, formatted as ref prefix/tab/full reference genome name

--trim will remove scaffolds shorter than the given size

To run the pipeline with fastq data of interest:
```
ml py-biopython/1.70_py27
ml bwa
ml samtools
ml R/3.6.1
python ~/scripts/sppIDer_sherlock.py --out ID_outputPrefix --ref
combinedRef_XXX.fasta --r1 ID.R1.fastq --r2 ID.R2.fastq
```

--out will be the prefix for all the produced files

--ref the combined reference genome produced above

--r1 read1 of individual

--r2 read2 of individual

--byGroup or --byBP two options for summarizing the depth (byGroup is faster, but dirtier)

### General use programs