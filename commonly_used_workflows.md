# Commonly used workflows

Please add more workflows and pipelines as you develop them!

## Contents
1.  [Parsing Illumina data - Tn5 libraries](#parsing-illumina-data---tn5-libraries)
    - 1.1 [i5 barcode parsing](#i5-barcode-parsing)
    - 1.2 [i7 barcode generation](#i7-barcode-generation)
    - 1.3 [parse by i5 and i7 barcodes](#parse-by-i5-and-i7-barcodes)
    - 1.4 [Quick quality checks after parsing](#quick-quality-checks-after-parsing)
2. [Parsing Illumina data - Quail libraries or similar](#parsing-illumina-data---quail-libraries-or-similar)
    - 2.1 [generate a file with the individual id and index](#generate-a-file-with-the-individual-id-and-index)
    - 2.2 [parse by i7 index](#parse-by-i7-index)
3. [Ancestry analysis with AncestryHMM](#ancestry-analysis-with-ancestryhmm)
    - 3.1 [For first time use](#for-first-time-use)
    - 3.2 [For subsequent use](#for-subsequent-use)
4. [Simulating hybrid genomes with Simulate_hybrid_genomes](#simulating-hybrid-genomes-with-simulate_hybrid_genomes)
5. [Ancestry tsv files to admixture mapping input](#ancestry-tsv-files-to-admixture-mapping-input)
6. [Ancestry tsv files to rQTL input](#ancestry-tsv-files-to-rqtl-input)
7. [Ancestry tsv files to average ancestry in windows](#ancestry-tsv-files-to-average-ancestry-in-windows)
8. [Thin Ancestry tsv files by genetic distance](#thin-ancestry-tsv-files-by-genetic-distance)
9. [Ancestry tsv to identifying ancestry transitions](#ancestry-tsv-to-identifying-ancestry-transitions)
10. [Ancestry tsv files to admixture LD decay](#ancestry-tsv-files-to-admixture-ld-decay)
11. [Merging multiple ancestry tsv files run separately](#merging-multiple-ancestry-tsv-files-run-separately)
12. [Identifying minor parent deserts and islands](#identifying-minor-parent-deserts-and-islands)
13. [Mapping and variant calling with GATK](#mapping-and-variant-calling-with-gatk)
14. [Case/Control GWAS from low coverage data](#casecontrol-gwas-from-low-coverage-data)
15. [Pseudohaploid calls for GWAS or population structure analysis from low coverage data](#pseudohaploid-calls-for-gwas-or-population-structure-analysis-from-low-coverage-data)
16. [Mapping Binary & Continuous Traits with PLINK](#mapping-binary--continuous-traits-with-plink)
17. [g.vcf files to Admixtools input](#gvcf-files-to-admixtools-input)
    - 17.1 [Batch script to generate insnp files](#batch-script-to-generate-insnp-files)
18. [g.vcf files to ped input for use with PLINK](#gvcf-files-to-ped-input-for-use-with-plink)
19. [D-statistic analysis with Admixtools](#d-statistic-analysis-with-admixtools)
20. [F4 ratio tests from fasta files]( #f4-ratio-tests-from-fasta-files)
21. [F4 Ratio test with Dsuite](#f4-ratio-test-with-dsuite)
    - 21.1 [Dtrios](#dtrios)
    - 21.2 [Fbranch](#fbranch)
22. [g.vcf files to pseudo-fasta files](#gvcf-files-to-pseudo-fasta-files)
23. [Mapping RNAseq data and transcript quantification](#mapping-rnaseq-data-and-transcript-quantification)
    - 23.1 [Mapping to a transcriptome](#mapping-to-a-transcriptome)
    - 23.2 [Mapping to a genome](#mapping-to-a-genome)
    - 23.3 [Transcript quantification](#transcript-quantification)
    - 23.4 [Statistical analysis with DESeq2](#statistical-analysis-with-deseq2)
24. [Allele specific expression from RNAseq data](#allele-specific-expression-from-rnaseq-data)
    - 24.1 [For first time use](#for-first-time-use-1)
    - 24.2 [For subsequent use](#for-subsequent-use-1)
25. [Genome assembly from PacBio HiFi data](#genome-assembly-from-pacbio-hifi-data)
26. [Genome assembly with supernova - 10x data](#genome-assembly-with-supernova---10x-data)
27. [Motif and binding site predictions](#motif-and-binding-site-predictions)
28. [Liftovers using cactus alignments](#liftovers-using-cactus-alignments)
    - 28.1 [Other useful basic commands with halTools](#other-useful-basic-commands-with-haltools)
29. [Demographic Inference with ABC](#demographic-inference-with-abc)
    - 29.1 [Simulations setup](#simulations-setup)
    - 29.2 [Post simulations population specific inference](#post-simulations-population-specific-inference)
30. [Demographic Inference with PSMC](#demographic-inference-with-psmc)
    - 30.1 [plotting PSMC output](#plotting-psmc-output)
31. [Building a recombination map with LDhelmet](#building-a-recombination-map-with-ldhelmet)

## Parsing Illumina data - Tn5 libraries
Note: this content is also in dropbox in the guide called "parsing_tn5_data.txt"
Before you can use Tn5 data, it needs to be parsed by both i5 and i7 barcode. Check out this helpful document from Patrick Reilly for more background File:IlluminaParsing
v1-2.pdf
### i5 barcode parsing

First, parse by i5 barcode:
1) generate a file with the plate id and i5 name \t the i5 sequence (see i5_indices.tsv or i5_indices_1-24_revcomp.tsv in Dropbox). For example:

```
i5-23_COACVI2018 GGCTCGAA
i5-2_CHAFV2018-ACUAV2018 TTGGCGTT
i5-3_ACUAV2018-CHAFI-II2018-ACUAVI2015 ATCAGCGC
i5-4_ACUAVI2015-CHAFXI2017-COACXI2017 TAAATAGG
```

- Note: whether you need the forward or reverse complement i5 sequences depends on which Illumina sequencer you used - currently NextSeq and HiSeq 4000 need the
reverse complement

2) Make sure to save this file as a windows formatted text file (or Tab delimited text)

3) copy the file to the server

4) convert the file to Unix text format

    `dos2unix [filename]`

Note: if you look at your file and see ^M characters, run the following as well:

    `mac2unix [filename]`

### i7 barcode generation
Note: you need to follow these steps once for every plate
1) open "Tn5_i7s_convert_plate_to_barcodefile.xls" from Dropbox
2) paste your plate layout in the space indicated
3) copy columns A & B which are automatically generated and paste using into a new excel document by selecting:
Edit->paste special->Values
4) change the file format to "windows formatted text". Save the file name as Platename_i7_barcodes
    - for example "COAC-V1_2018_CHAF-V-V1_2018_Tn5_data_July2018_i7_barcodes.txt"
5) spot check several cells spread throughout the document to ensure that the right sample name, i7 id, and i7 sequence have been matched up
6) copy the file to the server
7) convert the file to unix text format

    `dos2unix [filename]`

### parse by i5 and i7 barcodes
1) Submit a slurm job to run parsing by i5 barcode. The usage of this job is:
/home/groups/schumer/shared_bin/Lab_shared_scripts/divideConquerParser.sh [# read files] [List of FASTQs IN QUOTES]
[# cores to use] [i5 barcode file] [which position in FASTQ list is i5 index read]
example:

```
/home/groups/schumer/shared_bin/Lab_shared_scripts/divideConquerParser.sh 4 "TLMCI2015_COACVI2018-2017_CHAFV2018-
2017_ACUAV2018_ACUAVI2015_S0_R1_alllanes_combined.fastq.gz TLMCI2015_COACVI2018-2017_CHAFV2018-
2017_ACUAV2018_ACUAVI2015_S0_R2_alllanes_combined.fastq.gz TLMCI2015_COACVI2018-2017_CHAFV2018-
2017_ACUAV2018_ACUAVI2015_S0_I1_alllanes_combined.fastq.gz TLMCI2015_COACVI2018-2017_CHAFV2018-
2017_ACUAV2018_ACUAVI2015_S0_I2_alllanes_combined.fastq.gz" 10 i5_library_COACVI2018_CHAFV2018_ACUAV2018_ACUAVI2015
4 1
```

- Note: parsing by i5 is not necessary if you have only used a single i5 in the sequencing run, you can move directly to parsing by i7
- Note: make sure the number of cores in your slurm job matches the number requested in the command line. For example, for the above an appropriate slurm header
would be:

```
#!/bin/bash
#SBATCH --job-name=parse_i5
#SBATCH --time=24:00:00
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem=32000
```

2) After this job has finished, you can move on the parsing by i7 barcode. You will need to run one i7 job for *each* plate. For example, if the i5 barcode file had five lines, you will need to run five i7 parsing jobs.

    If you look at the output of the i5 parsing you will notice that each parsed set generated three files, so in this next round of parsing you will have three files instead of four but otherwise the usage is the same:

```
/home/groups/schumer/shared_bin/Lab_shared_scripts/divideConquerParser.sh [# read files] [List of FASTQs IN QUOTES]
[# cores to use] [i7 barcode file] [which position in FASTQ list is i7 index read]
```

```
/home/groups/schumer/shared_bin/Lab_shared_scripts/divideConquerParser.sh 3 "i5-3_ACUAV2018-CHAFI-II2018-
ACUAVI2015_read_1.fastq.gz i5-3_ACUAV2018-CHAFI-II2018-ACUAVI2015_read_2.fastq.gz i5-3_ACUAV2018-CHAFI-II2018-
ACUAVI2015_read_3.fastq.gz" 10 20180716_Tn5_library_ACUA_V_2018_CHAF_I_II_2018_ACUA_VI_2015_i7_barcodes.txt 3 1
```

### Quick quality checks after parsing
1) Check the Part*logs files
    - ensure that there are few reads in the control wells (<0.01% in a good library)
    - check that coverage is relatively even among individuals (few individuals <0.3% or >2%)
2) Check duplication levels
    - check duplication levels with picardtools for a few individuals per run
    - acceptable levels are <20% (0.2 in picard output)

```
module load biology
module load bwa

bwa mem -t 3 -M xma_washu_4.4.2-jhp_0.1_combined-unplaced-mito.fa COACVI1801_i7-13_read_1.fastq.gz COACVI1801_i7-
13_read_2.fastq.gz > COACVI1801.sam

java -jar /home/groups/schumer/shared_bin/SortSam.jar INPUT=COACVI1801.sam OUTPUT=COACVI1801.sorted.bam
SORT_ORDER=coordinate

java -jar /home/groups/schumer/shared_bin/BuildBamIndex.jar INPUT=COACVI1801.sorted.bam

java -jar /home/groups/schumer/shared_bin/MarkDuplicates.jar I=COACVI1801.sorted.bam O=COACVI1801.dedup.bam
CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=COACVI1801.metrics
```

## Parsing Illumina data - Quail libraries or similar
Quail libraries are parsed in a similar way, but only have a single index so only need a single matching index file.

#### generate a file with the individual id and index
This files needs to have the individual id and i7/FC-1 name \t the i7/FC-1 sequence (see Quail_FC1_index_sequence.xls). For example:

```
FC1-9_HUICXI17M01 GATCAG
FC1-10_HUICXI17M04 TAGCTT
FC1-12_HUICXI17JM06 CTTGTA
FC1-13_HUICXI17M02 AGTCAA
FC1-14_HUICXI17M03 AGTTCC
FC1-15_HUICXI17M05 ATGTCA
FC1-16_HUICXI17JM07 CCGTCC
```

#### parse by i7 index
Because these libraries are only parsed by i7 barcode, you can drop the second index read in the parsing
Usage:

```
/home/groups/schumer/shared_bin/Lab_shared_scripts/divideConquerParser.sh [# read files] [List of FASTQs IN QUOTES]
[# cores to use] [i7 barcode file] [which position in FASTQ list is i7 index read]
```

For example:

```
/home/groups/schumer/shared_bin/Lab_shared_scripts/divideConquerParser.sh 3
"Xcortezi_resequence_HUIC_sc_wt_November2018_S0_R1_allcombined.fastq.gz
Xcortezi_resequence_HUIC_sc_wt_November2018_S0_R2_allcombined.fastq.gz
Xcortezi_resequence_HUIC_sc_wt_November2018_S0_I1_allcombined.fastq.gz" 10 Xcortezi_spottedcaudal_November2018 3
```

### Ancestry analysis with AncestryHMM
Local ancestry inference is one of the most common workflows in the lab. The steps we use to run local ancestry inference for birchmanni x malinche are outlined below and
is also available in dropbox ("Shared_lab_resources/Common_commands_and_pipelines/Running_Ancestry_HMM_on_Sherlock.txt").
#### For first time use
1) Make a folder within your personal lab member folder to perform the analysis, e.g.:

```
mkdir Ancestry_HMM_runs
cd Ancestry_HMM_runs
```

2) link or copy the files you plan to use for analysis (genome assemblies and ancestry informative sites files). We keep up to date lab versions of these in shared_resources,
e.g:

```
cp /home/groups/schumer/shared_bin/shared_resources/ancestry_hmm/Base_ancestry_run_files_May2023/* ./
```

3) make a folder to link the raw reads files from OAK. For example:

```
mkdir reads_files
cd reads_files
ln -s /oak/stanford/groups/schumer/data/All_swordtail_low_coverage_Tn5_data/XcorXbir-F2-24-S255_XmalXbir-F2-23-24-
S210_XmalXbir-OMCA-V-23_XxipXcou-F2-F3-VI-X-23_Tn5_FEB2024/XmalXbir-F2* ./
```

4) make a reads list file that has read1 for each individual in column 1 and read2 for each individual in column 2. If you data has the typical format from admera you can do
this from your main analysis directory by running the following:

```
ls ./reads_files/*.R1.fastq.gz > read1_list
ls ./reads_files/*.R2.fastq.gz > read2_list
paste read1_list read2_list > combined_read_list
```

5) make a copy of the sample configuration file for your specific analysis

```
cp hmm_configuration_file.cfg hmm_configuration_file_runCHAF_samples2017.cfg
```

6) use a text editor or emacs to edit the parameters in this file (see parameter details in Dropbox guide, README, or at:
https://github.com/Schumerlab/Ancestry_HMM_pipeline)

7) load required packages and modules:

```
module load armadillo
module load biology
module load samtools
module load bcftools
module load py-pysam/0.14.1_py27
module load bwa
module load boost
module load R
export PATH="/home/groups/schumer/shared_bin/ngsutils/bin:$PATH"
export PATH="/home/groups/schumer/shared_bin/Ancestry_HMM/src:$PATH"
export PYTHONPATH=/home/groups/schumer/shared_bin:$PYTHONPATH
```

8) make a batch script for submission, e.g.:

```
#!/bin/bash
#SBATCH --job-name=run_hmm
#SBATCH --time=01:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32000
#SBATCH --mail-user=youremail@stanford.edu

perl /home/groups/schumer/shared_bin/ancestryinfer_July2022/Ancestry_HMM_parallel_v7.pl
hmm_configuration_file_runCHAF_samples2017.cfg
```

you can also include the module/export commands in your slurm script after the #SBATCH lines
See example here:
```
/home/groups/schumer/shared_bin/shared_resources/ancestry_hmm/Base_ancestry_run_files_May2023/run_hmm_example_new.sh
```

9) submit your job, e.g.:
`sbatch run_hmm.sh`

#### For subsequent use
1) edit your configuration files and job submit script as needed
2) export dependencies or make sure they are present in your slurm submission file

```
module load armadillo
module load biology
module load samtools
module load bcftools
module load py-pysam/0.14.1_py27
module load bwa
export PATH="/home/groups/schumer/shared_bin/ngsutils/bin:$PATH"
export PATH="/home/groups/schumer/shared_bin/Ancestry_HMM/src:$PATH"
export PYTHONPATH=/home/groups/schumer/shared_bin:$PYTHONPATH
```

3) submit your job, e.g.:
`sbatch run_hmm.sh`

## Simulating hybrid genomes with Simulate_hybrid_genomes
It's often useful to simulate admixed genomes to evaluate efficacy of ancestry calling under different scenarios.
1) For first time use:

```
module load perl
cpan Math::Random
```

Make sure the following programs are in your path:

```
fastahack
seqtk
wgsim
bedtools2
```
They are all available in:

`/home/groups/schumer/shared_bin`

You can either export these each time you use them:

`export PATH="/home/groups/schumer/shared_bin:$PATH"`
or put this export command in your .bashrc file
make a directory for the analysis:

```
mkdir Simulate_hybrid_genomes_runs
cd Simulate_hybrid_genomes_runs
ln -s /home/groups/schumer/shared_bin/Simulate_hybrid_genomes/* ./
cp hybrid_simulation_configuration.cfg myrun_hybrid_simulation_configuration.cfg
```

2) After step 1 or for subsequent use:

    Export dependencies:

    `export PATH="/home/groups/schumer/shared_bin:$PATH"`

    Edit the configuration file:

    `emacs myrun_hybrid_simulation_configuration.cfg`

3) Run the simulation
perl simulate_admixed_genomes_v2.pl myrun_hybrid_simulation_configuration.cfg
4) You can use these simulated reads as input into the Ancestry_HMM_pipeline
5) Afterwards you can summarize accuracy with the following script (part of the Simulate_hybrid_genomes GitHub):
perl post_hmm_accuracy_shell.pl ancestry-probs-par1_file ancestry-probs_par2_file path_to_simulation_reads_folder

## Ancestry tsv files to admixture mapping input
1) After running AncestryHMM, convert your ancestry tsv files to hard calls using parsetsv_to_genotypes_Dec2017_v2.pl

`perl parsetsv_to_genotypes_Dec2017_v2.pl ancestry-probs-par1_file ancestry-probs-par2_file
genotypes_outputfile_name`

- 1a) It's often a good idea to filter your tsv files to remove markers violating hardy-weingberg equilibrium. This script can be used to do that. This script is set up to run this
filter before you convert to genotypes format:

`perl determine_HWE_filter_markers.pl ancestry-par1 ancestry-par2 bonferonni_pval_thresh
path_to:transpose_nameout.pl`

*Writes a file named infile_deviating_markers with deviating markers from HWE, and filtered ancestry tsv files: ancestry-par1_HWE.tsv ancestry-par2_HWE.tsv*

- 1b) Optionally filter your genotypes file with a script for filtering genotypes file of redundant columns. The number in the command line corresponds to the number of markers that can differ between adjacent columns for the column to be retained:

`perl filter_identical_columns_threshold.pl genotypes_file num_markers_differentiation path_to:transpose_nameout.pl`

*Writes an output file with genotypes_file.identicalfilter.txt*

- 1c) Another option for filtering is to filter based on physical distance rather than redundancy. You may lose more information this way but it is unbiased with respect to
errors. Be careful to chose a reasonable threshold based on admixture LD decay in the specific population you are looking at.

`perl thin_genotypes_by_physical_distance.pl genotypes_file distance_thresh_bp path_to:transpose_nameout.pl`

*Writes an output file named genotypes_file_thinned_physical_dist.txt*

2) Generate a hybrid index file for these individuals with mixture proportions for each individual

`perl parsetsv_ancestry_v2.pl ancestry-probs-par1_file ancestry-probs-par2_file > hybrid_index_file`

3) Match your previously generated phenotypes file with genotypes and hybrid indices for each individual

`perl match_phenotypes_names_with_genotypes_and_index_file.pl phenotypes_file_name genotypes_file hybrid_index_file`

*This will generate output files appended with _matched_to_genotypes _matched_to_phenotypes _matched_to_hybrid_index*

4) Check file output carefully to make sure each has the appropriate number of lines and that individuals are listed in the same order in each file

5) By default the admixture mapping scripts test the likelihood of a model of phenotype~hybrid_index versus a model of phenotype~genotype_focal_site + hybrid_index. If you want to include more covariates (i.e. body size) you will need to modify the script.

For binomial traits (0/1) 

`Rscript perform_glm_admixture_mapping_v2_binomialtrait.R genotypes_file hybrid_index_file
phenotypes_file focal_column_number outfile_name_tag`

For continuous traits that follow a gaussian distribution

`Rscript perform_glm_admixture_mapping_v2_gaussian.R genotypes_file
hybrid_index_file phenotypes_file focal_column_number outfile_name_tag`

*The outfile name will be: genotypes_file_results_gaussian_v2_outfile_name_tag*

6) To plot your results (as long as you used the X. birchmanni 10x genome), convert the file format:

`perl convert_birchmanni10x_mapping_output_manhattan_plot_input.pl
genotypes_file_results_gaussian_v2_outfile_name_tag`

7) Then you can make a nice Rplot of these results using the adapted_qqman.R script
in R, source the adapted_qqman.R script:

`source("adapted_qqman.R")`

Load the data in R:

`data<-read.csv(file="genotypes_file_results_gaussian_v2_outfile_name_tag",sep="\t",head=FALSE)`

Reformat the data:
```
data_trim<-{}
data_trim$CHR<-data$chrom
data_trim$BP<-data$marker
data_trim$P<- data$likelihood.diff
data_trim<-as.data.frame(data_trim)
data_trim<-na.omit(data_trim)
```

Plot the data, e.g.:

`manhattan(data_trim,genomewideline=8,suggestiveline=6)`


## Ancestry tsv files to rQTL input
1) After running AncestryHMM, convert your ancestry tsv files to hard calls using parsetsv_to_genotypes_Dec2017_v2.pl

`perl parsetsv_to_genotypes_Dec2017_v2.pl ancestry-probs-par1_file ancestry-probs-par2_file
genotypes_outputfile_name`

- 1a) It's often a good idea to filter your tsv files to remove markers violating hardy-weingberg equilibrium. This script can be used to do that. This script is set up to run this
filter before you convert to genotypes format:

`perl determine_HWE_filter_markers.pl ancestry-par1 ancestry-par2 bonferonni_pval_thresh
path_to:transpose_nameout.pl`

*Writes a file named infile_deviating_markers with deviating markers from HWE, and filtered ancestry tsv files: ancestry-par1_HWE.tsv ancestry-par2_HWE.tsv*

2) Optionally filter your genotypes file with a script for filtering genotypes file of redundant columns. The number in the command line corresponds to the number of markers
that can differ between adjacent columns for the column to be retained:

`perl filter_identical_columns_threshold.pl genotypes_file num_markers_differentiation path_to:transpose_nameout.pl`

*Writes an output file with genotypes_file.identicalfilter.txt'*

- 2a) Another option for filtering is to filter based on physical distance rather than redundancy. You may lose more information this way but it is unbiased with respect to
errors:

`perl thin_genotypes_by_physical_distance.pl genotypes_file distance_thresh_bp path_to:transpose_nameout.pl`

*Writes an output file named genotypes_file_thinned_physical_dist.txt*

3) Generate a hybrid index file for these individuals with mixture proportions for each individual

`perl parsetsv_ancestry_v2.pl ancestry-probs-par1_file ancestry-probs-par2_file > hybrid_index_file`

4) Match your previously generated phenotypes file with genotypes and hybrid indices for each individual

`perl match_phenotypes_names_with_genotypes_and_index_file.pl phenotypes_file_name genotypes_file hybrid_index_file`

*This will generate output files appended with _matched_to_genotypes _matched_to_phenotypes _matched_to_hybrid_index*

5) Check file output carefully to make sure each has the appropriate number of lines and that individuals are listed in the same order in each file

6) Convert the genotypes to rQTL format

`perl genotypes_to_rqtl_April2019_v4.pl genotypes_file.identicalfilter.txt_matched_to_phenos
phenotypes_file_phenos_second_column`

*writes an output file named: infile.rqtl.csv*

7) Now you're ready to run rQTL! See rQTL documentation for details: http://www.rqtl.org/manual/qtl-manual.pdf
Use the following structure to read in this data format properly:

`data<-read.cross("csv",dir="./","infile.rqtl.csv",na.strings=c("NA"),estimate.map=FALSE)`

## Ancestry tsv files to average ancestry in windows
1) Generate an file with average ancestry at every site that has greater than a user-defined threshold number of individuals and user-defined posterior probability cutoff
(recommended >=0.9)

`perl calculate_avg_ancestry_by_site_tsv_files_v5.pl infilepar1 infilepar2 num_ind_thresh posterior_prob_thresh
path_to:transpose_nameout.pl > outfile`

- 1a) It's often a good idea to filter your tsv files to remove markers violating hardy-weingberg equilibrium. This script can be used to do that. This script is set up to run this
filter before you convert to genotypes format:

`perl determine_HWE_filter_markers.pl ancestry-par1 ancestry-par2 bonferonni_pval_thresh
path_to:transpose_nameout.pl`

*Writes a file named infile_deviating_markers with deviating markers from HWE, and filtered ancestry tsv files: ancestry-par1_HWE.tsv ancestry-par2_HWE.tsv*

2) Summarize average ancestry in windows based on the intervals in a bed file:

`Rscript average_ancestry_bybedintervals.R infile bins.bed chr_name`

*The script will write an output file called average_infile_ancestry_bins.bed*

- 2a) If you want to summarize aver ancestry in windows for all chromosomes in the bed file:

`Rscript average_ancestry_bybedintervals_wholeGenome.R infile bins.bed`

*The script will write an output file called average_infile_ancestry_bins_WG*

## Thin Ancestry tsv files by genetic distance
1) After running AncestryHMM, extract the marker names:

`head -n 1 ancestry-probs-par1_allchrs.tsv > my_markers`

transpose these markers and reformat:

```
perl transpose_nameout.pl my_markers
perl -pi -e 's/:/\t/g' my_markers_transposed
```

2) Thin the markers from the X. birchmanni 10x assembly to the desired genetic distance. It's a good idea to also impose a physical distance threshold because of
inaccuracies in the recombination map at a fine scale:

`Rscript calculate_genetic_distance_and_thin_adjacent_markers.R my_markers_transposed cMthresh physical_dist_thresh`

*This script will write an output file called infile_cMdistances_thinned_distancethresh*

**Note that there may be greater distance between markers if the region is marker poor, this will be recorded in the output file**

3) Select these markers from your tsv file

```
cut -f 1,2 infile_cMdistances_thinned_distancethresh | perl -p -e 's/\t/:/g' > my_thinned_markers_reformat

perl select_passing_markers_multi_geno_files_v2.pl my_thinned_markers_reformat ancestry-probs-par1_allchrs.tsv
path_to:transpose_nameout.pl outfile_name_par1

perl select_passing_markers_multi_geno_files_v2.pl my_thinned_markers_reformat ancestry-probs-par2_allchrs.tsv
path_to:transpose_nameout.pl outfile_name_par2
```


## Ancestry tsv to identifying ancestry transitions
1) After running AncestryHMM, convert your ancestry tsv files to hard calls using parsetsv_to_genotypes_Dec2017_v2.pl

`perl parsetsv_to_genotypes_Dec2017_v2.pl ancestry-probs-par1_file ancestry-probs-par2_file
genotypes_outputfile_name`

2) Run intervals script:

`Rscript identify_intervals_10x_genomes.R genotypes_outputfile_name path_to:transpose_nameout.pl`

*Writes output to: genotypes_outputfile_name_focal_chr, with one file for each chromosome*

## Ancestry tsv files to admixture LD decay
We often want to quantify the change in admixture LD over genetic or physical distance. This is important for dating hybrid populations, understanding what regions are
independent from each other in hybrid populations, selection, and much more.

To quantify admixture LD, first convert genotypes files to plink input files. It's a good idea to thin the data first either by uniqueness or by physical or genetic distance.
1) Convert genotypes files to plink input files

`perl /home/groups/schumer/shared_bin/Lab_shared_scripts/convert_msg_genotypes_to_plink.pl
genotypes_ACUA_2018.txt.identical_filter.txt`

- 1a) You can filter your genotypes files to generate a less redundant input dataset. Either the filter_identical_columns_threshold.pl or thin_genotypes_by_physical_distance.pl
will work. I recommend thinning by physical distance for this application:

`perl thin_genotypes_by_physical_distance.pl genotypes_file distance_thresh_bp path_to:transpose_nameout.pl`

*Writes an output file named genotypes_file_thinned_physical_dist.txt*

`perl filter_identical_columns_threshold.pl genotypes_file num_markers_differentiation path_to:transpose_nameout.pl`

*Writes an output file with genotypes_file.identicalfilter.txt'*

- 2a) You can filter your data with plink before quantifying LD. This command excludes AIMs that are missing in 50% or more of individuals and that are >98% one ancestry
or another

`/home/groups/schumer/shared_bin/plink --file genotypes_ACUA_2018.txt.identical_filter.txt --recode --tab --geno 0.5
--maf 0.02 --allow-extra-chr --out genotypes_ACUA_2018.txt.identical_filter.txt.trim`

- 2b) Run plink to quantify admixture LD. The downstream scripts assume that you are outputting R2 and D so modify accordingly if needed

`/home/groups/schumer/shared_bin/plink --file genotypes_ACUA_2018.txt.identical_filter.txt.trim --ld-window-kb 3000
--ld-window 20000 --r2 d --hardy --hwe 0.00001 --out ACUA --ld-window-r2 0 --allow-extra-chr`

*See plink documentation (https://web.archive.org/web/20240414053322/https://www.cog-genomics.org/plink/1.9/) for details on these parameters*

- 3a) Convert plink physical distance into genetic distance based on the X. birchmanni 10x genome (do not use this script for other genome assembly versions):

`Rscript /home/groups/schumer/shared_bin/Lab_shared_scripts/convert_plink_LD_genetic_distance_xbir10x.R ACUA.ld
ScyDAA6-2-HRSCAF-26`

*This will output a file called ACUA.ld_ScyDAA6-2-HRSCAF-26_cMdistances*

- 3b) If you want a genome-wide picture combine files from multiple chromosomes. e.g.:

`cat ACUA.ld_ScyDAA6-*_cMdistances | grep -v totalcM > ACUA.ld_allscaffolds_combined_cMdistances`

Note that to use the downstream scripts you'll need to add back in the original header

4) To quantify admixture decay over genetic distance in bins of N cM:

`Rscript /home/groups/schumer/shared_bin/Lab_shared_scripts/convert_cM_version_plink_admixture_decay.R
ACUA.ld_allscaffolds_combined_cMdistances 0.2`

*Will produce a file called: ACUA.ld_allscaffolds_cMdistances_admixture_ld_decay_cMdist*

5) Use the decay in admixture LD over genetic distance to estimate the time of initial admixture. Note that there are many assumptions that go into this estimate, including
large population sizes, a single pulse of migration, etc, so interpret with caution:

`Rscript /home/groups/schumer/shared_bin/Lab_shared_scripts/estimateage_plink_admixture_decay_v2.R
ACUA.ld_allscaffolds_cMdistances_admixture_ld_decay_cMdist 0.1`

*The last argument here is for excluding distances smaller than mincM - here 0.1 - since the recombination map is not as reliable at small scales*

*This command will produce a pdf plotting the admixture decay named ACUA.ld_allscaffolds_cMdistances_admixture_ld_decay_cMdist_lddecay.pdf and will print the
estimated mixture time and standard deviation to the screen*

*a is the coefficient of the model, b is the estimated admixture time, and c is the LD decay asymptote* 

## Merging multiple ancestry tsv files run separately
- Warning!!! it is only appropriate to do this if the runs were performed with the exact same parameters in the cfg file
- It is always better to run all individuals through the HMM together if possible, but sometimes this is not practically for very large numbers of individuals

The example shown here merges three ancestry tsv files, but can be done for any number:

1) generate lists markers shared between all files

```
head -n 1 ancestry-probs-par1_allchrs_batch1.tsv > marker_list_batch1
head -n 1 ancestry-probs-par1_allchrs_batch2.tsv > marker_list_batch2
head -n 1 ancestry-probs-par1_allchrs_batch3.tsv > marker_list_batch3
```

transpose files so there is one marker per line:
```
perl /home/groups/schumer/shared_bin/Lab_shared_scripts/transpose_nameout.pl marker_list_batch1
perl /home/groups/schumer/shared_bin/Lab_shared_scripts/transpose_nameout.pl marker_list_batch2
perl /home/groups/schumer/shared_bin/Lab_shared_scripts/transpose_nameout.pl marker_list_batch3
```

*this will output files named e.g. marker_list_batch1_transposed*

2) make a list of marker files:

`ls marker_list_batch*_transposed > all_markers_list`

3) identify markers that are covered in all tsv files:

`perl /home/groups/schumer/shared_bin/Lab_shared_scripts/generate_list_of_passing_markers_multi_geno_files.pl
all_markers_list`

*this will write an output file with passing_markers_FILENAME, i.e. for this example: passing_markers_all_markers_list*

4) select these markers from your tsv files:

```
perl /home/groups/schumer/shared_bin/Lab_shared_scripts/select_passing_markers_list_genotype_files_v2.pl
passing_markers_all_markers_list ancestry-probs-par1_allchrs_batch1.tsv
/home/groups/schumer/shared_bin/Lab_shared_scripts/ ancestry-probs-par1_allchrs_batch1_selectshared.tsv

perl /home/groups/schumer/shared_bin/Lab_shared_scripts/select_passing_markers_list_genotype_files_v2.pl
passing_markers_all_markers_list ancestry-probs-par2_allchrs_batch1.tsv
/home/groups/schumer/shared_bin/Lab_shared_scripts/ ancestry-probs-par2_allchrs_batch1_selectshared.tsv

perl /home/groups/schumer/shared_bin/Lab_shared_scripts/select_passing_markers_list_genotype_files_v2.pl
passing_markers_all_markers_list ancestry-probs-par1_allchrs_batch2.tsv
/home/groups/schumer/shared_bin/Lab_shared_scripts/ ancestry-probs-par1_allchrs_batch2_selectshared.tsv

perl /home/groups/schumer/shared_bin/Lab_shared_scripts/select_passing_markers_list_genotype_files_v2.pl
passing_markers_all_markers_list ancestry-probs-par2_allchrs_batch2.tsv
/home/groups/schumer/shared_bin/Lab_shared_scripts/ ancestry-probs-par2_allchrs_batch2_selectshared.tsv

perl /home/groups/schumer/shared_bin/Lab_shared_scripts/select_passing_markers_list_genotype_files_v2.pl
passing_markers_all_markers_list ancestry-probs-par1_allchrs_batch3.tsv
/home/groups/schumer/shared_bin/Lab_shared_scripts/ ancestry-probs-par1_allchrs_batch3_selectshared.tsv

perl /home/groups/schumer/shared_bin/Lab_shared_scripts/select_passing_markers_list_genotype_files_v2.pl
passing_markers_all_markers_list ancestry-probs-par2_allchrs_batch3.tsv
/home/groups/schumer/shared_bin/Lab_shared_scripts/ ancestry-probs-par2_allchrs_batch3_selectshared.tsv
```

5) merge these outfiles:
- *Warning!! carefully check all file names when doing this merging, it is easy to make mistakes such as swapping par1 and par2 files which will have serious consequences*

We need to retain the header from one file of each the parent 1 and parent 2 files, so we will leave these two files unmodified:

```
ancestry-probs-par1_allchrs_batch1_selectshared.tsv

ancestry-probs-par2_allchrs_batch1_selectshared.tsv
```

For the other files we need to trim the header:

```
tail -n +2 ancestry-probs-par1_allchrs_batch2_selectshared.tsv > ancestry-probs-
par1_allchrs_batch2_selectshared_noheader.tsv

tail -n +2 ancestry-probs-par1_allchrs_batch2_selectshared.tsv > ancestry-probs-
par2_allchrs_batch2_selectshared_noheader.tsv

tail -n +2 ancestry-probs-par1_allchrs_batch3_selectshared.tsv > ancestry-probs-
par1_allchrs_batch3_selectshared_noheader.tsv

tail -n +2 ancestry-probs-par1_allchrs_batch3_selectshared.tsv > ancestry-probs-
par2_allchrs_batch3_selectshared_noheader.tsv
```

Combine all files for parent 1 and parent 2:

```
cat ancestry-probs-par1_allchrs_batch1_selectshared.tsv ancestry-probs-
par1_allchrs_batch2_selectshared_noheader.tsv ancestry-probs-par1_allchrs_batch3_selectshared_noheader.tsv >
ancestry-probs-par1_allchrs_combinedindividuals_sharedmarkers.tsv

cat ancestry-probs-par2_allchrs_batch1_selectshared.tsv ancestry-probs-
par2_allchrs_batch2_selectshared_noheader.tsv ancestry-probs-par2_allchrs_batch3_selectshared_noheader.tsv >
ancestry-probs-par2_allchrs_combinedindividuals_sharedmarkers.tsv
```

## Identifying minor parent deserts and islands
This workflow will identify regions of high or low minor parent ancestry in the genome for a population of interest.
First convert the ancestry-tsv files to average ancestry by site (XXX is the 4-etter code for the population):

`perl /home/groups/schumer/shared_bin/Lab_shared_scripts/calculate_avg_ancestry_by_site_tsv_files_v5.pl ancestry-
probs-par1.tsv ancestry-probs-par2.tsv 10 0.9 /home/groups/schumer/shared_bin/Lab_shared_scripts/ >
average_ancestry_by_site_XXX.txt`

Then get the ancestry in 0.05cM windows (the 0.05cM windowed Xbir genome is here
(`/home/groups/schumer/shared_bin/shared_resources/xbir_genome_windowed/`):

`Rscript /home/groups/schumer/shared_bin/Lab_shared_scripts/average_ancestry_bybedintervals_wholeGenome.R
average_ancestry_by_site_XXX.txt ancestry_xbir_genome_0.05cM_windows.bed`

Then you can run:

`Rscript identifyDesertsIslands_bySite_2.5per_filter0.5cM_10SNPs_mergeShort_v2.R XXX`

This script will look for the two files produced previously, `average_ancestry_by_site_XXX` and
`average_average_ancestry_by_site_XXX.txt_ancestry_xbir_genome_0.05cM_windows.bed_WG`.

It will error out if these files are not in the directory you are running this script.
It produces raw desert and island files from just the site-wise data and final versions of these desert and island where regions where the 0.05cM window that contains the
regions midpoint is also an ancestry outlier and where regions with <10 SNPs are removed (`cMpass`). Then regions < 50,000bps apart are merged
(`cMpass_shortMerged`).

## Mapping and variant calling with GATK
One of the most common workflows in the lab is moving from fastq.gz files to vcf files. This involves a large number of steps and can differ between programs used for
mapping and variant calling. Below we outline the most commonly used steps in our lab and two shell scripts that can be used to run them in an automated way:

0) before you start, index your reference genome bwa index, gatk, and samtools. This only needs to be done the first time you run this analysis in a particular directory.

```
#load modules
module load biology
module load bwa
module load samtools
#index
bwa index ref_genome.fa
java -jar /home/groups/schumer/shared_bin/CreateSequenceDictionary.jar R=ref_genome.fa O=ref_genome.dict
samtools faidx ref_genome.fa
```

1) mapping reads to the reference genome:

```
module load biology
module load bwa
bwa mem -t 3 -M -R '@RG\tID:id\tSM:quail_libraryprep\tPL:illumina\tLB:lib1\tPU:illuminaHiSeq' ref_genome.fa
read_1.fq.gz read_2.fq.gz > file.sam
```

2) process and sort sam file

```
java -jar /home/groups/schumer/shared_bin/SortSam.jar INPUT=file.sam OUTPUT=file.sorted.bam SORT_ORDER=coordinate
java -jar /home/groups/schumer/shared_bin/BuildBamIndex.jar INPUT=file.sorted.bam
```

3) mark read duplicates in the bam file

```
java -jar /home/groups/schumer/shared_bin/MarkDuplicates.jar INPUT=file.sorted.bam OUTPUT=file.sorted.dedup.bam METRICS_FILE=file.sorted.metrics
java -jar /home/groups/schumer/shared_bin/BuildBamIndex.jar INPUT=file.sorted.dedup.bam
```

4) re-align reads around indels

```
java -jar /home/groups/schumer/shared_bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ref_genome.fa -I
file.sorted.dedup.bam -o file.sorted.dedup.bam.list
java -jar /home/groups/schumer/shared_bin/GenomeAnalysisTK.jar -T IndelRealigner -R ref_genome.fa -I
file.sorted.dedup.bam -targetIntervals file.sorted.dedup.bam.list -o file.sorted.dedup.realigned.bam
```

5) perform variant calling

```
java -jar /home/groups/schumer/shared_bin/GenomeAnalysisTK.jar -T HaplotypeCaller -R ref_genome.fa -I
file.sorted.dedup.realigned.bam --genotyping_mode DISCOVERY -L chromosome_targets.list -stand_emit_conf 10 -
stand_call_conf 30 -ERC GVCF -o file.sorted.dedup.realigned.bam.rawvariants.g.vcf
java -jar /home/groups/schumer/shared_bin/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ref_genome.fa --variant
file.sorted.dedup.realigned.bam.rawvariants.g.vcf --sample_ploidy 2 --max_alternate_alleles 4 --
includeNonVariantSites --standard_min_confidence_threshold_for_calling 30 -o file.sorted.dedup.realigned.bam.g.vcf
```

We have two shell scripts that will run these steps for lists of reads.

The first runs the mapping, sorting, de-duplicates, and realigns indels:


`perl /home/groups/schumer/shared_bin/Lab_shared_scripts/submit_bwa_jobs_list.pl list_of_fastq_files PE_or_SE
slurm_submission_file path_to_picard_tools_and_gatk_jars genome_to_map_to file_name_tag`

Example: `perl /home/groups/schumer/shared_bin/Lab_shared_scripts/submit_bwa_jobs_list.pl
combined_read_list_coac_ld_map PE /home/groups/schumer/shared_bin/Lab_shared_scripts/submit_slurm_bwa.sh
/home/groups/schumer/shared_bin xiphophorus_birchmanni_10x_12Sep2018_yDAA6.fasta xbir-10x`

The reads file should be in the format: `id1_read1.fq.gz\tid1_read2.fq.gz\n`


The second shell script, which can be run after the first finishes, runs the actual variant calling steps (see Submitting slurm jobs for tips). It requires a list of sorted, realigned
bam files and a chromosome targets list, among other parameters. This can include all chromosomes but splitting into several batches is helpful for running in parallel.

`perl /home/groups/schumer/shared_bin/Lab_shared_scripts/submit_gatk-indel-hc_jobs_list.pl
list_of_sorted_realigned_bam_files slurm_submission_file_to_use absolute_path_to_GATK_jars genome_assembly_used
vcf_target_chromosome_list output_tag`

Example: `perl /home/groups/schumer/shared_bin/Lab_shared_scripts/submit_gatk-indel-hc_jobs_list.pl
xmac_mapped_bam_list /home/groups/schumer/shared_bin/Lab_shared_scripts/submit_gatk-indel-hc.sh
/home/groups/schumer/shared_bin/ xma_washu_4.4.2-jhp_0.1_combined-unplaced-mito.fa xmac_targets.group1.list group1`

## Case/Control GWAS from low coverage data
Genome-wide association studies (GWAS) are a common way to identify the genetic basis of polymorphisms that exist within a single (unstructured) population. For simple
(ie binary) traits that can be classified as a "case" or as a "control" (eg spotted caudal vs wildtype), it is possible to identify SNPs associated with one trait or the other from
low coverage data, given a large enough sample size.

Because there are multiple steps involved, the lab has a number of shell scripts that simplify the process:

1) Getting set up
1a) make a reads list for input in the next step. The format should be the entire path to the file (get this with the "pwd" command) and the file name. For PE reads, this will
have two columns, one for read1 and one for read2. This file should be tab delimited:

```
tail -n 5 gwas_higher_coverage_indivs.txt
/path/to/file/COAC_98SC_J_combined_read_1.fastq.gz /path/to/file/COAC_98SC_J_combined_read_2.fastq.gz
/path/to/file/COAC_98wt_M_combined_read_1.fastq.gz /path/to/file/COAC_98wt_M_combined_read_2.fastq.gz
/path/to/file/COAC_99SC_J_combined_read_1.fastq.gz /path/to/file/COAC_99SC_J_combined_read_2.fastq.gz
/path/to/file/COAC_99wt_M_combined_read_1.fastq.gz /path/to/file/COAC_99wt_M_combined_read_2.fastq.gz
/path/to/file/COAC_9SC_J_combined_read_1.fastq.gz /path/to/file/COAC_9SC_J_combined_read_2.fastq.gz
```

generate a case control list which is simply a list of 0s and 1s, in the same order and corresponding to the reads list. For example, for the part of the reads list posted above,
the case control list for a SC vs WT GWAS would be:

```
tail -n 5 case_control_status.txt
1
0
1
0
1
```

- 1b) the relevant scrips are all stored in /home/groups/schumer/shared_bin/Lab_shared_scripts/

```
run_samtools_only_parental_v2.pl
run_mpileup_bcftools_GWAS.pl
print_alleles_depth_freq_chi_per_site_GWAS.pl
merge_files_using_two_columns_sharing_values.pl #optional: copy these to your working directory
```

- 1c) Load appropriate modules. This is a list of the one's you might need:

```
module load biology
module load bwa
module load samtools
module load bcftools
module load java
```

- 1d) If you haven't already done so, index your reference genome with bwa index, gatk, and samtools. This only needs to be done the first time you run this analysis in a
particular directory.

```
bwa index ref_genome.fa
java -jar /home/groups/schumer/shared_bin/CreateSequenceDictionary.jar R=ref_genome.fa O=ref_genome.dict
samtools faidx ref_genome.fa
```

2) generate sam file of reads mapped to reference genome (most data in lab is PE [paired end]). reformat sam file to bam file, sort file, remove duplicates, and realigns around
indels:

    usage: `perl run_samtools_only_parental_v2.pl read_list genome SE_or_PE`

    example: `perl run_samtools_only_parental_v2.pl gwas_higher_coverage_indivs.txt ref_genome.fasta PE`

3) run samtools mpileup by scaffold (increases speed of run by parallelizing by scaffold). 
    
    generates a .vcf file
    usage: `perl run_mpileup_bcftools_GWAS.pl read_list case_control_status genome1 samtools_legacy_path focal_scaff`

    example: `perl run_mpileup_bcftools_GWAS.pl gwas_higher_coverage_indivs.txt case_control_status.txt ref_genome.fasta
/home/groups/schumer/shared_bin/ ScyDAA6-1508-HRSCAF-1794`

4) reformats .vcf output for each chromosome into a more easily readable summary output
usage: `perl print_alleles_depth_freq_chi_per_site_GWAS.pl mpileup_file_scaff.vcf low_dp_threshold > mpileup_file_scaff.vcf.summary`

    example: `perl print_alleles_depth_freq_chi_per_site_GWAS_v3.pl gwas_higher_coverage_indivs.txtindivs.allindiv.ScyDAA6-
1508-HRSCAF-1794.vcf 30 > gwas_higher_coverage_indivs.txtindivs.allindiv.ScyDAA6-1508-HRSCAF-1794.vcf.summary`

- 4a) check if files are the right size. sometimes step 4a doesn't work: 

    `ls -sh *summary`

5) combine the summary file for each scaffold into single file containing all scaffolds. precise steps will depend based on analysis.
`gwas_higher_coverage_indivs.txtindivs.allindiv.allchroms.vcf.summary`

    As an example, here's what I did to combine all the individual summary files into one summary file + remove duplicate headers (-Gabe):
`cat *.summary > [PREFIX].allindiv.allchroms.vcf.summary`

    `awk '/group\tpos\tref_allele\talt_allele\tdepth\tMAF\tAF-group1\tAF-group2\tLR1\tLR2\tpchi/&&c++>0 {next} 1'
[PREFIX].allindiv.allchroms.vcf.summary > [PREFIX].allindiv.allchroms.vcf.summary.dedupheaders`

6) *optional step*: only keep known, high confidence SNPs
    
    usage: `perl merge_files_using_two_columns_sharing_values.pl file1 0 1 file2 0 1`

    example: `perl merge_files_using_two_columns_sharing_values.pl
COAC_population_confirmed_snps_Xbirchmanni_10X_genome_flipped_sorted.bed 0 1
gwas_higher_coverage_indivs.txtindivs.allindiv.allchroms.vcf.summary 0 1`

7) convert scaffold names to chromosome numbers. there are currently perl scripts for the xbir_10x, xbir_pacbio, and xvar_scaff genomes. output is a file named
*_formanhattan with numbers in column 1.


    usage: `perl convert_birchmanni10x_mapping_output_manhattan_plot_input.pl summary_file`

    example: `perl convert_birchmanni10x_mapping_output_manhattan_plot_input.pl gwas_higher_coverage_indivs.txtindivs.allindiv.allchroms.vcf.summary`

You are now ready to plot the output! The significance line will vary and can be determined using permutation, but for the X. birchmanni COAC population, it's in the
ballpark of 6.5.

## Pseudohaploid calls for GWAS or population structure analysis from low coverage data
It is useful for many analyses to generate 'pseudohaploid' calls, where a single read is selected at random at each covered site. To begin this process, generate bam files for
the samples of interest, as described above.

Once you have bams you can use bcftools to generate vcf files that retain information at each allele. I recommend using the targets function if your goal is to evaluate the
same SNPs in multiple samples so you can more easily combine those files.

(update GAP 20231108: I wrote a script (https://web.archive.org/web/20240414053322/https://github.com/Schumerlab/Lab_shared_scripts/blob/master/mpileup_out_to_plink_format.py) that takes the output of the "Case/Control GWAS from low coverage data" as an input and outputs the same data in the format below)

0) generate targets file. The targets file can be generated by making a list that looks like this:

```
chr-01 20478 G,T
chr-01 20542 G,C
chr-01 20572 T,C
chr-01 21199 C,T
chr-01 21402 A,G
```

and then bgzipping and indexing it, for example by:

```
module load htslib
bgzip -c all_chrs_targets.tsv > all_chrs_targets.tsv.gz && tabix -s1 -b2 -e2 all_chrs_targets.tsv.gz
```

1) Next you can run bcftools on only those target SNPs. For example:

    `bcftools mpileup -f xbir-COAC-16-VIII-22-M_v2023.fa -T all_chrs_targets.tsv.gz my.unique.realigned.bam | bcftools
call -mO z -i -T all_chrs_targets.tsv.gz -C alleles -i -o my.unique.realigned.vcf.gz`

2) Now that we have a vcf we can generate pseudohaploid calls. First unzip the file:

`gunzip my.unique.realigned.vcf.gz`

Next run the pseudo haploid script to randomly sample alleles at covered sites (NAs will be given to sites that are not covered):

`perl pseudo_haploid_calls_bcftools.pl my.unique.realigned.vcf`

*This will generate an output file called my.unique.realigned.vcf.pseudohap.txt*

3) If you are running multiple individuals, run all of these individuals and combine them (making sure the line numbers of the pseudohap.txt output files are the same -
sometimes INDELs can cause issues).

This data frame can be used for PCA (be aware of effects of asymmetry in missing data) or to do a pseudohaploid GWAS.

## Mapping Binary & Continuous Traits with PLINK
This workflow supposes that you have variant calls as generated in "Case/Control GWAS from low coverage data" that were converted to pseudohaploid calls using the
workflow "Pseudohaploid calls for GWAS or population structure analysis from low coverage data". I'm going to present an updated method of doing this here, just to reflect
exactly what I did when I was running this worflow. Also I think you could run most, if not all, of these steps on a sdev nodes (espsecially the PLINK steps, those go very quickly). I think the longest step here is the variant calling in step 1 which took a couple hours. -GAP

1) Instead of running: `perl pseudo_haploid_calls_bcftools.pl my.unique.realigned.vcf`
Run this instead:

```
module load biology samtools bcftools
perl /home/groups/schumer/shared_bin/Lab_shared_scripts/generate_pseudohaploid_calls_from_bams.pl bam_list.txt
fasta_file.fa chr_targets.tsv.gz
```

Where `bam_list.txt` is list of .bam files (and paths if necessary), `fasta_file.fa` is the reference genome (which should have relevant index files in the same directory), and
`chr_targets.tsv.gz` is the compressed set of reformatted target SNPs.

This will generate *.pseudohaploid.txt files for each .bam file

2) check that the *.pseudohaploid.txt files have the same line counts as chr_targets.tsv.gz (each *.pseudohaploid.txt file should have the same number of lines)

3) If all that looks good, concatenate all the *pseudohap.txt files:

`ls *pseudohap.txt > my_pseudohap_list`

Then run the following script:

`perl /home/groups/schumer/shared_bin/Lab_shared_scripts/combine_pseudohaploid_files_pop_structure_analysis.pl
my_pseudohap_list`

This should produce a file with a name like xmul_pseudohap_list_data_frame

4) Convert the combined pseudohap file to a .ped and .map files for PLINK. You will need your phenotype data for this step and it must be matched exactly to your
combined pseudohap file in number and order. For binary traits, use 1/2 for case/control. For continuous traits, just include the value of whatever you're measuring. The
workflow for both cases is the same until the actual GWAS step.

For this step, you will need your phenotype data in the following format:

```
TABQ-18-XI-22_Xmul_25_M.R1.fastq.pseudohap.txt_trim 1
TABQ-18-XI-22_Xmul_26_M.R1.fastq.pseudohap.txt_trim 1
TABQ-18-XI-22_Xmul_27_M.R1.fastq.pseudohap.txt_trim 1
TABQ-18-XI-22_Xmul_28_M.R1.fastq.pseudohap.txt_trim 2
TABQ-18-XI-22_Xmul_50_M.R1.fastq.pseudohap.txt_trim 2
TABQ-19-XI-22_Xmul_01_M.R1.fastq.pseudohap.txt_trim 1
TABQ-19-XI-22_Xmul_02_M.R1.fastq.pseudohap.txt_trim 1
TABQ-19-XI-22_Xmul_03_M.R1.fastq.pseudohap.txt_trim 1
```

Run the following script to generate the .map and .ped files:

`perl /home/groups/schumer/shared_bin/Lab_shared_scripts/convert_pseudohap_to_plink_pheno.pl
my_pseudohap_list_data_frame.txt my_phenotypes.txt`

At this point, you can run regular association tests with GWAS using PLINK (--assoc). However, to account for population structure, it is better to run a genomic PCA to
identify potential principal components reflecting population structure that you can include in your GWAS as covariates.

5) Run a genomic PCA with PLINK to look for population structure. PLINK will use "my_pseudohap_list_data_frame.txt" as a base name for the other input and output files
you need for PCA or GWAS or whatever other PLINK analysis you want to run. As long as you do everything in the same directory at this point, you'll be in good shape.

```
module load plink
plink --noweb --file my_pseudohap_list_data_frame.txt --pca --mind 0.75 --allow-extra-chr --allow-no-sex --out
my_pseudohap_list_data_frame.txt
```

This will produce 2 files: my_pseudohap_list_data_frame.txt.eigenval and my_pseudohap_list_data_frame.txt.eigenvec. The former contains the eigenvalues for the PCA,
while the later contains the actual PC values for each individual in your sample. At this point you should download both of these and plot them in R. Important: `--mind 0.75`
will filter out samples for missingness. Feel free to play with this parameter, but 0.75 has worked for me + seems to work for other swordtail species according to Molly.

Do the following separately in R:

For the eigenval file, you can calculate percent variance explained for each PC as PVE = eigenval/sum of all eigenvals * 100. Plotting the distribution of PVE for each PC
can tell you PCs are driving the variance.

For the eigenvec file, treat this like a regular dataframe and plot PC1 on the x-axis and PC2 on the y-axis to get your PCA plot.

Once you have identified the PCs contributing the most to the overall variation, you can move on to running a GWAS with covariates.

6) Extract the first n PCs (in my case it was 3 but check your PVE distribution):

`cat my_pseudohap_list_data_frame_.txt.eigenvec | perl -p -e 's/ +/\t/g' | cut -f 1-n | perl -p -e 's/\t/ /g' > my_pseudohap_list_data_frame_.txt.eigenvec.firstn`

At this point, you should make a "keeplist" file of the number of individuals that **didn't** get filtered out during the PCA as well as their "fam IDs" that PLINK assigns them
(PLINK output will tell you how many got filtered out). The reason for this is that if you try to use PCs as covariates with the full my_pseudohap_list_data_frame.txt, PLINK
will throw an error because the number of samples in the keeplist and the dataframe are different. You can do that with the following command:


`awk '{print $1"\t"$2"}' my_pseudohap_list_data_frame.txt.eigenvec.firstn > my_pca_keeplist.txt`


7) Run GWAS with PLINK. The following command will run a regular association test without covariates (--assoc) as well as a logistic test with covariates (--logistic). **If
you are trying to map a continuous trait, swap --logistic for --linear.**

`plink --file my_pseudohap_list_data_frame.txt --assoc --allow-extra-chr --allow-no-sex --logistic hide-covar --keep
my_pca_keeplist.txt --covar my_pseudohap_list_data_frame.txt.eigenvec.firstn --out my_pseudohap_GWAS_PC1-
n_covariates`

This will result in 2 files, my_pseudohap_GWAS_PC1-n_covariates.assoc and my_pseudohap_GWAS_PC1-n_covariates.assoc.logistic (or .qassoc if you did --linear). These
can be downloaded and plotted in R at this point in Manhattan plots. (Pro tip: read these files in with read.csv("filename", sep="")).

## g.vcf files to Admixtools input
Admixtools is a program developed by the Reich lab that can perform a large number of analyses [[1] (https://web.archive.org/web/20240414053322/https://github.com/DRe
ichLab/AdmixTools)], including calculating D-statistics for four populations.
The following is a workflow for converting from .g.vcf files generated by GATK to input for Admixtools.

1) for each .g.vcf file, generate an insnp file with masked and variant basepairs. See documentation of insnp_v9_gatk3.4_gvcf.py for hard-call filter options.

`python /home/groups/schumer/shared_bin/insnp_v9_gatk3.4_gvcf.py file.g.vcf file.g.vcf.insnp 20 5 40 10 10 4 -12.5
-8.0 5`

2) make a list of all insnp files you would like summarized for analysis using Admixtools, one on each line

3) using this file, run the following script:

`perl /home/groups/schumer/shared_bin/Lab_shared_scripts/join_insnp_files_to_admixtools_v5.pl insnp_list
outfile_name_tag 1 5 500 /home/groups/schumer/shared_bin/Lab_shared_scripts`


where insnp_list is the list generated in step 2, 1 indicates that the reference individual should be printed [0 to omit], 5 indicates that a site will be excluded if more than 5
individuals have missing/masked basepairs, 500 indicates the physical distance over which to thin adjacent markers, and /home/groups/schumer/shared_bin/ is the path to the
dependency script: overlap_list_retain_unmatched-for-join-insnp.pl

This will generate three output files that can be used as input for Admixtools:

`outfile_name_tag.geno outfile_name_tag.snp outfile_name_tag.ind`

4) edit outfile_name_tag.ind to contain the species names of interest in the third column

5) Admixtools complains if the individual names are too long and if the chromosome names are not chr1, chr2, etc. This might need modification depending on the analysis
being run. A quick fix is a perl regex find and replace:

`perl -pi -e 's/group1/chr1/g' outfile_name_tag.snp`

### Batch script to generate insnp files
If you have a lot of g.vcf files to process, it's easier to use a batch script we have in the shared folder to submit all of the inns jobs.

1) make a list of the gvcf files you want to process, if they aren't in the same folder you are working in make sure to include the path

2) run the shell script

usage is: `perl /home/groups/schumer/shared_bin/Lab_shared_scripts/submit_insnp_jobs_list.pl vcf_list
path_to_dependency_scripts submission_header.sh`

example: `perl /home/groups/schumer/shared_bin/Lab_shared_scripts/submit_insnp_jobs_list.pl COAC_vcf_list
/home/groups/schumer/shared_bin/Lab_shared_scripts
/home/groups/schumer/shared_bin/Lab_shared_scripts/submit_slurm_insnp.sh`

## g.vcf files to ped input for use with PLINK
This follows steps 1-5 of the previous section g.vcf files to Admixtools input, then uses the convertf package from Admixtools to convert to ped format:

6) make a parameter file for `convertf`

`cp /home/groups/schumer/shared_bin/AdmixTools/convertf/par.EIGENSTRAT.PED parLDmapbir_eigen_to_ped`

Use a text editor to edit the input parameters to match your file names (which are in eigenstrat format):

```
genotypename: example.eigenstratgeno
snpname: example.snp
indivname: example.ind
and specify the names of the output file for plink:
genotypeoutname: example.ped
snpoutname: example.pedsnp
indivoutname: example.pedind
```

7) Load the modules required for Admixtools and run the convertf program

```
module load openblas
module load gsl
/home/groups/schumer/shared_bin/AdmixTools/bin/convertf -p parLDmapbir_eigen_to_ped
```

## D-statistic analysis with Admixtools
Once you have .geno, .snp, and .ind files, you can run D-statistic analysis with Admixtools

1) Generate a parameter file for Admixtools by copying the example file:

`cp /home/groups/schumer/shared_bin/AdmixTools/examples/parqpDstat parameter_file_forDstat`

2) Edit the parameter file lines:

```
DIR: ../data/ #path to the directory where your files are, can be ./
SSS: allmap #if your files have a common prefix you can put it here
indivname: DIR/SSS.ind #the path and name of your .ind file, this string corresponds to a file ../data/allmap.ind
snpname: DIR/SSS.snp #the path and name of your snp file
genotypename: DIR/SSS.geno #the path and name of your geno file
poplistname: list_qpDstat #the list of populations you are looking at, one on each line
```

3) Load the required modules and run Admixtools:

```
module load openblas
module load gsl
/home/groups/schumer/shared_bin/AdmixTools/bin/qpDstat -p parameter_file_forDstat
```

## F4 ratio tests from fasta files
We have a shell script to convert a list of fasta files to input compatible with admixtools. These scripts and example input files can be found in the folder:

`F4_ratio_window_analysis_fasta_alignments_v2`

This folder contains all the scripts needed to run F4 ratios in windows based on aligned fasta. The shell script to run these steps is:

`run_eigenstrat_to_F4_fasta_alignment_v4.pl`

Usage of this script is:
`perl run_eigenstrat_to_F4_fasta_alignment_v4.pl aligned_fastas.fa individual_file_for_convertf snp_window_size`

See test.indiv for an example of the individual_file_for_convertf and test.fa for an example of the aligned_fastas.fa file

## F4 Ratio test with Dsuite

### Dtrios
Dtrios will take a joint .vcf file created from a list of .bam files for different species (including a .bam file for your outgroup) and calculate F4 statistics for each possible trio
in the .bam list.
1) Prepare .bam list
Fairly straightforward, just create a file with a list of the filenames of your .bam files

example: `ls *.bam > my_bamlist.txt`

2a) Perform joint variant calling with your .bam list to generate a .vcf file using bcftools. If you need a .bcf file for some reason, you can always generate that and convert it
to .vcf for Dsuite. You can also just directly go to .vcf by replacing -Ob to -Ov in the second half of the command (and also changing the outfile suffix to .vcf). Just to
demonstrate how to convert though, this is how you go from .bcf to .vcf.

```
module load biology
module load bcftools
module load vcftools
module load gcc
```

Get .bcf

`bcftools mpileup -Ou -f [reference_genome] -b [my_bamlist.txt] | bcftools call -mv -Ob -o [my_calls.bcf]`

Convert to .vcf

`bcftools convert [my_calls.bcf] -Ov -o [my_calls.vcf]`

2b) Filter .vcf file

`vcftools --maf 0.05 --max-maf 0.9 --vcf [my_calls.vcf] --remove-indels --recode --out
[my_calls_noindels_maf_filter]`

3) Prepare SETS.txt
This is another input file that matches filenames to species labels, with two tab-delimited columns containing the file prefix and the species. This needs to be done for every
species in the list of .bam files. IMPORTANT: for the outgroup, the file prefix is the same as the rest of the species but the species label MUST be "Outgroup" (and make
sure it is capitalized).

`[file_prefix]\t[species_label]`

Here's what it looks like for these 3 example .bam files (X. continens, X. multilineatus, outgroup = X. variatus):

```
Xcontinens_2_Gabe.R1.fastq.gz_Xmac.sorted.dedup.realigned.bam
Xmul-03-III-22_CourterM-MixedMeso_Quail.R1.fastq.gz_Xmac.sorted.dedup.realigned.bam
barcode_trimmed_10X_reads_variatus_R1_001_subset.fastq.gz_Xmac.sorted.dedup.realigned.bam
```

(brackets to make visualizing easier, don't include in actual SETS.txt file)

```
[Xcontinens_2_Gabe.R1]\t[Xcon]
[Xmul-03-III-22_CourterM-MixedMeso_Quail.R1]\t[Xmul]
[barcode_trimmed_10X_reads_variatus_R1_001_subset]\t[Outgroup]
```

4) Run Dtrios with tree

`/home/groups/schumer/shared_bin/Dsuite/./Build/Dsuite Dtrios [my_calls_noindels_maf_filter.vcf] [SETS.txt] -t
[newick_tree] -o [prefix_for_outfiles]`

Tree should look something like this:

`(((((Xbir:1.0,Xmal:1.0):1.0,Xcor:1.0):1.0,((Xcon:1.0,Xmon:1.0):1.0,Xnez:1.0):1.0):1.0,(Xmul:1.0,(Xpyg:1.0,Xnig:1.0):1.0):1.0),Xvar:1.0);`

You'll get 5 output files, the *_BBAA.txt probably being the most relevant (has D and F4 statistics) and the *_tree.txt file for input into Fbranch. For more info you can
check out the github[2] (https://web.archive.org/web/20240414053322/https://github.com/millanek/Dsuite).

### Fbranch
1) Run Fbranch:

`/home/groups/schumer/shared_bin/Dsuite/./Build/Dsuite Fbranch [newick_tree] [Dtrios_output_tree] >
f4_branch_results.txt`

2) Plot the Fbranch results to see pairwise gene flow signals using dtools.py. For --outgroup, use the same species label as that used in the newick tree (not just "Outgroup"
like in Dtrios).

`/home/groups/schumer/shared_bin/Dsuite/utils/dtools.py [f4_branch_results.txt] [newick_tree] --outgroup
[outgroup_species_label]`


## g.vcf files to pseudo-fasta files
To do this, we first want to generate .insnp files with the following script. This takes a list of your g.vcf filenames (that are in the working directory) as input. You do not need
to submit a separate job for this command since it uses a shared submission script (i.e. you can copy/paste/run this directly in the terminal). Note: be sure to use *.g.vcf files,
not *.rawvariants.g.vcf files otherwise the pseudo-fasta won't be updated.

`perl /home/groups/schumer/shared_bin/Lab_shared_scripts/submit_insnp_jobs_list.pl [g_vcf_list]
[path_to_insnp_and_coverage_scripts] [slurm_submission_script]`

example: `perl /home/groups/schumer/shared_bin/Lab_shared_scripts/submit_insnp_jobs_list.pl Xmul_g_vcf_list
/home/groups/schumer/shared_bin/Lab_shared_scripts
/home/groups/schumer/shared_bin/Lab_shared_scripts/submit_slurm_insnp.sh`

Now, we can generate the pseudo-fasta files using seqtk:

`/home/groups/schumer/shared_bin/seqtk [seqtk_command] [reference_genome] [.insnp_file_from_previous_step] >
[PREFIX].pseudoref.fasta`

example: `/home/groups/schumer/shared_bin/seqtk mutfa xiphophorus_birchmanni_10x_12Sep2018_yDAA6_mito.fasta Xmul-03-
III-22_CourterM-MixedMeso_Quail.R1.fastq.gz_xbir-10x.sorted.dedup.realigned.bam.xbir-10.g.vcf.cov-corrected.insnp >
Xmultilineatus_03-III-22_CourterM-MixedMeso_Quail_pseudoref_10x_12Sep2018_yDAA6.fasta`

The resulting file is the pseudo-fasta!

## Mapping RNAseq data and transcript quantification
There are two options for mapping RNAseq data:
1) transcriptome-based mapping
2) mapping to a genome with an exon-aware mapper

### Mapping to a transcriptome
For mapping to a transcriptome I recommend bwa or kallisto. Note that results from both programs will be dependent on the quality of the reference transcriptome.
For bwa:
1) format transcriptome

```
module load biology
module load bwa
module load samtools
bwa index transcriptome.fa
```

2) map to transcriptome

`bwa mem transcriptome.fa read1.fq.gz read2.fq.gz > out.sam`

3) convert to bam

`samtools sort out.sam -o out.bam`

For kallisto:
1) format transcriptome
```
module load kallisto
kallisto index -I my_name_tag my_transcriptome.fa
```

2) run kallisto for each individual. See kallisto documentation for parameter details, here is an example command line:

`kallisto quant -I my_name_tag -o outfile_name --single -l 300 -s 50 --rf-stranded --bootstrap-samples=50
mysample.read1.fq.gz mysample.read2.fq.gz`

3) output of kallisto should be analyzed with the sleuth bioconductor package

### Mapping to a genome

For mapping to a genome with exon-awareness I recommend using STAR:

1) format reference genome

`/home/groups/schumer/shared_bin/STAR/bin/Linux_x86_64_static/STAR --runMode genomeGenerate --genomeDir ./ --
genomeFastaFiles genome.fa --runThreadN 2 --sjdbGTFfile file.gtf --sjdbOverhang 99`

2) Map reads

`/home/groups/schumer/shared_bin/STAR/bin/Linux_x86_64_static/STAR --genomeDir ./ --readFilesIn reads1.fq.gz
reads2.fq.gz --readFilesCommand zcat --outFileNamePrefix name_prefix --runThreadN 10 --outSAMtype BAM
SortedByCoordinate`

## Transcript quantification

There are lots of options for expression quantification and testing for differential expression after you have generated bam files. Some recommended tools include:
- RSEM (https://web.archive.org/web/20240414053322/https://deweylab.github.io/RSEM/)
- sleuth (https://web.archive.org/web/20240414053322/https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html)
*Note: sleuth should be used with kallisto output*

## Statistical analysis with DESeq2
Whether you use an exon-aware mapper and reference genome or a pseudoaligner like kallisto, you will ultimately have some abundance measures that you would like to
statistically analyze. We typically do this analysis using DESeq2 in R. The analysis below is for data aligned with kallisto

First install the following packages in R:

## Allele specific expression from RNAseq data
Our lab pipeline for ASE is packaged into a parallelized pipeline: ncASE_pipeline
### For first time use
### For subsequent use

## Genome assembly from PacBio HiFi data
Long-read sequencing has changed the game! Shortreads (e.g. illumina), which are generally 150bp PE, are cheap and useful for SNP calling. But they have substantial
shortcomings in complex regions of the genome and assemblies based on this technology are fragmented and incomplete.

PacBio HiFi (highly accurate [>99% accuracy] long reads [~10-20kb]) is a leading technology for genome assembly. While such single molecule sequencing is generally
error prone (especially homopolymer INDELs), PacBio HiFi achieves >99% accuracy by reading the same molecule multiple times and taking the consensus. Oxford
Nanopore reads, by contrast, are less accurate (>96-98% with the latest technology) but can be substantially longer (100s of kb to a few kb). Each technology has pros and
cons, and are complementary rather than interchangeable, but we focus here on PacBio data.

Sequencing companies may return the raw data, CCS data, or HiFi data in either bam (unaligned) or fastq format. HiFi reads (after a little processing) is the input for
everything downstream. Raw data is massive and there is no real benefit to backing it up (after confirming the HiFi data are good).

### Resources
We generally use the hifi assembler hifiasm, which is the industry standard for highly accurate, contiguous, phased diploid assemblies. HiFiasm can also integrate, parental
short-read, HiC and ONT data collected from the same individual as well to increase contiguity and phasing performance.

HiFiasm Paper: Cheng, H., Concepcion, G.T., Feng, X., Zhang, H., Li H. (2021) Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. Nat
Methods, 18:170-175. https://doi.org/10.1038/s41592-020-01056-5

HiFiasm docs: https://hifiasm.readthedocs.io/en/latest/index.html

HiFiasm manual on sherlock: `man /home/groups/schumer/shared_bin/hifiasm/hifiasm.1`

### 0a. Getting set up
make a new directory and cd. I like to call this directory my {assembly_prefix}-PacBioHiFi (eg xcor-PTHC-VII-22-M-PacBioHiFi)

```
mkdir ${new_dir}
cd
mkdir raw q20 hifiasm_default nanoplot_filt_out
```

link your hifi reads bam and bam.pbi files in your raw folder

`ln -s /oak/path/to/raw/data/*hifi_reads* raw/`

### 0b. Remove CCS reads not meeting HiFi threshold Q>20 (not necessary for UW)
rationale: Hifi reads are defined as circular consensus (CCS) reads with mean read base accuracy > 99% (i.e. a quality score of >20). If the reads are called "hifi" then they
should meet this threshold but companies are inconsistent about this. If some reads do not meet the Q20 threshold, you should filter the reads (e.g using bamtools) to just
include HiFi reads.

programs: bamtools

resources: ~6-18 hours, 8 threads

input file: ccs.bam

output file: ccs.hifi.bam

```
ml biology bamtools
bamtools filter -in raw/ccs.bam -out raw/ccs.hifi.bam -tag "rq":">=0.99"
```

note: you can set different quality thresholds based on rq value. 0.99 is hifi and = q20 phred score. You can go more stringent but will lose read depth, so assembly might get
worse. q25 is has been another threshold we've explored, but q20 is more standard.

### 1. Convert .bam to .fq.gz (you might have to remove adapters, see hifiadapter filt section)

rationale: assemblers take fastq files, but sequencing companies give back bam files (with a long prefix that doesn't really mean anything). thankfully, it's quite easy to convert bam to fastq.

programs: pbindex bam2fastq

github: https://github.com/PacificBiosciences/pbbioconda

location: local miniconda installation

resources: ~1 hours, 4 threads

input file: bam

output file: fastq.gz

```
pbtools
pbindex raw/hifi.bam
bam2fastq raw/hifi.bam -o q20/hifi
```

2. Assess read quality (recommended)

rationale: It's a good idea to check how the sequencing went -- nanoplot calculates the important stats and visualizes the reads so you can check for irregularities. The main
things to check are your read length and read quality (QV). If you didn't filter the bam by quality score, you should also check that the sequencing company sent you HiFi
reads (CCS reads Q>20) and not CCS reads (CCS reads Q>0). HiFiasm requires HiFi reads but does not check that your reads are HiFi (not CCS).

Note: in the most recent update, fastq.gz read QV is no longer equivalent to bam RQ. We are still trying to get to the bottom of this, but for now are filtering based on rq>=30
in the bam.

programs: nanoplot

github: https://github.com/wdecoster/NanoPlot/

location: install locally on sherlock using pip or conda

NOTE: if program doesn't run try downgrading seaborn, see https://github.com/wdecoster/NanoPlot/issues/222

```
conda install seaborn==0.10.1
conda install NanoPlot
```

resources: ~1.5 hours, 4 threads

input file: q20/hifi.fastq.gz

output file: {PREFIX}.fastq.gz.NanoPlot-report.html {PREFIX}.fastq.gz.NanoStats.txt

```
mkdir nanoplot_filt_out
NanoPlot --fastq hifi.fastq.gz -p PREFIX.hifi.fastq.gz. -t4 -o nanoplot_filt_out
```

Confirm all your reads are Q>20 and check read length (you'll likely have to report this later)

### 3. Genome assembly with HiFiasm

rationale: There are many genome assembly programs for HiFi data, but HiFiasm usually performs the best and is quite fast. The default parameters provide several .gfa files
(like .fasta files, but are graphs with edges connecting sequences), which can easily be converted to .fasta files. The typical HiFiasm consensus fasta for Xiphophorus is
usually a little under chromosome-level, with contigs not spanning all satellite sequences. Playing around with different parameters might increase your contiguity, but
default parameters seem to work best for Xiphophorus. Heng Li suggests increasing -D and --max-kocc increase assembly contiguity at the expense of run time. From
testing, increasing -r and -k also tend to increase contiguity in reptiles, but sometimes these parameters in combination with -D and --max-kocc break contigs. These
parameters don't seem to work in a predictable way in Xiphophorus, so default is probably the best bet.

Note: for ${ASM_PREFIX} it is usually best to name according to sample convention species-pop-(snumber)-date.

programs: hifiasm


github: https://github.com/chhylp123/hifiasm

location: /home/groups/schumer/shared_bin/hifiasm/

resources: ~6 hours, 48 threads

input file: hifi.fastq.gz

output file: ${ASM_PREFIX}*gfa

NOTES: Depending on coverage and genome size, you might have to allocate more memory (for a 100x coverage Xnezzy, --mem=140GB worked). Also set ntasks to
number of threads in your command (in this case, --ntasks=48). The job will get preempted repeatedly unless you set -p schumer. Swordtail genomes tend to assemble in 6-
24 hours, although they might queue for 3-5 days.

```
export PATH=$PATH:/home/groups/schumer/shared_bin/hifiasm
hifiasm -o hifiasm_default/${ASM_PREFIX}.asm -t 48 q20/hifi.fastq.gz
```

note: set your ASM_PREFIX to your sample name and read quality (eg xcor-PTHC-VII-22-M.q20). This will help other people keep track of the files.
```
awk '/^S/{print ">"$2;print $3}' hifiasm_default/${ASM_PREFIX}.asm.bp.p_ctg.gfa >
hifiasm_default/${ASM_PREFIX}.asm.bp.p_ctg.fa #convert .gfa file to fasta
```

### 4. Visualizing output
Visualizing output: the program Bandage is a great way to visualize your phased diploid assembly (*p.utg.gfa). You can BLAST & visualize the results in Bandage as well.
Bandage is a desktop application that can be downloaded here: https://github.com/rrwick/Bandage

### 5. Scaffolding the contigs into chromosome models

rationale: In Xiphophorus, most chromosomes are still in 2-3 pieces (with breaks most common around centromeres or other satellite-rich regions. It's often useful to scaffold
these contigs into chromosome-level models (with the breaks represented with N's). If Hi-C or Omni-C data exist, these can be used to scaffold the genomes using programs
like YaHS. If not, RagTag is an easy way to scaffold to a reference genome. Note: hifiasm rarely (but sometimes) make false joins with contigs containing sequence from
multiple chromosomes. These should be broken before putting contigs into ragtag.

programs: RagTag

github: https://github.com/malonge/RagTag

location: installed locally

```
ml biology samtools python/3.9.0
pip3 install --user RagTag
```

resources: ~10 minutes, 8 threads

input files: Reference genome with only one chr-21 version (the more ancestral one) and contigs from last step
${ASM_PREFIX}.asm.bp.p_ctg.fa.

Note: It's best to use same species confirmed with HiC data. If this does not exist, best practice is to run RagTag 3 times, using 2 closely related species and an outgroup and
then compare the differences between the output files. It is also best to not use species with highly derived chromosomal architecture (e.g. not birchmanni or maculatus).
Hellerii and malinche are two good species to use.

for example, you can generate one for X. malinche like this:

make a file called chr-21.list

this file should have 2 lines, representing each sex chromosome:

```
chr-21-X
chr-21-Y

cat <(seqkit grep -v -f chr-21.list xmal-CHIC-XI-20-M_v2023.2.fa) xmal-CHIC-XI-20-M_v2023.2_chr-21-Y.fa | seqkit
sort -n -o xmal-CHIC-XI-20-M_v2023.2_chr-21-Yonly.fa
```

or one from X. multilineatus like this:

`seqkit grep -p chr-21-Y -v xmul-TABQ-V-22-M_v2023.2.fa -o xmul-TABQ-V-22-M_v2023.2_chr-21-Xonly.fa`

output file: ${ASM_PREFIX}.asm.bp.p_ctg.fa_RagTag.fa

programs: REF=xmal-CHIC-XI-20-M_v2023.2_chr-21-Yonly.fa programs: QUERY=${ASM_PREFIX}.asm.bp.p_ctg.fa

programs: samtools faidx $REF

programs: ragtag.py scaffold $REF $QUERY -t 4

Note: resulting scaffolded assemblies can be compared to each other using minimap2 and visualized with the pafr package. This is a good way to ensure there are no large

scale chromosome rearrangements (i.e. that scaffolding is robust).
### 6. Representing both sex chromosomes in the assembly

rationale: In Xiphophorus, the sex determining region is usually on the distal end of chr-21 and can range from 2-12Mb of unique sequence that contains several important
genes and differs between X and Y. hifiasm usually picks the longer haplotype to represent in the primary contigs, so the RagTag output will either have X or Y (usually Y)
represented.

hifiasm also produces partially phased haplotypes (partially phased because haplotype switches can occur in regions with little heterozygosity or profound structural
differences). Depending on the complexity/architecture of the sex chromosome, the entire region might be fully phased and present in the *p.utg.gfa. If not, it should
certainly be present in the hap1 or hap2 .gfa files (although be warned these might have a switch error and result in chimeric XY chroms).

To represent both haplotypes, map the *p.utg.fa and *hap1*.fa *hap2*.fa back to the ${ASM_PREFIX}.asm.bp.p_ctg.fa_RagTag.fa using minimap2.

programs: minimap2

location: /home/groups/schumer/shared_bin/minimap2

input file: 
`${ASM_PREFIX}.asm.bp.p_ctg.fa_RagTag.fa ${ASM_PREFIX}.asm.hap1.p_ctg.fa ${ASM_PREFIX}.asm.hap2.p_ctg.fa`

output file: 
```
${ASM_PREFIX}.asm.hap1.p_ctg.fa_2_${ASM_PREFIX}.asm.bp.p_ctg.fa_RagTag.fa.paf
${ASM_PREFIX}.asm.hap2.p_ctg.fa_2_${ASM_PREFIX}.asm.bp.p_ctg.fa_RagTag.fa.paf
```

```
/home/groups/schumer/shared_bin/minimap2 ${ASM_PREFIX}.asm.bp.p_ctg.fa_RagTag.fa ${ASM_PREFIX}.asm.hap1.p_ctg.fa >
${ASM_PREFIX}.asm.hap1.p_ctg.fa_2_${ASM_PREFIX}.asm.bp.p_ctg.fa_RagTag.fa.paf
/home/groups/schumer/shared_bin/minimap2 ${ASM_PREFIX}.asm.bp.p_ctg.fa_RagTag.fa ${ASM_PREFIX}.asm.hap2.p_ctg.fa >
${ASM_PREFIX}.asm.hap2.p_ctg.fa_2_${ASM_PREFIX}.asm.bp.p_ctg.fa_RagTag.fa.paf
```

Breakpoints can be identified from the paf file and the reference can be replaced with the haplotype (e.g. using seqkit or bedtools). X or Y can be identified based on
homology to other species. These scaffolds should be named chr-21-Y and chr-21-Y.
For most Xiphophorus genomes, we include the X and Y but mask the recombining parts of the Y in the final assembly. The unmasked X and unmasked Y sequence should
be also saved, since it is often useful to have them.

### 7. Misc.
When you're happy with the quality of the assembly, scaffolds can be renamed to remove the *RagTag suffix, eg using seqkit rename

#### 1a. Removing adapter contamination (non-pacbio data, eg Cantata)

rationale: sequencing companies are supposed to remove adapter sequences before they send us the HiFi reads back in .bam format, but sometimes a small proportion of
reads still contain adapter sequences. This affects assembly quality and sometimes leads NCBI to reject the genomes. The program HiFiAdapterFilt removes these reads
from a .bam or fq.gz file and converts it to a .fq.gz, which is the format you'll need going forward.

Note: HiFiAdapterFilt is a bit finicky and not very smart. Program should be run exclusively on schumer nodes. HiFiAdapterFilt makes lots of intermediate files and if it gets
preempted it will restart and try to filter the intermediate files, which leads to bad results. As written, the following code will run on every *bam and *fastq file in the
directory, so set up your directory with only the hifi.bam file you created in the previous step.

programs: HiFiAdapterFilt

github: https://github.com/sheinasim/HiFiAdapterFilt

location: /home/groups/schumer/shared_bin/HiFiAdapterFilt

resources: ~4 hours, 8-16 threads

input file: bam (or fq.gz)


output file: filt.fq.gz

```
ml biology bamtools ncbi-blast+
export PATH=$PATH:/home/groups/schumer/shared_bin/HiFiAdapterFilt/
export PATH=$PATH:/home/groups/schumer/shared_bin/HiFiAdapterFilt/DB
sh hifiadapterfilt.sh -t 8 -l 30
```

## Genome assembly with supernova - 10x data
Below is a quick guide but check out Patrick Reilly's guide (also in box) File:SupernovaAssembly_v1.pdf
1) Make sure supernova is in your path, e.g.:

`export PATH=/home/groups/schumer/shared_bin/supernova-2.1.1:$PATH`

It's a good idea to embed this in your supernova script

2) set up your slurm script for supernova

supernova has high memory and time requirements

recommended: at least 60 hours, 20 cpus, and high memory (~10G)

3) Example supernova commands:

`supernova run --id Xbirchmanni_10X_ref --fastqs /home/data/Xbirchmanni_10Xchromium_Hudsonalpha_July2018_raw_data/ -
-description run1 --maxreads 280000000
optimal number of reads targets 40-56X coverage
supernova mkoutput --
asmdir=/home/data/Xbirchmanni_10Xchromium_Hudsonalpha_July2018_raw_data/Xbirchmanni_10X_ref/outs/assembly --
outprefix=Xbirchmanni_10X_assembly --style=pseudohap`

## Motif and binding site predictions
1) Generate a position weight matrix for your zinc finger sequence at zf.princeton.edu and save as a text file

2) Convert this output into a file format compatible with MEME

`Rscript /home/groups/schumer/shared_bin/motiflogo.R cat_prdm9_PWM_zfp.txt`

3) Build a background model for the genome of interest

`/home/groups/schumer/shared_bin/meme-5.0.2/src/fasta-get-markov -m 2 -dna
GCF_000181335.3_Felis_catus_9.0_genomic.fna fcatus.bg`

4) Generate a MEME formatted position weight matrix

`/home/groups/schumer/shared_bin/meme-5.0.2/scripts/uniprobe2meme -numseqs 1 -bg fcatus.bg -pseudo 1
cat_prdm9_PWM_zfp.txt.PWM.txt > cat_prdm9_PWM_zfp.meme.txt`

5) Run FIMO on the genome and position weight matrix:

`/home/groups/schumer/shared_bin/meme-5.0.2/src/fimo -bgfile fcatus.bg -oc ./FIMO-results_cat_predictedPRDM9 -thresh
0.00001 --max-stored-scores 300000 -verbosity 5 cat_prdm9_PWM_zfp.meme.txt
GCF_000181335.3_Felis_catus_9.0_genomic.fna`

## Liftovers using cactus alignments
With the new PacBio genomes you may need to make new alignments.
Making alignments with cactus

If you are running liftovers using cactus alignments for the first time, make sure to make a local copy of Bernard's alignments:

`cp /home/groups/schumer/data/cactus_alignments/xip_complete.hal ./`

For example, to liftover between maculatus coordinates and birchmanni coordinates:
1) load required modules

`ml gcc`

2) run liftover:

`/home/groups/schumer/tools/hal/bin/halLiftover xip_complete.hal mac_coords.bed birchmanii bir_liftover_coords.bed`

### Other useful basic commands with halTools
See basic stats including names of the aligned sequences:

`/home/groups/schumer/tools/hal/bin/halStats xip_complete.hal`

Pulls out genomes in fasta format, some contigs might be named with NCBI names:

`/home/groups/schumer/tools/hal/bin/hal2fasta xip_complete.hal`

Generate synteny blocks, useful for looking at genomic rearrangements:

`/home/groups/schumer/tools/hal/bin/halSynteny xip_complete.hal`

## Demographic Inference with ABC

We can estimate hybridization timing, population size, ancestry proportion, and parental migration rates with Approximate Bayesian Computation (ABC) simulating a single
swordtail chromosome with SLiM and fitting population specific parameters of average ancestry, variance in ancestry, and median minor parent tract length.

### Simulations setup
To make csv of input parameters for 1 million simulations:

`Rscript /oak/stanford/groups/schumer/data/ABC_simulation_trees/simParameterInput.R`

Example to run a single SLiM simulation:
```
slim -d SEED=1 -d POPSIZE=2000 -d GEN=120 -d INIT_PROP=0.75 -d PAR1MIG=4.5e-06 -d PAR2MIG=0.002
/home/groups/schumer/shared_bin/Lab_shared_scripts/neutral_admixture_ABC_migration.slim
```

### Post simulations population specific inference
To make demographic inferences you will need population specific parameters: mean hybrid index, coefficient of variation in hybrid index, median minor parent tract
length. It's best to get these just from chromosome 2 since that is what's simulated with the SLiM script. Finally you will also need the number of individuals you used to
get these parameters.

For each simulated tree for each population you can create a population specific simulated tsv matching the number of individuals used to get population parameters.

`python3 /oak/stanford/groups/schumer/data/ABC_simulation_trees/slim_genetree_to_indiv_ancestries.py
dem_abc_sim*seed*.trees *number_of_individuals_in_population* dem_abc_sim_mig*seed*.tsv`

The output name can be change to whatever you want.
These tsvs can be summarized to get the simulation mean and coefficient of variation of ancestry and length of minor tracts. 

`Rscript /oak/stanford/groups/schumer/data/ABC_simulation_trees/abc_summary_stats_chr_win.R *seed* dem_abc_sim_mig*seed*.tsv >> summary_params.txt`

Once this has been done for a large set of simulations you can then determine which simulations fall within 5% of the empirical population values.

## Demographic Inference with PSMC

### plotting PSMC output
after running psmc_plot.pl with -R flag on your primary and concatenated bootstrapped files, you should have 1 *txt file for the primary and 100 *txt files for the bootstraps


```
awk '{print $0, "\t", "species name", "\t", "primary"}' primary_psmc.txt > species_nonbootstrap.txt

# add 2 columns noting your species (or population) and if it was a bootstrap or not

for i in `ls species_round1-100*txt`; do awk '{print $0, "\t" ,"species bootstrap","\t", FILENAME}' $i; done > species_bootstrap.txt

#for all your bootstrap files, add 2 columns noting your species (or population) + if it was a bootstrap and the name of your file. FILENAME is a variable that awk recognizes and can be left as is.

cat species_nonbootstrap.txt species_bootstrap.txt > PSMC_R_plotting_file.txt 

#concatenate your modified psmc output files. now ready for plotting in R
```

in R

```
PSMC_data <- read.csv("PSMC_R_plotting_file.txt", header=FALSE, sep = "\t")

ggplot() +
geom_step(data=subset(PSMC_data, PSMC_data$V6==" species bootstrap "), aes(x=V1, y=V2*10, group=V7), color="#69b9cd", alpha=0.5) +
geom_step(data=subset(PSMC_data, PSMC_data$V6==" species "), aes(x=V1, y=V2*10, group=V7), color="#092d44") +
scale_x_log10( #set limits, breaks, labels however you please
# limits=c(10000, 12000000),
# breaks = c(10000,100000,1000000, 10000000),
# labels = c("10000","100000","1000000","10000000")) +
scale_y_continuous(
# expand=c(0,0), limits = c(0,1250), breaks = c(seq(0,1250,250)) # set limits, breaks, labels however you please) +
theme_classic(base_size = 12) +
theme(axis.title.x = element_text(size = 14),
axis.title.y = element_text(size = 14),
axis.text.y = element_text(size = 12),
axis.text = element_text(size = 12),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
legend.position = "none") +
xlab("years in past") +
ylab(expression(paste("N"[e],"*", 10^{3},sep=" "))) 

# can change scaling so 10^4, 10^6 etc, just also modify the scaling factor in geom_step(). psmc plot automatically outputs 10^4 annotation_logticks(side = "b")
```

## Building a recombination map with LDhelmet

LDhelmet workflow