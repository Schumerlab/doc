# Commonly used scripts

Note: Please add any script that you write and have thoroughly tested to the lab git repository [1] (https://github.com/Schumerlab/Lab_shared_scripts) and document it here

On sherlock our lab shared bin with these scripts can be found here: `/home/groups/schumer/shared_bin/Lab_shared_scripts`

## Contents

1. Usage of commonly used scripts

- 1.1 [Scripts that convert between file formats](#Scripts-that-convert-between-file-formats)
- 1.2 [Scripts that merge files, filter files, or match file contents to lists](#Scripts-that-merge-files-filter-files-or-match-file-contents-to-lists)
- 1.3 [Scripts or commands that extract sequences from files](#Scripts-or-commands-that-extract-sequences-from-files)
- 1.4 [Scripts for managing fastq files](#Scripts-for-managing-fastq-files)
- 1.5 [Scripts for mapping or variant calling](#Scripts-for-mapping-or-variant-calling)
- 1.6 [Scripts that manipulate or summarize fasta files](#Scripts-that-manipulate-or-summarize-fasta-files)
- 1.7 [Scripts for QTL and admixture mapping](#Scripts-for-QTL-and-admixture-mapping)
- 1.8 [Scripts to calculate hybrid index and make ancestry plots post ancestryinfer](#Scripts-to-calculate-hybrid-index-and-make-ancestry-plots-post-ancestryinfer)
- 1.9 [Scripts to average by site or in windows](#Scripts-to-average-by-site-or-in-windows)
- 1.10 [Scripts to summarize features in windows like coding, conserved, and repetitive basepairs](#Scripts-to-summarize-features-in-windows-like-coding-conserved-and-repetitive-basepairs)
- 1.11 [Scripts to identify recombination transitions](#Scripts-to-identify-recombination-transitions)
- 1.12 [Scripts relating to ancestry tracts](#Scripts-relating-to-ancestry-tracts)
- 1.13 [Scripts to summarize genetic distance or admixture LD decay](#Scripts-to-summarize-genetic-distance-or-admixture-LD-decay)
- 1.14 [Scripts that generate or use insnp files](#Scripts-that-generate-or-use-insnp-files)
- 1.15 [Scripts for admix'em simulations](#Scripts-for-admixem-simulations)
- 1.16 [Scripts for SLiM simulations or SLiM output](#Scripts-for-SLiM-simulations-or-SLiM-output)
- 1.17 [Miscellaneous](#Miscellaneous)
- 1.18 [Need to be annotated](#Need-to-be-annotated)

## Usage of commonly used scripts

### Scripts that convert between file formats

Convert fasta to phylip format:

`perl Fasta2Phylip.pl infile.fa > outfile.phy`

Convert phylip to fasta:

`perl Phylip2Fasta.pl infile.phy outfile.fa`

Convert MSG genotypes files to plink file format:

`perl convert_msg_genotypes_to_plink.pl genotypes_file_name`

_Writes two output file: genotypes_file_name.ped; genotypes_file_name.map_

Convert a fasta file to a fastq file:

`perl fasta_to_fastq.pl infile.fasta > outfile.fq`

Convert the ancestry tsv output from AncestryHMM to a hard-calls genotypes files:

`perl parsetsv_to_genotypes_Dec2017_v2.pl ancestry-probs-par1.tsv ancestry-probs-par2.tsv genotypes_output_file_name`

Convert the 3-way ancestry tsv files from AncestryHMM to a hard-calls multi genotypes file:

`perl parse_3way_tsv_to_genotypes_file.pl ancestry-file-stem posterior_prob_threshold > genotypes_output_file_name`

example:

`perl parse_3way_tsv_to_genotypes_file.pl allchrs.tsv 0.9 > genotypes_CAPS_all3ways.txt`

Convert the MSG genotypes file to a format that is compatible with rQTL (still requires phenotypes added, see workflow documentation):

`perl genotypes_to_rqtl_Feb2018_v3.pl genotypes_file_name`

_Writes an output file: genotypes_file_name.rqtl.csv_

Convert a phylip file to an input file for treemix:

`perl phy_to_treemix.pl input.phy population_keys_list > outfile`

_The population keys list contains the number of individuals per population, e.g.: 2\n1\n3\n_

Convert a samtools vcf file to input counts for ASE pipeline:

`perl samtools_vcf_to_ASE_counts.pl infile.vcf`

_Writes an output file: infile.vcf\_ASE\_counts_

Convert a sequence in a clustal alignment to a fasta alignment:

`perl convert_clustal_alignment_to_fasta_alignment.pl clustal_alignment.txt seq_name_to_extract`

### Scripts that merge files, filter files, or match file contents to lists

- Note: there are many useful scripts for combining and filtering files on the FAS scriptome [2] (http://archive.sysbio.harvard.edu/CSB/resources/computational/scriptome/UNIX/)*

Combine read files for two different sequencing runs of the same individual. The two files to be combined *must* be in the same order:

`perl combine_reads_two_lists.pl list_full_path_to_file1_set list_full_path_to_file2_set`

_Writes a new file for each individual, using the name in list 1 appended with \_combined_

Script to identify markers in tsv files violating hardy-weinberg equilibrium and filter them:

`perl determine_HWE_filter_markers.pl ancestry-par1 ancestry-par2 bonferonni_pval_thresh path_to:transpose_nameout.pl`

_Writes a file named infile\_deviating\_markers with deviating markers from HWE, and filtered ancestry tsv files: ancestry-par1\_HWE.tsv ancestry-par2\_HWE.tsv_

Script for filtering genotypes file of redundant columns. The number in the command line corresponds to the number of markers that can differ between adjacent columns for the column to be retained:

`perl filter_identical_columns_threshold.pl genotypes_file num_markers_differentiation path_to:transpose_nameout.pl`

_Writes an output file with genotypes\_file.identicalfilter.txt. Note: columns that differ by NAs are treated as different from the previous column_

Script for grepping a list from a separate file with an option to use grep -w or regular grep

`perl grep_list.pl select_list file_to_grep_from outfile_name grep_w_1_or_0`

Similar script but automatically retains header line:

`perl grep_list_keep_header.pl select_list file_to_grep_from outfile_name grep_w_1_or_0`

Alternately, exclude a list from a file:

`perl grep_v_list.pl list_to_grep_v focal_file outfile_name grepw_0_1`

Takes a list of files, reads all of their entries, and writes out a combined file with duplicate lines removed:

`perl combine_lines_remove_duplicates_list.pl list_of_files outfile_name`

Filter a genotypes file based on an allowed proportion of missing markers. Depends on the script transpose_nameout.pl which is on Sherlock in `/home/groups/schumer/shared_bin/Lab_shared_scripts`

`perl filter_missing_markers_genotypes_file_thresh_v2.pl infile_genotypes path_to:transpose_nameout.pl prop_missing_allowed`

_Writes out an outfile named infile\_genotypes\_filtered\_at\_thresh_

Merge files based shared values in two columns. See FAS scriptome [3] (http://archive.sysbio.harvard.edu/CSB/resources/computational/scriptome/UNIX/) for other merging options based on a single column:

`perl merge_files_using_two_columns_sharing_values.pl file1 column1_file1 column2_file1 file2 column1_file2 column2_file2`

Subset an MSG genotypes file based on a list of markers. Depends on the script transpose_nameout.pl (on Sherlock /home/groups/schumer/shared_bin)

`perl select_passing_markers_multi_geno_files_v2.pl list_of_markers_to_select genotypes_file path_to:transpose_nameout.pl outfile_name`

Subset files based on shared values in the first two columns to exclude values in file1 that are also found in file2:

`perl exclude_shared_values_lists_based_on_first_two_columns.pl file_to_select_from file_to_exclude_from`

_Wites an output file called file\_to\_select\_from\_excluded\_overlap\_file\_to\_exclude\_from_

### Scripts or commands that extract sequences from files

Extract fasta sequences using entries from a bed-formatted list:

`perl extract_bed_seqs.pl list_of_seqs.bed fasta_to_extract outfile_tag_name`

_Writes output to a file named list\_of\_seqs.bed\_tag.fa_

Generate predicted transcript sequence using a gtf file containing the exon coordinates and strand information and the relevant fasta file for the _X. maculatus_ genome :

`perl extract_gtf_seqs_mergetranscript_xmac_genome.pl gene_of_interest.gtf fasta outfile_tag`

_Writes output to a file named gene\_of\_interest.gtf\_tag.fa_

This script is the same as above but prints to standard out for the X. maculatus genome :

`perl extract_gtf_seqs_mergetranscript_printstdout_xmac_genome.pl gene_of_interest.gtf fasta name_tag`

This script extracts the predicted transcript sequence using a gtf file containing the exon coordinates and strand information and the relevant fasta file for the _X. birchmanni_ and _X. malinche_ 10x genomes:

`perl extract_gtf_seqs_mergetranscript_10x_assembly.pl list_of_exons.gtf fasta_file outfile_tag`

This script is the same as above but prints to standard out for the _X. birchmanni_ and _X. malinche_ 10x genomes:

`perl extract_gtf_seqs_mergetranscript_printstdout_10x_assembly.pl list_of_exons.gtf fasta_file outfile_tag`

This script is similar to the above scripts, works on the 10x gtfs, but extracts exons without merging them. Each exon coordinate should be a line in the gtf file, the name tag will be appended to the exon name:

`perl extract_gtf_seqs_print_exons_10x_assembly.pl exons_single_gene.gtf fasta_file name_tag`

fastahack is also an incredibly useful program from extracting sequences from fasta files. It can be found on Sherlock in /home/groups/schumer/shared_bin

`/home/groups/schumer/shared_bin/fastahack genome.fa -r chr:start_coord..stop_coord > outfile.fa`

To extract the predicted cDNA sequences for a list of transcripts you can use the following scripts.

This script extracts predicted cDNA sequences for a list of genes in the 10x genome format, given a gtf file, and genome sequence fasta:

`perl generate_transcript_seqs_list.pl my_gene_list gtf_file genome_sequence outfile_name`

This script extracts predicted cDNA sequences including introns for a list of genes in the 10x genome format, given a gtf file, and genome sequence fasta:

`perl generate_transcript_seqs_list.pl my_gene_list gtf_file genome_sequence outfile_name`

To extract predicted cDNA sequences and run codeml to calculate dN/dS, you can use the following script. This script takes two genomes and gtf files and runs pairwise dN/dS analysis with a list of genes.

`perl extract_bir_mal_10x_seqs_run_codeml.pl list_of_genes_to_analyze gtf_species1.gtf gtf_species2.gtf species1.fa species2.fa codeml.ctl_file`

_Note: It can be used on two different species’ gtfs if you have carefully curated them to include single-copy genes of the same length between the species, but it tends to work better with a pseudo reference for one of the species. Here is an example of how you might use it to calculate dN/dS in X. birchmanni and X. malinche:_

`perl extract_bir_mal_10x_seqs_run_codeml.pl list_of_genes xiphophorus_birchmanni_10x_12Sep2018_yDAA6.gtf xiphophorus_birchmanni_10x_12Sep2018_yDAA6.gtf xiphophorus_birchmanni_10x_12Sep2018_yDAA6.fasta Xmalinche_backbone_10x_12Sep2018_yDAA6.fasta codeml.ctl`

### Scripts for managing fastq files

Scripts for parsing fastq files by i5 and i7 index (see details in: parsing_tn5_data.txt in Dropbox and on the wiki). This script must point to a directory with barcode_splitter.py so if you work with it on a different server edit the path

`divideConquerParser.sh 4 "reads_R1_allanes_combined.fastq.gz reads_R2_allanes_combined.fastq.gz reads_I1_allanes_combined.fastq.gz reads_I2_allanes_combined.fastq.gz" 10 i5_library 4 1`

Count the number of reads in a list of zipped fastq files:

`perl count_reads_fastq_list.pl list_of_fastqgz_files`

_Writes an outfile named list\_of\_fastqgz\_files\_counts_

Count the number of reads in a list of pair-end zipped fastq files:

`perl count_reads_fastq_list.pl list_of_PE_fastqgz_files`

Format of reads list should be:

`read1.fq.gz\tread2.fq.gz\n`

_Writes an outfile named list\_of\_PE\_fastqgz\_files\_counts_

Count the number of reads in the batch of negatives on a Tn5 plate. This script is set up to recognize any file that has NEG in the name, and should be run in the folder of interest.

`perl count_negative_reads_fastq_list_PE.pl library_name`

_Writes an outfile named library\_name\_negative\_counts\_results_

### Scripts for mapping or variant calling

This script takes a list of SE or PE reads and maps them to a particular reference genome and uses samtools to take them through the generation of sorted bam files:

`perl run_map_and_samtools_parental.pl read_list_file genome_to_map_to read_type_SE_or_PE`

_The read list should contain one column if the data is SE and two columns \t delimited if the data is paired end. The output file for each entry in the read list will be: samplename.read1.fq.sorted.unique.bam_

This script can be run after running run_map_and_samtools_parental.pl to run bcftools on the bams generated by the mapping/samtools script:

`perl run_bcftools_v2.pl reads_list genome_assembly`

_The reads\_list file is the same file used as input in run\_map\_and\_samtools\_parental.pl and the command should be executed in the same directory. For each bam file, the script will generate a pileup and vcf file._

This script can be run after executing a case-control GWAS with samtools legacy. This script takes the results of a case-control GWAS using the samtools legacy vcf and reformats them for plotting:

`perl print_alleles_depth_freq_chi_per_site_GWAS.pl infile.legacy.vcf > outfile`

### Scripts that manipulate or summarize fasta files

Take a fasta file with IUPAC ambiguity codes and perform a base by base coin flip for use in programs that do not accomodate polymorphism (e.g. RAxML):

`perl coinflip_fasta_poly.pl infile.fa > outfile_coinflip.fa`

Calculate divergence and polymorphism (dxy and pi) in two co-linear fasta sequences:

`perl fa2polydiv_window_list_v3.pl fasta_file1 fasta_file2 list_of_scaffolds_to_analyze window_size_bp`

Extract a scaffold by name from a fasta file:

`perl getScaffold_samtools.pl file.fa contig_name > outfile`

Summarize the N50 and length distribution of an assembly and optionally print scaffolds greater than a particular size:

`perl getScaffLengthDist.pl fasta_file min_contig_length print_scaffolds_greater_than_min_contig_length_0_1`

Summarize the N50 and length distribution of an assembly and optionally print scaffolds smaller than a particular size:

`perl getScaffLengthDist_printsmallerthan.pl fasta_file min_contig_length print_scaffolds_greater_than_min_contig_length_0_1`

Count polymorphic basepairs based on IUPAC ambiguity codes in the fasta file. Note: not all fasta files include this coding of polymorphic bases.

`perl count_IUPAC_poly_single_genome.pl fasta_file > outfile`

Count polymorphic basepairs based on ambiguity codes in a fasta file using a given window size:

`perl count_IUPAC_poly_single_genome_window.pl fasta_file window_size_inbp > outfile`

Summarize base composition in a list of windows:

`perl count_basecomp_window.pl list_of_windows_to_analyze.bed genome.fa path_to_fastahack > outfile`

Remove a list of sequences from a fasta file:

`perl remove_record_list_from_fasta.pl list_of_sequences_to_remove fasta_file`

_Writes an outfile with these sequences removed named fasta\_file.filtered_

Translate a fasta file in a give frame:

`python seqs_processor_and_translator_bin_V118_AGCT.py inputfile.fa outputfile_prefix DNA 1 1 NOBIN 20`

Generate reverse complement of a fasta sequence:

`perl rev_com_v2.pl fasta_name`

_Writes out a file appended with -revcom.fa_

Extract variable sites from an aligned multi-fasta file:

```
module load biology py-biopython

python Variable_sites_extractor.py -v myfasta.fa -o myvariablefa.fa
```

_For first time use of this script you may need to install progressbar2:_

`pip install progressbar2`

Split fasta files into shorter segments as specified by a bed file:

`perl split_fastas_list_by_bed_windows.pl fasta_file_list windows.bed`

_Split results will be written to individual files called focal\_file\_start\_stop.fa_

Extract sequences from a bed file and format for RAxML analysis:

`perl extract_chrom_to_raxml_list.pl list_of_fasta_files random_select_file.bed outfile_tag_name`

Format of fasta file:

file1.fa\n file2.fa\n

_Writes out a file called: outfile\_tag\_name\_concatenated.fasta_

### Scripts for QTL and admixture mapping

See also: QTL mapping (https://openwetware.org/wiki/Schumer_lab:_Commonly_used_workflows#Ancestry_tsv_files_to_rQTL_input) and Admixture mapping (https://openwetware.org/wiki/Schumer_lab:_Commonly_used_workflows#Ancestry_tsv_files_to_admixture_mapping_input)

Match a phenotypes file with a genotypes and hybrid index file:

`perl match_phenotypes_names_with_genotypes_and_index_file.pl phenotypes_file genotypes_file hybrid_index`

*Write output files appended with modified names indicating matching*

Scripts for admixture mapping with binomial traits, see workflows here also: Admixture mapping (Commonly_used_workflows: Ancestry_tsv_files_to_admixture_mapping_input)

`Rscript perform_glm_admixture_mapping_v2_binomialtrait.R genotypes_file hybrid_index_file phenotypes_file focal_column_number name_tag`

`Rscript perform_glm_admixture_mapping_v2_gaussian.R genotypes_file hybrid_index_file phenotypes_file focal_column_number name_tag`

Convert the output of the glm scripts to easy input for plotting in R:

`perl convert_birchmanni10x_mapping_output_manhattan_plot_input.pl mapping_results.txt`

### Scripts to calculate hybrid index and make ancestry plots post ancestryinfer
Generate hybrid index file for a pair of ancestry-tsv files:

`perl parsetsv_ancestry_v2.pl ancestry-probs-par1.tsv ancestry-probs-par2.tsv`

*Writes hybrid index and heterozygosity per individual*

Generate hybrid index file for 3 way ancestry HMM files (six files total):

`perl parse_3way_tsv_ancestry.pl file_name_stem posterior_probability_threshold`

Example:

`perl parse_3way_tsv_ancestry.pl allchrs.tsv 0.9`

*Writes the proportion of the genome from each parent and ancestry heterozygosity of each type per individual*

Generate plots of ancestry per individual:

First convert ancestry-tsv files to hard called genotypes:

`perl parsetsv_to_genotypes_Dec2017_v2.pl ancestry-probs-par.tsv ancestry-probs-par2.tsv genotypes_file_name.txt`

Then get list of individuals in genotypes file:

`cut -f 1 genotypes_file_name.txt > id_list`

Then create plots.

For a single chromosome:

`perl genotypes_to_ancestry_plot.R genotypes_file_name.txt scaffold_name id_list`

Or for several chromosomes that will all be in the same pdf:

`perl genotypes_to_ancestry_plot_wChrList.R genotypes_file_name.txt scaffold_list.txt id_list`

### Scripts to average by site or in windows

Write average ancestry per site for every site in two ancestry tsv files (par1 and par2 files):

`perl calculate_avg_ancestry_by_site_tsv_files_v5.pl infilepar1 infilepar2 num_ind_thresh posterior_prob_thresh path_to:transpose_nameout.pl > outfile`

*Write average ancestry per gene using the ancestry by site file generated with calculate_avg_ancestry_by_site_tsv_files_v5.pl:*

`perl average_ancestry_by_gene_list.pl ancestry_by_site_file gene_list gtf_file`

Writes output to average_ancestrybysitesfilename_ancestry_intervals.bed

Write average ancestry per site for every site in an MSG genotypes file:

`perl parse_genotypes_ancestry_bysite.pl genotypes_file path_to:transpose_nameout.pl`

*Writes average ancestry per site to an output file named after the input file and appended with ancestry_by_site*

Average ancestry by site files in windows specified in a bed file. See workflow here: (Commonly used workflows: Ancestry tsvs to average ancestry)

`Rscript average_ancestry_bybedintervals.R ancestry_by_site_file bins.bed chr_name`

Newer version runs the script across all chromosome:

`Rscript average_ancestry_bybedintervals_v2.R ancestry_by_site_results bins.bed`

*This script takes the output of calculate_avg_ancestry_by_site_tsv_files_v5.pl and averages them in bed intervals and outputs average ancestry per window. The outfiles will be generated by chromosome and will be called: average_[infile]ancestry[bedfile][chromosome]*

Average recombination rates output by LDhelmet in windows specified in a bed file, appropriately weighing per base pair recombination rates:

`Rscript average_LDhelmet_rates_res_file_window.R modified_LD_helmet_output chromosome_name window_size_bps chrom_length_file`

*Note the latest recombination maps can be found here:*
`/home/groups/schumer/shared_bin/Xbirchmanni_LD_recombination_map_10xgenome_March2019`

### Scripts to summarize features in windows (i.e. coding, conserved, repetitive basepairs)

Determine the number of coding, conserved, or repetitive basepairs in a set of bed-formatted windows, using features in another file (can be bed or gtf formatted):

`Rscript calculate_features_per_window.R bins.bed features.gtf_or_features.bed path_to_bedtools_bin/[0 for global install] outfile_name`

for example:

`Rscript calculate_features_per_window.R combined_0.05cM_results_TLMC_March2019.bed /home/groups/schumer/shared_bin/shared_resources/Xbirchmanni-10x_12Sep2018_yDAA6-annotation-output/xiphophorus_birchmanni_10x_12Sep2018_yDAA6.gtf 0 combined_TLMC_0.05cM_windows_addcoding`

This will calculate the number of coding basepairs per window based on the windows provided in combined_0.05cM_results_TLMC_March2019.bed

*For just conserved basepairs for all chromosomes across the bed file.*

`Rscript calculate_conservedBPs_per_window_v2.R windowBins conserved.bed`

for example:

`Rscript calculate_conservedBPs_per_window_v2.R xbir_allChrs_100kb_windows xma_mostconserved_bases_v2_xbir.bed`

*This counts the number of conserved base pairs in a windowed bed file and creates a version with wConservedBPs.bed as the suffix.*

To count synonymous and non-synonymous basepairs in windowed bins.

`Rscript calculate_dNdSsites_per_window.R windowBins ntChanges.bed`

*This creates a version with wSynNonsyn.bed as the suffix. Where the last column is the non-synonymous count and the 2nd to last column is the synonymous count.*

### Scripts to identify recombination transitions

Use a genotypes file produced from ancestry-par1 and ancestry-par2 tsv files to identify ancestry transitions in each individual and the interval over which they occur:

`Rscript identify_intervals_10x_genomes.R genotypes_file_name path_to:transpose_nameout.pl`

*Writes out outfiles for each chromosome named: genotypes_file_name_chr_intervals. Each line is a ancestry transition event with the start and stop of the interval and the individual it was detected in:*

```
1407648 1801361 2

3158148 3306561 2

12714653 12933987 2

18340742 18501540 2
```

*Note that this is coded for the X. birchmanni 10x genome and would need to be modified to provide chromosome names for other genome versions*

### Scripts relating to ancestry tracts

This script takes a genotypes file and outputs inferred ancestry tracts for a given chromosome. The individual id list allows you to match the tract to the individual it came from:

`Rscript genotypes_to_ancestry_tract_lengths.R genotypes_format_file.txt chromosome_name_to_run individual_id_list_with_header outfile_name`

*Note: the individual id list can be generated using the following steps: head -n 1 genotypes_format_file.txt > id_list_tmp perl transpose_nameout.pl id_list_tmp*

This script prints inferred recombination intervals based on a genotypes file. It will output one file per chromosome in the 10x genome, appended with the word “_intervals”:

`Rscript identify_intervals_10x_genomes_v2.R genotypes_format_file.txt path_to:transpose_nameout.pl`

This script takes the file output by genotypes_to_ancestry_tract_lengths.R and converts it to estimated length in cMs, given a modified LDhelmet bed outfile. These modified bed LDhelmet outfiles can be found here on Sherlock:
`/home/groups/schumer/shared_bin/shared_resources/Xbirchmanni_LD_recombination_map_10xgenome_March2019/block_penalty_`

`Rscript convert_ancestry_tract_lengths_to_cM_lengths.R tract_length_file_chr rec_file_chr`

*The output file is the input tract length file name with “cM_lengths” appended. Examples of file formats can be found at the top of the script*

### Scripts to summarize genetic distance or admixture LD decay

Generate bed intervals for windows of a certain genetic size in cMs for the X. birchmanni 10x recombination map:

`Rscript convert_LD_rate_to_cM_rate.R chromname-bir10x window_size_cM`

*Generates an output file called: cM_windows_window_size_xbirchmanni10x_chr*

Generate files based on plink output for approximate dating of time of admixture using the decay of admixture LD (see Commonly used workflows for details on generating these files and workflows):

Add distance between the two markers in cMs to a plink-generated LD file:

`Rscript convert_plink_LD_genetic_distance_xbir10x.R plink_ld_file focal_chromosome_name`

*Outputs a file called: plink_ld_file_chrom_cMdistances*

these results can then be summarized by running additional scripts. To average LD metrics in windows of a particular cM distances:

`Rscript convert_cM_version_plink_admixture_decay.R plink_ld_file_cMdistances window_size_cM`

*Outputs a file called: plink_ld_file_chrom_cMdistances_admixture_ld_decay_cMdist*

these results can then be used to estimate age of admixture and plot decay in D:

`Rscript estimateage_plink_admixture_decay.R plink_output_ld_decay_D.ld_admixture_ld_decay`

### Scripts that generate or use insnp files

Generate a masking insnp file based on a vcf file. For example, use GATK's variant selecting tools to identify variants that have high or low depth for masking GATK Variant Selector (https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php)

`perl gatk_vcf_to_masked_insnp.pl gatk.g.vcf`

*writes an outfile called infile.mask.insnp*

Generates an insnp file to update and mask a fasta file. Variants that pass the quality thresholds are marked for updating and variants that fail the quality thresholds as well as invariants of low coverage/quality are marked for masking:

`python insnp_v9_gatk3.4_gvcf.py file.g.vcf file.g.vcf.insnp gq_cutoff dpcutoff mq_cutoff qd_cutoff fs_cutoff sor_cutoff mqrs_cutoff readpos_cutoff indel_window`

All but the indel window parameter (which masks basepairs around indels) are GATK parameters so see more details at their website GATK hard-call parameters (https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set)

Currently suggested parameters for swordtails:

`20 10 40 10 10 4 -12.5 -8.0 5`

### Scripts for admix'em simulations

See this section of the wiki for all scripts related to admix'em simulations: Admixture simulations with admix'em (Important Simulation tools: Admixture simulations with admix'em)

### Scripts for SLiM simulations or SLiM output

This script takes a SLiM vcf and generages input files for running admixtools:

`perl SLiM_vcf_to_admixtools_format.pl SLiM_format.vcf`

*This script takes a SLiM vcf and generages input files for running admixtools. The file names will be the input vcf name appended with .geno, .snp, and .ind*

### Miscellaneous

Script to link locally to a list of files elsewhere on the server:

`perl generate_local_link_file_list.pl infile_list`

*The infile list should contain the full path to the files and have one file per line*

Script to remove a list of files, use with caution!

`perl rm_list.pl infile_list`

*The infile list should contain the full path to the files and have one file per line*

Cancel a list of slurm jobs. The list can either be of slurm job ids or slurm file:

`perl slurm_cancel_jobs_list.pl slurm_list`

Submit a slurm script recursively for every chromosome in the birchmanni 10x genome:

`perl submit_all_10x_chroms_shell.pl name_of_starting_chrom slurm.sh`

Retain only every n windows from a bed file:

`perl thin_to_every_n_windows.pl file.bed retain_every_n_windows[integer]`

*Writes an output file called file.bed_thinned_by_n.bed*

Script to rename X. birchmanni 10x scaffold names as their corresponding chromosome names:

`perl replace_bir10x_names_with_chrom_names.pl file_to_rename`

*Writes an output file named file_to_rename.chromnames*

Estimate genetic distance between markers in the X. birchmanni 10x genome and optionally thin by physical and genetic distance:

`Rscript calculate_genetic_distance_and_thin_adjacent_markers.R markers_to_thin cM_threshold physical_threshold`

For example:

`Rscript calculate_genetic_distance_and_thin_adjacent_markers.R ACUA_markers_to_thin 0.25 2000`

the input markers for thinning should be in this format:

```
Scaffold_name\tpos1

Scaffold_name\tpos2

Scaffold_name\tposn
```

This script can be run on two completely co-linear genomes to document all non-polymorphic and non-missing sites that differ between them:

`perl identify_AIMs_two_genomes.pl fasta1 fasta2 > outfile`

Calculate base composition in a fasta file. This script takes a fasta as input and prints out base composition for each contig/sequence:

`perl basecompl.pl fasta > outfile`

This script takes a PWM produced by zf.princeton.edu and converts to the PWM format used by MEME and FIMO:

`Rscript motiflogo.R input.PWM.txt species_name`

### Need to be annotated

`gvcf_to_pseudo_fasta.pl`

`random_flip_samtools_low_coverage_vcf_for_PCA.pl`