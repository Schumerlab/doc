# Schumer lab: Where can I find commonly used files

## Contents
1. [Where can I find shared files and resources?](#Where-can-I-find-shared-files-and-resources)
2. [Where can I find up-to-date files for local ancestry inference](#Where-can-I-find-up-to-date-files-for-local-ancestry-inference)
3. [Which genome assembly and annotation file should I use?](#Which-genome-assembly-and-annotation-file-should-I-use)
4. [Where can I find i5 and i7 barcode files?](#Where-can-I-find-i5-i7-barcode-files)
  - 4.1 [Tn5 libraries](#Tn5-libraries)
  - 4.2 [Quail libraries](#Quail-libraries)
  - 4.3 [RNAseq libraries](#RNAseq-libraries)
5. [Which recombination map should I use?](#Which-recombination-map-should-I-use)
6. [Where can I find previously run ancestry.tsv files?](#where-can-i-find-previously-run-ancestrytsv-files)
7. [Where can I find the X. birchmanni genome split into windows?](#Where-can-I-find-the-X-birchmanni-genome-split-into-windows)

## Where can I find shared files and resources?
Most genome assemblies, recombination maps, simulation tools, annotation files, etc can be found in our shared_resources folder on Sherlock (also available in Dropbox):

`/home/groups/schumer/shared_bin/shared_resources`

The most important resources are outlined below

## Where can I find up-to-date files for local ancestry inference
The most up-to-date files for running local ancestry inference can be found here:

`/home/groups/schumer/shared_bin/shared_resources/Base_ancestry_run_files_August2020`

This folder is set up so that you can copy all of its contents to the directory on slack where you are running local ancestry inference. For example:

`cp /home/groups/schumer/shared_bin/shared_resources/Base_ancestry_run_files_August2020/* ./`

## Which genome assembly and annotation file should I use?
For *X. birchmanni* the current assembly is:

`/home/groups/schumer/shared_bin/shared_resources/xiphophorus_birchmanni_10x_12Sep2018_yDAA6.fasta`

and the annotation files are:

`/home/groups/schumer/shared_bin/shared_resources/Xbirchmanni-10x_12Sep2018_yDAA6-annotation-output`

For *X. malinche* the current assembly is:

`/home/groups/schumer/shared_bin/shared_resources/xiphophorus_malinche_10x_14Nov2018_HyZ96.fasta`

and the annotation files are:

`/home/groups/schumer/shared_bin/shared_resources/xiphophorus_malinche_14Nov2018_HyZ96_annotation_output`

For *X. maculatus* the current assembly is:

`/home/groups/schumer/shared_bin/shared_resources/xma_washu_4.4.2-jhp_0.1_combined-unplaced-mito.fa`

and the annotation files are:

`/home/groups/schumer/shared_bin/shared_resources/Xiphophorus_maculatus_LG.Xipmac4.4.2.81.gtf`

## Where can I find i5 and i7 barcode files?

### Tn5 libraries
For Tn5 libraries, i5 barcode sequences can be found in box under this path:

`Schumer_lab_resources/Extractions_and_library_preps/Tn5_i5_indices.xls`

For Tn5 libraries, i7 barcode sequences can be found in box under this path:

`Schumer_lab_resources/Extractions_and_library_preps/Tn5_i7_indices.xls`

### Quail libraries

Quail libraries only have an i7 index, which can be found here:

`Schumer_lab_resources/Extractions_and_library_preps/Quail_FC1_index_sequence.xls`

## RNAseq libraries

i5 and i7 indices for kapa RNAseq libraries can be found here:

`Schumer_lab_resources/Extractions_and_library_preps/IDX Unique Dual Indexes_forkapa_RNAseq.xls`

## Which recombination map should I use?

- The current LD recombination map for X. birchmanni can be found here:

`/home/groups/schumer/shared_bin/shared_resources/Xbirchmanni_LD_recombination_map_10xgenome_March2019`

- The current crossover map for F2s can be found here:

`/home/groups/schumer/shared_bin/shared_resources/Xbirchmanni_Xmalinche_F2_map_March2019`

- The current crossover map for all artificial hybrids can be found here:

`/home/groups/schumer/shared_bin/shared_resources/Xbirchmanni_Xmalinche_all_artificial_hybrids_map_March2019`

## Where can I find previously run ancestry.tsv files?

Results of previous AncestryHMM runs can be found here:

`/oak/stanford/groups/schumer/data/Ancestry_tsv_results_files`

Make sure to put your matching configuration files here too!

## Where can I find the X. birchmanni genome split into windows?

Bed formatted files of the X. birchmanni reference genome split into windows and with recombination rate, number of coding and conserved basepairs, and number of synonymous and non-synonymous basepairs:

`/home/groups/schumer/shared_bin/shared_resources/xbir_genome_windowed/`

Window sizes are: 5kb, 10kb, 50kb, 100kb, 250kb, 500kb, 1Mb, 0.05cM, 0.1cM, 0.25cM, 0.5cM, 1cM

For basepair windowed data there are three columns: scaffold, start, end

Files annotated with extra information are `*recRate_codingFeats_wConservedBPs_wSynNonsyn.bed`

With eight columns: scaffold, start, end, mean recombination rate, number of coding bps, number of conserved bps, number of synonymous bps, and number of non-synonymous bps

For Centimorgan windowed data there are five columns: scaffold, start, end, mean recombination rate in that window, number of SNPs

Files annotated with extra information are `*codingFeats_wConservedBPs_wSynNonsyn.bed`

With nine columns: scaffold, start, end, number of SNPs, mean recombination rate, number of coding bps, number of conserved bps, number of synonymous bps, and number of non-synonymous bps
