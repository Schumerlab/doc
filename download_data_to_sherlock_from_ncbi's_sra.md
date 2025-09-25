# Download data to Sherlock from NCBI's SRA

## Get and convert data

NCBI's SRA is a huge raw data database that can be very useful. To download data directly to Sherlock:

1) Find the dataset you want to download on the [SRA https://www.ncbi.nlm.nih.gov/sra]
2) Find the run id under Runs
For example, in this [link https://www.ncbi.nlm.nih.gov/sra/SRX4394151] the run id is: SRR7525606
3) On Sherlock, navigate to the directory you want to download data to.
construct the appropriate wget command based on the SRA run id.
For example, for SRR7525606 the command would be:

`wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR752/SRR7525606/SRR7525606.sra`

and enter it into the command line

4) Convert the SRA file into fastq.gz files using fastq-dump from the sratoolkit

`/home/groups/schumer/shared_bin/fastq-dump SRR7525606.sra --gzip --split-
files`