# ONT sequencing and data acquisition

## Background

The lab does 3 primary workflows to achieve high quality, near T2T assemblies:

[Ligation Sequencing Kit](https://nanoporetech.com/document/genomic-dna-by-ligation-sqk-lsk114?device=PromethION): For standard ONT protocols. Low DNA input requirements (1ug in 48ul), moderate read lengths (20-40kb N50s), and high throughput. Library prep takes ½ day (requires several bead cleanups + time to resuspend) and 24 hours sequencing for a high coverage swordtail genome (~30-40Gb data). You can expect ~100-150ish Gb data off a flowcell if you do washes.

[Rapid Sequencing Kit](https://nanoporetech.com/document/rapid-sequencing-sqk-rad114?device=PromethION): For quick and dirty sequencing. Very low DNA input requirements (100ng DNA in 10ul), moderate read lengths (10-25kb N50s), and low throughput. Library prep takes 1hr and 24 hours sequencing for a high coverage swordtail genome. You can expect ~50-80ish Gb data off a flowcell if you do washes. Ligation is generally preferred because more economical use of flowcells.

[Ultralong Sequencing Kit](https://nanoporetech.com/document/ultra-long-dna-sequencing-kit-sqk-ulk114): For use cases involving very long reads. High DNA input requirements (10-40ug in 150+ ul), long read lengths (50-100kb N50), and low throughput. Library prep takes ½ day (plus 24 hours to resuspend) and 3 days of sequencing (with washes) for a high coverage swordtail genome. You can expect ~50-80ish Gb data off a flowcell if you do washes.

### Ordering reagents
Put requests in Quartzy like usual, but procurement is through HHMI, who need to generate a PO, which takes time. Flowcells come in packs of 4 and kits come with enough reagents for 6 runs. We could try doing half reactions at some point.


### Receiving reagents
Kits are stored at -20C. We put most reagents in plastic tub boxes in the small freezer but keep the flow cell buffers in the large freezer.
Flowcells are stored at 4C. Immediately after being received, flowcells should be checked. If they are below warranty (5000 pores), we can get a refund (this would save the lab $800). Immediately after being received, flowcells should also be logged on the spreadsheet (see #quick_links_lab_info).

## Sequencing
### Flowcell use
While ONT says flowcells expire after 3 months, we’ve successfully sequenced using flowcells >1 year expired. The number of pores decreases over time, so try to use older flowcells first. Generally, a flowcell with >2000-3000 pores should sequence well. 500-2000 pores can be used for testing and troubleshooting, or for protocols where not as much data are needed (e.g. amplicon sequencing).

### Flowcell reuse
Flowcells can be washed and reused using the [wash protocol](https://nanoporetech.com/document/flow-cell-wash-kit-exp-wsh004?device=PromethION). After wash + storage buffer, but before putting in fridge, run a flowcell check and add a sticker with the number of remaining pores. Also, fill out the flowcell spreadsheet. We save and reuse flowcells with > 500 pores. Flowcells should be stored in the storage bags in the 4C.

**Note**: theoretically, the DNAse in the washmix gets rid of 99.9% of old DNA but we have not tested this. Not all analyses will be sensitive to contaminant DNA, but if yours are (eg genome assembly), consider using a flowcell that previously sequenced a diverged species (eg northern swordtail vs platyfish). With this level of divergence (2%), it should be very easy identify contaminant reads. More closely related species can be used too, but we recommend not using the same species. A formal de-contamination pipeline is in the works. Using the barcoding kit is another way to avoid contamination.

### P2 Solo 
We have a P2 Solo, which can run 2 flowcells at once. This machine cost $25,000 so please be careful with it.

### Computer
The computer is a powerful/beautiful but finicky machine. This is mostly because of what we’re asking it to do. It runs some of the latest Nvidia GPUs, but is on the Ubuntu 22.04 OS. This causes some problems with the drivers and prevents live basecalling. We currently run the sequencer without basecalling (outputting pod5 raw data) and then basecall on Sherlock. This precludes analyses requiring adaptive sampling. We will update to Ubuntu 24.04 soon, which will hopefully solve these problems. We’ve been getting help from Eric Campbell from IT when the issues are too complex for us to address. The computer has a UPS which prevents restarts after short power interruption, which messes up the drivers, and also has a wired ethernet connection.

**Note**: if you turn on the P2, but it doesn't light up and minknow doesn't see it, there might be a problem with the thunderbolt controller being enabled (BIOS menu). Restarting the computer might fix the issue.

## Data acquisition and management

### File naming conventions
MinKnow is fairly user friendly. Please name your experiment with the following convention:

YEARMODY_INITIALS_KIT[LSK/RSK/ULK]_experiment (**eg 20250925_TOD_LSK_xbir-fiberseq-brain**)

Then name the SAMPLE with the standard swordtail naming convention, plus “_KIT” (e.g.
**xbir-COAC-S238-28-VIII-25-M01_LSK**)

You can join existing experiments too. This will allow you to sequence multiple samples/ individuals in the same project, which will keep OAK less cluttered.
The 4TB or 2TB data drives can be used for temporary data storage while sequencing. Generally, 100Gb of bases = 1TB of data (pod5 files are heavy).

### Data upload to oak
From the terminal on the computer, copy data to oak. This usually takes about an hour for a high coverage swordtail genome.
```
scp -r /media/4tb-data/minknow/data/20250925_TOD_LSK_xbir-fiberseq-brain/xbir-COAC-S238-28-VIII-25-M01_LSK/ <sunetID>@dtn.sherlock.stanford.edu:/oak/stanford/groups/schumer/data/Nanopore_data/
```
Eventually, we will figure out a system to back up all pod5 files to elm or some other long-term storage.

### Basecalling
On sherlock, we run the dorado basecaller in the super accuracy mode (with modified bases enabled if PCR-free). By default, we filter the reads by q10 and length 5000bp in the fastq.gz.

`dorado` is installed in shared_bin. Dependencies and other useful tools can be easily installed in a conda environment. Finally, download the superaccurate models.
```
conda create -n nanopore -c bioconda -c conda-forge pod5 samtools seqkit
conda activate nanopore

cd /path/to/scratch/dir/basecall/
dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.0.0
dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mCG_5hmCG@v3
```

After this occurs the basecalling (in parallel) and merging the combined outputs.
```
cd /path/to/scratch/dir/basecall/

ln -s /oak/stanford/groups/schumer/data/Nanopore_data/20250925_TOD_LSK_xbir-fiberseq-brain/xbir-COAC-S238-28-VIII-25-M01_LSK .

cp /home/groups/schumer/shared_bin/Lab_shared_scripts/sub_ont_basecall_dorado.sh .

sbatch sub_ont_basecall_dorado.sh xbir-COAC-S238-28-VIII-25-M01_LSK

cp /home/groups/schumer/shared_bin/Lab_shared_scripts/sub_ont_basecall_merge.sh .

sbatch sub_ont_basecall_merge.sh xbir-COAC-S238-28-VIII-25-M01_LSK
```
After you basecall, please put the combined ubam and filtered fastq files in the nanopore directory in OAK.
```
cp xbir-COAC-S238-28-VIII-25-M01_LSK.dorado9.1_sup@v5.0.0_CpG.q10.* /oak/stanford/groups/schumer/data/Nanopore_data/
```

### Analysis
For genome assembly, the combined `<sample>.dorado9.1_sup@v5.0.0_CpG.q10.l5000.fastq.gz` file from basecalling can be used as input to `hifiasm` (with the `--ont` option specified). The assembly workflow from this point on is basically the same as with hifi data.

If you want to map reads to a reference, you can use `dorado aligner` to map either the combined ubam or fastq.gz file to an existing assembly.
