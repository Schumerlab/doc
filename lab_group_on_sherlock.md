# Schumer lab: Lab group on Sherlock

Please remember to clean up files you will not use or keep most of your temporary files on /scratch! It is really easy to forget about this and use up all of our space so make a habit of routinely removing temporary files.

## Contents
1. [Lab organization on Sherlock](#Lab-organization-on-Sherlock)
- 1.1 [Overall organization](#Overall-organization)
- 1.2 [Where to find specific files or programs](#Where-to-find-specific-files-or-programs)
- 1.3 [Getting started](#Getting-started)
2. [Storing data and data backup](#Storing-data-and-data-backup)
3. [Interactive jobs](#Interactive-jobs)
4. [Common programs and dependencies and how to load them](#Common-programs-and-dependencies-and-how-to-load-them)
- 4.1 [R](#R)
- 4.2 [perl](#perl)

## Lab organization on Sherlock

### Overall organization

Files and folders containing active scripts or small data files for analyses can be kept in:

`$GROUP_HOME`

which is the equivalent path to:

`/home/groups/schumer`

This directory does not have much space (1 TB shared) so should not be used for big files of any kind. A better place for that is:

`$GROUP_SCRATCH`

which is the equivalent path to:

`/scratch/groups/schumer/Molly`

$GROUP_SCRATCH files are automatically deleted after being inactive for 6 months so make sure to backup on $GROUP_HOME or $OAK if you need files or scripts long-term

Long term data storage and large files should kept in our lab directory:

`$OAK`

which is the equivalent path to:

`/oak/stanford/groups/schumer`

### Where to find specific files or programs

Genome assemblies, annotation files, recombination maps and other shared resources can be found here:

`/home/groups/schumer/shared_bin/shared_resources`

Shared scripts and the lab git repository are here:

`/home/groups/schumer/shared_bin/Lab_shared_scripts`

Commonly used programs in the lab that are not available via module loading can be found here:

`/home/groups/schumer/shared_bin`

To make these programs globally available, you can export the path to this folder:

`export PATH="/home/groups/schumer/shared_bin:$PATH"`

**To avoid doing this every time you login to Sherlock you can edit you .bashrc file, which can be found here:**

`/home/users/$USER/\.bashrc`

and add the export command below the line that says:

`# User specific aliases and functions`

export PATH="/home/groups/schumer/shared_bin:$PATH"

In addition to the above command you might also want to automatically load certain modules, e.g.:

`module load perl`

`module load R`

### Getting started

Make your own person folder inside our lab group folder here:

`/home/groups/schumer/lab_member_folders`

## Storing data and data backup

Our lab raw data repository and repository for large files is on Oak, and can be accessed in three ways:

1) `cd /oak/stanford/groups/schumer/data`

2) `cd $OAK/data`

3) `I've put a link to this folder in our lab Sherlock directory: cd $GROUP_HOME/data`

- NOTE: Oak is not regularly backed up and cannot be seen as a permanent data repository. Other options include external hard drives and NCBI SRA. It is possible to set the release data on the SRA so that release can occur after publication if necessary, this can be edited any time.
- For directions on how to deposit on the SRA, see this file in dropbox: `depositing_data_to_NCBIs_SRA.txt`

## Interactive jobs

Before running any programs from the login node on Sherlock, remember to request an interactive job:

you can use our dedicated nodes:

`srun --pty -p schumer -t 0-2 bash`

or request general nodes:

`srun --pty -t 0-2 bash`

You can also submit jobs, this is covered elsewhere Submitting slurm jobs

## Common programs and dependencies and how to load them

Many programs on Sherlock are available as modules. To see the modules that are available, type:

`module avail`

Lots of biology specific modules have to be loaded after loading the general biology module:

`module load biology`

Any module can be loaded by typing:

`module load $MODULE_NAME`

after which the commands will be available globally

### R

You will probably want to load R interactively on Sherlock at some point, to do quick analyses of your data or run analysis/generate figures on a dataset too large to load quickly on your Desktop.

To access R, first load the module:

`module load R`

To run R from the command line simply type R then enter.

To quit, type: `q()`

### perl

Perl is globally available on Sherlock, but to install packages you need to load the perl module to use cpan:

`module load perl`

Many scripts we have in the lab depend on certain perl modules, they are listed below:

`Math::Random`

`List::Util`

`List::MoreUtils`

You can install perl modules by running the following from your terminal:

`cpan MODULE`

for example:

`cpan Math::Random`
