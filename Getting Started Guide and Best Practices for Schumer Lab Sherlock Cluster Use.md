# Getting Started Guide and Best Practices for Schumer Lab Sherlock Cluster Use

- 1.1 [Introduction](#Introduction)
- 1.2 [Sherlock](#Sherlock)
- 1.3 [Sherlock Partition and Node Specifications](#Sherlock-Partition-and-Node-Specifications)
- 2.1 [Storage](#Storage)
- 2.2 [Storage locations on sherlock](#Storage-locations-on-sherlock)
- 3.1 [Logging on](#Logging-on)
- 4.1 [Submitting a job](#Submitting-a-job)
- 4.2 [Making your jobs run faster](#Making-your-jobs-run-faster)
- 4.3 [Increasing the RAM](#Increasing-the-RAM)
- 4.4 [Increasing the cores used by your jobs](#Increasing-the-cores-used-by-your-jobs)
- 4.5 [Cores vs RAM](#Cores-vs-RAM)
- 5.1 [How do I actually start submitting jobs?](#How-do-I-actually-start-submitting-jobs?)
- 5.2 [Submission Scripts](#Submission-Scripts)
- 5.3 [Output](#Output)
- 5.4 [Using multiple threads](#Using-multiple-threads)
- 5.5 [Notifications](#Notifications)
- 5.6 [Understanding how long your jobs take to run](#Understanding-how-long-your-jobs-take-to-run)
- 5.7 [Loading software from modules](#Loading-software-from-modules)
- 5.8 [Change into the directory you want to be in](#Change-into-the-directory-you-want-to-be-in)
- 5.9 [Run the commands you want to run](#Run-the-commands-you-want-to-run)
- 5.10 [Example SLURM script](#Example-SLURM-script)
- 6.1 [Interacting with SLURM](#Interacting-with-SLURM)
- 6.2 [Sending submission scripts to the scheduler](#Sending-submission-scripts-to-the-scheduler)
- 6.3 [Monitoring your job while it runs](#Monitoring-your-job-while-it-runs)
- 6.4 [Need to delete a job (#Need-to-delete-a-job)
- 6.5 [Interactively log on to a node](#Interactively-log-on-to-a-node)
- 7.1 [Transferring Files](#Transferring-Files)
- 8.1 [Software Installation and Modules](#Software-Installation-and-Modules)
- 8.2 [Conda package installation](#Conda-package-installation)
- 8.3 [Accessing frequently used modules or packages](#Accessing-frequently-used-modules-or-packages)
- 8.4 [Installation Requests](#Installation-Requests)
- 9.1 [Common workflows](#Common-workflows)
- 10.1 [Editing Files](#Editing-Files)
- 11.1 [Reproducibility](#Reproducibility)
- 11.2 [File naming and links](#File-naming-and-links)
- 11.3 [General Reproducibility Practices](#General-Reproducibility-Practices)
- 12.1 [Schumer Lab Best Practices](#Schumer-Lab-Best-Practices)
- 13.1 [Troubleshooting Common Errors](#Troubleshooting-Common-Errors)

# Introduction 
Hello! This document is intended as a guide to help new users from the Schumer lab access the computing resources available at Stanford. Also included are some best practice guidelines to standardize computing practices, encourage reproducibility, and help maintain equal access to these resources for all users. This guide is a compilation of specific information about sherlock’s computing resources primarily gathered from the [Sherlock Documentation](https://www.sherlock.stanford.edu/docs/) as well as computing knowledge the author has gathered through the years. It is a living document and is thus subject to change. 

Additionally, users should note that this is not an introduction to bioinformatics or working in Unix. All users should become acquainted with the basics of these skills - I personally recommend [Practical Computing for Biologists](https://www.amazon.com/Practical-Computing-Biologists-Steven-Haddock/dp/0878933913) as it is a great tool to become familiar with navigating computing through a command line interface. This guide will seldom go into detail on specific commands so if the reader desires to know more about a specific command in the text, please turn first to Google or try accessing the command’s help page (e.g. man exampleCommand or exampleCommand -h).

## Sherlock 
There is one main computing cluster available to us through Stanford for our use: Sherlock. Sherlock, one of the top 500 computers IN THE WORLD!!, currently consists of approximately 2,100 compute nodes providing 75,000 CPU cores and 1,200 gpus. If you have been added as a member of the Schumer Lab on sherlock (please contact Molly Schumer to be added), then you have access to **normal** nodes, nodes dedicated to labs belonging to Humanities and Sciences (the **hns** partition), an exclusive set of nodes reserved for labs who have paid extra for dedicated nodes (the **owners** partition), as well as the Schumer lab dedicated nodes (the **schumer** partition). This entitlement means that our jobs have priority when running on the **schumer** partition over those that belong to people not in our user group. If other users not in our group are running jobs on the **schumer** partition, and you request the resources they are using, your jobs will preempt their jobs and their jobs will re-enter the SLURM scheduler to be re-assigned to a new node. Conversely, if you are running a job on the **owners** partition, your job will be preempted if a user from the owner’s lab request a job requiring the resources you are using. Your job will be cancelled and resubmitted to the SLURM scheduler to await the next open node that matches your job’s requirements. You are able to submit jobs to the **normal**, **hns**, **owners**, and **schumer** partitions (and possibly more! e.g., **gpu** and **bigmem**), but your priority access is limited to jobs on **schumer**. 

Below I have summarized the memory and node specifications for each of these partitions and their nodes. 

## Sherlock Partition and Node Specifications
To get some basic information on partition resources and job limits, you can run:

```
(base) [khunnicu@sh02-ln03 login ~]$ sh_part  
```  

| partiton    | total nodes | total CPU cores | total **gpu**s | max runtime | mem per node (GB) |
| ----------- | ----------- | --------------- | -------------- | ----------- | ----------------- |
| **normal**  | 248         | 6568            | 0              | 7d          | 128-384           |
| **bigmem**  | 11          | 824             | 0              | 1d          | 384-4096          |
| **gpu**     | 47          | 1444            | 192            | 2d          | 191-2048          |
| **hns**     | 144         | 5420            | 12             | 7d          | 128-1536          |
| **schumer** | 8           | 224             | 0              | 7d          | 191-256           |
| **owners**  | 1741        | 66272           | 980            | 2d          | 128-4096          |

# Storage 
There are a few places where you can store your data and associated files during analysis. In the following sections, I detail where you can store longer-term backed-up data, shorter-term (mid-analysis) non-backed up data, and miscellaneous small scripts.

**As a general rule, always backup your scripts used for your analyses as well as non-replaceable raw data (i.e. data you have generated but not yet stored in a public archive like [SRA](https://www.ncbi.nlm.nih.gov/sra)).** For raw data, a general rule of thumb is to have at least one additional copy of this data stored through your own personal storage and one copy on the Schumer lab $OAK directory (explained further below). Intermediate files (i.e. alignment files like .bams or .sams) seldom belong in long-term lab storage and more often belong in short-term storage, though you are free to backup any important intermediate bam files in your personal storage. On **sherlock**, we have access to a $GROUP_SCRATCH directory that is essentially limitless (100Tb shared across all lab members) short-term storage (90 days). Additionally, you have access to a small home directory $HOME where small miscellaneous scripts can be stored (~15Gb). Please note however that files in this location can only be accessed by you so avoid storing any files in your $HOME directory that other users (e.g. Molly) may need access to. An important consideration for this directory is that it is often a default location that many software packages will install themselves into. The most prevalent example of this is **miniconda** (a convenient and useful software manager that allows you to install packages that require specific dependencies). I would strongly recommend that users do not use their $HOME directories for this purpose as these packages (especially **conda**) can rapidly grow to exceed your 15Gb storage limit and having to reinstall **conda** and all associated environments elsewhere is not a fun time.

* What does this strange notation mean (e.g. $OAK)? If you have encountered this notation with a dollar sign followed by text characters when interacting with a command line environment then congratulations! You have interacted with your first variable. This syntax is actually representing a variable where a specific value you care about has been assigned to a named storage location. If you want to investigate this further, you can see that $OAK is actually short-hand for a specific directory on sherlock. To see what that directory is you can `echo` the variable name and the terminal will print out the directory it is pointing to!

```
(base) [SUNET_ID@sh03-ln04 login ~]$ echo $OAK
/oak/stanford/groups/schumer
```

To avoid causing storage problems, keep a close eye on what intermediate files you are keeping in short-term and long-term storage (when deleting think about how long it’s been since you’ve used this file, how easy would it be to regenerate, how likely are you to use it again in the next X months, etc). 

Also, for files you are interested in keeping on the cluster but are not actively using, consider whether they can be compressed so that your storage use is minimized. The most common way to compress a file is by zipping or gzipping (gzipping is more common for genome sequence data). Simply run `gzip yourFile` and if it is possible for the file to be compressed, then the space occupied by the file will decrease. Also some special file formats such as **.sam** files have unique compression formats. In this case, **.sam** files can be compressed into **.bam** files which take up much less space and can usually be directly read by downstream software. **Samtools** is a good software suite to achieve this compression type.

I recommend that individuals periodically check on how much storage space is available in our short-term and long-term storage locations and how much of the total space you are using so that we can all avoid causing problems for others. Much of this information is printed by default when you log into a new instance of **sherlock**. 

```
+---------------------------------------------------------------------------+
| Disk usage for SUNET_ID (group: schumer)                                  |
+---------------------------------------------------------------------------+
|    Directory |  volume /   limit [          use%] | inodes /  limit (use%)|
+---------------------------------------------------------------------------+
          HOME |   2.2GB /  15.0GB [|          14%] |      - /      - (  -%)
    GROUP_HOME | 949.7GB /   1.0TB [||||||||   92%] |      - /      - (  -%)
       SCRATCH |  13.5GB / 100.0TB [            0%] |   3.7K /  20.0M (  0%)
 GROUP_SCRATCH |  55.5TB / 100.0TB [|||||      55%] |   1.8M /  20.0M (  9%)
           OAK | 203.2TB / 250.0TB [||||||||   81%] |   7.6M /  37.5M ( 20%)
+---------------------------------------------------------------------------+
```

To manually assess how much space is left in a shared storage location, you’ll want to use either the `df -h /path/` command. 

```
(base) [SUNET_ID@sh03-ln04 login ~]$ df -h $GROUP_SCRATCH 
Filesystem                              Size  Used Avail Use% Mounted on
10.0.10.51@o2ib7:10.0.10.52@o2ib7:/fir  100T   56T   45T  56% /scratch
```
  
As seen above, as of the time of writing (25 Mar 2026), the short-term storage location ($OAK) on **sherlock** for the Schumer lab has 45TB still free. Be extremely wary of adding to this total when the Use% exceeds 90%, as a larger discussion about lab data storage is probably warranted. 

To determine how much storage you are using in a particular directory, navigate to the folder you are interested in measuring, then use the du -sh command to get the total storage used by all of your files and directories in that location. If you are interested in breaking this down by directory, use `du -sh *`.

```
(base) [SUNET_ID@sh03-ln04 login ~]$ du -sh
2.4G .

(base) [SUNET_ID@sh03-ln04 login ~]$ du -sh *
224K adapters
208K busco_downloads
1.1M perl5
546M scripts
```

As you can see, as of the time of writing, my home directory ($HOME; symbolized by the ~ ) currently contains 2.4 GB worth of files, and the second command allows us to see that the majority of this is taken up by my ~/scripts/ directory.

## Storage locations on sherlock

*For longer-term, backed-up storage:* As part of the Schumer group, you will have access to our $OAK directory. The entire lab $OAK directory has a limit of **250Tb**. This is where we store raw data as a lab (/oak/stanford/groups/schumer/data) and also where you can keep your longer-term files. You will need to create a folder for yourself within this directory at /oak/stanford/groups/schumer/data/lab_member_folders/yourName. To do this, use the `mkdir /oak/stanford/groups/schumer/data/lab_member_folders/yourName` command. The general lab convention for naming these folders is your first name in all lower case.

*For short-term (mid-analysis), not backed-up storage:* You will also be granted access to the lab $GROUP_SCRATCH directory (/scratch/groups/schumer/) when your account is created. As above, you will need to create a directory for yourself within this directory for your intermediate files (e.g. /scratch/groups/schumer/kelsie/). This location has a total **limit of 100Tb** of storage shared across all lab members. This is where you can keep your larger intermediate files. An additional storage location for your short-term storage needs is a personal scratch directory ($SCRATCH). By default, most lab sherlock users prefer to use $GROUP_SCRATCH over their personal $SCRATCH. As with $GROUP_SCRATCH, $SCRATCH has a storage limit of **100Tb** for 90 days although you do not share this storage with anyone. An important consideration here is that others will not be able to access files in your personal $SCRATCH so keep that in mind if others (e.g. Molly) will need to access these files. **IMPORTANT NOTE:** The sherlock system automatically deletes files within all scratch directories that exceed 90 days from time of creation (aka **the scratch purge**). It is important that you both (1) backup files you are likely to continue needing for future analyses and (2) that you do not intentionally modify the time stamps associated with files to prolong their duration on $GROUP_SCRATCH (if you do not know what I am referring to here, then you are unlikely to accidentally modify the time stamps of your files). Modification of time stamps on files can trigger the sherlock managers to flag our lab group as suspicious and potentially result in negative consequences for our lab’s use of sherlock. 

*Special cases of long-term storage:* As a group, we have access to one additional shared storage directory: $GROUP_HOME. This directory is unique in that it is capped at 1Tb total across all lab members, and there is no way for the lab to gain additional space beyond this 1Tb cap. This is not a lot of storage space! As a result, **we strongly discourage storing any but the most essential shared files in this location**. An example of a key directory is: /home/groups/schumer/shared_bin. This is a directory that houses many packages that lab members frequently use across projects. You will notice that within the $GROUP_HOME directory that there is a directory of lab member folders at /home/groups/schumer/lab_member_folders. This is not a space to store data or intermediate files given the restricted storage limits in this directory. We frequently hit our 1Tb storage limit due to users accidentally (or intentionally) writing intermediate files to this location which causes very bad problems for the lab! We recommend keeping an eye on any files you write to this location. **IMPORTANT NOTE:** There is a directory within $OAK which appears to be a safe place to write intermediate files:  /oak/stanford/groups/schumer/schumer/lab_member_folders. However, you have been deceived! This location is actually just a link to /home/groups/schumer/lab_member_folders, and any files written to this $OAK location will actually count against our $GROUP_HOME storage limit. Please pay close attention to the section on $OAK storage above for the appropriate path at which you can set up your own long-term storage folder.

# Logging on 
To log on to the sherlock, you will be using `ssh` through a command line interface of your choice (Mac users use terminal; Windows users will need to do some research to figure out what the appropriate equivalent is). The general format of the ssh command is `ssh SUNETID@login.sherlock.stanford.edu`. 

The usernames correspond to your SUNET ID, and by default, your password is your SUNET password. 

If you are using a GUI to interact with your files (such as Cyberduck or others; see **transferring files**), the address from your ssh command is what you will log in with **using SFTP** (cyberduck default is FTP). Your username and password are also the same.

Want to minimize the number of multi-factor authentication requests you have to complete when logging on to sherlock? Check out this [helpful guide](https://www.sherlock.stanford.edu/docs/advanced-topics/connection/?h=#avoiding-multiple-duo-prompts)!

# Submitting a job 
Congratulations! You have now logged on to the cluster(s) and are ready to run a job.

What do we mean by job? A job is a singular invocation of a command or script that you are interested in running on the cluster. For example, you want to align raw whole genome sequencing data from one cottontail individual, SYLV001, to the European rabbit genome. This action would be performed by invoking the following command (this is not a tutorial on how to use bwa-mem2 mem; I accept no liability if the following example breaks upon your own use): 

```
bwa-mem2 mem -t 8 -R '@RG\tID:sample\tSM:sample.1' \
/scratch/groups/schumer/kelsie/genome_links/ocun_genome.fna \
/scratch/groups/schumer/kelsie/post-bwa/SYLV001.fastq > SYLV001.sam
```

Sure, you could sit on one of the login nodes (the default nodes you log onto; e.g. SUNET_ID@sh03-ln04 login), move to your relevant directory, copy and paste the above, and it will run. HOWEVER, all users automatically log onto the login nodes, and running this command on the login nodes will slow down anything else that runs on the head node and result in you receiving warnings from the sherlock admin. In practical terms, this means that every other user logging on to the cluster will experience slow downs when changing directories or running basic commands. You also want to avoid doing this for selfish reasons as your job will in all likelihood run slower than you want it to. Additionally, you probably have several different cottontail individuals you need to align (read: you have multiple jobs that need to run). So how do we run our jobs so that they run as fast as possible and disrupt other users as little as possible?

## Making your jobs run faster 
This can generally be accomplished by either increasing the RAM available to your job and/or dispersing it across multiple nodes (read: cores/threads, I will use these interchangeably throughout the document). Note there are limitations associated with both of these actions which I describe below.

## Increasing the RAM 
You can accomplish this by requesting more RAM for your jobs on **schumer** nodes (up to XXX GB of RAM) or in extreme cases, by running high-memory jobs on nodes that intrinsically have more RAM available. For example, the **bigmem** partition on Sherlock has **384-4096Gb of RAM** available for individual jobs. How much memory a job takes will depend on what software you are using and what that software is doing. Are you calling SNPs across whole genome data for 20 bunny rabbits using clunky software? Then, this might be a job type that you’d want to run on a higher memory node. Are you aligning 60+ hamster testes transcriptomes to the dwarf hamster reference genome? Then you don’t need a super high memory node - each individual job in this case is very low memory and will even run on all the typical sherlock partitions with little issue. When in doubt, try out **8Gb** or **16Gb** of RAM and adjust upward if necessary. For project with low RAM needs then, what can you do to make these jobs go faster?

## Increasing the cores used by your jobs 
Got 60+ hamster testes transcriptomes to align? Then what you need in this case is not more memory for each of the 60 jobs, but rather more processing power for each job. This is done by “multi-threading” each of your jobs. Each job occurs on a particular node, but within each node there are several cores/threads available for your use. Essentially what multi-threading means is that the node is taking the job you want done, then splitting it into smaller sub-jobs (known as tasks) that can each run independently and concurrently, making your overall job run faster. The documentation of your software of interest will go into whether or not your software can take advantage of multiple cores and will have flags to specify how many threads you want your job to use. 

Make sure you reference the available number of cores per node of interest above. You cannot exceed the number of cores available per node or use the same core twice. This also means that once you start using threads on a node, no one else can use them either. Refer to the [Best Practices](#Best Practices) section below for proper etiquette on how to proceed when you have a lot of computing you need to get done all at once.

## Cores vs RAM
**Should I increase the available RAM or increase the number of cores used?** This is going to heavily depend on the software you are using and your overall goal. Could I use **schumer** to run the above example of 8 hamster alignment jobs requesting 4 threads for each job? **That would make the most sense to get a job to run faster, right (we’ve optimized RAM and number of threads)? Well, not necessarily.** As I discussed above, aligning RNASeq data is actually a fairly low memory activity, so giving the jobs all the RAM in the world might result in a partial speed up, but in this case, it’s not going to speed up my overall project as much as being able to run more individual hamster samples at once and have each job take advantage of more cores. Running 12 jobs at one time in this case more significantly speeds up the project than running 8 jobs at one time with more memory for each. Additionally, running the same 8 individual jobs but requesting 8 threads each instead of 4 threads each is also going to be more beneficial than upping the available RAM. 

**So does that mean I should always go the route of increasing the number of cores I use instead?** 
Again, not necessarily. Some software packages (e.g. STAR, an alternate RNASeq aligner requires up to 32GB of RAM for each job; or minimap2, a long read sequence aligner) are inherently high memory. The best course of action is to consult your software documentation or try googling around to see what others have done to make their jobs run the best with that particular software. 

# How do I actually start submitting jobs? 
So you’ve decided where you want to submit your jobs and how many threads you want to request for each of them. Now it’s time to actually tell the cluster what you want it to do. Sherlock has a built in system to accept your jobs and then tell the individual nodes what to do and when. This is accomplished through a scheduler that is analogous to a human secretary. You (and everyone else on the cluster) send in a request, and the scheduler put that request into action. It will make sure a job goes to a particular node (or one of a series of nodes) that you’ve requested, and if someone is already using that node, they will place your jobs in a queue while the scheduler waits for that resource to become available (recall from above that if you are using a node, no one else can use that node while you are using it). These schedulers are also what manage entitlement queues/priority users. Recall from above that on **schumer**, we have priority over anyone not in our lab, and that when we submit jobs to **schumer**, anything running on **schumer** that isn’t Schumer lab owned will be preempted. This is done by the scheduler. **Sherlock** uses **SLURM** (Simple Linux Utility for Resource Management), but you may encounter other schedulers in your future computing efforts (e.g. **SGE**, although this has been phased out of most modern clusters). 

## Submission Scripts 
You will have learned by now that we can use scripts to communicate with our computers (If not, refer back to “[Practical Computing for Biologists](https://www.amazon.com/Practical-Computing-Biologists-Steven-Haddock/dp/0878933913)” before proceeding). When communicating with schedulers, we are again using scripts to communicate with them. They vary in content based on which scheduler your cluster uses, but the overall structure is the same. All scripts of this type will start with a block of specifications (indicated by the starting # symbol) letting the scheduler know a few things about what you want. Following this block, you will put the actual commands you are interested in running. I will begin by walking through the basic specifications you will include in the beginning, and then end this section with example scripts for SLURM that you can build yours off of.


**Every time we invoke a script, we start by letting the computer know what language we are speaking**. **SLURM** “speaks” unix/bash, so when we make scripts to communicate with our schedulers, we always want to make sure we start with this line. 

`#!/bin/bash`
  
**Now we want to give the job a name so that when we communicate with the scheduler**, both of us know what job we are talking about. Note, when displaying job names, schedulers tend to only display the first ~8 characters, so it is in your best interest to find a way to abbreviate your job title so that you can easily tell which job is which. It is also important to consistently name your jobs the same format for automation purposes (see File Naming). For this example, I am running hisat2 (an RNASeq aligner) on the hamster individual SSSS-188-4M-SP to the Phodopus campbelli genome. A good short job name for this would be th_SS-188-4M-SP_pcam. I would be able to easily see from the abbreviated name which step in the process I’m at (alignment with hisat2=th) and what sample it is (SSSS-188-4M-SP = SS-188-4M-SP). It is possible to get the full job name when requesting more information about a running job (see Monitoring your job below) but it requires more effort than you want to expend for a quick check. 

`#SBATCH --job-name=JOBNAME`

**Next, we want to tell the schedulers where (read: what partition) we want out jobs to go.** You can specify one partition or multiple partitions here. If multiple, separate each partition with a comma and no space. 

`#SBATCH -p hns,owners,schumer`

## Output
**When you run any command, it will often have output that displays while the job is running.** This output is usually written to either STDOUT or STDERR (the difference between the two is beyond the scope of this tutorial). Below I show you how to tell the scheduler to save both types of output to files that you can access while the job is running or after it is completed.

```
#SBATCH --output=th_SS-188-4M-SP_pcam.out 
#SBATCH --error=th_SS-188-4M-SP_pcam.err
```

In this case, you are telling SLURM that you want the STDERR and STDOUT output saved as different files. You have specified here what you want those files to be named. Again, be consistent; future you will thank you for it.
  
## Using multiple threads 
This is how you tell the scheduler you want four threads used on this job. If you don’t need multi-threading, these specifications are not necessary and omitting them will result in your job being allocated 1 core/cpu by default.

```
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
```

## Notifications 
The schedulers can email you with notifications about the status of your jobs. You can specify exactly what types of notifications you want. I default to getting all types of notifications so if you want more specificity turn to google. I recommend as you advance to running multiple jobs in parallel that you create an email filter to move these out of your inbox as they quickly become overwhelming. I find that filters for messages containing “sherlock Slurm” catches them all for me.

```
#SBATCH --mail-user=SUNET_ID@stanford.edu
#SBATCH --mail-type=ALL
```
  
With this last block, we have completed the specification section of the submission script. Now we need to tell the nodes what commands to run.

## Understanding how long your jobs take to run 
This is purely personal preference, but I strongly recommend you begin and end the section containing your commands with the date command. This command prints your current date and time to the outfiles we specified above. Doing this for all your jobs gives you a record of how long a particular job took (take the start time and subtract it from the finish time). From this you can start to understand how long a particular software usually takes to run which helps you plan your work more efficiently in the future and can let you know if something fishy is happening with one of your jobs if it is running longer/shorter than you would expect.

```
(base) [SUNET_ID@sh03-ln04 login ~]$ date
Wed Mar 25 01:47:14 PDT 2026
```

Additionally, if you are interested in learning more runtime details about a job you’ve already run, you can use the following command:

```
sacct -j <jobid> --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss,MaxVMSize,nnodes,ncpus,nodelist
```

This should report back to you: your UserID, your JobID, the Jobname, the partition on which it ran, the job state (completed, cancelled, failed, etc), time you submitted this command, start time, end time, time elapsed, maximum resident set size of all tasks in job, maximum virtual memory size of all tasks in job (helpful if you are getting errors suggesting you are hitting memory caps), number nodes allocated to this job, number cpus allocated to this job, and the nodes the job ran on. See [here](https://slurm.schedmd.com/sacct.html) for more info on the different flags for **sacct**.

## Loading software from modules 
For all pieces of software that you want to use for your job, you need to tell the node to load them from the modules. More information about modules and software installation are found below in the [Software Installation and Modules](#[Software Installation and Modules]) section. The script is now just talking to the node we were interested in directly (i.e. not through the scheduler) so this is just basic unix/bash. For example, if I want to use hisat2 in my script, I need to run the command below before I can use hisat2: 

```
(base) [SUNET_ID@sh03-ln04 login ~]$ module load hisat2/2.1.0
```

Having trouble finding modules that already exist? Try module spider! This gives you information on how to load the module you are interested in as well as if any other modules are required as dependencies:  

```
(base) [SUNET_ID@sh03-ln04 login ~]$ module spider hisat2

-------------------------------------------------------------------------------------------  hisat2: hisat2/2.1.0  

-------------------------------------------------------------------------------------------    Description:

      HISAT2 is a fast and sensitive alignment program for mapping next-generation sequencing reads (both DNA and RNA) to a population of human genomes (as well as to a single reference genome).

    You will need to load all module(s) on any one of the lines below before the "hisat2/2.1.0" module is available to load.

      biology
```

## Change into the directory you want to be in 

Remember we have defaulted to the working directory being the location where you submitted your submission script from. Generally, you may not want to be in that location (which should be your home directory which is limited in storage space). You will likely run out of memory mid-job doing this since this home directory is capped at 100GB.

`cd /scratch/groups/schumer/kelsie/melanoma_RNAseq/`

## Run the commands you want to run 
Self explanatory. 

```
hisat2 -p 4 -k 100 -x /scratch/groups/schumer/kelsie/genome_links/xbir-COAC-16-VIII-22-M_v2023.1.fa   \

-1 /scratch/groups/schumer/kelsie/melanoma_RNAseq/post-trimm/paired/sample_1.paired.fq.gz \

-2 /scratch/groups/schumer/kelsie/melanoma_RNAseq/post-trimm/paired/sample_2.paired.fq.gz \

-S /scratch/groups/schumer/kelsie/melanoma_RNAseq/post-trimm/post-hisat/sample_xbir.sam
```

*What’s with the weird “\” at the end of each line? For REALLY long commands, you can increase the readability and edit-ability of your script by ending lines in bash with a “\”. Bash reads this as “keep going” and interprets your more aesthetically-pleasing multi-lined file as one single command (which it is) instead of artificially truncating the command.

**End the script by printing the date!** See above.

## Example SLURM script
```
#!/bin/bash
#SBATCH --job-name=example
#SBATCH --output=example.out
#SBATCH --error=example.err
#SBATCH -p hns,owners,schumer
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=SUNET_ID@stanford.edu 
#SBATCH --mail-type=ALL

date

module load softwareOfChoice/version

cd /scratch/groups/schumer/yourDirectory/fullPath/

softwareOfChoice -F flagsNeeded -i inputFiles -o outputDirectory

date
```
  
# Interacting with SLURM

## Sending submission scripts to the scheduler 
You’ve made your submission scripts. Now you need to send them to the scheduler. From the location where your script is and where you want your STDERR/STDOUT files to go, you will submit jobs using the following commands. In this example, you’ve named your submission script hisat2_testRun.submit

`sbatch hisat2_testRun.submit`  

## Monitoring your job while it runs 
Want to see the state of your job (is it queued and waiting for space, currently running, finished)? Alternatively, want to see more information about your job (confirm the jobname, check how long it has been running, etc.)?  Run the following from anywhere on the cluster.

`squeue`

```
squeue -u username lets you see only jobs associated with that user name
squeue jobid lets you see detailed information associated with your job of interest
squeue -p schumer lets you see how busy the **schumer** nodes are at a given moment
```

## Need to delete a job (or several) in the queue or kill a job already in progress?
`scancel [#####]`
Get the job ID from squeue

## Interactively log on to a node 
Need to debug a script and want somewhere to test it for a short time? Want to set up a file transfer that’s going to take awhile? You can log directly onto the **sdev** partition and run your commands directly from the command line (no submission script needed). Alternatively, you can request a **schumer** or **hns** node if your debugging is likely to exceed 1hr or require more memory. To log off interactive sessions, either type exit while still logged on. 

```
sdev  

srun --time 4:00:00 --mem=32G --ntasks 1 -p schumer,**hns** --pty bash -i
```

# Transferring Files 
If you are interested in moving files from one machine to the cluster (i.e. from Novogene to the cluster, from your personal laptop to the cluster, etc), there are several options available to you. First, you can use a GUI such as [Cyberduck](https://cyberduck.io/download/), Filezilla, etc. with the **SFTP protocol** and the addresses listed in the Logging On section, and drag and drop your desired files of interest. 

A more reliable way to transfer files is through the **rsync** (or **scp**) tool as this gives you the best chance of having your file transfer proceed without interruptions. A final option for large transfers is the **GLOBUS** GUI which is not discussed here. A guide for file transferring specific to sherlock can be found [here](https://srcc.stanford.edu/private/oak-gateways). In your desired destination, run the **rsync** command to initiate file transfer. It is recommended that you use the **md5** (mac) or **md5sum** (linux) command on your files to ensure file integrity has been maintained. Every instance of a particular file has a unique md5 code that should be maintained after file transfer. Run this command on each file before and after transfer and they should match if the transfer was successful. File size will not necessarily match before and after transfer so use this comparison as a last resort.

# Software Installation and Modules 
This section details how to access existing installed software on the cluster, how to install your own software through package managers, and how to submit your own requests for new software/updated versions. When you access the cluster, there will already be a subset of software already installed by previous user requests. We use the module system on sherlock which means that we access software by loading the appropriate module (as opposed to linking to specific software through your path which you may have read about in *Practical Computing*). To see a list of what is already installed, run **module avail** and it will list what software is installed and what the location of the software is. You will always need to load modules (module load /pathToSoftware/) within your submission scripts so that they load on the node your job runs on.  

Additional relevant software on Sherlock may be found in /home/groups/schumer/shared_bin/. Invoking software that hasn’t been loaded yet, will result in an error message. 

## Conda package installation
Want to easily install packages not already available in the module system? Then **conda** may be the solution for you! **Conda** (or **miniconda**/**anaconda**/etc) is a package and environment manager to allow for easy distribution and installation of software. **Conda** allows users to create a “virtual environment”. For computing systems like **sherlock**, it is not uncommon that multiple versions of packages will need to be installed that are used for different scenarios. A common example of this are different versions of python. Virtual environments can help solve this problem of potentially conflicting versions of software by creating an isolated environment that allows you to use specific versions of packages and their dependencies at a given moment in time. 

To install **conda**, follow the installation guide available through [Anaconda for Linux systems](https://www.anaconda.com/docs/getting-started/miniconda/install/linux-install). **IMPORTANT NOTE:** The sherlock admins [do not officially recommend using conda on sherlock](https://www.sherlock.stanford.edu/docs/software/using/anaconda/?h=conda#introduction). As a result of this, sherlock does not support current versions of conda. As of 2026, they recommend that if you insist upon using conda as a package manager, that you install an older version of miniconda as current versions will not work.

Once a working version of **conda** is installed, creating a new **conda** environment is simple:  
  
`conda create -n name_of_environment`

To activate this environment, simply run:  
  
`conda activate name_of_environment`

Once you have activated an environment, you will be able to install packages associated with this environment. An example is provided below, but I recommend googling “conda + the name of the package” to get the precise syntax needed for your package of interest.  
  
`conda install -c bioconda package_name`

You should now have access to the software of your choosing whenever you activate this conda environment! **IMPORTANT NOTE:** While conda activate name_of_environment will work for interactive jobs, a different syntax is needed to activate conda environments within a SLURM submission script:

`source activate /path_to_miniconda/miniconda3/envs/name_of_environment`

You can find this path to your installation of miniconda by running the command **which conda**:

```
(base) [SUNET_ID@sh02-09n06 ~]$ which conda
/path_to_miniconda/miniconda3/bin/conda
```

## Accessing frequently used modules or packages

If you find yourself frequently using packages in certain directories or frequently loading the same modules, you can add these packages to your bash profile to have them automatically load into every new sherlock instance. Your bash profile is a file located at ~/.bashrc. If you are new to sherlock, this file may not have yet been initialized. To get started creating your bash profile, begin by accessing this file with a text editor such as nano.

```
(base) [khunnicu@sh02-09n06 ~]$ nano ~/.bashrc
```

To have modules load automatically, you will add lines of this format to this file:

`module load perl`

To have packages within specific folders load automatically, you will add lines of this format to the file:

`export PATH="/home/groups/schumer/shared_bin:$PATH"`

Some package managers and individual packages will edit your bash profile for you as part of installation (e.g. **miniconda**). It is best to not modify these additions to ensure these packages continue to work. **IMPORTANT NOTE:** packages and modules loaded into your bash profile are loaded in the order in which they occur in the file. This means that if you have conflicting package versions listed in your bash profile, the one that is ultimately loaded upon login will be the last version listed in the file. Take note of this as this can result in errors in your job (see [Troubleshooting Common Errors](#Troubleshooting Common Errors) for more helpful tips).

Here is an example bash profile:

`(base) [SUNET_ID@sh02-09n06 ~]$ cat ~/.bashrc`

```
# User specific aliases and functions
export PATH="/home/groups/schumer/shared_bin/:$PATH"

module load perl
module load devel
module load gcc/10.1.0
module load R
module load boost/1.69.0
module load gsl
```

## Installation Requests 
After confirming that your software of interest is not installed and that you are unable to install it on your own, you can place a new request by emailing SRCC support at [srcc-support@stanford.edu](mailto:srcc-support@stanford.edu). SRCC manages several clusters, so when you send an email, specify:
- Your name and SUNet ID
- The cluster (Sherlock)
- Details about your problem, including error messages and your attempts at debugging
- Links to the documentation for your software (installation instructions/the homepage/etc) as well as a link to download and any information about what version you want
Make sure your requests are polite in tone and it’s always a good idea to say “Thank you” when submitting a request. 

# Common workflows
Our lab policy is to publish thoroughly tested useful scripts to our lab GitHub to encourage reproducibility. These scripts and common workflows can be found here: [https://github.com/Schumerlab/doc](https://github.com/Schumerlab/doc). As of time of writing, our lab GitHub includes the following documents to help guide your bioinformatics journey:
- Common_errors.md
- commonly_used_programs.md
- commonly_used_scripts.md
- commonly_used_workflows.md
- download_data_to_sherlock_from_admera.md
- download_data_to_sherlock_from_box_dropbox_or_googledocs.md
- download_data_to_sherlock_from_ncbi's_sra.md
- important_simulation_tools.md
- lab_github.md
- lab_group_on_sherlock.md
- ONT_sequencing_data_acquisition.md
- Submitting_slurm_jobs.md
- upload_data_from_sherlock_to_ncbi's_sra.md
- Useful Unic:perl:awk commands.md
- where_can_I_find_commonly_used_files.md

# Editing Files 
There are a few built in file editors on the cluster for editing files directly without copying/pasting from your local machine. The most straightforward editor that I’d recommend for most users is **nano**. **Vim** and emacs **are** also available, but I do not recommend them for beginners. All of these options are accessible by typing `nano/vim/emacs fileName`.

# Reproducibility 
An overview of some best practices to encourage reproducibility and transparency in our research.

## File naming and links 
In the interest of reproducibility, the Schumer lab has a few recommended guidelines for file naming. 
- Always use a consistent format for similar files. For example, an appropriate way of naming different alignment files for different individuals could be softwareUsed_sampleID_to_refGenomeName.bam. 
- When appropriate, include a date in the file name (very useful for your final analyses such as a filtered vcf file or a featureCounts file).
- **Do not name your files with spaces in the file name**. This makes it extremely difficult to interact with them on the command line. There are several work arounds to spaces such as hyphens (**the-file-name.txt**), underscores (**the_file_name.txt**), or capitalization of separate words (**theFileName.txt**). Pick a system that works best with your brain and stick to it. This is important not only for work in **unix** but also statistics software such as R which also does not allow spaces in filenames or variable names.
- **Avoid permanently renaming raw data files.** A good alternative to this is through links. You can make links in your scratch working directory that have consistent, reproducible names with `ln -s sourceFile myFile` . There are two types of links, soft links and hard links. Soft links just act as a pointer to the original file and do not take up significant storage space while hard links essentially act as a duplicate file and take up the same amount of space as the original file. Softlinks are not taken by scratch purges so are a good way to keep important file links in scratch directories.
- **Folder Structure:** When analyzing big data with multiple analysis steps, you can quickly be overwhelmed by the number of files you’ve generated.  Before beginning a project, think critically about what the best way to organize all these files will be for you. I personally prefer a recursive structure where each level of the analysis is nested within the previous level. For example, rawData -> analysis_1 -> analysis_2 -> analysis_3 -> finalSummaryFile.txt . But consider whether a parallel style (e.g. rawData -> analysis_1, rawData -> analysis_2, rawData -> analysis_3) makes more sense for your project.

## General Reproducibility Practices
- You should also create a **ReadMe.txt** file for each project that includes the necessary metadata associated with your raw data (what sequencers were they generated on, what type of data is it, metadata for the individuals, etc). Future you will appreciate it and it will help you when writing up your results. 
- **You should also keep a log** of all commands you run/jobs you submitted for the entirety of the project. You should always be able to theoretically recreate any file you generated at any step in your analysis pipeline. A standard tool in our lab for this purpose is **Benchling**.
- I do not recommend keeping code logs in **word** as autocorrect is not your friend when you are typing commands or file paths. I recommend installing a good text editor such as **TextWrangler** or **BBedit** or a markdown editor like **Obsidian** as these are very compatible with code and will not autocorrect your text.
- You should also log all commands run in downstream analysis such as in **R**. [Rmarkdown](https://rmarkdown.rstudio.com/) through Rstudio is a great tool for this. It is worth your time investment to learn how to create reproducible code documents.
- Check the Schumer Lab Dropbox for additional reproducibility guidelines.

# Schumer Lab Best Practices 
As discussed throughout, our computing resources are shared among us and other labs. Once a user requests a resource, other users are prevented from also accessing that resource. Because of the inherent nature of these limits, it is in all of our best interests to make sure that our own computing use does not prevent others from also making progress on their own projects. There are a few ways you as an individual can contribute to keeping the cluster running smoothly and operating fairly:
- Think both about how many nodes and how many cores each job uses when submitting multiple jobs. The runtime of your jobs is another important consideration. It is especially important to be conscious of your total use of the clusters if your jobs will be running for multiple days/weeks and would cause long term delays to other users. 
- If someone is using the server in a way that is hampering your own usage, just slack them and ask. **Be compassionate towards new users who may not realize some behaviors could be problematic.** The Schumer lab would like to avoid public shaming culture. 

# Troubleshooting Common Errors
So your workflow didn’t run as expected... bummer! Here are some examples of how to troubleshoot common errors (minimap2 used as an example), but similar lessons apply to other programs. These sorts of errors happen all the time (even for experienced programmers) are good questions to ask yourself before asking for help.
- **Did you check the slurm output?** The .err and .out files usually contain logs of what has gone wrong during your job run.
- **Is the program actually available to run?** Some common instances where this might not be the case is 1) if you forgot to load it as a module, 2) if you didn’t activate your conda environment, or 3) if the executable isn’t in your bin/working directory (and you need to provide the full path).  

To check you can type the program in the command line or in your script, e.g. **minimap2** or **which minimap2** to see if it runs. You might need to follow the above steps (in your script if it doesn’t run).

```
(base) [SUNET_ID@sh04-ln05 login ~]$ minimap2
Usage: minimap2 [options] <target.fa>|<target.idx> [query.fa] [...]

(base) [SUNET_ID@sh04-ln05 login ~]$ which minimap2
/home/groups/schumer/lab_member_folders/kelsie/packages/minimap2/minimap2
```

- **Are the files in your current working directory?** A very common issue is if your data do not actually exist in the directory that you’re working in. It’s often helpful to **ls** all the files that will go into your command.
- **Do your files have the correct stuff in them?** Especially when running pipelines, it’s possible that your file isn’t what you think it is. For instance, maybe a file exists but is empty (run **ls -lh** to check file size).
- **Have you read the program’s documentation?** Most programs have a help command **minimap2** -h or you can run **man minimap2** to read the full manual or you can look online for the documentation (which should exist for most programs, [https://github.com/lh3/minimap2](https://github.com/lh3/minimap2)). Incorrectly invoking a command (such as omitting required parameters) can often lead to errors that will terminate your job.
- **Have you tried to interpret the error?** Hopefully the program’s error messages are interpretable, but if not, try Google-ing the error, browsing online forums (e.g. StackOverflow), or asking ChatGPT/Claude. **IMPORTANT NOTE:** these AI tools are trained on StackOverflow etc, so you might get similar answers, although success rates vary.

Finally, be sure to join the Schumer Lab “bioinformatics_help” slack channel! This is a great resource for any computing issues you run into that weren’t addressed in this guide (as is google!) and a great place for communication about resource use between users. **Happy computing!**

By Kelsie E. Hunnicutt
Lasted updated: 30Mar26
