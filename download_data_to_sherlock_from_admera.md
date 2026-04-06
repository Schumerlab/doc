# Download data to Sherlock from Admera

## Download data with ftp
1) Make sure your computer or the lab computer is set to stay awake and running for ~24-48 hours (depends on amount of data to download and download speed (which you have absolutely no control of))

2) To download data to Sherlock with ftp, log on to sherlock and navigate to

`/oak/stanford/groups/schumer/data`

3) Make a new folder
```
mkdir Admera_tmp_"quote #" Ex. Admera_tmp_72
cd Admera_tmp_"quote #"
```

4) Start an interactive job and wait for screen to say "resources allocated". this sometimes takes a minute. you can change the time after -t based on how long you think this data set will take to download

`srun --pty -p hns,schumer -t 48:00:00 bash`

5) Load the ftp module:

`module load lftp`

6) Login using sftp to the Admera host, using the user name they gave you. For example:
`sftp -P 2225 moly.schumer@38.122.175.98`

7) Enter the password they gave you. It is included in the email about the data transfer.

8) Change directory into the project number they provided, for example:

`cd 19176-01`

9) Make sure you have all the files you expect (should be number of samples times two)

`ls *gz | wc -l`

10) Retrieve all the files you want. For example:

`mget *.fastq.gz`

11) Make sure there is no error message saying the job was killed! If the job was killed delete the last file that was downloaded and repeat steps 4-10. If there is no error message exit the lftp server then exit the interactive job.

```
exit
exit
```

12) confirm that you have the expected number of files again just to be sure

`ls *gz | wc -l`

13) Use the sample key provided by Admera to make a file called sample.key that contains a list of admera IDs and sample IDs separated by only at "@"

Example:

```
19176-01-01@sample1
19176-01-02@sample2
```

If the file is not named sample.key the following step will not work.

14) make sure there are no spaces or "/"s in the sample names! If there are "/"'s remove them manually! you can run
this command to remove spaces.

`sed 's/ //g' sample.key > sample.key`

15) run the rename script to rename files. This may take a few minutes.

`bash /home/groups/schumer/shared_bin/Lab_shared_scripts/rename_admera_reads.sh`

16) Make sure there are no error messages and that the numbers the script prints out are the total number of files
you expect.

17) Change permissions to the files you downloaded so that only you (the owner) have write permission and such
that all other users have read+execute permissions.

`chmod 755 *`

18) Make sure that the permissions are correct (-rwxr-xr-x+) by using `ls -lh`

19) Make a file called library_names using the full names of the samples.

Example:

```
Xnez-PTES-V-23-J260.R1.fastq.gz
Xnez-PTES-V-23-J260.R2.fastq.gz
Xnez-PTES-V-23-J266.R1.fastq.gz
Xnez-PTES-V-23-J266.R2.fastq.gz
```

20) Check Negative read counts and note any extremes.

`perl /home/groups/schumer/shared_bin/Lab_shared_scripts/count_negative_reads_fastq_list_PE.pl library_names`

21) Move reads to their appropriate folders. You can use Globus to do so. Let lab members know where the data is and enter the paths into the library inventory in DropBox as well as the Lab_shared_data_map google sheet.

## Download data with FileZilla (External hard-drive needed)
- You need an external hard-drive to perform these downloads. We are unable to connect FileZilla to OAK directly, so you need to download the files locally.

1) Download and install FTP software. Visit [FileZilla's official website](https://filezilla-project.org/download.php?type=client). Select the appropriate version, and install it.

2) Enter 38.122.175.98 for the ‘Host’ value.

3) Click on "Quickconnect” • For first-time users, check "always trust the host and add the key to the cache" when prompted.

4) The project directory will appear in the lower-right area. Double-click to open and view files.

5) In the software, the left side is your local computer, and the right side is Admera Health's FTP server. Select the local disk on the left, then right-click on the FTP server on the right to perform operations like "download."

6) To expedite downloads, consider using simultaneous transfer options. For instance, click ‘Edit’ in the top left corner. From the dropdown menu click ‘settings’. Under ‘Select Page’ click ‘Transfers’. In the left side menu, adjust ‘Maximum simultaneous transfers’ and set the value to 10.

7) Use Globus to transfer files from local location to OAK into the appropriate folder.
