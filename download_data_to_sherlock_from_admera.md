# Download data to Sherlock from Admera

## Download data with ftp
1) Make sure your computer or the lab computer is set to stay awake and running for ~4-16 hours (depends on

amount of data to download and download speed (which you have absolutely no control of))

2) To download data to Sherlock with ftp, log on to sherlock and navigate to

`/oak/stanford/groups/schumer/data`

3) Make a new folder
```
mkdir Admera_tmp
cd Admera_tmp
```

4)Start and interactive job and wait for screen to say "resources allocated". this sometimes takes a minute. you can change the time after -t based on how long you think this data set will take to download

`srun --pty -p hns,normal,schumer -t 12:00:00 bash`

5) Load the ftp module:

`module load lftp`

6) Login using sftp to the Admera host, using the user name they gave you. For example:
`sftp -P 2222 moly.schumer@bpssftp.admerahealth.com`

7) Enter the password they gave you.

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

`chmod a+rx *`

18) Make sure that the permissions are correct (-rwxr-xr-x+) by using `ls -lh`

19) Move reads to their appropriate folders. Let lab members know where the data is and enter the paths into the
library inventory in Box.