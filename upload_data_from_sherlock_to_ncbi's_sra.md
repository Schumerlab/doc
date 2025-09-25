# Upload data from Sherlock to NCBI's SRA

In general, all raw sequencing data should be uploaded to the NCBI Sequence Read Archive (SRA) by the time of
publication. Note that you only have to submit sequencing reads that have not been uploaded as part of previous
projects, so make sure to coordinate with other lab members and search the SRA database to understand which
sequences have or have not already been added.

These are the steps provided by the NCBI as of July 2021, with specific details edited to match our workflow on the
Sherlock cluster. Most information here is also available at https://submit.ncbi.nlm.nih.gov/subs/sra/ (https://web.archive.org/web/20230602213258/https://submit.ncbi.nlm.nih.gov/subs/sra/)

1. Log in to the SRA Submission Portal Wizard.
2. Create new SRA submission (click on the button New submission).
3. Register your project (Bioproject) and biological samples (Biosamples) if you did not register them before at
BioProject and BioSample databases, respectively, or in other Submission Portal Wizards. Please refer to the
SRA-specific guidelines here (https://web.archive.org/web/20230602213258/https://www.ncbi.nlm.nih.gov/sra/docs/submitbio/) and to an example BioSample attributes file from the lab here (https://web.archive.org/web/20230602213258/https://submit.ncbi.nlm.nih.gov/api/2.0/files/swty8c9y/mitonuc_sra_biosample_submission.xlsx/?format=attachment).
4. Submit SRA metadata - information that will link your project, samples/experiments and file names. Please
refer to the SRA Metadata Overview (https://web.archive.org/web/20230602213258/https://www.ncbi.nlm.nih.gov/sra/docs/submitmeta/) and to an example metadata file from the lab here (https://web.archive.org/web/20230602213258/https://submit.ncbi.nlm.nih.gov/api/2.0/files/eqs2q2bz/mitonuc_sra_metadata.tsv/?format=attachment).
5. Once you have uploaded metadata, you will receive options for uploading sequence data files.
    1. Choose the File Transfer Protocol (FTP) option, at which point you'll receive credentials that you'll need for the transfer: an address, a username, and a password, as well as an account folder of the form
`uploads/<NCBI_given_folder_name>`
    2. Log on to Sherlock and navigate to the source folder where all files for submission are located
    3. Set up an FTP using rclone:
        1. Type `rclone config`
        2. Type "n" to set up new remote
        3. Input a name for the remote
        4. Type 10 to choose FTP as the type of storage to configure
        5. Input the address provided by NCBI as the FTP host to connect to
        6. Input the username provided by NCBI as the FTP username
        7. Leave FTP port blank and press enter to use default FTP port
        8. Input the password provided by NCBI as the FTP password
4. Make a folder where the files will be uploaded using a command of the following form: 

    `rclone mkdir <FTPname>:/uploads/<NCBI_given_folder_name>/<desired_upload_folder_name>`
5. Copy the files using 

    `rclone copy -L -P <address of folder with files to be
    transferred> <FTPname>:/uploads/<NCBI_given_folder_name> <desired_upload_folder_name> --include "*.fastq.gz"`

This is assuming that the reads to be uploaded are in fastq.gz format, but you can adjust the `--include` option to send
other types of data.
You will be able to see when the files have been processed and approved on the SRA Submission portal