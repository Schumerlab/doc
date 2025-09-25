## Very useful resources

Here are links to some very useful resources for quick file manipulation or bioinformatic functions:

**1) FAS scriptome**

The FAS scriptome is an invaluable resource for all sorts of file manipulations:

[https://web.archive.org/web/20221019203818/http://archive.sysbio.harvard.edu/CSB/resources/computational/scriptome/UNIX/](https://web.archive.org/web/20221019203818/http://archive.sysbio.harvard.edu/CSB/resources/computational/scriptome/UNIX/)


**2) Stephen Turner's collection of one-liners**

[https://github.com/stephenturner/oneliners](https://github.com/stephenturner/oneliners)


**3) Heng Li's collection of one-liners**

[http://lh3lh3.users.sourceforge.net/biounix.shtml](http://lh3lh3.users.sourceforge.net/biounix.shtml)

## Useful general purpose one-liners

### Trim header off file

```tail -n +2 filename > newfile```

###Remove a particular column from your file

In this example, column 10:

```cut -f 10 --complement filename > newfile```

###Use awk to duplicate or modify a column

Duplicate the second column of a file:

```awk -v OFS='\t' '$2=$2"\t"$2' myfile > mynewfile```

Rewrite the first column of a file as first column _ second column:

```awk -v OFS='\t' '$1=$1"_"$2' myfile > mynewfile```

###Use awk to select rows

Select all rows where a particular column contains a certain word:

```awk -F"\t" '$3 == "transcript" { print}' myfile > mynewfile```

Select all rows where a particular column does not contain a certain word:

```awk -F"\t" '$2 !== "N" { print}' myfile > mynewfile```

###Use perl to find and replace

Find and replace an exact string match:

```perl -pi -e 's/find/replace/g' filename```

You can also find and replace with wild cards, etc:

```perl -pi -e 's/^[^_]*find/find/g' filename```

###Split or shuffle your file

Split your file into a set of files containing n lines (in this case 20):

```split -d -e -l 20 myfile_name myfile_name_split```

Randomly sample n lines from your file:

```shuf -n 10000 myfile_name > my_sampled_file```

###Count unique lines in a file

```cat myfile.txt | uniq | wc -l```

###Sort a file numerically

In this example, sort using the value in the second column:

```sort -nk 2 my_file.bed > my_sorted_file.bed```


## Useful commands for manipulating bam files

###Calculate the average coverage of a bam file

```module load biology
module load samtools
samtools depth my_bam_file | awk '{sum+=$3} END { print sum/NR}' 

```


###Calculate the number of mapped reads in a bam file

```samtools view -c -F 260 my_bam_file```