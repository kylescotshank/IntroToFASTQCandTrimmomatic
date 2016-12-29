# Introduction to Quality Assessment


<p align="center">
<kbd>
  <img src="Images/quality.jpg"/>
 </kbd>
 </p>

Quality control and filtering of sequencing reads is one of the most important steps in a sequence analysis pipeline. However, it is not always trivial to figure out which reads needs adjustment and which can be left untouched. In this tutorial, we explain the basics of the `FASTA` and `FASTQ` file formats, including the Phred score concept, an important quality metric used in a majority of quality control bioinformatics tools such as FASTQC. We also demonstrate how to understand and interpret these quality metrics. We then proceed to show how one can use `Trimmomatic`, a common tool used to remove read fragments when the situation is appropriate. 

  * [FASTA Files](#fasta-files)
  * [FASTQ Files](#fastq-files)
  * [Phred Scores](#phred-scores)
  * [FastQC](#fastqc)
   * [FastQC Results](#fastqc-results)
  * [Trimmomatic](#trimmomatic)

***

##FASTA Files

FASTA format is a text-based format for representing either nucleotide sequences or peptide sequences, in which nucleotides or amino acids are represented using single-letter codes. These codes are from [IUPAC](https://iupac.org). 

A sequence in FASTA format is represented as a series of lines, each of which should be no longer than 120 characters and usually do not exceed 80 characters. Some background: This restriction was probably put into place to allow for preallocation of fixed line sizes in software, as at the time most users relied on terminals which could display 80 or 132 characters per line. Most people preferred the bigger font in 80-character modes and so it became the recommended fashion to use 80 characters or less (often 70) in FASTA lines.

The first line in a FASTA file starts either with a ">" (greater-than) symbol followed by a unique sequence identifier. 

### Example

Below is an example of the Zebrafish (_Danio Rerio_) reference genome. 
__Note__: "1" is the sequence identifier for chromosome 1.

```
>1 dna:chromosome chromosome:GRCz10:1:1:58871917:1 REF
GATCTTAAACATTTATTCCCCCTGCAAACATTTTCAATCATTACATTGTCATTTCCCCTC
CAAATTAAATTTAGCCAGAGGCGCACAACATACGACCTCTAAAAAAGGTGCTGTAACATG
TACCTATATGCAGCACCACTATATGAGAGCGGCATAGCAGTGTTTAGTCACTTGGTTGCT
TTGTTTATATTAACTTGAAAGTGTGTTTTAGCTATTGAGTTTAAACAAAGGGAGCGGTTT
ACATTGAATTAAAGGCAACTACTGATGGGTTGTGTAATGTTTCAAAGAGCTGTTGCAGCA
TGAGTGGAAAATAAAACCGTATTAGTGCTGCCTGGCCCAGTTTGGCACAAAATGGAGCGA
TTCCATTAAGAGAACGATTCAGCATAAGTGGAACAGCTAAAGTTTATGAAAATTTTTAAT
CTGGATGTAGAGAATCTCATAACACAGAAACAGCACTCCTAAAGATGCATTTATACTTCT
GCATAGAGCACACAAGTATGCTTCAGCACAACCTGTGCATGGTCACATAGCCCTTGCTGT
```

### Codes/Documentation

Supported nucleic acid codes:

| Nucleic Acid Code |	Meaning	| Expanded Definition |
| ---|---|---------|---------------- |
| A	| A	| Adenine |
| C	| C	| Cytosine |
| G	| G	| Guanine |
| T	| T	| Thymine |
| U	| U	| Uracil |
| R	| A or G | Purine |
| Y	| C, T or U |	Pyrimidines |
| K	| G, T or U	| Bases which are ketones |
| M	| A or C	| Bases with aMino groups |
| S	| C or G	| Strong interaction |
| W	| A, T or U |	Weak interaction |
| B	| not A (i.e. C, G, T or U) |	 |
| D	| not C (i.e. A, G, T or U)	|  |
| H	| not G (i.e., A, C, T or U) |	|
| V	| neither T nor U (i.e. A, C or G) |	|
| N	| A C G T U	| Nucleic acid |
| -	| gap of indeterminate length	| |


***

Supported amino acid codes:

| Amino Acid Code |	Meaning |
| --- | --- |
| A | Alanine |
| B	| Aspartic acid (D) or Asparagine (N) |
| C	| Cysteine |
| D |	Aspartic acid |
| E |	Glutamic acid |
| F |	Phenylalanine |
| G	| Glycine |
| H	| Histidine |
| I	| Isoleucine |
| J	| Leucine (L) or Isoleucine (I) |
| K	| Lysine |
| L	| Leucine |
| M	| Methionine |
| N	| Asparagine |
| O	| Pyrrolysine |
| P	| Proline |
| Q	| Glutamine |
| R	| Arginine |
| S	| Serine |
| T	| Threonine |
| U	| Selenocysteine |
| V	| Valine |
| W	| Tryptophan |
| Y	| Tyrosine |
| Z	| Glutamic acid (E) or Glutamine (Q) |
| X	| any |
| *	| translation stop |
| -	| gap of indeterminate length |

###
***
But what about if we want to include measures of read *quality* in addition to our sequence information? For that, we'll need a new kind of format - the `FASTQ` format. 
***

##FASTQ Files

A FASTQ file is a FASTA file with extra lines that contain metadata related to the quality of the reads. 

A FASTQ file containing a single sequence might look like this:

```
@SEQ_ID
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
```

To further explain:

  * Line 1 begins with a '@' character and is followed by a sequence identifier and an optional description (like a FASTA description).
  * Line 2 is the raw sequence letters.
  * Line 3 begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again.
  * Line 4 encodes the quality values for the sequence in Line 2, and must contain the same number of symbols as letters in the sequence.

The character '!' represents the lowest quality while '~' is the highest. Here are the quality value characters in left-to-right increasing order of quality (ASCII):

```
 !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
```

Note that there are a variety of different *ways* that the above quality value characters are used, depending on the particular scoring metric that is applied. The most commonly seen encoding systems are Sanger “Q + 33” Shift,  Sanger “Q + 33” ASCII GLYPH,  Illumina 1.3+ “Q + 64” Shift, and  Illumina 1.3+ “Q + 64” ASCII GLYPH. You'll notice that the biggest difference between the systems is the "shift" - i.e., whether or not the scoring method is in base 33 or base 64. 

##Phred Scores

A Phred quality score is a measure of the quality of the identification of the nucleobases generated by automated DNA sequencing. Perhaps the most important use of Phred quality scores is the automatic determination of accurate, quality-based consensus sequences.

The phred quality score, or *Q* value, is defined as a property which is logarithmically related to the base-calling error probabilities P. This involves a bit more statistics than will be covered in this lesson, but the basic formula is:

Q = –10 log<sub>10</sub>(P)

A table of Q-scores (for two different encoding systems) looks like this:

<p align="center">
<kbd>
  <img src="Images/qscores.gif"/>
 </kbd>
 </p>

We now have the knowledge that we need to run (and, more importantly, understand!) FastQC.

##FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a quality control tool for high throughput sequence data created by the Babraham Institute. It's an incredibly useful open-source Java tool that is found in many sequence analysis pipelines. Though you can use an implementation of FastQC in Galaxy, we will focus on how to run it on the command line.

###Getting Started

We've taken the liberty of preparing a set of analysis scripts so that you can learn how to run common RNA-Seq analysis tools on our Linux cluster.  In order for the analysis to finish quickly, we will be using truncated data (only reads that map to chromosome 16). The automated scripts allow you to run the following workflow:

  * Examine data quality prior to trimming using FastQC
  * Trim reads using Trimmomatic
  * Examine data quality after trimming using FastQC

The scripts that you create with the following steps will *actually* allow you to do quite a bit more than this, but we will only focus on the actions above. 

####Step 1:

Login to terminal with assigned user name and password.
```
   ssh stu01@dirigo.mdibl.org
```

####Step 2:

Change directory to /nextgen3/class/Bioinfo
```
   cd /nextgen3/class/Bioinfo
```

####Step 3:

Load R (type R at command prompt)
```
   R
```

####Step 4:

Once in R, type: source(“applied_bio_airways_setup.R”)
```
   source("applied_bio_airways_setup.R")
```

####Step 5:

When prompted for username type it and press enter.

####Step 6:

Step 6.) At the confirmation prompt, type yes if it correct, or anything else to reenter the username.  Note: all files and directories will be setup at this point. 

####Step 7:

Quit out of R by typing q() and then entering n to not save the R session.

####Step 8:

Change directory to /nextgen3/class/Bioinfo/Results/<username>/Scripts. For stu01, you would type:
```
   cd /nextgen3/class/Bioinfo/Results/stu01/Scripts
```
####Step 9:

Once in the “Scripts” directory, the user can submit any of the bash scripts through ```qsub``` command.  The order for job submission is: fastqc_prior, trim, fastqc_post. Please run the jobs in this order. 

***

###Wait, what the heck is ```qsub```?

```qsub``` is the command used for job submission to a Linux cluster. It takes several command line arguments and can also use special directives found in the submission scripts or command file. 

***

###Okay, got it. Keep going.

Now, run the following command:

```
qsub fastqc_prior_trim_trunc_sample_1.txt
```

This will place several new files into your ```FASTQC_prior_trim_trunc/``` directory. The important one at the moment is [this one](http://applbio.mdibl.org/Results/stu00/FASTQC_prior_trim_trunc/sample_1_R1_fastqc.html). Open this file in your web browser. 

Some other commands of interest:

  * `qstat`: get information about the current status of the jobs
  * `qdel`: remove a job from the queue. 

##FastQC Results

<p align="center">
<kbd>
  <img src="Images/fastqc_res_1.tiff"/>
 </kbd>
 </p>

The FASTQC report is divided into several sections. For example, you can find general information about number and length of reads in the first section, called, “Basic statistics”. The list of potentially undeleted adapters will be in the section, titled, “overrepresented sequences”. Detailed explanation related to the meaning of, and potential problems with, each section of the FASTQC report can be found on this page. 

You can navigate to the corresponding parts of the report by clicking on the items in the table of contents placed on the left. Each item in the table of contents is marked by a green tick, an orange exclamation mark, or a red cross. A green tick means that this aspect of the data is of satisfactory quality, an orange exclamation mark means that something is slightly abnormal, and a red cross suggests that the corresponding parameter of the source data is very unusual and needs correction. As you can guess, the red crosses are of the greatest interest for us. To improve the quality of our sequencing data, we need to fix the parameters that are marked by red crosses.

We will leave the assessment of the data on the cluster to you, and instead will present below some obvious differences between *good* and *bad* examples of data quality.

###Good Data

<p align="center">
<kbd>
  <img src="Images/fastqc_res_good.tiff"/>
 </kbd>
 </p>

To summarize this graph: for each base position a "box-and-whisker" type plot is drawn. The elements of the plot are as follows:

  * The central red line is the median value
  * The yellow box represents the inter-quartile range (25-75%)
  * The upper and lower whiskers represent the 10% and 90% points
  * The blue line represents the mean quality

The y-axis on the graph shows the quality scores. The higher the score the better the base call. The background of the graph divides the y axis into very good quality calls (green), calls of reasonable quality (orange), and calls of poor quality (red). The quality of calls on most platforms will degrade as the run progresses, so it is common to see base calls falling into the orange area towards the end of a read.

It should be mentioned that there are number of different ways to encode a quality score in a FastQ file. FastQC attempts to automatically determine which encoding method was used, but in some very limited datasets it is possible that it will guess this incorrectly (ironically only when your data is universally very good!). The title of the graph will describe the encoding FastQC thinks your file used.

###Bad Data

<p align="center">
<kbd>
  <img src="Images/fastqc_res_bad.tiff"/>
 </kbd>
 </p>

 As you can see above - the sequence quality doesn't just degrade slightly at the end - it's bad through and through! This could be cause to reject this particular sample, once all other factors are taken into account.

 You can additionally use information from FastQC to decide whether or not one needs to "trim" the reads at a given base position. 

###Adapter Content

<p align="center">
<kbd>
  <img src="Images/fastqc_res_trim.tiff"/>
 </kbd>
 </p>

One obvious class of sequences which you might want to analyse are adapter sequences. It is useful to know if your library contains a significant amount of adapter in order to be able to assess whether you need to adapter trim or not. This module therefore does a specific search for a set of separately defined Kmers and will give you a view of the total proportion of your library which contain these Kmers. A results trace will always be generated for all of the sequences present in the adapter config file so you can see the adapter content of your library, even if it's low.

The plot itself shows a cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position. Once a sequence has been seen in a read it is counted as being present right through to the end of the read so the percentages you see will only increase as the read length goes on.

In situations like the above, we need to trim our reads - a task for which we will use `Trimmomatic`.

***

##Trimmomatic

<p align="center">
<kbd>
  <img src="Images/multiplex.jpg"/>
 </kbd>
 </p>

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is a powerful tool that performs a variety of useful trimming tasks for illumina paired-end and single ended data.

Use this line at the command line when you are ready:

```
qsub trimmomatic_trunc_sample_1.txt
```

This will trim your paired-end reads to remove adapter sequences, placing the new files into ```FASTQ_trimmed_trunc/``` directory. 

Now, if you re-run your FastQC pipeline with the following script:

```
qsub fastqc_post_trim_trunc_sample_1.txt
```

you'll be able to re-analyze your adapter content and see that it is no longer problematic. 

###A Closer Look at Trimmomatic Commands

Let's take a look at the different arguments that can be passed to ```Trimmomatic``` on the command line.

An example command may look something like this:

```{java} 
java -jar trimmomatic-0.35.jar PE -phred33 input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

This will perform the following:

  * Remove specific adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
  * Remove leading low quality or N bases (below quality 3) (LEADING:3)
  * Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
  * Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
  * Drop reads below the 36 bases long (MINLEN:36)

Here is a list of the arguments that ```trimmomatic``` can take, along with their meaning:

  * ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
  * SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
  * LEADING: Cut bases off the start of a read, if below a threshold quality
  * TRAILING: Cut bases off the end of a read, if below a threshold quality
  * CROP: Cut the read to a specified length
  * HEADCROP: Cut the specified number of bases from the start of the read
  * MINLEN: Drop the read if it is below a specified length
  * TOPHRED33: Convert quality scores to Phred-33
  * TOPHRED64: Convert quality scores to Phred-64


