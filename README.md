# Introduction to Quality Assessment

Quality control and filtering of sequencing reads is one of the most important steps in a bioinformatics analysis pipeline. However, it is not always trivial to figure out which reads needs adjustment and which can be left untouched. In this tutorial, we explain the basics of the ``FASTA` and `FASTQ` file formats, including the Phred score concept, an important quality metric used in a majority of quality control bioinformatics tools such as FASTQC. We also demonstrate how to understand and interpret these quality metrics. We then proceed to show how one can use `Trimmomatic`, a common tool used to remove read fragments when the situation is appropriate. 

An Outline is Below:

  * [FASTA Files](#fasta-files)
  * [FASTQ Files](#fastq-files)

***

## FASTA Files

FASTA format is a text-based format for representing either nucleotide sequences or peptide sequences, in which nucleotides or amino acids are represented using single-letter codes. 

A sequence in FASTA format is represented as a series of lines, each of which should be no longer than 120 characters and usually do not exceed 80 characters. Some background: This restriction was probably put into place to allow for preallocation of fixed line sizes in software, as at the time most users relied on terminals which could display 80 or 132 characters per line. Most people preferred the bigger font in 80-character modes and so it became the recommended fashion to use 80 characters or less (often 70) in FASTA lines.

The first line in a FASTA file starts either with a ">" (greater-than) symbol or, less frequently, a ";" (semicolon) and was taken as a comment. Subsequent lines starting with a semicolon would be ignored by software. Since the only comment used was the first, it quickly became used to hold a summary description of the sequence, often starting with a unique library accession number, and with time it has become commonplace use to always use ">" for the first line and to not use ";" comments (which would otherwise be ignored).

### Example

Below is an example of the Zebrafish (_Dario Rerio_) reference genome. 
__Note__: `>` (and less commonly, `;` denote comments - usually summary description of the sequence.)

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

Supported nucleaic acid codes:

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

But what about if we want to include measures of read *quality* in addition to our sequence information? For that, we'll need a new kind of format - the `FASTQ` format. 

## FASTQ Files