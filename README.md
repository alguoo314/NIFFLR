# NIFFLR: Novel IsoForm Finder using Long RNASeq reads

NIFFLR is a tool that identifies and quantifies both known and novel isoforms in samples sequences using long-read RNA sequencing data. NIFFLR recovers known transcripts and assembles novel transcripts present in the data by aligning exons from a reference annotation to the long reads. 

Reference: https://f1000research.com/articles/14-608

# Installation insructions

## Installation from BioConda
To install NIFFLR from BioConda, run:
```
mamba install -c bioconda -c conda-forge nifflr
```
Then you can run NIFFLR by executing nifflr.sh

## Installation from source
### Compilation Dependencies
To successfully compile NIFFLR, make sure the following development libraries are installed on your system:

Boost development libraries (e.g., libboost-all-dev on Ubuntu/Debian)

zlib development libraries (e.g., zlib1g-dev on Ubuntu/Debian)

For example, on Ubuntu/Debian, run:

sudo apt-get update

sudo apt-get install libboost-all-dev zlib1g-dev

### Installation
To install, first download the latest distribution from the github release page https://github.com/alguoo314/NIFFLR/releases. Then untar/unzip the package nifflr-X.X.X.tgz, cd to the resulting folder and run `./install.sh`.  The installation script will configure and make all necessary packages.  The nifflr.sh executable will appear under bin/

### Only for developers
You can clone the development tree, but then there are dependencies such as swig and yaggo (http://www.swig.org/ and https://github.com/gmarcais/yaggo) that must be available on the PATH:

```
$ git clone https://github.com/alguoo314/NIFFLR
$ cd NIFFLR
$ git submodule init
$ git submodule update
$ cd jellyfish/ && git checkout develop
$ cd ../PacBio && git checkout jf_aligner
$ cd ../ufasta && git checkout master
$ cd ..
```
after that you can run 'make' to compile the package.  The binaries will appear under build/inst/bin.  To create a distribution, run 'make install'.
Note that on some systems you may encounter a build error due to lack of xlocale.h file, because it was removed in glibc 2.26.  xlocale.h is used in Perl extension modules used by jf-aligner.  To fix/work around this error, you can upgrade the Perl extensions, or create a symlink for xlocale.h to /etc/local.h or /usr/include/locale.h, e.g.:
```
ln -s /usr/include/locale.h /usr/include/xlocale.h
```

# Usage
To run NIFFLR, run the main executable script nifflr.sh with options:
```
NIFFLR version 2.0.0
Usage: nifflr.sh [options]
Options:
Options (default value in (), *required):
-B, --bases double                      minimum percentage of exon bases in matching K-mers (35.0)
-m, --mer int                           alignment K-mer size (12)
-k, --known float                       minimum (must be > than) intron junction coverage for detection of known transcripts (0.0)
-n, --novel float                       minimum (must be > than) intron junction coverage for detection of novel transcripts (2.0)
-f, --fasta string                      *fasta/fastq file containing the reads, file can ge gzipped, multiple files should be listed in single quotes e.g. 'file1.fastq file2.fastq'
-r, --ref path                          *fasta file containing the genome sequence
-g, --gtf path                          *GTF file for the genome annotation
-p, --prefix string                     prefix of the output files (output)
-t, --threads int                       number of threads (16)
-e, --allowed_exon_gap_or_overlap int   maximum allowed gap or overlap between two adjacent aligned exons for building a valid transcript (15)
-k, --keep                              if set, all the intermediate files will be kept
-h, --help                              this message
-v, --verbose                           verbose mode (False)
--version  
```

# Outputs
NIFFLR produces the following output files:

<output_prefix>.quantify.tsv -- tab-separated four column file with the following columns:
transcript_id -- if of the transcript, either from the input reference or TCONS_* for novel transcripts
read_count -- number of reads assigned to the transcript
intron chain -- intron chain of the transcript
min_junction_count -- minimum number of reads spanning an intron junction in the transcript

<output_prefix>.transcripts.gtf -- GTF file of detected reference and novel transcripts.  Novel transcripts have "nifflr" in the 2nd column of the GTF file.






