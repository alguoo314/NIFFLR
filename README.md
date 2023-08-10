# NIFFLR: Novel IsoForm Finder using Long RNASeq reads

# Installation insructions

To install, first download the latest distribution from the github release page https://github.com/alekseyzimin/jf_aligner/releases.Then untar/unzip the package nifflr-X.X.X.tgz, cd to the resulting folder and run `./install.sh`.  The installation script will configure and make all necessary packages.  The nifflr.sh executable will appear under bin/

Only for developers:  you can clone the development tree, but then there are dependencies such as swig and yaggo (http://www.swig.org/ and https://github.com/gmarcais/yaggo) that must be available on the PATH:

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

# Usage:
```
Usage: nifflr.sh [options]
Options:
Options (default value in (), *required):
-B, --bases double      For jf_aligner, filter base on percent of bases matching (17.0)
-d, --discard           If supplied, all the intermediate files will be removed (False)
-f, --fasta string      *Path to the fasta/fastq file containing the reads, file can ge gzipped
-r, --ref path          *Path to the fasta file containing the genome sequence
-g, --gff path          *Path to the GFF file for the genome annotation
-m, --mer uint32        Mer size (15)
-p, --prefix string     Prefix of the output gtf files (output)
-q, --quantification    If supplied, niffler will assign the reads back to the reference transcripts based on coverages (False)
-t, --threads uint16    Number of threads (16)
-h, --help              This message
-v, --verbose           Verbose mode (False)
```
