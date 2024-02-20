#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
INPUT_GTF="input.gtf"
INPUT_READS="reads.fasta"
REF="reference.fasta"
OUTPUT_PREFIX="output"
DISCARD_INTERM=false
QUANT=false
JF_THREADS=16
BASES=17.0
MER=15
if tty -s < /dev/fd/1 2> /dev/null; then
    GC='\e[0;32m'
    RC='\e[0;31m'
    NC='\e[0m'
fi

trap abort 1 2 15
function abort {
log "Aborted"
kill -9 0
exit 1
}

log () {
    dddd=$(date)
    echo -e "${GC}[$dddd]${NC} $@"
}

function error_exit {
    dddd=$(date)
    echo -e "${RC}[$dddd]${NC} $1" >&2
    exit "${2:-1}"
}

function usage {
    echo "Usage: nifflr.sh [options]"
    echo "Options:"
    echo "Options (default value in (), *required):"
    echo "-B, --bases double      For jf_aligner, filter base on percent of bases matching (17.0)"
    echo "-d, --discard           If supplied, all the intermediate files will be removed (False)"
    echo "-f, --fasta string      *Path to the fasta/fastq file containing the reads, file can ge gzipped"
    echo "-r, --ref path          *Path to the fasta file containing the genome sequence"
    echo "-g, --gtf path          *Path to the GTF file for the genome annotation"
    echo "-m, --mer uint32        Mer size (15)"
    echo "-p, --prefix string     Prefix of the output files (output)"
    echo "-q, --quantification    If supplied, NIFFLR will assign the reads back to the reference transcripts based on coverages (False)"
    echo "-t, --threads uint16    Number of threads (16)"
    echo "-h, --help              This message"
    echo "-v, --verbose           Verbose mode (False)"
}

while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -g|--gtf)
            export INPUT_GTF="$2"
            shift
            ;;
        -f|--fasta)
            export INPUT_READS="$2"
            shift
            ;;
	-B|--bases) 
	    export BASES="$2"
	    shift
            ;;
	-m|--mer)
	    export MER="$2"
	    shift
	    ;;
	 -r|--ref)
            export REF="$2"
            shift
            ;;
	 -p|--prefix)
            export OUTPUT_PREFIX="$2"
            shift
            ;;
        -d|--discard)
	    export DISCARD_INTERM=true;
	    shift
            ;;
	-q|--quantification)
            export QUANT=true;
            shift
            ;;
	-t|--threads)
            export JF_THREADS="$2"
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        -h|--help|-u|--usage)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option $1"
            exit 1        # unknown option
            ;;
    esac
    shift
done

if [ ! -s $MYPATH/create_exon_fasta.py ];then
error_exit "create_exon_fasta.py not found in $MYPATH. It must be in the directory as this script"
fi

if [ ! -s $MYPATH/add_read_counts.py ];then
error_exit "add_read_counts.py not found in $MYPATH. It must be in the directory as this script"
fi

if [ ! -s $MYPATH/majority_vote.py ];then
error_exit "majority_vote.py not found in $MYPATH. It must be in the directory as this script"
fi

if [ ! -s $MYPATH/find_path.py ];then
error_exit "find_path.py not found in $MYPATH. It must be in the directory as this script"
fi

if [ ! -s $MYPATH/generate_gtf.py ];then
error_exit "generate_gtf.py not found in $MYPATH. It must be in the directory as this script"
fi

if [ "$QUANT" = true ] && [ ! -s $MYPATH/quantification.py ];then
error_exit "quantification.py not found in $MYPATH but the --quantification switch is provided by the user"
fi

if [ ! -s $INPUT_GTF ];then
error_exit "The input gff file does not exist. Please supply a valid gff file."
fi

if [ ! -s $REF ];then
error_exit "The reference file does not exist. Please supply a valid reference fasta file.."
fi

if [ ! -s $INPUT_READS ];then
error_exit "The input reads file does not exist. Please supply a valid fasta file containing the reads."
fi

if [ ! -e nifflr.exons_extraction.success ];then
  log "Extracting exons from the GFF file and putting them into a fasta file" && \
  log "All exons are listed as in the positive strand" && \
  python $MYPATH/create_exon_fasta.py -r $REF -g <(gffread -F $INPUT_GTF) -o $OUTPUT_PREFIX.exons.fna  && \
  rm -f nifflr.alignment.success && \
  touch nifflr.exons_extraction.success || error_exit "exon extraction failed"
fi

if [ ! -e nifflr.alignment.success ];then
  log "Running jf_aligner to align between the reads and the reference exons, folowed by finding the best path through the exons in each read" && \
  SIZE=$(grep -v ">" $OUTPUT_PREFIX.exons.fna | awk '{sum += length} END {print sum}') && \
  chmod +x $MYPATH/majority_vote.py && \
  chmod +x $MYPATH/find_path.py && \
  zcat -f $INPUT_READS | fastqToFasta.pl |jf_aligner -t $JF_THREADS -B $BASES -m $MER -s $SIZE -q /dev/stdin -r $OUTPUT_PREFIX.exons.fna --coords /dev/stdout | \
  $MYPATH/majority_vote.py | \
  $MYPATH/find_path.py -o $OUTPUT_PREFIX.best_paths.fasta.tmp && \
  mv $OUTPUT_PREFIX.best_paths.fasta.tmp $OUTPUT_PREFIX.best_paths.fasta && \
  rm -f nifflr.gtf_generation.success && \
  touch nifflr.alignment.success || error_exit "jf_aligner or majority voting or finding the best path failed. Please see the detailed error messages."
fi

if [ ! -e nifflr.gtf_generation.success ];then
  log "Generating the gtf file which converts the paths of exons to transcripts" && \
  python $MYPATH/generate_gtf.py -i $OUTPUT_PREFIX.best_paths.fasta -g $OUTPUT_PREFIX.good_output.gtf -b  $OUTPUT_PREFIX.bad_output.gtf  && \
  rm -f nifflr.count.success  && \
  touch nifflr.gtf_generation.success || error_exit "GTF generation failed"
fi

if [ "$QUANT" = true ] && [ ! -e nifflr.quantification.success ];then
  log "Performing quantification of reference transcripts"
  sort -S 10% -k1,1 -Vs $OUTPUT_PREFIX.good_output.gtf | gffread -F > $OUTPUT_PREFIX.sorted.good_output.gff && \
  sort -S 10% -k1,1 -Vs $INPUT_GTF | awk -F '\t' '{if($3=="mRNA" || $3=="exon" || $3=="transcript") print $0}' |gffread -F > $OUTPUT_PREFIX.sorted.ref.gff && \
  awk '!/^#/ && !seen[$1]++ {print $1}' $OUTPUT_PREFIX.sorted.ref.gff > chr_names.txt && \
  python $MYPATH/quantification.py -a $OUTPUT_PREFIX.sorted.good_output.gff -r $OUTPUT_PREFIX.sorted.ref.gff -o $OUTPUT_PREFIX.ref.reads.assigned.gff -c chr_names.txt && \
  rm $OUTPUT_PREFIX.sorted.ref.gff $OUTPUT_PREFIX.sorted.good_output.gff chr_names.txt && \
  log "Reference transcripts quantified in $OUTPUT_PREFIX.ref.reads.assigned.gff" && \
  touch nifflr.quantification.success || error_exit "Reference transcripts quantification failed"    
fi

if [ ! -e nifflr.count.success ];then
  log "Performing quantification of assembled transcripts" && \
  gffcompare -STC $OUTPUT_PREFIX.good_output.gtf -o $OUTPUT_PREFIX 1>gffcmp.out 2>&1 && \
  sort -S 10% -k1,1 -Vs $OUTPUT_PREFIX.combined.gtf | gffread -F > $OUTPUT_PREFIX.sorted.combined.gff && \
  sort -S 10% -k1,1 -Vs $OUTPUT_PREFIX.good_output.gtf | gffread -F > $OUTPUT_PREFIX.sorted.good_output.gff && \
  awk '!/^#/ && !seen[$1]++ {print $1}' $OUTPUT_PREFIX.sorted.combined.gff > chr_names.txt && \
  python $MYPATH/quantification.py -a $OUTPUT_PREFIX.sorted.good_output.gff -r $OUTPUT_PREFIX.sorted.combined.gff -o $OUTPUT_PREFIX.asm.reads.assigned.gff -c chr_names.txt && \
  rm $OUTPUT_PREFIX.sorted.combined.gff $OUTPUT_PREFIX.sorted.good_output.gff chr_names.txt && \
  filter_by_threshold.pl 0.02 < output.asm.reads.assigned.gff >  $OUTPUT_PREFIX.asm.reads.assigned.filtered.gff && \
  touch nifflr.count.success || error_exit "Assembled transcripts quantification failed"
fi

if [ -e nifflr.count.success ];then
  log "Assembled transcripts and quantifications are in $OUTPUT_PREFIX.asm.reads.assigned.filtered.gff"
fi


