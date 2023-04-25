#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
INPUT_GFF="input.gff"
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
    echo "Usage: wrapper_script.sh [options]"
    echo "Options:"
    echo "Options (default value in (), *required):"
    echo "-B, --bases double      For jf_aligner, filter base on percent of bases matching (17.0)"
    echo "-d, --discard           If supplied, all the intermediate files will be removed (False)"
    echo "-f, --fasta string      *Path to the fasta file containing the reads"
    echo "-r, --ref path          *Path to the fasta file containing the reference (often refseq)"
    echo "-g, --gff path          *Path to the reference GFF file"
    echo "-m, --mer uint32        Mer size (15)"
    echo "-p, --prefix string     Prefix of the output gtf files (output)"
    echo "-q, --quantification    If supplied, niffler will assign the reads back to the reference transcripts based on coverages (False)"
    echo "-t, --threads uint16    Number of threads (16)"
    echo "-h, --help              This message"
    echo "-v, --verbose           Verbose mode (False)"
}

while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -g|--gff)
            export INPUT_GFF="$2"
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

if [ ! -s $INPUT_GFF ];then
error_exit "The input gff file does not exist. Please supply a valid gff file."
fi

if [ ! -s $REF ];then
error_exit "The reference file does not exist. Please supply a valid reference fasta file.."
fi

if [ ! -s $INPUT_READS ];then
error_exit "The input reads file does not exist. Please supply a valid fasta file containing the reads."
fi

if [ ! -e niffler.exons_extraction.success ];then
  log "Extracting exons from the GFF file and putting them into a fasta file" && \
  log "All exons are listed as in the positive strand" && \
  python $MYPATH/create_exon_fasta.py -r $REF -g $INPUT_GFF -o $OUTPUT_PREFIX.exons.fna -n $OUTPUT_PREFIX.negative_direction_exons.csv  && \
  rm -f niffler.alignment.success && \
  touch niffler.exons_extraction.success || error_exit "exon extraction failed"
fi

if [ ! -e niffler.alignment.success ];then
  log "Running jf_aligner to align between the reads and the reference exons, folowed by finding the best path through the exons in each read" && \
  SIZE=$(grep -v ">" $OUTPUT_PREFIX.exons.fna | awk '{sum += length} END {print sum}') && \
  chmod +x $MYPATH/majority_vote.py && \
  chmod +x $MYPATH/find_path.py && \
  jf_aligner -t $JF_THREADS -B $BASES -m $MER -s $SIZE -p $INPUT_READS -r $OUTPUT_PREFIX.exons.fna --coords /dev/stdout | $MYPATH/majority_vote.py -n $OUTPUT_PREFIX.negative_direction_exons.csv | $MYPATH/find_path.py -o $OUTPUT_PREFIX.best_paths.fasta && \
  rm -f niffler.gtf_generation.success && \
  touch niffler.alignment.success || error_exit "jf_aligner or majority voting or finding the best path failed. Please see the detailed error messages."
fi


if [ ! -e niffler.gtf_generation.success ];then
  log "Generating the gtf file which converts the pathes of exons as transcripts" && \
  python $MYPATH/generate_gtf.py -i $OUTPUT_PREFIX.best_paths.fasta -g $OUTPUT_PREFIX.good_output.gtf -b  $OUTPUT_PREFIX.bad_output.gtf -n $OUTPUT_PREFIX.negative_direction_exons.csv  && \
  rm -f niffler.gfftools.success  && \
  touch niffler.gtf_generation.success || error_exit "GTF generation failed"
fi

if [ "$QUANT" = true ] && [ ! -e niffler.quantification.success ];then
  log "Performing reference transcripts quantification"
  sort -k1,1 -V -s $OUTPUT_PREFIX.good_output.gtf | gffread -F > $OUTPUT_PREFIX.sorted.good_output.gff && \
  python $MYPATH/quantification.py -a $OUTPUT_PREFIX.sorted.good_output.gff -r $INPUT_GFF -o $OUTPUT_PREFIX.reads.assigned.gff && \
  touch niffler.quantification.success || error_exit "Reference transcripts quantification failed"    
fi

if [ ! -e niffler.gfftools.success ];then
  log "Running gffread -T --cluster-only and gffcompare -r  to group transcripts into loci and compare with the reference exons" && \
  gffread -T --cluster-only $OUTPUT_PREFIX.good_output.gtf &> $OUTPUT_PREFIX.after_gffread.gtf && \
  gffcompare -T -r $INPUT_GFF $OUTPUT_PREFIX.after_gffread.gtf -o $OUTPUT_PREFIX && \
  rm -f niffler.count.success  && \
  touch niffler.gfftools.success || error_exit "gffread or gffcompare failed, please check the error messages for details"
fi

if [ ! -e niffler.count.success ];then
  log "Adding the number of reads corresponded to each transcript onto the gtf anno file produced in the above step" && \
  python $MYPATH/add_read_counts.py -a $OUTPUT_PREFIX.annotated.gtf -u $OUTPUT_PREFIX.good_output.gtf -o $OUTPUT_PREFIX.reads_num_added_annotated.gtf && \
  touch niffler.gfftools.success || error_exit "Adding read counts to gtf anno files failed"
fi
