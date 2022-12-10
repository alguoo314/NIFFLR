#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
INPUT_GFF="input.gff"
INPUT_READS="reads.fasta"
REF="reference.fasta"
OUTPUT_PREFIX="output"
DISCARD_INTERM=false
MIN_MATCH=15
MIN_CLUSTER=31
DELTA=false
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
    echo "-c,--mincluster=uint32  Sets the minimum length of a cluster of matches (31)"
    echo "-d. --discard           If supplied, all the intermediate files will be removed (False)"
    echo "-f, --fasta             *Path to the fasta file containing the reads"
    echo "-r, --ref               *Path to the fasta file containing the reference (often refseq)"
    echo "-g, --gff               *Path to the reference GFF file"
    echo "-l,--minmatch=uint32    Set the minimum length of a single exact match in nucmer (15)"
    echo "-n, --nucmer_delta      User provided nucmer file. If provided, the program will skip the nucmer process"
    echo "-p, --prefix            The prefix of the output gtf files (output)"
    echo "-h, --help              This message"
    echo "-v, --verbose           Output information (False)"
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
	-l|--minmatch)
            export MIN_MATCH="$2"
            shift
            ;;
	 -c|--mincluster)
            export MIN_CLUSTER="$2"
            shift
            ;;
	 -r|--ref)
            export REF="$2"
            shift
            ;;
	 -n|--nucmer_delta)
            export DELTA="$2"
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
python create_exon_fasta.py -r $REF -g $INPUT_GFF -o $OUTPUT_PREFIX.exons.fna -n $OUTPUT_PREFIX.negative_direction_exons.csv  && \
rm -f niffler.nucmer.success && \
touch niffler.exons_extraction.success || error_exit "exon extraction failed"
fi

if [ ! -e niffler.nucmer.success ];then
    if [[ "$DELTA" = false || ! -s $DELTA ]] ; then
	log "Nucmer delta file not provided or the path is invalid" && \
	log "Running nucmer to align between the exons and the reads" && \
	nucmer --batch 100000 -l $MIN_MATCH -c $MIN_CLUSTER -p $OUTPUT_PREFIX -t 32 $OUTPUT_PREFIX.exons.fna $INPUT_READS
    else
	log "Using existing nucmer file" && \
	cp $DELTA $OUTPUT_PREFIX.delta
    fi
    rm -f niffler.voting.success && \
    touch niffler.nucmer.success || error_exit "nucmer failed"
fi


if [ ! -e niffler.voting.success ];then
log "Perform majority voting such that for each read, only exons of the most-mapped gene to each read is kept under the read" && \
grep ">" -A 1 $OUTPUT_PREFIX.delta  > first_two_lines_only.delta && \
sed '3~3d' first_two_lines_only.delta > modified.first_two_lines_only.delta
sed 'N;s/\n/ /g' modified.first_two_lines_only.delta > one_line_per_match.txt && \
sort -k2,2 --parallel=32 --buffer-size=80% one_line_per_match.txt > $OUTPUT_PREFIX.sorted_one_line_per_match.txt
python majority_vote.py -i $OUTPUT_PREFIX.sorted_one_line_per_match.txt -o $OUTPUT_PREFIX.majority_voted.fasta -n $OUTPUT_PREFIX.negative_direction_exons.csv &&\
rm -f niffler.find_path.success && \    
touch niffler.voting.success || error_exit "Filtering by majority voting failed"    
if [ "$DISCARD_INTERM" = true ] ; then
    log "Removing intermediate files one_line_per_match.txt, $OUTPUT_PREFIX.sorted_one_line_per_match.txt, and first_two_lines_only.delta" &&\
    rm one_line_per_match.txt && \
    rm modified.first_two_lines_only.delta && \
    rm first_two_lines_only.delta
fi
fi

if [ ! -e niffler.find_path.success ];then
log "Finding best path through the exons in each read" && \
python find_path.py -i $OUTPUT_PREFIX.majority_voted.fasta -o $OUTPUT_PREFIX.best_paths.fasta  && \
rm -f niffler.gtf_generation.success  && \    
touch niffler.find_path.success || error_exit "Finding the best path failed"
fi

if [ ! -e niffler.gtf_generation.success ];then
log "Generating the gtf file which converts the pathes of exons as transcripts" && \
python generate_gtf.py -i $OUTPUT_PREFIX.best_paths.fasta -g $OUTPUT_PREFIX.good_output.gtf -b  $OUTPUT_PREFIX.bad_output.gtf -n $OUTPUT_PREFIX.negative_direction_exons.csv  && \
rm -f niffler.gfftools.success  && \ 
touch niffler.gtf_generation.success || error_exit "GTF generation failed"
fi

if [ ! -e niffler.gfftools.success ];then
log "Running gffread -T --cluster-only and gffcompare -r  to group transcripts into loci and compare with the reference exons" && \
gffread -T --cluster-only $OUTPUT_PREFIX.good_output.gtf &> $OUTPUT_PREFIX.after_gffread.gtf && \
gffcompare -r $INPUT_GFF $OUTPUT_PREFIX.after_gffread.gtf -o $OUTPUT_PREFIX && \
mv $OUTPUT_PREFIX.$OUTPUT_PREFIX.after_gffread.gtf.refmap $OUTPUT_PREFIX.after_gffread.gtf.refmap && \
mv $OUTPUT_PREFIX.$OUTPUT_PREFIX.after_gffread.gtf.tmap $OUTPUT_PREFIX.after_gffread.gtf.tmap &&\    
rm -f niffler.count.success  && \
touch niffler.gfftools.success || error_exit "gffread or gffcompare failed, please check the error messages for details"
fi

if [ ! -e niffler.count.success ];then
log "Adding the number of reads corresponded to each transcript onto the gtf anno file produced in the above step" && \
python add_read_counts.py -a $OUTPUT_PREFIX.annotated.gtf -u $OUTPUT_PREFIX.good_output.gtf -o $OUTPUT_PREFIX.reads_num_added_annotated.gtf && \
touch niffler.gfftools.success || error_exit "Adding read counts to gtf anno files failed"
fi


