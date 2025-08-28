#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
INPUT_GTF="input.gtf"
INPUT_READS="reads.fasta"
REF="reference.fasta"
OUTPUT_PREFIX="output"
KEEP_INTERM=0
QUANT=false
JF_THREADS=16
BASES=35
MER=12
GAP_OVERLAP_ALLOWANCE=15
MAX_AVG_OVERLAP=5
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
    echo "-B, --bases double      For jf_aligner, filter base on percent of bases matching (35.0)"
    echo "-k, --keep              If set, all the intermediate files will be kept"
    echo "-f, --fasta string      *Path to the fasta/fastq file containing the reads, file can ge gzipped, multiple files should be listed in single quotes e.g. 'file1.fastq file2.fastq'"
    echo "-r, --ref path          *Path to the fasta file containing the genome sequence"
    echo "-g, --gtf path          *Path to the GTF file for the genome annotation"
    echo "-m, --mer uint32        Mer size (12)"
    echo "-p, --prefix string     Prefix of the output files (output)"
    echo "-t, --threads uint16    Number of threads (16)"
    echo "-e, --allowed_exon_gap_or_overlap uint16   Threshold for the allowed bases of gaps or overlaps between two adjacent exons in mapped reads for building a valid  transcript (15)"
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
	-e|--allowed_exon_gap_or_overlap)
	    export GAP_OVERLAP_ALLOWANCE="$2"
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
        -k|--keep)
	    export KEEP_INTERM=1;
            echo "Will keep intermediate files"
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

if [ ! -e nifflr.exons_extraction.success ];then
  log "Extracting exons from the transcripts" && \
  python $MYPATH/create_exon_fasta.py -r $REF --gtf <(gffread -T $INPUT_GTF) -o $OUTPUT_PREFIX.exons.fna  && \
  rm -f nifflr.alignment.success && \
  touch nifflr.exons_extraction.success || error_exit "exon extraction failed"
fi

if [ ! -e nifflr.alignment.success ];then
  log "Aligning reads to exons" && \
  SIZE=$(grep -v ">" $OUTPUT_PREFIX.exons.fna | awk '{sum += length} END {print sum}') && \
  zcat -f $INPUT_READS | \
  $MYPATH/fastqToFasta.pl |\
  $MYPATH/psa_aligner -t $JF_THREADS -B $BASES -m $MER --psa-min $MER -s $SIZE -q /dev/stdin -r $OUTPUT_PREFIX.exons.fna --coords /dev/stdout | \
  $MYPATH/majority_vote.py | \
  $MYPATH/find_path.py -o $OUTPUT_PREFIX.best_paths.$MER.fasta.tmp && \
  mv $OUTPUT_PREFIX.best_paths.$MER.fasta.tmp $OUTPUT_PREFIX.best_paths.fasta && \
  rm -f nifflr.gtf_generation.success && \
  touch nifflr.alignment.success || error_exit "jf_aligner or majority voting or finding the best path failed. Please see the detailed error messages."
fi

if [ ! -e nifflr.gtf_generation.success ];then
  log "Converting alignments of exons to transcripts" && \
  python $MYPATH/generate_gtf.py -i $OUTPUT_PREFIX.best_paths.fasta -g $OUTPUT_PREFIX.all.gtf && \
  perl -F'\t' -ane '{
    if($F[2] eq "transcript" && $F[8] =~ /transcript_id\s"(\S+)";\ssource_reads\s"(\S+)";\slongest_mapped_read_len\s"(\S+)";\sbest_matched_reads_avg_penality_score\s"(\S+)";\sbest_matched_reads_max_penality_score\s"(\S+)"/){
      print STDERR "$1\t$3\t$4\t$5\n";
      if($4<='$MAX_AVG_OVERLAP' && $5<='$GAP_OVERLAP_ALLOWANCE'){
        $flag=1;
      }else{
        $flag=0;
      }
    }print if($flag);
  }' $OUTPUT_PREFIX.all.gtf 1>$OUTPUT_PREFIX.all.gtf.tmp 2>$OUTPUT_PREFIX.stats.txt.tmp && \
  mv $OUTPUT_PREFIX.all.gtf.tmp $OUTPUT_PREFIX.gtf && \
  mv $OUTPUT_PREFIX.stats.txt.tmp $OUTPUT_PREFIX.stats.txt && \
  rm -f nifflr.quantification.success  && \
  touch nifflr.gtf_generation.success || error_exit "GTF generation failed"
fi

if [ ! -e nifflr.quantification.success ] && [ -e nifflr.gtf_generation.success ];then
  log "Performing filtering and quantification of assembled transcripts" && \
  #first we figure out which reference transcripts are present
  #fixing junctions
  gffread --tlf $INPUT_GTF | fix_junctions.pl $OUTPUT_PREFIX.gtf |gffread -M -T > $OUTPUT_PREFIX.fix.gtf.tmp && \
  mv $OUTPUT_PREFIX.fix.gtf.tmp $OUTPUT_PREFIX.fix.gtf && \
  gffread --tlf $INPUT_GTF | fix_junctions.pl <(perl -F'\t' -ane '{if($F[8] =~ /^gene_id "(\S+)"; transcript_id "(\S+)"; source_reads "(\S+)"; longest_mapped_read_len "(\S+)"; best_matched_reads_avg_penality_score "(\S+)"; best_matched_reads_max_penality_score "(\S+)";/){$flag=($5<2 ||$6<5) ? 1 : 0;}print if($flag);}' $OUTPUT_PREFIX.gtf) |gffread -M -T > $OUTPUT_PREFIX.fix.filter.gtf.tmp && \
  mv $OUTPUT_PREFIX.fix.filter.gtf.tmp $OUTPUT_PREFIX.fix.filter.gtf && \
  trmap -c '=c' $INPUT_GTF $OUTPUT_PREFIX.fix.gtf | quantify.pl $OUTPUT_PREFIX.gtf  > $OUTPUT_PREFIX.quantify_ref.txt.tmp && \
  mv $OUTPUT_PREFIX.quantify_ref.txt.tmp $OUTPUT_PREFIX.quantify_ref.txt && \
  perl -ane '{$h{$F[1]}=1 if($F[7] > 0 || ($F[7] ==-1 && $F[5]>3)||$F[0] eq "unique_ref");}END{open(FILE,"gffread -T '$INPUT_GTF' | ");while($line=<FILE>){chomp($line);@f=split(/\t/,$line);if($f[2] eq "transcript"){if($f[8] =~ /transcript_id "(\S+)";/){$flag=defined($h{$1}) ? 1:0;}}print $line,"\n" if($flag);}}' $OUTPUT_PREFIX.quantify_ref.txt > $OUTPUT_PREFIX.known.gtf.tmp && \
  mv $OUTPUT_PREFIX.known.gtf.tmp $OUTPUT_PREFIX.known.gtf && \
  gffcompare -STC $OUTPUT_PREFIX.fix.filter.gtf -o ${OUTPUT_PREFIX}_combine 1>gffcmp.out 2>&1 && \
  rm -f ${OUTPUT_PREFIX}_combine.{redundant.gtf,stats,tracking} &&\
  trmap -c '=c' ${OUTPUT_PREFIX}_combine.combined.gtf $OUTPUT_PREFIX.fix.gtf | quantify.pl $OUTPUT_PREFIX.gtf > $OUTPUT_PREFIX.quantify_novel.txt.tmp && \
  mv $OUTPUT_PREFIX.quantify_novel.txt.tmp $OUTPUT_PREFIX.quantify_novel.txt && \
  perl -ane '{$h{$F[1]}=1 if($F[7] > 2);}END{open(FILE,"'${OUTPUT_PREFIX}'_combine.combined.gtf");while($line=<FILE>){chomp($line);@f=split(/\t/,$line);if($f[2] eq "transcript"){if($f[8] =~ /transcript_id "(\S+)";/){$flag=defined($h{$1}) ? 1:0;}}print $line,"\n" if($flag);}}' $OUTPUT_PREFIX.quantify_novel.txt > $OUTPUT_PREFIX.novel.gtf.tmp && \
  mv $OUTPUT_PREFIX.novel.gtf.tmp $OUTPUT_PREFIX.novel.gtf && \
  gffcompare -T -r $OUTPUT_PREFIX.known.gtf $OUTPUT_PREFIX.novel.gtf -o $OUTPUT_PREFIX.combine 1>gffcmp.out 2>&1 && \
  rm -f ${OUTPUT_PREFIX}.combine.{loci,tracking} &&\
  gffread -T $OUTPUT_PREFIX.known.gtf <(gffread --nids <(perl -F'\t' -ane '{print "$1\n" if($F[8] =~ /transcript_id "(\S+)";(.+) class_code "(c|=)";/);}' $OUTPUT_PREFIX.combine.annotated.gtf) $OUTPUT_PREFIX.novel.gtf ) > $OUTPUT_PREFIX.transcripts.gtf.tmp && \
  mv $OUTPUT_PREFIX.transcripts.gtf.tmp $OUTPUT_PREFIX.transcripts.gtf && \
  trmap -c '=c' ${OUTPUT_PREFIX}.transcripts.gtf $OUTPUT_PREFIX.fix.gtf | quantify.pl $OUTPUT_PREFIX.gtf | grep -v "^unique_ref" > $OUTPUT_PREFIX.quantify_transcripts.txt.tmp && \
  mv $OUTPUT_PREFIX.quantify_transcripts.txt.tmp $OUTPUT_PREFIX.quantify_transcripts.txt && \
  touch nifflr.quantification.success || error_exit "Reference transcripts quantification failed"
fi


if [ -e nifflr.quantification.success ];then
  log "Assembled transcripts are in $OUTPUT_PREFIX.transcripts.gtf, transcript read counts are in $OUTPUT_PREFIX.quantify_transcripts.txt" && \
  if [ $KEEP_INTERM -lt 1 ];then
    rm -f $OUTPUT_PREFIX.fix.gtf $OUTPUT_PREFIX.fix.filter.gtf $OUTPUT_PREFIX.quantify_ref.txt $OUTPUT_PREFIX.known.gtf $OUTPUT_PREFIX.quantify_novel.txt $OUTPUT_PREFIX.combine.annotated.gtf 
  fi
fi


