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
  gffcompare -STC $OUTPUT_PREFIX.gtf -o $OUTPUT_PREFIX 1>gffcmp.out 2>&1 && \
  rm -rf $OUTPUT_PREFIX.{redundant.gtf,stats,loci,tracking} && \
  sort -S 20% -k1,1 -k4,4V -Vs $OUTPUT_PREFIX.combined.gtf | gffread -F > $OUTPUT_PREFIX.sorted.combined.gff && \
  sort -S 20% -k1,1 -k4,4V -Vs $OUTPUT_PREFIX.all.gtf | \
  tee >(awk '!/^#/ && !seen[$1]++ {print $1}' > chr_names.txt) |\
  gffread -F > $OUTPUT_PREFIX.sorted.gff && \
  python $MYPATH/count_junction_coverage.py -i $OUTPUT_PREFIX.sorted.gff -s $OUTPUT_PREFIX.exon_junction_counts.csv -c $OUTPUT_PREFIX.full_exon_junction_counts.csv && \
  python $MYPATH/quantification.py -a $OUTPUT_PREFIX.sorted.gff -r $OUTPUT_PREFIX.sorted.combined.gff -o $OUTPUT_PREFIX.asm.reads.assigned.gff -c chr_names.txt --single_junction_coverage $OUTPUT_PREFIX.exon_junction_counts.csv --full_junction_coverage $OUTPUT_PREFIX.full_exon_junction_counts.csv && \
  $MYPATH/filter_by_threshold.pl 0.0025 < $OUTPUT_PREFIX.asm.reads.assigned.gff >  $OUTPUT_PREFIX.asm.reads.assigned.prelim.gff && \
  gffcompare -T -r $INPUT_GTF $OUTPUT_PREFIX.asm.reads.assigned.prelim.gff -o combine 1>gffcmp.out 2>&1 && \
  rm -f combine.{tracking,loci,stats} && \
  perl -F'\t' -ane '{
    if($F[2] eq "transcript"){
      if($F[8] =~/transcript_id\s"(\S+)";\sgene_id\s"(\S+)";\sgene_name\s"(\S+)";.*\scmp_ref\s"(\S+)";\sclass_code\s"(\S)";/){
        if($5 eq "=" || $5 eq "c"){
          print "$1\n$4\n";
        }
      }
    }
  }' < combine.annotated.gtf > $OUTPUT_PREFIX.transcripts_identified.txt && \
  cat <(gffread -T $INPUT_GTF | perl -F'\t' -ane 'BEGIN{
    open(FILE,"'$OUTPUT_PREFIX'.transcripts_identified.txt");
    while($line=<FILE>){
      chomp($line);
      $h{$line}=1;
    }
  }{
    if($F[2] eq "transcript"){
    $flag=0;
    if($F[8] =~/transcript_id\s"(\S+)";/){
      $flag=1 if(defined($h{$1}));
    }
  }print if($flag);
  }' | tee known.gtf ) \
  <( gffread -TF $OUTPUT_PREFIX.asm.reads.assigned.prelim.gff| \
  perl -F'\t' -ane 'BEGIN{
    open(FILE,"'$OUTPUT_PREFIX'.stats.txt");
    while($line=<FILE>){
      chomp($line);
      @f=split(/\t/,$line);
      $avg_gap{$f[0]}=$f[2];
      $max_gap{$f[0]}=$f[3];
    }
    open(FILE,"'$OUTPUT_PREFIX'.transcripts_identified.txt");
    while($line=<FILE>){
      chomp($line);
      $h{$line}=1;
    }
  }{
    if($F[2] eq "transcript"){
      $flag=0;
      if($F[8] =~/transcript_id\s"(\S+)";\sgene_id\s"(\S+)";.*\soId\s"(\S+)";.*read_num\s"(\S+)";.*full_chain_reads_coverage\s"(\S+)";/){
        $flag=1 if(not(defined($h{$1})) && ($5>3 || $avg_gap{$3}<2 || $max_gap{$3}<5));
      }
    }
    if($flag){
      @ff=split(/;/,$F[8]);
      print join("\t",@F[0..7]),"\t",join(";",@ff[0..3]),"\n";
    }
  }' | tee novel.gtf ) | \
  sort -S 20% -k1,1 -k4,4V -Vs | \
  awk -F '\t' '{if($3=="mRNA" || $3=="exon" || $3=="transcript") print $0}' |gffread > $OUTPUT_PREFIX.sorted.ref.gff && \
  python $MYPATH/quantification.py -a $OUTPUT_PREFIX.sorted.gff -r $OUTPUT_PREFIX.sorted.ref.gff -o /dev/stdout -c chr_names.txt --single_junction_coverage $OUTPUT_PREFIX.exon_junction_counts.csv --full_junction_coverage $OUTPUT_PREFIX.full_exon_junction_counts.csv |
  perl -F'\t' -ane '{
    if($F[2] eq "transcript"){
      $flag=1;
      if($F[8] =~/ID=(\S+);geneID=(\S+);gene_name=(\S+);read_num=(\S+);transcript_support=(\S+);least_junction_reads_coverage=(\S+);full_chain_reads_coverage=(\d+);covered_junctions=(\d+)\/(\d+)/){
        $flag=0 if($5<0.75 && $8==0 && $9>0);
      }
    }
    print if($flag);
  }'  | \
  sort -S 20% -k1,1 -k4,4V -Vs | \
  awk -F '\t' '{if($3=="mRNA" || $3=="exon" || $3=="transcript") print $0}' |gffread > $OUTPUT_PREFIX.sorted.ref2.gff && \
  python $MYPATH/quantification.py -a $OUTPUT_PREFIX.sorted.gff -r $OUTPUT_PREFIX.sorted.ref2.gff -o /dev/stdout -c chr_names.txt --single_junction_coverage $OUTPUT_PREFIX.exon_junction_counts.csv --full_junction_coverage $OUTPUT_PREFIX.full_exon_junction_counts.csv > $OUTPUT_PREFIX.transcripts.reads.assigned.gff.tmp && \
  mv $OUTPUT_PREFIX.transcripts.reads.assigned.gff.tmp $OUTPUT_PREFIX.transcripts.reads.assigned.gff && \
  awk -F '\t' '{if($3=="transcript"){n=split($9,a,";");for(i=1;i<=n;i++){if(a[i] ~ /^ID=/){id=substr(a[i],4);}else if(a[i] ~ /^read_num=/){if(substr(a[i],10)>=1){rn=substr(a[i],10);print id"\t"rn}}}}}' $OUTPUT_PREFIX.transcripts.reads.assigned.gff  > $OUTPUT_PREFIX.transcript_read_counts.txt.tmp  && \
  mv $OUTPUT_PREFIX.transcript_read_counts.txt.tmp $OUTPUT_PREFIX.transcript_read_counts.txt && \
  touch nifflr.quantification.success || error_exit "Reference transcripts quantification failed"
fi


if [ -e nifflr.quantification.success ];then
  log "Assembled transcripts with coverage and read count information are in $OUTPUT_PREFIX.transcripts.reads.assigned.gff, transcript read counts are in $OUTPUT_PREFIX.transcript_read_counts.txt" && \
  if [ $KEEP_INTERM -lt 1 ];then
    rm -f ${OUTPUT_PREFIX}_uniq.{combined.gtf,combined.both.gtf,loci} gffcmp.out ${OUTPUT_PREFIX}.{stats.txt,combined.gtf,sorted.combined.gff,sorted.gff,exon_junction_counts.csv,asm.reads.assigned.gff,asm.reads.assigned.prelim.gff,transcripts_identified.txt,sorted.ref.gff,sorted.ref2.gff,exon_junction_counts.csv,full_exon_junction_counts.csv} novel.gtf known.gtf combine.annotated.gtf chr_names.txt scores.csv 
  fi
fi


