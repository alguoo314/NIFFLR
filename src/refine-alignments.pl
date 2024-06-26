#!/usr/bin/env perl
#this code uses nucmer to refine alignments after psa_aligner
my $reads_file=$ARGV[0];
my $exons_file=$ARGV[1];

my $seq="";
my $rn="";
open(FILE,$reads_file);
while($line=<FILE>){
  chomp($line);
  if($line=~ /^>/){
    $read_seqs{$rn}=$seq if(not($rn eq ""));
    $seq="";
    @f=split(/\s+/,$line);
    $rn=substr($f[0],1);
  }else{
    $seq.=$line;
  } 
}   
$read_seqs{$rn}=$seq if(not($rn eq ""));

$seq="";
$rn="";
open(FILE,$exons_file);
while($line=<FILE>){
  chomp($line);
  if($line=~ /^>/){
    if(not($rn eq "")){
      @f=split(/-/,$rn);
      $gene=join("-",@f[0..$#f-3]);
      $orientation{$gene}=$f[-3];
      $exons_for_gene{$gene}.="$rn ";
      $exon_seqs{$rn}=$seq;
      $seq="";
    }
    @f=split(/\s+/,$line);
    $rn=substr($f[0],1);
  }else{
    $seq.=$line;
  } 
}    
if(not($rn eq "")){
  $exon_seqs{$rn}=$seq;
  @f=split(/-/,$rn);
  $gene=join("-",@f[0..$#f-3]);
  $exons_for_gene{$gene}.="$rn ";
  $orientation{$gene}=$f[-3];
}

#input on STDIN
while($line=<STDIN>){
  chomp($line);
  if($line=~ /^>/){
    @f=split(/\s+/,$line);
    $rn=substr($f[0],1);
    $read_dir{$rn}=$f[1];
  }else{
    if(not(defined($gene_for_read{$rn}))){
      @f=split(/\s+/,$line,3);
      @ff=split(/-/,$f[1]);
      $gene_name=join("-",@ff[0..$#ff-3]);
      $gene_for_read{$rn}=$gene_name;
      $reads_for_gene{$gene_name}.="$rn ";
    }
  }
}

#do the alignments, one per gene
foreach $g(keys %reads_for_gene){
  my @reads=split(/\s+/,$reads_for_gene{$g});
  my @exons=split(/\s+/,$exons_for_gene{$g});
  my $dir=$orientation{$g};
  next if(scalar(@reads)==0 || scalar(@exons)==0);
  open(FILE,">/dev/shm/$g.reads.fa");
  foreach $r(@reads){
    my $seq=$read_seqs{$r};
    #if($dir eq "R" && $read_dir{$r} eq "+" || $dir eq "F" && $read_dir{$r} eq "-"){
    if($read_dir{$r} eq "-"){
      $seq=~tr/ACGTNacgtn/TGCANtgcan/;
      $seq=reverse($seq);
    }
    print FILE ">$r $dir $read_dir{$r}\n$seq\n";
  }
  open(FILE,">/dev/shm/$g.exons.fa");
  foreach $e (@exons){
    print FILE ">$e\n$exon_seqs{$e}\n";
  }
  $cmd="nucmer --maxmatch --forward --delta /dev/stdout -l 11 -c 23 -t 16 /dev/shm/$g.exons.fa /dev/shm/$g.reads.fa |show-coords -lHq /dev/stdin |perl -ane \x27BEGIN{\$n=1;\$readname=\"\";}{push(\@readnames,\$F[-1]);\$line=\"\$F[-2] \$F[3] \$F[4] \".(\$F[3]-\$F[0]).\" \".(\$F[11]-\$F[1]+\$F[4]).\" \$F[11] \".int(\$F[6]*\$F[9]/100);push(\@lines,\$line);\$pair=\"\$F[-2] \$F[-1]\";push(\@pairs,\$pair);\$matchq=\$F[6]*\$F[9];if(not(defined(\$h{\$pair})) || \$h{\$pair}<\$matchq){\$h{\$pair}=\$matchq;\$se{\$pair}=\$line;}}END{for(\$i=0;\$i<=\$#lines;\$i++){if(\$se{\$pairs[\$i]} eq \$lines[\$i]){\$se{\$pairs[\$i]}=-1;if(not(\$readnames[\$i] eq \$readname)){print \">\$readnames[\$i]\\t+\\n\";\$readname=\$readnames[\$i];\$n=1;}print \"exon\$n \$lines[\$i]\\n\";\$n++;}}}\x27";
  #print "$cmd\n";
  system($cmd);
  unlink("/dev/shm/$g.reads.fa","/dev/shm/$g.exons.fa");
}
