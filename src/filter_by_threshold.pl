#!/usr/bin/env perl
my $threshold=$ARGV[0];
while($line=<STDIN>){
  chomp($line);
  push(@lines,$line);
}
my $current_count;
my $current_support;
my $minjcov;
my %maxjcov=();
for($j=0;$j<=$#lines;$j++){
  @F=split(/\t/,$lines[$j]);
  if($F[2] eq "transcript"){
    @f=split(/;/,$F[8]);
    for($i=0;$i<=$#f;$i++){
      if($f[$i] =~ /^geneID=/){
        @ff=split(/=/,$f[$i]);
        $geneID=$ff[1];
      }elsif($f[$i] =~ /^read_num=/){
        @ff=split(/=/,$f[$i]);
        $count[int($ff[1])]++;
        $current_count=$ff[1];
      }elsif($f[$i] =~ /^transcript_support=/){
        @ff=split(/=/,$f[$i]);
        $current_support=$ff[1];
      }elsif($f[$i] =~ /^full_junction_reads_coverage=/){
        @ff=split(/=/,$f[$i]);
        $current_full_cov=$ff[1];
      }
    }
  }
  $linefullcov[$j]=$current_full_cov;
  $linecount[$j]=$current_count;
  $linesupport[$j]=$current_support;
}
$thresh=0;
for($i=0;$i<$#count;$i++){
  $thresh+=$i*$count[$i];
}
$thresh*=$threshold;
$min_count=0;
$n=0;
for($i=0;$i<$#count;$i++){
  $n+=$i*$count[$i];
  if($n>$thresh && defined($count[$i])){
    $min_count=$i;
    $i=$#count;
  }
}
print "#gff\n#produced by NIFFLR\n#min read count = $min_count\n";
for($j=0;$j<=$#lines;$j++){
  #print $lines[$j],"\n" if($linecount[$j] > $min_count || $linesupport[$j] > 0.9);
  print $lines[$j],"\n" if($linecount[$j] > $min_count || $linesupport[$j] > 0.95 || $linefullcov[$j]>3);
}
