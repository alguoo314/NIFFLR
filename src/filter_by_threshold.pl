#!/usr/bin/env perl
my $threshold=$ARGV[0];
while($line=<STDIN>){
  chomp($line);
  push(@lines,$line);
}
my $current_count;
for($j=0;$j<=$#lines;$j++){
  @F=split(/\t/,$lines[$j]);
  if($F[2] eq "transcript"){
    @f=split(/;/,$F[8]);
    for($i=0;$i<=$#f;$i++){
      if($f[$i] =~ /^read_num=/){
        @ff=split(/=/,$f[$i]);
        $count[int($ff[1])]++;
        $current_count=$ff[1];
      }
    }
  }
  $linecount[$j]=$current_count;
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
print "min count = $min_count\n";
for($j=0;$j<=$#lines;$j++){
  print $lines[$j],"\n" if($linecount[$j]>$min_count);
}
