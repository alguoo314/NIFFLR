#!/usr/bin/env perl
#this code reads the number of reads in  each transcript from output.gtf and then quantifies each reference transcript in the output of trmap
#assumes input of "trmap -c '=c' output.combined.gtf output.gtf" on STDIN
my $all_gtf=$ARGV[0];
my %num_reads=();
my %junctions=();
my %read_counts=();
my $qry_transcript;
my @intron_chain=();
my $exons=();

open(FILE,$all_gtf);
while($line=<FILE>){
  chomp($line);
  @f=split(/\t/,$line);
  if($f[2] eq "transcript"){
    if($f[8] =~ /^gene_id "(\S+)"; transcript_id "(\S+)"; source_reads "(\S+)"/){
      $gene=$1;
      $transcript=$2;
      @reads=split(/,/,$3);
      $num_reads{$transcript}=scalar(@reads);
    }
  }
}

my @match_lines=();
while($line=<STDIN>){
  chomp($line);
  if($line=~/^>/){
    #print "$line\n";
    process_matches(@match_lines) if(scalar(@match_lines)>0);
    @f=split(/\s/,$line);
    $qry_transcript=substr($f[0],1);
    @intron_chain=split(/,/,$f[-1]);
    @match_lines=();
  }else{
    push(@match_lines,$line);
  }
}

foreach $t(keys %read_counts){
  my $min_count=1000000000;
  foreach $c(@{$junctions{$t}}){
    $min_count=$c if($c< $min_count);
  }
  print "transcript $t exons $exons{$t} reads $read_counts{$t} min_count $min_count junction_counts ",join(" ",@{$junctions{$t}}),"\n";
}

sub process_matches{
  my @matches=@_;
  my $num_matches=0;
  my %matched_ref=();
  #print "DEBUG processing matches for $qry_transcript \n",join(",",@intron_chain),"\n",join("\n",@matches),"\n";
  foreach $m (@matches){
    @f=split(/\t/,$m);
    $ref_name=$f[5];
    @ref_intron_chain=split(/,/,$f[6]);
    $exons{$ref_name}=$f[6] if(not(defined($exons{$ref_name})));
    my $ic_match=0;
    my $match_start=-1;
    for($i=0;$i<=$#ref_intron_chain-$#intron_chain;$i++){
      #print "DEBUG $i $ref_intron_chain[$i] $intron_chain[$i]\n";
      if($ref_intron_chain[$i] eq $intron_chain[0]){
	$ic_match++;
	$match_start=$i;
	for($j=1;$j<=$#intron_chain;$j++){
	  if($ref_intron_chain[$i+$j] eq $intron_chain[$j]){
	    $ic_match++;
	  }else{
	    $match_start=-1;
	    last;
	  }
	}
      }
    }
    if($ic_match==scalar(@intron_chain)){
      #print "DEBUG found match $ref_name\n";
      $matched_ref{$ref_name}=1;
      $num_matches++;
      if(scalar(@intron_chain)>1){
        if(defined($junctions{$ref_name})){
	  for($i=$match_start;$i<$match_start+$ic_match-1;$i++){
	    $junctions{$ref_name}->[$i]++;
	  }
        }else{
	  my @arr=();
	  for($i=0;$i<$#ref_intron_chain;$i++){
	    push(@arr,0);
	  }
	  for($i=$match_start;$i<$match_start+$ic_match-1;$i++){
	    $arr[$i]++;
	  }
	  $junctions{$ref_name}=\@arr;
        }
      }
    }
  }
  if($num_matches>0){
    foreach $m (@matches){
      @f=split(/\t/,$m);
      my $ref_name=$f[5];
      if(defined($matched_ref{$ref_name})){
        #print "DEBUG updating $ref_name $exons{$ref_name} ",join(" ",@{$junctions{$ref_name}}),"\n";
        $read_counts{$ref_name}+=$num_reads{$qry_transcript}/$num_matches;
      }
    }
  }else{
    print "unmatched $qry_transcript\n";
  }
}


