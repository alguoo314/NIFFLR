#!/usr/bin/env perl
#this script perfroms inexact matching of intron chains
#assumes input of "trmap -c '=cj' ref.gtf output.novel.gtf" on STDIN
my $slack=20;
my %matches=();
my $qry_transcript;
my @intron_chain=();
my $exons=();
my %unique_refs=();

while($line=<STDIN>){
  chomp($line);
  if($line=~/^>/){
    process_inexact_matches(@match_lines) if(scalar(@match_lines)>0);
    @f=split(/\s/,$line);
    $qry_transcript=substr($f[0],1);
    @intron_chain=split(/,/,$f[-1]);
    @match_lines=();
  }else{
    push(@match_lines,$line);
  }
}
process_inexact_matches(@match_lines) if(scalar(@match_lines)>0);

foreach $t(keys %matches){
  print "$t\n";
}

sub process_inexact_matches{
  my @matches=@_;
  my $num_matches=0;
  my %matched_ref=();
  my @cff=split(/-/,$intron_chain[0]);
  my $ic_match=scalar(@intron_chain);
  #print "DEBUG processing matches for $qry_transcript \n",join(",",@intron_chain),"\n",join("\n",@matches),"\n";
  foreach $m (@matches){
    @f=split(/\t/,$m);
    $ref_name=$f[5];
    if($f[0]=~/(c|=)/){
      $matches{$qry_transcript}=1;
      last;
    }
    @ref_intron_chain=split(/,/,$f[6]);
    my $match=0;
    for($i=0;$i<=$#ref_intron_chain-$#intron_chain;$i++){
      @cr=split(/-/,$ref_intron_chain[$i]);
      #find the first match
      if($cff[0]>=$cr[0] && ($cr[1]-$cff[1]<$slack && $cr[1]>=$cff[1])){#found match
	$match++;
	for($j=$i+1;$j<=$#ref_intron_chain;$j++){
          @cr=split(/-/,$ref_intron_chain[$j]);
          @cf=split(/-/,$intron_chain[$j-$i]);
          if(($cf[0]-$cr[0]<$slack && $cf[0]>=$cr[0]) && ($cr[1]-$cf[1]<$slack && $cr[1]>=$cf[1])){
            $match++;
          }else{
            last;
          }
        }
        $i=$#ref_intron_chain-$#intron_chain;
      }
    }
    if($match == $ic_match){
      $matches{$qry_transcript}=1;
      last;
    }
  }
}


