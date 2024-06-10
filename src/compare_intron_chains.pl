#!/usr/bin/env perl

my $ref_gtf=$ARGV[0];
my $qry_gtf=$ARGV[1];
my $chromosome;

my $intronchain="";
open(FILE,$ref_gtf);
while($line=<FILE>){
  chomp($line);
  my @F=split(/\t/,$line);
  $chromosome=$F[0];
  if($F[2] eq "transcript"){
    if($F[8] =~ /transcript_id\s"(\S+)";/){
      $intronchains_ref{$chromosome}.="$intronchain " if(length($intronchain)>0);
      $chromosome=$F[0];
      $intronchain="$1";
      $tid=$1;
    }
  }elsif($F[2] eq "exon"){
    $intronchain.=":$F[3]:$F[4]";
    $ref_len{$tid}+=$F[4]-$F[3]+1;
  }
}
$intronchains_ref{$chromosome}.=" $intronchain" if(length($intronchain)>0);

$intronchain="";
open(FILE,$qry_gtf);
while($line=<FILE>){
  chomp($line);
  my @F=split(/\t/,$line);
  if($F[2] eq "transcript"){
    if($F[8] =~ /transcript_id\s"(\S+)";.*;\sread_num\s"(\S+)";/){
      $intronchains_qry{$chromosome}.="$intronchain " if(length($intronchain)>0);
      $chromosome=$F[0];
      $intronchain = $1;
      $read_num{$1} = $2;
    }
  }elsif($F[2] eq "exon"){
    $intronchain.= ":$F[3]:$F[4]";
  }
}
$intronchains_qry{$chromosome}.=" $intronchain" if(length($intronchain)>0);


foreach $c(keys %intronchains_ref){
  my %ref_counts=();
  my %ref_assign=();
  my %assigned=();
  my %ref_trim_counts=();
  my %ref_trim_assign=();
  my @ref_ids=();
  my @ref_chains=();
  my @qry_ids=();
  my @qry_trim_chains=();
  my @qry_chains=();
  my @i_ref=split(/\s/,$intronchains_ref{$c});
  foreach $r(@i_ref){
    my ($transcript_id,$chain)=split(/:/,$r,2);
    #print "DEBUG chromosome $c transcript $transcript_id $chain\n";
    push(@ref_ids,$transcript_id);
    push(@ref_chains,$chain);
  }
  next if(not(defined($intronchains_qry{$c})));
  my @i_qry=split(/\s/,$intronchains_qry{$c});
  foreach $q(@i_qry){
    my @t=split(/:/,$q);
    push(@qry_ids,$t[0]);
    if($#t>2){
      push(@qry_trim_chains,join(":",@t[2..($#t-1)]));
    }else{
      push(@qry_trim_chains,"-");
    }
    push(@qry_chains,join(":",@t[1..$#t]));
    #print "DEBUG query chromosome $c transcript $t[0] $qry_chains[-1] $qry_trim_chains[-1]\n";
  }
  for(my $i=0;$i<=$#qry_ids;$i++){
    if($qry_trim_chains[$i] eq "-"){#single exon
    for(my $j=0;$j<=$#ref_chains;$j++){
        if(index($ref_chains[$j],$qry_chains[$i])>-1){
          $ref_counts{$ref_ids[$j]}++;
          $ref_assign{$qry_ids[$i]}.="$ref_ids[$j] ";
          $assigned{$qry_ids[$i]}=1;
        }
      }
    }else{#multi-exon, check both trim and not trim
      for(my $j=0;$j<=$#ref_chains;$j++){
        #print "DEBUG comparing $qry_ids[$i] $ref_ids[$j] $ref_chains[$j] $qry_chains[$i]\n";
        if(index($ref_chains[$j],$qry_chains[$i])>-1){
          #print "DEBUG success match\n";
          $ref_counts{$ref_ids[$j]}++;
          $ref_assign{$qry_ids[$i]}.="$ref_ids[$j] ";
          $assigned{$qry_ids[$i]}=1;
        }else{
          if(index($ref_chains[$j],$qry_trim_chains[$i])>-1){
            #print "DEBUG success trim\n";
            $ref_trim_counts{$ref_ids[$j]}++;
            $ref_trim_assign{$qry_ids[$i]}.="$ref_ids[$j] ";
            $assigned{$qry_ids[$i]}=1;
          }
        }
      }
    }
  }
  foreach $k(keys %assigned){
    my @full=();
    print "$k\n";
    if(defined($ref_assign{$k})){
      @full=split(/\s+/,$ref_assign{$k}) if(defined($ref_assign{$k}));
      my $min_len=1000000000;
      foreach $f (@full){
        $min_len=$ref_len{$f} if($ref_len{$f}<$min_len);
      }
      foreach $f (@full){
        print "$f full\n" if($ref_len{$f}==$min_len);
      }
    }else{
      my @trim=();
      @trim=split(/\s+/,$ref_trim_assign{$k}) if(defined($ref_trim_assign{$k}));
      my $min_len=1000000000;
      foreach $f (@trim){
        $min_len=$ref_len{$f} if($ref_len{$f}<$min_len);
      }
      foreach $f (@trim){
        print "$f trim\n" if($ref_len{$f}==$min_len);
      }
    }
  }
}



