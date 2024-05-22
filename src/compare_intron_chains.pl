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
      $intronchain="$1";
    }
  }elsif($F[2] eq "exon"){
    $intronchain.=":$F[3]:$F[4]";
  }
}
$intronchains_ref{$chromosome}.="$intronchain" if(length($intronchain)>0);

$intronchain="";
open(FILE,$qry_gtf);
while($line=<FILE>){
  chomp($line);
  my @F=split(/\t/,$line);
  $chromosome=$F[0];
  if($F[2] eq "transcript"){
    if($F[8] =~ /transcript_id\s"(\S+)";/){
      $intronchains_qry{$chromosome}.="$intronchain " if(length($intronchain)>0);
      $intronchain="$1";
    }
  }elsif($F[2] eq "exon"){
    $intronchain.=":$F[3]:$F[4]";
  }
}
$intronchains_qry{$chromosome}.="$intronchain" if(length($intronchain)>0);

foreach $c(keys %intronchains_ref){
  my %ref_counts=();
  my %ref_assign=();
  my @i_ref=split(/\s/,$intronchains_ref{$c});
  foreach $r(@i_ref){
    my ($transcript_id,$chain)=split(/:/,$r,2);
    push(@ref_ids,$transcript_id);
    push(@ref_chains,$chain);
  }
  next if(not(defined($intronchains_qry{$c})));
  my @i_qry=split(/\s/,$intronchains_qry{$c});
  for(my $i=0;$i<=$#i_qry;$i++){
    my @t=split(/:/,$i_qry[$i]);
    if($#t>2){
      $qry_trim_chain=join(":",@t[2..$#t]);
      for(my $j=0;$j<=$#ref_chains;$j++){
        if(index($ref_chains[$j],$qry_trim_chain)>-1){
          $ref_counts{$ref_ids[$j]}++;
          $ref_assign{$t[0]}.="$ref_ids[$j] ";
        }
      }
    }else{
      $qry_trim_chain=join(":",@t[1..2]);
      for(my $j=0;$j<=$#ref_chains;$j++){
        if($ref_chains[$j] eq $qry_trim_chain){
          $ref_counts{$ref_ids[$j]}++;
          $ref_assign{$t[0]}.="$ref_ids[$j] ";
        }
      }
    }
  }
  foreach $k(keys %ref_assign){
    my @f=split(/\s+/,$ref_assign{$k});
    my $max_count=0;
    foreach $t(@f){
      $max_count=$ref_counts{$t} if($ref_counts{$t}>$max_count);
    }
    print "$k\n";
    my $count=0;
    foreach $t(@f){
      $count++ if($ref_counts{$t}==$max_count);
    }
    if($count<3){
      foreach $t(@f){
        print "$t\n" if($ref_counts{$t}==$max_count); 
      }
    }
  }
}


