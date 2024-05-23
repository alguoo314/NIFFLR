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
  }
  for(my $i=0;$i<=$#qry_ids;$i++){
    if($qry_trim_chains[$i] eq "-"){#single exon, must equal
      for(my $j=0;$j<=$#ref_chains;$j++){
        if($ref_chains[$j] eq $qry_chains[$i]){
          $ref_counts{$ref_ids[$j]}++;
          $ref_assign{$qry_ids[$i]}.="$ref_ids[$j] ";
          $assigned{$qry_ids[$i]}=1;
        }
      }
    }else{#multi-exon, check both trim and not trim
      for(my $j=0;$j<=$#ref_chains;$j++){
        if(index($ref_chains[$j],$qry_chains[$i])>-1){
          $ref_counts{$ref_ids[$j]}++;
          $ref_assign{$qry_ids[$i]}.="$ref_ids[$j] ";
          $assigned{$qry_ids[$i]}=1;
        }else{
          if(index($ref_chains[$j],$qry_trim_chains[$i])>-1){
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
    if(defined($ref_assign{$k})){
      @full=split(/\s+/,$ref_assign{$k}) if(defined($ref_assign{$k}));
      if(scalar(@full)==1){
        print "$k\n$full[0]\n";
      }
    }else{
      my @trim=();
      @trim=split(/\s+/,$ref_trim_assign{$k}) if(defined($ref_trim_assign{$k}));
      if(scalar(@trim)==1){
          print "$k\n$trim[0]\n";
      }
    }
  }
}



