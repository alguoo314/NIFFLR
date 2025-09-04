#!/usr/bin/env perl
while($line=<STDIN>){
  chomp($line);
  if($line=~/^>/){
    my @f=split(/\s+/,$line);
    my $num_lines=substr($f[0],1);
    my @elines=();
    my %lines=();
    my %kmers=();
    my %max_offset=();
    print ">$f[1]\t";
    for(my $i=0;$i<$num_lines;$i++){
      $line=<STDIN>;
      #print $line;
      chomp($line);
      my @info=split(/\s+/,$line);
      my @exon_info=split(/\-/,$info[-1]);
      my $chrom_gene_ori=join("-",@exon_info[0..$#exon_info-2]);
      if($info[2]<$info[3]){
        $chrom_gene_ori.="+";
      }else{
        $chrom_gene_ori.="-";
      }
      if($max_offset{$chrom_gene_ori}<$info[1]){
        my $factor=$max_offset{$chrom_gene_ori}<=$info[0] ? 1 : ($info[1]-$max_offset{$chrom_gene_ori})/($info[1]-$info[0]);
        $kmers{$chrom_gene_ori}+=$info[4]*$factor;
        $max_offset{$chrom_gene_ori}=$info[1];
        if($info[2]<$info[3]){#forward match
          $exon_start= $info[0]-$info[2]+1;
          $exon_end=$info[1]+($info[10]-$info[3]);
        }else{
          $exon_start=$info[0]-($info[10]-$info[2]);
          $exon_end=$info[1]+$info[3]-1;
        }
        push(@{$lines{$chrom_gene_ori}},"$info[-1] $info[0] $info[1] $exon_start $exon_end $info[10] $info[4]");
      }
    }#done reading lines
    my $max_kmer=0;
    my $max_key="";
    foreach my $k (keys %kmers){
      if($kmers{$k}>$max_kmer){
        $max_kmer=$kmers{$k};
        $max_key=$k;
      }
    }
    my $num=1;
    my @lines_sorted=sort {(split(/\s/,$a))[1] <=> (split(/\s/,$b))[1]} @{$lines{$max_key}};
    print substr($max_key,-1),"\n";
    foreach my $l(@lines_sorted){
      print "exon$num $l\n";
      $num++;
    }
  }
}
