#!/usr/bin/env perl
#
#convert fastq on STDIN  to fasta on STDOUT, leave untouched if fasta

while($line=<STDIN>){
  if(substr($line,0,1) eq "@"){ #type is fastq
    print ">",substr($line,1); 
    $seq=""; 
    $nlines=0;
    while($line=<STDIN>){
      chomp($line); 
      last if(substr($line,0,1) eq "+");
      $seq.=$line;
      $nlines++;
    } 
    print $seq,"\n"; 
    $slines=0;
    while($line=<STDIN>){
      $slines++; 
      last if($slines==$nlines);
    }
  }elsif(substr($line,0,1) eq ">"){ #type is fasta, output unchanged
    print $line;
    while($line=<STDIN>){
      print $line;
    }
  }
}
