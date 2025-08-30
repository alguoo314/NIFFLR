#!/usr/bin/env perl
#pipe in output of gffread --tlf on the reference annotation
while($line=<STDIN>){
  chomp($line);
  @F=split(/\t/,$line);
  @f=split(/;/,$F[8]);
  for($i=0;$i<=$#f;$i++){
    if($f[$i]=~/^exons=/){
      @ff=split(/-/,$f[$i]);
      if($#ff>1){
        for($j=1;$j<$#ff;$j++){
          $junc{"$F[0] $F[6] $ff[$j]"}=1;
        }
      }
      last;
    }
  }
}
open(FILE,"gffread --tlf $ARGV[0] | ");
while($line=<FILE>){
  chomp($line);
  @F=split(/\t/,$line);
  @f=split(/;/,$F[8]);
  for($i=0;$i<=$#f;$i++){
    if($f[$i]=~/^exons=/){
      @ff=split(/-/,$f[$i]);
      if($#ff>1){
        for($j=1;$j<$#ff;$j++){
          if(not(defined($junc{"$F[0] $F[6] $ff[$j]"}))){#never seen this junction. try to fix +-5 bp
            ($jd,$ja)=split(/,/,$ff[$j]);
            for($k=1;$k<=10;$k++){
              if(defined($junc{"$F[0] $F[6] ".($jd+$k).",".($ja)})){
                $ff[$j]=($jd+$k).",".($ja);
                $k=11;
              }elsif(defined($junc{"$F[0] $F[6] ".($jd-$k).",".($ja)})){
                $ff[$j]=($jd-$k).",".($ja);
                $k=11;
              }
            }
            for($k=1;$k<=10;$k++){
              if(defined($junc{"$F[0] $F[6] ".($jd).",".($ja+$k)})){
                $ff[$j]=($jd).",".($ja+$k);
                $k=11;
              }elsif(defined($junc{"$F[0] $F[6] ".($jd).",".($ja-$k)})){
                $ff[$j]=($jd).",".($ja-$k);
                $k=11;
              }
            }
          }
        }
      }
    #print "DEBUG $f[$i] ",join("-",@ff),"\n";
    $f[$i]=join("-",@ff);
    last;
    }
  }
  print join("\t",@F[0..7]),"\t",join(";",@f),"\n";
}
  

