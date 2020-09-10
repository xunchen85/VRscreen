#!/usr/bin/perl
use strict;
#
# Author:       Xun Chen
# Email:        xunchen85@gmail.com
# Date:         09/10/2020
#
my $line="";
my @line=();
my $input_sampleID=$ARGV[0];

######## ARGV[1] define the type of QC, including: h (human), v (vector), and r (repeat)

my %vector=();
### open vector
if (($ARGV[1] && $ARGV[1] =~ "v") || !($ARGV[1])) {
 open VEC,"${input_sampleID}.vector";
 while (<VEC>){
  @line=split;
  $vector{$line[0]}=0;
 }
}
my %human=();
### open human
if (($ARGV[1] && $ARGV[1] =~ "h") || !($ARGV[1])) {
 open HUM,"${input_sampleID}.human";
 while (<HUM>){
  @line=split;
  $human{$line[0]}=0;
 }
}

my %repeats=();
### open repeat2
if (($ARGV[1] && $ARGV[1] =~ "r") || !($ARGV[1])) {
 open REP,"${input_sampleID}.repeat2";
 while (<REP>){
  @line=split;
  @{$repeats{$line[0]}}=@line;
 }
}

open IN,"${input_sampleID}.reads2";
while (<IN>){
 @line=split;

### filter vector and human
 if (exists ($vector{$line[4]}) || exists ($human{$line[4]})){next;}
### Get alignment information 
 my $start1=0;
 my $end=0;
 my $length=0;
 my $start_1=0;
 my $end_1=0;
 my $len_1=0;
 my $start_2=0;
 my $end_2=0;
 my $len_2=0;
 my $len_1_f=0;
 my $len_2_f=0;
# first end
  my @temp1=($line[9]=~/(\d+)/g);
  my @temp2=($line[9]=~/([A-Z])/g);
  for(my $i=0;$i<@temp1;$i++){
  if($temp2[$i] eq "S"){   
   if($i==0){$start1=$temp1[$i]+1;$end=$temp1[$i];$length=$temp1[$i];}
   if($i==@temp1-1){$length+=$temp1[$i];}
                       }
  elsif($temp2[$i] eq "D"){next;}
  elsif($temp2[$i] eq "M" || $temp2[$i] eq "I"){
   if($i==0){$start1=1;$end=$temp1[$i];$length=$temp1[$i];}
   elsif($i==@temp1-1){$length+=$temp1[$i];$end+=$temp1[$i];}
   else{$end+=$temp1[$i];$length+=$temp1[$i]}
  }
                              }
  $start_1=$start1;$end_1=$end;$len_1=$length;
  $len_1_f=$end_1-$start_1+1;
  my @line2=();
  if ($line[3] =~ "F"){
   if (exists($repeats{$line[4]."/1"})){ @line2=@{$repeats{$line[4]."/1"}};}
  } else {
   if (exists($repeats{$line[4]."/2"})){@line2=@{$repeats{$line[4]."/2"}};}
  }
    for(my $i=1;$i<@line2;$i+=4){
     if ($line2[$i]<=$start_1 && $line2[$i+1]>=$end_1){
      $len_1_f=$len_1_f-($end_1-$start_1+1);
     } elsif ($line2[$i]>=$start_1 && $line2[$i+1]<=$end_1){
      $len_1_f=$len_1_f-($line2[$i+1]-$line2[$i]+1);
     } elsif ($line2[$i]<=$start_1 && $line2[$i+1]>=$start_1){
      $len_1_f=$len_1_f-($line2[$i+1]-$start_1+1);
     } elsif ($line2[$i]<=$end_1 && $line2[$i+1]>=$end_1){
      $len_1_f=$len_1_f-($end_1-$line2[$i]+1);
     }
   }
# second end
 if ($line[21]){
  my @temp1=($line[21]=~/(\d+)/g);
  my @temp2=($line[21]=~/([A-Z])/g);
  for(my $i=0;$i<@temp1;$i++){
  if($temp2[$i] eq "S"){     
   if($i==0){$start1=$temp1[$i]+1;$end=$temp1[$i];$length=$temp1[$i];}
   if($i==@temp1-1){$length+=$temp1[$i];}
                       }  
  elsif($temp2[$i] eq "D"){next;}
  elsif($temp2[$i] eq "M" || $temp2[$i] eq "I"){
   if($i==0){$start1=1;$end=$temp1[$i];$length=$temp1[$i];}
   elsif($i==@temp1-1){$length+=$temp1[$i];$end+=$temp1[$i];}
   else{$end+=$temp1[$i];$length+=$temp1[$i]}  }
                              }
  $start_2=$start1;$end_2=$end;$len_2=$length;
  $len_2_f=$end_2-$start_2+1;
  my @line2=();
  if ($line[15] =~ "F"){
   if (exists($repeats{$line[4]."/1"})){@line2=@{$repeats{$line[4]."/1"}};}
  } else {
   if (exists($repeats{$line[4]."/2"})){@line2=@{$repeats{$line[4]."/2"}};}
  }
    for(my $i=1;$i<@line2;$i+=4){
     if ($line2[$i]<=$start_2 && $line2[$i+1]>=$end_2){
      $len_2_f=$len_2_f-($end_2-$start_2+1);
     } elsif ($line2[$i]>=$start_2 && $line2[$i+1]<=$end_2){
      $len_2_f=$len_2_f-($line2[$i+1]-$line2[$i]+1);
     } elsif ($line2[$i]<=$start_2 && $line2[$i+1]>=$start_2){
      $len_2_f=$len_2_f-($line2[$i+1]-$start_2+1);
     } elsif ($line2[$i]<=$end_2 && $line2[$i+1]>=$end_2){
      $len_2_f=$len_2_f-($end_2-$line2[$i]+1);
     }
   } 
  if ($len_2_f + $len_1_f <50) {next;}
  else {
   $line[13]=$line[7]."_".$start_1."_".$end_1."|".$line[19]."_".$start_2."_".$end_2;
  }
 } else {
  if ($len_1_f <50) {next;}
  else {$line[13]=$line[7]."_".$start_1."_".$end_1;}
 }
 print "@line\n";
}
