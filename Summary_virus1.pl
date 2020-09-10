#!usr/bin/perl
use strict;
#
# Author:       Xun Chen
# Email:        xunchen85@gmail.com
# Date:         09/10/2020
#
my $line="";
my @line=();
my %gi=();
my $candidate_virus="";
my $taxonomy=$ARGV[2];
my $virus_list=$ARGV[3];
my $len="";
my @depth=();
my $mapped_len=0;
my @temp1=();
######################################## find the size of each GI
my %virus_list=();
open VIRUS,"$taxonomy";
while(<VIRUS>){
 @line=split;
 $virus_list{$line[0]}=$line[3];
 }

######################################## find the virus name of each GI
my %virus_name=();
open VIRUS_NAME,"${virus_list}";
while(<VIRUS_NAME>){
 @line=split;
 $virus_name{$line[0]}=$line[1];
 }

######################################## Each GI;
my %Candidate_GI=();
open GI_NAME,"$ARGV[1]";
while(<GI_NAME>){
 @line=split;
 $len=$virus_list{$line[0]}-1;
 @depth=(0) x $len;
 for (my $i=8;$i<@line;$i+=8){
#  print "aa:$line[0] $line[$i+7]\n";
  if ($line[$i+4]=~"P"){
   my $temp1=$line[$i+7];
   $temp1=~s/\|/_/;
   @temp1=split /\_/, $temp1;
   for (my $j=$temp1[0];$j<$temp1[0]+($temp1[2]-$temp1[1]+1);$j++){
    if($depth[$j] eq 0){$mapped_len++;}
    $depth[$j]++;}
   for (my $j=$temp1[3];$j<$temp1[3]+($temp1[5]-$temp1[4]+1);$j++){
    if($depth[$j] eq 0){$mapped_len++;}
    $depth[$j]++;} 
  } else {
   my $temp1=$line[$i+7];
   $temp1=~s/\|/_/;
   @temp1=split /\_/, $temp1;
   for (my $j=$temp1[0];$j<$temp1[0]+($temp1[2]-$temp1[1]+1);$j++){
    if($depth[$j] eq 0){$mapped_len++;}
    $depth[$j]++;}
  }
 }
 @{$Candidate_GI{$line[0]}}=(@line[1..7],$mapped_len,join("_",@depth));
 $mapped_len=0;
 @depth=();
}

######################################## virus_v and virus_g
open F1,"$ARGV[0]";
while(<F1>){
 @line=split;
 $candidate_virus=$line[0];
 for(my $i=8;$i<@line;$i+=8){
  my @temp1=split /\|/, $line[$i+5];
  if($#temp1>1){                  ##### if one read map to multiple GIs with top AS;
   if(exists($gi{$temp1[0]})){$gi{$temp1[0]}+=1;}
   else{$gi{$temp1[0]}=1;}
   if(exists($gi{$temp1[1]})){$gi{$temp1[1]}+=1;}
   else{$gi{$temp1[1]}=1;}
               }
  else{                           ##### only one target GI;
   if(exists($gi{$temp1[0]})){$gi{$temp1[0]}+=1;}
   else{$gi{$temp1[0]}=1;}
      }
                            }
   my @temp2=sort{$gi{$b}<=>$gi{$a}} keys %gi;
   my @candidate1=();
   my @candidate2=();
   for(my $i=0;$i<@temp2;$i++){
   if(exists($virus_list{$temp2[$i]})){
    if(exists($Candidate_GI{$temp2[$i]})){
     print "$ARGV[0] $candidate_virus $temp2[$i] $gi{$temp2[$i]} $virus_list{$temp2[$i]} @line[1..7] @{$Candidate_GI{$temp2[$i]}}[0..8]\n";
     delete $Candidate_GI{$temp2[$i]};}
    else {print "$ARGV[0] $candidate_virus $temp2[$i] $gi{$temp2[$i]} $virus_list{$temp2[$i]} @line[1..7] 0 0 0 0 0 0 0\n";}
                                      }
                              }
   %gi=();
           }

####################################### virus_g left
my @temp1=keys %Candidate_GI;
for(my $i=0;$i<@temp1;$i++){
 if(exists($virus_name{$temp1[$i]}) && exists($virus_list{$temp1[$i]})){
  print"$ARGV[0] $virus_name{$temp1[$i]} $temp1[$i] 0 $virus_list{$temp1[$i]} 0 0 0 0 0 0 0 @{$Candidate_GI{$temp1[$i]}}[0..8]\n";
 }
}
