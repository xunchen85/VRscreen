#!usr/bin/perl
use strict;
#
## Author:       Xun Chen
## Email:        xunchen85@gmail.com
## Date:         09/10/2020
#
#
my $directory=$ARGV[1];
open S3,">$ARGV[0].fa";
open F1,"$ARGV[0]_unmapped_1.fq"||di("di too");
open F2,"$ARGV[0]_unmapped_2.fq"||di("di too");
my @line=split;
my $read_name;

while($read_name=<F1>){
   ######### forward
   my $read_1=$read_name;
   my @temp11=split /\s+/, $read_1;
   $read_name=$temp11[0];
   $read_name=~s/\/1//;
   $read_name=~s/^@//;
   my $temp01_1=<F1>;
   my $temp02_1=<F1>;
   my $temp03_1=<F1>;
   my $temp04_1=$read_1.$temp01_1.$temp02_1.$temp03_1;

   ######### reverse
   my $read_2=<F2>;
   my $temp01_2=<F2>;
   my $temp02_2=<F2>;
   my $temp03_2=<F2>;
   my $temp04_2=$read_2.$temp01_2.$temp02_2.$temp03_2;
   $read_1=~s/^@/>/;
   $read_2=~s/^@/>/;
   print S3 "$read_1$temp01_1$read_2$temp01_2";
}
