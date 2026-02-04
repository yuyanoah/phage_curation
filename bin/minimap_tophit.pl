#!/usr/bin/perl
#Date: 2019/08/29
#Developer: Yuya Kiguchi
#Description: minimapの結果から算出された各queryのaligment identityとquery coverageが記載されたfileから最もidentityが高い配列をtop hitとして出力する
#	      sort -k1,1 -k3,3nrでsortしておく

use strict;
use warnings;

my $first = 0;
my $id;
while(<>){
    chomp;
    my @array = split(/\t/, $_);
    if($first == 0){
        print "$_\n";
        $id = $array[0];
        $first++;
    }
    elsif($id eq $array[0]){
        next;
    }
    else{
        print "$_\n";
        $id = $array[0];
    }
}
