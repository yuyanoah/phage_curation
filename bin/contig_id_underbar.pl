#!/usr/bin/perl

use strict;
use warnings;

while(<>){
    chomp;
    my @array = split(/_/, $_);
    my $len = @array;
    my $id = "";
    for (my $i = 0; $i < $len - 1; $i++){
        $id = $id . "_$array[$i]";
    }
    $id =~ s/^_//g;
    print "$id\n";
}
