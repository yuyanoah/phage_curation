#!/usr/bin/perl
#Developer: Yuya Kiguchi
#Date: 2018/02/18
#Modified : 2021/2/25 by Steve
# Modified : 2021/11/18 by Yuya
#pafファイルをsort後(sort -k 1,1 -k 6,6 -k 3,3n)のファイルをインプットとして各queryのidentityとcoverageを算出する
#alignmentが1つのreferenceに対して複数箇所にされている場合にはidentityは各alignmentの平均を算出し、query coverageは全alignmentがカバーする領域から算出される
#alignment箇所が重複している場合には重複領域の配列長を除いてquery coverageが算出される
#usage: perl cal_id_cov_paf.pl input.paf > output
#output format (tab-delimited): 1. query id 2. reference id 3. identity 4. query coverage

use strict;
use warnings;

my @array = ();         #配列
my $count = 0;          #1つのreferenceに対するalignment数
my $q_id = "dummy";     #query id
my $ref_id = "dummy";   #reference id
my $q_length = 0;       #query長
my $align_beg = 0;      #query上のalignment begin position
my $align_end = 0;      #query上のalignment end position
my $align_sum = 0;      #合計alignment長
my $id_sum = 0;         #合計identity
my $overlap_sum = 0;    #重なってalignmentされた領域長の合計

######### @array param ##########
#   $array[0]: query id
#   $array[1]: query length
#   $array[2]: query start (0-based)
#   $array[3]: query end (0-based)
#   $array[4]: strand
#   $array[5]: reference id
#   $array[6]: reference length
#   $array[7]: alignment start position at reference (0-based)
#   $array[8]: alignment end position at reference (0-based)
#   $array[9]: matched base pair
#   $array[10]: alignment block length
#   $array[11]: mapping quality
#   $array[12]: mapping quality (0-255; 255 for missing)
#################################

sub max ($$) { $_[$_[0] < $_[1]] }
sub min ($$) { $_[$_[0] > $_[1]] }

while(<>){
    chomp;
    @array = split(/\t/, $_);
    if($count == 0){                                                    #1行目の処理
        $q_id = $array[0];
        $ref_id = $array[5];
        $q_length = $array[1];
        $align_beg = $array[2];
        $align_end = $array[3];
        $align_sum = $align_sum + $array[3] - $array[2];
        $id_sum = $id_sum + $array[9]/($array[3] - $array[2]);
        $count++;
        next;
    }
    if(eof){                                                         #最終行の処理
        if($q_id ne $array[0]){                                         #queryが異なる場合
            my $id_ave = $id_sum/$count;
            my $align_length = $align_sum - $overlap_sum;
            my $q_cov = $align_length/$q_length;
            print "$q_id\t$ref_id\t$id_ave\t$q_cov\n";
            my $id_last = $array[9]/($array[3] - $array[2]);
            my $q_cov_last = ($array[3] - $array[2])/$array[1];
            print "$array[0]\t$array[5]\t$id_last\t$q_cov_last\n";
            last;
        }
        elsif($q_id eq $array[0] and $ref_id ne $array[5]){             #queryが同じかつreferenceが異なる場合
            my $id_ave = $id_sum/$count;
            my $align_length = $align_sum - $overlap_sum;
            my $q_cov = $align_length/$q_length;
            print "$q_id\t$ref_id\t$id_ave\t$q_cov\n";
            my $id_last = $array[9]/($array[3] - $array[2]);
            my $q_cov_last = ($array[3] - $array[2])/$array[1];
            print "$array[0]\t$array[5]\t$id_last\t$q_cov_last\n";
            last;
        }
        else{                                                           #queryが同じかつreferenceが同じ場合
            $overlap_sum = $overlap_sum + max(0, min($align_end, $array[3]) - max($align_beg, $array[2]));
            $q_id = $array[0];
            $ref_id = $array[5];
            $align_beg = min($align_beg, $array[2]);
            $align_end = max($align_end, $array[3]);
            $q_length = $array[1];
            $align_sum = $align_sum + $array[3] - $array[2];
            $id_sum = $id_sum + $array[9]/($array[3] - $array[2]);
            $count++;
            my $id_ave = $id_sum/$count;
            my $align_length = $align_sum - $overlap_sum;
            my $q_cov = $align_length/$q_length;
            print "$q_id\t$ref_id\t$id_ave\t$q_cov\n";
            last;
        }
    }
    elsif($q_id eq $array[0] and $ref_id ne $array[5]){                 #queryが同じかつreferenceが異なる場合
        my $id_ave = $id_sum/$count;
        my $align_length = $align_sum - $overlap_sum;
        my $q_cov = $align_length/$q_length;
        $id_sum = 0;
        $align_sum = 0;
        $overlap_sum = 0;
        print "$q_id\t$ref_id\t$id_ave\t$q_cov\n";
        $q_id = $array[0];
        $ref_id = $array[5];
        $align_end = $array[3];
        $q_length = $array[1];
        $align_sum = $align_sum + $array[3] - $array[2];
        $id_sum = $id_sum + $array[9]/($array[3] - $array[2]);
        $count = 1;
    }
    elsif($q_id eq $array[0] and $ref_id eq $array[5]){                 #queryが同じかつreferenceが同じ場合
        $overlap_sum = $overlap_sum + max(0, min($align_end, $array[3]) - max($align_beg, $array[2]));
        $ref_id = $array[5];
        $align_beg = min($align_beg, $array[2]);
        $align_end = max($align_end, $array[3]);
        $q_length = $array[1];
        $align_sum = $align_sum + $array[3] - $array[2];
        $id_sum = $id_sum + $array[9]/($array[3] - $array[2]);
        $count++;
    }
    elsif($q_id ne $array[0]){                                          #queryが異なる場合
        my $id_ave = $id_sum/$count;
        my $align_length = $align_sum - $overlap_sum;
        my $q_cov = $align_length/$q_length;
        $id_sum = 0;
        $align_sum = 0;
        $overlap_sum = 0;
        print "$q_id\t$ref_id\t$id_ave\t$q_cov\n";
        $q_id = $array[0];
        $ref_id = $array[5];
        $align_end = $array[3];
        $q_length = $array[1];
        $align_sum = $align_sum + $array[3] - $array[2];
        $id_sum = $id_sum + $array[9]/($array[3] - $array[2]);
        $count = 1;
    }
}
