#!/bin/bash

for i in `ls *.paf`
do
# ソート(query ID順、reference ID順、start position順)
sort -k 1,1 -k 6,6 -k 3,3n $i > $i.sort

# alignmentごとのstrand, identityとquery coverage出力
perl bin/cal_id_cov_paf.pl $i.sort > $i.sort.idcov

# >=70% identity, >=10% query coverageを抽出
cat $i.sort.idcov | awk '$3 >= 0.70' | awk '$4 >= 0.10' > $i.sort.idcov.id70cov10

# pick tophit
sort -k1,1 -k4,4nr $i.sort.idcov.id70cov10 > $i.sort.idcov.id70cov10.sort
perl bin/minimap_tophit.pl $i.sort.idcov.id70cov10.sort > $i.sort.idcov.id70cov10.sort.tophit
awk '{print $1}' $i.sort.idcov.id70cov10.sort.tophit > $i.sort.idcov.id70cov10.sort.tophit.id
done

