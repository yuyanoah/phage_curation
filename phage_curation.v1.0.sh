#!/bin/bash
#Date: 2021/07/29
#Developer: Yuya Kiguchi
#Usage: phage_curation.v1.0.sh input.fasta threads

fasta=$1
input=${fasta//.fasta/}
th=$2

### Version ###
version=1.0
###############

### Database ###
plasmid=database/refseq_plasmid_release202.plasmid71.mmi               ### Plasmid database ###
plasmid_MGs=database/plasmid-enrich-Pfam.full_length_sequences.dmnd   ### Plasmid-enriched genes ###
################

echo "### Phage curation pipeline for short-read assembly ###"
echo "### $version ###"
echo "### input: $fasta ###"
echo -e "### No. of threads: $th ###\n"

# Detect circular contigs (CCs)
echo "[`date +'%Y/%m/%d %H:%M:%S'`] ### Detect circular contigs ###"
cmd="perl bin/check_circularity_illumina.pl \
--force $input.fasta $input.CCs \
> $input.CCs.log 2>&1"
echo "[`date +'%Y/%m/%d %H:%M:%S'`] $cmd"
eval $cmd
# Pick CCs
cat $input.CCs/*.fa.cut.fa > $input.CCs.fasta
# Pick linear contigs (LCs)
seqkit fx2tab -nl $input.CCs.fasta | awk '{print $1}' > $input.CCs.id
find $input.CCs -name "*.fa" | xargs cat | seqkit grep -v -f $input.CCs.id > $input.LCs.fasta

# CheckV for LCs
echo "[`date +'%Y/%m/%d %H:%M:%S'`] ### CheckV for linear contigs with >90% completeness and <10% contamination ###"
cmd="checkv end_to_end -t $th $input.LCs.fasta $input.LCs.checkV -d database/checkv-db-v1.5 > $input.LCs.checkV.log 2>&1"
echo "[`date +'%Y/%m/%d %H:%M:%S'`] $cmd"
eval $cmd
cp $input.LCs.checkV/quality_summary.tsv $input.LCs.checkV.quality_summary.tsv
awk '$10!="NA"' $input.LCs.checkV.quality_summary.tsv | awk '$10>=90' | awk '$13<10' > $input.LCs.checkV.quality_summary.comp90cont10.tsv
awk '{print $1}' $input.LCs.checkV.quality_summary.comp90cont10.tsv > $input.LCs.checkV.quality_summary.comp90cont10.id
seqkit grep -f $input.LCs.checkV.quality_summary.comp90cont10.id $input.LCs.fasta > $input.LCs.checkV.fasta

# CheckV for CCs
echo "[`date +'%Y/%m/%d %H:%M:%S'`] ### CheckV for circular contigs with <10% contamination ###"
cmd="checkv end_to_end -t $th $input.CCs.fasta $input.CCs.checkV -d database/checkv-db-v1.5 > $input.CCs.checkV.log 2>&1"
echo "[`date +'%Y/%m/%d %H:%M:%S'`] $cmd"
eval $cmd
cp $input.CCs.checkV/quality_summary.tsv $input.CCs.checkV.quality_summary.tsv
awk '$13<10' $input.CCs.checkV.quality_summary.tsv > $input.CCs.checkV.quality_summary.cont10.tsv
awk '{print $1}' $input.CCs.checkV.quality_summary.cont10.tsv > $input.CCs.checkV.quality_summary.cont10.id
seqkit grep -f $input.CCs.checkV.quality_summary.cont10.id $input.CCs.fasta > $input.CCs.checkV.fasta

# Combine CCs and LCs
cmd="cat $input.LCs.checkV.fasta $input.CCs.checkV.fasta > $input.CCs-LCs.checkV.fasta"
echo "[`date +'%Y/%m/%d %H:%M:%S'`] $cmd"
eval $cmd

# Remove contigs with bacterial marker genes
# 40 marker genes
echo "[`date +'%Y/%m/%d %H:%M:%S'`] ### Remove contigs with bacterial marker genes ###"
cmd="prodigal -p meta \
-i $input.CCs-LCs.checkV.fasta \
-a $input.CCs-LCs.checkV.prodigal.faa \
-d $input.CCs-LCs.checkV.prodigal.fna \
-o $input.CCs-LCs.checkV.prodigal \
-s $input.CCs-LCs.checkV.prodigal.s \
> $input.CCs-LCs.checkV.prodigal.log 2>&1"
echo "[`date +'%Y/%m/%d %H:%M:%S'`] $cmd"
eval $cmd

cmd="database/fetchMG/fetchMG.pl -m extraction -p $input.CCs-LCs.checkV.prodigal.faa -o $input.CCs-LCs.checkV.prodigal.fmg"
echo "[`date +'%Y/%m/%d %H:%M:%S'`] $cmd"
eval $cmd
grep -v "#" $input.CCs-LCs.checkV.prodigal.fmg/$input.CCs-LCs.checkV.prodigal.all.marker_genes_scores.table | awk '{print $1}' | perl bin/contig_id_underbar.pl | sort | uniq > $input.CCs-LCs.checkV.prodigal.fmg.id

# RNAmmer
cmd="rnammer \
-S bac \
-m lsu,ssu,tsu \
-f $input.CCs-LCs.checkV.rnammer.fasta \
-gff $input.CCs-LCs.checkV.rnammer.gff \
$input.CCs-LCs.checkV.fasta \
> $input.CCs-LCs.checkV.rnammer.log 2>&1"
echo "[`date +'%Y/%m/%d %H:%M:%S'`] $cmd"
eval $cmd
grep -v "#" $input.CCs-LCs.checkV.rnammer.gff | awk '{print $1}' | sort | uniq > $input.CCs-LCs.checkV.rnammer.id

# Alignment to plasmid database
cmd="minimap2 \
-t $th \
$plasmid \
$input.CCs-LCs.checkV.fasta \
> $input.CCs-LCs.checkV.map2plasmid.paf"
echo "[`date +'%Y/%m/%d %H:%M:%S'`] $cmd"
eval $cmd

cmd="bin/minimap_parse_id70cov10.sh"
echo "[`date +'%Y/%m/%d %H:%M:%S'`] $cmd"
eval $cmd

# Detection plasmid-enriched genes
cmd="diamond \
blastp \
-p $th \
--more-sensitive \
-q $input.CCs-LCs.checkV.prodigal.faa \
-d $plasmid_MGs \
-f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen -e 1e-10 \
-o $input.CCs-LCs.checkV.prodigal.plasmid_MGs.tab \
--tmpdir . \
> $input.CCs-LCs.checkV.prodigal.plasmid_MGs.log 2>&1"
echo "[`date +'%Y/%m/%d %H:%M:%S'`] $cmd"
eval $cmd
awk '$3 >= 50' $input.CCs-LCs.checkV.prodigal.plasmid_MGs.tab | awk '(($8-$7+1)/$13) >= 0.5' > $input.CCs-LCs.checkV.prodigal.plasmid_MGs.tab.id50cov50
awk '{print $1}' $input.CCs-LCs.checkV.prodigal.plasmid_MGs.tab.id50cov50 | perl bin/contig_id_underbar.pl | sort | uniq > $input.CCs-LCs.checkV.prodigal.plasmid_MGs.tab.id50cov50.id

# Remove viral candidates with bacterial and plasmid contigs
cmd="cat \
$input.CCs-LCs.checkV.prodigal.fmg.id \
$input.CCs-LCs.checkV.rnammer.id \
$input.CCs-LCs.checkV.map2plasmid.paf.sort.idcov.id70cov10.sort.tophit.id \
$input.CCs-LCs.checkV.prodigal.plasmid_MGs.tab.id50cov50.id \
| sort | uniq \
> $input.CCs-LCs.checkV.remove.id"
echo "[`date +'%Y/%m/%d %H:%M:%S'`] $cmd"
eval $cmd
seqkit grep -v -f $input.CCs-LCs.checkV.remove.id $input.CCs-LCs.checkV.fasta > $input.CCs-LCs.checkV.rm_bacteria_plasmid.fasta
# $input.CCs-LCs.checkV.map2gtdb.paf.sort.idcov.id70cov10.sort.tophit.id

echo "########### Final output #############"
echo "Curated phage genomes : $input.CCs-LCs.checkV.rm_bacteria_plasmid.fasta"
echo "######################################"
echo "[`date +'%Y/%m/%d %H:%M:%S'`] ### Finished ###"

