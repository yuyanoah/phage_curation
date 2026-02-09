# Phage contigs curation
This pipeline removes chromosomal or plasmid contigs from phage contigs obtained via metagenomic assembly (hybridSPAdes or metaSPAdes).
Input contigs must be predicted as potential phage contigs by phage prediction tools such as Virsorter2.
This pipeline was developed in the ellie environment of the Laboratory for Symbiotic Microbiome Sciences, in RIKEN.

# Dependencies
- seqkit (2.3.0)
- checkv (0.7.0)
- prodigal (V2.6.3)
- rnammer (1.2)
- minimap2 (2.24-r1122)
- diamond (2.1.21)

# Set up in RIKEN ellie environment
Clone git
```
git clone https://github.com/yuyanoah/phage_curation.git
cd phage_curation
```
Setup database
```
Setup fetchMG database
wget http://vm-lux.embl.de/~kultima/share/mOTU/fetchMGv1.0.tar.gz
mv fetchMGv1.0.tar.gz database/
tar -xzvf database/fetchMGv1.0.tar.gz -C database/

Download checkv database
checkv download_database database/

Please download plasmid and plasmid-enriched Pfam from here into "database" directory.
https://doi.org/10.5281/zenodo.18465132

Create plasmid-enriched-pfam database index
diamond makedb --in database/plasmid-enrich-Pfam.full_length_sequences.fasta.gz --db database/plasmid-enrich-Pfam.full_length_sequences.dmnd

Create plasmid database
minimap2 -d database/refseq_plasmid_release202.plasmid71.mmi database/refseq_plasmid_release202.plasmid71.fna.gz

```

# Usage
Run the below command.
```
phage_curation.v1.0.sh input.fasta num_of_threads
```

# Workflow of this pipeline
## STEP 1.
Detect circular contigs for detecting potentially complete phage genomes.

## STEP 2.
Check the contamination rate of circular phage condidate contigs using CheckV.
The cutoff is <10% contamination.

## STEP 3.
Cehck the completeness and contamination rate of linear phage condidate contigs using CheckV.
The cutoff is >90% completeness and <10% contamination.

## STEP 4.
Filter the potential bacterial chromosomal contigs.
Contigs having >0 bacterial marker genes, which are predicted using fetchMG, are removed as chromosomal contigs.
Contigs having >0 rRNA genes, which are predicted using RNAmmer, are removed as chromosomal contigs.

## STEP 5.
Filter the potential plasmid contigs.
Contigs aligned with plasmid database with >70% identity and >10% covered fraction are removed as plasmid contigs.

# Output
Curated phage contigs: ~.CCs-LCs.checkV.rm_bacteria_plasmid.fasta

# Citation
Yuya Kiguchi†, Daiki Takewaki†, et al., Ecological and functional characteristics of gut phageome in patients with multiple sclerosis

