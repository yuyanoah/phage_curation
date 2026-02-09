Setup database
```
Setup fetchMG database
wget http://vm-lux.embl.de/~kultima/share/mOTU/fetchMGv1.0.tar.gz
mv fetchMGv1.0.tar.gz database/
tar -xzvf database/fetchMGv1.0.tar.gz -C database/

Please download plasmid and plasmid-enriched Pfam from here into "database" directory.
https://doi.org/10.5281/zenodo.18465132

Create plasmid-enriched-pfam database index
diamond makedb --in database/plasmid-enrich-Pfam.full_length_sequences.fasta.gz --db database/plasmid-enrich-Pfam.full_length_sequences.dmnd

Create plasmid database
minimap2 -d database/refseq_plasmid_release202.plasmid71.mmi database/refseq_plasmid_release202.plasmid71.fna.gz

```

