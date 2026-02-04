Download database
After download these database, place the files in database directory and make index file.
```
Please download from here.
https://doi.org/10.5281/zenodo.18465132

Create plasmid-enriched-pfam database index
diamond makedb --in plasmid-enrich-Pfam.full_length_sequences.fasta.gz --db plasmid-enrich-Pfam.full_length_sequences.dmnd


Create plasmid database
minimap2 -d refseq_plasmid_release202.plasmid71.mmi refseq_plasmid_release202.plasmid71.fna.gz


Download checkv database
checkv download_database database/

```

