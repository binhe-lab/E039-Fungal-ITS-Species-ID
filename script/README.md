# script

This folder contains command-line scripts for the ITS identification pipeline.

Scripts:

- `00_setup_notes.sh`: documents required local software.
- - `01_build_fungal_ITS_blastdb.sh`: downloads the NCBI RefSeq Targeted Loci fungal ITS FASTA and builds a local BLAST database.
- - `02_process_seq_tarball_and_blast.sh`: processes a `.tar` archive of `.Seq` files, creates a trimmed FASTA, runs BLAST, and writes output.
