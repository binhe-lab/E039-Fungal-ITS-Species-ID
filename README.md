# Yeast ITS Identification Pipeline

## Purpose

Use ITS sequencing to validate species and strain identity.

## Overview

This repository contains a command-line pipeline for validating yeast species or strain identities from Sanger ITS sequencing results.

The pipeline is designed for `.Seq` files returned from Sanger sequencing. It:

1. accepts a `.tar` archive of individual `.Seq` files;
2. renames the archive with a date prefix;
3. extracts the sequences;
4. converts them into a combined FASTA file;
5. trims the first 50 nucleotides from each sequence;
6. keeps the next 800 bp;
7. searches the sequences against a local fungal ITS BLAST database;
8. writes ranked BLAST results to the `output/` folder.

The reference database is built from the NCBI RefSeq Targeted Loci fungal ITS FASTA file:

```text
https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.ITS.fna.gz
```

Parameters such renaming rules and trimming amount can be easily adjusted in the shell script code.

## Installation

Install the following command-line tools before running the pipeline:

```text
bash
tar
gzip
curl or wget
NCBI BLAST+
```

On macOS with Homebrew:

```bash
brew install blast curl
```

On Linux using conda/mamba:

```bash
conda create -n its-blast -c bioconda -c conda-forge python>=3.10 blast curl
conda activate its-blast
```

Check that BLAST+ is available:

```bash
makeblastdb -version
blastn -version
```

## Setup

From the repository root, run:

`bash script/01_build_fungal_ITS_blastdb.sh`

This will:

1. create data/blastdb/ if needed;
1. download the NCBI curated fungal ITS FASTA file;
1. decompress it;
1. build a nucleotide BLAST database using makeblastdb.

The database files will be stored in:

`data/blastdb/`

## Running a query at the command line

Place a Sanger sequencing .tar archive in:

`data/query/`

The tarball should contain individual plain-text sequence files ending in `.Seq`.

Then run:

```bash
bash script/02_process_seq_tarball_and_blast.sh data/query/YOUR_FILE.tar
```

Example:

```bash
bash script/02_process_seq_tarball_and_blast.sh data/query/yeast_ITS_results.tar
```

The script will:

1. rename the tarball by adding the current date prefix, for example:
    `20260504-yeast_ITS_results.tar`
2. extract the `.Seq` files;
3. make a combined FASTA: `data/query/20260504-yeast_ITS_results.fasta`
4. trim each sequence by removing the first 50 nt and keeping the next 800 bp;
5. run blastn;
6. write BLAST output to:
    `output/20260504-yeast_ITS_results.blast.tsv`

## Running the web app

The repository also includes a small Flask frontend. Lab members can either
upload a `.tar`, `.tar.gz`, or `.tgz` archive to run a new BLAST search, or open
an existing result already stored in `output/`.

Install the Python dependency:

```bash
python3 -m pip install -r requirements.txt
```

(Optional) Run unit testing
```bash
python3 -m pip install pytest # if you don't already have pytest
python3 -m pytest tests/ -v
```

Start the app from the repository root:

```bash
python3 app/app.py
```

If port 5000 is already in use, choose another port:

```bash
PORT=5001 python3 app/app.py
```

Then open:

```text
http://127.0.0.1:5000/
```

Or use the alternate port you selected, for example
`http://127.0.0.1:5001/`.

The web app uses the same command-line pipeline script and writes generated
files to `data/query/` and `output/`.

The home page has two workflows:

- `Run New BLAST`: upload a sequencing tarball and process it through the
  command-line pipeline.
- `Open Existing Result`: choose a saved result from `output/` and view it
  without uploading or rerunning BLAST.

When a run finishes, the results page shows both the plain text BLAST summary
and a query browser. The query browser lets users move through one query at a
time with a dropdown or previous/next buttons, and can display aligned query
and subject sequences with subject identities shown as dots.

Older saved results that only have `output/*.summary.txt` can still be opened
in the simple output view. Results with matching `output/*.blast.tsv` files also
support the summary table, per-query browser, and dotted alignment view.

BLAST+ must be available to the process running the web app. Either install
`blastn` on your default `PATH`, or provide the executable path in the web form.
On this machine, for example:

```text
/usr/local/ncbi/blast/bin/blastn
```

The command-line script also accepts this path through `BLASTN_BIN`:

```bash
BLASTN_BIN=/usr/local/ncbi/blast/bin/blastn bash script/02_process_seq_tarball_and_blast.sh data/query/YOUR_FILE.tar
```
