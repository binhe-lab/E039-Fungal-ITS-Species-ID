#!/usr/bin/env bash

# This script does not install software automatically.
# It documents the required software for this pipeline.

cat <<'EOF'
Required command-line tools:

  bash
  tar
  gzip
  curl or wget
  NCBI BLAST+

macOS/Homebrew:

  brew install blast curl

Linux/conda/mamba:

  mamba create -n its-blast -c bioconda -c conda-forge blast curl
  mamba activate its-blast

Check installation:

  makeblastdb -version
  blastn -version

EOF
