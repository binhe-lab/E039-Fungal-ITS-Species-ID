#!/usr/bin/env bash
# this script downloads the curated fungal ITS sequences from NCBI
# and generates a blast database. this only needs to be done once
# or when the database is updated
# written with help from ChatGPT
# 2026-05-09

set -euo pipefail

# the purpose of the following line is to get the parent directory of this script file
# so that subsequent commands can use the absolute path of the parent folder as the 
# basis 
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# directory to store the blastdb
DB_DIR="${REPO_ROOT}/data/blastdb"

# URL to download the ITS sequence file
URL="https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.ITS.fna.gz"

# path to the downloaded fasta file
FASTA_GZ="${DB_DIR}/fungi.ITS.fna.gz"
FASTA="${DB_DIR}/fungi.ITS.fna"

# blastdb prefix
DB_PREFIX="${DB_DIR}/fungi_ITS"

# create the blastdb folder if it doesn't already exist
mkdir -p "${DB_DIR}"


if [ ! -f "${FASTA_GZ}" ]; then
	echo "Downloading NCBI RefSeq Targeted Loci fungal ITS database..."
    wget -O "${FASTA_GZ}" "${URL}"
else
	echo "ITS sequence file already downloaded"
fi

echo "Decompressing FASTA..."
gzip -dc "${FASTA_GZ}" > "${FASTA}"

echo "Building BLAST database..."
echo ${FASTA}
echo ${DB_PREFIX}
makeblastdb -in ${FASTA} -dbtype nucl -parse_seqids -title "NCBI Targeted Loci Fungal ITS" -out "${DB_PREFIX}"

echo "Done."
echo "BLAST database prefix: ${DB_PREFIX}"

rm ${FASTA}
