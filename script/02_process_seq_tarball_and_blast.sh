#!/usr/bin/env bash

# ==============================================================================
# Script: 02_process_seq_tarball_and_blast.sh
# Description: Automates Sanger sequencing processing: extracts .seq files,
#              trims them, generates a FASTA, and runs a BLAST search.
# ==============================================================================

set -euo pipefail

# 1. Input Validation
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 data/query/FILE.tar" >&2
    exit 1
fi

INPUT_TAR="$1"

# Check to see if the input tarball file exists
[[ ! -f "${INPUT_TAR}" ]] && { echo "Error: Tarball not found." >&2; exit 1; }

# 2. Directory & Path Setup
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QUERY_DIR="${QUERY_DIR:-${REPO_ROOT}/data/query}"
OUTPUT_DIR="${OUTPUT_DIR:-${REPO_ROOT}/output}"
DB_PREFIX="${REPO_ROOT}/data/blastdb/fungi_ITS"
BLASTN_BIN="${BLASTN_BIN:-blastn}"

mkdir -p "${QUERY_DIR}" "${OUTPUT_DIR}"

# Naming: Use date prefix and strip common extensions (.tar, .tar.gz)
DATE_PREFIX="$(date +%Y%m%d)"
BASE_NAME="$(basename "${INPUT_TAR}")"
BASE_NO_EXT="${BASE_NAME%.tar*}"
# remove the date prefix if the filename already has it
BASE_NO_DATE="$(echo "${BASE_NO_EXT}" | sed -E 's/^[0-9]{8}-//')"
DATED_BASE="${DATE_PREFIX}-${BASE_NO_DATE}"

DATED_TAR="${QUERY_DIR}/${DATED_BASE}.tar"
EXTRACT_DIR="${QUERY_DIR}/${DATED_BASE}_seq_files"
FASTA="${QUERY_DIR}/${DATED_BASE}.fasta"
BLAST_OUT="${OUTPUT_DIR}/${DATED_BASE}.blast.tsv"

# 3. File Preparation
# Copy input to query folder with date prefix if it's not already there
[[ "${INPUT_TAR}" != "${DATED_TAR}" ]] && cp "${INPUT_TAR}" "${DATED_TAR}"

echo "Extracting sequences..."
mkdir -p "${EXTRACT_DIR}"
tar -xf "${DATED_TAR}" -C "${EXTRACT_DIR}"

# 4. Trim and Convert to FASTA
# 4.1 Collect and sort files (assumes lowercase .seq)
for f in $(ls "${EXTRACT_DIR}"/*.Seq | sort); do
    
    # 4.2 Get the filename without extension
    name=$(basename "$f" .Seq)
    
    # 4.3 Clean and Uppercase the sequence
    # Removes all whitespace and makes it readable for BLAST
    seq=$(tr -d '[:space:]' < "$f" | tr '[:lower:]' '[:upper:]')
    
    # 4.4 Trim (Skip 50, take next 800bp)
    trimmed="${seq:50:800}"
    
    # 4.5 Format and append to the FASTA file
    printf ">%s\n%s\n" "$name" "$trimmed"

done > "${FASTA}"

[[ ! -s "${FASTA}" ]] && { echo "Error: No sequences generated." >&2; exit 1; }

# 5. BLAST Search
# Verify the BLAST database exists before attempting search
if [[ ! -f "${DB_PREFIX}.nhr" && ! -f "${DB_PREFIX}.00.nhr" ]]; then
    echo "Error: BLAST DB not found at ${DB_PREFIX}" >&2
    exit 1
fi

if ! command -v "${BLASTN_BIN}" >/dev/null 2>&1; then
    echo "Error: blastn not found. Install BLAST+ on PATH or set BLASTN_BIN to the blastn executable path." >&2
    exit 1
fi

SUMMARY_OUT="${OUTPUT_DIR}/${DATED_BASE}.summary.txt"

echo "Running BLAST..."

"${BLASTN_BIN}" \
    -query "${FASTA}" \
    -db "${DB_PREFIX}" \
    -task blastn \
    -dust no \
    -max_target_seqs 10 \
	-max_hsps 1 \
	-outfmt "6 qseqid sseqid pident qcovs evalue bitscore stitle qstart qend sstart send qseq sseq" \
	| sort -k1,1 -k6,6gr \
	| awk -F '\t' '
BEGIN { OFS="\t" }
{
    query=$1

    if (++count[query] > 2)
        next

    print
}
' > "${BLAST_OUT}"

awk -F '\t' '
BEGIN { OFS="\t" }
{
	query=$1; subject_id=$2; pident=$3; qcov=$4; evalue=$5; bitscore=$6; title=$7;

    species_line = title
    #sub("^" subject_id "[[:space:]]+", "", species_line)

	n = split(species_line, words, " ")
    species = words[1] " " words[2]

    hits[query] = hits[query] \
        "  Hit " count[query] ":\n" \
        "    Species: " species "\n" \
        "    Full reference: " title "\n" \
        "    % identity: " pident "\n" \
        "    Query coverage: " qcov "%\n" \
        "    E-value: " evalue "\n" \
        "    Bitscore: " bitscore "\n"
}
END {
    for (query in hits) {
        print "Query: " query
        print hits[query]
        print ""
    }
}
' "${BLAST_OUT}" > "${SUMMARY_OUT}"

# 6. Optional: Cleanup intermediate extraction folder
rm -rf "${EXTRACT_DIR}"

echo "Done."
echo "Input tarball:  ${DATED_TAR}"
echo "Combined FASTA: ${FASTA}"
echo "BLAST table:    ${BLAST_OUT}"
echo "BLAST summary:  ${SUMMARY_OUT}"
