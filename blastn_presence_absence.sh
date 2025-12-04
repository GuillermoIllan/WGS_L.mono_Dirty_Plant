#!/usr/bin/env bash
# This script performs a BLASTN search of a set of loci/alleles against a genome database
# and converts the results into a binary presence/absence matrix (genes × genomes).
# It is intended for summarizing gene presence across a panel of genomes.
#
# Steps:
#  1) Run BLASTN with identity and query coverage filters.
#  2) From the BLAST output, extract unique (gene, sample) pairs indicating presence.
#  3) Build sorted lists of genes and samples.
#  4) Construct a 1/0 presence/absence matrix (genes vs. genomes).
#
# Assumptions:
#  - The BLAST database contains genomes indexed in such a way that the sample ID can be
#    extracted from the subject title (stitle) as the first field before a "|" character
#    and/or before the first whitespace.
#  - Query IDs (qseqid) may encode allele variants using a suffix "_<number>", and the
#    base gene name is obtained by removing that suffix (e.g. "geneA_3" → "geneA").
#
# Usage (see function `usage` below) and typical command:
#   ./script.sh -db db_genomes -q loci.fasta -o results -n mytag
#
# The final output matrix is written to:
#   $OUTDIR/matrix_<TAG>.tsv

set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 -db <DB_PREFIX> -q <GENES_FASTA> -o <OUT_DIR> -n <TAG> [-t 8] [--minid 90] [--minqcov 80] [--task megablast|dc-megablast]

  Required arguments:
    -db       Prefix of the BLAST database (e.g. 'db_genomes')
    -q        FASTA file containing ALL loci/alleles (BLAST query)
    -o        Output directory
    -n        Run tag (used to name output files)

  Optional arguments:
    -t        Number of threads (default: 8)
    --minid   Minimum percent identity for BLAST hits (default: 90)
    --minqcov Minimum query coverage percentage for HSPs (default: 80)
    --task    BLASTN task: megablast | dc-megablast (default: megablast)

This script runs BLASTN against a genome database, filters hits by identity and
query coverage, and converts the results into a binary gene presence/absence
matrix (genes × genomes).
EOF
  exit 1
}

# Initialize parameters with defaults
DB=""
QUERY=""
OUTDIR=""
TAG=""

THREADS=8
MINID=90
MINQCOV=80
TASK="megablast"

# -------------------------
# Parse command-line args
# -------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -db) DB="$2"; shift 2 ;;
    -q)  QUERY="$2"; shift 2 ;;
    -o)  OUTDIR="$2"; shift 2 ;;
    -n)  TAG="$2"; shift 2 ;;
    -t)  THREADS="$2"; shift 2 ;;
    --minid)   MINID="$2"; shift 2 ;;
    --minqcov) MINQCOV="$2"; shift 2 ;;
    --task)    TASK="$2"; shift 2 ;;
    -h|--help) usage ;;
    *)
      echo "Unknown option: $1"
      usage
      ;;
  esac
done

# Check that required arguments are provided
[[ -n "$DB" && -n "$QUERY" && -n "$OUTDIR" && -n "$TAG" ]] || usage
mkdir -p "$OUTDIR"

echo "### 1) Running BLASTN with identity and query coverage filters"
blastn -task "$TASK" \
  -query "$QUERY" \
  -db "$DB" \
  -evalue 1e-20 \
  -perc_identity "$MINID" \
  -qcov_hsp_perc "$MINQCOV" \
  -num_threads "$THREADS" \
  -out "$OUTDIR/all.${TAG}.blast6.tsv" \
  -outfmt "6 qseqid sseqid stitle pident length qlen qcovhsp bitscore evalue"

echo "### 2) Extracting (GENE, SAMPLE) pairs and collapsing to presence"
# We assume that:
#  - The sample ID is encoded in the subject title (stitle, column 3) as the
#    first token before '|' and/or before the first whitespace.
#  - qseqid may contain an allele suffix of the form *_<number>. If present,
#    that suffix is removed to obtain the base gene name.
#
# The AWK command reads each BLAST hit and:
#  - Derives 'st' (sample ID) from stitle.
#  - Derives 'gene' from qseqid (removing allele suffix).
#  - Prints unique (gene, sample) pairs.
awk 'BEGIN{FS=OFS="\t"}
{
  q  = $1      # qseqid (query ID)
  st = $3      # stitle (subject title, full description)

  # Extract sample ID from stitle:
  #   - remove everything after the first "|"
  #   - then remove everything after the first space
  sub(/\|.*/, "", st)
  sub(/ .*/,  "", st)

  # Extract base gene name from qseqid by removing trailing "_<number>", if present
  gene = q
  if (match(q, /_[0-9]+$/)) {
    gene = substr(q, 1, RSTART - 1)
  }

  print gene, st
}' "$OUTDIR/all.${TAG}.blast6.tsv" \
| sort -u > "$OUTDIR/${TAG}.pairs.tsv"

echo "### 3) Generating sorted lists of genes and samples"
cut -f1 "$OUTDIR/${TAG}.pairs.tsv" | sort -u > "$OUTDIR/${TAG}.genes.list"
cut -f2 "$OUTDIR/${TAG}.pairs.tsv" | sort -u > "$OUTDIR/${TAG}.samples.list"

echo "### 4) Building binary presence/absence matrix (genes × genomes; 1/0)"
# This AWK block:
#  - Reads the (gene, sample) pairs and stores presence in an associative array.
#  - Reads the ordered gene and sample lists.
#  - Prints a matrix with:
#       rows    = genes
#       columns = samples
#       cells   = 1 if (gene, sample) is present; 0 otherwise.
awk -v pairs="$OUTDIR/${TAG}.pairs.tsv" '
  BEGIN {
    FS  = OFS = "\t"
  }
  # First file: pairs of (gene, sample) → record presence.
  FILENAME == ARGV[1] {
    pres[$1, $2] = 1
    next
  }
  # Second file: list of genes
  FILENAME == ARGV[2] {
    genes[++ng] = $1
    next
  }
  # Third file: list of samples
  FILENAME == ARGV[3] {
    samples[++ns] = $1
    next
  }
  END {
    # Print header: "Gene" followed by all sample IDs
    printf "Gene"
    for (i = 1; i <= ns; i++) {
      printf OFS samples[i]
    }
    printf "\n"

    # Print one row per gene
    for (i = 1; i <= ng; i++) {
      g = genes[i]
      printf g
      for (j = 1; j <= ns; j++) {
        s = samples[j]
        # Check if this (gene, sample) combination was observed
        printf OFS (((g SUBSEP s) in pres) ? 1 : 0)
      }
      printf "\n"
    }
  }
' "$OUTDIR/${TAG}.pairs.tsv" "$OUTDIR/${TAG}.genes.list" "$OUTDIR/${TAG}.samples.list" \
> "$OUTDIR/matrix_${TAG}.tsv"

echo "✅ Presence/absence matrix written to: $OUTDIR/matrix_${TAG}.tsv"

