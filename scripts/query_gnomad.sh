#!/usr/bin/env bash
# Fetch gnomAD AF fields for a list of sites from a remote, tabix-indexed VCF.
# Nearby sites are merged into short intervals and queried in parallel, which is
# much faster than one HTTP seek per site. Only exact positions are kept.
# Usage: query_gnomad.sh <sites.tsv (CHROM POS)> <vcf url> <out.tsv> <threads>
set -euo pipefail
sites=$1; url=$2; out=$(realpath -m "$3"); threads=${4:-4}
work=$(mktemp -d "${out}.tmp.XXXX")
trap 'rm -rf "$work"' EXIT

sort -k1,1 -k2,2n -u "$sites" > "$work/sites"
echo -e "CHROM\tPOS\tREF\tALT\tFILTER\tAF\tAF_grpmax\tnhomalt" > "$out"
[ -s "$work/sites" ] || exit 0

# merge sites closer than 2 kb into one region string chr:start-end
awk -v pad=2000 '{ if ($1==c && $2-e<=pad) e=$2; else { if (c!="") print c":"s"-"e; c=$1; s=$2; e=$2 } }
     END { if (c!="") print c":"s"-"e }' "$work/sites" > "$work/regions"

# Download the index once: parallel jobs that each cache it in the working
# directory overwrite each other's copy and silently return nothing.
cd "$work"
curl -sSfL -o remote.tbi "${url}.tbi"

# one bcftools call per region, each writing its own part file
nl -ba regions | xargs -P "$threads" -L 1 sh -c \
  'bcftools query -r "$1" -f "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%AF\t%AF_grpmax\t%nhomalt\n" \
     "'"$url"'##idx##remote.tbi" > part_$0.tsv || exit 255'

# keep only records at the requested positions
cat part_*.tsv | awk 'NR==FNR { want[$1":"$2]=1; next } ($1":"$2) in want' sites - \
  | sort -k1,1 -k2,2n | uniq >> "$out"   # plain "sort -u -k" would drop other alleles at a site
