#!/usr/bin/env bash

set -euo pipefail

# Arguments
INPUT="$1"     # VCF (facets), .segments.txt (ascat), or .cns (cnvkit)
CALLER="$2"    # facets | ascat | cnvkit
OUTDIR="$3"    # output directory

# Directory for output .seg files
mkdir -p "$OUTDIR"

# Derive sample ID from file name
BASENAME=$(basename "$INPUT")
case "$CALLER" in
    facets)
        SAMPLE_ID=$(echo "$BASENAME" | sed -E 's/^[^_]+_[^_]+_([^\.]+).*/\1/')
        ;;
    ascat)
        SAMPLE_ID=$(echo "$BASENAME" | sed -E 's/\.segments\.txt$//')
        ;;
    cnvkit)
        SAMPLE_ID=$(echo "$BASENAME" | sed -E 's/\.cns$//')
        ;;
    *)
        echo "[ERROR] Unknown caller: $CALLER" >&2
        exit 1
        ;;
esac

SEG_FILE="${OUTDIR}/${SAMPLE_ID}.seg"

# Conversion based on caller
case "$CALLER" in
    facets)
        # Ensure index exists
        if [[ ! -f "${INPUT}.tbi" ]]; then
            echo "[INFO] Indexing $INPUT..." >&2
            tabix -p vcf "$INPUT"
        fi

        # Extract CNV segments from VCF
        {
            echo -e "Sample\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean"
            bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/NUM_MARK\t%INFO/CNLR_MEDIAN\n' "$INPUT" 2>/dev/null | \
                awk -v sample="$SAMPLE_ID" 'BEGIN {OFS="\t"} {print sample, $1, $2, $3, $4, $5}'
        } > "$SEG_FILE"
        ;;

    ascat)
        # Reformat ASCAT .segments.txt
        {
        echo -e "Sample\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean"
        awk 'NR>1 {
            sample=$1
            chr=$2
            start=$3
            end=$4
            nMajor=$5
            nMinor=$6
            cn=nMajor+nMinor
            # Avoid log2(0): if CN=0, set Segment_Mean to -Inf
            if (cn>0) {
                segmean=log(cn/2)/log(2)
            } else {
                segmean=-Inf
            }
            print sample, chr, start, end, "NA", segmean
        }' OFS="\t" "$INPUT"
        } > "$SEG_FILE"
        ;;

    cnvkit)
        # Reformat CNVkit .cns
        cnvkit.py segment "$INPUT" -o "${OUTDIR}/${SAMPLE_ID}.cns"
        CNS="${OUTDIR}/${SAMPLE_ID}.cns"
        {
        echo -e "Sample\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean"
        awk -v sample="$SAMPLE_ID" 'NR>1 {
            print sample, $1, $2, $3, $7, $6
        }' OFS="\t" "$INPUT"
        } > "$SEG_FILE"
        ;;

esac

echo "[INFO] Converted $INPUT -> $SEG_FILE" >&2