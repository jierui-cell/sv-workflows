#!/usr/bin/env bash
set -Eeuo pipefail
pair_id=$1
reference_fai=$2
tumor=$3
normal=$4
step=initialization
trap 'rc=$?; printf "FAILED step=%s line=%s exit=%s time=%s\n" "$step" "$LINENO" "$rc" "$(date -Is)" >&2; exit "$rc"' ERR
: "${SLURM_JOB_ID:?Validation must run on a compute node through Slurm}"
printf 'START pair=%s job=%s host=%s time=%s\n' "$pair_id" "$SLURM_JOB_ID" "$(hostname)" "$(date -Is)"
samtools --version
cut -f1,2 "$reference_fai" > reference.contigs.tsv
for role in tumor normal; do
    bam=$tumor
    if [[ $role == normal ]]; then bam=$normal; fi
    step="$role.quickcheck"
    samtools quickcheck -v "$bam"
    step="$role.header"
    samtools view -H "$bam" > "$role.header.sam"
    awk -F '\t' '$1=="@HD" {for(i=2;i<=NF;i++) if($i=="SO:coordinate") ok=1} END {exit !ok}' "$role.header.sam"
    awk -F '\t' '$1=="@SQ" {sn="";ln=""; for(i=2;i<=NF;i++){if($i~/^SN:/)sn=substr($i,4);if($i~/^LN:/)ln=substr($i,4)} print sn "\t" ln}' "$role.header.sam" > "$role.contigs.tsv"
    diff -u reference.contigs.tsv "$role.contigs.tsv"
    awk -F '\t' '$1=="@RG" {for(i=2;i<=NF;i++) if($i~/^SM:/) print substr($i,4)}' "$role.header.sam" | sort -u > "$role.samples.txt"
    test "$(wc -l < "$role.samples.txt")" -eq 1
    step="$role.indexed_query"
    samtools view -c "$bam" 'chr1:1000000-1001000' > "$role.interval_read_count.txt"
done
step=distinct_sample_names
if cmp -s tumor.samples.txt normal.samples.txt; then
    echo 'Tumor and normal have identical BAM SM names; review before calling.' >&2
    exit 1
fi
printf 'PASS pair=%s time=%s\n' "$pair_id" "$(date -Is)"
printf 'pair_id\tstatus\tscope\n%s\tPASS\theader_dictionary_EOF_and_small_indexed_query_only\n' "$pair_id" > "$pair_id.validation.tsv"
