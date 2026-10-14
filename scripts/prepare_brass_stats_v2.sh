#!/bin/bash -l
# Generic safe pair IDs; preserve historical prepare_brass_stats.sh for existing runs.
#SBATCH --time=06:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
set -Eeuo pipefail
: "${SV_JOB_ID:?Launch through scripts/workflow.py}"
cd "${SV_PROJECT_DIR:?}"
root=$PWD
pair=${1:?Expected pair ID}
[[ $pair =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]]
run_dir=${BRASS_STATS_RUN_DIR:-"$root/logs/brass-stats-$SV_JOB_ID"}
out=${BRASS_STATS_OUTDIR:-"$root/results/brass-preparation/$pair/stats-$SV_JOB_ID"}
mkdir -p "$run_dir" "$out"
exec > >(tee -a "$run_dir/driver.log") 2>&1
step=initialize
trap 'rc=$?; echo "FAILED step=$step line=$LINENO exit=$rc time=$(date -Is)"; exit "$rc"' ERR
echo "START job=$SV_JOB_ID host=$(hostname) time=$(date -Is)"
cp scripts/prepare_brass_stats_v2.sh "$run_dir/"
cp "${BRASS_MANIFEST:?}" "$run_dir/samplesheet.tsv"
image=${BRASS_IMAGE:?}
sha256sum "$image" > "$run_dir/image.sha256"
for role in tumor normal; do
    col=8
    [[ $role == normal ]] && col=13
    bam=$(python3 - "$run_dir/samplesheet.tsv" "$pair" "$role" <<'PY'
import csv,sys
rows=[r for r in csv.DictReader(open(sys.argv[1]),delimiter='\t') if r['patient_id']==sys.argv[2]]
assert len(rows)==1
print(rows[0][sys.argv[3]+'_bam'])
PY
)
    test -f "$bam"
    step="$role.bam_stats"
    bamdir=$(dirname "$bam")
    singularity exec --cleanenv --bind "$root:$root,$bamdir:$bamdir:ro" "$image" \
        bam_stats -i "$bam" -o "$out/$role.bam.bas" -@ "$SV_CPUS"
    test -s "$out/$role.bam.bas"
done
sha256sum "$out"/*.bas > "$run_dir/outputs.sha256"
echo "COMPLETE prerequisites_only=$out time=$(date -Is)"
