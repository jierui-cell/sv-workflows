# HPC Deployment

Set module names, account/partition and caller resources in your cohort YAML. The launcher loads the modules before starting Nextflow; each compute task loads the same list.

## Launch on Slurm

Submit from the repository root:

~~~bash
mkdir -p logs
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=7-00:00:00 --output=logs/cohort-%j.log \
  scripts/run.slurm.sh callers --config config/private-cohort.yaml
~~~

`run.slurm.sh` selects the repository’s `.venv` and starts `workflow.py`, which loads modules and invokes Nextflow. Nextflow submits preflight and caller jobs using the YAML account/partition. No separate `module load` step is needed.

The driver stays active for the cohort; each caller’s YAML `time` limits an individual task. Set durations within your site’s queue and maintenance limits.

Submitting from elsewhere? Set `SV_WORKFLOW_ROOT` to the checkout. `SV_PYTHON` selects a different launcher Python.

## Capacity and Concurrency

The cohort template permits eight simultaneous pairs **per caller** (`max_parallel_callers: 8`). The scheduler queues work until resources are available.

Nextflow’s default submission window is 24 outstanding tasks. To allow all five callers to reach eight tasks each, copy `config/site.nextflow.example.config` to `config/private-site.config` and set:

~~~groovy
executor.queueSize = 48
~~~

Add `--site-config config/private-site.config` to your launch command. Tune caller CPU/RAM in YAML and concurrency to match your allocation and storage throughput. Keep `work/` on persistent shared storage for resume.

## Other Schedulers and Local Execution

Inside your site’s driver allocation, activate `.venv` and choose the scheduler profile:

~~~bash
source .venv/bin/activate
python scripts/workflow.py callers --config config/private-cohort.yaml \
  --profile pbs --site-config config/private-site.config
~~~

Use `lsf` or `sge` instead of `pbs` as needed; put queue/accounting options in the site config. These profiles are site templates; Slurm is compute-node tested.

For a workstation or allocated compute node, use `--profile local --allow-local`. This runs on the current machine; it does not request a node.

## Resume and Troubleshooting

Find launcher logs under `logs/callers-TIMESTAMP/` and task errors in the reported work directory’s `.command.err`. Fix the cause, then [resume the same run](../README.md#4-optional-restart-an-interrupted-run). Keep `work/` and `.nextflow/`; matching completed tasks are reused.

The optional [entry-point tests](validation.md#entry-point-tests) exercise preflight, caller stubs and resume.
