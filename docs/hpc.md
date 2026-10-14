# HPC Deployment

Use the same YAML and scientific settings across `slurm`, `pbs`, `lsf`, `sge` or `local` profiles. Slurm has been exercised on compute nodes; the other scheduler profiles provide configurable site templates.

## Launch on Slurm

Activate Python and load tools/modules before submission. From the repository root:

~~~bash
mkdir -p logs
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  --cpus-per-task=1 --mem=4G --time=7-00:00:00 --output=logs/cohort-%j.log \
  scripts/run.slurm.sh callers --config config/private-cohort.yaml
~~~

The driver schedules per-pair preflight and caller tasks. Child-job account/partition come from your YAML. Its walltime covers the entire cohort, while each caller's `time` sets the individual task limit. Choose site-appropriate durations; upcoming maintenance can delay allocations that extend into the reserved window.

The driver resolves the checkout from `SLURM_SUBMIT_DIR`. To submit from another directory, export `SV_WORKFLOW_ROOT=/absolute/path/to/sv-workflows`. `SV_PYTHON` can select the launcher's Python executable.

## Capacity and Concurrency

The cohort template uses `max_parallel_callers: 8`: eight pairs for each enabled caller. Five callers permit 40 tasks. The scheduler fills available nodes and queues the remainder.

The default Nextflow submission window is 24 outstanding tasks, including preflight. To increase that window, put this in your private site config:

~~~groovy
executor.queueSize = 48
~~~

Then add `--site-config /path/to/private/site.config` to the launch command. Increase or decrease `max_parallel_callers` independently. Caller CPU/RAM requests are in YAML; shared filesystem throughput and scratch capacity also affect throughput.

Tasks use persistent shared work directories for resume. A site can configure node-local scratch through its Nextflow settings and container mounts; keep results/logs on persistent storage.

## Other Schedulers and Local Execution

Run the launcher inside your site's usual driver allocation:

~~~bash
python scripts/workflow.py callers --config config/private-cohort.yaml \
  --profile pbs --site-config /path/to/private/site.config
~~~

Replace `pbs` with `lsf` or `sge`. Queue/accounting/parallel-environment settings belong in the site config and `scheduler_options`. Preflight recognizes these scheduler job environments.

On a workstation or an already allocated compute node:

~~~bash
python scripts/workflow.py callers --config config/private-cohort.yaml \
  --profile local --allow-local
~~~

`--allow-local` acknowledges the current machine; it does not allocate a node. HPC caller computation belongs on compute nodes.

## Resume and Troubleshooting

Use the official Nextflow launcher, which supplies the JVM settings needed for persistent task caching.

If a task fails, already submitted tasks finish and Nextflow reports the failure. Find the cause in the driver log and task `.command.err`, then relaunch the original command with `--resume SESSION_UUID`. Keep the original work directory. A new scientific experiment uses a new output directory.

The optional [entry-point setup test](validation.md#entry-point-tests) checks preflight and cached resume. Module names can be set in YAML; leave them null when tools are already on PATH.
