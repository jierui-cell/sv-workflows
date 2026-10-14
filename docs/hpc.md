# HPC profiles

The launcher supports `--profile slurm`, `pbs`, `lsf`, `sge`, and `local`; all share scientific defaults. Slurm is the production-tested backend. Other scheduler profiles are configuration templates, not yet validated on those schedulers. Check the resolved scheduler command and perform a one-pair pilot at your site.

Use YAML for scientific parameters and a private Nextflow site config for scheduler-specific behavior. Example: `--site-config config/site.nextflow.example.config`. PBS/LSF/SGE accounting flags belong in `scheduler_options`; do not assume Slurm flags work elsewhere. Some sites require additional parallel-environment or memory directives.

For Slurm, submit `scripts/run.slurm.sh` using your own account, partition, time and driver memory flags. Child-task account/partition come from YAML. Never run real callers directly on a login node. Local execution requires `--allow-local`, an explicit acknowledgement that you are using an allocated node or workstation; it does not allocate resources itself. Numerical signature controls also belong on compute nodes.

Submit from the repository root. Because Slurm copies batch scripts into its spool, the driver resolves the repository using `SV_WORKFLOW_ROOT` when provided, otherwise `SLURM_SUBMIT_DIR`. It checks the repository files before launching. `SV_PYTHON` can specify an absolute environment interpreter. Preflight recognizes Slurm, PBS, LSF and SGE job IDs; local execution uses the explicit acknowledgement passed by the launcher. Environment-recognition tests are not evidence of actual execution on those other schedulers.

## Capacity and concurrency

`max_parallel_callers` caps each caller process independently. Start with 1–3; raise it only after a successful same-configuration pilot. `executor.queueSize` provides a separate outstanding-task bound. Eight tasks per caller is not necessarily the fastest configuration when shared BAM/reference I/O is saturated. RAM, scratch capacity, filesystem throughput, CPU allocation and queue availability are separate constraints.

Tasks default to persistent shared work directories for dependable resume and container access. Node-local scratch is optional and site-specific: verify its filesystem, available capacity, mounts and cleanup lifecycle first. Do not set a generic `/tmp` path for large WGS intermediates. Final results/logs must remain on persistent storage.

## Failures and resume

Use the official Nextflow launcher, not a minimal `java -jar` wrapper. The official launcher supplies JVM module-access flags needed by its libraries. The entry-point test verifies actual cached resume; successful task execution alone does not prove task records were persisted. Do not relax cache checks or delete previous work to hide a failed resume test.

The strategy is `finish`, not `ignore`: submitted tasks can finish, then the run exits unsuccessfully. There are no blanket automatic retries. Fix the logged cause and use `--resume SESSION_UUID`, with unchanged scientific inputs/settings and the original work directory. Resume only the intended session; concurrent launchers make an unqualified latest-session resume risky. For a parameter experiment choose a new output directory rather than resuming into old results.

Long cohorts need a driver allocation longer than the entire run, not only one task. Restarting the driver does not recreate deleted work artifacts. Configure optional site modules with YAML or your site initialization; no lab account, partition or module is hard-coded in the published launch path.
