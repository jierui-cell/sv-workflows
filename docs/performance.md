# Runtime and Memory Comparisons

All five callers use `measure_benchmark.py` with the same task boundary: after preflight, from caller preparation/configuration through filtering and final QC. `elapsed_seconds` excludes queue delay and shared reference preparation.

All callers save benchmarks at `CALLER/PAIR/output/benchmark.json`. Performance exports and signature-input discovery use this same location.

## What Is Measured

On a readable, verified Slurm job cgroup, the recorder samples resident memory excluding unmapped filesystem cache. Otherwise it records command-session process-tree RSS and a distinct `memory_method_id`. Process PSS is included where available. Sampling is once per second.

If cgroup access fails, the caller continues with process-tree measurement. Warnings and measurement status are saved; unavailable values remain missing. If sampling changes method mid-run, the primary reported peak uses the process-tree series throughout.

GNU time/wait4 stage measurements and Slurm accounting are supplementary. Slurm MaxRSS can include filesystem cache, so use the primary benchmark field for resident-memory comparison.

## Export a Comparison Table

~~~bash
python scripts/summarize_benchmarks.py --help
~~~

The exporter prepares balanced pair/caller tables and checks dataset/reference identity, method IDs, task scope, duplicates and settings. Compare one measurement method at a time. Report caller CPU/RAM allocations alongside elapsed time and peak resident memory.

Concurrency improves throughput; shared I/O can influence individual runtimes. For controlled benchmarks, hold inputs, CPUs and filesystem conditions consistent. Historical runs without this wrapper retain their original measurement scope.
