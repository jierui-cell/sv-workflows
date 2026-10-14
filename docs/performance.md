# Comparable runtime and memory

All integrated caller processes use the same `measure_benchmark.py` wrapper. `elapsed_seconds` measures the recorded post-preflight caller task boundary, including preparation/filter/QC within that boundary. Queue delay and shared upstream reference preparation are not included.

On the validated Slurm site, peak resident memory is sampled from the verified job cgroup, excluding unmapped page cache. On other backends, a job-isolated cgroup has not been demonstrated; the public recorder uses command-session process-tree RSS and records a different `memory_method_id`. This excludes unmapped cache but can double-count shared pages. Process PSS is also recorded when `/proc` permissions permit; missing values remain missing. One-second sampling can miss brief peaks.

GNU `time -v`/`wait4` stage maxima and Slurm accounting are supplementary measurements, not interchangeable with simultaneous process/job peak resident memory. In particular, this site's historical Slurm MaxRSS was dominated by filesystem cache and must not be interpreted as application RSS.

Use `scripts/summarize_benchmarks.py --help` to export balanced pair/caller tables. The exporter rejects mixed datasets, references, measurement methods, benchmark boundaries, duplicate outputs and inconsistent settings within a caller. Hold CPUs, scientific parameters, filesystem conditions and overlap policy constant; report algorithm-specific RAM allocations separately. Concurrent execution helps throughput but can bias runtime comparison through I/O contention.

Historical results lack uniform instrumentation and are not retroactively made comparable by this refactor. Measure new runs using the new wrapper and retain native logs/QC alongside summaries.
