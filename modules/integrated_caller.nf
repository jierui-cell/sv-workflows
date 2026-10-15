process INTEGRATED_CALLER {
    tag "${tool}:${meta.id}"
    cpus { settings.cpus as int }
    memory { "${settings.memory_gb} GB" }
    time { settings.time }
    maxForks (params.max_parallel_callers as int)
    scratch false
    module { ((params.modules ?: []) + (params.modules_calling ?: '').tokenize(':')).unique().join(':') }
    publishDir "${params.outdir}/${tool}/${meta.id}", mode: 'copy'
    input:
    tuple val(meta), val(tool), val(settings), path(manifest)
    path runner
    path meter
    path monitor
    path benchmarker
    output:
    path 'output/*'
    script:
    def spec = groovy.json.JsonOutput.toJson([tool:tool, id:meta.id,
        project_dir:projectDir.toString(), manifest:manifest.toString(), settings:settings,
        allow_local_execution:params.allow_local_execution])
    def benchmarkSpec = groovy.json.JsonOutput.toJson([caller:tool, pair_id:meta.id,
        dataset_id:params.dataset_id, reference:file(params.reference).toString(),
        allocated_cpus:task.cpus, requested_memory_gb:task.memory.toGiga(),
        benchmark_scope:'caller preparation through final output QC',
        output_class:tool == 'lumpy' ? 'research_candidates_unvalidated' : 'native_brass_bnd',
        settings:settings])
    """
    cat > task.json <<'TASK_SETTINGS'
    ${spec}
    TASK_SETTINGS
    cat > benchmark.meta.json <<'BENCHMARK_SETTINGS'
    ${benchmarkSpec}
    BENCHMARK_SETTINGS
    mkdir -p output
    python3 '${benchmarker}' output/benchmark.json benchmark.meta.json -- \
      python3 '${meter}' output/call.time.txt \
      python3 '${monitor}' output/resource-usage.json bash '${runner}' task.json
    """
    stub:
    """
    mkdir -p output
    echo 'STUB: no SV calls produced' > output/STUB.txt
    """
}
