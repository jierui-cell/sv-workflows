process CUE_PAIR {
    tag "${meta.id}"
    cpus { settings.cpus as int }
    memory { "${settings.memory_gb} GB" }
    time { settings.time }
    maxForks (params.max_parallel_callers as int)
    scratch false
    module { params.modules_calling }
    publishDir "${params.outdir}/cue/${meta.id}", mode: 'copy'
    input:
    tuple val(meta), val(settings), path(tumor), path(tumor_bai), path(normal), path(normal_bai)
    path reference_fai
    path runner
    path benchmarker
    output:
    path 'output/*'
    script:
    def spec = groovy.json.JsonOutput.toJson([id:meta.id, root:projectDir.toString(),
        tumor_bam:meta.tumor_bam, normal_bam:meta.normal_bam,
        fai:file(params.reference).toString()+'.fai', settings:settings])
    def benchmarkSpec = groovy.json.JsonOutput.toJson([caller:'cue', pair_id:meta.id,
        dataset_id:params.dataset_id, reference:file(params.reference).toString(),
        allocated_cpus:task.cpus, requested_memory_gb:task.memory.toGiga(),
        benchmark_scope:'sequential tumor and normal inference through paired BAM-screen QC',
        output_class:'research_only_unvalidated_somatic_candidates', settings:settings])
    """
    cat > task.json <<'TASK_SETTINGS'
    ${spec}
    TASK_SETTINGS
    cat > benchmark.meta.json <<'BENCHMARK_SETTINGS'
    ${benchmarkSpec}
    BENCHMARK_SETTINGS
    mkdir -p output
    python3 '${benchmarker}' output/benchmark.json benchmark.meta.json -- python3 '${runner}' task.json
    """
    stub:
    """
    mkdir -p output
    echo 'STUB: no Cue calls produced' > output/STUB.txt
    """
}
