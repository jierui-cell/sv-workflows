process CALL_SV {
    tag "${tool}:${meta.id}"
    cpus { settings.cpus as int }
    memory { "${settings.memory_gb} GB" }
    time { settings.time }
    maxForks (params.max_parallel_callers as int)
    scratch false
    module { params.modules_calling }
    publishDir "${params.outdir}/${tool}/${meta.id}", mode: 'copy', pattern: 'output/*',
        saveAs: { name -> params.compact_publish ? name.tokenize('/').last() : name }
    input:
    tuple val(meta), val(tool), val(settings), path(image), path(tumor), path(tumor_bai), path(normal), path(normal_bai)
    path reference_files
    path runner
    path benchmarker
    output:
    path 'output/*'
    script:
    def effective = settings + [cpus:task.cpus, memory_gb:task.memory.toGiga()]
    def spec = groovy.json.JsonOutput.toJson([tool:tool, settings:effective, id:meta.id,
        image:image.toString(), reference:reference_files[0].toString(),
        tumor:tumor.toString(), normal:normal.toString(), project_dir:projectDir.toString()])
    def benchmarkSpec = groovy.json.JsonOutput.toJson([caller:tool, pair_id:meta.id,
        dataset_id:params.dataset_id, reference:file(params.reference).toString(),
        allocated_cpus:task.cpus, requested_memory_gb:task.memory.toGiga(),
        benchmark_scope:'caller preparation through final output QC',
        output_class:tool == 'gridss' ? (settings.somatic_filter ? 'somatic_filtered_no_pon' : 'raw_joint_not_somatic') : 'somatic',
        settings:settings])
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
    def spec = groovy.json.JsonOutput.toJson([tool:tool, settings:settings, id:meta.id, project_dir:projectDir.toString()])
    """
    cat > task.json <<'TASK_SETTINGS'
    ${spec}
    TASK_SETTINGS
    python3 '${runner}' task.json --stub
    """
}
