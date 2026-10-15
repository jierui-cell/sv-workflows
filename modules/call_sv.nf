process CALL_SV {
    tag "${tool}:${meta.id}"
    cpus { settings.cpus as int }
    memory { "${settings.memory_gb} GB" }
    time { settings.time }
    maxForks (params.max_parallel_callers as int)
    scratch false
    module { ((params.modules ?: []) + (params.modules_calling ?: '').tokenize(':')).unique().join(':') }
    publishDir "${params.outdir}/${tool}/${meta.id}", mode: 'copy', pattern: 'output/*'
    input:
    tuple val(meta), val(tool), val(settings), path(image), path(tumor), path(tumor_bai), path(normal), path(normal_bai)
    path reference_files
    path runner
    path benchmarker
    path exclude_regions
    val exclude_sha256
    output:
    path 'output/*'
    script:
    def effective = settings + [cpus:task.cpus, memory_gb:task.memory.toGiga()]
    def spec = groovy.json.JsonOutput.toJson([tool:tool, settings:effective, id:meta.id,
        image:image.toString(), reference:reference_files[0].toString(),
        tumor:tumor.toString(), normal:normal.toString(), project_dir:projectDir.toString(),
        genome_build:params.genome_build, exclude_bed:exclude_regions ? exclude_regions[0].toString() : null])
    def benchmarkSpec = groovy.json.JsonOutput.toJson([caller:tool, pair_id:meta.id,
        dataset_id:params.dataset_id, reference:file(params.reference).toString(),
        allocated_cpus:task.cpus, requested_memory_gb:task.memory.toGiga(),
        benchmark_scope:'caller preparation through final output QC',
        output_class:tool == 'gridss' ? (settings.somatic_filter ? 'somatic_filtered_no_pon' : 'raw_joint_not_somatic') : 'somatic',
        settings:settings + [genome_build:params.genome_build, exclude_bed:params.exclude_bed,
                            exclude_bed_sha256:exclude_sha256]])
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
    def spec = groovy.json.JsonOutput.toJson([tool:tool, settings:settings, id:meta.id, project_dir:projectDir.toString(),
        reference:reference_files[0].toString(), exclude_bed:exclude_regions ? exclude_regions[0].toString() : null])
    """
    cat > task.json <<'TASK_SETTINGS'
    ${spec}
    TASK_SETTINGS
    python3 '${runner}' task.json --stub
    """
}
