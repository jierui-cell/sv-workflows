process JABBA_COVERAGE {
    tag "${meta.id}"
    cpus { params.jabba.coverage_cpus as int }
    memory { "${params.jabba.coverage_memory_gb} GB" }
    time { params.jabba.coverage_time }
    maxForks ((params.jabba?.max_parallel_coverage ?: 1) as int)
    scratch false
    publishDir "${params.outdir}/jabba_coverage/${meta.id}", mode: 'copy', pattern: 'output/*'
    input:
    tuple val(meta), path(junctions), path(tumor), path(tumor_bai), path(normal), path(normal_bai)
    path gcmap
    path frag_image
    path cbs_image
    path cbs_script
    path runner
    path benchmarker
    output:
    tuple val(meta), path(junctions), path('output/coverage.rds'), emit: prepared
    path 'output/*', emit: artifacts
    script:
    def spec = groovy.json.JsonOutput.toJson([pair_id:meta.id, project_dir:projectDir.toString(),
        tumor:tumor.toString(), normal:normal.toString(), gcmap:gcmap.toString(),
        frag_image:frag_image.toString(), cbs_image:cbs_image.toString(),
        cbs_script:cbs_script.toString(), settings:params.jabba])
    def benchmarkSpec = groovy.json.JsonOutput.toJson([caller:'jabba_coverage', pair_id:meta.id,
        dataset_id:params.dataset_id, reference:file(params.reference).toString(),
        allocated_cpus:task.cpus, requested_memory_gb:task.memory.toGiga(),
        benchmark_scope:'fragCounter GC/mappability-corrected tumor and normal coverage plus upstream CBS tumor-normal ratio/segmentation; excludes SvABA and JaBbA optimization',
        output_class:'coverage_ratio_and_segmentation',
        settings:[coverage_method:params.jabba.coverage_method,
                  coverage_options:params.jabba.coverage_options,
                  fragcounter_image:params.jabba.fragcounter_image,
                  cbs_image:params.jabba.cbs_image,
                  cbs_source_commit:params.jabba.cbs_source_commit,
                  cpus:task.cpus, memory_gb:task.memory.toGiga()]])
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
    touch output/coverage.rds
    echo '{"status":"stub_not_benchmarked"}' > output/benchmark.json
    """
}
