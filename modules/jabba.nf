process JABBA {
    tag "${meta.id}"
    cpus { params.jabba.cpus as int }
    memory { "${params.jabba.memory_gb} GB" }
    time { params.jabba.time }
    maxForks (params.max_parallel_callers as int)
    scratch false
    module { params.modules_calling }
    publishDir "${params.outdir}/jabba/${meta.id}", mode: 'copy', pattern: 'output/*'
    input:
    tuple val(meta), path(junctions), path(coverage)
    path image
    path runner
    path benchmarker
    output:
    path 'output/*'
    script:
    def options = params.jabba.options
    def spec = groovy.json.JsonOutput.toJson([pair_id:meta.id, junctions:junctions.toString(),
        coverage:coverage.toString(), image:image.toString(), project_dir:projectDir.toString(),
        purity:meta.purity, ploidy:meta.ploidy, settings:params.jabba])
    def benchmarkSpec = groovy.json.JsonOutput.toJson([caller:'jabba', pair_id:meta.id,
        dataset_id:params.dataset_id, reference:file(params.reference).toString(),
        allocated_cpus:task.cpus, requested_memory_gb:task.memory.toGiga(),
        benchmark_scope:'PASS BND filtering and JaBbA graph inference from precomputed SvABA junctions and corrected coverage through output QC; excludes upstream SvABA and coverage preparation',
        output_class:'graph-selected junctions, dependent on SvABA input',
        settings:[options:options, image:params.jabba.image,
                  solver:params.jabba.solver, cplex_dir:params.jabba.cplex_dir,
                  upstream_junction_caller:params.jabba.upstream_junction_caller,
                  coverage_method:params.jabba.coverage_method,
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
    echo '{"status":"stub_not_benchmarked"}' > output/benchmark.json
    echo 'STUB: no graph produced' > output/STUB.txt
    """
}
