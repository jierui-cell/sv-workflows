nextflow.enable.dsl = 2

params.manifest = "${projectDir}/input_manifest.tsv"
params.method = 'cosmic_sv32'
params.genome_build = 'GRCh37'
params.sequencing_type = 'WGS'
params.pcawg_mask_policy = 'local_cen_telo_exclusion'
params.outdir = ''
params.python = "${projectDir}/.venv/bin/python"
params.groups = 'all'
params.pair_id = ''
params.pair_ids = ''
params.pcawg_stage = 'signatures'
params.pcawg_min_size_bp = 1
params.pcawg_timing_source = 'none'  // Optional annotation: none | reference
params.pcawg_timing_reference = "${projectDir}/../resources/replication_timing/pcawg_mean_tertiles"
params.clustersv_source = "${projectDir}/../resources/clustersv/source"
params.pcawg_rscript = "${projectDir}/pcawg/Rscript.sh"
params.pcawg_fdr_first = 0.01
params.pcawg_fdr_final = 0.05
params.pcawg_footprint_fdr = 0.01
params.pcawg_rank_min = 2
params.pcawg_rank_max = 10
params.pcawg_rank = 0
params.pcawg_nrun = 100
params.pcawg_seed = 1234
params.pcawg_nmf_epsilon = 2.220446049250313e-16
params.pcawg_nmf_max_iter = 2000
params.min_size_bp = 1000
params.cosmic_version = 3.6
params.skip_fit = false
params.account = null
params.partition = null
params.max_parallel_groups = 4
params.group_cpus = 4
params.group_memory_gb = 32
params.group_time = '4h'

process SIGNATURE_GROUP {
    tag "${cohort}:${caller}"
    cpus { params.group_cpus as int }
    memory { "${params.group_memory_gb} GB" }
    time { params.group_time }
    maxForks (params.max_parallel_groups as int)
    publishDir {
        def root = file(params.outdir ?: "${projectDir}/../results/${params.method == 'cosmic_sv32' ? 'sv_signatures_corrected' : 'sv_signatures_pcawg'}").toAbsolutePath()
        params.pair_id ? "${root}/smoke/${cohort}/${caller}/${params.pair_id}" :
            (params.pair_ids ? "${root}/pilot/${cohort}/${caller}" : "${root}/${cohort}/${caller}")
    }, mode: 'copy'
    input:
    tuple val(cohort), val(caller), val(expected)
    path manifest
    path runner
    path helper_files
    output:
    tuple val(cohort), val(caller), path('group_output')
    script:
    def pythonExecutable = params.python == 'python3' ? 'python3' : file(params.python).toAbsolutePath().toString()
    def pilot = params.pair_id ? "--pair-ids ${params.pair_id}" : (params.pair_ids ? "--pair-ids ${params.pair_ids}" : '')
    def skip = params.skip_fit ? '--skip-fit' : ''
    def pcawg = params.method == 'pcawg_clustersv'
    def root = projectDir.parent.toAbsolutePath().toString()
    def methodArgs = pcawg ? """
        --clustersv-source '${file(params.clustersv_source).toAbsolutePath()}'
        --upstream-lock upstream.json --rscript '${file(params.pcawg_rscript).toAbsolutePath()}'
        --cluster-wrapper cluster.R --nmf-wrapper nmf.R
        --memory-monitor measure_process_tree.py
        --stage '${params.pcawg_stage}' --min-size-bp ${params.pcawg_min_size_bp}
        --genome-build '${params.genome_build}' --mask-policy '${params.pcawg_mask_policy}'
        --fdr-first ${params.pcawg_fdr_first} --fdr-final ${params.pcawg_fdr_final}
        --footprint-fdr ${params.pcawg_footprint_fdr} --rank-min ${params.pcawg_rank_min}
        --rank-max ${params.pcawg_rank_max} --rank ${params.pcawg_rank}
        --nrun ${params.pcawg_nrun} --seed ${params.pcawg_seed}
        --nmf-epsilon ${params.pcawg_nmf_epsilon} --nmf-max-iter ${params.pcawg_nmf_max_iter}
        ${params.pcawg_timing_source == 'reference' ? "--timing-reference '${file(params.pcawg_timing_reference).toAbsolutePath()}' --timing-wrapper replication_timing.py" : ''}
    """.stripIndent().replaceAll('\\s*\\n\\s*', ' ') : "--genome-build '${params.genome_build}' --min-size-bp ${params.min_size_bp} --cosmic-version ${params.cosmic_version} ${skip}"
    def control = pcawg ? "test_clustersv.py --pipeline-root '${root}' --clustering-only --clustersv-source '${file(params.clustersv_source).toAbsolutePath()}' --rscript '${file(params.pcawg_rscript).toAbsolutePath()}'" : 'test_clustering.py'
    def verification = pcawg ? "'${pythonExecutable}' verify_pcawg.py --group-output group_output --expected-count ${expected} --output group_output/verification.json" : ''
    """
    set -euo pipefail
    export MPLBACKEND=Agg
    export OMP_NUM_THREADS=${task.cpus}
    '${pythonExecutable}' ${control} > cluster-selftest.log 2>&1
    '${pythonExecutable}' '${runner}' --manifest '${manifest}' --cohort '${cohort}' --caller '${caller}' \
        --output group_output --expected-pairs ${expected} --cpus ${task.cpus} ${pilot} ${methodArgs}
    ${verification}
    cp cluster-selftest.log group_output/cluster-selftest.log
    """
}

workflow {
    if (params.genome_build != 'GRCh37') error 'Supported signature routes require GRCh37; other assemblies are not implemented'
    if (params.sequencing_type != 'WGS') error 'Genome-wide signature methods require WGS; capture assays need a capture-aware model'
    if (!(params.method in ['cosmic_sv32', 'pcawg_clustersv'])) error 'method must be cosmic_sv32 or pcawg_clustersv'
    if (!(params.pcawg_stage in ['clustering', 'features', 'signatures'])) error 'pcawg_stage must be clustering, features or signatures'
    if (!(params.pcawg_timing_source in ['none','reference'])) error 'pcawg_timing_source must be none or reference'
    def resolvedOutdir = file(params.outdir ?: "${projectDir}/../results/${params.method == 'cosmic_sv32' ? 'sv_signatures_corrected' : 'sv_signatures_pcawg'}").toAbsolutePath()
    def manifest = file(params.manifest, checkIfExists: true)
    def lines = manifest.readLines()
    def required = ['cohort', 'pair_id', 'caller', 'vcf', 'source_status', 'source_sha256']
    if (!lines || lines[0].trim().split('\\t', -1).toList() != required) error "Manifest header must be ${required.join(', ')}"
    def parsed = lines.drop(1).collect { it.split('\\t', -1) }
    if (!parsed || parsed.any { it.size() != 6 || !it[0] || !it[1] || !it[2] || !it[3] }) error 'Malformed six-column input manifest'
    if (parsed.collect { [it[0], it[1], it[2]] }.unique().size() != parsed.size()) error 'Duplicate cohort/pair/caller in manifest'
    if (parsed.any { !(it[0] ==~ /[A-Za-z0-9_-]+/) || !(it[1] ==~ /[A-Za-z0-9_-]+/) || !(it[2] in ['manta','delly','gridss','svaba','brass']) }) error 'Invalid sample, cohort or caller identifier'
    def allGroups = parsed.collect { "${it[0]}:${it[2]}".toString() }.unique()
    def wanted = params.groups == 'all' ? allGroups : params.groups.toString().tokenize(',')
    if (!wanted || wanted.size() != wanted.unique().size() || wanted.any { !(it in allGroups) })
        error "groups must be 'all' or unique cohort:caller names from ${allGroups.join(', ')}"
    if ((params.min_size_bp as int) < 1 || (params.group_cpus as int) < 1 ||
        (params.group_memory_gb as int) < 1 || (params.max_parallel_groups as int) < 1)
        error 'Invalid positive numeric SV32 resource or size parameter'
    if (params.pair_id && !(params.pair_id ==~ /[A-Za-z0-9_-]+/)) error 'Invalid pair_id'
    if (params.pair_id && params.pair_ids) error 'Use pair_id or pair_ids, not both'
    def ids = params.pair_id ? [params.pair_id.toString()] : params.pair_ids.toString().tokenize(',')
    if (ids.any { !(it ==~ /[A-Za-z0-9_-]+/) } || ids.unique().size() != ids.size()) error 'Invalid or duplicate pair_ids'
    def runner = file("${projectDir}/${params.method == 'cosmic_sv32' ? 'run_group.py' : 'run_clustersv.py'}", checkIfExists: true)
    def converter = file("${projectDir}/vcf_to_bedpe.py", checkIfExists: true)
    def normalizer = file("${projectDir}/normalize_exposures.py", checkIfExists: true)
    if (params.python != 'python3' && !file(params.python).isFile()) error "Pinned signature Python environment missing: ${params.python}"
    def helpers = [converter, normalizer, file("${projectDir}/../scripts/reference_policy.py", checkIfExists:true)]
    if (params.method == 'cosmic_sv32') helpers.add(file("${projectDir}/test_clustering.py", checkIfExists: true))
    else {
        helpers.addAll(['run_group.py','verify_pcawg.py','clustersv_features.py','pcawg/cluster.R','pcawg/nmf.R','pcawg/stable_nmf.R','pcawg/upstream.json','pcawg/test_clustersv.py','pcawg/replication_timing.py'].collect { file("${projectDir}/${it}", checkIfExists: true) })
        helpers.add(file("${projectDir}/../scripts/measure_process_tree.py", checkIfExists: true))
        if (!file(params.pcawg_rscript).isFile() || !file("${params.clustersv_source}/R/clustering_index.R").isFile()) error 'Missing ClusterSV resources; run pcawg setup first'
        if (params.pcawg_timing_source == 'reference' &&
            (!file("${params.pcawg_timing_reference}/reference.json").isFile() ||
             !file("${params.pcawg_timing_reference}/verification.json").isFile()))
            error 'Missing verified tertile timing reference; run pcawg/replication_run.sbatch first'
    }
    def provenance = file("${resolvedOutdir}/provenance/${workflow.runName}")
    provenance.mkdirs()
    def resolved = new LinkedHashMap(params)
    resolved.outdir = resolvedOutdir.toString()
    file("${provenance}/resolved-params.json").text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(resolved))
    manifest.copyTo(file("${provenance}/input_manifest.tsv"))
    def selected = wanted.collect { name ->
        def bits = name.split(':')
        def n = parsed.count { it[0] == bits[0] && it[2] == bits[1] && (!ids || it[1] in ids) }
        def expected = ids ? ids.size() : parsed.count { it[0] == bits[0] && it[2] == bits[1] }
        if (n != expected) error "${name}: expected ${expected} source VCFs, found ${n}"
        def target = params.pair_id ? "${resolvedOutdir}/smoke/${bits[0]}/${bits[1]}/${params.pair_id}/group_output" :
            (params.pair_ids ? "${resolvedOutdir}/pilot/${bits[0]}/${bits[1]}/group_output" : "${resolvedOutdir}/${bits[0]}/${bits[1]}/group_output")
        if (file(target).exists()) error "Output already exists: ${target}; choose a new outdir to preserve the earlier run"
        tuple(bits[0], bits[1], expected)
    }
    SIGNATURE_GROUP(Channel.fromList(selected), manifest, runner, helpers)
}
