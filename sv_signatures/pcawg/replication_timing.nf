nextflow.enable.dsl = 2

params.account = null
params.partition = null
params.timing_python = "${projectDir}/../../environments/replication-timing-python/bin/python"
params.timing_outdir = "${projectDir}/../../resources/replication_timing/pcawg_mean_tertiles"
params.timing_label_policy = 'tertiles'
params.timing_input_reference = ''
params.timing_bin_bp = 1000
params.timing_min_coverage = 1.0
params.timing_chunk_bins = 5000
params.timing_cpus = 1
params.timing_memory_gb = 8
params.timing_time = '1h'
params.timing_verify_only = false

process DOWNLOAD_TIMING {
    tag { cell_line }
    cpus 1
    memory '1 GB'
    time '30m'
    publishDir "${file(params.timing_outdir).toAbsolutePath()}/source", mode: 'copy', overwrite: false
    input:
    val cell_line
    path source_lock
    path helper
    output:
    path '*.bigWig'
    script:
    """
    set -euo pipefail
    '${file(params.timing_python).toAbsolutePath()}' '${helper}' download \
        --sources '${source_lock}' --cell-line '${cell_line}'
    """
}

process BUILD_TIMING {
    tag 'PCAWG_three_cell_line_mean'
    cpus { params.timing_cpus as int }
    memory { "${params.timing_memory_gb} GB" }
    time { params.timing_time }
    publishDir "${file(params.timing_outdir).toAbsolutePath()}", mode: 'copy', overwrite: false,
        saveAs: { name -> name.tokenize('/').last() }
    input:
    path tracks
    path source_lock
    path sizes
    path helper
    path tests
    path memory_monitor
    output:
    path 'reference/*'
    path 'build.memory.json'
    path 'tests.log'
    script:
    """
    set -euo pipefail
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
    '${file(params.timing_python).toAbsolutePath()}' '${tests}' > tests.log 2>&1
    '${file(params.timing_python).toAbsolutePath()}' '${memory_monitor}' build.memory.json \
        '${file(params.timing_python).toAbsolutePath()}' '${helper}' build \
        --sources '${source_lock}' --chrom-sizes '${sizes}' --tracks ${tracks} \
        --outdir reference --label-policy '${params.timing_label_policy}' \
        --bin-bp ${params.timing_bin_bp} --chunk-bins ${params.timing_chunk_bins} \
        --min-coverage ${params.timing_min_coverage}
    """
}

process VERIFY_TIMING {
    tag 'independent_reference_checks'
    cpus 1
    memory '2 GB'
    time '10m'
    publishDir "${file(params.timing_outdir).toAbsolutePath()}", mode: 'copy', overwrite: true
    input:
    path reference_files, stageAs: 'reference/*'
    path tracks, stageAs: 'reference/source/*'
    path helper
    path verifier
    output:
    path 'verification.json'
    script:
    """
    set -euo pipefail
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
    '${file(params.timing_python).toAbsolutePath()}' '${verifier}' \
        --reference reference --output verification.json
    """
}

process RELABEL_TIMING {
    tag 'reuse_verified_continuous_reference'
    cpus 1
    memory '2 GB'
    time '10m'
    publishDir "${file(params.timing_outdir).toAbsolutePath()}", mode: 'copy', overwrite: false,
        saveAs: { name -> name.replaceFirst('^reference/', '') }
    input:
    path base_reference
    path source_lock
    path helper
    path tests
    path monitor
    output:
    path 'reference/*'
    path 'label.memory.json'
    path 'tests.log'
    script:
    """
    set -euo pipefail
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
    '${file(params.timing_python).toAbsolutePath()}' '${tests}' > tests.log 2>&1
    '${file(params.timing_python).toAbsolutePath()}' '${monitor}' label.memory.json \
        '${file(params.timing_python).toAbsolutePath()}' '${helper}' relabel \
        --reference '${base_reference}' --outdir reference --label-policy '${params.timing_label_policy}'
    """
}

workflow {
    if (!(params.timing_label_policy in ['continuous','panprostate_quartiles','tertiles']))
        error 'timing_label_policy must be continuous, panprostate_quartiles or tertiles'
    if (!file(params.timing_python).isFile()) error 'Missing isolated timing Python environment; use replication_run.sbatch'
    def root = projectDir.parent.parent
    def lock = file("${projectDir}/replication_sources.json", checkIfExists: true)
    def helper = file("${projectDir}/replication_timing.py", checkIfExists: true)
    def sizes = file("${root}/resources/clustersv/source/references/hs37d5.chrom_sizes", checkIfExists: true)
    def tests = file("${projectDir}/test_replication_timing.py", checkIfExists: true)
    def monitor = file("${root}/scripts/measure_process_tree.py", checkIfExists: true)
    def verifier = file("${projectDir}/verify_replication_reference.py", checkIfExists: true)
    def outputRoot = file(params.timing_outdir).toAbsolutePath()
    def existing = file("${outputRoot}/reference.json")
    if (existing.exists()) {
        def settings = new groovy.json.JsonSlurper().parseText(existing.text)
        if (settings.label_policy != params.timing_label_policy || settings.bin_bp != (params.timing_bin_bp as int) ||
            settings.min_input_coverage != (params.timing_min_coverage as double))
            error 'Existing reference settings differ; choose a new timing_outdir to avoid mixing outputs'
        if (!workflow.resume && !params.timing_verify_only)
            error 'Reference already exists; use -resume with unchanged settings or timing_verify_only'
    }
    if (params.timing_verify_only) {
        def products = ['reference.json','replication_timing.npz','replication_timing.bw',
                        'replication_timing.bed.tsv.gz','chromosome_qc.tsv'].collect {
            file("${outputRoot}/${it}", checkIfExists: true)
        }
        def sourceInfo = new groovy.json.JsonSlurper().parseText(lock.text)
        def tracks = sourceInfo.sources.collect { file("${outputRoot}/source/${it.filename}", checkIfExists: true) }
        VERIFY_TIMING(Channel.value(products), Channel.value(tracks), helper, verifier)
    } else {
        def inputReference = params.timing_input_reference ? file(params.timing_input_reference,checkIfExists: true) :
            file("${root}/resources/replication_timing/pcawg_three_cell_line_mean")
        def canReuse = false
        if (inputReference.isDirectory() && file("${inputReference}/reference.json").exists()) {
            def base = new groovy.json.JsonSlurper().parseText(file("${inputReference}/reference.json").text)
            canReuse = base.bin_bp == (params.timing_bin_bp as int) &&
                base.min_input_coverage == (params.timing_min_coverage as double)
            if (!canReuse && params.timing_input_reference)
                error 'Prepared input reference has different bin/coverage settings; leave timing_input_reference empty to rebuild'
        }
        if (canReuse) {
            RELABEL_TIMING(inputReference,lock,helper,tests,monitor)
            // The reference directory is already populated with its source tracks.
            def products = RELABEL_TIMING.out[0].map { files -> files.findAll { it.fileName.toString() != 'source' } }
            def sourceInfo = new groovy.json.JsonSlurper().parseText(lock.text)
            def tracks = RELABEL_TIMING.out[0].map { files ->
                def dir = files.find { it.fileName.toString() == 'source' }
                sourceInfo.sources.collect { dir.resolve(it.filename) }
            }
            VERIFY_TIMING(products, tracks, helper, verifier)
        } else {
            if (params.timing_input_reference) error 'timing_input_reference is not a prepared reference directory'
            def downloads = DOWNLOAD_TIMING(Channel.of('NHEK','GM12878','IMR90'), lock, helper).collect()
            BUILD_TIMING(downloads, lock, sizes, helper, tests, monitor)
            VERIFY_TIMING(BUILD_TIMING.out[0], downloads, helper, verifier)
        }
    }
}
