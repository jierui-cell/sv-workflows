nextflow.enable.dsl = 2
include { CALL_SV as CALL_DELLY } from './modules/call_sv'
include { CALL_SV as CALL_GRIDSS } from './modules/call_sv'
include { CALL_SV as CALL_SVABA } from './modules/call_sv'
include { INTEGRATED_CALLER } from './modules/integrated_caller'
include { CUE_PAIR } from './modules/cue_pair'
include { JABBA } from './modules/jabba'
include { JABBA_COVERAGE } from './modules/jabba_coverage'

params.input = "${projectDir}/examples/samplesheet.tsv"
params.reference = ''
params.callers = 'preflight'
params.pair_ids = ''
params.outdir = "${projectDir}/results"

process VALIDATE_PAIR {
    tag "${meta.id}"
    cpus 1
    memory '2 GB'
    time '30m'
    maxForks (params.max_parallel_preflight as int)
    module { params.modules_samtools }
    publishDir "${params.outdir}/preflight/${meta.id}", mode: 'copy',
        saveAs: { name -> (name.endsWith('.bam') || name.endsWith('.bai')) ? null : name }
    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path reference_fai
    path validator
    output:
    val meta, emit: validated
    path '*.validation.tsv'
    path 'validation.log'
    path '*.header.sam'
    path '*.samples.txt'
    path '*.interval_read_count.txt'
    script:
    """
    set -euo pipefail
    bash '${validator}' '${meta.id}' '${reference_fai}' '${tumor_bam}' '${normal_bam}' 2>&1 | tee validation.log
    """
}

process MANTA {
    tag "${meta.id}"
    cpus { params.manta_cpus as int }
    memory { "${params.manta_memory_gb} GB" }
    time { params.manta_time }
    maxForks (params.max_parallel_callers as int)
    scratch false
    module { params.modules_bcftools }
    publishDir "${params.outdir}/manta/${meta.id}", mode: 'copy',
        saveAs: { name -> params.compact_publish ? name.tokenize('/').last() : name }
    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path reference
    path reference_fai
    path image
    path configurator
    path benchmarker
    output:
    path 'manta/results/variants/*', optional: true
    path 'somaticSV.stats.txt'
    path 'artifacts.sha256'
    path 'manta.*.ini'
    path 'manta.overrides.json'
    path 'manta.time.txt', optional: true
    path 'benchmark.json'
    script:
    def tumorDir = meta.tumor_dir
    def normalDir = meta.normal_dir
    def refDir = file(params.reference).parent
    def options = groovy.json.JsonOutput.toJson(params.manta_options)
    def benchmarkSpec = groovy.json.JsonOutput.toJson([caller:'manta', pair_id:meta.id,
        dataset_id:params.dataset_id, reference:file(params.reference).toString(),
        allocated_cpus:task.cpus, requested_memory_gb:task.memory.toGiga(),
        benchmark_scope:'caller configuration through final output QC', output_class:'somatic',
        settings:[options:params.manta_options, image:params.manta_image,
                  workflow_memory_gb:params.manta_workflow_memory_gb,
                  cpus:task.cpus, memory_gb:task.memory.toGiga()]])
    """
    cat > benchmark.meta.json <<'BENCHMARK_SETTINGS'
    ${benchmarkSpec}
    BENCHMARK_SETTINGS
    cat > manta-task.sh <<'MANTA_TASK'
    set -euo pipefail
    sha256sum '${image}' '${reference}' '${reference_fai}' > artifacts.sha256
    cat > manta.overrides.json <<'MANTA_OPTIONS'
    ${options}
    MANTA_OPTIONS
    singularity exec --cleanenv --bind '${projectDir}:${projectDir}' '${image}' python '${configurator}' manta.overrides.json --reference '${reference}'
    singularity exec --cleanenv --bind '${projectDir}:${projectDir},${refDir}:${refDir}:ro,${tumorDir}:${tumorDir}:ro,${normalDir}:${normalDir}:ro' '${image}' \
      configManta.py --normalBam '${normal_bam}' --tumorBam '${tumor_bam}' \
      --referenceFasta '${reference}' --config manta.effective.ini --runDir manta
    python3 '${projectDir}/scripts/measure_command.py' manta.time.txt singularity exec --cleanenv --bind '${projectDir}:${projectDir},${refDir}:${refDir}:ro,${tumorDir}:${tumorDir}:ro,${normalDir}:${normalDir}:ro' '${image}' \
      python manta/runWorkflow.py -m local -j ${task.cpus} -g ${params.manta_workflow_memory_gb}
    bcftools view -Ov -o /dev/null manta/results/variants/somaticSV.vcf.gz
    bcftools stats manta/results/variants/somaticSV.vcf.gz > somaticSV.stats.txt
    MANTA_TASK
    python3 '${benchmarker}' benchmark.json benchmark.meta.json -- bash manta-task.sh
    """
    stub:
    def options = groovy.json.JsonOutput.toJson(params.manta_options)
    """
    cat > manta.overrides.json <<'MANTA_OPTIONS'
    ${options}
    MANTA_OPTIONS
    singularity exec --cleanenv --bind '${projectDir}:${projectDir}' '${image}' python '${configurator}' manta.overrides.json --reference '${reference}'
    sha256sum manta.effective.ini > artifacts.sha256
    echo 'STUB: configuration tested; no SV calls produced' > somaticSV.stats.txt
    echo '{"status":"stub_not_benchmarked"}' > benchmark.json
    """
}

workflow {
    if (!params.reference) error 'Set reference to the exact FASTA used for BAM alignment.'
    if ((params.manta_cpus as int) < 1 || (params.manta_workflow_memory_gb as int) < 1 ||
        (params.manta_workflow_memory_gb as int) >= (params.manta_memory_gb as int))
        error 'Manta requires positive CPUs and workflow memory below allocated memory.'
    def provenance = file("${params.outdir}/provenance/${workflow.runName}")
    provenance.mkdirs()
    file("${provenance}/resolved-params.json").text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params))
    def requested = params.callers.tokenize(',')
    if ('lumpy' in requested) error 'LUMPY is discontinued and not supported in this release.'
    if (requested.any { !(it in ['preflight', 'manta', 'delly', 'svaba', 'gridss', 'brass', 'lumpy', 'cue', 'jabba']) }) {
        error 'Implemented Nextflow callers: manta, delly, svaba, gridss, brass, lumpy, cue, jabba.'
    }
    if (requested.any { it != 'preflight' } && !params.dataset_id)
        error 'Set dataset_id in YAML so benchmark results from different datasets cannot be mixed.'
    if (!file(params.input).exists()) {
        error "Samplesheet does not exist: ${params.input}"
    }
    file(params.input).copyTo(file("${provenance}/samplesheet.tsv"))
    // Header and every row are checked synchronously before any job channel.
    def sheetRows = Samplesheet.read(file(params.input))
    sheetRows.findAll { (!params.pair_ids || it.patient_id in params.pair_ids.tokenize(',')) && it.missing_metadata }.each { row ->
        log.warn "Samplesheet ${row.patient_id}: unavailable metadata: ${row.missing_metadata.join(', ')}; values will not be inferred."
    }
    def unknownIds = params.pair_ids ? params.pair_ids.tokenize(',') - sheetRows.collect { it.patient_id } : []
    if (unknownIds) error "Unknown pair_ids: ${unknownIds.join(', ')}"

    pair_ch = Channel.fromList(sheetRows)
        .filter { row -> !params.pair_ids || row.patient_id in params.pair_ids.tokenize(',') }
        .ifEmpty { error 'No pairs matched the requested samplesheet/pair_ids selection.' }
        .map { row ->
            def meta = [
                id: row.patient_id,
                tumor_dir: file(row.tumor_bam).parent.toString(),
                normal_dir: file(row.normal_bam).parent.toString(),
                tumor_bam: row.tumor_bam,
                tumor_bai: row.tumor_bai,
                normal_bam: row.normal_bam,
                normal_bai: row.normal_bai,
                sex: row.sex,
                phenotype: row.tumor_phenotype,
                purity: row.purity,
                ploidy: row.ploidy,
                missing_metadata: row.missing_metadata
            ]
            tuple(
                meta,
                file(row.tumor_bam, checkIfExists: true),
                file(row.tumor_bai, checkIfExists: true),
                file(row.normal_bam, checkIfExists: true),
                file(row.normal_bai, checkIfExists: true)
            )
        }

    VALIDATE_PAIR(pair_ch, file("${params.reference}.fai", checkIfExists: true),
                  file("${projectDir}/scripts/validate_pair.sh", checkIfExists: true))
    if ('manta' in requested) {
        ready_ch = VALIDATE_PAIR.out.validated.map { meta ->
            tuple(meta, file(meta.tumor_bam), file(meta.tumor_bai),
                  file(meta.normal_bam), file(meta.normal_bai))
        }
        MANTA(ready_ch,
            file(params.reference, checkIfExists: true),
            file("${params.reference}.fai", checkIfExists: true),
            file(params.manta_image, checkIfExists: true),
            file("${projectDir}/scripts/configure_manta.py", checkIfExists: true),
            file("${projectDir}/scripts/measure_benchmark.py", checkIfExists: true))
    } else if (requested == ['preflight']) {
        log.info 'PREFLIGHT ONLY: no SV callers requested.'
    }
    def expanded = requested.findAll { it in ['delly', 'svaba', 'gridss'] }
    if (expanded) {
        if ((params.max_parallel_callers as int) < 1) error 'max_parallel_callers must be positive'
        expanded.each { caller ->
            if (!(params[caller] instanceof Map)) error "Missing YAML configuration for ${caller}"
        }
        def suffixes = ['', '.fai']
        if (expanded.any { it in ['svaba', 'gridss'] }) suffixes += ['.amb', '.ann', '.bwt', '.pac', '.sa']
        expanded.each { caller ->
            def jobs = VALIDATE_PAIR.out.validated.map { meta ->
                tuple(meta, caller, params[caller], file(params[caller].image, checkIfExists:true),
                    file(meta.tumor_bam), file(meta.tumor_bai), file(meta.normal_bam), file(meta.normal_bai))
            }
            def referenceFiles = suffixes.collect { file("${params.reference}${it}", checkIfExists:true) }
            if (caller == 'gridss')
                referenceFiles.add(file("${file(params.reference).parent}/${file(params.reference).baseName}.dict", checkIfExists:true))
            def runner = file("${projectDir}/scripts/run_caller.py", checkIfExists:true)
            def benchmarker = file("${projectDir}/scripts/measure_benchmark.py", checkIfExists:true)
            if (caller == 'delly') CALL_DELLY(jobs, referenceFiles, runner, benchmarker)
            if (caller == 'gridss') CALL_GRIDSS(jobs, referenceFiles, runner, benchmarker)
            if (caller == 'svaba') CALL_SVABA(jobs, referenceFiles, runner, benchmarker)
        }
    }
    def integrated = requested.findAll { it in ['brass', 'lumpy'] }
    if ('brass' in requested && params.brass.assembly != 'GRCh37')
        error 'The current validated BRASS bundle adapter supports GRCh37 only.'
    if (integrated) {
        if ((params.max_parallel_callers as int) < 1) error 'max_parallel_callers must be positive'
        integrated.each { caller ->
            def cfg = params[caller]
            if (!(cfg instanceof Map)) error "Missing YAML configuration for ${caller}"
            if ((cfg.cpus as int) < 1 || (cfg.memory_gb as int) < 1 || !cfg.time)
                error "Invalid compute resources for ${caller}"
            if (caller == 'lumpy' && ((cfg.genotype_workers as int) < 1 ||
                    (cfg.genotype_workers as int) > (cfg.cpus as int)))
                error 'LUMPY genotype_workers must be between 1 and allocated CPUs.'
        }
        if ('lumpy' in integrated) {
            def allowed = params.lumpy.pilot_pair_ids ?: []
            def selected = sheetRows.findAll { !params.pair_ids || it.patient_id in params.pair_ids.tokenize(',') }
                    .collect { it.patient_id }
            def outside = selected - allowed
            if (outside) error "LUMPY integrated path is pilot-only; pairs outside lumpy.pilot_pair_ids: ${outside.join(', ')}"
        }
        def manifest = file(params.input, checkIfExists:true)
        integratedJobs = VALIDATE_PAIR.out.validated.flatMap { meta ->
            integrated.collect { caller -> tuple(meta, caller, params[caller] + [reference: file(params.reference).toAbsolutePath().toString()], manifest) }
        }
        INTEGRATED_CALLER(integratedJobs,
            file("${projectDir}/scripts/run_integrated_caller_v2.sh", checkIfExists:true),
            file("${projectDir}/scripts/measure_command.py", checkIfExists:true),
            file("${projectDir}/scripts/measure_process_tree.py", checkIfExists:true),
            file("${projectDir}/scripts/measure_benchmark.py", checkIfExists:true))
    }
    if ('cue' in requested) {
        if (!(params.cue instanceof Map)) error 'Missing YAML configuration for Cue'
        if ((params.cue.cpus as int) < 1 || (params.cue.memory_gb as int) < 1 || !params.cue.time)
            error 'Invalid Cue compute resources'
        if (params.cue.pair_execution != 'sequential')
            error 'Cue benchmark currently supports pair_execution: sequential only.'
        if (!file("${projectDir}/resources/cue/READY.sha256").exists())
            error 'Cue resources are not ready: prepare Cue first.'
        def cueJobs = VALIDATE_PAIR.out.validated.map { meta ->
            tuple(meta, params.cue, file(meta.tumor_bam), file(meta.tumor_bai),
                  file(meta.normal_bam), file(meta.normal_bai))
        }
        CUE_PAIR(cueJobs,
            file("${params.reference}.fai", checkIfExists:true),
            file("${projectDir}/scripts/run_cue_pair_task.py", checkIfExists:true),
            file("${projectDir}/scripts/measure_benchmark.py", checkIfExists:true))
    }
    if ('jabba' in requested) {
        if (!(params.jabba instanceof Map) || !params.jabba.inputs || !params.jabba.image)
            error 'JaBbA requires jabba.inputs and jabba.image in YAML.'
        if ((params.jabba.cpus as int) < 1 || (params.jabba.memory_gb as int) < 2 || !params.jabba.time)
            error 'Invalid JaBbA compute resources.'
        if (params.jabba.upstream_junction_caller != 'svaba' || !params.jabba.coverage_method)
            error 'JaBbA requires recorded SvABA junction source and coverage method.'
        def coverageMode = params.jabba.coverage_mode ?: 'provided'
        if (!(coverageMode in ['provided', 'generate']))
            error 'jabba.coverage_mode must be provided or generate.'
        def selected = sheetRows.findAll { !params.pair_ids || it.patient_id in params.pair_ids.tokenize(',') }
                    .collect { it.patient_id }
        def jabbaInputs = JabbaInputs.read(file(params.jabba.inputs, checkIfExists:true), projectDir,
                                           selected as Set, coverageMode == 'provided')
        def missing = selected - jabbaInputs.keySet()
        if (missing) error "JaBbA input TSV has no rows for: ${missing.join(', ')}"
        if (coverageMode == 'generate') {
            if (!params.jabba.fragcounter_image || !params.jabba.cbs_image || !params.jabba.gcmap_dir ||
                !params.jabba.cbs_source_commit || (params.jabba.coverage_cpus as int) < 1 ||
                (params.jabba.coverage_memory_gb as int) < 2 || !params.jabba.coverage_time ||
                (params.jabba.max_parallel_coverage as int) < 1)
                error 'JaBbA coverage generation requires images, GC/map tracks, pinned CBS source and resources.'
            def coverageJobs = VALIDATE_PAIR.out.validated.map { meta ->
                tuple(meta, file(jabbaInputs[meta.id].junction_vcf, checkIfExists:true),
                      file(meta.tumor_bam), file(meta.tumor_bai),
                      file(meta.normal_bam), file(meta.normal_bai))
            }
            JABBA_COVERAGE(coverageJobs,
                file(params.jabba.gcmap_dir, checkIfExists:true),
                file(params.jabba.fragcounter_image, checkIfExists:true),
                file(params.jabba.cbs_image, checkIfExists:true),
                file("${projectDir}/resources/jabba/nf-gos/bin/cbsFH_ZC.R", checkIfExists:true),
                file("${projectDir}/scripts/run_jabba_coverage.py", checkIfExists:true),
                file("${projectDir}/scripts/measure_benchmark.py", checkIfExists:true))
            JABBA(JABBA_COVERAGE.out.prepared,
                file(params.jabba.image, checkIfExists:true),
                file("${projectDir}/scripts/run_jabba.py", checkIfExists:true),
                file("${projectDir}/scripts/measure_benchmark.py", checkIfExists:true))
        } else {
            def jabbaJobs = VALIDATE_PAIR.out.validated.map { meta ->
                def row = jabbaInputs[meta.id]
                tuple(meta, file(row.junction_vcf, checkIfExists:true),
                      file(row.coverage, checkIfExists:true))
            }
            JABBA(jabbaJobs,
                file(params.jabba.image, checkIfExists:true),
                file("${projectDir}/scripts/run_jabba.py", checkIfExists:true),
                file("${projectDir}/scripts/measure_benchmark.py", checkIfExists:true))
        }
    }
}
