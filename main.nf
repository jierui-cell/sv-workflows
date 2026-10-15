nextflow.enable.dsl = 2
include { CALL_SV as CALL_DELLY } from './modules/call_sv'
include { CALL_SV as CALL_GRIDSS } from './modules/call_sv'
include { CALL_SV as CALL_SVABA } from './modules/call_sv'
include { INTEGRATED_CALLER } from './modules/integrated_caller'

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
    module { ((params.modules ?: []) + (params.modules_samtools ?: '').tokenize(':')).unique().join(':') }
    publishDir "${params.outdir}/preflight/${meta.id}", mode: 'copy',
        saveAs: { name -> (name.endsWith('.bam') || name.endsWith('.bai')) ? null : name }
    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path reference_fai
    path validator
    path reference_checker
    val is_stub_run
    path exclude_regions
    path target_regions
    output:
    val meta, emit: validated
    path '*.validation.tsv'
    path 'validation.log'
    path '*.header.sam'
    path '*.samples.txt'
    path '*.interval_read_count.txt'
    path 'reference.validation.json'
    path 'input-regions/*', optional: true
    script:
    def regionCheck = (exclude_regions || target_regions) ? "python3 '${projectDir}/scripts/regions.py' ${exclude_regions ? "--exclude '${exclude_regions[0]}'" : ''} ${target_regions ? "--targets '${target_regions[0]}'" : ''} --fai '${reference_fai}' --out input-regions" : ''
    """
    set -euo pipefail
    ${regionCheck}
    python3 '${reference_checker}' --fai '${reference_fai}' --build '${params.genome_build}' \
      --output reference.validation.json ${is_stub_run ? '--allow-synthetic' : ''}
    bash '${validator}' '${meta.id}' '${reference_fai}' '${tumor_bam}' '${normal_bam}' '${params.allow_local_execution ? '1' : '0'}' 2>&1 | tee validation.log
    """
}

process MANTA {
    tag "${meta.id}"
    cpus { params.manta_cpus as int }
    memory { "${params.manta_memory_gb} GB" }
    time { params.manta_time }
    maxForks (params.max_parallel_callers as int)
    scratch false
    module { ((params.modules ?: []) + (params.modules_bcftools ?: '').tokenize(':')).unique().join(':') }
    publishDir "${params.outdir}/manta/${meta.id}", mode: 'copy',
        saveAs: { name -> params.compact_publish ? name.tokenize('/').last() : name }
    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path reference
    path reference_fai
    path image
    path configurator
    path benchmarker
    path exclude_regions
    val exclude_sha256
    path target_regions
    val target_sha256
    output:
    path 'manta/results/variants/*', optional: true
    path 'somaticSV.stats.txt'
    path 'artifacts.sha256'
    path 'manta.*.ini'
    path 'manta.overrides.json'
    path 'manta.time.txt', optional: true
    path 'output/benchmark.json'
    path 'regions/*', optional: true
    script:
    def tumorDir = meta.tumor_dir
    def normalDir = meta.normal_dir
    def refDir = file(params.reference).toRealPath().parent
    def options = groovy.json.JsonOutput.toJson(params.manta_options)
    def regionPreparation = (exclude_regions || target_regions) ? "python3 '${projectDir}/scripts/regions.py' ${exclude_regions ? "--exclude '${exclude_regions[0]}'" : ''} ${target_regions ? "--targets '${target_regions[0]}'" : ''} --fai '${reference_fai}' --out regions --manta" : ''
    def regionFlags = (exclude_regions || target_regions) ? "--callRegions regions/manta.call-regions.bed.gz" : ''
    def assayFlags = params.sequencing_type == 'WGS' ? '' : '--exome'
    def binds = "${projectDir}:${projectDir},\$PWD:\$PWD,${refDir}:${refDir}:ro,${tumorDir}:${tumorDir}:ro,${normalDir}:${normalDir}:ro"
    def benchmarkSpec = groovy.json.JsonOutput.toJson([caller:'manta', pair_id:meta.id,
        dataset_id:params.dataset_id, reference:file(params.reference).toString(),
        allocated_cpus:task.cpus, requested_memory_gb:task.memory.toGiga(),
        benchmark_scope:'caller configuration through final output QC', output_class:'somatic',
        settings:[options:params.manta_options, image:params.manta_image,
                  workflow_memory_gb:params.manta_workflow_memory_gb,
                  cpus:task.cpus, memory_gb:task.memory.toGiga(),
                  genome_build:params.genome_build, exclude_bed:params.exclude_bed,
                  exclude_bed_sha256:exclude_sha256, sequencing_type:params.sequencing_type,
                  target_bed:params.target_bed, target_bed_sha256:target_sha256]])
    """
    cat > benchmark.meta.json <<'BENCHMARK_SETTINGS'
    ${benchmarkSpec}
    BENCHMARK_SETTINGS
    cat > manta-task.sh <<'MANTA_TASK'
    set -euo pipefail
    ${regionPreparation}
    sha256sum '${image}' '${reference}' '${reference_fai}' > artifacts.sha256
    cat > manta.overrides.json <<'MANTA_OPTIONS'
    ${options}
    MANTA_OPTIONS
    singularity exec --cleanenv --bind "${binds}" '${image}' python '${configurator}' manta.overrides.json --reference '${reference}'
    singularity exec --cleanenv --bind "${binds}" '${image}' \
      configManta.py --normalBam '${normal_bam}' --tumorBam '${tumor_bam}' \
      --referenceFasta '${reference}' --config manta.effective.ini --runDir manta ${assayFlags} ${regionFlags}
    python3 '${projectDir}/scripts/measure_command.py' manta.time.txt singularity exec --cleanenv --bind "${binds}" '${image}' \
      python manta/runWorkflow.py -m local -j ${task.cpus} -g ${params.manta_workflow_memory_gb}
    bcftools view -Ov -o /dev/null manta/results/variants/somaticSV.vcf.gz
    bcftools stats manta/results/variants/somaticSV.vcf.gz > somaticSV.stats.txt
    MANTA_TASK
    python3 '${benchmarker}' output/benchmark.json benchmark.meta.json -- bash manta-task.sh
    """
    stub:
    def options = groovy.json.JsonOutput.toJson(params.manta_options)
    def regionPreparation = (exclude_regions || target_regions) ? "python3 '${projectDir}/scripts/regions.py' ${exclude_regions ? "--exclude '${exclude_regions[0]}'" : ''} ${target_regions ? "--targets '${target_regions[0]}'" : ''} --fai '${reference_fai}' --out regions --manta" : ''
    def regionFlags = (exclude_regions || target_regions) ? '--callRegions regions/manta.call-regions.bed.gz' : ''
    def assayFlags = params.sequencing_type == 'WGS' ? '' : '--exome'
    def binds = "${projectDir}:${projectDir},\$PWD:\$PWD,${file(params.reference).toRealPath().parent}:${file(params.reference).toRealPath().parent}:ro,${meta.tumor_dir}:${meta.tumor_dir}:ro,${meta.normal_dir}:${meta.normal_dir}:ro"
    """
    ${regionPreparation}
    cat > manta.overrides.json <<'MANTA_OPTIONS'
    ${options}
    MANTA_OPTIONS
    singularity exec --cleanenv --bind "${binds}" '${image}' python '${configurator}' manta.overrides.json --reference '${reference}'
    singularity exec --cleanenv --bind "${binds}" '${image}' configManta.py \
      --normalBam '${normal_bam}' --tumorBam '${tumor_bam}' --referenceFasta '${reference}' \
      --config manta.effective.ini --runDir manta ${assayFlags} ${regionFlags}
    sha256sum manta.effective.ini > artifacts.sha256
    echo 'STUB: configuration tested; no SV calls produced' > somaticSV.stats.txt
    mkdir -p output
    echo '{"status":"stub_not_benchmarked"}' > output/benchmark.json
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
    if (!(params.sequencing_type in ['WGS','WES','targeted'])) error 'sequencing_type must be WGS, WES or targeted'
    if (params.sequencing_type != 'WGS' && requested.any { !(it in ['preflight','manta']) })
        error 'Capture-assay settings are implemented for Manta; select callers=manta'
    if (params.sequencing_type == 'targeted' && !params.target_bed) error 'targeted sequencing requires target_bed'
    if (params.target_bed && (params.sequencing_type == 'WGS' || requested.any { !(it in ['preflight','manta']) }))
        error 'target_bed is implemented for WES/targeted Manta runs'
    def targetRegions = params.target_bed ? [file(params.target_bed, checkIfExists:true).toRealPath()] : []
    def targetHash = params.target_bed ? java.security.MessageDigest.getInstance('SHA-256')
        .digest(file(params.target_bed).bytes).encodeHex().toString() : 'none'
    if (!(params.genome_build in ['GRCh37','GRCh38']) && !(workflow.stubRun && params.genome_build == 'synthetic'))
        error 'genome_build must be GRCh37 or GRCh38; synthetic is permitted only with -stub-run'
    if (params.exclude_bed && requested.any { !(it in ['preflight','manta','delly','gridss']) })
        error 'exclude_bed is implemented for Manta, DELLY and GRIDSS only; run other callers separately without this setting'
    def excludeRegions = params.exclude_bed ? [file(params.exclude_bed, checkIfExists:true).toRealPath()] : []
    def excludeHash = params.exclude_bed ? java.security.MessageDigest.getInstance('SHA-256')
        .digest(file(params.exclude_bed).bytes).encodeHex().toString() : 'none'
    if (requested.any { !(it in ['preflight', 'manta', 'delly', 'svaba', 'gridss', 'brass']) }) {
        error 'Supported callers: manta, delly, svaba, gridss, brass.'
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
        .map { original ->
            def row = new LinkedHashMap(original)
            ['tumor_bam','tumor_bai','normal_bam','normal_bai'].each { key ->
                row[key] = file(row[key], checkIfExists: true).toRealPath().toString()
            }
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
                  file("${projectDir}/scripts/validate_pair.sh", checkIfExists: true),
                  file("${projectDir}/scripts/reference_policy.py", checkIfExists: true), workflow.stubRun, excludeRegions, targetRegions)
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
            file("${projectDir}/scripts/measure_benchmark.py", checkIfExists: true), excludeRegions, excludeHash, targetRegions, targetHash)
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
            if (caller == 'delly') CALL_DELLY(jobs, referenceFiles, runner, benchmarker, excludeRegions, excludeHash)
            if (caller == 'gridss') CALL_GRIDSS(jobs, referenceFiles, runner, benchmarker, excludeRegions, excludeHash)
            if (caller == 'svaba') CALL_SVABA(jobs, referenceFiles, runner, benchmarker, excludeRegions, excludeHash)
        }
    }
    def integrated = requested.findAll { it == 'brass' }
    if ('brass' in requested && params.brass.assembly != 'GRCh37')
        error 'The current validated BRASS bundle adapter supports GRCh37 only.'
    if (integrated) {
        if ((params.max_parallel_callers as int) < 1) error 'max_parallel_callers must be positive'
        integrated.each { caller ->
            def cfg = params[caller]
            if (!(cfg instanceof Map)) error "Missing YAML configuration for ${caller}"
            if ((cfg.cpus as int) < 1 || (cfg.memory_gb as int) < 1 || !cfg.time)
                error "Invalid compute resources for ${caller}"
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
}
