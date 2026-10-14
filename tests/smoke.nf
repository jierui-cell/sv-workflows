nextflow.enable.dsl = 2
params.outdir = 'results/toy'
process SMOKE {
    publishDir params.outdir, mode: 'copy'
    cpus 1
    memory '1 GB'
    input:
    path converter
    path fixture
    output:
    path 'toy.bedpe.tsv'
    path 'conversion.json'
    path 'audit.tsv'
    path 'PASS.txt'
    script:
    """
    python3 '${converter}' --vcf '${fixture}' --caller manta --sample TOY_PAIR --outdir toy-converted
    cp toy-converted/TOY_PAIR.bedpe.tsv toy.bedpe.tsv
    cp toy-converted/event_audit.tsv audit.tsv
    cp toy-converted/summary.json conversion.json
    python3 -c 'import json; assert json.load(open("conversion.json"))["retained_events"] == 3'
    echo 'PASS: synthetic VCF conversion and Nextflow task execution; not a caller accuracy test' > PASS.txt
    """
}
workflow {
    SMOKE(file("${projectDir}/../sv_signatures/vcf_to_bedpe.py"), file("${projectDir}/../examples/toy.somatic.vcf"))
}
