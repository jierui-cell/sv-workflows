nextflow.enable.dsl = 2

workflow {
    def header = Samplesheet.REQUIRED.join('\t')
    def valid = 'TEST\t/tumor.bam\t/tumor.bam.bai\t/normal.bam\t/normal.bam.bai\tfemale\tmelanoma\t0.8\t2.0'
    assert Samplesheet.parse([header, valid]).size() == 1
    def absent = valid.replace('female', '').replace('0.8\t2.0', '\t')
    assert Samplesheet.parse([header, absent])[0].missing_metadata == ['sex', 'purity', 'ploidy']
    def cases = [
        [[header.replace('tumor_bam', 'tumour_bam'), valid], 'tumor_bam'],
        [[header + '\tsex', valid + '\tfemale'], 'duplicate'],
        [[header, valid.replace('/tumor.bam\t', '\t')], 'tumor_bam'],
        [[header, valid, valid], 'duplicate patient_id'],
        [[header, valid.replace('TEST', '../unsafe')], 'unsafe patient_id'],
        [[header, valid.replace('0.8', '1.2')], 'purity'],
        [[header, valid.replace('2.0', 'NaN')], 'ploidy'],
        [[header, valid.replace('/tumor.bam\t', '/normal.bam\t')], 'same file'],
        [[header, valid + '\textra'], 'expected 9 columns'],
        [[header], 'no data rows']
    ]
    cases.each { test ->
        def caught = false
        try { Samplesheet.parse(test[0]) }
        catch (IllegalArgumentException e) {
            assert e.message.contains(test[1]): e.message
            caught = true
        }
        assert caught: "Expected samplesheet rejection: ${test[1]}"
    }
    log.info 'PASS: samplesheet valid/optional metadata and 10 negative cases'
}
