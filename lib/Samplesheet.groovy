/** Strict local, unquoted TSV schema. Validate before constructing any channel. */
class Samplesheet {
    static final List REQUIRED = ['patient_id', 'tumor_bam', 'tumor_bai', 'normal_bam',
                                 'normal_bai', 'sex', 'tumor_phenotype', 'purity', 'ploidy']
    static final List NONEMPTY = REQUIRED.take(5)

    static List<Map> read(java.nio.file.Path path) {
        parse(java.nio.file.Files.readAllLines(path))
    }

    static List<Map> parse(List<String> lines) {
        if (!lines) throw new IllegalArgumentException('Samplesheet is empty')
        def header = lines[0].replaceFirst('^\uFEFF', '').split('\t', -1).toList()
        if (header.any { !it || it != it.trim() } || header.toSet().size() != header.size())
            throw new IllegalArgumentException('Samplesheet has empty, whitespace-padded or duplicate column names')
        def missing = REQUIRED - header
        if (missing) throw new IllegalArgumentException("Samplesheet missing required columns: ${missing.join(', ')}")
        def seen = [] as Set
        def rows = []
        lines.drop(1).eachWithIndex { line, index ->
            def lineNo = index + 2
            if (!line.trim()) throw new IllegalArgumentException("Samplesheet line ${lineNo}: blank row")
            def cells = line.split('\t', -1).toList()
            if (cells.size() != header.size())
                throw new IllegalArgumentException("Samplesheet line ${lineNo}: expected ${header.size()} columns, found ${cells.size()}")
            def row = [header, cells].transpose().collectEntries { [(it[0]): it[1]] }
            NONEMPTY.each { key ->
                if (!row[key]?.trim()) throw new IllegalArgumentException("Samplesheet line ${lineNo}: empty required value '${key}'")
            }
            if (!(row.patient_id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/))
                throw new IllegalArgumentException("Samplesheet line ${lineNo}: unsafe patient_id '${row.patient_id}'")
            if (!seen.add(row.patient_id)) throw new IllegalArgumentException("Samplesheet duplicate patient_id: ${row.patient_id}")
            if (row.sex && !(row.sex in ['female', 'male', 'unknown']))
                throw new IllegalArgumentException("Samplesheet ${row.patient_id}: sex must be female, male, unknown or blank")
            ['purity', 'ploidy'].each { key ->
                if (row[key]) {
                    double number
                    try { number = Double.parseDouble(row[key]) }
                    catch (NumberFormatException e) { throw new IllegalArgumentException("Samplesheet ${row.patient_id}: ${key} must be numeric or blank") }
                    if (!Double.isFinite(number) || number <= 0 || (key == 'purity' && number > 1))
                        throw new IllegalArgumentException("Samplesheet ${row.patient_id}: invalid ${key} '${row[key]}'")
                }
            }
            if (new File(row.tumor_bam).canonicalPath == new File(row.normal_bam).canonicalPath)
                throw new IllegalArgumentException("Samplesheet ${row.patient_id}: tumor_bam and normal_bam identify the same file; use an isolated, explicitly prepared control")
            row.missing_metadata = (REQUIRED - NONEMPTY).findAll { !row[it] }
            rows << row
        }
        if (!rows) throw new IllegalArgumentException('Samplesheet has no data rows')
        return rows
    }
}
