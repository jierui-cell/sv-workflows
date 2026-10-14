/** Separate, strict per-pair inputs because JaBbA consumes precomputed junctions and coverage. */
class JabbaInputs {
    static final List REQUIRED = ['pair_id', 'junction_vcf', 'coverage']

    static Map<String, Map> read(java.nio.file.Path path, java.nio.file.Path projectRoot,
                                 Set<String> selectedIds, boolean requireCoverage) {
        def lines = java.nio.file.Files.readAllLines(path)
        if (!lines) throw new IllegalArgumentException('JaBbA input TSV is empty')
        def header = lines[0].replaceFirst('^\uFEFF', '').split('\t', -1).toList()
        if (header.toSet().size() != header.size() || header.any { !it || it != it.trim() })
            throw new IllegalArgumentException('JaBbA input TSV has invalid or duplicate headers')
        def missing = REQUIRED - header
        if (missing) throw new IllegalArgumentException("JaBbA input TSV missing columns: ${missing.join(', ')}")
        def rows = [:]
        lines.drop(1).eachWithIndex { line, index ->
            def number = index + 2
            def cells = line.split('\t', -1).toList()
            if (cells.size() != header.size())
                throw new IllegalArgumentException("JaBbA input TSV line ${number}: expected ${header.size()} fields")
            def row = [header, cells].transpose().collectEntries { [(it[0]): it[1].trim()] }
            if (!row.pair_id) throw new IllegalArgumentException("JaBbA input TSV line ${number}: empty pair_id")
            if (!(row.pair_id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/))
                throw new IllegalArgumentException("JaBbA input TSV line ${number}: unsafe pair_id")
            if (rows.containsKey(row.pair_id))
                throw new IllegalArgumentException("JaBbA input TSV duplicate pair_id: ${row.pair_id}")
            if (!selectedIds.contains(row.pair_id)) {
                rows[row.pair_id] = row
                return
            }
            (requireCoverage ? ['junction_vcf', 'coverage'] : ['junction_vcf']).each { key ->
                if (!row[key]) throw new IllegalArgumentException("JaBbA input TSV line ${number}: empty ${key}")
                def supplied = new File(row[key])
                def resolved = (supplied.isAbsolute() ? supplied : projectRoot.resolve(row[key]).toFile()).canonicalFile
                if (!resolved.isFile() || resolved.length() == 0)
                    throw new IllegalArgumentException("JaBbA ${row.pair_id}: ${key} not found or empty: ${resolved}")
                row[key] = resolved.path
            }
            rows[row.pair_id] = row
        }
        if (!rows) throw new IllegalArgumentException('JaBbA input TSV has no pairs')
        rows
    }
}
