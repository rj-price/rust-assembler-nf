/*
 * Genome-size arithmetic shared between modules.
 *
 * Lives in lib/ because Nextflow puts that directory on the classpath for the whole
 * workflow, including module script blocks. A plain `def` in main.nf would not be visible
 * there, and duplicating the parsing in each module is how the two sizes drifted apart
 * before (COVERAGE_SUMMARY claimed to derive expected depth and then hard-coded it).
 *
 * All of this is evaluated at TASK time, not config-parse time, so `--genome_size` given on
 * the command line is always the value used.
 */
class GenomeSize {

    /*
     * '525m' / '1.05g' / '525000000' -> 525000000L
     *
     * Accepts the k/m/g suffixes the assemblers use, case-insensitively, and bare integers.
     */
    static long parse(def size) {
        if (size == null) {
            throw new IllegalArgumentException("genome size is null — set --genome_size")
        }
        def s = size.toString().trim().toLowerCase()
        def m = (s =~ /^([0-9]*\.?[0-9]+)\s*([kmg]?)b?$/)
        if (!m.matches()) {
            throw new IllegalArgumentException(
                "cannot parse genome size '${size}' — expected e.g. 525m, 1.05g or 525000000")
        }
        def value = m[0][1].toBigDecimal()
        def mult  = ['': 1L, 'k': 1_000L, 'm': 1_000_000L, 'g': 1_000_000_000L][m[0][2]]
        return (value * mult).toLong()
    }

    /*
     * Expected total size of a dikaryotic assembly that has RETAINED both nuclei.
     *
     * Two nuclei, so 2 x haploid, with a tolerance band either side. The band is not
     * symmetric on purpose: a little over is ordinary (residual duplication, uncollapsed
     * repeat), whereas landing under 2x is the collapse this pipeline exists to detect.
     */
    static long dikaryonMin(def haploid, def frac) {
        return (parse(haploid) * 2 * frac.toBigDecimal()).toLong()
    }

    static long dikaryonMax(def haploid, def frac) {
        return (parse(haploid) * 2 * frac.toBigDecimal()).toLong()
    }

    /*
     * Expected read depth per haplotype: yield spread over both nuclei.
     *
     * Collapsed/shared sequence draws reads from both nuclei and so sits at twice this.
     * Returns null when the yield is unknown, which the caller reports as "not compared"
     * rather than inventing a number.
     */
    static BigDecimal expectedHaplotypeDepth(def yieldBases, def haploid) {
        if (yieldBases == null) return null
        return (yieldBases.toBigDecimal() / (2.0 * parse(haploid)))
                   .setScale(1, java.math.RoundingMode.HALF_UP)
    }
}
