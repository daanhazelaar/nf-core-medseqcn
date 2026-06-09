include { SUBSAMPLE_BAM  } from '../../modules/local/subsamplebam.nf'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'

workflow EQUALIZE_COVERAGE {

    take:
    bam         // channel: [ meta, bam ]   (meta carries assay + methylation)
    bai         // channel: [ meta, bai ]
    flagstat    // channel: [ meta, flagstat_file ]

    main:

    // Join bam, bai, flagstat; parse mapped read count from flagstat.
    // Group samples by meta.patient to identify comparison pairs.
    bam
        .join(bai)
        .join(flagstat)
        .map { meta, bam, bai, flagstat_file ->
            // flagstat line format: "123456 + 0 mapped (80.00% : N/A)"
            def mapped_reads = flagstat_file.text
                .readLines()
                .find { it.contains(' mapped (') }
                ?.split(' ')[0]
                ?.toLong() ?: 0L
            // Use patient ID as group key; fall back to sample ID if patient not set
            def group_key = meta.patient ?: meta.id
            return [ group_key, meta, bam, bai, meta.assay, meta.methylation, mapped_reads ]
        }
        .groupTuple(by: 0)
        .flatMap { group_key, metas, bams, bais, assays, methylations, read_counts ->
            // Build a list of per-sample maps for easier manipulation
            def samples = [ metas, bams, bais, assays, methylations, read_counts ]
                .transpose()
                .collect { m, b, bi, a, meth, rc ->
                    [ meta: m, bam: b, bai: bi, assay: a, methylation: meth, reads: rc ]
                }

            def result  = []
            def handled = [] as Set

            // Defined comparison pairs along the (assay, methylation) axis:
            //   - quadf  (nonmeth) <-> swgs   (nonmeth)
            //   - quadm  (nonmeth) <-> medseq (nonmeth)
            //   - quadm  (meth)    <-> medseq (meth)
            // Within each pair the lower-coverage sample sets the target;
            // the higher-coverage sample is downsampled to match.
            def pairs = [
                [ [ 'quadf',  'nonmeth' ], [ 'swgs',   'nonmeth' ] ],
                [ [ 'quadm',  'nonmeth' ], [ 'medseq', 'nonmeth' ] ],
                [ [ 'quadm',  'meth'    ], [ 'medseq', 'meth'    ] ],
            ]
            pairs.each { key_a, key_b ->
                def s_a = samples.find { it.assay == key_a[0] && it.methylation == key_a[1] }
                def s_b = samples.find { it.assay == key_b[0] && it.methylation == key_b[1] }
                if (s_a) handled << "${key_a[0]}_${key_a[1]}".toString()
                if (s_b) handled << "${key_b[0]}_${key_b[1]}".toString()

                if (s_a && s_b) {
                    def min_reads = Math.min(s_a.reads, s_b.reads) as double
                    result << [ s_a.meta, s_a.bam, s_a.bai, min_reads / s_a.reads ]
                    result << [ s_b.meta, s_b.bam, s_b.bai, min_reads / s_b.reads ]
                } else {
                    // Unpaired sample: pass through unchanged
                    if (s_a) result << [ s_a.meta, s_a.bam, s_a.bai, 1.0 ]
                    if (s_b) result << [ s_b.meta, s_b.bam, s_b.bai, 1.0 ]
                }
            }

            // Any sample not covered by a defined pair passes through unchanged
            samples.findAll { !("${it.assay}_${it.methylation}".toString() in handled) }.each { s ->
                result << [ s.meta, s.bam, s.bai, 1.0 ]
            }

            return result
        }
        .set { ch_with_fractions }

    // All samples pass through SUBSAMPLE_BAM.
    // Fraction < 1.0 → samtools view -s (downsample).
    // Fraction = 1.0 → cp (passthrough, renamed to _equalized.bam for consistency).
    SUBSAMPLE_BAM ( ch_with_fractions )

    SAMTOOLS_INDEX ( SUBSAMPLE_BAM.out.bam )

    SUBSAMPLE_BAM.out.bam
        .join(SAMTOOLS_INDEX.out.bai)
        .set { out }

    emit:
    bam = out.map { meta, bam, bai -> [ meta, bam ] }
    bai = out.map { meta, bam, bai -> [ meta, bai ] }
}
