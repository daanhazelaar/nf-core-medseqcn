include { SUBSAMPLE_BAM  } from '../../modules/local/subsamplebam.nf'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'

workflow EQUALIZE_COVERAGE {

    take:
    bam         // channel: [ meta, bam ]
    bai         // channel: [ meta, bai ]
    flagstat    // channel: [ meta, flagstat_file ]
    samplesheet // channel: [ meta, fastqs, methylated_bam, assay, sex ]

    main:

    // Join bam, bai, flagstat and assay (from samplesheet); parse mapped read count
    // from flagstat. Group samples by meta.patient to identify comparison pairs.
    bam
        .join(bai)
        .join(flagstat)
        .join(samplesheet.map { meta, fastqs, methylated_bam, assay, sex -> [ meta, assay ] })
        .map { meta, bam, bai, flagstat_file, assay ->
            // flagstat line format: "123456 + 0 mapped (80.00% : N/A)"
            def mapped_reads = flagstat_file.text
                .readLines()
                .find { it.contains(' mapped (') }
                ?.split(' ')[0]
                ?.toLong() ?: 0L
            // Use patient ID as group key; fall back to sample ID if patient not set
            def group_key = meta.patient ?: meta.id
            return [ group_key, meta, bam, bai, assay, mapped_reads ]
        }
        .groupTuple(by: 0)
        .flatMap { group_key, metas, bams, bais, assays, read_counts ->
            // Build a list of per-sample maps for easier manipulation
            def samples = [ metas, bams, bais, assays, read_counts ]
                .transpose()
                .collect { m, b, bi, a, rc -> [ meta: m, bam: b, bai: bi, assay: a, reads: rc ] }

            def result  = []
            def handled = [] as Set

            // Defined comparison pairs: quadf<->swgs and quadm<->medseq.
            // Within each pair the lower-coverage sample sets the target;
            // the higher-coverage sample is downsampled to match.
            def pairs = [ [ 'quadf', 'swgs' ], [ 'quadm', 'medseq' ] ]
            pairs.each { assay_a, assay_b ->
                def s_a = samples.find { it.assay == assay_a }
                def s_b = samples.find { it.assay == assay_b }
                if (s_a) handled << assay_a
                if (s_b) handled << assay_b

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

            // Any assay not covered by a defined pair passes through unchanged
            samples.findAll { !(it.assay in handled) }.each { s ->
                result << [ s.meta, s.bam, s.bai, 1.0 ]
            }

            return result
        }
        .branch { meta, bam, bai, fraction ->
            subsample:   fraction < 1.0   // needs downsampling
            passthrough: true              // already at target coverage
        }
        .set { ch_branched }

    // Downsample the higher-coverage member of each pair
    SUBSAMPLE_BAM ( ch_branched.subsample )

    SAMTOOLS_INDEX ( SUBSAMPLE_BAM.out.bam )

    // Recombine downsampled and pass-through samples into a single output channel
    SUBSAMPLE_BAM.out.bam
        .join(SAMTOOLS_INDEX.out.bai)
        .mix(
            ch_branched.passthrough
                .map { meta, bam, bai, fraction -> [ meta, bam, bai ] }
        )
        .set { out }

    emit:
    bam = out.map { meta, bam, bai -> [ meta, bam ] }
    bai = out.map { meta, bam, bai -> [ meta, bai ] }
}
