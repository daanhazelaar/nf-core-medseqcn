include { SAMTOOLS_MERGE                              } from '../../modules/local/samtoolsmerge.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_COMBINED   } from '../../modules/nf-core/samtools/index/main'

workflow COMBINE_QUAD_FRACTIONS {

    take:
    bam   // channel: [ meta, bam ]   (final processed BAMs, equalized when enabled)
    bai   // channel: [ meta, bai ]

    main:

    // Group all samples of a patient together so the QUAD-M and QUAD-F fractions
    // can be merged into a single combined pseudo-sample for copy-number analysis.
    bam
        .join(bai)
        .map { meta, bam_file, bai_file ->
            def group_key = meta.patient ?: meta.id
            return [ group_key, meta, bam_file, bai_file, meta.assay, meta.methylation ]
        }
        .groupTuple(by: 0)
        .map { group_key, metas, bams, bais, assays, methylations ->
            def samples = [ metas, bams, bais, assays, methylations ]
                .transpose()
                .collect { m, b, bi, a, meth ->
                    [ meta: m, bam: b, bai: bi, assay: a, methylation: meth ]
                }

            // Combine only the non-methylated fractions of QUAD-M and QUAD-F.
            def s_quadf = samples.find { it.assay == "quadf" && it.methylation == "nonmeth" }
            def s_quadm = samples.find { it.assay == "quadm" && it.methylation == "nonmeth" }

            // Only emit a combined sample when both fractions are present for the patient.
            if (s_quadf && s_quadm) {
                // Synthesize the combined-sample meta: inherit shared fields (patient, sex,
                // ...) from the QUAD-F meta and overwrite the fraction-specific tags.
                def combined_meta = s_quadf.meta + [
                    id          : "${group_key}_quadcombined",
                    assay       : "quadcombined",
                    methylation : "nonmeth"
                ]
                return [ combined_meta, [ s_quadf.bam, s_quadm.bam ] ]
            }
            return null
        }
        .filter { it != null }
        .set { ch_to_merge }

    // MODULE: SAMTOOLS_MERGE — naive concatenation, no re-markdup (the two fractions
    // are separate libraries with distinct read names, so there are no name collisions).
    SAMTOOLS_MERGE ( ch_to_merge )

    SAMTOOLS_INDEX_COMBINED ( SAMTOOLS_MERGE.out.bam )

    SAMTOOLS_MERGE.out.bam
        .join(SAMTOOLS_INDEX_COMBINED.out.bai)
        .set { out }

    emit:
    bam = out.map { meta, bam_file, bai_file -> [ meta, bam_file ] }
    bai = out.map { meta, bam_file, bai_file -> [ meta, bai_file ] }
}
