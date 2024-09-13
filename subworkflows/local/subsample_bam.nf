
include { SAMTOOLS_VIEW as SUBSAMPLE     } from '../../modules/nf-core/samtools/view/main'
// include { SAMTOOLS_VIEW as SUBSAMPLE_2     } from '../../modules/nf-core/samtools/view/main'
// include { SAMTOOLS_VIEW as SUBSAMPLE_3     } from '../../modules/nf-core/samtools/view/main'
// include { SAMTOOLS_VIEW as SUBSAMPLE_4     } from '../../modules/nf-core/samtools/view/main'
// include { SAMTOOLS_VIEW as SUBSAMPLE_5     } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
// include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_2 } from '../../modules/nf-core/samtools/index/main'
// include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_3 } from '../../modules/nf-core/samtools/index/main'
// include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_4 } from '../../modules/nf-core/samtools/index/main'
// include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_5 } from '../../modules/nf-core/samtools/index/main'


workflow SUBSAMPLE_BAM {

    take:
    bam
    index
    fasta

    main:

    ch_versions = Channel.empty()

    // bam.view()

    bam
        .join(index)
        .combine(params.subsample_fractions)
        .flatMap {meta, bam, bai, subsample ->
            subsample.collect { value ->
                return [ meta + [subsample_fraction: value], bam, bai, subsample]
            }
        }
        // .flatten()
        .set {testing}

    // testing.map {meta, bam, bai, subsample -> return [ meta, bam, bai]}
    //     .view()

    // testing.map {meta, bam, bai, subsample -> return [ meta, bam, bai]}
    //     .map{meta, bam, bai, subsample -> return[ meta, methylated_bam ]}
    //     .view()

    SUBSAMPLE (
        testing.map {meta, bam, bai, subsample -> return [ meta, bam, bai]},
        fasta,
        []
    )

    SUBSAMPLE.out.bam.view()

    // SAMTOOLS_INDEX ( SUBSAMPLE.out.bam )

    // SUBSAMPLE_1.out.bam
    //     .map{meta, bam -> return [ meta + [subsample:0.75], bam ]}
    //     .view()

    // SUBSAMPLE_1.out.bam
    //     .map { meta, bam ->
    //         def values = params.subsample_fractions
    //         values.collect { value ->

    //             // def new_meta = [meta.id + "_${value}"]

    //             return [ meta + [subsamples: value], bam]  // Return a new tuple
    //         }
    // }
    //     .flatten() // Flatten the channel to merge the inner lists into one stream


    // testing
    //     .map {
    //         meta, bam, subsample ->
    //             return [meta, subsample]
    //     }
    //     .view()

    // SUBSAMPLE_1.out.bam.
    //     .map { it ->
    //         def values = [1,2,3,4,5]

    //         values.collect {value ->
    //             def new_meta = item.meta + "_${value}"
    //             return [new_meta, ]
    //         }

    //     }

    // SAMTOOLS_INDEX_1 ( SUBSAMPLE_1.out.bam )

    // SUBSAMPLE_2 (
    //     bam.join(index),
    //     fasta,
    //     []
    // )

    // SAMTOOLS_INDEX_2 ( SUBSAMPLE_2.out.bam )

    // SUBSAMPLE_3 (
    //     bam.join(index),
    //     fasta,
    //     []
    // )

    // SAMTOOLS_INDEX_3 ( SUBSAMPLE_3.out.bam )

    // SUBSAMPLE_4 (
    //     bam.join(index),
    //     fasta,
    //     []
    // )

    // SAMTOOLS_INDEX_4 ( SUBSAMPLE_4.out.bam )

    // SUBSAMPLE_5 (
    //     bam.join(index),
    //     fasta,
    //     []
    // )

    // SAMTOOLS_INDEX_5 ( SUBSAMPLE_5.out.bam )


    // ch_versions = ch_versions.mix(SAMTOOLS_INDEX_1.out.versions.first())

    emit:
    bam = bam
    // bam_1      = SUBSAMPLE_1.out.bam
    // bai_1      = SAMTOOLS_INDEX_1.out.bai
    // bam_2      = SUBSAMPLE_2.out.bam
    // bai_2      = SAMTOOLS_INDEX_2.out.bai
    // bam_3      = SUBSAMPLE_3.out.bam
    // bai_3      = SAMTOOLS_INDEX_3.out.bai
    // bam_4      = SUBSAMPLE_4.out.bam
    // bai_4      = SAMTOOLS_INDEX_4.out.bai
    // bam_5      = SUBSAMPLE_5.out.bam
    // bai_5      = SAMTOOLS_INDEX_5.out.bai

    // versions = ch_versions
}

