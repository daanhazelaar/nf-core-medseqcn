
include { SAMTOOLS_VIEW as SUBSAMPLE     } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_SUBSAMPLE } from '../../modules/nf-core/samtools/index/main'


workflow SUBSAMPLE_BAM {

    take:
    samplesheet
    bam
    index
    fasta

    main:

    ch_versions = Channel.empty()

    samplesheet
        .join(bam)
        .join(index)
        .combine(params.subsample_fractions)
        .flatMap {meta, fastqs, methylated_bam, assay, sex, bam, bai, subsample ->
            subsample.collect { value ->
                return [ meta + [subsample_fraction: value], bam, bai, subsample, assay, sex]
                // return [ [[id:meta.id + "_${value}"], [single_end:meta.single_end], [subsample_fraction: value]], bam, bai, subsample]
            }
        }
        .set {bam_subsample_split}

    // bam_subsample_split.view()

    SUBSAMPLE (
        bam_subsample_split.map {meta, bam, bai, subsample, assay, sex -> return [ meta, bam, bai]},
        fasta,
        []
    )

    SAMTOOLS_INDEX_SUBSAMPLE ( SUBSAMPLE.out.bam )

    bam_subsample_split
        .map {meta, bam, bai, subsample, assay, sex -> return [ meta, subsample, assay, sex]}
        .set {more_meta}


    SUBSAMPLE.out.bam
        .join(SAMTOOLS_INDEX_SUBSAMPLE.out.bai)
        .join(more_meta)
        .set { bam_out }

    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_SUBSAMPLE.out.versions.first())

    emit:
    bam           = bam_out
    // bam           = SUBSAMPLE.out.bam
    // bai           = SAMTOOLS_INDEX_SUBSAMPLE.out.bai

    // versions = ch_versions
}

