include { REMOVE_READ_LOW_MAPQ      } from '../../modules/local/removelowmapq.nf'
include { SAMTOOLS_INDEX            } from '../../modules/nf-core/samtools/index/main'

workflow FILTER_MAPQ {

    take:
    bam
    index
    min_mapq

    main:

    ch_versions = Channel.empty()

    REMOVE_READ_LOW_MAPQ (
        bam.join(index),
        params.min_mapq
    )

    SAMTOOLS_INDEX (
        REMOVE_READ_LOW_MAPQ.out.bam
    )

    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    bam      = REMOVE_READ_LOW_MAPQ.out.bam    // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
