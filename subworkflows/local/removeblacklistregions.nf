

include { BEDTOOLS_INTERSECT       } from '../../modules/nf-core/bedtools/intersect/main'
include { SAMTOOLS_INDEX           } from '../../modules/nf-core/samtools/index/main'

workflow REMOVE_BLACKLIST_REGIONS {

    take:
    bam
    index
    blacklist
    chrom_sizes

    main:

    ch_versions = Channel.empty()

    // MUDOLE: BEDTOOLS_INTERSECT
    BEDTOOLS_INTERSECT (
        bam.combine(blacklist),
        chrom_sizes.map{[ [:], it]}
    )

    SAMTOOLS_INDEX ( BEDTOOLS_INTERSECT.out.intersect )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    bam      = BEDTOOLS_INTERSECT.out.intersect     // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai               // channel: [ val(meta), [ bai ] ]

    versions = ch_versions                          // channel: [ versions.yml ]
}

