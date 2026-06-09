
include { GET_READ_IDS_FROM_BAM                                  } from '../../modules/local/getreadidsfrombam.nf'
include { SAMTOOLS_VIEW as SUBSET_BAM_READ_IDS                   } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_NONMETH               } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_METH                  } from '../../modules/nf-core/samtools/index/main'

workflow SUBSET_MEDSEQ_DATA {

    take:
    bam
    index
    samplesheet
    fasta

    main:

    // Split channel based on assay type. Methylated-read identification is performed for
    // medseq, quadm, and quadf (all LpnPI-based). swgs passes through unchanged.
    bam
        .join(index)
        .join(samplesheet)
        .branch { meta, bam, bai, fastq, methylated_bam, assay, sex ->
            medseq_based: assay == "medseq" || assay == "quadm" || assay == "quadf"
            swgs: assay == "swgs"
        }
        .set{ ch_reads_split_assay }

    // MODULE: GET_READ_IDS_FROM_BAM
    GET_READ_IDS_FROM_BAM (
        ch_reads_split_assay.medseq_based
            .map{meta, bam, bai, fastq, methylated_bam, assay, sex -> return[ meta, methylated_bam ]}
    )

    // MODULE: SUB_SET_BAM_READ_IDS
    // Emits both the methylated-matching reads (.bam) and the non-methylated remainder (.unselected).
    SUBSET_BAM_READ_IDS (
        ch_reads_split_assay.medseq_based
            .map{ meta, bam, bai, fastq, methylated_bam, assay, sex -> return[ meta, bam, bai ]}
            .join(GET_READ_IDS_FROM_BAM.out),
        fasta.map{[ [:], it]}
    )

    // Non-methylated stream: applies to medseq/quadm/quadf. Tag meta.methylation = "nonmeth"
    // and suffix meta.id so all downstream filenames disambiguate without changes to publish rules.
    SUBSET_BAM_READ_IDS.out.unselected
        .map { meta, bam ->
            def m = meta + [ methylation: "nonmeth", id: "${meta.id}_nonmeth" ]
            [ m, bam ]
        }
        .set { ch_nonmeth_bam }

    SAMTOOLS_INDEX_NONMETH ( ch_nonmeth_bam )

    ch_nonmeth_bam
        .join(SAMTOOLS_INDEX_NONMETH.out.bai)
        .set { ch_nonmeth_out }

    // Methylated stream: only emitted for medseq and quadm (per project scope).
    // quadf is dropped here — its methylated fraction is not analyzed downstream.
    SUBSET_BAM_READ_IDS.out.bam
        .filter { meta, bam -> meta.assay == "medseq" || meta.assay == "quadm" }
        .map { meta, bam ->
            def m = meta + [ methylation: "meth", id: "${meta.id}_meth" ]
            [ m, bam ]
        }
        .set { ch_meth_bam }

    SAMTOOLS_INDEX_METH ( ch_meth_bam )

    ch_meth_bam
        .join(SAMTOOLS_INDEX_METH.out.bai)
        .set { ch_meth_out }

    // swgs passthrough: tag as nonmeth so the methylation axis is consistent for downstream pairing.
    ch_reads_split_assay.swgs
        .map { meta, bam, bai, fastq, methylated_bam, assay, sex ->
            def m = meta + [ methylation: "nonmeth" ]
            [ m, bam, bai ]
        }
        .set { ch_swgs_out }

    // Mix all streams into the unified output channel.
    ch_nonmeth_out
        .mix(ch_meth_out)
        .mix(ch_swgs_out)
        .set { out }

    emit:
    bam = out.map{meta, bam, bai -> [ meta, bam ]}
    bai = out.map{meta, bam, bai -> [ meta, bai ]}

}
