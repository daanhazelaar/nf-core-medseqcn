
include { GET_READ_IDS_FROM_BAM     } from '../../modules/local/getreadidsfrombam.nf'
include { SUBSET_BAM_READ_IDS } from '../../modules/local/subsetbamreadids.nf'
include { SAMTOOLS_INDEX  } from '../../modules/nf-core/samtools/index/main'

workflow SUBSET_MEDSEQ_DATA {

    take:
    bam
    index
    samplesheet
    fasta

    main:

    // Split channel based on assay type medseq or swgs. Subsetting is only performend for medseq data.
    bam
        .join(index)
        .join(samplesheet)
        .branch { meta, bam, bai, fastq, methylated_bam, assay, sex ->
            medseq: assay == "medseq"
            swgs: assay == "swgs"
        }
        .set{ ch_reads_split_assay }

    // MODULE: GET_READ_IDS_FROM_BAM
    GET_READ_IDS_FROM_BAM (
        ch_reads_split_assay.medseq
            .map{meta, bam, bai, fastq, methylated_bam, assay, sex -> return[ meta, methylated_bam ]}
    )

    // MODULE: SUB_SET_BAM_READ_IDS
    SUBSET_BAM_READ_IDS (
        ch_reads_split_assay.medseq
            .map{ meta, bam, bai, fastq, methylated_bam, assay, sex -> return[ meta, bam, bai ]}
            .join(GET_READ_IDS_FROM_BAM.out),
        fasta.map{[ [:], it]}
    )

    // MODULE: SAMTOOLS_INDEX
    SAMTOOLS_INDEX (
        SUBSET_BAM_READ_IDS.out.unselected
    )

    // Combinde subset medseq data and swgs data in channel
    SUBSET_BAM_READ_IDS.out.unselected
        .join(SAMTOOLS_INDEX.out.bai)
        .join(samplesheet)
        .mix(ch_reads_split_assay.swgs)
        .set{out}

    emit:
    bam = out.map{meta, bam, index, fastq, methylated_bam, assay, sex -> [ meta, bam ]}
    bai = out.map{meta, bam, index, fastq, methylated_bam, assay, sex -> [ meta, index ]}

}

