/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                    } from '../modules/nf-core/fastqc/main'
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_medseqcn_pipeline'

include { PREPARE_REFERENCE_GENOME  } from '../subworkflows/local/prepare_reference_genome'

include { FASTP                     } from '../modules/nf-core/fastp/main'
include { BWA_MEM                   } from '../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_RAW        } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_RAW      } from '../modules/nf-core/samtools/index/main'

include { SUBSET_MEDSEQ_DATA        } from '../subworkflows/local/subsetmedseqdata.nf'
include { SIZE_SELECTION            } from '../modules/local/sizeselection.nf'
include { SAMBAMBA_MARKDUP          } from '../modules/nf-core/sambamba/markdup/main'
include { FILTER_MAPQ               } from '../subworkflows/local/filtermapq.nf'
include { REMOVE_BLACKLIST_REGIONS  } from '../subworkflows/local/removeblacklistregions.nf'

include { BAM_SORT_STATS_SAMTOOLS   } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { SAMTOOLS_COVERAGE         } from '../modules/nf-core/samtools/coverage/main'

include { HMMCOPY_READCOUNTER       } from '../modules/nf-core/hmmcopy/readcounter/main'
include { ICHORCNA_RUN_CUSTOM       } from '../modules/local/ichorcnaruncustom.nf'

include { SUBSAMPLE_BAM             } from '../subworkflows/local/subsample_bam.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MEDSEQCN {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // SUBWORKFLOW: PREPARE_REFERENCE_GENOME
    PREPARE_REFERENCE_GENOME ()

    // MODULE: Run FastQC
    FASTQC (
        ch_samplesheet.map{meta, fastqs, methylated_bam, assay, sex -> return[ meta, fastqs ]}
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    // MODULE: MULTIQC
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    // MODULE: FASTP
    discard_trimmed_pass = false
    save_trimmed_fail = false
    save_merged = false

    FASTP (
        ch_samplesheet.map{meta, fastqs, methylated_bam, assay, sex -> return[ meta, fastqs ]},
        [],
        discard_trimmed_pass,
        save_trimmed_fail,
        save_merged
    )

    // MODULE: BWA_MEM
    BWA_MEM (
        FASTP.out.reads,
        PREPARE_REFERENCE_GENOME.out.bwa_index,
        PREPARE_REFERENCE_GENOME.out.fasta.map{[ [:], it]},
        false
    )

    // MODULE: SAMTOOLS_SORT
    SAMTOOLS_SORT_RAW (
        BWA_MEM.out.bam,
        PREPARE_REFERENCE_GENOME.out.fasta.map{[ [:], it]}
    )

    // MODULE: SAMTOOLS_INDEX
    SAMTOOLS_INDEX_RAW (
        SAMTOOLS_SORT_RAW.out.bam
    )

    // SUBWORKFLOW: SUBSET_MEDSEQ_DATA
    SUBSET_MEDSEQ_DATA (
        SAMTOOLS_SORT_RAW.out.bam,
        SAMTOOLS_INDEX_RAW.out.bai,
        ch_samplesheet,
        PREPARE_REFERENCE_GENOME.out.fasta
    )

    // CUSTOM MUDOLE: SIZE_SELECTION
    if (params.insert_size_selection) {
        SIZE_SELECTION (
            SUBSET_MEDSEQ_DATA.out.bam.join(SUBSET_MEDSEQ_DATA.out.bai),
            params.min_size_selection,
            params.max_size_selection
        )
    }

    // MUDOLE: SAMBAMBA_MARKDUP
    SAMBAMBA_MARKDUP (
        (params.insert_size_selection ? SIZE_SELECTION.out.bam : SUBSET_MEDSEQ_DATA.out.bam)
    )

    // SUBWORKFLOW: FILTER_MAPQ
    if (params.remove_read_low_mapq) {
        FILTER_MAPQ (
            SAMBAMBA_MARKDUP.out.bam,
            SAMBAMBA_MARKDUP.out.bai,
            params.min_mapq
        )
    }

    // SUBWORKFLOW: REMOVE_BLACKLIST_REGIONS
    REMOVE_BLACKLIST_REGIONS (
        (params.remove_read_low_mapq ? FILTER_MAPQ.out.bam : SAMBAMBA_MARKDUP.out.bam),
        (params.remove_read_low_mapq ? FILTER_MAPQ.out.bai : SAMBAMBA_MARKDUP.out.bai),
        PREPARE_REFERENCE_GENOME.out.blacklist,
        PREPARE_REFERENCE_GENOME.out.chrom_sizes
    )

    // SUBWORKFLOW: SUBSAMPLE_BAM
    if (params.subsample) {
        SUBSAMPLE_BAM (
            ch_samplesheet,
            REMOVE_BLACKLIST_REGIONS.out.bam,
            REMOVE_BLACKLIST_REGIONS.out.bai,
            PREPARE_REFERENCE_GENOME.out.fasta.map{[ [:], it]}
        )
    }

    // SUBWORKFLOW: BAM_SORT_STATS_SAMTOOLS
    BAM_SORT_STATS_SAMTOOLS (
        (params.subsample ? SUBSAMPLE_BAM.out.bam.map{meta, bam, bai, subsample, assay, sex -> return [ meta, bam]} : REMOVE_BLACKLIST_REGIONS.out.bam),
        PREPARE_REFERENCE_GENOME.out.fasta.map{[ [:], it]}
    )

    // MUDOLE: SAMTOOLS_COVERAGE
    SAMTOOLS_COVERAGE (
        BAM_SORT_STATS_SAMTOOLS.out.bam.join(BAM_SORT_STATS_SAMTOOLS.out.bai),
        PREPARE_REFERENCE_GENOME.out.fasta.map{[ [:], it]},
        PREPARE_REFERENCE_GENOME.out.fai.map{[ [:], it]}
    )

    // MUDOLE: HMMCOPY_READCOUNTER
    HMMCOPY_READCOUNTER (
        BAM_SORT_STATS_SAMTOOLS.out.bam.join(BAM_SORT_STATS_SAMTOOLS.out.bai)
    )

    SUBSAMPLE_BAM.out.bam
        .map{meta, bam, bai, subsample, assay, sex -> return [ meta, subsample, assay, sex]}
        .join(HMMCOPY_READCOUNTER.out.wig)
        .set{input_ichorcna}

    ICHORCNA_RUN_CUSTOM (
        input_ichorcna.map{meta, subsample, assay, sex, wig -> return [ meta, wig, sex]},
        Channel.value(file(params.panel_of_normals_medseq)),
        Channel.value(file(params.gc_wig)),
        Channel.value(file(params.map_wig)),
        Channel.value(file(params.centromere)),
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
