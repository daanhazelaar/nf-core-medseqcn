//
// Uncompress and prepare reference genome files
// Most of this file I got from: https://github.com/nf-core/atacseq/blob/2.1.2/subworkflows/local/prepare_genome.nf

include { UNTAR                 } from '../../modules/nf-core/untar/main'
include { GUNZIP                } from '../../modules/nf-core/gunzip/main'
include { BWA_INDEX             } from '../../modules/nf-core/bwa/index/main'
include { GET_CHROMOMOSOME_SIZES} from '../../modules/local/getchromosomesizes.nf'


workflow PREPARE_REFERENCE_GENOME {
    take:
    // prepare_tool_index // string  : tool to prepare index for

    main:
    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    ch_fasta = Channel.empty()
    if (params.fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], params.fasta ] ).gunzip.map{ it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(params.fasta))
    }

    //
    // Create chromosome sizes file
    //
    GET_CHROMOMOSOME_SIZES ( ch_fasta.map { [ [:], it ] } )
    ch_chrom_sizes = GET_CHROMOMOSOME_SIZES.out.sizes.map { it[1] }
    ch_fai         = GET_CHROMOMOSOME_SIZES.out.fai.map{ it[1] }
    ch_versions    = ch_versions.mix(GET_CHROMOMOSOME_SIZES.out.versions)


    //
    // Uncompress blacklist file if required
    //
    ch_blacklist = Channel.empty()
    if (params.blacklist) {
        if (params.blacklist.endsWith('.gz')) {
            ch_blacklist = GUNZIP_BLACKLIST ( [ [:], params.blacklist ] ).gunzip.map{ it[1] }
            ch_versions  = ch_versions.mix(GUNZIP_BLACKLIST.out.versions)
        } else {
            ch_blacklist = Channel.value(file(params.blacklist))
        }
    }

    // Uncompress BWA index if required
    //
    prepare_tool_index ='bwa'

    ch_bwa_index = Channel.empty()
    if (prepare_tool_index == 'bwa') {
        if (params.bwa_index) {
            if (params.bwa_index.endsWith('.tar.gz')) {
                ch_bwa_index = UNTAR_BWA_INDEX ( [ [:], params.bwa_index ] ).untar
                ch_versions  = ch_versions.mix(UNTAR_BWA_INDEX.out.versions)
            } else {
                ch_bwa_index = [ [:], file(params.bwa_index) ]
            }
        } else {
            ch_bwa_index = BWA_INDEX ( ch_fasta.map { [ [:], it ] } ).index
            ch_versions  = ch_versions.mix(BWA_INDEX.out.versions)
        }
    }

    emit:
    fasta         = ch_fasta                      //    path: genome.fasta
    fai           = ch_fai                        //    path: genome.fai
    blacklist     = ch_blacklist                  //    path: encode blacklist file
    // gtf           = ch_gtf                        //    path: genome.gtf
    // gene_bed      = ch_gene_bed                   //    path: gene.bed
    // tss_bed       = ch_tss_bed                    //    path: tss.bed
    chrom_sizes   = ch_chrom_sizes                //    path: genome.sizes
    // filtered_bed  = ch_genome_filtered_bed        //    path: *.include_regions.bed
    bwa_index     = ch_bwa_index                  //    path: bwa/index/
    // bowtie2_index = ch_bowtie2_index              //    path: bowtie2/index/
    // chromap_index = ch_chromap_index              //    path: genome.index
    // star_index    = ch_star_index                 //    path: star/index/
    // autosomes     = ch_genome_autosomes           //    path: *.autosomes.txt
    // macs_gsize    = ch_macs_gsize                 // integer: MACS2 genome size

    versions      = ch_versions.ifEmpty(null)     // channel: [ versions.yml ]
}
