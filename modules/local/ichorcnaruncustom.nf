process ICHORCNA_RUN_CUSTOM {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ichorcna:0.3.2--pl5321r42hdfd78af_2' :
        'biocontainers/r-ichorcna:0.3.2--pl5321r42hdfd78af_2' }"

    input:
    tuple val(meta), path(wig)
    path gc_wig
    path map_wig
    path panel_of_normals
    path centromere
    tuple val(meta), val(sex)


    output:
    // tuple val(meta), path("${meta.id}.cna.seg")     , emit: seg_cna
    // tuple val(meta), path("${meta.id}.seg.txt")     , emit: seg_txt
    // tuple val(meta), path("${meta.id}.seg")         , emit: seg
    // tuple val(meta), path("${meta.id}.params.txt")  , emit: ichorcna_params
    path "versions.yml"                             , emit: versions
    path "${meta.id}_ichorCNA_files"                               , emit: outputDir
    // path "${meta.id}.RData"                         , emit: rdata_object
    // path "${meta.id}.ALL_RESULTS.RData"             , emit: all_results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pon = panel_of_normals ? "--normalPanel ${panel_of_normals}" : ''
    def centro = centromere ? "--centromere ${centromere}" : ''
    def sex = sex ? "--sex ${sex}" : ''
    def VERSION = '0.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    runIchorCNA_custom.R \\
        $args \\
        --WIG ${wig} \\
        --id ${prefix} \\
        --gcWig ${gc_wig} \\
        --mapWig ${map_wig} \\
        ${pon} \\
        ${centro} \\
        ${sex} \\
        --outDir .

    mkdir ./${meta.id}_ichorCNA_files
    cp -r ${meta.id}.ALL_RESULTS.RData ./${meta.id}_ichorCNA_files
    cp -r ${meta.id}.RData ./${meta.id}_ichorCNA_files
    cp -r ${meta.id}.seg ./${meta.id}_ichorCNA_files
    cp -r ${meta.id}/ ./${meta.id}_ichorCNA_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ichorcna: $VERSION
    END_VERSIONS
    """
}

