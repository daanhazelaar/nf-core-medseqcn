process SIZE_SELECTION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"


    input:
    tuple val(meta), path(bam), path(index)
    val(min)
    val(max)

    output:
    tuple val(meta), path("*.bam"),     emit: bam,  optional: true

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-fmt sam") ? "sam" :
                    args.contains("--output-fmt cram") ? "cram" :
                    "bam"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    samtools view -h --threads ${task.cpus-1} ${bam} |
    awk 'substr(\$0,1,1)=="@" || (\$9>=$min && \$9<=$max) || (\$9<=-$min && \$9>=-$max)' |
    samtools view -b > ${prefix}.${extension}
    """
}
