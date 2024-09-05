process GET_READ_IDS_FROM_BAM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"


    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.txt"),     emit: readIDs,  optional: true

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools view --threads ${task.cpus-1} ${bam} | cut -f 1 >  ${prefix}_methylated_reads.txt
    """
}


