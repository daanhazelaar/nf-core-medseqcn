process SUBSAMPLE_BAM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' :
        'biocontainers/samtools:1.20--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(bai), val(fraction)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml",            emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Groovy conditional: evaluated by Nextflow before the script is submitted,
    // so no shell arithmetic tools (bc, awk) are required inside the container.
    if ( fraction < 1.0 ) {
        // samtools view -s seed.fraction (e.g. 42.75 = seed 42, keep 75%)
        def subsample_seed_frac = 42 + fraction
        """
        samtools view \\
            -s ${subsample_seed_frac} \\
            -b \\
            -o ${prefix}_equalized.bam \\
            ${bam}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else {
        """
        cp ${bam} ${prefix}_equalized.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
}
