#!/usr/bin/env nextflow

params.conda = "$moduleDir/environment.yml"

process split_by_chr {
    label "high_resource"
    conda "${params.conda}"
    tag "${prefix}:${chrom}"

    input:
        tuple val(group), val(prefix), path(bam_file), path(bam_file_index), val(chrom)
    
    output:
        tuple val(group), val(chrom), path(name)

    script:
    name = "${prefix}.aligned.${chrom}.bam"
    """
    samtools view \
        -@ ${task.cpus} \
        -b \
        -o ${name} \
        ${bam_file} \
        ${chrom}
    """
}

process merge_bams {
    label "high_resource"
    conda "${params.conda}"
    tag "${group}:${chrom}"
    publishDir "${params.outdir}/merged/${group}"

    input:
        tuple val(group), val(chrom), path(bam_files)
    
    output:
        tuple val(group), val(chrom), path(name), path("${name}.bai")
    
    script:
    name = "${group}.${chrom}.aligned.bam"
    if bam_files.size() == 1: {
        """
        ln -s ${bam_files[0]} ${name}
        samtools index \
            -@ ${task.cpus} \
            ${name}
        """
    } else: {
        """
        samtools merge \
            -@ ${task.cpus} \
            -o ${name} \
            ${bam_files}

        samtools index \
            -@ ${task.cpus} \
            ${name}
        """
    }

}


workflow {
    chroms = Channel.fromPath(params.chrom_sizes)
        | splitCsv(header: false, sep: "\t")
        | map(row -> row[0])

    Channel.fromPath(params.samples_file)
        | splitCsv(header: true, sep: "\t")
        | map(row -> tuple(
                row.group,
                row.prefix,
                file(row.bam),
                file(row.bam_index)
            )
        )
        | combine(chroms)
        | split_by_chr
        | groupTuple([0, 1])
        | merge_bams
}