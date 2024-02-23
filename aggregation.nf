#!/usr/bin/env nextflow

params.conda = "$moduleDir/environment.yml"

def set_key_for_group_tuple(ch) {
  ch.groupTuple()
    | map{ it -> tuple(groupKey(it[0], it[1].size()), *it[1..(it.size()-1)]) }
    | transpose()
}

process split_by_chr {
    label "high_resource"
    conda "${params.conda}"
    tag "${prefix}:${chrom}"

    input:
        tuple val(group_key), val(group), val(prefix), path(bam_file), path(bam_file_index), val(chrom)
    
    output:
        tuple val(group_key), val(group), path(name)

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
    tag "${group_key}"
    publishDir "${params.outdir}/merged/${groups[0]}"

    input:
        tuple val(group_key), val(groups), path(bam_files)
    
    output:
        tuple val(group_key), path(name), path("${name}.bai")
    
    script:
    name = "${group_key}.aligned.bam"
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
                row.sample_id,
                file(row.bam),
                file(row.bam_index)
            )
        ) // group, sample, bam, bam_index
        | combine(chroms) // group, sample, bam, bam_index, chrom
        | map(it -> tuple("${it[0]}.${it[4]}", *it)) // new_id, group, sample, bam, bam_index, chrom
        | set_key_for_group_tuple
        | split_by_chr
        | groupTuple([0, 1])
        | merge_bams
}