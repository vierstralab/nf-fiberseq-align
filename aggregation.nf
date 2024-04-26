#!/usr/bin/env nextflow

params.conda = "$moduleDir/environment.yml"

def set_key_for_group_tuple(ch) {
  ch.groupTuple()
    | map{ it -> tuple(groupKey(it[0], it[1].size()), *it[1..(it.size()-1)]) }
    | transpose()
}

process merge_bams {
    label "high_resource"
    conda "${params.conda}"
    tag "${group_key}"
    publishDir "${params.outdir}/merged/${groups[0]}"

    input:
        tuple val(group_key), path(bam_files), path(bam_indices)
    
    output:
        tuple val(group_key), path(name), path("${name}.bai")
    
    script:
    name = "${group_key}.aligned.bam"
    if (bam_files.size() == 1) {
        """
        ln -s ${bam_files[0]} ${name}
        samtools index \
            -@ ${task.cpus} \
            ${name}
        """
    } else {
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
    Channel.fromPath(params.samples_file)
        | splitCsv(header: true, sep: "\t")
        | map(row -> tuple(
            row.group,
            file(row.bam),
            file(row.bam_index ?: "${row.bam}.bai")
        )) // group, sample, bam, bam_index
        | set_key_for_group_tuple
        | groupTuple()
        | merge_bams
}
