#!/usr/bin/env nextflow

params.conda = "$moduleDir/environment.yml"

process call_5mc {
    label "high_resource"
    tag "${prefix}"
    conda "${params.conda}"

    input:
        tuple val(prefix), path(bam)

    output:
        tuple val(prefix), path(name)

    script:
    name = "${prefix}.predicted_5mc.bam"
    """
    jasmine --keep-kinetics \
        --num-threads ${task.cpus} \
        --log-level INFO \
        ${bam} ${name}
    """
}
process call_m6a {
    label "high_resource"
    tag "${prefix}"
    module "fibertools-rs/0.3.2"
    publishDir "${params.outdir}/predicted_m6a"

    input:
        tuple val(prefix), path(bam)

    output:
        tuple val(prefix), path(name)

    script:
    name = "${prefix}.predicted_m6a.bam"
    """
    ft predict-m6a -t ${task.cpus} -v ${bam} ${name}
    """
}

process align_bam {
    label "high_resource"
    conda "${params.conda}"
    tag "${prefix}"
    publishDir "${params.outdir}/aligned_bam"

    input:
        tuple val(prefix), path(bam)

    output:
        tuple val(prefix), path(name), path("${name}.bai")

    script:
    name = "${prefix}.aligned.bam"
    """
    pbmm2 align \
        --preset CCS \
        --sort \
        --num-threads ${task.cpus} \
        --bam-index BAI \
        --sort-threads ${task.cpus} \
        ${params.genome_fasta_file} \
        ${bam} \
        ${name}
    """
}

process extract_signal {
    tag "${prefix}"
    publishDir "${params.outdir}/per_fiber_signal"
    conda "${params.conda}"
    scratch true
    cpus 8

    input:
        tuple val(prefix), path(bam), path(bam_index)
    output:
        tuple val(prefix), path(name)

    script:
    name = "${prefix}.signal.bed.gz"
    """
    ft extract ${bam} -t ${task.cpus} --all ${name}
    """
}

process to_ucsc_format {
    tag "${prefix}"
    publishDir "${params.outdir}/ucsc_format"
    conda "${params.conda}"

    input:
        tuple val(prefix), path(per_fiber_signal)
    output:
        tuple val(prefix), path(name), path("${name}.tbi")
    
    script:
    name = "${prefix}.ucsc.bed.gz"
    """
    zcat ${per_fiber_signal} \
        | tail -n +2 \
        | awk -v OFS='\t' -F'\t' \
            '{ \
                split(\$28, positions, ","); \
                for (i in positions) { \
                    if (positions[i] != "" && positions[i] != "-1" && positions[i] != ".") { \
                        print \$1, positions[i], positions[i] + 1, \$6, \$4  \
                    } \
                } \
            }'  \
        | sort-bed - > tmp.bed
    
    bedtools getfasta \
        -fi ${params.genome_fasta_file} \
        -bed tmp.bed \
        -bedOut \
        | bgzip -c > ${name}

    tabix ${name}
    """
}

workflow {
    data = Channel.fromPath(params.samples_file)
        | splitCsv(header: true, sep: "\t")
        | map(row -> tuple(row.sample_id, file(row.reads), row['Instrument']))
        | branch {
            with_5mc: it[2] == "Revio"
            no_5mc: true
        }
    data.no_5mc
        | map(it -> tuple(it[0], it[1]))
        | call_5mc
        | mix(
            data.with_5mc | map(it -> tuple(it[0], it[1]))
        )
        | call_m6a
        | align_bam
        | extract_signal
        | to_ucsc_format
}

workflow extractSignal {
    Channel.fromPath(params.samples_file)
        | splitCsv(header: true, sep: "\t")
        | map(row -> tuple(row.sample_id, file(row.bam), file(row.bam_index)))
        | extract_signal
        | to_ucsc_format
}
