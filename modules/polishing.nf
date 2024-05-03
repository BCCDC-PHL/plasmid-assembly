process medaka {

    tag { sample_id }

    input:
    tuple val(sample_id), path(reads), path(assembly)

    output:
    tuple val(sample_id), path("medaka_output/${sample_id}_plassembler_hybrid.fa")

    script:
    """
    medaka_consensus \
	-t ${task.cpus} \
	-m ${params.medaka_model} \
	-i ${reads[2]} \
	-d ${assembly} \
	-o medaka_output

    mv medaka_output/consensus.fasta medaka_output/${sample_id}_plassembler_hybrid.fa
    """

}

process polypolish {

    tag { sample_id }

    input:
    tuple val(sample_id), path(reads), path(assembly)

    output:
    tuple val(sample_id), path("polypolish_output/${sample_id}_plassembler_hybrid.fa")

    script:
    """
    bwa index ${assembly}

    bwa mem \
	-t ${task.cpus} \
	-a \
	${assembly} \
	${reads[0]} \
	> ${sample_id}_R1.sam

    bwa mem \
	-t ${task.cpus} \
	-a \
	${assembly} \
	${reads[1]} \
	> ${sample_id}_R2.sam

    polypolish filter \
	--in1 ${sample_id}_R1.sam \
	--in2 ${sample_id}_R2.sam \
	--out1 ${sample_id}_filtered_R1.sam \
	--out2 ${sample_id}_filtered_R2.sam

    mkdir -p polypolish_output

    polypolish polish \
	${assembly} \
	${sample_id}_filtered_R1.sam \
	${sample_id}_filtered_R2.sam \
	> polypolish_output/${sample_id}_plassembler_hybrid.fa
    """

}
