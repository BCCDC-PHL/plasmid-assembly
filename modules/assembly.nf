process flye {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_flye_hybrid_assembly_info.tsv", mode: 'copy'
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_flye_hybrid_assembly_completeness.csv*", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_flye_hybrid_chromosome.fa"),             emit: chromosome_assembly
    tuple val(sample_id), path("${sample_id}_flye_hybrid_assembly_info.csv"),         emit: assembly_info
    tuple val(sample_id), path("${sample_id}_flye_hybrid_assembly_completeness.csv"), emit: assembly_completeness
    tuple val(sample_id), path("${sample_id}_flye_output"),                           emit: flye_outdir

    script:
    """
    flye \
	--threads ${task.cpus} \
	--${params.flye_model} \
	${reads[2]} \
	--out-dir ${sample_id}_flye_output

    cp ${sample_id}_flye_output/assembly.fasta ${sample_id}_flye_hybrid.fa
    sed 's/#//g' ${sample_id}_flye_output/assembly_info.txt > ${sample_id}_flye_hybrid_assembly_info.tsv

    check_assembly_completion.py \
	--sample-id ${sample_id} \
	--min-chromosome-contig-length ${params.chromosome_length} \
	${sample_id}_flye_hybrid_assembly_info.tsv \
	> ${sample_id}_flye_hybrid_assembly_completeness.csv

    extract_chromosome.py \
	--sample-id ${sample_id} \
	--min-chromosome-contig-length ${params.chromosome_length} \
	--assembly-info ${sample_id}_flye_hybrid_assembly_info.tsv \
	--assembly-input ${sample_id}_flye_hybrid.fa \
	--output-chromosome-fasta ${sample_id}_flye_hybrid_chromosome.fa \
        --output-assembly-info ${sample_id}_flye_hybrid_assembly_info.csv
    """
}

process reorient_contigs {

    tag { sample_id }

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_flye_hybrid_chromosome.fa"), emit: reoriented_assembly

    script:
    """
    dnaapler all \
	--threads ${task.cpus} \
	-i ${assembly} \
	--autocomplete nearest \
	--db dnaa,repa \
	--prefix ${sample_id} \
    	-o ${sample_id}_dnaapler_output

    mv ${assembly} input_assembly.fa
    cp ${sample_id}_dnaapler_output/${sample_id}_reoriented.fasta ${sample_id}_flye_hybrid_chromosome.fa
    """
}

process plassembler {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_plassembler*", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), path(flye_output), path(database)

    output:
    tuple val(sample_id), path("${sample_id}_plassembler_hybrid_plasmids.fa"),                                 emit: plasmid_assembly
    tuple val(sample_id), path("${sample_id}_plassembler_summary.tsv"),                                        emit: plassembler_summary
    tuple val(sample_id), path("${sample_id}_plassembler.log"),                                                emit: plassembler_log
    tuple val(sample_id), path("${sample_id}_plassembler_provenance.yml"),                                     emit: provenance

    script:
    """
    printf -- "- process_name: plassembler\\n"                                                 >> ${sample_id}_plassembler_provenance.yml
    printf -- "  tools:\\n"                                                                    >> ${sample_id}_plassembler_provenance.yml
    printf -- "    - tool_name: plassembler\\n"                                                >> ${sample_id}_plassembler_provenance.yml
    printf -- "      tool_version: \$(plassembler --version | cut -d ' ' -f 3)\\n"             >> ${sample_id}_plassembler_provenance.yml
    printf -- "      parameters:\\n"                                                           >> ${sample_id}_plassembler_provenance.yml
    printf -- "        - parameter: database\\n"                                               >> ${sample_id}_plassembler_provenance.yml
    printf -- "          value: \$(readlink ${database})\\n"                                   >> ${sample_id}_plassembler_provenance.yml
    printf -- "        - parameter: chromosome_length\\n"                                      >> ${sample_id}_plassembler_provenance.yml
    printf -- "          value: ${params.chromosome_length}\\n"                                >> ${sample_id}_plassembler_provenance.yml
    printf -- "        - parameter: skip_qc\\n"                                                >> ${sample_id}_plassembler_provenance.yml
    printf -- "          value: null\\n"                                                       >> ${sample_id}_plassembler_provenance.yml
    printf -- "        - parameter: keep_chromosome\\n"                                        >> ${sample_id}_plassembler_provenance.yml
    printf -- "          value: null\\n"                                                       >> ${sample_id}_plassembler_provenance.yml

    plassembler run \
        --threads ${task.cpus} \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -l ${reads[2]} \
        -d ${database} \
	-c ${params.chromosome_length} \
	--skip_qc \
	--flye_directory ${flye_output} \
        --outdir ${sample_id}_plassembler_output

    cp ${sample_id}_plassembler_output/plassembler_plasmids.fasta ${sample_id}_plassembler_hybrid_plasmids.fa
    cp ${sample_id}_plassembler_output/plassembler_summary.tsv ${sample_id}_plassembler_summary.tsv
    cp ${sample_id}_plassembler_output/plassembler_*.log ${sample_id}_plassembler.log
    """
}


process combine_chromosome_and_plasmid_assemblies {

    tag { sample_id }

    executor 'local'

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_hybrid.fa", mode: 'copy'

    input:
    tuple val(sample_id), path(chromosome_assembly), path(plasmid_assembly)

    output:
    tuple val(sample_id), path("${sample_id}_hybrid.fa")

    script:
    """
    cat ${chromosome_assembly} ${plasmid_assembly} > ${sample_id}_hybrid.fa
    """
}
