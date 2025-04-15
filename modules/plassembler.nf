process plassembler {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_plassembler*", mode: 'copy'
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_flye*", mode: 'copy'
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_unicycler*", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), path(database)

    output:
    tuple val(sample_id), path("${sample_id}_plassembler_hybrid.fa"),                                          emit: plassembler_assembly
    tuple val(sample_id), path("${sample_id}_plassembler_hybrid_chromosome.fa"),                               emit: plassembler_chromosome
    tuple val(sample_id), path("${sample_id}_plassembler_hybrid_plasmids.fa"),                               emit: plassembler_plasmids
    tuple val(sample_id), path("${sample_id}_plassembler_summary.tsv"),                                        emit: plassembler_summary
    tuple val(sample_id), path("${sample_id}_plassembler.log"),                                                emit: plassembler_log
    tuple val(sample_id), path("${sample_id}_flye_long.fa"), path("${sample_id}_flye_long.gfa"),               emit: flye_assembly
    tuple val(sample_id), path("${sample_id}_unicycler_hybrid.fa"), path("${sample_id}_unicycler_hybrid.gfa"), emit: unicycler_assembly, optional: true
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

    plassembler run \
        --threads ${task.cpus} \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -l ${reads[2]} \
        -d ${database} \
	-c ${params.chromosome_length} \
	--keep_chromosome \
        --outdir ${sample_id}_assembly

    cat ${sample_id}_assembly/chromosome.fasta ${sample_id}_assembly/plassembler_plasmids.fasta > ${sample_id}_plassembler_hybrid.fa
    cp ${sample_id}_assembly/chromosome.fasta ${sample_id}_plassembler_hybrid_chromosome.fa
    cp ${sample_id}_assembly/plassembler_plasmids.fasta ${sample_id}_plassembler_hybrid_plasmids.fa
    cp ${sample_id}_assembly/plassembler_summary.tsv ${sample_id}_plassembler_summary.tsv
    cp ${sample_id}_assembly/plassembler_*.log ${sample_id}_plassembler.log

    cp ${sample_id}_assembly/flye_output/assembly.fasta ${sample_id}_flye_long.fa
    cp ${sample_id}_assembly/flye_output/assembly_graph.gfa ${sample_id}_flye_long.gfa

    if [ -f ${sample_id}_assembly/unicycler_output/assembly.fasta ]; then
        cp ${sample_id}_assembly/unicycler_output/assembly.fasta ${sample_id}_unicycler_hybrid.fa
        cp ${sample_id}_assembly/unicycler_output/assembly.gfa ${sample_id}_unicycler_hybrid.gfa
    fi
    """
}
