process prokka {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_plassembler_hybrid_prokka.{gbk,gff}", mode: 'copy'

    input:
      tuple val(sample_id), path(assembly)

    output:
      tuple val(sample_id), path("${sample_id}_plassembler_hybrid_prokka.gbk"),            emit: gbk
      tuple val(sample_id), path("${sample_id}_plassembler_hybrid_prokka.gff"),            emit: gff
      tuple val(sample_id), path("${sample_id}_plassembler_hybrid_prokka_provenance.yml"), emit: provenance

    script:
      """
      printf -- "- process_name: prokka\\n"                                          >> ${sample_id}_plassembler_hybrid_prokka_provenance.yml
      printf -- "  tools:\\n"                                                        >> ${sample_id}_plassembler_hybrid_prokka_provenance.yml
      printf -- "    - tool_name: prokka\\n"                                         >> ${sample_id}_plassembler_hybrid_prokka_provenance.yml  
      printf -- "      tool_version: \$(prokka --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_plassembler_hybrid_prokka_provenance.yml  
      printf -- "      parameters:\\n"                                               >> ${sample_id}_plassembler_hybrid_prokka_provenance.yml
      printf -- "        - parameter: --compliant\\n"                                >> ${sample_id}_plassembler_hybrid_prokka_provenance.yml
      printf -- "          value: null\\n"                                           >> ${sample_id}_plassembler_hybrid_prokka_provenance.yml

      prokka --cpus ${task.cpus} --compliant --locustag ${sample_id} --centre "BCCDC-PHL" --prefix "${sample_id}" ${assembly}

      cp ${sample_id}/${sample_id}.gbk ${sample_id}_plassembler_hybrid_prokka.gbk
      cp ${sample_id}/${sample_id}.gff ${sample_id}_plassembler_hybrid_prokka.gff
      """
}
