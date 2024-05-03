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


process bakta {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}*.{gbk,gff,json,log}", mode: 'copy'

    input:
      tuple val(sample_id), path(assembly)

    output:
      tuple val(sample_id), path("${sample_id}_plassembler_hybrid_bakta.gbk"),            emit: gbk
      tuple val(sample_id), path("${sample_id}_plassembler_hybrid_bakta.gff"),            emit: gff
      tuple val(sample_id), path("${sample_id}_plassembler_hybrid_bakta.json"),           emit: json
      tuple val(sample_id), path("${sample_id}_plassembler_hybrid_bakta.log"),            emit: log
      tuple val(sample_id), path("${sample_id}_plassembler_hybrid_bakta_provenance.yml"), emit: provenance

    script:
      """
      printf -- "- process_name: bakta\\n"                                     >> ${sample_id}_plassembler_hybrid_bakta_provenance.yml
      printf -- "  tools:\\n"                                                  >> ${sample_id}_plassembler_hybrid_bakta_provenance.yml
      printf -- "    - tool_name: bakta\\n"                                    >> ${sample_id}_plassembler_hybrid_bakta_provenance.yml
      printf -- "      tool_version: \$(bakta --version | cut -d ' ' -f 2)\\n" >> ${sample_id}_plassembler_hybrid_bakta_provenance.yml
      printf -- "      parameters:\\n"                                         >> ${sample_id}_plassembler_hybrid_bakta_provenance.yml
      printf -- "        - parameter: --db\\n"                                 >> ${sample_id}_plassembler_hybrid_bakta_provenance.yml
      printf -- "          value: ${params.bakta_db}\\n"                       >> ${sample_id}_plassembler_hybrid_bakta_provenance.yml
      printf -- "        - parameter: --compliant\\n"                          >> ${sample_id}_plassembler_hybrid_bakta_provenance.yml
      printf -- "          value: null\\n"                                     >> ${sample_id}_plassembler_hybrid_bakta_provenance.yml
      printf -- "        - parameter: --keep-contig-headers\\n"                >> ${sample_id}_plassembler_hybrid_bakta_provenance.yml
      printf -- "          value: null\\n"                                     >> ${sample_id}_plassembler_hybrid_bakta_provenance.yml

      mkdir tmp

      bakta \
        --force \
        --threads ${task.cpus} \
        --tmp-dir ./tmp \
        --debug \
        --db ${params.bakta_db} \
        --compliant \
        --keep-contig-headers \
        --locus-tag ${sample_id} \
        --prefix "${sample_id}" \
        ${assembly}

      cp ${sample_id}.gff3 ${sample_id}_plassembler_hybrid_bakta.gff
      cp ${sample_id}.gbff ${sample_id}_plassembler_hybrid_bakta.gbk
      cp ${sample_id}.json ${sample_id}_plassembler_hybrid_bakta.json
      cp ${sample_id}.log  ${sample_id}_plassembler_hybrid_bakta.log
      """
}
