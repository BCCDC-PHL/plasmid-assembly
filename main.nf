#!/usr/bin/env nextflow

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { hash_files }                 from './modules/hash_files.nf'
include { fastp }                      from './modules/fastp.nf'
include { fastp_json_to_csv }          from './modules/fastp.nf'
include { filtlong }                   from './modules/long_read_qc.nf'
include { nanoq as nanoq_pre_filter }  from './modules/long_read_qc.nf'
include { nanoq as nanoq_post_filter } from './modules/long_read_qc.nf'
include { merge_nanoq_reports }        from './modules/long_read_qc.nf'
include { plassembler }                from './modules/plassembler.nf'
include { prokka }                     from './modules/prokka.nf'
include { bakta }                      from './modules/bakta.nf'
include { quast }                      from './modules/quast.nf'
include { parse_quast_report }         from './modules/quast.nf'
include { bandage }                    from './modules/long_read_qc.nf'
include { pipeline_provenance }        from './modules/provenance.nf'
include { collect_provenance }         from './modules/provenance.nf'


workflow {
    ch_start_time = Channel.of(LocalDateTime.now())
    ch_pipeline_name = Channel.of(workflow.manifest.name)
    ch_pipeline_version = Channel.of(workflow.manifest.version)

    ch_pipeline_provenance = pipeline_provenance(ch_pipeline_name.combine(ch_pipeline_version).combine(ch_start_time))

    if (params.samplesheet_input != 'NO_FILE') {
	ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], [it['R1'], it['R2'], it['LONG']]] }
	ch_short_reads = ch_fastq.map{ it -> [it[0], [it[1][0], it[1][1]]] }
	ch_long_reads = ch_fastq.map{ it -> [it[0], it[1][2]] }    
    } else {
	ch_short_reads = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], [it[1], it[2]]] }.unique{ it -> it[0] }
	ch_long_reads = Channel.fromPath( params.long_reads_search_path ).map{ it -> [it.baseName.split("_")[0], [it]] }
	ch_fastq = ch_short_reads.join(ch_long_reads).map{ it -> [it[0], it[1] + it[2]] }
    }

    ch_db = Channel.fromPath(params.db)

    main:
    ch_provenance = ch_fastq.map{ it -> it[0] }

    hash_files(ch_fastq.combine(Channel.of("fastq-input")))

    fastp(ch_short_reads)
    fastp_json_to_csv(fastp.out.json)
    nanoq_pre_filter(ch_long_reads.combine(Channel.of("pre_filter")))
    filtlong(ch_long_reads)
    nanoq_post_filter(filtlong.out.filtered_reads.combine(Channel.of("post_filter")))
    merge_nanoq_reports(nanoq_pre_filter.out.report.join(nanoq_post_filter.out.report))
    plassembler(fastp.out.trimmed_reads.join(filtlong.out.filtered_reads).map{ it -> [it[0], [it[1], it[2], it[3]]] }.combine(ch_db))
    ch_assembly = plassembler.out.plassembler_assembly

    if (params.prokka) {
	prokka(ch_assembly)
    }

    if (params.bakta) {
	bakta(ch_assembly)
    }

    quast(ch_assembly)
    // bandage(plassembler.out.assembly_graph)

    parse_quast_report(quast.out.tsv)

    // Collect csv & tsv outputs from all samples into collection files
    if (params.collect_outputs) {
	fastp_json_to_csv.out.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "${params.collected_outputs_prefix}_fastp.csv", storeDir: "${params.outdir}")
	parse_quast_report.out.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "${params.collected_outputs_prefix}_quast.csv", storeDir: "${params.outdir}")
	if (params.hybrid || params.long_only) {
	    merge_nanoq_reports.out.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "${params.collected_outputs_prefix}_nanoq.csv", storeDir: "${params.outdir}")
	}
    }

    //
    // Provenance collection processes
    // The basic idea is to build up a channel with the following structure:
    // [sample_id, [provenance_file_1.yml, provenance_file_2.yml, provenance_file_3.yml...]]
    // ...and then concatenate them all together in the 'collect_provenance' process.
    ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it -> [it[0], [it[1]]] }
    ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

    ch_provenance = ch_provenance.join(nanoq_pre_filter.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(filtlong.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(nanoq_post_filter.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

    ch_provenance = ch_provenance.join(plassembler.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

    if (params.prokka) {
	ch_provenance = ch_provenance.join(prokka.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    }

    if (params.bakta) {
	ch_provenance = ch_provenance.join(bakta.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    }

    ch_provenance = ch_provenance.join(quast.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

    collect_provenance(ch_provenance)
}
