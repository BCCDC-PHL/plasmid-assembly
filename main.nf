#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { hash_files }                 from './modules/hash_files.nf'
include { fastp }                      from './modules/short_read_qc.nf'
include { fastp_json_to_csv }          from './modules/short_read_qc.nf'
include { filtlong }                   from './modules/long_read_qc.nf'
include { nanoq as nanoq_pre_filter }  from './modules/long_read_qc.nf'
include { nanoq as nanoq_post_filter } from './modules/long_read_qc.nf'
include { merge_nanoq_reports }        from './modules/long_read_qc.nf'
include { flye }                       from './modules/assembly.nf'
include { reorient_contigs }           from './modules/assembly.nf'
include { plassembler }                from './modules/assembly.nf'
include { combine_chromosome_and_plasmid_assemblies } from './modules/assembly.nf'
include { medaka as medaka_round_1 }   from './modules/polishing.nf'
include { medaka as medaka_round_2 }   from './modules/polishing.nf'
include { polypolish }                 from './modules/polishing.nf'
include { prokka }                     from './modules/annotation.nf'
include { bakta }                      from './modules/annotation.nf'
include { ale }                        from './modules/assembly_qc.nf'
include { quast }                      from './modules/assembly_qc.nf'
include { parse_quast_report }         from './modules/assembly_qc.nf'
include { bandage }                    from './modules/assembly_qc.nf'
include { pipeline_provenance }        from './modules/provenance.nf'
include { collect_provenance }         from './modules/provenance.nf'


workflow {

    ch_workflow_metadata = Channel.value([
	workflow.sessionId,
	workflow.runName,
	workflow.manifest.name,
	workflow.manifest.version,
	workflow.start,
    ])

    ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)

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

    ch_cleaned_reads = fastp.out.trimmed_reads.join(filtlong.out.filtered_reads).map{ it -> [it[0], [it[1], it[2], it[3]]] }

    flye(ch_cleaned_reads)

    ch_chromosome_assembly = flye.out.chromosome_assembly

    if (!params.skip_polishing && !params.skip_medaka) {
	ch_chromosome_assembly = medaka_round_1(ch_cleaned_reads.join(ch_chromosome_assembly))
    }

    reorient_contigs(ch_chromosome_assembly)

    ch_chromosome_assembly = reorient_contigs.out.reoriented_assembly

    if (!params.skip_polishing && !params.skip_medaka) {
	ch_chromosome_assembly = medaka_round_2(ch_cleaned_reads.join(ch_chromosome_assembly))
    }

    if (!params.skip_polishing && !params.skip_polypolish) {
	ch_chromosome_assembly = polypolish(ch_cleaned_reads.join(ch_chromosome_assembly))
    }

    ch_flye_outdir = flye.out.flye_outdir

    ale(ch_cleaned_reads.join(ch_chromosome_assembly))

    plassembler(ch_cleaned_reads.join(ch_flye_outdir).combine(ch_db))

    ch_plasmid_assembly = plassembler.out.plasmid_assembly

    ch_full_assembly = combine_chromosome_and_plasmid_assemblies(ch_chromosome_assembly.join(ch_plasmid_assembly))

    if (params.prokka) {
	prokka(ch_full_assembly)
    }

    if (params.bakta) {
	bakta(ch_full_assembly)
    }

    quast(ch_full_assembly)
    // bandage(plassembler.out.assembly_graph)

    parse_quast_report(quast.out.tsv)

    // Collect csv & tsv outputs from all samples into collection files
    if (params.collect_outputs) {
	fastp_json_to_csv.out.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "${params.collected_outputs_prefix}_fastp.csv", storeDir: "${params.outdir}")
	parse_quast_report.out.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "${params.collected_outputs_prefix}_quast.csv", storeDir: "${params.outdir}")
	if (params.hybrid || params.long_only) {
	    merge_nanoq_reports.out.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "${params.collected_outputs_prefix}_nanoq.csv", storeDir: "${params.outdir}")
	}
	flye.out.assembly_info.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "${params.collected_outputs_prefix}_flye_assembly_info.csv", storeDir: "${params.outdir}")
	flye.out.assembly_completeness.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "${params.collected_outputs_prefix}_flye_assembly_completeness.csv", storeDir: "${params.outdir}")
	plassembler.out.plassembler_summary.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "${params.collected_outputs_prefix}_plassembler_summary.tsv", storeDir: "${params.outdir}")
	ale.out.map{ it -> it[1] }.collectFile(keepHeader: true, sort: { it.text }, name: "${params.collected_outputs_prefix}_ale.csv", storeDir: "${params.outdir}")
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
