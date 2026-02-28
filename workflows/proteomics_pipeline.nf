#!/usr/bin/env nextflow

/*
 * OmniOmics AI - Proteomics Analysis Pipeline
 * Nextflow workflow for automated proteomics analysis
 */

nextflow.enable.dsl = 2

// Parameters
params.input = "${params.data_dir}/*.txt"
params.output_dir = "${params.results_dir}"
params.design_file = "${params.data_dir}/design.csv"
params.normalization = "median"
params.llm_enabled = true

// Include subworkflows
include { proteomics_qc } from './modules/proteomics_qc'
include { differential_expression } from './modules/differential_expression'
include { enrichment_analysis } from './modules/enrichment'
include { llm_interpretation } from './modules/llm_interpretation'

workflow {
    main:
        // Channel for input files
        input_ch = Channel.fromPath(params.input)
        
        // Step 1: Quality Control
        proteomics_qc(input_ch)
        
        // Step 2: Differential Expression
        differential_expression(
            proteomics_qc.out.filtered_data,
            params.design_file
        )
        
        // Step 3: Enrichment Analysis
        enrichment_analysis(differential_expression.out.de_genes)
        
        // Step 4: LLM Interpretation (optional)
        if (params.llm_enabled) {
            llm_interpretation(
                differential_expression.out.results,
                enrichment_analysis.out.pathways
            )
        }
        
    emit:
        de_results = differential_expression.out.results
        pathways = enrichment_analysis.out.pathways
        llm_results = llm_interpretation.out.interpretations
}
