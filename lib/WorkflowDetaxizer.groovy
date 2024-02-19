//
// This file holds several functions specific to the workflow/detaxizer.nf in the nf-core/detaxizer pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowDetaxizer {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {

        genomeExistsError(params, log)


        if (!params.fasta && !params.genome) {
            Nextflow.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file. If no fasta is specified a genome has to be specified, e.g. '--genome GRCh38' "
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Generate methods description for MultiQC
    //

    public static String toolCitationText(params) {

        def citation_text = [
                "Tools used in the workflow included:",
                "FastQC (Andrews 2010),",
                "MultiQC (Ewels et al. 2016),",
                "Kraken2 (Wood et al. 2019),",
                !params["skip_blastn"] ? "BLAST (Altschul et al. 1990)," : "",
                !params["skip_blastn"] | params["enable_filter"] ? "seqkit (Shen et al. 2016)," : "",
                "fastp (Chen et al. 2018)",
                "."
            ].join(' ').trim()

        return citation_text
    }

    public static String toolBibliographyText(params) {

        def reference_text = [
                "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
                "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>",
                "<li>Wood, D. E., Lu, J. & Langmead, B. (2019) Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257. doi: 10.1186/s13059-019-1891-0</li>",
                !params["skip_blastn"] ? "<li>Altschul, S. F., Gish, W., Miller, W., Myers, E. W. & Lipman, D. J. (1990) Basic local alignment search tool. Journal of Molecular Biology 215, 403–410. doi: 10.1016/s0022-2836(05)80360-2.</li>" : "",
                !params["skip_blastn"] | params["enable_filter"] ? "<li>Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. In Q. Zou (Ed.), PLOS ONE (Vol. 11, Issue 10, p. e0163962). Public Library of Science (PLoS). doi: 10.1371/journal.pone.0163962</li>" : "",
                "<li>Chen, S., Zhou, Y., Chen, Y. & Gu, J. (2018) fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34, i884–i890. doi: 10.1093/bioinformatics/bty560</li>"
            ].join(' ').trim()

        return reference_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml, params) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        // Pipeline DOI
        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        // Tool references
        meta["tool_citations"] = ""
        meta["tool_bibliography"] = ""

        meta["tool_citations"] = toolCitationText(params).replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
        meta["tool_bibliography"] = toolBibliographyText(params)


        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.error(error_string)
        }
    }
}
