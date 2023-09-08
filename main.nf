#!/usr/bin/env nextflow

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --bams sample.bam [Options]
    
    Inputs Options:
    --bam        Input file
    --bam_index     Input BAM Index file
    --fasta         Input genome ref
    --index         Input genome Index file
   
Resource Options:
    --max_cpus      Maximum number of CPUs (int)
                    (default: $params.max_cpus)  
    --max_memory    Maximum memory (memory unit)
                    (default: $params.max_memory)
    --max_time      Maximum time (time unit)
                    (default: $params.max_time)
    See here for more info: https://github.com/lifebit-ai/hla/blob/master/docs/usage.md
    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

// Define channels from repository files
projectDir = workflow.projectDir
params.outdir = "./results"

// Define Channels from input
Channel
    .fromPath(params.bam)
    .ifEmpty { exit 1, "Cannot find input BAM file : ${params.bam}" }
    .set { ch_input }

Channel
    .fromPath(params.bam_index)
    .ifEmpty { exit 1, "Cannot find input BAM file : ${params.bam}" }
    .set { ch_bindex }
    
Channel
    .fromPath(params.fasta)
    .ifEmpty { exit 1, "Cannot find input reference file : ${params.fasta}" }
    .set { ch_ref }
    
Channel
    .fromPath(params.index)
    .ifEmpty { exit 1, "Cannot find input reference file : ${params.index}" }
    .set { ch_ind}

// Define Process
process gatk {
    tag "$bam_file"
    label 'gatk'
    publishDir "${params.outdir}/gatk", mode: 'copy'
 
    input:
    file(bam_file) from ch_input
    file(bam_index) from ch_bindex
    file(fasta) from ch_ref
    file(index) from ch_ind
      
    output:
    file "output.vcf" into ch_out
    
    script:
    """    
    
    gatk --java-options "-Xmx4g" HaplotypeCaller \
     -R ${fasta} -I $bam_file -O output.vcf
     
     """
  }

ch_report_dir = Channel.value(file("${projectDir}/bin/report"))

process report {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file(report_dir) from ch_report_dir
    file(table) from ch_out
    
    output:
    file "multiqc_report.html" into ch_multiqc_report

    script:
    """
    cp -r ${report_dir}/* .
    Rscript -e "rmarkdown::render('report.Rmd',params = list(res_table='$table'))"
    mv report.html multiqc_report.html
    """
}
