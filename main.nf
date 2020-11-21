#!/usr/bin/env nextflow

params.outdir = "$baseDir/results"
params.help = ""

outdir = params.outdir



def helpMessage() {
    log.info"""
    ==================================================================
    ${workflow.manifest.name}  ~  version ${workflow.manifest.version}
    ==================================================================

    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run  -profile standard --outdir /output/path --ids ids_file.txt
    
    Options:
      --ids                         file containing the ids of the proteins, one per line
    Other options:
      --outdir                      The output directory where the results will be saved (default: $outdir)
      -w/--work-dir                 The temporary directory where intermediate data will be saved
      -profile                      Configuration profile to use. [standard, other_profiles] (default 'standard')
    """.stripIndent()
}

// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
}


