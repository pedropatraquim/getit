#!/usr/bin/env nextflow



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

id = Channel.fromPath(params.ids).
splitText(){it.replaceFirst(/\n/, "")}.
ifEmpty { error "Cannot find file matching: ${params.ids}" }

process read_file {
  tag "read the file"
  echo true
  publishDir "results/" 

  input:
    val id 

  output:
    file 'test_*.txt' into files

  script:
  """
  echo "pepe $id" > test_${id}.txt
  """
}

process copy_file {
  tag "copy files"

  input:
    file file from files
  

  script:
  """
   cp $file ~/${file}_2.txt
  """
}

// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
}


