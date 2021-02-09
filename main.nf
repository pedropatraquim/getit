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
      --db                          path to main db
      --libs                        file containing the name of each libs and its path, comma separated, one per line
    Other options:
      --outdir                      The output directory where the results will be saved (default: $outdir)
      -w/--work-dir                 The temporary directory where intermediate data will be saved
      -profile                      Configuration profile to use. [standard, other_profiles] (default 'standard')
    """.stripIndent()
}


id = Channel.fromPath(params.ids).
  splitText(){it.replaceFirst(/\n/, "")}.
  take( params.sample ). // for testing pourposes
  ifEmpty { error "No ids at  ${params.ids}" }

libs = Channel.fromPath(params.libs).
  splitCsv().
  //collate( 2 ).
  ifEmpty { error "No libs at ${params.libs}" }

//id_libs=id.combine(libs).flatten().collate(3)

process get_queries {
  tag "getting sequence for query $id"
  //echo true
  publishDir "$params.outdir/$id/"
  //errorStrategy 'ignore'
  cpus 1

  input:
    val id

  output:
    tuple val(id), path('*.fasta') optional true into seqsOut

  script:
  """
  samtools faidx $params.db $id -o ${id}.fasta
  """
}

process jackhammer {

  tag "perform jackhammers for $id at $lib"
  //echo true
  publishDir "$params.outdir/$id/$lib/"
  cpus 2


  input:
    tuple val(id), path('sequence.fasta'),val(lib),path("lib.fasta") from seqsOut.combine(libs)

  output:
    tuple val(id),val(lib),path ('lib.fasta'),file("*_doms.txt") into jackhmmerTables
    path "*_table.txt"
    path "*_out.txt"
  script:
  """
   jackhmmer  --seed 1 --cpu $task.cpus -o  $id'_'$lib'_out.txt' --tblout $id'_'$lib'_table.txt' --domtblout $id'_'$lib'_doms.txt' sequence.fasta lib.fasta
  """
}

process get_hits {
  tag "get hits from first jackhammer  for $id at $lib"
  //echo true
  publishDir "$params.outdir/$id/$lib/"
  //errorStrategy 'ignore'
  cpus 1

  input:
    tuple val(id),val(lib),path ('lib.fasta'),path("doms.txt") from jackhmmerTables
  output:
    tuple val(id),val(lib),path ('lib.fasta'),path("*_hits.fasta") optional true into fastas

  script:
  """
    grep -v "#" doms.txt  | cut -f1 -d" "| sort | uniq | grep . | cat > ${lib}_ids.txt
    if test -s ${lib}_ids.txt; then
      samtools faidx lib.fasta -r ${lib}_ids.txt -o ${lib}_hits.fasta;
    fi
  """
}



process reciprocal_jackhammers {
    cpus 2
    tag "run reciprocal jackhammers for $id at $lib "
    //echo true
    publishDir "$params.outdir/$id/$lib/hits/"
    //validExitStatus 0,1

    input:
      tuple val(id),val(lib),path ('lib.fasta'),path("hits.fasta") from fastas
    output:
      tuple val(id),val(lib),path ('lib.fasta'),path("hits_doms.txt") into reciprocal_jackhammers

  script:
  """
    jackhmmer  --seed 1 --cpu $task.cpus -o hits_out.txt --domtblout hits_doms.txt hits.fasta $params.db
  """
}

process find_reciprocal_hits{
  tag "find reciprocal hits for $id at $lib"
  //echo true
  publishDir "$params.outdir/$id/$lib/hits/"
  errorStrategy 'ignore'
  cpus 1


  input:
    tuple val(id),val(lib),val(lib_fasta),path("hits_doms.txt") from reciprocal_jackhammers
  output:
    path("reciprocal_hits.txt") into reciprocal_hits

  script:
  """
    grep $id hits_doms.txt | tr -s " " "\t" | cut -f 4 | sort | uniq | grep . | perl -ne 'chomp;print "\$_,${id},${lib},${lib_fasta}\n"' > reciprocal_hits.txt
  """
}

process get_reciprocal_seqs {
  tag "get sequences of reciprocal hits $hit for $id at $lib"
  //echo true
  publishDir "$params.outdir/$id/$lib/hits/"
  cpus 1

  input:
    tuple val(hit),val(id),val(lib),path ('lib.fasta') from reciprocal_hits.splitCsv()

  output:
    tuple val(hit),val(id),val(lib),path ("*.fasta") into reciprocal_hits_seqs

  script:
  """
    samtools faidx lib.fasta $hit -o ${hit}.fasta
  """
}


process alingments {
  tag "run alignments for $id with $hit at $lib"
  //echo true
  publishDir "$params.outdir/$id/$lib/hits/$hit"
  //errorStrategy 'finish'
  cpus 8

  input:
    tuple val(hit),val(id),val(lib),path (hit_fasta) from reciprocal_hits_seqs

  output:
    tuple val(hit),val(id),val(lib),path("*.aln") into alignments

  script:

  """
    cat $launchDir/$params.outdir/$id/${id}.fasta  $hit_fasta  > aln.fasta ;
    mafft --thread $task.cpus --clustalout --reorder --maxiterate 1000 --retree 1 --globalpair aln.fasta > ${hit}.aln
    rm aln.fasta;
  """
}



process score_alignments {

  tag "score alignments for $id and $hit from $lib"
  publishDir "$params.outdir/$id/$lib/hits/$hit"
  cpus 1

  input:
    tuple val(hit),val (id),val(lib),path(aln_file) from alignments
  output:
    tuple val(id), path( "*_score.txt") into scores
    path("*_lenght.txt")
  script:
  """
    #!/usr/bin/env perl

    open FILE,"$aln_file";
    \$seq_length=`grep $hit $aln_file | tr  -s " " "\\t"  | cut -f2  | tr -d "-" | tr -d "\\n" | wc -c`;
    open L,">${hit}_lenght.txt";
    print L "\$seq_length\n";
    open OUT,">${hit}_score.txt";
    while(<FILE>){
      next unless /[\\*:.]/;
      \$asterisk=  tr/\\*//;
      \$colon=  tr/://;
      \$dot=  tr/\\.//;
      \$score=(\$asterisk*100+\$colon*50+\$dot*40)/\$seq_length;
    }
    print OUT "$lib\t$hit\t\$score\\n";
  """
}

process collect_scores {
  tag "collect scores in single file for query $id"
  publishDir "$params.outdir/$id/"
  cpus 1

  input:
     tuple val(id),path(score_file) from scores.groupTuple()
  output:
     path "${id}_scores.txt"
  script:
  """
   cat $score_file >> ${id}_scores.txt
  """
}



// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
}
