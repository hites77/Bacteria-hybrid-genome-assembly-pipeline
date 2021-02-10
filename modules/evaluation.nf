nextflow.enable.dsl=2

include { getDirectory; mapToDirectory } from './commons.nf'

// output directories for each process relative to params.outdir
outdirs = {}

outdirs.evaluateChromosome = 'chromosome_eval/'
outdirs.evaluatePlasmid = 'plasmid_eval/'

outdirs.shortReadsCoverage = 'short_read_coverage/'
outdirs.longReadsCoverage = 'long_read_coverage/'
outdirs.prokkaAnnotate = 'prokka/'
outdirs.quastEvaluate = 'quast/'
outdirs.checkmEvaluate = 'checkm/'
outdirs.runKofamscan = 'kofamscan/'

process shortReadsCoverage {
    publishDir "${params.outdir}/${pubDirPrefix}/${outdirs.shortReadsCoverage}" , mode: 'copy', pattern: '{.command.sh,.command.log,.exitcode}', saveAs: { 'nextflow' + it }
    publishDir "${params.outdir}/${pubDirPrefix}/${outdirs.shortReadsCoverage}" , mode: 'copy'
    conda params.condaEnvsDir + '/urops-assembly'
    
    input:
    val pubDirPrefix
    path assemblyFa
    path illumina1Fq
    path illumina2Fq
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path 'stats.txt', emit: stats
    path 'histogram.txt', emit: histogram
    path 'bbmap_stderr.txt', emit: stderr

    script:
    """
    bbmap.sh in1=$illumina1Fq in2=$illumina2Fq ref=$assemblyFa nodisk \
                covstats=stats.txt covhist=histogram.txt \
            2> bbmap_stderr.txt
    """
}

process longReadsCoverage {
    publishDir "${params.outdir}/${pubDirPrefix}/${outdirs.longReadsCoverage}", mode: 'copy', pattern: '{.command.sh,.command.log,.exitcode}', saveAs: { 'nextflow' + it }
    publishDir "${params.outdir}/${pubDirPrefix}/${outdirs.longReadsCoverage}", mode: 'copy'
    conda params.condaEnvsDir + '/urops-assembly'
    
    input:
    val pubDirPrefix
    path assemblyFa
    path pacbioFq
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path 'stats.txt', emit: stats
    path 'histogram.txt', emit: histogram
    path 'pileup_stderr.txt', emit: stderr
    

    script:
    """
    minimap2 -a $assemblyFa $pacbioFq > mapping.sam
    pileup.sh in=mapping.sam out=stats.txt hist=histogram.txt 2> pileup_stderr.txt
    """
}

process prokkaAnnotate {
    publishDir "${params.outdir}/${pubDirPrefix}/${outdirs.prokkaAnnotate}", mode: 'copy', pattern: '{.command.sh,.command.log,.exitcode}', saveAs: { 'nextflow' + it }
    publishDir "${params.outdir}/${pubDirPrefix}", mode: 'copy'
    conda params.condaEnvsDir + '/urops-assembly'

    input:
    val pubDirPrefix
    path assemblyFa
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path outdirs.prokkaAnnotate + 'prokka.gff', emit: gff
    path outdirs.prokkaAnnotate + 'prokka.faa', emit: faa
    path outdirs.prokkaAnnotate + 'prokka.txt', emit: txt
    path outdirs.prokkaAnnotate + '*'

    script:
    """
    prokka --cpus $params.threads --outdir ${outdirs.prokkaAnnotate} --prefix prokka --addgenes --addmrna --compliant --rfam $assemblyFa
    """
}

process quastEvaluate {
    publishDir "${params.outdir}/${pubDirPrefix}/${outdirs.quastEvaluate}", mode: 'copy', pattern: '{.command.sh,.command.log,.exitcode}', saveAs: { 'nextflow' + it }
    publishDir "${params.outdir}/${pubDirPrefix}", mode: 'copy'
    conda params.condaEnvsDir + '/urops-assembly'
    
    input:
    val pubDirPrefix
    path assemblyFa
    path prokkaGff
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path outdirs.quastEvaluate + '*', emit: allFiles
    
    script:
    """
    quast $assemblyFa --circos -g $prokkaGff -t ${params.threads} --gene-finding --fragmented --conserved-genes-finding --rna-finding -o ${outdirs.quastEvaluate}
    """
}

process checkmEvaluate {
    publishDir "${params.outdir}/${pubDirPrefix}/${outdirs.checkmEvaluate}", mode: 'copy', pattern: '{.command.sh,.command.log,.exitcode}', saveAs: { 'nextflow' + it }
    publishDir "${params.outdir}/${pubDirPrefix}", mode: 'copy'
    conda params.condaEnvsDir + '/urops-checkm'
    
    input:
    val pubDirPrefix
    path assemblyFa
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path outdirs.checkmEvaluate + '*', emit: allFiles

    script:
    """
    mkdir input
    mv $assemblyFa input/assembly.fna
    checkm lineage_wf input ${outdirs.checkmEvaluate}
    """
}

process runKofamscan {
    conda params.condaEnvsDir + '/urops-kofamscan'
    publishDir "${params.outdir}/${pubDirPrefix}/${outdirs.runKofamscan}", mode: 'copy', pattern: '{.command.sh,.command.log,.exitcode}', saveAs: { 'nextflow' + it }
    publishDir "${params.outdir}/${pubDirPrefix}/${outdirs.runKofamscan}", mode: 'copy'


    input:
    val pubDirPrefix
    path prokkaFaa

    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path 'kofamscan_mapper.txt'
    path 'kofamscan_results.tsv'

    script:
    """
    # mapper format
    exec_annotation -f mapper -o kofamscan_mapper.txt -p ${params.kofam_profile} -k ${params.kofam_kolist} --cpu=${params.threads} $prokkaFaa

    # deatiled format (default)
    exec_annotation -o kofamscan_results.tsv -p ${params.kofam_profile} -k ${params.kofam_kolist} -f detail-tsv -E 0.01 --cpu=${params.threads} $prokkaFaa
    """
}

process makeChromosomeSummary {
    publishDir params.outdir, mode: 'copy'
    conda params.condaEnvsDir + '/urops-assembly'
    
    input:
    path shortReadsCoverageDir
    path longReadsCoverageDir
    path quastDir
    path prokkaTxt
    path assemblyFa
    path checkmDir
    
    output:
    path 'chromosome-summary.json'

    script:
    """
    chromosome_summary.py --short $shortReadsCoverageDir \
        --long $longReadsCoverageDir \
        --quast $quastDir \
        --prokka $prokkaTxt \
        --fasta $assemblyFa \
        --checkm $checkmDir \
        --out chromosome-summary.json
    """
}

process makePlasmidSummary {
    publishDir params.outdir, mode: 'copy'
    conda params.condaEnvsDir + '/urops-assembly'
    
    input:
    path shortReadsCoverageDir
    path longReadsCoverageDir
    path quastDir
    path prokkaTxt
    path plasmidFa
    
    output:
    path 'plasmid-summary.json'

    script:
    """
    plasmid_summary.py --short $shortReadsCoverageDir \
        --long $longReadsCoverageDir \
        --quast $quastDir \
        --prokka $prokkaTxt \
        --fasta $plasmidFa
        --out plasmid-summary.json
    """
}

workflow evaluateChromosome {
    take:
    pubDirPrefix
    assembly
    cleanedShortReads1
    cleanedShortReads2
    cleanedLongReads

    main:
    shortReadsCoverage(pubDirPrefix, assembly, cleanedShortReads1, cleanedShortReads2)
    longReadsCoverage(pubDirPrefix, assembly, cleanedLongReads)
    prokkaAnnotate(pubDirPrefix, assembly)
    quastEvaluate(pubDirPrefix, assembly, prokkaAnnotate.out.gff)
    checkmEvaluate(pubDirPrefix, assembly)
    makeChromosomeSummary(
        shortReadsCoverage.out.stats.map { file(it.parent) }, // HACK
        longReadsCoverage.out.stats.map { file(it.parent) }, // HACK
        mapToDirectory(quastEvaluate.out.allFiles),
        prokkaAnnotate.out.txt,
        assembly,
        mapToDirectory(checkmEvaluate.out.allFiles)
    )
}

workflow evaluatePlasmid {
    take:
    pubDirPrefix
    assembly
    cleanedShortReads1
    cleanedShortReads2
    cleanedLongReads

    main:
    shortReadsCoverage(pubDirPrefix, assembly, cleanedShortReads1, cleanedShortReads2)
    longReadsCoverage(pubDirPrefix, assembly, cleanedLongReads)
    prokkaAnnotate(pubDirPrefix, assembly)
    quastEvaluate(pubDirPrefix, assembly, prokkaAnnotate.out.gff)
    if (params.kofamscan) {
        runKofamscan(pubDirPrefix, prokkaAnnotate.out.faa)
    }
    makePlasmidSummary(
        shortReadsCoverage.out.stats.map { file(it.parent) }, // HACK
        longReadsCoverage.out.stats.map { file(it.parent) }, // HACK
        mapToDirectory(quastEvaluate.out.allFiles),
        prokkaAnnotate.out.txt,
        assembly
    )
}
