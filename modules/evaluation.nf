nextflow.enable.dsl=2

include { getDirectory; mapToDirectory } from './commons.nf'

params.outdir = 'genome-assembly/'

// output directories for each process relative to params.outdir
outdirs = {}

outdirs.evaluateChromosome = 'chromosome_eval/'
outdirs.evaluatePlasmid = 'plasmid_eval/'

outdirs.shortReadsCoverage = 'short_read_coverage/'
outdirs.longReadsCoverage = 'long_read_coverage/'
outdirs.prokkaAnnotate = 'prokka/'
outdirs.quastEvaluate = 'quast/'
outdirs.checkmEvaluate = 'checkm/'

process shortReadsCoverage {
    publishDir params.outdir + outdirs.shortReadsCoverage, mode: 'copy', pattern: '{.command.sh,.command.log,.exitcode}', saveAs: { 'nextflow' + it }
    publishDir "${params.outdir}/${pubDirPrefix}/${outdirs.shortReadsCoverage}" , mode: 'copy'
    conda params.condaEnvsDir + 'urops-assembly'
    
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
    publishDir params.outdir + outdirs.longReadsCoverage, mode: 'copy', pattern: '{.command.sh,.command.log,.exitcode}', saveAs: { 'nextflow' + it }
    publishDir "${params.outdir}/${pubDirPrefix}/${outdirs.longReadsCoverage}", mode: 'copy'
    conda params.condaEnvsDir + 'urops-assembly'
    
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
    publishDir params.outdir + outdirs.prokkaAnnotate, mode: 'copy', pattern: '{.command.sh,.command.log,.exitcode}', saveAs: { 'nextflow' + it }
    publishDir "${params.outdir}/${pubDirPrefix}/${outdirs.prokkaAnnotate}", mode: 'copy'
    conda params.condaEnvsDir + 'urops-assembly'

    input:
    val pubDirPrefix
    path assemblyFa
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path outdirs.prokkaAnnotate + 'prokka.gff', emit: gff
    path outdirs.prokkaAnnotate + 'prokka.txt', emit: txt
    path outdirs.prokkaAnnotate + '*', emit: allFiles

    script:
    """
    prokka --cpus 0 --outdir ${outdirs.prokkaAnnotate} --prefix prokka --addgenes --addmrna --compliant --rfam $assemblyFa
    """
}

process quastEvaluate {
    publishDir params.outdir + outdirs.quastEvaluate, mode: 'copy', pattern: '{.command.sh,.command.log,.exitcode}', saveAs: { 'nextflow' + it }
    publishDir "${params.outdir}/${pubDirPrefix}/${outdirs.quastEvaluate}", mode: 'copy'
    conda params.condaEnvsDir + 'urops-assembly'
    
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
    publishDir params.outdir + outdirs.checkmEvaluate, mode: 'copy', pattern: '{.command.sh,.command.log,.exitcode}', saveAs: { 'nextflow' + it }
    publishDir "${params.outdir}/${pubDirPrefix}/${outdirs.checkmEvaluate}", mode: 'copy'
    conda params.condaEnvsDir + 'urops-checkm'
    
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

process makeChromosomeSummary {
    publishDir params.outdir, mode: 'copy'
    conda params.condaEnvsDir + 'urops-assembly'
    
    input:
    path shortReadsCoverageDir
    path longReadsCoverageDir
    path quastDir
    path prokkaTxt
    path circularitySummary
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
        --circularity $circularitySummary \
        --assembly $assemblyFa \
        --checkm $checkmDir \
        --out chromosome-summary.json
    """
}

process makePlasmidSummary {
    publishDir params.outdir, mode: 'copy'
    conda params.condaEnvsDir + 'urops-assembly'
    
    input:
    path shortReadsCoverageDir
    path longReadsCoverageDir
    path quastDir
    path prokkaTxt
    path platonDir
    
    output:
    path 'plasmid-summary.json'

    script:
    """
    plasmid_summary.py --short $shortReadsCoverageDir \
        --long $longReadsCoverageDir \
        --quast $quastDir \
        --prokka $prokkaTxt \
        --platon $platonDir
    """
}

workflow evaluateChromosome {
    take:
    assembly
    cleanedShortReads1
    cleanedShortReads2
    cleanedLongReads
    circularitySummary

    main:
    shortReadsCoverage(outdirs.evaluateChromosome, assembly, cleanedShortReads1, cleanedShortReads2)
    longReadsCoverage(outdirs.evaluateChromosome, assembly, cleanedLongReads)
    prokkaAnnotate(outdirs.evaluateChromosome, assembly)
    quastEvaluate(outdirs.evaluateChromosome, assembly, prokkaAnnotate.out.gff)
    checkmEvaluate(outdirs.evaluateChromosome, assembly)
    makeChromosomeSummary(
        shortReadsCoverage.out.stats.map { file(it.parent) }, // HACK
        longReadsCoverage.out.stats.map { file(it.parent) }, // HACK
        mapToDirectory(quastEvaluate.out.allFiles),
        prokkaAnnotate.out.txt,
        circularitySummary,
        assembly,
        mapToDirectory(checkmEvaluate.out.allFiles)
    )
}

workflow evaluatePlasmid {
    take:
    assembly
    cleanedShortReads1
    cleanedShortReads2
    cleanedLongReads
    platonTsv

    main:
    shortReadsCoverage(outdirs.evaluatePlasmid, assembly, cleanedShortReads1, cleanedShortReads2)
    longReadsCoverage(outdirs.evaluatePlasmid, assembly, cleanedLongReads)
    prokkaAnnotate(outdirs.evaluatePlasmid, assembly)
    quastEvaluate(outdirs.evaluatePlasmid, assembly, prokkaAnnotate.out.gff)
    makePlasmidSummary(
        shortReadsCoverage.out.stats.map { file(it.parent) }, // HACK
        longReadsCoverage.out.stats.map { file(it.parent) }, // HACK
        mapToDirectory(quastEvaluate.out.allFiles),
        prokkaAnnotate.out.txt,
        platonTsv
    )
}
