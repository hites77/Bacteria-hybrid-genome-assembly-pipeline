nextflow.enable.dsl=2

include { getDirectory; mapToDirectory } from './modules/commons.nf'

/// PARAMS START HERE ///

/// Parameters for specific programs/scripts

// maximum number of iterations to run pilon for
params.pilonMaxIters = 6
// maximum number of iterations to run racon for
params.raconMaxIters = 4

// args for bin/bbduk_keep_percent.py
params.bbduk_keep_percent = 80
params.bbduk_start_trimq = 40

params.enablePublish = false
params.outdir = 'genome-assembly/'

/// PARAMS END HERE ///

// TODO validate params

// output directories for each process relative to params.outdir
outdirs = {}
outdirs.cleanShortReads = 'reads/short_cleaned/'
outdirs.cleanLongReads = 'reads/long_cleaned/'
outdirs.flyeAssembly = 'assembly/flye/'
outdirs.raconPolish = 'assembly/racon/'
outdirs.canuCorrect = 'reads/long_canu/'
outdirs.circlator = 'assembly/circlator/'
outdirs.pilonPolish = 'assembly/pilon/'
outdirs.separateChromosomesAndPlasmids = 'assembly/separate/'

outdirs.evaluateChromosome = 'chromosome_eval/'
outdirs.evaluatePlasmid = 'plasmid_eval/'

outdirs.shortReadsCoverage = 'short_read_coverage/'
outdirs.longReadsCoverage = 'long_read_coverage/'
outdirs.prokkaAnnotate = 'prokka/'
outdirs.quastEvaluate = 'quast/'
outdirs.checkmEvaluate = 'checkm/'

// TODO validate: all dirs end with a slash, no spaces
// TODO Extract out env names?

process cleanShortReads {
    publishDir params.outdir + outdirs.cleanShortReads, mode: 'copy', pattern: '{.command.sh,.command.log,.exitcode}', saveAs: { 'nextflow' + it }
    publishDir params.outdir, mode: 'copy'
    conda params.condaEnvsDir + 'urops-assembly'

    input:
    path illumina1Fq
    path illumina2Fq

    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path outdirs.cleanShortReads + 'illumina1.fq', emit: fq1
    path outdirs.cleanShortReads + 'illumina2.fq', emit: fq2
    path outdirs.cleanShortReads + 'trimq_used.txt'

    script:
    """
    mkdir -p ${outdirs.cleanShortReads}
    bbduk_keep_percent.py \
            --in1 $illumina1Fq --in2=$illumina2Fq \
            --out1 ${outdirs.cleanShortReads}/illumina1.fq --out2 ${outdirs.cleanShortReads}/illumina2.fq \
            --infodir ${outdirs.cleanShortReads} \
            --keep_percent ${params.bbduk_keep_percent} --start_trimq ${params.bbduk_start_trimq} --args qtrim=rl minlength=40
    """
}

process cleanLongReads {
    publishDir params.outdir + outdirs.cleanLongReads, mode: 'copy', pattern: '{.command.sh,.command.log,.exitcode}', saveAs: { 'nextflow' + it }
    publishDir params.outdir, mode: 'copy'
    conda params.condaEnvsDir + 'urops-assembly'

    input:
    path pacbioFq
    path illumina1Fq
    path illumina2Fq

    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path outdirs.cleanLongReads + 'pacbio.fq', emit: fq

    script:
    """
    mkdir -p ${outdirs.cleanLongReads}
    filtlong -1 $illumina1Fq -2 $illumina2Fq \
        --min_length 1000 --keep_percent 90 --trim --split 500 --mean_q_weight 10 \
        $pacbioFq > ${outdirs.cleanLongReads}/pacbio.fq
    """
}

process flyeAssembly {
    publishDir params.outdir + outdirs.flyeAssembly, mode: 'copy', pattern: '{.command.sh,.command.log,.exitcode}', saveAs: { 'nextflow' + it }
    publishDir params.outdir, mode: 'copy'
    conda params.condaEnvsDir + 'urops-assembly'

    input:
    path pacbioFq

    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path outdirs.flyeAssembly + 'assembly.fasta', emit: assemblyFa
    path outdirs.flyeAssembly + '*', emit: allFiles
    path 'any_circular_contigs.txt', emit: anyCircularFile

    script:
    """
    mkdir -p ${outdirs.flyeAssembly} # flye can only create 1 dir
    flye --plasmids --threads $params.threads --pacbio-raw $pacbioFq -o ${outdirs.flyeAssembly}
    flye_circularity.py ${outdirs.flyeAssembly} > 'any_circular_contigs.txt'
    """
}

process raconPolish {
    publishDir params.outdir + outdirs.raconPolish, mode: 'copy', pattern: '{.command.sh,.command.log,.exitcode}', saveAs: { 'nextflow' + it }
    publishDir params.outdir + outdirs.raconPolish, mode: 'copy'
    conda params.condaEnvsDir + 'urops-assembly'

    input:
    path assemblyFa
    path pacbioFq

    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path 'final_racon_assembly.fa', emit: assemblyFa
    path 'final_racon_log.tsv', emit: logFile

    script:
    """
    run_racon.py --in_assembly $assemblyFa --in_pacbio $pacbioFq --out_prefix final_racon --threads $params.threads --maxiters $params.raconMaxIters --args "-m 8 -x -6 -g -8 -w 500"
    """
}

process canuCorrect {
    publishDir params.outdir + outdirs.canuCorrect, mode: 'copy', pattern: '{.command.sh,.command.log,.exitcode}', saveAs: { 'nextflow' + it }
    publishDir params.outdir, mode: 'copy'
    conda params.condaEnvsDir + 'urops-assembly'
    
    input:
    path pacbioFq
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path outdirs.canuCorrect + 'canu.correctedReads.fasta.gz', emit: pacbioFa
    path outdirs.canuCorrect + 'canu*', emit: allFiles

    script:
    """
    canu -correct -p canu -d ${outdirs.canuCorrect} genomeSize=5m -pacbio $pacbioFq useGrid=false
    """
}

process circlator {
    publishDir params.outdir + outdirs.circlator, mode: 'copy', pattern: '{.command.sh,.command.log,.exitcode}', saveAs: { 'nextflow' + it }
    publishDir params.outdir, mode: 'copy'
    conda params.condaEnvsDir + 'urops-circlator'
    
    input:
    path assemblyFa
    path pacbioFa
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path outdirs.circlator + '06.fixstart.fasta', emit: assemblyFa
    path outdirs.circlator + '*', emit: allFiles
    path 'circularity_summary.json', emit: circularitySummary

    script:
    """
    # circlator can't handle nested directories
    circlator all $assemblyFa $pacbioFa circlator-temp
    circlator_circularity_summary.py circlator-temp > circularity_summary.json

    mkdir -p ${outdirs.circlator}
    mv circlator-temp/* ${outdirs.circlator}
    """
}

process pilonPolish {
    publishDir params.outdir + outdirs.pilonPolish, mode: 'copy', pattern: '{.command.sh,.command.log,.exitcode}', saveAs: { 'nextflow' + it }
    publishDir params.outdir + outdirs.pilonPolish, mode: 'copy'
    conda params.condaEnvsDir + 'urops-assembly'
    
    input:
    path assemblyFa
    path illumina1Fq
    path illumina2Fq
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path 'final_assembly.fasta', emit: assemblyFa
    path 'pilon*.changes'
    path 'pilon_info.tsv'

    script:
    """
    run_pilon.py --assembly $assemblyFa --reads1 $illumina1Fq --reads2 $illumina2Fq --out final_assembly.fasta \
                --maxiters $params.pilonMaxIters --threads $params.threads
    """
}

process separateChromosomesAndPlasmids {
    publishDir params.outdir + outdirs.separateChromosomesAndPlasmids, mode: 'copy', pattern: '{.command.sh,.command.log,.exitcode}', saveAs: { 'nextflow' + it }
    publishDir params.outdir + outdirs.separateChromosomesAndPlasmids, mode: 'copy'
    conda params.condaEnvsDir + 'urops-assembly'
    
    input:
    path assemblyFa
    path flyeDir
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path 'assembly.chromosome.fasta', emit: chromosomeFa
    path 'assembly.plasmid.fasta' ,optional: true, emit: plasmidFa
    path 'assembly.tsv' ,optional: true, emit: platonTsv
    path 'assembly.json' ,optional: true, emit: platonJson

    script:
    """
    if [ -s $flyeDir/22-plasmids/plasmids_raw.fasta ] ; then
        echo plasmids present
        platon -p assembly -t $params.threads $assemblyFa -d \$PLATONDB
    else
        echo plasmids absent
        ln -s $assemblyFa assembly.chromosome.fasta
    fi
    """
}

// evaluation

process shortReadsCoverage {
    publishDir params.outdir + outdirs.shortReadsCoverage, mode: 'copy', pattern: '{.command.sh,.command.log,.exitcode}', saveAs: { 'nextflow' + it }
    // publishDir params.outdir + "${pubDirPrefix}" + outdirs.shortReadsCoverage, mode: 'copy', saveAs: makeNextflowLogClosure("${pubDirPrefix}" + outdirs.shortReadsCoverage), enabled: params.enablePublish
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
    // publishDir params.outdir + "${pubDirPrefix}" + outdirs.longReadsCoverage, mode: 'copy', saveAs: makeNextflowLogClosure("${pubDirPrefix}" + outdirs.longReadsCoverage), enabled: params.enablePublish
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
    // publishDir params.outdir + "${pubDirPrefix}", mode: 'copy', saveAs: makeNextflowLogClosure("${pubDirPrefix}" + outdirs.prokkaAnnotate), enabled: params.enablePublish
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
    // publishDir params.outdir + "${pubDirPrefix}", mode: 'copy', saveAs: makeNextflowLogClosure("${pubDirPrefix}" + outdirs.quastEvaluate), enabled: params.enablePublish
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
    // publishDir params.outdir + "${pubDirPrefix}", mode: 'copy', saveAs: makeNextflowLogClosure("${pubDirPrefix}" + outdirs.checkmEvaluate), enabled: params.enablePublish
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

process splitByCircularity {
    input:
    path assemblyFa
    path anyCircularFile

    output:
    path 'circular-assembly.fa' ,optional: true, emit: circularAssembly
    path 'linear-assembly.fa' ,optional: true, emit: linearAssembly 

    script:
    """
    if [[ `cat $anyCircularFile` == yes ]]; then
        ln -s $assemblyFa circular-assembly.fa
    else
        ln -s $assemblyFa linear-assembly.fa
    fi
    """
}

workflow circulariseAssembly {
    take:
    assembly
    longReads
    
    main:
    canuCorrect(longReads)
    circlator(assembly, canuCorrect.out.pacbioFa)

    emit:
    assembly = circlator.out.assemblyFa
    circularitySummary = circlator.out.circularitySummary
}

process makeLinearAssemblyCircSummary {
    input:
    path assemblyFa

    output:
    path 'circularity_summary.json', emit: circularitySummary

    script:
    """echo \"all linear\" > circularity_summary.json"""
}

workflow assembleGenome {
    take:
    rawIllumina1Fq
    rawIllumina2Fq
    rawPacbioFq

    main:
    cleanShortReads(rawIllumina1Fq, rawIllumina2Fq)

    cleanedShort1 = cleanShortReads.out.fq1
    cleanedShort2 = cleanShortReads.out.fq2

    cleanLongReads(rawPacbioFq, cleanedShort1, cleanedShort2)

    cleanedLong = cleanLongReads.out.fq

    flyeAssembly(cleanedLong)
    raconPolish(flyeAssembly.out.assemblyFa, cleanedLong)
    splitByCircularity(raconPolish.out.assemblyFa, flyeAssembly.out.anyCircularFile)

    circulariseAssembly(splitByCircularity.out.circularAssembly, cleanedLong)
    makeLinearAssemblyCircSummary(splitByCircularity.out.linearAssembly)

    assembly = circulariseAssembly.out.assembly.mix(splitByCircularity.out.linearAssembly)
    circularitySummary = circulariseAssembly.out.circularitySummary.mix(makeLinearAssemblyCircSummary.out.circularitySummary)

    pilonPolish(assembly, cleanedShort1, cleanedShort2)
    separateChromosomesAndPlasmids(pilonPolish.out.assemblyFa, mapToDirectory(flyeAssembly.out.allFiles))

    emit:
    chromosomeFa = separateChromosomesAndPlasmids.out.chromosomeFa
    plasmidFa = separateChromosomesAndPlasmids.out.plasmidFa
    cleanedShortReads1 = cleanedShort1
    cleanedShortReads2 = cleanedShort2
    cleanedLongReads = cleanedLong
    circularitySummary = circularitySummary
    platonTsv = separateChromosomesAndPlasmids.out.platonTsv
    
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

workflow full {
    rawIllumina1Fq = params.rawIllumina1
    rawIllumina2Fq = params.rawIllumina2
    rawPacbioFq = params.rawPacbio

    assembleGenome(rawIllumina1Fq, rawIllumina2Fq, rawPacbioFq)
    evaluateChromosome(
        assembleGenome.out.chromosomeFa,
        assembleGenome.out.cleanedShortReads1,
        assembleGenome.out.cleanedShortReads2,
        assembleGenome.out.cleanedLongReads,
        assembleGenome.out.circularitySummary
    )
    evaluatePlasmid(
        assembleGenome.out.plasmidFa,
        assembleGenome.out.cleanedShortReads1,
        assembleGenome.out.cleanedShortReads2,
        assembleGenome.out.cleanedLongReads,
        assembleGenome.out.platonTsv
    )
}

workflow {
    full()
}
