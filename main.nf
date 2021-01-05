nextflow.enable.dsl=2

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

// TODO validate params

// output directories for each process relative to params.outdir
params.o = {}
params.o.cleanShortReads = 'reads/short_cleaned/'
params.o.cleanLongReads = 'reads/long_cleaned/'
params.o.flyeAssembly = 'assembly/flye/'
params.o.raconPolish = 'assembly/racon/'
params.o.canuCorrect = 'reads/long_canu/'
params.o.circularise = 'assembly/circlator/'
params.o.pilonPolish = 'assembly/pilon/'
params.o.separateChromosomesAndPlasmids = 'assembly/separate/'
params.o.shortReadsCoverage = 'assembly_eval/short_read_coverage/'
params.o.longReadsCoverage = 'assembly_eval/long_read_coverage/'
params.o.prokkaAnnotate = 'assembly_eval/prokka/'
params.o.quastEvaluate = 'assembly_eval/quast/'
params.o.checkmEvaluate = 'assembly_eval/checkm/'
params.o.makeSummary = 'summary/'

/// PARAMS END HERE ///

// TODO validate: all dirs end with a slash, no spaces
// TODO Extract out env names?

/**
 * Returns a closure to be used with publishDir's saveAs parameter which ensures
 * .command.sh, .command.log and .command.sh are be published to params.oudir + params.o_pubdir.
 *
 * @param o_pubdir: params.o.processName eg. process.o.cleanShortReads
 *
 */
def makeNextflowLogClosure(o_pubdir) {
    return { // it = file name
        if (it == '.exitcode' || it == '.command.log' || it == '.command.sh' ) {
            return params.outdir + o_pubdir + 'nextflow' + it
        } else {
            return it
        }
    }
}

/**
 * Get longest common directory of a list of files.
 */
def getDirectory(fileList) {
    // make paths absolute
    for (int i=0; i < fileList.size(); i++) {
        fileList[i] = fileList[i].toAbsolutePath()
    }

    // try to find longest common directory
    def directory = fileList[0].isDirectory() ? fileList[0] : file(fileList[0].parent)
    boolean continueFlag = false
    while (true) {
        continueFlag = false
        for (int i=0; i < fileList.size(); i++) {
            if (fileList[i] != directory) {
                continueFlag = true
                if (fileList[i].toString().length() >= directory.toString().length()) {
                    fileList[i] = file(fileList[i].parent)
                }

                if (fileList[i].toString().length() < directory.toString().length()) {
                    directory = fileList[i]
                }
            }
        }
        if (!continueFlag) {
            break
        }
    }

    return directory
}

/**
 * Transforms a channel of lists of files to a channel of directories by applying getDirectory to each list.
 */
def mapToDirectory(fileListChan) {
    return fileListChan.map { getDirectory(it) }
}

process cleanShortReads {
    publishDir params.outdir, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.cleanShortReads)
    conda params.condaEnvsDir + 'urops-assembly'

    input:
    path illumina1Fq
    path illumina2Fq

    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path params.o.cleanShortReads + 'illumina1.fq', emit: fq1
    path params.o.cleanShortReads + 'illumina2.fq', emit: fq2
    path params.o.cleanShortReads + 'trimq_used.txt'

    script:
    """
    mkdir -p ${params.o.cleanShortReads}
    bbduk_keep_percent.py \
            --in1 $illumina1Fq --in2=$illumina2Fq \
            --out1 ${params.o.cleanShortReads}/illumina1.fq --out2 ${params.o.cleanShortReads}/illumina2.fq \
            --infodir ${params.o.cleanShortReads} \
            --keep_percent ${params.bbduk_keep_percent} --start_trimq ${params.bbduk_start_trimq} --args qtrim=rl minlength=40
    """
}

process cleanLongReads {
    publishDir params.outdir, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.cleanLongReads)
    conda params.condaEnvsDir + 'urops-assembly'

    input:
    path pacbioFq
    path illumina1Fq
    path illumina2Fq

    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path params.o.cleanLongReads + 'pacbio.fq', emit: fq

    script:
    """
    mkdir -p ${params.o.cleanLongReads}
    filtlong -1 $illumina1Fq -2 $illumina2Fq \
        --min_length 1000 --keep_percent 90 --trim --split 500 --mean_q_weight 10 \
        $pacbioFq > ${params.o.cleanLongReads}/pacbio.fq
    """
}

process flyeAssembly {
    publishDir params.outdir, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.flyeAssembly)
    conda params.condaEnvsDir + 'urops-assembly'

    input:
    path pacbioFq

    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path params.o.flyeAssembly + 'assembly.fasta', emit: assemblyFa
    path params.o.flyeAssembly + '*', emit: allFiles

    script:
    """
    mkdir -p ${params.o.flyeAssembly} # flye can only create 1 dir
    flye --plasmids --threads $params.threads --pacbio-raw $pacbioFq -o ${params.o.flyeAssembly}
    """
}

process raconPolish {
    publishDir params.outdir + params.o.raconPolish, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.raconPolish)
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
    publishDir params.outdir, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.canuCorrect)
    conda params.condaEnvsDir + 'urops-assembly'
    
    input:
    path pacbioFq
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path params.o.canuCorrect + 'canu.correctedReads.fasta.gz', emit: pacbioFa
    path params.o.canuCorrect + 'canu*', emit: allFiles

    script:
    """
    canu -correct -p canu -d ${params.o.canuCorrect} genomeSize=5m -pacbio $pacbioFq useGrid=false
    """
}

process circularise {
    publishDir params.outdir, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.circularise)
    conda params.condaEnvsDir + 'urops-circlator'
    
    input:
    path assemblyFa
    path pacbioFa
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path params.o.circularise + '06.fixstart.fasta', emit: assemblyFa
    path params.o.circularise + '*', emit: allFiles

    script:
    """
    # circlator can't handle nested directories
    circlator all $assemblyFa $pacbioFa circlator-temp
    mkdir -p ${params.o.circularise}
    mv circlator-temp/* ${params.o.circularise}
    """
}

process pilonPolish {
    publishDir params.outdir + params.o.pilonPolish, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.pilonPolish)
    conda params.condaEnvsDir + 'urops-assembly'
    
    input:
    path assemblyFa
    path illumina1Fq
    path illumina2Fq
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path 'pilon*'
    path 'final_assembly.fasta', emit: assemblyFa
    path 'incomplete' optional true

    script:
    """
    run_pilon.py --assembly $assemblyFa --reads1 $illumina1Fq --reads2 $illumina2Fq --maxiters $params.pilonMaxIters --threads $params.threads
    """
}

process separateChromosomesAndPlasmids {
    publishDir params.outdir + params.o.separateChromosomesAndPlasmids, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.separateChromosomesAndPlasmids)
    conda params.condaEnvsDir + 'urops-assembly'
    
    input:
    path assemblyFa
    path flyeDir
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path 'assembly.chromosome.fasta', emit: chromosomeFa
    path 'assembly.plasmid.fasta' optional true, emit: plasmidFa
    path 'assembly.tsv' optional true, emit: platonTsv
    path 'assembly.json' optional true, emit: platonJson

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

process shortReadsCoverage {
    publishDir params.outdir + params.o.shortReadsCoverage, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.shortReadsCoverage), enabled: params.enablePublish
    conda params.condaEnvsDir + 'urops-assembly'
    
    input:
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
    publishDir params.outdir + params.o.longReadsCoverage, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.longReadsCoverage), enabled: params.enablePublish
    conda params.condaEnvsDir + 'urops-assembly'
    
    input:
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
    publishDir params.outdir, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.prokkaAnnotate), enabled: params.enablePublish
    conda params.condaEnvsDir + 'urops-assembly'

    input:
    path assemblyFa
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path params.o.prokkaAnnotate + 'prokka.gff', emit: gff
    path params.o.prokkaAnnotate + 'prokka.txt', emit: txt
    path params.o.prokkaAnnotate + '*', emit: allFiles

    script:
    """
    prokka --cpus 0 --outdir ${params.o.prokkaAnnotate} --prefix prokka --addgenes --addmrna --compliant --rfam $assemblyFa
    """
}

process quastEvaluate {
    publishDir params.outdir, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.quastEvaluate), enabled: params.enablePublish
    conda params.condaEnvsDir + 'urops-assembly'
    
    input:
    path assemblyFa
    path prokkaGff
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path params.o.quastEvaluate + '*', emit: allFiles
    
    script:
    """
    quast $assemblyFa --circos -g $prokkaGff -t ${params.threads} --gene-finding --fragmented --conserved-genes-finding --rna-finding -o ${params.o.quastEvaluate}
    """
}

process checkmEvaluate {
    publishDir params.outdir, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.checkmEvaluate), enabled: params.enablePublish
    conda params.condaEnvsDir + 'urops-checkm'
    
    input:
    path assemblyFa
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path params.o.checkmEvaluate + '*', emit: allFiles

    script:
    """
    mkdir input
    mv $assemblyFa input/assembly.fna
    checkm lineage_wf input ${params.o.checkmEvaluate}
    """
}

process makeChromosomeSummary {
    // publishDir params.outdir, mode: 'copy', saveAs: makeNextflowLogClosure()
    publishDir params.outdir, mode: 'copy'
    conda params.condaEnvsDir + 'urops-assembly'
    
    input:
    path shortReadsCoverageDir
    path longReadsCoverageDir
    path quastDir
    path prokkaTxt
    path circlatorDir
    path checkmDir
    
    output:
    // path '.command.sh'
    // path '.command.log'
    // path '.exitcode'
    path 'chromosome-summary.json'

    script:
    """
    chromosome_summary.py --short $shortReadsCoverageDir \
        --long $longReadsCoverageDir \
        --quast $quastDir \
        --prokka $prokkaTxt \
        --circlator $circlatorDir \
        --checkm $checkmDir \
        --out chromosome-summary.json
    """
}

process makePlasmidSummary {
    // publishDir params.outdir, mode: 'copy', saveAs: makeNextflowLogClosure()
    publishDir params.outdir, mode: 'copy'
    conda params.condaEnvsDir + 'urops-assembly'
    
    input:
    path shortReadsCoverageDir
    path longReadsCoverageDir
    path quastDir
    path prokkaTxt
    path platonDir
    
    output:
    // path '.command.sh'
    // path '.command.log'
    // path '.exitcode'
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

workflow assembleGenome {
    take:
    rawIllumina1Fq
    rawIllumina2Fq
    rawPacbioFq

    emit:
    chromosomeFa = separateChromosomesAndPlasmids.out.chromosomeFa
    plasmidFa = separateChromosomesAndPlasmids.out.plasmidFa
    cleanedShortReads1 = cleanedShort1
    cleanedShortReads2 = cleanedShort2
    cleanedLongReads = cleanedLong
    circlatorDir = mapToDirectory(circularise.out.allFiles)
    platonTsv = separateChromosomesAndPlasmids.out.platonTsv
    
    main:
    cleanShortReads(rawIllumina1Fq, rawIllumina2Fq)

    cleanedShort1 = cleanShortReads.out.fq1
    cleanedShort2 = cleanShortReads.out.fq2

    cleanLongReads(rawPacbioFq, cleanedShort1, cleanedShort2)

    cleanedLong = cleanLongReads.out.fq

    flyeAssembly(cleanedLong)
    raconPolish(flyeAssembly.out.assemblyFa, cleanedLong)
    canuCorrect(cleanedLong)
    circularise(raconPolish.out.assemblyFa, canuCorrect.out.pacbioFa)
    pilonPolish(circularise.out.assemblyFa, cleanedShort1, cleanedShort2)
    separateChromosomesAndPlasmids(pilonPolish.out.assemblyFa, mapToDirectory(flye.out.allFiles))
}

workflow evaluateChromosome {
    take:
    assembly
    cleanedShortReads1
    cleanedShortReads2
    cleanedLongReads
    circlatorDir

    main:
    shortReadsCoverage(assembly, cleanedShortReads1, cleanedShortReads2)
    longReadsCoverage(assembly, cleanedLongReads)
    prokkaAnnotate(assembly)
    quastEvaluate(assembly, prokkaAnnotate.out.gff)
    checkmEvaluate(assembly)
    makeChromosomeSummary(
        shortReadsCoverage.out.stats.map { file(it.parent) }, // HACK
        longReadsCoverage.out.stats.map { file(it.parent) }, // HACK
        mapToDirectory(quastEvaluate.out.allFiles),
        prokkaAnnotate.out.txt,
        circlatorDir,
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
    shortReadsCoverage(assembly, cleanedShortReads1, cleanedShortReads2)
    longReadsCoverage(assembly, cleanedLongReads)
    prokkaAnnotate(assembly)
    quastEvaluate(assembly, prokkaAnnotate.out.gff)
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

    assembleGenome(rawIllumina1F1, rawIllumina2Fq, rawPacbioFq)
    evaluateChromosome(
        assembleGenome.out.chromosomeFa,
        assembleGenome.out.cleanedShortReads1,
        assembleGenome.out.cleanedShortReads2,
        assembleGenome.out.cleanedLongReads,
        assembleGenome.out.circlatorDir
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
