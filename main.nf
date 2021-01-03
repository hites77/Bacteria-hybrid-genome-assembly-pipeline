nextflow.enable.dsl=2

/// PARAMS START HERE ///

/// Parameters for specific programs/scripts

// maximum number of iterations to run pilon for
params.pilonMaxIters = 6

// args for bin/bbduk_keep_percent.py
params.bbduk_keep_percent = 80
params.bbduk_start_trimq = 40


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
    run_racon.py --in_assembly $assemblyFa --in_pacbio $pacbioFq --out_prefix final_racon --threads $params.threads --maxiters 4 --args "-m 8 -x -6 -g -8 -w 500"
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

process shortReadsCoverage {
    publishDir params.outdir + params.o.shortReadsCoverage, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.shortReadsCoverage)
    conda params.condaEnvsDir + 'urops-assembly'
    
    input:
    path assemblyFa
    path illumina1Fq
    path illumina2Fq
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path 'stats.txt'
    path 'histogram.txt'
    path 'bbmap_stderr.txt'

    script:
    """
    bbmap.sh in1=$illumina1Fq in2=$illumina2Fq ref=$assemblyFa nodisk \
                covstats=stats.txt covhist=histogram.txt \
            2> bbmap_stderr.txt
    """
}

process longReadsCoverage {
    publishDir params.outdir + params.o.longReadsCoverage, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.longReadsCoverage)
    conda params.condaEnvsDir + 'urops-assembly'
    
    input:
    path assemblyFa
    path pacbioFq
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path 'stats.txt'
    path 'histogram.txt'
    path 'pileup_stderr.txt'
    

    script:
    """
    minimap2 -a $assemblyFa $pacbioFq > mapping.sam
    pileup.sh in=mapping.sam out=stats.txt hist=histogram.txt 2> pileup_stderr.txt
    """
}

process prokkaAnnotate {
    publishDir params.outdir, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.prokkaAnnotate)
    conda params.condaEnvsDir + 'urops-assembly'

    input:
    path assemblyFa
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path params.o.prokkaAnnotate + 'prokka.gff', emit: gff
    path params.o.prokkaAnnotate + '*', emit: allFiles

    script:
    """
    prokka --cpus 0 --outdir ${params.o.prokkaAnnotate} --prefix prokka --addgenes --addmrna --compliant --rfam $assemblyFa
    """
}

process quastEvaluate {
    publishDir params.outdir, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.quastEvaluate)
    conda params.condaEnvsDir + 'urops-assembly'
    
    input:
    path assemblyFa
    path prokkaGff
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path params.o.quastEvaluate + '*'
    
    script:
    """
    quast $assemblyFa --circos -g $prokkaGff -t ${params.threads} --gene-finding --fragmented --conserved-genes-finding --rna-finding -o ${params.o.quastEvaluate}
    """
}

process checkmEvaluate {
    publishDir params.outdir, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.checkmEvaluate)
    conda params.condaEnvsDir + 'urops-checkm'
    
    input:
    path assemblyFa
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path params.o.checkmEvaluate + '*'

    script:
    """
    mkdir input
    mv $assemblyFa input/assembly.fna
    checkm lineage_wf input ${params.o.checkmEvaluate}
    """
}

process makeSummary {
    publishDir params.outdir + params.o.makeSummary, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.makeSummary)
    conda params.condaEnvsDir + 'urops-assembly'

    // just use process.out for these
    input:
    val shortReadsDone
    val longReadsDone
    val quastDone
    val prokkaDone
    val circlatorDone
    val checkmDone

    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path 'summary.json'

    script:
    """
    assembly_summary.py summary.json \
        ${params.outdir + params.o.shortReadsCoverage} \
        ${params.outdir + params.o.longReadsCoverage} \
        ${params.outdir + params.o.quastEvaluate} \
        ${params.outdir + params.o.prokkaAnnotate} \
        ${params.outdir + params.o.circularise} \
        ${params.outdir + params.o.checkmEvaluate}
    """
}

workflow full {
    rawIllumina1Fq = params.rawIllumina1
    rawIllumina2Fq = params.rawIllumina2
    rawPacbioFq = params.rawPacbio
    
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

    finalAssembly = pilonPolish.out.assemblyFa

    shortReadsCoverage(finalAssembly, cleanedShort1, cleanedShort2)
    longReadsCoverage(finalAssembly, cleanedLong)
    prokkaAnnotate(finalAssembly)
    quastEvaluate(finalAssembly, prokkaAnnotate.out.gff)
    checkmEvaluate(finalAssembly)

    makeSummary(shortReadsCoverage.out[0], longReadsCoverage.out[0], quastEvaluate.out[0], prokkaAnnotate.out[0], circularise.out[0], checkmEvaluate.out[0])
}

workflow {
    full()
}

// workflow quickTest {
//     // without canu
    
//     rawIllumina1Fq = '/home/chloe/Documents/NUS/UROPS/server-data/S8E_3_1/reads/raw/illumina1.fq.gz'
//     rawIllumina2Fq = '/home/chloe/Documents/NUS/UROPS/server-data/S8E_3_1/reads/raw/illumina2.fq.gz'
//     rawPacbioFq = '/home/chloe/Documents/NUS/UROPS/server-data/S8E_3_1/reads/raw/pacbio.fq.gz'
    
//     cleanShortReads(rawIllumina1Fq, rawIllumina2Fq)

//     cleanedShort1 = cleanShortReads.out.fq1
//     cleanedShort2 = cleanShortReads.out.fq2

//     cleanLongReads(rawPacbioFq, cleanedShort1, cleanedShort2)

//     cleanedLong = cleanLongReads.out.fq

//     flyeAssembly(cleanedLong)
//     raconPolish(flyeAssembly.out.assemblyFa, cleanedLong)
//     circularise(raconPolish.out.assemblyFa, cleanedLong)
//     pilonPolish(circularise.out.assemblyFa, cleanedShort1, cleanedShort2)

//     finalAssembly = pilonPolish.out.assemblyFa

//     shortReadsCoverage(finalAssembly, cleanedShort1, cleanedShort2)
//     longReadsCoverage(finalAssembly, cleanedLong)
//     prokkaAnnotate(finalAssembly)
//     quastEvaluate(finalAssembly, prokkaAnnotate.out.gff)
//     checkmEvaluate(finalAssembly)

//     makeSummary(shortReadsCoverage.out[0], longReadsCoverage.out[0], quastEvaluate.out[0], prokkaAnnotate.out[0], circularise.out[0], checkmEvaluate.out[0])
// }
