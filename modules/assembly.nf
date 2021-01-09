nextflow.enable.dsl=2

include { getDirectory; mapToDirectory } from './commons.nf'

params.outdir = 'genome-assembly/'

// maximum number of iterations to run pilon for
params.pilonMaxIters = 6
// maximum number of iterations to run racon for
params.raconMaxIters = 4

// args for bin/bbduk_keep_percent.py
params.bbdukKeepPercent = 80
params.bbdukStartTrimq = 40
params.bbdukMinTrimq = 28
params.bbdukArgs = 'qtrim=rl minlength=40'

outdirs = {}
outdirs.cleanShortReads = 'reads/short_cleaned/'
outdirs.cleanLongReads = 'reads/long_cleaned/'
outdirs.flyeAssembly = 'assembly/flye/'
outdirs.raconPolish = 'assembly/racon/'
outdirs.canuCorrect = 'reads/long_canu/'
outdirs.circlator = 'assembly/circlator/'
outdirs.pilonPolish = 'assembly/pilon/'
outdirs.separateChromosomesAndPlasmids = 'assembly/separate/'

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
            --keep_percent $params.bbdukKeepPercent --start_trimq $params.bbdukStartTrimq \
            --min_trimq $params.bbdukMinTrimq --args $params.bbdukArgs
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
    flye_circularity.py ${outdirs.flyeAssembly} > ${outdirs.flyeAssembly}/any_circular_contigs.txt
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
    mkdir -p ${outdirs.circlator}

    # circlator can't handle nested directories
    circlator all $assemblyFa $pacbioFa circlator-temp
    circlator_circularity_summary.py circlator-temp > ${outdirs.circlator}/circularity_summary.json

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
