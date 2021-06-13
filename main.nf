nextflow.enable.dsl=2

include { assembleGenome } from './modules/assembly.nf'
include { evaluateChromosome } from './modules/evaluation.nf'
include { validateParams } from './assemble.nf'
include { getCondaEnv } from './commons.nf'

// parameter validation
validateParams()

// processes and workflows

process findReadyToEvaluate {
    conda "${getCondaEnv(params.mainEnv)}"

    input:
    path assemblyFa

    output:
    path 'assembly_for_eval.fa', optional: true, emit: ready
    path 'assembly_not_for_eval.fa', optional: true, emit: notReady

    script:
    """
    numberOfContigs=`number_of_contigs.py $assemblyFa`
    if [[ ! \$numberOfContigs =~ ^[0-9]+\$ ]]; then
        echo "Error in number_of_contigs.py:"
        echo \$numberOfContigs
        exit 1
    fi
    if [[ \$numberOfContigs == 1 ]]; then
        ln -s $assemblyFa assembly_for_eval.fa
    else
        ln -s $assemblyFa assembly_not_for_eval.fa
    fi
    """
}

process handleNotReadyAssembly {
    input:
    path notReadyAssembly

    exec:
    log.info """
Assembly contains multiple contigs. Assembly evaluation will not be carried out.

Decide how you want to split the contigs and perform evaluation manually using evaluateChromosome.nf and evaluatePlasmid.nf.
    """
}

workflow {
    rawIllumina1Fq = file(params.illumina1)
    rawIllumina2Fq = file(params.illumina2)
    rawLongReadsFq = file(params.pacbio == null ? params.nanopore : params.pacbio)
    longReadType = params.pacbio == null ? LongRead.NANOPORE : LongRead.PACBIO

    assembleGenome(rawIllumina1Fq, rawIllumina2Fq, rawLongReadsFq, longReadType)
    findReadyToEvaluate(assembleGenome.out.assembly)
    evaluateChromosome('evaluation',
                       findReadyToEvaluate.out.ready,
                       assembleGenome.out.cleanedShortReads1,
                       assembleGenome.out.cleanedShortReads2,
                       assembleGenome.out.cleanedLongReads)
    handleNotReadyAssembly(findReadyToEvaluate.out.notReady)
}
