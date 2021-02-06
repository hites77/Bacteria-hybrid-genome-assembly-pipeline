nextflow.enable.dsl=2

include { checkAllDependencies } from './modules/dependencyChecks.nf'
include { assembleGenome } from './modules/assembly.nf'
include { evaluateChromosome } from './modules/evaluation.nf'
include { validateParams } from './assemble.nf'

// parameter validation
validateParams()

// processes and workflows

process findReadyToEvaluate {
    conda params.condaEnvsDir + 'urops-assembly'

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
    def doneDepChecksChannel
    if (params.skipDepChecks) {
        doneDepChecksChannel = Channel.of(true)
    } else {
        checkAllDependencies()
        doneDepChecksChannel = checkAllDependencies.out
    }


    // ensure assembly only begins after dependency checks are done
    rawIllumina1Fq = doneDepChecksChannel.map({ params.illumina1 })
    rawIllumina2Fq = doneDepChecksChannel.map({ params.illumina2 })
    rawLongReadsFq = depChecksDone.map({ params.pacbio == null ? params.nanopore : params.pacbio })
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
