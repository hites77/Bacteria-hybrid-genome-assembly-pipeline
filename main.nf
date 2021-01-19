nextflow.enable.dsl=2

include { checkAllDependencies } from './modules/dependencyChecks.nf'
include { assembleGenome } from './modules/assembly.nf'
include { evaluateChromosome } from './modules/evaluation.nf'

// parameter validation

if (params.assembly != null) {
    log.info "The --assembly parameter will be ignored."
}

if (params.illumina1 == null || params.illumina2 == null || params.pacbio == null || params.outdir == null) {
    log.error "--illumina1, --illumina 2, --pacbio, and --outdir are required parameters."
    exit 1
}

if (params.forceCirclator && params.noCirclator) {
    log.error "--forceCirclator and --noCirclator are mutually exclusive. Pass only 1 flag, or pass neither of them."
    exit 1
}

params.outdir = params.outdir + "/"

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
    rawPacbioFq = doneDepChecksChannel.map({ params.pacbio })

    assembleGenome(rawIllumina1Fq, rawIllumina2Fq, rawPacbioFq)
    findReadyToEvaluate(assembleGenome.out.assembly)
    evaluateChromosome('evaluation',
                       findReadyToEvaluate.out.ready,
                       assembleGenome.out.cleanedShortReads1,
                       assembleGenome.out.cleanedShortReads2,
                       assembleGenome.out.cleanedLongReads)
    handleNotReadyAssembly(findReadyToEvaluate.out.notReady)
}
