nextflow.enable.dsl=2

include { assembleGenome } from './modules/assembly.nf'
include { evaluateChromosome; evaluatePlasmid } from './modules/evaluation.nf'
include { checkDependencies } from './modules/dependency_checks.nf'

// TODO validate params
// TODO validate: all dirs end with a slash, no spaces

workflow full {
    take:
    depChecksDone
    
    main:
    // HACK ensure assembly only start when the dependency checks are finished
    rawIllumina1Fq = depChecksDone.map({ params.illumina1 }) 
    rawIllumina2Fq = depChecksDone.map({ params.illumina2 })
    rawPacbioFq = depChecksDone.map({ params.pacbio })

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
    def doneChannel
    if (params.skipDepChecks) {
        doneChannel = Channel.of(true)
    } else {
        checkDependencies()
        doneChannel = checkDependencies.out
    }

    full(doneChannel)
}
