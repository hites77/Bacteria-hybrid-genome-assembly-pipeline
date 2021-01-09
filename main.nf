nextflow.enable.dsl=2

include { assembleGenome } from './modules/assembly.nf'
include { evaluateChromosome; evaluatePlasmid } from './modules/evaluation.nf'


/// Parameters for specific programs/scripts

// TODO validate params


// TODO validate: all dirs end with a slash, no spaces
// TODO Extract out env names?

workflow full {
    rawIllumina1Fq = params.illumina1
    rawIllumina2Fq = params.illumina2
    rawPacbioFq = params.pacbio

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
