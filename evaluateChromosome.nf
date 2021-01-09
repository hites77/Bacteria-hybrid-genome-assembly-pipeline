nextflow.enable.dsl=2

include { evaluateChromosome } from './modules/evaluation.nf'

workflow {
    evaluateChromosome(
        file(params.assembly),
        file(params.illumina1),
        file(params.illumina2),
        file(params.pacbio),
        file(params.circularity)
    )
}
