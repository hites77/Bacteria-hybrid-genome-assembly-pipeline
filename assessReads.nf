nextflow.enable.dsl=2

include { assessReads } from './modules/readQuality.nf'

workflow {
    assessReads(params.illumina1, params.illumina2, params.pacbio)
}
