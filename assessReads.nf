nextflow.enable.dsl=2

include { assessReads } from './modules/readQuality.nf'

workflow {
    assessReads(params.rawIllumina1, params.rawIllumina2, params.rawPacbio)
}
