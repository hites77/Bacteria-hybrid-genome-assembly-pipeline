nextflow.enable.dsl=2

include { cleanShortReads; cleanLongReads } from './modules/assembly.nf'
include { assessReads } from './modules/readQuality.nf'

workflow {
    cleanShortReads(params.rawIllumina1, params.rawIllumina2)
    cleanLongReads(params.rawPacbio, cleanShortReads.out.fq1, cleanShortReads.out.fq2)
    assessReads(cleanShortReads.out.fq1, cleanShortReads.out.fq2, cleanLongReads.out.fq)
}
