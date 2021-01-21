nextflow.enable.dsl=2

include { runKofamscan, prokkaAnnotate } from './modules/evaluation.nf'

workflow {
    contig = file(params.assembly)

    prokkaAnnotate('', contig)
    runKofamscan('', prokkaAnnotate.out.faa)
}
