nextflow.enable.dsl=2

import LongRead
include { assembleGenome } from './modules/assembly.nf'
include { evaluateChromosome; evaluatePlasmid } from './modules/evaluation.nf'

/**
* Ensures that:
* - --illumina1, --illumina2, --outdir and exactly one of --nanopore or --pacbio are present.
* - At most one of --forceCirclator and --noCirclator are present
* 
* Gives warning (but still continues to run) if:
* - --assembly is passed
*/
def validateParams() {
    if (params.assembly != null) {
        log.info "The --assembly parameter will be ignored."
    }

    // check required params are present
    if (params.illumina1 == null || params.illumina2 == null || params.outdir == null) {
        log.error "--illumina1, --illumina2, and --outdir are required parameters."
        exit 1
    }

    if (params.nanopore == null && params.pacbio == null) {
        log.error "Long reads are missing. Pass the long reads either using --nanopore or --pacbio."
        exit 1
    }

    if (params.nanopore != null && params.pacbio != null) {
        log.error "Cannot pass both --nanopore and --pacbio. Pass only one of them."
        exit 1
    }

    if (params.forceCirclator && params.noCirclator) {
        log.error "--forceCirclator and --noCirclator are mutually exclusive. Pass only 1 flag, or pass neither of them."
        exit 1
    }
}

validateParams()

workflow {
    rawIllumina1Fq = file(params.illumina1)
    rawIllumina2Fq = file(params.illumina2)
    rawLongReadsFq = file(params.pacbio == null ? params.nanopore : params.pacbio)
    longReadType = params.pacbio == null ? LongRead.NANOPORE : LongRead.PACBIO

    assembleGenome(rawIllumina1Fq, rawIllumina2Fq, rawLongReadsFq, longReadType)
}
