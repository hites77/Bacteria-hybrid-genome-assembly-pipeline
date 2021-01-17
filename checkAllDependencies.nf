nextflow.enable.dsl=2

include { checkAllDependencies } from './modules/dependencyChecks.nf'

if (params.assembly != null) {
    log.info "The --assembly parameter will be ignored."
}

if (params.illumina1 != null) {
    log.info "The --illumina1 parameter will be ignored."
}

if (params.illumina2 != null) {
    log.info "The --illumina2 parameter will be ignored."
}

if (params.pacbio != null) {
    log.info "The --pacbio parameter will be ignored."
}

if (params.outdir != null) {
    log.info "The --outdir parameter will be ignored."
}

workflow {
    checkAllDependencies()
}
