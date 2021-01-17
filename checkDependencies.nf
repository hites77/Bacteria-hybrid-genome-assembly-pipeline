nextflow.enable.dsl=2

include { checkDependencies } from './modules/dependencyChecks.nf'

workflow {
    checkDependencies()
}
