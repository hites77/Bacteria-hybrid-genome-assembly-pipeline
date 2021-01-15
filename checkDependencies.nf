nextflow.enable.dsl=2

include { checkDependencies } from './modules/dependency_checks.nf'

workflow {
    checkDependencies()
}
