nextflow.enable.dsl=2

include { checkAllDependencies } from './modules/dependencyChecks.nf'

workflow {
    checkAllDependencies()
}
