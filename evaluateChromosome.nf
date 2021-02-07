nextflow.enable.dsl=2

include { evaluateChromosome } from './modules/evaluation.nf'

if (params.assembly == null || params.illumina1 == null || params.illumina2 == null
    || params.longReads == null || params.outdir == null) {
    log.error "--assembly, --illumina1, --illumina2, --longReads, and --outdir are required parameters."
    exit 1
}

workflow {
    chromosome = file(params.assembly)
    illumina1 = file(params.illumina1)
    illumina2 = file(params.illumina2)
    longReads = file(params.longReads)

    evaluateChromosome('', chromosome, illumina1, illumina2, longReads)
}
