nextflow.enable.dsl=2

import LongRead
include { assembleGenome } from './modules/assembly.nf'
include { evaluateChromosome; evaluatePlasmid } from './modules/evaluation.nf'
include { testSamtools; testBwa; testBbduk; testFiltlong; testFlye;
         testCirclator; testRacon; testCanu; testPilon; testMinimap2;
         testPython_assemblyEnv; testPython_circlatorEnv } from './modules/dependencyChecks.nf'

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

workflow checkDependencies {
    main:
    def doneChannel
    if (params.skipDepChecks) {
        doneChannel = Channel.of(true)
    } else {
        testSamtools()
        testBwa()
        testBbduk()
        testFiltlong()
        testFlye()
        testCirclator()
        testRacon()
        testCanu()
        testPilon()
        testMinimap2()
        testPython_assemblyEnv()
        testPython_circlatorEnv()

        doneChannel = Channel.of(1)
            .mix(testSamtools.out[0], testBwa.out[0],
                 testBbduk.out[0], testFiltlong.out[0],
                 testFlye.out[0], testCirclator.out[0],
                 testRacon.out[0], testCanu.out[0],
                 testPilon.out[0], testMinimap2.out[0],
                 testPython_assemblyEnv.out[0], testPython_circlatorEnv.out[0])
            .toList()
            .map({ true })
    }

    emit:
    done = doneChannel
}

workflow assemble {
    take:
    depChecksDone
    
    main:
    // ensure assembly only starts when the dependency checks are finished
    rawIllumina1Fq = depChecksDone.map({ params.illumina1 }) 
    rawIllumina2Fq = depChecksDone.map({ params.illumina2 })
    rawLongReadsFq = depChecksDone.map({ params.pacbio == null ? params.nanopore : params.pacbio })
    longReadType = params.pacbio == null ? LongRead.NANOPORE : LongRead.PACBIO

    assembleGenome(rawIllumina1Fq, rawIllumina2Fq, rawLongReadsFq, longReadType)
}

workflow {
    checkDependencies()
    assemble(checkDependencies.out)
}
