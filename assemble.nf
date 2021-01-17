nextflow.enable.dsl=2

include { assembleGenome } from './modules/assembly.nf'
include { evaluateChromosome; evaluatePlasmid } from './modules/evaluation.nf'
include { testSamtools; testBwa; testBbduk; testFiltlong; testFlye;
         testCirclator; testRacon; testCanu; testPilon; testMinimap2;
         testPython_assemblyEnv; testPython_circlatorEnv } from './modules/dependencyChecks.nf'

// TODO validate params
// TODO validate: all dirs end with a slash, no spaces

workflow checkDepsIfNecessary {
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
    rawPacbioFq = depChecksDone.map({ params.pacbio })

    assembleGenome(rawIllumina1Fq, rawIllumina2Fq, rawPacbioFq)
}

workflow {
    checkDepsIfNecessary()
    assemble(checkDepsIfNecessary.out)
}
