nextflow.enable.dsl=2

include { evaluatePlasmid } from './modules/evaluation.nf'
include { testPileup; testBbmap; testProkka; testQuast;
         testMinimap2; testPython_assemblyEnv; testPython_checkmEnv } from './modules/dependencyChecks.nf'

workflow checkDepsIfNecessary {
    main:
    def doneChannel
    if (params.skipDepChecks) {
        doneChannel = Channel.of(true)
    } else {
        testPileup()
        testBbmap()
        testProkka()
        testQuast()
        testMinimap2()
        testPython_assemblyEnv()
        
        doneChannel = Channel.of(1)
            .mix(testPileup.out[0], testBbmap.out[0],
                 testProkka.out[0], testQuast.out[0],
                 testMinimap2.out[0], testPython_assemblyEnv.out[0])
            .toList()
            .map({ true })
    }

    emit:
    done = doneChannel
}

workflow {
    checkDepsIfNecessary()

    // ensure evaluation only starts when the dependency checks are finished
    plasmid = checkDepsIfNecessary.out.map({ file(params.plasmid) })
    illumina1 = checkDepsIfNecessary.out.map({ file(params.illumina1) })
    illumina2 = checkDepsIfNecessary.out.map({ file(params.illumina2) })
    pacbio = checkDepsIfNecessary.out.map({ file(params.pacbio) })

    evaluatePlasmid('', plasmid, illumina1, illumina2, pacbio)
}
