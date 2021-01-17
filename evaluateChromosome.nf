nextflow.enable.dsl=2

include { evaluateChromosome } from './modules/evaluation.nf'
include { testPileup; testBbmap; testProkka; testQuast;
         testMinimap2; testCheckm; testPython_assemblyEnv;
         testPython_checkmEnv } from './modules/dependencyChecks.nf'

if (params.assembly == null || params.illumina1 == null || params.illumina2 == null
    || params.pacbio == null || params.outdir == null) {
    log.error "--assembly, --illumina1, --illumina2, --pacbio, and --outdir are required parameters."
    exit 1
}

workflow checkDependencies {
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
        testCheckm()
        testPython_assemblyEnv()
        testPython_checkmEnv()
        
        doneChannel = Channel.of(1)
            .mix(testPileup.out[0], testBbmap.out[0],
                 testProkka.out[0], testQuast.out[0],
                 testMinimap2.out[0], testCheckm.out[0],
                 testPython_assemblyEnv.out[0], testPython_checkmEnv.out[0])
            .toList()
            .map({ true })
    }

    emit:
    done = doneChannel
}

workflow {
    checkDependencies()

    // ensure evaluation only starts when the dependency checks are finished
    chromosome = checkDependencies.out.map({ file(params.assembly) })
    illumina1 = checkDependencies.out.map({ file(params.illumina1) })
    illumina2 = checkDependencies.out.map({ file(params.illumina2) })
    pacbio = checkDependencies.out.map({ file(params.pacbio) })

    evaluateChromosome('', chromosome, illumina1, illumina2, pacbio)
}
