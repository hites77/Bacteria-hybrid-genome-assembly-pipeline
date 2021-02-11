nextflow.enable.dsl=2

testIllumina1 = file('test-data/illumina1.fq')
testIllumina2 = file('test-data/illumina2.fq')
testPacbio = file('test-data/pacbio.fq') 
testAssembly = file('test-data/assembly.fa')
testAsmPacbioSam = file('test-data/assembly-pacbio.sam')
testAsmIlluminaBam = file('test-data/assembly-illumina.bam')
testGff = file('test-data/assembly-prokka.gff')

process testSamtools {
    errorStrategy 'ignore'
    conda params.condaEnvsDir + '/urops-assembly'

    output:
    path "input.bam.bai"

    script:
    """
    cp $testAsmIlluminaBam input.bam
    samtools index input.bam
    """
}

process testBwa {
    errorStrategy 'ignore'
    conda params.condaEnvsDir + '/urops-assembly'
    
    output:
    path "input.fa.amb"
    path "input.fa.ann"
    path "input.fa.bwt"
    path "input.fa.pac"
    path "input.fa.sa"

    script:
    """
    cp $testAssembly input.fa
    bwa index input.fa
    """
}

process testBbduk {
    errorStrategy 'ignore'
    conda params.condaEnvsDir + '/urops-assembly'

    output:
    path 'cleaned_illumina1.fq'
    path 'cleaned_illumina2.fq'

    script:
    """
    bbduk.sh -Xmx1g in1=$testIllumina1 in2=$testIllumina2 out1=cleaned_illumina1.fq out2=cleaned_illumina2.fq trimq=35
    """
}

process testFiltlong {
    errorStrategy 'ignore'
    conda params.condaEnvsDir + '/urops-assembly'

    output:
    path 'cleaned_pacbio.fq'

    script:
    """
    filtlong -1 $testIllumina1 -2 $testIllumina2 $testPacbio --keep_percent 90 --trim > cleaned_pacbio.fq
    """
}

process testFlye {
    errorStrategy 'ignore'
    conda params.condaEnvsDir + '/urops-assembly'

    output:
    path 'flye/assembly.fasta'

    script:
    testPacbio2 = file('test-data/pacbio2.fa')
    """
    flye --plasmids --pacbio-raw $testPacbio2 -o flye -m 1000
    """
}

process testCirclator {
    errorStrategy 'ignore'
    conda params.condaEnvsDir + '/urops-circlator'

    output:
    path 'circlator-test/06.fixstart.fasta'

    script:
    testPacbio = file("test-data/pacbio.fq") 
    """
    circlator all $testAssembly $testPacbio circlator-test
    """

}

process testRacon {
    errorStrategy 'ignore'
    conda params.condaEnvsDir + '/urops-assembly'

    output:
    path 'test_racon_assembly.fa'

    script:
    """
    racon $testPacbio $testAsmPacbioSam $testAssembly > test_racon_assembly.fa
    """
}

process testCanu {
    errorStrategy 'ignore'
    conda params.condaEnvsDir + '/urops-assembly'
    
    output:
    path 'canu.correctedReads.fasta.gz'

    script:
    """
    canu -correct -p canu genomeSize=70k -pacbio $testPacbio useGrid=false
    """
}

process testPilon {
    errorStrategy 'ignore'
    conda params.condaEnvsDir + '/urops-assembly'

    output:
    path 'pilon.fasta'
    path 'pilon.changes'

    script:
    """
    java -Xmx2G -jar \$PILONJAR --genome $testAssembly --frags $testAsmIlluminaBam --changes --output pilon
    """
}
process testPlaton {
    errorStrategy 'ignore'
    conda params.condaEnvsDir + '/urops-assembly'

    output:
    path 'assembly.chromosome.fasta'
    path 'assembly.tsv'
    path 'assembly.json'

    script:
    """
    platon -p assembly $testAssembly -d \$PLATONDB    
    """
}

process testPileup {
    errorStrategy 'ignore'
    conda params.condaEnvsDir + '/urops-assembly'

    output:
    path 'stats.txt'
    path 'histogram.txt'
    
    script:
    """
    pileup.sh in=$testAsmPacbioSam out=stats.txt hist=histogram.txt
    """
}

process testBbmap {
    errorStrategy 'ignore'
    conda params.condaEnvsDir + '/urops-assembly'

    output:
    path 'stats.txt'
    path 'histogram.txt'
    
    script:
    """
    bbmap.sh in1=$testIllumina1 in2=$testIllumina2 ref=$testAssembly nodisk covstats=stats.txt covhist=histogram.txt
    """
}

process testProkka {
    errorStrategy 'ignore'
    conda params.condaEnvsDir + '/urops-assembly'

    output:
    path 'prokka-test/prokka.gff' // will be others, but just check 1 i guess

    script:
    """
    prokka --outdir prokka-test --prefix prokka --addgenes --addmrna --compliant --rfam $testAssembly
    """
}

process testQuast {
    errorStrategy 'ignore'
    conda params.condaEnvsDir + '/urops-assembly'

    output:
    path 'quast/transposed_report.tsv'
    
    script:
    """
    quast $testAssembly --circos -g $testGff --gene-finding --fragmented --conserved-genes-finding --rna-finding -o quast
    """

}

process testMinimap2 {
    errorStrategy 'ignore'
    conda params.condaEnvsDir + '/urops-assembly'

    output:
    path 'minimap-test.sam'

    script:
    """
    minimap2 -a $testAssembly $testAssembly > minimap-test.sam
    """
}

process testCheckm {
    errorStrategy 'ignore'
    conda params.condaEnvsDir + '/urops-checkm'

    output:
    path 'checkm-results'

    // there will be an error while running checkm if the data folder is empty
    script:
    """
    mkdir input
    cp  $testAssembly input/assembly.fna
    checkm lineage_wf input checkm-results
    """
}

process testPython_assemblyEnv {
    errorStrategy 'ignore'
    conda params.condaEnvsDir + '/urops-assembly'

    script:
    """
    test_python.py
    """
}

process testPython_circlatorEnv {
    errorStrategy 'ignore'
    conda params.condaEnvsDir + '/urops-circlator'

    script:
    """
    test_python.py
    """
}

process testPython_checkmEnv {
    errorStrategy 'ignore'
    conda params.condaEnvsDir + '/urops-checkm'

    script:
    """
    test_python.py
    """
}

process failTest {
    errorStrategy 'ignore'
    conda params.condaEnvsDir + '/urops-checkm'

    script:
    """
    run_nonsesnse
    """
}

// TODO is there a way to run everything to completion before erroring?
workflow checkAllDependencies {
    main:
    // urops-assembly environment
    testBbduk()
    testFiltlong()
    testFlye()
    testBwa()
    testPilon()
    testPlaton()
    testRacon()
    testCanu()
    testSamtools()
    testMinimap2()
    testBbmap()
    testPileup()
    testProkka()
    testQuast()
    testPython_assemblyEnv()

    // urops-circlator environment
    testCirclator()
    testPython_circlatorEnv()

    // urops-checkm environment
    testCheckm()
    testPython_checkmEnv()

    failTest()
}
