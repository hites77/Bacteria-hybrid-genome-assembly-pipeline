nextflow.enable.dsl=2

process assessShortReads {
    publishDir "${params.outdir}/fastqc", mode: 'move'
    conda params.condaEnvsDir + "urops-assembly"
    
    input:
    path illumina1Fq
    path illumina2Fq
    
    output:
    path '*_fastqc.html'
    path '*_fastqc.zip'

    script:
    """
    fastqc -t $params.threads $illumina1Fq $illumina2Fq
    """
}

process assessLongReads {
    publishDir "${params.outdir}/nanoplot", mode: 'move'
    conda params.condaEnvsDir + "urops-assembly"
    
    input:
    path pacbioFq
    
    output:
    path '*'

    script:
    """
    NanoPlot -t $params.threads --fastq $pacbioFq -f png --N50 --dpi 300 -o ./
    """
}

workflow assessReads {
    take:
    illumina1Fq
    illumina2Fq
    pacbioFq

    main:
    assessShortReads(illumina1Fq, illumina2Fq)
    assessLongReads(pacbioFq)
}
