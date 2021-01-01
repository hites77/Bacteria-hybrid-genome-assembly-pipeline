nextflow.enable.dsl=2

// shouldn't need to override

// output directory relative to outdir
params.o = {}
params.o.cleanShortReads = 'reads/short_cleaned/'
params.o.cleanLongReads = 'reads/long_cleaned/'
params.o.flyeAssembly = 'assembly/flye/'
params.o.raconPolish = 'assembly/racon/'
params.o.canuCorrect = 'assembly/canu/'
params.o.circularise = 'assembly/circlator/'
params.o.pilonPolish = 'assembly/pilon/'

/**
 * Returns a closure to be used with publishDir's saveAs parameter which ensures
 * .command.sh, .command.log and .command.sh are be published to params.oudir + params.o_pubdir.
 *
 * @param o_pubdir: params.o.processName eg. process.o.cleanShortReads
 *
 */
def makeNextflowLogClosure(o_pubdir) {
    return { // it = file name
        if (it == '.exitcode' || it == '.command.log' || it == '.command.sh' ) {
            return params.outdir + o_pubdir + 'nextflow' + it
        } else {
            return it
        }
    }
}

process cleanShortReads {
    publishDir params.outdir, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.cleanShortReads)
    conda params.condaEnvsDir + 'urops-assembly'

    input:
    path illumina1Fq
    path illumina2Fq

    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path params.o.cleanShortReads + 'illumina1.fq', emit: fq1
    path params.o.cleanShortReads + 'illumina2.fq', emit: fq2

    script:
    """
    mkdir -p ${params.o.cleanShortReads}
    bbduk.sh -Xmx1g \
            in1=$illumina1Fq in2=$illumina2Fq \
            out1=${params.o.cleanShortReads}/illumina1.fq out2=${params.o.cleanShortReads}/illumina2.fq \
            qtrim=rl trimq=34 minlength=40
    """
}

process cleanLongReads {
    publishDir params.outdir, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.cleanLongReads)
    conda params.condaEnvsDir + 'urops-assembly'

    input:
    path pacbioFq
    path illumina1Fq
    path illumina2Fq

    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path params.o.cleanLongReads + 'pacbio.fq', emit: fq

    script:
    """
    mkdir -p ${params.o.cleanLongReads}
    filtlong -1 $illumina1Fq -2 $illumina2Fq \
        --min_length 1000 --keep_percent 90 --trim --split 500 --mean_q_weight 10 \
        $pacbioFq > ${params.o.cleanLongReads}/pacbio.fq
    """
}

process flyeAssembly {
    publishDir params.outdir, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.flyeAssembly)
    conda params.condaEnvsDir + 'urops-assembly'

    input:
    path pacbioFq

    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path params.o.flyeAssembly + 'assembly.fasta', emit: assemblyFa
    path params.o.flyeAssembly + '*', emit: allFiles

    script:
    """
    mkdir -p ${params.o.flyeAssembly}
    flye --plasmids --threads $params.threads --pacbio-raw $pacbioFq -o ${params.o.flyeAssembly}
    """
}

process raconPolish {
    publishDir params.outdir, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.raconPolish)
    conda params.condaEnvsDir + 'urops-assembly'

    input:
    path assemblyFa
    path pacbioFq

    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path params.o.raconPolish + 'assembly.fasta', emit: assemblyFa

    script:
    """
    mkdir -p ${params.o.raconPolish}
    minimap2 -t ${params.threads} -ax map-pb $assemblyFa $pacbioFq > assembly_pacbio_alignment.sam
    racon -m 8 -x -6 -g -8 -w 500 -t ${params.threads} \
        $pacbioFq assembly_pacbio_alignment.sam $assemblyFa \
        > ${params.o.raconPolish}/assembly.fasta
    """
}

process canuCorrect {
    publishDir params.outdir, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.canuCorrect)
    conda params.condaEnvsDir + 'urops-assembly'
    
    input:
    path pacbioFq
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path params.o.canuCorrect + 'canu.correctedReads.fasta.gz', emit: pacbioFa
    path params.o.canuCorrect + 'canu*', emit: allFiles

    script:
    """
    canu -correct -p canu -d ${params.o.canuCorrect} genomeSize=5m -pacbio $pacbioFq # useGrid=false
    """
}

process circularise {
    publishDir params.outdir, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.circularise)
    conda params.condaEnvsDir + 'urops-circlator'
    
    input:
    path assemblyFa
    path pacbioFa
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path params.o.circularise + '06.fixstart.fasta', emit: assemblyFa
    path params.o.circularise + '*', emit: allFiles

    script:
    """
    circlator all $assemblyFa $pacbioFa ${params.o.circularise}
    """
}

process pilonPolish {
    publishDir params.outdir + params.o.pilonPolish, mode: 'copy', saveAs: makeNextflowLogClosure(params.o.pilonPolish)
    conda params.condaEnvsDir + 'urops-assembly'
    
    input:
    path assemblyFa
    path illumina1Fq
    path illumina2Fq
    
    output:
    path '.command.sh'
    path '.command.log'
    path '.exitcode'
    path 'pilon*'
    path 'incomplete' optional true

    script:
    """
#!/usr/bin/env python3
import subprocess
from pathlib import Path

def run(assemblyFasta, illumina1Fastq, illumina2Fastq, pilonJar, threads, maxIters):
    bamFile = preprocess(assemblyFasta, illumina1Fastq, illumina2Fastq, threads)
    run_pilon(assemblyFasta, bamFile, 1, pilonJar)
    for i in range(2, maxIters + 1):
        changes_file = Path(f"pilon{i-1}.changes")
        assemblyFasta = f"pilon{i-1}.fasta"
        if not changes_file.exists() or changes_file.read_text().strip() == '':
            break
        bamFile = preprocess(assemblyFasta, illumina1Fastq, illumina2Fastq, threads)
        run_pilon(assemblyFasta, bamFile, i, pilonJar)
    clean()

def preprocess(assemblyFasta, illumina1Fastq, illumina2Fastq, threads):
    bamFile = assemblyFasta[:-len('fasta')] + 'bam'
    subprocess.run(f"bwa index {assemblyFasta}", shell=True, check=True)
    subprocess.run(f"bwa mem -t {threads} {assemblyFasta} {illumina1Fastq} {illumina2Fastq} | samtools sort --threads {threads} > {bamFile}", shell=True, check=True)
    subprocess.run(f"samtools index {bamFile}", shell=True, check=True)
    return bamFile

def run_pilon(assemblyFasta, bamFile, index, pilonJar):
    subprocess.run(f"java -Xmx13G -jar \$PILONJAR --genome {assemblyFasta} --frags {bamFile} --changes --output pilon{index}", shell=True, check=True)

def clean():
    subprocess.run('rm *.bai', shell=True)
    subprocess.run('rm *.bwt', shell=True)
    subprocess.run('rm *.amb', shell=True)
    subprocess.run('rm *.ann', shell=True)
    subprocess.run('rm *.pac', shell=True)
    subprocess.run('rm *.sa', shell=True)

run('$assemblyFasta', '$illumina1Fq', '$illumina2Fq', $params.threads, 6)
    """
}

workflow {
    rawIllumina1Fq = '/home/chloe/Documents/NUS/UROPS/server-data/S8E_3_1/reads/raw/illumina1.fq.gz'
    rawIllumina2Fq = '/home/chloe/Documents/NUS/UROPS/server-data/S8E_3_1/reads/raw/illumina2.fq.gz'
    rawPacbioFq = '/home/chloe/Documents/NUS/UROPS/server-data/S8E_3_1/reads/raw/pacbio.fq.gz'
    
    cleanShortReads(rawIllumina1Fq, rawIllumina2Fq)
    cleanLongReads(rawPacbioFq, cleanShortReads.out.fq1, cleanShortReads.out.fq2)
    flyeAssembly(cleanLongReads.out.fq)
    raconPolish(flyeAssembly.out.assemblyFa, cleanLongReads.out.fq)
    canuCorrect(cleanLongReads.out.fq)
    circularise(raconPolish.out.assemblyFa, canuCorrect.out.pacbioFa)
    pilonPolish(circularise.out.assemblyFa, cleanShortReads.out.fq1, cleanShortReads.out.fq2)
}
