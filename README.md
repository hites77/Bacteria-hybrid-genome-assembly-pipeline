# Genome assembly pipeline

This is a hybrid genome assembly pipeline for bacterial genomes written in Nextflow.

## The pipeline

_Note: Name of Nextflow process is written in parenthesis._

### Assembly

**Inputs:**

- Paired end illumina short reads, in 2 separate fastq files (may be zipped).
- Pacbio long reads, as 1 fastq file (may be zipped).

**Steps:** 

1. Clean short reads using Bbduk (`cleanShortReads`)
1. Clean long reads using Filtlong (`cleanLongReads`)
1. Form an initial assembly using the long reads with Flye (`flyeAssembly`)
1. Polish the assembly using the long reads with Racon (`raconPolish`)
1. Circularise the assembly using Circlator (`circularise`)
   - In addition to the assembly, Circlator takes in the long reads as input. We use long reads corrected using Canu as input (`canuCorrect`).
1. Polish the assembly using Pilon (`pilonPolish`)

### Assembly evaluation

TODO

## Setup

1. Create the conda environments from yaml files in the `conda-envs` directory
   - Use only one of `assembly.yml` and `assembly-nscc.yml`. The difference is that `assembly.yml` includes samtools and bwa while `assembly-nscc.yml` does not.
   - For `assembly.yml`/`assembly-nscc.yml`, change the PILONJAR variable according. It is meant to be the path to Pilon's jar file. It is usually at `<path to conda directory>/pkgs/pilon-1.23-2/share/pilon-1.23-2/pilon-1.23.jar`.
1. Install [Nextflow](https://www.nextflow.io). It has been tested using v20.10 so far (latest stable release as of time of writing).

## Running the pipeline

### 1. Set variables

All parameters which can be adjusted can be found at the top of `main.nf` and in `nextflow.config`. Parameters which are likely to be device specific (eg. number of threads) are defined in profiles in `nextflow.config`.

Recommended variables to adjust: 
- Device specific parameters (under a profile in `nextflow.config`)

### 2. Start the pipeline

There are a few ways to run the pipeline:

**As a job on the NSCC Aspire 1 server:**

Submit the following job script using `qsub`:
``` sh
#!/bin/bash
#PBS -q normal
#PBS -P Personal
#PBS -j oe
#PBS -l select=1:ncpus=10:mpiprocs=10
#PBS -l walltime=04:00:00
cd $PBS_O_WORKDIR
module load java
module load openmpi
module load samtools/1.3
module load bwa/0.7.13

# cd to the directory containing this repo 
mpirun --pernode nextflow run main.nf -profile nscc -with-mpi
```

**As a local job:**

`nextflow run main.nf -profile nscc` (if running on the NSCC server) or `nextflow run main.nf -profile local` if running locally.

## Output files

Every process will create a separate folder which contain:

- The process's output eg. fasta files.
- `nextflow.command.sh`: the script run by Nextflow.
- `nextflow.command.log`: the scripts's stderr and stdout.
- `nextflow.exitcode`: the script's exit code.

Once the pipeline has finished running, several reports are also generated inside the directory that nextflow was started from (ie. the current directory when `nextflow run ...` was issued):
- `nextflow-report.html`: a comprehensive summary of the pipeline ([Example](https://www.nextflow.io/docs/latest/tracing.html#execution-report))
- `nextflow-timeline.html`: a timeline showing the duration of each process. ([Example](https://www.nextflow.io/docs/latest/tracing.html#timeline-report))
- `nextflow-trace.tsv`: a table of information about each process. ([Example](https://www.nextflow.io/docs/latest/tracing.html#trace-report))

Note: if there are existing reports in the current directory, the existing reports will get renamed to `nextflow-report.html.1`, `nextflow-report.html.2` etc.


TODO:
- Work dir, deleting work dir
- Debugging using report, `-resume`
