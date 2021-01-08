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
   - Racon is run up to 4 times or until there are no changes.
1. If the genome is circular according to Flye, then circularise the assembly using Circlator (`circlator`)
   - In addition to the assembly, Circlator requires long reads as input. We use long reads corrected using Canu (`canuCorrect`).
1. Polish the assembly using Pilon (`pilonPolish`)
   - Pilon is run up to 6 times or until there are no changes.

### Assembly evaluation

Evaluation is done separately for chromosomes and plasmids.

**Metrics recorded:**

- Coverage of short reads and long reads (`shortReadsCoverage`, `longReadsCoverage`)
- Annotations using Prokka (`prokkaAnnotate`)
- Statistics about the assembly itself eg. length, number of contigs, N50 using Quast (`quastEvaluate`)
- (chromosome only) Completeness and contamination using CheckM (`checkmEvaluate`)

**Summary document:**

The most important statistics are also saved to a single summary document, `chromosome-summary.json` or `plasmid-summary.json` (`makeChromosomeSummary`, `makePlasmidSummary`):

```
{
    "avg short reads coverage": 273.4589,
    "short reads mapped": {
        "total reads": 9428394,
        "mapped reads:": 9424508,
        "proportion mapped reads": 0.9995878407287604
    },
    "avg long reads coverage": 161.3932,
    "long reads mapped": {
        "total reads": 86563,
        "mapped reads:": 86523,
        "proportion mapped reads": 0.9995379088063029
    },
    "length": 4159544,
    "contigs": 1,
    "N50": 4159544,
    "GC": 0.6025,
    "CDS": 3868,
    "is circular": {
        "contig_1": true
    },
    "completeness": 1.0,
    "contamination": 0.0002787844995818232,
    "errors": []
}
```

Note that `errors` refers to any errors encountered while creating the summary document (does not include the errors when running quast/prokka etc.)

## Output files

Every process will create a separate folder which contain:

- The process's output eg. fasta files.
- `nextflow.command.sh`: the script run by Nextflow.
- `nextflow.command.log`: the scripts's stderr and stdout.
- `nextflow.exitcode`: the script's exit code.

The name and location of the folder for each process is determined by the `params.outdir` and `params.o...` (eg. `params.o.cleanShortReads`) parameters defined at the top of `main.nf`.

Several reports are also generated inside the directory that nextflow was started from (ie. the current directory when `nextflow run ...` was issued):
- `nextflow-report.html`: a comprehensive summary of the pipeline ([Example](https://www.nextflow.io/docs/latest/tracing.html#execution-report))
- `nextflow-timeline.html`: a timeline showing the duration of each process. ([Example](https://www.nextflow.io/docs/latest/tracing.html#timeline-report))
- `nextflow-trace.tsv`: a table of information about each process. ([Example](https://www.nextflow.io/docs/latest/tracing.html#trace-report))

Note: if there are existing reports in the current directory, the existing reports will get renamed to `nextflow-report.html.1`, `nextflow-report.html.2` etc.

**(TODO) Nextflow's working directory (workDir):**

Nextflow stores temporary data in a working directory.
This will basically includes all of the files mentioned above and possibly additional intermediate files.
This can be safely deleted after the pipeline finishes running.
See more at: https://www.nextflow.io/docs/latest/script.html?highlight=workdir#implicit-variables


## Setup

1. Create the conda environments from yaml files in the `conda-envs` directory
   - Use only one of `assembly.yml` and `assembly-nscc.yml`. The difference is that `assembly.yml` includes samtools and bwa while `assembly-nscc.yml` does not.
   - Before creating `assembly.yml`/`assembly-nscc.yml`, change the following variables:
        - Set `PILONJAR` to the path to the Pilon's jar file. It is usually at `<path to conda directory>/pkgs/pilon-1.23-2/share/pilon-1.23-2/pilon-1.23.jar`.
        - Set `PLATONDB` to the path to Platon's database, which will be donwloaded in the next step.
1. Download databases needed by various tools
   - Platon: download the database at this [link](https://zenodo.org/record/4066768/files/db.tar.gz?download=1) and extract it to `PLATONDB`
   - Checkm: 
       1. Download and unzip the database at this [link](https://data.ace.uq.edu.au/public/CheckM_databases/)
       1. Tell Checkm where the database is:
          ```sh
          conda activate urops-checkm
          checkm data setRoot <path_to_checkm_database>
          ```
1. Install [Nextflow](https://www.nextflow.io). It has been tested using v20.10 so far (latest stable release as of time of writing).

## Running the pipeline

There are a few ways to run the pipeline:

### (A) As a job on the NSCC Aspire 1 server:

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
mpirun --pernode nextflow run main.nf -profile nscc -with-mpi --rawIllumina1 <path> --rawIllumina2 <path> --rawPacbio <path> --outdir <path> --condaEnvsDir <path>
```

Appropriate paths need to be passed to `--rawIllumina1`, `--rawIllumina2`, `--rawPacbio` and `--condaEnvsDir`. See the section [Setting Parameters](#setting-parameters) for more details as well as other parameters which can be set.

### (B) As a local job:

`nextflow run main.nf -profile nscc --rawIllumina1 <path> --rawIllumina2 <path> --rawPacbio <path> --outdir <path> --condaEnvsDir <path>` (if running on the NSCC server) or `nextflow run main.nf -profile local --rawIllumina1 <path> --rawIllumina2 <path> --rawPacbio <path> --outdir <path> --condaEnvsDir` (if running on desktop).

Appropriate paths need to be passed to `--rawIllumina1`, `--rawIllumina2`, `--rawPacbio` and `--condaEnvsDir`. See the section [Setting Parameters](#setting-parameters) for more details as well as other parameters which can be set.

### Setting Parameters

Several parameters can be set as a command line flag, eg. `nextflow run ... --outdir /path/to/somewhere` (see the [Nextflow docs](https://www.nextflow.io/docs/latest/getstarted.html?highlight=params#pipeline-parameters) for more info on how this is implemented). They are listed below:

**Inputs and outputs:**

- `rawIllumina1`: (required) Path to one set of raw Illumina reads. May be fastq or gzipped fastq.
- `rawIllumina2`: (required) Path to the other set of raw Illumina reads. May be fastq or gzipped fastq.
- `rawPacbio`: (required) Path to the raw Pacbio reads. May be fastq or gzipped fastq.
- `outdir`: (required) Path to the output directory, which is where all output files will be stored.

**Execution related:**

- `condaEnvsDir`: (required) Path to where the conda environments are stored. Usually it's `~/.conda/envs/`.
- `threads`: number of threads a process should use.

**Pipeline parameters:**

- Short read filtering and cleaning using BBDuk:
    - `bbdukKeepPercent`: ensure at least X% of reads are kept. Default: 80. \*
    - `bbdukStartTrimq`: highest possible `trimq` value. Default: 40. \*
    - `bbdukMinTrimq`: lowest permissible `trimq` value. Default: 28. \*
    - `bbdukArgs`: arguments other than inputs, outputs and `trimq` to pass to bbduk. Default: `qtrim=rl minlength=40`.
    - \* : You may have noticed that these aren't arguments that are part of bbduk. They're part of a script which helps finds the best possible `trimq` value to use with bbduk such that at least X% of reads are kept, similar to Filtlong. ([here](https://github.com/chloelee767/assembly-pipeline/blob/master/bin/bbduk_keep_percent.py))
- `pilonMaxIters`: maximum number of iterations to run Pilon for. Default: 6.
- `raconMaxIters`: maximum number of iterations to run Racon for. Default: 4.


## TODO debugging
TODO:
- Debugging using report, `-resume`

## Understanding the effect of each assembly step

To understand the changes each step makes on the genome assembly, look at the following files in the process' output directory:

- **Racon:** 
    - `final_assembly_log.tsv` gives the difference between assembly before and after running racon, as well as the number of times that racon was run.
- **Circlator:**
    - The `.log` files: describes what circlator did in each step. See the [Circlator wiki](https://github.com/sanger-pathogens/circlator/wiki) for documentation on the contents of the log file for each task.
- **Pilon:**
    - The number of iterations of pilon can be inferred from the numbering of the pilon files.
    - The exact changes made in each iteration of pilon is given in `pilonX.changes`. See the [Pilon wiki](https://github.com/broadinstitute/pilon/wiki/Output-File-Descriptions#changes) for the format.
