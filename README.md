# Genome assembly pipeline

This is a hybrid genome assembly pipeline for bacterial genomes written in [Nextflow](https://www.nextflow.io/).

## Quickstart

### 1. Installation

1. Clone this repository
1. Install Conda and [Nextflow](https://www.nextflow.io). Make sure the `nextflow` command is available in your path.
1. Setup the conda environments from the YAML files in the `conda-envs` directory:

    a. Create the environments

   ```sh
   cd conda-envs

    # create the urops-assembly environment
    conda env create -f assembly.yml # or assembly-nscc.yml
    # create the urops-checkm environment
    conda env create -f checkm.yml
    # create the urops-circlator environment
    conda env create -f circlator.yml
   ```

   Use `assembly-nscc.yml` instead of `assembly.yml` if you are using the NSCC Aspire 1 server.
   The difference between `assembly.yml` and `assembly-nscc.yml` is that `assembly-nscc.yml` does not include BWA and Samtools, since they are already pre-installed on the server.
   If you are using the NSCC server, remember to load the `bwa/0.7.13` and `samtools/1.3` modules before running the pipeline.

   b. Set the `PILONJAR` and `PLATONDB` environment variables for the urops-assembly environment:

    ```sh
    conda activate urops-assembly
    conda env config vars set PILONJAR=<path to pilon jar file> # change this
    # example: conda env config vars set PILONJAR=~/.conda/pkgs/pilon-1.23-2/share/pilon-1.23-2/pilon-1.23.jar
    conda env config vars set PLATONDB=<path to platon database> # change this
    ```

    - `PILONJAR` is the path to the Pilon's jar file. It is usually at `<path to conda directory>/pkgs/pilon-1.23-2/share/pilon-1.23-2/pilon-1.23.jar`.

    - `PLATONDB` is the path to Platon's database, which will be downloaded in the next step.

1. Download the databases needed by various tools:
   - **Platon database:** download the database at this [link](https://zenodo.org/record/4066768/files/db.tar.gz?download=1) and extract it to the file path specified for `PLATONDB`

### 2. Run the pipeline

The several different operations which can be run, such as genome assembly, assembly evaluation, read quality assessment and dependency checking.
Each operation is located in a different `.nf` script in the root directory.

A script is launched using the `nextflow run` command from the root directory, eg. `nextflow run main.nf`.
Various command line arguments can be passed to a script (eg. `--illumina1`, `-work-dir`).

Things to take note of:
- Some flags begin with a single dash `-`, while some begin with a double dash `--`.
- The scripts need to be launched from the root directory of this repository, in order for the configurations in the `nextflow.config` file to be applied.

If you are submitting a job on the NSCC Aspire 1 server, here is a template for a job script which can be submitted using `qsub`:

``` sh
#!/bin/bash
#PBS -q normal
#PBS -P Personal
#PBS -j oe
#PBS -l select=1:ncpus=10:mpiprocs=10
#PBS -l walltime=04:00:00
module load java # nextflow is a JVM program
module load openmpi # for mpirun command

# used by pipeline
module load samtools/1.3
module load bwa/0.7.13

cd $PBS_O_WORKDIR # cd to the directory qsub was run from
cd assembly-pipeline/ # cd to the directory containing this repository
mpirun --pernode nextflow run main.nf -with-mpi -profile nscc [pipeline parameters]
```

We running nextflow via the `mpirun` command in order to take advantage of the [OpenMPI standard](https://www.open-mpi.org/) for improved performance.
You can read more about how Nextflow uses OpenMPI [here](https://www.nextflow.io/docs/latest/ignite.html?highlight=mpi#execution-with-mpi) and [here](https://www.nextflow.io/blog/2015/mpi-like-execution-with-nextflow.html).

`[pipeline parameters]` should be replaced with the command line parameters for the `main.nf` script. The available parameters can be found in the following sections:

- The [scripts section](#scripts): this includes the parameters for each script.
- [Execution related parameters](#execution-related-parameters): this describes the parameters that are common to all scripts (eg. number of threads to use). This includes a description of the `-profile` parameter.

Note that `main.nf` may be replaced with any other script, eg. `assemble.nf`. 

Otherwise, if you are running the pipeline on a desktop computer, or don't wish to use OpenMPI just use the `nextflow run` command, eg.

``` sh
nextflow run main.nf [pipeline parameters]
```

### 3. View output files

The output files generated by each script can be found in the [scripts section](#scripts).

In addition, Nextflow has also been configured to generate several reports detailing the runtime, resources used and other information:

- `nextflow-report-<timestamp>.html`: a comprehensive summary of the pipeline ([Example](https://www.nextflow.io/docs/latest/tracing.html#execution-report))
- `nextflow-timeline-<timestamp>.html`: a timeline showing the duration of each process. ([Example](https://www.nextflow.io/docs/latest/tracing.html#timeline-report))
- `nextflow-trace-<timestamp>.tsv`: a table of information about each process. ([Example](https://www.nextflow.io/docs/latest/tracing.html#trace-report))

These reports are saved to either the [outdir folder](#assemble-parameter-desc) or the directory that nextflow was launched from (if outdir is not set).

### 4. Delete Nextflow's working directory
<a id="work-dir"></a>

Nextflow stores temporary data in a working directory.
This will basically include a copy of all the output files plus other intermediate files.
This can be safely deleted after the pipeline finishes running, **unless** the pipeline failed and you want to [resume it](#resuming-a-script).

Read more at: 

- https://github.com/danrlu/Nextflow_cheatsheet#the-working-directory
- https://www.nextflow.io/docs/latest/script.html?highlight=workdir#implicit-variables

## Scripts

This section contains a description of each Nextflow script, including what it does, its parameters and output files.

_Note about format: In the detailed description section for each script, the name of Nextflow process is written in parenthesis. Eg. (`cleanShortReads`)_

### Genome assembly + evaluation: `main.nf`

Perform genome assembly using [`assemble.nf`](#genome-assembly-assemblenf) followed by assembly evaluation using [`evaluateChromosome.nf`](#assembly-evaluation-evaluatechromosomenf-evaluateplasmidnf), if the assembly contains exactly 1 contig. If there is more than 1 contig, then evaluation will not be carried out. It will be up to the user to inspect the assembly and decide how they want to split the contigs for evaluation. Evaluation can be carried out using the [`evaluateChromosome.nf` and `evaluatePlasmid.nf` scripts](#evaluation).

Command line parameters are the same as `assemble.nf`.

### Genome assembly: `assemble.nf`

Assemble a genome from raw Illumina and Pacbio reads. The final assembly will be in the file `assembly/pilon/final_pilon_assembly.fa`. Some basic statistics (number of contigs, size of contigs, and circularity of contigs) will be given in the file `assembly/assembly-summary.json`.

``` sh
nextflow run assemble.nf --illumina1 <path> --illumina2 <path> --pacbio <path> --outdir <path> \
    [--shortReadsKeepPercent <percent>] [--shortReadsStartTrimq <percent>] [--shortReadsMinTrimq <percent>] \
    [--bbdukArgs <args>] [--filtlongArgs <args>] [--filtlongCheckThreshold <threshold>] \
    [--flyeArgs <args>] \
    [--raconMaxIters <number>] [--raconArgs <args>] [--pilonMaxIters <number>] [--pilonArgs <args>] \
    [--canuGenomeSize <size>] [--canuArgs <args>] [--circlatorArgs <args>] [--forceCirclator | --noCirclator] \
    # execution-related params
    [--skipDepChecks] [--threads <number>] [-work-dir <path>] [--condaEnvsDir <path>] [-profile <profiles>]
```

#### Detailed description

1. Check that all programs are working (`checkDependencies`).
1. Clean short reads using Bbduk (`cleanShortReads`)
1. Clean long reads using Filtlong (`cleanLongReads`)
1. Form an initial assembly using the long reads with Flye (`flyeAssembly`)
1. Polish the assembly using the long reads with Racon (`raconPolish`)
   - Racon is run up to 4 times or until there are no changes.
1. If the assembly is potentially circular, then attempt circularisation using Circlator (`circlator`). Canu corrected long reads (`canuCorrect`) are given to Circlator as input.
   - An assembly is considered potentially circular if Flye identifies any circular contigs or there are multiple linear contigs.
   - You can force Circlator to always run / never run using `--forceCirclator` / `--noCirclator`.
1. Polish the assembly using Pilon (`pilonPolish`)
   - Pilon is run up to 6 times or until there are no changes.
1. Save a summary of some basic statistics about the assembly (number of contigs, size of contigs, and circularity of contigs) to the file `assembly-summary.json` (`summariseAssembly`)

#### Output files

- `<outdir>/`
    - `assembly/`
        - `flye/`
          - Files in the Flye root directory.
          - `22-plasmids/`
        - `racon/`
            - `final_racon_assembly.fa`: assembly after polishing with Racon for several iterations.
            - `final_racon_log.tsv`: whether Racon made a difference, and number of Racon iterations. (TODO)
        - `circlator/` (if Circlator was run)
            - select Circlator output files.
        - `pilon/`
            - **`final_pilon_assembly.fa`: the final assembly**
            - `pilonX.changes`: changes made during each Pilon iteration
            - `pilon_info.tsv`: number of Pilon iterations, and whether or not Pilon converged.
        - `assembly-summary.json`
    - `reads/`
        - `short_cleaned/`
            - `illumina1.fq`: 1st set of cleaned Illumina reads.
            - `illumina2.fq`: 2nd set of cleaned Illumina reads.
            - `trimq_used.txt`: trimq score used by Bbduk.
        - `long_cleaned/`
            - `pacbio.fq`: cleaned pacbio reads
            - `above_10kb_reads_removed.tsv`: list of long reads above 10kb in size which where removed by Filtlong.
        - `long_canu/` (if canu was run)
            - `canu.correctedReads.fasta.gz`
            - Select Canu output files.

<a id="other-outputs"></a>

In addition to the files created specifically by each process, the stdout, stderr, exit code and exact script run are also saved:
- `nextflow.command.sh`: the script run by Nextflow.
- `nextflow.command.log`: the scripts's stderr and stdout.
- `nextflow.exitcode`: the script's exit code.

#### Parameter descriptions

<a id="assemble-parameter-desc"></a>

**Required parameters:**

- `--illumina1 <path>`: Path to 1st file for raw paired end Illumina reads. May be fastq or gzipped fastq.
- `--illumina2 <path>`: Path to 2nd file for raw paired end Illumina reads. May be fastq or gzipped fastq.
- `--pacbio <path>`: Path to the raw Pacbio reads. May be fastq or gzipped fastq.
- `--outdir <path>`: Path to the output directory, which is where all output files will be stored.

**Optional parameters:**

- Short read filtering and cleaning:
    - `--shortReadsKeepPercent <percent>`: Ensure at least X% of reads are kept. Default: 80.
    - `--shortReadsStartTrimq <trimq>`: Highest possible `trimq` value for bbduk. Default: 40.
    - `--shortReadsMinTrimq <trimq>`: Lowest permissible `trimq` value for bbduk. Default: 28.
    - `--bbdukArgs <args>`: Arguments (other than inputs, outputs and `trimq`) to pass to bbduk. Default: `qtrim=rl minlength=40`.
- Long read filtering and cleaning:
    - `--filtlongArgs <args>`: Arguments (other than inputs and outputs) to pass to Filtlong. Default: `--min_length 1000 --keep_percent 90 --trim --split 500 --mean_q_weight 10`.
    - `--filtlongCheckThreshold <number>[k|m|g]`: Flag reads above the given length (eg. 10, 10k, 10m, 10g) which were removed by Filtlong. The number of these reads will be printed to stdout and the IDs and lengths of these reads will be saved to a TSV file. Default: 10k.
- Flye:
    - `--flyeArgs <args>`: Arguments (other than inputs, outputs, threads and `--plasmids`) to pass to Flye. Default: none.
- Racon:
    - `--raconMaxIters <number>`: Maximum number of iterations to run Racon for. Default: 4.
    - `--raconArgs <args>`: Arguments (other than inputs, outputs, and threads) to pass to Racon. Default: `-m 8 -x -6 -g -8 -w 500`.
- Pilon:
    - `--pilonMaxIters <number>`: Maximum number of iterations to run Pilon for. Default: 6.
    - `--pilonArgs <args>`: Arguments (other than inputs, outputs, and  `--changes`) to pass to Pilon. Default: none.
- Circularisation:
    - `--canuGenomeSize <genome size>`: When specified, force Canu to use this genome size. See the Canu documentation for genomeSize for valid values. Otherwise, calculate the genome size from the assembly. Default: not specified.
    - `--canuArgs <args>`: Arguments (other than inputs, outputs and genome size) to pass to Canu. Default: none.
    - `--circlatorArgs <args>`: Arguments (other than inputs and outputs) to pass to Circlator. Default: none.
    - `--forceCirclator`: Force Circlator to be used, regardless of the state of the Flye assembly.
    - `--noCirclator`: Do not use Circlator, regardless of the state of the Flye assembly.
- Also see [execution related parameters](#execution-related-parameters)

<a id="evaluation"></a>

### Assembly evaluation: `evaluateChromosome.nf`, `evaluatePlasmid.nf`

To evaluate a chromosome:

``` sh
nextflow run evaluateChromosome.nf --illumina1 <path> --illumina2 <path> --pacbio <path> --assembly <path> \
    --outdir <path> \
    # execution-related params
    [--skipDepChecks] [--threads <number>] [-work-dir <path>] [--condaEnvsDir <path>] [-profile <profiles>]
```

To evaluate a plasmid, replace `evaluationChromosome.nf` with `evaluatePlasmid.nf`.

The only difference between the chromosome and plasmid evaluation is CheckM is not run for plasmids.


#### Detailed description

Before performing the analysis, the script checks that all the programs needed are working (`checkDependencies`).

The following metrics about the given chromosome/plasmid assembly are recorded:

- Coverage of short reads and long reads (`shortReadsCoverage`, `longReadsCoverage`)
- Annotations using Prokka (`prokkaAnnotate`)
- Statistics about the assembly itself eg. length, number of contigs, N50 using Quast (`quastEvaluate`)
- (chromosome only) Completeness and contamination using CheckM (`checkmEvaluate`)

The most important statistics are also saved to a single summary document, `chromosome-summary.json` or `plasmid-summary.json` (`makeChromosomeSummary`, `makePlasmidSummary`):

```
{
    "avg short reads coverage": {
        "contig_1": 372.3411
    },
    "short reads mapped": {
        "total reads": 9397580,
        "mapped reads:": 207163,
        "proportion mapped": 0.022044292253963253
    },
    "avg long reads coverage": {
        "contig_1": 210.1854
    },
    "long reads mapped": {
        "total reads": 86864,
        "mapped reads:": 3390,
        "proportion mapped": 0.03902652422177197
    },
    "assembly length": 69354,
    "contigs": 1,
    "N50": 69354,
    "GC": 0.7051000000000001,
    "CDS": 79,
    "size by contig": {
        "contig_1": 69354
    },
    "errors": [],
    "completeness": 0.0,
    "contamination": 0.0
}
```

`errors` is a list of any errors encountered while creating the summary document (does not include the errors when running quast/prokka etc.)

#### Output files

- `<outdir>/`
  - `short_read_coverage/`
  - `long_read_coverage/`
  - `prokka/`
  - `quast/`
  - `checkm/` (for chromosome evaluation only)
  - `chromsome-summary.json` / `plasmid-summary.json`

As with other scripts, `nextflow.command.sh`, `nextflow.command.log` and `nextflow.exitcode` is saved for every process. Read more [here](#other-outputs).

#### Parameter descriptions

**Required parameters:**

- `--illumina1 <path>`: Path to 1st file for **cleaned** paired end Illumina reads. May be fastq or gzipped fastq.
- `--illumina2 <path>`: Path to 2nd file for **cleaned** paired end Illumina reads. May be fastq or gzipped fastq.
- `--pacbio <path>`: Path to the **cleaned** Pacbio reads. May be fastq or gzipped fastq.
- `--assembly <path>`: Path to fasta file of chromosome/plasmid assembly.
- `--outdir <path>`: Path to the output directory, which is where all output files will be stored.

**Optional parameters:**

- See [execution related parameters](#execution-related-parameters)


### Check all dependencies: `checkAllDependencies.nf`

If you want to manually check all the programs are working correctly:

``` sh
nextflow run checkAllDependencies.nf [--threads <number>] [-work-dir <path>] [--condaEnvsDir <path>] [-profile <profiles>]
```

See [execution related parameters](#execution-related-parameters) for descriptions of the optional parameters.

### Assess read quality: `assessReads.nf`

TODO document

## Execution related parameters

These are parameters which do not control what is run, only how the pipeline is run (eg. number of threads). For the most part, they are common to all scripts.

- `--threads <num threads>`: Number of threads a process should use. Default: 1.
- `--condaEnvsDir <path>`: Path to where the conda environments are stored. Default: `~/.conda/envs/`.
- `-work-dir`: (note that there is only a single `-` at the front) Path to Nextflow's working directory, which is where temporary files will be stored. Default: `./work/`.
    - Read more [here](#work-dir)
- `--skipDepChecks`: Skip pre-pipeline dependency checks. Does not apply to the `checkAllDependencies.nf` script.
- `-profile <profiles>`: (note that there is only a single `-` at the front) List of [Nextflow configuration profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles) to apply. Options:
    - `local`: Applies settings for working on a normal laptop/desktop computer: `--threads 4`.
    - `nscc`: Applies settings for working on the NSCC Aspire 1 sever: `-work-dir ~/scratch/work/ --threads 20`.


## Troubleshooting

### Resuming a script

If a script terminated due to an error, you can resume the script from where it left off after fixing the error using the `-resume` flag. Note that you must **not** have deleted the working directory in order for resume to work.

Read more at:
- https://www.nextflow.io/docs/latest/getstarted.html?highlight=resume#modify-and-resume
- https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html
- https://www.nextflow.io/blog/2019/troubleshooting-nextflow-resume.html
