# Genome assembly pipeline

- [Introduction](#introduction)
- [Tools used](#tools-used)
- [Running time](#running-time)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick usage](#quick-usage)
  * [1. Evaluate the reads](#1-evaluate-the-reads)
  * [2. Assemble and evaluate the reads](#2-assemble-and-evaluate-the-reads)
    + [Adjusting parameters](#adjusting-parameters)
- [Detailed documentation](#detailed-documentation)
  * [Documentation for Nextflow scripts](#documentation-for-nextflow-scripts)
    + [Genome assembly + evaluation: `main.nf`](#genome-assembly--evaluation-mainnf)
    + [Genome assembly: `assemble.nf`](#genome-assembly-assemblenf)
      - [Parameter descriptions](#parameter-descriptions)
      - [Detailed description](#detailed-description)
      - [Output files](#output-files)
    + [Assembly evaluation: `evaluateChromosome.nf`, `evaluatePlasmid.nf`](#assembly-evaluation-evaluatechromosomenf-evaluateplasmidnf)
      - [Parameter descriptions](#parameter-descriptions-1)
      - [Detailed description](#detailed-description-1)
      - [Output files](#output-files-1)
    + [Check dependencies: `checkDependencies.nf`](#check-dependencies-checkdependenciesnf)
  * [Execution related parameters](#execution-related-parameters)
- [Troubleshooting](#troubleshooting)
  * [Help! Where are all my files inside the `work/` directory?](#help-where-are-all-my-files-inside-the-work-directory)
  * [Resuming a script](#resuming-a-script)
  * [Finding plasmids](#finding-plasmids)
- [Tips](#tips)
  * [Periodically clear Nextflow's work directory](#periodically-clear-nextflows-work-directory)
  * [Running the pipeline on clusters supporting OpenMPI](#running-the-pipeline-on-clusters-supporting-openmpi)
  * [Viewing the assembly graph](#viewing-the-assembly-graph)
  * [Repeated configs](#repeated-configs)

## Introduction

This is a hybrid genome assembly pipeline for bacterial genomes written in Nextflow.
It requires (1) paired end Illumina short reads and (2) either Pacbio or Nanopore long reads as input.

Both genome assembly as well as assembly evaluation are performed.
Some major features are: stringent, but configurable read filtering criteria; detection of plasmids and possible contamination; automated invocation of assembly evaluation when possible.
The pipeline tries to set reasonable defaults as much as possible, but most parameters can be adjusted easily if need be via command line options.

## Tools used

- Read filtering
  - Long reads: Filtlong
  - Short reads: BBduk
- Read correction: Canu (only used to generate the input for Circlator)
- Genome assembly
  - Assembly: Flye
  - Polishing: Racon, Pilon
  - Circularisation: Circlator
  - Plasmid detection: Platon
  - Misc: Minimap2, samtools, BWA
- Assembly evaluation
  - Read coverage: 
    - Short reads: BBmap 
    - Long reads: Minimap2, Pileup (from BBtools)
  - Statistics: Quast
  - Annotation: Prokka
  - Completeness and contamination: Checkm

## Running time

- Genome assembly: 3-4 hours using `--threads 20` on a compute cluster with Intel Xeon E5-2690 v3 processors (2.60GHz-3.50GHz)
- Assembly evaluation: 10-30 minutes using `--threads 20` on a compute cluster with Intel Xeon E5-2690 v3 processors (2.60GHz-3.50GHz)

## Requirements

- Linux
- Java (use Java version required by Nextflow)
- [Nextflow](https://www.nextflow.io/) v20.10 or later
- Conda. Preferably, use the current version, however older versions which use `conda activate` (rather than `source activate`) ought to work as well.
- Git (optional, but recommended)
    - Git will be used to clone the respository. It is also possible to download a zip of the respository instead, but it is preferable to use git as it is easier to download updates for the pipeline in the future.

## Installation

1. Install the dependencies listed in the [requirements](#requirements) section.
2. Download/clone this repository.
   If you have git, run the command below to clone the repository to the directory `assembly-pipeline`:
``` sh
git clone https://github.com/chloelee767/assembly-pipeline.git
```
Otherwise, download a zip file of the code and extract it to your desired location

3. Create the conda environments from the YAML files in the `conda-envs` directory:

    a. Create the environments

   ```sh
   # cd to the assembly-pipeline directory containing the code
   cd path/to/assembly-pipeline # change this

   cd conda-envs

    # create the ap-main environment
    conda env create -f main.yml
    # create the ap-checkm environment
    conda env create -f checkm.yml
    # create the ap-circlator environment
    conda env create -f circlator.yml
   ```

   b. Set the `PILONJAR` and `PLATONDB` environment variables for the `ap-main` environment:
   
   Make sure to use the full path (eg. `/home/your_username/path/to/file` not `~/path/to/file`), as conda will set the environment variables verbatim without path name expansion. 

    ```sh
    conda activate ap-main
    conda env config vars set PILONJAR=<path to pilon jar file> # change this
    # example: conda env config vars set PILONJAR=/home/your_username/.conda/pkgs/pilon-1.23-2/share/pilon-1.23-2/pilon-1.23.jar
    conda env config vars set PLATONDB=<path to platon database> # change this
    ```

    - `PILONJAR` is the path to the Pilon's jar file. It is usually at `<path to conda directory>/pkgs/pilon-1.23-2/share/pilon-1.23-2/pilon-1.23.jar`.

    - `PLATONDB` is the path to Platon's database, which will be downloaded in the next step.

4. Download the databases needed by various tools:
   - **Platon database:** download the database at this [link](https://zenodo.org/record/4066768/files/db.tar.gz?download=1) and extract it to the file path specified for `PLATONDB`

5. Check that all programs have been installed correctly:
    ```sh
    nextflow run checkDependencies.nf
    ```
    If your conda environments are stored somewhere other than `~/.conda/envs`, add the argument `--condaEnvsDir <path>`, replacing `<path>` with the path to the conda environments. 
    
    Check the stdout to see which dependency checks have failed. Note that `checkDependencies.nf` script itself will **NOT** fail and give a non-zero exit code, it will only print the errors to the command line.
  
    Example of the output when a dependency check has failed:
    ```
    N E X T F L O W  ~  version 20.10.0
    Launching `checkDependencies.nf` [magical_brattain] - revision: 77e3ff929e
    executor >  local (20)
    [73/0e90af] process > checkAllDependencies:testBbduk               [100%] 1 of 1 ✔
    [2d/f9b5b9] process > checkAllDependencies:testFiltlong            [100%] 1 of 1 ✔
    [c5/2b477a] process > checkAllDependencies:testFlye                [100%] 1 of 1 ✔
    [f5/fc9d37] process > checkAllDependencies:testBwa                 [100%] 1 of 1 ✔
    [47/ce0101] process > checkAllDependencies:testPilon               [100%] 1 of 1 ✔
    [02/fd44c6] process > checkAllDependencies:testPlaton              [100%] 1 of 1 ✔
    [f9/24d6a7] process > checkAllDependencies:testRacon               [100%] 1 of 1 ✔
    [2d/ddfc29] process > checkAllDependencies:testCanu                [100%] 1 of 1 ✔
    [b1/742dfb] process > checkAllDependencies:testSamtools            [100%] 1 of 1 ✔
    [8b/9a6576] process > checkAllDependencies:testMinimap2            [100%] 1 of 1 ✔
    [2b/be540a] process > checkAllDependencies:testBbmap               [100%] 1 of 1 ✔
    [ec/b2e849] process > checkAllDependencies:testPileup              [100%] 1 of 1 ✔
    [3e/474dec] process > checkAllDependencies:testProkka              [100%] 1 of 1 ✔
    [5e/008de6] process > checkAllDependencies:testQuast               [100%] 1 of 1 ✔
    [50/8e6ff6] process > checkAllDependencies:testPython_mainEnv      [100%] 1 of 1 ✔
    [60/336736] process > checkAllDependencies:testCirclator           [100%] 1 of 1 ✔
    [d8/f9cb08] process > checkAllDependencies:testPython_circlatorEnv [100%] 1 of 1 ✔
    [8b/0825d0] process > checkAllDependencies:testCheckm              [100%] 1 of 1 ✔
    [25/18e33e] process > checkAllDependencies:testPython_checkmEnv    [100%] 1 of 1 ✔
    [68/fcef28] process > checkAllDependencies:dummyTest               [100%] 1 of 1, failed: 1 ✔
    Completed at: 09-Feb-2021 17:09:04
    Duration    : 2m 27s
    CPU hours   : 0.1 (0.3% failed)
    Succeeded   : 19
    Ignored     : 1
    Failed      : 1
    ```
    
  In this case, the dummyTest dependency check has failed. 
  To figure out what caused the dummyTest to fail, you can check what was printed to the command line when dummyTest was run by opening the file `work/68/fcef28.../.command.log` (if you changed the working directory, replace `work` with the path to your working directory). Note that `fcef28...` refers to a folder whose name starts with `fcef28`.
  

## Quick usage

### 1. Evaluate the reads

Use Fastqc and Nanoplot (or any programs) to evaluate the quality of short reads and long reads (respectively) to decide if the [default criteria used for read filtering](#assemble-parameter-desc) are appropriate.

Fastqc and Nanoplot are included in the `ap-main` conda environment, so just run `conda activate ap-main` and you will be able to run the `fastqc` and `NanoPlot` commands.

### 2. Assemble and evaluate the reads

The command below runs the full pipeline (genome assembly + if possible, evaluation of the chromosome + any plasmids) with the default settings:

``` sh
cd assembly-pipeline/ # cd to the directory containing this repository
nextflow run main.nf --illumina1 short_reads_1.fq.gz --illumina2 short_reads_2.fq.gz \ 
    --pacbio long_reads.fq.gz --outdir assembly-results
```

If you are using Nanopore long reads instead, replace `--pacbio` with `--nanopore`.

If there are multiple contigs in the assembly, evaluation will not be carried out automatically and a message will be printed to the command line. Multiple contigs may mean that there are multiple bacterial species and hence multiple chromosomes and/or there are plasmids. To evaluate the assembly, you will have to:

1. Separate out each chromosome and plasmid from the final assembly file into separate files. You may use Platon, which is included in one of the conda environments, to help identify plasmids. See [this section](#finding-plasmids) for more details.

2. Invoke the evaluation script manually for each chromosome and plasmid using the [`evaluateChromosome.nf` and `evaluatePlasmid.nf`](#evaluation) scripts.

It is also possible to run just genome assembly, without evaluation, using the [`assemble.nf` script](#genome-assembly-assemblenf).

#### Adjusting parameters

Additional command line flags may be passed to adjust the criteria used for read filtering, based on step 1. 
Check out the documentation [here](#assemble-parameter-desc).

Depending on the machine you are using, some settings that commonly need to be adjusted are:
- **Path to the directory containing the conda environments** (default: `~/.conda/envs/`). To change this, add this command line flag: `--condaEnvsDir <path>` eg. `--condaEnvsDir /opt/conda/envs/`.
- **Memory to allocate to the JVM when running Pilon** (default: 13GB). To change this, add this command line flag: `--pilonMemory <memory>`. See the documentation for the `java -Xmx` command for valid values for `<memory>`.
- **Number of threads** to use for multithreaded programs (default: 1). To change this, add this command line flag: `--threads <number>`.
- **Directory where temporary files are stored** (default: `work/`). To change this, add this command line flag: `-work-dir <path>`. If you are running the pipeline on a cluster with limited permanent storage, you may want to store temporary files in the scratch/temp directory instead. 

Some things to take note of:
- Some flags begin with a single dash `-`, while some begin with a double dash `--`.
- The scripts need to be launched from the root directory of this repository, in order for the configurations in the `nextflow.config` file to be applied.

## Detailed documentation

### Documentation for Nextflow scripts

This section contains a description of each Nextflow script, including what it does, its parameters and output files.

_Note about format: In the detailed description section for each script, the name of Nextflow process is written in parenthesis. Eg. (`cleanShortReads`)_

#### Genome assembly + evaluation: `main.nf`

Perform genome assembly using [`assemble.nf`](#genome-assembly-assemblenf) followed by assembly evaluation using [`evaluateChromosome.nf`](#assembly-evaluation-evaluatechromosomenf-evaluateplasmidnf), if the assembly contains exactly 1 contig. If there is more than 1 contig, then evaluation will not be carried out. It will be up to the user to inspect the assembly and decide how they want to split the contigs for evaluation. Evaluation can be carried out using the [`evaluateChromosome.nf` and `evaluatePlasmid.nf` scripts](#evaluation).

**Nextflow Reports:**
<a id="nextflow-reports"></a>

In addition to the output files from `assemble.nf` and `evaluateChromsome.nf`, Nextflow also generates several reports summarising the pipeline execution, ie. what exactly it did when you run `nextflow run main.nf ...`  and useful information such as resource usage.

These files will be found in the `<outdir>`. They are listed below:

- `nextflow-report.html`: A document summarising the pipeline execution. Includes the exact command run, as well as resource usage, runtime and other information for each process. See [here](https://www.nextflow.io/docs/latest/tracing.html#execution-report) for a sample report.
- `nextflow-timeline.html`: A chart showing when each process ran and how long they took.
- `nextflow-trace.tsv`: Table of runtime, resource usage among other information, for each process. This is included in `nextflow-report.html`.


Command line parameters are the same as `assemble.nf`.

#### Genome assembly: `assemble.nf`

Assemble a genome from raw Illumina and Pacbio reads. The final assembly will be in the file `assembly/pilon/final_pilon_assembly.fa`. Some basic statistics (number of contigs, size of contigs, and circularity of contigs) will be given in the file `assembly/assembly-summary.json`.

``` sh
nextflow run assemble.nf --illumina1 <path> --illumina2 <path> (--pacbio | --nanopore) <path> --outdir <path> \
    [--shortReadsKeepPercent <percent>] [--shortReadsStartTrimq <percent>] [--shortReadsMinTrimq <percent>] \
    [--bbdukArgs <args>] [--filtlongArgs <args>] [--filtlongCheckThreshold <threshold>] \
    [--flyeArgs <args>] \
    [--raconMaxIters <number>] [--raconArgs <args>] [--pilonMaxIters <number>] [--pilonArgs <args>] \
    [--canuGenomeSize <size>] [--canuArgs <args>] [--circlatorArgs <args>] [--forceCirclator | --noCirclator] \
    # execution-related params
    [--threads <number>] [-work-dir <path>] [--condaEnvsDir <path>]
```

##### Parameter descriptions

<a id="assemble-parameter-desc"></a>

**Required parameters:**

- `--illumina1 <path>`: Path to 1st file for raw paired end Illumina reads. May be fastq or gzipped fastq.
- `--illumina2 <path>`: Path to 2nd file for raw paired end Illumina reads. May be fastq or gzipped fastq.
- Long reads. Must be one of:
  - `--pacbio <path>`: Path to the raw Pacbio reads. May be fastq or gzipped fastq.
  - `--nanopore <path>`: Path to the raw Nanopore reads. May be fastq or gzipped fastq.
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
    - `--pilonMemory <memory>`: Amount of memory to allocate to the JVM (via the `java -Xmx` flag) when running Pilon. See the documentation for the `java -Xmx` flag for valid values for `<memory>`. (default: `13G`)
- Circularisation:
    - `--canuGenomeSize <genome size>`: When specified, force Canu to use this genome size. See the Canu documentation for genomeSize for valid values. Otherwise, calculate the genome size from the assembly. Default: not specified.
    - `--canuArgs <args>`: Arguments (other than inputs, outputs and genome size) to pass to Canu. Default: none.
    - `--circlatorArgs <args>`: Arguments (other than inputs and outputs) to pass to Circlator. Default: none.
    - `--forceCirclator`: Force Circlator to be used, regardless of the state of the Flye assembly.
    - `--noCirclator`: Do not use Circlator, regardless of the state of the Flye assembly.
- Also see [execution related parameters](#execution-related-parameters)

##### Detailed description

1. Clean short reads using Bbduk (`cleanShortReads`)
   - Note that we have extended Bbduk so that it can find the best the trimq value to use such that at least X % of the reads are kept.
1. Clean long reads using Filtlong (`cleanLongReads`)
1. Form an initial assembly using the long reads with Flye (`flyeAssembly`)
1. Polish the assembly using the long reads with Racon (`raconPolish`)
   - Racon is run up to 4 times or until there are no changes.
1. If the assembly is potentially circular, then attempt circularisation using Circlator (`circlator`). Canu corrected long reads (`canuCorrect`) are given to Circlator as input.
   - An assembly is considered potentially circular if Flye identifies any circular contigs or there are multiple linear contigs.
   - You can force Circlator to always run / never run using `--forceCirclator` / `--noCirclator`.
1. Polish the assembly using Pilon (`pilonPolish`)
   - Pilon is run up to 6 times or until there are no changes.
1. Save a summary of some basic statistics about the assembly (number of contigs, size of contigs, and circularity of contigs, whether Flye _may_ have found contigs) to the file `assembly-summary.json` (`summariseAssembly`)

##### Output files

- `<outdir>/`
    - `assembly-summary.json`
    - `assembly/`
        - `flye/`
            - Files in the Flye root directory.
            - `22-plasmids/`
            - \* Nextflow files
        - `racon/`
            - `final_racon_assembly.fa`: assembly after polishing with Racon for several iterations.
            - `final_racon_log.tsv`: whether Racon made a difference, and number of Racon iterations.
            - \* Nextflow files
        - `circlator/` (if Circlator was run)
            - select Circlator output files
            - \* Nextflow files
        - `pilon/`
            - **`final_pilon_assembly.fa`: the final assembly**
            - `pilonX.changes`: changes made during each Pilon iteration
            - `pilon_info.tsv`: number of Pilon iterations, and whether or not Pilon converged.
            - \* Nextflow files
    - `reads/`
        - `short_cleaned/`
            - `illumina1.fq`: 1st set of cleaned Illumina reads.
            - `illumina2.fq`: 2nd set of cleaned Illumina reads.
            - `trimq_used.txt`: trimq score used by Bbduk.
            - \* Nextflow files
        - `long_cleaned/`
            - `pacbio.fq`: cleaned pacbio reads
            - `above_10kb_reads_removed.tsv`: list of long reads above 10kb in size which where removed by Filtlong.
            - \* Nextflow files
        - `long_canu/` (if canu was run)
            - `canu.correctedReads.fasta.gz`
            - Select Canu output files.
            - \* Nextflow files

<a id="other-outputs"></a>

**\* Nextflow files:**

In addition to the files created specifically by each process, the stdout, stderr, exit code and exact script run are also saved:
- `nextflow.command.sh`: the script run by Nextflow.
- `nextflow.command.log`: the scripts's stderr and stdout.
- `nextflow.exitcode`: the script's exit code.

**Nextflow Reports:**

As with `main.nf`, if `assemble.nf` was run alone (ie. you did `nextflow run assemble.nf ...`), the usaul [Nextflow reports](#nextflow-reports) will be created. 

<a id="evaluation"></a>

#### Assembly evaluation: `evaluateChromosome.nf`, `evaluatePlasmid.nf`

To evaluate a chromosome:

``` sh
nextflow run evaluateChromosome.nf --illumina1 <path> --illumina2 <path> --longReads <path> --assembly <path> \
    --outdir <path> \
    # execution-related params
    [--threads <number>] [-work-dir <path>] [--condaEnvsDir <path>]
```

To evaluate a plasmid, replace `evaluationChromosome.nf` with `evaluatePlasmid.nf`.

The only difference between the chromosome and plasmid evaluation is CheckM is not run for plasmids.

##### Parameter descriptions

**Required parameters:**

- `--illumina1 <path>`: Path to 1st file for **cleaned** paired end Illumina reads. May be fastq or gzipped fastq.
- `--illumina2 <path>`: Path to 2nd file for **cleaned** paired end Illumina reads. May be fastq or gzipped fastq.
- `--longReads <path>`: Path to the **cleaned** long reads. May be fastq or gzipped fastq.
- `--assembly <path>`: Path to fasta file of chromosome/plasmid assembly.
- `--outdir <path>`: Path to the output directory, which is where all output files will be stored.

**Optional parameters:**

- See [execution related parameters](#execution-related-parameters)

##### Detailed description

The following metrics about the given chromosome/plasmid assembly are recorded:

- Coverage of short reads and long reads (`shortReadsCoverage`, `longReadsCoverage`)
- Annotations using Prokka (`prokkaAnnotate`)
- Statistics about the assembly itself eg. length, number of contigs, N50 using Quast (`quastEvaluate`)
- (chromosome only) Completeness and contamination using CheckM (`checkmEvaluate`)

The most important statistics are also saved to a single summary document, `chromosome-summary.json` or `plasmid-summary.json` (`makeChromosomeSummary`, `makePlasmidSummary`):

```
{
    "avg short reads coverage": {
        "contig_1": 171.4541
    },
    "short reads mapped": {
        "total reads": 9397580,
        "mapped reads:": 6462316,
        "proportion mapped": 0.6876574607505337
    },
    "avg long reads coverage": {
        "contig_1": 109.7114
    },
    "long reads mapped": {
        "total reads": 102343,
        "mapped reads:": 73142,
        "proportion mapped": 0.7146751609782789
    },
    "assembly length": 4642116,
    "contigs": 1,
    "N50": 4642116,
    "GC": 0.741,
    "CDS": 4071,
    "size by contig": {
        "contig_1": 4642116
    },
    "errors": [],
    "completeness": 1.0,
    "contamination": 0.015414258188824661
}
```

`errors` is a list of any errors encountered while creating the summary document (does not include the errors when running quast/prokka etc.)

##### Output files

- `<outdir>/`
  - `short_read_coverage/`
  - `long_read_coverage/`
  - `prokka/`
  - `quast/`
  - `checkm/` (for chromosome evaluation only)
  - `chromsome-summary.json` / `plasmid-summary.json`

As with other scripts, `nextflow.command.sh`, `nextflow.command.log` and `nextflow.exitcode` is saved for every process (read more [here](#other-outputs)) and the [Nextflow reports](#nextflow-reports) will be generated if `evaluateChromosome.nf`/`evaluatePlasmid.nf` was run alone.

#### Check dependencies: `checkDependencies.nf`

``` sh
nextflow run checkDependencies.nf [-work-dir <path>] [--condaEnvsDir <path>]
```

See [execution related parameters](#execution-related-parameters) for descriptions of the optional parameters.

### Execution related parameters

These are parameters which do not control what is run, only how the pipeline is run (eg. number of threads). For the most part, they are common to all scripts.

- `--threads <num threads>`: Number of threads a process should use. Default: 1.
- `--condaEnvsDir <path>`: Path to where the conda environments are stored. Default: `~/.conda/envs/`.
- `-work-dir`: (note that there is only a single `-` at the front) Path to Nextflow's working directory, which is where temporary files will be stored. Default: `./work/`.
    - Read more [here](#work-dir)

## Troubleshooting

### Help! Where are all my files inside the `work/` directory?

TODO

### Resuming a script

If a script terminated due to an error, you can resume the script from where it left off after fixing the error using the `-resume` flag. Note that you must **not** have deleted the working directory in order for resume to work.

Read more at:
- https://www.nextflow.io/docs/latest/getstarted.html?highlight=resume#modify-and-resume
- https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html
- https://www.nextflow.io/blog/2019/troubleshooting-nextflow-resume.html

### Finding plasmids

If there are multiple contigs in the assembly, it could mean that either mena that there are plasmid(s) in addition to the chromosome, and/or there are multiple chromosomes. The program [Platon](https://github.com/oschwengers/platon)  can take in the final assembly as input, and output likely plasmids. We have included Platon in the `ap-main` environment, so to use it, run `conda activate ap-main` first. See [Platon's documentation](https://github.com/oschwengers/platon) for instructions on how to use it and how to interpret the results. Bear in mind that Platon may miss out plasmids sometimes.

## Tips 

### Periodically clear Nextflow's work directory

The work directory is where Nextflow stores temporary files. By default, this will be `work/`, inside the directory where Nextflow was launched from, however, it may be changed via the `-work-dir` flag.

Some of the temporary files from this pipeline are rather large (eg. reads), so it is a good idea to periodically delete the contents of the work directory to save space on your computer/compute cluster. You should clear the work directory only when **there are no pipelines running** and you **do not intend to resume a previously run pipeline**.

### Running the pipeline on clusters supporting OpenMPI

`mpirun --pernode nextflow run main.nf -with-mpi [pipeline parameters]`

You can nextflow via the `mpirun` command to take advantage of the [OpenMPI standard](https://www.open-mpi.org/) for improved performance.
You can read more about how Nextflow uses OpenMPI [here](https://www.nextflow.io/docs/latest/ignite.html?highlight=mpi#execution-with-mpi) and [here](https://www.nextflow.io/blog/2015/mpi-like-execution-with-nextflow.html).


### Viewing the assembly graph

Flye generates an assembly graph which can be visualised using programs such as [Bandage](https://github.com/rrwick/Bandage) or [AGB](https://github.com/almiheenko/AGB). The graph file is saved at `<outdir>/assembly/flye/assembly_graph.gfa`.

### Repeated configs 

TODO
(profiles, -c config)
