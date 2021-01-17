# Genome assembly pipeline

This is a hybrid genome assembly pipeline for bacterial genomes written in [Nextflow](https://www.nextflow.io/).

## Quickstart

### Installation

1. Create the conda environments from yaml files in the `conda-envs` directory
   - Use only one of `assembly.yml` and `assembly-nscc.yml`. The difference is that `assembly.yml` includes samtools and bwa while `assembly-nscc.yml` does not.
   - Before creating `assembly.yml`/`assembly-nscc.yml`, change the following variables:
        - Set `PILONJAR` to the path to the Pilon's jar file. It is usually at `<path to conda directory>/pkgs/pilon-1.23-2/share/pilon-1.23-2/pilon-1.23.jar`.
        - Set `PLATONDB` to the path to Platon's database, which will be donwloaded in the next step.
1. Download databases needed by various tools
   - Platon: download the database at this [link](https://zenodo.org/record/4066768/files/db.tar.gz?download=1) and extract it to `PLATONDB`
1. Install [Nextflow](https://www.nextflow.io) and make sure the `nextflow` command is available in your path.


### Running the pipeline

The several different operations which can be run, such as genome assembly, assembly evaluation, read quality assessment and dependency checking. 
Each operation is located in a different `.nf` script in the root directory. 

A script is launched using the `nextflow run` command from the root directory, eg. `nextflow run assemble.nf`.
Various command line arguments can be passed to a script (eg. `--illumina1`, `-work-dir`). 
See the [scripts section](#scripts) for descriptions of the scripts available and their parameters.

**Things to take note of:**

- Do note that some flags begin with a single dash `-`, while some begin with a double dash `--`.
- The scripts need to be launched from the root directory of this repository, in order for the configurations in the `nextflow.config` file to be applied.

**(TODO) Nextflow's working directory (workDir):**

Nextflow stores temporary data in a working directory.
This will basically includes all of the files mentioned above and possibly additional intermediate files.
This can be safely deleted after the pipeline finishes running.

See more at: https://www.nextflow.io/docs/latest/script.html?highlight=workdir#implicit-variables, https://github.com/danrlu/Nextflow_cheatsheet#the-working-directory

TODO document `-resume`

#### Recommended workflow:

1. Run `assemble.nf` to assemble the genome.
1. Check whether the resulting assembly has multiple contigs. If there are multiple contigs, decide how to separate them (chromosomes and plasmids, possibly belonging to multiple bacteria). It may be helpful to use Platon (included in the `urops-assembly` environment) and/or BLAST.
1. Evaluate chromsomes using `evaluateChromosome.nf`, and plasmids using `evaluatePlasmid.nf`.
1. Delete the work directory if you don't need any of the other files.


#### TIP: for compute clusters with OpenMPI

If you are running on the pipeline on a compute cluster which uses the Open MPI standard (eg. NSCC's Aspire 1 server), you can run the pipeline as if it were an MPI application:

``` sh
mpirun --pernode nextflow run <script name> -with-mpi [pipeline parameters]
```

Eg. `mpirun --pernode nextflow run assemble.nf -with-mpi --illumina1 <path> --illumina2 <path> --pacbio <path> --outdir <path>`

Read more at: https://www.nextflow.io/docs/latest/ignite.html?highlight=mpi#execution-with-mpi , https://www.nextflow.io/blog/2015/mpi-like-execution-with-nextflow.html

#### TIP: Minimal PBSPro job for the NSCC Aspire 1 server

Here is a minimal template for a job script which can be submitted using the `qsub` command:

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

cd assembly-pipeline/ # cd to the directory containing this repository
mpirun --pernode nextflow run assemble.nf -profile nscc -with-mpi [pipeline parameters]
```

## Scripts

_Note: Name of Nextflow process is written in parenthesis. Eg. (`cleanShortReads`)_

### Genome assembly: `assemble.nf`

Assemble a genome from raw Illumina and Pacbio reads. Some basic statistics (number of contigs, size of contigs, and circularity of contigs) are given. More detailed evaluation can be done using [`evaluateChromosome.nf` and `evaluatePlasmid.nf`](#evaluation).

``` sh
nextflow run assemble.nf --illumina1 <path> --illumina2 <path> --pacbio <path> --outdir <path> \
    [--shortReadsKeepPercent <percent>] [--shortReadsStartTrimq <percent>] [--shortReadsMinTrimq <percent>] \
    [--bbdukArgs <args>] [--pilonMaxIters <number>] [--raconMaxIters <number>] \ 
    [--canuGenomeSize <size>] [--forceCirclator | --noCirclator] \
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
   - (Coming soon) ~You can force Circlator to always run / never run using `--forceCirclator` / `--noCirclator`.~
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
    - `--bbdukArgs <args>`: Arguments other than inputs, outputs and `trimq` to pass to bbduk. Default: `qtrim=rl minlength=40`.
- `--pilonMaxIters <number>`: Maximum number of iterations to run Pilon for. Default: 6.
- `--raconMaxIters <number>`: Maximum number of iterations to run Racon for. Default: 4.
- Circularisation:
    - `--canuGenomeSize <genome size>`: When specified, force Canu to use this genome size. See the Canu documentation for genomeSize for valid values. Otherwise, calculate the genome size from the assembly. Default: not specified.
    - `--forceCirclator`: Force Circlator to be used, regardless of the state of the Flye assembly.
    - `--noCirclator`: Do not use Circlator, regardless of the state of the Flye assembly.
- Also see [execution related parameters](#execution-related-parameters)

<a id="evaluation"></a>
    
### Assembly evaluation: `evaluateChromosome.nf`, `evaluatePlasmid.nf`

To evaluate a chromosome:

``` sh
nextflow run evaluateChromosome.nf --illumina1 <path> --illumina2 <path> --pacbio <path> --chromosome <path> \
    --outdir <path> \
    # execution-related params
    [--skipDepChecks] [--threads <number>] [-work-dir <path>] [--condaEnvsDir <path>] [-profile <profiles>]
```

To evaluate a plasmid:

``` sh
nextflow run evaluatePlasmid.nf --illumina1 <path> --illumina2 <path> --pacbio <path> --plasmid <path> \
    --outdir <path> \
    # execution-related params
    [--skipDepChecks] [--threads <number>] [-work-dir <path>] [--condaEnvsDir <path>] [-profile <profiles>]
```


The only difference between the chromosome and plasmid evaluation is CheckM is not run for plasmids.


#### Detailed description

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
  - `checkm/` (for chromsome evaluation only)
  - `chromsome-summary.json` / `plasmid-summary.json`
  
As with other scripts, `nextflow.command.sh`, `nextflow.command.log` and `nextflow.exitcode` is saved for every process. See [here](#other-outputs).

#### Parameter descriptions

**Required parameters:**

- `--illumina1 <path>`: Path to 1st file for **cleaned** paired end Illumina reads. May be fastq or gzipped fastq.
- `--illumina2 <path>`: Path to 2nd file for **cleaned** paired end Illumina reads. May be fastq or gzipped fastq.
- `--pacbio <path>`: Path to the **cleaned** Pacbio reads. May be fastq or gzipped fastq.
- `--plasmid <path>`/`--chromosome <path>`: Path to fasta file of chromosome/plasmid assembly.
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

These are parameters which do not control what is run, only how the pipeline is run (eg. execution speed).

- `--threads <num threads>`: Number of threads a process should use. Default: 1.
- `--condaEnvsDir <path>`: Path to where the conda environments are stored. Default: `~/.conda/envs/`.
- `-work-dir`: (note that there is only a single `-` at the front) Path to Nextflow's working directory, which is where temporary files will be stored. Default: `./work/`.
    - When running on the NSCC Aspire 1 server, you can set this to a folder in the scratch directory to save space in your home directory.
- `--skipDepChecks`: Skip pre-pipeline dependency checks. Does not apply to the `checkAllDependencies.nf` script.
- `-profile <profiles>`: (note that there is only a single `-` at the front) List of [Nextflow profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles) to apply. Options: 
    - `local`: Applies settings for working on a laptop: `--threads 4`.
    - `nscc`: Applies settings for working on the NSCC Aspire 1 sever: `-work-dir ~/scratch/work/ --threads 20`.

