** Note, this is forked from [https://github.com/indrops/indrops](https://github.com/indrops/indrops) ** and was modified for our use.

## Project YAML file.

An inDrops project is composed of a series of sequencing runs, each including one (or several) indrops libraries within it. A sequencing run can further be split into several parts (effectively arbitrary chunks) to parallelize the analysis. Give example of a chunk

The project yaml file contains the details of all sequencing runs and libraries within a project. 

The same project can contain runs from different versions of the inDrops platform. 

A project will be aligned against the same reference genome with the same alignment parameters. 

We have provided the yaml files used for this paper above, though you will need to modify them for your own use (i.e. change paths). They are named "Early_med_bud.yaml" and this yaml file should be used for fastq #'s SRR8147022-SRR8147025. yaml file wound_healing_med_bud.yaml should be used with SRR8147026-SRR8147029. Finally, intact_and_contralateral.yaml should be used for SRR8147030-SRR8147033. These fastq's can all be found here: https://www.ncbi.nlm.nih.gov/sra?term=SRP167700. Once you have downloaded the files, it's important to rename each file so that it is compatible with the yaml files. The fastq's should have a structure like (note you don't need to call it Undetermined):

Undetermined_S0_{split_affix}_{read}_001.fastq.gz

You must deteremine which read the fastq file is associated with and then rename files adding the appropriate read information into {read}. Manual inspection of the SRR downloads should allow for confirmation of read. Read structure is as follows:

1. 61bp Read 1             transcript
2. 8bp Index Read 1 (i7)   8bp part 1 single cell barcode
3. 8bp Index Read 2 (i5)   library index
4. 14bp Read 2             part 2 cell barcode (8bp) / UMI (6bp) / Poly T

A table has been provided above further describing where to find each of the biological replicates within the fastq's deposited on SRA. It is important to note that each SRR number does NOT correlate with one sample. For example, SRR8147025 has reads from early and medium bud blastemas (see sample.table.txt above), but only contains ONE read. Therefore, you will need all fastq's (SRR8147022-SRR8147025) to obatin all 4 reads for running the pipeline. When running the yaml file these will be demultiplexed.

## Supported library versions
   - v1 : original design where R2 is the biological read and R1 is the metadata read. 
   - v2 : inversion of v1 where R1 is the biological read and R2 is the metadata read.
   - v3 : summer 2016 redesign requiring manual demultiplexing. R1 is the biological read.
          R2 carries the first half of the gel barcode, R3 carries the library index and R4
          the second half of the gel barcode, the UMI and a fraction of the polyA tail.

## Installation
The package requires
  - Python 2.7 (with the packages numpy, scipy, matplotlib, pysam>0.9.0, pyyaml, pyfasta). [See Appendix 2]
  - RSEM (1.2.16+)
  - Bowtie (1.1.1+)
  - samtools (1.3.1+) [See Appendix 3] *This specific version is needed to account for a BAM-format oddity in RSEM output.
  - Java 
The path to the directories containing these executables should be set in the project YAML.
If these executables can be found in the PATH variables, this project YAML paths can be left empty, or not specified.

### March 7th Notice -- PySAM version.

Previous installation instructions install PySAM version 0.6.0. To install the correct PySAM version, use the following commands:

  conda remove pysam
  conda install pip
  pip install pysam==0.9.1


## Project YAML file

An example YAML file is provided in `test/test_project.yaml`. It should contain the following information:

    project_name : "project_name"
    project_dir : "/path/to/project/dir"  #This dir should be user-owned and writable, all output will go into this dir.
    paths : 
      bowtie_index : "/path/to/index" #This index will be built automatically
      # The paths below can be omitted if the relevant directories are already on $PATH
      bowtie_dir : "/path/to/bowtie/dir/"
      python_dir : "/path/to/env/bins/"
      java_dir: "/path/to/java/dir/"
      rsem_dir: "/path/to/rsem/dir/"
      samtools_dir: "/path/to/samtools-1.3.1/bin/" #This needs to be version 1.3.1, 1.3 is not good enough!

    sequencing_runs : 
      # A list of sequencing runs which form the project. 
      # Each run should have:
      - name : "MyRun" # The name of the run will be used as a prefix in filenames, so keep it sane.
        version : "vN" # Can be 'v1', 'v2' or 'v3'

      # For a run with a single 'part', and a single library
        dir : "/path/to/run_files/"
        fastq_path : "{read}.fastq.gz" # Read with be replaced by R1, R2, R3, R4 as appropriate.
        library_name : "my_library"

        # This will expect to find the files:
        #    /path/to/run_files/R1.fastq.gz (and R2...)

       # For a run with several parts, but a single library
        dir : "/path/to/run_files/"
        fastq_path : "{split_affix}_{read}.fastq.gz" # Read with be replaced by R1, R2, R3, R4 as appropriate.
        split_affixes : ["L001", "L002"]
        library_name : "my_library"

        # This will expect to find the files:
        #    /path/to/run_files/L001_R1.fastq.gz (and R2...)
        #    /path/to/run_files/L002_R1.fastq.gz (and R2...)

       # For a run with several parts, several libraries, that have already been demultiplexed
        dir : "/path/to/run_files/"
        fastq_path : "{library_prefix}_{split_affix}_{read}.fastq.gz" # Read with be replaced by R1, R2, R3, R4 as appropriate.
        split_affixes : ["L001", "L002"]
        libraries : 
          - {library_name: "test_lib1", library_prefix: "lib1"}
          - {library_name: "test_lib2", library_prefix: "lib2"}

        # This will expect to find the files:
        #    /path/to/run_files/lib1_L001_R1.fastq.gz (and R2...)
        #    /path/to/run_files/lib1_L002_R1.fastq.gz (and R2...)
        #    /path/to/run_files/lib2_L001_R1.fastq.gz (and R2...)
        #    /path/to/run_files/lib2_L002_R1.fastq.gz (and R2...)

       # For a V3 run with several parts, with several libraries that are not already demultiplexed
        dir : "/path/to/run_files/"
        fastq_path : "{library_prefix}_{split_affix}_{read}.fastq.gz" # Read with be replaced by R1, R2, R3, R4 as appropriate.
        split_affixes : ["L001", "L002", "L003", "L004"]
        libraries :  # The library index is what the expected index read sequence (on a NextSeq, this is the reverse complement of the index sequence)
          - {library_name: "test_lib3", library_index: "ATAGAG"}
          - {library_name: "test_lib4", library_index: "AGAGGA"}

        # This will expect to find the files:
        #    /path/to/run_files/lib1_L001_R1.fastq.gz (and R2, R3, R4...)
        #    /path/to/run_files/lib1_L002_R1.fastq.gz (and R2, R3, R4...)
        #    /path/to/run_files/lib1_L003_R1.fastq.gz (and R2, R3, R4...)
        #    /path/to/run_files/lib1_L004_R1.fastq.gz (and R2, R3, R4...)

#### Note about v3 runs. 
The raw BCL files are needed for manual demultiplexing. Move the raw BCL files to a run directory, then use the following command to extract the R1,R2,R3 and R4 files.

    cd /run/dir/
    bcl2fastq --use-bases-mask y*,y*,y*,y* --mask-short-adapter-reads 0 --minimum-trimmed-read-length 0
    # The 'dir' used in the project YAML file should then be:
    #   /run/dir/Data/Intensities/BaseCalls/

## Analysis steps

### 0. Generate bowtie index from axolotl transcriptome
The version of the transcriptome we used can be downloaded here:https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE121737&format=file&file=GSE121737%5FAxolotl%2ETrinity%2ECellReports2017%2Efasta%2Egz

To generate a bowtie index with this transcriptome.
   ```
   bowtie-build Axolotl.Trinity.CellReports2017.fasta Axolotl.Trinity.CellReports2017.fasta
   
   ```
   

### 1. Filter
This iterates over sequencing run parts, optionally filtered by a list of sequencing parts, and a list of libraries of interest.

    python indrops.py project.yaml filter [--total-workers 1] [--worker-index 0]
          [-r --runs RUNS ] [-l --libraries LIBRARIES ]

    # --runs comma-separated list of runs :               If specified, step will be restricted to run parts coming #                                                     from runs in the list
    # --libraries comma-separated list of libraries :      If specified, step will be restricted to run parts that 
    #                                                     contain reads from a library in the list
    # 
    # Resulting workload (a list of run parts), will be split among N --total-workers,
    # where worker with --worker-index i will do steps (i, N+i, 2N+i, ...)

This step reads the raw FastQ files as input and filters them:
  - For every raw read, it determines if the read has the expected structure (depending on library version). 
  - For reads with correct structure, it runs Trimmomatic.
  - For reads surviving Trimmomatic, it finds and trims the polyA tail a maximum length of 4, and checks if the reads are still above MIN_LEN.
  - For surviving reads, it determines which fraction of the read is composed of runs of the same base (only considering runs of 5 or more). 
    It rejects reads whose fraction is greater than `low_complexity_filter_arguments:max_low_complexity_fraction`.

As output, for every input run part, this produces a filtered FastQ file for every library contained in that run. These files are referred to as 'parts of libraries'.

A log is created detailing what happened to every input read. An index is created that lists the number of reads found for every barcode. 

### 2. Identify abundant barcodes
This iterates over libraries, optionally filtered by a list. 

    python indrops.py project.yaml identify_abundant_barcodes [--total-workers 1] [--worker-index 0]
          [-l --libraries LIBRARIES]
  
    # --libraries comma-separated list of librares : If specified, step will be restricted to libraries in this list.
    # 
    # Resulting workload (a list of libraries), will be split among N --total-workers,
    # where worker with --worker-index i will do steps (i, N+i, 2N+i, ...)
    # 
    #    *Note* This step is fast, it does not need to be dispatched to several workers.

For each library, this collates the results of filtering all the sequencing run parts that have reads related to this library. It then outputs,
  - Histogram of the distribution barcode abundances
  - Summary table of filtering for that library
  - An index to be used by `sort`. 

### This is where we diverge from the original inDrops processing pipeline!

### 3. Sort reads according to their barcode of origin.
This iterates over parts of libraries, optionally filtered by a list.

    python extract_barcoded_reads.py project.yaml sort [--total-workers 1] [--worker-index 0]
          [-l --libraries LIBRARIES]

    # --libraries comma-separated list of libraries :    If specified, step will be restricted to library-run-parts
    #                                                   that contain reads from a library in the list
    # 
    # Resulting workload (a list of library-run-parts), will be split among N --total-workers,
    # where worker with --worker-index i will do steps (i, N+i, 2N+i, ...)
    #
    #    *Note* this step is currently memory intensive, as it loads the entire 'library-run-part' in memory. 

This sorts the reads according to the name of their barcode of origin. Barcodes with less than 250 total reads (across all library-run-parts) are ignored, and placed at the end of the file.

As output, this creates a gzipped FastQ file and an index of the byte offsets for every barcode with more than 250 reads.

### 4. Quantify expression
This iterates over a list of barcodes, from a list of optionally filtered libraries. 

    python extract_barcoded_reads.py project.yaml quantify --no-bam > barcoded_reads.fastq 
            [--total-workers 1] [--worker-index 0]
            [-l --libraries LIBRARIES] [-r --runs RUNS ]
            [--min-reads 750] [--min-counts 0]
            [--analysis prefix '']
            [--no-bam]

    # --min-reads INT :                                 Ignore barcodes with less than specified number of reads.
    # --min-counts INT :                                Ignore output for barcodes with less than the specified number
    #                                                   of UMIFM counts. This significantly speeds up
    #                                                   downstream processing.
    # --analysis-prefix STR :                           Prefix output data files with the specified prefix.
    #                                                   (filename --> prefix.filename)
    # --no-bam :                                        If specified, do not output and process BAM files. 
    # 
    # --libraries comma-separated list of libraries      If specified, step will be restricted to libraries
    #                                                   in this list.
    # --runs comma-separated list of runs               If specified, only align reads coming from these runs
    #                                                   [This is an uncommon use case.]
    # 
    # 
    # The resulting list of barcodes will be split among --total-workers, with worker identified by --worker-index.
    #    *Note* This step requires ~2Gb of memory. 

This step is resumable. If the same --analysis-prefix/--total-workers/--worker-index was previously running, another run will only quantify barcodes that were not previously quantified (or whose data was lost). To force requantification, delete files in /project_dir/library_dir/quant_dir/[prefix.]worker\*_[total_workers]\*

## note, must make the read name unique:
```
encode_read_number_in_fastq.pl barcoded_reads.fastq > barcoded_reads.adj.fq
```

# align reads to target transcripts using bowtie:
```
bowtie target.fasta.bowtie -q -p 20 -a --best --strata --chunkmbs 1000 --sam -m 200 -n 1 -l 15 -e 100 barcoded_reads.adj.fq  >  target.bowtie.adj.sam
```
# generate count matrix:
note: bam_to_count_matrix.pl can be modified to output more or less cells to the sc.counts.matrix file. This can be done by changing the value for max_top_cells. We expect about 3000 cells per library. 

```
bam_to_count_matrix.pl --bam target.bowtie.adj.sam > sc.counts.matrix
```
# annotate count matrix
Trinity (https://github.com/trinityrnaseq/trinityrnaseq/wiki) has a nice perl script to take care of this for us, so download and install Trinity for these next few steps (we used Trinity v2.5.1). You should run these next two commands from within the Trinity codebase. The script we will use for annotation can be found here: https://github.com/trinityrnaseq/trinityrnaseq/blob/master/Analysis/DifferentialExpression/rename_matrix_feature_identifiers.pl 

To perform the annotation, we need two things 1) the sc.counts.matrix and 2) an annotation mapping file for the transcriptome of interest. The annotation mapping file for the the Axolotl.Trinity.CellReports2017.fasta transcriptome is called Axo.Mar2014.Trinotate.xls.annot_mapping and is available above.

```
rename_matrix_feature_identifiers.pl sc.counts.matrix Axo.Mar2014.Trinotate.xls.annot_mapping > sc.counts.annotated.matrix
```

# create representative gene count matrix
Since the transcriptome we used has a lot of isoforms, we collapsed these down to one representative isoform. Perl script can be found here: https://github.com/trinityrnaseq/trinityrnaseq/blob/master/util/misc/trinity_trans_matrix_to_rep_trans_gene_matrix.pl and again run from the Trinity codebase. 

```
trinity_trans_matrix_to_rep_trans_gene_matrix.pl sc.counts.annotated.matrix > sc.counts.annotated.matrix.repGene
```

This count matrix is now ready to go into an analysis software. We used Seurat (https://satijalab.org/seurat/) for our manuscript! 


## Appendix 1: Parallelizing the analysis steps

### Analyzing only parts of the project.

Most parts of the analysis can be filtered by specifying a list of sequencing runs,
a list of sequencing libraries, or both. When a filter is provided, the analysis will
only be carried out on data matching the filter.

Every part of the analysis can be filtered based on both libraries and sequencing runs.

    # Will filter all parts from runs Run5 and Run6:
    python indrops.py test_project.yaml filter --runs Run5,Run6

    # Will sort all parts from all runs of libraries test_lib3 and test_lib4:
    python indrops.py test_project.yaml sort --libraries test_lib3,test_lib4

### Dividing the analysis between jobs

Most parts of the analysis can easily be divided for concurrent processing in different jobs,
by specifying the total number of jobs (--total-workers) and the index of the current worker (--worker-index). 

    # Submitting the 20 commands below would filter all run parts within the project in 20 different parts.
    python indrops.py test_project.yaml filter --total-workers 20 --worker-index [0-19]

## Appendix 2: Using a custom Python environment on the a cluster

### How to install conda and create a new environment

Download Miniconda (the anaconda package manager, without all the packages)

    mkdir -pv /user_owned/path
    cd  /user_owned/path
    wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh

Install Miniconda

    bash Miniconda-latest-Linux-x86_64.sh
    # Agree to license with “yes”, and choose to install in a directory that is user owned.
    # I installed it in: /groups/klein/adrian/miniconda

Create a new Python environment (in this example, in /groups/klein/adrian/pyndrops)
install Python2.7, Numpy, Scipy, Pandas, Matplotlib, PyYaml, PySAM

    conda create -p /groups/klein/adrian/pyndrops python numpy scipy pandas pyyaml matplotlib pip
    source activate /groups/klein/adrian/pyndrops
    pip install pyfasta pysam==0.9.1

## Appendix 3: Installing Samtools 1.3.1

    mkdir -pv SAMTOOLS_DIR
    cd SAMTOOLS_DIR
    wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
    tar xvfj samtools-1.3.1.tar.bz2
    cd samtools-1.3.1
    make
    make prefix=. install

Now add `SAMTOOLS_DIR/samtools-1.3.1/bin/` as the `samtools_dir` in your project YAML file.
