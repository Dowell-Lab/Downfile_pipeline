# Downfile_pipeline
Nextflow pipeline for making bedgraphs, bigwigs, and tdfs from CRAMs/BAMs

This pipeline is intended to make mapped read files of different formats for visualization or downstream bidirectional calling purposes.

# Installation

 `$ git clone https://github.com/Dowell-Lab/Downfile_pipeline`

# Requirements

Check each tool for configuration requirements.

- Nextflow (version >= 19.10.0)

- SAMtools

- BEDtools

- Python3

- IGVtools

## Configuration Files

Each run of the Downfile_pipeline requires a configuration file which specifies genome-specific variables for CRAM/BAM and bedgraph manipulation. Edit the configuration file to reflect the correct cluster paths for genome and chromosome size files.

## Loading requirements on SLURM

  Below is a summary of all FIJI modules needed to run Downfile_pipeline.

  ```
  module load samtools/1.8
  module load bedtools/2.28.0
  module load igvtools/2.3.75
  module load python/3.6.3
  ```

# Running Downfile_pipeline

## Usage

   The typical command for running the pipeline is as follows:

   `$ nextflow run main.nf -profile hg38 --crams '/project/*.sorted.cram' --workdir '/project/tempfiles' --outdir '/project/' --saveall`

   Below are the arguments:

   ```
    Required arguments:
         -profile                      Configuration profile to use. <genome_user>
         --crams                       Directory pattern for cram files: /project/*.sorted.cram (Required if --bams not specified).
         --bams                        Directory pattern for bam files: /project/*.sorted.bam (Required if --crams not specified).
         --workdir                     Nextflow working directory where all intermediate files are saved.

    Input File options:
        --singleEnd                    Specifies that the input files are not paired reads (default is paired-end).
        --r1_five_prime                (for dREG bigwig prep only) If input file is paired, specifies if read 1 has the 5 prime end
                                       (default R2 is five prime, must be manually determined)

    Save options:
        --outdir                       Specifies where to save the output from the nextflow run.
        --savebg                       Saves normalized bedgraphs.
        --savetfitbg                   Saves non-normalized and multimapped read filtered bedgraphs for Tfit/FStitch
        --savebw                       Saves normalized bigwigs for visualization.
        --savedregbw                   Saves 5' bigwigs for input into dREG.
        --savetdf                      Saves TDFs for visualization.
        --saveall                      Saves all output files.
   ```

## Credits

* Lynn Sanford ([@lynn-sanford](https://github.com/lynn-sanford))
