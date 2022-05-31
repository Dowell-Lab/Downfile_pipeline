#!/usr/bin/env nextflow
/*
========================================================================================
                     Downfile Nextflow Pipeline
========================================================================================
 Started 2022-05-18.
 #### Homepage / Documentation
 https://github.com/Dowell-Lab/Downfile_pipeline
 #### Authors
 Lynn Sanford <lynn.sanford@colorado.edu>
========================================================================================
========================================================================================
Pipeline steps:
    1. SAMtools -- convert CRAM file to BAM and extract millions mapped value
    2. BEDtools -- make normalized and non-normalized bedgraphs
    3. BEDtools -- make normalized bigwigs for genome browser
    4. BEDtools -- make 5' bigwigs for dREG
    5. IGV Tools -- make tdfs for genome browser
*/


def helpMessage() {
    log.info"""
    =======================================
     Downfile Pipeline v${params.version}
    =======================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf -profile slurm --crams '/project/*.cram' --outdir '/project/'
    Required arguments:
         -profile                      Configuration profile to use. <base, slurm>
         --crams                       Directory pattern for CRAM files: /project/*.cram (Required if --bams not specified)
         --bams                        Directory pattern for BAM files: /project/*.bam (Required if --crams not specified)
         --workdir                     Nextflow working directory where all intermediate files are saved.
         
    Input File options:
        --singleEnd                    Specifies that the input files are not paired reads (default is paired-end).
        --r1_five_prime                Specifies that for PE reads, the 5' end is on read 1 (default is read 2)
        
    Save options:
        --outdir                       Specifies where to save the output from the nextflow run.
        --savebg                       Saves normalized bedgraphs.
        --savetfitbg                   Saves non-normalized and multimapped read filtered bedgraphs for Tfit/FStitch
        --savebw                       Saves normalized bigwigs for visualization.
        --savedregbw                   Saves 5' bigwigs for input into dREG.
        --savetdf                      Saves TDFs for visualization.
        --saveall                      Saves all output files.
        
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.bedGraphToBigWig = "$baseDir/bin/bedGraphToBigWig"
params.rcc = "$baseDir/bin/rcc.py"
params.workdir = "./nextflowTemp"
output_docs = file("$baseDir/docs/output.md")

import java.text.SimpleDateFormat
def date = new java.util.Date()
def sdf = new SimpleDateFormat("yyMMdd")
output_date =  sdf.format(date)

// Validate inputs

if ( params.genome ){
    genome = file(params.genome)
    if( !genome.exists() ) exit 1, "Genome directory not found: ${params.genome}"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

software_versions = Channel.create()


// Header log info
log.info """========================================
Downfile Pipeline v${params.version}"
========================================"""
def summary = [:]
summary['Pipeline Name']    = 'Downfile Pipeline'
summary['Help Message']     = params.help
summary['Pipeline Version'] = params.version
summary['Run Name']         = custom_runName ?: workflow.runName
if(params.crams) summary['CRAMs']            = params.crams
if(params.bams) summary['BAMs']              = params.bams
summary['Genome Ref']       = params.genome
summary['Chromosome sizes'] = params.chrom_sizes
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Save bedgraphs']   = params.savebg ? 'YES' : params.saveall ? 'YES' : 'NO'
summary['Save tfit bedgraphs']               = params.savetfitbg ? 'YES' : params.saveall ? 'YES' : 'NO'
summary['Save vis bigwigs'] = params.savebw ? 'YES' : params.saveall ? 'YES' :'NO'
summary['Save dreg bigwigs']                 = params.savedregbw ? 'YES' : params.saveall ? 'YES' :'NO'
summary['Save tdfs']        = params.savetdf ? 'YES' : params.saveall ? 'YES' : 'NO'
summary['Max Memory']       = params.max_memory
summary['Max CPUs']         = params.max_cpus
summary['Max Time']         = params.max_time
summary['Output dir']       = params.outdir
summary['Working dir']      = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']     = "$HOME"
summary['Current user']     = "$USER"
summary['Current path']     = "$PWD"
summary['Output dir']       = params.outdir
summary['Script dir']       = workflow.projectDir
summary['Config Profile']   = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "======================================================="

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

/*
 * Parse software version numbers
 */
process get_software_versions {
    time '1h'

    output:
    stdout into software_versions

    script:
    """
    printf "downfile_version: %s\n" ${params.version}
    printf "nextflow_version: %s\n" ${workflow.nextflow.version}
    printf "samtools_version: %s\n" \$(samtools --version | head -1 | awk '{print \$NF}')
    printf "java_version: %s\n" \$(java -version 2>&1 | head -1 | awk -F '"' '{print \$2}')
    printf "bedtools_version: %s\n" \$(bedtools --version | head -1 | awk -F " v" '{print \$2}')
    printf "igvtools_version: %s\n" \$(igvtools version | head -1 | awk '{print \$3}')
    printf "gcc_version: %s\n" \$(gcc --version | head -1 | awk '{print \$NF}')
    printf "python_version: %s\n" \$(python3 --version | head -1 | awk '{print \$NF}')
    printf "pipeline_hash: %s\n" ${workflow.scriptId}
    """
}

software_versions.collectFile(name: "software_versions_downfile_${output_date}_${workflow.runName}.yaml", storeDir: "${params.outdir}/pipeline_info")

/*
 * STEP 1 - Make input file channels, convert to BAM format, sort and collect mapped reads values
 */

if (params.crams) {
    println "[Log 1]: Converting CRAM files to BAM files"
    println "[Log 1]: Genome file being used ..... $params.genome "
    println "[Log 1]: Cram file directory ........ $params.crams"
    println "[Log 1]: Working directory ... $params.workdir"
    println "[Log 1]: Output directory ... $params.outdir"
    cramfiles = Channel
                 .fromPath(params.crams)
                 .map { file -> tuple((file.simpleName), file)}

    process cram_to_bam {
        cpus 16
        queue 'short'
        memory '5 GB'
        time '1h30m'
        tag "$prefix"

        input:
        tuple val(prefix), file(cram) from cramfiles

        output:
        tuple val(prefix), file("${prefix}.bam") into bamfiles

        script:
        """
        samtools view -@ 16 -b -1 -T ${params.genome} ${cram} > ${prefix}.bam
        """
    }
} else {
    bamfiles = Channel
                .fromPath(params.bams)
                .map { file -> tuple(file.simpleName, file) }
}

process samtools {
    println "[Log 1]: Sorting BAM files"
    println "[Log 1]: Collecting mapstats"
    tag "$prefix"
    memory '60 GB'
    cpus 16
    publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
             if (filename.indexOf("flagstat") > 0)             "mapstats/$filename"
        else if (filename.indexOf("millionsmapped") > 0)       "mapstats/$filename"
        else null            
    }

    input:
    tuple val(prefix), file(bamfile) from bamfiles

    output:
    tuple val(prefix), file("${prefix}.downfilepipe.sorted.bam"), file("${prefix}.millionsmapped") into bams_for_bedgraph
    tuple val(prefix), file("${prefix}.downfilepipe.sorted.bam") into bams_for_tfit_conv, bams_for_dreg
    tuple val(prefix), file("${prefix}.flagstat") into flagstats

    script:
    """
    samtools sort -@ 16 ${bamfile} > ${prefix}.downfilepipe.sorted.bam
    samtools flagstat ${prefix}.downfilepipe.sorted.bam > ${prefix}.flagstat
    samtools view -@ 16 -F 0x904 -c ${prefix}.downfilepipe.sorted.bam > ${prefix}.millionsmapped
    """
}

process bam_conversion_tfit {
    println "[Log 1]: Filtering multimapped reads for Tfit"
    cpus 16
    queue 'short'
    memory '5 GB'
    time '2h'
    tag "$prefix"

    when:
    params.savetfitbg || params.saveall

    input:
    tuple val(prefix), file(bam) from bams_for_tfit_conv

    output:
    tuple val(prefix), file("${prefix}.mmfilt.sorted.bam") into bams_for_tfit

    script:
    """
    samtools view -@ 16 -h -q 1 ${bam} | \
        grep -P '(NH:i:1|^@)' | \
        samtools view -h -b > ${prefix}.mmfilt.sorted.bam
    """
}

/*
 *STEP 2 - Create non-normalized and normalized bedGraphs
 */

process tfit_bedgraphs {
    println "[Log 2]: Generating FStitch/Tfit bedgraphs"
    println "[Log 2]: Genome information ..... $params.genome "
    println "[Log 2]: Chromosome Sizes ....... $params.chrom_sizes"
    tag "$prefix"
    queue 'short'
    memory '40 GB'
    time '4h'

    publishDir "${params.outdir}/mapped/bedgraphs_tfit", mode: 'copy', pattern: "${prefix}.bedGraph"
    publishDir "${params.outdir}/mapped/bedgraphs_fstitch", mode: 'copy', pattern: "${prefix}.pos.bedGraph"
    publishDir "${params.outdir}/mapped/bedgraphs_fstitch", mode: 'copy', pattern: "${prefix}.neg.bedGraph"

    when:
    params.saveall || params.savetfitbg

    input:
    tuple val(prefix), file(bamfile) from bams_for_tfit

    output:
    tuple val(prefix), file("${prefix}.bedGraph"), file("${prefix}.pos.bedGraph"), file("${prefix}.neg.bedGraph") into non_normalized_bedgraphs

    script:
    if (params.singleEnd) {
    """    
    genomeCoverageBed \
        -bg \
        -strand + \
        -g ${params.chrom_sizes} \
        -ibam ${bamfile} \
        > ${prefix}.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -strand - \
        -g ${params.chrom_sizes} \
        -ibam ${bamfile} \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        > ${prefix}.neg.bedGraph

    cat ${prefix}.pos.bedGraph \
        ${prefix}.neg.bedGraph \
        > ${prefix}.unsorted.bedGraph
        
    sortBed \
        -i ${prefix}.unsorted.bedGraph \
        > ${prefix}.bedGraph
    """
    } else {
    """   
    samtools view \
        -h -b -f 0x0040 \
        ${bamfile} \
        > ${prefix}.first_pair.bam
        
    samtools view \
        -h -b -f 0x0080 \
        ${bamfile} \
        > ${prefix}.second_pair.bam
        
    genomeCoverageBed \
        -bg \
        -split \
        -strand - \
        -g ${params.chrom_sizes} \
        -ibam ${prefix}.first_pair.bam \
        | sortBed \
        > ${prefix}.first_pair.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -split \
        -strand + \
        -g ${params.chrom_sizes} \
        -ibam ${prefix}.first_pair.bam \
        | sortBed \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        > ${prefix}.first_pair.neg.bedGraph
                     
    genomeCoverageBed \
        -bg \
        -split \
        -strand + \
        -g ${params.chrom_sizes} \
        -ibam ${prefix}.second_pair.bam \
        | sortBed \
        > ${prefix}.second_pair.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -split \
        -strand - \
        -g ${params.chrom_sizes} \
        -ibam ${prefix}.second_pair.bam \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        | sortBed \
        > ${prefix}.second_pair.neg.bedGraph
                     
    unionBedGraphs \
        -i ${prefix}.first_pair.pos.bedGraph ${prefix}.second_pair.pos.bedGraph \
        | awk -F '\t' {'print \$1"\t"\$2"\t"\$3"\t"(\$4+\$5)'} \
        > ${prefix}.pos.bedGraph 
        
    unionBedGraphs \
        -i ${prefix}.first_pair.neg.bedGraph ${prefix}.second_pair.neg.bedGraph \
        | awk -F '\t' {'print \$1"\t"\$2"\t"\$3"\t"(\$4+\$5)'} \
        > ${prefix}.neg.bedGraph 
    
    cat ${prefix}.pos.bedGraph \
        ${prefix}.neg.bedGraph \
        > ${prefix}.unsorted.bedGraph
        
    sortBed \
        -i ${prefix}.unsorted.bedGraph \
        > ${prefix}.bedGraph
    """
    }
}

process norm_bedgraphs {
    println "[Log 2]: Generating normalized bedgraphs"
    println "[Log 2]: Genome information ..... $params.genome "
    println "[Log 2]: Chromosome Sizes ....... $params.chrom_sizes"
    tag "$prefix"
    queue 'short'
    memory '40 GB'
    time '4h'

    publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
        if ((params.saveall || params.savebg) && (filename == "${prefix}.rcc.bedGraph"))    "mapped/bedgraphs_norm/$filename"
        else null
    }

//    publishDir path: {params.saveall ? "${params.outdir}/mapped/bedgraphs_norm" : null}, mode: 'copy', pattern: "${prefix}.rcc.bedGraph"
//    publishDir path: {params.savebg ? "${params.outdir}/mapped/bedgraphs_norm" : null}, mode: 'copy', pattern: "${prefix}.rcc.bedGraph"

    when:
    params.saveall || params.savebg || params.savebw || params.savetdf

    input:
    tuple val(prefix), file(bamfile), file(millionsmapped) from bams_for_bedgraph

    output:
    tuple val(prefix), file("${prefix}.rcc.bedGraph") into bedgraph_tdf
    tuple val(prefix), file("${prefix}.pos.rcc.bedGraph"), file("${prefix}.neg.rcc.bedGraph") into bedgraph_bigwig

    script:
    if (params.singleEnd) {
    """
    genomeCoverageBed \
        -bg \
        -strand + \
        -g ${params.chrom_sizes} \
        -ibam ${bamfile} \
        > ${prefix}.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -strand - \
        -g ${params.chrom_sizes} \
        -ibam ${bamfile} \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        > ${prefix}.neg.bedGraph

    cat ${prefix}.pos.bedGraph \
        ${prefix}.neg.bedGraph \
        > ${prefix}.unsorted.bedGraph

    sortBed \
        -i ${prefix}.unsorted.bedGraph \
        > ${prefix}.bedGraph

    python ${params.rcc} \
        ${prefix}.bedGraph \
        ${millionsmapped} \
        ${prefix}.rcc.bedGraph

    python ${params.rcc} \
        ${prefix}.pos.bedGraph \
        ${millionsmapped} \
        ${prefix}.unsorted.pos.rcc.bedGraph
    sortBed -i ${prefix}.unsorted.pos.rcc.bedGraph > ${prefix}.pos.rcc.bedGraph

    python ${params.rcc} \
        ${prefix}.neg.bedGraph \
        ${millionsmapped} \
        ${prefix}.unsorted.neg.rcc.bedGraph
    sortBed -i ${prefix}.unsorted.neg.rcc.bedGraph > ${prefix}.neg.rcc.bedGraph
    """
    } else {
    """
    samtools view \
        -h -b -f 0x0040 \
        ${bamfile} \
        > ${prefix}.first_pair.bam

    samtools view \
        -h -b -f 0x0080 \
        ${bamfile} \
        > ${prefix}.second_pair.bam

    genomeCoverageBed \
        -bg \
        -split \
        -strand - \
        -g ${params.chrom_sizes} \
        -ibam ${prefix}.first_pair.bam \
        | sortBed \
        > ${prefix}.first_pair.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -split \
        -strand + \
        -g ${params.chrom_sizes} \
        -ibam ${prefix}.first_pair.bam \
        | sortBed \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        > ${prefix}.first_pair.neg.bedGraph

    genomeCoverageBed \
        -bg \
        -split \
        -strand + \
        -g ${params.chrom_sizes} \
        -ibam ${prefix}.second_pair.bam \
        | sortBed \
        > ${prefix}.second_pair.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -split \
        -strand - \
        -g ${params.chrom_sizes} \
        -ibam ${prefix}.second_pair.bam \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        | sortBed \
        > ${prefix}.second_pair.neg.bedGraph

    unionBedGraphs \
        -i ${prefix}.first_pair.pos.bedGraph ${prefix}.second_pair.pos.bedGraph \
        | awk -F '\t' {'print \$1"\t"\$2"\t"\$3"\t"(\$4+\$5)'} \
        > ${prefix}.pos.bedGraph

    unionBedGraphs \
        -i ${prefix}.first_pair.neg.bedGraph ${prefix}.second_pair.neg.bedGraph \
        | awk -F '\t' {'print \$1"\t"\$2"\t"\$3"\t"(\$4+\$5)'} \
        > ${prefix}.neg.bedGraph

    cat ${prefix}.pos.bedGraph \
        ${prefix}.neg.bedGraph \
        > ${prefix}.unsorted.bedGraph

    sortBed \
        -i ${prefix}.unsorted.bedGraph \
        > ${prefix}.bedGraph

    python ${params.rcc} \
        ${prefix}.bedGraph \
        ${millionsmapped} \
        ${prefix}.rcc.bedGraph

    python ${params.rcc} \
        ${prefix}.pos.bedGraph \
        ${millionsmapped} \
        ${prefix}.unsorted.pos.rcc.bedGraph
    sortBed -i ${prefix}.unsorted.pos.rcc.bedGraph > ${prefix}.pos.rcc.bedGraph

    python ${params.rcc} \
        ${prefix}.neg.bedGraph \
        ${millionsmapped} \
        ${prefix}.unsorted.neg.rcc.bedGraph
    sortBed -i ${prefix}.unsorted.neg.rcc.bedGraph > ${prefix}.neg.rcc.bedGraph
    """
    }
}

/*
 *STEP 3 - Generate normalized bigWigs for visualization
 */

process normalized_bigwigs {
    println "[Log 3]: Generating normalized bigwig files for visualization"
    tag "$prefix"
    memory '10 GB'
    queue 'short'

    publishDir "${params.outdir}/mapped/bigwigs_norm", mode: 'copy', pattern: "*.bw"

    when:
    params.saveall || params.savebw

    input:
    tuple val(prefix), file(pos_bedgraph), file(neg_bedgraph) from bedgraph_bigwig

    output:
    tuple val(prefix), file("*pos.rcc.bw"), file("*neg.rcc.bw") into normalized_bigwig

    script:
    """
    ${params.bedGraphToBigWig} ${pos_bedgraph} ${params.chrom_sizes} ${prefix}.pos.rcc.bw
    ${params.bedGraphToBigWig} ${neg_bedgraph} ${params.chrom_sizes} ${prefix}.neg.rcc.bw
    """
}


/*
 *STEP 4 - dREG prep : generate bigwigs for input into dREG
 */

process dreg_prep {
    println "[Log 4]: Generating bigwig files for dREG"

    tag "$prefix"
    memory '60 GB'
    cpus 16
    queue 'short'

    publishDir "${params.outdir}/mapped/bigwigs_dreg", mode: 'copy', pattern: "*.bw"
    publishDir "${params.outdir}/mapped/bedgraphs_dreg", mode: 'copy', pattern: "${prefix}.bedGraph"

    when:
    params.saveall || params.savedregbw

    input:
    tuple val(prefix), file(bamfile) from bams_for_dreg

    output:
    tuple val(prefix), file("${prefix}.pos.bw"), file("${prefix}.neg.bw") into dreg_bigwig
    tuple val(prefix), file("${prefix}.bedGraph") into dreg_bg

    script:
    if (params.singleEnd) {
        """
        bamToBed -i ${bamfile} | awk 'BEGIN{OFS="\t"} (\$5 > 0){print \$0}' | \
        awk 'BEGIN{OFS="\t"} (\$6 == "+") {print \$1,\$2,\$2+1,\$4,\$5,\$6}; (\$6 == "-") {print \$1, \$3-1,\$3,\$4,\$5,\$6}' \
        > ${prefix}.dreg.bed
        sortBed -i ${prefix}.dreg.bed > ${prefix}.dreg.sort.bed

        bedtools genomecov \
                -bg \
                -i ${prefix}.dreg.sort.bed \
                -g ${params.chrom_sizes} \
                -strand + \
                > ${prefix}.pos.bedGraph

        sortBed \
                -i ${prefix}.pos.bedGraph \
                > ${prefix}.pos.sort.bedGraph

        bedtools genomecov \
                -bg \
                -i ${prefix}.dreg.sort.bed \
                -g ${params.chrom_sizes} \
                -strand - \
                | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' > ${prefix}.neg.bedGraph

        sortBed \
                -i ${prefix}.neg.bedGraph \
                > ${prefix}.neg.sort.bedGraph

        ${params.bedGraphToBigWig} ${prefix}.pos.sort.bedGraph ${params.chrom_sizes} ${prefix}.pos.bw
        ${params.bedGraphToBigWig} ${prefix}.neg.sort.bedGraph ${params.chrom_sizes} ${prefix}.neg.bw

        cat ${prefix}.pos.bedGraph \
        ${prefix}.neg.bedGraph \
        > ${prefix}.unsorted.bedGraph

        sortBed \
        -i ${prefix}.unsorted.bedGraph \
        > ${prefix}.bedGraph
        """
    } else {
        if (params.r1_five_prime) {
            """
            samtools view -@ 16 -bf 0x2 ${bamfile} | samtools sort -n -@ 16 \
            > ${prefix}.dreg.bam

            bedtools bamtobed -bedpe -mate1 -i ${prefix}.dreg.bam \
              | awk 'BEGIN{OFS="\t"} (\$9 == "+") {print \$1,\$2,\$2+1,\$7,\$8,\$9}; (\$9 == "-") {print \$1,\$3-1,\$3,\$7,\$8,\$9}' \
              | sort -k 1,1 -k 2,2n > ${prefix}.dreg.sort.bed

            bedtools genomecov -bg \
              -i ${prefix}.dreg.sort.bed \
              -g ${params.chrom_sizes} \
              -strand + \
              > ${prefix}.pos.bedGraph

            bedtools genomecov -bg \
              -i ${prefix}.dreg.sort.bed \
              -g ${params.chrom_sizes} \
              -strand - \
              > ${prefix}.neg.noinv.bedGraph

            cat ${prefix}.neg.noinv.bedGraph \
              | awk 'BEGIN{OFS="\t"} {print \$1,\$2,\$3,-1*\$4}' \
              > ${prefix}.neg.bedGraph

            ${params.bedGraphToBigWig} ${prefix}.pos.bedGraph \
              ${params.chrom_sizes} ${prefix}.pos.bw

            ${params.bedGraphToBigWig} ${prefix}.neg.bedGraph \
              ${params.chrom_sizes} ${prefix}.neg.bw

            cat ${prefix}.pos.bedGraph \
            ${prefix}.neg.bedGraph \
            > ${prefix}.unsorted.bedGraph

            sortBed \
            -i ${prefix}.unsorted.bedGraph \
            > ${prefix}.bedGraph
            """
        } else {
            """
            samtools view -@ 16 -bf 0x2 ${bamfile} | samtools sort -n -@ 16 \
              > ${prefix}.dreg.bam

            bedtools bamtobed -bedpe -mate1 -i ${prefix}.dreg.bam \
              | awk 'BEGIN{OFS="\t"} (\$10 == "+") {print \$1,\$5,\$5+1,\$7,\$8,\$10}; (\$10 == "-") {print \$1,\$6-1,\$6,\$7,\$8,\$10}' \
              | sort -k 1,1 -k 2,2n > ${prefix}.dreg.sort.bed

            bedtools genomecov -bg \
              -i ${prefix}.dreg.sort.bed \
              -g ${params.chrom_sizes} \
              -strand + \
              > ${prefix}.pos.bedGraph

            bedtools genomecov -bg \
              -i ${prefix}.dreg.sort.bed \
              -g ${params.chrom_sizes} \
              -strand - \
              > ${prefix}.neg.noinv.bedGraph

            cat ${prefix}.neg.noinv.bedGraph \
              | awk 'BEGIN{OFS="\t"} {print \$1,\$2,\$3,-1*\$4}' \
              > ${prefix}.neg.bedGraph

            ${params.bedGraphToBigWig} ${prefix}.pos.bedGraph \
              ${params.chrom_sizes} ${prefix}.pos.bw

            ${params.bedGraphToBigWig} ${prefix}.neg.bedGraph \
              ${params.chrom_sizes} ${prefix}.neg.bw

            cat ${prefix}.pos.bedGraph \
            ${prefix}.neg.bedGraph \
            > ${prefix}.unsorted.bedGraph

            sortBed \
            -i ${prefix}.unsorted.bedGraph \
            > ${prefix}.bedGraph
            """
        }
    }
}

/*
 *STEP 5 - IGV Tools : generate tdfs for optimal visualization in Integrative Genomics Viewer (IGV)
 */

process igvtools {
    println "Log[5]: Generating TDFs"
    tag "$prefix"
    memory '30 GB'
    time '1h'
    queue 'short'
    errorStrategy 'ignore'

    // This often blows up due to a ceiling in memory usage, so we can ignore
    // and re-run later as it's non-essential.

    publishDir "${params.outdir}/mapped/tdfs", mode: 'copy', pattern: "*.tdf"

    when:
    params.saveall | params.savetdf

    input:
    tuple val(prefix), file(bg) from bedgraph_tdf

    output:
    tuple val(prefix), file("*.tdf") into tiled_data_ch

    script:
    """
    igvtools toTDF ${bg} ${prefix}.rcc.tdf ${params.chrom_sizes}
    """
 }


/*
 * Completion report
 */
workflow.onComplete {

    def report_fields = [:]
    report_fields['version'] = params.version
    report_fields['runName'] = custom_runName ?: workflow.runName
    report_fields['success'] = workflow.success
    report_fields['dateComplete'] = workflow.complete
    report_fields['duration'] = workflow.duration
    report_fields['exitStatus'] = workflow.exitStatus
    report_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    report_fields['errorReport'] = (workflow.errorReport ?: 'None')
    report_fields['commandLine'] = workflow.commandLine
    report_fields['projectDir'] = workflow.projectDir
    report_fields['summary'] = summary
    report_fields['summary']['Date Started'] = workflow.start
    report_fields['summary']['Date Completed'] = workflow.complete
    report_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    report_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    report_fields['summary']['Pipeline repository Git URL'] = workflow.repository ?: 'Not stored'
    report_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId ?: 'See hash'
    report_fields['summary']['Pipeline Git branch/tag'] = workflow.revision ?: 'See hash'
    report_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    report_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    report_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/report_template.txt")
    def txt_template = engine.createTemplate(tf).make(report_fields)
    def report_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/report_template.html")
    def html_template = engine.createTemplate(hf).make(report_fields)
    def report_html = html_template.toString()

    // Write summary HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report_downfile_${workflow.runName}.html" )
    output_hf.withWriter { w -> w << report_html }
    def output_tf = new File( output_d, "pipeline_report_downfile_${workflow.runName}.txt" )
    output_tf.withWriter { w -> w << report_txt }

    log.info "Downfile Pipeline Complete"

}
