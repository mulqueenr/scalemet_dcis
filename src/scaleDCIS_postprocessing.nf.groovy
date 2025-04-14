//Nextflow pipeline for processing Navin lab spatial curio+multiome libraries//

// Declare syntax version
nextflow.enable.dsl=2

// Script parameters
params.projDir="/data/rmulqueen/projects/scalebio_dcis" //Project directory for processing
params.runDir="/data/rmulqueen/projects/scalebio_dcis/data/240202_prelim1/scale_data" //Output from scalemethyl pipeline
params.outputPrefix="scale" //Output from scalemethyl pipeline
params.src = "/home/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/" //src directory for cloned repo
params.min_cnv_count= "100000" //minimum readcount to try copykit per cell

params.max_cpus=50
params.max_forks=100


log.info """

		================================================
		             scalebio DCIS Pipeline v1.0
		================================================
		Project Dir : ${params.projDir}
		Run Dir : ${params.runDir}
		Src Dir : ${params.src}
		NF Working Dir : ${workflow.launchDir}

		Max cpus : ${params.max_cpus}
        Max forks : ${params.max_forks}
		================================================

""".stripIndent()


//COPYKIT PROCESSING//
process COUNT_READS { 
    //Count reads to split bams by minimum read counts
	cpus "${params.max_cpus}"
    publishDir "${params.runDir}/cnv/read_counts", mode: 'copy', overwrite:true, pattern: "*tsv"
    cpus "${params.max_cpus}"

	input:
		path runDir
	output:
		path("cells_pf.tsv"), emit: cells_pf
        path("cells_pf.tsv"), emit: uniq_read_counts

    script:
		"""
        #use a function to parallelize
        count_reads() { 
        samtools view \$1 | awk -v b=\$1 '{split(\$1,a,":"); print a[8],b}' | sort | uniq -c | sort -k1,1n
        }
        export -f count_reads

        parallel -j ${task.cpus} count_reads ::: \$(find ${runDir}/alignments -maxdepth 5 -name '*bam') | sort -k1,1n > unique_read_counts.tsv
        awk '\$1>${params.min_cnv_count} {print \$1,\$2}' unique_read_counts.tsv > cells_pf.tsv
		"""
}


process SPLIT_BAMS {
	//Take output of COUNT_READS with each line being a tuple, then split from source bam to single cell bam
    publishDir "${params.runDir}/cnv/sc_bam", mode: 'copy', overwrite:true, pattern: "*bam"

	input:
		tuple val(read_count),val(idx),path(bam)
	output:
		path("*bam")
	script:
	"""
        outprefix=${bam.Name}
        outprefix=\$(echo \$outprefix | sed -e 's/.dedup.bam//g' -)
        ((samtools view -H $bam) && (samtools view $bam | awk -v i=$idx '{split(\$1,a,":"); if(a[8]==i); print \$0}')) | samtools view -bS > ${outprefix}.${idx}.bam
	"""
}


process COPYKIT { 
	//Generate cell level Fastq Files from BCL Files and generated white list
	//TODO This container should be updated to be in the SIF and not local run
	cpus "${params.max_cpus}"
	containerOptions "--bind ${params.src}:/src/,${params.projDir}"
    publishDir "${params.runDir}/cnv/copykit", mode: 'copy', overwrite:true, pattern: "*{tsv,rds}"
    publishDir "${params.runDir}/cnv/copykit_plots", mode: 'copy', overwrite:true, pattern: "*pdf"

	label 'copykit'
	input:
		path bams
	output:
		path("*.tsv"), emit: copykit_tsv
        path("*.rds"), emit: copykit_rds
        path("*.pdf"), emit: cnv_plots
    script:
		"""
        Rscript /src/copykit_cnvcalling.R \
            --input_dir . \
            --output_dir . \
            --output_prefix ${params.outputPrefix} \
            --task_cpus ${task.cpus}
		"""
}

// TRIM, ALIGN, and DEDUPLICATE READS
process AMETHYST_INIT {
	//TRIM READS OF ADAPTERS AND KNOWN METHYLATED REGIONS (GAP FILLS)
	cpus "${params.max_cpus}"
	publishDir "${params.outdir}/reports/adapter_trim", mode: 'copy', overwrite: true, pattern: "*.log"
	containerOptions "--bind ${params.src}:/src/,${params.projDir}"
    
	label 'amethyst'

	input:
		path runDir, stageAs: 'rundir/'
        path copykit_clones

	output:
		tuple val(cellid),path("*.R1_001.trim.fastq.gz"), path("*.R2_001.trim.fastq.gz"), emit: fqs
		path("*.trim_report.log"), emit: trim_log
	script:
		"""
        Rscript /src/amethyst_initial_processing.R \
        --input_dir rundir/. \
        --output_prefix ${params.outputPrefix} \
        --task_cpus ${task.cpus}
		"""
}


workflow {
    
    //Split bams and run copykit
    indir= channel.fromPath( params.runDir , type: 'dir')

    COUNT_READS(indir) 
    COUNT_READS.out.cells_pf.view()

    cells_pf = COUNT_READS.out.cells_pf
    .splitCsv( sep: "\t", header: false)
    .map { row -> tuple(row[0], row[1], row[2]) }
    .view()

    SPLIT_BAMS(cells_pf) \
    | collect \
    | COPYKIT

    copykit_clones = 
        COPYKIT.out.copykit_tsv 
        .collectFile(name: 'copykit_clones.txt', newLine: true)

    //Initiate amethyst object per sample    
    AMETHYST_INIT(indir, copykit_clones)
}