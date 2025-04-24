//Nextflow pipeline for processing Navin lab spatial curio+multiome libraries//

// Declare syntax version
nextflow.enable.dsl=2

// Script parameters
params.projDir="/data/rmulqueen/projects/scalebio_dcis" //Project directory for processing
params.runDir="/data/rmulqueen/projects/scalebio_dcis/data/240202_prelim1/scale_data" //Output from scalemethyl pipeline
params.outputPrefix="scale" //Output from scalemethyl pipeline
params.src = "/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/" //src directory for cloned repo
params.min_cnv_count= "100000" //minimum readcount to try copykit per cell


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

	input:
		path runDir
	output:
		path("cells_pf.tsv"), emit: cells_pf
        path("unique_read_counts.tsv"), emit: uniq_read_counts

    script:
		"""
        #use a function to parallelize
        count_reads() { 
        bam_in=\$1
        samtools view \$1 | awk -v b=\$bam_in '{split(\$1,a,":"); print a[8],b}' | sort | uniq -c | sort -k1,1n
        }
        export -f count_reads

        parallel -j ${task.cpus} count_reads ::: \$(find ${params.runDir}/alignments -maxdepth 5 -name '*bam') | sort -k1,1n > unique_read_counts.tsv
        awk 'OFS="," {if(\$1>${params.min_cnv_count}) print \$1,\$2,\$3}' unique_read_counts.tsv > cells_pf.tsv
		"""
}


process SPLIT_BAMS {
	//Take output of COUNT_READS with each line being a tuple, then split from source bam to single cell bam
    publishDir "${params.runDir}/cnv/sc_bam", mode: 'copy', overwrite:true, pattern: "*bam"
	maxForks "${params.max_cpus}"
	input:
		tuple val(read_count),val(idx),path(bam)
	output:
		path("*bam")
	script:
	"""
        outprefix="${bam.Name}"
        outprefix=\$(echo \$outprefix | sed -e 's/.dedup.bam//g' -)
        ((samtools view -H $bam) && (samtools view $bam | awk -v i=$idx '{split(\$1,a,":"); if(a[8]==i); print \$0}')) | samtools view -bS > \${outprefix}.${idx}.bam
	"""
}


process COPYKIT { 
	//Run COPYKIT on single cell bam files
	cpus "${params.max_cpus}"
	containerOptions "--bind ${params.src}:/src/,${params.projDir}"
    publishDir "${params.runDir}/cnv/copykit", mode: 'copy', overwrite:true, pattern: "*{tsv,rds}"
    publishDir "${params.runDir}/cnv/copykit_plots", mode: 'copy', overwrite:true, pattern: "*pdf"

	label 'copykit'
	input:
		path bams
	output:
		path("*.scCNA.tsv"), emit: copykit_tsv
        path("*.rds"), emit: copykit_rds
        path("*.pdf"), emit: cnv_plots
		path("${params.outputPrefix}.merged.cnv.tsv"), emit: copykit_clones
    script:
		"""
        Rscript /src/copykit_cnvcalling.R \
            --input_dir . \
            --output_dir . \
            --output_prefix ${params.outputPrefix} \
            --task_cpus ${task.cpus}

		cat *.scCNA.tsv > ${params.outputPrefix}.merged.cnv.tsv
		"""
}

//Initiate amethyst object with added CNV clones
process AMETHYST_INIT {
	//TRIM READS OF ADAPTERS AND KNOWN METHYLATED REGIONS (GAP FILLS)
	cpus "${params.max_cpus}"
    publishDir "${params.runDir}/amethyst/amethyst_plots", mode: 'copy', overwrite:true, pattern: "*pdf"
	publishDir "${params.runDir}/amethyst/", mode: 'copy', overwrite:true, pattern: "*rds"

	containerOptions "--bind ${params.src}:/src/,${params.runDir},${params.projDir},/home/rmulqueen/R/x86_64-conda-linux-gnu-library/4.4"
    //last container option only necessary because i didnt add optparse to SIF

	label 'amethyst'

	input:
		path copykit_clones

	output:
		path("*.pdf"), emit: amethyst_plots
		path("*.rds"), emit: amethyst_obj

	script:
		"""
        Rscript /src/amethyst_initial_processing.R \
        --input_dir ${params.runDir} \
        --output_prefix ${params.outputPrefix} \
		--cnv_clones ${copykit_clones} \
        --task_cpus ${task.cpus}
		"""
}


workflow {
    
    //Split bams and run copykit
    indir= channel.fromPath( params.runDir , type: 'dir')

    COUNT_READS(indir) 

    cells_pf = COUNT_READS.out.cells_pf
    .splitCsv( header: false)
    .map { row -> tuple(row[0], row[1], row[2]) }

    SPLIT_BAMS(cells_pf) \
    | collect \
    | COPYKIT        

    //Initiate amethyst object per sample    
    AMETHYST_INIT(COPYKIT.out.copykit_clones)
}

//Downstream:
//Merge amethyst objects, identify cell types
//Output amethyst for methyltree processing