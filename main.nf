#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*-----------------------------------------------------------------------------------------

	Please put your raw paired-end FASTQ files into the folder "reads/" and ensure the 
		that the files end with the suffix ".fastq.gz". The sample name should precede
		the read pair  designation. e.g., the files for SampleA would be 
		"SampleA_1.fastq.gz" and "SampleA_2.fastq.gz".
	
	Usage:
		nextflow run --genome <file> --transcriptome <file> main.nf
	
	Inputs:
	
		(Required)
		
		--genome:			Specifies the reference genome used for
							generating the index decoys.
						
		--transcriptome:	Specifies the reference transcriptome.	
			
		(Optional)
		
		--kmer:				Specifies the k-mer length used for
							indexing the reference transcriptome.
								Defaults to "31"
							
	Outputs:
		results/    
		|
		|---fastqc_raw/
		|
		|---fastqc_trimmed/
		|
		|---quantified_transcripts/
		|
		`---trimmed_reads/
		
-----------------------------------------------------------------------------------------*/		
 
 
/*---------------------------
	Define the parameters
---------------------------*/

params.genome = false
params.transcriptome = false
params.kmer = 31


/*-------------------------
	Define the workflow
-------------------------*/


workflow {

	if (params.genome == false | params.transcriptome == false) {
		log.info """
		
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ERROR! Please define the reference genome and transcriptome FASTA files.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-------------------------------------------------------------------------
Usage:
  	nextflow run --genome <file> --transcriptome <file> main.nf
  
Inputs:

	(Required)
	
	--genome:		Specifies the reference genome used for \
					generating the index decoys.
					
	--transcriptome:	Specifies the reference transcriptome.	
		
	(Optional)
  	
	--kmer:			Specifies the k-mer length used for \
						indexing the reference transcriptome.
					Defaults to "31"
-------------------------------------------------------------------------

"""
		exit 0
	}			
	else {
	reads_ch = Channel.fromFilePairs('reads/*_{1,2}.fastq.gz')
	FASTQC(reads_ch)
	trim_ch=TRIM(reads_ch)
	FASTQC_TRIM(trim_ch)
	genome_ch = Channel.fromPath( params.genome )
	transcripts_ch = Channel.fromPath( params.transcriptome )
	index_ch=INDEX(genome_ch, transcripts_ch)
	QUANT(index_ch,trim_ch)	
	}
}


/*-------------------------------------------------
	Define the processes called in the workflow
-------------------------------------------------*/


process FASTQC {
	tag "Processing $sampleId"
	
    publishDir './results/fastqc_raw', mode:'copy'

	input:
		tuple val(sampleId), file(reads_ch)
	
	output:
		tuple val(sampleId), path("${sampleId}_{1,2}_fastqc.html")
		tuple val(sampleId), path("${sampleId}_{1,2}_fastqc.zip")
	
    script:
		"""
		fastqc ${reads_ch} -q -t $task.cpus
		"""
}


process TRIM {
	tag "Processing $sampleId"
	
	publishDir './results/trimmed_reads/', pattern: '*_val_*.fq.gz', mode:'copy'

	input:
		tuple val(sampleId), file(reads_ch)
	
	output:
		tuple val(sampleId), path("${sampleId}_{1,2}_val_{1,2}.fq.gz"), emit: trimmed_reads
	
    script:
		"""
		trim_galore --max_n 4 --cores $task.cpus --paired --gzip ${reads_ch}
		"""
}


process FASTQC_TRIM {
	tag "Processing $sampleId"
	
    publishDir './results/fastqc_trimmed', mode:'copy'

	input:
		tuple val(sampleId), path(trimmed_reads)
	
	output:
		tuple val(sampleId), path("${sampleId}_{1,2}_val_{1,2}_fastqc.html")
		tuple val(sampleId), path("${sampleId}_{1,2}_val_{1,2}_fastqc.zip")
	
    script:
		"""
		fastqc ${trimmed_reads[0]} ${trimmed_reads[1]} -q -t $task.cpus
		"""
}


process INDEX {
	tag "Indexing reference genome"
    
	input:
		path genome_ch
		path transcripts_ch
	
    output:
		path index
	
    script:
		"""
		zcat ${genome_ch} | grep "^>" | cut -d " " -f 1 | sed 's/>//g' > decoys
		zcat ${transcripts_ch} ${genome_ch} > gentrome
		salmon index --threads $task.cpus -t gentrome -d decoys -i index -k $params.kmer
		"""
}


process QUANT {
	tag "Quantifying transcripts for $sampleId"
	
	publishDir "results/quantified_transcripts", mode:'copy'
	
    input:
		each index
		tuple val(sampleId), path(trimmed_reads)
	
	output:
		path "$sampleId/quant.sf"
	
    script:
		"""
		salmon quant \
			--threads $task.cpus --libType A -i $index \
			-1 ${trimmed_reads[0]} -2 ${trimmed_reads[1]} -o $sampleId
		"""
}