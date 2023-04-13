# RNAseq (v1.0.0)

This Nextflow pipeline quantifies transcript abundance from raw sequencing data.

## Contents
- [Workflow](#workflow)
- [Usage](#usage)
- [Dependencies](#dependencies)
- [Example](#example)
- [License](#license)
- [Citations](#citations)

Workflow
--------
The user provides raw paired-end sequencing reads. [**Trim Galore v0.6.10**](https://github.com/FelixKrueger/TrimGalore) is used to perform adapter and quality trimming. [**FastQC v0.11.9**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is used to provide some quality control information for both the raw and trimmed reads. [**Salmon v1.9.0**](https://salmon.readthedocs.io/en/latest/) is used to index the reference transcriptome and quantify transcript abundance.

Usage
-----
	Please put your raw paired-end FASTQ files into a new folder titled "reads/" and ensure 
		that the files end with the suffix ".fastq.gz". The sample name should precede
		the read pair  designation. e.g., the files for SampleA would be 
		"SampleA_1.fastq.gz" and "SampleA_2.fastq.gz".
	
	Usage:
		nextflow run --genome <file> --transcriptome <file> main.nf
	
	Inputs:
	
		(Required)
		
		--genome:		Specifies the reference genome used for
					generating the index decoys.
						
		--transcriptome:	Specifies the reference transcriptome.	
			
		(Optional)
		
		--kmer:			Specifies the k-mer length used for
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
    
Dependencies 
------------
The following programs must be installed on your computer:
* [**Nextflow**](https://github.com/nextflow-io/nextflow) (v22.10.6 or higher)
* [**Singularity**](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) (v3.8.7 or higher)
 
An easy way to install these is with [**Mamba**](https://github.com/mamba-org/mamba). This doesn't require the user to be root so can help to manage software packages on HPC systems.

With Mamba installed, one can fetch the dependencies using the following commands:
````
```
mamba install -c bioconda nextflow
mamba install -c conda-forge singularity
```
````

Containers for software used in the workflow are hosted on [**Docker Hub**](https://hub.docker.com/u/bryoinformatics) and will be downloaded automatically.

Example
------------
A working example can be performed using the sample data distributed with this package

Step 1: Download the repository and change directory.
````
```
git clone https://github.com/bryoinformatics/RNAseq.git
cd RNAseq
```
````

Step 2: Make a new directory called "reads/" and move the sample FASTQ files into it.
````
```
mkdir reads
cp sample_data/*.fastq.gz reads/
```
````

Step 3: Run the pipeline.
````
```
nextflow run --genome sample_data/genome.fa.gz --transcriptome sample_data/transcripts.fa.gz main.nf
```
````

License
-------
This workflow is released under a GNU General Public License (v3.0).

Citations
---------
* Andrews, S. (2019). FastQC: A quality control tool for high throughput sequence data. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
* Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), Article 4. https://doi.org/10.1038/nbt.3820
* Krueger, F. (2023). Babraham Bioinformatics—Trim Galore! (v0.6.10). https://zenodo.org/badge/latestdoi/62039322
* Kurtzer, G. M., Sochat, V., & Bauer, M. W. (2017). Singularity: Scientific containers for mobility of compute. PLOS ONE, 12(5), e0177459. https://doi.org/10.1371/journal.pone.0177459
* Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods, 14(4), Article 4. https://doi.org/10.1038/nmeth.4197


