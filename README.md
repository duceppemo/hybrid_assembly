## Notes

Current version needs the pacbio image files to be preprocess, meaning that the "filtered_subreads.fastq" file is available. It also requires a set of paired-end Illumina file for polishing. It's also aimed at assembling bacteria genomes.

## Dependencies

Software list:
* Java
* FastQC
* Centrifuge
* Krona tools
* Firefox
* Canu
* BBtools
* BWA
* Samtools
* Racon
* Minimap
* Pilon
* Picard tools
* Mauve
* Circlator
* Prokka
* BLAST+
* GNU parallel
* Entrez Direct (NCBI)

Other scripts:
* formatFasta.pl
* http_post.pl
* get_assemblies.sh
* findClosest.sh

## Usage

Youy need to update the "User Defined" section of the script prior running it. It will use all the available CPUs available.

## Workflow

1. Pacbio read quality control (QC) and taxonomic diversity analysis
2. De nove assembly of pacbio reads.
3. Polishing of assembly
4. Ordering of contig(s)
5. QC of assembly
6. Genome annotation



