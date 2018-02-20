#!/bin/bash

#script version
version="0.1"

# Process and assemble PacBio data from the 
# Includes stats, qc and de novo assembly


# ### IMPORTANT ###

# # you must be loging in as super upser "smrtportal" to run this script!
# # su -l smrtanalysis
# if [ $(whoami) != "smrtanalysis" ]; then
#     echo -e "You must loged in as \"smrtanalysis\" to run this script.\nUse: su -l smrtanalysis"
#     # exit 1
# fi


######################
#                    #
#    User Defined    #
#                    #
######################


# Assembly name
export prefix="MBWGS070"

#Annotation
kingdom="Bacteria"
genus="Mycobacterium"
species="bovis"
strain="MBWGS070"
gram="pos"
locus_tag="XXX"
centre="OLF"

# Analysis folder
baseDir=""${HOME}"/Desktop/Mbovis_canu/"$prefix""

# Pacbio reads
filtered_subreads="/media/6tb_raid10/data/Mbovis/canada/pacbio/MBWGS070/MBWGS070.filtered_subreads.fastq.gz"
ccs_reads="/media/6tb_raid10/data/Mbovis/canada/pacbio/MBWGS070/MBWGS070.ccs.fastq.gz"

# Illumina paired-end data
r1="/media/6tb_raid10/data/Mbovis/canada/illumina/MBWGS070_R1.fastq.gz"
r2="/media/6tb_raid10/data/Mbovis/canada/illumina/MBWGS070_R2.fastq.gz"

# Database to use for metagomic analysis of raw data (contamination)
# db="/media/6tb_raid10/db/centrifuge/p_compressed+h+v"
# db="/media/6tb_raid10/db/centrifuge/nt"
db="/media/6tb_raid10/db/centrifuge/2017-10-12_bact_vir_h"

# Where programs are installed
prog=""${HOME}"/prog"
scripts=""${HOME}"/scripts"

#genome size (g|m|k)
size="4350k"

# For SPAdes
kmer="21,33,55,77,99,127"

# For assembly trimming
smallest_contig=1000

#BUSCO
# busco_db="/media/6tb_raid10/db/busco/bacteria_odb9.tar.gz"
# busco_species="Escherichia coli"


######################
#                    #
#     Resources      #
#                    #
######################


# Computer performance
cpu=$(nproc) #total number of cores
mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB
memJava="-Xmx"$mem"g"


#######################
#                     #
#   Data Stucture     #
#                     #
#######################


# Folder structure
logs=""${baseDir}"/logs"
qc=""${baseDir}"/qc"
fastq=""${baseDir}"/fastq"
assemblies=""${baseDir}"/assemblies"
polished=""${baseDir}"/polished"
ordered=""${baseDir}"/ordered"
trimmed=""${baseDir}"/trimmed"
corrected=""${baseDir}"/corrected"
normalized=""${baseDir}"/normalized"
merged=""${baseDir}"/merged"
aligned=""${baseDir}"/aligned"
export annotation=""${baseDir}"/annotation"

# Create folders if do not exist
[ -d "$baseDir" ] || mkdir -p "$baseDir"
[ -d "$logs" ] || mkdir -p "$logs"
[ -d "$qc" ] || mkdir -p "$qc"
[ -d "$fastq" ] || mkdir -p "$fastq"
[ -d "$assemblies" ] || mkdir -p "$assemblies"
[ -d "$polished" ] || mkdir -p "$polished"
[ -d "$trimmed" ] || mkdir -p "$trimmed"
[ -d "$corrected" ] || mkdir -p "$corrected"
[ -d "$merged" ] || mkdir -p "$merged"
[ -d "$normalized" ] || mkdir -p "$normalized"
[ -d "$aligned" ] || mkdir -p "$aligned"
[ -d "$ordered" ] || mkdir -p "$ordered"
[ -d "$annotation" ] || mkdir -p "$annotation"


################
#              #
#     Log      #
#              #
################


# Date
echo -e "$(date)\n" | tee "${logs}"/log.txt
echo -e "User: $(whoami)" | tee -a "${logs}"/log.txt
echo -e "Processors: "$cpu"" | tee -a "${logs}"/log.txt
echo -e "Memory: "$mem"G" | tee -a "${logs}"/log.txt

#script version
echo -e "\npacbio_assembly.sh version "$version"\n" | tee -a "${logs}"/log.txt  # $0


####################
#                  #
#   Dependencies   #
#                  #
####################


# Check for dependencies and log versions
echo "Checking for dependencies..."

# java
if hash java 2>/dev/null; then 
    java -version 2>&1 1>/dev/null | grep "java version" | tr -d '"' | tee -a "${logs}"/log.txt
else
    echo >&2 "java was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# FastQC
if hash fastqc 2>/dev/null; then 
    fastqc -v | tee -a "${logs}"/log.txt
else
    echo >&2 "fastQC was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Centrifuge
if hash centrifuge 2>/dev/null; then 
    v=$(centrifuge --version | grep -F "version")
    echo "centrifuge $v" | tee -a "${logs}"/log.txt
else
    echo >&2 "centrifuge was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Check Centrifuge datase
if [ -s "${db}".1.cf ]; then
    echo -e "Centrifuge database: $(basename "$db")" | tee -a "${logs}"/log.txt
else
    echo "Could not find the centrifuge database"
    exit 1
fi

# Canu
if hash canu 2>/dev/null; then
    canu --version | tee -a "${logs}"/log.txt
else
    echo >&2 "canu was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# BWA
if hash bwa 2>/dev/null; then
    v=$(bwa 2>&1 1>/dev/null | grep -F "Version" | cut -d " " -f 2)
    echo "bwa $v" | tee -a "${logs}"/log.txt
else
    echo >&2 "bwa was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Samtools
if hash samtools 2>/dev/null; then
    samtools --version | grep -F 'samtools' | tee -a "${logs}"/log.txt
else
    echo >&2 "samtools was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Racon
if hash racon 2>/dev/null; then
    echo "racon v2017-04-18" | tee -a "${logs}"/log.txt
else
    echo >&2 "racon was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# minimap
if hash minimap 2>/dev/null; then
    v=$(minimap -V)
    echo "minimap v"${v}"" | tee -a "${logs}"/log.txt
else
    echo >&2 "minimap was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Pilon
java -jar "${prog}"/pilon/pilon-dev.jar &> /dev/null
if [ $? -ne 0 ] ; then  # "$?" = last execution status. 0 = all good
    echo >&2 "pilon was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
else
    java -jar "${prog}"/pilon/pilon-dev.jar | grep -F 'version' | tee -a "${logs}"/log.txt
fi

#BBtools (bbduk, tadpole, bbnorm, bbmerge)


#Circlator
if hash circlator; then
    version=$(circlator version)
    echo "Circlator "$version"" | tee -a "${logs}"/log.txt
else
    echo >&2 "circlator was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi 

#Prokka
if hash prokka 2>/dev/null; then
    prokka --version 2>&1 1>/dev/null | tee -a "${logs}"/log.txt
else
    echo >&2 "prokka was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

#Mauve
echo "Mauve version \"Snapshot 2015-02-13\"" | tee -a "${logs}"/log.txt


#############
#           #
#   FASTQ   #
#           #
#############


# Populate fastq files
# Reorder for better compression and speed
clumpify.sh "$memJava" \
    in="$r1" \
    in2="$r2" \
    out="${fastq}"/"${prefix}"_R1.fastq.gz \
    out2="${fastq}"/"${prefix}"_R2.fastq.gz \
    reorder \
    ziplevel=9 \
    2> >(tee -a "${logs}"/clumpify_illumina.txt)

#Do it for pacbio reads too?
clumpify.sh "$memJava" \
    in="$ccs_reads" \
    out="${fastq}"/"${prefix}"_ccs.fastq.gz \
    reorder \
    ziplevel=9 \
    2> >(tee -a "${logs}"/clumpify_ccs.txt)

clumpify.sh "$memJava" \
    in="$filtered_subreads" \
    out="${fastq}"/"${prefix}"_subreads.fastq.gz \
    reorder \
    ziplevel=9 \
    2> >(tee -a "${logs}"/clumpify_subreads.txt)


############
#          #
#   HGAP   #
#          #
############


# Run HGAP pipeline from the command line (instead of using the SMRT portal)
# TODO


################
#              #
#   QC - Raw   #
#              #
################


function run_fastqc()
{
    #first argument is output folder
    #all other are the fastq files to process
    args=("$@")  #store arguments in array
    output="${args[0]}"
    input=("${args[@]:1}")

    [ -d "$output" ] || mkdir -p "$output"

    fastqc \
        --o  "$output" \
        --noextract \
        --threads "$cpu" \
        ${input[@]}

    #Merge all FastQC reports together
    multiqc \
        -o "$output" \
        -n merged_reports.html \
        "$output"
}

# Create folder to store report
[ -d "${qc}"/fastqc/pacbio/raw ] || mkdir -p "${qc}"/fastqc/pacbio/raw
[ -d "${qc}"/fastqc/illumina/raw ] || mkdir -p "${qc}"/fastqc/illumina/raw

echo "Running FastQC on PacBio reads..."
run_fastqc \
    "${qc}"/fastqc/pacbio/raw \
    "${fastq}"/"${prefix}"_ccs.fastq.gz \
    "${fastq}"/"${prefix}"_subreads.fastq.gz

echo "Running FastQC on Illumina reads..."
run_fastqc \
    "${qc}"/fastqc/illumina/raw \
    "${fastq}"/"${prefix}"_R1.fastq.gz \
    "${fastq}"/"${prefix}"_R2.fastq.gz


######################
#                    #
#   Centrifuge raw   #
#                    #
######################


function run_centrifuge()
{
    #first argument is output file
    #all other are the fastq files to process
    args=("$@")  #store arguments in array
    output="${args[0]}"
    input=("${args[@]:1}")

    outFolder=$(dirname "$output")
    [ -d "$outFolder" ] || mkdir -p "$outFolder"

    # echo $input
    # echo $output
    # echo "${input[@]}"
    # echo "${#input[@]}"

    if [ "${#input[@]}" -eq 1 ]; then  #single-end
        option="-U "${input[0]}""
    elif [ "${#input[@]}" -eq 2 ]; then  #paired-end
        option="-1 "${input[0]}" -2 "${input[1]}""
    else  #error
        echo "Centrifuge can only process a single fastq file (unpaired) or one set of paired-end reads at the time"
        exit 1
    fi
    
    #build the command
    com="centrifuge \
        -p "$cpu" \
        -t \
        --seed "$RANDOM" \
        -x "$db" \
        "$option" \
        > "$output""

    # run the command
    eval "$com"

    if [ $? -eq 0 ]; then  #if no error from centrifuge
        # Prepare centrifuge result for display with Krona
        cat "$output" | \
            cut -f 1,3 | \
            ktImportTaxonomy /dev/stdin -o "${output%.tsv}".html

        # Visualize the resutls in web browser of choice
        firefox file://"${output%.tsv}".html &
    else
        echo "An error occured while processing "$input""
        exit 1
    fi
}

run_centrifuge \
    "${qc}"/centrifuge/"${prefix}"_subreads.tsv \
    "$filtered_subreads"

run_centrifuge \
    "${qc}"/centrifuge/"${prefix}"_ccs.tsv \
    "$ccs_reads"

run_centrifuge \
    "${qc}"/centrifuge/"${prefix}"_pe.tsv \
    "${fastq}"/"${prefix}"_R1.fastq.gz \
    "${fastq}"/"${prefix}"_R2.fastq.gz 


###########################
#                         #
#   Read Pre-Processing   #
#                         #
###########################


# Illumina

# #normalize reads
# #NOTE:  All depth parameters control kmer depth, not read depth.
# bbnorm.sh "$memJava" \
#     in="${fastq}"/"${prefix}"_R1.fastq.gz \
#     in2="${fastq}"/"${prefix}"_R2.fastq.gz \
#     out="${trimmed}"/"${prefix}"_Normalized_1P.fastq.gz \
#     out2="${trimmed}"/"${prefix}"_Normalized_2P.fastq.gz \
#     target=100 \
#     mindepth=2 \
#     2> >(tee "${logs}"/illumina_tile_filtering.txt)

#Removing low-quality regions
filterbytile.sh "$memJava" \
    in="${fastq}"/"${prefix}"_R1.fastq.gz \
    in2="${fastq}"/"${prefix}"_R2.fastq.gz \
    out="${trimmed}"/"${prefix}"_Filtered_1P.fastq.gz \
    out2="${trimmed}"/"${prefix}"_Filtered_2P.fastq.gz \
    2> >(tee "${logs}"/illumina_tile_filtering.txt)

#QC
[ -d "${qc}"/fastqc/illumina/tile ] || mkdir -p "${qc}"/fastqc/illumina/tile
run_fastqc \
    "${qc}"/fastqc/illumina/tile \
    "${trimmed}"/"${prefix}"_Filtered_1P.fastq.gz \
    "${trimmed}"/"${prefix}"_Filtered_2P.fastq.gz

#Quality trimming and adapter trimming
bbduk.sh "$memJava" \
    in="${trimmed}"/"${prefix}"_Filtered_1P.fastq.gz \
    in2="${trimmed}"/"${prefix}"_Filtered_2P.fastq.gz \
    ref="${prog}"/bbmap/resources/adapters.fa \
    ktrim=r k=23 mink=11 hdist=1 tbo tpe \
    qtrim=lr trimq=10 \
    minlen=64 \
    out="${trimmed}"/"${prefix}"_Trimmed_1P.fastq.gz \
    out2="${trimmed}"/"${prefix}"_Trimmed_2P.fastq.gz \
    ziplevel=9 \
    ordered=t \
    2> >(tee "${logs}"/illumina_trimming.txt)

#QC
[ -d "${qc}"/fastqc/illumina/trimmed ] || mkdir -p "${qc}"/fastqc/illumina/trimmed
run_fastqc \
    "${qc}"/fastqc/illumina/trimmed \
    "${trimmed}"/"${prefix}"_Trimmed_1P.fastq.gz \
    "${trimmed}"/"${prefix}"_Trimmed_2P.fastq.gz

#Removing synthetic artifacts and spike-ins
#Add to "ref" any contaminants found with centrifuge
bbduk.sh "$memJava" \
    in="${trimmed}"/"${prefix}"_Trimmed_1P.fastq.gz \
    in2="${trimmed}"/"${prefix}"_Trimmed_2P.fastq.gz \
    ref="${prog}"/bbmap/resources/sequencing_artifacts.fa.gz,"${prog}"/bbmap/resources/phix174_ill.ref.fa.gz \
    k=31 \
    out="${trimmed}"/"${prefix}"_Cleaned_1P.fastq.gz \
    out2="${trimmed}"/"${prefix}"_Cleaned_2P.fastq.gz \
    ziplevel=9 \
    ordered=t \
    2> >(tee "${logs}"/illumina_cleaning.txt)

#QC
[ -d "${qc}"/fastqc/illumina/cleaned ] || mkdir -p "${qc}"/fastqc/illumina/cleaned
run_fastqc \
    "${qc}"/fastqc/illumina/cleaned \
    "${trimmed}"/"${prefix}"_Cleaned_1P.fastq.gz \
    "${trimmed}"/"${prefix}"_Cleaned_2P.fastq.gz

#Correcting Illumina paired-end reads
#Phase 1
bbmerge.sh "$memJava" \
    in="${trimmed}"/"${prefix}"_Cleaned_1P.fastq.gz \
    in2="${trimmed}"/"${prefix}"_Cleaned_2P.fastq.gz \
    out="${corrected}"/"${prefix}"_Cor1_1P.fastq.gz \
    out2="${corrected}"/"${prefix}"_Cor1_2P.fastq.gz \
    ecco=t \
    mix=t \
    verystrict=t \
    ordered=t \
    ziplevel=9 \
    ordered ihist="${logs}"/illumina_ihist_corr_merge.txt \
    2> >(tee "${logs}"/illumina_corr_merging.txt)

#Phase2
clumpify.sh "$memJava" \
    in="${corrected}"/"${prefix}"_Cor1_1P.fastq.gz \
    in2="${corrected}"/"${prefix}"_Cor1_2P.fastq.gz \
    out="${corrected}"/"${prefix}"_Cor2_1P.fastq.gz \
    out2="${corrected}"/"${prefix}"_Cor2_2P.fastq.gz \
    ecc=t \
    passes=4 \
    reorder=t \
    ziplevel=9 \
    2> >(tee "${logs}"/illumina_corr_clumpify.txt)

#Phase3
tadpole.sh "$memJava" \
    in="${corrected}"/"${prefix}"_Cor2_1P.fastq.gz \
    in2="${corrected}"/"${prefix}"_Cor2_2P.fastq.gz \
    out="${corrected}"/"${prefix}"_Cor3_1P.fastq.gz \
    out2="${corrected}"/"${prefix}"_Cor3_2P.fastq.gz \
    ecc=t \
    k=62 \
    threads="$cpu" \
    mode=correct \
    2> >(tee "${logs}"/illumina_corr_tadpol.txt)

#Merging Illumina overlapping paried-end reads
bbmerge-auto.sh "$memJava" \
    in="${corrected}"/"${prefix}"_Cor3_1P.fastq.gz \
    in2="${corrected}"/"${prefix}"_Cor3_2P.fastq.gz \
    out="${merged}"/"${prefix}"_merged.fastq.gz \
    outu="${merged}"/"${prefix}"_unmerged_1P.fastq.gz \
    outu2="${merged}"/"${prefix}"_unmerged_2P.fastq.gz \
    strict=t \
    k=93 \
    extend2=80 \
    rem=t \
    ordered=t \
    ziplevel=9 \
    ihist="${logs}"/illumina_ihist_merge.txt \
    2> >(tee "${logs}"/illumina_merge.txt)

#QC
[ -d "${qc}"/fastqc/illumina/merged ] || mkdir -p "${qc}"/fastqc/illumina/merged
run_fastqc \
    "${qc}"/fastqc/illumina/merged \
    "${merged}"/"${prefix}"_merged.fastq.gz \
    "${merged}"/"${prefix}"_unmerged_1P.fastq.gz \
    "${merged}"/"${prefix}"_unmerged_2P.fastq.gz

#PacBio

#Remove SMRT adapter and split chimera
removesmartbell.sh "$memJava" \
    in="$ccs_reads" \
    out="${trimmed}"/"$prefix"_ccs_trimmed.fastq.gz \
    split=t \
    &> >(tee "${logs}"/ccs_trimmed.txt)

removesmartbell.sh "$memJava" \
    in="$filtered_subreads" \
    out="${trimmed}"/"$prefix"_subreads_trimmed.fastq.gz \
    split=t \
    &> >(tee "${logs}"/subreads_trimmed.txt)

#QC
[ -d "${qc}"/fastqc/pacbio/adapter ] || mkdir -p "${qc}"/fastqc/pacbio/adapter
run_fastqc \
    "${qc}"/fastqc/pacbio/adapter \
    "${trimmed}"/"$prefix"_ccs_trimmed.fastq.gz \
    "${trimmed}"/"$prefix"_subreads_trimmed.fastq.gz

#Remove internal control from pacbio data
#Add to "ref" any contaminants found with centrifuge
bbduk.sh "$memJava" \
    in="${trimmed}"/"$prefix"_ccs_trimmed.fastq.gz \
    ref=/media/6tb_raid10/ref/MG495226.1.fasta \
    k=31 \
    out="${trimmed}"/"${prefix}"_ccs_Cleaned.fastq.gz \
    ziplevel=9 \
    ordered=t \
    2> >(tee "${logs}"/pacbio_ccs_cleaning.txt)

bbduk.sh "$memJava" \
    in="${trimmed}"/"$prefix"_subreads_trimmed.fastq.gz \
    ref=/media/6tb_raid10/ref/MG495226.1.fasta \
    k=31 \
    out="${trimmed}"/"$prefix"_subreads_Cleaned.fastq.gz \
    ziplevel=9 \
    ordered=t \
    2> >(tee "${logs}"/pacbio_subreads_cleaning.txt)

#QC
[ -d "${qc}"/fastqc/pacbio/phix ] || mkdir -p "${qc}"/fastqc/pacbio/phix
run_fastqc \
    "${qc}"/fastqc/pacbio/phix \
    "${trimmed}"/"$prefix"_ccs_Cleaned.fastq.gz \
    "${trimmed}"/"$prefix"_subreads_Cleaned.fastq.gz

#self-correction
canu -correct \
    -p "$prefix" \
    -d "$corrected" \
    genomeSize="$size" \
    -pacbio-raw "${trimmed}"/"${prefix}"_subreads_Cleaned.fastq.gz
    # output = "${corrected}"/"${prefix}".correctedReads.fasta.gz

#trim
canu -trim \
    -p "$prefix" \
    -d "${corrected}"/trimmed \
    genomeSize="$size" \
    -pacbio-corrected "${corrected}"/"${prefix}".correctedReads.fasta.gz
    # output = "${corrected}"/trimmed/"${prefix}".trimmedReads.fasta.gz

function plot_read_length_distribution()
{
    name=$(basename "$1")
    #from bbmap
    readlength.sh \
        in="$1" \
        out="${qc}"/"${name}".length.txt \
        bin=500

    #plot to file
    gnuplot -e "set term png; \
        set output '"${qc}"/"${name}".length.png'; \
        set xlabel 'Read Length (bp)'; \
        set ylabel 'Count'; \
        plot '"${qc}"/"${name}".length.txt' using 1:2 with linespoints" 
}

#create plots of pacbio read length distribution befor and after correction
plot_read_length_distribution "$filtered_subreads"  # before correction
plot_read_length_distribution "${corrected}"/trimmed/"${prefix}".trimmedReads.fasta.gz  # after correction


################
#              #
#   Assembly   #
#              #
################


#de novo assembly
echo "De novo assembly of PacBio reads using canu..."
canu -assemble \
    -p "$prefix" \
    -d "$assemblies" \
    genomeSize="$size" \
    -pacbio-corrected "${corrected}"/trimmed/"${prefix}".trimmedReads.fasta.gz

mv "${assemblies}"/"${prefix}".contigs.fasta \
    "${assemblies}"/"${prefix}"_canu.fasta

genome=""${assemblies}"/"${prefix}"_canu.fasta"


#################
#               #
#   Contig ID   #
#               #
#################


function run_blast()
{
    blastn \
    -db nt \
    -query "$1" \
    -out "${1%.fasta}".blastn.tsv \
    -evalue "1e-10" \
    -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
    -num_threads "$cpu" \
    -max_target_seqs 1 \
    -max_hsps 1

    echo -e "qseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" \
        > "${1%.fasta}".blastn.tsv.tmp
    cat "${1%.fasta}".blastn.tsv >> "${1%.fasta}".blastn.tsv.tmp
    mv "${1%.fasta}".blastn.tsv.tmp "${1%.fasta}".blastn.tsv
}

run_blast "$genome"

#inspect manually blast output
#Remove any contigs from contaminant
#Remove low coverage or short contigs 


#######################
#                     #
#      Polishing      #
#                     #
#######################


# Illumina reads are used to polish the PacBio de novo assembly
# They are then used for a new de novo assembly with SPAdes,
# using the polished PacBio assembly as template for scaffolding (--trusted-contigs)
# PacBio subreads are also used for scaffolding
# PacBio Circular Consensus Sequences (CCS) are used as single-end reads for the assembly


# Introduces lots of indels in the assembly if using uncorrected subreads
# Helps reduce size of assembly
# Polish assembly with racon
echo "Mapping PacBio reads onto canu contigs using minimap..."
minimap2 \
    -t "$cpu" \
    "$genome" \
    "${trimmed}"/"$prefix"_subreads_Cleaned.fastq.gz \
    > "${polished}"/"${prefix}"_minimap_overlaps.paf

echo "Correcting contigs with PacBio reads using racon..."
racon \
    -t "$cpu" \
    "${trimmed}"/"$prefix"_subreads_Cleaned.fastq.gz \
    "${polished}"/"${prefix}"_minimap_overlaps.paf \
    "$genome" \
    "${polished}"/"${prefix}"_racon.fasta

genome="${polished}"/"${prefix}"_racon.fasta

# #convert fastq to fasta
# zcat "${corrected}"/"${prefix}"_subreads_Corrected.fastq.gz \
#     | sed -n '1~4s/^@/>/p;2~4p' \
#     > "${corrected}"/"${prefix}"_subreads_Corrected.fasta

# #map pacbio reads back to assembly
# blasr \
#    "${corrected}"/"${prefix}"_subreads_Corrected.fasta \
#    "$genome" \
#    --nproc "$cpu" \
#    --out "${polished}"/"${prefix}"_subreads_Corrected.sam \
#    --bam

# samtools view \
#     -@ "$cpu" \
#     -b -h \
#     -T "$genome" \
#     "${polished}"/"${prefix}"_subreads_Corrected.sam \
#     > "${polished}"/"${prefix}"_subreads_Corrected.bam

# bamtools sort \
#     -in "${polished}"/"${prefix}"_subreads_Corrected.bam \
#     -out "${polished}"/"${prefix}"_subreads_sorted.bam

# bamtools index \
#     -in "${polished}"/"${prefix}"_subreads_sorted.bam

# source activate quiver

# quiver \
#     -j "$cpu" \
#     -r "$genome" \
#     -o "${polished}"/"${prefix}"_quiver.fasta \
#     "${polished}"/"${prefix}"_subreads_sorted.bam

# source deactivate


function Polish()
{
    # Index genome
    bwa index "$1"

    #Align merged
    rg_pe_merged="@RG\tID:"${prefix}"_merged\tCN:"${centre}"\tLB:NexteraXT\tPL:ILLUMINA\tSM:"${prefix}""

    bwa mem -t "$cpu" -M -R "$rg_pe_merged" "$1" "$m" | \
        samtools view -@ "$cpu" -b -h -F 4 - | \
        samtools sort -@ "$cpu" -m 10G -o "${polished}"/"${prefix}"_merged.bam -

    samtools rmdup \
        -s \
        "${polished}"/"${prefix}"_merged.bam \
        "${polished}"/"${prefix}"_merged_nodup.bam

    samtools sort -@ "$cpu" -m 10G \
        -o "${polished}"/"${prefix}"_merged_sorted.bam \
        "${polished}"/"${prefix}"_merged_nodup.bam

    rm "${polished}"/"${prefix}"_merged.bam
    rm "${polished}"/"${prefix}"_merged_nodup.bam

    samtools index "${polished}"/"${prefix}"_merged_sorted.bam


    #Align unmerged
    rg_pe_unmerged="@RG\tID:"${prefix}"_unmerged\tCN:"${centre}"\tLB:NexteraXT\tPL:ILLUMINA\tSM:"${prefix}""

    bwa mem -t "$cpu" -M -R "$rg_pe_unmerged" "$1" "$mu1" "$mu2" | \
        samtools view -@ "$cpu" -b -h -F 4 - | \
        samtools sort -@ "$cpu" -m 10G -o "${polished}"/"${prefix}"_unmerged.bam -

    samtools rmdup \
        "${polished}"/"${prefix}"_unmerged.bam \
        "${polished}"/"${prefix}"_unmerged_nodup.bam

    samtools sort -@ "$cpu" -m 10G \
        -o "${polished}"/"${prefix}"_unmerged_sorted.bam \
        "${polished}"/"${prefix}"_unmerged_nodup.bam

    rm "${polished}"/"${prefix}"_unmerged.bam
    rm "${polished}"/"${prefix}"_unmerged_nodup.bam

    samtools index "${polished}"/"${prefix}"_unmerged_sorted.bam


    #Align pacbio CCS reads
    rg_ccs="@RG\tID:"${prefix}"_ccs\tCN:GQC\tLB:SMRTBELL\tPL:PACBIO\tSM:"${prefix}""

    bwa mem -x pacbio -t "$cpu" -M -R "$rg_ccs" "$1" "${trimmed}"/"${prefix}"_ccs_Cleaned.fastq.gz | \
        samtools view -@ "$cpu" -b -h -F 4 - | \
        samtools sort -@ "$cpu" -m 10G -o "${polished}"/"${prefix}"_ccs.bam -

    samtools rmdup \
        -s "${polished}"/"${prefix}"_ccs.bam \
        "${polished}"/"${prefix}"_ccs_nodup.bam

    samtools sort -@ "$cpu" -m 10G \
        -o "${polished}"/"${prefix}"_ccs_sorted.bam \
        "${polished}"/"${prefix}"_ccs_nodup.bam

    rm "${polished}"/"${prefix}"_ccs.bam
    rm "${polished}"/"${prefix}"_ccs_nodup.bam

    samtools index "${polished}"/"${prefix}"_ccs_sorted.bam

    #Align pacbio corrected subreads
    rg_subreads="@RG\tID:"${prefix}"_subreads\tCN:GQC\tLB:SMRTBELL\tPL:PACBIO\tSM:"${prefix}""

    bwa mem -x pacbio -t "$cpu" -M -R "$rg_subreads" "$1" "${corrected}"/trimmed/"${prefix}".trimmedReads.fasta.gz | \
        samtools view -@ "$cpu" -b -h -F 4 - | \
        samtools sort -@ "$cpu" -m 10G -o "${polished}"/"${prefix}"_subreads_Corrected.bam -

    samtools rmdup \
        -s "${polished}"/"${prefix}"_subreads_Corrected.bam \
        "${polished}"/"${prefix}"_subreads_Corrected_nodup.bam

    samtools sort -@ "$cpu" -m 10G \
        -o "${polished}"/"${prefix}"_subreads_Corrected_sorted.bam \
        "${polished}"/"${prefix}"_subreads_Corrected_nodup.bam

    rm "${polished}"/"${prefix}"_subreads_Corrected.bam
    rm "${polished}"/"${prefix}"_subreads_Corrected_nodup.bam

    samtools index "${polished}"/"${prefix}"_subreads_Corrected_sorted.bam

    #Merge all bam files
    samtools merge \
        -@ "$cpu" \
        - \
        "${polished}"/"${prefix}"_merged_sorted.bam \
        "${polished}"/"${prefix}"_unmerged_sorted.bam \
        "${polished}"/"${prefix}"_ccs_sorted.bam \
        "${polished}"/"${prefix}"_subreads_Corrected_sorted.bam | \
    samtools sort -@ "$cpu" -m 10G -o "${polished}"/"${prefix}"_all.bam
    
    samtools index "${polished}"/"${prefix}"_all.bam

    g1=$(basename "$genome")
    g2="${g1%.fasta}"

    if [ "$2" != "map" ]; then
        #Correct contigs using pilon
        java "$memJava" -jar "${prog}"/pilon/pilon-dev.jar \
            --threads "$cpu" \
            --genome "$1" \
            --unpaired "${polished}"/"${prefix}"_merged_sorted.bam \
            --unpaired "${polished}"/"${prefix}"_ccs_sorted.bam \
            --unpaired "${polished}"/"${prefix}"_subreads_Corrected_sorted.bam \
            --frags "${polished}"/"${prefix}"_unmerged_sorted.bam \
            --outdir "$polished" \
            --output "$2" \
            --changes \
            &> >(tee "${logs}"/"${prefix}".pilon"${i}")
    fi
}

# Illumina reads to use for polishing
# mu1=$(find "$merged" -type f -name "*_1P.fastq.gz") #unmerged R1
# mu2=$(sed "s/_1P/_2P/" <<< "$mu1")  #unmerged R2
# m=$(find "$merged" -type f -name "*merged.fastq.gz")  # merged
m="${merged}"/"${prefix}"_merged.fastq.gz
mu1="${merged}"/"${prefix}"_unmerged_1P.fastq.gz
mu2="${merged}"/"${prefix}"_unmerged_2P.fastq.gz

#do up to 10 rounds of polishing
#stop when no more changes are detected
for i in {1..10}; do
    if [ "$i" -gt 1 ]; then
        #check if previous polishing round made any changes
        if [ ! -s  "${polished}"/"${prefix}"_pilon"${i}".changes ]; then
            break  #Exit loop, no more polishing required
        fi
    fi
    
    #Correct contigs using pilon
    Polish "$genome"
done

#Check if polishing was completed after 10 rounds
if [ "$i" -eq 10 ] && [ -s "${polished}"/"${prefix}"_pilon"${i}".changes ]; then
    echo -e "\nGenome can probably be polished more\n" | tee -a "${logs}"/log.txt
fi


# Unicycler hybrid assembly
# output = "${assemblies}"/unicycler/assembly.fasta

#use pacbio ccs as single-end reads
# python3 "${prog}"/Unicycler/unicycler-runner.py \
#     -1 "${corrected}"/"${prefix}"_Cor3_1P.fastq.gz \
#     -2 "${corrected}"/"${prefix}"_Cor3_2P.fastq.gz \
#     -l "${corrected}"/trimmed/"${prefix}".trimmedReads.fasta.gz \
#     -s "${trimmed}"/"${prefix}"_ccs_Cleaned.fastq.gz \
#     -o "${assemblies}"/unicycler \
#     -t "$cpu" \
#     --keep 2 \
#     --no_correct \
#     --mode conservative \
#     --pilon_path "${prog}"/pilon/pilon-dev.jar

#Concatenate merged Illumina paired-end reads with PacBio CCS to feed as single-end reads to Unicycler
cat "$m" "${trimmed}"/"${prefix}"_ccs_Cleaned.fastq.gz \
    > "${corrected}"/"${prefix}"_single.fastq.gz

#use merged paired-end as single-end reads
# --keep 2 seems not to work properly. The alignement files are getting deleted.
python3 "${prog}"/Unicycler/unicycler-runner.py \
    -1 "$mu1" \
    -2 "$mu2" \
    -s "${corrected}"/"${prefix}"_single.fastq.gz \
    -l "${corrected}"/trimmed/"${prefix}".trimmedReads.fasta.gz \
    -o "${assemblies}"/unicycler \
    -t "$cpu" \
    --keep 3 \
    --no_correct \
    --verbosity 2 \
    --mode normal \
    --pilon_path "${prog}"/pilon/pilon-dev.jar

cp "${assemblies}"/unicycler/assembly.fasta \
    "${assemblies}"/"${prefix}"_unicycler.fasta

#Furthur polish
# python3 "${prog}"/Unicycler/unicycler_polish-runner.py \
#     -a "${assemblies}"/unicycler/assembly.fasta \
#     -1 "$mu1" \
#     -2 "$mu2" \2
#     --long_reads "${corrected}"/trimmed/"${prefix}".trimmedReads.fasta.gz \
#     --threads "$cpu" \
#     --pilon "${prog}"/pilon/pilon-dev.jar

run_blast "${assemblies}"/"${prefix}"_unicycler.fasta

#inspect manually blast output
#Remove low coverage or short contigs

#SPAdes hybrid assembly
#Canu assembly is less accurate than Unicycler assembly
spades.py \
    --only-assembler \
    -t "$cpu" \
    -m "$mem" \
    -k "$kmer" \
    --careful \
    --pacbio "${corrected}"/trimmed/"${prefix}".trimmedReads.fasta.gz \
    --trusted-contigs "$genome" \
    --trusted-contigs "${assemblies}"/"${prefix}"_unicycler.fasta \
    --s1 "$m" \
    --pe1-1 "$mu1" \
    --pe1-2 "$mu2" \
    --s2 "${trimmed}"/"${prefix}"_ccs_Cleaned.fastq.gz \
    -o "${assemblies}"/spades

#copy and rename assembly
cp "${assemblies}"/spades/scaffolds.fasta \
    "${assemblies}"/"${prefix}"_spades.fasta

genome="${assemblies}"/"${prefix}"_spades.fasta

run_blast "$genome"

#inspect manually blast output
#Remove low coverage or short contigs

#Polish with pilon
Polish "$genome"


###################
#                 #
#   Circularize   #
#                 #
###################


# Folder to store refseq genomes
[ -d "${ordered}"/refseq ] || mkdir -p "${ordered}"/refseq

# Download all assemblies from species
bash "${scripts}"/get_assemblies.sh \
    -q ""$genus" "$species"" \
    -t 'refseq' \
    -a 'Complete Genome' \
    -o ""${ordered}"/refseq"

# Use mash to find closest genome from refseq
bash "${scripts}"/findClosest.sh \
    "${ordered}"/refseq \
    "$genome" \
    "$cpu"

genome_basename=$(basename "$genome")
genome_name="${genome_basename%.fasta}"
distance=""${ordered}"/refseq/"${genome_name}".distance.tsv"
closest=$(cat "$distance" | head -n 1 | cut -f 1)  # The closest is the top one
closest_acc=$(zcat "$closest" | head -n 1 | cut -d " " -f 1 | tr -d ">")
closest_ID=$(zcat "$closest" | head -n 1 | cut -d " " -f 2- | cut -d "," -f 1)
score=$(cat "$distance" | head -n 1 | cut -f 6)  #its score out of 1000
acc=$(cut -d "_" -f 1,2 <<< $(basename "$closest"))  # accession number of the closest match

echo ""$prefix" closest genome is \""${closest_ID}" ("${closest_acc}")\" with a score of "${score}"/1000" | tee -a "${logs}"/log.txt  #Add information to log

# get all CDS from closest match
esearch -db nucleotide -query "$acc" \
    | efetch -format fasta_cds_na \
    > "${ordered}"/"${acc}"_cds.fasta

#extract only the dnaA entry
cat "${ordered}"/"${acc}"_cds.fasta \
    | awk 'BEGIN {RS = ">"} /dnaA|DnaA|hromosomal replication initiator protein/ {print ">" $0}' \
    > "${ordered}"/"${acc}"_dnaA.fasta

#Check if one and only one dnaA sequence
if [ $(cat "${ordered}"/"${acc}"_dnaA.fasta | grep -Ec "^>") -eq 0 ]; then
    echo "No dnaA gene was found in closest genome from refseq: "$acc"" | tee -a "${logs}"/log.txt
    # exit 1
elif [ $(cat "${ordered}"/"${acc}"_dnaA.fasta | grep -Ec "^>") -gt 1 ]; then
    echo "More than one dnaA gene was found in closest genome from refseq: "$acc"" | tee -a "${logs}"/log.txt
    echo "Keeping th efirst entry only"

    cat "${ordered}"/"${acc}"_dnaA.fasta | \
        awk 'BEGIN {RS = ">"} NR > 1 {print RS $0;exit;}' \
        > "${ordered}"/tmp.txt
    mv "${ordered}"/tmp.txt "${ordered}"/"${acc}"_dnaA.fasta
fi

#Cleanup
rm "${ordered}"/"${acc}"_cds.fasta

#check if one contig or more
if [ $(cat "$genome" | grep -Ec "^>") -eq 1 ]; then
    #just run fixstart
    circlator fixstart \
        --verbose \
        "$genome" \
        "${polished}"/"${prefix}"_circlator
    # circlator fixstart \
    #     --verbose \
    #     --genes_fa "${ordered}"/"${acc}"_dnaA.fasta \
    #     "$genome" \
    #     "${polished}"/"${prefix}"_circlator
else
    #--assemble_spades_k 127 \ #use this option to speed up
    circlator all \
        --verbose \
        --threads "$cpu" \
        --assembler "spades" \
        --assemble_spades_k 127 \
        --data_type "pacbio-corrected" \
        --bwa_opts "-x pacbio" \
        --b2r_discard_unmapped \
        --genes_fa "${ordered}"/"${acc}"_dnaA.fasta \
        "$genome" \
        "${corrected}"/trimmed/"${prefix}".trimmedReads.fasta.gz \
        "${polished}"/circlator

    cp "${polished}"/circlator/06.fixstart.fasta \
        "${polished}"/"${prefix}"_circlator.fasta
fi

genome="${polished}"/"${prefix}"_circlator.fasta
Polish "$genome"


#######################
#                     #
#   Contig ordering   #
#                     #
#######################


function order_contigs ()
{
    # Use Mauve in batch mode to order contigs with closest genome
    java "$memJava" -cp "${prog}"/mauve_snapshot_2015-02-13/Mauve.jar \
        org.gel.mauve.contigs.ContigOrderer \
        -output "${ordered}"/mauve \
        -ref "${closest%.gz}" \
        -draft "$1"

    #fix formating of Mauve output
    #find out how many "alignment" folder there is
    n=$(find "${ordered}"/mauve -maxdepth 1 -type d | sed '1d' | wc -l)

    # reformat with no sorting
    perl "${scripts}"/formatFasta.pl \
        -i "${ordered}"/mauve/alignment"${n}"/"$(basename "$1")" \
        -o "${ordered}"/"${prefix}"_ordered.fasta \
        -w 80

    genome="${ordered}"/"${prefix}"_ordered.fasta

    #align with progessiveMauve
    [ -d "${qc}"/mauve ] || mkdir -p "${qc}"/mauve

    "${prog}"/mauve_snapshot_2015-02-13/linux-x64/./progressiveMauve \
        --output="${qc}"/mauve/"${prefix}"_ordered.xmfa \
        "${closest%.gz}" \
        "$genome"

    # View alignemnt with mauve
    "${prog}"/mauve_snapshot_2015-02-13/./Mauve \
        "${qc}"/mauve/"${prefix}"_ordered.xmfa &

    rm -rf "${ordered}"/mauve
}

[ -s "${closest%.gz}" ] || pigz -d -k "$closest" # decompress if not present

# If more then one contig
if [ $(cat "$genome" | grep -Ec "^>") -gt 1 ];then
    remove_small_contigs "$genome"
    
    #if there is still more than one contig
    if [ $(cat "$genome" | grep -Ec "^>") -gt 1 ];then
        order_contigs "$genome"
    fi
fi

#cleanup
rm -rf "${ordered}"/refseq
find "$ordered" -type f -name "*.sslist" -exec rm {} \;


######################
#                    #
#   Assembly Stats   #
#                    #
######################


# With QUAST
quast.py \
    -o "${qc}"/quast \
    -t "$cpu" \
    "$genome"


# # With BUSCO
# python "${prog}"/busco/scripts/run_BUSCO.py \
#     -i "$genome" \
#     -c "$cpu" \
#     -o "${qc}"/busco \
#     -e 1e-06 \
#     -m "geno" \
#     -l "$busco_db" \
#     -sp "$busco_species" \
#     -t /tmp \
#     -z


###################
#                 #
#   Annotation    #
#                 #
###################


#Prokka
prokka  --outdir "$annotation" \
        --force \
        --prefix "$prefix" \
        --kingdom "$kingdom" \
        --genus "$genus" \
        --species "$species" \
        --strain "$strain" \
        --gram "$gram" \
        --locustag "$locustag" \
        --compliant \
        --centre "$centre" \
        --cpus "$cpu" \
        --rfam \
        "$genome"

#extract hypothetical proteins
cat "${annotation}"/"${prefix}".faa | \
    awk '{if(substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' | \
    grep --no-group-separator -A 1 -F "hypothetical protein" \
    > "${annotation}"/"${prefix}"_hypoth.faa

echo -e "Number of hypothetical proteins found by Prokka: $(cat "${annotation}"/"${prefix}"_hypoth.faa | grep -ic "hypothetical")" \
    | tee -a "${logs}"/log.txt

echo "Blasting hypothetical proteins on NR..."

#make sure $BLASTDB is set in environment variables
# export BLASTDB=/media/3tb_hdd/db/nr:/media/3tb_hdd/db/nt
blastp -query "${annotation}"/"${prefix}"_hypoth.faa \
    -db nr \
    -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
    -evalue 1e-30 \
    -max_target_seqs 1 \
    -max_hsps 1 \
    -num_threads "$cpu" \
    > "${annotation}"/"${prefix}"_hypoth.blastp


# cat "${annotation}"/"${prefix}"_hypoth.faa \
#     | parallel \
#         -j "$cpu" \
#         --block $((mem*1000/48)) \
#         --recstart '>' \
#         --env annotation \
#         --env prefix \
#         --pipe \
#             'blastp \
#                 -query - \
#                 -db nr \
#                 -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
#                 -evalue 1e-30 \
#                 -max_target_seqs 1 \
#                 -max_hsps 1 \
#                 -num_threads 1 \
#                 > "${annotation}"/"${prefix}"_hypoth.blastp'


#Fetch the fasta entry of the hits that do not contain "hypothetical"
#Re-filter for evalues
cat "${annotation}"/"${prefix}"_hypoth.blastp | \
    grep -viF "hypothetical" | \
    awk '{if($12 < 1e-30) {print}}' | \
    cut -f 2 | \
    cut -d "|" -f4 \
    > "${annotation}"/accession.list

#Download the sequences
perl "${scripts}"/http_post.pl \
    "${annotation}"/accession.list \
    "${annotation}"/extra_hits.fasta
# esearch -db protein -query "$(cat "${spadesOut}"/annotation/accession.list)" | \
#     efetch -db protein -format fasta \
#     > "${spadesOut}"/annotation/extra_hits.fasta

#TODO -> Cleanup sequence titles. E.g. "MULTISPECIES: "
# Make first letter uppercase
cat "${annotation}"/extra_hits.fasta | \
    sed 's/MULTISPECIES: //' | \
    sed 's/ ./\U&/'
    > "${annotation}"/extra_hits_renamed.fasta

# remove duplicate entries
cat "${annotation}"/"${sample}"/extra_hits_renamed.fasta \
    | awk 'BEGIN {RS=">"} NR>1 {sub("\n","\t"); gsub("\n",""); print RS$0}' \
    | sort -uk1,1 \
    | tr "\t" "\n" \
    > "${annotation}"/"${sample}"/extra_hits_nodup.fasta

#relaunch Prokka annotation with the new positive blast hit fasta file as reference
prokka  --outdir "$annotation" \
        --force \
        --prefix "$prefix" \
        --kingdom "$kingdom" \
        --genus "$genus" \
        --species "$species" \
        --strain "$strain" \
        --gram "$gram" \
        --locustag "$locustag" \
        --compliant \
        --centre "$centre" \
        --cpus "$cpu" \
        --rfam \
        --proteins "${annotation}"/extra_hits_nodup.fasta \
        "$genome"

echo -e "Number of hypothetical proteins remaining after the BLAST (1e-30): $(cat "${annotation}"/"${prefix}".faa | grep -ic "hypothetical")" \
    | tee -a "${logs}"/log.txt
