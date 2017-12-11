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
prefix="16-389"
genus="Mycobacterium"
species="bovis"
gram="pos"

# Analysis folder
baseDir=""${HOME}"/Desktop/Mbovis_canu/"$prefix""

# Pacbio reads
filtered_subreads="/media/6tb_raid10/data/Mbovis/canada/pacbio/Culture_16_389/Culture_16_389.filtered_subreads.longest.fastq.gz"
ccs_reads="/media/6tb_raid10/data/Mbovis/canada/pacbio/Culture_16_389/Culture_16_389.ccs.fastq.gz"

# Illumina paired-end data
r1="/media/6tb_raid10/data/Mbovis/canada/illumina/2016_outbreak/16-389_R1.fastq.gz"
r2="/media/6tb_raid10/data/Mbovis/canada/illumina/2016_outbreak/16-389_R2.fastq.gz"

# Database to use for metagomic analysis of raw data (contamination)
# db="/media/6tb_raid10/db/centrifuge/p_compressed+h+v"
db="/media/6tb_raid10/db/centrifuge/nt"

# Where programs are installed
prog=""${HOME}"/prog"
scripts=""${HOME}"/scripts"

#genome size (g|m|k)
size="4350k"

# For SPAdes
kmer="21,33,55,77,99,127"

# For assembly trimming
smallest_contig=1000


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
annotation=""${baseDir}"/annotation"

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
    2> >(tee -a "${logs}"/clumpify.txt)

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


# Create folder to store report
[ -d "${qc}"/fastqc/pacbio ] || mkdir -p "${qc}"/fastqc/pacbio
[ -d "${qc}"/fastqc/illumina ] || mkdir -p "${qc}"/fastqc/illumina

echo "Running FastQC on PacBio reads..."

# Run fastQC on pacbio reads
fastqc \
    --o  "${qc}"/fastqc/pacbio \
    --noextract \
    --threads "$cpu" \
    "${fastq}"/"${prefix}"_ccs.fastq.gz \
    "${fastq}"/"${prefix}"_subreads.fastq.gz

echo "Running FastQC on Illumina reads..."

# Run fastQC on fastq reads
fastqc \
    --o  "${qc}"/fastqc/illumina \
    --noextract \
    --threads "$cpu" \
    "${fastq}"/"${prefix}"_R1.fastq.gz \
    "${fastq}"/"${prefix}"_R2.fastq.gz


######################
#                    #
#   Centrifuge raw   #
#                    #
######################


# Create folder
[ -d "${qc}"/centrifuge ] || mkdir -p "${qc}"/centrifuge

echo "Running centrifuge on PacBio reads..."

# Run centrifuge
centrifuge \
    -p "$cpu" \
    -t \
    --seed "$RANDOM" \
    -x "$db" \
    -U "$filtered_subreads" \
    --report-file "${qc}"/centrifuge/"${prefix}"_report_pacbio.tsv \
    > "${qc}"/centrifuge/"${prefix}"_pacbio.tsv

# Prepare result for display with Krona
cat "${qc}"/centrifuge/"${prefix}"_pacbio.tsv | \
    cut -f 1,3 | \
    ktImportTaxonomy /dev/stdin -o  "${qc}"/centrifuge/"${prefix}"_pacbio.html

# Visualize the resutls in Firefow browser
firefox file://"${qc}"/centrifuge/"${prefix}"_pacbio.html &


echo "Running centrifuge on Illumina reads..."

# Run centrifuge
centrifuge \
    -p "$cpu" \
    -t \
    --seed "$RANDOM" \
    -x "$db" \
    -1 "${fastq}"/"${prefix}"_R1.fastq.gz \
    -2 "${fastq}"/"${prefix}"_R2.fastq.gz \
    --report-file "${qc}"/centrifuge/"${prefix}"_report.tsv \
    > "${qc}"/centrifuge/"${prefix}".tsv

# Prepare result for display with Krona
cat "${qc}"/centrifuge/"${prefix}".tsv | \
    cut -f 1,3 | \
    ktImportTaxonomy /dev/stdin -o  "${qc}"/centrifuge/"${prefix}".html

# Visualize the resutls in Firefow browser
firefox file://"${qc}"/centrifuge/"${prefix}".html &


###########################
#                         #
#   Read Pre-Processing   #
#                         #
###########################


# Illumina

#Removing low-quality regions
filterbytile.sh "$memJava" \
    in="${fastq}"/"${prefix}"_R1.fastq.gz \
    in2="${fastq}"/"${prefix}"_R2.fastq.gz \
    out="${trimmed}"/"${prefix}"_Filtered_1P.fastq.gz \
    out2="${trimmed}"/"${prefix}"_Filtered_2P.fastq.gz \
    2> >(tee "${logs}"/illumina_tile_filtering.txt)

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

#Removing synthetic artifacts and spike-ins
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


#PacBio

#CCS trimming
removesmartbell.sh "$memJava" \
    in="$ccs_reads" \
    out="${trimmed}"/"$prefix"_ccs_trimmed.fastq.gz \
    split=t \
    &> >(tee "${logs}"/ccs_trimming.txt)

#subread trimming
removesmartbell.sh "$memJava" \
    in="$filtered_subreads" \
    out="${trimmed}"/"$prefix"_subreads_trimmed.fastq.gz \
    split=t \
    &> >(tee "${logs}"/subreads_trimming.txt)

#Remove internal control from pacbio data
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

#subread hybrid error correction
pigz -dk "${merged}"/"${prefix}"_merged.fastq.gz \
    "${merged}"/"${prefix}"_unmerged_1P.fastq.gz \
    "${merged}"/"${prefix}"_unmerged_2P.fastq.gz

#7H30 to run
perl "${prog}"/proovread/bin/proovread \
    -l "${trimmed}"/"$prefix"_subreads_Cleaned.fastq.gz \
    -s "${merged}"/"${prefix}"_merged.fastq \
    -s "${merged}"/"${prefix}"_unmerged_1P.fastq \
    -s "${merged}"/"${prefix}"_unmerged_2P.fastq \
    -t "$cpu" \
    -p "${corrected}"/"${prefix}"_subreads_Corrected

rm "${merged}"/"${prefix}"_merged.fastq \
    "${merged}"/"${prefix}"_unmerged_1P.fastq \
    "${merged}"/"${prefix}"_unmerged_2P.fastq

mv "${corrected}"/"${prefix}"_subreads_Corrected/"${prefix}"_subreads_Corrected.trimmed.fq \
    "${corrected}"/"${prefix}"_subreads_Corrected.fastq

pigz "${corrected}"/"${prefix}"_subreads_Corrected.fastq

rm -rf "${corrected}"/"${prefix}"_subreads_Corrected


#self-correction?



################
#              #
#   Assembly   #
#              #
################


#de novo assembly
echo "De novo assembly of PacBio reads using canu..."
canu \
    -p "$prefix" \
    -d "$assemblies" \
    genomeSize="$size" \
    -pacbio-corrected "${corrected}"/"${prefix}"_subreads_Corrected.fastq.gz


#######################
#                     #
#      Polishing      #
#                     #
#######################


# Illumina reads are used to polish the PacBio de novo assembly
# They are then used for a new de novo assembly with SPAdes,
# using the polished PacBio assembly as template for scaffolding (--trusted-contigs)
# PacBio longest subreads are also used for scaffolding
# PacBio Circular Consensus Sequences (CCS) are used as single-end reads for the assembly


# Disabled because it introduces too many indels in the assembly
# # Polish assembly with racon
# echo "Mapping PacBio reads onto canu contigs using minimap..."
# minimap \
#     -t "$cpu" \
#     "${assemblies}"/"${prefix}".unitigs.fasta \
#     "${corrected}"/"${prefix}"_subreads_Corrected.fastq.gz \
#     > "${polished}"/"${prefix}"_minimap_overlaps.paf

# echo "Correcting contigs with PacBio reads using racon..."
# racon \
#     -t "$cpu" \
#     "${corrected}"/"${prefix}"_subreads_Corrected.fastq.gz \
#     "${polished}"/"${prefix}"_minimap_overlaps.paf \
#     "${assemblies}"/"${prefix}".unitigs.fasta \
#     "${polished}"/"${prefix}"_racon.fasta


# #convert fastq to fasta
# zcat "${corrected}"/"${prefix}"_subreads_Corrected.fastq.gz \
#     | sed -n '1~4s/^@/>/p;2~4p' \
#     > "${corrected}"/"${prefix}"_subreads_Corrected.fasta

# #map pacbio reads back to assembly
# blasr \
#    "${corrected}"/"${prefix}"_subreads_Corrected.fasta \
#    "${assemblies}"/"${prefix}".unitigs.fasta \
#    --nproc "$cpu" \
#    --out "${polished}"/"${prefix}"_subreads_Corrected.sam \
#    --bam

# samtools view \
#     -@ "$cpu" \
#     -b -h \
#     -T "${assemblies}"/"${prefix}".unitigs.fasta \
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
#     -r "${assemblies}"/"${prefix}".unitigs.fasta \
#     -o "${polished}"/"${prefix}"_quiver.fasta \
#     "${polished}"/"${prefix}"_subreads_sorted.bam

# source deactivate


function Polish()
{
    # Index genome
    bwa index "$genome"

    #Align merged
    bwa mem -t "$cpu" -M "$1" "$m" | \
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
    bwa mem -t "$cpu" -M "$1" "$mu1" "$mu2" | \
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
    bwa mem -x pacbio -t "$cpu" -M "$1" "${trimmed}"/"${prefix}"_ccs_Cleaned.fastq.gz | \
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

    #Align pacbio corrected longest subreads
    bwa mem -x pacbio -t "$cpu" -M "$1" "${corrected}"/"${prefix}"_subreads_Corrected.fastq.gz | \
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

    if [ "$3" != "map" ]; then
        #Correct contigs using pilon based on the Illumina reads
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
            &> >(tee "${logs}"/"${g2}".pilon)
    fi
}

# Illumina reads to use for polishing
mu1=$(find "$merged" -type f -name "*_1P.fastq.gz") #unmerged R1
mu2=$(sed "s/_1P/_2P/" <<< "$mu1")  #unmerged R2
m=$(find "$merged" -type f -name "*merged.fastq.gz")  # merged

# 1st round of polishing
genome="${assemblies}"/"${prefix}".unitigs.fasta  # Canu assembly
Polish "$genome" "${prefix}"_pilon1
genome="${polished}"/"${prefix}"_pilon1.fasta

#SPAdes hybrid assembly
spades.py \
    --only-assembler \
    -t "$cpu" \
    -m "$mem" \
    -k "$kmer" \
    --careful \
    --pacbio "${corrected}"/"${prefix}"_subreads_Corrected.fastq.gz \
    --untrusted-contigs "$genome" \
    --s1 "$m" \
    --pe1-1 "$mu1" \
    --pe1-2 "$mu2" \
    --s2 "${trimmed}"/"${prefix}"_ccs_Cleaned.fastq.gz \
    -o "${polished}"/spades

#copy and rename assembly
cp "${polished}"/spades/scaffolds.fasta "${polished}"/"${prefix}"_spades.fasta

#Polish with pilon
#Second round of Pilon correction
genome="${polished}"/"${prefix}"_spades.fasta
Polish "$genome" "${prefix}"_pilon2
genome="${polished}"/"${prefix}"_pilon2.fasta


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

distance=""${ordered}"/refseq/"${prefix}"_pilon2.distance.tsv"
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
    | awk 'BEGIN {RS = ">"} /dnaA/ {print ">" $0}' \
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
    "${corrected}"/"${prefix}"_subreads_Corrected.fastq.gz \
    "${polished}"/circlator

cp "${polished}"/circlator/06.fixstart.fasta \
    "${polished}"/"${prefix}"_circlator.fasta

genome="${polished}"/"${prefix}"_circlator.fasta
Polish "$genome" "${prefix}"_pilon3
genome="${polished}"/"${prefix}"_pilon3.fasta


#######################
#                     #
#   Contig ordering   #
#                     #
#######################


function remove_small_contigs ()
{
    # Remove contigs smaller than X bp
    perl "${prog}"/phage_typing/removesmallscontigs.pl \
        "$smallest_contig" \
        "$genome" \
        > "${genome%.fasta}"_"${smallest_contig}".fasta

    genome="${genome%.fasta}"_"${smallest_contig}".fasta
}

function order_contigs ()
{
    # Use Mauve in batch mode to order contigs with closest genome
    java "$memJava" -cp "${prog}"/mauve_snapshot_2015-02-13/Mauve.jar \
        org.gel.mauve.contigs.ContigOrderer \
        -output "${ordered}"/mauve \
        -ref "${closest%.gz}" \
        -draft "$genome"

    #fix formating of Mauve output
    #find out how many "alignment" folder there is
    n=$(find "${ordered}"/mauve -maxdepth 1 -type d | sed '1d' | wc -l)

    # reformat with no sorting
    perl "${scripts}"/formatFasta.pl \
        -i ""${ordered}"/mauve/alignment"${n}"/"${prefix}"_pilon3.fasta" \
        -o "${ordered}"/"${prefix}"_ordered.fasta \
        -w 80

    genome="${ordered}"/"${prefix}"_ordered.fasta
}

[ -s "${closest%.gz}" ] || pigz -d -k "$closest" # decompress if not present

# If more then one contig
if [ $(cat "$genome" | grep -Ec "^>") -gt 1 ];then
    remove_small_contigs
    
    #if there is still more than one contig
    if [ $(cat "$genome" | grep -Ec "^>") -gt 1 ];then
        order_contigs
    fi
fi

#align with progessiveMauve
[ -d "${qc}"/mauve ] || mkdir -p "${qc}"/mauve

"${prog}"/mauve_snapshot_2015-02-13/linux-x64/./progressiveMauve \
    --output="${qc}"/mauve/"${prefix}"_ordered.xmfa \
    "${closest%.gz}" \
    "$genome"

# View alignemnt with mauve
"${prog}"/mauve_snapshot_2015-02-13/./Mauve \
    "${qc}"/mauve/"${prefix}"_ordered.xmfa &

#cleanup
rm -rf "${ordered}"/mauve
rm -rf "${ordered}"/refseq
find "$ordered" -type f -name "*.sslist" -exec rm {} \;


#################
#               #
#   Contig ID   #
#               #
#################


blastn \
    -db nt \
    -query "$genome" \
    -out "${genome%.fasta}".blastn.tsv \
    -evalue "1e-10" \
    -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
    -num_threads "$cpu" \
    -max_target_seqs 1

echo -e "qseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" \
    > "${polished}"/tmp.txt
cat "${genome%.fasta}".blastn.tsv >> "${polished}"/tmp.txt
mv "${polished}"/tmp.txt "${genome%.fasta}".blastn.tsv


######################
#                    #
#   Assembly Stats   #
#                    #
######################


quast.py \
    -o "${qc}"/quast \
    -t "$cpu" \
    "$genome"


###################
#                 #
#   Annotation    #
#                 #
###################


#Prokka
prokka  --outdir "$annotation" \
        --force \
        --prefix "$prefix" \
        --centre "OLF" \
        --kingdom "Bacteria" \
        --genus "$genus" \
        --species "$species" \
        --strain "$prefix" \
        --gram "$gram" \
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
    -num_threads "$cpu" \
    > "${annotation}"/"${prefix}"_hypoth.blastp

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

#relaunch Prokka annotation with the new positive blast hit fasta file as reference
prokka  --outdir "$annotation" \
        --force \
        --prefix "$prefix" \
        --centre "OLF" \
        --kingdom "Bacteria" \
        --genus "$genus" \
        --species "$species" \
        --strain "$prefix" \
        --gram "neg" \
        --cpus "$cpu" \
        --rfam \
        --proteins "${annotation}"/extra_hits.fasta \
        "$genome"

echo -e "Number of hypothetical proteins remaining after the BLAST (1e-30): $(cat "${annotation}"/"${prefix}".faa | grep -ic "hypothetical")" \
    | tee -a "${logs}"/log.txt
