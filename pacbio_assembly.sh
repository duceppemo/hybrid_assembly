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
prefix="CFF_09A980"
genus="Campylobacter"
species="fetus"

# Analysis folder
baseDir=""${HOME}"/analyses/pacbio_test/"$prefix""

# Pacbio raw image files location
# Point to where the ".metada.xml" file is located
# A subdirectory named "Analysis_Results" containing the associated ".bas.h5" and ".bax.h5" file must be present
# raw_data="/home/smrtanalysis/pacbio_data/CFV/08A948"  # for hgap
filtered_subreads="/media/6tb_raid10/data/campylobacter/raw_reads/olf/pacbio/CFF/CFF_09A980.filtered_subreads.fastq.gz"

# Illumina paired-end data
r1="/media/6tb_raid10/data/campylobacter/raw_reads/olf/illumina/CFF/09A980_combined_R1.fastq.gz"
r2="/media/6tb_raid10/data/campylobacter/raw_reads/olf/illumina/CFF/09A980_combined_R2.fastq.gz"

# Database to use for metagomic analysis of raw data (contamination)
# db="/media/6tb_raid10/db/centrifuge/p_compressed+h+v"
db="/media/6tb_raid10/db/centrifuge/nt"

# Where programs are installed
prog=""${HOME}"/prog"
scripts=""${HOME}"/scripts"
picard=""${prog}"/picard-tools/picard.jar"

# Maximum number of cores used per sample for parallel processing
# A highier value reduces the memory footprint.
# maxProc=6

#genome size (g|m|k)
size="2m"


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
merged=""${baseDir}"/merged"
aligned=""${baseDir}"/aligned"
annotation=""${baseDir}"/annotation"

# Create folders if do not exist
# "||" if test is false
# "&&" if test is true
[ -d "$baseDir" ] || mkdir -p "$baseDir"
[ -d "$logs" ] || mkdir -p "$logs"
[ -d "$qc" ] || mkdir -p "$qc"
[ -d "$fastq" ] || mkdir -p "$fastq"
[ -d "$assemblies" ] || mkdir -p "$assemblies"
[ -d "$polished" ] || mkdir -p "$polished"
[ -d "$trimmed" ] || mkdir -p "$trimmed"
[ -d "$merged" ] || mkdir -p "$merged"
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


# Check Centrifuge datase
if [ -s "${db}".1.cf ]; then
    echo -e "\nCentrifuge database: $(basename "$db")" | tee -a "${logs}"/log.txt
else
    echo "Could not find the centrifuge database"
    exit 1
fi

#Picard tool

#Circlator

#Prokka

#Mauve


#############
#           #
#   FASTQ   #
#           #
#############


# Populate fastq files
# TODO


############
#          #
#   HGAP   #
#          #
############


# Run HGAP pipeline from the command line (instead of useing the SMRT portal)
# TODO


####################
#                  #
#   FastQc - Raw   #
#                  #
####################


# Create folder to store report
[ -d "${qc}"/fastqc/pacbio ] || mkdir -p "${qc}"/fastqc/pacbio
[ -d "${qc}"/fastqc/fastq ] || mkdir -p "${qc}"/fastqc/fastq

echo "Running FastQC on PacBio reads..."

# Run fastQC on pacbio reads
fastqc \
    --o  "${qc}"/fastqc/pacbio \
    --noextract \
    --threads "$cpu" \
    "$filtered_subreads"

echo "Running FastQC on Illumina reads..."

# Run fastQC on fastq reads
fastqc \
    --o  "${qc}"/fastqc/fastq \
    --noextract \
    --threads "$cpu" \
    "$r1" \
    "$r2"


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
    -1 "$r1" \
    -2 "$r2" \
    --report-file "${qc}"/centrifuge/"${prefix}"_report.tsv \
    > "${qc}"/centrifuge/"${prefix}".tsv

# Prepare result for display with Krona
cat "${qc}"/centrifuge/"${prefix}".tsv | \
    cut -f 1,3 | \
    ktImportTaxonomy /dev/stdin -o  "${qc}"/centrifuge/"${prefix}".html

# Visualize the resutls in Firefow browser
firefox file://"${qc}"/centrifuge/"${prefix}".html &


########################
#                      #
#   de novo assembly   #
#                      #
########################


echo "De novo assembly of PacBio reads using canu..."

# canu
# assembly file name is "${prefix}".unitigs.fasta
canu \
    -p "$prefix" \
    -d "$assemblies" \
    genomeSize="$size" \
    -pacbio-raw "$filtered_subreads"

# Polish assembly with racon

echo "Mapping pacbio reads onto canu contigs using minimap..."

minimap \
    -t "$cpu" \
    "${assemblies}"/"${prefix}".unitigs.fasta \
    "$filtered_subreads" \
    > "${polished}"/"${prefix}"_minimap_overlaps.paf


echo "Correcting contigs using racon..."
racon \
    -t "$cpu" \
    "$filtered_subreads" \
    "${polished}"/"${prefix}"_minimap_overlaps.paf \
    "${assemblies}"/"${prefix}".unitigs.fasta \
    "${polished}"/"${prefix}"_racon.fasta

genome=""${polished}"/"${prefix}"_racon.fasta"

# Index contigs, map Illumina reads to contigs by BWA, sorting and indexing using samtools:

echo "Indexing genome for bwa..."

# Index genome assembly polished by racon
bwa index "$genome"


####################
#                  #
#     Trimming     #
#                  #
####################


echo -e "Quality trimming Illumina paired-end reads..."

[ -e "${logs}"/trimming.txt ] && rm "${logs}"/trimming.txt

#Trimming
bbduk.sh "$memJava" \
    in1="$r1" \
    in2="$r2" \
    ref="${prog}"/bbmap/resources/nextera.fa.gz \
    ktrim=r k=23 mink=11 hdist=1 tbo tpe \
    qtrim=lr trimq=10 \
    minlen=64 \
    out1="${trimmed}"/"${prefix}"_Trimmed_1P.fastq.gz \
    out2="${trimmed}"/"${prefix}"_Trimmed_2P.fastq.gz \
    pigz=t \
    unpigz=t \
    2> >(tee -a "${logs}"/trimming.txt)


###################
#                 #
#     Merging     #
#                 #
###################


echo -e "Merging Illumina paried-end reads..."

t1=$(find "$trimmed" -type f | grep -F "fastq.gz" | grep -F "_1P")
t2=$(sed "s/_1P/_2P/" <<< "$t1")

[ -e "${logs}"/merging.txt ] && rm "${logs}"/merging.txt

bbmerge.sh "$memJava" \
    in1="$t1" \
    in2="$t2" \
    out="${merged}"/"${prefix}"_merged.fastq.gz \
    outu1="${merged}"/"${prefix}"_unmerged_1P.fastq.gz \
    outu2="${merged}"/"${prefix}"_unmerged_2P.fastq.gz \
    pigz=t \
    ziplevel=5 \
    unpigz=t \
    2> >(tee -a "${logs}"/merging.txt)


#######################
#                     #
#      Alignment      #
#                     #
####################### 


echo "Mapping Illumina reads on PacBio corrected assembly..."

#unmerged
mu1=$(find "$merged" -type f | grep -F "fastq.gz" | grep -F "_1P")
mu2=$(sed "s/_1P/_2P/" <<< "$mu1")

#merged
m=$(find "$merged" -type f | grep -F "fastq.gz" | grep -F "_merged")

#Align merged
bwa mem -x intractg -t "$cpu" -r 1 -a -M "$genome" "$m" | \
    samtools view -@ "$cpu" -b -h -F 4 - | \
    samtools sort -@ "$cpu" -m 10G -o "${polished}"/"${prefix}"_merged.bam -

# remove duplicates for sinle-end reads
java "$memJava" -jar "$picard" MarkDuplicates \
    INPUT="${polished}"/"${prefix}"_merged.bam \
    OUTPUT="${polished}"/"${prefix}"_merged_nodup.bam \
    METRICS_FILE="${polished}"/"${prefix}"_merged_duplicates.txt \
    ASSUME_SORTED=true \
    REMOVE_DUPLICATES=true

# Merge and sort all merged files
samtools sort -@ "$cpu" -m 10G \
    -o "${polished}"/"${prefix}"_merged_sorted.bam \
    "${polished}"/"${prefix}"_merged_nodup.bam

# index bam file
samtools index "${polished}"/"${prefix}"_merged_sorted.bam

#Align unmerged
bwa mem -x intractg -t "$cpu" -r 1 -a -M "${polished}"/"${prefix}"_racon.fasta "$mu1" "$mu2" | \
    samtools view -@ "$cpu" -b -h -F 4 - | \
    samtools sort -@ "$cpu"  -m 10G -o "${polished}"/"${prefix}"_unmerged.bam -

echo "Removing duplicates in bam file..."

samtools rmdup \
    "${polished}"/"${prefix}"_unmerged.bam \
    "${polished}"/"${prefix}"_unmerged_nodup.bam

# Merge and sort all merged files
samtools sort -@ "$cpu" -m 10G \
    -o "${polished}"/"${prefix}"_unmerged_sorted.bam \
    "${polished}"/"${prefix}"_unmerged_nodup.bam

samtools index "${polished}"/"${prefix}"_unmerged_sorted.bam

#Correct contigs using pilon based on the Illumina reads
# java "$memJava" -XX:+UseConcMarkSweepGC -XX:-UseGCOverheadLimit \
java "$memJava" -jar "${prog}"/pilon/pilon-dev.jar \
    --threads "$cpu" \
    --genome "${polished}"/"${prefix}"_racon.fasta \
    --unpaired "${polished}"/"${prefix}"_merged_sorted.bam \
    --jumps "${polished}"/"${prefix}"_unmerged_sorted.bam \
    --outdir "$polished" \
    --output "${prefix}"_pilon \
    --changes

#Cleanup
find "$polished" -type f ! -name "*pilon*" -exec rm {} \;


#Blast on nt
blastn \
    -db nt \
    -query "${polished}"/"${prefix}"_pilon.fasta \
    -out "${baseDir}"/"${prefix}".blastn \
    -evalue "1e-30" \
    -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
    -num_threads "$cpu" \
    -max_target_seqs 1

#add header to blast file
echo -e "qseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" \
    > "${baseDir}"/tmp.txt
cat "${baseDir}"/"${prefix}".blastn >> "${baseDir}"/tmp.txt
mv "${baseDir}"/tmp.txt "${baseDir}"/"${prefix}".blastn


#map pacbio reads back to assembly
bwa index "$polished"/"${prefix}"_pilon.fasta

bwa mem -x pacbio -t "$cpu" -r 1 -a -M "$polished"/"${prefix}"_pilon.fasta "$filtered_subreads" | \
    samtools view -@ "$cpu" -b -h - | \
    samtools sort -@ "$cpu" -m 10G -o "${aligned}"/"${prefix}"_pacbio_merged.bam -

# remove duplicates for sinle-end reads
java "$memJava" -jar "$picard" MarkDuplicates \
    INPUT="${aligned}"/"${prefix}"_pacbio_merged.bam \
    OUTPUT="${aligned}"/"${prefix}"_pacbio_merged_nodup.bam \
    METRICS_FILE="${aligned}"/"${prefix}"_pacbio_merged_duplicates.txt \
    ASSUME_SORTED=true \
    REMOVE_DUPLICATES=true

#remove bam with duplicates
rm "${aligned}"/"${prefix}"_pacbio_merged.bam

#index bam file
samtools index "${aligned}"/"${prefix}"_pacbio_merged_nodup.bam

#Visualize mapping in IGV (or genome browser of choice)


#######################
#                     #
#   Contig ordering   #
#                     #
#######################


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
    "$polished"/"${prefix}"_pilon.fasta \
    "$cpu"

# The closest is the top one
closest=$(cat "${ordered}"/refseq/"${prefix}"_pilon"${smallest_contig}".distance.tsv | head -n 1 | cut -f 1)

#its score out of 1000
score=$(cat "${ordered}"/refseq/"${prefix}"_pilon"${smallest_contig}".distance.tsv | head -n 1 | cut -f 6)

#Add information to log
echo ""$prefix" closest genome is "$(basename "${closest%.*}")" with a score of "$score"/1000" | tee -a "${logs}"/log.txt

# Use Mauve in batch mode to order contigs with closest genome
[ -s "${closest%.gz}" ] || pigz -d -k "$closest"  # decompress if not present

#Run both mauve and circlator?

if [ $(cat "$polished"/"${prefix}"_pilon.fasta | grep -Ec "^>") -gt 1 ] && [[ -z $(cat "${baseDir}"/"${prefix}".blastn | grep -i "plasmid") ]]; then
    java "$memJava" -cp "${prog}"/mauve_snapshot_2015-02-13/Mauve.jar \
        org.gel.mauve.contigs.ContigOrderer \
        -output "${ordered}"/mauve \
        -ref "${closest%.gz}" \
        -draft "$polished"/"${prefix}"_pilon.fasta

    #fix formating of Mauve output
    #find out how many "alignment" folder there is
    n=$(find "${ordered}"/mauve -maxdepth 1 -type d | sed '1d' | wc -l)

    # reformat with no sorting
    perl "${scripts}"/formatFasta.pl \
        -i "${ordered}"/mauve/alignment"${n}"/"${prefix}"_pilon.fasta \
        -o "${ordered}"/"${prefix}"_ordered.fasta \
        -w 80

    [ -d "${qc}"/mauve ] || mkdir -p "${qc}"/mauve

    #align with progessiveMauve
    "${prog}"/mauve_snapshot_2015-02-13/linux-x64/./progressiveMauve \
        --output="${qc}"/mauve/"${prefix}"_ordered.xmfa \
        "${closest%.gz}" \
        "${ordered}"/"${prefix}"_ordered.fasta

    # View alignemnt with mauve
    # "${prog}"/mauve_snapshot_2015-02-13/./Mauve

else
    # accession number of the closest match
    acc=$(cut -d "_" -f 1,2 <<< $(basename "$closest"))
    
    # get all CDS from closest match
    esearch -db nucleotide -query "$acc" \
        | efetch -format fasta_cds_na \
        > "${ordered}"/"${acc}"_cds.fasta

    #extract only the dnaA entry
    cat "${ordered}"/"${acc}"_cds.fasta \
        | awk 'BEGIN {RS = ">"} /dnaA/ {print ">" $0}' \
        > "${ordered}"/"${acc}"_dnaA.fasta

    rm "${ordered}"/"${acc}"_cds.fasta

    #circlator
    circlator fixstart \
        --genes_fa "${ordered}"/"${acc}"_dnaA.fasta \
        "${polished}"/"${prefix}"_pilon.fasta \
        "${ordered}"/"${prefix}"_ordered
fi

#cleanup
rm -rf "${ordered}"/mauve
rm -rf "${ordered}"/refseq
find "$ordered" -type f -name "*.sslist" -exec rm {} \;


#############
#           #
#   QUAST   #
#           #
#############


quast.py \
    -o "${qc}"/quast \
    -t "$cpu" \
    "${ordered}"/"${prefix}"_ordered.fasta


###################
#                 #
#   Annotation   #
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
        --gram "neg" \
        --cpus "$cpu" \
        --rfam \
        "${ordered}"/"${prefix}"_ordered.fasta

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
        "${ordered}"/"${prefix}"_ordered.fasta

echo -e "Number of hypothetical proteins remaining after the BLAST (1e-30): $(cat "${annotation}"/"${prefix}".faa | grep -ic "hypothetical")" \
    | tee -a "${logs}"/log.txt
