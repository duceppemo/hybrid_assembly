#!/bin/bash


##############################
#
#   Get assemblies or refseq sequences from GenBank
#   using their command line tool "edirect" (uses E-utilities)
#
#   more info on edirect: https://www.ncbi.nlm.nih.gov/books/NBK179288/
#   more info on E-utilities: https://www.ncbi.nlm.nih.gov/books/NBK25500/
#
##############################


function print_help() {
    echo "\
Usage: bash get_assemblies.sh [option]

Mandatory flags:

    -q | -l     Query text (e.g. "Mycobacterim bovis") OR list file with accession numbers

    -t          Type of sequece to retrieve ("refseq" or "assemblies")

    -o          Output folder

Optional flags:

    -h          Print this help

    -r          Rename files according to sequence header
            
    -a          Assembly level restriction (Chromosome, Complete Genome, Contig, Scaffold)
                Choose only one. Not using this option will download all assembly levels

    -u          Uncompress all downloaded files\
"
exit 1
}


#http://wiki.bash-hackers.org/howto/getopts_tutorial
# Default values
query=""
out=""
db=""
list=""
rename=0
level="all"
uncompress=0

options=':q:o:t:rul:a:h'

while getopts "$options" opt; do
    case "$opt" in
        q) export query="$OPTARG";;
        t) export db="$OPTARG";;
        o) export out="$OPTARG";;
        l) export list="$OPTARG";;
        a) level="$OPTARG";;
        r) rename=1;;
        u) uncompress=1;;
        h) print_help;;
        \?) echo "Invalid option: -"$OPTARG"" >&2; print_help;;
        :) echo "Option -"$OPTARG" requires an argument." >&2; print_help;;
    esac
done

# Get frags and arguments
shift $((OPTIND - 1))

# Exit if argument is missing
if [[ ( -z "$list" && -z "$query") ]] || [ -z "$out" ] || [ -z "$db" ]; then 
    echo "All \"-q\", \"-t\" and \"-o\" options are mandatory"
    echo "Query: "$query", List: "$list", Out: "$out", DB: "$db", Level: "$level""
    print_help
fi

# if [ -z "$fasta" ] && [ -z "$fastq" ] || [ -z "$primers" ] || [ -z "$output" ]; then 
#     echo "All \"-f or -q\", \"-p\" and \"-o\" options are mandatory"
#     print_help
# fi

# Check if only one of "-q" or "-l" is being used
if [ -n "$query" ] && [ -n "$list" ]; then
    echo "\
You can either submit a query term with the \"-q\" option flag
or submit a list of accessions with the \"-l\" option flag"
    print_help
fi

# Check if valid list file
if [ -n "$list" ]; then
    if [ ! -s "$list" ]; then
        echo "The list file provided is no valid"
        print_help
    fi
fi

# Check if good database type entered
if [ "$db" == "refseq" ] || [ "$db" == "assemblies" ]; then
    echo "Downloading "$query" "$level" sequences from "$db"..."
else
    echo "\"-t\" option is madatory and should be \"refseq\" or \"assemblies\""
    print_help
fi

# <FtpPath_GenBank>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/009/445/GCA_000009445.1_ASM944v1</FtpPath_GenBank>
# <FtpPath_RefSeq>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/445/GCF_000009445.1_ASM944v1</FtpPath_RefSeq>
# <FtpPath_Assembly_rpt>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/445/GCF_000009445.1_ASM944v1/GCF_000009445.1_ASM944v1_assembly_report.txt</FtpPath_Assembly_rpt>

# Set database type
if [ "$db" == "refseq" ]; then
    export db="FtpPath_RefSeq"
else  # if [ "$db" == "assemblies" ]; then
    export db="FtpPath_GenBank"
fi

#Create Directory if does not exist
[ -d "$out" ] || mkdir -p "$out"

# Filtering
export assembly_level=""

#check if valid level used
if [ -n "$level" ]; then
    case "$level" in  # indirect variable expansion
        Chromosome)         export assembly_level="$level";;
        'Complete Genome')  export assembly_level="$level";;
        Contig)             export assembly_level="$level";;
        Scaffold)           export assembly_level="$level";;
        all)                export assembly_level="$level";;
    esac

    if [ -z "$assembly_level" ]; then
        echo "Invalid assembly level term \""$level"\""
        print_help
    fi
fi

# <AssemblyStatus>Complete Genome</AssemblyStatus>
# <AssemblyStatus>Contig</AssemblyStatus>
# <AssemblyStatus>Scaffold</AssemblyStatus>
# <AssemblyStatus>Chromosome</AssemblyStatus>

if [ "$assembly_level" != "all" ]; then  #Get only selected assembly levels
    function search_all ()
    {
        esearch -db assembly -query "$1" \
            | efetch -format docsum \
            | xtract -pattern DocumentSummary -element "$db" AssemblyStatus \
            | grep -F "$assembly_level" \
            | cut -f 1 \
            | awk -F"/" '{print $0"/"$NF"_genomic.fna.gz"}'
    }

    if [ -z "$list" ]; then
        search_all "$query" \
            | parallel \
                --bar \
                --env out \
                --delay 0.3 \
                'wget -P $out {} 2> /dev/null'
    else
        # make function available to parallel
        export -f search_all

        cat "$list" \
            | parallel \
                --delay 0.3 \
                --env search_all \
                --env db \
                --env assembly_level \
                'search_all {}' \
            | parallel \
                --bar \
                --env out \
                --delay 0.3 \
                'wget -P $out {} 2> /dev/null'
    fi
else  # Get all sequences
    function search_some ()
    {
        esearch -db assembly -query "$1" \
            | efetch -format docsum \
            | xtract -pattern DocumentSummary -element "$db" \
            | awk -F"/" '{print $0"/"$NF"_genomic.fna.gz"}'
    }

    if [ -z "$list" ]; then
        search_some "$query" \
            | parallel \
                --bar \
                --env out \
                --delay 0.3 \
                'wget -P $out {} 2> /dev/null'
    else
        # make function available to parallel
        export -f search_some

        cat "$list" \
            | parallel \
                --delay 0.3 \
                --env search_some \
                --env db \
                --env assembly_level \
                'search_some {}' \
            | parallel \
                --bar \
                --env out \
                --delay 0.3 \
                'wget -P $out {} 2> /dev/null'
    fi
            
fi

# Rename files according to fasta header
if [ "$rename" -eq 1 ]; then
    echo "Renaming downloaded files according to fasta header..."

    function rename()
    {
        path=$(dirname "$1")

        new_name=$(zcat "$1" | grep -E "^>" | head -n 1 | cut -d " " -f 2- | cut -d "," -f 1 \
            | tr " " "_" | tr "/" "_" | tr -d "(" | tr -d ")" | tr -d "." \
            | sed   -e 's/contig.*//' \
                    -e 's/Contig.*//' \
                    -e 's/cont.*//' \
                    -e 's/genomic.*//' \
                    -e 's/scaffold.*//' \
                    -e 's/chrom.*//' \
                    -e 's/Chrom.*//' \
                    -e 's/_$//' \
                    -e 's/_=.*//' \
                    -e 's/_NODE.*//' \
                    -e 's/_genome_assembly//' \
                    -e 's/_complete_genome//')

        mv "$1" "${path}"/"${new_name}".fna.gz
    }

    export -f rename

    find "$out" -type f -name "*_genomic.fna.gz" \
        | parallel --env rename 'rename {}'
fi

# | sed   -e 's/genome assembly//' \
#                 -e 's/genomic scaffold//'

#Rename assemblies
# for i in $(find "$out" -type f -name "*_genomic.fna.gz"); do
#     a=$(zcat "$i" | grep -E "^>" | head -n 1 | cut -d " " -f 1 | tr -d ">")

#     if [[ -n $(echo "$a" | grep -F "01000001.1") ]]; then
#         mv "$i" "${out}"/"${a:0:4}00000000.fna.gz"
#     else
#         mv "$i" "${out}"/"${a}".fna.gz
#     fi
# done

#uncompress
if [ "$uncompress" -eq 1 ]; then
    echo "Decompressing files..."
    find "$out" -type f -name "*fna.gz" \
        | parallel --bar 'pigz -d {}'
fi


