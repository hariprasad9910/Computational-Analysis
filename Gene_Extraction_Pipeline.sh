# Gene Extraction Pipeline

This repository contains a comprehensive pipeline for downloading SRA data, processing FASTQ files, performing read mapping, and extracting specific gene sequences (like rbcL) from the reference genome.

## Pipeline Overview

The pipeline performs the following steps:
1. Download SRA data for rice samples from different IRGC accessions
2. Convert SRA to FASTQ format and split paired-end reads
3. Perform quality control on the FASTQ files
4. Map reads to the reference genome using STAR
5. Extract specific genes (e.g., rbcL) from the reference genome

## Prerequisites

- SRA Toolkit (prefetch, fastq-dump)
- FastQC
- STAR aligner
- Samtools
- BEDTools
- Bash

# Make the script executable
chmod +x gene_extraction.sh
```

## Script Contents

```bash
#!/bin/bash

# Gene Extraction Pipeline
# This script downloads SRA data, processes it to FASTQ, performs QC,
# maps reads to a reference genome, and extracts specific genes

# Configuration - Edit these variables as needed
BASE_DIR="/home/hari/gene_extraction"
GENOME_DIR="${BASE_DIR}/mapping/cpar"
GTF_FILE="${BASE_DIR}/mapping/GCF_034140825.1_ASM3414082v1_genomic.gtf"
GFF_FILE="${BASE_DIR}/mapping/GCF_034140825.1_ASM3414082v1_genomic.gff"
GENOME_FASTA="${BASE_DIR}/mapping/GCF_034140825.1_ASM3414082v1_genomic.fna"
OUTPUT_BASE_DIR="${BASE_DIR}/alignment_results"
FASTQC_OUTPUT="${BASE_DIR}/fastqc"
THREADS=2

# Create necessary directories
mkdir -p "${BASE_DIR}"
mkdir -p "${OUTPUT_BASE_DIR}"
mkdir -p "${FASTQC_OUTPUT}"

# Sample information
declare -A samples=(
    ["IRGC_127150"]="ERR605303" 
    ["IRGC_127340"]="ERR612137" 
    ["IRGC_127205"]="ERR635084" 
    ["IRGC_127135"]="ERR622382"
)

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check dependencies
check_dependencies() {
    local missing_deps=()
    
    for cmd in prefetch fastq-dump fastqc STAR samtools bedtools; do
        if ! command_exists "$cmd"; then
            missing_deps+=("$cmd")
        fi
    done
    
    if [ ${#missing_deps[@]} -gt 0 ]; then
        echo "❌ Error: Missing dependencies: ${missing_deps[*]}"
        echo "Please install the required tools before running this script."
        exit 1
    fi
    
    echo "✅ All dependencies are installed."
}

# Step 1: Download SRA data
download_sra_data() {
    echo "📥 Downloading SRA data..."
    
    for accession in "${!samples[@]}"; do
        sra_id="${samples[$accession]}"
        accession_dir="${BASE_DIR}/${accession}"
        mkdir -p "$accession_dir"
        
        # Remove any lock files if they exist
        rm -f "${accession_dir}/${sra_id}/${sra_id}.sra.lock" 2>/dev/null
        
        echo "  • Downloading ${sra_id} for ${accession}..."
        prefetch "$sra_id" -O "$accession_dir"
    done
    
    echo "✅ SRA data download complete."
}

# Step 2: Convert SRA to FASTQ and split paired-end reads
convert_to_fastq() {
    echo "🔄 Converting SRA to FASTQ format..."
    
    for accession in "${!samples[@]}"; do
        sra_id="${samples[$accession]}"
        accession_dir="${BASE_DIR}/${accession}"
        sra_file="${accession_dir}/${sra_id}/${sra_id}.sra"
        
        if [ -f "$sra_file" ]; then
            echo "  • Processing ${sra_id} for ${accession}..."
            fastq-dump --split-files --gzip "$sra_file" -O "$accession_dir"
        else
            echo "  ❌ SRA file not found: $sra_file"
        fi
    done
    
    echo "✅ FASTQ conversion complete."
}

# Step 3: Check quality of reads
check_read_quality() {
    echo "🔍 Checking read quality with FastQC..."
    
    for accession in "${!samples[@]}"; do
        accession_dir="${BASE_DIR}/${accession}"
        
        if ls "${accession_dir}"/*fastq.gz 1> /dev/null 2>&1; then
            echo "  • Running FastQC on ${accession} reads..."
            fastqc "${accession_dir}"/*fastq.gz -o "${FASTQC_OUTPUT}" -t "$THREADS"
        else
            echo "  ❌ No FASTQ files found in: $accession_dir"
        fi
    done
    
    echo "✅ FastQC analysis complete. Results saved in ${FASTQC_OUTPUT}"
}

# Step 4: Generate sample IDs list
generate_sample_ids() {
    echo "📋 Generating sample IDs list..."
    
    find "${BASE_DIR}" -type f -name "*.gz" | awk -F'[_/]' '{print $(NF-1)}' | sort | uniq > "${BASE_DIR}/sample_ids.txt"
    
    echo "✅ Sample IDs list generated at ${BASE_DIR}/sample_ids.txt"
}

# Step 5: Read mapping with STAR
map_reads() {
    echo "🧬 Mapping reads to reference genome using STAR..."
    
    # Create output directory
    mkdir -p "$OUTPUT_BASE_DIR"
    
    # Check if genome directory exists
    if [ ! -d "$GENOME_DIR" ]; then
        echo "❌ Error: Genome directory not found: $GENOME_DIR"
        return 1
    fi
    
    # Check if GTF file exists
    if [ ! -f "$GTF_FILE" ]; then
        echo "❌ Error: GTF file not found: $GTF_FILE"
        return 1
    fi
    
    # Read each sample ID from the file
    while read -r sample_id; do
        # Find the IRGC_* directory containing the sample
        SAMPLE_DIR=$(find "$BASE_DIR" -type d -name "IRGC_*" -exec test -e {}/"${sample_id}_1.fastq.gz" \; -print | head -n 1)
        
        # Check if the sample directory was found
        if [[ -n "$SAMPLE_DIR" ]]; then
            # Define file paths for paired-end reads
            READ1="$SAMPLE_DIR/${sample_id}_1.fastq.gz"
            READ2="$SAMPLE_DIR/${sample_id}_2.fastq.gz"
            
            # Define output directory for STAR results
            OUTPUT_DIR="$OUTPUT_BASE_DIR/${sample_id}/"
            mkdir -p "$OUTPUT_DIR"
            
            # Run STAR alignment
            echo "  🚀 Processing sample: $sample_id (Found in $SAMPLE_DIR)"
            STAR --runThreadN "$THREADS" \
                 --genomeDir "$GENOME_DIR" \
                 --sjdbGTFfile "$GTF_FILE" \
                 --readFilesIn "$READ1" "$READ2" \
                 --readFilesCommand zcat \
                 --outFileNamePrefix "${OUTPUT_DIR}/${sample_id}_" \
                 --outSAMtype BAM SortedByCoordinate \
                 --limitBAMsortRAM 3000000000 \
                 --quantMode GeneCounts
            
            echo "  ✅ Alignment completed for: $sample_id (Results saved in $OUTPUT_DIR)"
        else
            echo "  ❌ ERROR: FASTQ files not found for sample: $sample_id"
        fi
    done < "${BASE_DIR}/sample_ids.txt"
    
    echo "✅ Read mapping complete."
}

# Step 6: Prepare the Reference Genome
prepare_reference() {
    echo "🧪 Preparing reference genome..."
    
    if [ -f "$GENOME_FASTA" ]; then
        samtools faidx "$GENOME_FASTA"
        echo "✅ Reference genome indexed."
    else
        echo "❌ Error: Reference genome FASTA not found: $GENOME_FASTA"
        return 1
    fi
}

# Step 7: Extract Gene Coordinates
extract_gene_coordinates() {
    echo "🔎 Extracting rbcL gene coordinates..."
    
    if [ -f "$GFF_FILE" ]; then
        grep -w "rbcL" "$GFF_FILE"
        
        # Create a BED file for rbcL gene
        # Note: These coordinates should be verified and updated based on your reference genome
        echo -e "NC_001320.1\t54095\t55528" > "${BASE_DIR}/rbcL.bed"
        echo "✅ rbcL gene coordinates extracted to ${BASE_DIR}/rbcL.bed"
    else
        echo "❌ Error: GFF file not found: $GFF_FILE"
        return 1
    fi
}

# Step 8: Extract the rbcL Gene Sequence
extract_gene_sequence() {
    echo "📝 Extracting rbcL gene sequence..."
    
    if [ -f "${BASE_DIR}/rbcL.bed" ] && [ -f "$GENOME_FASTA" ]; then
        bedtools getfasta -fi "$GENOME_FASTA" -bed "${BASE_DIR}/rbcL.bed" -fo "${BASE_DIR}/rbcL.fasta"
        echo "✅ rbcL gene sequence extracted to ${BASE_DIR}/rbcL.fasta"
    else
        echo "❌ Error: Missing required files for gene extraction."
        return 1
    fi
}

# Main execution
main() {
    echo "======================================"
    echo "      Gene Extraction Pipeline        "
    echo "======================================"
    echo "Starting pipeline at $(date)"
    echo ""
    
    # Check if required tools are installed
    check_dependencies
    
    # Uncomment the function calls you want to run
    # download_sra_data
    # convert_to_fastq
    # check_read_quality
    # generate_sample_ids
    # map_reads
    # prepare_reference
    # extract_gene_coordinates
    # extract_gene_sequence
    
    echo ""
    echo "Pipeline completed at $(date)"
    echo "======================================"
}

# Run the main function
main
```

## Usage

1. Edit the configuration variables at the top of the script to match your system and data paths.
2. Uncomment the function calls in the main() function that you want to run.
3. Run the script:

```bash
./gene_extraction.sh
```

## Step-by-Step Execution

You can run individual steps of the pipeline by uncommenting the corresponding function calls in the main() function:

1. **download_sra_data**: Downloads SRA data for the specified accessions
2. **convert_to_fastq**: Converts SRA files to FASTQ format and splits paired-end reads
3. **check_read_quality**: Performs quality control using FastQC
4. **generate_sample_ids**: Creates a list of sample IDs for mapping
5. **map_reads**: Maps reads to the reference genome using STAR
6. **prepare_reference**: Indexes the reference genome using Samtools
7. **extract_gene_coordinates**: Extracts coordinates for the rbcL gene
8. **extract_gene_sequence**: Extracts the rbcL gene sequence from the reference genome

## Important Notes

- Verify the rbcL gene coordinates in the `extract_gene_coordinates` function to ensure they match your reference genome.
- The pipeline is designed for paired-end reads. If you have single-end reads, modify the mapping commands accordingly.
- Adjust the `THREADS` variable to match your system's capabilities.

## Troubleshooting

- If SRA downloads fail, check your internet connection and ensure the SRA IDs are correct.
- If mapping fails, verify that the reference genome and GTF file are properly formatted.
- If gene extraction fails, check that the gene coordinates in the BED file are correct.



## Author

[Hariprasad T]
