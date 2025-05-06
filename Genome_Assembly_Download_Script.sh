# Genome Assembly Download Script

This script automates the downloading of genomic data from NCBI for multiple accession IDs using the `datasets` command-line tool.

## Features

- Downloads multiple genomic data formats (GFF3, RNA, CDS, protein, genome, sequence reports)
- Handles errors gracefully
- Checks for prerequisites
- Skips previously downloaded accessions
- Shows progress information

## Prerequisites

- NCBI's `datasets` command-line tool must be installed
- Bash shell environment

## Usage

1. Create a file named `accession.txt` with one NCBI accession ID per line
2. Run the script: `bash assembly.sh`

## Installation

```bash
# Make the script executable
chmod +x assembly.sh
```

## The Script

```bash
#!/bin/bash

# Genome Assembly Download Script
# Downloads genomic data for a list of NCBI accessions

# Configuration
ACCESSIONS_FILE="accession.txt"
OUTPUT_DIR="genome_data"
LOG_FILE="download_log.txt"

# Function to check if prerequisites are met
check_prerequisites() {
    # Check if datasets command is available
    if ! command -v datasets &> /dev/null; then
        echo "Error: 'datasets' command not found."
        echo "Please install NCBI datasets tool: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/"
        exit 1
    fi
    
    # Check if accession file exists
    if [ ! -f "$ACCESSIONS_FILE" ]; then
        echo "Error: Accession file '$ACCESSIONS_FILE' not found."
        echo "Please create a file named '$ACCESSIONS_FILE' with one accession ID per line."
        exit 1
    fi
    
    # Create output directory if it doesn't exist
    mkdir -p "$OUTPUT_DIR"
}

# Function to download data for a single accession
download_accession() {
    local accession=$1
    local output_file="${OUTPUT_DIR}/${accession}.zip"
    
    # Skip if already downloaded
    if [ -f "$output_file" ]; then
        echo "Skipping $accession - already downloaded"
        return 0
    fi
    
    echo "Downloading data for $accession ..."
    
    # Attempt to download with error handling
    if datasets download genome accession "$accession" \
        --include gff3,rna,cds,protein,genome,seq-report \
        --filename "$output_file" 2>> "$LOG_FILE"; then
        echo "✓ Successfully downloaded $accession"
        return 0
    else
        echo "✗ Failed to download $accession. See $LOG_FILE for details."
        return 1
    fi
}

# Main execution
main() {
    echo "=== Genome Assembly Download Script ==="
    echo "Starting at $(date)"
    echo "-----------------------------------"
    
    # Initialize the log file
    echo "Download log - $(date)" > "$LOG_FILE"
    
    # Check prerequisites
    check_prerequisites
    
    # Count total accessions
    total_accessions=$(grep -c "[^[:space:]]" "$ACCESSIONS_FILE")
    echo "Found $total_accessions accession IDs to process"
    
    # Process each accession
    count=0
    success=0
    failed=0
    
    while read -r accession; do
        # Skip empty lines
        if [[ -z "$accession" ]]; then
            continue
        fi
        
        # Increment counter
        ((count++))
        
        # Show progress
        echo "[$count/$total_accessions] Processing $accession"
        
        # Download the accession
        if download_accession "$accession"; then
            ((success++))
        else
            ((failed++))
        fi
    done < "$ACCESSIONS_FILE"
    
    # Summary
    echo "-----------------------------------"
    echo "Download summary:"
    echo "- Total accessions: $total_accessions"
    echo "- Successfully downloaded: $success"
    echo "- Failed: $failed"
    echo "- Output directory: $OUTPUT_DIR"
    echo "- Log file: $LOG_FILE"
    echo "Completed at $(date)"
}

# Run the script
main
```

## How to Test

1. Create a file named `accession.txt` with a few accession IDs, for example:
   ```
   GCF_000001405.40
   GCF_000001635.27
   GCF_000002035.6
   ```

2. Run the script:
   ```bash
   ./assembly.sh
   ```

3. Check the `genome_data` directory for downloaded files

## Error Codes

- Exit code 1: Missing prerequisites (datasets tool or accession file)
