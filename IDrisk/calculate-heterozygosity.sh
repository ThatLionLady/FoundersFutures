#!/bin/bash

#!/bin/sh
set -euo pipefail

########################################
# Setup and Logging
########################################

# Save log in directory command was launched from
LAUNCH_DIR=$(pwd)
SCRIPT=$(basename "$0")
SCRIPT=${SCRIPT%.*}  # Removes the file extension (whatever it is)
LOG="${LAUNCH_DIR}/$(date '+%Y%m%d_%H%M%S')_${SCRIPT}.log"

# Log Header
{ 
    echo "Script: $0"
    echo "Started: $(date)"
    echo "Host: $(hostname)"
    echo
} > ${LOG}

########################################
# Argument Check
########################################


# Check if at least one command-line argument is provided
if [ "$#" -lt 2 ]; then
    echo ""    
    echo "Usage: $0 {loop|parallel|solo sample} <INPUT VCF>"
    echo ""
    exit 1
fi

########################################
# Environment Check
########################################

REQUIRED_TOOLS=("vcftools")

for tool in "${REQUIRED_TOOLS[@]}"; do
    if ! command -v "$tool" &> /dev/null; then
        echo "[ERROR] Required tool '$tool' not found in PATH." | tee -a "$LOG"
        exit 1
    fi
done


########################################
# Functions
########################################

end_function(){

    # At the beginning of what is to be tracked:
    # Record the start time for end_function
    # START_TIME=$(date +%s)

    # Record the end time
    END_TIME=$(date +%s)

    # Calculate the total execution time in seconds
    EXECUTION_TIME=$((END_TIME - START_TIME))

    # Calculate hours, minutes, and remaining seconds
    HOURS=$((EXECUTION_TIME / 3600))
    MINUTES=$((EXECUTION_TIME % 3600 / 60))
    SECONDS=$((EXECUTION_TIME % 60))

    # Print the execution time
    echo "[$(date '+%F %T')]: $1 completed in ${HOURS}:${MINUTES}:${SECONDS}"
    echo "----------------------------------------"

}

START_TIME=$(date +%s)

########################################
# Redirect all output after header to log
########################################

# This sends all stdout + stderr to both the terminal and the log.
exec > >(tee -a "$LOG") 2>&1

########################################
# Error & Signal Traps
########################################

# Log error messages automatically
trap 'echo "[ERROR] $(date): Script failed at line $LINENO (exit code: $?)" ' ERR

# Clean exit on Ctrl+C or job termination
cleanup() {
    echo "[INFO] $(date): Received termination signal. Cleaning up..."
    # Example cleanup steps: rm -f temp_*
    end_function "Script (terminated early)"
    exit 1
}
trap cleanup SIGINT SIGTERM

##########################################

VCF=$2
export VCF

# Function to process each sample
process_sample() {
    SAMPLE=$1

    vcftools --gzvcf ${VCF} --indv ${SAMPLE} --exclude-bed IDrisk/BED/roh_${SAMPLE}.bed --recode --stdout | gzip -c > IDrisk/HET/noROH_indiv_${SAMPLE}.vcf.gz
    vcftools --gzvcf IDrisk/HET/noROH_indiv_${SAMPLE}.vcf.gz --het --out IDrisk/HET/noROH_indiv_${SAMPLE}

}

export -f process_sample

# Function to run the loop
loop_task() {
    # Read samples into an array
    readarray -t samples < <(bcftools query -l "${VCF}")

    for SAMPLE in "${samples[@]}"; do
        process_sample ${SAMPLE}
    done
}

# Function to run the samples in parallel
parallel_task(){
    # Read samples into an array
    readarray -t samples < <(bcftools query -l "${VCF}")
    parallel process_sample ::: "${samples[@]}"
}

# Function to run the solo task
solo_task() {
    SAMPLE=$1
    process_sample ${SAMPLE}
}

# Main script execution based on the input parameter
case "$1" in
  loop)
    loop_task
    ;;
  parallel)
    parallel_task
    ;;
  solo)
    if [ -z "$3" ]; then
        echo "Usage: $0 solo <sample> <INPUT VCF>"
        exit 1
    fi
    solo_task $3
    ;;
  *)
    echo "Usage: $0 {loop|parallel|solo sample} <INPUT VCF>"
    exit 1
    ;;
esac