#!/bin/bash
set -euo pipefail

########################################
# Setup and Logging
########################################

# Save log in directory command was launched from
LAUNCH_DIR=$(pwd)
mkdir -p "${LAUNCH_DIR}/Logs/"
SCRIPT=$(basename "$0")
SCRIPT=${SCRIPT%.*}  # Removes the file extension (whatever it is)
TIME=$(date '+%Y%m%d_%H%M%S')
LOG="${LAUNCH_DIR}/Logs/${TIME}_$(hostname)_${SCRIPT}.log"

# Log Header
{ 
    echo "Script: $0"
    echo "Started: $(date)"
    echo "Host: $(hostname)"
    echo
} > ${LOG}

########################################
# Environment Check
########################################

REQUIRED_TOOLS=("slim" "loop-progress")

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

########################################
# Script
########################################

MODELS_DIR=ZDC/SLiM/Scripts/Models/
INPUT_DIR=ZDC/SLiM/Input/
FAI=References/LionRefs/LionRef.fasta.fai

CHR=A1

echo "Running SLiM on Chromosome ${CHR}"

LENGTH=$(awk -v chr="${CHR}" '$1 == chr {print $2}' "${FAI}")

OUT=${LAUNCH_DIR}/Output/
mkdir -p ${OUT}

# Check that required input files exist
TRUE_FOUNDERS="${INPUT_DIR}/${CHR}.TrueFounders.vcf"
REINTRODUCED="${INPUT_DIR}/${CHR}.Reintroduced.vcf"
SUPPLEMENT1="${INPUT_DIR}/${CHR}.SAB-supplement.vcf"
SUPPLEMENT2="${INPUT_DIR}/${CHR}.LMP-supplement.vcf"

# Reintro = Reintroduced (N = 18)
if [[ ! -f "${REINTRODUCED}" ]]; then
    echo "Missing file: ${REINTRODUCED}, skipping chromosome ${CHR}"
    continue
fi    

# noSupp = True Founders -> no Supplementation (N = 19)
if [[ ! -f "${TRUE_FOUNDERS}" ]]; then
    echo "Missing file: ${TRUE_FOUNDERS}, skipping chromosome ${CHR}"
    continue
fi

# wSupp (N = 19 + 2) = Supplementation (N = 2)
if [[ ! -f "${SUPPLEMENT1}" ]]; then
    echo "Missing file: ${SUPPLEMENT1}, skipping chromosome ${CHR}"
    continue
fi

# wSupp (N = 19 + 2 + 2) = Supplementation (N = 2 + 2)
if [[ ! -f "${SUPPLEMENT2}" ]]; then
    echo "Missing file: ${SUPPLEMENT2}, skipping chromosome ${CHR}"
    continue
fi


# Run Models
for MODEL in w1Supp duoSupp noSupp Reintro; do

    echo "${CHR} --- Running SLiM for 10 iterations of ${MODEL}Sim20_vcf.slim"
    mkdir -p "${OUT}/${MODEL}/"
    cd "${OUT}/${MODEL}/"

    for ITERATION in {0..9}; do

        if [[ -f "${OUT}/${MODEL}/Run${ITERATION}/20_slimout.txt" ]]; then
            echo "  - Run${ITERATION} already complete, skipping."
            continue
        fi

        RUNLOG=${MODEL}_A1_${ITERATION}.log

    	START_TIME=$(date +%s)
    	echo "  - Run${ITERATION} START: $(date)" 
    	slim -d I=${ITERATION} "${MODELS_DIR}${MODEL}Sim20_vcf.slim" > ${RUNLOG} 2>&1

        if [[ $? -ne 0 ]]; then
            echo "  - SLiM failed on ${CHR} Run${ITERATION} (wSupp)"
        fi

    	end_function "RUN${ITERATION}" 

    done

done

echo "--- ${MODEL}Sim20_vcf.slim Complete" 
echo 

echo "${CHR} COMPLETE"
echo

rm progress.tmp
echo "TIME TO PLOT! Good luck!"