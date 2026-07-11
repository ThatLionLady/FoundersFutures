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

total=0

# Output header
echo -e "Simulation\tRun\tCycle\tN_samples" > slim_output_sample_count.txt

# Loop through all ROH files
# Loop through specified directories
for SIM_DIR in Output/noSupp Output/w1Supp Output/duoSupp Output/Reintro; do
    if [[ -d "${SIM_DIR}" ]]; then
        SIM_NAME=$(basename "${SIM_DIR}" | cut -c1 | tr '[:lower:]' '[:upper:]')

        echo "STARTING $(date): A1_${SIM_NAME}"
        
        # Loop through Run directories inside each sim directory
        for RUN_DIR in $SIM_DIR/Run*/; do
            RUN_NAME=$(basename "$RUN_DIR" | sed 's/[^0-9]//g')  # Extract only the numeric part

            echo "STARTING $(date): A1_${SIM_NAME}${RUN_NAME}"

            for CYCLE in $(seq 5 5 20); do
                VCF=${RUN_DIR}${CYCLE}_SLiM_output.T.m0.vcf # extract CYCLE from VCF file

                N=$(bcftools query -l $VCF | wc -l)

                # Output result
                echo -e "${SIM_NAME}\t${RUN_NAME}\t${CYCLE}\t${N}" >> slim_output_sample_count.txt

                total=$((total + N))
            done
        done
    fi
done

echo "Total samples across VCF output: $total"