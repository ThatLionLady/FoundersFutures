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

# Output file
OUTPUT="compiled_data.csv"

# Print header to output file
echo "Sim,Run,cycle,N(t-1),N(t),Ne_heterozygosity,heterozygosity,inbreeding_load" > "$OUTPUT"

# Loop through specified directories
for SIM_DIR in Output/noSupp Output/Reintro Output/w1Supp Output/duoSupp; do
    if [[ -d "$SIM_DIR" ]]; then
        SIM_NAME=$(basename "$SIM_DIR" | cut -c1 | tr '[:lower:]' '[:upper:]')
        
        # Loop through Run directories inside each sim directory
        for RUN_DIR in "$SIM_DIR"/Run*/; do
            RUN_NAME=$(basename "$RUN_DIR" | sed 's/[^0-9]//g')  # Extract only the numeric part

            # Set the correct data file
            DATA_FILE="${RUN_DIR}Ne_log.csv"

            if [[ -f "$DATA_FILE" ]]; then
                awk -v sim="$SIM_NAME" -v run="$RUN_NAME" 'NR>1 {print sim "," run "," $0}' "$DATA_FILE" >> "$OUTPUT"
            else
                echo "Warning: No Ne_log.csv file found in $RUN_DIR"
            fi
        done
    else
        echo "Warning: Directory $SIM_DIR not found"
    fi

done


echo "Data compilation complete: $OUTPUT"
