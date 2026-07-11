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
echo -e "Chromosome\tSimulation\tRun\tCycle\tN_samples" > counts4IDrisk.txt

# Loop through all ROH files
for file in ROH/ROH_*.hom; do
    if [[ -f "$file" ]]; then
        # Get the filename without path
        filename=$(basename "$file")

        # Remove extension and split parts
        name=${filename%.hom}  # remove .hom
        # Example: ROH_A1_N0_10 → parts=(ROH A1 N0 10)
        IFS='_' read -r prefix chr sim_run cycle <<< "$name"

        # Extract sim and run (e.g., sim_run=N0 → sim=N, run=0)
        sim="${sim_run:0:1}"
        run="${sim_run:1}"

        # Count unique IDs
        n=$(tail -n +2 "$file" | awk '{print $2}' | sort -u | wc -l)

        # Output result
        echo -e "${chr}\t${sim}\t${run}\t${cycle}\t${n}" >> counts4IDrisk.txt

        total=$((total + n))
    fi
done

echo "Total FROH calculations across all files: $total"