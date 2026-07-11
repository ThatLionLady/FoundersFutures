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
export -f end_function

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

mkdir -p VCFtools/
mkdir -p HET/
mkdir -p Lists/

# Function to run the samples in parallel
parallel_task(){
    SAMPLE=$1
    CHR=$2
    SIM_NAME=$3
    RUN_NAME=$4
    CYCLE=$5
    VCF=$6

    BFILE=A1_${SIM_NAME}${RUN_NAME}_${CYCLE}_${SAMPLE}
    vcftools_LOG=Logs/vcftools_A1_${SIM_NAME}${RUN_NAME}_${CYCLE}.log
    > "${vcftools_LOG}"

    BED="BED/roh_${BFILE}.bed"
    ROH_VCF="VCFtools/noROH_indiv_${BFILE}.vcf.gz"
    OUT="HET/noROH_indiv_${BFILE}"

    START_TIME=$(date +%s)
    echo "$(date) START: ${SAMPLE} (${BFILE})" >> "${vcftools_LOG}"

	if [[ ! -f "${ROH_VCF}" ]]; then

        echo "INFO: ${ROH_VCF} does not exist"

		if [[ ! -f ${BED} ]]; then

	    	echo "$(date) WARNING: BED file missing for ${BFILE}, no sites excluded from VCF" >> "${vcftools_LOG}"
	    	echo "$(date) WARNING: ${SAMPLE} (${BFILE}) has no BED file, ie no ROH > 2Mb. VCF does not exclude sites."
	    	# Run vcftools, and catch failures
	    	vcftools --gzvcf ${VCF} --indv ${SAMPLE} --recode --stdout | gzip -c > ${ROH_VCF}
		
		else

			# Run vcftools, and catch failures
	    	vcftools --gzvcf ${VCF} --indv ${SAMPLE} --exclude-bed ${BED} --recode --stdout | gzip -c > ${ROH_VCF}

	    fi

	else
		echo "$(date) SKIPPING: ${SAMPLE} (${BFILE}) VCF already exists."
	fi

	if [[ -f ${ROH_VCF} ]]; then

        echo "$(date) INFO: ${ROH_VCF} exists"

        if [[ ! -f "${OUT}.het" ]]; then

            echo "$(date) INFO: Running VCFtools for ${OUT}"
            vcftools --gzvcf ${ROH_VCF} --het --out "${OUT}"

        else

            echo "$(date): INFO: ${OUT}.het exists"

        fi

        if [[ -f "${OUT}.het" ]]; then

            echo "$(date) SUCCESS: ${SAMPLE} (${BFILE})"
            echo "$(date) SUCCESS: ${SAMPLE} (${BFILE})" >> "${vcftools_LOG}"

        else
            
            echo "$(date) ERROR: Failed het stats for ${SAMPLE} (${BFILE})"
            echo "$(date) ERROR: Failed het stats for ${SAMPLE} (${BFILE})" >> "${vcftools_LOG}"
        
        fi

	else

        echo "$(date) ERROR: Failed to extract noROH VCF for ${SAMPLE} (${BFILE})"
    
    fi

    end_function ${BFILE}

}

export -f parallel_task

# Record environment (only if ~/.parallel/ignored_vars doesn't exist)
[[ -f ~/.parallel/ignored_vars ]] || parallel --record-env

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

            	if [[ -f "${VCF}" ]]; then
            		echo "STARTING $(date): A1_${SIM_NAME}${RUN_NAME}_${CYCLE} (${VCF})"

            		readarray -t samples < <(bcftools query -l "${VCF}")
            		printf "%s\n" "${samples[@]}" > "Lists/samples_A1_${SIM_NAME}${RUN_NAME}_${CYCLE}.txt"
            		parallel -j $(nproc) parallel_task ::: "${samples[@]}" ::: "A1" ::: "${SIM_NAME}" ::: "${RUN_NAME}" ::: "${CYCLE}" ::: "${VCF}"

				fi
			done
			echo "COMPLETE $(date): A1_${SIM_NAME}${RUN_NAME}"
			echo "---"
		done
		echo "COMPLETE $(date): A1_${SIM_NAME}"
		echo "----"
	fi
done
echo "COMPLETE $(date): A1"
echo "-----"
echo

echo "ALL DONE. Proceed to SLiM_nonROH-het-per-kb.R"


