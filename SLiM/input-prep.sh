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

REQUIRED_TOOLS=("bcftools" "loop-progress")

for tool in "${REQUIRED_TOOLS[@]}"; do
    if ! command -v "$tool" &> /dev/null; then
        echo "[ERROR] Required tool '$tool' not found in PATH." | tee -a "$LOG"
        exit 1
    fi
done

# Specific version of bcftools

REQUIRED_BCFTOOLS_VERSION="1.14"

BCFTOOLS_VERSION=$(
    bcftools --version |
    head -n1 |
    awk '{print $2}'
)

if [ "$(printf '%s\n' "$REQUIRED_BCFTOOLS_VERSION" "$BCFTOOLS_VERSION" | sort -V | head -n1)" != "$REQUIRED_BCFTOOLS_VERSION" ]; then
    echo "[ERROR] bcftools >= ${REQUIRED_BCFTOOLS_VERSION} required, found ${BCFTOOLS_VERSION}"
    exit 1
fi

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

VCF=ZDC/LASTfc.split.setGT-2fdp40.fixed.4pixy-m25idp800.vcf.gz
OUT=ZDC/SLiM/Input
REF=References/LionRefs/LionRef.fasta

echo "##################################################"
echo "$(date): Make SLiM Input VCFs"
echo

if [[ ! -f ${OUT}/A1.all.vcf.gz ]]; then
	# Reduce LN21 VCF to chromosome A1
	bcftools view -r A1 ${VCF} -Oz -o ${OUT}/A1.all.vcf.gz
	tabix ${OUT}/A1.all.vcf.gz
fi

if [[ ! -f ${OUT}/A1.all.multi.vcf.gz ]]; then
	# Make L21 VCF multiallelic (single site)
	bcftools norm -m + --threads 10 -Oz -o ${OUT}/A1.all.multi.vcf.gz ${OUT}/A1.all.vcf.gz
	tabix ${OUT}/A1.all.vcf.gz
	#bcftools query -f '%CHROM\t%POS\n' ${OUT}/A1.all.multi.vcf.gz > ${OUT}/site_list.txt
	#awk '{print $1"\t"$2-1"\t"$2}' ${OUT}/site_list.txt > ${OUT}/site_list.bed
fi

echo "-------------------------------------------------"
echo "$(date): Make LMP VCFs"
echo

LMP1="ZDC/VCFs/LMP-M001-B03.AS.vcf.gz"
LMP2="ZDC/VCFs/LMP-M002-B03.AS.vcf.gz"

if [[ ! -f ${OUT}/A1.LMP-M001.filtered.vcf.gz ]] && [[ ! -f ${OUT}/A1.LMP-M002.filtered.vcf.gz ]]; then

	for LMP in $LMP1 $LMP2; do

		START_TIME=$(date +%s)

		echo "--> ${LMP}"
		NAME=$(basename $LMP .AS.vcf.gz)
		ID=${NAME:0:8}
		OLD=$(bcftools query -l ${LMP})
		echo "--> ID = $ID"	

		if [[ ! -f ${OUT}/A1.${ID}.full.vcf.gz ]]; then

			# Make A1 VCF
			bcftools view -r A1 --threads 10 ${LMP} -Oz -o ${OUT}/A1.${ID}.full.vcf.gz
			tabix ${OUT}/A1.${ID}.full.vcf.gz

		fi

		if [[ ! -f ${OUT}/A1.${ID}.reheader.vcf.gz ]]; then
			
			# Reheader VCF
			bcftools reheader -s <(echo "${ID}-B03 ${ID}") ${OUT}/A1.${ID}.full.vcf.gz -o ${OUT}/A1.${ID}.reheader.vcf.gz
			tabix ${OUT}/A1.${ID}.reheader.vcf.gz
			
			echo "# Check that missing sites are truely missing from the input VCFs (site-investigation.sh)"		

		fi

		if [[ ! -f ${OUT}/A1.${ID}.filtered.vcf.gz ]]; then

		    # BY VIEW INCLUDE
		    # Filter out indels and keep sites with DP between 2 and 40
		    bcftools view -V indels --threads 10 ${OUT}/A1.${ID}.reheader.vcf.gz -Oz -o ${OUT}/A1.${ID}.filtered.vcf.gz -i 'FORMAT/DP>=2 & FORMAT/DP<=40'
		    tabix ${OUT}/A1.${ID}.filtered.vcf.gz

		fi

		end_function "$ID VCF"

	done

fi

echo "-------------------------------------------------"
echo "$(date): Merge LMP VCFs to make Supplement2"
echo

if [[ -f ${OUT}/A1.LMP-M001.filtered.vcf.gz ]] && [[ -f ${OUT}/A1.LMP-M002.filtered.vcf.gz ]]; then

	if [[ ! -f ${OUT}/A1.complete.vcf.gz ]]; then

		START_TIME=$(date +%s)
		# Merge VCFs
		bcftools merge --threads 16 -Oz -o ${OUT}/A1.complete.vcf.gz \
			${OUT}/A1.all.multi.vcf.gz ${OUT}/A1.LMP-M001.filtered.vcf.gz ${OUT}/A1.LMP-M002.filtered.vcf.gz
		tabix ${OUT}/A1.complete.vcf.gz
		echo "A1.complete.vcf.gz = $(bcftools query -l ${OUT}/A1.complete.vcf.gz | wc -l)"
		end_function "Merge"

	fi

# FYI: 
# at every site there currently are 17-21 GT at every position from the original file.
# AFTER MERGING:
# where LMP adds sites there will be 21-22 GT missing.
# if both LMP GT are missing, there will be 17-21 still.
# if one LMP GT is missng, there will be 18-22 GT.
# if both LMP GT present, there will be 19-23 GT.

	if [[ ! -f ${OUT}/A1.complete.tagged.vcf.gz ]]; then	
		START_TIME=$(date +%s)
		# Fix REF/ALT and fill F_MISSING tags
		bcftools view -a --threads 5 ${OUT}/A1.complete.vcf.gz | \
		bcftools +fill-tags --threads 5 -Oz -o ${OUT}/A1.complete.tagged.vcf.gz -- -t F_MISSING 
		tabix ${OUT}/A1.complete.tagged.vcf.gz
		end_function "Tagging"
	fi

	if [[ ! -f ${OUT}/A1.complete.filtered.vcf.gz ]]; then	
		START_TIME=$(date +%s)
		# Filter for F_MISSING removing sites with <50% (sites introduced by LMPs)
		bcftools view --threads 10 -i "F_MISSING<=0.5" ${OUT}/A1.complete.tagged.vcf.gz -Oz -o ${OUT}/A1.complete.filtered.vcf.gz
		tabix ${OUT}/A1.complete.filtered.vcf.gz
		end_function "Filtering"
	fi

	if [[ ! -f ${OUT}/A1.complete.split.vcf.gz ]]; then

		START_TIME=$(date +%s)
		# Normalize (aka split) multiallelic sites
		bcftools norm -m-any -f ${REF} --threads 16 -Oz -o ${OUT}/A1.complete.split.vcf.gz ${OUT}/A1.complete.filtered.vcf.gz
		tabix ${OUT}/A1.complete.split.vcf.gz
		end_function "Normalize"

	fi

	if [[ ! -f ${OUT}/A1.complete.split.VO.vcf.gz ]]; then

		START_TIME=$(date +%s)
		# Normalize (aka split) multiallelic sites
		bcftools view -v snps --threads 16 -Oz -o ${OUT}/A1.complete.split.VO.vcf.gz ${OUT}/A1.complete.split.vcf.gz
		tabix ${OUT}/A1.complete.split.VO.vcf.gz
		end_function "Variants"

	fi


	if [[ ! -f ${OUT}/A1.complete.split.VO.set.vcf.gz ]]; then

		START_TIME=$(date +%s)
		# BY SETTING GENOTYPES
		bcftools +setGT --threads 10 ${OUT}/A1.complete.split.VO.vcf.gz -Oz -o ${OUT}/A1.complete.split.VO.set.vcf.gz -- -t q -n 0 -i 'GT="./."'
		tabix ${OUT}/A1.complete.split.VO.set.vcf.gz
		end_function "Clear Missing GTs"

	fi

	echo "-------------------------------------------------"
	echo "$(date): Make Individual Group VCFs"
	echo

	if [[ ! -f ${OUT}/A1.LMP-supplement.vcf ]]; then	
		START_TIME=$(date +%s)
		echo "--> Supplement 2 (i.e. Limpopo males)"
		bcftools view -s LMP-M001,LMP-M002 ${OUT}/A1.complete.split.VO.set.vcf.gz | \
		bcftools view -g ^miss -o ${OUT}/A1.LMP-supplement.vcf
		echo "A1.LMP-supplement.vcf = $(bcftools query -l ${OUT}/A1.LMP-supplement.vcf | wc -l)"
		end_function "LMP-supplement"
	fi

	if [[ ! -f ${OUT}/A1.SAB-supplement.vcf ]]; then
		START_TIME=$(date +%s)
		echo "--> Supplement 1 (i.e. Sabi males)"
		bcftools view -s SAB-M001,SAB-M002 ${OUT}/A1.complete.split.VO.set.vcf.gz | \
		bcftools view -g ^miss -o ${OUT}/A1.SAB-supplement.vcf
		echo "A1.SAB-supplement.vcf = $(bcftools query -l ${OUT}/A1.SAB-supplement.vcf | wc -l)"
		end_function "SAB-supplement"
	fi

	if [[ ! -f ${OUT}/A1.TrueFounders.vcf ]]; then
		START_TIME=$(date +%s)
		echo "--> # True founders (i.e. 18 reintroduced + MOZ male)"
		bcftools view -s KAR-F001,KAR-F002,KHA-F001,KHA-F002,MAK-F001,MAK-F002,MAK-F003,MAK-F004,MAK-F005,MAK-F006,TEM-F001,TEM-F002,TEM-F003,TEM-F004,TSW-F004,TSW-F005,MOZ-M001,TSW-M001,TSW-M002 ${OUT}/A1.complete.split.VO.set.vcf.gz | \
		bcftools view -g ^miss -o ${OUT}/A1.TrueFounders.vcf
		echo "A1.TrueFounders.vcf = $(bcftools query -l ${OUT}/A1.TrueFounders.vcf | wc -l)"
		end_function "TrueFounders"
	fi

	if [[ ! -f ${OUT}/A1.Reintroduced.vcf ]]; then
		START_TIME=$(date +%s)
		echo "--> Reintroduced (i.e. 18 South African lions)"
		bcftools view -s KAR-F001,KAR-F002,KHA-F001,KHA-F002,MAK-F001,MAK-F002,MAK-F003,MAK-F004,MAK-F005,MAK-F006,TEM-F001,TEM-F002,TEM-F003,TEM-F004,TSW-F004,TSW-F005,TSW-M001,TSW-M002 ${OUT}/A1.complete.split.VO.set.vcf.gz | \
		bcftools view -g ^miss -o ${OUT}/A1.Reintroduced.vcf
		echo "A1.Reintroduced.vcf = $(bcftools query -l ${OUT}/A1.Reintroduced.vcf | wc -l)"
		end_function "Reintroduced"
	fi

else

	echo "WARNING: LMP VCFs not found. A1.LMP-supplement.vcf"

fi

