- [Founders \& Futures Analysis Methods](#founders--futures-analysis-methods)
  - [Mapping](#mapping)
  - [Calling](#calling)
- [Baseline Population Assessment](#baseline-population-assessment)
  - [Nucleotide Diversity (PIXY)](#nucleotide-diversity-pixy)
  - [Runs of Homozygosity](#runs-of-homozygosity)
  - [SNP Heterozygosity](#snp-heterozygosity)
  - [Population Structure (VCF2PCACluster)](#population-structure-vcf2pcacluster)
  - [Private Alleles](#private-alleles)
  - [IDRISK](#idrisk)
- [Viability Assessment via Forward Simulations](#viability-assessment-via-forward-simulations)

# Founders & Futures Analysis Methods

## Mapping

- [PALEOMIX](https://paleomix.readthedocs.io/en/stable/) pipeline:
  - BWA-MEM v0.7.17 mapped reads to the lion reference genome.
  - Samtools (v1.18) performed fixmate, sorting, merging, and calmd
  - Picard MarkDuplicates (v2.27.4) marked PCR duplicates.

The final BAM was generated with `samtools view -bF 3328`, filtering out unmapped, secondary, supplementary, QC-fail, and duplicate reads.

## Calling

- **Per-sample variant calling**: 
  - All sites were called independently for each individual using bcftools mpileup and call, 
  - with genotype quality, read depth, and allele balance (`'FMT/AB:1=float(ssum(FMT/AD[*:1]) / FMT/DP)'`) annotations added. 
  - Low-confidence variants (QUAL < 50) were removed prior to downstream analyses.
- **Multi-sample VCF generation**: 
  - Individual VCFs were merged into a single multi-sample VCF and multiallelic variants were decomposed into biallelic records (norm).
- **Genotype filtering**: 
  - Genotypes with sequencing depth <2× or >40× were masked as missing, and variant summary statistics (e.g., depth and missingness) were recalculated.
- **Variant filtering**: 
  - Indels and sex chromosome variants were excluded, and variants were filtered by missingness and total read depth. 
  - Analysis-specific VCFs were generated, including datasets permitting up to 80% missing data (for PIXY) and complete datasets containing only autosomal SNPs with no missing genotypes for downstream population genetic analyses.


```
bcftools mpileup -f References/LionRef.fasta -a DP,AD | bcftools call -m -f GQ,GP -Oz -o KAR-F001-B02.AS.vcf.gz
```
> *Date=Fri Jul 25 14:04:17 2025*

```
bcftools fill-tags -Oz -o KAR-F001.ASwAB.vcf.gz -- KAR-F001.AS.vcf.gz -t 'FMT/AB:1=float(ssum(FMT/AD[*:1]) / FMT/DP)'
```
> *Date=Wed Dec 31 10:55:49 2025*


```
bcftools view -i QUAL>=50 -Oz -o KAR-F001.ASwAB.qual50.vcf.gz KAR-F001.ASwAB.vcf.gz
```
> *Date=Thu Feb  5 17:14:58 2026*

```
bcftools merge -o OMGts.qual50.m.vcf.gz -Oz *.ASwAB.qual50.vcf.gz
```
> *Date=Thu Feb  5 21:40:45 2026*

```
bcftools norm -m-any -f References/LionRefs/LionRef.fasta --threads 30 -Oz -o LASTfc.split.vcf.gz OMGts.qual50.m.r.vcf.gz
```
> *Date=Thu Feb 19 18:06:52 2026*

```
bcftools setGT --threads 30 -Oz -o LASTfc.split.setGT-2fdp40.vcf.gz -- LASTfc.split.vcf.gz -t q -n . -i 'FMT/DP!="." & (FMT/DP<2 | FMT/DP>40)'
```
> *Date=Thu Feb 26 11:52:16 2026*

```
bcftools view -a --threads 10 LASTfc.split.setGT-2fdp40.vcf.gz
```
> *Date=Thu Feb 26 15:47:48 2026*

```
bcftools fill-tags --threads 10 -- -t all
```
> *Date=Thu Feb 26 15:47:48 2026*

```
bcftools fill-tags --threads 10 -Oz -o LASTfc.split.setGT-2fdp40.fixed.vcf.gz -- -t DP:1=int(sum(FORMAT/DP))
```
> *Date=Thu Feb 26 15:47:48 2026*

```
bcftools view --threads 10 -V indels -e 'CHROM=="Y" || CHROM=="X"' LASTfc.split.setGT-2fdp40.fixed.vcf.gz
```
> *Date=Thu Feb 26 21:02:53 2026*

```
bcftools view -i F_MISSING==0 -Oz -o LASTfc.split.setGT-2fdp40.fixed.4pixy-m00idp800.vcf.gz LASTfc.split.setGT-2fdp40.fixed.4pixy-m80idp800.vcf.gz
```
> *Date=Mon Mar  9 16:46:32 2026*

```
bcftools view -v snps -Oz -o LASTfc.split.setGT-2fdp40.fixed.4pixy-m00idp800.VO.vcf.gz LASTfc.split.setGT-2fdp40.fixed.4pixy-m00idp800.vcf.gz
```
> *Date=Tue Mar 10 11:49:11 2026*

# Baseline Population Assessment

## Nucleotide Diversity ([PIXY](https://pixy.readthedocs.io/en/latest/))

See [.sh]().

## Runs of Homozygosity

```sh
plink --allow-extra-chr \
      --homozyg \
      --homozyg-density 150 \
      --homozyg-gap 2000 \
      --homozyg-kb 200 \
      --homozyg-snp 30 \
      --homozyg-window-het 3 \
      --homozyg-window-missing 20 \
      --homozyg-window-snp 30 \
      --homozyg-window-threshold 0.05 \
      --out ${DIR}/ROH/${NAME} \
      --bfile ${DIR}/${NAME}
```
*Modified parameters from https://doi.org/10.1186/s12864-018-4489-0*  

## SNP Heterozygosity

```sh
vcftools --gzvcf ${VO} --het --out ${DIR}/HET/vcftools_${NAME}
```

## Population Structure ([VCF2PCACluster](https://github.com/hewm2008/VCF2PCACluster))

```sh
VCF2PCACluster --InVCF ${VCF} --OutPut PCAs/${NAME} -InSampleGroup pop.info 
```

## Private Alleles

See [.R]().

## ID<sub>RISK</sub>

See [IDrisk folder]() for scripts. 

# Viability Assessment via Forward Simulations

Genetically explicit forward simulations were done using [SLiM v4](https://messerlab.org/slim/)

See [SLiM folder]() for the lion-specific models and scripts to predict population size and heterozygosity of the population over a span of 20 years (cycles) in four scenarios: (1) the South African reintroduced lions only, (2) the population as it was at reintroduction (including the resident Mozambican male), (3) with supplementation, adding the Sabie males in cycle three, and (4) supplementing with an additional two males, LMP, in cycle six from additional data generated for this secondary supplementation event.
